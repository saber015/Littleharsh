! *************************************************************************************************************** !
!
! Last modified by C 12/05/2016
!   - Set up for SHS (fingers crossed)
!
! *************************************************************************************************************** !
!
! Contains:   start - start file - called from littleharsh : calls ygrid, def_k, proc_lims, getbounds
!                     - reads input file, sets up domain, allocates memory, shares info between procs
!             ygrid - Called from start
!                     - Creates (streched) y grid points (seperate for u/w and v points)
!                       - Different grid streching in channel and immersed boundaries
!                     - Calculates dtheta/dy transformation
!                     - Calculates dy2i for all points (2nd derrivative) 
!             def_k - Called from start
!                     - Defines the x and z wavenumbers and their squares
!                       - Also calculates modified wavenumber but not used...
!         getbounds - Called from start
!                     - Calls boundary file and sets up gemomerty for immersed boundaries
! proc_lims_columns - Called from start
!                     - Calculates column lists for each proc and jlims
!                     - Also calculates dk
!                     - Calculates SHS boundary conditions
!  proc_lims_planes - Called from start
!                     - Defiens plane limits to balance between procs
!                     - Defines jgal
!            getini - Called from littleharsh : 
!                     - Grabs the initial conditions from file
!                     - )Or sets up parabolic profile)
!        init_stats - Called from getini : 
!                     - Initialises statistics
!
!
! *************************************************************************************************************** !

subroutine start(myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!    START    !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'            ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer iband,iproc,inputInt(15),k
  integer i
  real(8) inputR(20)           ! Message passing array containing input parameters

  iproc = 0
  do while (pnodes<np)
    iproc  = iproc+1
    pnodes = 2**iproc
  end do

! The master proc is in charge of: 
!  - reading input file
!  - initialise N, Ngal
!  - initialise the y-grid
!  and send this data to the other procs.
  if (myid==0) then
    open(40,file='input.in',form='formatted')
    do i = 1,9
      read(40,10)             ! Input file header
    end do
    read(40,10) Re            ! Reynolds number based on the bulk velocity and the channel half-height 
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,10) alp           ! The streamwise periodicity Lx of the box is Lx=2*pi/alp
    read(40,10) bet           ! The spanwise periodicity Lz of the box is Lz=2*pi/bet
    Lx = 2d0*pi/alp
    Lz = 2d0*pi/bet
    Ly = 2d0

    ! Ngal is the physical space matrix
    ! N is the Fourier space matrix
    ! They contain the same information except that N is 2/3Ngal in the number of x- and z-grid points,
    !  this is done to avoid aliasing.
    !            0   ----nb of bands----   4
    !        x:[ 0 , -nb of grid points- ,-2 ]
    ! Ngal = z:[ 0 , -nb of grid points- , 0 ]
    !       yu:[ a ,   b   ,  c  ,   d   , 0 ]
    !	yv:[ a ,   b   ,  c  ,   d   , 0 ]
    !   a: index of the wall. a+1 is the first point in the flow
    !   b: index of the last element of the first  band. b+1 is the first element of the next band
    !   c: index of the last element of the second band. c+1 is the symmetric point of b
    !   d: index of the last element of the third  band. d+1 is the wall and the symmetric point of a 
    allocate(Ngal(4,0:nband+1))
    Ngal = 0
    ! grid points in the 'x' direction
    read(40,20) Ngal(1,botband)
    read(40,20) Ngal(1,midband)
    Ngal(1,topband) = Ngal(1,botband)
    Ngal(1,nband+1) = -2
    ! grid points in the 'z' direction
    read(40,20) Ngal(2,botband)
    read(40,20) Ngal(2,midband)
    Ngal(2,topband) = Ngal(2,botband)
    Ngal(2,nband+1) = 0
    ! grid points in the 'y' direction from 'y=-1' to 'y=1'
    read(40,20) Ngal(3,nband)
    Ngal(3,nband) = Ngal(3,nband)-2

    allocate(N(4,0:nband+1))
    allocate(Ny(3,0:nband+1))
    N = Ngal
    do iband = 1,nband
      N(1,iband) = 2*(Ngal(1,iband)/3)
      N(2,iband) = 2*(Ngal(2,iband)/3)
    end do
    read(40,10) dyq              ! dymax/dymin
    read(40,20) ppp              ! y polinomial exponent
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank

    ! Initial conditions
    read(40,20) flag_init     ! specifies initial conditions: 1 new, 2 continuing, 3 parabolic
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,30) fnameimb  ! initial conditions file path, if flag_init=1
    read(40,10)               ! - Blank
    read(40,30) dirin         ! directory of the initial conditions
    read(40,30) dirout        ! directory where output files are saved
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank

    ! Parameters of the numerical method I
    read(40,10) CFL           ! CFL
    read(40,20) nwrite        ! frecuency of output in terminal, meassured in time steps
    read(40,21) flag_ctpress  ! sets constant flow rate (0) or pressure gradient (1)
    read(40,10)               ! - Blank
    if (flag_ctpress==1) then
      read(40,10) mpgx      ! PRESCRIBED MEAN PRESSURE GRADIENT
    else if (flag_ctpress==0) then
      read(40,10)
    else
      write(*,*) 'Wrong presure flag'
      stop
    end if
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank

    ! Parameters of the numerical method II
    read(40,10) Kib       ! immersed boundaries forcing parameter
    read(40,10) maxerr    ! maximum error for convergence in pressure
    read(40,10) maxt      ! maximum simulation time 't'
    read(40,10) maxA      ! boundary condition stability precission limiter
    maxA = -maxA
    read(40,20) physlim_bot      ! Defines the last physical space plane (bottom band) for tridiag solve
    read(40,20) physlim_top      ! Defines the last physical space plane (top band) for tridiag solve
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank

    ! Surface geometry
    read(40,20) geometry_type ! 0 smooth, 1 hydrophobic, 2 roughness
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,10)               ! - Blank
    read(40,20) ntilex        ! number of tiles in the periodic box
                              ! subroutine boundary checks that 'mod(Ngal(1,1)/ntilex) = 0'
    read(40,20) ntilez        ! number of tiles in the periodic box
                              ! subroutine boundary checks that 'mod(Ngal(2,1)/ntilez) = 0'
    if (geometry_type /= 0) then
      if (ntilex <= 0 .OR. ntilez <= 0) then
        write(*,*) 'Tiles must be positive. Moron.'
        stop
      end if 
    end if 
    read(40,10) posth         ! riblet height to width ratio
    read(40,20) npeakx        ! number of points representing the blade tip
    read(40,20) npeakz        ! number of points representing the blade tip
    read(40,10) Lfracx        ! number of points representing the blade tip
    read(40,10) Lfracz        ! number of points representing the blade tip
    read(40,10) shift_stag    ! fraction of spacing that 2nd row is shifted 
    read(40,20) dsty          ! number of points representing canopy in y 
    read(40,20) dstx          ! number of points representing canopy in x
    read(40,20) dstz          ! number of points representing canopy in z
    
    if (geometry_type == 2) then
      if (npeakx <= 0 .OR. npeakz <= 0) then
        write(*,*) 'The number of points per peak must be positive'
        stop
      end if 
    end if 

    if (geometry_type == 3) then
      if (dstx <= 0 .OR. dstz <= 0 .OR. dsty <= 0) then
        write(*,*) 'The number of points per peak must be positive'
        stop
      end if 
    end if 
    ! Height of the fine meshed band over the riblet's tips. 
    ! h_ny(1) = 1.5*posth means that the lower band extends up to 1.5 times the height of the riblets.
    allocate(h_ny(botband))
    if (geometry_type == 0) then
      post_spacing = 0d0
      posth        =   0d0
      h_ny(1)      = .25d0
    elseif (geometry_type == 1) then
      post_spacing = Lz/ntilez
      h_ny(1)      = 1.5d0*post_spacing
h_ny(1)      = 0.75d0*post_spacing
print *, "h_ny(1)      = 0.75d0*post_spacing"
! h_ny(1)      = .25d0
! print *, "h_ny(1)      = .25d0"
if(Lfracx==0.and.Lfracz==0)then
h_ny(1)=0.25d0
endif

    elseif (geometry_type == 3) then
      post_spacing = Lz/ntilez
      !posth        = posth
      h_ny(1)      = 2.0d0*post_spacing

    else
      post_spacing = Lz/ntilez
      posth        = posth*post_spacing
      h_ny(1)      = 1.5d0*posth
    end if
    close(40)

    ! Creates the geometry in the y-direction
    call y_grid_canopy

    allocate(u11(0:nn+2))

    write(ext1,'(i4.4)') N(1,1)       ! Writes the number of point in the x-direction in the first band into ext1
    write(ext2,'(i4.4)') N(2,1)       ! Writes the number of point in the z-direction in the first band into ext2
    write(ext3,'(i4.4)') nn+2         ! Writes the number of point in the y-direction                   into ext3
    write(ext4,'(i5.5)') int(10d0*t)  ! Writes the time                                                 into ext4

    filout     = 'boundary_'
    boundfname = trim(dirout)//trim(filout)//ext1//'x'//ext2//'x'//ext3//'.dat'     ! Example: boundary_192x360x182.dat
                                                                                    !  index('string','t') returns 2.
                                                                                    !  index(filout,' ') returns
                                                                                    !  the value of the first empty cell.
    filout     = 'hist_'
    fnameima   = trim(dirout)//trim(filout)//ext1//'x'//ext2//'x'//ext3//'.dat'     ! Example:     hist_0192x0360x0182.dat

    !open(30,file=fnameima,form='unformatted')
    !write(30) t,Re,alp,bet,mpgx,nband
    !write(30) N

inquire(file=fnameima, exist=exist_file_hist)
if (exist_file_hist) then
open (30,file=fnameima,form='unformatted',status="old", position="append", action="write")
else
open (30,file=fnameima,form='unformatted',status="new", action="write")
write(30) t,Re,alp,bet,mpgx,nband
write(30) N
end if
    

    
    filout     = 'c_four_'
    fnameima   = trim(dirout)//trim(filout)//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat' ! Example:   c_four_0192x0360x0182_t01554.dat

! The master thread send to the others all the information they'll need
    inputR( 1) = Re
    inputR( 2) = alp
    inputR( 3) = bet
    inputR( 4) = Lx
    inputR( 5) = Ly
    inputR( 6) = Lz
    inputR( 7) = Kib
    inputR( 8) = dtheta
    inputR( 9) = dthetai
    inputR(10) = CFL
    inputR(11) = maxt
    inputR(12) = mpgx
    inputR(13) = ddthetavi
    inputR(14) = Lfracx
    inputR(15) = Lfracz
    inputR(16) = posth
    inputR(17) = shift_stag
    inputInt(1) = flag_init
    inputInt(2) = nwrite
    inputInt(3) = nn
    inputInt(4) = ntilex
    inputInt(5) = ntilez
    inputInt(6) = npeakx
    inputInt(7) = npeakz
    inputInt(8) = flag_ctpress
    inputInt(9) = geometry_type
    inputInt(10)= physlim_bot
    inputInt(11)= physlim_top
    inputInt(12)= dsty
    inputInt(13)= dstx
    inputInt(14)= dstz

    do iproc=1,np-1
      call MPI_SEND(inputR  ,                     20,MPI_REAL8  ,iproc,       iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(inputInt,                     15,MPI_INTEGER,iproc,  1000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(N       ,4*(nband+2)            ,MPI_INTEGER,iproc,  2000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(Ngal    ,4*(nband+2)            ,MPI_INTEGER,iproc,  3000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(Ny      ,3*(nband+2)            ,MPI_INTEGER,iproc,  4000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(yu      ,   N(4,nband)-N(4,0)+2 ,MPI_REAL8  ,iproc,  7000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(dthdyu  ,   N(4,nband)-N(4,0)+2 ,MPI_REAL8  ,iproc,  8000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(dyu2i   ,3*(N(4,nband)-N(4,0)+2),MPI_REAL8  ,iproc,  9000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(yv      ,   N(3,nband)-N(3,0)+2 ,MPI_REAL8  ,iproc, 10000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(dthdyv  ,   N(3,nband)-N(3,0)+2 ,MPI_REAL8  ,iproc, 11000+iproc,MPI_COMM_WORLD,ierr)
      call MPI_SEND(dyv2i   ,3*(N(3,nband)-N(3,0)+2),MPI_REAL8  ,iproc, 12000+iproc,MPI_COMM_WORLD,ierr)
    end do

  else

! All the procs, except the master, receive and store the data
! Every proc knows the physical parameters and the whole geometry describing the problem (including N and Ngal)
    call MPI_RECV(inputR  ,20,MPI_REAL8  ,0,     myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(inputInt,15,MPI_INTEGER,0,1000+myid,MPI_COMM_WORLD,status,ierr)
    Re        = inputR( 1)
    alp       = inputR( 2)
    bet       = inputR( 3)
    Lx        = inputR( 4)
    Ly        = inputR( 5)
    Lz        = inputR( 6)
    Kib       = inputR( 7)
    dtheta    = inputR( 8)
    dthetai   = inputR( 9)
    CFL       = inputR(10)
    maxt      = inputR(11)
    mpgx      = inputR(12)
    ddthetavi = inputR(13)
    Lfracx    = inputR(14)
    Lfracz    = inputR(15)
    posth     = inputR(16)
    shift_stag= inputR(17)
    flag_init     = inputInt(1)
    nwrite        = inputInt(2)
    nn            = inputInt(3)
    ntilex        = inputInt(4)
    ntilez        = inputInt(5)
    npeakx        = inputInt(6)
    npeakz        = inputInt(7)
    flag_ctpress  = inputInt(8)
    geometry_type = inputInt(9)
    physlim_bot = inputInt(10)
    physlim_top = inputInt(11)
    dsty        = inputInt(12)
    dstx        = inputInt(13)
    dstz        = inputInt(14)
 

    allocate(N   (4,0:nband+1))
    allocate(Ngal(4,0:nband+1))
    allocate(Ny  (3,0:nband+1))
    call MPI_RECV(N      ,4*(nband+2)            ,MPI_INTEGER,0, 2000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(Ngal   ,4*(nband+2)            ,MPI_INTEGER,0, 3000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(Ny     ,3*(nband+2)            ,MPI_INTEGER,0, 4000+myid,MPI_COMM_WORLD,status,ierr)
    allocate(yu    (  N(4,0):N(4,nband)+1))
    allocate(dthdyu(  N(4,0):N(4,nband)+1))
    allocate(dyu2i (3,N(4,0):N(4,nband)+1))
    allocate(yv    (  N(3,0):N(3,nband)+1))
    allocate(dthdyv(  N(3,0):N(3,nband)+1))
    allocate(dyv2i (3,N(3,0):N(3,nband)+1))
    call MPI_RECV(yu     ,   N(4,nband)-N(4,0)+2 ,MPI_REAL8  ,0, 7000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(dthdyu ,   N(4,nband)-N(4,0)+2 ,MPI_REAL8  ,0, 8000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(dyu2i  ,3*(N(4,nband)-N(4,0)+2),MPI_REAL8  ,0, 9000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(yv     ,   N(3,nband)-N(3,0)+2 ,MPI_REAL8  ,0,10000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(dthdyv ,   N(3,nband)-N(3,0)+2 ,MPI_REAL8  ,0,11000+myid,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(dyv2i  ,3*(N(3,nband)-N(3,0)+2),MPI_REAL8  ,0,12000+myid,MPI_COMM_WORLD,status,ierr)
  end if

! From here on, every proc do everything which follows.

  nstat = min(20,nwrite)

  ! Initialise the FFT
  allocate(buffR_x(nband),buffC_z(nband))
  allocate(buffRal_x(nband),buffCal_z(nband))
  do iband = 1,nband
    allocate(buffR_x  (iband)%b( 2*N   (1,iband)+18), &
&            buffC_z  (iband)%b(10*N   (2,iband)+19), &
&            buffRal_x(iband)%b( 2*Ngal(1,iband)+18), &
&            buffCal_z(iband)%b(10*Ngal(2,iband)+19))
  end do
  
  do iband = 1,nband
    call rfti(N   (1,iband),buffR_x  (iband)%b)
    call cfti(N   (2,iband),buffC_z  (iband)%b)
    call rfti(Ngal(1,iband),buffRal_x(iband)%b)
    call cfti(Ngal(2,iband),buffCal_z(iband)%b)
  end do

  ! Define wavenumbers k1F_x and k1F_z
  call def_k

  ! Runge-Kutta coefficients (Le&Moin, 1991)
  aRK(1) =   4d0/15d0    ! aRK = alpha
  aRK(2) =   1d0/15d0    ! aRK = alpha
  aRK(3) =   1d0/ 6d0    ! aRK = alpha
  bRK(1) =   4d0/15d0    ! bRK = beta
  bRK(2) =   1d0/15d0    ! bRK = beta 
  bRK(3) =   1d0/ 6d0    ! bRK = beta
  gRK(1) =   8d0/15d0    ! gRK = gamma
  gRK(2) =   5d0/12d0    ! gRK = gamma
  gRK(3) =   3d0/ 4d0    ! gRK = gamma
  cRK(1) =   0d0         ! cRK = zeta
  cRK(2) = -17d0/60d0    ! cRK = zeta
  cRK(3) = - 5d0/12d0    ! cRK = zeta
  dRK(1) =   8d0/15d0    ! dRK = gamma + zeta
  dRK(2) =   2d0/15d0    ! dRK = gamma + zeta
  dRK(3) =   1d0/ 3d0    ! dRK = gamma + zeta

  ! Runge-Kutta coefficients (Spalart, 1991)
  ! We need to check what is the right coefficient for the pressure
  !aRK(1) =  29d0/ 96d0    ! aRK = alpha
  !aRK(2) = - 3d0/ 40d0    ! aRK = alpha
  !aRK(3) =   1d0/  6d0    ! aRK = alpha
  !bRK(1) =  37d0/160d0    ! bRK = beta
  !bRK(2) =   5d0/ 24d0    ! bRK = beta 
  !bRK(3) =   1d0/  6d0    ! bRK = beta
  !gRK(1) =   8d0/ 15d0    ! gRK = gamma
  !gRK(2) =   5d0/ 12d0    ! gRK = gamma
  !gRK(3) =   3d0/  4d0    ! gRK = gamma
  !cRK(1) =   0d0          ! cRK = zeta
  !cRK(2) = -17d0/ 60d0    ! cRK = zeta
  !cRK(3) = - 5d0/ 12d0    ! cRK = zeta
  !dRK(1) =   8d0/ 15d0    ! dRK = gamma + zeta
  !dRK(2) =   2d0/ 15d0    ! dRK = gamma + zeta
  !dRK(3) =   1d0/  3d0    ! dRK = gamma + zeta

  err = 10d0*maxerr
  !dtv = Re/(-k2F_x(N(1,1)/2)-k2F_z(N(2,1)/2)+4d0/(yu(N(4,0)+1)-yu(N(4,0)))**2) !E! yv or yu??? Smaller timestep?....
!print *, "dtv yv"
  dtv = Re/(-k2F_x(N(1,1)/2)-k2F_z(N(2,1)/2)+4d0/(yv(N(3,0)+1)-yv(N(3,0)))**2) !E! yv or yu??? Smaller timestep?....

10 FORMAT(7X,D10.1)
20 FORMAT(7X,I10)
21 FORMAT(12X,I5)
30 FORMAT(7X,A90)

  ! Calculates the different variables that define planes (phys) and pencils (spectral) for every proc 
  ! Computes planelim, bandPL, crossband and dk 
  ! Computes sband and eband
  !call proc_lims(myid)
  call proc_lims_columns(myid)
  call proc_lims_planes (myid)

  ! This function send to every proc a list containing points (and weights) of the immersed boundary.
  !  If the planes a proc has to solve are in the middle of the channel, the don't have to deal with the
  !  immersed boundaries, and then the list is empty.
  ! list_ib      stores the 'solid' points at the immersed boundary
  ! A_ib         stores the weights
  ! nyIB1, nyIB2 store  the planes that contain 'solid' point for a certain proc 
  call getbounds(myid,status,ierr)

  ! dnz = nb of points per tile
  if (geometry_type /= 0) then
    dnx = Ngal(1,1)/ntilex
    dnz = Ngal(2,1)/ntilez
!     dnx = N(1,1)/ntilex
!     dnz = N(2,1)/ntilez
!     print *, "dnx/z N"
  else
    dnx = 0
    dnz = 0
!     !Add Conditional for smooth wall?
!     dnx = Ngal(1,1)/ntilex
!     dnz = Ngal(2,1)/ntilez
  end if

  ! Creates a vector of indices with an extra tile and the peoriodicity.
  ! Ex: Ngal = 10, tilez = 3 -> indkor(0:11) = [10,1,2,...,10,1]
  ! It's only used in stats.f90 and only in bands 1 and 3 (phys)
  allocate(indkor(-dnz/2+1:Ngal(2,1)+dnz/2))
  do k = -dnz/2+1,0
    indkor(k) = k + Ngal(2,1)
  end do
  do k = 1,Ngal(2,1)
    indkor(k) = k
  end do
  do k = Ngal(2,1)+1,Ngal(2,1)+dnz/2
    indkor(k) = k - Ngal(2,1)
  end do
  allocate(indior(-dnx/2+1:Ngal(1,1)+dnx/2))
  do i = -dnx/2+1,0
    indior(i) = i + Ngal(1,1)
  end do
  do i = 1,Ngal(1,1)
    indior(i) = i
  end do
  do i = Ngal(1,1)+1,Ngal(1,1)+dnx/2
    indior(i) = i - Ngal(1,1)
  end do

  ! Grid weighting for extrapolating and interpolating the ghost points at the wall
  allocate(gridweighting(nband,2))
  gridweighting(botband,1) =-(yv(N(3,0)  )-yu(N(4,0)  ))/(yv(N(3,0)  )-yu(N(4,0)+1)) !band 1 bottom
  gridweighting(botband,2) = (yu(N(4,1)+1)-yv(N(3,1)+1))/(yv(N(3,1)+1)-yu(N(4,1)  )) !band 1 top
  gridweighting(midband,1) =-(yv(N(3,0)  )-yu(N(4,0)  ))/(yv(N(3,0)  )-yu(N(4,0)+1)) !band 2 bottom
  gridweighting(midband,2) = (yu(N(4,3)+1)-yv(N(3,3)+1))/(yv(N(3,3)+1)-yu(N(4,3)  )) !band 2 top
  gridweighting(topband,1) =-(yv(N(3,2)  )-yu(N(4,2)  ))/(yv(N(3,2)  )-yu(N(4,2)+1)) !band 3 bottom
  gridweighting(topband,2) = (yu(N(4,3)+1)-yv(N(3,3)+1))/(yv(N(3,3)+1)-yu(N(4,3)  )) !band 3 top
  ! Used in interp_v and v_corr
  allocate(gridweighting_interp(nband,2))
  gridweighting_interp(botband,1) = (yu(N(4,0)  )-yu(N(4,0)+1))/(yv(N(3,0)  )-yu(N(4,0)+1))
  gridweighting_interp(midband,1) = (yu(N(4,0)  )-yu(N(4,0)+1))/(yv(N(3,0)  )-yu(N(4,0)+1))
  gridweighting_interp(topband,1) = (yu(N(4,2)  )-yu(N(4,2)+1))/(yv(N(3,2)  )-yu(N(4,2)+1))
  gridweighting_interp(botband,2) = (yu(N(4,1)+1)-yu(N(4,1)  ))/(yv(N(3,1)+1)-yu(N(4,1)  ))
  gridweighting_interp(midband,2) = (yu(N(4,3)+1)-yu(N(4,3)  ))/(yv(N(3,3)+1)-yu(N(4,3)  ))
  gridweighting_interp(topband,2) = (yu(N(4,3)+1)-yu(N(4,3)  ))/(yv(N(3,3)+1)-yu(N(4,3)  ))

  ! PL = PLane
  ! By default variables are stored in columns

!  allocate( div_cPL         (igal,kgal,jgal(pgrid,1)-1:jgal(pgrid,2)+1))
!  allocate( div_outPL       (igal,kgal,jgal(pgrid,1)-1:jgal(pgrid,2)+1))
!  allocate( u_outPL         (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
!  allocate( v_outPL         (igal,kgal,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
!  allocate( w_outPL         (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))


  allocate( u1PL       (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate( u2PL       (igal,kgal,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate( u3PL       (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate( u1PL_itp   (igal,kgal,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate( u2PL_itp   (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate( u3PL_itp   (igal,kgal,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate( ppPL       (igal,kgal,jgal(pgrid,1)-1:jgal(pgrid,2)+1))
  allocate(Nu1PL       (igal,kgal,jgal(ugrid,1)  :jgal(ugrid,2)  ))
  allocate(Nu2PL       (igal,kgal,jgal(vgrid,1)  :jgal(vgrid,2)  ))
  allocate(Nu3PL       (igal,kgal,jgal(ugrid,1)  :jgal(ugrid,2)  ))
  allocate(Nu1PL_dy    (igal,kgal,jgal(ugrid,1)  :jgal(ugrid,2)  ))
  allocate(Nu2PL_dy    (igal,kgal,jgal(vgrid,1)  :jgal(vgrid,2)  ))
  allocate(Nu3PL_dy    (igal,kgal,jgal(ugrid,1)  :jgal(ugrid,2)  ))
  allocate(uu_cPL      (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate(uv_fPL      (igal,kgal,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(uw_cPL      (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate(vv_cPL      (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate(vw_fPL      (igal,kgal,jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(ww_cPL      (igal,kgal,jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  
  allocate(du1dy_planes(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(du2dy_planes(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(du3dy_planes(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  
  allocate(du1dy_planes2(Ngal(1,bandPL(myid))+2,Ngal(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(du2dy_planes2(Ngal(1,bandPL(myid))+2,Ngal(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(du3dy_planes2(Ngal(1,bandPL(myid))+2,Ngal(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  
  !allocate(Qcrit(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  !allocate(Qcrit(N(1,2)+2,N(2,2),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  
  allocate(Qcrit(Ngal(1,2)+2,Ngal(2,2),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  
  allocate( u1PLN(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate( u2PLN(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate( u3PLN(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate( u1PL_itpN(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate( u2PL_itpN(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  allocate( u3PL_itpN(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate( ppPLN(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(pgrid,1)-1:jgal(pgrid,2)+1))
  
!  allocate(div_out(sband:eband))
!  allocate(u_out(sband:eband))
!  allocate(v_out(sband:eband))
!  allocate(w_out(sband:eband))
  
  allocate(u1_itp(sband:eband))
  allocate(u2_itp(sband:eband))
  allocate(u3_itp(sband:eband))
  allocate(Nu1_dy(sband:eband))
  allocate(Nu2_dy(sband:eband))
  allocate(Nu3_dy(sband:eband))
  allocate(uv_f  (sband:eband))
  allocate(vw_f  (sband:eband))
  allocate(vv_c  (sband:eband))
  allocate(DG    (sband:eband))
  allocate(du1dy_columns(sband:eband))
  allocate(du2dy_columns(sband:eband))
  allocate(du3dy_columns(sband:eband))
  
  do iband = sband,eband
!    allocate(div_out     (iband)%f   (  jlim(1,pgrid,iband)  :jlim(2,pgrid,iband)  ,columns_num(iband,myid)))
!    allocate( u_out      (iband)%f   (  jlim(1,ugrid,iband)  :jlim(2,ugrid,iband)  ,columns_num(iband,myid)))
!    allocate( v_out      (iband)%f   (  jlim(1,vgrid,iband)  :jlim(2,vgrid,iband)  ,columns_num(iband,myid)))
!    allocate( w_out      (iband)%f   (  jlim(1,ugrid,iband)  :jlim(2,ugrid,iband)  ,columns_num(iband,myid)))

    allocate( u1_itp      (iband)%f   (  jlim(1,vgrid,iband)  :jlim(2,vgrid,iband)  ,columns_num(iband,myid)))
    allocate( u2_itp      (iband)%f   (  jlim(1,ugrid,iband)  :jlim(2,ugrid,iband)  ,columns_num(iband,myid)))
    allocate( u3_itp      (iband)%f   (  jlim(1,vgrid,iband)  :jlim(2,vgrid,iband)  ,columns_num(iband,myid)))
    allocate(Nu1_dy       (iband)%f   (  jlim(1,ugrid,iband)+1:jlim(2,ugrid,iband)-1,columns_num(iband,myid)))
    allocate(Nu2_dy       (iband)%f   (  jlim(1,vgrid,iband)+1:jlim(2,vgrid,iband)-1,columns_num(iband,myid)))
    allocate(Nu3_dy       (iband)%f   (  jlim(1,ugrid,iband)+1:jlim(2,ugrid,iband)-1,columns_num(iband,myid)))
    allocate( uv_f        (iband)%f   (  jlim(1,vgrid,iband)  :jlim(2,vgrid,iband)  ,columns_num(iband,myid)))
    allocate( vv_c        (iband)%f   (  jlim(1,ugrid,iband)  :jlim(2,ugrid,iband)  ,columns_num(iband,myid)))
    allocate( vw_f        (iband)%f   (  jlim(1,vgrid,iband)  :jlim(2,vgrid,iband)  ,columns_num(iband,myid)))
    allocate( DG          (iband)%f_dg(3,jlim(1,pgrid,iband)  :jlim(2,pgrid,iband)  ,columns_num(iband,myid)))
    allocate(du1dy_columns(iband)%f   (  jlim(1,vgrid,iband)  :jlim(2,vgrid,iband)  ,columns_num(iband,myid)))
    allocate(du2dy_columns(iband)%f   (  jlim(1,vgrid,iband)  :jlim(2,vgrid,iband)  ,columns_num(iband,myid)))
    allocate(du3dy_columns(iband)%f   (  jlim(1,vgrid,iband)  :jlim(2,vgrid,iband)  ,columns_num(iband,myid)))
  end do

  do iband = sband,eband
!   div_out(iband)%f       = 0d0
!    u_out(iband)%f        = 0d0
!    v_out(iband)%f        = 0d0
!    w_out(iband)%f        = 0d0
 
    u1_itp(iband)%f        = 0d0
    u2_itp(iband)%f        = 0d0
    u3_itp(iband)%f        = 0d0
    Nu1_dy(iband)%f        = 0d0
    Nu2_dy(iband)%f        = 0d0
    Nu3_dy(iband)%f        = 0d0
    uv_f  (iband)%f        = 0d0
    vw_f  (iband)%f        = 0d0
    vv_c  (iband)%f        = 0d0
    DG    (iband)%f_dg     = 0d0
    du1dy_columns(iband)%f = 0d0
    du3dy_columns(iband)%f = 0d0
  end do

  !C! Matrix fir the laplacian of the pressure
  !Build (and decompose) matirx - only needs to be done once
  

  do iband = sband,eband
    call LU_buildP(jlim(1,pgrid,iband),jlim(2,pgrid,iband),myid,iband,DG)
  enddo
  
  
   allocate(du1PL(igal,kgal,nyuIB1(myid):nyuIB2(myid)))
   allocate(du2PL(igal,kgal,nyvIB1(myid):nyvIB2(myid)))
   allocate(du3PL(igal,kgal,nyuIB1(myid):nyuIB2(myid)))

  allocate(   wx(igal,kgal,jgal(vgrid,1)-1  :jgal(vgrid,2)+1  ))
  !allocate(   wx(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(vgrid,1)-1  :jgal(vgrid,2)+1  ))

  iband = midband
  allocate(spU (jlim(1,ugrid,iband):jlim(2,ugrid,iband),columns_num(iband,myid)))
  allocate(spV (jlim(1,vgrid,iband):jlim(2,vgrid,iband),columns_num(iband,myid)))
  allocate(spW (jlim(1,ugrid,iband):jlim(2,ugrid,iband),columns_num(iband,myid)))
  allocate(spUV(jlim(1,ugrid,iband):jlim(2,ugrid,iband),columns_num(iband,myid)))
  allocate(spP (jlim(1,pgrid,iband):jlim(2,pgrid,iband),columns_num(iband,myid)))
  spU  = 0d0
  spV  = 0d0
  spW  = 0d0
  spUV = 0d0
  spP  = 0d0

end subroutine

subroutine ygrid
! TODO check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      YGRID     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Creates the geometry in the y-direction
! y(j) = a(j-q) + b(j-q)**p

! y is the physical coordinate
! theta is the mapping coordinate (constant increments: dtheta)
!  Since (1) dtheta is constant and (2) the mapping between y and theta is analytical,
!  the second order of the centered finite difference is preserved.

! n   is the total number of point from wall to wall
! q   is the index of the centerline (n/2) (y(q) = 0)
! dyq is dymax/dymin (max stretching)
! p   polinomial exponent
! a and b are obtained such that dy(0)/dy(q)=dyq and y(0)=-1

  use declaration
  implicit none
  integer j,iband,sumj,dny
  real(8) aaa,bbb,qqq,dy0,ddy,y0

  nn  = N(3,nband)+1 !v points (as collocated points)
  qqq = nn/2d0
  bbb = (dyq-1d0)/(qqq*(1d0-qqq**(ppp-1)+dyq*((1d0-qqq)*qqq**(ppp-1)-(1d0-qqq)**ppp)))
  aaa = (1d0-bbb*qqq**ppp)/qqq

  dtheta    = 1d0
  dthetai   = 1d0/dtheta
  !ddthetai=.5d0/dtheta
  ddthetavi = 1.0d0/dtheta

  dy0 = aaa*(1d0-qqq)+bbb*(1d0-qqq)**ppp+1d0

  if (geometry_type==2) then
    dny  = ceiling(posth/dy0)
    sumj = 0
    do j = -dny,-1
      sumj = sumj+j
    end do
    ! Lineal (almost constant) stretching in the immersed boundary region (below y=-1)
    ddy = (posth-dy0*dny)/sumj
  else
    dny = 0
    ddy = 0d0
  end if

  allocate(yu    (  -dny:nn+dny+1))
  allocate(dyu2i (3,-dny:nn+dny+1))
  allocate(dthdyu(  -dny:nn+dny+1))
  
  allocate(yv    (  -dny:nn+dny  ))
  allocate(dyv2i (3,-dny:nn+dny  ))
  allocate(dthdyv(  -dny:nn+dny  ))

!   !'Superhydrophobic' grid ------------------------------------------------------------------------
   yv(-dny)=-1d0-posth
    yu(-dny)=yv(-dny)-dy0*.5d0-ddy*(-dny-.5d0)*0.5d0
    do j=-dny-1,-dny,-1             ! Below the valleys
      yv(j)=yv(j+1)-dy0-ddy*(j  )
      yu(j)=yu(j+1)-dy0-ddy*(j-.5d0  )
    end do
    do j=-dny+1,-1                         ! Immersed boundary region
      yv(j)=yv(j-1)+dy0+ddy*(j-1)
      yu(j)=yu(j-1)+dy0+ddy*(j-1.5d0)
    end do
    
    if(dny>0)then !yu(0) point
    yu(0)=yu(-1)+dy0+ddy*(-1.5d0) !If immersed boundaries part of immersed boundary grid
    else
    yu(0)=-2-(aaa*(.5d0-qqq)+bbb*(.5d0-qqq)**ppp) !Else grid point reflected about -1
    endif
   do j=0,nn                              ! Flow
     yv(j)=aaa*(j-qqq)+bbb*(j-qqq)**ppp
     yu(j+1)=aaa*(j+.5d0-qqq)+bbb*(j+.5d0-qqq)**ppp  
   end do
    do j=nn+1,nn+dny                       ! Immersed boundary region and below valleys
      yv(j)=yv(j-1)+dy0-ddy*(j-nn)
      yu(j)=-yu(nn-j+1)
    end do
    yu(nn+dny+1)=-yu(-dny)
 

!  !C! Reflect immersed regions around y=-1/1
!  do j=0,-dny,-1 !E!2
!    yu(j)=-1-(yu(-j+1)+1) !E!2
!  enddo !E!2
!  do j=-1,-dny,-1 !E!2
!    yv(j)=-1-(yv(-j)+1) !E!2
!  enddo !E!2
!  
!  do j=nn+1,nn+dny        !E!2
!    yv(j)=-yv(nn-j) !E!2
!    yu(j)=-yu(nn-j+1) !E!2
!  enddo !E!2
!  yu(nn+dny+1)=-yu(-dny) !E!2
 
 
 
 
 ! Analytical expression for dtheta/dy in the lineal regions, and the polinomial one.
    do j=-dny,-1
      dthdyv(j)=dtheta/(dy0+ddy*(j-.5d0)) !E! yv(j) = dy0*j + ddy*0.5*(j+1)*j --> dy/dj = dy0 + ddy*(j+0.5)
      dthdyu(j)=dtheta/(dy0+ddy*(j-1d0)) !E! yuw(0) = yv(0) + 0.5*(dy0+ddy); yuw(1) = yv(1) + 0.5*(dy0+2*ddy) --> dy/dj = dy0 + ddy*(y+1)
    end do
    if(dny>0)then !dthdyu(0) point
    dthdyu(0)=dtheta/(dy0+ddy*(-1d0)) !E! If immersed boundaries part of immersed boundary grid
    else
   dthdyu(0)=1d0/(aaa+ppp*bbb*(.5d0-qqq)**(ppp-1)) !E! Else ghost point
    endif
 
   do j=0,nn
     dthdyv(j)=1d0/(aaa+ppp*bbb*(j-qqq)**(ppp-1))
     dthdyu(j+1)=1d0/(aaa+ppp*bbb*(j+.5d0-qqq)**(ppp-1))
   end do
    do j=nn+1,nn+dny
      dthdyv(j)=dtheta/(dy0-ddy*(j-nn+.5d0))
      dthdyu(j)=dthdyu(nn-j+1)
    end do
    dthdyu(nn+dny+1)=dthdyu(-dny)
    
   do j=0,-dny,-1 !E!2
   dthdyu(j)=dthdyu(-j+1) !E!2
   enddo !E!2
   do j = -1,-dny,-1 !E!2
   dthdyv(j)=dthdyv(-j) !E!2
   enddo !E!2
   
   do j=nn+1,nn+dny        !E!2
   dthdyu(j)=dthdyu(nn-j+1) !E!2
   dthdyv(j)=dthdyv(nn-j) !E!2
   enddo
   dthdyu(nn+dny       +1)=dthdyu(-dny)

 !   !'Wall' grid (just smooth channel) ------------------------------------------------------------------------
    
    yu(0)=aaa*(-1d0+.5d0-qqq)+bbb*(-1d0+.5d0-qqq)**ppp  
    do j=0,nn                              ! Flow
      yv(j)=aaa*(j-qqq)+bbb*(j-qqq)**ppp
      yu(j+1)=aaa*(j+.5d0-qqq)+bbb*(j+.5d0-qqq)**ppp  
    end do
  
  ! Analytical expression for dtheta/dy in the lineal regions, and the polinomial one.
    dthdyu(0)=1d0/(aaa+ppp*bbb*(-1d0+.5d0-qqq)**(ppp-1))
    do j=0,nn
      dthdyv(j)=1d0/(aaa+ppp*bbb*(j-qqq)**(ppp-1))
      dthdyu(j+1)=1d0/(aaa+ppp*bbb*(j+.5d0-qqq)**(ppp-1))
    end do
 
!    do j = -dny,nn+dny+1
!   print *, j, yu(j), yv(j), dthdyu(j), dthdyv(j)
!   enddo
!   stop
  
  
  
!    !'Original' grid ------------------------------------------------------------------------
!    
!   yv(-dny)=-1d0-posth
!   yu(-dny)=yv(-dny)-dy0*.5d0-ddy*(-dny-.5d0)*0.5d0
!   do j=-dny-1,-dny       ,-1             ! Below the valleys
!     yv(j)=yv(j+1)-dy0-ddy*(j  )
!     yu(j)=yu(j+1)-dy0-ddy*(j-.5d0  )
!   end do
!   do j=-dny+1,-1                         ! Immersed boundary region
!     yv(j)=yv(j-1)+dy0+ddy*(j-1)
!     yu(j)=yu(j-1)+dy0+ddy*(j-1.5d0)
!   end do
!   
!   if(dny>0)then !yu(0)
!   yu(0)=yu(-1)+dy0+ddy*(-1.5d0) !C! Part of immersed boundaries
!   else
!   yu(0)=-2-(aaa*(.5d0-qqq)+bbb*(.5d0-qqq)**ppp) !C! host point
!   endif
!   
!   do j=0,nn                              ! Flow
!     yv(j)=aaa*(j-qqq)+bbb*(j-qqq)**ppp
!     yu(j+1)=aaa*(j+.5d0-qqq)+bbb*(j+.5d0-qqq)**ppp  
!   end do
!   do j=nn+1,nn+dny                       ! Immersed boundary region and below valleys
!     yv(j)=yv(j-1)+dy0-ddy*(j-nn)
!     yu(j)=-yu(nn-j+1)
!   end do
!   yu(nn+dny+1)=-yu(-dny)
!   
!  
! ! Analytical expression for dtheta/dy in the lineal regions, and the polinomial one.
!   do j=-dny,-1
!     dthdyv(j)=dtheta/(dy0+ddy*(j-.5d0)) !E! yv(j) = dy0*j + ddy*0.5*(j+1)*j --> dy/dj = dy0 + ddy*(j+0.5)
!     dthdyu(j)=dtheta/(dy0+ddy*(j-1d0)) !E! yuw(0) = yv(0) + 0.5*(dy0+ddy); yuw(1) = yv(1) + 0.5*(dy0+2*ddy) --> dy/dj = dy0 + ddy*(y+1)
!   end do
!    if(dny>0)then !dthdyu(0) point
!    dthdyu(0)=dtheta/(dy0+ddy*(-1d0)) !E! If immersed boundaries part of immersed boundary grid
!    else
!    dthdyu(0)=1d0/(aaa+ppp*bbb*(.5d0-qqq)**(ppp-1)) !E! Else ghost point
!    endif
!   do j=0,nn
!     dthdyv(j)=1d0/(aaa+ppp*bbb*(j-qqq)**(ppp-1))
!     dthdyu(j+1)=1d0/(aaa+ppp*bbb*(j+.5d0-qqq)**(ppp-1))
!   end do
!   do j=nn+1,nn+dny
!     dthdyv(j)=dtheta/(dy0-ddy*(j-nn+.5d0))
!     dthdyu(j)=dthdyu(nn-j+1)
!   end do
!   dthdyu(nn+dny+1)=dthdyu(-dny) 
!   


  

! Second-order finite difference coefficient for 1. u(j-1), 2. u(j), 3. u(j+1)
  j=-dny
  dyu2i(1,j)= 0d0
  dyu2i(2,j)=-2d0/(yu(j+1)-yu(j))**2
  dyu2i(3,j)= 1d0/(yu(j+1)-yu(j))**2
  
  dyv2i(1,j)= 0d0
  dyv2i(2,j)=-2d0/(yv(j+1)-yv(j))**2
  dyv2i(3,j)= 1d0/(yv(j+1)-yv(j))**2
  do j=-dny+1,nn+dny-1
    dyu2i(1,j)=2d0/((yu(j)-yu(j-1))*(yu(j+1)-yu(j-1)))
    dyu2i(2,j)=-2d0/((yu(j+1)-yu(j))*(yu(j)-yu(j-1)))
    dyu2i(3,j)=2d0/((yu(j+1)-yu(j))*(yu(j+1)-yu(j-1)))
    
    dyv2i(1,j)=2d0/((yv(j)-yv(j-1))*(yv(j+1)-yv(j-1)))
    dyv2i(2,j)=-2d0/((yv(j+1)-yv(j))*(yv(j)-yv(j-1)))
    dyv2i(3,j)=2d0/((yv(j+1)-yv(j))*(yv(j+1)-yv(j-1)))
  end do
  j=nn+dny
  dyu2i(1,j)=2d0/((yu(j)-yu(j-1))*(yu(j+1)-yu(j-1)))
  dyu2i(2,j)=-2d0/((yu(j+1)-yu(j))*(yu(j)-yu(j-1)))
  dyu2i(3,j)=2d0/((yu(j+1)-yu(j))*(yu(j+1)-yu(j-1)))
  
  dyv2i(1,j)= 1d0/(yv(j)-yv(j-1))**2
  dyv2i(2,j)=-2d0/(yv(j)-yv(j-1))**2
  dyv2i(3,j)= 0d0
  
  j=nn+dny+1
  dyu2i(1,j)= 1d0/(yu(j)-yu(j-1))**2
  dyu2i(2,j)=-2d0/(yu(j)-yu(j-1))**2
  dyu2i(3,j)= 0d0
  
  nn=N(3,nband)

  N(3,0)=-dny
  N(4,0)=-dny
  do iband = 1,botband
    j = N(3,iband-1)
    h_ny(iband) = yu(0) + h_ny(iband)
    do while (yu(j)<h_ny(iband))
      j = j+1
    end do
    N(3,iband) = j
    N(3,nband-iband) = N(3,nband)-j
  end do
  N(3,nband)=nn+dny

  N(4,1) = N(3,1)+1 !Band 1-2
  N(4,2) = N(3,2) !Band 2-3 
  N(4,3) = N(3,3)+1 !Band 3-end 
  
  do iband=0,nband+1
    Ngal(3,iband)=N(3,iband)
    Ngal(4,iband)=N(4,iband)
  end do
  
  Ngal(4,1) = Ngal(3,1)+1 !Band 1-2
  Ngal(4,2) = Ngal(3,2) !Band 2-3 
  Ngal(4,3) = Ngal(3,3)+1 !Band 3-end 
    
  !C! Ny - temp variable, needs removing...
  !Ny(1,:) - v grid points
  !Ny(2,:) - u/w grid points
  !Ny(3,:) - p grid points
  Ny(ugrid,:  ) = N(4,:)
  Ny(vgrid,:  ) = N(3,:)
  Ny(pgrid,0  ) = N(4,0)+1
  Ny(pgrid,1  ) = N(4,1)-1
  Ny(pgrid,2  ) = N(4,2)+1
  Ny(pgrid,3  ) = N(4,3)-1
  
end subroutine

subroutine def_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    DEFINE K    !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define the wavenumbers k1 and their squares k2 = k1**2.
! They are defined for Fourier k1F, which gives an spectral accuracy,
!  which is in principle good. However, the highest frequencies are
!  not filtered and some wiggling may appear.
! In order to filter high frequencies a Differential wavenumber
!  is also defined, k1D, which reproduces the effective k of the
!  central difference scheme.

! ATTENTION: k1D is implemented in the code, but here, after being
!  defined, they are reasigned the same value as k1F.
!  Hence, when k1D is found in the code, it is really using k1F.

! Notes:
! alp = 2*pi/Lx (Lx being the streamwise periodicity of the box)
! bet = 2*pi/Lz (Lx being the spanwise   periodicity of the box)

  use declaration
  implicit none
  integer i,k,iband
  real(8) dzN   ,dzNi   ,dzN2i
  real(8) dxN   ,dxNi   ,dxN2i
!  real(8) dzNgal,dzNgali

  allocate(k1F_x  (0:N   (1,1)/2))
  allocate(k2F_x  (0:N   (1,1)/2))
  allocate(k1F_z  (1:N   (2,1)  ))
  allocate(k2F_z  (1:N   (2,1)  ))
 
  !!!!!!   define differential operator eigenvalues !!!!!!
  do i = 0,N(1,1)/2
    k1F_x(i) =  im* alp*i
    k2F_x(i) = -   (alp*i)**2
   
  end do
  do k = 1,N(2,1)/2
    k1F_z(k) =  im* bet*(k-1)
    k2F_z(k) = -   (bet*(k-1))**2

  end do
!    k1F_z(N(2,1)/2+1) = 0d0 !These are now included
!    k2F_z(N(2,1)/2+1) = 0d0 !These are now included
  do k = N(2,1)/2+1,N(2,1)
    k1F_z(k) = -im* bet*(N(2,1)-k+1)
    k2F_z(k) = -   (bet*(N(2,1)-k+1))**2
  end do
 
  
  
! For modified wavenumbers (2nd order centered difference)
 !  dzN   = Lz/N(2,1)
 !  dzNi  = 1d0/dzN
 !  dzN2i = dzNi**2
 !  dxN   = Lx/(N(1,1))
 !  dxNi  = 1d0/dxN
 !  dxN2i = dxNi**2
!! 
 !  do i = 0,N(1,1)/2
 !    !Modifed k
 !    k1F_x(i)= im*  sin(alp*(i)*dxN     )*dxNi
 !    k2F_x(i)= 2d0*(cos(alp*(i)*dxN)-1d0)*dxN2i
 !  end do
!!     
 !  do k = 1,N(2,1)/2
 !    !Modifed k
 !    k1F_z(k)= im*  sin(bet*(k-1)*dzN     )*dzNi
 !     k2F_z(k)= 2d0*(cos(bet*(k-1)*dzN)-1d0)*dzN2i
 !  end do
 !  do k = N(2,1)/2+1,N(2,1)
 !    !Modifed k
 !    k1F_z(k)=-im*  sin(bet*(N(2,1)-k+1)*dzN)     *dzNi
 !    k2F_z(k)= 2d0*(cos(bet*(N(2,1)-k+1)*dzN)-1d0)*dzN2i
 !  end do

  
  !print *, k1F_z(N(2,1))
  
end subroutine

subroutine getbounds(myid,status,ierr)
! TODO check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    GETBOUNDS   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This function send to every proc a list containing points (and weights) of the immersed boundary.
!  If the planes a proc has to solve are in the middle of the channel, the don't have to deal with the
!  immersed boundaries, and then the list is empty.


  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer j,iproc,ilist,ilist2
  integer nlist_ib_s_bot, nlist_ib_f_bot
  integer, allocatable:: inputNu(:), inputNv(:)                 ! Message passing array containing input parameters
  integer, allocatable:: list_ib_s_bot(:,:,:), list_ib_s_top(:,:,:)
  integer, allocatable:: list_ib_f_bot(:,:,:),list_ib_f_top(:,:,:)
  real(8), allocatable:: list_ib_f_w_bot(:,:,:),list_ib_f_w_top(:,:,:)
  integer, allocatable:: list_pointer_s(:,:), list_pointer_f(:,:)

  allocate(nlist_ib_s(2))
  allocate(nlist_ib_f(2))
  allocate(nyuIB1(0:np-1))
  allocate(nyuIB2(0:np-1))
  allocate(nyvIB1(0:np-1))
  allocate(nyvIB2(0:np-1))

  ! if nribs==0, smooth channel
  if (geometry_type/=3) then
    nlist_ib_s = 0
    nlist_ib_f = 0
    nyuIB1   = -88
    nyuIB2   = -99
    nyvIB1   = -88
    nyvIB2   = -99
  else
    allocate(inputNu(2*np+1+5))
    allocate(inputNv(2*np+1+5))

! MASTER

    if (myid==0) then
      ! Creates the geometry. 
      !  list_ib1 and list_ib2 store the points of the immersed boundaries for the lower and upper bands, resp.
      !  ny11, ny21, ny12 and ny22 stores the y-limits of the imm boundaries. From bottom to top.
      !  A_ib1 and A_ib2 stores the weights of the points in the imm boundaries.
      write(*,*) 'Creating geometry'
      call boundary_circ_rough
      write(*,*) 'Attention, nlist_ib(ugrid) must be equal to nlist_ib(vgrid)'
  
      open(10,file=boundfname,form='unformatted',access='stream')
      read(10) Lx,Ly,Lz
      read(10) Ngal,nlist_ib_s_bot,nlist_ib_f_bot,nyu11,nyu21,nyu12,nyu22,nyv11,nyv21,nyv12,nyv22
      allocate(list_ib_s_bot(  3,nlist_ib_s_bot,2))
      allocate(list_ib_f_bot(  9,nlist_ib_f_bot,2))
      allocate(list_ib_f_w_bot(2,nlist_ib_f_bot,2))
      allocate(list_ib_s_top(  3,nlist_ib_s_bot,2))
      allocate(list_ib_f_top(  9,nlist_ib_f_bot,2))
      allocate(list_ib_f_w_top(2,nlist_ib_f_bot,2))
      read(10) list_ib_s_bot, list_ib_s_top  
      read(10) list_ib_f_bot, list_ib_f_top
      read(10) list_ib_f_w_bot, list_ib_f_w_top 
      close(10)

      ! Tells every proc how many of their planes contain points of the immersed boundaries
      do iproc = 0,np-1
!TODO change these statements using pandPL or something similar
        if (iproc<np/2) then                                ! Lower riblets
          nyuIB1(iproc) = max(nyu11,planelim(ugrid,1,iproc))
          nyuIB2(iproc) = min(nyu21,planelim(ugrid,2,iproc))
          nyvIB1(iproc) = max(nyv11,planelim(vgrid,1,iproc))
          nyvIB2(iproc) = min(nyv21,planelim(vgrid,2,iproc))
        else                                                ! Upper riblets
          nyuIB1(iproc) = max(nyu12,planelim(ugrid,1,iproc))
          nyuIB2(iproc) = min(nyu22,planelim(ugrid,2,iproc))
          nyvIB1(iproc) = max(nyv12,planelim(vgrid,1,iproc))
          nyvIB2(iproc) = min(nyv22,planelim(vgrid,2,iproc))
        end if
      end do

      ! Prepares the information of the LOWER riblets to be sent to half of the procs
      do iproc = 1,np/2-1
        nlist_ib_s = 0
        nlist_ib_f = 0
        allocate(list_pointer_s(nlist_ib_s_bot,2))
        allocate(list_pointer_f(nlist_ib_f_bot,2))
        ! Checks the points of the imm boundaries that belong to planes of the proc 'iproc'
        ! If the point is in one of those planes, it is added to a list
        ! ugrid solid points bottom  
        do ilist = 1,nlist_ib_s_bot
          j = list_ib_s_bot(3,ilist,ugrid)                          ! 1:i 2:k 3:j
          if (j>=nyuIB1(iproc) .and. j<=nyuIB2(iproc)) then
            nlist_ib_s(ugrid) = nlist_ib_s(ugrid)+1
            list_pointer_s(nlist_ib_s(ugrid),ugrid) = ilist
          end if
        end do
        ! vgrid solid points bottom 
        do ilist = 1,nlist_ib_s_bot
          j = list_ib_s_bot(3,ilist,vgrid)                          ! 1:i 2:k 3:j
          if (j>=nyvIB1(iproc) .and. j<=nyvIB2(iproc)) then
            nlist_ib_s(vgrid) = nlist_ib_s(vgrid)+1
            list_pointer_s(nlist_ib_s(vgrid),vgrid) = ilist
          end if
        end do
        ! ugrid forcing points bottom  
        do ilist = 1,nlist_ib_f_bot
          j = list_ib_f_bot(3,ilist,ugrid)                          ! 1:i 2:k 3:j
          if (j>=nyuIB1(iproc) .and. j<=nyuIB2(iproc)) then
            nlist_ib_f(ugrid) = nlist_ib_f(ugrid)+1
            list_pointer_f(nlist_ib_f(ugrid),ugrid) = ilist
          end if
        end do
        ! vgrid forcing points bottom 
        do ilist = 1,nlist_ib_f_bot
          j = list_ib_f_bot(3,ilist,vgrid)                          ! 1:i 2:k 3:j
          if (j>=nyvIB1(iproc) .and. j<=nyvIB2(iproc)) then
            nlist_ib_f(vgrid) = nlist_ib_f(vgrid)+1
            list_pointer_f(nlist_ib_f(vgrid),vgrid) = ilist
          end if
        end do
        ! The list is reordered and stored in a local list
        allocate(s_list_ib(3,nlist_ib_s(ugrid),2))
        allocate(f_list_ib(9,nlist_ib_f(ugrid),2))
        allocate(w_list_ib(2,nlist_ib_f(ugrid),2))
        ! ugrid solid points bottom 
        do ilist2 = 1,nlist_ib_s(ugrid)
          ilist = list_pointer_s(ilist2,ugrid)
          s_list_ib(1,ilist2,ugrid) = list_ib_s_bot(1,ilist,ugrid) ! i
          s_list_ib(2,ilist2,ugrid) = list_ib_s_bot(2,ilist,ugrid) ! k
          s_list_ib(3,ilist2,ugrid) = list_ib_s_bot(3,ilist,ugrid) ! j
        end do
        ! vgrid solid points bottom 
        do ilist2 = 1,nlist_ib_s(vgrid)
          ilist = list_pointer_s(ilist2,vgrid)
          s_list_ib(1,ilist2,vgrid) = list_ib_s_bot(1,ilist,vgrid) ! i
          s_list_ib(2,ilist2,vgrid) = list_ib_s_bot(2,ilist,vgrid) ! k
          s_list_ib(3,ilist2,vgrid) = list_ib_s_bot(3,ilist,vgrid) ! j
        end do
        deallocate(list_pointer_s)
        ! ugrid forcing points bottom 
        do ilist2 = 1,nlist_ib_f(ugrid)
          ilist = list_pointer_f(ilist2,ugrid)
          f_list_ib(1,ilist2,ugrid) = list_ib_f_bot(1,ilist,ugrid) ! i
          f_list_ib(2,ilist2,ugrid) = list_ib_f_bot(2,ilist,ugrid) ! k
          f_list_ib(3,ilist2,ugrid) = list_ib_f_bot(3,ilist,ugrid) ! j
          f_list_ib(4,ilist2,ugrid) = list_ib_f_bot(4,ilist,ugrid) ! i2
          f_list_ib(5,ilist2,ugrid) = list_ib_f_bot(5,ilist,ugrid) ! k2
          f_list_ib(6,ilist2,ugrid) = list_ib_f_bot(6,ilist,ugrid) ! j2
          f_list_ib(7,ilist2,ugrid) = list_ib_f_bot(7,ilist,ugrid) ! i3
          f_list_ib(8,ilist2,ugrid) = list_ib_f_bot(8,ilist,ugrid) ! k3
          f_list_ib(9,ilist2,ugrid) = list_ib_f_bot(9,ilist,ugrid) ! j3

          w_list_ib(1,ilist2,ugrid) = list_ib_f_w_bot(1,ilist,ugrid) ! w2
          w_list_ib(2,ilist2,ugrid) = list_ib_f_w_bot(2,ilist,ugrid) ! w3
        end do
        ! vgrid forcing points bottom 
        do ilist2 = 1,nlist_ib_f(vgrid)
          ilist = list_pointer_f(ilist2,vgrid)
          f_list_ib(1,ilist2,vgrid) = list_ib_f_bot(1,ilist,vgrid) ! i
          f_list_ib(2,ilist2,vgrid) = list_ib_f_bot(2,ilist,vgrid) ! k
          f_list_ib(3,ilist2,vgrid) = list_ib_f_bot(3,ilist,vgrid) ! j
          f_list_ib(4,ilist2,vgrid) = list_ib_f_bot(4,ilist,vgrid) ! i2
          f_list_ib(5,ilist2,vgrid) = list_ib_f_bot(5,ilist,vgrid) ! k2
          f_list_ib(6,ilist2,vgrid) = list_ib_f_bot(6,ilist,vgrid) ! j2
          f_list_ib(7,ilist2,vgrid) = list_ib_f_bot(7,ilist,vgrid) ! i3
          f_list_ib(8,ilist2,vgrid) = list_ib_f_bot(8,ilist,vgrid) ! k3
          f_list_ib(9,ilist2,vgrid) = list_ib_f_bot(9,ilist,vgrid) ! j3

          w_list_ib(1,ilist2,vgrid) = list_ib_f_w_bot(1,ilist,vgrid) ! w2
          w_list_ib(2,ilist2,vgrid) = list_ib_f_w_bot(2,ilist,vgrid) ! w3
        end do
        deallocate(list_pointer_f)
        ! Send info to the corresponding procs
        ! The list with 'imm boundary points you have to handle with' is sent to every proc
        inputNu(1) = nlist_ib_s(ugrid)
        inputNv(1) = nlist_ib_s(vgrid)
        inputNu(2*np+1+5) = nlist_ib_f(ugrid)
        inputNv(2*np+1+5) = nlist_ib_f(vgrid)
        do j = 0,np-1
          inputNu(j+2   ) = nyuIB1(j)
          inputNu(j+2+np) = nyuIB2(j)
          inputNv(j+2   ) = nyvIB1(j)
          inputNv(j+2+np) = nyvIB2(j)
        end do
        inputNu(2*np+1+1) = nyu11
        inputNu(2*np+1+2) = nyu21
        inputNu(2*np+1+3) = nyu12
        inputNu(2*np+1+4) = nyu22
        inputNv(2*np+1+1) = nyv11
        inputNv(2*np+1+2) = nyv21
        inputNv(2*np+1+3) = nyv12
        inputNv(2*np+1+4) = nyv22

        call MPI_SEND(inputNu,2*np+1+5           ,MPI_INTEGER,iproc,20000+iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(inputNv,2*np+1+5           ,MPI_INTEGER,iproc,21000+iproc,MPI_COMM_WORLD,ierr)
        ! ATTENTION nlist_ib ugrid and vgrid are equal, otherwise it wont work
        call MPI_SEND(s_list_ib,3*nlist_ib_s(ugrid)*2,MPI_INTEGER,iproc,22000+iproc,MPI_COMM_WORLD,ierr) 
        call MPI_SEND(f_list_ib,9*nlist_ib_f(ugrid)*2,MPI_INTEGER,iproc,23000+iproc,MPI_COMM_WORLD,ierr) 
        call MPI_SEND(w_list_ib,2*nlist_ib_f(ugrid)*2,MPI_REAL8,iproc,24000+iproc,MPI_COMM_WORLD,ierr) 
        deallocate(s_list_ib, f_list_ib, w_list_ib)
      end do

      ! Prepares the information of the UPPER riblets to be sent to half of the procs
      do iproc = np/2,np-1
        nlist_ib_s = 0
        nlist_ib_f = 0
        allocate(list_pointer_s(nlist_ib_s_bot,2))
        allocate(list_pointer_f(nlist_ib_f_bot,2))
        ! Checks the points of the imm boundaries that belong to planes of the proc 'iproc'
        ! If the point is in one of those planes, it is added to a list
        ! ugrid solid points top
        do ilist = 1,nlist_ib_s_bot
          j = list_ib_s_top(3,ilist,ugrid)
          if (j>=nyuIB1(iproc) .and. j<=nyuIB2(iproc)) then
            nlist_ib_s(ugrid) = nlist_ib_s(ugrid)+1
            list_pointer_s(nlist_ib_s(ugrid),ugrid) = ilist
          end if
        end do
        ! vgrid solid points top
        do ilist = 1,nlist_ib_s_bot
          j = list_ib_s_top(3,ilist,vgrid)
          if (j>=nyvIB1(iproc) .and. j<=nyvIB2(iproc)) then
            nlist_ib_s(vgrid) = nlist_ib_s(vgrid)+1
            list_pointer_s(nlist_ib_s(vgrid),vgrid) = ilist
          end if
        end do
        ! ugrid forcing points top
        do ilist = 1,nlist_ib_f_bot
          j = list_ib_f_top(3,ilist,ugrid)
          if (j>=nyuIB1(iproc) .and. j<=nyuIB2(iproc)) then
            nlist_ib_f(ugrid) = nlist_ib_f(ugrid)+1
            list_pointer_f(nlist_ib_f(ugrid),ugrid) = ilist
          end if
        end do
        ! vgrid forcing points top
        do ilist = 1,nlist_ib_f_bot
          j = list_ib_f_top(3,ilist,vgrid)
          if (j>=nyvIB1(iproc) .and. j<=nyvIB2(iproc)) then
            nlist_ib_f(vgrid) = nlist_ib_f(vgrid)+1
            list_pointer_f(nlist_ib_f(vgrid),vgrid) = ilist
          end if
        end do
        ! The list is reordered and stored in a local list
        allocate(s_list_ib(3,nlist_ib_s(ugrid),2))
        allocate(f_list_ib(9,nlist_ib_f(ugrid),2))
        allocate(w_list_ib(2,nlist_ib_f(ugrid),2))
        ! ugrid solid points top
        do ilist2 = 1,nlist_ib_s(ugrid)
          ilist = list_pointer_s(ilist2,ugrid)
          s_list_ib(1,ilist2,ugrid) = list_ib_s_top(1,ilist,ugrid)  ! i
          s_list_ib(2,ilist2,ugrid) = list_ib_s_top(2,ilist,ugrid)  ! k
          s_list_ib(3,ilist2,ugrid) = list_ib_s_top(3,ilist,ugrid)  ! j
        end do
        ! vgrid solid points top
        do ilist2 = 1,nlist_ib_s(vgrid)
          ilist = list_pointer_s(ilist2,vgrid)
          s_list_ib(1,ilist2,vgrid) = list_ib_s_top(1,ilist,vgrid)  ! i
          s_list_ib(2,ilist2,vgrid) = list_ib_s_top(2,ilist,vgrid)  ! k
          s_list_ib(3,ilist2,vgrid) = list_ib_s_top(3,ilist,vgrid)  ! j
        end do
        deallocate(list_pointer_s)
        ! ugrid forcing points top
        do ilist2 = 1,nlist_ib_f(ugrid)
          ilist = list_pointer_f(ilist2,ugrid)
          f_list_ib(1,ilist2,ugrid) = list_ib_f_top(1,ilist,ugrid)  ! i
          f_list_ib(2,ilist2,ugrid) = list_ib_f_top(2,ilist,ugrid)  ! k
          f_list_ib(3,ilist2,ugrid) = list_ib_f_top(3,ilist,ugrid)  ! j
          f_list_ib(4,ilist2,ugrid) = list_ib_f_top(4,ilist,ugrid)  ! i2
          f_list_ib(5,ilist2,ugrid) = list_ib_f_top(5,ilist,ugrid)  ! k2
          f_list_ib(6,ilist2,ugrid) = list_ib_f_top(6,ilist,ugrid)  ! j2
          f_list_ib(7,ilist2,ugrid) = list_ib_f_top(7,ilist,ugrid)  ! i3
          f_list_ib(8,ilist2,ugrid) = list_ib_f_top(8,ilist,ugrid)  ! k3
          f_list_ib(9,ilist2,ugrid) = list_ib_f_top(9,ilist,ugrid)  ! j3

          w_list_ib(1,ilist2,ugrid) = list_ib_f_w_top(1,ilist,ugrid)  ! w2
          w_list_ib(2,ilist2,ugrid) = list_ib_f_w_top(2,ilist,ugrid)  ! w3
        end do
        ! vgrid forcing points top
        do ilist2 = 1,nlist_ib_f(vgrid)
          ilist = list_pointer_f(ilist2,vgrid)
          f_list_ib(1,ilist2,vgrid) = list_ib_f_top(1,ilist,vgrid)  ! i
          f_list_ib(2,ilist2,vgrid) = list_ib_f_top(2,ilist,vgrid)  ! k
          f_list_ib(3,ilist2,vgrid) = list_ib_f_top(3,ilist,vgrid)  ! j
          f_list_ib(4,ilist2,vgrid) = list_ib_f_top(4,ilist,vgrid)  ! i2
          f_list_ib(5,ilist2,vgrid) = list_ib_f_top(5,ilist,vgrid)  ! k2
          f_list_ib(6,ilist2,vgrid) = list_ib_f_top(6,ilist,vgrid)  ! j2
          f_list_ib(7,ilist2,vgrid) = list_ib_f_top(7,ilist,vgrid)  ! i3
          f_list_ib(8,ilist2,vgrid) = list_ib_f_top(8,ilist,vgrid)  ! k3
          f_list_ib(9,ilist2,vgrid) = list_ib_f_top(9,ilist,vgrid)  ! j3

          w_list_ib(1,ilist2,vgrid) = list_ib_f_w_top(1,ilist,vgrid)  ! w2
          w_list_ib(2,ilist2,vgrid) = list_ib_f_w_top(2,ilist,vgrid)  ! w3
        end do
        deallocate(list_pointer_f)
        ! The list with 'imm boundary points you have to handle with' is sent to every proc
        inputNu(1) = nlist_ib_s(ugrid)
        inputNv(1) = nlist_ib_s(vgrid)
        inputNu(2*np+1+5) = nlist_ib_f(ugrid)
        inputNv(2*np+1+5) = nlist_ib_f(vgrid)
        do j = 0,np-1
          inputNu(j+2   ) = nyuIB1(j)
          inputNu(j+2+np) = nyuIB2(j)
          inputNv(j+2   ) = nyvIB1(j)
          inputNv(j+2+np) = nyvIB2(j)
        end do
        inputNu(2*np+1+1) = nyu11
        inputNu(2*np+1+2) = nyu21
        inputNu(2*np+1+3) = nyu12
        inputNu(2*np+1+4) = nyu22
        inputNv(2*np+1+1) = nyv11
        inputNv(2*np+1+2) = nyv21
        inputNv(2*np+1+3) = nyv12
        inputNv(2*np+1+4) = nyv22
        call MPI_SEND(inputNu,2*np+1+5           ,MPI_INTEGER,iproc,20000+iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(inputNv,2*np+1+5           ,MPI_INTEGER,iproc,21000+iproc,MPI_COMM_WORLD,ierr)
        ! ATTENTION nlist_ib ugrid and vgrid are equal, otherwise it wont work
        call MPI_SEND(s_list_ib,3*nlist_ib_s(ugrid)*2,MPI_INTEGER,iproc,22000+iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(f_list_ib,9*nlist_ib_f(ugrid)*2,MPI_INTEGER,iproc,23000+iproc,MPI_COMM_WORLD,ierr)
        call MPI_SEND(w_list_ib,2*nlist_ib_f(ugrid)*2,MPI_REAL8  ,iproc,24000+iproc,MPI_COMM_WORLD,ierr)

!        write(*,*) 'myid', myid
!        write(*,*) 'iproc', iproc
!        write(*,*) 'nlist', nlist_ib_f(2)
!        write(*,*)  'j', f_list_ib(3,100,2)
!        write(*,*)  'w', w_list_ib(1,100,2)
!        write(*,*)
        deallocate(s_list_ib, f_list_ib, w_list_ib)
      end do

      ! Prepares the information to be sent to the master proc 
      iproc = 0
      nlist_ib_s = 0
      nlist_ib_f = 0
      allocate(list_pointer_s(nlist_ib_s_bot,2))
      allocate(list_pointer_f(nlist_ib_f_bot,2))
      ! ugrid solid points bottom 
      do ilist = 1,nlist_ib_s_bot
        j = list_ib_s_bot(3,ilist,ugrid)
        if (j>=nyuIB1(iproc) .and. j<=nyuIB2(iproc)) then
          nlist_ib_s(ugrid) = nlist_ib_s(ugrid)+1
          list_pointer_s(nlist_ib_s(ugrid),ugrid) = ilist
        end if
      end do
      ! vgrid solid points bottom   
      do ilist = 1,nlist_ib_s_bot
        j = list_ib_s_bot(3,ilist,vgrid)
        if (j>=nyvIB1(iproc) .and. j<=nyvIB2(iproc)) then
          nlist_ib_s(vgrid) = nlist_ib_s(vgrid)+1
          list_pointer_s(nlist_ib_s(vgrid),vgrid) = ilist
        end if
      end do
      ! ugrid forcing points bottom 
      do ilist = 1,nlist_ib_f_bot
        j = list_ib_f_bot(3,ilist,ugrid)
        if (j>=nyuIB1(iproc) .and. j<=nyuIB2(iproc)) then
          nlist_ib_f(ugrid) = nlist_ib_f(ugrid)+1
          list_pointer_f(nlist_ib_f(ugrid),ugrid) = ilist
        end if
      end do
      ! vgrid forcing points bottom   
      do ilist = 1,nlist_ib_f_bot
        j = list_ib_f_bot(3,ilist,vgrid)
        if (j>=nyvIB1(iproc) .and. j<=nyvIB2(iproc)) then
          nlist_ib_f(vgrid) = nlist_ib_f(vgrid)+1
          list_pointer_f(nlist_ib_f(vgrid),vgrid) = ilist
        end if
      end do

      allocate(s_list_ib(3,nlist_ib_s(ugrid),2))
      allocate(f_list_ib(9,nlist_ib_f(ugrid),2))
      allocate(w_list_ib(2,nlist_ib_f(ugrid),2))
      ! ugrid solid points bottom  
      do ilist2 = 1,nlist_ib_s(ugrid)
        ilist = list_pointer_s(ilist2,ugrid)
        s_list_ib(1,ilist2,ugrid) = list_ib_s_bot(1,ilist,ugrid)  ! i
        s_list_ib(2,ilist2,ugrid) = list_ib_s_bot(2,ilist,ugrid)  ! k
        s_list_ib(3,ilist2,ugrid) = list_ib_s_bot(3,ilist,ugrid)  ! j
      end do
      ! vgrid solid points bottom
      do ilist2 = 1,nlist_ib_s(vgrid)
        ilist = list_pointer_s(ilist2,vgrid)
        s_list_ib(1,ilist2,vgrid) = list_ib_s_bot(1,ilist,vgrid)  ! i
        s_list_ib(2,ilist2,vgrid) = list_ib_s_bot(2,ilist,vgrid)  ! k
        s_list_ib(3,ilist2,vgrid) = list_ib_s_bot(3,ilist,vgrid)  ! j
      end do
      ! ugrid forcing points bottom  
      do ilist2 = 1,nlist_ib_f(ugrid)
        ilist = list_pointer_f(ilist2,ugrid)
        f_list_ib(1,ilist2,ugrid) = list_ib_f_bot(1,ilist,ugrid)  ! i
        f_list_ib(2,ilist2,ugrid) = list_ib_f_bot(2,ilist,ugrid)  ! k
        f_list_ib(3,ilist2,ugrid) = list_ib_f_bot(3,ilist,ugrid)  ! j
        f_list_ib(4,ilist2,ugrid) = list_ib_f_bot(4,ilist,ugrid)  ! i2
        f_list_ib(5,ilist2,ugrid) = list_ib_f_bot(5,ilist,ugrid)  ! k2
        f_list_ib(6,ilist2,ugrid) = list_ib_f_bot(6,ilist,ugrid)  ! j2
        f_list_ib(7,ilist2,ugrid) = list_ib_f_bot(7,ilist,ugrid)  ! i3
        f_list_ib(8,ilist2,ugrid) = list_ib_f_bot(8,ilist,ugrid)  ! k3
        f_list_ib(9,ilist2,ugrid) = list_ib_f_bot(9,ilist,ugrid)  ! j3

        w_list_ib(1,ilist2,ugrid) = list_ib_f_w_bot(1,ilist,ugrid)  ! w2
        w_list_ib(2,ilist2,ugrid) = list_ib_f_w_bot(2,ilist,ugrid)  ! w3
      end do
      ! vgrid forcing points bottom
      do ilist2 = 1,nlist_ib_f(vgrid)
        ilist = list_pointer_f(ilist2,vgrid)
        f_list_ib(1,ilist2,vgrid) = list_ib_f_bot(1,ilist,vgrid)  ! i
        f_list_ib(2,ilist2,vgrid) = list_ib_f_bot(2,ilist,vgrid)  ! k
        f_list_ib(3,ilist2,vgrid) = list_ib_f_bot(3,ilist,vgrid)  ! j
        f_list_ib(4,ilist2,vgrid) = list_ib_f_bot(4,ilist,vgrid)  ! i2
        f_list_ib(5,ilist2,vgrid) = list_ib_f_bot(5,ilist,vgrid)  ! k2
        f_list_ib(6,ilist2,vgrid) = list_ib_f_bot(6,ilist,vgrid)  ! j2
        f_list_ib(7,ilist2,vgrid) = list_ib_f_bot(7,ilist,vgrid)  ! i3
        f_list_ib(8,ilist2,vgrid) = list_ib_f_bot(8,ilist,vgrid)  ! k3
        f_list_ib(9,ilist2,vgrid) = list_ib_f_bot(9,ilist,vgrid)  ! j3

        w_list_ib(1,ilist2,vgrid) = list_ib_f_w_bot(1,ilist,vgrid)  ! w2
        w_list_ib(2,ilist2,vgrid) = list_ib_f_w_bot(2,ilist,vgrid)  ! w3
      end do
      deallocate(list_pointer_s)
      deallocate(list_pointer_f)
      deallocate(list_ib_s_bot,list_ib_s_top)

! SLAVES

    else
      ! The information is received by every proc
      call MPI_RECV  (inputNu  ,2*np+1+5   ,MPI_INTEGER,0,20000+myid,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV  (inputNv  ,2*np+1+5   ,MPI_INTEGER,0,21000+myid,MPI_COMM_WORLD,status,ierr)
      nlist_ib_s(ugrid) = inputNu(1)
      nlist_ib_s(vgrid) = inputNv(1)
      nlist_ib_f(ugrid) = inputNu(2*np+1+5)
      nlist_ib_f(vgrid) = inputNv(2*np+1+5)
      do iproc = 0,np-1
        nyuIB1(iproc) = inputNu(iproc+2   )
        nyuIB2(iproc) = inputNu(iproc+2+np)
        nyvIB1(iproc) = inputNv(iproc+2   )
        nyvIB2(iproc) = inputNv(iproc+2+np)
      end do
      nyu11 = inputNu(2*np+1+1)
      nyu21 = inputNu(2*np+1+2)
      nyu12 = inputNu(2*np+1+3)
      nyu22 = inputNu(2*np+1+4)
      nyv11 = inputNv(2*np+1+1)
      nyv21 = inputNv(2*np+1+2)
      nyv12 = inputNv(2*np+1+3)
      nyv22 = inputNv(2*np+1+4)

      allocate(s_list_ib(3,nlist_ib_s(ugrid),2))
      allocate(f_list_ib(9,nlist_ib_f(ugrid),2))
      allocate(w_list_ib(2,nlist_ib_f(ugrid),2))
      ! ATTENTION nlist_ib ugrid and vgrid are equal, otherwise it wont work
      call MPI_RECV  (s_list_ib,3*nlist_ib_s(ugrid)*2,MPI_INTEGER,0,22000+myid,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV  (f_list_ib,9*nlist_ib_f(ugrid)*2,MPI_INTEGER,0,23000+myid,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV  (w_list_ib,2*nlist_ib_f(ugrid)*2,MPI_REAL8  ,0,24000+myid,MPI_COMM_WORLD,status,ierr)
!if (myid == 7) then
!  write(*,*) 'myid', myid
!  write(*,*) 'nlist', nlist_ib_f(2)
!  write(*,*)  'j', f_list_ib(3,100,2)
!  write(*,*)  'w', w_list_ib(1,100,2)
!end if
    end if
    deallocate(inputNu)
    deallocate(inputNv)
  end if
!  write(*,*) 'myid', myid
!  write(*,*) 'nlist', nlist_ib_f(2)
!  write(*,*)  'j', f_list_ib(3,100,2)
!  write(*,*)  'w', w_list_ib(1,100,2)
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  stop

end subroutine

subroutine proc_lims_columns(myid)
  ! Calculates the different variables that define columns (spectral) for every proc 
  ! Creates columns_num, columns_i, columns_k, and jlim
  ! Creating list of columns with x and z coordinates:
  !   columns_i(column_index,iband,myid), columns_k(column_index,iband,myid) and columns_num(iband,myid)
  ! The distribution for x and z is done in spectral space (Ring geometry)
  !  'band'  refers to sections of the ring 
  !    1. Bottom center (short pencils)     | 2 | 3 | 2 |
  !    2. Sides         (long  pencils)   y | 2 |   | 2 |
  !    3. Top    center (short pencils)     | 2 | 1 | 2 |
  !                                               z      

  use declaration
  implicit none

  include 'mpif.h'                         
  integer :: ierr,status

  integer :: iblock, iband, iproc, accum, column, i, k, myid
  integer :: columns_num_total, columns_num_long, columns_num_short
  integer :: max_columns_num
  integer :: columns_num_short_proc, columns_num_short_proc_rem
  integer :: columns_num_long_proc,  columns_num_long_proc_rem
  integer, allocatable :: columns_short_list_i(:), columns_short_list_k(:)
  integer, allocatable :: columns_long_list_i(:) , columns_long_list_k(:)
  
  
  integer :: postpointsx,postpointsz,texturepointsx,texturepointsz
  integer :: ip,kp
  
  ! Creating list of columns with x and z coordinates:
  !   columns_i(column_index,iband,myid), columns_k(column_index,iband,myid) and columns_num(iband,myid)

  ! NOTE 
  !   The columns corresponding to the last Fourier modes are also included.
  !   In the original code they are skipped in z. The * indicates places which need changes in order to skip some columns.

  ! Number of columns to be divided
  columns_num_total = (N(1,botband)/2 + 1) * (N(2,botband)) !*
  columns_num_long  = (N(1,midband)/2 + 1) * (N(2,midband)) !*
  columns_num_short = columns_num_total - columns_num_long    !*
  ! Distributing the columns among procs (and extra columns/remainder)
  columns_num_short_proc     = columns_num_short/np
  columns_num_short_proc_rem = columns_num_short - columns_num_short_proc*np
  columns_num_long_proc      = columns_num_long /np
  columns_num_long_proc_rem  = columns_num_long  - columns_num_long_proc *np
   
  ! Array storing the number of columns per proc
  ! Ex. columns_num(myid) contains the number/amount of columns of proc numb "myid"
  allocate(columns_num(nband,0:np-1))
  do iproc = 0,columns_num_short_proc_rem-1
    columns_num(botband,iproc) = columns_num_short_proc +1
    columns_num(topband,iproc) = columns_num_short_proc +1
  end do
  do iproc = 0,columns_num_long_proc_rem -1
    columns_num(midband,iproc) = columns_num_long_proc  +1
  end do
  do iproc = columns_num_short_proc_rem,np-1
    columns_num(botband,iproc) = columns_num_short_proc
    columns_num(topband,iproc) = columns_num_short_proc
  end do
  do iproc = columns_num_long_proc_rem, np-1
    columns_num(midband,iproc) = columns_num_long_proc
  end do
  
  ! Array with the coordinates of all points
  ! Short columns
  allocate(columns_short_list_i(columns_num_short))
  allocate(columns_short_list_k(columns_num_short))
  column = 0
  do   k = 1               , N(2,midband)/2    !do   k =  0               , N(2,midband)/2-1
    do i = N(1,midband)/2+1, N(1,botband)/2 
      column = column +1
      columns_short_list_i(column) = i
      columns_short_list_k(column) = k
    end do
  end do
  do   k = N(2,midband)/2+1, N(2,botband)/2    !do   k =  N(2,midband)/2  , N(2,botband)/2-1 
    do i = 0               , N(1,botband)/2 
      column = column +1
      columns_short_list_i(column) = i
      columns_short_list_k(column) = k
    end do
  end do
  do   k = N(2,botband)/2+1, N(2,botband) - N(2,midband)/2 + 0   !do   k = -N(2,botband)/2+1,-N(2,midband)/2
    do i = 0               , N(1,botband)/2 
      column = column +1
      columns_short_list_i(column) = i
      columns_short_list_k(column) = k
    end do
  end do
  
  do   k = N(2,botband) - N(2,midband)/2 + 1, N(2,botband)       !do   k = -N(2,midband)/2+1,-1
    do i = N(1,midband)/2+1, N(1,botband)/2 
      column = column +1
      columns_short_list_i(column) = i
      columns_short_list_k(column) = k
    end do
  end do

  
  ! Array with the coordinates of all points
  ! Long columns
  allocate(columns_long_list_i(columns_num_long))
  allocate(columns_long_list_k(columns_num_long))
  column = 0
  do   k =  1               ,N(2,midband)/2   !do   k =  0               ,N(2,midband)/2-1  
    do i =  0               ,N(1,midband)/2
      column = column +1
      columns_long_list_i(column) = i
      columns_long_list_k(column) = k
    end do
  end do
  do   k = N(2,botband) - N(2,midband)/2 + 1, N(2,botband)  !do   k = -N(2,midband)/2+1,-1
    do i =  0               ,N(1,midband)/2
      column = column +1
      columns_long_list_i(column) = i
      columns_long_list_k(column) = k
    end do
  end do
  
  ! Creating the lists with column coordinates
  ! columns_i(ncolumn,iband,myid) has the x coordiante of column number ncolumn in band iband.
  max_columns_num = maxval(columns_num)
  allocate(columns_i(max_columns_num,nband,0:np-1))
  allocate(columns_k(max_columns_num,nband,0:np-1))

  accum = 0
  do iproc = 0,np-1
    do column = 1,columns_num(botband,iproc)
      columns_i(column,botband,iproc) = columns_short_list_i(column + accum)
      columns_k(column,botband,iproc) = columns_short_list_k(column + accum)
      columns_i(column,topband,iproc) = columns_short_list_i(column + accum)
      columns_k(column,topband,iproc) = columns_short_list_k(column + accum)
    end do
    accum = accum + columns_num(botband,iproc)
  end do
  
  accum = 0
  do iproc = 0,np-1
    do column = 1,columns_num(midband,iproc)
      columns_i(column,midband,iproc) = columns_long_list_i(column + accum)
      columns_k(column,midband,iproc) = columns_long_list_k(column + accum)
    end do
    accum = accum + columns_num(midband,iproc)
  end do
  
  deallocate(columns_short_list_i)
  deallocate(columns_short_list_k)
  deallocate(columns_long_list_i)
  deallocate(columns_long_list_k)

  ! Building jlim. jlim stores the height of column in a certain band
  ! jlim(bottom/top,grid,band)
  allocate(jlim(2,3,nband)) ! jlim(bottom(1)/top(2),grid,iband) 
  ! u and w: Include ghosts points
  jlim(1,ugrid,1) = N(4,0)
  jlim(2,ugrid,1) = N(4,1)+1 
  jlim(1,ugrid,2) = N(4,0)
  jlim(2,ugrid,2) = N(4,3)+1
  jlim(1,ugrid,3) = N(4,2)
  jlim(2,ugrid,3) = N(4,3)+1
  ! v
  jlim(1,vgrid,1) = N(3,0)
  jlim(2,vgrid,1) = N(3,1)+1 
  jlim(1,vgrid,2) = N(3,0)
  jlim(2,vgrid,2) = N(3,3)+1
  jlim(1,vgrid,3) = N(3,2)
  jlim(2,vgrid,3) = N(3,3)+1
  ! p: Ghost points not included
  jlim(1,pgrid,1) = N(4,0)+1
  jlim(2,pgrid,1) = N(4,1)
  jlim(1,pgrid,2) = N(4,0)+1
  jlim(2,pgrid,2) = N(4,3)
  jlim(1,pgrid,3) = N(4,2)+1
  jlim(2,pgrid,3) = N(4,3)

  ! Used in FOU3D and stats
  ! It's a shift in z, used to align modes in different bands
  allocate(dk(max_columns_num,nband,0:np-1,sband:eband))
  dk = 0
  do iproc = 0,np-1
    do iblock = sband,eband
      do iband = sband,eband
        do column = 1,columns_num(iblock,iproc)
          if (columns_k(column,iblock,iproc) > N(2,1)/2) then
            dk(column,iblock,iproc,iband) = N(2,1)-Ngal(2,iband)
          end if
        end do
      end do
    end do
  end do
  
  allocate(dk_phys(max_columns_num,nband,0:np-1,sband:eband))
  dk_phys = 0
  do iproc = 0,np-1
    do iblock = sband,eband
      do iband = 2,2
        do column = 1,columns_num(iblock,iproc)
          if (columns_k(column,iblock,iproc) > N(2,1)/2) then
            dk_phys(column,iblock,iproc,iband) = N(2,1)-N(2,iband)
          end if
        end do
      end do
    end do
  end do
  
  !Define Boundary condition lists for SHS
  if(geometry_type==1)then
    !1 - No-slip
    !2 - Free-shear
  
    allocate(planeBC(0:(N(1,1)+2),N(2,1)))
  
    planeBC = 2 !2 - Free-shear

    texturepointsx = (N(1,1))/ntilex
    texturepointsz = N(2,1)/ntilez
    
    if(Lfracx==0.and.Lfracz==0)then
      postpointsx = texturepointsx
      postpointsz = texturepointsz
    else
  
      if(mod((N(1,1)/2d0/ntilex/Lfracx)+1,1d0)/=0)then
	print *, "Post x problem...", N(1,1)/2d0/ntilex/Lfracx+1,N(1,1)/2d0,ntilex,Lfracx
	stop
      endif
  
      if(mod((N(2,1)/ntilez/Lfracz)+1,1d0)/=0)then
	print *, "Post z problem...", N(2,1)/ntilez/Lfracz+1,N(1,1),ntilez,Lfracz
	stop
      endif
  
    postpointsx = (N(1,1)/ntilex/Lfracx)+1
    postpointsz = (N(2,1)/ntilez/Lfracz)+1
    
    endif
  
    if(postpointsz>N(2,1))then
      postpointsz=N(2,1) !overflow check...
    endif
  
    !Create geometry for 1 post
    do i=1,postpointsx
      do k=1,postpointsz
	planeBC(i,k) = 1 !1 - No-slip
      enddo
    enddo 
  
    !Duplicate for all posts
    do ip=1,ntilex
      do kp=1,ntilez
	do i=1,postpointsx
	  do k=1,postpointsz
	    planeBC(i+(ip-1)*texturepointsx,k+(kp-1)*texturepointsz) = planeBC(i,k)
	  enddo        
	enddo
      enddo
    enddo

  endif
  
!planeBC = 1 !2 - Free-shear
!bslip = (0.535d0*6d0)/(180d0)
  
end subroutine

subroutine proc_lims_planes(myid)
  ! Calculates the different variables that define planes (phys) for every proc 
  ! Defines i-k-,j-gal, planelim and bandPL
  ! First computes roughly the proportional load that each proc should handle.
  !  Then it tries to share all the planes, allocating more procs for the finer planes.
  !  Once balanced, it establishes the limit planes (planelim) for the procs.

  use declaration
  implicit none

  integer :: myid
  integer :: nplanes, iplanes
  integer :: iband, iproc, rem
  real(8) :: loadT, loadP, maxload
  real(8), allocatable:: load_band(:),i_load(:)
  
  integer :: proc_planes
  integer :: rem_planes

  allocate(load_band(nband),i_load(nband))
  allocate(procs  (0:nband))
  allocate(procs_b(0:nband))
  allocate(procs_c(0:nband))
  allocate(planelim  (3,2,0:np-1))
  allocate(limPL_incw(3,2,0:np-1))
  allocate(limPL_excw(3,2,0:np-1))
  allocate(limPL_FFT(3,2,0:np-1))
  allocate(bandPL(0:np-1))
  allocate(bandPL_FFT(0:np-1))
  
  limPL_FFT = 9999

  ! Computing mean load per proc
  loadT = 0d0
  do iband = 1,nband
    ! Roughly the load of a band is equal to the number of points
    load_band(iband) = (Ngal(1,iband)+2)*Ngal(2,iband)*(Ngal(4,iband)-N(4,iband-1)) ! Number of points (nx+2)*nz*ny(band)
    loadT            = loadT + load_band(iband)
  end do
  loadP = loadT/np ! Load per proc
  ! Assigns one proc per plane to the fine-grid bands
  ! Ussually we'll run out of procs. This will be corrected in a following step
  procs(0) = 0  ! Procs per band
  rem      = 0
  nplanes         = N(4,botband) - N(4,0)                     ! Nb of planes in a band (phys)
  procs(botband)  = nplanes                                   ! One proc per plane (too greedy)
  procs(topband)  = nplanes                                   !
  rem             = rem + procs(botband) + procs(topband)
  i_load(botband) = load_band(botband)/nplanes                ! Computational load per plane
  i_load(topband) = load_band(topband)/nplanes                ! Computational load per plane
  ! Assigns the remaining procs to the central band (which are coarser and thus less computing demanding)
  nplanes         = N(4,midband) - N(4,botband)
  procs(midband)  = np - rem
  iplanes         = ceiling(nplanes*1d0/procs(midband))
  i_load(midband) = load_band(midband)/nplanes*iplanes
  ! Balance the load of procs.
    ! Iterates till the central band has a positive number of procs
    iplanes = 1
    do while (procs(midband)<=0)
      iplanes = iplanes + 1
      rem = 0
      nplanes = N(4,botband) - N(4,0)
      procs(botband)  = ceiling(nplanes*1d0/iplanes)
      procs(topband)  = procs(botband)
      rem             = rem + procs(botband) + procs(topband)
      i_load(botband) = load_band(botband)/nplanes*iplanes
      i_load(topband) = load_band(topband)/nplanes*iplanes
      nplanes         = N(4,midband) - N(4,botband)
      procs(midband)  = np - rem
      i_load(midband) = load_band(midband)/nplanes*ceiling(nplanes*1d0/procs(midband))
    end do
    ! Few fine grid planes per proc, and more coarse grid planes such that the computing time is more or less the same
    maxload = 2d0*maxval(i_load)
    do while (maxval(i_load)<maxload)
      iplanes = iplanes+1
      maxload = maxval(i_load)
      rem = 0
      nplanes         = N(4,botband)-N(4,0)
      procs(botband)  = ceiling(nplanes*1d0/iplanes)
      procs(topband)  = procs(botband)
      rem             = rem + procs(botband) + procs(topband)
      i_load(botband) = load_band(botband)/nplanes*iplanes
      i_load(topband) = load_band(topband)/nplanes*iplanes
      nplanes         = N(4,midband) - N(4,botband)
      procs(midband)  = np - rem
      i_load(midband) = load_band(midband)/nplanes*ceiling(nplanes*1d0/procs(midband))
    end do
    ! Once the load has been balanced, we go one step backwards to reduce a bit the load of the 'fine' procs.
    iplanes = iplanes-1
    rem = 0
    nplanes         = N(4,botband) - N(4,0)
    procs(botband)  = ceiling(nplanes*1d0/iplanes)
    procs(topband)  = procs(botband)
    rem             = rem + procs(botband) + procs(topband)
    i_load(botband) = load_band(botband)/nplanes*iplanes
    i_load(topband) = load_band(topband)/nplanes*iplanes
    nplanes         = N(4,midband) - N(4,botband)
    procs(midband)  = np - rem
    i_load(midband) = load_band(midband)/nplanes*ceiling(nplanes*1d0/procs(midband))

  ! TODO
  ! Once we know the number of procs per band ('procs') we asign a particular set of planes to each proc
  !  planelim does NOT include the walls (first and last points in ugrid, vgrid).
  !  For the pgrid it does not include the first and last points for consistency. Although they actually matter
  !    this function is never called for pgrid.
  !   planelim(1,:,:) = planelim(vgrid,bottom(1)/top(2),myid) v grid points
  !   planelim(2,:,:) = planelim(ugrid,bottom(1)/top(2),myid) u/w grid points
  !   planelim(3,:,:) = planelim(pgrid,bottom(1)/top(2),myid) p grid points
  procs_b = procs
  procs_c = procs
  ! u and w grid
  do iband = 1,nband
    ! Compute the load per band and print it to screen
    nplanes = N(4,iband) - N(4,iband-1)
    iplanes = nplanes/procs_b(iband)
    rem     = nplanes - iplanes*procs_b(iband)
    if (rem==0) then
      i_load(iband) = load_band(iband)/nplanes*iplanes
    else
      i_load(iband) = load_band(iband)/nplanes*(iplanes+1)
    end if
    if (myid==0) then
      write(*,27) 'iband',iband,'planes',nplanes,'procs',procs(iband),'    proc_load',i_load(iband)
      27 format(a5,i2,a7,i4,a6,i4,a14,f14.2)
    end if
    ! Assign to each procs the planes they have to handle
    do iproc = 1,rem
      planelim(ugrid,1,iproc-1+procs_b(iband-1)) = (iplanes+1)*(iproc-1) +1   ! From
      planelim(ugrid,2,iproc-1+procs_b(iband-1)) = (iplanes+1)* iproc         ! To
    end do
    do iproc = rem+1,procs_b(iband)
      planelim(ugrid,1,iproc-1+procs_b(iband-1)) = iplanes*(iproc-rem-1) + (iplanes+1)*rem +1
      planelim(ugrid,2,iproc-1+procs_b(iband-1)) = iplanes*(iproc-rem  ) + (iplanes+1)*rem
    end do
    do iproc = 1,procs_b(iband)
      planelim(ugrid,1,iproc-1+procs_b(iband-1)) = planelim(ugrid,1,iproc-1+procs_b(iband-1)) + N(4,iband-1)
      planelim(ugrid,2,iproc-1+procs_b(iband-1)) = planelim(ugrid,2,iproc-1+procs_b(iband-1)) + N(4,iband-1)
    end do   
    procs_b(iband) = procs_b(iband-1) + procs_b(iband)
  end do
  ! v grid
  do iband = 1,nband
    nplanes = N(3,iband) - N(3,iband-1)
    iplanes = nplanes/procs(iband)
    rem     = nplanes - iplanes*procs(iband)
    do iproc = 1,rem
      planelim(vgrid,1,iproc-1+procs(iband-1)) = (iplanes+1)*(iproc-1) +1    ! From
      planelim(vgrid,2,iproc-1+procs(iband-1)) = (iplanes+1)* iproc          ! To
    end do
    do iproc = rem+1,procs(iband)
      planelim(vgrid,1,iproc-1+procs(iband-1)) = iplanes*(iproc-1-rem) + (iplanes+1)*rem +1
      planelim(vgrid,2,iproc-1+procs(iband-1)) = iplanes*(iproc  -rem) + (iplanes+1)*rem 
    end do
    do iproc = 1,procs(iband)
      planelim(vgrid,1,iproc-1+procs(iband-1)) = planelim(vgrid,1,iproc-1+procs(iband-1)) + N(3,iband-1)
      planelim(vgrid,2,iproc-1+procs(iband-1)) = planelim(vgrid,2,iproc-1+procs(iband-1)) + N(3,iband-1)
    end do   
    procs(iband) = procs(iband-1) + procs(iband)
  end do
  ! p grid
  do iband = botband,botband
    nplanes = N(4,iband) - N(4,iband-1) - 2 ! No ghost points
    iplanes = nplanes/procs_c(iband)
    rem     = nplanes - iplanes*procs_c(iband)
    do iproc = 1,rem
      planelim(pgrid,1,iproc-1+procs_c(iband-1)) = (iplanes+1)*(iproc-1) +1   ! From
      planelim(pgrid,2,iproc-1+procs_c(iband-1)) = (iplanes+1)* iproc         ! To
    end do
    do iproc = rem+1,procs_c(iband)
      planelim(pgrid,1,iproc-1+procs_c(iband-1)) = iplanes*(iproc-rem-1) + (iplanes+1)*rem +1
      planelim(pgrid,2,iproc-1+procs_c(iband-1)) = iplanes*(iproc-rem  ) + (iplanes+1)*rem
    end do
    do iproc = 1,procs_c(iband)
      planelim(pgrid,1,iproc-1+procs_c(iband-1)) = planelim(pgrid,1,iproc-1+procs_c(iband-1)) + N(4,iband-1) +1 !E!first point shifted
      planelim(pgrid,2,iproc-1+procs_c(iband-1)) = planelim(pgrid,2,iproc-1+procs_c(iband-1)) + N(4,iband-1) +1 !E!first point shifted
    end do   
    procs_c(iband) = procs_c(iband-1) + procs_c(iband)
  end do
  do iband = midband,midband
    nplanes = N(4,iband) - N(4,iband-1) + 2
    iplanes = nplanes/procs_c(iband)
    rem     = nplanes - iplanes*procs_c(iband)
    do iproc = 1,rem
      planelim(pgrid,1,iproc-1+procs_c(iband-1)) = (iplanes+1)*(iproc-1) +1-1   ! From
      planelim(pgrid,2,iproc-1+procs_c(iband-1)) = (iplanes+1)* iproc      -1   ! To
    end do
    do iproc = rem+1,procs_c(iband)
      planelim(pgrid,1,iproc-1+procs_c(iband-1)) = iplanes*(iproc-rem-1)+(iplanes+1)*rem +1 -1
      planelim(pgrid,2,iproc-1+procs_c(iband-1)) = iplanes*(iproc-rem  )+(iplanes+1)*rem    -1
    end do
    do iproc = 1,procs_c(iband)
      planelim(pgrid,1,iproc-1+procs_c(iband-1)) = planelim(pgrid,1,iproc-1+procs_c(iband-1)) + N(4,iband-1)
      planelim(pgrid,2,iproc-1+procs_c(iband-1)) = planelim(pgrid,2,iproc-1+procs_c(iband-1)) + N(4,iband-1)
    end do   
    procs_c(iband)=procs_c(iband-1)+procs_c(iband)
  end do
  do iband = topband,topband
    nplanes = N(4,iband) - N(4,iband-1) - 2 ! No ghost point2
    iplanes = nplanes/procs_c(iband)
    rem     = nplanes - iplanes*procs_c(iband)
    do iproc = 1,rem
      planelim(pgrid,1,iproc-1+procs_c(iband-1)) = (iplanes+1)*(iproc-1) +1   ! From
      planelim(pgrid,2,iproc-1+procs_c(iband-1)) = (iplanes+1)* iproc         ! To
    end do
    do iproc = rem+1,procs_c(iband)
      planelim(pgrid,1,iproc-1+procs_c(iband-1)) = iplanes*(iproc-rem-1)+(iplanes+1)*rem +1
      planelim(pgrid,2,iproc-1+procs_c(iband-1)) = iplanes*(iproc-rem  )+(iplanes+1)*rem
    end do
    do iproc = 1,procs_c(iband)
      planelim(pgrid,1,iproc-1+procs_c(iband-1)) = planelim(pgrid,1,iproc-1+procs_c(iband-1)) + N(4,iband-1) +1
      planelim(pgrid,2,iproc-1+procs_c(iband-1)) = planelim(pgrid,2,iproc-1+procs_c(iband-1)) + N(4,iband-1) +1
    end do   
    procs_c(iband) = procs_c(iband-1) + procs_c(iband)
  end do

  ! TODO
  ! limPL_incw is like planelim but including first and last points that were previously removed.
  !  It is only used in planes_to_modes_UVP, modes_to_planes_UVP, record_out and stats
  limPL_incw = planelim
  limPL_incw(:,1,0   ) = planelim(:,1,0   ) -1
  limPL_incw(:,2,np-1) = planelim(:,2,np-1) +1

  ! TODO
  ! limPL_excw is like planelim (i.e. EXCluding Walls) but the pgrid includes all points (there is no reason to take them out in first place)
  ! This function is used in ops_in_planes, a few functions when computing the advective term
  limPL_excw = planelim
  !limPL_excw(pgrid,1,0   ) = planelim(pgrid,1,0   ) -1
  !limPL_excw(pgrid,2,np-1) = planelim(pgrid,2,np-1) +1

  if (myid==0) then
    write(*,*) ''
  end if

  ! bandPL(iproc) returns the band the proc works on.
  ! In the physical space, procs only act on a single band
  do iproc = 0,np-1
    do iband = nband,1,-1
      if (iproc<procs(iband)) then
        bandPL(iproc) = iband
      end if
    end do
  end do

  ! Limits of the 'planes' in x and z: A certain thickness in y, and the whole box in x and z
  do iband = nband,1,-1
    if (myid<procs(iband)) then
      igal = Ngal(1,iband)+2
      kgal = Ngal(2,iband)
    end if
  end do

  jgal(ugrid,1) = limPL_excw(ugrid,1,myid)
  jgal(ugrid,2) = limPL_excw(ugrid,2,myid)
  jgal(vgrid,1) = limPL_excw(vgrid,1,myid)
  jgal(vgrid,2) = limPL_excw(vgrid,2,myid)
  jgal(pgrid,1) = limPL_excw(pgrid,1,myid)
  jgal(pgrid,2) = limPL_excw(pgrid,2,myid)
  

  !For FFT planes

  !Bottom band
  proc_planes = floor(1d0*(physlim_bot-jlim(1,ugrid,2)+1)/(np/2))
  rem_planes = (physlim_bot-jlim(1,ugrid,2)+1)-(np/2)*proc_planes
  iplanes = jlim(1,ugrid,2)-1

    do iproc = 0,np/2-1
      limPL_FFT(:,1,iproc) = iplanes + 1
      
      if(iproc<rem_planes)then
        iplanes = iplanes + proc_planes + 1
      else
        iplanes = iplanes + proc_planes
      endif
      limPL_FFT(:,2,iproc) = iplanes
      bandPL_FFT(iproc)=1

    enddo

  !Top band
  proc_planes = floor(1d0*(jlim(2,ugrid,2)-physlim_top+1)/(np/2))
  rem_planes = (jlim(2,ugrid,2)-physlim_top+1)-(np/2)*proc_planes
  iplanes = physlim_top-1

    do iproc = np/2,np-1
      limPL_FFT(:,1,iproc) = iplanes + 1
      
      if(iproc-np/2<rem_planes)then
        iplanes = iplanes + proc_planes + 1
      else
        iplanes = iplanes + proc_planes
      endif
      limPL_FFT(:,2,iproc) = iplanes
      bandPL_FFT(iproc)=3
    enddo
    
    limPL_FFT(vgrid,2,np-1)=limPL_FFT(ugrid,2,np-1)-1


  deallocate(load_band,i_load)

end subroutine

subroutine getini(u1,u2,u3,p,div,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!     GETINI     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Grab the initial condition and initialize some variables

  use declaration
  implicit none

  include 'mpif.h'                                  ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid         ! MPI variables
  integer iband
  type(cfield)  u1(sband:eband),u2(sband:eband),u3(sband:eband)
  type(cfield)  p (sband:eband)
  type(cfield) div(sband:eband)

  if (myid==0) then
    write(*,*) 'Launching...'
  end if

  if (flag_init==1) then       ! Initial conditions borrowed from another 
    if (myid==0) then
      write(*,*) 'starting from multi-block, flat channel?'
      write(*,*) 'start file: ',trim(dirin)//trim(fnameimb)
    end if
    call mblock_ini(u1,u2,u3,p,myid,status,ierr)         ! Initializes u1, u2, u3, p, u1PL, u2PL, u3PL, ppPL
    if (myid==0) then
      call flowrateIm(Qx,u1(midband)%f(N(4,0),1))        ! Column 1 of proc 0 is mode (0,1) [the meeeean]
      if (flag_ctpress/=1) then
        QxT = Qx
      end if
    end if
    iter0 = 0
    mpgz  = 0d0
    t     = 0d0
  else if (flag_init==2) then  ! Continuing simulation
    if (myid==0) then
      write(*,*) 'continuing simulation'
      write(*,*) 'start file:',fnameimb
    end if
    call mblock_ini(u1,u2,u3,p,myid,status,ierr)         ! Initializes u1, u2, u3, p, u1PL, u2PL, u3PL, ppPL
    if (myid==0) then
      call flowrateIm(Qx,u1(midband)%f(N(4,0),1))        ! Column 1 of proc 0 is mode (0,1) [the meeeean]
      ! TODO change this condition for 'flag_ctpress == 0'
      if (flag_ctpress/=1) then  ! Constant flow rate (flag_ctpress == 1, constant pressure gradient)
        QxT = Qx
      end if
    end if
    mpgz = 0d0
  else if (flag_init==3) then  ! Parabolic profile
    if (myid==0) then
      write(*,*) 'parabolic profile'
    end if
    call mblock_ini_parabolic_profile(u1,u2,u3,p,myid,status,ierr)         ! Initializes u1, u2, u3, p, u1PL, u2PL, u3PL, ppPL
    if (myid==0) then
      call flowrateIm(Qx,u1(midband)%f(N(4,0),1))        ! Column 1 of proc 0 is mode (0,1) [the meeeean]
      ! TODO change this condition for 'flag_ctpress == 0'
      if (flag_ctpress/=1) then  ! Constant flow rate (flag_ctpress == 1, constant pressure gradient)
        QxT = Qx
      end if
    end if
    iter0 = 0
    mpgz  = 0d0
    t     = 0d0
  else 
    write(*,*) 'INITIAL CONDITIONS NOT IMPLEMENTED YET'
    call MPI_FINALIZE(ierr)
    stop
  end if

  ! Broadcast the time step and time.
  ! If it's initialized from a previous simulation this value is already known, otherwise it's set to 0
  call MPI_BCAST(iter0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(t    ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
  iter   = iter0
  !  iter0=iter-nstat
  iwrite = iter

  ! 'Probably' this is used to initialize the divergence in the case of a new simulation
  do iband = sband,eband
    call divergence(div(iband)%f,u1(iband)%f,u2(iband)%f,u3(iband)%f,iband,myid)
  end do

  call init_stats(myid)

  if (myid==0) then
    write(*,*) 'Lx    ',Lx
    write(*,*) 'Ly    ',Ly
    write(*,*) 'Lz    ',Lz
    write(*,*) 'Nx    ',N(1,1:nband)
    write(*,*) 'Nz    ',N(2,1:nband)
    write(*,*) 'Nyv   ',N(3,0:nband)
    write(*,*) 'Nyu   ',N(4,0:nband)
    write(*,*) 'Ngalx ',Ngal(1,1:nband)
    write(*,*) 'Ngalz ',Ngal(2,1:nband)
    write(*,*) 'Ngaly ',Ngal(3,1:nband)-Ngal(3,0:nband-1)
    write(*,*) 'Ngalyu',Ngal(4,1:nband)-Ngal(4,0:nband-1)
    write(*,*) ''
    write(*,*) 'dymin ',yu(1)-yu(0)
    write(*,*) 'dymax ',yu((N(4,nband)+1)/2+1)-yu((N(4,nband)+1)/2)
    write(*,*) ''
    write(*,*) 'Re ',Re
    write(*,*) 't  ',t
  end if

end subroutine

subroutine mblock_ini(u1,u2,u3,p,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!   MBLOCK INI   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'                                  ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid         ! MPI variables

  integer iband,j,i,k,column
  real(8) sigmaz,z,sigmax,x,fact
  type(cfield) u1(sband:eband),u2(sband:eband),u3(sband:eband),p(sband:eband)

  u1PL = 0d0
  u2PL = 0d0
  u3PL = 0d0
  ppPL = 0d0
  do iband = sband,eband
    u1(iband)%f = 0d0
    u2(iband)%f = 0d0
    u3(iband)%f = 0d0
    p (iband)%f = 0d0
  end do

  call read_in(myid)
   
  
  call planes_to_modes_UVP(u1,u1PL,2,myid,status,ierr)
  call planes_to_modes_UVP(u2,u2PL,1,myid,status,ierr)
  call planes_to_modes_UVP(u3,u3PL,2,myid,status,ierr)
  call planes_to_modes_UVP(p ,ppPL,3,myid,status,ierr)
  
  ! Chris' trick
  !Pressure reset
!  print *, "pressure reset"
!  do iband=sband,eband
!   p (iband)%f=0d0
!  end do
 
!  if(myid==0)then
!  u1(2)%f(:,1)=u1(2)%f(:,1)+0.243d0!0.607692307d0
!  endif


end subroutine

subroutine mblock_ini_parabolic_profile(u1,u2,u3,p,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! MBLOCK INI PARABOLIC PROFILE !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! TODO Needs to be checked

  use declaration
  implicit none

  include 'mpif.h'                                  ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid         ! MPI variables

  integer iband
  type(cfield) u1(sband:eband),u2(sband:eband),u3(sband:eband),p(sband:eband)
  integer nx,nz
  integer j,j1,j2,ju1,ju2,jv1,jv2,jp1,jp2
  integer, allocatable:: nxxu(:),nzzu(:),nxxv(:),nzzv(:)
  real(8), allocatable:: buffSR(:,:)

  integer column,i,k

  u1PL = 0d0
  u2PL = 0d0
  u3PL = 0d0
  ppPL = 0d0
  do iband = sband,eband
    u1(iband)%f = 0d0
    u2(iband)%f = 0d0
    u3(iband)%f = 0d0
    p (iband)%f = 0d0
  end do

  ! Function to create parabolic profile
  allocate(nxxu(N(4,0):N(4,nband)+1),nzzu(N(4,0):N(4,nband)+1))
  allocate(nxxv(N(3,0):N(3,nband)+1),nzzv(N(3,0):N(3,nband)+1))
  nxxu(N(4,0))=N(1,1)+2
  nzzu(N(4,0))=N(2,1)
  nxxv(N(3,0))=N(1,1)+2
  nzzv(N(3,0))=N(2,1)
  do iband=1,nband
    do j=N(4,iband-1)+1,N(4,iband)
      nxxu(j)=N(1,iband)+2
      nzzu(j)=N(2,iband)
    end do
    do j=N(3,iband-1)+1,N(3,iband)
      nxxv(j)=N(1,iband)+2
      nzzv(j)=N(2,iband)
    end do
  end do
  nxxv(N(3,nband)+1)=N(1,nband)+2
  nzzv(N(3,nband)+1)=N(2,nband)
  nxxu(N(4,nband)+1)=N(1,nband)+2
  nzzu(N(4,nband)+1)=N(2,nband)
  if (myid==0) then
   ju1=jgal(2,1)-1
   ju2=jgal(2,2)
   ju1=max(ju1,N(4,0))
  else
   ju1=jgal(2,1)
   ju2=jgal(2,2)
   if (jgal(2,2)==N(4,nband)) then
     ju2=jgal(2,2)+1
   end if
   ju1=max(ju1,N(4,0))
   ju2=min(ju2,N(4,nband)+1)
  end if
  if (myid==0) then
   jv1=jgal(1,1)-1
   jv2=jgal(1,2)
   jv1=max(jv1,N(3,0))
  else
   jv1=jgal(1,1)
   jv2=jgal(1,2)
   if (jgal(1,2)==N(3,nband)) then
     jv2=jgal(1,2)+1
   end if
   jv1=max(jv1,N(3,0))
   jv2=min(jv2,N(3,nband)+1)
  end if
  if (myid==0) then
   jp1=jgal(3,1)-1
   jp2=jgal(3,2)
   jp1=max(jp1,N(4,0)+1)
  else
   jp1=jgal(3,1)
   jp2=jgal(3,2)
   if (jgal(3,2)==N(4,nband)-1) then
     jp2=jgal(3,2)+1
   end if
   jp1=max(jp1,N(4,0)+1)
   jp2=min(jp2,N(4,nband)+1-1)
  end if
  !!!!!!!!!!!!!!    u1    !!!!!!!!!!!!!!
  do j=ju1,ju2
    nx=nxxu(j)
    nz=nzzu(j)
    allocate(buffSR(nx,nz))
    buffSR     =0d0
    buffSR(1,1)=(0.5d0)*Re*mpgx*(yu(j)**2-1)                    ! Parabolic profile: 1/2*Re*dp/dx*(y^2-1) (Reynolds bulk)
    call buff_to_u(u1PL(1,1,j),buffSR,nx,nz,igal,kgal)
    deallocate(buffSR)
  end do
  !!!!!!!!!!!!!!    u2    !!!!!!!!!!!!!!
  do j=jv1,jv2
    nx=nxxv(j)
    nz=nzzv(j)
    allocate(buffSR(nx,nz))
    buffSR     =0
    buffSR(1,1)=1e-6                                          ! Some noise
    call buff_to_u(u2PL(1,1,j),buffSR,nx,nz,igal,kgal)
    deallocate(buffSR)
  end do
  !!!!!!!!!!!!!!    u3    !!!!!!!!!!!!!!
  do j=ju1,ju2
    nx=nxxu(j)
    nz=nzzu(j)
    allocate(buffSR(nx,nz))
    buffSR     =0d0
    buffSR(1,1)=1e-5                                          ! Some noise 
    call buff_to_u(u3PL(1,1,j),buffSR,nx,nz,igal,kgal)
    deallocate(buffSR)
  end do
  !!!!!!!!!!!!!!    p     !!!!!!!!!!!!!!
  do j=jp1,jp2
    nx=nxxu(j)
    nz=nzzu(j)
    allocate(buffSR(nx,nz))
    buffSR     =0d0
    buffSR(1,1)=1d0                                           ! Some noise 
    call buff_to_u(ppPL(1,1,j),buffSR,nx,nz,igal,kgal)
    deallocate(buffSR)
  end do
  deallocate(nxxu,nzzu,nxxv,nzzv)

  call planes_to_modes_UVP(u1,u1PL,ugrid,myid,status,ierr)
  call planes_to_modes_UVP(u2,u2PL,vgrid,myid,status,ierr)
  call planes_to_modes_UVP(u3,u3PL,ugrid,myid,status,ierr)
  call planes_to_modes_UVP(p ,ppPL,pgrid,myid,status,ierr)

end subroutine

subroutine read_in(myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    READ IN     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reads the file with the inital condition

! The variables u1PL, u2PL, u3PL and ppPL are sent to every proc
! The procs only stores the planes they have to compute
!  (jgal is the local name of planelim for j)

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid
  integer nx,nz,nxin,nzin,nband2,iband
  integer j,jin,iproc,dummI,ju1,ju2,jv1,jv2,jp1,jp2
  real(8) dummRe,Re2,alp2,bet2
  real(8), allocatable:: buffSR(:,:),dumm_y(:)
  integer, allocatable:: dummint(:),N2(:,:)
  integer, allocatable:: nxxu(:),nzzu(:),nxxv(:),nzzv(:),nxxp(:),nzzp(:)

  filout   = fnameimb(3:index(fnameimb,' ')-1)

  if (myid==0) then
    ! Read the file with the initial conditions
    fnameimb = trim(dirin)//'/u1'//filout
    open(10,file=fnameimb,form='unformatted')
    allocate(dummint(88))
    read(10) t,Re2,alp2,bet2,dummRe,nband2,iter0,dummint
!    if (flag_ctpress==1) then
!      mpgx=dummRe
!    end if
    deallocate(dummint)

    ! Checks that the size of the variables in the file matches with the setup of actual problem
    ! Checks nband, alp, bet, N
    if (nband2/=nband) then
      write(*,*) ''
      write(*,*) 'wrong startfile?'
      write(*,*) 'nband,nband_old=',nband,nband2
      stop
    else if (alp2/=alp .or. bet2/=bet) then
      write(*,*) ''
      write(*,*) 'wrong startfile?'
      write(*,*) 'alp,alp_old=',alp,alp2
      write(*,*) 'bet,bet_old=',bet,bet2
      write(*,*) ''
    end if
    write(*,*) 'Re_old,Re',Re2,Re
    allocate(N2(4,0:nband2+1))
    read(10) N2
    close(10)

    dummI=0
    do iband=1,nband
      if (N(1,iband)/=N2(1,iband)) then
        dummI=2
      else if (N(2,iband)/=N2(2,iband)) then
        dummI=2
      else if (N(3,iband)/=N2(3,iband)) then
        dummI=2
      end if
    end do
    if (N(3,0)/=N2(3,0)) then
      dummI=2
    else if (N(3,nband)/=N2(3,nband)) then
      dummI=2
    end if
    if (dummI==1) then
      write(*,*) ''
      write(*,*) 'wrong startfile?',dummI
      write(*,*) 'N_new='
      write(*,*) N(1,0:nband+1)
      write(*,*) N(2,0:nband+1)
      write(*,*) N(3,0:nband+1)
      write(*,*) N(4,0:nband+1)
      write(*,*) ''
      write(*,*) 'N_old='
      write(*,*) N2(1,0:nband+1)
      write(*,*) N2(2,0:nband+1)
      write(*,*) N2(3,0:nband+1)
      write(*,*) N2(4,0:nband+1)
      stop
    else if (dummI==2) then
      write(*,*) ''
      write(*,*) 'WARNING startfile size',dummI
      write(*,*) 'N_new='
      write(*,*) N(1,0:nband+1)
      write(*,*) N(2,0:nband+1)
      write(*,*) N(3,0:nband+1)
      write(*,*) N(4,0:nband+1)
      write(*,*) ''
      write(*,*) 'N_old='
      write(*,*) N2(1,0:nband+1)
      write(*,*) N2(2,0:nband+1)
      write(*,*) N2(3,0:nband+1)
      write(*,*) N2(4,0:nband+1)
      write(*,*) ''
    end if

    ! Despite possible unmatches, it sends N to the procs
    do iproc=1,np-1
      call MPI_SEND(N2,4*(nband2+2),MPI_INTEGER,iproc,121*iproc,MPI_COMM_WORLD,ierr)
    end do

    ! nxx(j) and nzz(j) store the number of points at a plane j
    allocate(nxxu(N2(4,0):N2(4,nband)+1),nzzu(N2(4,0):N2(4,nband)+1))
    allocate(nxxv(N2(3,0):N2(3,nband)+1),nzzv(N2(3,0):N2(3,nband)+1))
    allocate(nxxp(N2(4,0)+1:N2(4,nband)+1-1),nzzp(N2(4,0)+1:N2(4,nband)+1-1))
    nxxu(N2(4,0))=N2(1,1)+2
    nzzu(N2(4,0))=N2(2,1)
    nxxv(N2(3,0))=N2(1,1)+2
    nzzv(N2(3,0))=N2(2,1)
    nxxp(N2(4,0)+1)=N2(1,1)+2
    nzzp(N2(4,0)+1)=N2(2,1)
    do iband=1,nband
      do j=N2(4,iband-1)+1,N2(4,iband)
        nxxu(j)=N2(1,iband)+2
        nzzu(j)=N2(2,iband)
      end do
      do j=N2(3,iband-1)+1,N2(3,iband)
        nxxv(j)=N2(1,iband)+2
        nzzv(j)=N2(2,iband)
      end do
      do j=max(N2(4,iband-1),N2(4,0)+1),min(N2(4,iband),N2(4,nband)-1)
        nxxp(j)=N2(1,iband)+2
        nzzp(j)=N2(2,iband)
      end do
      nxxp(N2(4,1))=nxxp(N2(4,1)+1)
      nzzp(N2(4,1))=nzzp(N2(4,1)+1)
      nxxp(N2(4,2))=nxxp(N2(4,2)-1)
      nzzp(N2(4,2))=nzzp(N2(4,2)-1)
      nxxp(N2(4,2)+1)=nxxp(N2(4,2)-1)
      nzzp(N2(4,2)+1)=nzzp(N2(4,2)-1)
    end do
    nxxu(N2(4,nband)+1)=N2(1,nband)+2
    nzzu(N2(4,nband)+1)=N2(2,nband)
    nxxv(N2(3,nband)+1)=N2(1,nband)+2
    nzzv(N2(3,nband)+1)=N2(2,nband)
    nxxp(N2(4,nband)+1-1)=N2(1,nband)+2
    nzzp(N2(4,nband)+1-1)=N2(2,nband)
    !filout=fnameimb(3:index(fnameimb,' ')-1)
    !filout=fnameimb(10:index(fnameimb,' ')-1) !22 as now in subfolder

    write(*,*) 'getting'

    ! Reads u1, u2, u3 and p,  and send them to the procs

    !!!!!!!!!!!!!!    u1    !!!!!!!!!!!!!!
    fnameimb = trim(dirin)//'/u1'//filout
    write(*,*) 'u1 from file ',trim(fnameimb),','
    open(10,file=fnameimb,form='unformatted')
    read(10)
    read(10)
    read(10)
    ju1=jgal(2,1)-1
    ju2=jgal(2,2)
    ju1=max(ju1,N2(4,0))
    do j=N2(4,0),ju1-1
      read(10)
    end do
    do j=ju1,ju2
      nx=nxxu(j)
      nz=nzzu(j)
      allocate(buffSR(nx,nz))
      read(10) jin,dummI,nxin,nzin,dummRe,buffSR
      call buff_to_u(u1PL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
      if (nx/=nxin .or. nz/=nzin) then
        write(*,*) 'WARNING!: unexpected size of plane',j
      end if
    end do
    do iproc=1,np-1
      ju1=planelim(2,1,iproc)
      ju2=planelim(2,2,iproc)
      if (planelim(2,2,iproc)==N(4,nband)) then
        ju2=planelim(2,2,iproc)+1
      end if
      ju1=max(ju1,N2(4,0))
      ju2=min(ju2,N2(4,nband)+1)
      do j=ju1,ju2
        nx=nxxu(j)
        nz=nzzu(j)
        allocate(buffSR(nx,nz))
        read(10) jin,dummI,nxin,nzin,dummRe,buffSR
        call MPI_SEND(buffSR,nx*nz,MPI_REAL8,iproc,123*iproc,MPI_COMM_WORLD,ierr)
        deallocate(buffSR)
        if (nx/=nxin .or. nz/=nzin) then
          write(*,*) 'WARNING!: unexpected size of plane',j
        end if
      end do
    end do
    close(10)

    !!!!!!!!!!!!!!    u2    !!!!!!!!!!!!!!
    fnameimb = trim(dirin)//'/u2'//filout
    write(*,*) 'u2 from file ',trim(fnameimb),','
    open(10,file=fnameimb,form='unformatted')
    read(10)
    read(10)
    read(10)
    jv1=jgal(1,1)-1
    jv2=jgal(1,2)
    jv1=max(jv1,N2(3,0))
    do j=N2(3,0),jv1-1
      read(10)
    end do
    do j=jv1,jv2
      nx=nxxv(j)
      nz=nzzv(j)
      allocate(buffSR(nx,nz))
      read(10) jin,dummI,nxin,nzin,dummRe,buffSR
      call buff_to_u(u2PL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
      if (nx/=nxin .or. nz/=nzin) then
        write(*,*) 'WARNING!: unexpected size of plane',j
      end if
    end do
    do iproc=1,np-1
      jv1=planelim(1,1,iproc)
      jv2=planelim(1,2,iproc)
      if (planelim(1,2,iproc)==N(3,nband).and.iproc==np-1) then
        jv2=planelim(1,2,iproc)+1
      end if
      jv1=max(jv1,N2(3,0))
      jv2=min(jv2,N2(3,nband)+1)
!jv2=min(jv2,N2(3,nband))
      do j=jv1,jv2
        nx=nxxv(j)
        nz=nzzv(j)
        allocate(buffSR(nx,nz))
        read(10) jin,dummI,nxin,nzin,dummRe,buffSR
        call MPI_SEND(buffSR,nx*nz,MPI_REAL8,iproc,124*iproc,MPI_COMM_WORLD,ierr)
        deallocate(buffSR)
        if (nx/=nxin .or. nz/=nzin) then
          write(*,*) 'WARNING!: unexpected size of plane',j
        end if
      end do
    end do
    close(10)

    !!!!!!!!!!!!!!    u3    !!!!!!!!!!!!!!
    fnameimb = trim(dirin)//'/u3'//filout
    write(*,*) 'u3 from file ',trim(fnameimb),','
    open(10,file=fnameimb,form='unformatted')
    read(10)
    read(10)
    read(10)
    ju1=jgal(2,1)-1
    ju2=jgal(2,2)
    ju1=max(ju1,N2(4,0))
    do j=N2(4,0),ju1-1
      read(10)
    end do
    do j=ju1,ju2
      nx=nxxu(j)
      nz=nzzu(j)
      allocate(buffSR(nx,nz))
      read(10) jin,dummI,nxin,nzin,dummRe,buffSR
      call buff_to_u(u3PL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
      if (nx/=nxin .or. nz/=nzin) then
        write(*,*) 'WARNING!: unexpected size of plane',j
      end if
    end do
    do iproc=1,np-1
      ju1=planelim(2,1,iproc)
      ju2=planelim(2,2,iproc)
      if (planelim(2,2,iproc)==N(4,nband)) then
        ju2=planelim(2,2,iproc)+1
      end if
      ju1=max(ju1,N2(4,0))
      ju2=min(ju2,N2(4,nband)+1)
      do j=ju1,ju2
        nx=nxxu(j)
        nz=nzzu(j)
        allocate(buffSR(nx,nz))
        read(10) jin,dummI,nxin,nzin,dummRe,buffSR
        call MPI_SEND(buffSR,nx*nz,MPI_REAL8,iproc,125*iproc,MPI_COMM_WORLD,ierr)
        deallocate(buffSR)
        if (nx/=nxin .or. nz/=nzin) then
          write(*,*) 'WARNING!: unexpected size of plane',j
        end if
      end do
    end do
    close(10)

    !!!!!!!!!!!!!!    p     !!!!!!!!!!!!!!

    fnameimb = trim(dirin)//'/p'//filout
    write(*,*) 'and p from file ',trim(fnameimb),','
    open(10,file=fnameimb,form='unformatted')
    read(10)
    read(10)
    read(10)
    jp1=jgal(3,1)-1
    jp2=jgal(3,2)
    jp1=max(jp1,N2(4,0)+1)
    do j=N2(4,0)+1,jp1-1
      read(10)
    end do
    do j=jp1,jp2
      nx=nxxp(j)
      nz=nzzp(j)
      allocate(buffSR(nx,nz))
      read(10) jin,dummI,nxin,nzin,dummRe,buffSR
      call buff_to_u(ppPL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
      if (nx/=nxin .or. nz/=nzin) then
        write(*,*) 'WARNING!: unexpected size of plane',j
      end if
    end do
    do iproc=1,np-1
      jp1=planelim(3,1,iproc)
      jp2=planelim(3,2,iproc)
      if (planelim(3,2,iproc)==N(4,nband)-1.and.iproc==np-1) then
        jp2=planelim(3,2,iproc)+1
      end if
      jp1=max(jp1,N2(4,0)+1)
      jp2=min(jp2,N2(4,nband)+1-1)
!jp2=min(jp2,N2(4,nband)+1-1-1)
      do j=jp1,jp2
        nx=nxxp(j)
        nz=nzzp(j)
        allocate(buffSR(nx,nz))
        read(10) jin,dummI,nxin,nzin,dummRe,buffSR
        call MPI_SEND(buffSR,nx*nz,MPI_REAL8,iproc,126*iproc,MPI_COMM_WORLD,ierr)
        deallocate(buffSR)
        if (nx/=nxin .or. nz/=nzin) then
          write(*,*) 'WARNING!: unexpected size of plane',j
        end if
      end do
    end do
    close(10)
    write(*,*) ''

    deallocate(nxxu,nzzu,nxxv,nzzv,nxxp,nzzp)
    deallocate(N2)

  else

    ! The procs receive u1, u2, u3 and p
    ! Those variables are stored in u1PL, u2PL, u3PL and ppPL
    ! The procs only stores the planes they have to compute
    !  (jgal is the local name of planelim for j)

    allocate(N2(4,0:nband+1))
    call MPI_RECV(N2,4*(nband+2),MPI_INTEGER,0,121*myid,MPI_COMM_WORLD,status,ierr)
    allocate(nxxu(N2(4,0):N2(4,nband)+1),nzzu(N2(4,0):N2(4,nband)+1))
    allocate(nxxv(N2(3,0):N2(3,nband)+1),nzzv(N2(3,0):N2(3,nband)+1))
    allocate(nxxp(N2(4,0)+1:N2(4,nband)+1-1),nzzp(N2(4,0)+1:N2(4,nband)+1-1))
    nxxu(N2(4,0))=N2(1,1)+2
    nzzu(N2(4,0))=N2(2,1)
    nxxv(N2(3,0))=N2(1,1)+2
    nzzv(N2(3,0))=N2(2,1)
    nxxp(N2(4,0)+1)=N2(1,1)+2
    nzzp(N2(4,0)+1)=N2(2,1)
    do iband=1,nband
      do j=N2(4,iband-1)+1,N2(4,iband)
        nxxu(j)=N2(1,iband)+2
        nzzu(j)=N2(2,iband)
      end do
      do j=N2(3,iband-1)+1,N2(3,iband)
        nxxv(j)=N2(1,iband)+2
        nzzv(j)=N2(2,iband)
      end do
      do j=max(N2(4,iband-1),N2(4,0)+1),min(N2(4,iband),N2(4,nband)-1)
        nxxp(j)=N2(1,iband)+2
        nzzp(j)=N2(2,iband)
      end do
      nxxp(N2(4,1))=nxxp(N2(4,1)+1)
      nzzp(N2(4,1))=nzzp(N2(4,1)+1)
      nxxp(N2(4,2))=nxxp(N2(4,2)-1)
      nzzp(N2(4,2))=nzzp(N2(4,2)-1)
      nxxp(N2(4,2)+1)=nxxp(N2(4,2)-1)
      nzzp(N2(4,2)+1)=nzzp(N2(4,2)-1)
    end do
    nxxu(N2(4,nband)+1)=N2(1,nband)+2
    nzzu(N2(4,nband)+1)=N2(2,nband)
    nxxv(N2(3,nband)+1)=N2(1,nband)+2
    nzzv(N2(3,nband)+1)=N2(2,nband)
    nxxp(N2(4,nband)+1-1)=N2(1,nband)+2
    nzzp(N2(4,nband)+1-1)=N2(2,nband)
   
    ju1=jgal(2,1)
    ju2=jgal(2,2)
    if (jgal(2,2)==N(4,nband)) then
      ju2=jgal(2,2)+1
    end if
    ju1=max(ju1,N2(4,0))
    ju2=min(ju2,N2(4,nband)+1)
    jv1=jgal(1,1)
    jv2=jgal(1,2)
    if (jgal(1,2)==N(3,nband).and.myid==np-1) then
      jv2=jgal(1,2)+1
    end if
    jv1=max(jv1,N2(3,0))
    jv2=min(jv2,N2(3,nband)+1)
!jv2=min(jv2,N2(3,nband)+1-1)
    jp1=jgal(3,1)
    jp2=jgal(3,2)
    if (jgal(3,2)==N(4,nband)-1.and.myid==np-1) then
      jp2=jgal(3,2)+1
    end if
    jp1=max(jp1,N2(4,0)+1)
    jp2=min(jp2,N2(4,nband)+1-1)
!jp2=min(jp2,N2(4,nband)+1-1-1)
    !!!!!!!!!!!!!!    u1    !!!!!!!!!!!!!!
    do j=ju1,ju2
      nx=nxxu(j)
      nz=nzzu(j)
      allocate(buffSR(nx,nz))
      call MPI_RECV(buffSR,nx*nz,MPI_REAL8,0,123*myid,MPI_COMM_WORLD,status,ierr)
      call buff_to_u(u1PL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
    end do
    !!!!!!!!!!!!!!    u2    !!!!!!!!!!!!!!
    do j=jv1,jv2
      nx=nxxv(j)
      nz=nzzv(j)
      allocate(buffSR(nx,nz))
      call MPI_RECV(buffSR,nx*nz,MPI_REAL8,0,124*myid,MPI_COMM_WORLD,status,ierr)
      call buff_to_u(u2PL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
    end do
    !!!!!!!!!!!!!!    u3    !!!!!!!!!!!!!!
    do j=ju1,ju2
      nx=nxxu(j)
      nz=nzzu(j)
      allocate(buffSR(nx,nz))
      call MPI_RECV(buffSR,nx*nz,MPI_REAL8,0,125*myid,MPI_COMM_WORLD,status,ierr)
      call buff_to_u(u3PL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
    end do
    !!!!!!!!!!!!!!    p     !!!!!!!!!!!!!!
    do j=jp1,jp2
      nx=nxxp(j)
      nz=nzzp(j)
      allocate(buffSR(nx,nz))
      call MPI_RECV(buffSR,nx*nz,MPI_REAL8,0,126*myid,MPI_COMM_WORLD,status,ierr)
      call buff_to_u(ppPL(1,1,j),buffSR,nx,nz,igal,kgal)
      deallocate(buffSR)
    end do
    deallocate(nxxu,nzzu,nxxv,nzzv,nxxp,nzzp)
    deallocate(N2)
  end if

end subroutine

subroutine init_stats(myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    init stats  !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer ju1,ju2,jv1,jv2,jp1,jp2,myid

  allocate(Um (limPL_incw(ugrid,1,myid):limPL_incw(ugrid,2,myid)),U2m (limPL_incw(ugrid,1,myid):limPL_incw(ugrid,2,myid)))
  allocate(Vm (limPL_incw(vgrid,1,myid):limPL_incw(vgrid,2,myid)),V2m (limPL_incw(vgrid,1,myid):limPL_incw(vgrid,2,myid)))
  allocate(Wm (limPL_incw(ugrid,1,myid):limPL_incw(ugrid,2,myid)),W2m (limPL_incw(ugrid,1,myid):limPL_incw(ugrid,2,myid)))
  allocate(Pm (limPL_incw(pgrid,1,myid):limPL_incw(pgrid,2,myid)),P2m (limPL_incw(pgrid,1,myid):limPL_incw(pgrid,2,myid)))
  allocate(wxm(limPL_incw(vgrid,1,myid):limPL_incw(vgrid,2,myid)),wx2m(limPL_incw(vgrid,1,myid):limPL_incw(vgrid,2,myid)))
  allocate(UVm(limPL_incw(ugrid,1,myid):limPL_incw(ugrid,2,myid)))
  allocate(UWm(limPL_incw(ugrid,1,myid):limPL_incw(ugrid,2,myid)))
  allocate(VWm(limPL_incw(vgrid,1,myid):limPL_incw(vgrid,2,myid)))
  
  if (bandPL(myid)==1) then
    ju1 = limPL_incw(ugrid,1,myid)-N(4,0)
    ju2 = limPL_incw(ugrid,2,myid)-N(4,0)
    jv1 = limPL_incw(vgrid,1,myid)-N(3,0)
    jv2 = limPL_incw(vgrid,2,myid)-N(3,0)
    jp1 = limPL_incw(pgrid,1,myid)-(N(4,0)+1)
    jp2 = limPL_incw(pgrid,2,myid)-(N(4,0)+1)
    allocate(UmC (dnx,dnz,ju1:ju2),U2mC (dnx,dnz,ju1:ju2))
    allocate(VmC (dnx,dnz,jv1:jv2),V2mC (dnx,dnz,jv1:jv2))
    allocate(WmC (dnx,dnz,ju1:ju2),W2mC (dnx,dnz,ju1:ju2))
    allocate(PmC (dnx,dnz,jp1:jp2),P2mC (dnx,dnz,jp1:jp2))
    allocate(wxmC(dnx,dnz,jv1:jv2),wx2mC(dnx,dnz,jv1:jv2))
    allocate(UVmC(dnx,dnz,ju1:ju2))
    allocate(UWmC(dnx,dnz,ju1:ju2))
    allocate(VWmC(dnx,dnz,jv1:jv2))  
    UmC   = 0d0
    U2mC  = 0d0
    VmC   = 0d0
    V2mC  = 0d0
    WmC   = 0d0
    W2mC  = 0d0
    PmC   = 0d0
    P2mC  = 0d0
    UVmC  = 0d0
    UWmC  = 0d0
    VWmC  = 0d0
    wxmC  = 0d0
    wx2mC = 0d0
  else if (bandPL(myid)==nband) then
    ju1 = -limPL_incw(ugrid,2,myid)+N(4,nband)+1
    ju2 = -limPL_incw(ugrid,1,myid)+N(4,nband)+1
    jv1 = -limPL_incw(vgrid,2,myid)+N(3,nband)+1
    jv2 = -limPL_incw(vgrid,1,myid)+N(3,nband)+1
    jp1 = -limPL_incw(pgrid,2,myid)+N(4,nband)
    jp2 = -limPL_incw(pgrid,1,myid)+N(4,nband)
    allocate(UmC (dnx,dnz,ju1:ju2),U2mC (dnx,dnz,ju1:ju2))
    allocate(VmC (dnx,dnz,jv1:jv2),V2mC (dnx,dnz,jv1:jv2))
    allocate(WmC (dnx,dnz,ju1:ju2),W2mC (dnx,dnz,ju1:ju2))
    allocate(PmC (dnx,dnz,jp1:jp2),P2mC (dnx,dnz,jp1:jp2))
    allocate(wxmC(dnx,dnz,jv1:jv2),wx2mC(dnx,dnz,jv1:jv2))
    allocate(UVmC(dnx,dnz,ju1:ju2))
    allocate(UWmC(dnx,dnz,ju1:ju2))
    allocate(VWmC(dnx,dnz,jv1:jv2))
    UmC   = 0d0
    U2mC  = 0d0
    VmC   = 0d0
    V2mC  = 0d0
    WmC   = 0d0
    W2mC  = 0d0
    PmC   = 0d0
    P2mC  = 0d0
    UVmC  = 0d0
    UWmC  = 0d0
    VWmC  = 0d0
    wxmC  = 0d0
    wx2mC = 0d0
  end if

  Um    = 0d0
  U2m   = 0d0
  Vm    = 0d0
  V2m   = 0d0
  Wm    = 0d0
  W2m   = 0d0
  Pm    = 0d0
  P2m   = 0d0
  UVm   = 0d0
  UWm   = 0d0
  VWm   = 0d0
  wxm   = 0d0
  wx2m  = 0d0
  istat = 0

 if (myid==0) then
   write(ext4,'(i5.5)') int(t)
   fnameimb = 'stats_'//ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
 end if

end subroutine

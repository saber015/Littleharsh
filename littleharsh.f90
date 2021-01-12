! *************************************************************************************************************** !
!
! Last modified by C 12/05/2016
!   - Set up for SHS
!
!  TODO
!    - Pressure correction of ghost points
!    - flowrate_corr LUsol0 boundary conditions if not no-slip?
!
! *************************************************************************************************************** !
!
! Contains:
!       littleharsh - Main program : Calls start, getini, RHS0_u1, RHS0_u2, RHS0_u3, nonlinear, solveU, solveV, meanflow_ctU, SolveP, v_corr, meanflow_ctP, finalize
!                     - Does the main loop
!        divergence - Calculates the divergence
!                     - Called from a few places (inc. v_corr, getini, solveP)
!       laplacian_U - Calculates the laplacian of u/w
!                     - Called from RHS0_u1 and RHS0_u3
!       laplacian_V - Calculates the laplacian of v
!                     - Called from RHS0_u2
!             error - Calculates the error
!                     - Called in record out and finalize
!                     - maximum of the divergance
!                       - Should be machine round-off
!          finalize - Records out last step
!                     - Called at last step from littleharsh
!                     - Should be renamed finalise
!        flowrateIm - Calculates the mass flow rate (for real u?)
!                     - Called from flowrate_corr and meanflow_ctP
!        flowrateRe - Calculates the mass flow rate
!                     - Called from flowrate_corr
!     flowrate_corr - Calculates correction to mean pressure gradient and u for constant mass flow
!                     - Called from meanflow_ctU
!                     - Calls from flowrateRe, flowrateIm, LUsol0
!            maxvel - Calculates the mximum velocitiy
!                     - Called from record_out
!      meanflow_ctP - Calculates instantainous mass flow rate if constant pressure gradient used 
!                     - Called form littleharsh
!                     - Calls flowrateIm
!      meanflow_ctU - Calculates corrections to ensure constant mass flow rate
!                     - Called form littleharsh
!                     - Calls flowrate_corr
!            solveP - Solves the pressure 
!                     - Called form littleharsh
!                     - Calls LUsolP
!           RHS0_u1 - Solves the RHS of u1 (excluding lastest nonlinear term) includes mpg
!                     - Called form littleharsh
!           RHS0_u2 - Solves the RHS of u2 (excluding lastest nonlinear term)
!                     - Called form littleharsh
!           RHS0_u3 - Solves the RHS of u3 (excluding lastest nonlinear term)
!                     - Called form littleharsh
!            solveU - Solves the u/w velocities 
!                     - Called form littleharsh
!                     - Calls LUsol
!            solveV - Solves the v velocity
!                     - Called form littleharsh
!                     - Calls LUsol
!            v_corr - Pressure correction step
!                     - Called form littleharsh
!
! *************************************************************************************************************** !

program littleharsh

  use declaration
  implicit none

  include 'mpif.h'                         
  integer status(MPI_STATUS_SIZE),ierr,myid 

  integer iband
  type(cfield), allocatable:: u1(:),u2(:),u3(:),p(:),psi(:)
  type(cfield), allocatable:: div(:)
  type(cfield), allocatable:: Nu1(:),Nu2(:),Nu3(:)
  type(cfield), allocatable:: du1(:),du2(:),du3(:)
  type(cfield), allocatable:: Lu1(:),Lu2(:),Lu3(:)
  
  real(8):: z,sigmaz
  integer column,i,k,j

  real(8):: flagslstat

  call MPI_INIT(ierr)
  if (ierr.ne.MPI_SUCCESS) pause 'ierr'
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)      ! np   = Number of processes in the group of comm
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)    ! myid = rank of the calling process

  if (myid==0) then
    write(*,*) ''
    write(*,*) '---------------LITTLE HARSH---------------'
    write(*,*) "         'That's a little harsh'          "
    write(*,*) '------------------------------------------'
    write(*,*) ''
  end if

  ! Initialise
  call start(myid,status,ierr)

itersl=iter0
nstatsl = nwrite
  
  ! Allocate memory for main variables
  allocate( u1  (sband:eband), u2 (sband:eband), u3 (sband:eband))
  allocate( Lu1 (sband:eband), Lu2(sband:eband), Lu3(sband:eband))
  allocate(du1  (sband:eband), du2(sband:eband), du3(sband:eband))
  allocate(Nu1  (sband:eband), Nu2(sband:eband), Nu3(sband:eband))
  allocate( p  (sband:eband))
  allocate( psi(sband:eband))
  allocate( div(sband:eband))
  do iband = sband,eband
    allocate( u1  (iband)%f(jlim(1,ugrid,iband)  :jlim(2,ugrid,iband)  ,columns_num(iband,myid)))
    allocate( u2  (iband)%f(jlim(1,vgrid,iband)  :jlim(2,vgrid,iband)  ,columns_num(iband,myid)))
    allocate( u3  (iband)%f(jlim(1,ugrid,iband)  :jlim(2,ugrid,iband)  ,columns_num(iband,myid)))
    allocate(Lu1  (iband)%f(jlim(1,ugrid,iband)  :jlim(2,ugrid,iband)  ,columns_num(iband,myid)))
    allocate(Lu2  (iband)%f(jlim(1,vgrid,iband)  :jlim(2,vgrid,iband)  ,columns_num(iband,myid)))
    allocate(Lu3  (iband)%f(jlim(1,ugrid,iband)  :jlim(2,ugrid,iband)  ,columns_num(iband,myid)))
    allocate(du1  (iband)%f(jlim(1,ugrid,iband)  :jlim(2,ugrid,iband)  ,columns_num(iband,myid)))
    allocate(du2  (iband)%f(jlim(1,vgrid,iband)  :jlim(2,vgrid,iband)  ,columns_num(iband,myid)))
    allocate(du3  (iband)%f(jlim(1,ugrid,iband)  :jlim(2,ugrid,iband)  ,columns_num(iband,myid)))
    allocate(Nu1  (iband)%f(jlim(1,ugrid,iband)+1:jlim(2,ugrid,iband)-1,columns_num(iband,myid)))
    allocate(Nu2  (iband)%f(jlim(1,vgrid,iband)+1:jlim(2,vgrid,iband)-1,columns_num(iband,myid)))
    allocate(Nu3  (iband)%f(jlim(1,ugrid,iband)+1:jlim(2,ugrid,iband)-1,columns_num(iband,myid)))
    allocate( p   (iband)%f(jlim(1,pgrid,iband)  :jlim(2,pgrid,iband)  ,columns_num(iband,myid)))
    allocate( psi (iband)%f(jlim(1,pgrid,iband)  :jlim(2,pgrid,iband)  ,columns_num(iband,myid)))
    allocate( div (iband)%f(jlim(1,pgrid,iband)  :jlim(2,pgrid,iband)  ,columns_num(iband,myid)))
  end do

  do iband = sband,eband
     u1  (iband)%f = 0d0
     u2  (iband)%f = 0d0
     u3  (iband)%f = 0d0
    Lu1  (iband)%f = 0d0
    Lu2  (iband)%f = 0d0
    Lu3  (iband)%f = 0d0
    du1  (iband)%f = 0d0
    du2  (iband)%f = 0d0
    du3  (iband)%f = 0d0
    Nu1  (iband)%f = 0d0
    Nu2  (iband)%f = 0d0
    Nu3  (iband)%f = 0d0
     p   (iband)%f = 0d0
     psi (iband)%f = 0d0
     div (iband)%f = 0d0
  end do

  ! Get the initial conditions
  call getini(u1,u2,u3,p,div,myid,status,ierr)

nextqt = floor(t*10d0)/10d0+0.1d0
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (myid==0) then
    write(*,*) ''
    write(*,*) 'Initialised'
  end if
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! MAIN LOOP 
  do while (t<maxt) ! This is the original condition
  !do while (t<maxt .AND. iter < 200)
    ! Runge-Kutta substeps
    do kRK = 1,3
!    do kRK = 1,1
      ! Build     linear terms of right-hand-side of Navier-Stokes equation
      call RHS0_u1(du1,u1,Nu1,p,Lu1,myid)
      call RHS0_u2(du2,u2,Nu2,p,Lu2,myid)
      call RHS0_u3(du3,u3,Nu3,p,Lu3,myid)
      
      ! Build non-linear terms of right-hand-side of Navier-Stokes equation
      call nonlinear(Nu1,Nu2,Nu3,u1,u2,u3,du1,du2,du3,p,div,myid,status,ierr)
      ! Resolve the matricial system
      call solveU(u1,du1,Lu1,myid)
      call solveV(u2,du2,Lu2,myid)
      call solveU(u3,du3,Lu3,myid)
      
     ! Compute the pressure gradient if constant mass flow condition is set
     if (flag_ctpress==0) then
       if (myid == 0) then
         call meanflow_ctU(u1,myid) ! Computes the convenient mpgx
       endif
     end if
      
      ! Solves pressure matrix
      call solveP(p,psi,u1,u2,u3,div,myid)

      ! Correct velocity for incompressibility
      call v_corr(u1,u2,u3,psi,div,myid,status,ierr)
    
    end do
    
    ! Evaluate mass flow if constant pressure gradient is set
    if (flag_ctpress==1) then
      if (myid == 0) then
        call meanflow_ctP(u1,myid)
      end if
    end if

    iter = iter+1

  end do 

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! Deallocate memory and shut communications
  call finalize(u1,u2,u3,p,div,myid,status,ierr)

end program

subroutine divergence(div,u1,u2,u3,iband,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   DIVERGENCE   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Divergence in Fourier space. (Pencils)
! Get the result of the divergence at the centers (pgrid)
! Div(u) =  du/dx +  dw/dz + dv/dy
!        = i*kx*u + i*kz*w + FinDiff(v) 
!        = i*kx*u + i*kz*w + (v(i) - v(i-1))/(DeltaTheta)*(dtheta/dyu)

  use declaration
  implicit none

  integer i,k,j,iband,myid,column
  complex(8) div(jlim(1,pgrid,iband):jlim(2,pgrid,iband),columns_num(iband,myid))
  complex(8)  u1(jlim(1,ugrid,iband):jlim(2,ugrid,iband),columns_num(iband,myid))
  complex(8)  u2(jlim(1,vgrid,iband):jlim(2,vgrid,iband),columns_num(iband,myid))
  complex(8)  u3(jlim(1,ugrid,iband):jlim(2,ugrid,iband),columns_num(iband,myid))
  complex(8) kx,kzF

  do column = 1,columns_num(iband,myid)
    i   = columns_i(column,iband,myid)
    k   = columns_k(column,iband,myid)
    kx  = k1F_x(i)  ! im*kx
    kzF = k1F_z(k)  ! im*kz
    !!!!!!!!      First band, finite diffs      !!!!!!!!
    do j = jlim(1,pgrid,iband),jlim(2,pgrid,iband)
      div(j,column)= kx                 *  u1(j  ,column) &
&                   +kzF                *  u3(j  ,column) &
&                   +ddthetavi*dthdyu(j)*( u2(j  ,column) &
&                                         -u2(j-1,column))
    end do
  end do

end subroutine

subroutine laplacian_U(Lu,u,Luy,iband,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    LAPLACIAN   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Second-order Central Difference Scheme laplacian using 
!   a 3-point stencil (Ferziger pp. 50).
!   Notes:
!     dy2i = 2/dy2; dy2(1,j) = [y(j  )-y(j-1)]*[y(j+1)-y(j-1)]  Upwind
!                   dy2(2,j) = [y(j+1)-y(j  )]*[y(j  )-y(j-1)]  Centered
!                   dy2(3,j) = [y(j+1)-y(j  )]*[y(j+1)-y(j-1)]  Backward

  use declaration
  implicit none

  integer i,k,j,iband,column,myid
  complex(8) Lu (jlim(1,ugrid,iband):jlim(2,ugrid,iband),columns_num(iband,myid))
  complex(8) Luy(jlim(1,ugrid,iband):jlim(2,ugrid,iband),columns_num(iband,myid))
  complex(8)  u(jlim(1,ugrid,iband):jlim(2,ugrid,iband),columns_num(iband,myid))
  real(8) k2x,k2z

  do column = 1,columns_num(iband,myid)
    i = columns_i(column,iband,myid)
    k = columns_k(column,iband,myid)
    k2x = k2F_x(i)
    k2z = k2F_z(k)
    do j = jlim(1,ugrid,iband)+1,jlim(2,ugrid,iband)-1
      Lu(j,column) =         dyu2i(1,j) *u(j-1,column) &
&                  +(k2x+k2z+dyu2i(2,j))*u(j  ,column) &
&                  +         dyu2i(3,j) *u(j+1,column)
     Luy(j,column) =         dyu2i(1,j) *u(j-1,column) &
&                  +(        dyu2i(2,j))*u(j  ,column) &
&                  +         dyu2i(3,j) *u(j+1,column)
    end do
  end do

end subroutine

subroutine laplacian_V(Lu,u,Luy,iband,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    LAPLACIAN   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Second-order Central Difference Scheme laplacian using 
!   a 3-point stencil (Ferziger pp. 50).
!   Notes:
!     dy2i = 2/dy2; dy2(1,j) = [y(j  )-y(j-1)]*[y(j+1)-y(j-1)]  Upwind
!                   dy2(2,j) = [y(j+1)-y(j  )]*[y(j  )-y(j-1)]  Centered
!                   dy2(3,j) = [y(j+1)-y(j  )]*[y(j+1)-y(j-1)]  Backward

  use declaration
  implicit none

  integer i,k,j,iband,column,myid
  complex(8) Lu (jlim(1,vgrid,iband):jlim(2,vgrid,iband),columns_num(iband,myid))
  complex(8) Luy(jlim(1,vgrid,iband):jlim(2,vgrid,iband),columns_num(iband,myid))
  complex(8)  u(jlim(1,vgrid,iband):jlim(2,vgrid,iband),columns_num(iband,myid))
  real(8) k2x,k2z

  do column = 1,columns_num(iband,myid)
    i = columns_i(column,iband,myid)
    k = columns_k(column,iband,myid)
    k2x = k2F_x(i)
    k2z = k2F_z(k)
    do j = jlim(1,vgrid,iband)+1,jlim(2,vgrid,iband)-1
      Lu(j,column) =         dyv2i(1,j) *u(j-1,column) &
&                  +(k2x+k2z+dyv2i(2,j))*u(j  ,column) &
&                  +         dyv2i(3,j) *u(j+1,column)
     Luy(j,column) =         dyv2i(1,j) *u(j-1,column) &
&                  +(        dyv2i(2,j))*u(j  ,column) &
&                  +         dyv2i(3,j) *u(j+1,column)
    end do
  end do

end subroutine

subroutine error(A,myid,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!     ERROR      !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer ierr

  integer j,iband,column,myid
  type(cfield) A(sband:eband)
  real(8) erri,errband

  erri = 0d0
  do iband = sband,eband
    errband = 0d0
    do column = 1,columns_num(iband,myid)
      do j = jlim(1,pgrid,iband),jlim(2,pgrid,iband)
        err     = abs(A(iband)%f(j,column))
        errband = max(errband,err)
      end do
    end do
    erri = max(erri,errband)
  end do

  call MPI_ALLREDUCE(erri,err,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

end subroutine

subroutine finalize(u1,u2,u3,p,div,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    FINALIZE    !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer iproc,iband,j
  real(8) val
  type(cfield)  u1(sband:eband), u2(sband:eband), u3(sband:eband)
  type(cfield)   p(sband:eband),div(sband:eband)

  iwrite = iwrite+nwrite
  call error(div,ierr)

  u1PL = 0d0
  u2PL = 0d0
  u3PL = 0d0
  ppPL = 0d0
  call modes_to_planes_UVP(u1PL,u1,2,myid,status,ierr)
  call modes_to_planes_UVP(u2PL,u2,1,myid,status,ierr)
  call modes_to_planes_UVP(u3PL,u3,2,myid,status,ierr)
  call modes_to_planes_UVP(ppPL, p,3,myid,status,ierr)
 
  iband = bandPL(myid)

  call record_out(u1,myid)

  if (myid/=0) then
    call MPI_SEND(1d0,1,MPI_REAL8,0,101+myid,MPI_COMM_WORLD,ierr)
  else
    do iproc=1,np-1
      call MPI_RECV(val,1,MPI_REAL8,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
      write(*,*) 'process',iproc,'over'
      call MPI_SEND(1d0,1,MPI_REAL8,iproc,100,MPI_COMM_WORLD,ierr)
    end do
  end if
  if (myid/=0) then
    call MPI_RECV(val,1,MPI_REAL8,0,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
  else
     write(*,*) 'process',0,'over'
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (myid==0) then
    write(*,*) ''
    write(*,*) '---------------LITTLE HARSH---------------'
    write(*,*) "         'That wasn't that harsh'         "
    write(*,*) '------------------------------------------'
    write(*,*) ''
  end if
  call MPI_FINALIZE(ierr)

end subroutine

subroutine flowrateIm(Qu,u)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!         FLOW RATE          !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Flow rate: integrating 1st Fourier mode over the whole channel

  use declaration
  implicit none
  integer j
  real(8) Qu
  complex(8) u(N(4,0):N(4,nband)+1)

  Qu = 0d0
  do j = 1,nn+1
    !Qu = Qu + .5d0*(yu(j+1) - yu(j-1))*dreal(u(j) !Old collocated version
    Qu = Qu + (yv(j) - yv(j-1))*dreal(u(j)) !Mixed ugrid vgrid
  end do

end subroutine

subroutine flowrateRe(Qu,u)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!        FLOW RATE 0         !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer j
  real(8) Qu
  real(8) u(1:nn+2)

  Qu = 0d0
  do j = 1,nn+1
    !Qu = Qu + .5d0*(yu(j+1)-yu(j-1))*u(j)
    Qu = Qu + (yv(j) - yv(j-1))*u(j) !Mixed ugrid vgrid
  end do

end subroutine

subroutine flowrate_corr(u,mpg,g)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!        FLOW RATE 0         !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer j
  real(8) Qx0,Qx1,bdt,g,mpg
  complex(8) u(N(4,0):N(4,nband)+1)

  !!!!!!!!!   dummy u11:  !!!!!!!!!
  ! Super meeeean!!
  bdt = -dRK(kRK)*dt
  u11(0) = 0d0
  do j = 1,nn+1
    u11(j) = bdt
  end do
  u11(nn+2) = 0d0
  call LUsol0(u11,0,nn+2)
  !!!!!!!!! u11 flowrate: !!!!!!!!!
  call flowrateRe(Qx1,u11)
  !!!!! uncorrected flowrate: !!!!!
  call flowrateIm(Qx0,u  )
  !!!!!  pressure correction: !!!!!
  g   = (QxT-Qx0)/Qx1
  mpg = mpg+g
  !!!!!  velocity correction: !!!!!
  do j = 0,nn+2
    u(j) = u(j)+g*u11(j)
  end do

end subroutine

subroutine maxvel(u)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!          MAX VEL           !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  complex(8) u(N(4,0):N(4,nband)+1)
  integer j

  Umax = 0d0
  do j = 0,nn+2
    Umax = max(Umax,dreal(u(j)))
  end do

end subroutine

subroutine meanflow_ctP(u1,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! MEAN FLOW WITH CONSTANT MPG  !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer      :: myid
  type(cfield) :: u1(sband:eband)

  ! Mean (mode 0,1) is in proc 0 column 1
  if (myid == 0) then 
    call flowrateIm(Qx,u1(midband)%f(jlim(1,ugrid,midband),1)) 
  end if

end subroutine

subroutine meanflow_ctU(u1,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!       MEAN FLOW CORR       !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer i,k,myid
  type(cfield) u1(sband:eband)

  ! Mean (mode 0,1) is in proc 0 column 1
  if (myid == 0) then 
    call flowrate_corr(u1(midband)%f(jlim(1,ugrid,midband),1),mpgx,dgx)
  end if

end subroutine

!subroutine meanpressgrad(mpg,u)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!   MEAN PRESSURE GRADIENT   !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  use declaration
!  implicit none
!  real(8) mpg,du0,du1
!  complex(8) u(N(4,0):N(4,nband)+1)

!  du0= dreal(u(1) )/(yu(1   )-yu(0 ))
!  du1=-dreal(u(nn))/(yu(nn+2)-yu(nn+1))

!  mpg=(du1-du0)/(Ly*Re)

!end subroutine

subroutine solveP(p,psi,u1,u2,u3,div,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!     SOLVE P    !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer :: i,k,j,iband,column,myid
  type(cfield) ::  u1(sband:eband),u2  (sband:eband),u3 (sband:eband)
  type(cfield) ::  p (sband:eband),psi (sband:eband)
  type(cfield) :: div(sband:eband)
  real(8) :: dtirk
 
  dtirk = dti/(aRK(kRK)+bRK(kRK))
 
  do iband = sband,eband

    call divergence(div(iband)%f,u1(iband)%f,u2(iband)%f,u3(iband)%f,iband,myid)

!    div_out(iband)%f = div(iband)%f  

     
    do column = 1,columns_num(iband,myid)
      do j = jlim(1,pgrid,iband),jlim(2,pgrid,iband)
        psi(iband)%f(j,column) = dtirk*div(iband)%f(j,column) 
      end do
    end do 
    
    !Boundary condition for pressure
    if (myid==0) then
      if (iband==midband) then
        psi(iband)%f(jlim(1,pgrid,iband),1) = 0d0 !C! Psi (mean) at bottom of channel = 0
      end if
    end if

!! For modified wavenumbers, need BC on last pressure mode (as k2x=k2z=0)
! do column = 1,columns_num(iband,myid)
! i = columns_i(column,iband,myid)
! k = columns_k(column,iband,myid)
! if(abs(k1F_x(i)*k1F_x(i))<10e-10.and.abs(k1F_z(k)*k1F_z(k))<10e-10)then
! psi(iband)%f(jlim(1,pgrid,iband),column) = 0d0
! endif
! enddo
    

    call LUsolP(psi(iband)%f,myid,iband,jlim(1,pgrid,iband),jlim(2,pgrid,iband))

    do column = 1,columns_num(iband,myid) 
      do j = jlim(1,pgrid,iband),jlim(2,pgrid,iband)
        p(iband)%f(j,column) = p(iband)%f(j,column)+psi(iband)%f(j,column)
      enddo
    enddo
   
  end do

end subroutine

subroutine RHS0_u1(du1,u1,Nu1,p,Lu1,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!     RHS U1     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer i,k,j,iband,column,myid
  type(cfield)  u1(sband:eband)
  type(cfield) Lu1(sband:eband)
  type(cfield) du1(sband:eband)
  type(cfield) Nu1(sband:eband)
  type(cfield)   p(sband:eband)
  real(8) C1,C2
  complex(8) C3

  !C1 = (aRK(kRK)+bRK(kRK))/Re !For solving for du
  C1 = aRK(kRK)/Re !For solving for u
  C2 = -cRK(kRK)

  do iband = sband,eband
    call laplacian_U(du1(iband)%f,u1(iband)%f,Lu1(iband)%f,iband,myid)
    do column = 1,columns_num(iband,myid)
      i = columns_i(column,iband,myid)
      C3 = -dRK(kRK)*k1F_x(i)
      do j = jlim(1,ugrid,iband)+1,jlim(2,ugrid,iband)-1
        du1(iband)%f(j,column) = C1*du1(iband)%f(j,column) &
&                              + C2*Nu1(iband)%f(j,column) &
&                              + C3* p (iband)%f(j,column)
!        Lu1(iband)%f(j,column) = du1(iband)%f(j,column)
      end do
    end do
  end do

  !Apply mean pressure gradient int x
  ! mode (0,1)
  ! TODO TRICK: We know that mode(0,1) is the first column of proc 0
  if (myid==0) then
    C1 = -dRK(kRK)*mpgx
    do j=N(4,0)+1,N(4,nband)
    !do j=1,nn+1              !       In case you want not to apply pressure gradient in the immersed boundary
      du1(midband)%f(j,1) = du1(midband)%f(j,1) + C1
    end do
  end if


end subroutine

subroutine RHS0_u2(du2,u2,Nu2,p,Lu2,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!     RHS  U2    !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer i,k,j,iband,column,myid
  type(cfield)  u2(sband:eband)
  type(cfield) Lu2(sband:eband)
  type(cfield) du2(sband:eband)
  type(cfield) Nu2(sband:eband)
  type(cfield)   p(sband:eband)
  real(8) C1,C2,C3

  !C1 = (aRK(kRK)+bRK(kRK))/Re !For solving for du
  C1 = aRK(kRK)/Re !For solving for u
  C2 = -cRK(kRK)

  do iband = sband,eband
    call laplacian_V(du2(iband)%f,u2(iband)%f,Lu2(iband)%f,iband,myid)
    do column = 1,columns_num(iband,myid)
      do j = jlim(1,vgrid,iband)+1,jlim(2,vgrid,iband)-1
        C3 = -dRK(kRK)*ddthetavi*dthdyv(j)
        du2(iband)%f(j,column) = C1* du2(iband)%f(j  ,column) &
&                              + C2* Nu2(iband)%f(j  ,column) &
&                              + C3*( p (iband)%f(j+1,column) & ! centeres to faces --> (j+1)-(j)
&                                    -p (iband)%f(j,  column))
!        Lu2(iband)%f(j,column) =  du2(iband)%f(j,  column)
      end do
    end do
  end do
  
end subroutine

subroutine RHS0_u3(du3,u3,Nu3,p,Lu3,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!     RHS  U3    !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer i,k,j,iband,column,myid
  type(cfield)  u3(sband:eband)
  type(cfield) Lu3(sband:eband)
  type(cfield) du3(sband:eband)
  type(cfield) Nu3(sband:eband)
  type(cfield)   p(sband:eband)
  real(8) C1,C2
  complex(8) C3

  !C1 = (aRK(kRK)+bRK(kRK))/Re !For solving for du
  C1 = aRK(kRK)/Re !For solving for u
  C2 = -cRK(kRK)

  do iband = sband,eband
    call laplacian_U(du3(iband)%f,u3(iband)%f,Lu3(iband)%f,iband,myid)
    do column = 1,columns_num(iband,myid)
      k = columns_k(column,iband,myid)
      C3 = -dRK(kRK)*k1F_z(k)
      do j = jlim(1,ugrid,iband)+1,jlim(2,ugrid,iband)-1
        du3(iband)%f(j,column) = C1*du3(iband)%f(j,column) &
&                              + C2*Nu3(iband)%f(j,column) &
&                              + C3* p (iband)%f(j,column)
!       Lu3(iband)%f(j,column) = C1*du3(iband)%f(j,column)
      end do
    end do
  end do

end subroutine

subroutine solveU(u,du,Lu,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!    SOLVE U1    !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer i,k,j,iband,column,myid
  type(cfield)  u(sband:eband)
  type(cfield) Lu(sband:eband)
  type(cfield) du(sband:eband)

  do iband = sband,eband
    do column = 1,columns_num(iband,myid)
      i = columns_i(column,iband,myid)
      k = columns_k(column,iband,myid)
      du(iband)%f(jlim(1,ugrid,iband),column) = 0d0
      do j = jlim(1,ugrid,iband)+1,jlim(2,ugrid,iband)-1
!         du(iband)%f(j,column) = dt*(du(iband)%f(j,column)) !For solving for du
        du(iband)%f(j,column) = u(iband)%f(j,column)+dt*(du(iband)%f(j,column)) !For solving for u
      end do
      du(iband)%f(jlim(2,ugrid,iband),column) = 0d0
    end do
  end do

  call LUsolU(u,du,Lu,ugrid,myid)
!call LUsolV(du,ugrid,myid) !Can use for smooth channel

!  do iband = sband,eband
!    do column = 1,columns_num(iband,myid)
!      do j = jlim(1,ugrid,iband),jlim(2,ugrid,iband)
!        u(iband)%f(j,column) = du(iband)%f(j,column) !For solving for u
!      enddo
!    enddo
!  end do

!if (myid ==0 ) then
!  write(*,*) 'du'
!  do j = jlim(1,ugrid,2),jlim(2,ugrid,2)
!    write(*,*) REAL(du(2)%f(j,1)), AIMAG(du(2)%f(j,1))
!  end do 
!end if 

!if (myid ==0 ) then
!  write(*,*) 'u'
!  do j = jlim(1,ugrid,2),jlim(2,ugrid,2)
!    write(*,*) REAL(u(2)%f(j,1))!, AIMAG(u(2)%f(j,1))
!  end do
!end if

end subroutine

subroutine solveV(u,du,Lu,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!    SOLVE U1    !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  use declaration
  implicit none
  integer i,k,j,iband,column,myid
  type(cfield)  u(sband:eband)
  type(cfield) Lu(sband:eband)
  type(cfield) du(sband:eband)


  do iband = sband,eband
    do column = 1,columns_num(iband,myid)
      i = columns_i(column,iband,myid)
      k = columns_k(column,iband,myid)
      du(iband)%f(jlim(1,vgrid,iband),column) = 0d0
      do j = jlim(1,vgrid,iband)+1,jlim(2,vgrid,iband)-1
!         du(iband)%f(j,column) = dt*(du(iband)%f(j,column)) !For solving for du
        du(iband)%f(j,column) = u(iband)%f(j,column)+dt*(du(iband)%f(j,column)) !For solving for u
      end do
      du(iband)%f(jlim(2,vgrid,iband),column) = 0d0
    end do
  end do

  call LUsolV(u,du,Lu,vgrid,myid)

end subroutine

subroutine v_corr(u1,u2,u3,psi,div,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!     V_CORR     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  
  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,iband,column
  complex(8) kx,kz,dtrk,dtrk_u2
  type(cfield) u1(sband:eband),u2(sband:eband),u3(sband:eband)
  type(cfield) psi(sband:eband),div(sband:eband)
  real (8),pointer :: vcorrPL(:,:,:)
  
  real(8) weighting
  

 
  allocate(vcorrPL(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(ugrid,1)-1:jgal(ugrid,2)+1))
  !!!!!!!!!      u1:      !!!!!!!!!
  dtrk = dt*(aRK(kRK)+bRK(kRK))
 
  do iband = sband,eband
    do column = 1,columns_num(iband,myid)
      i  = columns_i(column,iband,myid)
      kx = dtrk*k1F_x(i)
      
      do j = jlim(1,ugrid,iband)+1,jlim(2,ugrid,iband)-1
!        u_out(iband)%f(j,column) = u1(iband)%f(j,column)
       
        u1(iband)%f(j,column) = u1(iband)%f(j,column)-kx*psi(iband)%f(j,column)
      end do

    end do
  end do
  
  !Enforce BC's
  do iband = 1,2
    j = jlim(1,ugrid,iband)
    do column = 1,columns_num(iband,myid)
      u1(iband)%f(j,column) = -gridweighting(1,1)*u1(iband)%f(j+1,column)
     enddo
  enddo
  
  do iband = 2,3
    j = jlim(2,ugrid,iband)
    do column = 1,columns_num(iband,myid)
      u1(iband)%f(j,column) = -gridweighting(3,2)*u1(iband)%f(j-1,column)
    enddo
  enddo

   
!!!!!!!!!      u2:      !!!!!!!!!
  dtrk_u2 = dtrk*ddthetavi
  do iband = sband,eband
    do column = 1,columns_num(iband,myid)
      
      do j = jlim(1,vgrid,iband)+1,jlim(2,vgrid,iband)-1
!        v_out(iband)%f(j,column) = u2(iband)%f(j,column)

        u2(iband)%f(j,column) = u2(iband)%f(j,column)-dtrk_u2*dthdyv(j)*(psi(iband)%f(j+1,column)-psi(iband)%f(j,column))
      end do
!       u2(iband)%f(jlim(1,vgrid,iband),column) = 0d0
!       u2(iband)%f(jlim(2,vgrid,iband),column) = 0d0

 u2(iband)%f(jlim(1,vgrid,iband),column) = 0d0
 u2(iband)%f(jlim(2,vgrid,iband),column) = 0d0

    end do
  end do

  !!!!!!!!!      u3:      !!!!!!!!!
  do iband = sband,eband
    do column = 1,columns_num(iband,myid)
      k  = columns_k(column,iband,myid)
      kz = dtrk*k1F_z(k)
      
      do j=jlim(1,ugrid,iband)+1,jlim(2,ugrid,iband)-1
!        w_out(iband)%f(j,column) = u3(iband)%f(j,column)

        u3(iband)%f(j,column) = u3(iband)%f(j,column)-kz*psi(iband)%f(j,column)
      end do

    end do
  end do


  !Enforce BC's
  do iband = 1,2
    j = jlim(1,ugrid,iband)
    do column = 1,columns_num(iband,myid)
      u3(iband)%f(j,column) = -gridweighting(1,1)*u3(iband)%f(j+1,column)
     enddo
  enddo
  
  do iband = 2,3
    j = jlim(2,ugrid,iband)
    do column = 1,columns_num(iband,myid)
      u3(iband)%f(j,column) = -gridweighting(3,2)*u3(iband)%f(j-1,column)
    enddo
  enddo
  
  
  !!!!!!!!!  divergence:  !!!!!!!!!
  do iband = sband,eband
    call divergence(div(iband)%f,u1(iband)%f,u2(iband)%f,u3(iband)%f,iband,myid)
  end do
deallocate(vcorrPL)
end subroutine

subroutine boundary_fibers
!-----------------------------------------------------------------------------------------
! !!! Limits of the boundaries at ugrid 

! ATTENTION, it raises an error if the number of points in x (Ngal(1,1)) and z (Ngal(2,1))
!   is not a multiple of the number of tiles (ntilex, ntilez)

! It is called by 'getbounds' in 'start.f90'


!   ------------------------ ny22
!   \  /\  /\  /\  /\  /\  /
!    \/  \/  \/  \/  \/  \/
!   ------------------------ ny12
!
!
!
!   ------------------------ ny21
!    /\  /\  /\  /\  /\  /\
!   /  \/  \/  \/  \/  \/  \
!   ------------------------ ny11
!------------------------------------------------------------------------------------------
! This subroutine creates the geometry: fibres aligned in streamwise direciton
! It defines and locates immersed boundary points(points used in IB calc) in physical space:
!  Stored in lists:
!       Solid points         -- deep inside the solid
!       Ghost points         -- where expected velocity is imposed
!       Surrounding points   -- four(p1, p2, p3, p4) neighour points of the probe point,
!                                only the index of p1 is stored
!       Weighting coefs      -- four, in a separated list
!  Not stored in lists(will not go into flow solver) :
!       Corresponding points -- each one corresponds to a ghost point, 
!                                they are symmetric with respect to boundary
!       Boundary points      -- Aligned with virtual boundary, in the centre 
!                                between Surr and Corr points, where v=0
! See for details: R. Mittal et al., JCP, 2008; J. Picot, S.Glockner, JCP, 2018 
!
! Procedure:
!  1. One element: Define points and calculate weightings for one fibre in zy plane, creat
!     base profile.
!  2. One column: Copy the base profile to each layer (in y direction), 
!     creat lists (in y and x directions).
!  3. Whole substrate: Repulicate the lists column by column (in z direction).
!
! Call:
!  interp -- input the shape function of the solid boundary, coords of ghost  points,
!             calculate the shortest distance from one ghost to the boundary,
!             locate the corresponding, surrounding and boundary points,
!             calculate the weighting coefficients,
!
!                                                                              25.2.2020  
!------------------------------------------------------------------------------------------       


  use declaration
  implicit none
  
  integer :: i,k,j,i1,j1,k1,ix,iz,dny,shift,grid
  integer :: ilist, glist ! ilist: solid points glist: forcing points 10.5.2020
  integer :: ilistu_yz, ilistv_yz, glistu_yz, glistv_yz ! 19.2.2020
  integer :: inxyz_u, inxyz_v, gnxyz_u, gnxyz_v ! 19.2.2020
  integer :: j_ele, k_ele, j_layer
  integer :: k_left, k_right
  integer :: nlist_ibv, nlist_ibu, nlist_gv, nlist_gu ! 4.2.2020
  integer :: ncolumnv, ncolumnu, gcolumnv, gcolumnu   ! 4.2.2020
  integer :: nxyz_solid, nxyz_solid1, points_lol, shift_x, shift_z
  integer :: rest, nyz_solid, nsolid1
  real(8) :: sizesolid, maxheight, minheight
  ! Lists of solid points
  integer, allocatable :: list_ibu1(:,:), list_ibu2(:,:)
  integer, allocatable :: list_ibv1(:,:), list_ibv2(:,:)
  ! Lists of forcing points 20.2.2020
  integer, allocatable :: list_gu1(:,:), list_gu2(:,:)
  integer, allocatable :: list_gv1(:,:), list_gv2(:,:)
  ! Lists of weighting coefficients 24.2.2020
  real(8), allocatable :: list_wu1(:,:), list_wu2(:,:)
  real(8), allocatable :: list_wv1(:,:), list_wv2(:,:)
  ! Lists of boundary points and corresponding points 
  ! Only for geometry, not taking into flow solver 24.2.2020
  real(8), allocatable :: list_geomu1(:,:), list_geomu2(:,:)
  real(8), allocatable :: list_geomv1(:,:), list_geomv2(:,:)
! Lists of limits 20.2.2020
  integer, allocatable :: nyu1_list(:,:), nyu2_list(:,:)
  integer, allocatable :: nyv1_list(:,:), nyv2_list(:,:)
! The base profile in yz plane 24.2.2020
  real(8), pointer :: baseu(:,:,:), basev(:,:,:)   ! Base arrays of distance 
  real(8), pointer :: worku(:,:,:), workv(:,:,:)   ! Working arrays of distance 
! Geometry function
  real(8) :: f
!  real(8) :: radius, dc, dcy, dcz   ! Distance to the centre of circular section,
                                    !  defined in declaration.f90    
  logical :: ll, lr, lt, lb ! four neighbour points of the ghost point
  integer :: intp1(4)
  real(8) :: weicoef(3)
  real(8) :: zbound, ybound, zcorr, ycorr ! coordinates of boundary and corresponding points
  integer :: btol ! Shift from base profile to each layer, should be an integer 3.3.2020

  
  ! Check that the number of grid points in x and z is a multiple of the number of tiles 
  if (ntilex.eq.0) then
    write(*,*) 'ERROR: ntilex equals 0'
    stop
  end if
  if (ntilez.eq.0) then
    write(*,*) 'ERROR: ntilez equals 0'
    stop
  end if
  if (mod(Ngal(1,1),ntilex)/=0) then
    write(*,*) 'ERROR: nx not multiple of ntilex'
    stop
  end if
  if (mod(Ngal(2,1),ntilez)/=0) then
    write(*,*) 'ERROR: nz not multiple of ntilez'
    stop
  end if
  
  dnx = Ngal(1,1)/ntilex  ! Width:  Number of points per tile the in streamwise direction
  dnz = Ngal(2,1)/ntilez  ! Width:  Number of points per tile the in spanwise   direction
  dny = dsty + qp 
  write(*,*) 'dny, dsty, qp ', dny, dsty, qp
  
  shift_x = 0!0.5d0*( dstx)
  shift_z = 2 !shift away from left bound !0 !dstz-1 !0.5d0*( dstz) ! 23.1.2020
  
  nyu11 = -dny+1
  nyu21 = 0
  
  nyv11 = -dny
  nyv21 = 0-1            !-1 To make the boundary align with u_grid

  nyu12 = Ngal(4,3)-dny+1
  nyu22 = Ngal(4,3)

  nyv12 = Ngal(3,3)+1-dny+1  !+1 To make the boundary align with u_grid
  nyv22 = Ngal(3,3)+1

! Above parameters are reserved to confine the circle in the square
! Should have more straightfoward ways 5.12.2019

  !print*, "NGAL", Ngal

  print*, "nyv11, nyv21, nyv12, nyv22:", nyv11, nyv21, nyv12, nyv22
  print*, "nyu11, nyu21, nyu12, nyu22:", nyu11, nyu21, nyu12, nyu22

  if ( ((nyv21-nyv11) /= (nyu21-nyu11)) .or. ((nyv22-nyv12) /= (nyu22-nyu12))) then
	write(*,*) 'Attention, nlist_ib(ugrid) must be equal to nlist_ib(vgrid)'
  end if

!----------------------------------------------------------------------
! CREATE THE PATTERN. 
! if (geometry_flag = 4) then
!    Permeable substrate: collocated circle-cross-scetion fibres
!           
!           O O O O
!	    O O O O
!
!    For staggered fibres, simply change 'shift_z' above
!
!           O O O O
!          O O O O O
!----------------------------------------------------------------------
! The solid boundary (circle) is defined with u-gird points,
!  applies to v-grid mesh. Hence two sets of lists of immersed boundary
!  points are required for u- and v-gird respectively.
!----------------------------------------------------------------------
  !--------------------------------------------------------------------
  ! Calculate nb of solid points in y in the first solid element
  !--------------------------------------------------------------------
  ! Index lists coefficient, defined as the ratio of circular area and square area 18.11.2019
  ! Cancelled. No stretching in the fibrous coating. 9.12.2019
!  acoef = pi / 4d0
  !
!  sizesolid = nsolid * dyp
!  write(*,*) 'Size of the solid blocks = ', sizesolid 
!  do i = -1, -qp, -1
!	maxheight = abs(yu(i) - yu(0))
!	if (maxheight > sizesolid) then
!		if (abs(maxheight-sizesolid) < abs(minheight-sizesolid)) then
!			nsolid1 = -i !nb of points in the first solid element in vgrid
!			write(*,*) 'Size of the first solid block: ', maxheight
!		else
!			nsolid1 = -i + 1
!			write(*,*) 'Size of the first solid block: ', minheight
!		end if
!		exit
!	end if
!	minheight = abs(yu(i) - yu(0))
!  end do
!  write(*,*) 'nb of points in the first solid element in vgrid', nsolid1
!  write(*,*) 'nb of points in the first fluid region in vgrid', qp - nsolid1

  !--------------------------------------------------------------------
  ! Calculate the number of layers in substrates
  !--------------------------------------------------------------------
  !nb of elements in y:
  nelements_y = int(dsty / (nsolid + nfluid)) + 1 
  write(*,*) 'nelements_y = ', nelements_y
  rest = dsty - (nelements_y - 1) * (nsolid + nfluid)
  write(*,*) 'rest = ', rest

  if (rest == nsolid) then
	nelements_y = nelements_y + 1
	write(*,*) 'Solid obstacles near the top and bottom walls.'

  elseif (rest == 0) then
	write(*,*) 'Fluid near the top and bottom walls.'

  elseif (rest < nsolid) then
	write(*,*) 'WARNING! the solid elements near the walls incomplete. Check dsty'
	nelements_y = nelements_y + 1

  elseif (rest > nsolid) then
	write(*,*) 'WARNING! the number of fluid points near to the walls fewer than through the substrate. Check dsty'
	nelements_y = nelements_y + 1

  end if
  write(*,*) 'nelements_y = ', nelements_y
!  write(*,*) 'ncolumnv, ncolumnu = ', ncolumnv, ncolumnu
!  write(*,*) 'nlist_ibv', nlist_ibv
!  write(*,*) 'nlist_ibu', nlist_ibu

!----------------------------------------------------------------------
! Define the function of circular section in j-k plane
! In u-grid (as grid points in j and k always equal?)
!
!    j
!    ^
!    |      ___          1+nsolid
!    |   _-     -_
!    |  -         -
!    | *     *     *     1+nsolid/2
!    |  _         _ 
!    |   -       -
! 1  |      ---          1
! 0 -|-----------------> k
!    0 1  centre  dstz
!
! centre point (1+nsolid/2, 1+nsolid/2)
  radius = nsolid/2d0
!
! for nsolid = even .OR. odd
!
! (j - (nsolid/2 + 1))^2 + (k - (nsolid/2 + 1))^2 - (nsolid/2)^2 = 0
!
!----------------------------------------------------------------------
! Base profile
! Calculate nb of ghost and solid points in one fibre in zy-plane
! 
! u-grid---------------------------------------------------------------

  allocate(baseu(10,0:dstz+1,0:nsolid+1+1)) ! Remember to deallocate 17.2.2020
  
! Sweeping the zy plane of one element, calculating distance to centre.
! Storing the porfile into the base array.
! A row of cells is added across each domain boundary to handle the domain 
!  boundary condition.
  do j_ele = 0, nsolid+1+1 
     do k_ele = 0, dstz+1

!        dcy = abs(j_ele - (1d0+nsolid/2d0))
!        dcz = abs(k_ele - (1d0+nsolid/2d0))
        baseu(1 ,k_ele,j_ele) = f(k_ele*1d0,j_ele*1d0,radius)!sqrt(dcy**2 + dcz**2)

     end do
  end do

! Assigning points into lists
  ilistu_yz = 0
  glistu_yz = 0

  do j_ele = 1, nsolid+1
     do k_ele = 1, dstz
        
        ll = baseu(1, k_ele-1,j_ele  ) <= 0d0!radius
        lr = baseu(1, k_ele+1,j_ele  ) <= 0d0!radius
        lt = baseu(1, k_ele  ,j_ele+1) <= 0d0!radius
        lb = baseu(1, k_ele  ,j_ele-1) <= 0d0!radius

        if (baseu(1, k_ele,j_ele) .gt. 0d0) then
           if (ll .or. lr .or. lt .or. lb) then
              ! This is a immersed boundary point
              glistu_yz = glistu_yz + 1
              ! Calculating IBP's corresponding points on the solid boundary
              ! Calculating the weighting coefficients for interpolation
              zcorr = 0d0
              ycorr = 0d0
              zbound = 0d0
              ybound = 0d0
              intp1 = 0d0
              weicoef = 0d0
              call interp(k_ele*1d0,j_ele*1d0,intp1,weicoef,zcorr,ycorr,zbound,ybound,ugrid)
              ! Store them into the base profile
              baseu(2 ,k_ele,j_ele) = intp1(1)   ! k of fluid point 2
              baseu(3 ,k_ele,j_ele) = intp1(2)   ! j of fluid point 2
              baseu(4 ,k_ele,j_ele) = intp1(3)   ! k of fluid point 3
              baseu(5 ,k_ele,j_ele) = intp1(4)   ! j of fluid point 3
              baseu(6 ,k_ele,j_ele) = weicoef(1) ! weighting of boundary point 1
              baseu(7 ,k_ele,j_ele) = weicoef(2) ! weighting of fluid point 2
              baseu(8 ,k_ele,j_ele) = weicoef(3) ! weighting of fluid point 3
              baseu(9 ,k_ele,j_ele) = zbound     ! z coord of boundary point 1
              baseu(10,k_ele,j_ele) = ybound     ! y coord of boundary point 1
!write(*,*) 'forcing k, j =', k_ele, j_ele
!write(*,*) 'fluid 2, 3   =', intp1(1),intp1(2),intp1(3),intp1(4)
!write(*,*) 'weighting    =', weicoef(1),weicoef(2),weicoef(3)
           end if
        else 
           ! This is a solid point
           ilistu_yz = ilistu_yz + 1
        end if
     end do
  end do


! v-grid---------------------------------------------------------------

  allocate(basev(10,0:dstz+1,-1:nsolid+1+1)) ! Remember to deallocate 19.2.2020

! Sweeping the zy plane of one element, calculating distance to centre.
! Storing the profile into the base array.
! Two rows of cells is added across each domain boundary to handle the domain 
!  boundary condition.
  do j_ele = -1, nsolid+1+1
     do k_ele = 0, dstz+1

!        dcy = abs(j_ele + 0.5d0 - (1d0+nsolid/2d0))
!        dcz = abs(k_ele - (1d0+nsolid/2d0))
        basev(1, k_ele,j_ele) = f(k_ele*1d0,j_ele+0.5d0,radius)!sqrt(dcy**2 + dcz**2)

     end do
  end do

! Assigning points into lists
  ilistv_yz = 0
  glistv_yz = 0

  do j_ele = 0, nsolid+1
     do k_ele = 1, dstz

        ll = basev(1, k_ele-1,j_ele  ) <= 0d0!radius
        lr = basev(1, k_ele+1,j_ele  ) <= 0d0!radius
        lt = basev(1, k_ele  ,j_ele+1) <= 0d0!radius
        lb = basev(1, k_ele  ,j_ele-1) <= 0d0!radius

        if (basev(1, k_ele,j_ele) .gt. 0d0) then
           if (ll .or. lr .or. lt .or. lb) then
              glistv_yz = glistv_yz + 1      
              ! Calculating ghost's corresponding points in flow domain 
              ! Calculating the weighting coefficients for interpolation
              zcorr = 0d0
              ycorr = 0d0
              zbound = 0d0
              ybound = 0d0
              intp1 = 0d0
              weicoef = 0d0
              call interp(k_ele*1d0,j_ele*1d0,intp1,weicoef,zcorr,ycorr,zbound,ybound,vgrid)
              ! Store them into the base profile
              basev(2 ,k_ele,j_ele) = intp1(1)   ! k of fluid point 2
              basev(3 ,k_ele,j_ele) = intp1(2)   ! j of fluid point 2
              basev(4 ,k_ele,j_ele) = intp1(3)   ! k of fluid point 3
              basev(5 ,k_ele,j_ele) = intp1(4)   ! j of fluid point 3
              basev(6 ,k_ele,j_ele) = weicoef(1) ! weighting of boundary point 1
              basev(7 ,k_ele,j_ele) = weicoef(2) ! weighting of fluid point 2
              basev(8 ,k_ele,j_ele) = weicoef(3) ! weighting of fluid point 3
              basev(9 ,k_ele,j_ele) = zbound     ! z coord of boundary point 1
              basev(10,k_ele,j_ele) = ybound     ! y coord of boundary point 1
           end if
        else 
           ilistv_yz = ilistv_yz + 1
        end if
     end do
  end do

!----------------------------------------------------------------------
! Calculate nb of IB points for the first column of fibres 19.2.2020
!----------------------------------------------------------------------
  ! These nbs of points are for ghost and solid points 
  !  in u- and v-grid respectively.

  inxyz_u = dstx * ilistu_yz   ! Solid points in 1 fibre in u-grid
  inxyz_v = dstx * ilistv_yz   ! Solid points in 1 fibre in v-grid
  gnxyz_u = dstx * glistu_yz   ! Ghost points in 1 fibre in u-grid
  gnxyz_v = dstx * glistv_yz   ! Ghost points in 1 fibre in v-grid
write(*,*) 'gnxyz_v =', gnxyz_v

  if (rest == nsolid) then ! To be modified 19.2.2020
!        nelements_y = nelements_y + 1
        ncolumnv = (nelements_y-1) * nxyz_solid + nxyz_solid1  ! nb of IB pointsin the one elements column in vgrid
        ncolumnu = (nelements_y-1) * (nxyz_solid + nyz_solid) - nyz_solid &
                &+ (nxyz_solid1 + nyz_solid)                         ! in ugrid 
        nlist_ibv  = ncolumnv * ntilex * ntilez  ! Number of IB points in allelements in vgrid.
        nlist_ibu  = ncolumnu * ntilex * ntilez  ! Number of IB points in allelements in ugrid. 
  elseif (rest == 0) then  ! Currently working
        ncolumnv = nelements_y * inxyz_v ! nb of solid points in the first element column in vgrid
        ncolumnu = nelements_y * inxyz_u                                                ! in ugrid
        gcolumnv = nelements_y * gnxyz_v ! nb of ghost points in the first element column in vgrid 
        gcolumnu = nelements_y * gnxyz_u                                                ! in ugrid 
        nlist_ibv  = ncolumnv * ntilex * ntilez  ! Number of solid points in all elements in vgrid.
        nlist_ibu  = ncolumnu * ntilex * ntilez  ! Number of solid points in all elements in ugrid. 
        nlist_gv = gcolumnv * ntilex * ntilez  ! Number of ghost points in all elements in v gird
        nlist_gu = gcolumnu * ntilex * ntilez  ! Number of ghost points in all elements in u gird
  elseif (rest < nsolid) then ! To be modified 19.2.2020
        ncolumnv  = ncolumnv + (nelements_y - 1 ) * nxyz_solid 
        ncolumnu  = ncolumnu + (nelements_y - 1 ) * (nxyz_solid + (dstz-2) * dstx)
        nlist_ibv   = ncolumnv * ntilex * ntilez  ! Number of IB points in allelements in vgrid.
        nlist_ibu   = ncolumnu * ntilex * ntilez  ! Number of IB points in allelements in ugrid. 
  elseif (rest > nsolid) then ! To be modified 19.2.2020
!        write(*,*) 'WARNING! the number of fluid points near to the walls fewerthan through the substrate. Check dsty'
!        nelements_y = nelements_y + 1
        ncolumnv  = (nelements_y-1) * nxyz_solid + nxyz_solid1 ! nb of IB pointsin the first element column in vgrid
        ncolumnu  = (nelements_y-1) * (nxyz_solid + nyz_solid) &
                &+ (nxyz_solid1 + nyz_solid)                         ! in ugrid
        nlist_ibv   = ncolumnv * ntilex * ntilez  ! Number of IB points in allelements in vgrid.
        nlist_ibu   = ncolumnu * ntilex * ntilez  ! Number of IB points in allelements in ugrid. 
  end if
!  write(*,*) 'nelements_y = ', nelements_y
  write(*,*) 'ncolumnv, ncolumnu = ', ncolumnv, ncolumnu
  write(*,*) 'nlist_ibv', nlist_ibv
  write(*,*) 'nlist_ibu', nlist_ibu
  write(*,*) 'gcolumnv, gcolumnu =', gcolumnv, gcolumnu
  write(*,*) 'nlist_gv', nlist_gv
  write(*,*) 'nlist_gu', nlist_gu

!----------------------------------------------------------------------
! Lists of limits to the geometry 
!
!        _____________ 2: upper limit
!     _-     -_
!    -         -
!   *     *     *----- 3: centre limit
!    _         _ 
!     -       -
!        ------------- 1: lower limit
!
! 20.2.2020: Only work for 'rest = 0', to be extended to general cases
!
!----------------------------------------------------------------------
  allocate(nyu1_list(2,nelements_y), nyu2_list(2,nelements_y))
  allocate(nyv1_list(2,nelements_y), nyv2_list(2,nelements_y)) ! Remember to deallocate 20.2.2020
  
  ! Bottom substrate
  nyu1_list = 0d0
  nyv1_list = 0d0        

  do j_layer = 1, nelements_y

     nyu1_list(1,j_layer) = nyu21 - (j_layer - 1) * (nsolid + nfluid) - nsolid
     nyu1_list(2,j_layer) = nyu21 - (j_layer - 1) * (nsolid + nfluid)
!     nyu1_list(3,j_layer) = (nyu1_list(1,j_layer) + nyu1_list(2,j_layer)) / 2d0 

     nyv1_list(1,j_layer) = nyv21 - (j_layer - 1) * (nsolid + nfluid) - nsolid 
     nyv1_list(2,j_layer) = nyv21 - (j_layer - 1) * (nsolid + nfluid) + 1
!     nyv1_list(3,j_layer) = nyu1_list(3,j_layer) ! The circle is defined in u-gird.

  end do

  ! Top substrate
  nyu2_list = 0d0
  nyv2_list = 0d0

  do j_layer = 1, nelements_y

     nyu2_list(1,j_layer) = nyu12 + (j_layer - 1) * (nsolid + nfluid) 
     nyu2_list(2,j_layer) = nyu12 + (j_layer - 1) * (nsolid + nfluid) + nsolid
!     nyu2_list(3,j_layer) = (nyu2_list(1,j_layer) + nyu2_list(2,j_layer)) / 2d0

     nyv2_list(1,j_layer) = nyv12 + (j_layer - 1) * (nsolid + nfluid) - 1
     nyv2_list(2,j_layer) = nyv12 + (j_layer - 1) * (nsolid + nfluid) + nsolid 
!     nyv2_list(3,j_layer) = nyu2_list(3,j_layer) ! The circle is defined in u-grid.

  end do

!------------------------------------------------------------------
! First column of obstacles, bottom boundary
!------------------------------------------------------------------  
  allocate(list_ibu1(3,nlist_ibu), list_ibv1(3,nlist_ibv))
  allocate(list_gu1(9,nlist_gu), list_gv1(9,nlist_gv))
  allocate(list_wu1(3,nlist_gu), list_wv1(3,nlist_gv))
  allocate(list_geomu1(3,nlist_gu), list_geomv1(3,nlist_gv))
!   write(*,*) 'u limits -------------------'
!   write(*,*) 'nyu11_list', nyu11_list
!   write(*,*) 'nyu21_list', nyu21_list
!   write(*,*) 'v limits -------------------'
!   write(*,*) 'nyv11_list', nyv11_list
!   write(*,*) 'nyv21_list', nyv21_list

  ! indexes for first column of elements in vgrid
  ilist = 0
  ! index for ghost points
  glist = 0
  
  do j_layer = nelements_y, 1, -1

     ! Copy the distance profile to the working array
     allocate(workv(10,0:dstz+1, nyv1_list(1,j_layer)-1:nyv1_list(2,j_layer)+1))
     workv = basev

     ! Shift from base profile to each layer (shift downward)
     btol = -(1+nsolid)-(j_layer-1)*(nsolid+nfluid) !-0.5d0

     do j_ele = nyv1_list(1,j_layer), nyv1_list(2,j_layer)
        j1 = j_ele
  
        do k_ele = 1, dstz
           ! Collocated or staggered
!           if(mod(j_layer,2) .eq. 0) then
           k1 = k_ele + shift_z
!           else 
!           k1 = k_ele
!           end if
        
           ll = workv(1, k_ele-1,j_ele  ) <= 0d0!radius
           lr = workv(1, k_ele+1,j_ele  ) <= 0d0!radius
           lt = workv(1, k_ele  ,j_ele+1) <= 0d0!radius
           lb = workv(1, k_ele  ,j_ele-1) <= 0d0!radius     
 
           if (workv(1, k_ele,j_ele) .gt. 0d0) then
              if (ll .or. lr .or. lt .or. lb) then
                 ! Forcing points 
                 do i = 1, dstx
                    i1 = i + shift_x
                    glist = glist + 1
   
                    list_gv1(1,glist) = i1                              ! i of forcing point
                    list_gv1(2,glist) = k1                              ! k of forcing point
                    list_gv1(3,glist) = j1                              ! j of forcing point
                    list_gv1(4,glist) = i1                              ! i of fluid point 2
                    list_gv1(5,glist) = workv(2,k_ele,j_ele) + shift_z  ! k of fluid point 2
                    list_gv1(6,glist) = workv(3,k_ele,j_ele) + btol     ! j of fluid point 2
                    list_gv1(7,glist) = i1                              ! i of fluid point 3
                    list_gv1(8,glist) = workv(4,k_ele,j_ele) + shift_z  ! k of fluid point 3
                    list_gv1(9,glist) = workv(5,k_ele,j_ele) + btol     ! j of fluid point 3

                    list_wv1(1,glist) = workv(6,k_ele,j_ele)            ! weighting of boundary point 1
                    list_wv1(2,glist) = workv(7,k_ele,j_ele)            ! weighting of    fluid point 2
                    list_wv1(3,glist) = workv(8,k_ele,j_ele)            ! weighting of    fluid point 3
                    
                    list_geomv1(1,glist) = i1                                   ! x coordinate of the boundary point
                    list_geomv1(2,glist) = workv(9 ,k_ele,j_ele) + shift_z      ! z coordinate of the boundary point
                    list_geomv1(3,glist) = workv(10,k_ele,j_ele) + btol-0.5d0   ! y coordinate of the boundary point
                 end do
              end if
           else
              ! Solid points
              do i = 1, dstx
                 i1 = i + shift_x
                 ilist = ilist + 1

                 list_ibv1( 1,ilist) = i1     ! i of the obstacle
                 list_ibv1( 2,ilist) = k1     ! k of the obstacle
                 list_ibv1( 3,ilist) = j1     ! j of the obstacle 
              end do
           end if
        end do
     end do

     deallocate(workv)

  end do
   
  write(*,*) 'ilist in vgrid', ilist
  write(*,*) 'glist in vgrid', glist

  if ( ilist /= ncolumnv ) then
      write(*,*) 'ERROR: in vgrid ilist is not equal to ncolumnv'
      stop
  end if 
  if ( glist .ne. gcolumnv ) then
      write(*,*) 'ERROR: in vgrid glist is not equal to gcolumnv'
      stop
  end if

  
  ! indexes for first column of elements in ugrid
  ilist = 0
  ! index for ghost points
  glist = 0

  do j_layer = nelements_y, 1, -1

     ! Copy the distance profile
     allocate(worku(10,0:dstz+1, nyu1_list(1,j_layer)-1:nyu1_list(2,j_layer)+1))
     worku = baseu
     ! Shift from the base profile to each layer (shift downward)
     btol = -(1+nsolid)-(j_layer-1)*(nsolid+nfluid)

     do j_ele = nyu1_list(1,j_layer), nyu1_list(2,j_layer)
        j1 = j_ele
  
        do k_ele = 1, dstz
           ! Collocated or staggered
!           if(mod(j_ele,2) .eq. 0) then
           k1 = k_ele + shift_z
!           else
!           k1 = k_ele
!           end if
        
           ll = worku(1, k_ele-1,j_ele  ) <= 0d0!radius
           lr = worku(1, k_ele+1,j_ele  ) <= 0d0!radius
           lt = worku(1, k_ele  ,j_ele+1) <= 0d0!radius
           lb = worku(1, k_ele  ,j_ele-1) <= 0d0!radius

           if (worku(1, k_ele,j_ele) .gt. 0d0) then
              if (ll .or. lr .or. lt .or. lb) then
                 ! Forcing points
                 do i = 1, dstx
                    i1 = i + shift_x
                    glist = glist + 1

                    list_gu1(1,glist) = i1                              ! i of forcing point
                    list_gu1(2,glist) = k1                              ! k of forcing point
                    list_gu1(3,glist) = j1                              ! j of forcing point
                    list_gu1(4,glist) = i1                              ! i of fluid point 2
                    list_gu1(5,glist) = worku(2,k_ele,j_ele) + shift_z  ! k of fluid point 2
                    list_gu1(6,glist) = worku(3,k_ele,j_ele) + btol     ! j of fluid point 2
                    list_gu1(7,glist) = i1                              ! i of fluid point 3
                    list_gu1(8,glist) = worku(4,k_ele,j_ele) + shift_z  ! k of fluid point 3
                    list_gu1(9,glist) = worku(5,k_ele,j_ele) + btol     ! j of fluid point 3

                    list_wu1(1,glist) = worku(6,k_ele,j_ele)
                    list_wu1(2,glist) = worku(7,k_ele,j_ele)
                    list_wu1(3,glist) = worku(8,k_ele,j_ele)

                    list_geomu1(1,glist) = i1                             ! x coordinate of the boundary point 
                    list_geomu1(2,glist) = worku(9 ,k_ele,j_ele) + shift_z! z coordinate of the boundary point 
                    list_geomu1(3,glist) = worku(10,k_ele,j_ele) + btol   ! y coordinate of the boundary point
                 end do
              end if
           else
              ! Solid points                
              do i = 1, dstx
                 i1 = i + shift_x
                 ilist = ilist + 1

                 list_ibu1( 1,ilist) = i1   ! i coordinate of the obstacle
                 list_ibu1( 2,ilist) = k1   ! k coordinate of the obstacle
                 list_ibu1( 3,ilist) = j1   ! j coordinate of the obstacle 	
              end do
           end if
        end do
     end do

     deallocate(worku)

  end do

  write(*,*) 'ilist in ugrid', ilist
  write(*,*) 'glist in ugrid', glist

  if ( ilist /= ncolumnu ) then
      write(*,*) 'ERROR: in ugrid ilist is not equal to ncolumnu'
      stop
  end if
  if ( glist .ne. gcolumnu ) then
      write(*,*) 'ERROR: in ugrid glist is not equal to gcolumnu'
      stop
  end if
  
!------------------------------------------------------------------
! First column of obstacles top boundary
!------------------------------------------------------------------
  allocate(list_ibu2(3,nlist_ibu), list_ibv2(3,nlist_ibv)) 
  allocate(list_gu2(9,nlist_gu), list_gv2(9,nlist_gv))
  allocate(list_wu2(3,nlist_gu), list_wv2(3,nlist_gv))
  allocate(list_geomu2(3,nlist_gu), list_geomv2(3,nlist_gv))

  ! indexes for first column of elements in vgrid
  ilist = 0
  ! index for ghost points
  glist = 0

  do j_layer = nelements_y, 1, -1

     ! Copy the distance profile to the working array
     allocate(workv(10,0:dstz+1, nyv2_list(1,j_layer)-1:nyv2_list(2,j_layer)+1))
     workv = basev
     ! Shift the base profile to each layer (shift upward)
     btol = +(nyv12-1)+(j_layer-1)*(nsolid+nfluid)

     do j_ele = nyv2_list(1,j_layer), nyv2_list(2,j_layer)
        j1 = j_ele
  
        do k_ele = 1, dstz
           ! Collocated or staggered
!           if(mod(j_ele,2) .eq. 0) then
           k1 = k_ele + shift_z
!           else
!           k1 = k_ele
!           end if

           ll = workv(1 ,k_ele-1,j_ele  ) <= 0d0!radius
           lr = workv(1 ,k_ele+1,j_ele  ) <= 0d0!radius
           lt = workv(1 ,k_ele  ,j_ele+1) <= 0d0!radius
           lb = workv(1 ,k_ele  ,j_ele-1) <= 0d0!radius

           if (workv(1 ,k_ele,j_ele) .gt. 0d0) then
              if (ll .or. lr .or. lt .or. lb) then
                 ! Forcing points
                 do i = 1, dstx
                     i1 = i + shift_x
                     glist = glist + 1

                     list_gv2(1,glist) = i1                             ! i of forcing point
                     list_gv2(2,glist) = k1                             ! k of forcing point
                     list_gv2(3,glist) = j1                             ! j of forcing point
                     list_gv2(4,glist) = i1                             ! i of fluid point 2
                     list_gv2(5,glist) = workv(2,k_ele,j_ele) + shift_z ! k of fluid point 2
                     list_gv2(6,glist) = workv(3,k_ele,j_ele) + btol    ! j of fluid point 2
                     list_gv2(7,glist) = i1                             ! i of fluid point 3
                     list_gv2(8,glist) = workv(4,k_ele,j_ele) + shift_z ! k of fluid point 3
                     list_gv2(9,glist) = workv(5,k_ele,j_ele) + btol    ! j of fluid point 3

                     list_wv2(1,glist) = workv(6,k_ele,j_ele)
                     list_wv2(2,glist) = workv(7,k_ele,j_ele)
                     list_wv2(3,glist) = workv(8,k_ele,j_ele)

                     list_geomv2(1,glist) = i1                                      ! x coordinate of the boundary point 
                     list_geomv2(2,glist) = workv(9 ,k_ele,j_ele) + shift_z         ! z coordinate of the boundary point 
                     list_geomv2(3,glist) = workv(10,k_ele,j_ele) + btol -1+0.5d0   ! y coordinate of the boundary point
                 end do
              end if
           else
              ! Solid points
              do i = 1, dstx
                 i1 = i + shift_x
                 ilist = ilist + 1

                 list_ibv2( 1,ilist) = i1   ! i coordinate of the obstacle
                 list_ibv2( 2,ilist) = k1   ! k coordinate of the obstacle
                 list_ibv2( 3,ilist) = j1   ! j coordinate of the obstacle 
              end do
           end if
        end do
     end do 
     
     deallocate(workv)

  end do

  if ( ilist /= ncolumnv ) then 
      write(*,*) 'ERROR: in vgrid ilist is not equal to ncolumnv'
      stop
  end if  
  if ( glist .ne. gcolumnv ) then
      write(*,*) 'ERROR: in vgrid glist is not equal to gcolumnv'
      stop
  end if

  ! indexes for first column of elements in ugrid
  ilist = 0
  ! index for ghost points
  glist = 0 

  do j_layer = nelements_y, 1, -1

     ! Copy the distance profile to the working array
     allocate(worku(10,0:dstz+1, nyu2_list(1,j_layer)-1:nyu2_list(2,j_layer)+1))
     worku = baseu
     ! Shift the base profile to each layer (shift upward)
     btol = +(nyu12-1)+(j_layer-1)*(nsolid+nfluid)

     do j_ele = nyu2_list(1,j_layer), nyu2_list(2,j_layer)
        j1 = j_ele
  
        do k_ele = 1, dstz
           ! Collocated or staggered
!           if(mod(j_ele,2) .eq. 0) then
           k1 = k_ele + shift_z
!           else
!           k1 = k_ele
!           end if
 
           ll = worku(1 ,k_ele-1,j_ele  ) <= 0d0!radius
           lr = worku(1 ,k_ele+1,j_ele  ) <= 0d0!radius
           lt = worku(1 ,k_ele  ,j_ele+1) <= 0d0!radius
           lb = worku(1 ,k_ele  ,j_ele-1) <= 0d0!radius

           if (worku(1 ,k_ele,j_ele) .gt. 0d0) then
              if (ll .or. lr .or. lt .or. lb) then
                  ! Foring points 
                  do i = 1, dstx
                     i1 = i + shift_x
                     glist = glist + 1

                     list_gu2(1,glist) = i1                             ! i of forcing point
                     list_gu2(2,glist) = k1                             ! k of forcing point
                     list_gu2(3,glist) = j1                             ! j of forcing point
                     list_gu2(4,glist) = i1                             ! i of fluid point 2
                     list_gu2(5,glist) = worku(2,k_ele,j_ele) + shift_z ! k of fluid point 2
                     list_gu2(6,glist) = worku(3,k_ele,j_ele) + btol    ! j of fluid point 2
                     list_gu2(7,glist) = i1                             ! i of fluid point 3
                     list_gu2(8,glist) = worku(4,k_ele,j_ele) + shift_z ! k of fluid point 3
                     list_gu2(9,glist) = worku(5,k_ele,j_ele) + btol    ! j of fluid point 3

                     list_wu2(1,glist) = worku(6,k_ele,j_ele)
                     list_wu2(2,glist) = worku(7,k_ele,j_ele)
                     list_wu2(3,glist) = worku(8,k_ele,j_ele)

                     list_geomu2(1,glist) = i1                             ! x coordinate of the boundary point 
                     list_geomu2(2,glist) = worku(9 ,k_ele,j_ele) + shift_z! z coordinate of the boundary point 
                     list_geomu2(3,glist) = worku(10,k_ele,j_ele) + btol   ! y coordinate of the boundary point
                  end do
               end if
            else
               ! Solid points
               do i = 1, dstx
                  i1 = i + shift_x
                  ilist = ilist + 1

                  list_ibu2( 1,ilist) = i1   ! i coordinate of the obstacle
                  list_ibu2( 2,ilist) = k1   ! k coordinate of the obstacle
                  list_ibu2( 3,ilist) = j1   ! j coordinate of the obstacle 
               end do       
           end if
        end do
     end do
     
     deallocate(worku)

  end do
 
  if ( ilist /= ncolumnu ) then
      write(*,*) 'ERROR: in ugrid ilist is not equal to ncolumnu'
      stop
  end if 
  if ( glist .ne. gcolumnu ) then
      write(*,*) 'ERROR: in ugrid glist is not equal to gcolumnu'
      stop
  end if

!-----------------------------------------------------------------  
! REPLICATE THE PATTERN in x and z
!-----------------------------------------------------------------
!   This section should be common for all geometries 

  
  do ix = 1, ntilex
     do iz = 1, ntilez!-1 ! 9.1.2020
        shift = ncolumnv * ( iz - 1 ) + ncolumnv * ntilez * ( ix - 1 ) 
        do ilist = 1, ncolumnv
           list_ibv1(1,ilist+shift) = list_ibv1(1,ilist) + dnx*(ix - 1)   ! i coordinate
           list_ibv1(2,ilist+shift) = list_ibv1(2,ilist) + dnz*(iz - 1)   ! k coordinate
           list_ibv1(3,ilist+shift) = list_ibv1(3,ilist)                  ! j coordinate

           list_ibv2(1,ilist+shift) = list_ibv2(1,ilist) + dnx*(ix - 1)   ! i coordinate
           list_ibv2(2,ilist+shift) = list_ibv2(2,ilist) + dnz*(iz - 1)   ! k coordinate
           list_ibv2(3,ilist+shift) = list_ibv2(3,ilist)                  ! j coordinate
        end do
                
        shift = gcolumnv * ( iz - 1) + gcolumnv * ntilez * ( ix - 1 )          ! 4.2.2020
        do glist = 1, gcolumnv
           list_gv1(1,glist+shift) = list_gv1(1,glist) + dnx*(ix - 1)     ! i forcing
           list_gv1(2,glist+shift) = list_gv1(2,glist) + dnz*(iz - 1)     ! k forcing
           list_gv1(3,glist+shift) = list_gv1(3,glist)                    ! j forcing
           list_gv1(4,glist+shift) = list_gv1(4,glist) + dnx*(ix - 1)     ! i fluid 2
           list_gv1(5,glist+shift) = list_gv1(5,glist) + dnz*(iz - 1)     ! k fluid 2
           list_gv1(6,glist+shift) = list_gv1(6,glist)                    ! j fluid 2
           list_gv1(7,glist+shift) = list_gv1(7,glist) + dnx*(ix - 1)     ! i fluid 3
           list_gv1(8,glist+shift) = list_gv1(8,glist) + dnz*(iz - 1)     ! k fluid 3
           list_gv1(9,glist+shift) = list_gv1(9,glist)                    ! j fluid 3

           list_wv1(1,glist+shift) = list_wv1(1,glist)                    ! weighting of boundary 1
           list_wv1(2,glist+shift) = list_wv1(2,glist)                    ! weighting of    fluid 2
           list_wv1(3,glist+shift) = list_wv1(3,glist)                    ! weighting of    fluid 3

           list_geomv1(1,glist+shift) = list_geomv1(1,glist) + dnx*(ix - 1)
           list_geomv1(2,glist+shift) = list_geomv1(2,glist) + dnz*(iz - 1)
           list_geomv1(3,glist+shift) = list_geomv1(3,glist)

           list_gv2(1,glist+shift) = list_gv2(1,glist) + dnx*(ix - 1)     ! i forcing
           list_gv2(2,glist+shift) = list_gv2(2,glist) + dnz*(iz - 1)     ! k forcing
           list_gv2(3,glist+shift) = list_gv2(3,glist)                    ! j forcing
           list_gv2(4,glist+shift) = list_gv2(4,glist) + dnx*(ix - 1)     ! i fluid 2
           list_gv2(5,glist+shift) = list_gv2(5,glist) + dnz*(iz - 1)     ! k fluid 2
           list_gv2(6,glist+shift) = list_gv2(6,glist)                    ! j fluid 2
           list_gv2(7,glist+shift) = list_gv2(7,glist) + dnx*(ix - 1)     ! i fluid 3
           list_gv2(8,glist+shift) = list_gv2(8,glist) + dnz*(iz - 1)     ! k fluid 3
           list_gv2(9,glist+shift) = list_gv2(9,glist)                    ! j fluid 3

           list_wv2(1,glist+shift) = list_wv2(1,glist)                    ! weighting of boudnary 1
           list_wv2(2,glist+shift) = list_wv2(2,glist)                    ! weighting of    fluid 2
           list_wv2(3,glist+shift) = list_wv2(3,glist)                    ! weighting of    fluid 3

           list_geomv2(1,glist+shift) = list_geomv2(1,glist) + dnx*(ix - 1)
           list_geomv2(2,glist+shift) = list_geomv2(2,glist) + dnz*(iz - 1)
           list_geomv2(3,glist+shift) = list_geomv2(3,glist)
        end do
     end do
  end do

  do ix = 1, ntilex
     do iz = 1, ntilez!-1 ! 9.1.2020
        shift = ncolumnu * ( iz - 1 ) + ncolumnu * ntilez * ( ix - 1 ) 
        do ilist = 1, ncolumnu
           list_ibu1(1,ilist+shift) = list_ibu1(1,ilist) + dnx*(ix - 1)   ! i coordinate
           list_ibu1(2,ilist+shift) = list_ibu1(2,ilist) + dnz*(iz - 1)   ! k coordinate
           list_ibu1(3,ilist+shift) = list_ibu1(3,ilist)                  ! j coordinate

           list_ibu2(1,ilist+shift) = list_ibu2(1,ilist) + dnx*(ix - 1)   ! i coordinate
           list_ibu2(2,ilist+shift) = list_ibu2(2,ilist) + dnz*(iz - 1)   ! k coordinate
           list_ibu2(3,ilist+shift) = list_ibu2(3,ilist)                  ! j coordinate
        end do

        shift = gcolumnu * ( iz - 1 ) + gcolumnu * ntilez * ( ix - 1 )         ! 4.2.2020
        do glist = 1, gcolumnu
           list_gu1(1,glist+shift) = list_gu1(1,glist) + dnx*(ix - 1)     ! i forcing
           list_gu1(2,glist+shift) = list_gu1(2,glist) + dnz*(iz - 1)     ! k
           list_gu1(3,glist+shift) = list_gu1(3,glist)                    ! j
           list_gu1(4,glist+shift) = list_gu1(4,glist) + dnx*(ix - 1)     ! i fluid 2
           list_gu1(5,glist+shift) = list_gu1(5,glist) + dnz*(iz - 1)     ! k
           list_gu1(6,glist+shift) = list_gu1(6,glist)                    ! j
           list_gu1(7,glist+shift) = list_gu1(7,glist) + dnx*(ix - 1)     ! i fluid 3
           list_gu1(8,glist+shift) = list_gu1(8,glist) + dnz*(iz - 1)     ! k
           list_gu1(9,glist+shift) = list_gu1(9,glist)                    ! j

           list_wu1(1,glist+shift) = list_wu1(1,glist)                    ! weighting
           list_wu1(2,glist+shift) = list_wu1(2,glist)
           list_wu1(3,glist+shift) = list_wu1(3,glist)

           list_geomu1(1,glist+shift) = list_geomu1(1,glist) + dnx*(ix - 1)
           list_geomu1(2,glist+shift) = list_geomu1(2,glist) + dnz*(iz - 1)
           list_geomu1(3,glist+shift) = list_geomu1(3,glist)

           list_gu2(1,glist+shift) = list_gu2(1,glist) + dnx*(ix - 1)     ! i forcing
           list_gu2(2,glist+shift) = list_gu2(2,glist) + dnz*(iz - 1)     ! k
           list_gu2(3,glist+shift) = list_gu2(3,glist)                    ! j
           list_gu2(4,glist+shift) = list_gu2(4,glist) + dnx*(ix - 1)     ! i fluid 2
           list_gu2(5,glist+shift) = list_gu2(5,glist) + dnz*(iz - 1)     ! k
           list_gu2(6,glist+shift) = list_gu2(6,glist)                    ! j
           list_gu2(7,glist+shift) = list_gu2(7,glist) + dnx*(ix - 1)     ! i fluid 3
           list_gu2(8,glist+shift) = list_gu2(8,glist) + dnz*(iz - 1)     ! k
           list_gu2(9,glist+shift) = list_gu2(9,glist)                    ! j

           list_wu2(1,glist+shift) = list_wu2(1,glist)                    ! weighting
           list_wu2(2,glist+shift) = list_wu2(2,glist)
           list_wu2(3,glist+shift) = list_wu2(3,glist)

           list_geomu2(1,glist+shift) = list_geomu2(1,glist) + dnx*(ix - 1)
           list_geomu2(2,glist+shift) = list_geomu2(2,glist) + dnz*(iz - 1)
           list_geomu2(3,glist+shift) = list_geomu2(3,glist)
        end do
     end do
  end do
  


! Save the lists into a file
  open(10,file=trim(dirout)//'boundary_'//ext1//'x'//ext2//'x'//ext3//'.dat',form='unformatted',access='stream')
  write(10) Lx, Ly, Lz
  write(10) Ngal, nlist_ibu, nlist_ibv, nlist_gu, nlist_gv, ncolumnv, ncolumnu, &  ! 4.2.2020
  & nyu11, nyu21, nyu12, nyu22, nyv11, nyv21, nyv12, nyv22
  write(10) nelements_y, dstz
!, nyu11_list, nyu21_list, nyu12_list, nyu22_list, &
!  & nyv11_list, nyv21_list, nyv12_list, nyv22_list 
  write(10) list_ibu1, list_ibu2
  write(10) list_ibv1, list_ibv2
  write(10) list_gu1, list_gu2 ! 3.2.2020
  write(10) list_gv1, list_gv2 ! 3.2.2020
  write(10) list_wu1, list_wu2 ! 24.2.2020
  write(10) list_wv1, list_wv2 ! 24.2.2020
  write(10) list_geomu1, list_geomu2 ! 24.2.2020 for geometry check
  write(10) list_geomv1, list_geomv2 ! 24.2.2020 for geometry check
  write(10) nyu1_list, nyu2_list, nyv1_list, nyv2_list ! 27.2.2020 for geometry check
  close(10)
  print*, "Geom complete"
!   open(15,file=trim(dirout)//'boundary_'//ext1//'x'//ext2//'x'//ext3// 'nlist1.dat',form='unformatted')
!    write(15) list_ib1
!   close(15)
!     open(17,file=trim(dirout)//'boundary_'//ext1//'x'//ext2//'x'//ext3// 'nlist2.dat',form='unformatted')
!    write(17) list_ib2
!   close(17)
!   
!   write(*,*) list_ib1
!   stop
!
!  deallocate(l11_ratio, l12_ratio, l21_ratio, l22_ratio) ! 7.2.2020
!  write(*,*) '?????????????????????'
  deallocate(list_ibu1,list_ibu2, list_ibv1, list_ibv2)
  deallocate(list_gu1, list_gu2, list_gv1, list_gv2) ! 3.2.2020
  deallocate(baseu, basev) ! 20.2.2020
  deallocate(nyu1_list, nyu2_list, nyv1_list, nyv2_list) ! 20.2.2020
  deallocate(list_wu1, list_wu2, list_wv1, list_wv2) ! 24.2.2020
  deallocate(list_geomu1, list_geomu2, list_geomv1, list_geomv2) ! 24.2.2020

end subroutine 

real(8) function f(z, y, r)
! Define the function of the cross-section, determine whether points are inside
!  or outside of the solid boundary ,return logical values of each point in the domain.
! Circular cross section, z, y, radius
  implicit none

  real(8) :: z, y
  real(8) :: r

  f = (z - (r + 1))**2 + (y - (r + 1))**2 - r**2 

end function

subroutine interp(k_ele,j_ele,intp1,weicoef,zcorr,ycorr,zbound,ybound,grid)
!----------------------------------------------------------------------
! This routine idnetifies the interpolation stencil and calculates the 
!  weighting coefficients, based on the position of the input forcing
!  point and a given shape function.
!
!----------------------------------------------------------------------
! Linear interpolation
! Ref. RANS solvers with adaptive structured boundary non-conforming 
!  grids. S. Majumdar et al. 2001
!
! The desired velocity on IB points is linearly interpolated by three
!  points nearby: boundary point and two closest fluid points, which 
!  are supposed to always be the two neighbour vertices of one mesh 
!  unit, i.e. not the diagonal vertex. 
! 
! Linear reconstruction of the velocity field:
! Define velocity field u(x,y) = c1 * x + c2 * y + c3 
! Identify three points with known velocity
!  u1(x1,y1) velo on boundary point 1
!  u2(x2,y2) velo on fluid point 2
!  u3(x3,y3) velo on fluid point 3 (clockwise)
! Rewrite above equation is terms of known velocity
!  u(x,y) =   w1 * u1   weighting * velo on boundary point 1
!           + w2 * u2   weighting * velo on fluid point 2
!           + w3 * u3   weighting * velo on fluid point 3 (clockwise)
! Do linear algebra, yield
!    w1         x1 y1 1            x
!  [ w2 ] = [ [ x2 y2 1 ]^-1 ]^T [ y ]
!    w3         x3 y3              1
!
!----------------------------------------------------------------------
!    ^ y (j)
!    |
! y2 |  *p4         *p3
!    |        
! y  |        *p
!    |
! y1 |  *p1(intp1)  *p2
!   -|--------------------> z (k)
!        z1     z     z2
!
! Bilinear interpolation
!  hc = (z - z1)/(z2 - z1)
!  vc = (y - y1)/(y2 - y1)
!  
!  p = (1-h)(1-v)p1 + h(1-v)p2 + hvp3 + (1-h)vp4
!       c1             c2        c3      c4
! Reference gives to W. H. Press et al. Numerical recipes in Fortran.
!----------------------------------------------------------------------

  use declaration
  implicit none

  integer :: i, j
  integer :: z1, z2, y1, y2
  integer :: loc(1)
  integer :: aa, grid
  integer :: n_z1, n_y1
  real(8) :: k_ele, j_ele, bb
  real(8) :: zbound, ybound, zcorr, ycorr, n_z, n_y
  real(8) :: zflow1, yflow1, zflow2, yflow2
  real(8) :: zforc, yforc
  real(8) :: hc, vc 
  real(8) :: z
  integer :: segm ! The curve boundary is devided into segments
  real(8) :: y, mindis
  real(8), allocatable :: dist(:)
  integer :: intp1(4)
  real(8) :: weicoef(3)
  real(8) :: denominator
 
  if (grid .eq. vgrid) then
     bb = 0.5d0
  elseif (grid .eq. ugrid) then
     bb = 0d0
  end if

  ! Find the shortest distance from a IB point to the boundary
  segm = 10000      ! Increase this parameter for better accuracy 
  allocate(dist(segm)) 


  if (j_ele+bb .ge. 1+radius) then
     aa = 1d0  ! Upper half circle
  else 
     aa = -1d0 ! Lower half circle
  end if

  do i = 1, segm
     z = (dstz-1d0)/(segm-1d0)*(i-1d0)+1
     y = aa*sqrt(radius**2-(z-(radius+1d0))**2)+(radius+1d0) ! From shape function
     dist(i) = sqrt((z-k_ele)**2+(y-(j_ele+bb))**2)
  end do
  
  mindis = minval(dist)

  ! Locate the boundary and correspoinding points
  loc = minloc(dist)
  zbound = (dstz-1d0)/(segm-1d0)*(loc(1)-1d0)+1d0
  ybound = aa*sqrt(radius**2-(zbound-(radius+1d0))**2)+(radius+1d0) ! From shape function

    ! Unit normal vector of the boundary at the above point
    ! Vector(n) == grad(f) / |grad(f)|
    ! Pointing outward
    n_z = (zbound-(radius+1d0))/sqrt((zbound-(radius+1d0))**2+(ybound-(radius+1d0))**2)
    n_y = (ybound-(radius+1d0))/sqrt((zbound-(radius+1d0))**2+(ybound-(radius+1d0))**2)
    ! Only take their signs
    n_z1 = n_z / abs(n_z)
    n_y1 = n_y / abs(n_y)

  ! Find the coordinate of the corresponding point
  ! This point is onlt used to locate the interpolation stencil
  ! Forcing point is in the midway between boundary and corresponding point
  zcorr = k_ele + 1d0*mindis*n_z
  ycorr = j_ele+bb + 1d0*mindis*n_y 

  ! Locating the fluid points
  ! Indices, not coordinates
  if (n_z .eq. 0) then
     ! The top or bottom point of a circle
     z1 = k_ele - 1
     y1 = j_ele + n_y1
     z2 = k_ele + 1
     y2 = j_ele + n_y1
  elseif (n_y .eq. 0) then
     ! The left or right point of a circle
     z1 = k_ele + n_z1
     y1 = j_ele - 1
     z2 = k_ele + n_z1
     y2 = j_ele + 1
  else
     z1 = k_ele + n_z1 !floor(zcorr)
     y1 = j_ele        !floor(ycorr - bb)
     z2 = k_ele        !z1 + 1
     y2 = j_ele + n_y1 !y1 + 1 
  end if

  intp1(1) = z1
  intp1(2) = y1
  intp1(3) = z2
  intp1(4) = y2

  ! Calculating the weighting coefficient  
  ! Coordinates, not indices
  zflow1 = z1 * 1d0
  yflow1 = y1 + bb
  zflow2 = z2 * 1d0
  yflow2 = y2 + bb
  zforc  = k_ele
  yforc  = j_ele + bb

  denominator=  (ybound*zflow1 - yflow1*zbound - ybound*zflow2 + yflow2*zbound + yflow1*zflow2 - yflow2*zflow1) 
  weicoef(1) =  (yflow1*zflow2 - yflow2*zflow1 - yflow1*zforc + yforc*zflow1 + yflow2*zforc - yforc*zflow2)/denominator
  weicoef(2) = -(ybound*zflow2 - yflow2*zbound - ybound*zforc + yforc*zbound + yflow2*zforc - yforc*zflow2)/denominator
  weicoef(3) =  (ybound*zflow1 - yflow1*zbound - ybound*zforc + yforc*zbound + yflow1*zforc - yforc*zflow1)/denominator

  deallocate(dist)

  return
end subroutine

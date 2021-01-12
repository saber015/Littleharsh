subroutine boundary_circ_rough
! This function creates the geometry- tall, sparse filaments.
! It defines lists with the points at the forcing points (list_ib1),
!   

! ATTENTION, it raises an error if the number of points in x (Ngal(1,1)) and z (Ngal(2,1))
!   is not a multiple of the number of tiles (ntilex, ntilez)

! It is called by 'getbounds' in 'start.f90'

!                          |- Periodic unit -|
!    dloly    ***          ***               ***                         ***           ***           ***                    
!     .       ***          ***               ***                         ***           ***           ***                    
!    dsty      *            *                 *                           *             *             *                     
! Y   2        *            *                 *                           *             *             *                      
!     1        *            *                 *                           *             *             *                     
!   wall  *    *     *      *         *       *      *      *      *      *      *      *      *      *      *     *  
!                           1      2      3                           7      8     9(1)  10(2)  11(3)
!                                                      X,Z


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



  use declaration
  implicit none
  integer k_ele, i_ele, j,ilist,ilist_xz,ix,iz,shift,grid
  integer nlist_s_xz, nlist_f_xz
  integer nlist_ib_s_bot, nlist_ib_f_bot
  integer points_stem_s, points_stem_f, shift_x, shift_z
  integer, allocatable:: list_ib_s_bot(:,:,:), list_ib_s_top(:,:,:)
  integer, allocatable:: list_ib_f_bot(:,:,:), list_ib_f_top(:,:,:)
  real(8), allocatable:: list_ib_f_w_bot(:,:,:), list_ib_f_w_top(:,:,:)
  integer, allocatable:: list_ib_s_xz(:,:)
  real(8) :: radius
  
  real(8), allocatable :: test_circ(:,:), list_ib_f_xz(:,:) 
  logical :: rig, lef, top, bot
  integer :: interpp(4)
  real(8) :: weicoef(2), zbound, xbound, circ

  
  ! Check that the number of grid points in x and z is a multiple of the number of tiles 
  if (ntilex.eq.0) then
    write(*,*) 'ERROR: ntilex equals 0'
    stop
  end if
  if (ntilez.eq.0) then
    write(*,*) 'ERROR: ntilez equals 0'
    stop
  end if
  if (mod(ntilex,2).ne.0) then
    write(*,*) 'ERROR: ntilex is not even'
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
  if (mod(dstx,2).eq.0) then
    write(*,*) 'WARNING: would be better if dstx was odd'
    stop
  end if

  dnx = Ngal(1,1)/ntilex  ! Width:  Number of points per tile the in streamwise direction
  dnz = Ngal(2,1)/ntilez  ! Width:  Number of points per tile the in spanwise   direction
  shift_x = 2!0.5d0*( dstx)
  shift_z = 2!0.5d0*( dstz)
  nyu11 = -dsty+1
  nyu21 = 0
  nyv11 = -dsty
  nyv21 = 0-1            !-1 To make the boundary align with u_grid

  nyu12 = Ngal(4,3)-dsty+1
  nyu22 = Ngal(4,3)

  nyv12 = Ngal(3,3)+1-dsty+1  !+1 To make the boundary align with u_grid
  nyv22 = Ngal(3,3)+1

!print*, "NGAL", Ngal

print*, "nyv11, 21, 12, 22", nyv11, nyv21, nyv12, nyv22
print*, "nyu11, 21, 12, 22", nyu11, nyu21, nyu12, nyu22

! CREATE CIRCULAR CROSS SECTION IN I-K PLANE MEL
!
!    i
!    ^
!    |      ___          dstx
!    |   _-     -_
!    |  -         -
!    | *     *     *     1+(dtsx-1)/2
!    |  _         _
!    |   -       -
! 1  |      ---          1
! 0 -|-----------------> k
!    0 1  centre  dstx
!
! centre point (1+(dstx-1)/2, 1+(dstx-1)/2)
  radius = (dstx-1)/2d0

  allocate(test_circ(0:dstx+1,0:dstx+1))
  allocate(list_ib_s_xz(2, dstx*dstx)) ! Too big but doesnt matter MEL
  allocate(list_ib_f_xz(10,dstx*dstx))

  do k_ele = 0, dstx+1
    do i_ele = 0, dstx+1
      test_circ(i_ele, k_ele) = circ(i_ele*1d0, k_ele*1d0, radius)
    end do
  end do

  nlist_s_xz = 0 !num of points for one stem in one xz plane
  nlist_f_xz = 0
  list_ib_s_xz = 0
  list_ib_f_xz = 0

  do k_ele = 1, dstx
    do i_ele = 1, dstx

      rig = test_circ(i_ele-1,k_ele  ) <= 0d0
      lef = test_circ(i_ele+1,k_ele  ) <= 0d0
      bot = test_circ(i_ele  ,k_ele+1) <= 0d0
      top = test_circ(i_ele  ,k_ele-1) <= 0d0 
      
      if (test_circ(i_ele, k_ele) .gt. 0d0) then
        if (rig .or. lef .or. bot .or. top) then
          ! This is a immersed boundary point
          nlist_f_xz = nlist_f_xz + 1
          ! Calculating IBP's corresponding points on the solid boundary
          ! Calculating the weighting coefficients for interpolation
          interpp = 0
          weicoef = 0d0
          zbound = 0d0
          xbound = 0d0
          call lin_interp(i_ele*1d0, k_ele*1d0, interpp, weicoef, zbound, xbound, radius, grid)
          ! Store them into the base profile
          list_ib_f_xz(1, nlist_f_xz) = i_ele
          list_ib_f_xz(2, nlist_f_xz) = k_ele
          list_ib_f_xz(3, nlist_f_xz) = interpp(2)   ! i of fluid point 2
          list_ib_f_xz(4, nlist_f_xz) = interpp(1)   ! k of fluid point 2
          list_ib_f_xz(5, nlist_f_xz) = interpp(4)   ! i of fluid point 3
          list_ib_f_xz(6, nlist_f_xz) = interpp(3)   ! k of fluid point 3
!          list_ib_f_xz(11, nlist_f_xz) = weicoef(3) ! weighting of boundary point 1 if v is non-zero
          list_ib_f_xz(7, nlist_f_xz) = weicoef(1) ! weighting of fluid point 2
          list_ib_f_xz(8, nlist_f_xz) = weicoef(2) ! weighting of fluid point 3
          list_ib_f_xz(9, nlist_f_xz) = xbound ! i of circular boundary point
          list_ib_f_xz(10, nlist_f_xz) = zbound ! k of circular boundary point

!write(*,*) 'forcing i, j =', i_ele, k_ele
!write(*,*) 'fluid2  i, j =', interpp(2),interpp(1)
!write(*,*) 'fluid3  i, j =', interpp(4),interpp(3)
!write(*,*) 'xbound,zbound=', xbound,zbound
!write(*,*) 'weighting    =', weicoef(1),weicoef(2),weicoef(3)
!write(*,*)
        end if
      else
           ! This is a solid point
           nlist_s_xz = nlist_s_xz + 1
           list_ib_s_xz(1, nlist_s_xz) = i_ele 
           list_ib_s_xz(2, nlist_s_xz) = k_ele
      end if
    end do
  end do

! CREATE THE PATTERN. 
!  This part can be generalise calling different functions
!  Here we could implement a switch to use different geometries
! if (geometry_flag = 17) then
!   call boundary_blabla
!            _______
!           |       |
!           |__   __|
!              | |
!              | |
!              | |
!              |_|

! CREATE FIRST STEM MEL

  points_stem_s = dsty*nlist_s_xz   ! Total number of forcing points in one stem  
  points_stem_f = dsty*nlist_f_xz   ! Total number of forcing points in one stem 
!  if(shift_stag .eq. 0d0) then
  nlist_ib_s_bot = ntilex*ntilez*points_stem_s    ! Number of points in all stems.    
  nlist_ib_f_bot = ntilex*ntilez*points_stem_f    ! Number of points in all stems.    
!  else 
!  write(*,*) 'shift_stag ', shift_stag
!  nlist_ib_s_bot = (ntilex*(ntilez-1)+ntilex/2)*points_stem_s    ! Number of points in all stems
!  nlist_ib_f_bot = (ntilex*(ntilez-1)+ntilex/2)*points_stem_f    ! Number of points in all stems
!  nlist_ib_f_bot = ntilex*ntilez*points_stem_f    ! Number of points in all stems.    
!  end if

  allocate(list_ib_s_bot(  3,nlist_ib_s_bot,2))
  allocate(list_ib_f_bot(  9,nlist_ib_f_bot,2))
  allocate(list_ib_f_w_bot(2,nlist_ib_f_bot,2))
  allocate(list_ib_s_top(  3,nlist_ib_s_bot,2))
  allocate(list_ib_f_top(  9,nlist_ib_f_bot,2))
  allocate(list_ib_f_w_top(2,nlist_ib_f_bot,2))

! FIRST STEM BOTTOM VGRID
  ilist = 0
  do j = nyv11,nyv21
    do ilist_xz = 1, nlist_s_xz
  
      ilist = ilist + 1
	  
      list_ib_s_bot(1, ilist, vgrid) = list_ib_s_xz(1, ilist_xz) + shift_x ! i coordinate of the forcing point
      list_ib_s_bot(2, ilist, vgrid) = list_ib_s_xz(2, ilist_xz) + shift_z ! k coordinate of the forcing point
      list_ib_s_bot(3, ilist, vgrid) = j                                ! j coordinate of the forcing point 
	
    end do
  end do

  if (ilist/=points_stem_s) then
    write(*,*) 'ERROR: ilist is not equal to points_stem_s 1', ilist, points_stem_s
    stop
  end if 

  ilist = 0
  do j = nyv11, nyv21
    do ilist_xz = 1, nlist_f_xz

      ilist = ilist + 1

      list_ib_f_bot(1, ilist, vgrid) =  list_ib_f_xz(1, ilist_xz) +shift_x
      list_ib_f_bot(2, ilist, vgrid) =  list_ib_f_xz(2, ilist_xz) +shift_z
      list_ib_f_bot(3, ilist, vgrid) =  j

      list_ib_f_bot(4, ilist, vgrid) =  list_ib_f_xz(3, ilist_xz) +shift_x
      list_ib_f_bot(5, ilist, vgrid) =  list_ib_f_xz(4, ilist_xz) +shift_z
      list_ib_f_bot(6, ilist, vgrid) =  j

      list_ib_f_bot(7, ilist, vgrid) =  list_ib_f_xz(5, ilist_xz) +shift_x
      list_ib_f_bot(8, ilist, vgrid) =  list_ib_f_xz(6, ilist_xz) +shift_z
      list_ib_f_bot(9, ilist, vgrid) =  j

      list_ib_f_w_bot(1, ilist, vgrid) =  list_ib_f_xz(7, ilist_xz) 
      list_ib_f_w_bot(2, ilist, vgrid) =  list_ib_f_xz(8, ilist_xz) 

    end do
  end do

  if (ilist/=points_stem_f) then
    write(*,*) 'ERROR: ilist is not equal to points_stem_f 1'
    stop
  end if 
  
! FIRST STEM BOTTOM UGRID
  ilist = 0
  do j = nyu11,nyu21
    do ilist_xz = 1, nlist_s_xz
  
      ilist = ilist + 1
	  
      list_ib_s_bot(1, ilist, ugrid) = list_ib_s_xz(1, ilist_xz) + shift_x ! i coordinate of the forcing point
      list_ib_s_bot(2, ilist, ugrid) = list_ib_s_xz(2, ilist_xz) + shift_z ! k coordinate of the forcing point
      list_ib_s_bot(3, ilist, ugrid) = j                                ! j coordinate of the forcing point 
	
    end do
  end do

  if (ilist/=points_stem_s) then
    write(*,*) 'ERROR: ilist is not equal to points_stem_s 2'
    stop
  end if 

  ilist = 0
  do j = nyu11, nyu21
    do ilist_xz = 1, nlist_f_xz

      ilist = ilist + 1

      list_ib_f_bot(1, ilist, ugrid) =  list_ib_f_xz(1, ilist_xz) +shift_x
      list_ib_f_bot(2, ilist, ugrid) =  list_ib_f_xz(2, ilist_xz) +shift_z
      list_ib_f_bot(3, ilist, ugrid) =  j

      list_ib_f_bot(4, ilist, ugrid) =  list_ib_f_xz(3, ilist_xz) +shift_x
      list_ib_f_bot(5, ilist, ugrid) =  list_ib_f_xz(4, ilist_xz) +shift_z
      list_ib_f_bot(6, ilist, ugrid) =  j

      list_ib_f_bot(7, ilist, ugrid) =  list_ib_f_xz(5, ilist_xz) +shift_x
      list_ib_f_bot(8, ilist, ugrid) =  list_ib_f_xz(6, ilist_xz) +shift_z
      list_ib_f_bot(9, ilist, ugrid) =  j

      list_ib_f_w_bot(1, ilist, ugrid) =  list_ib_f_xz(7, ilist_xz) 
      list_ib_f_w_bot(2, ilist, ugrid) =  list_ib_f_xz(8, ilist_xz) 

    end do
  end do

  if (ilist/=points_stem_f) then
    write(*,*) 'ERROR: ilist is not equal to points_stem_f 2'
    stop
  end if 

! FIRST STEM TOP VGRID
  ilist = 0
  do j = nyv12,nyv22
    do ilist_xz = 1, nlist_s_xz
  
      ilist = ilist + 1
	  
      list_ib_s_top(1, ilist, vgrid) = list_ib_s_xz(1, ilist_xz) + shift_x ! i coordinate of the forcing point
      list_ib_s_top(2, ilist, vgrid) = list_ib_s_xz(2, ilist_xz) + shift_z ! k coordinate of the forcing point
      list_ib_s_top(3, ilist, vgrid) = j                                ! j coordinate of the forcing point 
	
    end do
  end do

  if (ilist/=points_stem_s) then
    write(*,*) 'ERROR: ilist is not equal to points_stem_s 3'
    stop
  end if 

  ilist = 0
  do j = nyv12, nyv22
    do ilist_xz = 1, nlist_f_xz

      ilist = ilist + 1

      list_ib_f_top(1, ilist, vgrid) =  list_ib_f_xz(1, ilist_xz) +shift_x
      list_ib_f_top(2, ilist, vgrid) =  list_ib_f_xz(2, ilist_xz) +shift_z
      list_ib_f_top(3, ilist, vgrid) =  j

      list_ib_f_top(4, ilist, vgrid) =  list_ib_f_xz(3, ilist_xz) +shift_x
      list_ib_f_top(5, ilist, vgrid) =  list_ib_f_xz(4, ilist_xz) +shift_z
      list_ib_f_top(6, ilist, vgrid) =  j

      list_ib_f_top(7, ilist, vgrid) =  list_ib_f_xz(5, ilist_xz) +shift_x
      list_ib_f_top(8, ilist, vgrid) =  list_ib_f_xz(6, ilist_xz) +shift_z
      list_ib_f_top(9, ilist, vgrid) =  j

      list_ib_f_w_top(1, ilist, vgrid) =  list_ib_f_xz(7, ilist_xz) 
      list_ib_f_w_top(2, ilist, vgrid) =  list_ib_f_xz(8, ilist_xz) 

    end do
  end do

  if (ilist/=points_stem_f) then
    write(*,*) 'ERROR: ilist is not equal to points_stem_f 3'
    stop
  end if
 
! FIRST STEM TOP UGRID
  ilist = 0
  do j = nyu12,nyu22
    do ilist_xz = 1, nlist_s_xz
  
      ilist = ilist + 1
	  
      list_ib_s_top(1, ilist, ugrid) = list_ib_s_xz(1, ilist_xz) + shift_x ! i coordinate of the forcing point
      list_ib_s_top(2, ilist, ugrid) = list_ib_s_xz(2, ilist_xz) + shift_z ! k coordinate of the forcing point
      list_ib_s_top(3, ilist, ugrid) = j                                ! j coordinate of the forcing point 
	
    end do
  end do

  if (ilist/=points_stem_s) then
    write(*,*) 'ERROR: ilist is not equal to points_stem_s 4'
    stop
  end if 

  ilist = 0
  do j = nyu12, nyu22
    do ilist_xz = 1, nlist_f_xz

      ilist = ilist + 1

      list_ib_f_top(1, ilist, ugrid) =  list_ib_f_xz(1, ilist_xz) +shift_x
      list_ib_f_top(2, ilist, ugrid) =  list_ib_f_xz(2, ilist_xz) +shift_z
      list_ib_f_top(3, ilist, ugrid) =  j

      list_ib_f_top(4, ilist, ugrid) =  list_ib_f_xz(3, ilist_xz) +shift_x
      list_ib_f_top(5, ilist, ugrid) =  list_ib_f_xz(4, ilist_xz) +shift_z
      list_ib_f_top(6, ilist, ugrid) =  j

      list_ib_f_top(7, ilist, ugrid) =  list_ib_f_xz(5, ilist_xz) +shift_x
      list_ib_f_top(8, ilist, ugrid) =  list_ib_f_xz(6, ilist_xz) +shift_z
      list_ib_f_top(9, ilist, ugrid) =  j

      list_ib_f_w_top(1, ilist, ugrid) =  list_ib_f_xz(7, ilist_xz) 
      list_ib_f_w_top(2, ilist, ugrid) =  list_ib_f_xz(8, ilist_xz) 

    end do
  end do

  if (ilist/=points_stem_f) then
    write(*,*) 'ERROR: ilist is not equal to points_stem_f 4'
    stop
  end if 

! REPLICATE THE PATTERN- STEM
! This section should be common for all geometries
  do grid = 1,2
    do ix = 1,ntilex
      do iz = 1,ntilez
!      do iz = 1,ntilez-mod(ix+1,2)
        shift = points_stem_s*(iz-1) + points_stem_s*ntilez*(ix-1) 
!        shift = points_stem_s*(iz-1 + floor(ix/2d0)) + points_stem_s*(ntilez-1)*(ix-1) 
        do ilist = 1,points_stem_s
          list_ib_s_bot(1,ilist+shift,grid) = list_ib_s_bot(1,ilist,grid) + dnx*(ix-1)  ! i 
          list_ib_s_bot(2,ilist+shift,grid) = 1+mod(list_ib_s_bot(2,ilist,grid) + dnz*(iz-1) &! k
&                                           + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))  
          list_ib_s_bot(3,ilist+shift,grid) = list_ib_s_bot(3,ilist,grid)               ! j 

          list_ib_s_top(1,ilist+shift,grid) = list_ib_s_top(1,ilist,grid) + dnx*(ix-1)  ! i 
          list_ib_s_top(2,ilist+shift,grid) = 1+mod(list_ib_s_top(2,ilist,grid) + dnz*(iz-1) &! k
&                                           + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))  
          list_ib_s_top(3,ilist+shift,grid) = list_ib_s_top(3,ilist,grid)               ! j
        end do
      end do
    end do
  end do
  
  do grid = 1,2
    do ix = 1,ntilex
      do iz = 1,ntilez
!      do iz = 1,ntilez-mod(ix+1,2)
!        shift = points_stem_f*(iz-1 + floor(ix/2d0)) + points_stem_f*(ntilez-1)*(ix-1) 
        shift = points_stem_f*(iz-1) + points_stem_f*ntilez*(ix-1) 
        do ilist = 1,points_stem_f
          list_ib_f_bot(1,ilist+shift,grid) = list_ib_f_bot(1,ilist,grid) + dnx*(ix-1)  ! i
!          list_ib_f_bot(2,ilist+shift,grid) = list_ib_f_bot(2,ilist,grid) + dnz*(iz-1)  ! k 
          list_ib_f_bot(2,ilist+shift,grid) = 1+mod(list_ib_f_bot(2,ilist,grid) + dnz*(iz-1) &! k
&                                           + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))  
          list_ib_f_bot(3,ilist+shift,grid) = list_ib_f_bot(3,ilist,grid)               ! j 
          list_ib_f_bot(4,ilist+shift,grid) = list_ib_f_bot(4,ilist,grid) + dnx*(ix-1)  ! i
!          list_ib_f_bot(5,ilist+shift,grid) = list_ib_f_bot(5,ilist,grid) + dnz*(iz-1)  ! k 
          list_ib_f_bot(5,ilist+shift,grid) = 1+mod(list_ib_f_bot(5,ilist,grid) + dnz*(iz-1) &! k
&                                           + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))  
          list_ib_f_bot(6,ilist+shift,grid) = list_ib_f_bot(6,ilist,grid)               ! j 
          list_ib_f_bot(7,ilist+shift,grid) = list_ib_f_bot(7,ilist,grid) + dnx*(ix-1)  ! i
!          list_ib_f_bot(8,ilist+shift,grid) = list_ib_f_bot(8,ilist,grid) + dnz*(iz-1)  ! k 
          list_ib_f_bot(8,ilist+shift,grid) = 1+mod(list_ib_f_bot(8,ilist,grid) + dnz*(iz-1) &! k
&                                           + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))  
          list_ib_f_bot(9,ilist+shift,grid) = list_ib_f_bot(9,ilist,grid)               ! j 

          list_ib_f_w_bot(1,ilist+shift,grid) = list_ib_f_w_bot(1,ilist,grid)
          list_ib_f_w_bot(2,ilist+shift,grid) = list_ib_f_w_bot(2,ilist,grid)
 
          list_ib_f_top(1,ilist+shift,grid) = list_ib_f_top(1,ilist,grid) + dnx*(ix-1)  ! i 
!          list_ib_f_top(2,ilist+shift,grid) = list_ib_f_top(2,ilist,grid) + dnz*(iz-1)  ! k 
          list_ib_f_top(2,ilist+shift,grid) = 1+mod(list_ib_f_top(2,ilist,grid) + dnz*(iz-1) &! k
&                                           + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))  
          list_ib_f_top(3,ilist+shift,grid) = list_ib_f_top(3,ilist,grid)               ! j 
          list_ib_f_top(4,ilist+shift,grid) = list_ib_f_top(4,ilist,grid) + dnx*(ix-1)  ! i 
!          list_ib_f_top(5,ilist+shift,grid) = list_ib_f_top(5,ilist,grid) + dnz*(iz-1)  ! k 
          list_ib_f_top(5,ilist+shift,grid) = 1+mod(list_ib_f_top(5,ilist,grid) + dnz*(iz-1) &! k
&                                           + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))  
          list_ib_f_top(6,ilist+shift,grid) = list_ib_f_top(6,ilist,grid)               ! j 
          list_ib_f_top(7,ilist+shift,grid) = list_ib_f_top(7,ilist,grid) + dnx*(ix-1)  ! i 
!          list_ib_f_top(8,ilist+shift,grid) = list_ib_f_top(8,ilist,grid) + dnz*(iz-1)  ! k 
          list_ib_f_top(8,ilist+shift,grid) = 1+mod(list_ib_f_top(8,ilist,grid) + dnz*(iz-1) &! k
&                                           + mod(ix+1,2)*dnz*shift_stag-1d0,1d0*Ngal(2,3))  
          list_ib_f_top(9,ilist+shift,grid) = list_ib_f_top(9,ilist,grid)               ! j 


          list_ib_f_w_top(1,ilist+shift,grid) = list_ib_f_w_top(1,ilist,grid)                
          list_ib_f_w_top(2,ilist+shift,grid) = list_ib_f_w_top(2,ilist,grid)                
        end do
      end do
    end do
  end do
! Save the lists into a file
  open(10,file=trim(dirout)//'boundary_'//ext1//'x'//ext2//'x'//ext3//'.dat',form='unformatted',access='stream')
  write(10) Lx,Ly,Lz
  write(10) Ngal,nlist_ib_s_bot,nlist_ib_f_bot,nyu11,nyu21,nyu12,nyu22,nyv11,nyv21,nyv12,nyv22
  write(10) list_ib_s_bot, list_ib_s_top
  write(10) list_ib_f_bot, list_ib_f_top
  write(10) list_ib_f_w_bot, list_ib_f_w_top
  write(10) list_ib_s_xz, list_ib_f_xz
  close(10)
print*, "Geom complete"
!  write(*,*) 'list_ib_f_w_bot', list_ib_f_w_bot 
  deallocate(test_circ)
  deallocate(list_ib_s_xz, list_ib_f_xz)
  deallocate(list_ib_s_bot, list_ib_s_top)
  deallocate(list_ib_f_bot, list_ib_f_top)
  deallocate(list_ib_f_w_bot, list_ib_f_w_top)
 
end subroutine 

real(8) function circ(x, z, r)
  
  implicit none 

  real(8) :: x, z, r

  circ = (x - (r + 1))**2 + (z - (r + 1))**2 - r**2

end function

subroutine lin_interp(i_ele, k_ele, interpp, weicoef, zbound, xbound, radius, grid)

  use declaration
  implicit none

  real(8) :: i_ele, k_ele, x, z
  integer :: z2, x2, z3, x3, interpp(4)
  integer :: n_z_sign, n_x_sign, aa, grid 
  integer :: n_segm, segm, min_loc(1)
  real(8) :: radius, min_dist
  real(8) :: zbound, xbound, n_z, n_x
  real(8) :: weicoef(2), denominator 

  real(8), allocatable :: dist(:)

  n_segm = 10000

  allocate(dist(n_segm))

  if (i_ele .ge. 1d0+radius) then
    aa = 1d0
  else
    aa = -1d0
  end if

!  write(*,*)
!  write(*,*) 'i_ele', i_ele 
!  write(*,*) 'k_ele', k_ele

  do segm = 1, n_segm
    z = (segm-1d0)*2d0*radius/(n_segm-1d0) +1d0
    x = aa*sqrt(radius**2-(z-(radius+1d0))**2) + radius + 1d0 
    dist(segm) = sqrt((x-i_ele)**2 + (z-k_ele)**2) 
  end do

  min_dist = minval(dist)
  min_loc = minloc(dist)

  zbound = (min_loc(1)-1d0)*2d0*radius/(n_segm-1d0) + 1d0
  xbound = aa*sqrt(radius**2-(zbound-(radius+1d0))**2) + radius + 1d0

!  write(*,*) 'min_dist', min_dist 
!  write(*,*) 'min_loc(1)', min_loc(1)
!  write(*,*) 'xbound', xbound, 'zbound', zbound
!  write(*,*)

  n_z = (zbound-(radius+1d0))/sqrt((zbound-(radius+1d0))**2+(xbound-(radius+1d0))**2)
  n_x = (xbound-(radius+1d0))/sqrt((zbound-(radius+1d0))**2+(xbound-(radius+1d0))**2)

  n_z_sign = n_z / abs(n_z)
  n_x_sign = n_x / abs(n_x)

!write(*,*) 'i_ele, k_ele', i_ele, k_ele
!write(*,*) 'n_x', n_x
!write(*,*) 'n_z', n_z
!write(*,*) 'n_x_sign', n_x_sign
!write(*,*) 'n_z_sign', n_z_sign
!write(*,*) 'i_ele -r', abs(i_ele-radius-1d0)
!write(*,*) 'k_ele -r', abs(k_ele-radius-1d0)
!write(*,*) 

  if (abs(i_ele-radius-1d0) .eq. abs(k_ele-radius-1d0)) then
!  if (abs(n_x) .eq. abs(n_z)) then !only accurate to 3dp so doesnt work
  ! Diagonal point 
    z2 = k_ele + n_z_sign
    x2 = i_ele  
    z3 = k_ele 
    x3 = i_ele + n_x_sign 

  elseif (abs(n_x) .lt. abs(n_z)) then 
  ! Between 10.5 o clock and 1.5 o clock  
  ! Between 4.5 o clock and 7.5 o clock 
    z2 = k_ele + n_z_sign
    x2 = i_ele   
    z3 = k_ele + n_z_sign
    x3 = i_ele + n_x_sign

  elseif (abs(n_x) .gt. abs(n_z)) then 
  ! Between 1.5 o clock and 4.5 o clock 
  ! Between 1.5 o clock and 4.5 o clock 
    z2 = k_ele 
    x2 = i_ele + n_x_sign 
    z3 = k_ele + n_z_sign
    x3 = i_ele + n_x_sign

  else
 
  write(*,*) 'Mel fucked up'
  stop
 
  end if

  interpp(1) = z2 
  interpp(2) = x2 
  interpp(3) = z3 
  interpp(4) = x3 

  denominator=  (xbound*z2 - x2*zbound - xbound*z3 + x3*zbound + x2*z3 - x3*z2)

!  weicoef(1) =  (x2*z3 - x3*z2 - x2*k_ele + i_ele*z2 + x3*k_ele - i_ele*z3)/denominator ! weighting for boundary point if v is non-zero

  weicoef(1) = -(xbound*z3 - x3*zbound - xbound*k_ele + i_ele*zbound + x3*k_ele - i_ele*z3)/denominator

  weicoef(2) =  (xbound*z2 - x2*zbound - xbound*k_ele + i_ele*zbound + x2*k_ele - i_ele*z2)/denominator 

  deallocate(dist)

end subroutine


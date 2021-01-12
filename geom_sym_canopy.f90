subroutine boundary_canopies
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
   integer i,k,j,i1,j1,k1,ilist,ix,iz,dny,shift,grid
   integer nlist_ib1, nlist_ib2
   integer points_stem, points_lol, shift_x, shift_z
   integer, allocatable:: list_ib1(:,:,:), list_ib2(:,:,:)
  

  
   ! Check that the number of grid points in x and z is a multiple of the number of tiles
   if (ntilex.eq.0) then
      write(*,*) 'ERROR: ntilex equals 0'
      stop
   end if
   if (ntilez.eq.0) then
      write(*,*) 'ERROR: ntilez equals 0'
      stop
   end if
   if (mod(Ngal(1,1),ntilez)/=0) then
      write(*,*) 'ERROR: nx not multiple of ntilex'
      stop
   end if
   if (mod(Ngal(2,1),ntilez)/=0) then
      write(*,*) 'ERROR: nz not multiple of ntilez'
      stop
   end if

   dnx = Ngal(1,1)/ntilex  ! Width:  Number of points per tile the in streamwise direction
   dnz = Ngal(2,1)/ntilez  ! Width:  Number of points per tile the in spanwise   direction
   shift_x = 0!0.5d0*( dstx)
   shift_z = 0!0.5d0*( dstz)
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

   ! First stem
   points_stem = (dsty)*dstx*dstz   ! Total number of forcing points in one stem
   nlist_ib1 = ntilex*ntilez*points_stem    ! Number of points in all stems.
   allocate(list_ib1(3,nlist_ib1,2))

   ilist = 0
   do j = nyv11,nyv21 !
      j1 = j
      do i = 1,dstx
         i1 = i+shift_x
         do k = 1,dstz
            k1 = k+shift_z
            ilist = ilist + 1
	  
            !	  list_ib1( 1,ilist,ugrid) = i1     ! i coordinate of the forcing point
            !	  list_ib1( 2,ilist,ugrid) = k1     ! k coordinate of the forcing point
            !	  list_ib1( 3,ilist,ugrid) = j1   ! j coordinate of the forcing point
	  
            list_ib1( 1,ilist,vgrid) = i1     ! i coordinate of the forcing point
            list_ib1( 2,ilist,vgrid) = k1     ! k coordinate of the forcing point
            list_ib1( 3,ilist,vgrid) = j1     ! j coordinate of the forcing point
	
         end do
      end do
   end do

  
   if (ilist/=points_stem) then
      write(*,*) 'ERROR: ilist is not equal to points_stem'
      stop
   end if
  
   ilist = 0
   do j = nyu11,nyu21 !
      j1 = j
      do i = 1,dstx
         i1 = i+shift_x
         do k = 1,dstz
            k1 = k+shift_z
            ilist = ilist + 1
	  
            list_ib1( 1,ilist,ugrid) = i1     ! i coordinate of the forcing point
            list_ib1( 2,ilist,ugrid) = k1     ! k coordinate of the forcing point
            list_ib1( 3,ilist,ugrid) = j1   ! j coordinate of the forcing point
	  
         !	  list_ib1( 1,ilist,vgrid) = i1     ! i coordinate of the forcing point
         !	  list_ib1( 2,ilist,vgrid) = k1     ! k coordinate of the forcing point
         !	  list_ib1( 3,ilist,vgrid) = j1     ! j coordinate of the forcing point
	
         end do
      end do
   end do



   if (ilist/=points_stem) then
      write(*,*) 'ERROR: ilist is not equal to points_stem'
      stop
   end if

   nlist_ib2 = ntilex*ntilez*points_stem    ! Number of points in all stems.
  
   print*, "Nlist ", nlist_ib1, nlist_ib2

   allocate(list_ib2(3,nlist_ib1,2))

   ilist = 0
   do j = nyv12,nyv22 !
      j1 = j
      do i = 1,dstx
         i1 = i+shift_x
         do k = 1,dstz
            k1 = k+shift_z
            ilist = ilist + 1
	  
            !	  list_ib2( 1,ilist,ugrid) = i1     ! i coordinate of the forcing point
            !	  list_ib2( 2,ilist,ugrid) = k1     ! k coordinate of the forcing point
            !	  list_ib2( 3,ilist,ugrid) = j1     ! j coordinate of the forcing point
	  
            list_ib2( 1,ilist,vgrid) = i1     ! i coordinate of the forcing point
            list_ib2( 2,ilist,vgrid) = k1     ! k coordinate of the forcing point
            list_ib2( 3,ilist,vgrid) = j1-1     ! j coordinate of the forcing point
	
         end do
      end do
   end do

 
   if (ilist/=points_stem) then
      write(*,*) 'ERROR: ilist is not equal to points_stem'
      stop
   end if

   ilist = 0
   do j = nyu12,nyu22 !
      j1 = j
      do i = 1,dstx
         i1 = i+shift_x
         do k = 1,dstz
            k1 = k+shift_z
            ilist = ilist + 1
	  
            list_ib2( 1,ilist,ugrid) = i1     ! i coordinate of the forcing point
            list_ib2( 2,ilist,ugrid) = k1     ! k coordinate of the forcing point
            list_ib2( 3,ilist,ugrid) = j1     ! j coordinate of the forcing point
	  
         !	  list_ib2( 1,ilist,vgrid) = i1     ! i coordinate of the forcing point
         !	  list_ib2( 2,ilist,vgrid) = k1     ! k coordinate of the forcing point
         !	  list_ib2( 3,ilist,vgrid) = j1-1     ! j coordinate of the forcing point
	
         end do
      end do
   end do
 
   if (ilist/=points_stem) then
      write(*,*) 'ERROR: ilist is not equal to points_stem'
      stop
   end if

  
   ! REPLICATE THE PATTERN- STEM
   !   This section should be common for all geometries
   do grid = 1,2
      do ix = 1,ntilex
         do iz = 1,ntilez
            shift = points_stem*(iz-1) + points_stem*ntilez*(ix-1)
            do ilist = 1,points_stem
               list_ib1(1,ilist+shift,grid) = list_ib1(1,ilist,grid) + dnx*(ix-1)     ! i coordinate
               list_ib1(2,ilist+shift,grid) = list_ib1(2,ilist,grid) + dnz*(iz-1)     ! k coordinate
               list_ib1(3,ilist+shift,grid) = list_ib1(3,ilist,grid)                  ! j coordinate

               list_ib2(1,ilist+shift,grid) = list_ib2(1,ilist,grid) + dnx*(ix-1)     ! i coordinate
               list_ib2(2,ilist+shift,grid) = list_ib2(2,ilist,grid) + dnz*(iz-1)     ! k coordinate
               list_ib2(3,ilist+shift,grid) = list_ib2(3,ilist,grid)                  ! j coordinate
            end do
         end do
      end do
   end do
  


   ! Save the lists into a file
   open(10,file=trim(dirout)//'boundary_'//ext1//'x'//ext2//'x'//ext3//'.dat',form='unformatted',access='stream')
   write(10) Lx,Ly,Lz
   write(10) Ngal,nlist_ib1,nlist_ib2, nyu11, nyu21, nyu12, nyu22, nyv11, nyv21, nyv12, nyv22
   write(10) list_ib1,list_ib2
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
   deallocate(list_ib1,list_ib2)
end subroutine 

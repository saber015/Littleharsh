subroutine y_grid_substrate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!   YGRID for substrates  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Creates the geometry in the y-direction

!  y is the physical coordinate
!  theta is the mapping coordinate (constant increments: dtheta)
!  Since (1) dtheta is constant and (2) the mapping between y and
!  theta is analytical,
!  the second order of the centered finite difference is preserved.
 
!  Three regions:
!  1.- Channel: y(j) = a/5(j-q)**p + b/3(j-q)**(p-2) + c(j-q) + d 
!      nn   is the total number of points within the channel (in vgrid), j in vgrid (0, nn)
!      q    is the index of the centerline (nn/2) (y(q) = 0)
!      dyq  is dymax/dymin (max stretching) within the channel
!      p    degree of polinomial for y grid within the channel. Here 5
!      a, b, c and d are obtained such that dy(q)/dy(0)=dyq, y(0)=-1, y(nn)=1, d2y(0)=0 
!
!      For non-stretching y gird                  dy(0)=dyp, y(0)=-1, y(nn)=1, d2y(0)=0
!
!  2.- First tile: 
!      y(j) = a/5(j-q)**p + b/4(j-q)**(p-1) + c/3(j-q)**(p-2) + d/2(j-q)**(p-3) + e/3(j-q)**(p) + f (eqiv. upper subs)
!      posth        height of permeable substrate in outer units
!      nsolid       nb of points in solid element in vgrid (in ugrid, 1 more point), in dy-constant region
!      nfluid       nb of points in fluid element in vgrid (in ugrid, 1 fewer point), in dy-constant region
!      dsty         nb of points in y in dy-constant region
!      dyp          Deltay in the constant-dy zone (found with posth, nsolid, nfluid, dsty)
!      sizep        size of the first tile in y in outer units (found with nsolid, nfluid, dyp)
!      qp           nb of points in the first tile (found iteratively)

!      a, b, c, d, e and f are obtained such that y(0)=-1, dy(0)=dy(0)chan, d2y(0)=d2y(0)chan,
!                   y(-qp)=-1, dy(-qp)=dyp, d2y(-qp)=0    
!
!  3.- Deep within the substrate: y(j) = dyp*j - (1+sizep) + dyp*qp (eqiv. upper subs)

   use declaration
   implicit none
   
   integer :: j,iband, dny
   real(8) :: aaa, bbb, ccc, qqq
   real(8) :: aaain1, bbbin1, cccin1, eeein1, fffin1
   real(8) :: aaain1_t, bbbin1_t, cccin1_t, eeein1_t, fffin1_t 
   real(8) :: aat, aab, dy1 
   real(8) :: jqp, dymin, dymax, sizep, dyv0, yu_ghost, dyinterface
!   integer :: qp
   real(8), allocatable :: yuin1(:), yvin1(:), dyuin1(:), dyvin1(:)
  
  ! v points in y go from 0 at the bottom interface to nn at the top interface
  ! Attention: the value of nn changes along the subroutine. Bad practice, this should be changed
   nn  = N(3,nband)+1
   print*, "nn = ", nn
   qqq = nn/2d0
   
   !ygrid coefficients within the channel
   !aaa = 15d0 * (dyq - 1d0) / ( (8d0*dyq + 7d0)*(qqq**ppp) )
   !bbb =-30d0 * (dyq - 1d0) / ( (8d0*dyq + 7d0)*(qqq**(ppp-2)) )
   !ccc = 15d0 * dyq / ( (8d0*dyq + 7d0)*qqq )
   
   !--------------------------------------------------------------------
   ! Find qp (nb of points in the first tile)
   !--------------------------------------------------------------------
   ! dyv0: dv at y=0 at vgrid using the grid within the channel
   !dyv0 = aaa*(0-qqq)**(ppp-1) + bbb*(0-qqq)**(ppp-3) + ccc
   ! The small offset between the interface yv(0)=-1 and where the solid elements start yu(0) (ghost point)
   yu_ghost = -1d0-pi/(ntilez*(nsolid+nfluid)*1d0)/0.5d0 ! 12.3.2020 !aaa*(-.5d0-qqq)**ppp/5d0 + bbb*(-.5d0-qqq)**(ppp-2)/3d0 + ccc*(-.5d0-qqq) 
   dyinterface = -(yu_ghost - (-1d0)) ! In reality, not exactly this, but better than nothing
   write(*,*) 'yu_ghost, dyinterface ', yu_ghost, dyinterface
   
   !nsolid = 3       ! nb of points in solid element in vgrid (in ugrid, 1 more point)
   !nfluid   = 3       ! nb of points in the fluid in vgrid (in ugrid, 1 fewer point)
   !dyp = 1d0/180d0       ! Deltay in the constant-Deltay zone
   !sizep  = (nsolid + nfluid) * dyp   ! size of the first tile in y
   ! posth = sizep + (dsty-1)*dyp + dyp/2, where sizep = (nsolid+nfluid)*dyp + dyinterface
   dyp = (posth - dyinterface) / (nfluid + nsolid + dsty -0.5d0) 
   write(*,*) 'posth, dyp, nsolid, nfluid ', posth, dyp, nsolid, nfluid 
   sizep = (nsolid + nfluid)*dyp
   write(*,*) 'sizep, dyp ', sizep, dyp

   qp = nsolid + nfluid ! 12.3.2020 ! ceiling(sizep/dyv0) ! max nb of points in the first tile in vgrid
   write(*,*) 'sizep, dyv0', sizep, dyv0
   write(*,*) 'qp', qp
   write(*,*) 'dyp', dyp
   dymin = 0d0;
   dymax = dyp + 1d0; 

   !ygrid coefficients within the channel ! 12.3.2020
   aaa = 5d0*(dyp*qqq-1d0)*(ppp-2d0) /(2d0*(-qqq)**ppp*(ppp-1d0))
   bbb =-3d0*ppp*qqq**2*(dyp*qqq-1d0)/(2d0*(-qqq)**ppp*(ppp-3d0))
   ccc =-(2d0*ppp-3d0*dyp*qqq-ppp**2+2d0*dyp*ppp*qqq)/(qqq*(ppp**2-4d0*ppp+3d0))
   
   ! Bottom substrate --------------------------
   !do while ( (dymin<dyv0) .or. (dymax>dyp) ) 
	!	qp = qp - 1
	!	jqp = qp + 0.5d0 ! The limit for this region is in ugrid, i.e. yv(-qp-0.5)
!		if (qp /= ceiling(sizep*dyv0)-1) then
!			deallocate(yvin1, yuin1)
!			deallocate(dthdyv_in1, dthdyu_in1)
!		end if
	!	
	!	allocate(yvin1(-qp:-1), yuin1(-qp:0))
	!	allocate(dyvin1(-qp:-1), dyuin1(-qp:0))

	!	eeein1 = -aaa*qqq**4 + ccc
	!	fffin1 = -1d00
	!	aaain1 = -15d0/jqp**4 *(eeein1 + dyp) + 30d0*fffin1/jqp**5 + 30d0*(1d0 + sizep) / jqp**5
	!	bbbin1 = 2d0*aaain1*jqp - 2d0*eeein1/jqp**3 + 2d0*dyp/jqp**3
	!	cccin1 = 0.5d0 * (-4d0*aaain1*jqp**2 + 3d0*bbbin1*jqp)
	
	!	yuin1(-qp)= aaain1*(-qp - .5)**ppp/5.0 + bbbin1*(-qp - .5)**(ppp-1)/4.0 &
	!			&+ cccin1*(-qp - .5)**(ppp-2)/3.0 + eeein1*(-qp - .5)  + fffin1
	!	dyuin1(-qp) = aaain1*(-qp - .5)**(ppp-1) + bbbin1*(-qp - .5)**(ppp-2) &
	!			&+ cccin1*(-qp - .5)**(ppp-3) + eeein1
	!	do j = -qp, -1
	!		yvin1(j)  = aaain1*(j     )**ppp/5.0 + bbbin1*(j     )**(ppp-1)/4.0 &
	!			&+ cccin1*(j     )**(ppp-2)/3.0 + eeein1*(j     )  + fffin1
	!		yuin1(j+1)= aaain1*(j + .5)**ppp/5.0 + bbbin1*(j + .5)**(ppp-1)/4.0 &
	!			&+ cccin1*(j + .5)**(ppp-2)/3.0 + eeein1*(j + .5)  + fffin1
   
	!		dyvin1(j)   = aaain1*(j     )**(ppp-1) + bbbin1*(j     )**(ppp-2) &
	!			&+ cccin1*(j     )**(ppp-3) + eeein1
	!		dyuin1(j+1) = aaain1*(j + .5)**(ppp-1) + bbbin1*(j + .5)**(ppp-2) &
	!			&+ cccin1*(j + .5)**(ppp-3) + eeein1
	!	end do
	!	dymin = minval(dyvin1)
	!	write(*,*) 'while loop: qp', qp
	!	write(*,*) 'min dv, dyvch', dymin, dyv0
	!	!write(*,*) 'dyvin1 ', dyvin1
	!	dymax = maxval(dyvin1)
	!		
	!	deallocate(yvin1, yuin1)
	!	deallocate(dyvin1, dyuin1)
   !end do
   !qp = nsolid + nfluid   
   write(*,*) 'final qp is', qp
   
   !-------------------------------------------------------------------
   
   ! dny: nb of points in y within the substrate (vgrid)
   dny = dsty + qp;
   allocate(yu    (  -dny:nn+dny+1))
   allocate(dyu2i (3,-dny:nn+dny+1))
   allocate(dthdyu(  -dny:nn+dny+1))
  
   allocate(yv    (  -dny:nn+dny  ))
   allocate(dyv2i (3,-dny:nn+dny  ))
   allocate(dthdyv(  -dny:nn+dny  ))
   
   !--------------------------------------------------------------------
   ! 1. Within the channel
   !		vgrid --> j \in [0, nn]
   !		ugrid --> j \in [1, nn]
   !--------------------------------------------------------------------
   ! coefficients aaa, bbb, ccc above

   !rem dy1 = aaa*(-qqq)**4 + bbb*(-qqq)**2 + ccc
   
   !rem aab = (posth-dy1*dny)/(dny**3)
  
   dtheta     = 1d0
   dthetai    = 1d0/dtheta
   ddthetavi  = 1.0d0/dtheta

   do j = 0, nn-1
     yv(j)   = aaa*(j     -qqq)**ppp/5d0 + bbb*(j     -qqq)**(ppp-2)/3d0 + ccc*(j     -qqq)
     yu(j+1) = aaa*(j+.5d0-qqq)**ppp/5d0 + bbb*(j+.5d0-qqq)**(ppp-2)/3d0 + ccc*(j+.5d0-qqq)
   
     dthdyv(j)   = 1d0 / ( aaa*(j     -qqq)**(ppp-1) + bbb*(j     -qqq)**(ppp-3) + ccc )
     dthdyu(j+1) = 1d0 / ( aaa*(j+.5d0-qqq)**(ppp-1) + bbb*(j+.5d0-qqq)**(ppp-3) + ccc )
   end do
   yv(nn) = aaa*(nn-qqq)**ppp/5d0 + bbb*(nn-qqq)**(ppp-2)/3d0 + ccc*(nn-qqq)
   dthdyv(nn) = 1d0/(aaa*(nn-qqq)**(ppp-1) + bbb*(nn-qqq)**(ppp-3) + ccc)


!   !--------------------------------------------------------------------
!   ! 2. First tile
!   !        Bottom substrate
!   !		vgrid --> j \in [-qp, -1]
!   !		ugrid --> j \in [-qp, 0]
!   !        Top substrate
!   !		vgrid --> j \in [nn+1, nn+qp]
!   !		ugrid --> j \in [nn+1, nn+qp+1]
!   !--------------------------------------------------------------------
!   
!   ! Bottom substrate --------------------------------------------------
!   jqp = qp + 0.5d0
!   eeein1 = 0d0 !-aaa*qqq**4 + ccc
!   fffin1 = 0d0 !-1d0
!   aaain1 = 0d0 !-15d0/jqp**4 *(eeein1 + dyp) + 30d0*fffin1/jqp**5 + 30d0*(1d0 + sizep)/jqp**5
!   bbbin1 = 0d0 !2d0*aaain1*jqp - 2d0*eeein1/jqp**3 + 2d0*dyp/jqp**3
!   cccin1 = 0d0 !0.5d0 * (-4d0*aaain1*jqp**2 + 3d0*bbbin1*jqp)
!
!   yu(-qp)= dyp * (-qp - .5d0) - (1d0 + sizep) + dyp * (qp + 0.5d0) !aaain1*(-qp - .5)**ppp/5.0 + bbbin1*(-qp - .5)**(ppp-1)/4.0 &
!		!&+ cccin1*(-qp - .5)**(ppp-2)/3.0 + eeein1*(-qp - .5)  + fffin1
!   dthdyu(-qp) = 1d0 / dyp !( aaain1*(-qp - .5d0)**(ppp-1) + bbbin1*(-qp - .5d0)**(ppp-2) &
!		!&+ cccin1*(-qp - .5)**(ppp-3) + eeein1 )
!   do j = -qp, -1
!		yv(j)   = dyp * (j       ) - (1d0 + sizep) + dyp * (qp + 0.5d0) !aaain1*(j     )**ppp/5d0 + bbbin1*(j     )**(ppp-1)/4d0 &
!			!&+ cccin1*(j     )**(ppp-2)/3d0 + eeein1*(j     )  + fffin1
!		yu(j+1) = dyp * (j + .5d0) - (1d0 + sizep) + dyp * (qp + 0.5d0) !aaain1*(j + .5d0)**ppp/5d0 + bbbin1*(j + .5d0)**(ppp-1)/4d0 &
!			!&+ cccin1*(j + .5d0)**(ppp-2)/3d0 + eeein1*(j + .5d0)  + fffin1
!  
!		dthdyv(j)   = 1d0 / dyp !( aaain1*(j       )**(ppp-1) + bbbin1*(j       )**(ppp-2) &
!			!&+ cccin1*(j     )**(ppp-3) + eeein1 )
!		dthdyu(j+1) = 1d0 / dyp !( aaain1*(j + .5d0)**(ppp-1) + bbbin1*(j + .5d0)**(ppp-2) &
!			!&+ cccin1*(j + .5)**(ppp-3) + eeein1 )
!   end do
! 
!   ! Top substrate -----------------------------------------------------
!	eeein1_t = 0d0 !-aaa * qqq**4 + ccc
!	fffin1_t = 0d0 !+1d0
!	aaain1_t = 0d0 !-15d0/jqp**4 * ( eeein1_t + dyp) - 30d0 * fffin1_t/jqp**5 + 30d0*(1 + sizep)/jqp**5
!	bbbin1_t = 0d0 !-2d0*aaain1_t*jqp + 2d0*eeein1_t/jqp**3 - 2d0*dyp/jqp**3
!	cccin1_t = 0d0 !0.5d0*( -4d0*aaain1_t*jqp**2 - 3d0*bbbin1_t*jqp )
!
!	do j = 1, qp
!		yv(nn+j) = dyp * (nn + j       ) + (1d0 + sizep) - dyp * (nn + qp + 0.5d0) !aaain1_t*(j       )**ppp/5d0 + bbbin1_t*(j       )**(ppp-1)/4d0 &
!			!&+ cccin1_t*(j       )**(ppp-2)/3.0 + eeein1_t*(j       )  + fffin1_t 
!		yu(nn+j) = dyp * (nn + j - .5d0) + (1d0 + sizep) - dyp * (nn + qp + 0.5d0) !aaain1_t*(j - .5d0)**ppp/5d0 + bbbin1_t*(j - .5d0)**(ppp-1)/4d0 &
!			!&+ cccin1_t*(j - .5d0)**(ppp-2)/3d0 + eeein1_t*(j - .5d0)  + fffin1_t 
!   
!		dthdyv(nn+j) = 1d0 / dyp !( aaain1_t*(j       )**(ppp-1) + bbbin1_t*(j       )**(ppp-2) &
!			!&+ cccin1_t*(j       )**(ppp-3) + eeein1_t )
!		dthdyu(nn+j) = 1d0 / dyp !( aaain1_t*(j - .5d0)**(ppp-1) + bbbin1_t*(j - .5d0)**(ppp-2) &
!			!&+ cccin1_t*(j - .5d0)**(ppp-3) + eeein1_t )
!	end do
!    yu(nn+qp+1) = dyp * (nn + qp + .5d0) + (1d0 + sizep) - dyp * (nn + qp + 0.5d0) !aaain1_t*(qp + .5d0)**ppp/5d0 + bbbin1_t*(qp + .5d0)**(ppp-1)/4d0 &
!			!&+ cccin1_t*(qp + .5d0)**(ppp-2)/3d0 + eeein1_t*(qp + .5d0)  + fffin1_t
!	dthdyu(nn+qp+1) = 1d0 / dyp !( aaain1_t*(qp + .5d0)**(ppp-1) + bbbin1_t*(qp + .5d0)**(ppp-2) &
!			!&+ cccin1_t*(qp + .5d0)**(ppp-3) + eeein1_t );
!   !--------------------------------------------------------------------
!   ! 3. Deep within the substrate
!   !        Bottom substrate
!   !		vgrid --> j \in [-dsty-qp, -qp-1]
!   !		ugrid --> j \in [-dsty-qp, -qp-1]
!   !        Top substrate
!   !		vgrid --> j \in [nn+qp+1, nn+qp+dsty]
!   !		ugrid --> j \in [nn+qp+1, nn+qp+dsty+1]
!   !--------------------------------------------------------------------
!   
!   ! Bottom substrate ----------------------------------------
!	do j = -dsty-qp, -qp-1 ! 4.12.2019
!		yv(j) = dyp * (j       ) - (1d0 + sizep) + dyp * (qp + 0.5d0)
!		yu(j) = dyp * (j - .5d0) - (1d0 + sizep) + dyp * (qp + 0.5d0) 
!    
!		dthdyv(j) = 1d0 / dyp
!		dthdyu(j) = 1d0 / dyp
!	end do
!   
!   ! Top substrate -----------------------------------------
!   do j = nn+qp+1, nn+dsty+qp ! 4.12.2019
!		yv(j) = dyp * (j       ) + (1d0 + sizep) - dyp * (nn + qp + 0.5d0)
!		yu(j+1) = dyp * (j + .5d0) + (1d0 + sizep) - dyp * (nn + qp + 0.5d0)
!    
!		dthdyv(j) = 1d0 / dyp
!		dthdyu(j+1) = 1d0 / dyp
!   end do
!
!   print*, "dthdyu = ", dthdyu(0), " ", dthdyu(nn+1)
!   print*, "dthdyv = ", dthdyv(0), " ", dthdyv(nn)
!   print*, "Shape u, v", shape(dthdyu), ", ", shape(dthdyv)

   !--------------------------------------------------------------------
   ! 2. Within the substrate
   !        Bottom substrate
   !    vgrid --> j \in [-dsty-qp, -1]
   !    ugrid --> j \in [-dsty-qp,  0]
   !        Top substrate
   !    vgird --> j \in [nn+1, nn+qp+dsty  ]
   !    ugird --> j \in [nn+1, nn+qp+dsty+1]
   !--------------------------------------------------------------------

   ! Bottom substrate --------------------------------------------------
   do j = -dsty-qp, -1
      yv(j) = dyp * (j       ) - 1d0 !- (1d0 + sizep) + dyp * (qp + 0.5d0) 
      yu(j) = dyp * (j - .5d0) - 1d0 !- (1d0 + sizep) + dyp * (qp + 0.5d0)

      dthdyv(j) = 1d0 / dyp
      dthdyu(j) = 1d0 / dyp
   end do

   yu(0) = dyp * (0d0 - .5d0) - 1d0 !- (1d0 + sizep) + dyp * (qp + 0.5d0)
   dthdyu(0) = 1d0 / dyp

   ! Top substrate -----------------------------------------------------
   do j = nn+1, nn+qp+dsty
      yv(j  ) = dyp * (j       ) + 1d0 - nn * dyp !+ (1d0 + sizep) - dyp * (nn + qp + 0.5d0)
      yu(j+1) = dyp * (j + .5d0) + 1d0 - nn * dyp !+ (1d0 + sizep) - dyp * (nn + qp + 0.5d0)

      dthdyv(j) = 1d0 / dyp
      dthdyu(j+1) = 1d0 / dyp
   end do

   yu(nn+1) = dyp * (nn + .5d0) + 1d0 - nn * dyp !+ (1d0 + sizep) - dyp * (nn + qp + 0.5d0)
   dthdyu(nn+1) = 1d0 / dyp

   print*, "dthdyu = ", dthdyu(0), " ", dthdyu(nn+1)
   print*, "dthdyv = ", dthdyv(0), " ", dthdyv(nn)
   print*, "Shape u, v", shape(dthdyu), ", ", shape(dthdyv)

!-----------------------------------------------------------------------
! Second-order finite difference coefficient for 1. u(j-1), 2. u(j), 3. u(j+1)
!-----------------------------------------------------------------------

  !NOT USED? do we have access to nyu21... now?
!  dyub2 = ((yu(nyu21+1)-yu(nyu21))*(yu(nyu21+1)-yu(nyu21-1)))
!  dyut2 = ((yu(nyu12)-yu(nyu12-1))*(yu(nyu12+1)-yu(nyu12-1)))

!  dyvb2 = ((yv(nyv21+1)-yv(nyv21))*(yv(nyv21+1)-yv(nyv21-1)))
!  dyvt2 = ((yv(nyv12)-yv(nyv12-1))*(yv(nyv12+1)-yv(nyv12-1)))
  !------

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

  ! Attention! nn changes
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
  N(4,2) = N(3,2)   !Band 2-3 
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


open(10,file=trim(dirout)//'y_grid_substrate.dat',form='unformatted',access='stream')
write(*,*) 'printing ', trim(dirout)//'y_grid_substrate.dat'
write(10) N, dny
write(10) yu, yv, dthdyu, dthdyv
write(10) dyu2i, dyv2i
close(10)


end subroutine

subroutine y_grid_canopy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!   YGRID for canopies  !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Creates the geometry in the y-direction
!  ! y(j) = a(j-q) + b(j-q)**p
  
!  y is the physical coordinate
!  theta is the mapping coordinate (constant increments: dtheta)
!  Since (1) dtheta is constant and (2) the mapping between y and
!  theta is analytical,
!  the second order of the centered finite difference is preserved.
  
!  n   is the total number of point from wall to wall
!  q   is the index of the centerline (n/2) (y(q) = 0)
!  dyq is dymax/dymin (max stretching)
!  p   polynomial exponent
!  a and b are obtained such that dy(nb)/dy(q)=dyq and y(0)=-h_canopy

   use declaration
   implicit none
   integer j,iband, dny
   real(8) aaa, bbb, aat, aab, dy1, dye, ccc, qqq
  
   dny = dsty;
   ! v points in y go from 0 at the bottom interface to nn at the top interface
   !Attention: the value of nn changes along the subroutine. This should be changed
   nn  = N(3,nband)+1
   print*, "nn = ", nn
   qqq = nn/2d0
   
   aaa = 15d0*(dyq - 1d0)/((8d0*dyq + 7d0)*(qqq**ppp))
   bbb =-30d0*(dyq - 1d0)/((8d0*dyq + 7d0)*(qqq**(ppp-2)))
   ccc = 15d0*dyq/((8d0*dyq + 7d0)*qqq)
   
   dy1 = aaa*(-qqq)**4 + bbb*(-qqq)**2 + ccc
   
   aab = (posth-dy1*dny)/(dny**3)
  
   dtheta     = 1d0
   dthetai    = 1d0/dtheta
   ddthetavi  = 1.0d0/dtheta

  
   allocate(yu    (  -dny:nn+dny+1))
   allocate(dyu2i (3,-dny:nn+dny+1))
   allocate(dthdyu(  -dny:nn+dny+1))
  
   allocate(yv    (  -dny:nn+dny  ))
   allocate(dyv2i (3,-dny:nn+dny  ))
   allocate(dthdyv(  -dny:nn+dny  ))


 !!! Grid between tips
   !Q: shouldn't be it nn? yv(nn)? 
   !Q: what about yu(0)?
   do j=0,nn-1
     yv(j)= aaa*(j       -qqq)**ppp/5d0 + bbb*(j     -qqq)**(ppp-2)/3d0 + ccc*(j     -qqq)
     yu(j+1)= aaa*(j+.5d0-qqq)**ppp/5d0 + bbb*(j+.5d0-qqq)**(ppp-2)/3d0 + ccc*(j+.5d0-qqq)
   
     dthdyv(j) = 1d0/(aaa*(j       -qqq)**(ppp-1)+bbb*(j     -qqq)**(ppp-3)+ccc)
     dthdyu(j+1) = 1d0/(aaa*(j+.5d0-qqq)**(ppp-1)+bbb*(j+.5d0-qqq)**(ppp-3)+ccc)

   end do

   !Q: and last point? yv(nn) and yu(nn+1)
   if (dny == 0) then

      yu(0) = aaa*(-0.5d0-qqq)**ppp/5d0+bbb*(-0.5d0-qqq)**(ppp-2)/3d0+ccc*(-0.5d0-qqq)
      dthdyu(0) = 1d0/(aaa*(-.5d0-qqq)**(ppp-1)+bbb*(-.5d0-qqq)**(ppp-3)+ccc)
   
   else

!!!! Grid inside canopies (CHANGE FOR POROUS)
 
      yu(-dny)    = aab*(-dny-.5d0)**3 + dy1*(-dny-0.5d0) -1d0
      dthdyu(-dny)= 1d0/(3d0*aab*(-dny-.5d0)**2 + dy1) 
      do j=-dny,-1
        yv(j)  = aab*(j     )**3 + dy1*(j     ) - 1d0
        yu(j+1)= aab*(j+.5d0)**3 + dy1*(j+.5d0) - 1d0

        dthdyv(j)  = 1d0/(3d0*aab*(j     )**2 + dy1) 
        dthdyu(j+1)= 1d0/(3d0*aab*(j+.5d0)**2 + dy1)
      end do
         
	  !Q: yv(nn) and yu(nn+1) are here
      do j=0,dny
        yv(j+nn)  = aab*(j     )**3 + dy1*(j     ) +1d0
        yu(j+nn+1)= aab*(j+.5d0)**3 + dy1*(j+.5d0) +1d0

        dthdyv(j+nn)  = 1d0/(3d0*aab*(j     )**2 + dy1) 
        dthdyu(j+nn+1)= 1d0/(3d0*aab*(j+.5d0)**2 + dy1)
      end do
   end if

   print*, "dthdyu = ", dthdyu(0), " ", dthdyu(nn+1)
   print*, "Shape u, v", shape(dthdyu), ", ", shape(dthdyv)


! Second-order finite difference coefficient for 1. u(j-1), 2. u(j), 3. u(j+1)

  !Q: USED? do we have access to nyu21... now?
  dyub2 = ((yu(nyu21+1)-yu(nyu21))*(yu(nyu21+1)-yu(nyu21-1)))
  dyut2 = ((yu(nyu12)-yu(nyu12-1))*(yu(nyu12+1)-yu(nyu12-1)))
        write(*,*) 'dyub2, dyut2 =', dyub2, dyut2 

  dyvb2 = ((yv(nyv21+1)-yv(nyv21))*(yv(nyv21+1)-yv(nyv21-1)))
  dyvt2 = ((yv(nyv12)-yv(nyv12-1))*(yv(nyv12+1)-yv(nyv12-1)))
        write(*,*) 'dyvb2, dyvt2 =', dyvb2, dyvt2
  !------

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
  N(4,2) = N(3,2)   !Band 2-3 
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


open(10,file=trim(dirout)//'y_grid_canopy.dat',form='unformatted',access='stream')
write(*,*) 'printing ', trim(dirout)//'y_grid_canopy.dat'
write(10) N, dny
write(10) yu, yv, dthdyu, dthdyv
write(10) dyu2i, dyv2i
close(10)


end subroutine


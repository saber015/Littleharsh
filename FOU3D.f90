! *************************************************************************************************************** !
!
! Last modified by C 12/05/2016
!   - Set up for SHS
!
! TODO
!   - Linear interpolation for interp?...
!
! *************************************************************************************************************** !
!
! Contains: nonlinear  - Called from newribs : Several called functions within file
!                        - Solves the nonlinear advective term
!                        - Records everything to file
!                        - Implementation of immersed boundaries
!           interp_u/v - Called from nonlinear
!                        - Interperlates the grids onto each other to calculate the nonlinear term
!                          - (linear interpolarion)
!           der_x/z    - Called from nonlinear
!                        - Calcualtes the x/z derrivative
!           der_yu/v_h - Called from nonlinear
!                        - Calcualtes the y derrivative
!           dtc_calc   - Called form nonlinear
!                        - Calculates the convective(?) timestep
!           imm_bounds - Called from nonlinear
!                        - Immersed boundaries
!           record_out - Called from nonlinear
!                        - Records to file
!
! Also includes four_to_phys_.. for IFFT
!   phys_to_four_.. for FFT
!   modes_to_planes_.. for modes to planes
!   planes_to_modes_.. for planes to modes
!
!
! *************************************************************************************************************** !

! ****** Modified C 31/08/2015 ****** !
! Advective terms now calculated in conservation form
!	Might not actually be quicker when using immersed boundaries

subroutine nonlinear(Nu1,Nu2,Nu3,u1,u2,u3,du1,du2,du3,p,div,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!   NONLINEAR TERMS  !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid
  
  integer iband,flagst,flagwr,flagslinst,flagqwr,j,i,k,column
  real(8) C1
  type(cfield)  u1(sband:eband), u2(sband:eband), u3(sband:eband)
  type(cfield) du1(sband:eband),du2(sband:eband),du3(sband:eband)
  type(cfield) Nu1(sband:eband),Nu2(sband:eband),Nu3(sband:eband)
  type(cfield)  p (sband:eband),div(sband:eband)
  
  if (iter-iter0>=nstat .and. kRK==1) then
    flagst = 1
    iter0  = iter
  else
    flagst = 0
  end if

  if (iter>=iwrite .and. kRK==1) then
    flagwr = 1
    iwrite = iwrite+nwrite
  else
    flagwr = 0
  end if
  
!   if (t>=nextqt) then
!     flagqwr = 1
!     nextqt = nextqt+10.0d0
!   else
!     flagqwr = 0
!   end if
  
  u1PL_itp = 0d0
  u2PL_itp = 0d0
  u3PL_itp = 0d0

  !C! Interpolate the grid velocities to the other grid points
  call interp_u(u1_itp,u1,myid)
  call interp_v(u2_itp,u2,myid)
  call interp_u(u3_itp,u3,myid) 

  u1PL  = 0d0
  u2PL  = 0d0
  u3PL  = 0d0
  Nu1PL = 0d0
  Nu2PL = 0d0
  Nu3PL = 0d0
  du1PL = 0d0
  du2PL = 0d0
  du3PL = 0d0

  !C! Shift 6 velocity fields into planes
  !!!!!!!!!!  modes to planes: !!!!!!!!!!
  call modes_to_planes_UVP ( u1PL,    u1,    ugrid,myid,status,ierr)
  call modes_to_planes_UVP ( u2PL,    u2,    vgrid,myid,status,ierr)
  call modes_to_planes_UVP ( u3PL,    u3,    ugrid,myid,status,ierr)
  call modes_to_planes_UVP ( u1PL_itp,u1_itp,vgrid,myid,status,ierr)
  call modes_to_planes_UVP ( u2PL_itp,u2_itp,ugrid,myid,status,ierr)
  call modes_to_planes_UVP ( u3PL_itp,u3_itp,vgrid,myid,status,ierr)

  !!!!!!!!!!!!!   spectra:  !!!!!!!!!!!!!
  if (flagst==1) then
    call spectra(u1,u2,u2_itp,u3,p,myid)
  end if
  
!   !!!!!!!!! Q criterion  !!!!!!!!!
!   if (flagqwr==1) then !Q criterion output every t
!     
!     do iband = sband,eband
!       call der_yv_h_wx(du1dy_columns(iband)%f,u1(iband)%f,iband,myid)
!       call der_yv_h_wx(du2dy_columns(iband)%f,u2_itp(iband)%f,iband,myid)
!       call der_yv_h_wx(du3dy_columns(iband)%f,u3(iband)%f,iband,myid)
!     end do
!     
!     call modes_to_planes_UVP(du1dy_planes2,du1dy_columns,vgrid,myid,status,ierr)
!     call modes_to_planes_UVP(du2dy_planes2,du2dy_columns,vgrid,myid,status,ierr)
!     call modes_to_planes_UVP(du3dy_planes2,du3dy_columns,vgrid,myid,status,ierr)
! 
!     call Qcriterion(myid)  
!     
!   endif

  !!!!!!!!!!!!! record out: !!!!!!!!!!!!!
  if (flagwr==1) then
    call error(div,myid,ierr)
    ppPL = 0d0
!    div_outPL = 0d0
    call modes_to_planes_UVP(ppPL,p,3,myid,status,ierr)
!    call modes_to_planes_UVP(div_outPL,div_out,3,myid,status,ierr)
!    call modes_to_planes_UVP(div_cPL,  div,  3,myid,status,ierr)
!    call modes_to_planes_UVP(u_outPL,u_out,ugrid,myid,status,ierr)
!    call modes_to_planes_UVP(v_outPL,v_out,vgrid,myid,status,ierr)
!    call modes_to_planes_UVP(w_outPL,w_out,ugrid,myid,status,ierr)
    !call modes_to_planes_UVP(ppPL,div,3,myid,status,ierr) !Output divergence for checking
    call record_out(u1,myid)
  end if

  !!!!!!!!! four to ops: !!!!!!!!!
  call ops_in_planes(myid,flagst) !C! ops in planes to compute velocity products and x/z derriatives

  !C! Shift y derrivative products to modes
  call planes_to_modes_UVP(uv_f,uv_fPL,vgrid,myid,status,ierr)
  call planes_to_modes_UVP(vv_c,vv_cPL,ugrid,myid,status,ierr)
  call planes_to_modes_UVP(vw_f,vw_fPL,vgrid,myid,status,ierr)

  !C! Shift x/z derrivatives to modes
  call planes_to_modes_NUVP(Nu1,Nu1PL,ugrid,myid,status,ierr)
  call planes_to_modes_NUVP(Nu2,Nu2PL,vgrid,myid,status,ierr)
  call planes_to_modes_NUVP(Nu3,Nu3PL,ugrid,myid,status,ierr)

  !C! Calculate y derrivatives
   do iband = sband,eband
     call der_yu_h(Nu1_dy(iband)%f,uv_f(iband)%f,iband,myid)
     call der_yv_h(Nu2_dy(iband)%f,vv_c(iband)%f,iband,myid)
     call der_yu_h(Nu3_dy(iband)%f,vw_f(iband)%f,iband,myid)

   !C! Calculate final advective term
     do column = 1,columns_num(iband,myid)
       do j = jlim(1,ugrid,iband)+1,jlim(2,ugrid,iband)-1
         Nu1(iband)%f(j,column) = Nu1(iband)%f(j,column)+Nu1_dy(iband)%f(j,column)
         Nu3(iband)%f(j,column) = Nu3(iband)%f(j,column)+Nu3_dy(iband)%f(j,column)
       end do
       do j = jlim(1,vgrid,iband)+1,jlim(2,vgrid,iband)-1
         Nu2(iband)%f(j,column) = Nu2(iband)%f(j,column)+Nu2_dy(iband)%f(j,column)
       enddo
     enddo
   enddo   

!   ! Low pass filter function ! 20.4.2020
!   do iband = sband, eband
!      do column = 1, columns_num(iband,myid)
!         K_lpf = exp(xi2*(columns_i(column,iband,myid)*1d0/(N(1,iband)/2))**xi1) &
!             & * exp(xi2*(columns_k(column,iband,myid)*1d0/maxval(columns_k(:,iband,:)))**xi1)
!         do j = jlim(1,ugrid,iband)+1, jlim(2,ugrid,iband)-1
!            Nu1(iband)%f(j,column) = K_lpf * Nu1(iband)%f(j,column)
!            Nu3(iband)%f(j,column) = K_lpf * Nu3(iband)%f(j,column)
!         end do
!         do j = jlim(1,vgrid,iband)+1, jlim(2,vgrid,iband)-1
!            Nu2(iband)%f(j,column) = K_lpf * Nu2(iband)%f(j,column)
!         end do
!      end do
!   end do

  !!!!!!!!!!!  CFL and stats: !!!!!!!!!!!
  if (kRK==1) then
  ! Calculating timestep
    call dtc_calc(u1PL,u2PL,u3PL,bandPL(myid),myid)
    call MPI_ALLREDUCE(dt,dtc,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
    dt  = min(2.5d0*dtv,CFL*dtc)
    dti = 1d0/dt
    t   = t+dt
    ! Calculating and writing stats and spectra
    if (flagst==1) then
    
    !C! Calculate Omega_x
    
    du3dy_planes2 = 0d0
    call modes_to_planes_UVP (du3dy_planes2,u2,vgrid,myid,status,ierr)   

    do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      call der_z(du3dy_planes2(1,1,j),wx(:,:,j),k1F_z,bandPL(myid)) !C! u2PL in Fourier space
    end do
    do iband = sband,eband
      call der_yv_h_wx(du3dy_columns(iband)%f,u3(iband)%f,iband,myid)
    end do
    
    call modes_to_planes_UVP(du3dy_planes2,du3dy_columns,vgrid,myid,status,ierr)
    
    do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      do i = 1,igal
        do k = 1,kgal
          wx(i,k,j) = -wx(i,k,j)+du3dy_planes2(i,k,j) !omega_x at faces!
        end do
      end do
      call four_to_phys_du(wx(1,1,j),bandPL(myid))
    end do

    
      ppPL = 0d0
      call modes_to_planes_UVP (ppPL,p,pgrid,myid,status,ierr)   
      do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
	call four_to_phys_du(ppPL(1,1,j),bandPL(myid))
      enddo

      call stats(myid,status,ierr) 

      u2PL=0d0

    end if
    
    if (flagwr==1) then
      call write_spect(myid,status,ierr) 
      call write_stats(myid,status,ierr) 
    end if    

        
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Update RHS with new advective term. du now full RHS
  C1 = -gRK(kRK)
  do iband = sband,eband
    do column = 1,columns_num(iband,myid)
      do j = jlim(1,ugrid,iband)+1,jlim(2,ugrid,iband)-1
        du1(iband)%f(j,column) = du1(iband)%f(j,column)+C1*Nu1(iband)%f(j,column)
        du3(iband)%f(j,column) = du3(iband)%f(j,column)+C1*Nu3(iband)%f(j,column)
      enddo
      do j = jlim(1,vgrid,iband)+1,jlim(2,vgrid,iband)-1
        du2(iband)%f(j,column) = du2(iband)%f(j,column)+C1*Nu2(iband)%f(j,column)
      enddo
    enddo
  enddo

  !C! Calculate immersed boundaries
!   call modes_to_planes_dU(du1PL,du1,myid,status,ierr)
!   call modes_to_planes_dV(du2PL,du2,myid,status,ierr)
!   call modes_to_planes_dU(du3PL,du3,myid,status,ierr)
!   
!   call imm_bounds_u(du1PL,u1PL,Nu1PL,myid,status,ierr)
!   call imm_bounds_v(du2PL,u2PL,Nu2PL,myid,status,ierr)
!   call imm_bounds_u(du3PL,u3PL,Nu3PL,myid,status,ierr)
! 
!   call planes_to_modes_dU(du1,du1PL,myid,status,ierr)
!   call planes_to_modes_dV(du2,du2PL,myid,status,ierr)
!   call planes_to_modes_dU(du3,du3PL,myid,status,ierr)

end subroutine

subroutine interp_u(u_itp,u,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!  interp u/w grid to v grid  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use declaration
  implicit none
  
  integer j,iband,column,myid
  type(cfield)  u(sband:eband)
  type(cfield)  u_itp(sband:eband)
  
 do iband = sband,eband
   do column = 1,columns_num(iband,myid)
     ! We interpolate everything. vgrid has got one less point than ugrid
     do j = jlim(1,vgrid,iband),jlim(2,vgrid,iband)  
       u_itp(iband)%f(j,column) = ((yv(j)-yu(j))*u(iband)%f(j+1,column)+(yu(j+1)-yv(j))*u(iband)%f(j,column))/(yu(j+1)-yu(j))
     end do
   end do
 end do
  
end subroutine

subroutine interp_v(u_itp,u,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!  interp v grid to u/w grid  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use declaration
  implicit none
  
  integer j,iband,column,myid
  type(cfield)  u(sband:eband)
  type(cfield)  u_itp(sband:eband)
  
 do iband = sband,eband
   do column = 1,columns_num(iband,myid)
     do j = jlim(1,ugrid,iband)+1,jlim(2,ugrid,iband)-1
       u_itp(iband)%f(j,column) = ((yu(j)-yv(j-1))*u(iband)%f(j,column)+(yv(j)-yu(j))*u(iband)%f(j-1,column))/(yv(j)-yv(j-1))
     end do
     u_itp(iband)%f(jlim(1,ugrid,iband),column) = &
&      gridweighting_interp(iband,1)*(u(iband)%f(jlim(1,vgrid,iband),column)-u_itp(iband)%f(jlim(1,ugrid,iband)+1,column)) &
&      + u_itp(iband)%f(jlim(1,ugrid,iband)+1,column)
     u_itp(iband)%f(jlim(2,ugrid,iband),column) = &
&      gridweighting_interp(iband,2)*(u(iband)%f(jlim(2,vgrid,iband),column)-u_itp(iband)%f(jlim(2,ugrid,iband)-1,column)) &
&      + u_itp(iband)%f(jlim(2,ugrid,iband)-1,column)
   end do
 end do
                  
end subroutine

subroutine der_x(u,dudx,kx,iband)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!  NONLINEAR DER TERMS  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! iband is 'bandPL(myid)' (physical space), however the variables are in Fourier space

  use declaration
  implicit none

  integer i,k,iband
  complex(8) u   (0:Ngal(1,iband)/2,Ngal(2,iband))
  complex(8) dudx(0:Ngal(1,iband)/2,Ngal(2,iband))
  complex(8) kx(0:N(1,1)/2)

  do k = 1,N(2,iband)/2
    do i = 0,N(1,iband)/2
      dudx(i,k) = kx(i)*u(i,k)
    end do
    do i = N(1,iband)/2+1,Ngal(1,iband)/2           !!!!!!!!!!!!  Zeros for the antialiasing region (x-dir)
      dudx(i,k) = 0d0                               !!!!!!!!!!!!  must be explicitly specified
    end do
  end do
  !Zeros includes mode Nz/2+1 mode (zero for advection derrivatives)
  do k = N(2,iband)/2+1,Ngal(2,iband)-N(2,iband)/2+1!!!!!!!!!!!!  Zeros for the antialiasing region (z-dir)
    do i = 0,Ngal(1,iband)/2                        !!!!!!!!!!!!
      dudx(i,k) = 0d0                               !!!!!!!!!!!!  must be explicitly specified
    end do
  end do
  do k = Ngal(2,iband)-N(2,iband)/2+2,Ngal(2,iband)
    do i = 0,N(1,iband)/2
      dudx(i,k) = kx(i )*u(i,k)
    end do
    do i = N(1,iband)/2+1,Ngal(1,iband)/2           !!!!!!!!!!!!  Zeros for the antialiasing region (x-dir)
      dudx(i,k) = 0d0                               !!!!!!!!!!!!  must be explicitly specified
    end do
  end do

end subroutine

subroutine der_x_N(u,dudx,kx,iband) !For N
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!  NONLINEAR DER TERMS  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! iband is 'bandPL(myid)' (physical space), however the variables are in Fourier space

  use declaration
  implicit none

  integer i,k,iband,dk2,k2
  complex(8) u   (0:N(1,iband)/2,N(2,iband))
  complex(8) dudx(0:N(1,iband)/2,N(2,iband))
  complex(8) kx(0:N(1,1)/2)

  do k = 1,N(2,iband)
    do i = 0,N(1,iband)/2
      dudx(i,k) = kx(i)*u(i,k)
    end do
  end do
  
end subroutine

subroutine der_z(u,dudz,kz,iband)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!  NONLINEAR DER TERMS  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! iband is 'bandPL(myid)' (physical space), however the variables are in Fourier space

  use declaration
  implicit none

  integer i,k,iband,dk2,k2
  complex(8) u   (0:Ngal(1,iband)/2,Ngal(2,iband))
  complex(8) dudz(0:Ngal(1,iband)/2,Ngal(2,iband))
  complex(8) kz(1:N(2,1))

  dk2=N(2,1)-Ngal(2,iband)

  do k = 1,N(2,iband)/2
    do i = 0,N(1,iband)/2
      dudz(i,k) = kz(k)*u(i,k)
    end do
    do i = N(1,iband)/2+1,Ngal(1,iband)/2           !!!!!!!!!!!!  Zeros for the antialiasing region (x-dir)
      dudz(i,k) = 0d0
    end do
  end do
  !Zeros includes mode Nz/2+1 mode (zero for advection derrivatives)
  do k = N(2,iband)/2+1,Ngal(2,iband)-N(2,iband)/2+1!!!!!!!!!!!!  Zeros for the antialiasing region (z-dir)
    do i = 0,Ngal(1,iband)/2
      dudz(i,k) = 0d0
    end do
  end do
  do k = Ngal(2,iband)-N(2,iband)/2+2,Ngal(2,iband)
    k2 = k+dk2
    do i = 0,N(1,iband)/2
      dudz(i,k) = kz(k2)*u(i,k)
    end do
    do i = N(1,iband)/2+1,Ngal(1,iband)/2           !!!!!!!!!!!!  Zeros for the antialiasing region (x-dir)
      dudz(i,k) = 0d0
    end do
  end do

end subroutine

subroutine der_z_N(u,dudz,kz,iband)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!  NONLINEAR DER TERMS  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! iband is 'bandPL(myid)' (physical space), however the variables are in Fourier space

  use declaration
  implicit none

  integer i,k,iband,dk2,k2
  complex(8) u   (0:N(1,iband)/2,N(2,iband))
  complex(8) dudz(0:N(1,iband)/2,N(2,iband))
  complex(8) kz(1:N(2,1))

  if(iband/=2)then
    do k = 1,N(2,iband)
      do i = 0,N(1,iband)/2
	dudz(i,k) = kz(k)*u(i,k)
      end do
    end do
  else
    do k = 1,N(2,iband)/2
      do i = 0,N(1,iband)/2
	dudz(i,k) = kz(k)*u(i,k)
      end do
    end do
    dk2=N(2,1)-N(2,2)
    do k = N(2,iband)/2+1,N(2,iband)
      do i = 0,N(1,iband)/2
	dudz(i,k) = kz(k+dk2)*u(i,k)
      end do
    end do
  endif

end subroutine

subroutine der_yu_h(dudy,u,iband,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      DER Y     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      f-->c     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer j,column,iband,myid
  complex(8)  u   (jlim(1,vgrid,iband)  :jlim(2,vgrid,iband)  ,columns_num(iband,myid))
  complex(8) dudy (jlim(1,ugrid,iband)+1:jlim(2,ugrid,iband)-1,columns_num(iband,myid))

  do column = 1,columns_num(iband,myid)
    do j = jlim(1,ugrid,iband)+1,jlim(2,ugrid,iband)-1
      dudy(j,column) = (u(j,column)-u(j-1,column))*dthdyu(j)*ddthetavi
    end do
  end do

end subroutine

subroutine der_yu_h_write(dudy,u,iband,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      DER Y     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      f-->c     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer j,column,iband,myid
  complex(8)  u   (jlim(1,vgrid,iband)  :jlim(2,vgrid,iband)  ,columns_num(iband,myid))
  complex(8) dudy (jlim(1,ugrid,iband)+1:jlim(2,ugrid,iband)-1,columns_num(iband,myid))

if (myid .eq. 0 .and. iband .eq. sband) then
open(unit=9,file='debug_der_yu_h.dat')
rewind 9
write(9,*) 'dudy_jlim =', jlim(1,ugrid,iband)+1, jlim(2,ugrid,iband)-1

  do column = 1,columns_num(iband,myid)
    do j = jlim(1,ugrid,iband)+1,jlim(2,ugrid,iband)-1
      dudy(j,column) = (u(j,column)-u(j-1,column))*dthdyu(j)*ddthetavi
write(9,*) 'column, j =', column, j
write(9,*) '   u(j-1) =', u(j-1,column)
write(9,*) '     dudy =', dudy(j,column)
    end do
  end do

close(unit=9)
end if

end subroutine

subroutine der_yv_h(dudy,u,iband,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      DER Y     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      c-->f     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer j,column,iband,myid
  complex(8)  u   (jlim(1,ugrid,iband)  :jlim(2,ugrid,iband)  ,columns_num(iband,myid))
  complex(8) dudy (jlim(1,vgrid,iband)+1:jlim(2,vgrid,iband)-1,columns_num(iband,myid))

  do column = 1,columns_num(iband,myid)
    do j = jlim(1,vgrid,iband)+1,jlim(2,vgrid,iband)-1
      dudy(j,column) = (u(j+1,column)-u(j,column))*dthdyv(j)*ddthetavi
    end do
  end do

end subroutine

subroutine der_yv_h_wx(dudy,u,iband,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      DER Y     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!      c-->f     !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  integer j,column,iband,myid
  complex(8)  u   (jlim(1,ugrid,iband):jlim(2,ugrid,iband),columns_num(iband,myid))
  complex(8) dudy (jlim(1,vgrid,iband):jlim(2,vgrid,iband),columns_num(iband,myid))

  do column = 1,columns_num(iband,myid)
    do j = jlim(1,vgrid,iband),jlim(2,vgrid,iband)
      dudy(j,column) = (u(j+1,column)-u(j,column))*dthdyv(j)*ddthetavi
    end do
  end do

end subroutine

subroutine dtc_calc(u1,u2,u3,iband,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!       CFL      !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer i,k,j,iband,myid
  real(8) u1mloc,u2mloc,u3mloc
  real(8) u1(igal,kgal,jgal(2,1)-1:jgal(2,2)+1)
  real(8) u2(igal,kgal,jgal(1,1)-1:jgal(1,2)+1)
  real(8) u3(igal,kgal,jgal(2,1)-1:jgal(2,2)+1)

  u1mloc = maxval(abs(u1))
  u2mloc = 1d-6

  do j = jgal(1,1),jgal(1,2)
    do k = 1,kgal
      do i = 1,igal-2
        u2mloc = max(u2mloc,abs(u2(i,k,j))*dthdyv(j))
      end do
    end do
  end do

  u3mloc = maxval(abs(u3))
  
  dt = min(1d0/(alp*(N(1,iband)/2)*u1mloc),1d0/(bet*(N(2,iband)/2)*u3mloc),1d0/(u2mloc*dthetavi)) ! 9.3.2020

end subroutine

subroutine four_to_phys_u(u1,u2,u3,iband)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!      FOUR TO PHYS     !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Transforms from Fourier Space to Physical Space u and its derivatives
! Its only used in ops in planes at FOU3D.f90

  use declaration
  implicit none
  integer iband
  real(8) u1 (Ngal(1,iband)+2,Ngal(2,iband))
  real(8) u2 (Ngal(1,iband)+2,Ngal(2,iband))
  real(8) u3 (Ngal(1,iband)+2,Ngal(2,iband))
  
  u1(:,Ngal(2,iband)/2+1)=0d0
  u2(:,Ngal(2,iband)/2+1)=0d0
  u3(:,Ngal(2,iband)/2+1)=0d0
  
  call cft(u1,Ngal(1,iband)+2,2,(N(1,iband)+2)/2,1,buffCal_z(iband)%b)
  call rft(u1,Ngal(1,iband)+2,Ngal(2,iband),1,buffRal_x(iband)%b)
  call cft(u2,Ngal(1,iband)+2,2,(N(1,iband)+2)/2,1,buffCal_z(iband)%b)
  call rft(u2,Ngal(1,iband)+2,Ngal(2,iband),1,buffRal_x(iband)%b)
  call cft(u3,Ngal(1,iband)+2,2,(N(1,iband)+2)/2,1,buffCal_z(iband)%b)
  call rft(u3,Ngal(1,iband)+2,Ngal(2,iband),1,buffRal_x(iband)%b)
  
end subroutine

subroutine four_to_phys_du(du,iband)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!      FOUR TO PHYS     !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Transforms from Fourier Space to Physical Space u and its derivatives
! Its only used in ops in planes at FOU3D.f90

  use declaration
  implicit none
  integer iband
  real(8) du  (Ngal(1,iband)+2,Ngal(2,iband))

  du(:,Ngal(2,iband)/2+1)=0d0
  
  call cft(du   ,Ngal(1,iband)+2,2,(N(1,iband)+2)/2,1,buffCal_z(iband)%b)
  call rft(du   ,Ngal(1,iband)+2,Ngal(2,iband),1,buffRal_x(iband)%b)

end subroutine

subroutine four_to_phys_N(du,iband)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!      FOUR TO PHYS     !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Transforms from Fourier Space to Physical Space u and its derivatives
! Its only used in ops in planes at FOU3D.f90

  use declaration
  implicit none
  integer iband
  real(8) du  (N(1,iband)+2,N(2,iband))

  !du(:,N(2,iband)/2+1)=0d0
  
  call cft(du   ,N(1,iband)+2,2,(N(1,iband)+2)/2,1,buffC_z(iband)%b)
  call rft(du   ,N(1,iband)+2,N(2,iband),1,buffR_x(iband)%b)

end subroutine


subroutine modes_to_planes_UVP (xPL,x,grid,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! MODES TO PLANES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,dki,plband,grid
  integer column
  integer iband,jband
  integer inode,yourid
  integer msizeR,msizeS
  type(cfield) x  (sband:eband)
  real(8)      xPL(igal,kgal,jgal(grid,1)-1:jgal(grid,2)+1)
  complex(8), allocatable :: buffS(:,:),buffR(:,:)

  ! Loop for itself
  ! Transpose the cube that it already owns
  plband = bandPL(myid)
  yourid = myid
  do iband = sband,eband
    jband = iband
    jminR = max(planelim(grid,1,  myid),jlim(1,grid,jband)+1)
    jmaxR = min(planelim(grid,2,  myid),jlim(2,grid,jband)-1)
    if (jminR==Ny(grid,plband-1)+1 .and. jmaxR>=jminR) then
      jminR = jminR-1
    end if
    if (jmaxR==Ny(grid,plband  )   .and. jmaxR>=jminR) then
      jmaxR = jmaxR+1
    end if
    do j = jminR,jmaxR
      do column = 1,columns_num(jband,yourid)
        i = columns_i(column,jband,yourid)
        k = columns_k(column,jband,yourid) - dk(column,jband,yourid,bandPL(myid))
        xPL(2*i+1,k,j) = dreal(x(jband)%f(j,column))
        xPL(2*i+2,k,j) = dimag(x(jband)%f(j,column))
      end do
    end do
  end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)
    if (yourid<np) then
      do iband = sband,eband
        !jband=crossband(iband,yourid)
        jband = iband
        jminS = max(planelim(grid,1,yourid),jlim(1,grid,iband)+1)
        jmaxS = min(planelim(grid,2,yourid),jlim(2,grid,iband)-1)
        jminR = max(planelim(grid,1,  myid),jlim(1,grid,jband)+1)
        jmaxR = min(planelim(grid,2,  myid),jlim(2,grid,jband)-1)
        if (jminS==Ny(grid,0)+1  ) then
          jminS = jminS-1
        end if
        if (jmaxS==Ny(grid,nband)) then
          jmaxS = jmaxS+1
        end if
        if (jminR==Ny(grid,0)+1  ) then
          jminR = jminR-1
        end if
        if (jmaxR==Ny(grid,nband)) then
          jmaxR = jmaxR+1
        end if
        allocate(buffS(jminS:jmaxS,columns_num(iband,  myid)))
        allocate(buffR(jminR:jmaxR,columns_num(jband,yourid)))
        msizeS = 2*(columns_num(iband,  myid)*(jmaxS-jminS+1))  ! 2 times because it's complex
        msizeR = 2*(columns_num(jband,yourid)*(jmaxR-jminR+1))
        msizeS = max(msizeS,0)
        msizeR = max(msizeR,0)

        do j = jminS,jmaxS
          do column = 1,columns_num(iband,myid)
            buffS(j,column) = x(iband)%f(j,column)

          end do
        end do

        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*iband+11*jband, &
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*iband+7*jband, &
&                         MPI_COMM_WORLD,status,ierr)
        do j = jminR,jmaxR
          do column = 1,columns_num(jband,yourid)
            i = columns_i(column,jband,yourid)
            k = columns_k(column,jband,yourid) - dk(column,jband,yourid,bandPL(myid))
            xPL(2*i+1,k,j) = dreal(buffR(j,column))
            xPL(2*i+2,k,j) = dimag(buffR(j,column))
          end do
        end do

        deallocate(buffR,buffS)

      end do
    end if
  end do

end subroutine

subroutine modes_to_planes_phys (xPL,x,grid,myid,myiband,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! MODES TO PLANES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,dki,plband,grid,myiband
  integer column
  integer iband,jband
  integer inode,yourid
  integer msizeR,msizeS
  type(cfield) x  (sband:eband)
  real(8)      xPL(N(1,myiband)+2,N(2,myiband),jgal(grid,1)-1:jgal(grid,2)+1)
  complex(8), allocatable :: buffS(:,:),buffR(:,:)

  plband = bandPL(myid)
  yourid = myid
  do iband = sband,eband
    jband = iband
    jminR = max(planelim(grid,1,  myid),jlim(1,grid,jband)+1)
    jmaxR = min(planelim(grid,2,  myid),jlim(2,grid,jband)-1)
    if (jminR==Ny(grid,0)+1  ) then
          jminR = jminR-1
        end if
        if (jmaxR==Ny(grid,nband)) then
          jmaxR = jmaxR+1
        end if
    do j = jminR,jmaxR
      do column = 1,columns_num(jband,yourid)
        i = columns_i(column,jband,yourid)
        k = columns_k(column,jband,yourid) - dk_phys(column,jband,yourid,bandPL(myid))
        xPL(2*i+1,k,j) = dreal(x(jband)%f(j,column))
        xPL(2*i+2,k,j) = dimag(x(jband)%f(j,column))
      end do
    end do
  end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)
    if (yourid<np) then
      do iband = sband,eband
        !jband=crossband(iband,yourid)
        jband = iband
        jminS = max(planelim(grid,1,yourid),jlim(1,grid,iband)+1)
        jmaxS = min(planelim(grid,2,yourid),jlim(2,grid,iband)-1)
        jminR = max(planelim(grid,1,  myid),jlim(1,grid,jband)+1)
        jmaxR = min(planelim(grid,2,  myid),jlim(2,grid,jband)-1)
        if (jminS==Ny(grid,0)+1  ) then
          jminS = jminS-1
        end if
        if (jmaxS==Ny(grid,nband)) then
          jmaxS = jmaxS+1
        end if
        if (jminR==Ny(grid,0)+1  ) then
          jminR = jminR-1
        end if
        if (jmaxR==Ny(grid,nband)) then
          jmaxR = jmaxR+1
        end if
        allocate(buffS(jminS:jmaxS,columns_num(iband,  myid)))
        allocate(buffR(jminR:jmaxR,columns_num(jband,yourid)))
        msizeS = 2*(columns_num(iband,  myid)*(jmaxS-jminS+1))  ! 2 times because it's complex
        msizeR = 2*(columns_num(jband,yourid)*(jmaxR-jminR+1))
        msizeS = max(msizeS,0)
        msizeR = max(msizeR,0)

        do j = jminS,jmaxS
          do column = 1,columns_num(iband,myid)
            buffS(j,column) = x(iband)%f(j,column)

          end do
        end do

        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*iband+11*jband, &
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*iband+7*jband, &
&                         MPI_COMM_WORLD,status,ierr)
        do j = jminR,jmaxR
          do column = 1,columns_num(jband,yourid)
            i = columns_i(column,jband,yourid)
            k = columns_k(column,jband,yourid) - dk_phys(column,jband,yourid,bandPL(myid))
            xPL(2*i+1,k,j) = dreal(buffR(j,column))
            xPL(2*i+2,k,j) = dimag(buffR(j,column))
          end do
        end do

        deallocate(buffR,buffS)

      end do
    end if
  end do

end subroutine

subroutine modes_to_planes_phys_lims (xPL,x,nystart,nyend,grid,myid,myiband,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! MODES TO PLANES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,dki,plband,grid,myiband
  integer column,nystart,nyend
  integer iband,jband
  integer inode,yourid
  integer msizeR,msizeS
  type(cfield) x  (sband:eband)
  real(8)      xPL(N(1,myiband)+2,N(2,myiband),limPL_FFT(grid,1,myid):limPL_FFT(grid,2,myid))
  complex(8), allocatable :: buffS(:,:),buffR(:,:)

  plband = bandPL_FFT(myid)
  yourid = myid
  do iband = sband,eband
    jband = iband
!     jminR = max(max(limPL_FFT(grid,1,  myid),jlim(1,grid,jband)+1),nystart)
!     jmaxR = min(min(limPL_FFT(grid,2,  myid),jlim(2,grid,jband)-1),nyend)
jminR = max(max(limPL_FFT(grid,1,  myid),jlim(1,grid,jband)),nystart)
jmaxR = min(min(limPL_FFT(grid,2,  myid),jlim(2,grid,jband)),nyend)
!     if (jminR==Ny(grid,0)+1  ) then
!           jminR = jminR-1
!         end if
!         if (jmaxR==Ny(grid,nband)) then
!           jmaxR = jmaxR+1
!         end if
    do j = jminR,jmaxR
      do column = 1,columns_num(jband,yourid)
        i = columns_i(column,jband,yourid)
        k = columns_k(column,jband,yourid) - dk_phys(column,jband,yourid,bandPL_FFT(myid))
        xPL(2*i+1,k,j) = dreal(x(jband)%f(j,column))
        xPL(2*i+2,k,j) = dimag(x(jband)%f(j,column))
      end do
    end do
  end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)
    if (yourid<np) then
      do iband = sband,eband
        !jband=crossband(iband,yourid)
        jband = iband
!         jminS = max(max(limPL_FFT(grid,1,yourid),jlim(1,grid,iband)+1),nystart)
!         jmaxS = min(min(limPL_FFT(grid,2,yourid),jlim(2,grid,iband)-1),nyend)
!         jminR = max(max(limPL_FFT(grid,1,  myid),jlim(1,grid,jband)+1),nystart)
!         jmaxR = min(min(limPL_FFT(grid,2,  myid),jlim(2,grid,jband)-1),nyend)
jminS = max(max(limPL_FFT(grid,1,yourid),jlim(1,grid,iband)),nystart)
jmaxS = min(min(limPL_FFT(grid,2,yourid),jlim(2,grid,iband)),nyend)
jminR = max(max(limPL_FFT(grid,1,  myid),jlim(1,grid,jband)),nystart)
jmaxR = min(min(limPL_FFT(grid,2,  myid),jlim(2,grid,jband)),nyend)
!         if (jminS==Ny(grid,0)+1  ) then
!           jminS = jminS-1
!         end if
!         if (jmaxS==Ny(grid,nband)) then
!           jmaxS = jmaxS+1
!         end if
!         if (jminR==Ny(grid,0)+1  ) then
!           jminR = jminR-1
!         end if
!         if (jmaxR==Ny(grid,nband)) then
!           jmaxR = jmaxR+1
!         end if
        allocate(buffS(jminS:jmaxS,columns_num(iband,  myid)))
        allocate(buffR(jminR:jmaxR,columns_num(jband,yourid)))
        msizeS = 2*(columns_num(iband,  myid)*(jmaxS-jminS+1))  ! 2 times because it's complex
        msizeR = 2*(columns_num(jband,yourid)*(jmaxR-jminR+1))
        msizeS = max(msizeS,0)
        msizeR = max(msizeR,0)

        do j = jminS,jmaxS
          do column = 1,columns_num(iband,myid)
            buffS(j,column) = x(iband)%f(j,column)

          end do
        end do

        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*iband+11*jband, &
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*iband+7*jband, &
&                         MPI_COMM_WORLD,status,ierr)
        do j = jminR,jmaxR
          do column = 1,columns_num(jband,yourid)
            i = columns_i(column,jband,yourid)
            k = columns_k(column,jband,yourid) - dk_phys(column,jband,yourid,bandPL_FFT(myid))
            xPL(2*i+1,k,j) = dreal(buffR(j,column))
            xPL(2*i+2,k,j) = dimag(buffR(j,column))
          end do
        end do

        deallocate(buffR,buffS)

      end do
    end if
  end do

end subroutine

subroutine modes_to_planes_phys_lims_2 (xPL,x,nystart,nyend,grid,myid,myiband,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! MODES TO PLANES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,dki,plband,grid,myiband
  integer column,nystart,nyend
  integer iband,jband
  integer inode,yourid
  integer msizeR,msizeS
  type(cfield) x  (sband:eband)
  real(8)      xPL(N(1,myiband)+2,N(2,myiband),jgal(grid,1)-1:jgal(grid,2)+1)
  complex(8), allocatable :: buffS(:,:),buffR(:,:)
  
  plband = bandPL(myid)
  yourid = myid
  do iband = sband,eband
    jband = iband
    jminR = max(max(limPL_incw(grid,1,  myid),jlim(1,grid,jband)+1),nystart)
    jmaxR = min(min(limPL_incw(grid,2,  myid),jlim(2,grid,jband)-1),nyend)
    if (jminR==Ny(grid,0)+1  ) then
          jminR = jminR-1
        end if
        if (jmaxR==Ny(grid,nband)) then
          jmaxR = jmaxR+1
        end if
        
    do j = jminR,jmaxR
      do column = 1,columns_num(jband,yourid)
        i = columns_i(column,jband,yourid)
        k = columns_k(column,jband,yourid) - dk_phys(column,jband,yourid,bandPL(myid))
        xPL(2*i+1,k,j) = dreal(x(jband)%f(j,column))
        xPL(2*i+2,k,j) = dimag(x(jband)%f(j,column))
      end do
    end do
  end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)
    if (yourid<np) then
      do iband = sband,eband
        !jband=crossband(iband,yourid)
        jband = iband
        jminS = max(max(limPL_incw(grid,1,yourid),jlim(1,grid,iband)+1),nystart)
        jmaxS = min(min(limPL_incw(grid,2,yourid),jlim(2,grid,iband)-1),nyend)
        jminR = max(max(limPL_incw(grid,1,  myid),jlim(1,grid,jband)+1),nystart)
        jmaxR = min(min(limPL_incw(grid,2,  myid),jlim(2,grid,jband)-1),nyend)
        if (jminS==Ny(grid,0)+1  ) then
          jminS = jminS-1
        end if
        if (jmaxS==Ny(grid,nband)) then
          jmaxS = jmaxS+1
        end if
        if (jminR==Ny(grid,0)+1  ) then
          jminR = jminR-1
        end if
        if (jmaxR==Ny(grid,nband)) then
          jmaxR = jmaxR+1
        end if
                
        allocate(buffS(jminS:jmaxS,columns_num(iband,  myid)))
        allocate(buffR(jminR:jmaxR,columns_num(jband,yourid)))
        msizeS = 2*(columns_num(iband,  myid)*(jmaxS-jminS+1))  ! 2 times because it's complex
        msizeR = 2*(columns_num(jband,yourid)*(jmaxR-jminR+1))
        msizeS = max(msizeS,0)
        msizeR = max(msizeR,0)

        do j = jminS,jmaxS
          do column = 1,columns_num(iband,myid)
            buffS(j,column) = x(iband)%f(j,column)
          end do
        end do

        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*iband+11*jband, &
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*iband+7*jband, &
&                         MPI_COMM_WORLD,status,ierr)

        do j = jminR,jmaxR
          do column = 1,columns_num(jband,yourid)
            i = columns_i(column,jband,yourid)
            k = columns_k(column,jband,yourid) - dk_phys(column,jband,yourid,bandPL(myid))
            xPL(2*i+1,k,j) = dreal(buffR(j,column))
            xPL(2*i+2,k,j) = dimag(buffR(j,column))
          end do
        end do

        deallocate(buffR,buffS)

      end do
    end if
  end do

end subroutine

subroutine modes_to_planes_lims_2 (xPL,x,nystart,nyend,grid,myid,myiband,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! MODES TO PLANES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,dki,plband,grid,myiband
  integer column,nystart,nyend
  integer iband,jband
  integer inode,yourid
  integer msizeR,msizeS
  type(cfield) x  (sband:eband)
  real(8)      xPL(Ngal(1,myiband)+2,Ngal(2,myiband),jgal(grid,1)-1:jgal(grid,2)+1)
  complex(8), allocatable :: buffS(:,:),buffR(:,:)

  plband = bandPL(myid)
  yourid = myid
  do iband = sband,eband
    jband = iband
    jminR = max(max(limPL_incw(grid,1,  myid),jlim(1,grid,jband)),nystart)
    jmaxR = min(min(limPL_incw(grid,2,  myid),jlim(2,grid,jband)),nyend  )
    
    do j = jminR,jmaxR
      do column = 1,columns_num(jband,yourid)
        i = columns_i(column,jband,yourid)
        k = columns_k(column,jband,yourid) - dk(column,jband,yourid,bandPL(myid))
        xPL(2*i+1,k,j) = dreal(x(jband)%f(j,column))
        xPL(2*i+2,k,j) = dimag(x(jband)%f(j,column))
      end do
    end do
  end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)
    if (yourid<np) then
      do iband = sband,eband
        !jband=crossband(iband,yourid)
        jband = iband
        jminS = max(max(limPL_incw(grid,1,yourid),jlim(1,grid,iband)),nystart)
        jmaxS = min(min(limPL_incw(grid,2,yourid),jlim(2,grid,iband)),nyend  )
        jminR = max(max(limPL_incw(grid,1,  myid),jlim(1,grid,jband)),nystart)
        jmaxR = min(min(limPL_incw(grid,2,  myid),jlim(2,grid,jband)),nyend  )
        
        allocate(buffS(jminS:jmaxS,columns_num(iband,  myid)))
        allocate(buffR(jminR:jmaxR,columns_num(jband,yourid)))
        msizeS = 2*(columns_num(iband,  myid)*(jmaxS-jminS+1))  ! 2 times because it's complex
        msizeR = 2*(columns_num(jband,yourid)*(jmaxR-jminR+1))
        msizeS = max(msizeS,0)
        msizeR = max(msizeR,0)

        do j = jminS,jmaxS
          do column = 1,columns_num(iband,myid)
            buffS(j,column) = x(iband)%f(j,column)

          end do
        end do

        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*iband+11*jband, &
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*iband+7*jband, &
&                         MPI_COMM_WORLD,status,ierr)
        do j = jminR,jmaxR
          do column = 1,columns_num(jband,yourid)
            i = columns_i(column,jband,yourid)
            k = columns_k(column,jband,yourid) - dk(column,jband,yourid,bandPL(myid))
            xPL(2*i+1,k,j) = dreal(buffR(j,column))
            xPL(2*i+2,k,j) = dimag(buffR(j,column))
          end do
        end do

        deallocate(buffR,buffS)

      end do
    end if
  end do

end subroutine

subroutine modes_to_planes_dU(xPL,x,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! MODES TO PLANES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid,grid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,dki,plband
  integer column
  integer iband,jband
  integer inode,yourid
  integer msizeR,msizeS
  type(cfield) x  (sband:eband)
  real(8)      xPL(igal,kgal,nyuIB1(myid):nyuIB2(myid))
  complex(8), allocatable :: buffS(:,:),buffR(:,:)

  plband = bandPL(myid)

  yourid = myid
  do iband = sband,eband
    !jband=crossband(iband,yourid)
    jband = iband
    jminR = max(nyuIB1(myid),jlim(1,ugrid,jband)+1)
    jmaxR = min(nyuIB2(myid),jlim(2,ugrid,jband)-1)
    do j = jminR,jmaxR
      do column = 1,columns_num(jband,yourid)
        i = columns_i(column,jband,yourid)
        k = columns_k(column,jband,yourid) - dk(column,jband,yourid,bandPL(myid))
        xPL(2*i+1,k,j) = dreal(x(jband)%f(j,column))
        xPL(2*i+2,k,j) = dimag(x(jband)%f(j,column))
      end do
    end do
  end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)
    if (yourid<np) then
      do iband = sband,eband
        !jband=crossband(iband,yourid)
        jband = iband
        jminS = max(nyuIB1(yourid),jlim(1,ugrid,iband)+1)
        jmaxS = min(nyuIB2(yourid),jlim(2,ugrid,iband)-1)
        jminR = max(nyuIB1(myid  ),jlim(1,ugrid,jband)+1)
        jmaxR = min(nyuIB2(myid  ),jlim(2,ugrid,jband)-1)
        allocate(buffS(jminS:jmaxS,columns_num(iband,  myid)))
        allocate(buffR(jminR:jmaxR,columns_num(jband,yourid)))
        msizeS = 2*(columns_num(iband,  myid)*(jmaxS-jminS+1))  ! 2 times because it's complex
        msizeR = 2*(columns_num(jband,yourid)*(jmaxR-jminR+1))
        msizeS = max(msizeS,0)
        msizeR = max(msizeR,0)

        do j = jminS,jmaxS
          do column = 1,columns_num(iband,myid)
            buffS(j,column) = x(iband)%f(j,column)
          end do
        end do
        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*iband+11*jband, &
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*iband+7*jband, &
&                         MPI_COMM_WORLD,status,ierr)
        do j = jminR,jmaxR
          do column = 1,columns_num(jband,yourid)
            i = columns_i(column,jband,yourid)
            k = columns_k(column,jband,yourid) - dk(column,jband,yourid,bandPL(myid))
            xPL(2*i+1,k,j) = dreal(buffR(j,column))
            xPL(2*i+2,k,j) = dimag(buffR(j,column))
          end do
        end do

        deallocate(buffR,buffS)

      end do
    end if
  end do

end subroutine

subroutine modes_to_planes_dV(xPL,x,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! MODES TO PLANES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid,grid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,dki,plband
  integer column
  integer iband,jband
  integer inode,yourid
  integer msizeR,msizeS
  type(cfield) x  (sband:eband)
  real(8)      xPL(igal,kgal,nyvIB1(myid):nyvIB2(myid))
  complex(8), allocatable :: buffS(:,:),buffR(:,:)

  plband = bandPL(myid)

  yourid = myid
  do iband = sband,eband
    !jband=crossband(iband,yourid)
    jband = iband
    jminR = max(nyvIB1(myid),jlim(1,vgrid,jband)+1)
    jmaxR = min(nyvIB2(myid),jlim(2,vgrid,jband)-1)
    do j = jminR,jmaxR
      do column = 1,columns_num(jband,yourid)
        i = columns_i(column,jband,yourid)
        k = columns_k(column,jband,yourid) - dk(column,jband,yourid,bandPL(myid))
        xPL(2*i+1,k,j) = dreal(x(jband)%f(j,column))
        xPL(2*i+2,k,j) = dimag(x(jband)%f(j,column))
      end do
    end do
  end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)
    if (yourid<np) then
      do iband = sband,eband
        !jband=crossband(iband,yourid)
        jband = iband
        jminS = max(nyvIB1(yourid),jlim(1,vgrid,iband)+1)
        jmaxS = min(nyvIB2(yourid),jlim(2,vgrid,iband)-1)
        jminR = max(nyvIB1(myid  ),jlim(1,vgrid,jband)+1)
        jmaxR = min(nyvIB2(myid  ),jlim(2,vgrid,jband)-1)
        allocate(buffS(jminS:jmaxS,columns_num(iband,  myid)))
        allocate(buffR(jminR:jmaxR,columns_num(jband,yourid)))
        msizeS = 2*(columns_num(iband,  myid)*(jmaxS-jminS+1))  ! 2 times because it's complex
        msizeR = 2*(columns_num(jband,yourid)*(jmaxR-jminR+1))
        msizeS = max(msizeS,0)
        msizeR = max(msizeR,0)

        do j = jminS,jmaxS
          do column = 1,columns_num(iband,myid)
            buffS(j,column) = x(iband)%f(j,column)
          end do
        end do
        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*iband+11*jband, &
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*iband+7*jband, &
&                         MPI_COMM_WORLD,status,ierr)
        do j = jminR,jmaxR
          do column = 1,columns_num(jband,yourid)
            i = columns_i(column,jband,yourid)
            k = columns_k(column,jband,yourid) - dk(column,jband,yourid,bandPL(myid))
            xPL(2*i+1,k,j) = dreal(buffR(j,column))
            xPL(2*i+2,k,j) = dimag(buffR(j,column))
          end do
        end do

        deallocate(buffR,buffS)

      end do
    end if
  end do

end subroutine

subroutine Qcriterion(myid)

  use declaration
  implicit none
  
  integer i,k,j,myid,flagst,temp
  real(8) ddyi
  real(8), allocatable:: du1dx(:,:),du1dz(:,:),du2dx(:,:),du2dz(:,:),du3dx(:,:),du3dz(:,:)
  real(8), allocatable:: du1dy(:,:),du2dy(:,:),du3dy(:,:)
  
  real(8), allocatable:: utmp(:,:)
    
  allocate(du1dx(Ngal(1,2)+2,Ngal(2,2)),du1dz(Ngal(1,2)+2,Ngal(2,2)))
  allocate(du2dx(Ngal(1,2)+2,Ngal(2,2)),du2dz(Ngal(1,2)+2,Ngal(2,2)))
  allocate(du3dx(Ngal(1,2)+2,Ngal(2,2)),du3dz(Ngal(1,2)+2,Ngal(2,2)))
  
  allocate(du1dy(Ngal(1,2)+2,Ngal(2,2)))
  allocate(du2dy(Ngal(1,2)+2,Ngal(2,2)))
  allocate(du3dy(Ngal(1,2)+2,Ngal(2,2)))
  
  allocate(utmp(Ngal(1,2)+2,Ngal(2,2)))
  
  do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
        
    utmp(:,1:Ngal(2,2)/2)           = u1PL_itp(1:Ngal(1,2)+2,1:Ngal(2,2)/2,j)
    utmp(:,Ngal(2,2)/2+1:Ngal(2,2)) = u1PL_itp(1:Ngal(1,2)+2,Ngal(2,bandPL(myid))-Ngal(2,2)/2+1:Ngal(2,bandPL(myid)),j)
    call der_x(utmp(1,1),du1dx,k1F_x,2)
    call der_z(utmp(1,1),du1dz,k1F_z,2)
    utmp(:,1:Ngal(2,2)/2)           = u2PL(1:Ngal(1,2)+2,1:Ngal(2,2)/2,j)
    utmp(:,Ngal(2,2)/2+1:Ngal(2,2)) = u2PL(1:Ngal(1,2)+2,Ngal(2,bandPL(myid))-Ngal(2,2)/2+1:Ngal(2,bandPL(myid)),j)
    call der_x(utmp(1,1),du2dx,k1F_x,2)
    call der_z(utmp(1,1),du2dz,k1F_z,2)
    utmp(:,1:Ngal(2,2)/2)           = u3PL_itp(1:Ngal(1,2)+2,1:Ngal(2,2)/2,j)
    utmp(:,Ngal(2,2)/2+1:Ngal(2,2)) = u3PL_itp(1:Ngal(1,2)+2,Ngal(2,bandPL(myid))-Ngal(2,2)/2+1:Ngal(2,bandPL(myid)),j)
    call der_x(utmp(1,1),du3dx,k1F_x,2)
    call der_z(utmp(1,1),du3dz,k1F_z,2)
       
    call four_to_phys_du(du1dx(1,1),2)
    call four_to_phys_du(du1dz(1,1),2)
    call four_to_phys_du(du2dx(1,1),2)
    call four_to_phys_du(du2dz(1,1),2)
    call four_to_phys_du(du3dx(1,1),2)
    call four_to_phys_du(du3dz(1,1),2)
    
    du1dy(:,1:Ngal(2,2)/2)           = du1dy_planes2(1:Ngal(1,2)+2,1:Ngal(2,2)/2,j)
    du1dy(:,Ngal(2,2)/2+1:Ngal(2,2)) = du1dy_planes2(1:Ngal(1,2)+2,Ngal(2,bandPL(myid))-Ngal(2,2)/2+1:Ngal(2,bandPL(myid)),j)
    du2dy(:,1:Ngal(2,2)/2)           = du2dy_planes2(1:Ngal(1,2)+2,1:Ngal(2,2)/2,j)
    du2dy(:,Ngal(2,2)/2+1:Ngal(2,2)) = du2dy_planes2(1:Ngal(1,2)+2,Ngal(2,bandPL(myid))-Ngal(2,2)/2+1:Ngal(2,bandPL(myid)),j)
    du3dy(:,1:Ngal(2,2)/2)           = du3dy_planes2(1:Ngal(1,2)+2,1:Ngal(2,2)/2,j)
    du3dy(:,Ngal(2,2)/2+1:Ngal(2,2)) = du3dy_planes2(1:Ngal(1,2)+2,Ngal(2,bandPL(myid))-Ngal(2,2)/2+1:Ngal(2,bandPL(myid)),j)
    
    call four_to_phys_du(du1dy(1,1),2)
    call four_to_phys_du(du2dy(1,1),2) 
    call four_to_phys_du(du3dy(1,1),2)
    
    do i = 1,Ngal(1,2)
      do k = 1,Ngal(2,2)
	Qcrit(i,k,j) = -0.5d0*(du1dx(i,k)**2+du2dy(i,k)**2+du3dz(i,k)**2 + &
&2d0*du1dy(i,k)*du2dx(i,k)+2d0*du1dz(i,k)*du3dx(i,k)+2d0*du2dz(i,k)*du3dy(i,k))
      enddo
    enddo
  enddo
  
  call write_Qcrit(myid)
  
  deallocate(du1dx,du1dz)
  deallocate(du2dx,du2dz)
  deallocate(du3dx,du3dz)
  
  deallocate(du1dy)
  deallocate(du2dy)
  deallocate(du3dy)
  
  deallocate(utmp)
  
endsubroutine

subroutine ops_in_planes(myid,flagst)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! OPS IN PLANES !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr

  
  integer i,k,j,myid,flagst,temp
  real(8) ddyi
  real(8), allocatable:: du1dx(:,:),du1dz(:,:),du2dx(:,:),du2dz(:,:),du3dx(:,:),du3dz(:,:)

  allocate(du1dx(igal,kgal),du1dz(igal,kgal))
  allocate(du2dx(igal,kgal),du2dz(igal,kgal))
  allocate(du3dx(igal,kgal),du3dz(igal,kgal))
  
  du1dx = 0d0
  du1dz = 0d0
  du2dx = 0d0
  du2dz = 0d0
  du3dx = 0d0
  du3dz = 0d0

  do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid) ! 10.3.2020 excw --> incw

    call four_to_phys_u(u1PL(1,1,j),u2PL_itp(1,1,j),u3PL(1,1,j),bandPL(myid))

    do k = 1,Ngal(2,bandPL(myid))
      do i = 1,Ngal(1,bandPL(myid))
        uu_cPL(i,k,j) = u1PL    (i,k,j)*u1PL    (i,k,j)
        uw_cPL(i,k,j) = u1PL    (i,k,j)*u3PL    (i,k,j)
        vv_cPL(i,k,j) = u2PL_itp(i,k,j)*u2PL_itp(i,k,j)
        ww_cPL(i,k,j) = u3PL    (i,k,j)*u3PL    (i,k,j)
      end do
    end do
    
    call phys_to_four_du(uu_cPL(1,1,j),bandPL(myid))
    call phys_to_four_du(uw_cPL(1,1,j),bandPL(myid))    
    call phys_to_four_du(vv_cPL(1,1,j),bandPL(myid))    
    call phys_to_four_du(ww_cPL(1,1,j),bandPL(myid))    
  
    call der_x(uu_cPL(1,1,j),du1dx,k1F_x,bandPL(myid))
    call der_z(uw_cPL(1,1,j),du1dz,k1F_z,bandPL(myid))
    call der_x(uw_cPL(1,1,j),du3dx,k1F_x,bandPL(myid))
    call der_z(ww_cPL(1,1,j),du3dz,k1F_z,bandPL(myid))
   
    do k = 1,Ngal(2,bandPL(myid))
      do i = 1,Ngal(1,bandPL(myid))
        Nu1PL(i,k,j) = du1dx(i,k)+du1dz(i,k)
        Nu3PL(i,k,j) = du3dx(i,k)+du3dz(i,k)
      end do
    end do

  end do
  
!! Includes wall for output (stats) - needed as walls now not in N?
!  if(myid==0)then
!    j = limPL_incw(ugrid,1,myid)
!    call four_to_phys_u(u1PL(1,1,j),u2PL_itp(1,1,j),u3PL(1,1,j),bandPL(myid))
!  elseif(myid==np-1)then
!    j = limPL_incw(ugrid,2,myid)
!    call four_to_phys_u(u1PL(1,1,j),u2PL_itp(1,1,j),u3PL(1,1,j),bandPL(myid))
!  endif
    
  

  do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid) ! 10.3.2020 excw --> incw
    
    call four_to_phys_u(u1PL_itp(1,1,j),u2PL(1,1,j),u3PL_itp(1,1,j),bandPL(myid))
    
    do k = 1,Ngal(2,bandPL(myid))
      do i = 1,Ngal(1,bandPL(myid))
        uv_fPL(i,k,j) = u1PL_itp(i,k,j)*u2PL(i,k,j)
        vw_fPL(i,k,j) = u3PL_itp(i,k,j)*u2PL(i,k,j)
      end do
    end do
    
    call phys_to_four_du(uv_fPL(1,1,j),bandPL(myid))
    call phys_to_four_du(vw_fPL(1,1,j),bandPL(myid)) 

    call der_x(uv_fPL(1,1,j),du2dx,k1F_x,bandPL(myid))
    call der_z(vw_fPL(1,1,j),du2dz,k1F_z,bandPL(myid))

    do k = 1,Ngal(2,bandPL(myid))
      do i = 1,Ngal(1,bandPL(myid))
        Nu2PL(i,k,j) = du2dx(i,k)+du2dz(i,k)   
      end do
    end do
    
  end do


!! Includes wall for output (stats) - needed as walls now not in N?
!  if(myid==0)then
!    j = limPL_incw(vgrid,1,myid)
!    call four_to_phys_u(u1PL_itp(1,1,j),u2PL(1,1,j),u3PL_itp(1,1,j),bandPL(myid))
!  elseif(myid==np-1)then
!    j = limPL_incw(vgrid,2,myid)
!    call four_to_phys_u(u1PL_itp(1,1,j),u2PL(1,1,j),u3PL_itp(1,1,j),bandPL(myid))
!  endif
  
  
  

  deallocate(du1dx,du1dz,du2dx,du2dz,du3dx,du3dz)

end subroutine

subroutine phys_to_four_du(duPL,iband)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!      FOUR TO PHYS     !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Transforms from Physical Space to Fourier Space u 
! Its only used in ops in planes at FOU3D.f90

  use declaration
  implicit none

  integer iband
  real(8) duPL(Ngal(1,iband)+2,Ngal(2,iband))

  call rft(duPL,Ngal(1,iband)+2,Ngal(2,iband),-1,buffRal_x(iband)%b)
  call cft(duPL,Ngal(1,iband)+2,2,(N(1,iband)+2)/2,-1,buffCal_z(iband)%b)
  
  duPL(:,Ngal(2,iband)/2+1)=0d0 !oddball advective term = 0

end subroutine

subroutine phys_to_four_N(duPL,iband)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!      FOUR TO PHYS     !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Transforms from Physical Space to Fourier Space u 
! Its only used in ops in planes at FOU3D.f90

  use declaration
  implicit none

  integer iband
  real(8) duPL(N(1,iband)+2,N(2,iband))

  call rft(duPL,N(1,iband)+2,N(2,iband),-1,buffR_x(iband)%b)
  call cft(duPL,N(1,iband)+2,2,(N(1,iband)+2)/2,-1,buffC_z(iband)%b)
  
  !duPL(:,N(2,iband)/2+1)=0d0

end subroutine

subroutine planes_to_modes_UVP (x,xPL,grid,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! PLANES TO MODES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prepare the vectors for the Fourier transform
! The procs broadcast the data they have of a plane, and receive the data of a pencil,
!  they also transpose the data for the Fourier transform XZY -> YXZ

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,plband,grid
  integer column
  integer iband,jband
  integer inode,yourid
  integer msizeR,msizeS
  type(cfield) x  (sband:eband)
  real(8)      xPL(igal,kgal,jgal(grid,1)-1:jgal(grid,2)+1) ! 6.3.2020
  complex(8), allocatable:: buffS(:,:),buffR(:,:)

  ! Loop for itself
  ! Transpose the cube that it already owns
  plband = bandPL(myid) ! Return the band (phys) the proc works at
  do iband = sband,eband
    jminR = max(limPL_excw(grid,1,myid),jlim(1,grid,iband)+1)  ! Select the planes to transpose 
    jmaxR = min(limPL_excw(grid,2,myid),jlim(2,grid,iband)-1)
    if (jminR==Ny(grid,plband-1)+1 .and. jmaxR>=jminR) then   ! Special cases: interfaces
      jminR = jminR-1
    end if
    if (jmaxR==Ny(grid,plband  )   .and. jmaxR>=jminR) then
      jmaxR = jmaxR+1
    end if
    do j = jminR,jmaxR
      do column = 1,columns_num(iband,myid)
        i = columns_i(column,iband,myid)
        k = columns_k(column,iband,myid) - dk(column,iband,myid,bandPL(myid))
        x(iband)%f(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j)) ! Transposition: Reordering from XZY to YC
      end do
    end do
  end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)   ! XOR. It's used to pair procs 1-to-1
    if (yourid<np) then
      do iband = sband,eband
        !jband = crossband(iband,yourid)
        jband = iband
        jminS = max(limPL_excw(grid,1,  myid),jlim(1,grid,jband)+1)  ! Select the planes to be SENT.
        jmaxS = min(limPL_excw(grid,2,  myid),jlim(2,grid,jband)-1)  ! max and min because maybe this proc needs less planes that the other proc has
        jminR = max(limPL_excw(grid,1,yourid),jlim(1,grid,iband)+1)  ! Select the planes to be RECEIVED
        jmaxR = min(limPL_excw(grid,2,yourid),jlim(2,grid,iband)-1)
        ! Adding the walls (but not the interfaces)
        if (jminS==Ny(grid,0)+1  ) then
          jminS=jminS-1
        end if
        if (jmaxS==Ny(grid,nband)) then
          jmaxS=jmaxS+1
        end if
        if (jminR==Ny(grid,0)+1  ) then
          jminR=jminR-1
        end if
        if (jmaxR==Ny(grid,nband)) then
          jmaxR=jmaxR+1
        end if
        allocate(buffS(jminS:jmaxS,columns_num(jband,yourid)))
        allocate(buffR(jminR:jmaxR,columns_num(iband,  myid)))
        msizeS = 2*(columns_num(jband,yourid)*(jmaxS-jminS+1))     ! Size of the data to be SENDER (times 2, because it is complex)
        msizeR = 2*(columns_num(iband,  myid)*(jmaxR-jminR+1))     ! Size of the data to be RECEIVED
        msizeS = max(msizeS,0)                                     ! The size has to be 0 or positive. 
        msizeR = max(msizeR,0)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do j=jminS,jmaxS
          do column = 1,columns_num(jband,yourid)
            i = columns_i(column,jband,yourid)
            k = columns_k(column,jband,yourid) - dk(column,jband,yourid,bandPL(myid))
            buffS(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j))     ! The data is transposed and stored in a buffer
          end do
        end do
        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*iband+11*jband, &   ! SEND_RECV so it can send and receive at the same time
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*iband+7*jband, &
&                         MPI_COMM_WORLD,status,ierr)
        do j=jminR,jmaxR
          do column = 1,columns_num(iband,myid)
            x(iband)%f(j,column) = buffR(j,column)                         ! Store the data received
          end do
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(buffR,buffS)
      end do
    end if
  end do

end subroutine

subroutine planes_to_modes_phys_lims (x,xPL,nystart,nyend,grid,myid,myiband,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! PLANES TO MODES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prepare the vectors for the Fourier transform
! The procs broadcast the data they have of a plane, and receive the data of a pencil,
!  they also transpose the data for the Fourier transform XZY -> YXZ

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,plband,grid,myiband
  integer column,nystart,nyend
  integer iband,jband
  integer inode,yourid
  integer msizeR,msizeS
  type(cfield) x  (sband:eband)
  real(8)      xPL(N(1,myiband)+2,N(2,myiband),limPL_FFT(grid,1,myid):limPL_FFT(grid,2,myid))
  complex(8), allocatable:: buffS(:,:),buffR(:,:)

  ! Loop for itself
  ! Transpose the cube that it already owns
  plband = bandPL_FFT(myid) ! Return the band (phys) the proc works at
  do iband = sband,eband
    jminR = max(max(limPL_FFT(grid,1,myid),jlim(1,grid,iband)),nystart)  ! Select the planes to transpose 
    jmaxR = min(min(limPL_FFT(grid,2,myid),jlim(2,grid,iband)),nyend  )
    
    do j = jminR,jmaxR
      do column = 1,columns_num(iband,myid)
        i = columns_i(column,iband,myid)
        k = columns_k(column,iband,myid) - dk_phys(column,iband,myid,bandPL_FFT(myid))
        x(iband)%f(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j)) ! Transposition: Reordering from XZY to YC
      end do
    end do
  end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)   ! XOR. It's used to pair procs 1-to-1
    if (yourid<np) then
      do iband = sband,eband
        !jband = crossband(iband,yourid)
        jband = iband
        jminS = max(max(limPL_FFT(grid,1,  myid),jlim(1,grid,jband)),nystart)  ! Select the planes to be SENT.
        jmaxS = min(min(limPL_FFT(grid,2,  myid),jlim(2,grid,jband)),nyend  )  ! max and min because maybe this proc needs less planes that the other proc has
        jminR = max(max(limPL_FFT(grid,1,yourid),jlim(1,grid,iband)),nystart)  ! Select the planes to be RECEIVED
        jmaxR = min(min(limPL_FFT(grid,2,yourid),jlim(2,grid,iband)),nyend  )

        allocate(buffS(jminS:jmaxS,columns_num(jband,yourid)))
        allocate(buffR(jminR:jmaxR,columns_num(iband,  myid)))
        msizeS = 2*(columns_num(jband,yourid)*(jmaxS-jminS+1))     ! Size of the data to be SENDER (times 2, because it is complex)
        msizeR = 2*(columns_num(iband,  myid)*(jmaxR-jminR+1))     ! Size of the data to be RECEIVED
        msizeS = max(msizeS,0)                                     ! The size has to be 0 or positive. 
        msizeR = max(msizeR,0)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do j=jminS,jmaxS
          do column = 1,columns_num(jband,yourid)
            i = columns_i(column,jband,yourid)
            k = columns_k(column,jband,yourid) - dk_phys(column,jband,yourid,bandPL_FFT(myid))
            buffS(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j))     ! The data is transposed and stored in a buffer
          end do
        end do
        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*iband+11*jband, &   ! SEND_RECV so it can send and receive at the same time
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*iband+7*jband, &
&                         MPI_COMM_WORLD,status,ierr)
        do j=jminR,jmaxR
          do column = 1,columns_num(iband,myid)
            x(iband)%f(j,column) = buffR(j,column)                         ! Store the data received
          end do
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(buffR,buffS)
      end do
    end if
  end do
end subroutine

subroutine planes_to_modes_phys_lims_2 (x,xPL,nystart,nyend,grid,myid,myiband,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! PLANES TO MODES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prepare the vectors for the Fourier transform
! The procs broadcast the data they have of a plane, and receive the data of a pencil,
!  they also transpose the data for the Fourier transform XZY -> YXZ

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,plband,grid,myiband
  integer column,nystart,nyend
  integer iband,jband
  integer inode,yourid
  integer msizeR,msizeS
  type(cfield) x  (sband:eband)
  real(8)      xPL(N(1,myiband)+2,N(2,myiband),jgal(grid,1)-1:jgal(grid,2)+1)
  complex(8), allocatable:: buffS(:,:),buffR(:,:)

  ! Loop for itself
  ! Transpose the cube that it already owns
  plband = bandPL(myid) ! Return the band (phys) the proc works at
  do iband = sband,eband
    jminR = max(max(limPL_incw(grid,1,myid),jlim(1,grid,iband)),nystart)  ! Select the planes to transpose 
    jmaxR = min(min(limPL_incw(grid,2,myid),jlim(2,grid,iband)),nyend  )
    
    do j = jminR,jmaxR
      do column = 1,columns_num(iband,myid)
        i = columns_i(column,iband,myid)
        k = columns_k(column,iband,myid) - dk_phys(column,iband,myid,bandPL(myid))
        x(iband)%f(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j)) ! Transposition: Reordering from XZY to YC
      end do
    end do
  end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)   ! XOR. It's used to pair procs 1-to-1
    if (yourid<np) then
      do iband = sband,eband
        !jband = crossband(iband,yourid)
        jband = iband
        jminS = max(max(limPL_incw(grid,1,  myid),jlim(1,grid,jband)),nystart)  ! Select the planes to be SENT.
        jmaxS = min(min(limPL_incw(grid,2,  myid),jlim(2,grid,jband)),nyend  )  ! max and min because maybe this proc needs less planes that the other proc has
        jminR = max(max(limPL_incw(grid,1,yourid),jlim(1,grid,iband)),nystart)  ! Select the planes to be RECEIVED
        jmaxR = min(min(limPL_incw(grid,2,yourid),jlim(2,grid,iband)),nyend  )

        allocate(buffS(jminS:jmaxS,columns_num(jband,yourid)))
        allocate(buffR(jminR:jmaxR,columns_num(iband,  myid)))
        msizeS = 2*(columns_num(jband,yourid)*(jmaxS-jminS+1))     ! Size of the data to be SENDER (times 2, because it is complex)
        msizeR = 2*(columns_num(iband,  myid)*(jmaxR-jminR+1))     ! Size of the data to be RECEIVED
        msizeS = max(msizeS,0)                                     ! The size has to be 0 or positive. 
        msizeR = max(msizeR,0)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do j=jminS,jmaxS
          do column = 1,columns_num(jband,yourid)
            i = columns_i(column,jband,yourid)
            k = columns_k(column,jband,yourid) - dk_phys(column,jband,yourid,bandPL(myid))
            buffS(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j))     ! The data is transposed and stored in a buffer
          end do
        end do
        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*iband+11*jband, &   ! SEND_RECV so it can send and receive at the same time
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*iband+7*jband, &
&                         MPI_COMM_WORLD,status,ierr)
        do j=jminR,jmaxR
          do column = 1,columns_num(iband,myid)
            x(iband)%f(j,column) = buffR(j,column)                         ! Store the data received
          end do
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(buffR,buffS)
      end do
    end if
  end do
end subroutine

subroutine planes_to_modes_dU (x,xPL,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! PLANES TO MODES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prepare the vectors for the Fourier transform
! The procs broadcast the data they have of a plane, and receive the data of a pencil,
!  they also transpose the data for the Fourier transform XZY -> YXZ


  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,plband,grid
  integer iband,jband
  integer column
  integer inode,yourid
  integer msizeR,msizeS
  type(cfield) x  (sband:eband)
  real(8)      xPL(igal,kgal,nyuIB1(myid):nyuIB2(myid))
  complex(8), allocatable:: buffS(:,:),buffR(:,:)

! DELETE  plband=bandPL(myid)  

  ! Loop for itself
  ! Transpose the cube that it already owns
  plband = bandPL(myid) ! Return the band (phys) the proc works at
  do iband = sband,eband
    jminR = max(nyuIB1(myid),jlim(1,ugrid,iband)+1)  ! Select the planes to transpose 
    jmaxR = min(nyuIB2(myid),jlim(2,ugrid,iband)-1)
    do j = jminR,jmaxR
      do column = 1,columns_num(iband,myid)
        i = columns_i(column,iband,myid)
        k = columns_k(column,iband,myid) - dk(column,iband,myid,bandPL(myid))
        x(iband)%f(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j)) ! Transposition: Reordering from XZY to YC
      end do
    end do
  end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)   ! XOR. It's used to pair procs 1-to-1
    if (yourid<np) then
      do iband = sband,eband
        jband = iband
        jminS = max(nyuIB1(  myid),jlim(1,ugrid,jband)+1)  ! Select the planes to be SENT.
        jmaxS = min(nyuIB2(  myid),jlim(2,ugrid,jband)-1) 
        jminR = max(nyuIB1(yourid),jlim(1,ugrid,iband)+1)  ! Select the planes to be RECEIVED
        jmaxR = min(nyuIB2(yourid),jlim(2,ugrid,iband)-1)
        ! Special cases: interfaces
        allocate(buffS(jminS:jmaxS,columns_num(jband,yourid)))
        allocate(buffR(jminR:jmaxR,columns_num(iband,  myid)))
        msizeS = 2*(columns_num(jband,yourid)*(jmaxS-jminS+1))     ! Size of the data to be SENDER (times 2, because it is complex)
        msizeR = 2*(columns_num(iband,  myid)*(jmaxR-jminR+1))     ! Size of the data to be RECEIVED
        msizeS = max(msizeS,0)                                     ! The size has to be 0 or positive. 
        msizeR = max(msizeR,0)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do j=jminS,jmaxS
          do column = 1,columns_num(jband,yourid)
            i = columns_i(column,jband,yourid)
            k = columns_k(column,jband,yourid) - dk(column,jband,yourid,bandPL(myid))
            buffS(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j))     ! The data is transposed and stored in a buffer
          end do
        end do
        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*iband+11*jband, &   ! SEND_RECV so it can send and receive at the same time
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*iband+7*jband, &
&                         MPI_COMM_WORLD,status,ierr)
        do j=jminR,jmaxR
          do column = 1,columns_num(iband,myid)
            x(iband)%f(j,column) = buffR(j,column)                         ! Store the data received
          end do
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(buffR,buffS)
      end do
    end if
  end do

end subroutine

subroutine planes_to_modes_dV (x,xPL,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! PLANES TO MODES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prepare the vectors for the Fourier transform
! The procs broadcast the data they have of a plane, and receive the data of a pencil,
!  they also transpose the data for the Fourier transform XZY -> YXZ


  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,plband,grid
  integer iband,jband
  integer column
  integer inode,yourid
  integer msizeR,msizeS
  type(cfield) x  (sband:eband)
  real(8)      xPL(igal,kgal,nyvIB1(myid):nyvIB2(myid))
  complex(8), allocatable:: buffS(:,:),buffR(:,:)

! DELETE  plband=bandPL(myid)  

  ! Loop for itself
  ! Transpose the cube that it already owns
  plband = bandPL(myid) ! Return the band (phys) the proc works at
  do iband = sband,eband
    jminR = max(nyvIB1(myid),jlim(1,vgrid,iband)+1)  ! Select the planes to transpose 
    jmaxR = min(nyvIB2(myid),jlim(2,vgrid,iband)-1)
    do j = jminR,jmaxR
      do column = 1,columns_num(iband,myid)
        i = columns_i(column,iband,myid)
        k = columns_k(column,iband,myid) - dk(column,iband,myid,bandPL(myid))
        x(iband)%f(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j)) ! Transposition: Reordering from XZY to YC
      end do
    end do
  end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)   ! XOR. It's used to pair procs 1-to-1
    if (yourid<np) then
      !do iband=bandit(1),bandit(2),bandit(3)   ! Iterate first over the short pencils and then the long ones.
      do iband = sband,eband
        !jband = crossband(iband,yourid)
        jband = iband
        jminS = max(nyvIB1(  myid),jlim(1,vgrid,jband)+1)  ! Select the planes to be SENT.
        jmaxS = min(nyvIB2(  myid),jlim(2,vgrid,jband)-1) 
        jminR = max(nyvIB1(yourid),jlim(1,vgrid,iband)+1)  ! Select the planes to be RECEIVED
        jmaxR = min(nyvIB2(yourid),jlim(2,vgrid,iband)-1)
        ! Special cases: interfaces
        allocate(buffS(jminS:jmaxS,columns_num(jband,yourid)))
        allocate(buffR(jminR:jmaxR,columns_num(iband,  myid)))
        msizeS = 2*(columns_num(jband,yourid)*(jmaxS-jminS+1))     ! Size of the data to be SENDER (times 2, because it is complex)
        msizeR = 2*(columns_num(iband,  myid)*(jmaxR-jminR+1))     ! Size of the data to be RECEIVED
        msizeS = max(msizeS,0)                                     ! The size has to be 0 or positive. 
        msizeR = max(msizeR,0)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do j=jminS,jmaxS
          do column = 1,columns_num(jband,yourid)
            i = columns_i(column,jband,yourid)
            k = columns_k(column,jband,yourid) - dk(column,jband,yourid,bandPL(myid))
            buffS(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j))     ! The data is transposed and stored in a buffer
          end do
        end do
        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*iband+11*jband, &   ! SEND_RECV so it can send and receive at the same time
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*iband+7*jband, &
&                         MPI_COMM_WORLD,status,ierr)
        do j=jminR,jmaxR
          do column = 1,columns_num(iband,myid)
            x(iband)%f(j,column) = buffR(j,column)                         ! Store the data received
          end do
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(buffR,buffS)
      end do
    end if
  end do

end subroutine

subroutine planes_to_modes_NUVP(x,xPL,grid,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! PLANES TO MODES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Prepare the vectors for the Fourier transform
! The procs broadcast the data they have of a plane, and receive the data of a pencil,
!  they also transpose the data for the Fourier transform XZY -> YXZ


  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer i,k,j,jminS,jmaxS,jminR,jmaxR,plband,grid
  integer iband,jband
  integer column
  integer inode,yourid
  integer msizeR,msizeS
  type(cfield) x  (sband:eband)
  real(8)      xPL(igal,kgal,jgal(grid,1)-1:jgal(grid,2)+1)
  complex(8), allocatable:: buffS(:,:),buffR(:,:)

! DELETE  plband=bandPL(myid)  

  ! Loop for itself
  ! Transpose the cube that it already owns
  plband = bandPL(myid) ! Return the band (phys) the proc works at
  do iband = sband,eband
    jminR = max(limPL_excw(grid,1,myid),jlim(1,grid,iband)+1)  ! Select the planes to transpose 
    jmaxR = min(limPL_excw(grid,2,myid),jlim(2,grid,iband)-1)
    do j = jminR,jmaxR
      do column = 1,columns_num(iband,myid)
        i = columns_i(column,iband,myid)
        k = columns_k(column,iband,myid) - dk(column,iband,myid,bandPL(myid))
        x(iband)%f(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j)) ! Transposition: Reordering from XZY to YC
      end do
    end do
  end do

  do inode = 1,pnodes-1
    yourid = ieor(myid,inode)   ! XOR. It's used to pair procs 1-to-1
    if (yourid<np) then
      !do iband=bandit(1),bandit(2),bandit(3)   ! Iterate first over the short pencils and then the long ones.
      do iband = sband,eband
        !jband = crossband(iband,yourid)
        jband = iband
        jminS = max(limPL_excw(grid,1,  myid),jlim(1,grid,jband)+1)  ! Select the planes to be SENT.
        jmaxS = min(limPL_excw(grid,2,  myid),jlim(2,grid,jband)-1)  ! max and min because maybe this proc needs less planes that the other proc has
        jminR = max(limPL_excw(grid,1,yourid),jlim(1,grid,iband)+1)  ! Select the planes to be RECEIVED
        jmaxR = min(limPL_excw(grid,2,yourid),jlim(2,grid,iband)-1)
        ! Special cases: interfaces
        allocate(buffS(jminS:jmaxS,columns_num(jband,yourid)))
        allocate(buffR(jminR:jmaxR,columns_num(iband,  myid)))
        msizeS = 2*(columns_num(jband,yourid)*(jmaxS-jminS+1))     ! Size of the data to be SENDER (times 2, because it is complex)
        msizeR = 2*(columns_num(iband,  myid)*(jmaxR-jminR+1))     ! Size of the data to be RECEIVED
        msizeS = max(msizeS,0)                                     ! The size has to be 0 or positive. 
        msizeR = max(msizeR,0)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do j=jminS,jmaxS
          do column = 1,columns_num(jband,yourid)
            i = columns_i(column,jband,yourid)
            k = columns_k(column,jband,yourid) - dk(column,jband,yourid,bandPL(myid))
            buffS(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j))     ! The data is transposed and stored in a buffer
          end do
        end do
        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid+7*iband+11*jband, &   ! SEND_RECV so it can send and receive at the same time
&                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid+11*iband+7*jband, &
&                         MPI_COMM_WORLD,status,ierr)
        do j=jminR,jmaxR
          do column = 1,columns_num(iband,myid)
            x(iband)%f(j,column) = buffR(j,column)                         ! Store the data received
          end do
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(buffR,buffS)
      end do
    end if
  end do

end subroutine

subroutine record_out(u1,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   RECORD OUT   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  type(cfield) u1(sband:eband)
  integer nx,nz
  integer j,jmax,iproc
  real(8), allocatable:: buffSR(:,:)
  integer, allocatable:: dummint(:)
  
  real(8) Uslip  

  if (myid/=0) then
    nx = N(1,bandPL(myid))+2
    nz = N(2,bandPL(myid))
    allocate(buffSR(nx,nz))
    do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
      call u_to_buff(buffSR,u1PL(1,1,j),nx,nz,igal,kgal)
      call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,123*myid,MPI_COMM_WORLD,ierr)
    end do
    do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      call u_to_buff(buffSR,u2PL(1,1,j),nx,nz,igal,kgal)
      call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,124*myid,MPI_COMM_WORLD,ierr)
    end do
    do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
      call u_to_buff(buffSR,u3PL(1,1,j),nx,nz,igal,kgal)
      call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,125*myid,MPI_COMM_WORLD,ierr)
    end do
    do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
      call u_to_buff(buffSR,ppPL(1,1,j),nx,nz,igal,kgal)
      call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,126*myid,MPI_COMM_WORLD,ierr)
    end do
!    do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
!      call u_to_buff(buffSR,div_outPL(1,1,j),nx,nz,igal,kgal)
!      call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,127*myid,MPI_COMM_WORLD,ierr)
!    end do
!    do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
!      call u_to_buff(buffSR,u_outPL(1,1,j),nx,nz,igal,kgal)
!      call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,128*myid,MPI_COMM_WORLD,ierr)
!    end do
!    do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
!      call u_to_buff(buffSR,v_outPL(1,1,j),nx,nz,igal,kgal)
!      call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,129*myid,MPI_COMM_WORLD,ierr)
!    end do
!    do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
!      call u_to_buff(buffSR,w_outPL(1,1,j),nx,nz,igal,kgal)
!      call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,130*myid,MPI_COMM_WORLD,ierr)
!    end do
!    do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
!      call u_to_buff(buffSR,div_cPL(1,1,j),nx,nz,igal,kgal)
!      call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,131*myid,MPI_COMM_WORLD,ierr)
!    end do
!     nx = Ngal(1,2)+2
!     nz = Ngal(2,2)
!     do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
!       call MPI_SEND(Qcrit(:,:,j),nx*nz,MPI_REAL8,0,127*myid,MPI_COMM_WORLD,ierr)
!     end do
    deallocate(buffSR)
  else
    write(ext4,'(i5.5)') int(10d0*t)!int(t)!
    allocate(dummint(88))
    dummint = 0
    !!!!!!!!!!!!!    u1    !!!!!!!!!!!!!
    fnameima = trim(dirout)//'/u1_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
    open(10,file=fnameima,form='unformatted')
    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint 
    write(10) N
    write(10) yu,dthetavi,dthdyu
    nx = N(1,bandPL(myid))+2
    nz = N(2,bandPL(myid))
    allocate(buffSR(nx,nz))
    do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
       call u_to_buff(buffSR,u1PL(1,1,j),nx,nz,igal,kgal)
       write(10) j,1,nx,nz,yu(j),buffSR
    end do
    deallocate(buffSR)
    do iproc = 1,np-1
      nx = N(1,bandPL(iproc))+2
      nz = N(2,bandPL(iproc))
      allocate(buffSR(nx,nz))
      do j = limPL_incw(ugrid,1,iproc),limPL_incw(ugrid,2,iproc)
        call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,123*iproc,MPI_COMM_WORLD,status,ierr)
        write(10) j,1,nx,nz,yu(j),buffSR
      end do
      deallocate(buffSR)
    end do
    close(10)
    !!!!!!!!!!!!!    u2    !!!!!!!!!!!!!
    fnameima=trim(dirout)//'/u2_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
    open(10,file=fnameima,form='unformatted')
    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
    write(10) N
    write(10) yv,dthetavi,dthdyv
    nx = N(1,bandPL(myid))+2
    nz = N(2,bandPL(myid))
    allocate(buffSR(nx,nz))
    do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      call u_to_buff(buffSR,u2PL(1,1,j),nx,nz,igal,kgal)
      write(10) j,2,nx,nz,yv(j),buffSR
    end do
    deallocate(buffSR)
    do iproc = 1,np-1
      nx = N(1,bandPL(iproc))+2
      nz = N(2,bandPL(iproc))
      allocate(buffSR(nx,nz))
      do j = limPL_incw(vgrid,1,iproc),limPL_incw(vgrid,2,iproc)
        call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,124*iproc,MPI_COMM_WORLD,status,ierr)
        write(10) j,2,nx,nz,yv(j),buffSR
      end do
      deallocate(buffSR)
    end do
    close(10)
    !!!!!!!!!!!!!    u3    !!!!!!!!!!!!!
    fnameima = trim(dirout)//'/u3_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
    open(10,file=fnameima,form='unformatted')
    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
    write(10) N
    write(10) yu,dthetavi,dthdyu
    nx = N(1,bandPL(myid))+2
    nz = N(2,bandPL(myid))
    allocate(buffSR(nx,nz))
    do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
      call u_to_buff(buffSR,u3PL(1,1,j),nx,nz,igal,kgal)
      write(10) j,3,nx,nz,yu(j),buffSR
    end do
    deallocate(buffSR)
    do iproc = 1,np-1
      nx = N(1,bandPL(iproc))+2
      nz = N(2,bandPL(iproc))
      allocate(buffSR(nx,nz))
      do j = limPL_incw(ugrid,1,iproc),limPL_incw(ugrid,2,iproc)
        call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,125*iproc,MPI_COMM_WORLD,status,ierr)
        write(10) j,3,nx,nz,yu(j),buffSR
      end do
      deallocate(buffSR)
    end do
    close(10)
    !!!!!!!!!!!!!    p     !!!!!!!!!!!!!
    fnameima = trim(dirout)//'/p_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
    open(10,file=fnameima,form='unformatted')
    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
    write(10) N
    write(10) yu,dthetavi,dthdyu
    nx = N(1,bandPL(myid))+2
    nz = N(2,bandPL(myid))
    allocate(buffSR(nx,nz))
    do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
      call u_to_buff(buffSR,ppPL(1,1,j),nx,nz,igal,kgal)
      write(10) j,4,nx,nz,yu(j),buffSR
    end do
    deallocate(buffSR)
    do iproc = 1,np-1
      nx = N(1,bandPL(iproc))+2
      nz = N(2,bandPL(iproc))
      allocate(buffSR(nx,nz))
      do j = limPL_incw(pgrid,1,iproc),limPL_incw(pgrid,2,iproc)
        call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,126*iproc,MPI_COMM_WORLD,status,ierr)
        write(10) j,4,nx,nz,yu(j),buffSR
      end do
      deallocate(buffSR)
    end do
    close(10)


!    !!!!!!!!!!!!!    div     !!!!!!!!!!!!!
!    fnameima = trim(dirout)//'/div_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
!    open(10,file=fnameima,form='unformatted')
!    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
!    write(10) N
!    write(10) yu,dthetavi,dthdyu
!    nx = N(1,bandPL(myid))+2
!    nz = N(2,bandPL(myid))
!    allocate(buffSR(nx,nz))
!    do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
!      call u_to_buff(buffSR,div_outPL(1,1,j),nx,nz,igal,kgal)
!      write(10) j,4,nx,nz,yu(j),buffSR
!    end do
!    deallocate(buffSR)
!    do iproc = 1,np-1
!      nx = N(1,bandPL(iproc))+2
!      nz = N(2,bandPL(iproc))
!      allocate(buffSR(nx,nz))
!      do j = limPL_incw(pgrid,1,iproc),limPL_incw(pgrid,2,iproc)
!        call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,127*iproc,MPI_COMM_WORLD,status,ierr)
!        write(10) j,4,nx,nz,yu(j),buffSR
!      end do
!      deallocate(buffSR)
!    end do
!    close(10)




!    !!!!!!!!!!!!!    u1    !!!!!!!!!!!!!
!    fnameima = trim(dirout)//'/u1_precorr_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
!    open(10,file=fnameima,form='unformatted')
!    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint 
!    write(10) N
!    write(10) yu,dthetavi,dthdyu
!    nx = N(1,bandPL(myid))+2
!    nz = N(2,bandPL(myid))
!    allocate(buffSR(nx,nz))
!    do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
!       call u_to_buff(buffSR,u_outPL(1,1,j),nx,nz,igal,kgal)
!       write(10) j,1,nx,nz,yu(j),buffSR
!    end do
!    deallocate(buffSR)
!    do iproc = 1,np-1
!      nx = N(1,bandPL(iproc))+2
!      nz = N(2,bandPL(iproc))
!      allocate(buffSR(nx,nz))
!      do j = limPL_incw(ugrid,1,iproc),limPL_incw(ugrid,2,iproc)
!        call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,128*iproc,MPI_COMM_WORLD,status,ierr)
!        write(10) j,1,nx,nz,yu(j),buffSR
!      end do
!      deallocate(buffSR)
!    end do
!    close(10)
!    !!!!!!!!!!!!!    u2    !!!!!!!!!!!!!
!    fnameima=trim(dirout)//'/u2_precorr_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
!    open(10,file=fnameima,form='unformatted')
!    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
!    write(10) N
!    write(10) yv,dthetavi,dthdyv
!    nx = N(1,bandPL(myid))+2
!    nz = N(2,bandPL(myid))
!    allocate(buffSR(nx,nz))
!    do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
!      call u_to_buff(buffSR,v_outPL(1,1,j),nx,nz,igal,kgal)
!      write(10) j,2,nx,nz,yv(j),buffSR
!    end do
!    deallocate(buffSR)
!    do iproc = 1,np-1
!      nx = N(1,bandPL(iproc))+2
!      nz = N(2,bandPL(iproc))
!      allocate(buffSR(nx,nz))
!      do j = limPL_incw(vgrid,1,iproc),limPL_incw(vgrid,2,iproc)
!        call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,129*iproc,MPI_COMM_WORLD,status,ierr)
!        write(10) j,2,nx,nz,yv(j),buffSR
!      end do
!      deallocate(buffSR)
!    end do
!!    close(10)
!    !!!!!!!!!!!!!    u3    !!!!!!!!!!!!!
!    fnameima = trim(dirout)//'/u3_precorr_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
!    open(10,file=fnameima,form='unformatted')
!    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
!    write(10) N
!    write(10) yu,dthetavi,dthdyu
!    nx = N(1,bandPL(myid))+2
!    nz = N(2,bandPL(myid))
!    allocate(buffSR(nx,nz))
!    do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
!      call u_to_buff(buffSR,w_outPL(1,1,j),nx,nz,igal,kgal)
!      write(10) j,3,nx,nz,yu(j),buffSR
!    end do
!    deallocate(buffSR)
!    do iproc = 1,np-1
!      nx = N(1,bandPL(iproc))+2
!      nz = N(2,bandPL(iproc))
!      allocate(buffSR(nx,nz))
!      do j = limPL_incw(ugrid,1,iproc),limPL_incw(ugrid,2,iproc)
!        call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,130*iproc,MPI_COMM_WORLD,status,ierr)
!        write(10) j,3,nx,nz,yu(j),buffSR
!      end do
!      deallocate(buffSR)
!    end do
!    close(10)
!
!    !!!!!!!!!!!!!    div     !!!!!!!!!!!!!
!    fnameima = trim(dirout)//'/div_c_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
!    open(10,file=fnameima,form='unformatted')
!    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
!    write(10) N
!    write(10) yu,dthetavi,dthdyu
!    nx = N(1,bandPL(myid))+2
!    nz = N(2,bandPL(myid))
!    allocate(buffSR(nx,nz))
!  
!    do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
!      call u_to_buff(buffSR,div_cPL(1,1,j),nx,nz,igal,kgal)
!      write(10) j,4,nx,nz,yu(j),buffSR
!    end do
!    deallocate(buffSR)
!    do iproc = 1,np-1
!      nx = N(1,bandPL(iproc))+2
!      nz = N(2,bandPL(iproc))
!      allocate(buffSR(nx,nz))
!      do j = limPL_incw(pgrid,1,iproc),limPL_incw(pgrid,2,iproc)
!        call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,131*iproc,MPI_COMM_WORLD,status,ierr)
!        write(10) j,4,nx,nz,yu(j),buffSR
!      end do
!      deallocate(buffSR)
!    end do
!    close(10)

!     !!!!!!!!!!!!!    Qcrit    !!!!!!!!!!!!!
!     fnameima='output/Qcrit_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
!     open(10,file=fnameima,form='unformatted')
!     write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
!     write(10) N
!     write(10) yv,dthetavi,dthdyv
!     nx = Ngal(1,2)+2
!     nz = Ngal(2,2)
!     !allocate(buffSR(nx,nz))
!     do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
!       !call u_to_buff(buffSR,u2PL(1,1,j),nx,nz,igal,kgal)
!       write(10) j,2,nx,nz,yv(j),Qcrit(:,:,j)
!     end do
!     !deallocate(buffSR)
!     do iproc = 1,np-1
!       nx = Ngal(1,2)+2
!       nz = Ngal(2,2)
!       allocate(buffSR(nx,nz))
!       do j = limPL_incw(vgrid,1,iproc),limPL_incw(vgrid,2,iproc)
!         call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,127*iproc,MPI_COMM_WORLD,status,ierr)
!         write(10) j,2,nx,nz,yv(j),buffSR
!       end do
!       deallocate(buffSR)
!     end do
!     close(10)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate(dummint)
  end if

  if (myid==0) then

    Uslip = ((u1(2)%f(1,1))*(yv(0)-yu(0))+(u1(2)%f(0,1))*(yu(1)-yv(0)))/(yu(1)-yu(0))

    call flowrateIm(Qx,u1(midband)%f(N(4,0),1))
    call maxvel(u1(midband)%f(N(4,0),1))
    write(*,*) ''
    write(*,*) 'iter',iter
    write(*,*) 't   ',t
    write(*,*) 'dtv ',dtv
    write(*,*) 'dtc ',dtc
    write(*,*) 'dt  ',dt
    write(*,*) 'err ',err
    write(*,*) 'Qx  ',Qx
    if (flag_ctpress==0) then
      write(*,*) 'QxT ',QxT
      write(*,*) 'mpgx',mpgx
      write(*,*) 'dpgx',dgx
    else
      write(*,*) 'mpgx',mpgx
    end if
    write(*,*) 'Umax',Umax
    write(*,*) 'Uslp',Uslip
!    write(*,*) 'utau',utau
! Save to history file
    if (flag_ctpress==0) then
      write(30) flag_ctpress,iter,t,dtv,dtc,dt,err,Qx,QxT,mpgx,dgx,Umax,Uslip
    else
      write(30) flag_ctpress,iter,t,dtv,dtc,dt,err,Qx,    mpgx,    Umax,Uslip
    end if
    flush(30)
  end if

end subroutine


subroutine write_Qcrit(myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   RECORD OUT   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  integer nx,nz
  integer j,jmax,iproc
  real(8), allocatable:: buffSR(:,:)
  integer, allocatable:: dummint(:)


  if (myid/=0) then
    nx = Ngal(1,2)+2
    nz = Ngal(2,2)
    do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      call MPI_SEND(Qcrit(:,:,j),nx*nz,MPI_REAL8,0,127*myid,MPI_COMM_WORLD,ierr)
    end do
  else
    write(ext4,'(i5.5)') int(10d0*t)!int(t)!
    allocate(dummint(88))
    dummint = 0
    !!!!!!!!!!!!!    Qcrit    !!!!!!!!!!!!!
    fnameima=trim(dirout)//'/Qcrit_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
    open(10,file=fnameima,form='unformatted')
    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
    write(10) N
    write(10) yv,dthetavi,dthdyv
    nx = Ngal(1,2)+2
    nz = Ngal(2,2)
    !allocate(buffSR(nx,nz))
    do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      !call u_to_buff(buffSR,u2PL(1,1,j),nx,nz,igal,kgal)
      write(10) j,2,nx,nz,yv(j),Qcrit(:,:,j)
    end do
    !deallocate(buffSR)
    do iproc = 1,np-1
      nx = Ngal(1,2)+2
      nz = Ngal(2,2)
      allocate(buffSR(nx,nz))
      do j = limPL_incw(vgrid,1,iproc),limPL_incw(vgrid,2,iproc)
        call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,127*iproc,MPI_COMM_WORLD,status,ierr)
        write(10) j,2,nx,nz,yv(j),buffSR
      end do
      deallocate(buffSR)
    end do
    close(10)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate(dummint)
  end if

end subroutine

subroutine u_to_buff(buffSR,u,nx,nz,igal,kgal)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    U to BUFF   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Rearrange z-modes before writing and after reading

  implicit none

  integer nx,nz,igal,kgal
  real(8) u(igal,kgal)
  real(8) buffSR(nx,nz)
  integer i,k,dkk

   do k = 1,nz/2
     do i = 1,nx
       buffSR(i,k) = u(i,k    )
     end do
   end do
   do i = 1,nx
     !buffSR(i,nz/2+1) = 0d0
   end do
   do k = nz/2+1,nz
     dkk = kgal-nz
     do i = 1,nx
       buffSR(i,k) = u(i,k+dkk)
     end do
   end do

end subroutine

subroutine buff_to_u(u,buffSR,nx,nz,igal,kgal)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    BUFF to U   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Rearrange z-modes before writing and after reading

  implicit none
  integer nx,nz,igal,kgal
  real(8) u(igal,kgal)
  real(8) buffSR(nx,nz)
  integer i,k,dkk

  do k = 1,min(nz/2,kgal/2)
    do i = 1,min(nx,igal)
      u(i,k    ) = buffSR(i,k)
    end do
  end do
  dkk = kgal-nz
  do k = nz-min(nz/2,kgal/2)+1,nz
    do i = 1,min(nx,igal)
      u(i,k+dkk) = buffSR(i,k)
    end do
  end do

end subroutine

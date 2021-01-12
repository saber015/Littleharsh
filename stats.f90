subroutine stats(myid,status,ierr)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!      STATS     !!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use declaration
   implicit none

   include 'mpif.h'             ! MPI variables
   integer status(MPI_STATUS_SIZE),ierr,myid

   integer i,k,j

   istat = istat+1


   !!!!!!!!    U stats:    !!!!!!!!
   do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
      do k = 1,kgal
         do i = 1,igal-2
            Um (j) = Um (j)+u1PL(i,k,j)
            U2m(j) = U2m(j)+u1PL(i,k,j)**2
         end do
      end do
   end do
   !!!!!!!!    V stats:    !!!!!!!!
   do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      do k = 1,kgal
         do i = 1,igal-2
            Vm (j) = Vm (j)+u2PL(i,k,j)
            V2m(j) = V2m(j)+u2PL(i,k,j)**2
         end do
      end do
   end do
   !!!!!!!!    W stats:    !!!!!!!!
   do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
      do k = 1,kgal
         do i = 1,igal-2
            Wm (j) = Wm (j)+u3PL(i,k,j)
            W2m(j) = W2m(j)+u3PL(i,k,j)**2
         end do
      end do
   end do
   !!!!!!!!    P stats:    !!!!!!!!
   do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
      do k = 1,kgal
         do i = 1,igal-2
            Pm (j) = Pm (j)+ppPL(i,k,j)
            P2m(j) = P2m(j)+ppPL(i,k,j)**2
         end do
      end do
   end do
   !!!!!!!!    wx stats:   !!!!!!!!
   do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      do k = 1,N(2,bandPL(myid))
         do i = 1,N(1,bandPL(myid))
            wxm (j) = wxm (j)+wx(i,k,j)
            wx2m(j) = wx2m(j)+wx(i,k,j)**2
         end do
      end do
   end do
   !!!!!!!!    UV stats:   !!!!!!!!
   do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
      do k = 1,kgal
         do i = 1,igal-2
            UVm(j) = UVm(j)+u1PL(i,k,j)*u2PL_itp(i,k,j)
         end do
      end do
   end do
   !!!!!!!!    UW stats:   !!!!!!!!
   do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
      do k = 1,kgal
         do i = 1,igal-2
            UWm(j) = UWm(j)+u1PL(i,k,j)*u3PL(i,k,j)
         end do
      end do
   end do
   !!!!!!!!    VW stats:   !!!!!!!!
   do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
      do k = 1,kgal
         do i = 1,igal-2
            VWm(j) = VWm(j)+u2PL(i,k,j)*u3PL_itp(i,k,j)
         end do
      end do
   end do

  
   if (bandPL(myid)==1) then
      call conddown1(UmC,U2mC,u1PL(1,1,limPL_incw(ugrid,1,myid)),ugrid,myid)
      call conddown1(VmC,V2mC,u2PL(1,1,limPL_incw(vgrid,1,myid)),vgrid,myid)
      call conddown1(WmC,W2mC,u3PL(1,1,limPL_incw(ugrid,1,myid)),ugrid,myid)
      call conddown1(PmC,P2mC,ppPL(1,1,limPL_incw(pgrid,1,myid)),pgrid,myid)
      call conddown1(wxmC,wx2mC,wx,vgrid,myid)
      call conddown2(UVmC,u1PL(1,1,limPL_incw(ugrid,1,myid)),u2PL_itp(1,1,limPL_incw(ugrid,1,myid)),ugrid,myid)
      call conddown2(UWmC,u1PL(1,1,limPL_incw(ugrid,1,myid)),u3PL    (1,1,limPL_incw(ugrid,1,myid)),ugrid,myid)
      call conddown2(VWmC,u2PL(1,1,limPL_incw(vgrid,1,myid)),u3PL_itp(1,1,limPL_incw(vgrid,1,myid)),vgrid,myid)
   else if (bandPL(myid)==nband) then
      call condup1  (UmC,U2mC,u1PL(1,1,limPL_incw(ugrid,1,myid)),1d0,ugrid,myid)
      call condup1  (VmC,V2mC,u2PL(1,1,limPL_incw(vgrid,1,myid)),-1d0,vgrid,myid)
      call condup1  (WmC,W2mC,u3PL(1,1,limPL_incw(ugrid,1,myid)),1d0,ugrid,myid)
      call condup1  (PmC,P2mC,ppPL(1,1,limPL_incw(pgrid,1,myid)),1d0,pgrid,myid)
      call condup1  (wxmC,wx2mC,wx,1d0,vgrid,myid)
      call condup2  (UVmC,u1PL(1,1,limPL_incw(ugrid,1,myid)),u2PL_itp(1,1,limPL_incw(ugrid,1,myid)),1d0,-1d0,ugrid,myid)
      call condup2  (UWmC,u1PL(1,1,limPL_incw(ugrid,1,myid)),u3PL    (1,1,limPL_incw(ugrid,1,myid)),1d0,1d0,ugrid,myid)
      call condup2  (VWmC,u2PL(1,1,limPL_incw(vgrid,1,myid)),u3PL_itp(1,1,limPL_incw(vgrid,1,myid)),-1d0,1d0,vgrid,myid)
   end if

end subroutine

subroutine conddown1(YmC,Y2mC,x,grid,myid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!    CONDDOWN1   !!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use declaration
   implicit none

   integer i,i1,k,k1,j,j1,ix,iz,grid,myid
   real(8) x   (igal,kgal,limPL_incw(grid,1,myid)           :limPL_incw(grid,2,myid))
   !   real(8) x   (N(1,bandPL(myid))+2,N(2,bandPL(myid)),limPL_incw(grid,1,myid)           :limPL_incw(grid,2,myid))
   real(8) YmC (dnx ,dnz ,limPL_incw(grid,1,myid)-Ny(grid,0):limPL_incw(grid,2,myid)-Ny(grid,0))
   real(8) Y2mC(dnx ,dnz ,limPL_incw(grid,1,myid)-Ny(grid,0):limPL_incw(grid,2,myid)-Ny(grid,0))

   do j = limPL_incw(grid,1,myid),limPL_incw(grid,2,myid)
      j1 = j-Ny(grid,0)
      do k = 1,dnz
         do i = 1,dnx
            do iz = 1,ntilez
               !k1 = indkor(k+dnz*(iz-1))
               k1 = k+dnz*(iz-1)
               do ix = 1,ntilex
                  !i1 = indior(i+dnx*(ix-1))
                  i1 = i+dnx*(ix-1)
                  YmC (i,k,j1) = YmC (i,k,j1)+x(i1,k1,j)
                  Y2mC(i,k,j1) = Y2mC(i,k,j1)+x(i1,k1,j)**2
               end do
            end do
         end do
      end do
   end do

end subroutine

subroutine conddown2(YmC,x1,x2,grid,myid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!    CONDDOWN2   !!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use declaration
   implicit none

   integer i,i1,k,k1,j,j1,ix,iz,grid,myid
   real(8) x1  (igal,kgal,limPL_incw(grid,1,myid)           :limPL_incw(grid,2,myid))
   real(8) x2  (igal,kgal,limPL_incw(grid,1,myid)           :limPL_incw(grid,2,myid))
   !   real(8) x1  (N(1,bandPL(myid))+2,N(2,bandPL(myid)),limPL_incw(grid,1,myid)           :limPL_incw(grid,2,myid))
   !   real(8) x2  (N(1,bandPL(myid))+2,N(2,bandPL(myid)),limPL_incw(grid,1,myid)           :limPL_incw(grid,2,myid))
   real(8) YmC (dnx ,dnz ,limPL_incw(grid,1,myid)-Ny(grid,0):limPL_incw(grid,2,myid)-Ny(grid,0))

   do j = limPL_incw(grid,1,myid),limPL_incw(grid,2,myid)
      j1 = j-Ny(grid,0)
      do k = 1,dnz
         do i = 1,dnx
            do iz = 1,ntilez
               !k1 = indkor(k+dnz*(iz-1))
               k1 = k+dnz*(iz-1)
               do ix = 1,ntilex
                  !i1 = indior(i+dnx*(ix-1))
                  i1 = i+dnx*(ix-1)
                  YmC(i,k,j1) = YmC(i,k,j1)+ x1(i1,k1,j) &
                     &                                     *x2(i1,k1,j)
               end do
            end do
         end do
      end do
   end do

end subroutine

subroutine condup1(YmC,Y2mC,x,xsign,grid,myid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!    CONDDOWN2   !!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use declaration
   implicit none

   integer i,i1,k,k1,j,j1,ix,iz,grid,myid
   real(8) x   (igal,kgal, limPL_incw(grid,1,myid)                 : limPL_incw(grid,2,myid))
   !real(8) x   (N(1,bandPL(myid))+2,N(2,bandPL(myid)),limPL_incw(grid,1,myid)           :limPL_incw(grid,2,myid))
   real(8) YmC (dnx ,dnz ,-limPL_incw(grid,2,myid)+Ny(grid,nband)+1:-limPL_incw(grid,1,myid)+Ny(grid,nband)+1)
   real(8) Y2mC(dnx ,dnz ,-limPL_incw(grid,2,myid)+Ny(grid,nband)+1:-limPL_incw(grid,1,myid)+Ny(grid,nband)+1)
   real(8) xsign

   do j = limPL_incw(grid,1,myid),limPL_incw(grid,2,myid)
      j1 = -j+Ny(grid,nband)+1
      do k = 1,dnz
         do i = 1,dnx
            do iz = 1,ntilez
               !k1 = indkor(k+dnz*(iz-1))
               !k1 = indkor(dnz*iz+2-k)                   ! Opposite direction to bottom wall
               k1 = k+dnz*(iz-1)
               do ix = 1,ntilex
                  !i1 = indior(i+dnx*(ix-1))
                  !i1 = indior(dnx*ix+2-i)                 ! Opposite direction to bottom wall
                  i1 = i+dnx*(ix-1)
                  YmC (i,k,j1) = YmC (i,k,j1)+x(i1,k1,j)*xsign
                  Y2mC(i,k,j1) = Y2mC(i,k,j1)+x(i1,k1,j)**2
               end do
            end do
         end do
      end do
   end do

end subroutine

subroutine condup2(YmC,x1,x2,xsign1,xsign2,grid,myid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!    CONDDOWN2   !!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use declaration
   implicit none

   integer i,i1,k,k1,j,j1,ix,iz,grid,myid
   real(8) x1  (igal,kgal, limPL_incw(grid,1,myid)                 : limPL_incw(grid,2,myid))
   real(8) x2  (igal,kgal, limPL_incw(grid,1,myid)                 : limPL_incw(grid,2,myid))
   !   real(8) x1  (N(1,bandPL(myid))+2,N(2,bandPL(myid)),limPL_incw(grid,1,myid)           :limPL_incw(grid,2,myid))
   !   real(8) x2  (N(1,bandPL(myid))+2,N(2,bandPL(myid)),limPL_incw(grid,1,myid)           :limPL_incw(grid,2,myid))
   real(8) YmC (dnx ,dnz ,-limPL_incw(grid,2,myid)+Ny(grid,nband)+1:-limPL_incw(grid,1,myid)+Ny(grid,nband)+1)
   real(8) xsign1
   real(8) xsign2

   do j = limPL_incw(grid,1,myid),limPL_incw(grid,2,myid)
      j1 = -j+Ny(grid,nband)+1
      do k = 1,dnz
         do i = 1,dnx
            do iz = 1,ntilez
               !k1 = indkor(k+dnz*(iz-1))
               !k1 = indkor(dnz*iz+2-k)                   ! Opposite direction to bottom wall
               k1 = k+dnz*(iz-1)
               do ix = 1,ntilex
                  !i1 = indior(i+dnx*(ix-1))
                  !i1 = indior(dnx*ix+2-i)                 ! Opposite direction to bottom wall
                  i1 = i+dnx*(ix-1)
                  YmC(i,k,j1) = YmC(i,k,j1)+ x1(i1,k1,j)*xsign1 &
                     &                                     *x2(i1,k1,j)*xsign2
               end do
            end do
         end do
      end do
   end do

end subroutine

subroutine write_stats(myid,status,ierr)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!! WRITE STATS !!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use declaration
   implicit none

   include 'mpif.h'             ! MPI variables
   integer status(MPI_STATUS_SIZE),ierr,myid

   real(8), allocatable:: xu(:),xv(:),xp(:),xCu(:,:,:),xCv(:,:,:),xCp(:,:,:),buffCu(:,:,:),buffCv(:,:,:),buffCp(:,:,:)
   integer, allocatable:: dummint(:)
   integer msizeu,msizev,msizep
   integer grid

   if (myid/=0) then   ! SLAVES
      call MPI_SEND(Um  ,limPL_incw(ugrid,2,myid)-limPL_incw(ugrid,1,myid)+1,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(U2m ,limPL_incw(ugrid,2,myid)-limPL_incw(ugrid,1,myid)+1,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(Vm  ,limPL_incw(vgrid,2,myid)-limPL_incw(vgrid,1,myid)+1,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(V2m ,limPL_incw(vgrid,2,myid)-limPL_incw(vgrid,1,myid)+1,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(Wm  ,limPL_incw(ugrid,2,myid)-limPL_incw(ugrid,1,myid)+1,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(W2m ,limPL_incw(ugrid,2,myid)-limPL_incw(ugrid,1,myid)+1,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(Pm  ,limPL_incw(pgrid,2,myid)-limPL_incw(pgrid,1,myid)+1,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(P2m ,limPL_incw(pgrid,2,myid)-limPL_incw(pgrid,1,myid)+1,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(UVm ,limPL_incw(ugrid,2,myid)-limPL_incw(ugrid,1,myid)+1,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(UWm ,limPL_incw(ugrid,2,myid)-limPL_incw(ugrid,1,myid)+1,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(VWm ,limPL_incw(vgrid,2,myid)-limPL_incw(vgrid,1,myid)+1,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(wxm ,limPL_incw(vgrid,2,myid)-limPL_incw(vgrid,1,myid)+1,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      call MPI_SEND(wx2m,limPL_incw(vgrid,2,myid)-limPL_incw(vgrid,1,myid)+1,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
      if (bandPL(myid)==1) then
         ! Conditional statistics
         msizeu = dnx*dnz*(limPL_incw(ugrid,2,myid)-limPL_incw(ugrid,1,myid)+1)
         msizev = dnx*dnz*(limPL_incw(vgrid,2,myid)-limPL_incw(vgrid,1,myid)+1)
         msizep = dnx*dnz*(limPL_incw(pgrid,2,myid)-limPL_incw(pgrid,1,myid)+1)
         if(msizeu>0)then
      
            if(myid<Ny(ugrid,1)-Ny(ugrid,0))then
               call MPI_SEND(UmC  ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
               call MPI_SEND(U2mC ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            endif
            if(myid<Ny(vgrid,1)-Ny(vgrid,0))then
               call MPI_SEND(VmC  ,msizev,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
               call MPI_SEND(V2mC ,msizev,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            endif
            if(myid<Ny(ugrid,1)-Ny(ugrid,0))then
               call MPI_SEND(WmC  ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
               call MPI_SEND(W2mC ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            endif
            if(myid<Ny(pgrid,1)-Ny(pgrid,0))then
               call MPI_SEND(PmC  ,msizep,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
               call MPI_SEND(P2mC ,msizep,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            endif
            if(myid<Ny(ugrid,1)-Ny(ugrid,0))then
               call MPI_SEND(UVmC ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
               call MPI_SEND(UWmC ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            endif
            if(myid<Ny(vgrid,1)-Ny(vgrid,0))then
               call MPI_SEND(VWmC ,msizev,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
               call MPI_SEND(wxmC ,msizev,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
               call MPI_SEND(wx2mC,msizev,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            endif

         endif
      elseif(bandPL(myid)==nband)then
         ! Conditional statistics
         msizeu = dnx*dnz*(limPL_incw(ugrid,2,myid)-limPL_incw(ugrid,1,myid)+1)
         msizev = dnx*dnz*(limPL_incw(vgrid,2,myid)-limPL_incw(vgrid,1,myid)+1)
         msizep = dnx*dnz*(limPL_incw(pgrid,2,myid)-limPL_incw(pgrid,1,myid)+1)
         if(msizeu>0)then
	
            !if(myid<Ny(ugrid,1)-Ny(ugrid,0))then
            call MPI_SEND(UmC  ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            call MPI_SEND(U2mC ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            !endif
            !if(myid<Ny(vgrid,1)-Ny(vgrid,0))then
            call MPI_SEND(VmC  ,msizev,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            call MPI_SEND(V2mC ,msizev,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            !endif
            !if(myid<Ny(ugrid,1)-Ny(ugrid,0))then
            call MPI_SEND(WmC  ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            call MPI_SEND(W2mC ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            !endif
            !if(myid<Ny(pgrid,1)-Ny(pgrid,0))then
            call MPI_SEND(PmC  ,msizep,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            call MPI_SEND(P2mC ,msizep,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            !endif
            !if(myid<Ny(ugrid,1)-Ny(ugrid,0))then
            call MPI_SEND(UVmC ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            call MPI_SEND(UWmC ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            !endif
            !if(myid<Ny(vgrid,1)-Ny(vgrid,0))then
            call MPI_SEND(VWmC ,msizev,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            call MPI_SEND(wxmC ,msizev,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            call MPI_SEND(wx2mC,msizev,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
         !endif
      
         endif
	
      end if
    
   else                ! MASTER
      fnameimb = trim(dirout)//'stats_'//ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
      open(10,file=fnameimb,form='unformatted')
      allocate(dummint(88))
      dummint = 0
      write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
      deallocate(dummint)
      write(10) Ngal
      write(10) yu,dthetai,dthdyu
      write(10) yv,dthetai,dthdyv
      write(10) istat
      allocate(xu(Ny(ugrid,0):Ny(ugrid,nband)+1))
      allocate(xv(Ny(vgrid,0):Ny(vgrid,nband)+1))
      allocate(xp(Ny(pgrid,0):Ny(pgrid,nband)+1))
      call recvX(xu,Um  ,ugrid,status,ierr)                       ! All procs
      call recvX(xu,U2m ,ugrid,status,ierr)
      call recvX(xv,Vm  ,vgrid,status,ierr)
      call recvX(xv,V2m ,vgrid,status,ierr)
      call recvX(xu,Wm  ,ugrid,status,ierr)
      call recvX(xu,W2m ,ugrid,status,ierr)
      call recvX(xp,Pm  ,pgrid,status,ierr)
      call recvX(xp,P2m ,pgrid,status,ierr)
      call recvX(xu,UVm ,ugrid,status,ierr)
      call recvX(xu,UWm ,ugrid,status,ierr)
      call recvX(xv,VWm ,vgrid,status,ierr)
      call recvX(xv,wxm ,vgrid,status,ierr)
      call recvX(xv,wx2m,vgrid,status,ierr)
      deallocate(xu)
      deallocate(xv)
      deallocate(xp)
      msizeu = dnx*dnz*(limPL_incw(ugrid,2,myid)-limPL_incw(ugrid,1,myid)+1)
      if(msizeu>0)then
         allocate(xCu  (dnx,dnz,0:Ny(ugrid,1)-Ny(ugrid,0)),buffCu(dnx,dnz,0:Ny(ugrid,1)-Ny(ugrid,0)))
         allocate(xCv  (dnx,dnz,0:Ny(vgrid,1)-Ny(vgrid,0)),buffCv(dnx,dnz,0:Ny(vgrid,1)-Ny(vgrid,0)))
         allocate(xCp  (dnx,dnz,0:Ny(pgrid,1)-Ny(pgrid,0)),buffCp(dnx,dnz,0:Ny(pgrid,1)-Ny(pgrid,0)))
         call recvXC(xCu ,UmC  ,buffCu ,dnx,dnz,ugrid,status,ierr)
         call recvXC(xCu ,U2mC ,buffCu ,dnx,dnz,ugrid,status,ierr)
         call recvXC(xCv ,VmC  ,buffCv ,dnx,dnz,vgrid,status,ierr)
         call recvXC(xCv ,V2mC ,buffCv ,dnx,dnz,vgrid,status,ierr)
         call recvXC(xCu ,WmC  ,buffCu ,dnx,dnz,ugrid,status,ierr)
         call recvXC(xCu ,W2mC ,buffCu ,dnx,dnz,ugrid,status,ierr)
         call recvXC(xCp ,PmC  ,buffCp ,dnx,dnz,pgrid,status,ierr)
         call recvXC(xCp ,P2mC ,buffCp ,dnx,dnz,pgrid,status,ierr)
         call recvXC(xCu ,UVmC ,buffCu ,dnx,dnz,ugrid,status,ierr)
         call recvXC(xCu ,UWmC ,buffCu ,dnx,dnz,ugrid,status,ierr)
         call recvXC(xCv ,VWmC ,buffCv ,dnx,dnz,vgrid,status,ierr)
         call recvXC(xCv ,wxmC ,buffCv ,dnx,dnz,vgrid,status,ierr)
         call recvXC(xCv ,wx2mC,buffCv ,dnx,dnz,vgrid,status,ierr)
         deallocate(xCu,buffCu)
         deallocate(xCv,buffCv)
         deallocate(xCp,buffCp)
      endif
      close(10)
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

   if (bandPL(myid)==1 .or. bandPL(myid)==nband) then
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

end subroutine

subroutine recvX(x,Xm,grid,status,ierr)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!    recvX    !!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use declaration
   implicit none

   include 'mpif.h'             ! MPI variables
   integer status(MPI_STATUS_SIZE),ierr

   integer j,iproc,grid
   integer msize
   real(8) x(Ny(grid,0):Ny(grid,nband)+1)
   real(8) Xm(limPL_incw(grid,1,0):limPL_incw(grid,2,0))

   x = 0d0
   do j = limPL_incw(grid,1,0),limPL_incw(grid,2,0)
      x(j) = Xm(j)
   end do
   do iproc = 1,np-1
      msize = limPL_incw(grid,2,iproc)-limPL_incw(grid,1,iproc)+1
      call MPI_RECV(x(limPL_incw(grid,1,iproc)),msize,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,status,ierr)
   end do

   write(10) x

end subroutine

subroutine recvXC(xC,XCm,buffC,mx,mz,grid,status,ierr)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!    recvX    !!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use declaration
   implicit none

   include 'mpif.h'             ! MPI variables
   integer status(MPI_STATUS_SIZE),ierr

   integer i,k,j,iproc,grid
   integer msize,mx,mz
   real(8) xC   (mx,mz,0:Ny(grid,1)-Ny(grid,0))
   real(8) buffC(mx,mz,0:Ny(grid,1)-Ny(grid,0))
   real(8) XCm  (mx,mz,limPL_incw(grid,1,0)-Ny(grid,0):limPL_incw(grid,2,0)-Ny(grid,0))
  
  
   xC = 0d0
   buffC = 0d0
   do j = limPL_incw(grid,1,0)-Ny(grid,0),limPL_incw(grid,2,0)-Ny(grid,0)
      do k = 1,mz
         do i = 1,mx
            xC(i,k,j) = XCm(i,k,j)
         end do
      end do
   end do
 
   do iproc = 1,min(procs(1)-1,Ny(grid,1)-Ny(grid,0)-1)
      msize = mx*mz*(limPL_incw(grid,2,iproc)-limPL_incw(grid,1,iproc)+1)
      call MPI_RECV(xC   (1,1, limPL_incw(grid,1,iproc)-Ny(grid,0)      ),msize,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,status,ierr)
   end do
   do iproc = procs(nband-1),np-1
      msize = mx*mz*(limPL_incw(grid,2,iproc)-limPL_incw(grid,1,iproc)+1)
      call MPI_RECV(buffC(1,1,-limPL_incw(grid,2,iproc)+Ny(grid,nband)+1),msize,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,status,ierr)
   end do
   xC = xC+buffC
   write(10) xC

end subroutine

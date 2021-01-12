subroutine spectra(u1,u2,u2_itp2,u3,p,myid)

  use declaration
  implicit none

  integer myid
  type(cfield)  u1(sband:eband), u2(sband:eband), u2_itp2(sband:eband), u3(sband:eband)
  type(cfield)  p (sband:eband)

  call buildsp  (spU ,u1(midband)%f,                   ugrid,myid)
  call buildsp  (spV ,u2(midband)%f,                   vgrid,myid)
  call buildsp  (spW ,u3(midband)%f,                   ugrid,myid)
  call buildspUV(spUV,u1(midband)%f,u2_itp2(midband)%f,      myid)
  call buildsp  (spP ,p (midband)%f,                   pgrid,myid)

end subroutine

subroutine buildsp(sp,u,grid,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!    buildD   !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer j,column,grid,myid
  real(8)    sp(jlim(1,grid,midband):jlim(2,grid,midband),columns_num(midband,myid))
  complex(8) u (jlim(1,grid,midband):jlim(2,grid,midband),columns_num(midband,myid))

  do column = 1,columns_num(midband,myid)
    do j = jlim(1,grid,midband),jlim(2,grid,midband)
      sp(j,column) = sp(j,column)+dreal(u(j,column))**2+dimag(u(j,column))**2
    end do
  end do

end subroutine

subroutine buildspUV(sp,u,v,myid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!    buildD   !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none
  integer j,column,grid,myid
  real(8)    sp(jlim(1,ugrid,midband):jlim(2,ugrid,midband),columns_num(midband,myid))
  complex(8) u (jlim(1,ugrid,midband):jlim(2,ugrid,midband),columns_num(midband,myid))
  complex(8) v (jlim(1,ugrid,midband):jlim(2,ugrid,midband),columns_num(midband,myid))

  do column = 1,columns_num(midband,myid)
    do j = jlim(1,ugrid,midband),jlim(2,ugrid,midband)
        sp(j,column) = sp(j,column)+dreal(u(j,column))*dreal(v(j,column)) &
&                                  +dimag(u(j,column))*dimag(v(j,column))
    end do
  end do

end subroutine

subroutine write_spect(myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! WRITE SPECT !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid

  real(8), allocatable:: buffSPu(:,:,:),buffSPv(:,:,:),buffSPp(:,:,:)
  integer msizeu,msizev,msizep

  if (myid/=0) then
    msizeu = (jlim(2,ugrid,midband)-jlim(1,ugrid,midband)+1)*columns_num(midband,myid)
    msizev = (jlim(2,vgrid,midband)-jlim(1,vgrid,midband)+1)*columns_num(midband,myid)
    msizep = (jlim(2,pgrid,midband)-jlim(1,pgrid,midband)+1)*columns_num(midband,myid)
    call MPI_SEND(spU ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(spV ,msizev,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(spW ,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(spUV,msizeu,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(spP ,msizep,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
  else
    allocate(buffSPu(0:N(1,midband)/2,1:N(2,midband)/2+1,jlim(1,ugrid,midband):jlim(2,ugrid,midband)))
    allocate(buffSPv(0:N(1,midband)/2,1:N(2,midband)/2+1,jlim(1,vgrid,midband):jlim(2,vgrid,midband)))
    allocate(buffSPp(0:N(1,midband)/2,1:N(2,midband)/2+1,jlim(1,pgrid,midband):jlim(2,pgrid,midband)))

    fnameimb=trim(dirout)//'/spU_' //ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
    call recvwrspec(spU ,buffSPu,ugrid,myid,status,ierr)
    fnameimb=trim(dirout)//'/spV_' //ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
    call recvwrspec(spV ,buffSPv,vgrid,myid,status,ierr)
    fnameimb=trim(dirout)//'/spW_' //ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
    call recvwrspec(spW ,buffSPu,ugrid,myid,status,ierr)
    fnameimb=trim(dirout)//'/spUV_'//ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
    call recvwrspec(spUV,buffSPu,ugrid,myid,status,ierr)
    fnameimb=trim(dirout)//'/spP_' //ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
    call recvwrspec(spP ,buffSPp,pgrid,myid,status,ierr)
    deallocate(buffSPu)
    deallocate(buffSPv)
    deallocate(buffSPp)
  end if

  spU  = 0d0
  spV  = 0d0
  spW  = 0d0
  spUV = 0d0
  spP  = 0d0

end subroutine

subroutine recvwrspec(spX,buffSP,grid,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!    recv+write spec    !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr

  integer i,k,kk,j,iproc,column,grid,myid
  integer msize
  real(8)    spX(jlim(1,grid,midband):jlim(2,grid,midband),columns_num(midband,myid))
  !The following +1 in 1:N(2,midband)/2+1 shouldn't be there but is consistent with the postprocessing....
  real(8) buffSP(0:N(1,midband)/2,1:N(2,midband)/2+1,jlim(1,grid,midband):jlim(2,grid,midband))
  real(8), allocatable:: buffrecv(:,:)
  integer, allocatable:: dummint(:)

  buffSP = 0d0

  do j = jlim(1,grid,midband),jlim(2,grid,midband)
    do column = 1,columns_num(midband,myid)
      i = columns_i(column,midband,myid)
      k = columns_k(column,midband,myid)
      if (k > N(2,1)/2) then
        kk = N(2,1)+2-k
      else
        kk = k
      end if
      buffSP(i,kk,j) = spX(j,column)
    end do
  end do

  do iproc = 1,np-1
    msize = (jlim(2,grid,midband)-jlim(1,grid,midband)+1)*columns_num(midband,iproc)
    allocate(buffrecv(jlim(1,grid,midband):jlim(2,grid,midband),columns_num(midband,iproc)))
    call MPI_RECV(buffrecv,msize,MPI_REAL8,iproc,iproc,MPI_COMM_WORLD,status,ierr)
    do j = jlim(1,grid,midband),jlim(2,grid,midband)
      do column = 1,columns_num(midband,iproc)
        i = columns_i(column,midband,iproc)
        k = columns_k(column,midband,iproc)
        if (k > N(2,1)/2) then
          kk = N(2,1)+2-k
        else
          kk = k
        end if
        buffSP(i,kk,j) = buffSP(i,kk,j)+buffrecv(j,column)
      end do
    end do
    deallocate(buffrecv)
  end do

  do j = jlim(1,grid,midband),jlim(2,grid,midband)
    do k = 1,N(2,midband)/2+1
      do i = 1,N(1,midband)/2
        buffSP(i,k,j) = 2d0*buffSP(i,k,j)
      end do
    end do
  end do

  open(10,file=fnameimb,form='unformatted')
  allocate(dummint(88))
  dummint = 0
  write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
  deallocate(dummint)
  write(10) N
  write(10) yu,dthetai,dthdyu
  write(10) yv,dthetai,dthdyv
  write(10) istat
  write(10) buffSP
  close(10)

end subroutine

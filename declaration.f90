module declaration

  ! Static variables
  
  real(8), parameter :: pi = 4d0*datan(1d0)
  complex, parameter :: im = dcmplx(0d0,1d0)
  
  integer, parameter :: nband = 3
  integer, parameter :: sband = 1
  integer, parameter :: eband = nband
  integer, parameter :: botband = 1           ! 1
  integer, parameter :: midband = (nband+1)/2 ! 2
  integer, parameter :: topband = nband       ! 3

  integer, parameter :: ugrid = 2 ! centres     including ghost points: u and w
  integer, parameter :: vgrid = 1 ! faces: v
  integer, parameter :: pgrid = 3 ! centres not including ghost points: p

  
  ! All variables

  integer bandit(3)
  integer np,pnodes
  integer,pointer:: procs(:),procs_b(:),procs_c(:)
  integer,pointer:: N(:,:),Ngal(:,:),Ny(:,:)
  integer,pointer:: jlim(:,:,:) ! Change name ffs
  integer,pointer:: planelim(:,:,:)
  integer,pointer:: limPL_incw(:,:,:),limPL_excw(:,:,:)
  integer,pointer:: limPL_FFT(:,:,:)
  integer,pointer:: bandPL(:)
  integer,pointer:: bandPL_FFT(:)
  integer           jgal(3,2),igal,kgal !TODO get rid of jgal
  integer,pointer:: dk(:,:,:,:)
  integer,pointer:: columns_i(:,:,:)
  integer,pointer:: columns_k(:,:,:)
  integer,pointer:: columns_num(:,:)
  integer,pointer:: dk_phys(:,:,:,:)
  
  integer,pointer:: planeBC(:,:)
    
  integer physlim_bot
  integer physlim_top
  
  integer nn
  integer iter,iter0,nwrite,iwrite,itersl,nstatsl
  integer nstat,istat
  integer nstat_sl,istat_sl
  integer flag_init,flag_ctpress
  real(8) nextqt
  integer geometry_type
  integer kRK
  real(8) t,dt,CFL,maxt,dtv,dtc,dti,Re

  integer nribs,npeak !TODO remove
  integer dnx,dnz
  integer ntilex,ntilez,dsty,dstx,dstz
  integer npeakx,npeakz
  real(8) Lfracx,Lfracz
  real(8) dyub2,dyut2,dyvb2,dyvt2
  real(8) posth
  real(8) post_spacing
  real(8) shift_stag

  real(8) alp,bet             ! Wavelengths
  real(8) Lx,Ly,Lz            ! Size of the computational box
  real(8),pointer:: h_ny(:)
  real(8) Kib

  ! Grid
  real(8),pointer:: yu(:),dyu2i(:,:),dthdyu(:)
  real(8),pointer:: yv(:),dyv2i(:,:),dthdyv(:)
  real(8) dtheta,dthetai
  real(8) dthetavi,ddthetavi
  real(8),pointer:: gridweighting(:,:)
  real(8),pointer:: gridweighting_interp(:,:)
  integer ppp
  real(8) dyq

  ! Variables in planes
  real(8),pointer::  u1PL(:,:,:), u2PL(:,:,:), u3PL(:,:,:)
  real(8),pointer::  u1PL_itp(:,:,:), u2PL_itp(:,:,:), u3PL_itp(:,:,:)
  real(8),pointer:: Nu1PL(:,:,:),Nu2PL(:,:,:),Nu3PL(:,:,:)
  real(8),pointer:: du1PL(:,:,:),du2PL(:,:,:),du3PL(:,:,:)
  real(8),pointer::    wx(:,:,:), ppPL(:,:,:)
  real(8),pointer:: div_outPL(:,:,:),u_outPL(:,:,:),v_outPL(:,:,:),w_outPL(:,:,:)
  real(8),pointer:: div_cPL(:,:,:) 
 
  real(8),pointer:: Qcrit(:,:,:)

  real(8),pointer::  u1PLN(:,:,:), u2PLN(:,:,:), u3PLN(:,:,:), ppPLN(:,:,:)
  real(8),pointer::  u1PL_itpN(:,:,:), u2PL_itpN(:,:,:), u3PL_itpN(:,:,:)

  ! Cross products in planes
  real(8),pointer::  uu_cPL (:,:,:), uv_fPL (:,:,:), uw_cPL (:,:,:)
  real(8),pointer::  vv_cPL (:,:,:), vw_fPL (:,:,:), ww_cPL (:,:,:)
  real(8),pointer:: Nu1PL_dy(:,:,:),Nu2PL_dy(:,:,:),Nu3PL_dy(:,:,:)

  ! Spectra
  real(8),pointer::  spU(:,:), spV(:,:), spW(:,:)
  real(8),pointer:: spUV(:,:), spP(:,:)

  ! Statistics. Mean and conditional statistics
  real(8),pointer:: Um(:),U2m(:),UmC(:,:,:),U2mC(:,:,:)
  real(8),pointer:: Vm(:),V2m(:),VmC(:,:,:),V2mC(:,:,:)
  real(8),pointer:: Wm(:),W2m(:),WmC(:,:,:),W2mC(:,:,:)
  real(8),pointer:: Pm(:),P2m(:),PmC(:,:,:),P2mC(:,:,:)
  real(8),pointer:: UVm(:),UVmC(:,:,:)
  real(8),pointer:: UWm(:),UWmC(:,:,:)
  real(8),pointer:: VWm(:),VWmC(:,:,:)
  real(8),pointer:: wxm(:),wx2m(:),wxmC(:,:,:),wx2mC(:,:,:)
  integer,pointer:: indkor(:),indior(:)

  real(8),pointer:: u11(:)
  
  complex(8),pointer:: k1F_x(:),k1F_z(:)
  real(8)   ,pointer:: k2F_x(:),k2F_z(:)

  ! Runge-Kutta coefficients
  real(8) aRK(3),bRK(3),gRK(3),cRK(3),dRK(3)

  real(8) err,maxerr,maxA!,lambda,lambdaQ,Re_div,iRediv
  real(8) mpgx,mpgz,dgx,dgz,QxT,Qx,Qz,Umax,utau
  character*120 fnameima,fnameimb,fnameimc,boundfname,filout,directory
  logical exist_file_hist
  character*120 dirin,dirout
  character*4 ext1,ext2,ext3
  character*5 ext4

  integer, pointer :: nlist_ib_s(:), nlist_ib_f(:), nlist_ib(:)
  integer, pointer :: s_list_ib(:,:,:), f_list_ib(:,:,:), list_ib(:,:,:)
  real(8), pointer :: w_list_ib(:,:,:)
  integer, pointer :: nyuIB1(:),nyuIB2(:),nyvIB1(:),nyvIB2(:)
  integer          :: nyu11,nyu21,nyu12,nyu22,nyv11,nyv21,nyv12,nyv22
  real(8), pointer :: A_ib(:,:,:)

  type cfield
    complex(8),pointer:: f(:,:)
  end type cfield
  
  type rfield_dg
    real(8),pointer:: f_dg(:,:,:)
  end type rfield_dg

  type array
    real(8),pointer:: b(:)
  end type array

  type(array), pointer:: buffR_x(:)
  type(array), pointer:: buffRal_x(:)
  type(array), pointer:: buffC_z(:)
  type(array), pointer:: buffCal_z(:)
  
  
  type(cfield), allocatable:: u1_itp(:),u2_itp(:),u3_itp(:)
  type(cfield), allocatable:: Nu1_dy(:),Nu2_dy(:),Nu3_dy(:)
  type(cfield), allocatable:: uv_f(:), vw_f(:), vv_c(:)

  type(cfield), allocatable:: div_c(:),div_out(:),u_out(:),v_out(:),w_out(:) 
 
  ! Omega x
  real(8),      allocatable :: du1dy_planes(:,:,:)
  real(8),      allocatable :: du2dy_planes(:,:,:)
  real(8),      allocatable :: du3dy_planes(:,:,:)
  
  real(8),      allocatable :: du1dy_planes2(:,:,:)
  real(8),      allocatable :: du2dy_planes2(:,:,:)
  real(8),      allocatable :: du3dy_planes2(:,:,:)
  
  type(cfield), allocatable :: du1dy_columns(:)
  type(cfield), allocatable :: du2dy_columns(:)
  type(cfield), allocatable :: du3dy_columns(:)
  
  type(rfield_dg), allocatable:: DG(:)
  
  real(8) bslip
  
end module

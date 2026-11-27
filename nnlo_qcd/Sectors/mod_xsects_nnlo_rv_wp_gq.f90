module mod_xsects_nnlo_rv_wp_gq
  use mod_types; use common_def; use mod_consts_dp
#if (_DebugFeatures == 1)
  use mod_consts_debug
#endif
  use mod_process
  use mod_kinematics_gen; use mod_kinematicsLO; use mod_kinematicsNLO
  use mod_histograms
  use mod_lumi
  use mod_auxlumi_wp
  use mod_auxfunctions; use mod_procfunctions
  use mod_amplitudes_tree_wp
  use mod_amplitudes_loop_wp
  implicit none
  integer, parameter :: namp_lo   = 2

  private

  !-- notation: xsect_nnlo_rv_proc_ch
  public :: xsect_nnlo_rv_wp_gq

contains

  function xsect_nnlo_rv_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    use mod_limits_rv_dy_qg
    integer    :: xsect_nnlo_rv_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,SoftLim,Col1Lim,Col2Lim
    real(dp)   :: xx(kNLO_max_full)
    real(dp)   :: tau,ylab,jac,mu,muR,muF,muRsq,weightarray(nobs)
    real(dp)   :: xi1,xi2,f1(-6:6),f2(-6:6)
    real(dp)   :: res(namp_lo),res_loop(namp_lo),res_loop_j(2),respdf,respdf_gq,respdf_loop
    real(dp)   :: z,eta41,eta42,s12,s14,s24,E4oE1,E4oE2
    real(dp)   :: limits_rv_soft(2),limits_rv_col1(2),limits_rv_col2(2),limits_rv_softcol1(2),limits_rv_softcol2(2)
    real(dp)   :: Fint_rv(2),kin(2)
    logical    :: makecut
    real(dp)   :: pref
    integer    :: i
    real(dp)   :: my_lumi,q2,Lmu

    callcount = callcount + 1

    xsect_nnlo_rv_wp_gq = 0
    xx = zero
    
    HardProc=>MyKinSet%AllKinConfig(1)
    SoftLim=>MyKinSet%AllKinConfig(2)
    Col1Lim=>MyKinSet%AllKinConfig(3)
    Col2Lim=>MyKinSet%AllKinConfig(4)
 
   xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    call kinematics_nlo(tau,ylab,xx(kNLO_min:kNLO_max_full),HardProc,SoftLim,Col1Lim,Col2Lim)
      
! hard emission 
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)

    call cut_histo_w(xi1,xi2,5,HardProc%AmpMom,makecut,binarray,weightarray,mu)
    if (makecut .or. HardProc%flag) then
       kin(1)  = zero    
       Fint_rv(1) = zero
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       muRsq = muR**2
       call res_1loop_j_gq(HardProc%AmpMom,muRsq,res_loop_j)
       respdf_loop = gq_lumi_rv(res_loop_j,f1,f2)
       Fint_rv(1) = respdf_loop * jac * pref * HardProc%wgt
       kin(1) = Fint_rv(1)
    endif

    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) 

! col 1
    xi1 = Col1Lim%PartFrac(1)
    xi2 = Col1Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,Col1Lim%AmpMom,makecut,binarray,weightarray,mu)
    if (makecut .or. Col1Lim%flag) then
       kin(2) = zero    
       Fint_rv(2) = zero
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       muRsq = muR**2
       call res_1loop(Col1Lim%AmpMom,muR**2,res,res_loop)
       respdf_gq = gq_lumi_rv_lo(res,f1,f2)
       respdf_loop = gq_lumi_rv_lo(res_loop,f1,f2)
       z = Col1Lim%KinInv(1)
       s14 = Col1Lim%KinInv(2)
       eta41 = Col1Lim%KinInv(2)/four/Col1Lim%LimMom(1,1)/Col1Lim%LimMom(1,3)*z
       q2 = Col1Lim%KinInv(3)

       Lmu = log(muR**2/q2)
       limits_rv_col1=col_rv(Lmu,z,s14,eta41)

       Fint_rv(2) = -( limits_rv_col1(1)*respdf_gq + limits_rv_col1(2)*respdf_loop)* jac * pref * Col1Lim%wgt
       kin(2) = Fint_rv(2)
       

    endif
    call HistoFill(10,binarray,weightarray,kin(2),vegasweight,iternumber)

    ff(1) = sum(kin)

#if (_DebugFeatures == 1)
    Fint_c41only_rv_CheckLim = Fint_rv 
#endif
    
    return

  end function xsect_nnlo_rv_wp_gq

end module mod_xsects_nnlo_rv_wp_gq

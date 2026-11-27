module mod_xsects_nnlo_vv_wp_ns
  use mod_types; use common_def; use mod_consts_dp
  use mod_kinematics_gen; use mod_kinematicsLO
  use mod_process
  use mod_auxfunctions
  use mod_procfunctions
  use mod_histograms
  use mod_lumi
  use mod_auxlumi_wp
  use mod_amplitudes_loop_wp
#if (_DebugFeatures == 1) 
  use mod_consts_debug
#endif
  implicit none
  integer, parameter :: namp_lo = 4
  
  private

  public :: xsect_nnlo_vv_wp_ns

contains

  function xsect_nnlo_vv_wp_ns(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_vv_wp_ns,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: LOProc
    real(dp) :: xx(kLO_max_full)
    real(dp) :: tau,ylab,jac,xi1,xi2,mu,muR,muF,muRsq,weightarray(nobs)
    real(dp) :: f1(-6:6),f2(-6:6),shat
    real(dp) :: res(namp_lo),res0(namp_lo),res1sq(namp_lo),res2(namp_lo),respdf,FintVV
    logical :: makecut
    real(dp) :: pref
    
#if (_DebugFeatures)
    callcount = callcount + 1
#endif

    xsect_nnlo_vv_wp_ns = 0
    ff(1) = zero

    LOProc=>MyKinSet%AllKinConfig(1)

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    
    !-- FintVV below
    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    xi1 = sqrt(tau)*exp(+ylab)
    xi2 = sqrt(tau)*exp(-ylab)
    call kinematics_lo(xi1,xi2,xx(kLO_min:kLO_max_full),LOProc) !-- CHANGED HERE

    !-- cut-histo
    call cut_histo_w(xi1,xi2,4,LOProc%AmpMom,makecut,binarray,weightarray,mu)
    if (makecut .or. LOProc%flag) return
    call setscale(mu,muF,muR)
    muRsq = muR**2

    call get_pdf(xi1,xi2,muF,f1,f2)
    call get_prefactor(muR,2,pref)

    call res_1loop_sq(LOProc%AmpMom,muRsq,res0,res1sq)
    call res_2loop(LOProc%AmpMom,muRsq,res0,res2)
    res = res1sq + res2
    respdf = ns_lumi_lo(res,f1,f2)

    FintVV = respdf*LOProc%wgt*jac*pref

    call HistoFill(0,binarray,weightarray,FintVV,vegasweight,iternumber)
    call HistoFill(10,binarray,weightarray,zero,vegasweight,iternumber)

    ff(1) = FintVV

    return

  end function xsect_nnlo_vv_wp_ns

end module mod_xsects_nnlo_vv_wp_ns

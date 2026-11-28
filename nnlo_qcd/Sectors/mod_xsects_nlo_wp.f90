module mod_xsects_nlo
  use mod_types; use common_def; use mod_consts_dp
#if (_DebugFeatures == 1)
  use mod_consts_debug
#endif
  use mod_kinematics_gen; use mod_kinematicsLO;  use mod_kinematicsNLO
  use mod_process
  use mod_auxfunctions; use mod_procfunctions
  use mod_histograms
  use mod_lumi
  use mod_auxlumi_wp
  use mod_amplitudes_tree_wp
  use mod_amplitudes_loop_wp
  use mod_limits_dy_ns
  use mod_splitting
  use mod_splitting_aux
  implicit none
  integer, parameter :: namp_lo   = 2
  integer, parameter :: namp_nlo  = 6
 
  private

  public :: xsect_nlo_r, xsect_nlo_v, xsect_nlo_s

contains

  function xsect_nlo_r(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nlo_r,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,SoftLim,Col1Lim,Col2Lim
    real(dp) :: xx(kNLO_max_full)
    real(dp) :: tau,ylab,jac,mu,muR,muF,weightarray(nobs)
    real(dp) :: xi1,xi2,f1(-6:6),f2(-6:6),res(namp_lo),res_j(namp_nlo),respdf,FintNLO(4)
    real(dp) :: respdf_qqb,respdf_qg,respdf_gq
    complex(dp) ::  res_sping_j(6,-1:1,-1:1)
    integer  :: i,j
    real(dp) :: reseik,z,s14,s24
    logical :: makecut
    real(dp) :: pref

    callcount = callcount + 1

    xsect_nlo_r = 0
    ff(1) = zero
    xx = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    SoftLim=>MyKinSet%AllKinConfig(2)
    Col1Lim=>MyKinSet%AllKinConfig(3)
    Col2Lim=>MyKinSet%AllKinConfig(4)

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    call kinematics_nlo(tau,ylab,xx(kNLO_min:kNLO_max_full),HardProc,SoftLim,Col1Lim,Col2Lim)

    !-- regular piece
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut .or. HardProc%flag) then

       FintNLO(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,1,pref)

       call restree_j(HardProc%AmpMom,res_sping_j,res_j)

       if (do_nlo_qqb.eq.1) then; 
          respdf_qqb = ns_lumi_nlo(res_j,f1,f2)
       else 
          respdf_qqb = zero
       endif
       
       if (do_nlo_qg.eq.1) then; 
          respdf_qg = qg_lumi_nlo(res_j,f1,f2)
          respdf_gq = gq_lumi_nlo(res_j,f1,f2)
       else 
          respdf_qg = zero
          respdf_gq = zero
       endif

       respdf = respdf_qqb + respdf_qg + respdf_gq
       
       FintNLO(1) = respdf*HardProc%wgt*jac * pref

    endif

    call HistoFill(0,binarray,weightarray,FintNLO(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    
    !-- soft piece, not needed, but just for fun
    xi1 = SoftLim%PartFrac(1)
    xi2 = SoftLim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,SoftLim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut .or. SoftLim%flag .or. do_nlo_qqb.eq.0) then

       FintNLO(2) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,1,pref)

       call restree(SoftLim%AmpMom,res)
       respdf = ns_lumi_lo(res,f1,f2)

       reseik = Eik_s4(SoftLim%LimMom)

       FintNLO(2) = -respdf*reseik*SoftLim%wgt*jac * pref

    endif

    !-- collinear 1 piece
    xi1 = Col1Lim%PartFrac(1)
    xi2 = Col1Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,Col1Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut .or. Col1Lim%flag) then

       FintNLO(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,1,pref)

       call restree(Col1Lim%AmpMom,res)

       z = Col1Lim%KinInv(1)
       s14 = Col1Lim%KinInv(2)

       if (do_nlo_qqb.eq.1) then; 
          respdf_qqb = ns_lumi_lo(res,f1,f2)
          respdf_qqb = respdf_qqb * (two/s14) * Pqq_i(z)
       else 
          respdf_qqb = zero
       endif
       
       if (do_nlo_qg.eq.1) then; 
          respdf_gq = gq_lumi_lo(res,f1,f2)
          respdf_gq = respdf_gq * (two/s14) * Pqg_i(z)
       else 
          respdf_gq = zero

       endif

       respdf = respdf_qqb + respdf_gq 

       FintNLO(3) = -respdf*Col1Lim%wgt*jac * pref

       call HistoFill(1,binarray,weightarray,FintNLO(3),vegasweight,iternumber)
    endif

    !-- collinear 2 piece
    xi1 = Col2Lim%PartFrac(1)
    xi2 = Col2Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,Col2Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut .or. Col2Lim%flag) then

       FintNLO(4) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,1,pref)

       call restree(Col2Lim%AmpMom,res)

       z = Col2Lim%KinInv(1)
       s24 = Col2Lim%KinInv(2)

       if (do_nlo_qqb.eq.1) then; 
          respdf_qqb = ns_lumi_lo(res,f1,f2)
          respdf_qqb = respdf_qqb * (two/s24) * Pqq_i(z)
       else 
          respdf_qqb = zero
       endif
       
       if (do_nlo_qg.eq.1) then; 
          respdf_qg = qg_lumi_lo(res,f1,f2)
          respdf_qg = respdf_qg * (two/s24) * Pqg_i(z)
       else 
          respdf_qg = zero
       endif

       respdf = respdf_qqb + respdf_qg

       FintNLO(4) = -respdf*Col2Lim%wgt*jac*pref

    endif

    call HistoFill(10,binarray,weightarray,FintNLO(4),vegasweight,iternumber)

    ff(1) = FintNLO(1) + FintNLO(3) + FintNLO(4)

#if (_DebugFeatures == 1)
    FintNLO_CheckLim = FintNLO
#endif
    
!     print *, 'FintNLO(1)',FintNLO(1)
!     print *, 'FintNLO(2)',FintNLO(2)
!     print *, 'FintNLO(3)',FintNLO(3)
!     print *, 'FintNLO(4)',FintNLO(4)
!     print *, 'ff(1)',ff(1); pause

    return

  end function xsect_nlo_r

  
  function xsect_nlo_v(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nlo_v,ndim,ncomp,userdata,iternumber,binarray(nobs)
    type(KinConfig),pointer :: VProc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    real(dp) :: xx(kLO_max_full)
    real(dp) :: tau,ylab,jac,mu,muR,muF,muRsq,weightarray(nobs)
    real(dp) :: xi1,xi2,f1(-6:6),f2(-6:6),shat,res0(4),res(4),respdf,FintV
    logical :: makecut
    real(dp) :: pref

    callcount = callcount + 1

    xsect_nlo_v = 0
    ff(1) = zero
    xx = zero

    if (do_nlo_qqb .eq. 0) then

       FintV = zero

    else
       
       VProc=>MyKinSet%AllKinConfig(1)

       xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
       call random_number(xx(kLO_max_full))

       !-- FintV below
       call get_tauy(xx(1:tauy_max),tau,ylab,jac)
       xi1 = sqrt(tau)*exp(+ylab)
       xi2 = sqrt(tau)*exp(-ylab)
       call kinematics_lo(xi1,xi2,xx(kLO_min:kLO_max_full),VProc)

       call cut_histo_w(xi1,xi2,4,VProc%AmpMom,makecut,binarray,weightarray,mu)
       if (makecut) return

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,1,pref)

       muRsq = muR**2
       call res_1loop(VProc%AmpMom,mursq,res0,res)
       respdf = ns_lumi_lo(res,f1,f2)

       FintV = respdf*VProc%wgt*jac*pref

    endif

    call HistoFill(0,binarray,weightarray,FintV,vegasweight,iternumber)
    call HistoFill(10,binarray,weightarray,zero,vegasweight,iternumber)

    ff(1) = FintV

    !print *, 'ff(1)',ff(1); pause !-- note different mapping wrt old code

    return

  end function xsect_nlo_v

 
  function xsect_nlo_s(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nlo_s,ndim,ncomp,userdata,iternumber,binarray(nobs)
    type(KinConfig),pointer :: SubProc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    real(dp) :: xx(kLO_max_full),z
    real(dp) :: tau,ylab,jac,mu,muR,muF,muFsq,weightarray(nobs)
    real(dp) :: xi1,xi2,f1(-6:6),f2(-6:6),res(4),my_lumi,respdf,FintSub(3)
    logical :: makecut
    real(dp) :: pref,q2
    real(dp) :: respdf_qqb,respdf_qg,respdf_gq
    real(dp) :: subres_qqb(-1:1),subres_qg(-1:1)
    integer :: i

    callcount = callcount + 1

    xsect_nlo_s = 0
    ff(1) = zero
    xx = zero

    SubProc=>MyKinSet%AllKinConfig(1)

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))

    z=buff+onet*real(yRnd(kLO_max_full),dp)
    call get_tauyxi1xi2(xx(1:tauy_max),one,one,tau,ylab,xi1,xi2,jac) !-- tau = Minv^2/s_had
    q2 = tau*sh

    call kinematics_lo(xi1,xi2,xx(kLO_min:kLO_max_full),SubProc)

    call restree(SubProc%AmpMom,res)

    call cut_histo_w(xi1,xi2,4,SubProc%AmpMom,makecut,binarray,weightarray,mu)
    if (makecut) return

    call setscale(mu,muF,muR)
    call get_prefactor(muR,1,pref)
    muFsq = muF**2
    if (do_nlo_qqb .eq. 1) then
       subres_qqb = -calPp_qq(z) - PqqAP(z) * log(muFsq/q2)   
       subres_qqb(-1) = subres_qqb(-1) + one/three * Cf * pisq !-- remainder from virtual definition. 
                                                             !-- Half w.r.t. notes because multiplied by #legs below
    endif

    if (do_nlo_qg .eq. 1) then
       subres_qg = -calPp_qg(z) - PqgAP(z) * log(muFsq/q2)
    endif

    !-- 12

    if (do_nlo_qqb .eq. 1) then

       call get_pdf(xi1,xi2,muF,f1,f2)
       respdf = ns_lumi_lo(res,f1,f2)

       FintSub(1) = respdf*SubProc%wgt*jac*pref * (subres_qqb(-1) &
            -subres_qqb(+1)) !-- plus term 
       FintSub(1) = FintSub(1) * two !-- two legs

    else

       FintSub(1) = zero
       
    endif
    call HistoFill(0,binarray,weightarray,FintSub(1),vegasweight,iternumber)
           
    !-- (z,2)
    call get_tauyxi1xi2(xx(1:tauy_max),z,one,tau,ylab,xi1,xi2,jac) !-- tau = Minv^2/s_had
    call get_pdf(xi1,xi2,muF,f1,f2)

    if (do_nlo_qqb .eq. 1) then

       respdf_qqb = ns_lumi_lo(res,f1,f2)
       respdf_qqb = respdf_qqb * subres_qqb(0)

    else

       respdf_qqb = zero

    endif

    if (do_nlo_qg .eq. 1) then

       respdf_gq = gq_lumi_lo(res,f1,f2)
       respdf_gq = respdf_gq * subres_qg(0)

    else

       respdf_gq = zero

    endif

    respdf = respdf_qqb + respdf_gq

    FintSub(2) = respdf*SubProc%wgt*jac*pref
    call HistoFill(1,binarray,weightarray,FintSub(2),vegasweight,iternumber)


    !-- (1,z)
    call get_tauyxi1xi2(xx(1:tauy_max),one,z,tau,ylab,xi1,xi2,jac) !-- tau = Minv^2/s_had
    call get_pdf(xi1,xi2,muF,f1,f2)

    if (do_nlo_qqb .eq. 1) then

       respdf_qqb = ns_lumi_lo(res,f1,f2)
       respdf_qqb = respdf_qqb * subres_qqb(0) 
       
    else

       respdf_qqb = zero

    endif

    if (do_nlo_qg .eq. 1) then

       respdf_qg = qg_lumi_lo(res,f1,f2)
       respdf_qg = respdf_qg * subres_qg(0)

    else

       respdf_qg = zero

    endif

    respdf = respdf_qqb + respdf_qg

    FintSub(3) = respdf*SubProc%wgt*jac*pref
    call HistoFill(10,binarray,weightarray,FintSub(3),vegasweight,iternumber)

    ff(1) = FintSub(1) + FintSub(2) + FintSub(3)

!     print *, 'FintSub(1)', FintSub(1)
!     print *, 'FintSub(2)', FintSub(2)
!     print *, 'FintSub(3)', FintSub(3)
!     if (FintSub(3).ne.zero)then
!        print *, 'ff(1)', ff(1); pause
!     endif
!     print *, 'ff(1)', callcount, ff(1)

    return

  end function xsect_nlo_s

end module mod_xsects_nlo

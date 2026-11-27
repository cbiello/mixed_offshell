module mod_xsects_nnlo_sub_wp_gg
  use mod_types; use common_def; use mod_consts_dp
#if (_DebugFeatures == 1)
  use mod_consts_debug
#endif
  use mod_kinematics_gen; use mod_kinematicsNLOz;  use mod_kinematicsLO
  use mod_process
  use mod_histograms
  use mod_lumi
  use mod_auxlumi_wp
  use mod_amplitudes_tree_wp
  use mod_subtrfn_nnlo_dy_gg
  use mod_auxfunctions; use mod_procfunctions
  use mod_splitting
  implicit none

  private

  public :: xsect_nnlo_sub124_wp_gg,xsect_nnlo_sub12_wp_gg

contains


  function xsect_nnlo_sub124_wp_gg(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer      :: xsect_nnlo_sub124_wp_gg,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15)   :: yRnd(30),ff(1),vegasweight
    type(KinConfig), pointer :: HardProc_z124,Col42Lim_z124
    type(KinConfig), pointer :: HardProc_1z24,Col41Lim_1z24
    real(dp)     :: xx(8),z
    real(dp)     :: tau,ylab,jac,mu,muR,muF,weightarray(nobs)
    real(dp)     :: xi1,xi2,f1(-6:6),f2(-6:6),res(2),res_j(6),respdf
    complex(dp)  :: res_sping_j(6,-1:1,-1:1)
    real(dp)     :: M2,shat,z4,s24,s14
    real(dp)     :: subint_z124_hard,subint_1z24_hard,subint_z124_col42,subint_1z24_col41,Fint_gg_sub_z124(2),Fint_gg_sub_1z24(2)
    logical      :: makecut
    real(dp)     :: pref

    callcount = callcount + 1

    xsect_nnlo_sub124_wp_gg = 0
    ff(1) = zero

        
    HardProc_z124=>MyKinSet%AllKinConfig(1)
    Col42Lim_z124=>MyKinSet%AllKinConfig(2)
    !
    HardProc_1z24=>MyKinSet%AllKinConfig(3)
    Col41Lim_1z24=>MyKinSet%AllKinConfig(4)

    Fint_gg_sub_z124 = zero
    Fint_gg_sub_1z24 = zero

    ! initialize
    xx = zero
    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

    z=buff+onet*real(yRnd(kNLO_max_full),dp)
    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = sh*tau
       

    !-----------------------------------------------------------!
    !   construct the function ONLOqg * FF * FLMqg(z*1,2,4)/z --!
    !-----------------------------------------------------------!

    call kinematics_c42only_nloz(M2,ylab,z,one,xx(kNLO_min:kNLO_max_full),HardProc_z124,Col42Lim_z124)

    !-- regular piece
    xi1 = HardProc_z124%PartFrac(1)
    xi2 = HardProc_z124%PartFrac(2)

    call cut_histo_w(xi1,xi2,5,HardProc_z124%AmpMom,makecut,binarray,weightarray,mu)
    if ( makecut.or. HardProc_z124%flag ) then
       Fint_gg_sub_z124(1) = zero
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree_j(HardProc_z124%AmpMom,res_sping_j,res_j)
       respdf = gg_lumi_qxg(res_j,f1,f2) + gg_lumi_qg(res_j,f1,f2) 

       shat = sh*xi1*xi2
       subint_z124_hard = gg_sub_qgz124(shat,muF**2,z)
       Fint_gg_sub_z124(1) = subint_z124_hard*HardProc_z124%wgt*respdf*jac*pref
    endif
    call HistoFill(0,binarray,weightarray,Fint_gg_sub_z124(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    !-- collinear limit
    xi1 = Col42Lim_z124%PartFrac(1)
    xi2 = Col42Lim_z124%PartFrac(2)


    call cut_histo_w(xi1,xi2,4,Col42Lim_z124%AmpMom,makecut,binarray,weightarray,mu)

    if ( makecut.or. Col42Lim_z124%flag ) then
       Fint_gg_sub_z124(2) = zero
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(Col42Lim_z124%AmpMom,res)
       respdf = gg_lumi_qxq(res,f1,f2) + gg_lumi_qqx(res,f1,f2)

       z4 = Col42Lim_z124%KinInv(1)
       s24 = Col42Lim_z124%KinInv(2)
       shat = sh*xi1*xi2

       subint_z124_col42 = gg_sub_qgz124(shat,muF**2,z)
       Fint_gg_sub_z124(2) = -(two/s24)*Pqg_i(z4)*subint_z124_col42*Col42Lim_z124%wgt*respdf*jac*pref

       call HistoFill(1,binarray,weightarray,Fint_gg_sub_z124(2),vegasweight,iternumber) 

    endif


    !-----------------------------------------------------------!
    !   construct the function ONLOgq * FF * FLMqg(1,z*2,4)/z --!
    !-----------------------------------------------------------!

    call kinematics_C41only_nloz(M2,ylab,one,z,xx(kNLO_min:kNLO_max_full),HardProc_1z24,Col41Lim_1z24)

    !-- regular piece
    xi1 = HardProc_1z24%PartFrac(1)
    xi2 = HardProc_1z24%PartFrac(2)

    call cut_histo_w(xi1,xi2,5,HardProc_1z24%AmpMom,makecut,binarray,weightarray,mu)

    if ( makecut.or. HardProc_1z24%flag ) then
       Fint_gg_sub_1z24(1) = zero
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree_j(HardProc_1z24%AmpMom,res_sping_j,res_j)
       respdf = gg_lumi_gqx(res_j,f1,f2) + gg_lumi_gq(res_j,f1,f2) 

       shat = sh*xi1*xi2
       subint_1z24_hard = gg_sub_qgz124(shat,muF**2,z)
       Fint_gg_sub_1z24(1) = subint_1z24_hard*HardProc_1z24%wgt*respdf*jac*pref
       call HistoFill(3,binarray,weightarray,Fint_gg_sub_1z24(1),vegasweight,iternumber)
    endif


    !-- collinear limit
    xi1 = Col41Lim_1z24%PartFrac(1)
    xi2 = Col41Lim_1z24%PartFrac(2)

    call cut_histo_w(xi1,xi2,4,Col41Lim_1z24%AmpMom,makecut,binarray,weightarray,mu)

    if ( makecut.or. Col41Lim_1z24%flag ) then
       Fint_gg_sub_1z24(2) = zero
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(Col41Lim_1z24%AmpMom,res)
       respdf = gg_lumi_qxq(res,f1,f2) + gg_lumi_qqx(res,f1,f2)

       z4 = Col41Lim_1z24%KinInv(1)
       s14 = Col41Lim_1z24%KinInv(2)
       shat = sh*xi1*xi2

       subint_1z24_col41 = gg_sub_qgz124(shat,muF**2,z)
       Fint_gg_sub_1z24(2) = -(two/s14)*Pqg_i(z4)*subint_1z24_col41*Col41Lim_1z24%wgt*respdf*jac*pref

    endif
    call HistoFill(10,binarray,weightarray,Fint_gg_sub_1z24(2),vegasweight,iternumber) 
          
    ff(1) = sum(Fint_gg_sub_z124(1:2))  + sum(Fint_gg_sub_1z24(1:2))
    

#if (_DebugFeatures == 1)
    Fint_gg_sub_z124_CheckLim = Fint_gg_sub_z124
    Fint_gg_sub_1z24_CheckLim = Fint_gg_sub_1z24
#endif

    return

  end function xsect_nnlo_sub124_wp_gg



  

  function xsect_nnlo_sub12_wp_gg(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    implicit none
    integer          :: xsect_nnlo_sub12_wp_gg,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15)       :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: SubProc
    real(dp)         :: xx(5)
    real(dp)         :: tau,ylab,jac,mu,muR,muF,shat
    real(dp)         :: f1(-6:6),f2(-6:6),xi1,xi2,res(2),respdf
    real(dp)         :: z,zbar
    real(dp)         :: subint,FintSub,my_lumi
    logical          :: makecut
    real(dp)     :: pref
    real(dp)         :: weightarray(nobs)

    callcount = callcount + 1

    xsect_nnlo_sub12_wp_gg = 0
    ff(1) = zero
    xx = zero

    SubProc=>MyKinSet%AllKinConfig(1)


    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    
    z     = buff + onet*real(yRnd(kLO_max_full),dp)
    zbar  = buff + onet*real(yRnd(kLO_max_full+1),dp)

    !Only z1,zbar2 kinematics in this channel
    !Call this function with arguments z,zbar 
    call  get_tauyxi1xi2(xx(1:tauy_max),z,zbar,tau,ylab,xi1,xi2,jac)
    shat = tau*sh/z/zbar

    ! rescale xi1 & xi2 by z & zbar in these next two functions
    call kinematics_lo(z*xi1,zbar*xi2,xx(kLO_min:kLO_max_full),SubProc)
    call cut_histo_w(z*xi1,zbar*xi2,4,SubProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut) return

    call setscale(mu,muF,muR)
    call get_pdf(xi1,xi2,muF,f1,f2)
    call get_prefactor(muR,2,pref)   
    call restree(SubProc%AmpMom,res)
    respdf = gg_lumi_qxq(res,f1,f2) + gg_lumi_qqx(res,f1,f2)

    subint = gg_sub_qqbz1zb2(shat,muF**2,z,zbar)
    FintSub = respdf * SubProc%wgt * jac*pref * subint
       
    call HistoFill(0,binarray,weightarray,FintSub,vegasweight,iternumber)
    call HistoFill(10,binarray,weightarray,zero,vegasweight,iternumber)
       
    ff(1) = FintSub


  end function xsect_nnlo_sub12_wp_gg

end module mod_xsects_nnlo_sub_wp_gg

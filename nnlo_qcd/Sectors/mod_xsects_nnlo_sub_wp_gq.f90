module mod_xsects_nnlo_sub_wp_gq
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
  use mod_amplitudes_loop_wp
  use mod_subtrfn_nnlo_dy_gq
  use mod_subtrfn_nlo_dy
  use mod_auxfunctions; use mod_procfunctions
  use mod_limits_dy_qg
  use mod_splitting

  implicit none

  private

  public :: xsect_nnlo_sub124_wp_gq, xsect_nnlo_sub12_wp_gq, xsect_nnlo_subv_wp_gq

contains

  function xsect_nnlo_sub124_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer      :: xsect_nnlo_sub124_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15)   :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc_GQQ_z124,Col1Lim_GQQ_z124,HardProc_GQQ_1z24,Col1Lim_GQQ_1z24
    type(KinConfig),pointer :: HardProc_GQQ_124,Col1Lim_GQQ_124
    type(KinConfig),pointer :: HardProc_QQBG_z124,SoftLim_QQBG_z124,Col1Lim_QQBG_z124,Col2Lim_QQBG_z124,SCol1Lim_QQBG_z124,SCol2Lim_QQBG_z124
    real(dp) :: xx(8),z
    real(dp) :: tau,ylab,jac,mu,muR,muF,weightarray(nobs)
    real(dp) :: xi1,xi2,f1(-6:6),f2(-6:6),res(2),res_j(6),respdf,respdf_gq,respdf_qqb
    real(dp)   :: x4,s14,s24,z41,z42,eta41,eta42,M2,shat,E4
    complex(dp) :: res_sping_j(6,-1:1,-1:1),hell_mat(-1:1,-1:1),respdf_sping(-1:1,-1:1)
    real(dp)   :: subint_gqq_z124_hard,subint_gqq_z124_col1,subint_gqq_1z24_hard,subint_gqq_1z24_col1, &
         subint_gqq_124_hard, subint_gqq_124_col1,subint_qqbg_1z24_hard,&
         subint_qqbg_z124_hard,subint_qqbg_z124_soft,subint_qqbg_z124_col1,subint_qqbg_z124_col2,&
         subint_qqbg_z124_scol1,subint_qqbg_z124_scol2, subint_gqq_124_hard_pls, subint_gqq_124_col1_pls
    real(dp)   :: Fint_GQQ_sub_z124(2),Fint_GQQ_sub_1z24(2),Fint_GQQ_sub_124(2),Fint_QQBG_sub_z124(6), Fint_QQBG_sub_1z24
    real(dp)   :: Fint_GQQ_z124, Fint_GQQ_1z24, Fint_GQQ_124, Fint_QQBG_z124, Fint_QQBG_1z24
    logical    :: makecut
    real(dp)   :: pref,my_lumi
    integer    :: i,j

    callcount = callcount + 1
        
    xsect_nnlo_sub124_wp_gq = 0
    ff(1) = zero

    ! for processes with FLMgq(z1,2,4)    
    HardProc_GQQ_z124 => MyKinSet%AllKinConfig(1)
    Col1Lim_GQQ_z124  => MyKinSet%AllKinConfig(2)

    ! for processes with FLMgq(1,z2,4)    
    HardProc_GQQ_1z24 => MyKinSet%AllKinConfig(3)
    Col1Lim_GQQ_1z24  => MyKinSet%AllKinConfig(4)

    ! for processes with FLMgq(1,z,4)    
    HardProc_GQQ_124 => MyKinSet%AllKinConfig(5)
    Col1Lim_GQQ_124  => MyKinSet%AllKinConfig(6)

    ! for processes with FLMqqb(z1,2,4)    
    HardProc_QQBG_z124 => MyKinSet%AllKinConfig(7)
    SoftLim_QQBG_z124  => MyKinSet%AllKinConfig(8)
    Col1Lim_QQBG_z124  => MyKinSet%AllKinConfig(9)
    Col2Lim_QQBG_z124  => MyKinSet%AllKinConfig(10)
    SCol1Lim_QQBG_z124 => MyKinSet%AllKinConfig(11)
    SCol2Lim_QQBG_z124 => MyKinSet%AllKinConfig(12)


    Fint_GQQ_z124  = zero
    Fint_GQQ_1z24  = zero
    Fint_GQQ_124   = zero
    Fint_QQBG_z124  = zero
    Fint_QQBG_z124 = zero



        ! initialize
    xx = zero
    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

    z=buff+onet*real(yRnd(kNLO_max_full),dp)
    x4 = xx(xRHO_nlo)
    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = sh*tau

    !-------------------------------------------------------------------------------------------------
    !   construct the function ONLOgq * FF_gq * FLMgq(z*1,2,4)/z 
    !-------------------------------------------------------------------------------------------------

    call kinematics_c41only_nloz(M2,ylab,z,one,xx(kNLO_min:kNLO_max_full),HardProc_GQQ_z124,Col1Lim_GQQ_z124)

    !---- regular piece 
    
    xi1 = HardProc_GQQ_z124%PartFrac(1) 
    xi2 = HardProc_GQQ_z124%PartFrac(2) 
    call cut_histo_w(xi1,xi2,5,HardProc_GQQ_z124%AmpMom,makecut,binarray,weightarray,mu)
    
    if ( makecut.or. HardProc_GQQ_z124%flag ) then
       Fint_GQQ_sub_z124(1) = zero     
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree_j(HardProc_GQQ_z124%AmpMom,res_sping_j,res_j)
       respdf_gq = gq_lumi_nlo(res_j,f1,f2)

       eta41 = x4
       eta42 = one-x4
       shat = sh*xi1*xi2
       subint_gqq_z124_hard = gq_qg_sub_gqz124(shat,muF**2,z,eta41,eta42,eta41)
       Fint_GQQ_sub_z124(1) = subint_gqq_z124_hard * HardProc_GQQ_z124%wgt *respdf_gq *jac*pref    ! pls+reg
    endif

    call HistoFill(0,binarray,weightarray,Fint_GQQ_sub_z124(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not


    !---   col 1 piece  
    xi1 = Col1Lim_GQQ_z124%PartFrac(1) 
    xi2 = Col1Lim_GQQ_z124%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,Col1Lim_GQQ_z124%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.Col1Lim_GQQ_z124%flag) then
       Fint_GQQ_sub_z124(2) = zero 
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(Col1Lim_GQQ_z124%AmpMom,res)
       respdf_gq = gq_lumi_lo(res,f1,f2)
       z41 = Col1Lim_GQQ_z124%KinInv(1)
       s14 = Col1Lim_GQQ_z124%KinInv(2)

       eta41 = zero       
       eta42 = one
       shat = sh*xi1*xi2

       subint_gqq_z124_col1= gq_qg_sub_gqz124(shat,muF**2,z,eta41,eta42,x4)
       Fint_GQQ_sub_z124(2) =  - (subint_gqq_z124_col1) * ( two/s14*Pqg_i(z41)*respdf_gq ) * jac*pref*Col1Lim_GQQ_z124%wgt
       call HistoFill(1,binarray,weightarray,Fint_GQQ_sub_z124(2),vegasweight,iternumber) 
    endif
    Fint_GQQ_z124 = sum(Fint_GQQ_sub_z124(1:2))


    !----------------------------------------------------
    !   construct the function ONLOqqb * FF_qqb * FLMqqb(z*1,2,4)/z
    !------------------------------------------------------
    call kinematics_fullONLO_nloz(M2,ylab,z,one,xx(kNLO_min:kNLO_max_full),HardProc_QQBG_z124,SoftLim_QQBG_z124,Col1Lim_QQBG_z124,Col2Lim_QQBG_z124,SCol1Lim_QQBG_z124,SCol2Lim_QQBG_z124)

!---- regular piece 

    xi1 = HardProc_QQBG_z124%PartFrac(1) 
    xi2 = HardProc_QQBG_z124%PartFrac(2) 

    call cut_histo_w(xi1,xi2,5,HardProc_QQBG_z124%AmpMom,makecut,binarray,weightarray,mu)
    
    if ( makecut.or. HardProc_QQBG_z124%flag ) then
       Fint_QQBG_sub_z124(1) = zero     
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree_j(HardProc_QQBG_z124%AmpMom,res_sping_j,res_j)
       respdf_qqb = gq_lumi_ns(res_j,f1,f2)

       eta41 = x4
       eta42 = one-x4
       shat = sh*xi1*xi2
       
       subint_qqbg_z124_hard = gq_qg_sub_qqbz124(shat,muF**2,z,eta41,eta42,eta41)
       Fint_QQBG_sub_z124(1) = subint_qqbg_z124_hard * HardProc_QQBG_z124%wgt *respdf_qqb *jac*pref
       
       call HistoFill(2,binarray,weightarray,Fint_QQBG_sub_z124(1),vegasweight,iternumber) 
    endif

    !--- soft piece 
    xi1 = SoftLim_QQBG_z124%PartFrac(1)
    xi2 = SoftLim_QQBG_z124%PartFrac(2)

    call cut_histo_w(xi1,xi2,4,SoftLim_QQBG_z124%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.SoftLim_QQBG_z124%flag) then
       Fint_QQBG_sub_z124(2) = zero
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(SoftLim_QQBG_z124%AmpMom,res)
       respdf_qqb = gq_lumi_lo(res,f1,f2)

       eta41 = x4
       eta42 = one-x4
       shat = sh*xi1*xi2

       subint_qqbg_z124_soft = gq_qg_sub_qqbz124(shat,muF**2,z,eta41,eta42,eta41)
       Fint_QQBG_sub_z124(2) = -subint_qqbg_z124_soft * Eik_s4(SoftLim_QQBG_z124%LimMom)*SoftLim_QQBG_z124%wgt*respdf_qqb*jac*pref

       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_z124(2),vegasweight,iternumber)
    endif
    

!---   col 1 piece 
    xi1 = Col1Lim_QQBG_z124%PartFrac(1)
    xi2 = Col1Lim_QQBG_z124%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,Col1Lim_QQBG_z124%AmpMom,makecut,binarray,weightarray,mu)

     if (makecut.or.Col1Lim_QQBG_z124%flag) then
       Fint_QQBG_sub_z124(3) = zero
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(Col1Lim_QQBG_z124%AmpMom,res)
       respdf_qqb = gq_lumi_lo(res,f1,f2)

       z41 = Col1Lim_QQBG_z124%KinInv(1)
       s14 = Col1Lim_QQBG_z124%KinInv(2)

       eta41 = zero
       eta42 = one
       shat = sh*xi1*xi2

       subint_qqbg_z124_col1 = gq_qg_sub_qqbz124(shat,muF**2,z,eta41,eta42,x4)
       Fint_QQBG_sub_z124(3) = -subint_qqbg_z124_col1 * two/s14*Pqq_i(z41)*Col1Lim_QQBG_z124%wgt*respdf_qqb*jac*pref 
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_z124(3),vegasweight,iternumber)
    endif

    
!---   col 2 piece  
    xi1 = Col2Lim_QQBG_z124%PartFrac(1) 
    xi2 = Col2Lim_QQBG_z124%PartFrac(2) 
    call cut_histo_w(xi1,xi2,4,Col2Lim_QQBG_z124%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.Col2Lim_QQBG_z124%flag) then
       Fint_QQBG_sub_z124(4) = zero 
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(Col2Lim_QQBG_z124%AmpMom,res)
       respdf_qqb = gq_lumi_lo(res,f1,f2)

       z42 = Col2Lim_QQBG_z124%KinInv(1)
       s24 = Col2Lim_QQBG_z124%KinInv(2)
       eta41 = one
       eta42 = zero
       shat = sh*xi1*xi2
       
       subint_qqbg_z124_col2 = gq_qg_sub_qqbz124(shat,muF**2,z,eta41,eta42,eta41)
       Fint_QQBG_sub_z124(4) = -subint_qqbg_z124_col2 * two/s24*Pqq_i(z42)*Col2Lim_QQBG_z124%wgt*respdf_qqb*jac*pref

       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_z124(4),vegasweight,iternumber) 
    endif

    !---   soft + col 1 piece
    xi1 = SCol1Lim_QQBG_z124%PartFrac(1)
    xi2 = SCol1Lim_QQBG_z124%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,SCol1Lim_QQBG_z124%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.SCol1Lim_QQBG_z124%flag) then
       Fint_QQBG_sub_z124(5) = zero
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(SCol1Lim_QQBG_z124%AmpMom,res)
       respdf_qqb = gq_lumi_lo(res,f1,f2)

       eta41 = zero
       eta42 = one
       shat = sh*xi1*xi2
       s14 = SCol1Lim_QQBG_z124%KinInv(1)
       
       subint_qqbg_z124_scol1 = gq_qg_sub_qqbz124(shat,muF**2,z,eta41,eta42,x4)
       Fint_QQBG_sub_z124(5) = subint_qqbg_z124_scol1 * one/s14 * four*Cf*SCol1Lim_QQBG_z124%LimMom(1,1)/SCol1Lim_QQBG_z124%LimMom(1,3) *SCol1Lim_QQBG_z124%wgt*respdf_qqb*jac*pref
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_z124(5),vegasweight,iternumber)
    endif

!---   soft + col 2 piece 
    xi1 = SCol2Lim_QQBG_z124%PartFrac(1)
    xi2 = SCol2Lim_QQBG_z124%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,SCol2Lim_QQBG_z124%AmpMom,makecut,binarray,weightarray,mu)
    if (makecut.or. SCol2Lim_QQBG_z124%flag) then
       Fint_QQBG_sub_z124(6) = zero
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(SCol2Lim_QQBG_z124%AmpMom,res)
       respdf_qqb = gq_lumi_lo(res,f1,f2)
       eta41 = one
       eta42 = zero
       shat = sh*xi1*xi2
       s24 = SCol2Lim_QQBG_z124%KinInv(1)

       subint_qqbg_z124_scol2 = gq_qg_sub_qqbz124(shat,muF**2,z,eta41,eta42,eta41)
       Fint_QQBG_sub_z124(6) = subint_qqbg_z124_scol2 * one/s24*four*Cf*SCol2Lim_QQBG_z124%LimMom(1,2)/SCol2Lim_QQBG_z124%LimMom(1,3)*SCOl2Lim_QQBG_z124%wgt*respdf_qqb*jac*pref
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_z124(6),vegasweight,iternumber)
    endif

    Fint_QQBG_z124 = sum(Fint_QQBG_sub_z124(1:6))





    !-------------------------------------------------------------------------------------------------
    !   construct the function ONLOgq * FF_gq * FLMgq(1,z*2,4)/z 
    !-------------------------------------------------------------------------------------------------

    call kinematics_C41only_nloz(M2,ylab,one,z,xx(kNLO_min:kNLO_max_full),HardProc_GQQ_1z24,Col1Lim_GQQ_1z24)

    !---- regular piece 
    xi1 = HardProc_GQQ_1z24%PartFrac(1) 
    xi2 = HardProc_GQQ_1z24%PartFrac(2) 
    call cut_histo_w(xi1,xi2,5,HardProc_GQQ_1z24%AmpMom,makecut,binarray,weightarray,mu)
    
    if ( makecut.or. HardProc_GQQ_1z24%flag ) then
       Fint_GQQ_sub_1z24(1) = zero     
       Fint_QQBG_sub_1z24 = zero
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree_j(HardProc_GQQ_1z24%AmpMom,res_sping_j,res_j)


       respdf_gq = gq_lumi_nlo(res_j,f1,f2)
    
       eta41 = x4
       eta42 = one-x4
       shat = sh*xi1*xi2

       subint_gqq_1z24_hard = gq_qg_sub_gq1z24(shat,muF**2,z,eta41,eta42)
       Fint_GQQ_sub_1z24(1) = subint_gqq_1z24_hard * HardProc_GQQ_1z24%wgt *respdf_gq *jac*pref

       call HistoFill(3,binarray,weightarray,Fint_GQQ_sub_1z24(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    endif


    Fint_QQBG_1z24 = Fint_QQBG_sub_1z24 

    !---   col 1 piece  
    xi1 = Col1Lim_GQQ_1z24%PartFrac(1) 
    xi2 = Col1Lim_GQQ_1z24%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,Col1Lim_GQQ_1z24%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.Col1Lim_GQQ_1z24%flag) then
       Fint_GQQ_sub_1z24(2) = zero 
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(Col1Lim_GQQ_1z24%AmpMom,res)
       respdf_gq = gq_lumi_lo(res,f1,f2)

       z41 = Col1Lim_GQQ_1z24%KinInv(1)
       s14 = Col1Lim_GQQ_1z24%KinInv(2)
       eta41 = zero
       eta42 = one 
       shat = sh*xi1*xi2

       subint_gqq_1z24_col1= gq_qg_sub_gq1z24(shat,muF**2,z,eta41,eta42)

       Fint_GQQ_sub_1z24(2) =  - subint_gqq_1z24_col1 * ( two/s14*Pqg_i(z41)*respdf_gq ) * jac*pref*Col1Lim_GQQ_1z24%wgt
       call HistoFill(1,binarray,weightarray,Fint_GQQ_sub_1z24(2),vegasweight,iternumber) 
    endif
    Fint_GQQ_1z24 = sum(Fint_GQQ_sub_1z24(1:2))



    !----------------------------------------------------
    !   construct the function ONLOgq * FF_gq * FLMgq(1,2,4)+ (pls + del) pieces from above
    !----------------------------------------------------

    call kinematics_c41only_nloz(M2,ylab,one,one,xx(kNLO_min:kNLO_max_full),HardProc_GQQ_124,Col1Lim_GQQ_124)
    
!---- regular piece 

    xi1 = HardProc_GQQ_124%PartFrac(1) 
    xi2 = HardProc_GQQ_124%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,HardProc_GQQ_124%AmpMom,makecut,binarray,weightarray,mu)
    if ( makecut.or. HardProc_GQQ_124%flag ) then
       Fint_GQQ_sub_124(1) = zero     
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree_j(HardProc_GQQ_124%AmpMom,res_sping_j,res_j)
       respdf_gq = gq_lumi_nlo(res_j,f1,f2)

       eta41 = x4
       eta42 = one-x4
       E4 = HardProc_GQQ_124%LimMom(1,3)
       shat = sh*xi1*xi2

       subint_gqq_124_hard = gq_qg_sub_gq124(shat,muF**2,E4,eta41,eta42,eta41)
       subint_gqq_124_hard_pls = gq_qg_sub_gqz124_pls(shat,muF**2,z,eta41,eta42,eta41) + gq_qg_sub_gq1z24_pls(shat,muF**2,z,eta41,eta42)

       Fint_GQQ_sub_124(1) = (subint_gqq_124_hard - subint_gqq_124_hard_pls) &
            * HardProc_GQQ_124%wgt *respdf_gq *jac*pref
       call HistoFill(3,binarray,weightarray,Fint_GQQ_sub_124(1),vegasweight,iternumber) 
    endif

        
!---   col 1 piece  

    xi1 = Col1Lim_GQQ_124%PartFrac(1) 
    xi2 = Col1Lim_GQQ_124%PartFrac(2) 
    call cut_histo_w(xi1,xi2,4,Col1Lim_GQQ_124%AmpMom,makecut,binarray,weightarray,mu)
    if (makecut.or.Col1Lim_GQQ_124%flag) then
       Fint_GQQ_sub_124(2) = zero 
    else
       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(Col1Lim_GQQ_124%AmpMom,res)
       respdf_gq = gq_lumi_lo(res,f1,f2)

       z41 = Col1Lim_GQQ_124%KinInv(1)
       s14 = Col1Lim_GQQ_124%KinInv(2)
       eta42 = one
       eta41 = zero
       shat = sh*xi1*xi2
       E4 = Col1Lim_GQQ_124%LimMom(1,3)

       subint_gqq_124_col1 = gq_qg_sub_gq124(shat,muF**2,E4,eta41,eta42,x4)
       subint_gqq_124_col1_pls = gq_qg_sub_gqz124_pls(shat,muF**2,z,eta41,eta42,x4) + gq_qg_sub_gq1z24_pls(shat,muF**2,z,eta41,eta42)

       Fint_GQQ_sub_124(2) =  - (subint_gqq_124_col1 - subint_gqq_124_col1_pls ) &
            * ( one/s14*two*Pqg_i(z41)*respdf_gq ) * jac*pref*Col1Lim_GQQ_124%wgt
    endif
    call HistoFill(10,binarray,weightarray,Fint_GQQ_sub_124(2),vegasweight,iternumber) 

    Fint_GQQ_124 = sum(Fint_GQQ_sub_124(1:2))

    ff(1) = Fint_GQQ_z124 + Fint_GQQ_1z24 + Fint_GQQ_124 + Fint_QQBG_z124
    

!    print *, "xx",xx
!    print *, "z",z
!!    print *, "",Fint_GQQ_sub_z124(1) + Fint_GQQ_sub_1z24(1) + Fint_GQQ_sub_124(1) + Fint_QQBG_sub_z124(1)
!!     print *, "qqb", Fint_QQBG_sub_z124
!     print *, "gq z124", Fint_GQQ_sub_z124
!     print *, "gq 1z24", Fint_GQQ_sub_1z24
!     print *, "gq 124", Fint_GQQ_sub_124
!    print *, "ff",ff(1)
!    pause

#if (_DebugFeatures == 1)
    Fint_gq_qg_sub_GQQ_z124_CheckLim = Fint_GQQ_sub_z124
    Fint_gq_qg_sub_GQQ_1z24_CheckLim = Fint_GQQ_sub_1z24
    Fint_gq_qg_sub_GQQ_124_CheckLim = Fint_GQQ_sub_124
    Fint_gq_qg_sub_QQBG_z124_CheckLim = Fint_QQBG_sub_z124

#endif


    return
    
  end function xsect_nnlo_sub124_wp_gq

 
 
   function xsect_nnlo_sub12_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
     implicit none
     integer          :: xsect_nnlo_sub12_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs)
     real(dp15)       :: yRnd(30),ff(1),vegasweight
     type(KinConfig),pointer :: SubProc
     real(dp)         :: xx(5)
     real(dp)         :: tau,ylab,jac,mu,muR,muF,muRsq,muFsq,shat
     real(dp)         :: f1(-6:6),f2(-6:6),xi1,xi2,res(2),respdf_gq
     real(dp)         :: z,zbar,yy(1:4)
     real(dp)         :: intsub_zzbar(0:1),intsub_z
     real(dp)         :: FintSub_zzbar,FintSub_z,FintSub
     logical          :: makecut
     real(dp)         :: pref,my_lumi
     integer          :: i,j
#if (_DebugFeatures == 1)
     logical          :: verbose = .false.
#endif
     real(dp)               :: weightarray(nobs)
 
     callcount = callcount + 1
 
     xsect_nnlo_sub12_wp_gq = 0
     ff(1) = zero
     xx = zero 

     SubProc=>MyKinSet%AllKinConfig(1)

     xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
     call random_number(xx(kLO_max_full))

     z     = buff + onet*real(yRnd(kLO_max_full),dp)
     zbar  = buff + onet*real(yRnd(kLO_max_full+1),dp)
     yy(:)  = buff + onet*real(yRnd(kLO_max_full+2:kLO_max_full+5),dp)  !-- for integration of triple collinear subtraction terms, with max
     FintSub_zzbar = zero
     FintSub_z = zero
     
     !use 1,2 to generate shat
     call  get_tauyxi1xi2(xx(1:tauy_max),one,one,tau,ylab,xi1,xi2,jac)
     shat = tau*sh
 
 ! kinematics and amplitudes are the same for [1,zbar2],[z1,zbar]
    call kinematics_lo(xi1,xi2,xx(kLO_min:kLO_max_full),SubProc)
 
     call cut_histo_w(xi1,xi2,4,SubProc%AmpMom,makecut,binarray,weightarray,mu)
     if (makecut) return
 
     call setscale(mu,muF,muR)
     call get_prefactor(muR,2,pref)
     muFsq = muF**2
     muRsq = muR**2
 
     ! get integrated splitting functions
     intsub_zzbar = gq_qg_sub_qqbz1zbar2(shat,muFsq,z,zbar)
     intsub_z  = gq_qg_sub_qqbz12(shat,muFsq,z,yy)

     call restree(SubProc%AmpMom,res)
 
 ! z1,2
     call  get_tauyxi1xi2(xx(1:tauy_max),z,one,tau,ylab,xi1,xi2,jac)
     if (xi1 .lt. one .and. xi2 .lt. one) then
        call get_pdf(xi1,xi2,muF,f1,f2)
        respdf_gq = gq_lumi_lo(res,f1,f2)
        FintSub_z = respdf_gq * SubProc%wgt * jac*pref * (intsub_z + intsub_zzbar(1) )  !  reg[1,z2] + pls
     endif
 
 ! z,zbar
     call  get_tauyxi1xi2(xx(1:tauy_max),z,zbar,tau,ylab,xi1,xi2,jac)    
     if (xi1 .lt. one .and. xi2 .lt. one) then
        call get_pdf(xi1,xi2,muF,f1,f2)
        respdf_gq = gq_lumi_lo(res,f1,f2)
        FintSub_zzbar = respdf_gq * SubProc%wgt * jac*pref * intsub_zzbar(0)  ! reg + pls
     endif
 
     FintSub =  FintSub_zzbar + FintSub_z
 
        
     call HistoFill(0,binarray,weightarray,FintSub,vegasweight,iternumber)
     call HistoFill(10,binarray,weightarray,zero,vegasweight,iternumber)
        
     ff(1) = FintSub
     
     if ( (ff(1) .ne. ff(1))) then
        print *, "xx",xx
        print *, "z",z
        print *, "zbar",zbar
        print *, "yy",yy
        print *, "FintSub_z",FintSub_z
        print *, "FintSub_zzbar",FintSub_zzbar
        print *, "ff",ff(1)
     endif
     if (ff(1) .ne. ff(1)) ff(1) = zero
 
 
   end  function xsect_nnlo_sub12_wp_gq
 
 
 
   function xsect_nnlo_subv_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
     implicit none
     integer          :: xsect_nnlo_subv_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs)
     real(dp15)       :: yRnd(30),ff(1),vegasweight
     type(KinConfig),pointer :: SubProc
     real(dp)         :: xx(5)
     real(dp)         :: tau,ylab,jac,mu,muR,muF,muRsq,muFsq,q2
     real(dp)         :: f1(-6:6),f2(-6:6),xi1,xi2,res0(2),res(2),respdf_gq
     real(dp)         :: z,intsub(-1:1),Fintsub
     logical          :: makecut
     real(dp)         :: pref,my_lumi
     integer          :: i,j
#if (_DebugFeatures == 1)
     logical          :: verbose = .false.
#endif
     real(dp)               :: weightarray(nobs)
     

 
     callcount = callcount + 1
 
     xsect_nnlo_subv_wp_gq = 0
     ff(1) = zero
     xx = zero
     Fintsub=zero

     SubProc=>MyKinSet%AllKinConfig(1)


    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    z = buff + onet*real(yRnd(kLO_max_full),dp)

    call  get_tauyxi1xi2(xx(1:tauy_max),one,one,tau,ylab,xi1,xi2,jac)
    q2 = tau*sh

     call kinematics_lo(xi1,xi2,xx(kLO_min:kLO_max_full),SubProc)
 
     call cut_histo_w(xi1,xi2,4,SubProc%AmpMom,makecut,binarray,weightarray,mu)
     if (makecut) return
 
     call setscale(mu,muF,muR)
     muFsq = muF**2
     muRsq = muR**2

     
     call get_prefactor(muR,2,pref)

     call res_1loop(SubProc%AmpMom,muFsq,res0,res)
     ! get integrated splitting functions
     intsub    = gq_q_sub_qqbz12(q2,muFsq,z)
 
     
     ! z1,2
     call  get_tauyxi1xi2(xx(1:tauy_max),z,one,tau,ylab,xi1,xi2,jac)
     if (xi1 .lt. one .and. xi2 .lt. one) then
        call get_pdf(xi1,xi2,muF,f1,f2)
        respdf_gq = gq_lumi_lo(res,f1,f2)
        FintSub = respdf_gq * SubProc%wgt * jac*pref * intsub(0)
     endif
        
     call HistoFill(0,binarray,weightarray,FintSub,vegasweight,iternumber)
     call HistoFill(10,binarray,weightarray,zero,vegasweight,iternumber)
        
     ff(1) = FintSub
     
     if (ff(1) .gt. 1d15 .or. (ff(1) .ne. ff(1))) then
        print *, "xx",xx
        print *, "z",z
        print *, "ff",ff(1)
     endif
     if (ff(1) .ne. ff(1)) ff(1) = zero
 
 
   end  function xsect_nnlo_subv_wp_gq




! This is in limits_dy_ns and not in limits_dy_gq.
! I have to call limits_dy_gq for the Pggq_s5 limit
! but then I can't call limits_dy_ns because
! the name eik_s5 conflicts (even those these are different functions)
! So I jus put this limit here and only call limits_dy_gq
  function Eik_s4(p)
    real(dp), intent(in) :: p(4,5)
    real(dp) :: Eik_s4
    real(dp) :: p1(4),p2(4),p4(4)

    p1 = p(:,1)
    p2 = p(:,2)
    p4 = p(:,3)

    Eik_s4 = two*Cf*scr(p1,p2)/scr(p1,p4)/scr(p2,p4)

    return 

  end function Eik_s4


end module mod_xsects_nnlo_sub_wp_gq

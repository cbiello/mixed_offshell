module mod_xsects_nnlo_sub_wp_ns
! does all qqb, qbq, qq, qbqb channels!
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
  use mod_subtrfn_nnlo_dy_ns
  use mod_subtrfn_nnlo_dy_s
  use mod_subtrfn_nlo_dy
  use mod_auxfunctions; use mod_procfunctions
  use mod_limits_dy_ns
  use mod_splitting
  implicit none
  real(dp), parameter :: termNOSPING=one
  real(dp), parameter :: termSPING=one


  private

  public :: xsect_nnlo_sub124_wp_ns, xsect_nnlo_sub12_wp_ns, xsect_nnlo_subv_wp_ns

contains
  



  
  function xsect_nnlo_sub124_wp_ns(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    implicit none
    integer    :: xsect_nnlo_sub124_wp_ns,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc_z12,SoftLim_z12,Col1Lim_z12,Col2Lim_z12,SCol1Lim_z12,SCol2Lim_z12,&
         HardProc_1z2,SoftLim_1z2,Col1Lim_1z2,Col2Lim_1z2,SCol1Lim_1z2,SCol2Lim_1z2,&
         HardProc_12,SoftLim_12,Col1Lim_12,Col2Lim_12,SCol1Lim_12,SCol2Lim_12
    real(dp) :: xx(8)
    real(dp)   :: tau,ylab,jac,mu,muR,muF,muRsq,muFsq
    real(dp)   :: f1(-6:6),f2(-6:6),xi1,xi2,res(2),res_j(6),respdf,respdf_qqb,res_spincorr_hard_pdf,respdf_qg,respdf_gq
    real(dp)   :: z,s14,s24,z41,z42,eta41,eta42,M2,shat,x4,E4
    real(dp)   :: subint_qqbg_z124_hard ,subint_qqbg_1z24_hard , subint_qqbg_124_hard 
    real(dp)   :: subint_qqbg_z124_soft ,subint_qqbg_1z24_soft , subint_qqbg_124_soft 
    real(dp)   :: subint_qqbg_z124_col1 ,subint_qqbg_1z24_col1 , subint_qqbg_124_col1 
    real(dp)   :: subint_qqbg_z124_col2 ,subint_qqbg_1z24_col2 , subint_qqbg_124_col2 
    real(dp)   :: subint_qqbg_z124_scol1,subint_qqbg_1z24_scol1, subint_qqbg_124_scol1
    real(dp)   :: subint_qqbg_z124_scol2,subint_qqbg_1z24_scol2, subint_qqbg_124_scol2
    real(dp)   :: subint_gqq_z124_hard, subint_gqq_z124_col1,subint_qgq_1z24_hard, subint_qgq_1z24_col2
    real(dp)   :: Fint_QQBG_sub_z124(6),Fint_QQBG_sub_1z24(6),Fint_QQBG_sub_124(6)
    real(dp)   :: Fint_GQQ_sub_z124(1:2),Fint_GQQ_z124, Fint_QGQ_sub_1z24(1:2),Fint_QGQ_1z24
    real(dp)   :: Fint_QQBG_z124,Fint_QQBG_1z24,Fint_QQBG_124
    complex(dp) :: res_sping_j(6,-1:1,-1:1),hell_mat(-1:1,-1:1),respdf_sping(-1:1,-1:1)
    real(dp) :: ncoll(3),ndir(3)
    logical    :: makecut
    real(dp)   :: pref,my_lumi
    real(dp)   :: weightarray(nobs)
    integer    :: i,j

    callcount = callcount + 1

    xsect_nnlo_sub124_wp_ns = 0
    ff(1) = 0

    ! for processes with FLM[z1,2,4]    
    HardProc_z12=>MyKinSet%AllKinConfig(1)
    SoftLim_z12=>MyKinSet%AllKinConfig(2)
    Col1Lim_z12=>MyKinSet%AllKinConfig(3)
    Col2Lim_z12=>MyKinSet%AllKinConfig(4)
    SCol1Lim_z12=>MyKinSet%AllKinConfig(5)
    SCol2Lim_z12=>MyKinSet%AllKinConfig(6)

    ! for processes with FLM[1,z2,4]    
    HardProc_1z2=>MyKinSet%AllKinConfig(7)
    SoftLim_1z2=>MyKinSet%AllKinConfig(8)
    Col1Lim_1z2=>MyKinSet%AllKinConfig(9)
    Col2Lim_1z2=>MyKinSet%AllKinConfig(10)
    SCol1Lim_1z2=>MyKinSet%AllKinConfig(11)
    SCol2Lim_1z2=>MyKinSet%AllKinConfig(12)

    ! for processes with FLM[1,2,4]    
    HardProc_12=>MyKinSet%AllKinConfig(13)
    SoftLim_12=>MyKinSet%AllKinConfig(14)
    Col1Lim_12=>MyKinSet%AllKinConfig(15)
    Col2Lim_12=>MyKinSet%AllKinConfig(16)
    SCol1Lim_12=>MyKinSet%AllKinConfig(17)
    SCol2Lim_12=>MyKinSet%AllKinConfig(18)

    Fint_QQBG_sub_z124 = zero
    Fint_QQBG_sub_1z24 = zero
    Fint_QQBG_sub_124  = zero
    Fint_GQQ_sub_z124  = zero
    Fint_QGQ_sub_1z24  = zero


    
    ! initialize
    xx = zero
    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

    z=buff+onet*real(yRnd(kNLO_max_full),dp)
    x4 = xx(xRHO_nlo)
    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = sh*tau
    !----------------------------------------------------
    !   construct the function ONLO FF * FLM(z*1,2,4)/z
    !------------------------------------------------------
    call kinematics_fullONLO_nloz(M2,ylab,z,one,xx(kNLO_min:kNLO_max_full),HardProc_z12,SoftLim_z12,Col1Lim_z12,Col2Lim_z12,SCol1Lim_z12,SCol2Lim_z12)
    call kinematics_c41only_nloz(M2,ylab,z,one,xx(kNLO_min:kNLO_max_full),HardProc_z12,Col1Lim_z12)

!---- regular piece 

    xi1 = HardProc_z12%PartFrac(1) 
    xi2 = HardProc_z12%PartFrac(2) 

    call cut_histo_w(xi1,xi2,5,HardProc_z12%AmpMom,makecut,binarray,weightarray,mu)
    if ( makecut.or. HardProc_z12%flag ) then
       Fint_QQBG_sub_z124(1) = zero
       Fint_GQQ_sub_z124(1) = zero
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)     
       call get_prefactor(muR,2,pref)
       call restree_j(HardProc_z12%AmpMom,res_sping_j,res_j)
       respdf_qqb = ns_lumi_nlo(res_j,f1,f2)
       respdf_gq = s_lumi_gq(res_j,f1,f2)*tag_s_s + ns_lumi_gq(res_j,f1,f2)*tag_ns_s
       
       ! ideally these should be in KinInv or something
       eta41 = x4 
       eta42 = one-x4
       shat = sh*xi1*xi2

       subint_qqbg_z124_hard  = qqb_ns_sub_qqbz124(shat,muFsq,z,eta41,eta42,eta41)
       subint_gqq_z124_hard   =   qqb_s_sub_gqz124(shat,muFsq,z,eta41,eta42,eta41) 

       Fint_QQBG_sub_z124(1) =  HardProc_z12%wgt * subint_qqbg_z124_hard*respdf_qqb *jac*pref * termNOSPING
       Fint_GQQ_sub_z124(1)  =  HardProc_z12%wgt * subint_gqq_z124_hard *respdf_gq  *jac*pref * termNOSPING
    endif

    call HistoFill(0,binarray,weightarray,Fint_QQBG_sub_z124(1)+Fint_GQQ_sub_z124(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    
    !--- soft piece    

    xi1 = SoftLim_z12%PartFrac(1) 
    xi2 = SoftLim_z12%PartFrac(2) 

    call cut_histo_w(xi1,xi2,4,SoftLim_z12%AmpMom,makecut,binarray,weightarray,mu)

    if ( makecut.or.SoftLim_z12%flag ) then
       Fint_QQBG_sub_z124(2) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(SoftLim_z12%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)

       eta41 = x4 
       eta42 = one-x4
       shat = sh*xi1*xi2
       subint_qqbg_z124_soft = qqb_ns_sub_qqbz124(shat,muFsq,z,eta41,eta42,eta41)
       Fint_QQBG_sub_z124(2) = -subint_qqbg_z124_soft * Eik_s4(SoftLim_z12%LimMom) * SoftLim_z12%wgt*respdf_qqb*jac*pref  * termNOSPING
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_z124(2),vegasweight,iternumber) 
    endif

    
!---   col 1 piece  

    xi1 = Col1Lim_z12%PartFrac(1) 
    xi2 = Col1Lim_z12%PartFrac(2) 

     call cut_histo_w(xi1,xi2,4,Col1Lim_z12%AmpMom,makecut,binarray,weightarray,mu)

     if (makecut.or. Col1Lim_z12%flag ) then
       Fint_QQBG_sub_z124(3) = zero 
       Fint_GQQ_sub_z124(2) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(Col1Lim_z12%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)
! ??
       respdf_gq = s_lumi_xj(res,f1,f2)*tag_s_s + ns_lumi_lo(res,f1,f2)*tag_ns_s

       z41 = Col1Lim_z12%KinInv(1)
       s14 = Col1Lim_z12%KinInv(2)
       eta41 = zero
       eta42 = one
       shat = sh*xi1*xi2
       subint_qqbg_z124_col1 =  qqb_ns_sub_qqbz124(shat,muFsq,z,eta41,eta42,x4)
       subint_gqq_z124_col1   =   qqb_s_sub_gqz124(shat,muFsq,z,eta41,eta42,x4) 

       Fint_QQBG_sub_z124(3) = -subint_qqbg_z124_col1 * two/s14*Pqq_i(z41) *Col1Lim_z12%wgt *respdf_qqb *jac*pref  * termNOSPING
       Fint_GQQ_sub_z124(2)  = -subint_gqq_z124_col1  * two/s14*Pqg_i(z41) *Col1Lim_z12%wgt *respdf_gq  *jac*pref  * termNOSPING
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_z124(3)+Fint_GQQ_sub_z124(2),vegasweight,iternumber) 
    endif


!---   col 2 piece  

    xi1 = Col2Lim_z12%PartFrac(1) 
    xi2 = Col2Lim_z12%PartFrac(2) 

    call cut_histo_w(xi1,xi2,4,Col2Lim_z12%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.Col2Lim_z12%flag) then
       Fint_QQBG_sub_z124(4) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(Col2Lim_z12%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)

       
       z42 = Col2Lim_z12%KinInv(1)
       s24 = Col2Lim_z12%KinInv(2)
       eta41 = one
       eta42 = zero
       shat = sh*xi1*xi2
       subint_qqbg_z124_col2 = -qqb_ns_sub_qqbz124(shat,muFsq,z,eta41,eta42,eta41)
       Fint_QQBG_sub_z124(4) = subint_qqbg_z124_col2 * two/s24*Pqq_i(z42)*Col2Lim_z12%wgt*respdf_qqb*jac*pref  * termNOSPING
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_z124(4),vegasweight,iternumber) 
    endif


!---   soft + col 1 piece  

    xi1 = SCol1Lim_z12%PartFrac(1) 
    xi2 = SCol1Lim_z12%PartFrac(2) 

    call cut_histo_w(xi1,xi2,4,SCol1Lim_z12%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.SCol1Lim_z12%flag) then
       Fint_QQBG_sub_z124(5) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(SCol1Lim_z12%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)

       eta41 = zero
       eta42 = one
       shat = sh*xi1*xi2
       subint_qqbg_z124_scol1 =  qqb_ns_sub_qqbz124(shat,muFsq,z,eta41,eta42,x4)
       s14 = SCol1Lim_z12%KinInv(1)

       
       Fint_QQBG_sub_z124(5) = subint_qqbg_z124_scol1 * one/s14 * four*Cf*SCol1Lim_z12%LimMom(1,1)/SCol1Lim_z12%LimMom(1,3) &
            *SCol1Lim_z12%wgt*respdf_qqb*jac*pref  * termNOSPING
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_z124(5),vegasweight,iternumber) 
    endif


!---   soft + col 2 piece  

    xi1 = SCol2Lim_z12%PartFrac(1) 
    xi2 = SCol2Lim_z12%PartFrac(2) 

    call cut_histo_w(xi1,xi2,4,SCol2Lim_z12%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or. SCol2Lim_z12%flag) then
       Fint_QQBG_sub_z124(6) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(SCol2Lim_z12%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)
       
       eta41 = one
       eta42 = zero
       shat = sh*xi1*xi2

       s24 = SCol2Lim_z12%KinInv(1)
       subint_qqbg_z124_scol2 = qqb_ns_sub_qqbz124(shat,muFsq,z,eta41,eta42,eta41)

       Fint_QQBG_sub_z124(6) = subint_qqbg_z124_scol2 * one/s24*four*Cf*SCol2Lim_z12%LimMom(1,2)/SCol2Lim_z12%LimMom(1,3)&
            *SCol2Lim_z12%wgt*respdf_qqb*jac*pref * termNOSPING
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_z124(6),vegasweight,iternumber) 
    endif

    Fint_QQBG_z124 = sum(Fint_QQBG_sub_z124(1:6))
    Fint_GQQ_z124  = sum(Fint_GQQ_sub_z124(1:2))

    !----------------------------------------------------
    !   construct the function ONLO FF * FLM(1,z*2,4)/z
    !------------------------------------------------------

    call kinematics_fullONLO_nloz(M2,ylab,one,z,xx(kNLO_min:kNLO_max_full),HardProc_1z2,SoftLim_1z2,Col1Lim_1z2,Col2Lim_1z2,SCol1Lim_1z2,SCol2Lim_1z2)

    !---- regular piece 

   xi1 = HardProc_1z2%PartFrac(1) 
   xi2 = HardProc_1z2%PartFrac(2) 


   call cut_histo_w(xi1,xi2,5,HardProc_1z2%AmpMom,makecut,binarray,weightarray,mu)
    if (makecut.or.HardProc_1z2%flag ) then
       Fint_QQBG_sub_1z24(1) = zero
       Fint_QGQ_sub_1z24(1) = zero
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree_j(HardProc_1z2%AmpMom,res_sping_j,res_j)
       respdf_qqb = ns_lumi_nlo(res_j,f1,f2)
       respdf_qg = s_lumi_qg(res_j,f1,f2)*tag_s_s + ns_lumi_qg(res_j,f1,f2)*tag_ns_s
       
       eta41 = x4 
       eta42 = one-x4
       shat = sh*xi1*xi2

       subint_qqbg_1z24_hard  = qqb_ns_sub_qqbz124(shat,muFsq,z,eta42,eta41,eta42)
       subint_qgq_1z24_hard   =   qqb_s_sub_gqz124(shat,muFsq,z,eta42,eta41,eta42) 


       Fint_QQBG_sub_1z24(1) = subint_qqbg_1z24_hard * HardProc_1z2%wgt *respdf_qqb *jac*pref  * termNOSPING
       Fint_QGQ_sub_1z24(1)  = subint_qgq_1z24_hard  * HardProc_1z2%wgt *respdf_qg *jac*pref  * termNOSPING
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_1z24(1)+Fint_QGQ_sub_1z24(1),vegasweight,iternumber) 
    endif

       
!--- soft piece    

    xi1 = SoftLim_1z2%PartFrac(1) 
    xi2 = SoftLim_1z2%PartFrac(2) 

    call cut_histo_w(xi1,xi2,4,SoftLim_1z2%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.SoftLim_1z2%flag) then
       Fint_QQBG_sub_1z24(2) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(SoftLim_1z2%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)
       eta41 = x4 
       eta42 = one-x4
       shat = sh*xi1*xi2

       subint_qqbg_1z24_soft = qqb_ns_sub_qqbz124(shat,muFsq,z,eta42,eta41,eta42)
       Fint_QQBG_sub_1z24(2) = -subint_qqbg_1z24_soft * Eik_s4(SoftLim_1z2%LimMom)*SoftLim_1z2%wgt*respdf_qqb*jac*pref  * termNOSPING
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_1z24(2),vegasweight,iternumber) 
    endif


    
!---   col 1 piece  

    xi1 = Col1Lim_1z2%PartFrac(1) 
    xi2 = Col1Lim_1z2%PartFrac(2) 

     call cut_histo_w(xi1,xi2,4,Col1Lim_1z2%AmpMom,makecut,binarray,weightarray,mu)

     if (makecut.or.Col1Lim_1z2%flag) then
       Fint_QQBG_sub_1z24(3) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(Col1Lim_1z2%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)
       
       z41 = Col1Lim_1z2%KinInv(1)
       s14 = Col1Lim_1z2%KinInv(2)
       eta41 = zero
       eta42 = one
       shat = sh*xi1*xi2

       subint_qqbg_1z24_col1 =  qqb_ns_sub_qqbz124(shat,muFsq,z,eta42,eta41,eta42)


       Fint_QQBG_sub_1z24(3) = -subint_qqbg_1z24_col1 * two/s14*Pqq_i(z41)*Col1Lim_1z2%wgt*respdf_qqb*jac*pref  * termNOSPING

       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_1z24(3),vegasweight,iternumber) 
    endif

    
!---   col 2 piece  

    xi1 = Col2Lim_1z2%PartFrac(1) 
    xi2 = Col2Lim_1z2%PartFrac(2) 

    call cut_histo_w(xi1,xi2,4,Col2Lim_1z2%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.Col2Lim_1z2%flag) then
       Fint_QQBG_sub_1z24(4) = zero 
       Fint_QGQ_sub_1z24(2) = zero
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(Col2Lim_1z2%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)
       respdf_qg = s_lumi_ix(res,f1,f2)*tag_s_s + ns_lumi_lo(res,f1,f2)*tag_ns_s


       z42 = Col2Lim_1z2%KinInv(1)
       s24 = Col2Lim_1z2%KinInv(2)
       eta41 = one
       eta42 = zero
       shat = sh*xi1*xi2

       subint_qqbg_1z24_col2 = qqb_ns_sub_qqbz124(shat,muFsq,z,eta42,eta41,one-x4)
       subint_qgq_1z24_col2   =  qqb_s_sub_gqz124(shat,muFsq,z,eta42,eta41,one-x4)

       Fint_QQBG_sub_1z24(4) = -subint_qqbg_1z24_col2 * two/s24*Pqq_i(z42) * Col2Lim_1z2%wgt* respdf_qqb *jac*pref  * termNOSPING
       Fint_QGQ_sub_1z24(2) =  -subint_qgq_1z24_col2  * two/s24*Pqg_i(z42) * Col2Lim_1z2%wgt* respdf_qg *jac*pref  * termNOSPING
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_1z24(4)+Fint_QGQ_sub_1z24(2),vegasweight,iternumber) 
    endif

    !---   soft + col 1 piece  

    xi1 = SCol1Lim_1z2%PartFrac(1) 
    xi2 = SCol1Lim_1z2%PartFrac(2) 

    call cut_histo_w(xi1,xi2,4,SCol1Lim_1z2%AmpMom,makecut,binarray,weightarray,mu)

     if (makecut.or.SCol1Lim_1z2%flag) then
       Fint_QQBG_sub_1z24(5) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(SCol1Lim_1z2%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)

       s14 = SCol1Lim_1z2%KinInv(1)
       eta41 = zero
       eta42 = one
       shat = sh*xi1*xi2

       subint_qqbg_1z24_scol1 =  qqb_ns_sub_qqbz124(shat,muFsq,z,eta42,eta41,eta42)


       Fint_QQBG_sub_1z24(5) = subint_qqbg_1z24_scol1 * one/s14 * four*Cf*SCol1Lim_1z2%LimMom(1,1)/SCol1Lim_1z2%LimMom(1,3) &
            *SCol1Lim_1z2%wgt*respdf_qqb*jac*pref  * termNOSPING
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_1z24(5),vegasweight,iternumber) 
    endif


!---   soft + col 2 piece  

    xi1 = SCol2Lim_1z2%PartFrac(1) 
    xi2 = SCol2Lim_1z2%PartFrac(2) 

    call cut_histo_w(xi1,xi2,4,SCol2Lim_1z2%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.SCol2Lim_1z2%flag) then
       Fint_QQBG_sub_1z24(6) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(SCol2Lim_1z2%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)

       s24 = SCol2Lim_1z2%KinInv(1)
       eta41 = one
       eta42 = zero
       shat = sh*xi1*xi2

       subint_qqbg_1z24_scol2 = qqb_ns_sub_qqbz124(shat,muFsq,z,eta42,eta41,one-x4)

       Fint_QQBG_sub_1z24(6) = subint_qqbg_1z24_scol2 * one/s24*four*Cf*SCol2Lim_1z2%LimMom(1,2)/SCol2Lim_1z2%LimMom(1,3)&
            *SCol2Lim_1z2%wgt*respdf_qqb*jac*pref * termNOSPING
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_1z24(6),vegasweight,iternumber) 
    endif

    Fint_QQBG_1z24 = sum(Fint_QQBG_sub_1z24(1:6))
    Fint_QGQ_1z24 = sum(Fint_QGQ_sub_1z24(1:2))

    !----------------------------------------------------
    !   construct the function ONLO FF * FLM(1,2,4)+ pls pieces from above
    !------------------------------------------------------

    
    
    call kinematics_fullONLO_nloz(M2,ylab,one,one,xx(kNLO_min:kNLO_max_full),HardProc_12,SoftLim_12,Col1Lim_12,Col2Lim_12,SCol1Lim_12,SCol2Lim_12)


    !---- regular piece 

   xi1 = HardProc_12%PartFrac(1) 
   xi2 = HardProc_12%PartFrac(2) 

   call cut_histo_w(xi1,xi2,5,HardProc_12%AmpMom,makecut,binarray,weightarray,mu)

   if (makecut.or.HardProc_12%flag ) then
       Fint_QQBG_sub_124(1) = zero
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree_j(HardProc_12%AmpMom,res_sping_j,res_j)
       respdf_qqb = ns_lumi_nlo(res_j,f1,f2)
       respdf_sping = ns_lumi_nlo_sping(res_sping_j,f1,f2)

       eta41 = x4 
       eta42 = one-x4
       shat = sh*xi1*xi2

       E4 = HardProc_12%LimMom(1,3)

       subint_qqbg_124_hard = qqb_ns_sub_qqb124(shat,muFsq,eta41,eta42,eta41,eta42,E4)   ! 124 coeff
       subint_qqbg_124_hard = subint_qqbg_124_hard  - qqb_ns_sub_qqbz124_pls(shat,muFsq,z,eta41,eta42,eta41) - qqb_ns_sub_qqbz124_pls(shat,muFsq,z,eta42,eta41,eta42)   ! -pls(z=1) from [z1,2,4] and [1,z2,4]


       ! spin-correlated part
       ncoll(:) = HardProc_12%LimMom(2:4,4)
       ndir(:) =  HardProc_12%LimMom(2:4,5)
       hell_mat = comp_hell_mat(ncoll,ndir)
       call zipp(respdf_sping,hell_mat,res_spincorr_hard_pdf)

       Fint_QQBG_sub_124(1) = (subint_qqbg_124_hard*respdf_qqb*termNOSPING -(Ca*tag_ns_gg - two*Nf*TR*tag_ns_nf)/three*res_spincorr_hard_pdf*termSPING)*HardProc_12%wgt*jac*pref
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_124(1),vegasweight,iternumber)    
    endif



!--- soft piece    

    xi1 = SoftLim_12%PartFrac(1) 
    xi2 = SoftLim_12%PartFrac(2) 
    call cut_histo_w(xi1,xi2,4,SoftLim_12%AmpMom,makecut,binarray,weightarray,mu)
    if (makecut.or.SoftLim_12%flag) then
       Fint_QQBG_sub_124(2) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(SoftLim_12%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)

       eta41 = x4 
       eta42 = one-x4
       shat = sh*xi1*xi2

       E4 = SoftLim_12%LimMom(1,3)
       subint_qqbg_124_soft = qqb_ns_sub_qqb124(shat,muFsq,eta41,eta42,eta41,eta42,E4)   ! 124 coeff      

       subint_qqbg_124_soft = subint_qqbg_124_soft - qqb_ns_sub_qqbz124_pls(shat,muFsq,z,eta41,eta42,eta41) - qqb_ns_sub_qqbz124_pls(shat,muFsq,z,eta42,eta41,eta42)   ! -pls(z=1) from [z1,2,4] and [1,z2,4]
       Fint_QQBG_sub_124(2) = -subint_qqbg_124_soft * Eik_s4(SoftLim_12%LimMom)*SoftLim_12%wgt*respdf_qqb*jac*pref*termNOSPING
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_124(2),vegasweight,iternumber) 
    endif


    

!---   col 1 piece  

    xi1 = Col1Lim_12%PartFrac(1) 
    xi2 = Col1Lim_12%PartFrac(2) 

    call cut_histo_w(xi1,xi2,4,Col1Lim_12%AmpMom,makecut,binarray,weightarray,mu)

     if (makecut.or.Col1Lim_12%flag) then
       Fint_QQBG_sub_124(3) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(Col1Lim_12%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)

       z41 = Col1Lim_12%KinInv(1)
       s14 = Col1Lim_12%KinInv(2)      
       eta41 = zero
       eta42 = one
       shat = sh*xi1*xi2
       E4 = Col1Lim_12%LimMom(1,3)

       subint_qqbg_124_col1 = qqb_ns_sub_qqb124(shat,muFsq,eta41,eta42,x4,eta42,E4)   ! 124 coeff
       subint_qqbg_124_col1 = subint_qqbg_124_col1  - ( qqb_ns_sub_qqbz124_pls(shat,muFsq,z,eta41,eta42,x4) + qqb_ns_sub_qqbz124_pls(shat,muFsq,z,eta42,eta41,eta42))   ! -pls(z=1) from [z1,2,4] and [1,z2,4]

       Fint_QQBG_sub_124(3) =  - (+subint_qqbg_124_col1*Pqq_i(z41)*termNOSPING - (Ca*tag_ns_gg - two*Nf*TR*tag_ns_nf)/three*rPqq_i(z41)*termSPING ) * two/s14*Col1Lim_12%wgt*respdf_qqb*jac*pref
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_124(3),vegasweight,iternumber) 
    endif


!---   col 2 piece  

    xi1 = Col2Lim_12%PartFrac(1) 
    xi2 = Col2Lim_12%PartFrac(2) 

    call cut_histo_w(xi1,xi2,4,Col2Lim_12%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.Col2Lim_12%flag) then
       Fint_QQBG_sub_124(4) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(Col2Lim_12%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)

       z42 = Col2Lim_12%KinInv(1)
       s24 = Col2Lim_12%KinInv(2)
       eta41 = one
       eta42 = zero
       shat = sh*xi1*xi2
       E4 = Col2Lim_12%LimMom(1,3)

       subint_qqbg_124_col2 = qqb_ns_sub_qqb124(shat,muFsq,eta41,eta42,eta41,one-x4,E4)   ! 124 coeff
       subint_qqbg_124_col2 = subint_qqbg_124_col2  - qqb_ns_sub_qqbz124_pls(shat,muFsq,z,eta41,eta42,eta41) - qqb_ns_sub_qqbz124_pls(shat,muFsq,z,eta42,eta41,one-x4)   ! -pls(z=1) from [z1,2,4] and [1,z2,4]
       Fint_QQBG_sub_124(4) =  - (+subint_qqbg_124_col2*Pqq_i(z42)*termNOSPING - (Ca*tag_ns_gg - two*Nf*TR*tag_ns_nf)/three*rPqq_i(z42)*termSPING)*two/s24*Col2Lim_12%wgt*respdf_qqb*jac*pref
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_124(4),vegasweight,iternumber) 
    endif


    !---   soft + col 1 piece  

    xi1 = SCol1Lim_12%PartFrac(1) 
    xi2 = SCol1Lim_12%PartFrac(2) 

    call cut_histo_w(xi1,xi2,4,SCol1Lim_12%AmpMom,makecut,binarray,weightarray,mu)
    if (makecut.or.SCol1Lim_12%flag) then
       Fint_QQBG_sub_124(5) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(SCol1Lim_12%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)

       eta41 = zero
       eta42 = one
       shat = sh*xi1*xi2

       E4 = SCol1Lim_12%LimMom(1,3)
       subint_qqbg_124_scol1 = qqb_ns_sub_qqb124(shat,muFsq,eta41,eta42,x4,eta42,E4)   ! 124 coeff
       subint_qqbg_124_scol1 = subint_qqbg_124_scol1  - qqb_ns_sub_qqbz124_pls(shat,muFsq,z,eta41,eta42,x4) - qqb_ns_sub_qqbz124_pls(shat,muFsq,z,eta42,eta41,eta42)   ! -pls(z=1) from [z1,2,4] and [1,z2,4]
       s14 = SCol1Lim_12%KinInv(1)
       Fint_QQBG_sub_124(5) = subint_qqbg_124_scol1 * one/s14 * four*Cf*SCol1Lim_12%LimMom(1,1)/SCol1Lim_12%LimMom(1,3) *SCol1Lim_12%wgt*respdf_qqb*jac*pref* termNOSPING
       call HistoFill(1,binarray,weightarray,Fint_QQBG_sub_124(5),vegasweight,iternumber) 
    endif


!---   soft + col 2 piece 

    xi1 = SCol2Lim_12%PartFrac(1) 
    xi2 = SCol2Lim_12%PartFrac(2) 

    call cut_histo_w(xi1,xi2,4,SCol2Lim_12%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or. SCol2Lim_12%flag) then
       Fint_QQBG_sub_124(6) = zero 
    else
       call setscale(mu,muF,muR)
       muFsq = muF**2
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree(SCol2Lim_12%AmpMom,res)
       respdf_qqb = ns_lumi_lo(res,f1,f2)

       eta41 = one
       eta42 = zero
       shat = sh*xi1*xi2

       E4 = SCol2Lim_12%LimMom(1,3)
       subint_qqbg_124_scol2 = qqb_ns_sub_qqb124(shat,muFsq,eta41,eta42,eta41,one-x4,E4)   ! 124 coeff
       subint_qqbg_124_scol2 = subint_qqbg_124_scol2  - qqb_ns_sub_qqbz124_pls(shat,muFsq,z,eta41,eta42,eta41) - qqb_ns_sub_qqbz124_pls(shat,muFsq,z,eta42,eta41,one-x4)   ! -pls(z=1) from [z1,2,4] and [1,z2,4]

       s24 = SCol2Lim_12%KinInv(1)
       Fint_QQBG_sub_124(6) = subint_qqbg_124_scol2 * one/s24*four*Cf*SCol2Lim_12%LimMom(1,2)/SCol2Lim_12%LimMom(1,3)*SCol2Lim_12%wgt*respdf_qqb*jac*pref* termNOSPING
    endif

    call HistoFill(10,binarray,weightarray,Fint_QQBG_sub_124(6),vegasweight,iternumber)  !-- close histogram, no matter what

    Fint_QQBG_124 = sum(Fint_QQBG_sub_124(1:6))

 

    ff(1) = Fint_QQBG_z124 + Fint_QQBG_1z24 + Fint_QQBG_124 + Fint_GQQ_z124 + Fint_QGQ_1z24

        


#if (_DebugFeatures == 1)
    Fint_qqb_ns_sub_QQBG_z124_CheckLim = Fint_QQBG_sub_z124
    Fint_qqb_ns_sub_QQBG_1z24_CheckLim = Fint_QQBG_sub_1z24
    Fint_qqb_ns_sub_QQBG_124_CheckLim = Fint_QQBG_sub_124
    Fint_qqb_ns_sub_GQQ_z124_CheckLim = Fint_GQQ_sub_z124
    Fint_qqb_ns_sub_QGQ_1z24_CheckLim = Fint_QGQ_sub_1z24
#endif

    return

  end function xsect_nnlo_sub124_wp_ns

    

  function xsect_nnlo_sub12_wp_ns(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    implicit none
    integer    :: xsect_nnlo_sub12_wp_ns,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: SubProc
    real(dp)   :: xx(5)
    real(dp)   :: tau,ylab,jac,mu,muF,muR,muFsq,muRsq
    real(dp)   :: f1(-6:6),f2(-6:6),xi1,xi2,res(2),respdf_s,respdf_ns,respdf_v,shat
    real(dp)   :: z,zbar,yy(1:4)
    real(dp)   :: intsub_ns_noz,intsub_ns_z(0:1),intsub_ns_zbar(0:1),intsub_ns_zzbar(1:4)
    real(dp)   :: intsub_s_z,intsub_v_z,intsub_s_zbar,intsub_v_zbar
    real(dp)   :: FintSub_zzbar,FintSub_z,FintSub_zbar,FintSub_12,FintSub
    logical    :: makecut
    real(dp)   :: pref
    real(dp)   :: weightarray(nobs)

    callcount = callcount + 1

    xsect_nnlo_sub12_wp_ns = 0
    ff(1) = zero
    xx=zero

    SubProc=>MyKinSet%AllKinConfig(1)

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))

    z     = buff + onet*real(yRnd(kLO_max_full),dp)
    zbar  = buff + onet*real(yRnd(kLO_max_full+1),dp)
    yy(:)  = buff + onet*real(yRnd(kLO_max_full+2:kLO_max_full+5),dp)  !-- for integration of triple collinear subtraction terms, with max

    FintSub_12 = zero
    FintSub_z = zero
    FintSub_zbar = zero
    FintSub_zzbar = zero

    call  get_tauyxi1xi2(xx(1:tauy_max),one,one,tau,ylab,xi1,xi2,jac)
    shat = tau*sh   

! kinematics and amplitudes are the same for [1,2], [z1,2], [1,zbar2],[z1,zbar]
    call kinematics_lo(xi1,xi2,xx(kLO_min:kLO_max_full),SubProc)

    call cut_histo_w(xi1,xi2,4,SubProc%AmpMom,makecut,binarray,weightarray,mu)
    if (makecut) return

    call setscale(mu,muF,muR)
    call get_prefactor(muR,2,pref)
    muFsq = muF**2
    muRsq = muR**2

    ! get integrated splitting functions
    intsub_ns_noz   = qqb_ns_sub_qqb12(shat,muFsq)
    intsub_ns_z     = qqb_ns_sub_qqbz12(shat,muFsq,z,yy) 
    intsub_ns_zbar  = qqb_ns_sub_qqbz12(shat,muFsq,zbar,yy) 
    intsub_ns_zzbar = qqb_ns_sub_qqbz1zbar2(shat,muFsq,z,zbar)
    
    intsub_s_z    = qqb_s_sub_qqbz12(shat,muFsq,z,yy)
    intsub_s_zbar = qqb_s_sub_qqbz12(shat,muFsq,zbar,yy)

    intsub_v_z    = qqb_v_sub_qqbz12(shat,muFsq,z,yy)
    intsub_v_zbar = qqb_v_sub_qqbz12(shat,muFsq,zbar,yy)

    !print *, 'z',z,zbar
    !print *, intsub_ns_z(0),intsub_ns_zbar(0) !-- reg+pls
    !print *, intsub_ns_z(1),intsub_ns_zbar(1) !-- pls
    !print *, intsub_ns_z(0)*two,intsub_ns_z(1)*two
    !print *, intsub_s_z*tag_s_s,intsub_s_zbar*tag_s_s
    !print *, intsub_v_z,intsub_v_zbar
    !pause

    call restree(SubProc%AmpMom,res)

! 12
    if (xi1 .lt. one .and. xi2 .lt. one) then
       call get_pdf(xi1,xi2,muF,f1,f2)
       respdf_ns = ns_lumi_lo(res,f1,f2)
       FintSub_12    = respdf_ns * SubProc%Wgt * jac * pref* (intsub_ns_noz - intsub_ns_z(1) - intsub_ns_zbar(1) + intsub_ns_zzbar(4) )
    endif

! z1,2
    call  get_tauyxi1xi2(xx(1:tauy_max),z,one,tau,ylab,xi1,xi2,jac)
    if (xi1 .lt. one .and. xi2 .lt. one) then
       call get_pdf(xi1,xi2,muF,f1,f2)
       respdf_ns = ns_lumi_lo(res,f1,f2)
       respdf_s = s_lumi_xj(res,f1,f2)
       respdf_v = s_lumi_xj_id(res,f1,f2)
       FintSub_z = SubProc%Wgt * jac * pref* (respdf_ns *(intsub_ns_z(0) + intsub_ns_zzbar(2) + intsub_s_z * tag_ns_s) + &            ! NS, incl. contamination of singlet
                                              respdf_s  *intsub_s_z*tag_s_s + respdf_v * intsub_v_z*tag_s_v   )                            ! S
    endif

! 1,z2
    call  get_tauyxi1xi2(xx(1:tauy_max),one,zbar,tau,ylab,xi1,xi2,jac)
    if (xi1 .lt. one .and. xi2 .lt. one) then
       call get_pdf(xi1,xi2,muF,f1,f2)
       respdf_ns = ns_lumi_lo(res,f1,f2)
       respdf_s = s_lumi_ix(res,f1,f2)
       respdf_v = s_lumi_ix_id(res,f1,f2)
       FintSub_zbar  = SubProc%Wgt * jac * pref*(respdf_ns *  (intsub_ns_zbar(0) + intsub_ns_zzbar(3)+ intsub_s_zbar*tag_ns_s) + &          ! NS, incl. contamination of singlet
                                                 respdf_s  * intsub_s_zbar*tag_s_s + respdf_v * intsub_v_zbar*tag_s_v  )                           ! S
            endif

! zbar1,z2
    call  get_tauyxi1xi2(xx(1:tauy_max),z,zbar,tau,ylab,xi1,xi2,jac)
    if (xi1 .lt. one .and. xi2 .lt. one) then
       call get_pdf(xi1,xi2,muF,f1,f2)
       respdf_ns = ns_lumi_lo(res,f1,f2)
       FintSub_zzbar = respdf_ns * SubProc%Wgt * jac * pref* intsub_ns_zzbar(1)
    endif

    FintSub = FintSub_12 + FintSub_z + FintSub_zbar + FintSub_zzbar
    
    call HistoFill(0,binarray,weightarray,FintSub,vegasweight,iternumber)
    call HistoFill(10,binarray,weightarray,zero,vegasweight,iternumber)
       
    ff(1) = FintSub
    
!    print *, "ff",ff(1)
!    stop

    if (ff(1) .ne. ff(1)) then
       print *, "xx",xx
       print *, "z",z
       print *, "zbar",zbar
       print *, "yy",yy
       print *, "sub 12",FintSub_12
       print *, "sub 12",FintSub_z
       print *, "sub 12",FintSub_zbar
       print *, "sub 12",FintSub_zzbar
    endif
    if (ff(1) .ne. ff(1)) ff(1) = zero

  end function xsect_nnlo_sub12_wp_ns





  function xsect_nnlo_subv_wp_ns(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    implicit none
    integer    :: xsect_nnlo_subv_wp_ns,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: SubProc
    real(dp)   :: xx(5)
    real(dp)   :: tau,ylab,jac,mu,muR,muF,q2,muRsq,muFsq
    real(dp)   :: f1(-6:6),f2(-6:6),xi1,xi2,res(2),res0(2),respdf
    real(dp)   :: z
    real(dp)   :: intsub(-1:1)
    real(dp)   :: FintSub_z12,FintSub_1z2,FintSub_12, FintSub
    logical    :: makecut
    real(dp)   :: pref, my_lumi
    real(dp)   :: weightarray(nobs)
    integer    :: i,j


    callcount = callcount + 1

    xsect_nnlo_subv_wp_ns = 0
    ff(1) = zero
    xx = zero

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

    call res_1loop(SubProc%AmpMom,muRsq,res0,res)

    ! get integrated splitting function
    intsub = qqb_g_sub_qqbz12(q2,muFsq,z)

    ! z1,2
    call get_tauyxi1xi2(xx(1:tauy_max),z,one,tau,ylab,xi1,xi2,jac) 
    if (xi1 .lt. one .and. xi2 .lt. one) then
       call get_pdf(xi1,xi2,muF,f1,f2)
       respdf = ns_lumi_lo(res,f1,f2)
       FintSub_z12   = respdf * SubProc%Wgt * jac * pref* intsub(0)        ! reg+pls
    else 
       FintSub_z12 = zero
    endif


    ! 1,z2
    call get_tauyxi1xi2(xx(1:tauy_max),one,z,tau,ylab,xi1,xi2,jac) 
    if (xi1 .lt. one .and. xi2 .lt. one) then
       call get_pdf(xi1,xi2,muF,f1,f2)
       respdf = ns_lumi_lo(res,f1,f2)
       FintSub_1z2   = respdf * SubProc%Wgt * jac * pref* intsub(0)        ! reg+pls
    else 
       FintSub_1z2 = zero
    endif


    ! 1,2
    call get_tauyxi1xi2(xx(1:tauy_max),one,one,tau,ylab,xi1,xi2,jac) 
    if (xi1 .lt. one .and. xi2 .lt. one) then
       call get_pdf(xi1,xi2,muF,f1,f2)
       respdf = ns_lumi_lo(res,f1,f2)
       FintSub_12   = respdf * SubProc%Wgt * jac * pref* (intsub(-1)-intsub(+1))*two        ! del-pls, two legs
    else 
       FintSub_12 = zero
    endif

    FintSub = FintSub_12 + FintSub_z12 + FintSub_1z2

    call HistoFill(0,binarray,weightarray,FintSub,vegasweight,iternumber)
    call HistoFill(10,binarray,weightarray,zero,vegasweight,iternumber)

       
    ff(1) = FintSub  * tag_ns_gg

!    print *, "ff",ff(1)
!    stop


  end function xsect_nnlo_subv_wp_ns


! this should go in the splitting routine
  function rPqq_i(z)
    real(dp) :: z,rPqq_i

    rPqq_i = Cf * half * (one+z)**2/(one-z)
    rPqq_i = rPqq_i/z

    return

  end function rPqq_i


    
end module mod_xsects_nnlo_sub_wp_ns

    
    
    

module mod_xsects_nnlo_s_aq_raoul
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_kinematics_lo
  use mod_kinematics_lo_z_raoul
  use mod_kinematics_nlo_z_raoul
!  use mod_kinematics_nlo
  use mod_process
  use mod_auxfunctions
  use mod_aux_sectors
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_amplitudes_loop_ppll
  use mod_eikonals
  use mod_partitions
  use mod_splittings_bare
  use mod_subtrfn_nnlo_z_aq
  use mod_subtrfn_nlo_z

  implicit none
  integer, parameter :: nFint  = 18 !-- 6 if full, 3 if use S(C1+C2)-S = 0
  integer, parameter :: nkin   = 18 !-- 4 if full, 3 if use S(C1+C2)-S = 0
  logical, parameter :: fncheck = .false.     !-- print out values for int. sub. functions and their arguments as a check    

  real(dp), parameter :: charges_ns(4,4) = reshape(& 
       [-Qdn,   Qdn,  -Qup,   Qup,   &
         Qdn,  -Qdn,   Qup,  -Qup,   &
         Q_lep, Q_lep, Q_lep, Q_lep, &
        -Q_lep,-Q_lep,-Q_lep,-Q_lep]  &
        , [4,4])


#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_onloqcd_aq(nFint)
  real(dp), public, save :: FintNNLO_onloewk_aq(nFint)
#endif

  private

  public :: xsect_nnlo_sub12_aq_raoul
  public :: xsect_nnlo_subvqcd_aq_raoul
  public :: xsect_nnlo_onloqcd_aq_raoul
  public :: xsect_nnlo_onloewk_aq_raoul


contains


  !--------------------------------------------------
  !-- subtraction counterterms
  !--------------------------------------------------


  function xsect_nnlo_onloqcd_aq_raoul(yRnd,ff,vegasweight)
    ! The coefficients of ONLO^g FLM[z1_q,2_q,3,4,|5_g]
    integer :: xsect_nnlo_onloqcd_aq_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc_z12,SLim_z12,C1Lim_z12,C2Lim_z12,SC1Lim_z12,SC2Lim_z12
    !--                                                                                                                                                                                                      
    real(dp)    :: xx(kNLO_max_full),z
    real(dp)    :: respdf(ipdf),kin(6)
    real(dp)    :: res_lo(2,2),res_lo_tmp(2,2),res_nlo(2,2),res_nlo_tmp(2,2)
    real(dp)    :: Ec, e5, z5i, s5i,eik,eta51,eta52,eta51lim
    real(dp)    :: subint_z124_hard_u,subint_z124_hard_d,subint_z124_slim_u,subint_z124_slim_d
    real(dp)    :: subint_z124_c1lim_u,subint_z124_c1lim_d,subint_z124_c2lim_u,subint_z124_c2lim_d
    real(dp)    :: subint_z124_sc1lim_u,subint_z124_sc1lim_d,subint_z124_sc2lim_u,subint_z124_sc2lim_d

    xsect_nnlo_onloqcd_aq_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z=buff+onet*real(yRnd(kNLO_max_full),dp)

!!#if (_withchecks == 1)
!!    if (override) then
!!       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
!!       print *, 'overriding input'
!!    endif
!!#endif

    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_r) .or. ((one-z) .lt. buff_z) ) then
       failed_points = failed_points + 1
       return
    endif

!    xx(xE)   = 1E-8_dp
!    xx(xrho) = 1E-8_dp
!    z = one - 1E-8_dp

    call open_histo()

    !! ----------------------------- !!
    !! FLM[z1,2,3,4 | 5g]            !!
    !! ----------------------------- !!

    call kinematics_nlo_z_is(yr=xx,z1=z,z2=one,HardProc=HardProc_z12,C1Lim=C1Lim_z12,C2Lim=C2Lim_z12,SLim=SLim_z12,SC1Lim=SC1Lim_z12,SC2Lim=SC2Lim_z12,compute_etas=.true.)
    HardProc_z12%ids(1:5) = [0,0,id_el,-id_el,id_g]
    C1Lim_z12%ids(1:4)    = [0,0,id_el,-id_el]
    SLim_z12%ids(1:4)     = [0,0,id_el,-id_el]
    C2Lim_z12%ids(1:4)    = [0,0,id_el,-id_el]
    SC1Lim_z12%ids(1:4)   = [0,0,id_el,-id_el]
    sC2Lim_z12%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard

    call cut_histo(HardProc_z12)
    if (HardProc_z12%makecut.or.HardProc_z12%flag) then

       kin(1) = zero
       FintNNLO_onloqcd_aq(1) = zero

    else

       call res_tree_g_qqb(HardProc_z12%AmpMom,res_nlo)
       Ec = HardProc_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = HardProc_z12%Lim_etaij(1,5)
       eta51lim = eta51

       subint_z124_hard_u = aqb_qbg_sub_qqbg_zi(Ec,HardProc_z12%muf(1),z,eta51,eta51lim,Qup)
       subint_z124_hard_d = aqb_qbg_sub_qqbg_zi(Ec,HardProc_z12%muf(1),z,eta51,eta51lim,Qdn)

       if (fncheck) then
          print *, "Ec->",Ec,",mu->",Hardproc_z12%muf(1),",z->",z,",eta51->",eta51
          print *, "z124 u",subint_z124_hard_u
          print *, "z124 d",subint_z124_hard_d
          pause
       endif

       res_nlo_tmp(:,1) = res_nlo(:,1) * subint_z124_hard_d
       res_nlo_tmp(:,2) = res_nlo(:,2) * subint_z124_hard_u
       call get_respdf(aq_lumi,1,1,HardProc_z12,res_nlo_tmp,respdf)

       respdf = respdf*HardProc_z12%wgt

       kin(1) = respdf(1)
       FintNNLO_onloqcd_aq(1) = kin(1)
       call fill_histo(respdf,vegasweight)

    endif

    !-- C1

    call cut_histo(C1Lim_z12)
    if (C1Lim_z12%makecut.or.C1Lim_z12%flag) then

       kin(2) = zero
       FintNNLO_onloqcd_aq(2) = zero

    else

       call res_tree_qqb(C1Lim_z12%AmpMom,res_lo)

       Ec = C1Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = zero
       eta51lim = C1Lim_z12%Lim_etaij(1,5)

       subint_z124_c1lim_u = aqb_qbg_sub_qqbg_zi(Ec,C1Lim_z12%muf(1),z,eta51,eta51lim,Qup)
       subint_z124_c1lim_d = aqb_qbg_sub_qqbg_zi(Ec,C1Lim_z12%muf(1),z,eta51,eta51lim,Qdn)

       res_lo_tmp(:,1) = res_lo(:,1) * subint_z124_c1lim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_z124_c1lim_u


       call get_respdf(aq_lumi,1,1,C1Lim_z12,res_lo_tmp,respdf)

       z5i = C1Lim_z12%Lim_KinInv(1)
       s5i = C1Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z5i)/(one-z5i))&
            * C1Lim_z12%wgt

       FintNNLO_onloqcd_aq(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C2

    call cut_histo(C2Lim_z12)
    if (C2Lim_z12%makecut.or.C2Lim_z12%flag) then

       kin(3) = zero
       FintNNLO_onloqcd_aq(3) = zero

    else
       call res_tree_qqb(C2Lim_z12%AmpMom,res_lo)

       Ec = C2Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = one
       eta51lim = eta51

       subint_z124_c2lim_u = aqb_qbg_sub_qqbg_zi(Ec,C2Lim_z12%muf(1),z,eta51,eta51lim,Qup)
       subint_z124_c2lim_d = aqb_qbg_sub_qqbg_zi(Ec,C2Lim_z12%muf(1),z,eta51,eta51lim,Qdn)

       res_lo_tmp(:,1) = res_lo(:,1) * subint_z124_c2lim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_z124_c2lim_u


       call get_respdf(aq_lumi,1,1,C2Lim_z12,res_lo_tmp,respdf)

       z5i = C2Lim_z12%Lim_KinInv(1)
       s5i = C2Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z5i)/(one-z5i))&
            * C2Lim_z12%wgt

       FintNNLO_onloqcd_aq(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- S

    call cut_histo(SLim_z12)
    if (SLim_z12%makecut.or.SLim_z12%flag) then

       kin(4) = zero
       FintNNLO_onloqcd_aq(4) = zero

    else

       call res_tree_qqb(SLim_z12%AmpMom,res_lo)
       call get_qcd_eik(Cf,SLim_z12%Lim_etaij(1:5,1:5),[1,2],5,eik)

       Ec = SLim_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       
       eta51 = SLim_z12%Lim_etaij(1,5)
       eta51lim = eta51

       subint_z124_slim_u = aqb_qbg_sub_qqbg_zi(Ec,SLim_z12%muf(1),z,eta51,eta51lim,Qup)
       subint_z124_slim_d = aqb_qbg_sub_qqbg_zi(Ec,SLim_z12%muf(1),z,eta51,eta51lim,Qdn)

       res_lo_tmp(:,1) = res_lo(:,1) * subint_z124_slim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_z124_slim_u
       call get_respdf(aq_lumi,1,1,SLim_z12,res_lo_tmp,respdf)

       e5 = SLim_z12%Lim_Ei(5)

       respdf = -respdf*SLim_z12%wgt*eik/e5**2

       kin(4) = respdf(1)
       FintNNLO_onloqcd_aq(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif



    !-- SC1

    call cut_histo(SC1Lim_z12)
    if (SC1Lim_z12%makecut.or.SC1Lim_z12%flag) then

       kin(5) = zero
       FintNNLO_onloqcd_aq(5) = zero

    else

       call res_tree_qqb(SC1Lim_z12%AmpMom,res_lo)

       Ec = SC1Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = zero
       eta51lim = SC1Lim_z12%Lim_etaij(1,5)

       subint_z124_sc1lim_u = aqb_qbg_sub_qqbg_zi(Ec,SC1lim_z12%muf(1),z,eta51,eta51lim,Qup)
       subint_z124_sc1lim_d = aqb_qbg_sub_qqbg_zi(Ec,SC1lim_z12%muf(1),z,eta51,eta51lim,Qdn)

       res_lo_tmp(:,1) = res_lo(:,1) * subint_z124_sc1lim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_z124_sc1lim_u


       call get_respdf(aq_lumi,1,1,SC1Lim_z12,res_lo_tmp,respdf)

       e5 = SC1Lim_z12%Lim_KinInv(1)
       eta51 = SC1Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = respdf*Cf/e5**2/eta51* SC1Lim_z12%wgt

       FintNNLO_onloqcd_aq(5) = respdf(1)
       kin(5) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- SC2

    call cut_histo(SC2Lim_z12)
    if (SC2Lim_z12%makecut.or.SC2Lim_z12%flag) then

       kin(6) = zero
       FintNNLO_onloqcd_aq(6) = zero

    else

       call res_tree_qqb(SC2Lim_z12%AmpMom,res_lo)

       Ec = SC2Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = one
       eta51lim = eta51

       subint_z124_sc2lim_u = aqb_qbg_sub_qqbg_zi(Ec,SC2lim_z12%muf(1),z,eta51,eta51lim,Qup)
       subint_z124_sc2lim_d = aqb_qbg_sub_qqbg_zi(Ec,SC2lim_z12%muf(1),z,eta51,eta51lim,Qdn)


       res_lo_tmp(:,1) = res_lo(:,1) * subint_z124_sc2lim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_z124_sc2lim_u


       call get_respdf(aq_lumi,1,1,SC2Lim_z12,res_lo_tmp,respdf)

       e5 = SC2Lim_z12%Lim_KinInv(1)
       eta52 = SC2Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]

       respdf = respdf*Cf/e5**2/eta52* SC2Lim_z12%wgt

       FintNNLO_onloqcd_aq(6) = respdf(1)
       kin(6) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif


    if (is_nan(kin)) then
       print *, "xx",xx
       print *, "z",z
       print *, "kin",kin(1:6)
       stop
    endif


    ff(1) = sum(kin(1:6))

!!    if (kin(1)*kin(2)*kin(3)*kin(4)*kin(5)*kin(6) .ne. zero) then
!!
!!       print *, "kin in onlo z12", kin(1:6)
!!       print *, "ff",ff(1)
!!       pause
!!   endif


    call close_histo()


  end function xsect_nnlo_onloqcd_aq_raoul



 function xsect_nnlo_onloewk_aq_raoul(yRnd,ff,vegasweight)
    ! The coefficients of ONLO^gamma FLM[1_gamma,z2_q,3,4,|5_q]
    integer :: xsect_nnlo_onloewk_aq_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc_1z2,C1Lim_1z2,C2Lim_1z2,HardProc_12,C1Lim_12,C2Lim_12
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: kin(6)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(2,2),res_lo(2,2),res_tmp(2,2),res_loAA
    real(dp)    :: z,z5i,s5i,e5,eta51,eta52,Ec,eta51lim,eta52lim
    real(dp)    :: subint_1z24_hard,subint_1z24_c1lim,subint_1z24_c2lim
    real(dp)    :: subint_124_hard,subint_124_c1lim,subint_124_c2lim

    xsect_nnlo_onloewk_aq_raoul = 0

    ff(1) = zero
    kin = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z=buff+onet*real(yRnd(kNLO_max_full),dp)

 
#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if ( (xx(xRHO).lt.buff_r) .or. ((one-xx(xRHO)).lt.buff_r) .or. ((one-z) .lt. buff_z) ) then
       failed_points = failed_points + 1
       return
    endif

!    xx(xrho) = 1E-10_dp
!    z = one - 1E-8_dp



    call open_histo()

    !!!!------------------------!!!!!
    !!!! FLM[1_a,z.2_q,3,4|5_q] !!!!!
    !!!!------------------------!!!!!

    
    call kinematics_nlo_z_is(yr=xx,z1=one,z2=z,HardProc=HardProc_1z2,&
         C1Lim=C1Lim_1z2,C2Lim=C2Lim_1z2,compute_etas=.true.)
    HardProc_1z2%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim_1z2%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim_1z2%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard                                                                                                                                                                                                 
    call cut_histo(HardProc_1z2)
    if (HardProc_1z2%makecut.or.HardProc_1z2%flag) then

       kin(1) = zero
       FintNNLO_onloewk_aq(1) = zero
 else

       call res_tree_a_aq(HardProc_1z2%AmpMom,res_nlo)
       call get_respdf(aq_lumi,1,1,HardProc_1z2,res_nlo,respdf)

       Ec = HardProc_1z2%Lim_Ei(1)                      ! Ec = E1 = E2

       eta52 = HardProc_1z2%Lim_etaij(2,5)
       eta52lim = eta52

       subint_1z24_hard = aqb_qbg_sub_aqbq_zi(Ec,HardProc_1z2%muf(1),z,eta52,eta52lim)       
       if (fncheck) then
          print *, "Ec->",Ec,",mu->",Hardproc_1z2%muf(1),",z->",z,",eta52->",eta52
          print *, "1z24 u",subint_1z24_hard
          pause
       endif


       respdf = respdf*HardProc_1z2%wgt * subint_1z24_hard


       kin(1) = respdf(1)
       FintNNLO_onloewk_aq(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1                                                                                                                                                                                                   
    call cut_histo(C1Lim_1z2)

    if (C1Lim_1z2%makecut.or.C1Lim_1z2%flag) then

       kin(2) = zero
       FintNNLO_onloewk_aq(2) = zero

    else

       call res_tree_qqb(C1Lim_1z2%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       call get_respdf(aq_lumi,1,1,C1Lim_1z2,res_lo,respdf)

       z5i   = C1Lim_1z2%Lim_KinInv(1)
       s5i = C1Lim_1z2%Lim_KinInv(2)

       Ec = C1Lim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2

       eta52 = one             
       eta52lim = eta52

       subint_1z24_c1lim = aqb_qbg_sub_aqbq_zi(Ec,C1Lim_1z2%muf(1),z,eta52,eta52lim)       


       respdf = (-one)*respdf*(two/s5i)*(xn*Pgq_spav(z5i)/(one-z5i))&
            * C1Lim_1z2%wgt * subint_1z24_c1lim

       FintNNLO_onloewk_aq(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C2                                                                                                                                                                                                   
    call cut_histo(C2Lim_1z2)

    if (C2Lim_1z2%makecut.or.C2Lim_1z2%flag) then

       kin(3) = zero
       FintNNLO_onloewk_aq(3) = zero

    else

       z5i   = C2Lim_1z2%Lim_KinInv(1)
       s5i = C2Lim_1z2%Lim_KinInv(2)

       Ec = C2Lim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       eta52 = zero
       eta52lim = C2Lim_1z2%Lim_etaij(2,5)

       subint_1z24_c2lim = aqb_qbg_sub_aqbq_zi(Ec,C2Lim_1z2%muf(1),z,eta52,eta52lim)       


      !-- here we use spin-averages since there are no spin correlations for the aa -> e-e+ amplitude
       call res_treeAA_aa(C2Lim_1z2%AmpMom,res_loAA)
       !-- 1=xn*aveqa/aveaa                                                                                                      
       res_loAA = res_loAA * Q_lep2 * Pqq(z5i)/(one-z5i) !-- use Pqq/(1-z) because of definition of z  

       res_tmp(:,1) = Qdn2 * res_loAA
       res_tmp(:,2) = Qup2 * res_loAA

       call get_respdf(aq_lumi,1,1,C2Lim_1z2,res_tmp,respdf)

       respdf = (-one)*respdf*(two/s5i)&
            * C2Lim_1z2%wgt * subint_1z24_c2lim

       FintNNLO_onloewk_aq(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif



    !!!!------------------------!!!!!
    !!!! FLM[1_a,2_q,3,4|5_q]   !!!!!
    !!!!------------------------!!!!!

    
    call kinematics_nlo_z_is(yr=xx,z1=one,z2=one,HardProc=HardProc_12,&
         C1Lim=C1Lim_12,C2Lim=C2Lim_12,compute_etas=.true.)
    HardProc_12%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim_12%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim_12%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard                                                                                                                                                                                                 
    call cut_histo(HardProc_12)
    if (HardProc_12%makecut.or.HardProc_12%flag) then

       kin(4) = zero
       FintNNLO_onloewk_aq(4) = zero
 else

       call res_tree_a_aq(HardProc_12%AmpMom,res_nlo)
       call get_respdf(aq_lumi,1,1,HardProc_12,res_nlo,respdf)

       Ec = HardProc_12%Lim_Ei(1)                      ! Ec = E1 = E2
       E5 = HardProc_12%Lim_Ei(5)                      ! Ec = E1 = E2
       eta51 = HardProc_12%Lim_etaij(1,5)
       eta52 = HardProc_12%Lim_etaij(2,5)
       eta51lim = eta51
       eta52lim = eta52

       subint_124_hard = aqb_qbg_sub_aqbq(Ec,HardProc_12%muf(1),E5,eta52,eta51lim,eta52lim) - &
            aqb_qbg_sub_aqbq_zi_pls(Ec,HardProc_12%muf(1),z,eta52,eta52lim)       

       if (fncheck) then
          print *, "Ec->",Ec,",mu->",Hardproc_12%muf(1),",z->",z,",eta52->",eta52,",E5->",E5,",eta51->",eta51
          print *, "124",subint_124_hard
          pause
       endif


       respdf = respdf*HardProc_12%wgt * subint_124_hard


       kin(4) = respdf(1)
       FintNNLO_onloewk_aq(4) = kin(4)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1                                                                                                                                                                                                   
    call cut_histo(C1Lim_12)

    if (C1Lim_12%makecut.or.C1Lim_12%flag) then

       kin(5) = zero
       FintNNLO_onloewk_aq(5) = zero

    else

       call res_tree_qqb(C1Lim_12%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       call get_respdf(aq_lumi,1,1,C1Lim_12,res_lo,respdf)

       z5i   = C1Lim_12%Lim_KinInv(1)
       s5i = C1Lim_12%Lim_KinInv(2)

       Ec = C1Lim_12%Lim_Ei(1)       ! Ec = E1 = E2
       E5 = C1Lim_12%Lim_Ei(5)       ! Ec = E1 = E2

       eta51 = zero
       eta51lim = C1Lim_12%Lim_etaij(1,5)
       eta52 = one
       eta52lim = eta52

       subint_124_c1lim = aqb_qbg_sub_aqbq(Ec,C1Lim_12%muf(1),E5,eta52,eta51lim,eta52lim) - &
            aqb_qbg_sub_aqbq_zi_pls(Ec,C1Lim_12%muf(1),z,eta52,eta52lim)       

       respdf = (-one)*respdf*(two/s5i)*(xn*Pgq_spav(z5i)/(one-z5i))&
            * C1Lim_12%wgt * subint_124_c1lim

       FintNNLO_onloewk_aq(5) = respdf(1)
       kin(5) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C2                                                                                                                                                                                                   
    call cut_histo(C2Lim_12)

    if (C2Lim_12%makecut.or.C2Lim_12%flag) then

       kin(6) = zero
       FintNNLO_onloewk_aq(6) = zero

    else

       z5i   = C2Lim_12%Lim_KinInv(1)
       s5i = C2Lim_12%Lim_KinInv(2)

       Ec = C2Lim_12%Lim_Ei(1)                      ! Ec = E1 = E2
       E5 = C2Lim_12%Lim_Ei(5)                      ! Ec = E1 = E2

       eta51 = one
       eta51lim = eta51
       eta52 = zero
       eta52lim = C2Lim_12%Lim_etaij(2,5)

       subint_124_c2lim =  aqb_qbg_sub_aqbq(Ec,C2Lim_12%muf(1),E5,eta52,eta51lim,eta52lim) - &
            aqb_qbg_sub_aqbq_zi_pls(Ec,C2Lim_12%muf(1),z,eta52,eta52lim)       



      !-- here we use spin-averages since there are no spin correlations for the aa -> e-e+ amplitude
       call res_treeAA_aa(C2Lim_12%AmpMom,res_loAA)
       !-- 1=xn*aveqa/aveaa                                                                                                      
       res_loAA = res_loAA * Q_lep2 * Pqq(z5i)/(one-z5i) !-- use Pqq/(1-z) because of definition of z  

       res_tmp(:,1) = Qdn2 * res_loAA
       res_tmp(:,2) = Qup2 * res_loAA

       call get_respdf(aq_lumi,1,1,C2Lim_12,res_tmp,respdf)

       respdf = (-one)*respdf*(two/s5i)&
            * C2Lim_12%wgt * subint_124_c2lim

       FintNNLO_onloewk_aq(6) = respdf(1)
       kin(6) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

!!    if (kin(1)*kin(2)*kin(3)*kin(4)*kin(5)*kin(6) .ne. zero) then
!!
!!       print *, "kin in onlo z12", kin(1:6)
!!       print *, "ff",ff(1),ff(1)/kin(1)
!!       pause
!!   endif

    call check_ff(ff,xx,FintNNLO_onloewk_aq)

#if(_withchecks == 1)
!!    FintNLO_ew(1:3) = FintNNLO_onloewk_aq
#endif

  end function xsect_nnlo_onloewk_aq_raoul
  


  function xsect_nnlo_sub12_aq_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_sub12_aq_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--                                                                                                                                            
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintLO_s(3),kin(3)
    real(dp) :: res_lo(2,2), res_tmp(2,2),res_aa,respdf(ipdf)
    real(dp) :: Ec,mu,z,zbar,xis(1:2)
    real(dp) :: intsub_qqb_zzb_u(0:1), intsub_qqb_zzb_d(0:1),intsub_qqb_z12_u, intsub_qqb_z12_d
    real(dp) :: intsub_aa_u, intsub_aa_d

    xsect_nnlo_sub12_aq_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    z   =buff+onet*real(yRnd(kLO_max_full),dp)
    zbar=buff+onet*real(yRnd(kLO_max_full+1),dp)

    if ( ((one-z) .lt. buff_z) .or. ((one-zbar) .lt. buff_z) ) then
       failed_points = failed_points + 1
       return
    endif


#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()

    ! no eta-dependence in the int. subtr. functions
    ! so we can use the usual trick of putting the boosts in the pdfs    
    call kinematics_lo(xx,LOProc)

    Ec = LOProc%AmpMom(1,1)
    LOProc%ids(1:4) = [0,0,id_el,-id_el]

    xis(1:2) = LOProc%PartFrac(1:2)

    call cut_histo(LOProc)
    mu = LOProc%muf(1)

    if (LOProc%makecut.or.LOProc%flag) then

       kin = zero
       FintLO_s = zero

    else

       call res_tree_qqb(LOProc%AmpMom,res_lo)
       call res_treeAA_aa(LOProc%AmpMom,res_aa)

       intsub_qqb_zzb_u = aqb_qbg_sub_qqbz1zbar2(mu,Ec,z,zbar,Qup)
       intsub_qqb_zzb_d = aqb_qbg_sub_qqbz1zbar2(mu,Ec,z,zbar,Qdn)

       intsub_qqb_z12_u = aqb_qbg_sub_qqbz12(mu,Ec,z,Qup)
       intsub_qqb_z12_d = aqb_qbg_sub_qqbz12(mu,Ec,z,Qdn)

       intsub_aa_u = aqb_qbg_sub_aa1z2(mu,Ec,zbar,Qup)
       intsub_aa_d = aqb_qbg_sub_aa1z2(mu,Ec,zbar,Qdn)

       if (fncheck) then
          print *,"mu->",mu,",EC->",Ec,",z->",z,",zb->",zbar
          print *, "zzb u",intsub_qqb_zzb_u(0)
          print *, "zzb d",intsub_qqb_zzb_d(0)
          print *, "z12 u", intsub_qqb_z12_u - intsub_qqb_zzb_u(1)
          print *, "z12 d", intsub_qqb_z12_d - intsub_qqb_zzb_d(1)
          
          print *, "aa u", intsub_aa_u
          print *, "aa d", intsub_aa_d
          pause
       endif

       ! FLM[z1_q,z2_qb]
       res_tmp(:,1) = res_lo(:,1)*intsub_qqb_zzb_d(0)         ! reg+pls
       res_tmp(:,2) = res_lo(:,2)*intsub_qqb_zzb_u(0)

       LOProc%PartFrac(1) = xis(1)/z
       LOProc%PartFrac(2) = xis(2)/zbar

       call get_respdf(aq_lumi,1,1,LOProc,res_tmp,respdf)
       respdf = respdf*LOProc%wgt/z/zbar

       kin(1) = respdf(1)

       !print *, '[z,zb] 1',kin(1)
       
       call fill_histo(respdf,vegasweight)

       ! FLM[z1_q,2_qb]
       res_tmp(:,1) = res_lo(:,1)*( intsub_qqb_z12_d - intsub_qqb_zzb_d(1))
       res_tmp(:,2) = res_lo(:,2)*( intsub_qqb_z12_u - intsub_qqb_zzb_u(1)) 

       LOProc%PartFrac(1) = xis(1)/z
       LOProc%PartFrac(2) = xis(2)

       call get_respdf(aq_lumi,1,1,LOProc,res_tmp,respdf)
       respdf = respdf*LOProc%wgt/z

       kin(2) = respdf(1)

       !print *, '[z,2 ] 2',kin(2)
       
       call fill_histo(respdf,vegasweight)


       ! FLM[1_a,zbar 2_a]
       res_tmp(:,1) = res_aa*intsub_aa_d
       res_tmp(:,2) = res_aa*intsub_aa_u

       LOProc%PartFrac(1) = xis(1)
       LOProc%PartFrac(2) = xis(2)/zbar

       call get_respdf(aq_lumi,1,1,LOProc,res_tmp,respdf)
       respdf = respdf*LOProc%wgt/zbar

       kin(3) = respdf(1)
       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 3',kin(3)
       !pause
       
    endif

    ff(1) = sum(kin(1:3))

    call close_histo()


  end function xsect_nnlo_sub12_aq_raoul



  function xsect_nnlo_subvqcd_aq_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_subvqcd_aq_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--                                                                                                                                            
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintLO_s(1),kin(1)
    real(dp) :: res_tree(2,2), res_tmp(2,2),res_loop(2,2),respdf(ipdf)
    real(dp) :: q2,mu,z,xis(1:2)
    real(dp) :: intsub_qqb(-1:1,nintsub)


    xsect_nnlo_subvqcd_aq_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    z   =buff+onet*real(yRnd(kLO_max_full),dp)

    if ( ((one-z) .lt. buff_z) ) then
       failed_points = failed_points + 1
       return
    endif


#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()

    ! no eta-dependence in the int. subtr. functions
    ! so we can use the usual trick of putting the boosts in the pdfs    
    call kinematics_lo(xx,LOProc)

    LOProc%ids(1:4) = [0,0,id_el,-id_el]

    xis(1:2) = LOProc%PartFrac(1:2)
    q2 = two*scr(LOProc%AmpMom(1:4,1),LOProc%AmpMom(1:4,2))

    call cut_histo(LOProc)
    mu = LOProc%muf(1)

    if (LOProc%makecut.or.LOProc%flag) then

       kin = zero
       FintLO_s = zero

    else
       
       call res_qcdloop_qqb(LOProc%AmpMom,res_tree,res_loop)

       res_tmp(:,1) = Qdn**2 * res_loop(:,1)
       res_tmp(:,2) = Qup**2 * res_loop(:,2)

       intsub_qqb = sub_a_aq_qqb_z(q2,LOProc%muf**2,z)

       ! FLV[z1_q,2_qb]

       LOProc%PartFrac(1) = xis(1)/z
       LOProc%PartFrac(2) = xis(2)

       call get_respdf(aq_lumi,1,1,LOProc,res_tmp,respdf)
       respdf = respdf*LOProc%wgt/z *intsub_qqb(0,1)

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = kin(1)

    call close_histo()


  end function xsect_nnlo_subvqcd_aq_raoul





  




  
end module mod_xsects_nnlo_s_aq_raoul
  

module mod_xsects_nnlo_s_gq_raoul
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
  use mod_subtrfn_nnlo_z_gq

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
  real(dp), public, save :: FintNNLO_onloqcd_gq(4)
  real(dp), public, save :: FintNNLO_onloewk_is_gq(6)
  real(dp), public, save :: FintNNLO_onloewk_5i_gq(6)
#endif

  private

  public :: xsect_nnlo_sub12_gq_raoul
  public :: xsect_nnlo_onloqcd_gq_raoul
  public :: xsect_nnlo_onloewk_is_gq_raoul
  public :: xsect_nnlo_onloewk_fs_53_gq_raoul, xsect_nnlo_onloewk_fs_54_gq_raoul
  public :: xsect_nnlo_subvewk_gq_raoul


contains


  !--------------------------------------------------
  !-- subtraction counterterms
  !--------------------------------------------------


 function xsect_nnlo_onloqcd_gq_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_onloqcd_gq_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc_1z2,C1Lim_1z2,HardProc_12,C1Lim_12
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: kin(4)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(2,2),res_lo(2,2)
    real(dp)    :: z,z5i,s5i,Ec,E1,E2,E3,E4,E5,eta32,eta34,eta42,eta51,eta52,eta53,eta45,eta51log
    real(dp)    :: subint_1z24_hard_u,subint_1z24_c1lim_u,subint_1z24_hard_d,subint_1z24_c1lim_d
    real(dp)    :: subint_124_hard_u,subint_124_c1lim_u,subint_124_hard_ub,subint_124_c1lim_ub
    real(dp)    :: subint_124_hard_d,subint_124_c1lim_d,subint_124_hard_db,subint_124_c1lim_db

    xsect_nnlo_onloqcd_gq_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z=buff+onet*real(yRnd(kNLO_max_full),dp)

!!    xx(1) =   0.10858601476543928_dp    
!!    xx(2)  =   0.36772289726337964_dp    
!!    xx(3)  =    5.4784999255415655E-002_dp
!!    xx(4)  =    5.6489252055865755E-003_dp
!!    xx(5)  =   0.27148247080504168_dp    
!!    xx(6)  =   0.49417661293139803_dp    
!!    xx(7)  =    5.3751152034694251E-002_dp
!!    xx(8)  =   0.31692953172783789_dp    
!!    z      =   0.46703123509286615_dp     


#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if ( (xx(xRHO).lt.buff_r) .or. ((one-z) .lt. buff_z) ) then
       failed_points = failed_points + 1
       return
    endif

!    xx(xrho) = 1E-10_dp
!    z = one - 1E-8_dp

    call open_histo()

    !!!!------------------------!!!!!
    !!!! FLM[1_g,z.2_q,3,4|5_q] !!!!!
    !!!!------------------------!!!!!                                                                                                                                                                        
    call kinematics_nlo_z_is(yr=xx,z1=one,z2=z,HardProc=HardProc_1z2,&
         C1Lim=C1Lim_1z2,compute_etas=.true.)
    HardProc_1z2%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim_1z2%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard                                                                                                                                                                                                 
    call cut_histo(HardProc_1z2)
    if (HardProc_1z2%makecut.or.HardProc_1z2%flag) then

       kin(1) = zero
       FintNNLO_onloqcd_gq(1) = zero

    else

       call res_tree_g_gq(HardProc_1z2%AmpMom,res_nlo)

       Ec = HardProc_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       subint_1z24_hard_u = gqb_qba_sub_gqbq_zi(Ec,HardProc_1z2%muf(1),z,Qup)
       subint_1z24_hard_d = gqb_qba_sub_gqbq_zi(Ec,HardProc_1z2%muf(1),z,Qdn)
       if (fncheck) then       
          print *, "Ec->",Ec,",mu->",HardProc_1z2%muf(1),",z->",z
          print *, "1z24 u",subint_1z24_hard_u
          print *, "1z24 d",subint_1z24_hard_d
          pause
       endif

       res_nlo(:,1) = res_nlo(:,1)*subint_1z24_hard_d
       res_nlo(:,2) = res_nlo(:,2)*subint_1z24_hard_u

       call get_respdf(gq_lumi,1,1,HardProc_1z2,res_nlo,respdf)

       respdf = respdf*HardProc_1z2%wgt

       kin(1) = respdf(1)
       FintNNLO_onloqcd_gq(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

   !-- C1                                                                                                                                                                                                   
    call cut_histo(C1Lim_1z2)

    if (C1Lim_1z2%makecut.or.C1Lim_1z2%flag) then

       kin(2) = zero
       FintNNLO_onloqcd_gq(2) = zero

    else

       call res_tree_qqb(C1Lim_1z2%AmpMom,res_lo)


       z5i   = C1Lim_1z2%Lim_KinInv(1)
       s5i = C1Lim_1z2%Lim_KinInv(2)

       Ec = C1Lim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2

       subint_1z24_c1lim_u = gqb_qba_sub_gqbq_zi(Ec,C1Lim_1z2%muf(1),z,Qup)
       subint_1z24_c1lim_d = gqb_qba_sub_gqbq_zi(Ec,C1Lim_1z2%muf(1),z,Qdn)

       res_lo(:,1) = res_lo(:,1)*subint_1z24_c1lim_d
       res_lo(:,2) = res_lo(:,2)*subint_1z24_c1lim_u

       call get_respdf(gq_lumi,1,1,C1Lim_1z2,res_lo,respdf)

       !-- use Pgq so that in the soft limit we do not have 1/[1-(1-x)]
       !-- Tr = Cf * aveqg/aveqq
       respdf = (-one)*respdf*(two/s5i)*(tr*Pgq_spav(z5i)/(one-z5i)) &
            * C1Lim_1z2%wgt

       FintNNLO_onloqcd_gq(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif


    !!!!------------------------!!!!!
    !!!! FLM[1_g,2_q,3,4|5_q] !!!!!
    !!!!------------------------!!!!!                                                                                                                                                                        
    call kinematics_nlo_z_is(yr=xx,z1=one,z2=one,HardProc=HardProc_12,&
         C1Lim=C1Lim_12,compute_etas=.true.)
    HardProc_12%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim_12%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard                                                                                                                                                                                                 
    call cut_histo(HardProc_12)
    if (HardProc_12%makecut.or.HardProc_12%flag) then

       kin(3) = zero
       FintNNLO_onloqcd_gq(3) = zero

    else

       call res_tree_g_gq(HardProc_12%AmpMom,res_nlo)

       Ec = HardProc_12%Lim_Ei(1)                      ! Ec = E1 = E2
       E1 = HardProc_12%Lim_Ei(1)                      
       E2 = HardProc_12%Lim_Ei(2)                      
       E3 = HardProc_12%Lim_Ei(3)                      
       E4 = HardProc_12%Lim_Ei(4)                      
       E5 = HardProc_12%Lim_Ei(5)                      
       eta32 = HardProc_12%Lim_etaij(2,3)
       eta34 = HardProc_12%Lim_etaij(3,4)
       eta42 = HardProc_12%Lim_etaij(2,4)
       eta51 = HardProc_12%Lim_etaij(1,5)
       eta52 = HardProc_12%Lim_etaij(2,5)
       eta53 = HardProc_12%Lim_etaij(3,5)
       eta45 = HardProc_12%Lim_etaij(4,5)
       eta51log = eta51

       subint_124_hard_u = gqb_qba_sub_gqbq(Ec,HardProc_12%muf(1),E2,E3,E4,E5,eta32,eta34,eta42,eta51,eta52,eta53,eta45,eta51log,-Qup,Q_lep) - &
            gqb_qba_sub_gqbq_zi_pls(Ec,HardProc_12%muf(1),z,Qup)
       subint_124_hard_ub = gqb_qba_sub_gqbq(Ec,HardProc_12%muf(1),E2,E3,E4,E5,eta32,eta34,eta42,eta51,eta52,eta53,eta45,eta51log,Qup,Q_lep) - &
            gqb_qba_sub_gqbq_zi_pls(Ec,HardProc_12%muf(1),z,-Qup)
       subint_124_hard_d = gqb_qba_sub_gqbq(Ec,HardProc_12%muf(1),E2,E3,E4,E5,eta32,eta34,eta42,eta51,eta52,eta53,eta45,eta51log,-Qdn,Q_lep) - &
            gqb_qba_sub_gqbq_zi_pls(Ec,HardProc_12%muf(1),z,Qdn)

       subint_124_hard_db = gqb_qba_sub_gqbq(Ec,HardProc_12%muf(1),E2,E3,E4,E5,eta32,eta34,eta42,eta51,eta52,eta53,eta45,eta51log,Qdn,Q_lep) - &
            gqb_qba_sub_gqbq_zi_pls(Ec,HardProc_12%muf(1),z,-Qdn)

       if (fncheck) then       
          print *, "Ec->",Ec,",mu->",HardProc_12%muf(1),",z->",z,",E2->",E2,",E3->",E3,",E4->",E4,",E5->",E5,",eta32->",eta32,",eta34->",eta34,",eta42->",eta42,",eta51->",eta51,",eta52->",eta52,",eta53->",eta53,",eta45->",eta45
          print *, "124 u",subint_124_hard_u
          print *, "124 d",subint_124_hard_d
          print *, "124 ub",subint_124_hard_ub
          print *, "124 db",subint_124_hard_db

          pause
       endif


       res_nlo(1,1) = res_nlo(1,1)*subint_124_hard_db
       res_nlo(2,1) = res_nlo(2,1)*subint_124_hard_d
       res_nlo(1,2) = res_nlo(1,2)*subint_124_hard_ub
       res_nlo(2,2) = res_nlo(2,2)*subint_124_hard_u

       call get_respdf(gq_lumi,1,1,HardProc_12,res_nlo,respdf)

       respdf = respdf*HardProc_12%wgt

       kin(3) = respdf(1)
       FintNNLO_onloqcd_gq(3) = kin(3)

       call fill_histo(respdf,vegasweight)

    endif

   !-- C1                                                                                                                                                                                                   
    call cut_histo(C1Lim_12)

    if (C1Lim_12%makecut.or.C1Lim_12%flag) then

       kin(4) = zero
       FintNNLO_onloqcd_gq(4) = zero

    else

       call res_tree_qqb(C1Lim_12%AmpMom,res_lo)


       z5i   = C1Lim_12%Lim_KinInv(1)
       s5i = C1Lim_12%Lim_KinInv(2)

       Ec = C1Lim_12%Lim_Ei(1)                      ! Ec = E1 = E2
       E1 = C1Lim_12%Lim_Ei(1)                      
       E2 = C1Lim_12%Lim_Ei(2)                      
       E3 = C1Lim_12%Lim_Ei(3)                      
       E4 = C1Lim_12%Lim_Ei(4)                      
       E5 = C1Lim_12%Lim_Ei(5)                      
       eta32 = C1Lim_12%Lim_etaij(2,3)
       eta34 = C1Lim_12%Lim_etaij(3,4)
       eta42 = C1Lim_12%Lim_etaij(2,4)
       eta51 = zero 
       eta52 = one                               ! C51 eta52 = eta12 = 1
       eta53 = C1Lim_12%Lim_etaij(1,3)           ! C51 eta53 = eta13
       eta45 = C1Lim_12%Lim_etaij(1,4)           ! C51 eta54 = eta14
       eta51log = C1Lim_12%Lim_etaij(1,5)
       

       subint_124_c1lim_u = gqb_qba_sub_gqbq(Ec,C1Lim_12%muf(1),E2,E3,E4,E5,eta32,eta34,eta42,eta51,eta52,eta53,eta45,eta51log,-Qup,Q_lep) - &
            gqb_qba_sub_gqbq_zi_pls(Ec,C1Lim_12%muf(1),z,Qup)
       subint_124_c1lim_ub = gqb_qba_sub_gqbq(Ec,C1Lim_12%muf(1),E2,E3,E4,E5,eta32,eta34,eta42,eta51,eta52,eta53,eta45,eta51log,Qup,Q_lep) - &
            gqb_qba_sub_gqbq_zi_pls(Ec,C1Lim_12%muf(1),z,-Qup)
       subint_124_c1lim_d = gqb_qba_sub_gqbq(Ec,C1Lim_12%muf(1),E2,E3,E4,E5,eta32,eta34,eta42,eta51,eta52,eta53,eta45,eta51log,-Qdn,Q_lep) - &
            gqb_qba_sub_gqbq_zi_pls(Ec,C1Lim_12%muf(1),z,Qdn)
       subint_124_c1lim_db = gqb_qba_sub_gqbq(Ec,C1Lim_12%muf(1),E2,E3,E4,E5,eta32,eta34,eta42,eta51,eta52,eta53,eta45,eta51log,Qdn,Q_lep) - &
            gqb_qba_sub_gqbq_zi_pls(Ec,C1Lim_12%muf(1),z,-Qdn)


       res_lo(1,1) = res_lo(1,1)*subint_124_c1lim_db
       res_lo(2,1) = res_lo(2,1)*subint_124_c1lim_d
       res_lo(1,2) = res_lo(1,2)*subint_124_c1lim_ub
       res_lo(2,2) = res_lo(2,2)*subint_124_c1lim_u


       call get_respdf(gq_lumi,1,1,C1Lim_12,res_lo,respdf)

       !-- use Pgq so that in the soft limit we do not have 1/[1-(1-x)]
       !-- Tr = Cf * aveqg/aveqq
       respdf = (-one)*respdf*(two/s5i)*(tr*Pgq_spav(z5i)/(one-z5i)) &
            * C1Lim_12%wgt

       FintNNLO_onloqcd_gq(4) = respdf(1)
       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    

    ff(1) = sum(kin)

    call close_histo()

!!    if (kin(1)*kin(2)*kin(3) .ne. zero) then
!!       print *, "kin in onlo z12", kin(1:4)
!!       print *, "ff",ff(1)
!!       pause
!!    endif




  end function xsect_nnlo_onloqcd_gq_raoul



  function xsect_nnlo_onloewk_is_gq_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_onloewk_is_gq_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc_z12,SLim_z12,C1Lim_z12,C2Lim_z12,SC1Lim_z12,SC2Lim_z12
    !--
    real(dp)    :: xx(kNLO_max_full),z,damp
    real(dp)    :: respdf(ipdf),kin(18)
    real(dp)    :: res_lo(2,2),res_lo_tmp(2,2),res_nlo(2,2)
    real(dp)    :: Ec, e5, z5i, s5i,eik(4),eta51,eta52,eta53,eta54,eta51lim
    real(dp)    :: subint_z124_hard,subint_z124_slim,subint_z124_c1lim,subint_z124_c2lim,subint_z124_sc1lim,subint_z124_sc2lim

    xsect_nnlo_onloewk_is_gq_raoul = 0

    ff(1) = zero
    kin   = zero

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

    call open_histo()
    
    !! ----------------------------- !!
    !! FLM[z1,2,3,4 | 5gammma]       !!
    !! ----------------------------- !!

    call kinematics_nlo_z_is(yr=xx,z1=z,z2=one,HardProc=HardProc_z12,C1Lim=C1Lim_z12,C2Lim=C2Lim_z12,SLim=SLim_z12,SC1Lim=SC1Lim_z12,SC2Lim=SC2Lim_z12,compute_etas=.true.)
    HardProc_z12%ids(1:5) = [0,0,id_el,-id_el,id_a]
    C1Lim_z12%ids(1:4)    = [0,0,id_el,-id_el]
    SLim_z12%ids(1:4)     = [0,0,id_el,-id_el]
    C2Lim_z12%ids(1:4)    = [0,0,id_el,-id_el]
    SC1Lim_z12%ids(1:4)   = [0,0,id_el,-id_el]
    SC2Lim_z12%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard

    call cut_histo(HardProc_z12)
    if (HardProc_z12%makecut.or.HardProc_z12%flag) then

       kin(1) = zero
       FintNNLO_onloewk_is_gq(1) = zero

    else

       call res_tree_a_qqb(HardProc_z12%AmpMom,res_nlo)
       call partition_nlo_qed(HardProc_z12,damp,12)

       Ec = HardProc_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = HardProc_z12%Lim_etaij(1,5)
       eta52 = HardProc_z12%Lim_etaij(2,5)
       eta53 = HardProc_z12%Lim_etaij(3,5)
       eta54 = HardProc_z12%Lim_etaij(4,5)
       eta51lim = eta51
       subint_z124_hard = gqb_qba_sub_qqba_zi(Ec,HardProc_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)

       if (fncheck) then       
          print *, "Ec->",Ec,",mu->",HardProc_z12%muf(1),",z->",z,",eta51->",eta51,",eta52->",eta52,",eta53->",eta53,",eta45->",eta54
          print *, "z124 ",subint_z124_hard
          pause
       endif

       call get_respdf(gq_lumi,1,1,HardProc_z12,res_nlo,respdf)

       respdf = respdf*HardProc_z12%wgt*subint_z124_hard * damp

       kin(1) = respdf(1)
       FintNNLO_onloewk_is_gq(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1

    call cut_histo(C1Lim_z12)
    if (C1Lim_z12%makecut.or.C1Lim_z12%flag) then

       kin(2) = zero
       FintNNLO_onloewk_is_gq(2) = zero

    else

       call res_tree_qqb(C1Lim_z12%AmpMom,res_lo)

       Ec = C1Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51lim = C1Lim_z12%Lim_etaij(1,5)

       eta53 = C1Lim_z12%Lim_etaij(1,3)              ! C51 eta53 = eta13
       eta54 = C1Lim_z12%Lim_etaij(1,4)              ! C51 eta43 = eta14
       
       eta51 = zero
       eta52 = one
       subint_z124_c1lim = gqb_qba_sub_qqba_zi(Ec,C1Lim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)


       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       call get_respdf(gq_lumi,1,1,C1Lim_z12,res_lo,respdf)

       z5i = C1Lim_z12%Lim_KinInv(1)
       s5i = C1Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z5i)/(one-z5i))&
            * C1Lim_z12%wgt * subint_z124_c1lim

       FintNNLO_onloewk_is_gq(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C2

    call cut_histo(C2Lim_z12)
    if (C2Lim_z12%makecut.or.C2Lim_z12%flag) then

       kin(3) = zero
       FintNNLO_onloewk_is_gq(3) = zero

    else
       call res_tree_qqb(C2Lim_z12%AmpMom,res_lo)

       Ec = C2Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta53 = C2Lim_z12%Lim_etaij(2,3)              ! C52 eta53 = eta23
       eta54 = C2Lim_z12%Lim_etaij(2,4)              ! C52 eta54 = eta24

       eta51 = one
       eta51lim = eta51
       eta52 = zero

       subint_z124_c2lim = gqb_qba_sub_qqba_zi(Ec,C2Lim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)

       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       call get_respdf(gq_lumi,1,1,C2Lim_z12,res_lo,respdf)

       z5i = C2Lim_z12%Lim_KinInv(1)
       s5i = C2Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z5i)/(one-z5i))&
            * C2Lim_z12%wgt * subint_z124_c2lim


       FintNNLO_onloewk_is_gq(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- S

    call cut_histo(SLim_z12)
    if (SLim_z12%makecut.or.SLim_z12%flag) then

       kin(4) = zero
       FintNNLO_onloewk_is_gq(4) = zero

    else

       call res_tree_qqb(SLim_z12%AmpMom,res_lo)
       call partition_nlo_qed(SLim_z12,damp,12)


       Ec = SLim_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = SLim_z12%Lim_etaij(1,5)
       eta52 = SLim_z12%Lim_etaij(2,5)
       eta53 = SLim_z12%Lim_etaij(3,5)
       eta54 = SLim_z12%Lim_etaij(4,5)

       subint_z124_slim = gqb_qba_sub_qqba_zi(Ec,SLim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51)

       call get_qed_eik(charges_ns,SLim_z12%Lim_etaij,[1,2,3,4],5,eik)
       res_lo_tmp(:,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_lo_tmp(:,2) = res_lo(:,2) * eik(3:4) !-- up
       call get_respdf(gq_lumi,1,1,SLim_z12,res_lo_tmp,respdf)


       e5 = SLim_z12%Lim_Ei(5)

       respdf = -respdf*SLim_z12%wgt*damp/e5**2 * subint_z124_slim

       kin(4) = respdf(1)
       FintNNLO_onloewk_is_gq(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif



    !-- SC1

    call cut_histo(SC1Lim_z12)
    if (SC1Lim_z12%makecut.or.SC1Lim_z12%flag) then

       kin(5) = zero
       FintNNLO_onloewk_is_gq(5) = zero

    else

       call res_tree_qqb(SC1Lim_z12%AmpMom,res_lo)

       Ec = SC1Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = zero
       eta52 = one
       eta51lim = SC1Lim_z12%Lim_etaij(1,5)

       eta53 = SC1Lim_z12%Lim_etaij(1,3)              ! C51 eta53 = eta13
       eta54 = SC1Lim_z12%Lim_etaij(1,4)              ! C51 eta43 = eta14

       subint_z124_sc1lim = gqb_qba_sub_qqba_zi(Ec,SC1Lim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)

       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)



       call get_respdf(gq_lumi,1,1,SC1Lim_z12,res_lo,respdf)

       e5 = SC1Lim_z12%Lim_KinInv(1)
       eta51 = SC1Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = respdf/e5**2/eta51* SC1Lim_z12%wgt * subint_z124_sc1lim

       FintNNLO_onloewk_is_gq(5) = respdf(1)
       kin(5) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- SC2

    call cut_histo(SC2Lim_z12)
    if (SC2Lim_z12%makecut.or.SC2Lim_z12%flag) then

       kin(6) = zero
       FintNNLO_onloewk_is_gq(6) = zero

    else

       call res_tree_qqb(SC2Lim_z12%AmpMom,res_lo)

       Ec = SC2Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta53 = SC2Lim_z12%Lim_etaij(2,3)              ! C52 eta53 = eta23
       eta54 = SC2Lim_z12%Lim_etaij(2,4)              ! C52 eta54 = eta24

       eta51 = one
       eta51lim = eta51
       eta52 = zero

       subint_z124_sc2lim = gqb_qba_sub_qqba_zi(Ec,SC2Lim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)

       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)



       call get_respdf(gq_lumi,1,1,SC2Lim_z12,res_lo,respdf)

       e5 = SC2Lim_z12%Lim_KinInv(1)
       eta52 = SC2Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]

       respdf = respdf/e5**2/eta52* SC2Lim_z12%wgt * subint_z124_sc2lim

       FintNNLO_onloewk_is_gq(6) = respdf(1)
       kin(6) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif




    if (is_nan(kin)) then
       print *, "kin",kin(1:6)
       stop
    endif

    ff(1) = sum(kin(1:6))


!!    if (kin(1)*kin(2)*kin(3)*kin(4)*kin(5)*kin(6) .ne. zero) then
!!       print *, "kin in onlo z12", kin(1:6)
!!       print *, "ff",ff(1)
!!       pause
!!    endif


    call close_histo()



  end function xsect_nnlo_onloewk_is_gq_raoul



  function xsect_nnlo_onloewk_fs_53_gq_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_onloewk_fs_53_gq_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nnlo_onloewk_fs_53_gq_raoul = xsect_nnlo_onloewk_5i_gq_raoul(yRnd,ff,vegasweight,3,4)
    
  end function xsect_nnlo_onloewk_fs_53_gq_raoul

  function xsect_nnlo_onloewk_fs_54_gq_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_onloewk_fs_54_gq_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nnlo_onloewk_fs_54_gq_raoul = xsect_nnlo_onloewk_5i_gq_raoul(yRnd,ff,vegasweight,4,3)
    
  end function xsect_nnlo_onloewk_fs_54_gq_raoul




  function xsect_nnlo_onloewk_5i_gq_raoul(yRnd,ff,vegasweight,icoll,jother)
    integer :: xsect_nnlo_onloewk_5i_gq_raoul,icoll,jother
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc_z12,SLim_z12,CLim_z12,SCLim_z12
    !-- 
    real(dp)    :: xx(kNLO_max_full),z,damp
    real(dp)    :: respdf(ipdf),kin(12)
    real(dp)    :: res_lo(2,2),res_lo_tmp(2,2),res_nlo(2,2)
    real(dp)    :: Ec, e5, z5i, s5i,eik(4),eta51,eta52,eta53,eta54,eta51lim
    real(dp)    :: subint_z124_hard,subint_z124_slim,subint_z124_clim,subint_z124_sclim

    xsect_nnlo_onloewk_5i_gq_raoul = 0

    ff(1) = zero
    kin   = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z=buff+onet*real(yRnd(kNLO_max_full),dp)

!!    xx = (/0.24639744183397252_dp,       0.28115138657488475_dp,        4.6186198037180770E-002_dp,&
!!         5.9983770547065468E-002_dp,   4.1578741242558249E-004_dp,  0.59597880312345652_dp,       0.41114427370245565_dp,       0.89061826774764818_dp/)


!!#if (_withchecks == 1)
!!    if (override) then
!!       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
!!       print *, 'overriding input'
!!    endif
!!#endif

    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. ((one-z) .lt. buff_z) ) then
       failed_points = failed_points + 1
       return
    endif
!    xx(xE)   = 1E-8_dp
!    xx(xrho) = 1E-8_dp
    
    call open_histo()

    !! ----------------------------- !!
    !! FLM[z1,2,3,4 | 5gammma]       !!
    !! ----------------------------- !!

    call kinematics_nlo_z_fs(xx,z,one,icoll,jother,HardProc_z12,CLim_z12,SCLim_z12,SLim_z12)
    HardProc_z12%ids(1:5) = [0,0,id_el,-id_el,id_a]
    CLim_z12%ids(1:4)    = [0,0,id_el,-id_el]
    SLim_z12%ids(1:4)     = [0,0,id_el,-id_el]
    SCLim_z12%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard

    call cut_histo(HardProc_z12)
    if (HardProc_z12%makecut.or.HardProc_z12%flag) then

       kin(1) = zero
       FintNNLO_onloewk_5i_gq(1) = zero

    else

       call res_tree_a_qqb(HardProc_z12%AmpMom,res_nlo)
       call partition_nlo_qed(HardProc_z12,damp,icoll)

       Ec = HardProc_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = HardProc_z12%Lim_etaij(1,5)
       eta52 = HardProc_z12%Lim_etaij(2,5)
       eta53 = HardProc_z12%Lim_etaij(3,5)
       eta54 = HardProc_z12%Lim_etaij(4,5)
       eta51lim = eta51

       subint_z124_hard = gqb_qba_sub_qqba_zi(Ec,HardProc_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)

       call get_respdf(gq_lumi,1,1,HardProc_z12,res_nlo,respdf)
       respdf = respdf*HardProc_z12%wgt*subint_z124_hard * damp

       kin(1) = respdf(1)
       FintNNLO_onloewk_5i_gq(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C

    call cut_histo(CLim_z12)
    if (CLim_z12%makecut.or.CLim_z12%flag) then

       kin(2) = zero
       FintNNLO_onloewk_5i_gq(2) = zero

    else

       call res_tree_qqb(CLim_z12%AmpMom,res_lo)

       Ec = CLim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       if (icoll .eq. 3 .and. jother .eq. 4) then
          eta51 = CLim_z12%Lim_etaij(1,3)
          eta52 = CLim_z12%Lim_etaij(2,3)
          eta53 = zero
          eta54 = CLim_z12%Lim_etaij(3,4)
       elseif (icoll .eq. 4 .and. jother .eq. 3) then
          eta51 = CLim_z12%Lim_etaij(1,4)
          eta52 = CLim_z12%Lim_etaij(2,4)
          eta53 = CLim_z12%Lim_etaij(3,4)
          eta54 = zero
       else
          print *, "incorrect icoll and jother", icoll, jother
          stop
       endif
       eta51lim = eta51


       subint_z124_clim = gqb_qba_sub_qqba_zi(Ec,CLim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)

       call get_respdf(gq_lumi,1,1,CLim_z12,res_lo,respdf)

       z5i = CLim_z12%Lim_KinInv(1)
       s5i = CLim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]

       respdf = (-one)*respdf*(two/s5i)*(Q_lep2*Pqg(z5i))&
            * CLim_z12%wgt * subint_z124_clim

       FintNNLO_onloewk_5i_gq(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif


    !-- S and SC

    call cut_histo(SLim_z12)
    if (SLim_z12%makecut.or.SLim_z12%flag) then

       kin(3) = zero
       FintNNLO_onloewk_5i_gq(3) = zero

    else

       call res_tree_qqb(SLim_z12%AmpMom,res_lo)
       call partition_nlo_qed(SLim_z12,damp,icoll)


       ! soft
       Ec = SLim_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = SLim_z12%Lim_etaij(1,5)
       eta52 = SLim_z12%Lim_etaij(2,5)
       eta53 = SLim_z12%Lim_etaij(3,5)
       eta54 = SLim_z12%Lim_etaij(4,5)

       subint_z124_slim = gqb_qba_sub_qqba_zi(Ec,SLim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51)

       call get_qed_eik(charges_ns,SLim_z12%Lim_etaij,[1,2,3,4],5,eik)
       res_lo_tmp(:,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_lo_tmp(:,2) = res_lo(:,2) * eik(3:4) !-- up
       call get_respdf(gq_lumi,1,1,SLim_z12,res_lo_tmp,respdf)

       e5 = SLim_z12%Lim_Ei(5)

       respdf = -respdf*SLim_z12%wgt*damp/e5**2 * subint_z124_slim       

       kin(3) = respdf(1)
       FintNNLO_onloewk_5i_gq(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif



    !-- SC

    call cut_histo(SCLim_z12)
    if (SCLim_z12%makecut.or.SCLim_z12%flag) then

       kin(4) = zero
       FintNNLO_onloewk_5i_gq(4) = zero

    else

       call res_tree_qqb(SCLim_z12%AmpMom,res_lo)

       Ec = SCLim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       if (icoll .eq. 3 .and. jother .eq. 4) then
          eta51 = SCLim_z12%Lim_etaij(1,3)
          eta52 = SCLim_z12%Lim_etaij(2,3)
          eta53 = zero
          eta54 = SCLim_z12%Lim_etaij(3,4)
       elseif (icoll .eq. 4 .and. jother .eq. 3) then
          eta51 = SCLim_z12%Lim_etaij(1,4)
          eta52 = SCLim_z12%Lim_etaij(2,4)
          eta53 = SCLim_z12%Lim_etaij(3,4)
          eta54 = zero
       else
          print *, "incorrect icoll and jother", icoll, jother
          stop
       endif
       eta51lim = eta51

       subint_z124_sclim = gqb_qba_sub_qqba_zi(Ec,SCLim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)

       call get_respdf(gq_lumi,1,1,SCLim_z12,res_lo,respdf)

       e5 = SCLim_z12%Lim_KinInv(1)
       eta51 = SCLim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = respdf/e5**2/eta51* SCLim_z12%wgt * subint_z124_sclim

       FintNNLO_onloewk_5i_gq(4) = respdf(1)
       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif


    ff(1) = sum(kin(1:4))


!!    if (kin(1)*kin(2)*kin(3)*kin(4)*kin(5)*kin(6) .ne. zero) then
!!       print *, "kin in onlo z12", kin(1:6)
!!       print *, "ff",ff(1)
!!       pause
!!    endif


    call close_histo()




  end function xsect_nnlo_onloewk_5i_gq_raoul
  



  function xsect_nnlo_sub12_gq_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_sub12_gq_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc_z12,LOProc_zzb
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintLO_s(4),kin(2)
    real(dp) :: res_lo(2,2), res_tmp(2,2),respdf(ipdf)
    real(dp) :: Ec,E3,E4,eta31,eta32,eta41,eta42,eta34
    real(dp) :: intsub_zzb_pls_u,intsub_zzb_pls_ub,intsub_zzb_pls_d,intsub_zzb_pls_db
    real(dp) :: intsub_zzb_u,intsub_zzb_ub,intsub_zzb_d,intsub_zzb_db
    real(dp) :: intsub_z_u,intsub_z_ub,intsub_z_d,intsub_z_db
    real(dp) :: z,zbar,mu,ometa34

    xsect_nnlo_sub12_gq_raoul = 0

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


!    xx =[0.75383452378766036_dp,       0.42334610385648064_dp,       0.53206811774833285_dp,        2.0772771618513354E-002_dp,  0.95530224700033917_dp]
!    z = 0.99999996779763245_dp
!    zbar = 0.44012375727571074_dp

    ! debug
!!    z = 0.99999999_dp
!!    zbar=0.9999999_dp
    !!


#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif
        
    call open_histo()

    call kinematics_lo_zzb(xx,z,  one, LOProc_z12)   
    call kinematics_lo_zzb(xx,z,  zbar,LOProc_zzb)  

    LOProc_z12%ids(1:4) = [0,0,id_el,-id_el]
    LOProc_zzb%ids(1:4) = [0,0,id_el,-id_el]

    !! FLM[z1,2]
    call cut_histo(LOProc_z12)

    if (LOProc_z12%makecut.or.LOProc_z12%flag) then

       kin(1) = zero
       FintLO_s(1) = zero

    else

       call res_tree_qqb(LOProc_z12%AmpMom,res_lo)
       
       mu = LOProc_z12%muf(1)

       E3 = LOProc_z12%Lim_Ei(3)
       E4 = LOProc_z12%Lim_Ei(4)
       
       eta31 = LOProc_z12%Lim_etaij(3,1)
       eta32 = LOProc_z12%Lim_etaij(3,2)
       eta41 = LOProc_z12%Lim_etaij(4,1)
       eta42 = LOProc_z12%Lim_etaij(4,2)
       eta34 = LOProc_z12%Lim_etaij(3,4)
       ometa34 = LOProc_z12%Lim_KinInv(6)

       Ec = LOProc_z12%Lim_Ei(1)    ! = E1

       ! get the relevant int. subtr. functions
       intsub_zzb_pls_u  = gqb_qba_sub_qqbz1zbar2_pls(mu,z,zbar,Ec,Qup)
       intsub_zzb_pls_d  = gqb_qba_sub_qqbz1zbar2_pls(mu,z,zbar,Ec,Qdn)
       intsub_zzb_pls_ub = intsub_zzb_pls_u                  ! ~ Qq^2
       intsub_zzb_pls_db = intsub_zzb_pls_d                  ! ~ Qq^2

       intsub_z_d   = gqb_qba_sub_qqbz12(mu,z,Ec,E3,E4,eta31,eta32,eta41,eta42,eta34,ometa34,Qdn,Q_lep)
       intsub_z_u   = gqb_qba_sub_qqbz12(mu,z,Ec,E3,E4,eta31,eta32,eta41,eta42,eta34,ometa34,Qup,Q_lep)
       intsub_z_db  = gqb_qba_sub_qqbz12(mu,z,Ec,E3,E4,eta31,eta32,eta41,eta42,eta34,ometa34,-Qdn,Q_lep)
       intsub_z_ub  = gqb_qba_sub_qqbz12(mu,z,Ec,E3,E4,eta31,eta32,eta41,eta42,eta34,ometa34,-Qup,Q_lep)
       if (fncheck) then
          print *, "eta31->",eta31,", eta41->",eta41,", eta32->",eta32,", eta42->",eta42,", eta34->",eta34,", mu->",mu, ", E3->",E3,", E4->",E4, ", Ec->",Ec,",z->",z,",zb->",zbar
          print *, "z12 u", (intsub_zzb_pls_u + intsub_z_u)   
          print *, "z12 d", (intsub_zzb_pls_d + intsub_z_d)   
          print *, "z12 ub",( intsub_zzb_pls_ub + intsub_z_ub)
          print *, "z12 db",( intsub_zzb_pls_db + intsub_z_db)
          pause
       endif

       res_tmp(1,2) = res_lo(1,2)*(intsub_zzb_pls_u + intsub_z_u)
       res_tmp(1,1) = res_lo(1,1)*(intsub_zzb_pls_d + intsub_z_d)

       res_tmp(2,2) = res_lo(2,2)*( intsub_zzb_pls_ub + intsub_z_ub)
       res_tmp(2,1) = res_lo(2,1)*( intsub_zzb_pls_db + intsub_z_db)

!       print *, "12 coeff",intsub_zzb_u(4), - intsub_z_u(1), - intsub_zb_u(1)
       
       call get_respdf(gq_lumi,1,1,LOProc_z12,res_tmp,respdf)

       ! -- int. subtraction term
       respdf = respdf*LOProc_z12%wgt 
       kin(1) = respdf(1)
       call fill_histo(respdf,vegasweight)

    endif


    !! FLM[z1,zb2]
    call cut_histo(LOProc_zzb)
    if (LOProc_zzb%makecut.or.LOProc_zzb%flag) then

       kin(2) = zero
       FintLO_s(2) = zero

    else

       call res_tree_qqb(LOProc_zzb%AmpMom,res_lo)
       
       mu = LOProc_zzb%muf(1)

       Ec = LOProc_zzb%Lim_Ei(1)    ! = E1

       ! get the relevant int. subtr. functions
       intsub_zzb_u  = gqb_qba_sub_qqbz1zbar2(mu,z,zbar,Ec,Qup)
       intsub_zzb_d  = gqb_qba_sub_qqbz1zbar2(mu,z,zbar,Ec,Qdn)
       intsub_zzb_ub = intsub_zzb_u                  ! ~ Qq^2
       intsub_zzb_db = intsub_zzb_d                  ! ~ Qq^2
       if (fncheck) then
          print *,"mu->",mu,",EC->",Ec,",z->",z,",zb->",zbar
          print *, "zzb u",intsub_zzb_u
          print *, "zzb d",intsub_zzb_d
          print *, "zzb ub",intsub_zzb_ub
          print *, "zzb db",intsub_zzb_db
          pause
       endif

       res_tmp(1,2) = res_lo(1,2)*(intsub_zzb_u )
       res_tmp(1,1) = res_lo(1,1)*(intsub_zzb_d )

       res_tmp(2,2) = res_lo(2,2)*( intsub_zzb_ub )
       res_tmp(2,1) = res_lo(2,1)*( intsub_zzb_db )

       call get_respdf(gq_lumi,1,1,LOProc_zzb,res_tmp,respdf)

       ! -- int. subtraction term
       respdf = respdf*LOProc_zzb%wgt 
       kin(2) = respdf(1)
       call fill_histo(respdf,vegasweight)

    endif


    ff(1) = sum(kin(1:2))

    if (ff(1) .ne. ff(1)) then
       print *, "xx",xx
       print *, "z,zbar",z,zbar
       print *, "kin",kin
       stop
    endif
    
!!    if (ff(1) .ne. zero) then
!!       print *, "xx",xx,z,zbar
!!       print *, "kin", kin(1:2)
!!       print *, "ff",ff(1)
!!       stop
!!    endif


!!    print *, "z,zbar",z,zbar
!!    print *, "ff",ff(1)
!!    print *, "kin",kin(1:4,1)


    call close_histo()

!    call check_ff(ff,kLO_max_full+1,3,(/xx,z/),FintNLO_s)
    
#if(_withchecks == 1)
!    FintNLO_qcd_vs = FintNLO_s
#endif
    
  end function xsect_nnlo_sub12_gq_raoul






  function xsect_nnlo_subvewk_gq_raoul(yRnd,ff,vegasweight)
    use mod_subtrfn_nlo_z
    integer :: xsect_nnlo_subvewk_gq_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--                                                                                                                                                                                                      
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintNLO_s(1),kin(1)
    real(dp) :: res_tree(2,3), res_loop(2,3)
    real(dp) :: respdf_gq(ipdf)
    real(dp) :: z,xis(2),q2,intsub_gq(-1:1,nintsub)

    xsect_nnlo_subvewk_gq_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    z=buff+onet*real(yRnd(kLO_max_full),dp)

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

    call kinematics_lo(xx,LOProc)

    LOProc%ids(1:4) = [0,0,id_el,-id_el]
  
    ! save the unboosted xi1 and xi2
    xis(1:2) = LOProc%PartFrac(1:2)
    q2 = two*scr(LOProc%AmpMom(1:4,1),LOProc%AmpMom(1:4,2))

    call cut_histo(LOProc)
    if (LOProc%makecut.or.LOProc%flag) then

       kin = zero
       FintNLO_s = zero

    else

       call res_ewkloop_qqb(LOProc%AmpMom,res_tree,res_loop)

       ! get the relevant int. subtr. functions
       intsub_gq = sub_g_gq_qqb_z(q2,LOProc%muf**2,z)

       !-- [z,2]
       ! change xi1 and recompute pdfs

       LOProc%PartFrac(1) = xis(1)/z
       LOProc%PartFrac(2) = xis(2)

       call get_respdf(gq_lumi_splitb,1,1,LOProc,res_loop,respdf_gq)   ! use splitb lumi as in RV
       respdf_gq = respdf_gq*LOProc%wgt/z * intsub_gq(0,1)

       call fill_histo(respdf_gq,vegasweight)
       kin(1) = respdf_gq(1)

    endif


    ff(1) = kin(1)

!    if (ff(1) .ne. zero) then
!       print *, "xx",xx,z
!       print *, "kin", kin(1:3,1)
!       print *, "ff",ff(1)
!       stop                                                                                                                                                                                                 
!    endif                                                                                                                                                                                                   



    call close_histo()

!    call check_ff(ff,kLO_max_full+1,3,(/xx,z/),FintNLO_s)                                                                                                                                                   

#if(_withchecks == 1)
!    FintNLO_qcd_vs = FintNLO_s
#endif

  end function xsect_nnlo_subvewk_gq_raoul



  
end module mod_xsects_nnlo_s_gq_raoul
  

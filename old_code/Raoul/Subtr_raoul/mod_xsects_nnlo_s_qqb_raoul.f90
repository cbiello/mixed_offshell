module mod_xsects_nnlo_s_qqb_raoul
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
  use mod_subtrfn_nnlo_z_qqb

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
  real(dp), public, save :: FintNNLO_onloqcd_ns(nFint)
  real(dp), public, save :: FintNNLO_onloewk_is_ns(nFint)
  real(dp), public, save :: FintNNLO_onloewk_5i_ns(nFint)
#endif

  private
  public :: xsect_nnlo_sub12_ns_ga_raoul
  public :: xsect_nnlo_onloqcd_ns_ga_raoul
  public :: xsect_nnlo_onloewk_is_ns_ga_raoul
  public :: xsect_nnlo_onloewk_fs_53_ns_ga_raoul, xsect_nnlo_onloewk_fs_54_ns_ga_raoul
  public :: xsect_nnlo_subvqcd_ns_ga_raoul, xsect_nnlo_subvewk_ns_ga_raoul
  public :: xsect_nnlo_sub12_ns_qqb_raoul, xsect_nnlo_sub12_ns_qq_raoul


contains


  !--------------------------------------------------
  !-- subtraction counterterms
  !--------------------------------------------------


  function xsect_nnlo_onloqcd_ns_ga_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_onloqcd_ns_ga_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc_12,SLim_12,C1Lim_12,C2Lim_12,SC1Lim_12,SC2Lim_12
    type(KinConfig) :: HardProc_z12,SLim_z12,C1Lim_z12,C2Lim_z12,SC1Lim_z12,SC2Lim_z12
    type(KinConfig) :: HardProc_1z2,SLim_1z2,C1Lim_1z2,C2Lim_1z2,SC1Lim_1z2,SC2Lim_1z2
    !--                                                                                                                                                                                                      
    real(dp)    :: xx(kNLO_max_full),z
    real(dp)    :: respdf(ipdf),kin(nkin)
    real(dp)    :: res_lo(2,2),res_lo_tmp(2,2),res_nlo(2,2),res_nlo_tmp(2,2)
    real(dp)    :: Ec, e5, z5i, s5i,eik,eta51,eta52, E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34
    real(dp)    :: subint_z124_hard_u,subint_z124_hard_d,subint_z124_slim_u,subint_z124_slim_d
    real(dp)    :: subint_z124_c1lim_u,subint_z124_c1lim_d,subint_z124_c2lim_u,subint_z124_c2lim_d
    real(dp)    :: subint_z124_sc1lim_u,subint_z124_sc1lim_d,subint_z124_sc2lim_u,subint_z124_sc2lim_d
    real(dp)    :: subint_1z24_hard_u,subint_1z24_hard_d,subint_1z24_slim_u,subint_1z24_slim_d
    real(dp)    :: subint_1z24_c1lim_u,subint_1z24_c1lim_d,subint_1z24_c2lim_u,subint_1z24_c2lim_d
    real(dp)    :: subint_1z24_sc1lim_u,subint_1z24_sc1lim_d,subint_1z24_sc2lim_u,subint_1z24_sc2lim_d
    real(dp)    :: subint_124_hard_u,subint_124_hard_d,subint_124_slim_u,subint_124_slim_d
    real(dp)    :: subint_124_c1lim_u,subint_124_c1lim_d,subint_124_c2lim_u,subint_124_c2lim_d
    real(dp)    :: subint_124_sc1lim_u,subint_124_sc1lim_d,subint_124_sc2lim_u,subint_124_sc2lim_d

    real(dp)    :: subint_124_hard_ub,subint_124_hard_db,subint_124_slim_ub,subint_124_slim_db
    real(dp)    :: subint_124_c1lim_ub,subint_124_c1lim_db,subint_124_c2lim_ub,subint_124_c2lim_db
    real(dp)    :: subint_124_sc1lim_ub,subint_124_sc1lim_db,subint_124_sc2lim_ub,subint_124_sc2lim_db

    xsect_nnlo_onloqcd_ns_ga_raoul = 0

    ff(1) = zero
    kin = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z=buff+onet*real(yRnd(kNLO_max_full),dp)

!!#if (_withchecks == 1)
!!    if (override) then
!!       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
!!       print *, 'overriding input'
!!    endif
!!#endif

    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_r) .or. (one-z) .lt. buff_z) then
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
       FintNNLO_onloqcd_ns(1) = zero

    else

       call res_tree_g_qqb(HardProc_z12%AmpMom,res_nlo)
       Ec = HardProc_z12%Lim_Ei(1)                      ! Ec = E1 = E2


       subint_z124_hard_u = qqb_ga_sub_qqbg_zi(Ec,HardProc_z12%muf(1),z,HardProc_z12%Lim_etaij(1,5),HardProc_z12%Lim_etaij(1,5),Qup)
       subint_z124_hard_d = qqb_ga_sub_qqbg_zi(Ec,HardProc_z12%muf(1),z,HardProc_z12%Lim_etaij(1,5),HardProc_z12%Lim_etaij(1,5),Qdn)

       if (fncheck) then
          print *, "EC->",Ec,",mu->",HardProc_z12%muf(1),",eta51->",HardProc_z12%Lim_etaij(1,5),",z->",z
          print *,"z125 u",subint_z124_hard_u
          print *,"z125 d",subint_z124_hard_d
          pause
       endif

       res_nlo_tmp(:,1) = res_nlo(:,1) * subint_z124_hard_d
       res_nlo_tmp(:,2) = res_nlo(:,2) * subint_z124_hard_u
       call get_respdf(ns_lumi,1,1,HardProc_z12,res_nlo_tmp,respdf)

       respdf = respdf*HardProc_z12%wgt

       kin(1) = respdf(1)
       FintNNLO_onloqcd_ns(1) = kin(1)
       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 1', FintNNLO_onloqcd_ns(1)
       
    endif

    !-- C1

    call cut_histo(C1Lim_z12)
    if (C1Lim_z12%makecut.or.C1Lim_z12%flag) then

       kin(2) = zero
       FintNNLO_onloqcd_ns(2) = zero

    else

       call res_tree_qqb(C1Lim_z12%AmpMom,res_lo)

       Ec = C1Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       subint_z124_c1lim_u = qqb_ga_sub_qqbg_zi(Ec,C1Lim_z12%muf(1),z,zero,C1Lim_z12%Lim_etaij(1,5),Qup)
       subint_z124_c1lim_d = qqb_ga_sub_qqbg_zi(Ec,C1Lim_z12%muf(1),z,zero,C1Lim_z12%Lim_etaij(1,5),Qdn)

       res_lo_tmp(:,1) = res_lo(:,1) * subint_z124_c1lim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_z124_c1lim_u


       call get_respdf(ns_lumi,1,1,C1Lim_z12,res_lo_tmp,respdf)

       z5i = C1Lim_z12%Lim_KinInv(1)
       s5i = C1Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z5i)/(one-z5i))&
            * C1Lim_z12%wgt

       FintNNLO_onloqcd_ns(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 2', FintNNLO_onloqcd_ns(2)
       
    endif

    !-- C2

    call cut_histo(C2Lim_z12)
    if (C2Lim_z12%makecut.or.C2Lim_z12%flag) then

       kin(3) = zero
       FintNNLO_onloqcd_ns(2) = zero

    else
       call res_tree_qqb(C2Lim_z12%AmpMom,res_lo)

       Ec = C2Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       subint_z124_c2lim_u = qqb_ga_sub_qqbg_zi(Ec,C2Lim_z12%muf(1),z,one,one,Qup)
       subint_z124_c2lim_d = qqb_ga_sub_qqbg_zi(Ec,C2Lim_z12%muf(1),z,one,one,Qdn)

       res_lo_tmp(:,1) = res_lo(:,1) * subint_z124_c2lim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_z124_c2lim_u


       call get_respdf(ns_lumi,1,1,C2Lim_z12,res_lo_tmp,respdf)

       z5i = C2Lim_z12%Lim_KinInv(1)
       s5i = C2Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z5i)/(one-z5i))&
            * C2Lim_z12%wgt

       FintNNLO_onloqcd_ns(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 3', FintNNLO_onloqcd_ns(3)

    endif

    !-- S

    call cut_histo(SLim_z12)
    if (SLim_z12%makecut.or.SLim_z12%flag) then

       kin(4) = zero
       FintNNLO_onloqcd_ns(4) = zero

    else

       call res_tree_qqb(SLim_z12%AmpMom,res_lo)
       call get_qcd_eik(Cf,SLim_z12%Lim_etaij(1:5,1:5),[1,2],5,eik)

       Ec = SLim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       subint_z124_slim_u = qqb_ga_sub_qqbg_zi(Ec,SLim_z12%muf(1),z,SLim_z12%Lim_etaij(1,5),SLim_z12%Lim_etaij(1,5),Qup)
       subint_z124_slim_d = qqb_ga_sub_qqbg_zi(Ec,SLim_z12%muf(1),z,SLim_z12%Lim_etaij(1,5),SLim_z12%Lim_etaij(1,5),Qdn)


       res_lo_tmp(:,1) = res_lo(:,1) * subint_z124_slim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_z124_slim_u
       call get_respdf(ns_lumi,1,1,SLim_z12,res_lo_tmp,respdf)

       e5 = SLim_z12%Lim_Ei(5)

       respdf = -respdf*SLim_z12%wgt*eik/e5**2



       kin(4) = respdf(1)
       FintNNLO_onloqcd_ns(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 4', FintNNLO_onloqcd_ns(4)

    endif



    !-- SC1

    call cut_histo(SC1Lim_z12)
    if (SC1Lim_z12%makecut.or.SC1Lim_z12%flag) then

       kin(5) = zero
       FintNNLO_onloqcd_ns(5) = zero

    else

       call res_tree_qqb(SC1Lim_z12%AmpMom,res_lo)

       Ec = SC1Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       subint_z124_sc1lim_u = qqb_ga_sub_qqbg_zi(Ec,SC1Lim_z12%muf(1),z,zero,SC1Lim_z12%Lim_etaij(1,5),Qup)
       subint_z124_sc1lim_d = qqb_ga_sub_qqbg_zi(Ec,SC1Lim_z12%muf(1),z,zero,SC1Lim_z12%Lim_etaij(1,5),Qdn)

       res_lo_tmp(:,1) = res_lo(:,1) * subint_z124_sc1lim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_z124_sc1lim_u


       call get_respdf(ns_lumi,1,1,SC1Lim_z12,res_lo_tmp,respdf)

       e5 = SC1Lim_z12%Lim_KinInv(1)
       eta51 = SC1Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = respdf*Cf/e5**2/eta51* SC1Lim_z12%wgt

       FintNNLO_onloqcd_ns(5) = respdf(1)
       kin(5) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 5', FintNNLO_onloqcd_ns(5)

    endif

    !-- SC2

    call cut_histo(SC2Lim_z12)
    if (SC2Lim_z12%makecut.or.SC2Lim_z12%flag) then

       kin(6) = zero
       FintNNLO_onloqcd_ns(6) = zero

    else

       call res_tree_qqb(SC2Lim_z12%AmpMom,res_lo)

       Ec = SC2Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       subint_z124_sc2lim_u = qqb_ga_sub_qqbg_zi(Ec,SC2Lim_z12%muf(1),z,one,one,Qup)
       subint_z124_sc2lim_d = qqb_ga_sub_qqbg_zi(Ec,SC2Lim_z12%muf(1),z,one,one,Qdn)

       res_lo_tmp(:,1) = res_lo(:,1) * subint_z124_sc2lim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_z124_sc2lim_u


       call get_respdf(ns_lumi,1,1,SC2Lim_z12,res_lo_tmp,respdf)

       e5 = SC2Lim_z12%Lim_KinInv(1)
       eta52 = SC2Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]

       respdf = respdf*Cf/e5**2/eta52* SC2Lim_z12%wgt

       FintNNLO_onloqcd_ns(6) = respdf(1)
       kin(6) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 6', FintNNLO_onloqcd_ns(6)
       !print *, '[]'

    endif




    !! ----------------------------- !!
    !! FLM[1,z2,3,4 | 5g]            !!
    !! ----------------------------- !!


    call kinematics_nlo_z_is(yr=xx,z1=one,z2=z,HardProc=HardProc_1z2,C1Lim=C1Lim_1z2,C2Lim=C2Lim_1z2,SLim=SLim_1z2,SC1Lim=SC1Lim_1z2,SC2Lim=SC2Lim_1z2,compute_etas=.true.)
    HardProc_1z2%ids(1:5) = [0,0,id_el,-id_el,id_g]
    C1Lim_1z2%ids(1:4)    = [0,0,id_el,-id_el]
    SLim_1z2%ids(1:4)     = [0,0,id_el,-id_el]
    C2Lim_1z2%ids(1:4)    = [0,0,id_el,-id_el]
    SC1Lim_1z2%ids(1:4)   = [0,0,id_el,-id_el]
    sC2Lim_1z2%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard

    call cut_histo(HardProc_1z2)
    if (HardProc_1z2%makecut.or.HardProc_1z2%flag) then

       kin(7) = zero
       FintNNLO_onloqcd_ns(7) = zero

    else

       call res_tree_g_qqb(HardProc_1z2%AmpMom,res_nlo)
       Ec = HardProc_1z2%Lim_Ei(1)                      ! Ec = E1 = E2

       !print *, "Lmu",log(HardProc_1z2%muf(1)/two/Ec)
       !print *, "shat",HardProc_1z2%PartFrac(1)*HardProc_1z2%PartFrac(2)*sh
       !print *, "E1,E2",HardProc_1z2%lim_Ei(1),HardProc_1z2%lim_Ei(2)
       !print *, "Lmu?",half*log(HardProc_1z2%muf(1)**2/(HardProc_1z2%PartFrac(1)*HardProc_1z2%PartFrac(2)*sh))
       !pause



       subint_1z24_hard_u = qqb_ga_sub_qqbg_zi(Ec,HardProc_1z2%muf(1),z,HardProc_1z2%Lim_etaij(2,5),HardProc_1z2%Lim_etaij(2,5),Qup)
       subint_1z24_hard_d = qqb_ga_sub_qqbg_zi(Ec,HardProc_1z2%muf(1),z,HardProc_1z2%Lim_etaij(2,5),HardProc_1z2%Lim_etaij(2,5),Qdn)

       res_nlo_tmp(:,1) = res_nlo(:,1) * subint_1z24_hard_d
       res_nlo_tmp(:,2) = res_nlo(:,2) * subint_1z24_hard_u
       call get_respdf(ns_lumi,1,1,HardProc_1z2,res_nlo_tmp,respdf)

       respdf = respdf*HardProc_1z2%wgt

       kin(7) = respdf(1)
       FintNNLO_onloqcd_ns(7) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 7', FintNNLO_onloqcd_ns(7)

    endif

    !-- C1

    call cut_histo(C1Lim_1z2)
    if (C1Lim_1z2%makecut.or.C1Lim_1z2%flag) then

       kin(8) = zero
       FintNNLO_onloqcd_ns(8) = zero

    else

       call res_tree_qqb(C1Lim_1z2%AmpMom,res_lo)

       Ec = C1Lim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       subint_1z24_c1lim_u = qqb_ga_sub_qqbg_zi(Ec,C1Lim_1z2%muf(1),z,one,one,Qup)
       subint_1z24_c1lim_d = qqb_ga_sub_qqbg_zi(Ec,C1Lim_1z2%muf(1),z,one,one,Qdn)

       res_lo_tmp(:,1) = res_lo(:,1) * subint_1z24_c1lim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_1z24_c1lim_u


       call get_respdf(ns_lumi,1,1,C1Lim_1z2,res_lo_tmp,respdf)

       z5i = C1Lim_1z2%Lim_KinInv(1)
       s5i = C1Lim_1z2%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z5i)/(one-z5i))&
            * C1Lim_1z2%wgt

       FintNNLO_onloqcd_ns(8) = respdf(1)
       kin(8) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 8', FintNNLO_onloqcd_ns(8)

    endif

    !-- C2

    call cut_histo(C2Lim_1z2)
    if (C2Lim_1z2%makecut.or.C2Lim_1z2%flag) then

       kin(9) = zero
       FintNNLO_onloqcd_ns(9) = zero

    else
       call res_tree_qqb(C2Lim_1z2%AmpMom,res_lo)

       Ec = C2Lim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       subint_1z24_c2lim_u = qqb_ga_sub_qqbg_zi(Ec,C2Lim_1z2%muf(1),z,zero,C2Lim_1z2%Lim_etaij(2,5),Qup)
       subint_1z24_c2lim_d = qqb_ga_sub_qqbg_zi(Ec,C2Lim_1z2%muf(1),z,zero,C2Lim_1z2%Lim_etaij(2,5),Qdn)

       res_lo_tmp(:,1) = res_lo(:,1) * subint_1z24_c2lim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_1z24_c2lim_u


       call get_respdf(ns_lumi,1,1,C2Lim_1z2,res_lo_tmp,respdf)

       z5i = C2Lim_1z2%Lim_KinInv(1)
       s5i = C2Lim_1z2%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z5i)/(one-z5i))&
            * C2Lim_1z2%wgt

       FintNNLO_onloqcd_ns(9) = respdf(1)
       kin(9) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 9', FintNNLO_onloqcd_ns(9)

    endif

    !-- S

    call cut_histo(SLim_1z2)
    if (SLim_1z2%makecut.or.SLim_1z2%flag) then

       kin(10) = zero
       FintNNLO_onloqcd_ns(10) = zero

    else

       call res_tree_qqb(SLim_1z2%AmpMom,res_lo)
       call get_qcd_eik(Cf,SLim_1z2%Lim_etaij(1:5,1:5),[1,2],5,eik)

       Ec = SLim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2

       subint_1z24_slim_u = qqb_ga_sub_qqbg_zi(Ec,SLim_1z2%muf(1),z,SLim_1z2%Lim_etaij(2,5),SLim_1z2%Lim_etaij(2,5),Qup)
       subint_1z24_slim_d = qqb_ga_sub_qqbg_zi(Ec,SLim_1z2%muf(1),z,SLim_1z2%Lim_etaij(2,5),SLim_1z2%Lim_etaij(2,5),Qdn)


       res_lo_tmp(:,1) = res_lo(:,1) * subint_1z24_slim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_1z24_slim_u
       call get_respdf(ns_lumi,1,1,SLim_1z2,res_lo_tmp,respdf)

       e5 = SLim_1z2%Lim_Ei(5)

       respdf = -respdf*SLim_1z2%wgt*eik/e5**2



       kin(10) = respdf(1)
       FintNNLO_onloqcd_ns(10) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 10', FintNNLO_onloqcd_ns(10)

    endif



    !-- SC1

    call cut_histo(SC1Lim_1z2)
    if (SC1Lim_1z2%makecut.or.SC1Lim_1z2%flag) then

       kin(11) = zero
       FintNNLO_onloqcd_ns(11) = zero

    else

       call res_tree_qqb(SC1Lim_1z2%AmpMom,res_lo)

       Ec = SC1Lim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       subint_1z24_sc1lim_u = qqb_ga_sub_qqbg_zi(Ec,SC1Lim_1z2%muf(1),z,one,one,Qup)
       subint_1z24_sc1lim_d = qqb_ga_sub_qqbg_zi(Ec,SC1Lim_1z2%muf(1),z,one,one,Qdn)

       res_lo_tmp(:,1) = res_lo(:,1) * subint_1z24_sc1lim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_1z24_sc1lim_u


       call get_respdf(ns_lumi,1,1,SC1Lim_1z2,res_lo_tmp,respdf)

       e5 = SC1Lim_1z2%Lim_KinInv(1)
       eta51 = SC1Lim_1z2%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = respdf*Cf/e5**2/eta51* SC1Lim_1z2%wgt

       FintNNLO_onloqcd_ns(11) = respdf(1)
       kin(11) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 11', FintNNLO_onloqcd_ns(11)

    endif

    !-- SC2

    call cut_histo(SC2Lim_1z2)
    if (SC2Lim_1z2%makecut.or.SC2Lim_1z2%flag) then

       kin(12) = zero
       FintNNLO_onloqcd_ns(12) = zero

    else

       call res_tree_qqb(SC2Lim_1z2%AmpMom,res_lo)

       Ec = SC2Lim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       subint_1z24_sc2lim_u = qqb_ga_sub_qqbg_zi(Ec,SC2Lim_1z2%muf(1),z,zero,SC2Lim_1z2%Lim_etaij(2,5),Qup)
       subint_1z24_sc2lim_d = qqb_ga_sub_qqbg_zi(Ec,SC2Lim_1z2%muf(1),z,zero,SC2Lim_1z2%Lim_etaij(2,5),Qdn)

       res_lo_tmp(:,1) = res_lo(:,1) * subint_1z24_sc2lim_d
       res_lo_tmp(:,2) = res_lo(:,2) * subint_1z24_sc2lim_u


       call get_respdf(ns_lumi,1,1,SC2Lim_1z2,res_lo_tmp,respdf)

       e5 = SC2Lim_1z2%Lim_KinInv(1)
       eta52 = SC2Lim_1z2%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]

       respdf = respdf*Cf/e5**2/eta52* SC2Lim_1z2%wgt

       FintNNLO_onloqcd_ns(12) = respdf(1)
       kin(12) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 12', FintNNLO_onloqcd_ns(12)
       !print *, '[]'

    endif




    !! ----------------------------- !!
    !! FLM[1,2,3,4 | 5g]            !!
    !! ----------------------------- !!



    call kinematics_nlo_z_is(yr=xx,z1=one,z2=one,HardProc=HardProc_12,C1Lim=C1Lim_12,C2Lim=C2Lim_12,SLim=SLim_12,SC1Lim=SC1Lim_12,SC2Lim=SC2Lim_12,compute_etas=.true.)
    HardProc_12%ids(1:5) = [0,0,id_el,-id_el,id_g]
    C1Lim_12%ids(1:4)    = [0,0,id_el,-id_el]
    SLim_12%ids(1:4)     = [0,0,id_el,-id_el]
    C2Lim_12%ids(1:4)    = [0,0,id_el,-id_el]
    SC1Lim_12%ids(1:4)   = [0,0,id_el,-id_el]
    sC2Lim_12%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard

    call cut_histo(HardProc_12)
    if (HardProc_12%makecut.or.HardProc_12%flag) then

       kin(13) = zero
       FintNNLO_onloqcd_ns(13) = zero

    else

       call res_tree_g_qqb(HardProc_12%AmpMom,res_nlo)
       Ec = HardProc_12%Lim_Ei(1)                      ! Ec = E1 = E2

       E1    = HardProc_12%Lim_Ei(1)
       E2    = HardProc_12%Lim_Ei(2)
       E3    = HardProc_12%Lim_Ei(3)
       E4    = HardProc_12%Lim_Ei(4)
       eta31 = HardProc_12%Lim_etaij(1,3)
       eta41 = HardProc_12%Lim_etaij(1,4)
       eta32 = HardProc_12%Lim_etaij(2,3)
       eta42 = HardProc_12%Lim_etaij(2,4)
       eta34 = HardProc_12%Lim_etaij(3,4)

       subint_124_hard_u = qqb_ga_sub_qqbg(Ec,HardProc_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qup,Q_lep)
       subint_124_hard_d = qqb_ga_sub_qqbg(Ec,HardProc_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qdn,Q_lep)

       subint_124_hard_ub = qqb_ga_sub_qqbg(Ec,HardProc_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qup,Q_lep)
       subint_124_hard_db = qqb_ga_sub_qqbg(Ec,HardProc_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qdn,Q_lep)

       subint_124_hard_u = subint_124_hard_u - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,HardProc_12%muf(1),z,HardProc_12%Lim_etaij(1,5),HardProc_12%Lim_etaij(1,5),Qup) + & 
              qqb_ga_sub_qqbg_zi_pls(Ec,HardProc_12%muf(1),z,HardProc_12%Lim_etaij(2,5),HardProc_12%Lim_etaij(2,5),Qup))
       subint_124_hard_d = subint_124_hard_d -&
            ( qqb_ga_sub_qqbg_zi_pls(Ec,HardProc_12%muf(1),z,HardProc_12%Lim_etaij(1,5),HardProc_12%Lim_etaij(1,5),Qdn) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,HardProc_12%muf(1),z,HardProc_12%Lim_etaij(2,5),HardProc_12%Lim_etaij(2,5),Qdn) )

       subint_124_hard_ub = subint_124_hard_ub - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,HardProc_12%muf(1),z,HardProc_12%Lim_etaij(1,5),HardProc_12%Lim_etaij(1,5),-Qup) + & 
              qqb_ga_sub_qqbg_zi_pls(Ec,HardProc_12%muf(1),z,HardProc_12%Lim_etaij(2,5),HardProc_12%Lim_etaij(2,5),-Qup))
       subint_124_hard_db = subint_124_hard_db -&
            ( qqb_ga_sub_qqbg_zi_pls(Ec,HardProc_12%muf(1),z,HardProc_12%Lim_etaij(1,5),HardProc_12%Lim_etaij(1,5),-Qdn) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,HardProc_12%muf(1),z,HardProc_12%Lim_etaij(2,5),HardProc_12%Lim_etaij(2,5),-Qdn) )

       if (fncheck) then
          print *, "E1->",E1,",E2->",E2,",E3->",E3,",E4->",E4,",EC->",Ec,",mu->",HardProc_12%muf(1),"eta13->",eta31,",eta23->",eta32,",eta14->",eta41,",eta24->",eta42,"eta34->",eta34,",eta51->",HardProc_12%Lim_etaij(1,5),",z->",z
          print *,"125 u",subint_124_hard_u
          print *,"125 d",subint_124_hard_d
          print *,"125 ub",subint_124_hard_ub
          print *,"125 db",subint_124_hard_db
          pause
       endif

       res_nlo_tmp(1,1) = res_nlo(1,1) * subint_124_hard_d
       res_nlo_tmp(1,2) = res_nlo(1,2) * subint_124_hard_u

       res_nlo_tmp(2,1) = res_nlo(2,1) * subint_124_hard_db
       res_nlo_tmp(2,2) = res_nlo(2,2) * subint_124_hard_ub


       call get_respdf(ns_lumi,1,1,HardProc_12,res_nlo_tmp,respdf)

       respdf = respdf*HardProc_12%wgt

       kin(13) = respdf(1)
       FintNNLO_onloqcd_ns(13) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 13', FintNNLO_onloqcd_ns(13)

    endif

    !-- C1

    call cut_histo(C1Lim_12)
    if (C1Lim_12%makecut.or.C1Lim_12%flag) then

       kin(14) = zero
       FintNNLO_onloqcd_ns(14) = zero

    else

       call res_tree_qqb(C1Lim_12%AmpMom,res_lo)

       Ec = C1Lim_12%Lim_Ei(1)                      ! Ec = E1 = E2

       E1    = C1Lim_12%Lim_Ei(1)
       E2    = C1Lim_12%Lim_Ei(2)
       E3    = C1Lim_12%Lim_Ei(3)
       E4    = C1Lim_12%Lim_Ei(4)
       eta31 = C1Lim_12%Lim_etaij(1,3)
       eta41 = C1Lim_12%Lim_etaij(1,4)
       eta32 = C1Lim_12%Lim_etaij(2,3)
       eta42 = C1Lim_12%Lim_etaij(2,4)
       eta34 = C1Lim_12%Lim_etaij(3,4)

       subint_124_c1lim_u = qqb_ga_sub_qqbg(Ec,C1Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qup,Q_lep)
       subint_124_c1lim_d = qqb_ga_sub_qqbg(Ec,C1Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qdn,Q_lep)

       subint_124_c1lim_ub = qqb_ga_sub_qqbg(Ec,C1Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qup,Q_lep)
       subint_124_c1lim_db = qqb_ga_sub_qqbg(Ec,C1Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qdn,Q_lep)


       subint_124_c1lim_u = subint_124_c1lim_u - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,C1Lim_12%muf(1),z,zero,C1Lim_12%Lim_etaij(1,5),Qup) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,C1Lim_12%muf(1),z,one, one,Qup) )
       subint_124_c1lim_d = subint_124_c1lim_d - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,C1Lim_12%muf(1),z,zero,C1Lim_12%Lim_etaij(1,5),Qdn) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,C1Lim_12%muf(1),z,one, one, Qdn) )

       subint_124_c1lim_ub = subint_124_c1lim_ub - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,C1Lim_12%muf(1),z,zero,C1Lim_12%Lim_etaij(1,5),-Qup) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,C1Lim_12%muf(1),z,one, one, -Qup) )
       subint_124_c1lim_db = subint_124_c1lim_db - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,C1Lim_12%muf(1),z,zero,C1Lim_12%Lim_etaij(1,5),-Qdn) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,C1Lim_12%muf(1),z,one, one, -Qdn) )

       res_lo_tmp(1,1) = res_lo(1,1) * subint_124_c1lim_d
       res_lo_tmp(1,2) = res_lo(1,2) * subint_124_c1lim_u


       res_lo_tmp(2,1) = res_lo(2,1) * subint_124_c1lim_db
       res_lo_tmp(2,2) = res_lo(2,2) * subint_124_c1lim_ub


       call get_respdf(ns_lumi,1,1,C1Lim_12,res_lo_tmp,respdf)

       z5i = C1Lim_12%Lim_KinInv(1)
       s5i = C1Lim_12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z5i)/(one-z5i))&
            * C1Lim_12%wgt

       FintNNLO_onloqcd_ns(14) = respdf(1)
       kin(14) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 14', FintNNLO_onloqcd_ns(14)

    endif

    !-- C2

    call cut_histo(C2Lim_12)
    if (C2Lim_12%makecut.or.C2Lim_12%flag) then

       kin(15) = zero
       FintNNLO_onloqcd_ns(15) = zero

    else
       call res_tree_qqb(C2Lim_12%AmpMom,res_lo)

       Ec = C2Lim_12%Lim_Ei(1)                      ! Ec = E1 = E2

       E1    = C2Lim_12%Lim_Ei(1)
       E2    = C2Lim_12%Lim_Ei(2)
       E3    = C2Lim_12%Lim_Ei(3)
       E4    = C2Lim_12%Lim_Ei(4)
       eta31 = C2Lim_12%Lim_etaij(1,3)
       eta41 = C2Lim_12%Lim_etaij(1,4)
       eta32 = C2Lim_12%Lim_etaij(2,3)
       eta42 = C2Lim_12%Lim_etaij(2,4)
       eta34 = C2Lim_12%Lim_etaij(3,4)

       subint_124_c2lim_u = qqb_ga_sub_qqbg(Ec,C2Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qup,Q_lep)
       subint_124_c2lim_d = qqb_ga_sub_qqbg(Ec,C2Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qdn,Q_lep)

       subint_124_c2lim_ub = qqb_ga_sub_qqbg(Ec,C2Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qup,Q_lep)
       subint_124_c2lim_db = qqb_ga_sub_qqbg(Ec,C2Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qdn,Q_lep)


       subint_124_c2lim_u = subint_124_c2lim_u - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,C2Lim_12%muf(1),z,one, one, Qup) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,C2Lim_12%muf(1),z,zero,C2Lim_12%Lim_etaij(2,5),Qup) )
       subint_124_c2lim_d = subint_124_c2lim_d - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,C2Lim_12%muf(1),z,one, one, Qdn) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,C2Lim_12%muf(1),z,zero,C2Lim_12%Lim_etaij(2,5),Qdn) )

       subint_124_c2lim_ub = subint_124_c2lim_ub - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,C2Lim_12%muf(1),z,one, one,-Qup) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,C2Lim_12%muf(1),z,zero,C2Lim_12%Lim_etaij(2,5),-Qup) )
       subint_124_c2lim_db = subint_124_c2lim_db - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,C2Lim_12%muf(1),z,one, one,-Qdn) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,C2Lim_12%muf(1),z,zero,C2Lim_12%Lim_etaij(2,5),-Qdn) )

       res_lo_tmp(1,1) = res_lo(1,1) * subint_124_c2lim_d
       res_lo_tmp(1,2) = res_lo(1,2) * subint_124_c2lim_u


       res_lo_tmp(2,1) = res_lo(2,1) * subint_124_c2lim_db
       res_lo_tmp(2,2) = res_lo(2,2) * subint_124_c2lim_ub


       call get_respdf(ns_lumi,1,1,C2Lim_12,res_lo_tmp,respdf)

       z5i = C2Lim_12%Lim_KinInv(1)
       s5i = C2Lim_12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z5i)/(one-z5i))&
            * C2Lim_12%wgt

       FintNNLO_onloqcd_ns(15) = respdf(1)
       kin(15) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 15', FintNNLO_onloqcd_ns(15)

    endif

    !-- S

    call cut_histo(SLim_12)
    if (SLim_12%makecut.or.SLim_12%flag) then

       kin(16) = zero
       FintNNLO_onloqcd_ns(16) = zero

    else

       call res_tree_qqb(SLim_12%AmpMom,res_lo)
       call get_qcd_eik(Cf,SLim_12%Lim_etaij(1:5,1:5),[1,2],5,eik)

       E1    = SLim_12%Lim_Ei(1)
       E2    = SLim_12%Lim_Ei(2)
       E3    = SLim_12%Lim_Ei(3)
       E4    = SLim_12%Lim_Ei(4)
       eta31 = SLim_12%Lim_etaij(1,3)
       eta41 = SLim_12%Lim_etaij(1,4)
       eta32 = SLim_12%Lim_etaij(2,3)
       eta42 = SLim_12%Lim_etaij(2,4)
       eta34 = SLim_12%Lim_etaij(3,4)

       Ec = E1

       subint_124_slim_u = qqb_ga_sub_qqbg(Ec,SLim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qup,Q_lep)
       subint_124_slim_d = qqb_ga_sub_qqbg(Ec,SLim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qdn,Q_lep)

       subint_124_slim_ub = qqb_ga_sub_qqbg(Ec,SLim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qup,Q_lep)
       subint_124_slim_db = qqb_ga_sub_qqbg(Ec,SLim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qdn,Q_lep)

       subint_124_slim_u = subint_124_slim_u - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,SLim_12%muf(1),z,SLim_12%Lim_etaij(1,5),SLim_12%Lim_etaij(1,5),Qup)  + &
              qqb_ga_sub_qqbg_zi_pls(Ec,SLim_12%muf(1),z,SLim_12%Lim_etaij(2,5),SLim_12%Lim_etaij(2,5),Qup) )
       subint_124_slim_d = subint_124_slim_d - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,SLim_12%muf(1),z,SLim_12%Lim_etaij(1,5),SLim_12%Lim_etaij(1,5),Qdn) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,SLim_12%muf(1),z,SLim_12%Lim_etaij(2,5),SLim_12%Lim_etaij(2,5),Qdn) )

       subint_124_slim_ub = subint_124_slim_ub - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,SLim_12%muf(1),z,SLim_12%Lim_etaij(1,5),SLim_12%Lim_etaij(1,5),-Qup)  + &
              qqb_ga_sub_qqbg_zi_pls(Ec,SLim_12%muf(1),z,SLim_12%Lim_etaij(2,5),SLim_12%Lim_etaij(2,5),-Qup) )
       subint_124_slim_db = subint_124_slim_db - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,SLim_12%muf(1),z,SLim_12%Lim_etaij(1,5),SLim_12%Lim_etaij(1,5),-Qdn) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,SLim_12%muf(1),z,SLim_12%Lim_etaij(2,5),SLim_12%Lim_etaij(2,5),-Qdn) )


       res_lo_tmp(1,1) = res_lo(1,1) * subint_124_slim_d
       res_lo_tmp(1,2) = res_lo(1,2) * subint_124_slim_u


       res_lo_tmp(2,1) = res_lo(2,1) * subint_124_slim_db
       res_lo_tmp(2,2) = res_lo(2,2) * subint_124_slim_ub

       call get_respdf(ns_lumi,1,1,SLim_12,res_lo_tmp,respdf)

       e5 = SLim_12%Lim_Ei(5)

       respdf = -respdf*SLim_12%wgt*eik/e5**2



       kin(16) = respdf(1)
       FintNNLO_onloqcd_ns(16) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 16', FintNNLO_onloqcd_ns(16)

    endif



    !-- SC1

    call cut_histo(SC1Lim_12)
    if (SC1Lim_12%makecut.or.SC1Lim_12%flag) then

       kin(17) = zero
       FintNNLO_onloqcd_ns(17) = zero

    else

       call res_tree_qqb(SC1Lim_12%AmpMom,res_lo)

       Ec = SC1Lim_12%Lim_Ei(1)                      ! Ec = E1 = E2

       E1    = SC1Lim_12%Lim_Ei(1)
       E2    = SC1Lim_12%Lim_Ei(2)
       E3    = SC1Lim_12%Lim_Ei(3)
       E4    = SC1Lim_12%Lim_Ei(4)
       eta31 = SC1Lim_12%Lim_etaij(1,3)
       eta41 = SC1Lim_12%Lim_etaij(1,4)
       eta32 = SC1Lim_12%Lim_etaij(2,3)
       eta42 = SC1Lim_12%Lim_etaij(2,4)
       eta34 = SC1Lim_12%Lim_etaij(3,4)

       subint_124_sc1lim_u = qqb_ga_sub_qqbg(Ec,SC1Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qup,Q_lep)
       subint_124_sc1lim_d = qqb_ga_sub_qqbg(Ec,SC1Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qdn,Q_lep)

       subint_124_sc1lim_ub = qqb_ga_sub_qqbg(Ec,SC1Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qup,Q_lep)
       subint_124_sc1lim_db = qqb_ga_sub_qqbg(Ec,SC1Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qdn,Q_lep)


       subint_124_sc1lim_u = subint_124_sc1lim_u - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,SC1Lim_12%muf(1),z,zero,SC1Lim_12%Lim_etaij(1,5),Qup) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,SC1Lim_12%muf(1),z,one, one, Qup))
       subint_124_sc1lim_d = subint_124_sc1lim_d - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,SC1Lim_12%muf(1),z,zero,SC1Lim_12%Lim_etaij(1,5),Qdn) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,SC1Lim_12%muf(1),z,one, one,Qdn))

       subint_124_sc1lim_ub = subint_124_sc1lim_ub - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,SC1Lim_12%muf(1),z,zero,SC1Lim_12%Lim_etaij(1,5),-Qup) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,SC1Lim_12%muf(1),z,one, one,-Qup))
       subint_124_sc1lim_db = subint_124_sc1lim_db - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,SC1Lim_12%muf(1),z,zero,SC1Lim_12%Lim_etaij(1,5),-Qdn) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,SC1Lim_12%muf(1),z,one, one, -Qdn))

       res_lo_tmp(1,1) = res_lo(1,1) * subint_124_sc1lim_d
       res_lo_tmp(1,2) = res_lo(1,2) * subint_124_sc1lim_u

       res_lo_tmp(2,1) = res_lo(2,1) * subint_124_sc1lim_db
       res_lo_tmp(2,2) = res_lo(2,2) * subint_124_sc1lim_ub


       call get_respdf(ns_lumi,1,1,SC1Lim_12,res_lo_tmp,respdf)

       e5 = SC1Lim_12%Lim_KinInv(1)
       eta51 = SC1Lim_12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = respdf*Cf/e5**2/eta51* SC1Lim_12%wgt

       FintNNLO_onloqcd_ns(17) = respdf(1)
       kin(17) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 17', FintNNLO_onloqcd_ns(17)

    endif

    !-- SC2

    call cut_histo(SC2Lim_12)
    if (SC2Lim_12%makecut.or.SC2Lim_12%flag) then

       kin(18) = zero
       FintNNLO_onloqcd_ns(18) = zero

    else

       call res_tree_qqb(SC2Lim_12%AmpMom,res_lo)

       Ec = SC2Lim_12%Lim_Ei(1)                      ! Ec = E1 = E2

       E1    = SC2Lim_12%Lim_Ei(1)
       E2    = SC2Lim_12%Lim_Ei(2)
       E3    = SC2Lim_12%Lim_Ei(3)
       E4    = SC2Lim_12%Lim_Ei(4)
       eta31 = SC2Lim_12%Lim_etaij(1,3)
       eta41 = SC2Lim_12%Lim_etaij(1,4)
       eta32 = SC2Lim_12%Lim_etaij(2,3)
       eta42 = SC2Lim_12%Lim_etaij(2,4)
       eta34 = SC2Lim_12%Lim_etaij(3,4)

       subint_124_sc2lim_u = qqb_ga_sub_qqbg(Ec,SC2Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qup,Q_lep)
       subint_124_sc2lim_d = qqb_ga_sub_qqbg(Ec,SC2Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qdn,Q_lep)

       subint_124_sc2lim_ub = qqb_ga_sub_qqbg(Ec,SC2Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qup,Q_lep)
       subint_124_sc2lim_db = qqb_ga_sub_qqbg(Ec,SC2Lim_12%muf(1),E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qdn,Q_lep)


       subint_124_sc2lim_u = subint_124_sc2lim_u - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,SC2Lim_12%muf(1),z,one, one, Qup) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,SC2Lim_12%muf(1),z,zero,SC2Lim_12%Lim_etaij(2,5),Qup))
       subint_124_sc2lim_d = subint_124_sc2lim_d -&
            ( qqb_ga_sub_qqbg_zi_pls(Ec,SC2Lim_12%muf(1),z,one, one,Qdn) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,SC2Lim_12%muf(1),z,zero,SC2Lim_12%Lim_etaij(2,5),Qdn))


       subint_124_sc2lim_ub = subint_124_sc2lim_ub - &
            ( qqb_ga_sub_qqbg_zi_pls(Ec,SC2Lim_12%muf(1),z,one, one, -Qup) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,SC2Lim_12%muf(1),z,zero,SC2Lim_12%Lim_etaij(2,5),-Qup))
       subint_124_sc2lim_db = subint_124_sc2lim_db -&
            ( qqb_ga_sub_qqbg_zi_pls(Ec,SC2Lim_12%muf(1),z,one, one, -Qdn) + &
              qqb_ga_sub_qqbg_zi_pls(Ec,SC2Lim_12%muf(1),z,zero,SC2Lim_12%Lim_etaij(2,5),-Qdn))

       res_lo_tmp(1,1) = res_lo(1,1) * subint_124_sc2lim_d
       res_lo_tmp(1,2) = res_lo(1,2) * subint_124_sc2lim_u


       res_lo_tmp(2,1) = res_lo(2,1) * subint_124_sc2lim_db
       res_lo_tmp(2,2) = res_lo(2,2) * subint_124_sc2lim_ub


       call get_respdf(ns_lumi,1,1,SC2Lim_12,res_lo_tmp,respdf)

       e5 = SC2Lim_12%Lim_KinInv(1)
       eta52 = SC2Lim_12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]

       respdf = respdf*Cf/e5**2/eta52* SC2Lim_12%wgt

       FintNNLO_onloqcd_ns(18) = respdf(1)
       kin(18) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 18', FintNNLO_onloqcd_ns(18)
       !pause

    endif


    if (is_nan(kin)) then
       print *, "kin",kin(1:18)
       stop
    endif


    ff(1) = sum(kin(1:18))


!!    if (kin(1)*kin(2)*kin(3)*kin(4)*kin(5)*kin(6) .ne. zero) then
!!
!!       print *, "kin in onlo z12", kin(1:6)
!!       print *, "kin in onlo 1z2", kin(7:12)
!!       print *, "kin in onlo 12", kin(13:18)
!!       print *, kin(1)+kin(7)+kin(13)
!!       print *, kin(2)+kin(8)+kin(14)
!!       print *, kin(3)+kin(9)+kin(15)
!!       print *, kin(4)+kin(10)+kin(16)
!!       print *, kin(5)+kin(11)+kin(17)
!!       print *, kin(6)+kin(12)+kin(18)
!!       print *, "ff",ff(1)
!!       pause
!!    endif


    call close_histo()


  end function xsect_nnlo_onloqcd_ns_ga_raoul



  function xsect_nnlo_onloewk_is_ns_ga_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_onloewk_is_ns_ga_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc_12,SLim_12,C1Lim_12,C2Lim_12,SC1Lim_12,SC2Lim_12
    type(KinConfig) :: HardProc_z12,SLim_z12,C1Lim_z12,C2Lim_z12,SC1Lim_z12,SC2Lim_z12
    type(KinConfig) :: HardProc_1z2,SLim_1z2,C1Lim_1z2,C2Lim_1z2,SC1Lim_1z2,SC2Lim_1z2
    !--                                                                                                                                                                                                      
    real(dp)    :: xx(kNLO_max_full),z,damp
    real(dp)    :: respdf(ipdf),kin(18)
    real(dp)    :: res_lo(2,2),res_lo_tmp(2,2),res_nlo(2,2)
    real(dp)    :: Ec, e5, z5i, s5i,eik(4),eta51,eta52,eta53,eta54,eta51lim,eta52lim
    real(dp)    :: subint_z124_hard,subint_1z24_hard,subint_124_hard
    real(dp)    :: subint_z124_slim,subint_1z24_slim,subint_124_slim
    real(dp)    :: subint_z124_c1lim,subint_1z24_c1lim,subint_124_c1lim
    real(dp)    :: subint_z124_c2lim,subint_1z24_c2lim,subint_124_c2lim
    real(dp)    :: subint_z124_sc1lim,subint_1z24_sc1lim,subint_124_sc1lim
    real(dp)    :: subint_z124_sc2lim,subint_1z24_sc2lim,subint_124_sc2lim

    xsect_nnlo_onloewk_is_ns_ga_raoul = 0

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

    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_r) .or. (one-z) .lt. buff_z) then
       failed_points = failed_points + 1
       return
    endif
!    xx(xE)   = 1E-8_dp
!    xx(xrho) = 1E-8_dp
!    z = one - 1E-10_dp

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
    sC2Lim_z12%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard

    call cut_histo(HardProc_z12)
    if (HardProc_z12%makecut.or.HardProc_z12%flag) then

       kin(1) = zero
       FintNNLO_onloewk_is_ns(1) = zero

    else

       call res_tree_a_qqb(HardProc_z12%AmpMom,res_nlo)
       call partition_nlo_qed(HardProc_z12,damp,12)

       Ec = HardProc_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = HardProc_z12%Lim_etaij(1,5)
       eta52 = HardProc_z12%Lim_etaij(2,5)
       eta53 = HardProc_z12%Lim_etaij(5,3)
       eta54 = HardProc_z12%Lim_etaij(5,4)
       eta51lim = eta51
       subint_z124_hard = qqb_ga_sub_qqba_zi(Ec,HardProc_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)

       if (fncheck) then
          print *, "mu->",HardProc_z12%muf(1),",z->",z,",eta51->",eta51,"eta52->",eta52,",eta53->",eta53,",eta54->",eta54,",Ec->",Ec
          print *,"subint z124", subint_z124_hard
          pause
       endif

       call get_respdf(ns_lumi,1,1,HardProc_z12,res_nlo,respdf)

       respdf = respdf*HardProc_z12%wgt*subint_z124_hard * damp

       kin(1) = respdf(1)
       FintNNLO_onloewk_is_ns(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1

    call cut_histo(C1Lim_z12)
    if (C1Lim_z12%makecut.or.C1Lim_z12%flag) then

       kin(2) = zero
       FintNNLO_onloewk_is_ns(2) = zero

    else

       call res_tree_qqb(C1Lim_z12%AmpMom,res_lo)

       Ec = C1Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51lim = C1Lim_z12%Lim_etaij(1,5)
       
       eta51 = zero
       eta52 = one
       eta53 = C1Lim_z12%Lim_etaij(1,3)              ! C51 eta53 = eta13
       eta54 = C1Lim_z12%Lim_etaij(1,4)              ! C51 eta54 = eta14         
       subint_z124_c1lim = qqb_ga_sub_qqba_zi(Ec,C1Lim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)


       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       call get_respdf(ns_lumi,1,1,C1Lim_z12,res_lo,respdf)

       z5i = C1Lim_z12%Lim_KinInv(1)
       s5i = C1Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z5i)/(one-z5i))&
            * C1Lim_z12%wgt * subint_z124_c1lim

       FintNNLO_onloewk_is_ns(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C2

    call cut_histo(C2Lim_z12)
    if (C2Lim_z12%makecut.or.C2Lim_z12%flag) then

       kin(3) = zero
       FintNNLO_onloewk_is_ns(3) = zero

    else
       call res_tree_qqb(C2Lim_z12%AmpMom,res_lo)

       Ec = C2Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = one
       eta51lim = eta51
       eta52 = zero
       eta53 = C2Lim_z12%Lim_etaij(2,3)              ! C52 eta53 = eta23
       eta54 = C2Lim_z12%Lim_etaij(2,4)              ! C52 eta54 = eta24         
       subint_z124_c2lim = qqb_ga_sub_qqba_zi(Ec,C2Lim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)

       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       call get_respdf(ns_lumi,1,1,C2Lim_z12,res_lo,respdf)

       z5i = C2Lim_z12%Lim_KinInv(1)
       s5i = C2Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z5i)/(one-z5i))&
            * C2Lim_z12%wgt * subint_z124_c2lim


       FintNNLO_onloewk_is_ns(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- S

    call cut_histo(SLim_z12)
    if (SLim_z12%makecut.or.SLim_z12%flag) then

       kin(4) = zero
       FintNNLO_onloewk_is_ns(4) = zero

    else

       call res_tree_qqb(SLim_z12%AmpMom,res_lo)
       call partition_nlo_qed(SLim_z12,damp,12)


       Ec = SLim_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = SLim_z12%Lim_etaij(1,5)
       eta52 = SLim_z12%Lim_etaij(2,5)
       eta53 = SLim_z12%Lim_etaij(5,3)
       eta54 = SLim_z12%Lim_etaij(5,4)
       subint_z124_slim = qqb_ga_sub_qqba_zi(Ec,SLim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51)

       call get_qed_eik(charges_ns,SLim_z12%Lim_etaij,[1,2,3,4],5,eik)
       res_lo_tmp(:,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_lo_tmp(:,2) = res_lo(:,2) * eik(3:4) !-- up
       call get_respdf(ns_lumi,1,1,SLim_z12,res_lo_tmp,respdf)


       e5 = SLim_z12%Lim_Ei(5)

       respdf = -respdf*SLim_z12%wgt*damp/e5**2 * subint_z124_slim

       kin(4) = respdf(1)
       FintNNLO_onloewk_is_ns(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif



    !-- SC1

    call cut_histo(SC1Lim_z12)
    if (SC1Lim_z12%makecut.or.SC1Lim_z12%flag) then

       kin(5) = zero
       FintNNLO_onloewk_is_ns(5) = zero

    else

       call res_tree_qqb(SC1Lim_z12%AmpMom,res_lo)

       Ec = SC1Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = zero
       eta52 = one
       eta53 = SC1Lim_z12%Lim_etaij(1,3)              ! C51 eta53 = eta13
       eta54 = SC1Lim_z12%Lim_etaij(1,4)              ! C51 eta54 = eta14         
       eta51lim = SC1Lim_z12%Lim_etaij(1,5)

       subint_z124_sc1lim = qqb_ga_sub_qqba_zi(Ec,SC1Lim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)

       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)



       call get_respdf(ns_lumi,1,1,SC1Lim_z12,res_lo,respdf)

       e5 = SC1Lim_z12%Lim_KinInv(1)
       eta51 = SC1Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = respdf/e5**2/eta51* SC1Lim_z12%wgt * subint_z124_sc1lim

       FintNNLO_onloewk_is_ns(5) = respdf(1)
       kin(5) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- SC2

    call cut_histo(SC2Lim_z12)
    if (SC2Lim_z12%makecut.or.SC2Lim_z12%flag) then

       kin(6) = zero
       FintNNLO_onloewk_is_ns(6) = zero

    else

       call res_tree_qqb(SC2Lim_z12%AmpMom,res_lo)

       Ec = SC2Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = one
       eta51lim = eta51
       eta52 = zero
       eta53 = SC2Lim_z12%Lim_etaij(2,3)              ! C52 eta53 = eta23
       eta54 = SC2Lim_z12%Lim_etaij(2,4)              ! C52 eta54 = eta24         
       subint_z124_sc2lim = qqb_ga_sub_qqba_zi(Ec,SC2Lim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)

       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)



       call get_respdf(ns_lumi,1,1,SC2Lim_z12,res_lo,respdf)

       e5 = SC2Lim_z12%Lim_KinInv(1)
       eta52 = SC2Lim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]

       respdf = respdf/e5**2/eta52* SC2Lim_z12%wgt * subint_z124_sc2lim

       FintNNLO_onloewk_is_ns(6) = respdf(1)
       kin(6) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif




    !! ----------------------------- !!
    !! FLM[1,z2,3,4 | 5gammma]       !!
    !! ----------------------------- !!

    call kinematics_nlo_z_is(yr=xx,z1=one,z2=z,HardProc=HardProc_1z2,C1Lim=C1Lim_1z2,C2Lim=C2Lim_1z2,SLim=SLim_1z2,SC1Lim=SC1Lim_1z2,SC2Lim=SC2Lim_1z2,compute_etas=.true.)
    HardProc_1z2%ids(1:5) = [0,0,id_el,-id_el,id_a]
    C1Lim_1z2%ids(1:4)    = [0,0,id_el,-id_el]
    SLim_1z2%ids(1:4)     = [0,0,id_el,-id_el]
    C2Lim_1z2%ids(1:4)    = [0,0,id_el,-id_el]
    SC1Lim_1z2%ids(1:4)   = [0,0,id_el,-id_el]
    sC2Lim_1z2%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard

    call cut_histo(HardProc_1z2)
    if (HardProc_1z2%makecut.or.HardProc_1z2%flag) then

       kin(7) = zero
       FintNNLO_onloewk_is_ns(7) = zero

    else

       call res_tree_a_qqb(HardProc_1z2%AmpMom,res_nlo)
       call partition_nlo_qed(HardProc_1z2,damp,12)

       Ec = HardProc_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = HardProc_1z2%Lim_etaij(1,5)
       eta52 = HardProc_1z2%Lim_etaij(2,5)
       eta53 = HardProc_1z2%Lim_etaij(5,3)
       eta54 = HardProc_1z2%Lim_etaij(5,4)
       eta52lim = eta52
       subint_1z24_hard = qqb_ga_sub_qqba_zi(Ec,HardProc_1z2%muf(1),z,eta52,eta51,eta53,eta54,eta52lim)

       call get_respdf(ns_lumi,1,1,HardProc_1z2,res_nlo,respdf)

       respdf = respdf*HardProc_1z2%wgt*subint_1z24_hard * damp

       kin(7) = respdf(1)
       FintNNLO_onloewk_is_ns(7) = kin(7)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1

    call cut_histo(C1Lim_1z2)
    if (C1Lim_1z2%makecut.or.C1Lim_1z2%flag) then

       kin(8) = zero
       FintNNLO_onloewk_is_ns(8) = zero

    else

       call res_tree_qqb(C1Lim_1z2%AmpMom,res_lo)

       Ec = C1Lim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       
       eta51 = zero
       eta52 = one
       eta52lim = eta52
       eta53 = C1Lim_1z2%Lim_etaij(1,3)              ! C51 eta53 = eta13
       eta54 = C1Lim_1z2%Lim_etaij(1,4)              ! C51 eta54 = eta14         
       subint_1z24_c1lim = qqb_ga_sub_qqba_zi(Ec,C1Lim_1z2%muf(1),z,eta52,eta51,eta53,eta54,eta52lim)


       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       call get_respdf(ns_lumi,1,1,C1Lim_1z2,res_lo,respdf)

       z5i = C1Lim_1z2%Lim_KinInv(1)
       s5i = C1Lim_1z2%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z5i)/(one-z5i))&
            * C1Lim_1z2%wgt * subint_1z24_c1lim

       FintNNLO_onloewk_is_ns(8) = respdf(1)
       kin(8) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C2

    call cut_histo(C2Lim_1z2)
    if (C2Lim_1z2%makecut.or.C2Lim_1z2%flag) then

       kin(9) = zero
       FintNNLO_onloewk_is_ns(9) = zero

    else
       call res_tree_qqb(C2Lim_1z2%AmpMom,res_lo)

       Ec = C2Lim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       eta52lim = C2Lim_1z2%Lim_etaij(2,5)
       

       eta51 = one
       eta52 = zero
       eta53 = C2Lim_1z2%Lim_etaij(2,3)              ! C52 eta53 = eta23
       eta54 = C2Lim_1z2%Lim_etaij(2,4)              ! C52 eta54 = eta24         
       subint_1z24_c2lim = qqb_ga_sub_qqba_zi(Ec,C2Lim_1z2%muf(1),z,eta52,eta51,eta53,eta54,eta52lim)

       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       call get_respdf(ns_lumi,1,1,C2Lim_1z2,res_lo,respdf)

       z5i = C2Lim_1z2%Lim_KinInv(1)
       s5i = C2Lim_1z2%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z5i)/(one-z5i))&
            * C2Lim_1z2%wgt * subint_1z24_c2lim


       FintNNLO_onloewk_is_ns(9) = respdf(1)
       kin(9) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- S

    call cut_histo(SLim_1z2)
    if (SLim_1z2%makecut.or.SLim_1z2%flag) then

       kin(10) = zero
       FintNNLO_onloewk_is_ns(10) = zero

    else

       call res_tree_qqb(SLim_1z2%AmpMom,res_lo)
       call partition_nlo_qed(SLim_1z2,damp,12)


       Ec = SLim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = SLim_1z2%Lim_etaij(1,5)
       eta52 = SLim_1z2%Lim_etaij(2,5)
       eta53 = SLim_1z2%Lim_etaij(5,3)
       eta54 = SLim_1z2%Lim_etaij(5,4)
       subint_1z24_slim = qqb_ga_sub_qqba_zi(Ec,SLim_1z2%muf(1),z,eta52,eta51,eta53,eta54,eta52)

       call get_qed_eik(charges_ns,SLim_1z2%Lim_etaij,[1,2,3,4],5,eik)
       res_lo_tmp(:,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_lo_tmp(:,2) = res_lo(:,2) * eik(3:4) !-- up
       call get_respdf(ns_lumi,1,1,SLim_1z2,res_lo_tmp,respdf)


       e5 = SLim_1z2%Lim_Ei(5)

       respdf = -respdf*SLim_1z2%wgt*damp/e5**2 * subint_1z24_slim

       kin(10) = respdf(1)
       FintNNLO_onloewk_is_ns(10) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif



    !-- SC1

    call cut_histo(SC1Lim_1z2)
    if (SC1Lim_1z2%makecut.or.SC1Lim_1z2%flag) then

       kin(11) = zero
       FintNNLO_onloewk_is_ns(11) = zero

    else

       call res_tree_qqb(SC1Lim_1z2%AmpMom,res_lo)

       Ec = SC1Lim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = zero
       eta52 = one
       eta52lim = eta52
       eta53 = SC1Lim_1z2%Lim_etaij(1,3)              ! C51 eta53 = eta13
       eta54 = SC1Lim_1z2%Lim_etaij(1,4)              ! C51 eta54 = eta14         

       subint_1z24_sc1lim = qqb_ga_sub_qqba_zi(Ec,SC1Lim_1z2%muf(1),z,eta52,eta51,eta53,eta54,eta52lim)

       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)



       call get_respdf(ns_lumi,1,1,SC1Lim_1z2,res_lo,respdf)

       e5 = SC1Lim_1z2%Lim_KinInv(1)
       eta51 = SC1Lim_1z2%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = respdf/e5**2/eta51* SC1Lim_1z2%wgt * subint_1z24_sc1lim

       FintNNLO_onloewk_is_ns(11) = respdf(1)
       kin(11) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- SC2

    call cut_histo(SC2Lim_1z2)
    if (SC2Lim_1z2%makecut.or.SC2Lim_1z2%flag) then

       kin(12) = zero
       FintNNLO_onloewk_is_ns(12) = zero

    else

       call res_tree_qqb(SC2Lim_1z2%AmpMom,res_lo)

       Ec = SC2Lim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       eta52lim = SC2Lim_1z2%Lim_etaij(2,5)

       eta51 = one
       eta51lim = eta51
       eta52 = zero
       eta53 = SC2Lim_1z2%Lim_etaij(2,3)              ! C52 eta53 = eta23
       eta54 = SC2Lim_1z2%Lim_etaij(2,4)              ! C52 eta54 = eta24         
       subint_1z24_sc2lim = qqb_ga_sub_qqba_zi(Ec,SC2Lim_1z2%muf(1),z,eta52,eta51,eta53,eta54,eta52lim)

       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)



       call get_respdf(ns_lumi,1,1,SC2Lim_1z2,res_lo,respdf)

       e5 = SC2Lim_1z2%Lim_KinInv(1)
       eta52 = SC2Lim_1z2%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]

       respdf = respdf/e5**2/eta52* SC2Lim_1z2%wgt * subint_1z24_sc2lim

       FintNNLO_onloewk_is_ns(12) = respdf(1)
       kin(12) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif




    if (is_nan(kin)) then
       print *, "kin",kin(1:18)
       stop
    endif




    !! ----------------------------- !!
    !! FLM[1,2,3,4 | 5gammma]       !!
    !! ----------------------------- !!

    call kinematics_nlo_z_is(yr=xx,z1=one,z2=one,HardProc=HardProc_12,C1Lim=C1Lim_12,C2Lim=C2Lim_12,SLim=SLim_12,SC1Lim=SC1Lim_12,SC2Lim=SC2Lim_12,compute_etas=.true.)
    HardProc_12%ids(1:5) = [0,0,id_el,-id_el,id_a]
    C1Lim_12%ids(1:4)    = [0,0,id_el,-id_el]
    SLim_12%ids(1:4)     = [0,0,id_el,-id_el]
    C2Lim_12%ids(1:4)    = [0,0,id_el,-id_el]
    SC1Lim_12%ids(1:4)   = [0,0,id_el,-id_el]
    sC2Lim_12%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard

    call cut_histo(HardProc_12)
    if (HardProc_12%makecut.or.HardProc_12%flag) then

       kin(13) = zero
       FintNNLO_onloewk_is_ns(13) = zero

    else

       call res_tree_a_qqb(HardProc_12%AmpMom,res_nlo)
       call partition_nlo_qed(HardProc_12,damp,12)

       Ec = HardProc_12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = HardProc_12%Lim_etaij(1,5)
       eta52 = HardProc_12%Lim_etaij(2,5)
       eta53 = HardProc_12%Lim_etaij(5,3)
       eta54 = HardProc_12%Lim_etaij(5,4)
       eta51lim = eta51
       eta52lim = eta52
       subint_124_hard = qqb_ga_sub_qqba(Ec,HardProc_12%muf(1)) -&
            ( qqb_ga_sub_qqba_zi_pls(Ec,HardProc_12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim) + &
              qqb_ga_sub_qqba_zi_pls(Ec,HardProc_12%muf(1),z,eta52,eta51,eta53,eta54,eta52lim))

       if (fncheck) then
          print *, "mu->",HardProc_12%muf(1),",z->",z,",eta51->",eta51,"eta52->",eta52,",eta53->",eta53,",eta54->",eta54,",Ec->",Ec
          print *,"subint 124", subint_124_hard
          pause
       endif
            

       call get_respdf(ns_lumi,1,1,HardProc_12,res_nlo,respdf)

       respdf = respdf*HardProc_12%wgt*subint_124_hard * damp

       kin(13) = respdf(1)
       FintNNLO_onloewk_is_ns(13) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1

    call cut_histo(C1Lim_12)
    if (C1Lim_12%makecut.or.C1Lim_12%flag) then

       kin(14) = zero
       FintNNLO_onloewk_is_ns(14) = zero

    else

       call res_tree_qqb(C1Lim_12%AmpMom,res_lo)

       Ec = C1Lim_12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51lim = C1Lim_12%Lim_etaij(1,5)
       
       eta51 = zero
       eta52 = one
       eta53 = C1Lim_12%Lim_etaij(1,3)              ! C51 eta53 = eta13
       eta54 = C1Lim_12%Lim_etaij(1,4)              ! C51 eta54 = eta14         

       eta51lim = C1Lim_12%Lim_etaij(1,5)
       eta52lim = eta52

       subint_124_c1lim = qqb_ga_sub_qqba(Ec,C1Lim_12%muf(1)) -&
            (qqb_ga_sub_qqba_zi_pls(Ec,C1Lim_12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim) + &
             qqb_ga_sub_qqba_zi_pls(Ec,C1Lim_12%muf(1),z,eta52,eta51,eta53,eta54,eta52lim))


       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       call get_respdf(ns_lumi,1,1,C1Lim_12,res_lo,respdf)

       z5i = C1Lim_12%Lim_KinInv(1)
       s5i = C1Lim_12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z5i)/(one-z5i))&
            * C1Lim_12%wgt * subint_124_c1lim

       FintNNLO_onloewk_is_ns(14) = respdf(1)
       kin(14) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C2

    call cut_histo(C2Lim_12)
    if (C2Lim_12%makecut.or.C2Lim_12%flag) then

       kin(15) = zero
       FintNNLO_onloewk_is_ns(15) = zero

    else
       call res_tree_qqb(C2Lim_12%AmpMom,res_lo)

       Ec = C2Lim_12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = one
       eta51lim = eta51
       eta52 = zero
       eta53 = C2Lim_12%Lim_etaij(2,3)              ! C52 eta53 = eta23
       eta54 = C2Lim_12%Lim_etaij(2,4)              ! C52 eta54 = eta24         
       eta52lim = C2Lim_12%Lim_etaij(2,5)

       subint_124_c2lim = qqb_ga_sub_qqba(Ec,C2Lim_12%muf(1)) - &
            ( qqb_ga_sub_qqba_zi_pls(Ec,C2Lim_12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim) + &
              qqb_ga_sub_qqba_zi_pls(Ec,C2Lim_12%muf(1),z,eta52,eta51,eta53,eta54,eta52lim) )

       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       call get_respdf(ns_lumi,1,1,C2Lim_12,res_lo,respdf)

       z5i = C2Lim_12%Lim_KinInv(1)
       s5i = C2Lim_12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z5i)/(one-z5i))&
            * C2Lim_12%wgt * subint_124_c2lim


       FintNNLO_onloewk_is_ns(15) = respdf(1)
       kin(15) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- S

    call cut_histo(SLim_12)
    if (SLim_12%makecut.or.SLim_12%flag) then

       kin(16) = zero
       FintNNLO_onloewk_is_ns(16) = zero

    else

       call res_tree_qqb(SLim_12%AmpMom,res_lo)
       call partition_nlo_qed(SLim_12,damp,12)


       Ec = SLim_12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = SLim_12%Lim_etaij(1,5)
       eta52 = SLim_12%Lim_etaij(2,5)
       eta53 = SLim_12%Lim_etaij(5,3)
       eta54 = SLim_12%Lim_etaij(5,4)
       subint_124_slim = qqb_ga_sub_qqba(Ec,SLim_12%muf(1)) -&
            ( qqb_ga_sub_qqba_zi_pls(Ec,SLim_12%muf(1),z,eta51,eta52,eta53,eta54,eta51) + &
              qqb_ga_sub_qqba_zi_pls(Ec,SLim_12%muf(1),z,eta52,eta51,eta53,eta54,eta52) )


       call get_qed_eik(charges_ns,SLim_12%Lim_etaij,[1,2,3,4],5,eik)
       res_lo_tmp(:,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_lo_tmp(:,2) = res_lo(:,2) * eik(3:4) !-- up
       call get_respdf(ns_lumi,1,1,SLim_12,res_lo_tmp,respdf)


       e5 = SLim_12%Lim_Ei(5)

       respdf = -respdf*SLim_12%wgt*damp/e5**2 * subint_124_slim

       kin(16) = respdf(1)
       FintNNLO_onloewk_is_ns(16) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif



    !-- SC1

    call cut_histo(SC1Lim_12)
    if (SC1Lim_12%makecut.or.SC1Lim_12%flag) then

       kin(17) = zero
       FintNNLO_onloewk_is_ns(17) = zero

    else

       call res_tree_qqb(SC1Lim_12%AmpMom,res_lo)

       Ec = SC1Lim_12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = zero
       eta52 = one
       eta53 = SC1Lim_12%Lim_etaij(1,3)              ! C51 eta53 = eta13
       eta54 = SC1Lim_12%Lim_etaij(1,4)              ! C51 eta54 = eta14         
       eta51lim = SC1Lim_12%Lim_etaij(1,5)
       eta52lim = eta52

       subint_124_sc1lim = qqb_ga_sub_qqba(Ec,SC1Lim_12%muf(1)) - &
            ( qqb_ga_sub_qqba_zi_pls(Ec,SC1Lim_12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim) + &
              qqb_ga_sub_qqba_zi_pls(Ec,SC1Lim_12%muf(1),z,eta52,eta51,eta53,eta54,eta52lim) )

       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)



       call get_respdf(ns_lumi,1,1,SC1Lim_12,res_lo,respdf)

       e5 = SC1Lim_12%Lim_KinInv(1)
       eta51 = SC1Lim_12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = respdf/e5**2/eta51* SC1Lim_12%wgt * subint_124_sc1lim

       FintNNLO_onloewk_is_ns(17) = respdf(1)
       kin(17) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- SC2

    call cut_histo(SC2Lim_12)
    if (SC2Lim_12%makecut.or.SC2Lim_12%flag) then

       kin(18) = zero
       FintNNLO_onloewk_is_ns(18) = zero

    else

       call res_tree_qqb(SC2Lim_12%AmpMom,res_lo)

       Ec = SC2Lim_12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = one
       eta51lim = eta51
       eta52 = zero
       eta53 = SC2Lim_12%Lim_etaij(2,3)              ! C52 eta53 = eta23
       eta54 = SC2Lim_12%Lim_etaij(2,4)              ! C52 eta54 = eta24         
       eta52lim = SC2Lim_12%Lim_etaij(2,5)

       subint_124_sc2lim = qqb_ga_sub_qqba(Ec,SC2Lim_12%muf(1)) -&
            ( qqb_ga_sub_qqba_zi_pls(Ec,SC2Lim_12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim) + &
              qqb_ga_sub_qqba_zi_pls(Ec,SC2Lim_12%muf(1),z,eta52,eta51,eta53,eta54,eta52lim) )

       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)



       call get_respdf(ns_lumi,1,1,SC2Lim_12,res_lo,respdf)

       e5 = SC2Lim_12%Lim_KinInv(1)
       eta52 = SC2Lim_12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]

       respdf = respdf/e5**2/eta52* SC2Lim_12%wgt * subint_124_sc2lim

       FintNNLO_onloewk_is_ns(18) = respdf(1)
       kin(18) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif




    if (is_nan(kin)) then
       print *, "kin",kin(1:18)
       stop
    endif


    ff(1) = sum(kin(1:18))


!!    if (kin(1)*kin(2)*kin(3)*kin(4)*kin(5)*kin(6) .ne. zero) then
!!!!
!!       print *, "kin in onlo z12", kin(1:6)
!!       print *, "kin in onlo 1z2", kin(7:12)
!!       print *, "kin in onlo 12", kin(13:18)
!!       print *, kin(1)+kin(7)+kin(13)
!!       print *, kin(2)+kin(8)+kin(14)
!!       print *, kin(3)+kin(9)+kin(15)
!!       print *, kin(4)+kin(10)+kin(16)
!!       print *, kin(5)+kin(11)+kin(17)
!!       print *, kin(6)+kin(12)+kin(18)
!!       print *, "ff",ff(1)
!!       pause
!!    endif


    call close_histo()



  end function xsect_nnlo_onloewk_is_ns_ga_raoul



  function xsect_nnlo_onloewk_fs_53_ns_ga_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_onloewk_fs_53_ns_ga_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nnlo_onloewk_fs_53_ns_ga_raoul = xsect_nnlo_onloewk_5i_ns_ga_raoul(yRnd,ff,vegasweight,3,4)
    
  end function xsect_nnlo_onloewk_fs_53_ns_ga_raoul

  function xsect_nnlo_onloewk_fs_54_ns_ga_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_onloewk_fs_54_ns_ga_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nnlo_onloewk_fs_54_ns_ga_raoul = xsect_nnlo_onloewk_5i_ns_ga_raoul(yRnd,ff,vegasweight,4,3)
    
  end function xsect_nnlo_onloewk_fs_54_ns_ga_raoul




  function xsect_nnlo_onloewk_5i_ns_ga_raoul(yRnd,ff,vegasweight,icoll,jother)
    integer :: xsect_nnlo_onloewk_5i_ns_ga_raoul,icoll,jother
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc_z12,SLim_z12,CLim_z12,SCLim_z12
    type(KinConfig) :: HardProc_1z2,SLim_1z2,CLim_1z2,SCLim_1z2
    type(KinConfig) :: HardProc_12,SLim_12,CLim_12,SCLim_12
    !--                                                                                                                                                                                                      
    real(dp)    :: xx(kNLO_max_full),z,damp
    real(dp)    :: respdf(ipdf),kin(12)
    real(dp)    :: res_lo(2,2),res_lo_tmp(2,2),res_nlo(2,2)
    real(dp)    :: Ec, e5, z5i, s5i,eik(4),eta51,eta52,eta53,eta54,eta51lim,eta52lim
    real(dp)    :: subint_z124_hard,subint_1z24_hard,subint_124_hard
    real(dp)    :: subint_z124_slim,subint_1z24_slim,subint_124_slim
    real(dp)    :: subint_z124_clim,subint_1z24_clim,subint_124_clim
    real(dp)    :: subint_z124_sclim,subint_1z24_sclim,subint_124_sclim

    xsect_nnlo_onloewk_5i_ns_ga_raoul = 0

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

    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. ((one-z) .lt. buff_z)) then
       failed_points = failed_points + 1
       return
    endif
!    xx(xE)   = 1E-8_dp
!    xx(xrho) = 1E-8_dp
!    z = one - 1E-10_dp
    
    call open_histo()


    !! ----------------------------- !!
    !! FLM[z1,2,3,4 | 5gammma]       !!
    !! ----------------------------- !!

!    call kinematics_nlo_z_fs(yr=xx,z1=z,z2=one,icoll,jother,HardProc=HardProc_z12,CLim=CLim_z12,SLim=SLim_z12,SCLim=SCLim_z12)
    call kinematics_nlo_z_fs(xx,z,one,icoll,jother,HardProc=HardProc_z12,CLim=CLim_z12,SCLim=SCLim_z12,SLim=SLim_z12)
    HardProc_z12%ids(1:5) = [0,0,id_el,-id_el,id_a]
    CLim_z12%ids(1:4)    = [0,0,id_el,-id_el]
    SLim_z12%ids(1:4)     = [0,0,id_el,-id_el]
    SCLim_z12%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard

    call cut_histo(HardProc_z12)
    if (HardProc_z12%makecut.or.HardProc_z12%flag) then

       kin(1) = zero
       FintNNLO_onloewk_5i_ns(1) = zero

    else

       call res_tree_a_qqb(HardProc_z12%AmpMom,res_nlo)
       call partition_nlo_qed(HardProc_z12,damp,icoll)

       Ec = HardProc_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = HardProc_z12%Lim_etaij(1,5)
       eta52 = HardProc_z12%Lim_etaij(2,5)
       eta53 = HardProc_z12%Lim_etaij(3,5)
       eta54 = HardProc_z12%Lim_etaij(4,5)
       eta51lim = eta51
       subint_z124_hard = qqb_ga_sub_qqba_zi(Ec,HardProc_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51lim)

       call get_respdf(ns_lumi,1,1,HardProc_z12,res_nlo,respdf)
       respdf = respdf*HardProc_z12%wgt*subint_z124_hard * damp

       kin(1) = respdf(1)
       FintNNLO_onloewk_5i_ns(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C

    call cut_histo(CLim_z12)
    if (CLim_z12%makecut.or.CLim_z12%flag) then

       kin(2) = zero
       FintNNLO_onloewk_5i_ns(2) = zero

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

       subint_z124_clim = qqb_ga_sub_qqba_zi(Ec,CLim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51)


       call get_respdf(ns_lumi,1,1,CLim_z12,res_lo,respdf)

       z5i = CLim_z12%Lim_KinInv(1)
       s5i = CLim_z12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]

       respdf = (-one)*respdf*(two/s5i)*(Q_lep2*Pqg(z5i))&
            * CLim_z12%wgt * subint_z124_clim

       FintNNLO_onloewk_5i_ns(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif


    !-- S and SC

    call cut_histo(SLim_z12)
    if (SLim_z12%makecut.or.SLim_z12%flag) then

       kin(3) = zero
       FintNNLO_onloewk_5i_ns(3) = zero

    else

       call res_tree_qqb(SLim_z12%AmpMom,res_lo)
       call partition_nlo_qed(SLim_z12,damp,icoll)


       ! soft
       Ec = SLim_z12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = SLim_z12%Lim_etaij(1,5)
       eta52 = SLim_z12%Lim_etaij(2,5)
       eta53 = SLim_z12%Lim_etaij(3,5)
       eta54 = SLim_z12%Lim_etaij(4,5)
       subint_z124_slim = qqb_ga_sub_qqba_zi(Ec,SLim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51)

       call get_qed_eik(charges_ns,SLim_z12%Lim_etaij,[1,2,3,4],5,eik)
       res_lo_tmp(:,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_lo_tmp(:,2) = res_lo(:,2) * eik(3:4) !-- up
       call get_respdf(ns_lumi,1,1,SLim_z12,res_lo_tmp,respdf)

       e5 = SLim_z12%Lim_Ei(5)

       respdf = -respdf*SLim_z12%wgt*damp/e5**2 * subint_z124_slim       

       kin(3) = respdf(1)
       FintNNLO_onloewk_5i_ns(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif



    !-- SC

    call cut_histo(SCLim_z12)
    if (SCLim_z12%makecut.or.SCLim_z12%flag) then

       kin(4) = zero
       FintNNLO_onloewk_5i_ns(4) = zero

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
          eta52 = sCLim_z12%Lim_etaij(2,4)
          eta53 = SCLim_z12%Lim_etaij(3,4)
          eta54 = zero
       else
          print *, "incorrect icoll and jother", icoll, jother
          stop
       endif

       subint_z124_sclim = qqb_ga_sub_qqba_zi(Ec,SCLim_z12%muf(1),z,eta51,eta52,eta53,eta54,eta51)

       call get_respdf(ns_lumi,1,1,SCLim_z12,res_lo,respdf)

       e5 = SCLim_z12%Lim_KinInv(1)
       eta51 = SCLim_z12%Lim_KinInv(2)


       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = respdf/e5**2/eta51* SCLim_z12%wgt * subint_z124_sclim

       FintNNLO_onloewk_5i_ns(4) = respdf(1)
       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !! ----------------------------- !!
    !! FLM[1,z2,3,4 | 5gammma]       !!
    !! ----------------------------- !!

!    call kinematics_nlo_z_fs(yr=xx,z1=z,z2=one,icoll,jother,HardProc=HardProc_z12,CLim=CLim_z12,SLim=SLim_z12,SCLim=SCLim_z12)
    call kinematics_nlo_z_fs(xx,one,z,icoll,jother,HardProc_1z2,CLim_1z2,SCLim_1z2,SLim_1z2)
    HardProc_1z2%ids(1:5) = [0,0,id_el,-id_el,id_a]
    CLim_1z2%ids(1:4)    = [0,0,id_el,-id_el]
    SLim_1z2%ids(1:4)     = [0,0,id_el,-id_el]
    SCLim_1z2%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard

    call cut_histo(HardProc_1z2)
    if (HardProc_1z2%makecut.or.HardProc_1z2%flag) then

       kin(5) = zero
       FintNNLO_onloewk_5i_ns(5) = zero

    else

       call res_tree_a_qqb(HardProc_1z2%AmpMom,res_nlo)
       call partition_nlo_qed(HardProc_1z2,damp,icoll)

       Ec = HardProc_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = HardProc_1z2%Lim_etaij(1,5)
       eta52 = HardProc_1z2%Lim_etaij(2,5)
       eta53 = HardProc_1z2%Lim_etaij(3,5)
       eta54 = HardProc_1z2%Lim_etaij(4,5)
       eta52lim = eta52
       subint_1z24_hard = qqb_ga_sub_qqba_zi(Ec,HardProc_1z2%muf(1),z,eta52,eta51,eta53,eta54,eta52lim)

       call get_respdf(ns_lumi,1,1,HardProc_1z2,res_nlo,respdf)

       respdf = respdf*HardProc_1z2%wgt*subint_1z24_hard * damp

       kin(5) = respdf(1)
       FintNNLO_onloewk_5i_ns(5) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C

    call cut_histo(CLim_1z2)
    if (CLim_1z2%makecut.or.CLim_1z2%flag) then

       kin(6) = zero
       FintNNLO_onloewk_5i_ns(6) = zero

    else

       call res_tree_qqb(CLim_1z2%AmpMom,res_lo)

       Ec = CLim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2

       if (icoll .eq. 3 .and. jother .eq. 4) then
          eta51 = CLim_1z2%Lim_etaij(1,3)
          eta52 = CLim_1z2%Lim_etaij(2,3)
          eta53 = zero
          eta54 = CLim_1z2%Lim_etaij(3,4)
       elseif (icoll .eq. 4 .and. jother .eq. 3) then
          eta51 = CLim_1z2%Lim_etaij(1,4)
          eta52 = CLim_1z2%Lim_etaij(2,4)
          eta53 = CLim_1z2%Lim_etaij(3,4)
          eta54 = zero
       else
          print *, "incorrect icoll and jother", icoll, jother
          stop
       endif
       
       subint_1z24_clim = qqb_ga_sub_qqba_zi(Ec,CLim_1z2%muf(1),z,eta52,eta51,eta53,eta54,eta52)


       call get_respdf(ns_lumi,1,1,CLim_1z2,res_lo,respdf)

       z5i = CLim_1z2%Lim_KinInv(1)
       s5i = CLim_1z2%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]

       respdf = (-one)*respdf*(two/s5i)*(Q_lep2*Pqg(z5i))&
            * CLim_1z2%wgt * subint_1z24_clim

       FintNNLO_onloewk_5i_ns(6) = respdf(1)
       kin(6) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif


    !-- S and SC

    call cut_histo(SLim_1z2)
    if (SLim_1z2%makecut.or.SLim_1z2%flag) then

       kin(7) = zero
       FintNNLO_onloewk_5i_ns(7) = zero

    else

       call res_tree_qqb(SLim_1z2%AmpMom,res_lo)
       call partition_nlo_qed(SLim_1z2,damp,icoll)


       ! soft
       Ec = SLim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = SLim_1z2%Lim_etaij(1,5)
       eta52 = SLim_1z2%Lim_etaij(2,5)
       eta53 = SLim_1z2%Lim_etaij(3,5)
       eta54 = SLim_1z2%Lim_etaij(4,5)
       subint_1z24_slim = qqb_ga_sub_qqba_zi(Ec,SLim_1z2%muf(1),z,eta52,eta51,eta53,eta54,eta52)

       call get_qed_eik(charges_ns,SLim_1z2%Lim_etaij,[1,2,3,4],5,eik)
       res_lo_tmp(:,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_lo_tmp(:,2) = res_lo(:,2) * eik(3:4) !-- up
       call get_respdf(ns_lumi,1,1,SLim_1z2,res_lo_tmp,respdf)

       e5 = SLim_1z2%Lim_Ei(5)

       respdf = -respdf*SLim_1z2%wgt*damp/e5**2 * subint_1z24_slim       

       kin(7) = respdf(1)
       FintNNLO_onloewk_5i_ns(7) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif



    !-- SC

    call cut_histo(SCLim_1z2)
    if (SCLim_1z2%makecut.or.SCLim_1z2%flag) then

       kin(8) = zero
       FintNNLO_onloewk_5i_ns(8) = zero

    else

       call res_tree_qqb(SCLim_1z2%AmpMom,res_lo)

       Ec = SCLim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2

       if (icoll .eq. 3 .and. jother .eq. 4) then
          eta51 = SCLim_1z2%Lim_etaij(1,3)
          eta52 = SCLim_1z2%Lim_etaij(2,3)
          eta53 = zero
          eta54 = SCLim_1z2%Lim_etaij(3,4)
       elseif (icoll .eq. 4 .and. jother .eq. 3) then
          eta51 = SCLim_1z2%Lim_etaij(1,4)
          eta52 = SCLim_1z2%Lim_etaij(2,4)
          eta53 = SCLim_1z2%Lim_etaij(3,4)
          eta54 = zero
       else
          print *, "incorrect icoll and jother", icoll, jother
          stop
       endif

       subint_1z24_sclim = qqb_ga_sub_qqba_zi(Ec,SCLim_1z2%muf(1),z,eta52,eta51,eta53,eta54,eta52)

       call get_respdf(ns_lumi,1,1,SCLim_1z2,res_lo,respdf)

       e5 = SCLim_1z2%Lim_KinInv(1)
       eta51 = SCLim_1z2%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = respdf/e5**2/eta51* SCLim_1z2%wgt * subint_1z24_sclim

       FintNNLO_onloewk_5i_ns(8) = respdf(1)
       kin(8) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif



    !! ----------------------------- !!
    !! FLM[1,2,3,4 | 5gammma]       !!
    !! ----------------------------- !!

    call kinematics_nlo_z_fs(xx,one,one,icoll,jother,HardProc_12,CLim_12,SCLim_12,SLim_12)


    HardProc_12%ids(1:5) = [0,0,id_el,-id_el,id_a]
    CLim_12%ids(1:4)    = [0,0,id_el,-id_el]
    SLim_12%ids(1:4)     = [0,0,id_el,-id_el]
    SCLim_12%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard

    call cut_histo(HardProc_12)
    if (HardProc_12%makecut.or.HardProc_12%flag) then

       kin(9) = zero
       FintNNLO_onloewk_5i_ns(9) = zero

    else

       call res_tree_a_qqb(HardProc_12%AmpMom,res_nlo)
       call partition_nlo_qed(HardProc_12,damp,icoll)

       Ec = HardProc_12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = HardProc_12%Lim_etaij(1,5)
       eta52 = HardProc_12%Lim_etaij(2,5)
       eta53 = HardProc_12%Lim_etaij(3,5)
       eta54 = HardProc_12%Lim_etaij(4,5)
       eta52lim = eta52
       subint_124_hard = qqb_ga_sub_qqba(Ec,HardProc_12%muf(1)) - &
            (qqb_ga_sub_qqba_zi_pls(Ec,HardProc_12%muf(1),z,eta52,eta51,eta53,eta54,eta52) + &
             qqb_ga_sub_qqba_zi_pls(Ec,HardProc_12%muf(1),z,eta51,eta52,eta53,eta54,eta51) )

       call get_respdf(ns_lumi,1,1,HardProc_12,res_nlo,respdf)

       respdf = respdf*HardProc_12%wgt*subint_124_hard * damp

       kin(9) = respdf(1)
       FintNNLO_onloewk_5i_ns(9) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C

    call cut_histo(CLim_12)
    if (CLim_12%makecut.or.CLim_12%flag) then

       kin(10) = zero
       FintNNLO_onloewk_5i_ns(10) = zero

    else

       call res_tree_qqb(CLim_12%AmpMom,res_lo)

       Ec = CLim_12%Lim_Ei(1)                      ! Ec = E1 = E2

       if (icoll .eq. 3 .and. jother .eq. 4) then
          eta51 = CLim_12%Lim_etaij(1,3)
          eta52 = CLim_12%Lim_etaij(2,3)
          eta53 = zero
          eta54 = CLim_12%Lim_etaij(3,4)
       elseif (icoll .eq. 4 .and. jother .eq. 3) then
          eta51 = CLim_12%Lim_etaij(1,4)
          eta52 = CLim_12%Lim_etaij(2,4)
          eta53 = CLim_12%Lim_etaij(3,4)
          eta54 = zero
       else
          print *, "incorrect icoll and jother", icoll, jother
          stop
       endif
       
       subint_124_clim = qqb_ga_sub_qqba(Ec,CLim_12%muf(1)) - &
            ( qqb_ga_sub_qqba_zi_pls(Ec,CLim_12%muf(1),z,eta52,eta51,eta53,eta54,eta52) + &
              qqb_ga_sub_qqba_zi_pls(Ec,CLim_12%muf(1),z,eta51,eta52,eta53,eta54,eta51) )


       call get_respdf(ns_lumi,1,1,CLim_12,res_lo,respdf)

       z5i = CLim_12%Lim_KinInv(1)
       s5i = CLim_12%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]

       respdf = (-one)*respdf*(two/s5i)*(Q_lep2*Pqg(z5i))&
            * CLim_12%wgt * subint_124_clim

       FintNNLO_onloewk_5i_ns(10) = respdf(1)
       kin(10) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif


    !-- S and SC

    call cut_histo(SLim_12)
    if (SLim_12%makecut.or.SLim_12%flag) then

       kin(11) = zero
       FintNNLO_onloewk_5i_ns(11) = zero

    else

       call res_tree_qqb(SLim_12%AmpMom,res_lo)
       call partition_nlo_qed(SLim_12,damp,icoll)


       ! soft
       Ec = SLim_12%Lim_Ei(1)                      ! Ec = E1 = E2
       eta51 = SLim_12%Lim_etaij(1,5)
       eta52 = SLim_12%Lim_etaij(2,5)
       eta53 = SLim_12%Lim_etaij(3,5)
       eta54 = SLim_12%Lim_etaij(4,5)

       subint_124_slim = qqb_ga_sub_qqba(Ec,SLim_12%muf(1)) - &
            ( qqb_ga_sub_qqba_zi_pls(Ec,SLim_12%muf(1),z,eta52,eta51,eta53,eta54,eta52) + &
              qqb_ga_sub_qqba_zi_pls(Ec,SLim_12%muf(1),z,eta51,eta52,eta53,eta54,eta51) )

 
       call get_qed_eik(charges_ns,SLim_12%Lim_etaij,[1,2,3,4],5,eik)
       res_lo_tmp(:,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_lo_tmp(:,2) = res_lo(:,2) * eik(3:4) !-- up
       call get_respdf(ns_lumi,1,1,SLim_12,res_lo_tmp,respdf)

       e5 = SLim_12%Lim_Ei(5)

       respdf = -respdf*SLim_12%wgt*damp/e5**2 * subint_124_slim       

       kin(11) = respdf(1)
       FintNNLO_onloewk_5i_ns(11) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif



    !-- SC

    call cut_histo(SCLim_12)
    if (SCLim_12%makecut.or.SCLim_12%flag) then

       kin(12) = zero
       FintNNLO_onloewk_5i_ns(12) = zero

    else

       call res_tree_qqb(SCLim_12%AmpMom,res_lo)

       Ec = SCLim_12%Lim_Ei(1)                      ! Ec = E1 = E2

       eta51 = SCLim_12%Lim_etaij(1,5)
       eta52 = SCLim_12%Lim_etaij(2,5)

       if (icoll .eq. 3 .and. jother .eq. 4) then
          eta51 = SCLim_12%Lim_etaij(1,3)
          eta52 = SCLim_12%Lim_etaij(2,3)
          eta53 = zero
          eta54 = SCLim_12%Lim_etaij(3,4)
       elseif (icoll .eq. 4 .and. jother .eq. 3) then
          eta51 = SCLim_12%Lim_etaij(1,4)
          eta52 = SCLim_12%Lim_etaij(2,4)
          eta53 = SCLim_12%Lim_etaij(3,4)
          eta54 = zero
       else
          print *, "incorrect icoll and jother", icoll, jother
          stop
       endif

       subint_124_sclim = qqb_ga_sub_qqba(Ec,SCLim_12%muf(1)) - &
            ( qqb_ga_sub_qqba_zi_pls(Ec,SCLim_12%muf(1),z,eta52,eta51,eta53,eta54,eta52) + &
              qqb_ga_sub_qqba_zi_pls(Ec,SCLim_12%muf(1),z,eta51,eta52,eta53,eta54,eta51) )

       call get_respdf(ns_lumi,1,1,SCLim_12,res_lo,respdf)

       e5 = SCLim_12%Lim_KinInv(1)
       eta51 = SCLim_12%Lim_KinInv(2)


       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = respdf/e5**2/eta51* SCLim_12%wgt * subint_124_sclim



       FintNNLO_onloewk_5i_ns(12) = respdf(1)
       kin(12) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    if (is_nan(kin)) then
       print *, "kin",kin(1:12)
       stop
    endif


    ff(1) = sum(kin(1:12))


!!    if (kin(1)*kin(2)*kin(3)*kin(4) .ne. zero) then
!!!!
!!       print *, "kin in onlo z12", kin(1:4)
!!       print *, "kin in onlo 1z2", kin(5:8)
!!       print *, "kin in onlo 12", kin(9:12)
!!       print *, kin(1)+kin(5)+kin(9)
!!       print *, kin(2)+kin(6)+kin(10)
!!       print *, kin(3)+kin(7)+kin(11)
!!       print *, kin(4)+kin(8)+kin(12)
!!       print *, "ff",ff(1)
!!       pause
!!    endif


    call close_histo()




  end function xsect_nnlo_onloewk_5i_ns_ga_raoul
  



  function xsect_nnlo_sub12_ns_ga_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_sub12_ns_ga_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc_12,LOProc_z12,LOProc_1zb2,LOProc_zzb
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintLO_s(4),kin(4)
    real(dp) :: res_lo(2,2), res_tmp(2,2),respdf(ipdf)
    real(dp) :: Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,s13,s23,s14,s24,s34
    real(dp) :: intsub_u, intsub_d, intsub_zzb_u(1:4), intsub_zzb_d(1:4)
    real(dp) :: intsub_z_u(0:1),intsub_z_d(0:1),intsub_zb_u(0:1),intsub_zb_d(0:1)
    real(dp) :: intsub_ub, intsub_db, intsub_zzb_ub(1:4), intsub_zzb_db(1:4)
    real(dp) :: intsub_z_ub(0:1),intsub_z_db(0:1),intsub_zb_ub(0:1),intsub_zb_db(0:1)
    real(dp) :: z,zbar,mu

    xsect_nnlo_sub12_ns_ga_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))               
    z   =buff+onet*real(yRnd(kLO_max_full),dp)
    zbar=buff+onet*real(yRnd(kLO_max_full+1),dp)

    if ((one-z) .lt. buff_z .or. (one-zbar) .lt. buff_z) then
       failed_points = failed_points + 1
       return
    endif

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

    call kinematics_lo_zzb(xx,one,one, LOProc_12)   
    call kinematics_lo_zzb(xx,z,  one, LOProc_z12)   
    call kinematics_lo_zzb(xx,one,zbar,LOProc_1zb2) 
    call kinematics_lo_zzb(xx,z,  zbar,LOProc_zzb)  

    LOProc_12%ids(1:4) = [0,0,id_el,-id_el]
    LOProc_z12%ids(1:4) = [0,0,id_el,-id_el]
    LOProc_1zb2%ids(1:4) = [0,0,id_el,-id_el]
    LOProc_zzb%ids(1:4) = [0,0,id_el,-id_el]

    !! FLM[1,2]

    call cut_histo(LOProc_12)
    if (LOProc_12%makecut.or.LOProc_12%flag) then

       kin(1) = zero
       FintLO_s(1) = zero

    else

       call res_tree_qqb(LOProc_12%AmpMom,res_lo)
       
       mu = LOProc_12%muf(1)

       E1 = LOProc_12%Lim_Ei(1)
       E2 = LOProc_12%Lim_Ei(2)
       E3 = LOProc_12%Lim_Ei(3)
       E4 = LOProc_12%Lim_Ei(4)
       
       s13 = LOProc_12%Lim_KinInv(1)
       s23 = LOProc_12%Lim_KinInv(2)
       s14 = LOProc_12%Lim_KinInv(3)
       s24 = LOProc_12%Lim_KinInv(4)
       s34 = LOProc_12%Lim_KinInv(5)
       
       eta31 = LOProc_12%Lim_etaij(3,1)
       eta32 = LOProc_12%Lim_etaij(3,2)
       eta41 = LOProc_12%Lim_etaij(4,1)
       eta42 = LOProc_12%Lim_etaij(4,2)
       eta34 = LOProc_12%Lim_etaij(3,4)

       Ec = E1 


       ! get the relevant int. subtr. functions
       intsub_zzb_u = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,Qup)
       intsub_zzb_d = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,Qdn)
       intsub_z_u   = qqb_ga_sub_qqbz12(mu,z,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qup,Q_lep)
       intsub_z_d   = qqb_ga_sub_qqbz12(mu,z,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qdn,Q_lep)
       intsub_zb_u  = qqb_ga_sub_qqbz12(mu,zbar,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qup,Q_lep)
       intsub_zb_d  = qqb_ga_sub_qqbz12(mu,zbar,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qdn,Q_lep)
       intsub_u     = qqb_ga_sub_qqb12(mu,Ec,E3,E4,eta31,eta32,eta41,eta42,eta34,s13,s14,s23,s24,Qup,Q_lep)
       intsub_d     = qqb_ga_sub_qqb12(mu,Ec,E3,E4,eta31,eta32,eta41,eta42,eta34,s13,s14,s23,s24,Qdn,Q_lep)

       intsub_zzb_ub = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,-Qup)
       intsub_zzb_db = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,-Qdn)
       intsub_z_ub   = qqb_ga_sub_qqbz12(mu,z,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qup,Q_lep)
       intsub_z_db   = qqb_ga_sub_qqbz12(mu,z,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qdn,Q_lep)
       intsub_zb_ub  = qqb_ga_sub_qqbz12(mu,zbar,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qup,Q_lep)
       intsub_zb_db  = qqb_ga_sub_qqbz12(mu,zbar,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qdn,Q_lep)
       intsub_ub     = qqb_ga_sub_qqb12(mu,Ec,E3,E4,eta31,eta32,eta41,eta42,eta34,s13,s14,s23,s24,-Qup,Q_lep)
       intsub_db     = qqb_ga_sub_qqb12(mu,Ec,E3,E4,eta31,eta32,eta41,eta42,eta34,s13,s14,s23,s24,-Qdn,Q_lep)


       if (fncheck) then
          print *, "eta31->",eta31,", eta41->",eta41,", eta32->",eta32,", eta42->",eta42,", eta34->",eta34,", mu->",mu, ", E1->",E1,", E2->",E2,", E3->",E3,", E4->",E4, ", Ec->",Ec,",z->",z,",zb->",zbar
          print *, "12 u",(intsub_u + intsub_zzb_u(4) - intsub_z_u(1) - intsub_zb_u(1))
          print *, "12 d",(intsub_d + intsub_zzb_d(4) - intsub_z_d(1) - intsub_zb_d(1))
          print *, "12 ub",(intsub_ub + intsub_zzb_ub(4) - intsub_z_ub(1) - intsub_zb_ub(1))
          print *, "12 db",(intsub_db + intsub_zzb_db(4) - intsub_z_db(1) - intsub_zb_db(1))
          pause
       endif

       res_tmp(1,2) = res_lo(1,2)*(intsub_u + intsub_zzb_u(4) - intsub_z_u(1) - intsub_zb_u(1))
       res_tmp(1,1) = res_lo(1,1)*(intsub_d + intsub_zzb_d(4) - intsub_z_d(1) - intsub_zb_d(1))

       res_tmp(2,2) = res_lo(2,2)*(intsub_ub + intsub_zzb_ub(4) - intsub_z_ub(1) - intsub_zb_ub(1))
       res_tmp(2,1) = res_lo(2,1)*(intsub_db + intsub_zzb_db(4) - intsub_z_db(1) - intsub_zb_db(1))

       call get_respdf(ns_lumi,1,1,LOProc_12,res_tmp,respdf)

       ! -- int. subtraction term
       respdf = respdf*LOProc_12%wgt 
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif


    !! FLM[z1,2]
    call cut_histo(LOProc_z12)
    if (LOProc_z12%makecut.or.LOProc_z12%flag) then

       kin(2) = zero
       FintLO_s(2) = zero

    else

       call res_tree_qqb(LOProc_z12%AmpMom,res_lo)
       
       mu = LOProc_z12%muf(1)

       E1 = LOProc_z12%Lim_Ei(1)
       E2 = LOProc_z12%Lim_Ei(2)
       E3 = LOProc_z12%Lim_Ei(3)
       E4 = LOProc_z12%Lim_Ei(4)
       
       s13 = LOProc_z12%Lim_KinInv(1)
       s23 = LOProc_z12%Lim_KinInv(2)
       s14 = LOProc_z12%Lim_KinInv(3)
       s24 = LOProc_z12%Lim_KinInv(4)
       s34 = LOProc_z12%Lim_KinInv(5)
       
       eta31 = LOProc_z12%Lim_etaij(3,1)
       eta32 = LOProc_z12%Lim_etaij(3,2)
       eta41 = LOProc_z12%Lim_etaij(4,1)
       eta42 = LOProc_z12%Lim_etaij(4,2)
       eta34 = LOProc_z12%Lim_etaij(3,4)

       Ec = E1 

       ! get the relevant int. subtr. functions
       intsub_zzb_u = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,Qup)
       intsub_zzb_d = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,Qdn)

       intsub_z_u   = qqb_ga_sub_qqbz12(mu,z,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qup,Q_lep)
       intsub_z_d   = qqb_ga_sub_qqbz12(mu,z,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qdn,Q_lep)

       intsub_zzb_ub = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,-Qup)
       intsub_zzb_db = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,-Qdn)
       intsub_z_ub   = qqb_ga_sub_qqbz12(mu,z,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qup,Q_lep)
       intsub_z_db   = qqb_ga_sub_qqbz12(mu,z,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qdn,Q_lep)

       if (fncheck) then
          print *, "eta31->",eta31,", eta41->",eta41,", eta32->",eta32,", eta42->",eta42,", eta34->",eta34,", mu->",mu, ", E1->",E1,", E2->",E2,", E3->",E3,", E4->",E4, ", Ec->",Ec,",z->",z,",zb->",zbar
          print *, "z12 u",(intsub_zzb_u(2) + intsub_z_u(0))
          print *, "z12 d",(intsub_zzb_d(2) + intsub_z_d(0))
          print *, "z12 ub", ( intsub_zzb_ub(2) + intsub_z_ub(0))
          print *, "z12 db", ( intsub_zzb_db(2) + intsub_z_db(0))
          pause
       endif

       res_tmp(1,2) = res_lo(1,2)*(intsub_zzb_u(2) + intsub_z_u(0))
       res_tmp(1,1) = res_lo(1,1)*(intsub_zzb_d(2) + intsub_z_d(0))

       res_tmp(2,2) = res_lo(2,2)*( intsub_zzb_ub(2) + intsub_z_ub(0))
       res_tmp(2,1) = res_lo(2,1)*( intsub_zzb_db(2) + intsub_z_db(0))

!       print *, "12 coeff",intsub_zzb_u(4), - intsub_z_u(1), - intsub_zb_u(1)
       
       call get_respdf(ns_lumi,1,1,LOProc_z12,res_tmp,respdf)

       ! -- int. subtraction term
       respdf = respdf*LOProc_z12%wgt 
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif


    !! FLM[1,zb2]
    call cut_histo(LOProc_1zb2)
    if (LOProc_1zb2%makecut.or.LOProc_1zb2%flag) then

       kin(3) = zero
       FintLO_s(3) = zero

    else

       call res_tree_qqb(LOProc_1zb2%AmpMom,res_lo)
       
       mu = LOProc_1zb2%muf(1)

       E1 = LOProc_1zb2%Lim_Ei(1)
       E2 = LOProc_1zb2%Lim_Ei(2)
       E3 = LOProc_1zb2%Lim_Ei(3)
       E4 = LOProc_1zb2%Lim_Ei(4)
       
       s13 = LOProc_1zb2%Lim_KinInv(1)
       s23 = LOProc_1zb2%Lim_KinInv(2)
       s14 = LOProc_1zb2%Lim_KinInv(3)
       s24 = LOProc_1zb2%Lim_KinInv(4)
       s34 = LOProc_1zb2%Lim_KinInv(5)
       
       eta31 = LOProc_1zb2%Lim_etaij(3,1)
       eta32 = LOProc_1zb2%Lim_etaij(3,2)
       eta41 = LOProc_1zb2%Lim_etaij(4,1)
       eta42 = LOProc_1zb2%Lim_etaij(4,2)
       eta34 = LOProc_1zb2%Lim_etaij(3,4)

       Ec = E1 

       ! get the relevant int. subtr. functions
       intsub_zzb_u = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,Qup)
       intsub_zzb_d = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,Qdn)
       intsub_zb_u   = qqb_ga_sub_qqbz12(mu,zbar,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qup,Q_lep)
       intsub_zb_d   = qqb_ga_sub_qqbz12(mu,zbar,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qdn,Q_lep)

       intsub_zzb_ub = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,-Qup)
       intsub_zzb_db = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,-Qdn)
       intsub_zb_ub   = qqb_ga_sub_qqbz12(mu,zbar,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qup,Q_lep)
       intsub_zb_db   = qqb_ga_sub_qqbz12(mu,zbar,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,-Qdn,Q_lep)

       if (fncheck) then
          print *, "eta31->",eta31,", eta41->",eta41,", eta32->",eta32,", eta42->",eta42,", eta34->",eta34,", mu->",mu, ", E1->",E1,", E2->",E2,", E3->",E3,", E4->",E4, ", Ec->",Ec,",zb->",zbar,",z->",z
          print *, "1zb2 u",(intsub_zzb_u(3) + intsub_zb_u(0))
          print *, "1zb2 d",(intsub_zzb_d(3) + intsub_zb_d(0))
          print *, "1zb2 ub",( intsub_zzb_ub(3) + intsub_zb_ub(0))
          print *, "1zb2 db",( intsub_zzb_db(3) + intsub_zb_db(0))
          pause
       endif

       res_tmp(1,2) = res_lo(1,2)*(intsub_zzb_u(3) + intsub_zb_u(0))
       res_tmp(1,1) = res_lo(1,1)*(intsub_zzb_d(3) + intsub_zb_d(0))

       res_tmp(2,2) = res_lo(2,2)*( intsub_zzb_ub(3) + intsub_zb_ub(0))
       res_tmp(2,1) = res_lo(2,1)*( intsub_zzb_db(3) + intsub_zb_db(0))

!       print *, "12 coeff",intsub_zzb_u(4), - intsub_z_u(1), - intsub_zb_u(1)
       
       call get_respdf(ns_lumi,1,1,LOProc_1zb2,res_tmp,respdf)

       ! -- int. subtraction term
       respdf = respdf*LOProc_1zb2%wgt 
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif


    !! FLM[z1,zb2]
    call cut_histo(LOProc_zzb)
    if (LOProc_zzb%makecut.or.LOProc_zzb%flag) then

       kin(4) = zero
       FintLO_s(4) = zero

    else

       call res_tree_qqb(LOProc_zzb%AmpMom,res_lo)
       
       mu = LOProc_zzb%muf(1)

       E1 = LOProc_zzb%Lim_Ei(1)
       Ec = E1 

       ! get the relevant int. subtr. functions
       intsub_zzb_u = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,Qup)
       intsub_zzb_d = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,Qdn)
       intsub_zzb_ub = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,-Qup)
       intsub_zzb_db = qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,-Qdn)

       if (fncheck) then
          print *, "mu->",mu,",EC->",Ec,",z->",z,",zb->",zbar
          
          print *, "zzb u",intsub_zzb_u(1)
          print *, "zzb d",intsub_zzb_d(1)
          print *, "zzb ub",intsub_zzb_ub(1)
          print *, "zzb db",intsub_zzb_db(1)
          pause
       endif

       res_tmp(1,2) = res_lo(1,2)*(intsub_zzb_u(1) )
       res_tmp(1,1) = res_lo(1,1)*(intsub_zzb_d(1) )

       res_tmp(2,2) = res_lo(2,2)*( intsub_zzb_ub(1) )
       res_tmp(2,1) = res_lo(2,1)*( intsub_zzb_db(1) )

       call get_respdf(ns_lumi,1,1,LOProc_zzb,res_tmp,respdf)

       ! -- int. subtraction term
       respdf = respdf*LOProc_zzb%wgt 
       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin(1:4))

    

!!    if (ff(1) .ne. zero) then
!!       print *, "xx",xx,z
!!       print *, "kin", kin(1:4)
!!       print *, "ff",ff(1)
!!       pause
!!    endif


!!    print *, "z,zbar",z,zbar
!!    print *, "ff",ff(1)
!!    print *, "kin",kin(1:4,1)


    call close_histo()

!    call check_ff(ff,kLO_max_full+1,3,(/xx,z/),FintNLO_s)
    
#if(_withchecks == 1)
!    FintNLO_qcd_vs = FintNLO_s
#endif
    
  end function xsect_nnlo_sub12_ns_ga_raoul






  function xsect_nnlo_subvewk_ns_ga_raoul(yRnd,ff,vegasweight)
    use mod_subtrfn_nlo_z
    integer :: xsect_nnlo_subvewk_ns_ga_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--                                                                                                                                                                                                      
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintNLO_s(3),kin(3)
    real(dp) :: res_tree(2,3), res_loop(2,3)
    real(dp) :: respdf_qqb(ipdf)
    real(dp) :: z,xis(2),q2,intsub_qqb(-1:1,nintsub)

    xsect_nnlo_subvewk_ns_ga_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    z=buff+onet*real(yRnd(kLO_max_full),dp)

    if ((one-z) .lt. buff_z) then
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
       intsub_qqb = sub_g_qqb_qqb_z(q2,LOProc%muf**2,z)

       !-- [1,2]
       call get_respdf(ns_lumi_splitb,1,1,LOProc,res_loop,respdf_qqb)

       ! -- int. subtraction term
       respdf_qqb = respdf_qqb*LOProc%wgt * (intsub_qqb(-1,1) - intsub_qqb(1,1))*two     ! delta-pls, two legs
       kin(1) = respdf_qqb(1)

       call fill_histo(respdf_qqb,vegasweight)

       !-- [z,2]
       ! change xi1 and recompute pdfs

       LOProc%PartFrac(1) = xis(1)/z
       LOProc%PartFrac(2) = xis(2)

       call get_respdf(ns_lumi_splitb,1,1,LOProc,res_loop,respdf_qqb)
       respdf_qqb = respdf_qqb*LOProc%wgt/z * intsub_qqb(0,1)

       kin(2) = respdf_qqb(1)
       call fill_histo(respdf_qqb,vegasweight)

       !-- [1,z]
       ! change xi1 and recompute pdfs
       LOProc%PartFrac(1) = xis(1)
       LOProc%PartFrac(2) = xis(2)/z


       call get_respdf(ns_lumi_splitb,1,1,LOProc,res_loop,respdf_qqb)
       respdf_qqb = respdf_qqb*LOProc%wgt/z * intsub_qqb(0,1)
       kin(3) = respdf_qqb(1)

       call fill_histo(respdf_qqb,vegasweight)

    endif


    ff(1) = sum(kin(1:3))

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

  end function xsect_nnlo_subvewk_ns_ga_raoul


 function xsect_nnlo_subvqcd_ns_ga_raoul(yRnd,ff,vegasweight)
    use mod_subtrfn_nlo_z
    integer :: xsect_nnlo_subvqcd_ns_ga_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--                                                                                                                                                                                                      
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintNLO_s(3),kin(3)
    real(dp) :: res_tree(2,2),res_loop(2,2), res_tmp(2,2)
    real(dp) :: respdf_qqb(ipdf)
    real(dp) :: z,xis(2),q2,intsub_u(-1:1,nintsub),intsub_d(-1:1,nintsub),intsub_ub(-1:1,nintsub),intsub_db(-1:1,nintsub)
    real(dp) :: eta31, eta41, eta32, eta42
    xsect_nnlo_subvqcd_ns_ga_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    z=buff+onet*real(yRnd(kLO_max_full),dp)

    if ((one-z) .lt. buff_z) then
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

    ! values of etaij for int.  subtraction functions

    eta31 =  scr(LOProc%AmpMom(1:4,1),LOProc%AmpMom(1:4,3))/two/LOProc%AmpMom(1,1)/LOProc%AmpMom(1,3)
    eta41 =  scr(LOProc%AmpMom(1:4,1),LOProc%AmpMom(1:4,4))/two/LOProc%AmpMom(1,1)/LOProc%AmpMom(1,4)
    eta32 =  scr(LOProc%AmpMom(1:4,2),LOProc%AmpMom(1:4,3))/two/LOProc%AmpMom(1,2)/LOProc%AmpMom(1,3)
    eta42 =  scr(LOProc%AmpMom(1:4,2),LOProc%AmpMom(1:4,4))/two/LOProc%AmpMom(1,2)/LOProc%AmpMom(1,4)
    LOProc%ids(1:4) = [0,0,id_el,-id_el]                                                                                                                                                                    

    ! save the unboosted xi1 and xi2
    xis(1:2) = LOProc%PartFrac(1:2)
    q2 = two*scr(LOProc%AmpMom(1:4,1),LOProc%AmpMom(1:4,2))

    call cut_histo(LOProc)
    if (LOProc%makecut.or.LOProc%flag) then

       kin = zero
       FintNLO_s = zero

    else

       call res_qcdloop_qqb(LOProc%AmpMom,res_tree,res_loop)

       ! get the relevant int. subtr. functions for different quark charges
       intsub_u = sub_a_qqb_qqb_z(q2,LOProc%muf**2,z,eta31,eta32,eta41,eta42,Qup,Q_lep)
       intsub_d = sub_a_qqb_qqb_z(q2,LOProc%muf**2,z,eta31,eta32,eta41,eta42,Qdn,Q_lep)
       intsub_ub = sub_a_qqb_qqb_z(q2,LOProc%muf**2,z,eta31,eta32,eta41,eta42,-Qup,Q_lep)
       intsub_db = sub_a_qqb_qqb_z(q2,LOProc%muf**2,z,eta31,eta32,eta41,eta42,-Qdn,Q_lep)


       !-- [1,2]
       res_tmp(1,1)  = res_loop(1,1)*(intsub_d(-1,1)  - intsub_d(+1,1))*two ! delta-pls, two legs
       res_tmp(1,2)  = res_loop(1,2)*(intsub_u(-1,1)  - intsub_u(+1,1))*two ! delta-pls, two legs
       res_tmp(2,1)  = res_loop(2,1)*(intsub_db(-1,1)  - intsub_db(+1,1))*two ! delta-pls, two legs
       res_tmp(2,2)  = res_loop(2,2)*(intsub_ub(-1,1)  - intsub_ub(+1,1))*two ! delta-pls, two legs

       call get_respdf(ns_lumi,1,1,LOProc,res_tmp,respdf_qqb)

       respdf_qqb = respdf_qqb*LOProc%wgt
       call fill_histo(respdf_qqb,vegasweight)
       kin(1) = respdf_qqb(1)

       !-- [z,2]
       ! change xi1 and recompute pdfs
       LOProc%PartFrac(1) = xis(1)/z
       LOProc%PartFrac(2) = xis(2)

       res_tmp(1,1)  = res_loop(1,1)*intsub_d(0,1)      ! reg 
       res_tmp(1,2)  = res_loop(1,2)*intsub_u(0,1)      ! reg
       res_tmp(2,1)  = res_loop(2,1)*intsub_db(0,1)      ! reg
       res_tmp(2,2)  = res_loop(2,2)*intsub_ub(0,1)      ! reg

       call get_respdf(ns_lumi,1,1,LOProc,res_tmp,respdf_qqb)

       respdf_qqb = respdf_qqb(1)*LOProc%wgt/z
       call fill_histo(respdf_qqb,vegasweight)
       kin(2) = respdf_qqb(1)

       !-- [1,z]
       ! change xi1 and recompute pdfs
       LOProc%PartFrac(1) = xis(1)
       LOProc%PartFrac(2) = xis(2)/z

       res_tmp(1,1)  = res_loop(1,1)*intsub_d(0,1)      ! reg
       res_tmp(1,2)  = res_loop(1,2)*intsub_u(0,1)      ! reg
       res_tmp(2,1)  = res_loop(2,1)*intsub_db(0,1)      ! reg
       res_tmp(2,2)  = res_loop(2,2)*intsub_ub(0,1)      ! reg

       call get_respdf(ns_lumi,1,1,LOProc,res_tmp,respdf_qqb)
       respdf_qqb = respdf_qqb(1)*LOProc%wgt/z
       kin(3) = respdf_qqb(1)


       call fill_histo(respdf_qqb,vegasweight)

    endif

    ff(1) = sum(kin(1:3))

!    if (ff(1) .ne. zero) then
!       print *, "xx",xx,z
!       print *, "kin", kin(1:3)
!       print *, "ff",ff(1)
!       stop
!    endif                                                                                                                                                                                                   



    call close_histo()


  end function xsect_nnlo_subvqcd_ns_ga_raoul



  function xsect_nnlo_sub12_ns_qqb_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_sub12_ns_qqb_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    type(KinConfig) :: LOProc_z12, LOProc_1z2
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintLO_s(4),kin(4)
    real(dp) :: res_lo(2,2), res_tmp(2,2),respdf(ipdf),respdf_tmp(ipdf)
    real(dp) :: Ec,z,mu,intsub_u, intsub_d,xis(1:2)

    xsect_nnlo_sub12_ns_qqb_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))               
    z   =buff+onet*real(yRnd(kLO_max_full),dp)

    if ((one-z) .lt. buff_z) then
       failed_points = failed_points + 1
       return
    endif


!    print *, "z",z

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

       kin(1) = zero
       FintLO_s(1) = zero

    else

       call res_tree_qqb(LOProc%AmpMom,res_lo)
!       print *, "12 wgt /z", LOProc%wgt/z
!       print *, "12 Ec",Ec
!       print *, "12 mu",mu
       
       ! get the relevant int. subtr. functions
       intsub_u = qqb_qqb_sub_qqbz12(mu,Ec,z,Qup)
       intsub_d = qqb_qqb_sub_qqbz12(mu,Ec,z,Qdn)

       if (fncheck) then
          print *, "mu->",mu,",EC->",Ec,",z->",z
          print *, "intsub_u",intsub_u
          print *, "intsub_d",intsub_d
          pause 
       endif

       res_tmp(:,1) = res_lo(:,1)*intsub_d
       res_tmp(:,2) = res_lo(:,2)*intsub_u

       !! FLM[z1,2]
       LOProc%PartFrac(1) = xis(1)/z
       LOProc%PartFrac(2) = xis(2)
       call get_respdf(ns_lumi,1,1,LOProc,res_tmp,respdf)

       respdf_tmp = respdf*LOProc%wgt / z
       kin(1) = respdf(1)

       !! FLM[1,z2]
       LOProc%PartFrac(1) = xis(1)
       LOProc%PartFrac(2) = xis(2)/z
       call get_respdf(ns_lumi,1,1,LOProc,res_tmp,respdf)

       respdf = respdf_tmp + respdf*LOProc%wgt / z
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)


    endif

    ff(1) = kin(1)

!!--check--!!! new checking code starts here
!!--check--!!    call kinematics_lo_zzb(xx,z,  one, LOProc_z12)
!!--check--!!    call kinematics_lo_zzb(xx,one,z,LOProc_1z2)
!!--check--!!
!!--check--!!    
!!--check--!!    LOProc_z12%ids(1:4) = [0,0,id_el,-id_el]
!!--check--!!    LOProc_1z2%ids(1:4) = [0,0,id_el,-id_el]
!!--check--!!
!!--check--!!    call cut_histo(LOProc_z12)
!!--check--!!    if (LOProc_z12%makecut.or.LOProc_z12%flag) then
!!--check--!!
!!--check--!!       kin(1) = zero
!!--check--!!       FintLO_s(1) = zero
!!--check--!!
!!--check--!!    else
!!--check--!!
!!--check--!!       call res_tree_qqb(LOProc_z12%AmpMom,res_lo)
!!--check--!!       
!!--check--!!       ! get the relevant int. subtr. functions
!!--check--!!       mu = LOProc_z12%muf(1)
!!--check--!!       Ec = LOProc_z12%Lim_Ei(1)
!!--check--!!!       print *, "z12 wgt", LOProc_z12%wgt
!!--check--!!!    
!!--check--!!!       print *, "z12 Ec",Ec
!!--check--!!!       print *, "z12 mu",mu
!!--check--!!!       pause
!!--check--!!
!!--check--!!       intsub_u = qqb_qqb_sub_qqbz12(mu,Ec,z,Qup)
!!--check--!!       intsub_d = qqb_qqb_sub_qqbz12(mu,Ec,z,Qdn)
!!--check--!!
!!--check--!!       res_tmp(:,1) = res_lo(:,1)*intsub_d
!!--check--!!       res_tmp(:,2) = res_lo(:,2)*intsub_u
!!--check--!!
!!--check--!!       call get_respdf(ns_lumi,1,1,LOProc_z12,res_tmp,respdf)
!!--check--!!
!!--check--!!       ! -- int. subtraction term       
!!--check--!!       respdf = respdf*LOProc_z12%wgt
!!--check--!!       kin(1) = respdf(1)
!!--check--!!
!!--check--!!       call fill_histo(respdf,vegasweight)
!!--check--!!
!!--check--!!    endif
!!--check--!!
!!--check--!!    call cut_histo(LOProc_1z2)
!!--check--!!    if (LOProc_1z2%makecut.or.LOProc_1z2%flag) then
!!--check--!!
!!--check--!!       kin(2) = zero
!!--check--!!       FintLO_s(2) = zero
!!--check--!!
!!--check--!!    else
!!--check--!!
!!--check--!!       call res_tree_qqb(LOProc_1z2%AmpMom,res_lo)
!!--check--!!       
!!--check--!!       ! get the relevant int. subtr. functions
!!--check--!!       mu = LOProc_1z2%muf(1)
!!--check--!!       Ec = LOProc_1z2%Lim_Ei(2)
!!--check--!!
!!--check--!!       intsub_u = qqb_qqb_sub_qqbz12(mu,Ec,z,Qup)
!!--check--!!       intsub_d = qqb_qqb_sub_qqbz12(mu,Ec,z,Qdn)
!!--check--!!
!!--check--!!       res_tmp(:,1) = res_lo(:,1)*intsub_d
!!--check--!!       res_tmp(:,2) = res_lo(:,2)*intsub_u
!!--check--!!
!!--check--!!       call get_respdf(ns_lumi,1,1,LOProc_1z2,res_tmp,respdf)
!!--check--!!
!!--check--!!       ! -- int. subtraction term       
!!--check--!!       respdf = respdf*LOProc_1z2%wgt
!!--check--!!       kin(2) = respdf(1)
!!--check--!!
!!--check--!!       call fill_histo(respdf,vegasweight)
!!--check--!!
!!--check--!!    endif
!!--check--!!
!!--check--!!!    kin(2) = zero
!!--check--!!    
!!--check--!!    ff(1) = sum(kin(1:2))

!!    if (ff(1) .ne. zero) then
!!       print *, "xx",xx,z
!!       print *, "kin", kin(1:4)
!!       print *, "ff",ff(1)
!!       pause
!!    endif


!!    print *, "z,zbar",z,zbar
!!    print *, "ff",ff(1)
!!    print *, "kin",kin(1:4,1)


    call close_histo()

!    call check_ff(ff,kLO_max_full+1,3,(/xx,z/),FintNLO_s)
    
#if(_withchecks == 1)
!    FintNLO_qcd_vs = FintNLO_s
#endif
    
  end function xsect_nnlo_sub12_ns_qqb_raoul


  function xsect_nnlo_sub12_ns_qq_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_sub12_ns_qq_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintLO_s(4),kin(4)
    real(dp) :: res_lo(2,2), res_tmp(2,2),res_swap(2,2),respdf(ipdf),respdf_tmp(ipdf)
    real(dp) :: Ec,z,mu,intsub_u, intsub_d,xis(1:2)

    xsect_nnlo_sub12_ns_qq_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))               
    z = buff+onet*real(yRnd(kLO_max_full),dp)

    if ((one-z) .lt. buff_z) then
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

       kin(1) = zero
       FintLO_s(1) = zero

    else

       call res_tree_qqb(LOProc%AmpMom,res_lo)
       
       ! get the relevant int. subtr. functions
       intsub_u = qq_qq_sub_qqb1z2(mu,Ec,z,Qup)
       intsub_d = qq_qq_sub_qqb1z2(mu,Ec,z,Qdn)

       if (fncheck) then
          print *, "qq->qq"
          print *, "mu->",mu,",EC->",Ec,",z->",z
          print *, "intsub_u",intsub_u
          print *, "intsub_d",intsub_d
          pause
       endif

       res_tmp(:,1) = res_lo(:,1)*intsub_d
       res_tmp(:,2) = res_lo(:,2)*intsub_u


       !! FLM[z1,2]
       LOProc%PartFrac(1) = xis(1)/z
       LOProc%PartFrac(2) = xis(2)
       res_swap(1,:)  = res_tmp(2,:)    ! TC1 acting on qq gives qbq initial state
       res_swap(2,:)  = res_tmp(1,:)
       call get_respdf(qq_lumi,1,1,LOProc,res_swap,respdf)

       respdf_tmp = respdf*LOProc%wgt / z
       kin(1) = respdf(1)

       !! FLM[1,z2]
       LOProc%PartFrac(1) = xis(1)
       LOProc%PartFrac(2) = xis(2)/z
       call get_respdf(qq_lumi,1,1,LOProc,res_tmp,respdf)   ! TC256 acting on qq gives qqb initial state

       respdf = respdf_tmp + respdf*LOProc%wgt / z
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = kin(1)
    

!!    if (ff(1) .ne. zero) then
!!       print *, "xx",xx,z
!!       print *, "kin", kin(1:4)
!!       print *, "ff",ff(1)
!!       pause
!!    endif


!!    print *, "z,zbar",z,zbar
!!    print *, "ff",ff(1)
!!    print *, "kin",kin(1:4,1)


    call close_histo()

!    call check_ff(ff,kLO_max_full+1,3,(/xx,z/),FintNLO_s)
    
#if(_withchecks == 1)
!    FintNLO_qcd_vs = FintNLO_s
#endif
    
  end function xsect_nnlo_sub12_ns_qq_raoul



  




  
end module mod_xsects_nnlo_s_qqb_raoul
  

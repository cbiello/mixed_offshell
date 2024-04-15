module mod_xsects_nnlo_s_aq
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_aux_sectors
  use mod_process
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_amplitudes_loop_ppll
  use mod_hoppet_tools
  use mod_hoppet_nnlo
  use mod_int_sub_nnlo
  use mod_splittings_bare
  use mod_eikonals
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_s_aq(1),FintNNLO_s_aq_onlo(6)
#endif

  public :: xsect_nnlo_s_aq_vqcd,xsect_nnlo_s_aq_s12
  public :: xsect_nnlo_s_aq_oqcd,xsect_nnlo_s_aq_oewk

  private

contains

  !!*************************************************************************!!

  function xsect_nnlo_s_aq_vqcd(yRnd,ff,vegasweight)
    use mod_kinematics_lo
    integer :: xsect_nnlo_s_aq_vqcd
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: respdf(ipdf)
    real(dp) :: res_tree(2,2),res_nlo(2,2),kin(1)

    xsect_nnlo_s_aq_vqcd = 0

    ff(1) = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    
#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()
    call kinematics_lo(xx,LOProc)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(LOProc)

    if (LOProc%makecut) then

       kin(1) = zero

    else

       call res_qcdloop_qqb(LOProc%AmpMom,res_tree,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       !-- F[z.p1,p2]
       call get_respdf_hoppet(xPij,PDFs,ns_lumi,1,1,LOProc,res_nlo,respdf,myPDFs1_Lmu=[xPij_Lmu])
       respdf = xn*respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_aq = kin
#endif

  end function xsect_nnlo_s_aq_vqcd

  function xsect_nnlo_s_aq_s12(yRnd,ff,vegasweight)
    use mod_kinematics_lo
    implicit none
    integer :: xsect_nnlo_s_aq_s12
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: kin(1),FintNNLO_s_aq_s12(1)
    real(dp) :: respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf),respdf_3(ipdf),respdf_tmp(ipdf,4)
    real(dp) :: res_lo_1(1,1),res_lo_2(2,2)

    xsect_nnlo_s_aq_s12 = 0

    ff(1) = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    
#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()
    call kinematics_lo(xx,LOProc)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(LOProc)

    if (LOProc%makecut) then

       kin(1) = zero
       FintNNLO_s_aq_s12(1) = zero

    else

       !-- single boost
       !-- aa
       call res_treeAA_aa(LOProc%AmpMom,res_lo_1(1,1))
       call get_respdf_hoppet(PDFs,xPij_a,aa_lumi,1,1,LOProc,res_lo_1,respdf_1,&
         myPDFs2_Lmu=[xPij_a_Lmu,xPij_a_Lmu2])

       !-- electric charges already inside the pdfs
       respdf_1 = CF * respdf_1 * LOProc%wgt

       !-- qqb
       call res_tree_qqb(LOProc%AmpMom,res_lo_2)
       res_lo_2(:,1) = Qdn2 * res_lo_2(:,1)
       res_lo_2(:,2) = Qup2 * res_lo_2(:,2)
       call get_respdf_hoppet(xPij,PDFs,ns_lumi,1,1,LOProc,res_lo_2,respdf_2,&
         myPDFs1_Lmu=[xPij_Lmu,xPij_Lmu2])
       respdf_2 = CF * xn * respdf_2 * LOProc%wgt

       !-- double boost
       call get_respdf_hoppet(xPij_A_1,xPij_A_2,ns_lumi,1,1,LOProc,res_lo_2,respdf_tmp(:,1),&
          myPDFs1_Lmu=[xPij_A_1_Lmu],myPDFs2_Lmu=[xPij_A_2_Lmu])

       call get_respdf_hoppet(xPij_B_1,xPij_B_2,ns_lumi,1,1,LOProc,res_lo_2,respdf_tmp(:,2),&
          myPDFs1_Lmu=[xPij_B_1_Lmu])

       call get_respdf_hoppet(xPij_C_1,xPij_C_2,ns_lumi,1,1,LOProc,res_lo_2,respdf_tmp(:,3),&
          myPDFs2_Lmu=[xPij_C_2_Lmu])

       call get_respdf_hoppet(xPij_D_1,xPij_D_2,ns_lumi,1,1,LOProc,res_lo_2,respdf_tmp(:,4))

       respdf_3 = sum(respdf_tmp,2)
       respdf_3 = CF * xn * respdf_3 * LOProc%wgt

       !-- total
       respdf = respdf_1 + respdf_2 + respdf_3

       FintNNLO_s_aq_s12 = respdf(1)
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_aq = FintNNLO_s_aq_s12
#endif

  end function xsect_nnlo_s_aq_s12

  function xsect_nnlo_s_aq_oqcd(yRnd,ff,vegasweight)
    use mod_kinematics_nlo_z
    implicit none
    integer :: xsect_nnlo_s_aq_oqcd
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim,S5Lim,S5C1Lim,S5C2Lim
    !--
    real(dp) :: xx(kNLO_max_full),z
    real(dp) :: kin(4),FintNNLO_s_aq_oqcd(6)
    real(dp) :: EC,damp,eta51,eta52,s15,s25,z5,eik_qcd,E5
    real(dp) :: respdf(ipdf),intsub(ipdf),sv_logs(ipdf),PqgNLO_v(ipdf)
    real(dp) :: respdf_bak(ipdf),respdf_1(ipdf),respdf_2(ipdf),respdf_3(ipdf)
    real(dp) :: Pqg0,PqgNLO,res_nlo(2,2),res_lo(2,2)

    xsect_nnlo_s_aq_oqcd = 0

    ff(1) = zero
    kin = 0

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z = buff+onet*real(yRnd(kNLO_max_full),dp)

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_r)) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()
    call kinematics_nlo_z_is(xx,z,one,HardProc=HardProc,&
      C1Lim=C1Lim,C2Lim=C2Lim,SLim=S5Lim,SC1Lim=S5C1Lim,SC2Lim=S5C2Lim,&
      compute_etas=.false.)

    Pqg0 = Pgq_spav(z)
    PqgNLO = 2*((one-z)*z + Pqg0*log(one-z))
    PqgNLO_v = PqgNLO

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_g]
    call cut_histo(HardProc)

    if (HardProc%makecut) then

       FintNNLO_s_aq_oqcd(1) = zero
       kin(1) = zero

    else

       call res_tree_g_qqb(HardProc%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(aq_lumi,1,1,HardProc,res_nlo,respdf)

       EC    = HardProc%Lim_Ei(1)
       eta51 = HardProc%Lim_etaij(1,5)
       eta52 = HardProc%Lim_etaij(2,5)
       damp = eta52 * (one + eta51)

       call fill_sv_logs(HardProc%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqgNLO_v + Pqg0 * log(eta51/2) * damp
       intsub = intsub - Pqg0*sv_logs
       intsub = xn * intsub

       respdf = intsub * respdf * HardProc%wgt

       FintNNLO_s_aq_oqcd(1) = respdf(1)
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 51                           !!
    !!-----------------------------------------------------------------------!!
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C1Lim)

    if (C1Lim%makecut) then

       FintNNLO_s_aq_oqcd(2) = zero
       kin(2) = zero

    else

       call res_tree_qqb(C1Lim%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       call get_respdf(aq_lumi,1,1,C1Lim,res_lo,respdf)

       !-- use etas from hard process
       eta51 = C1Lim%Lim_KinInv(4)
       EC    = C1Lim%Lim_Ei(1)
       s15   = C1Lim%Lim_sij(1,5)
       z5    = C1Lim%Lim_z(1)

       call fill_sv_logs(C1Lim%muf(1)**2,4*EC**2,sv_logs)

       !-- damp = one
       intsub = PqgNLO_v + Pqg0 * log(eta51/2)
       intsub = intsub - Pqg0*sv_logs
       intsub = xn * intsub

       respdf = intsub * respdf &
              * CF * 2/s15 * Pqg(z5)/(one-z5) &
              * C1Lim%wgt
       respdf = - respdf

       FintNNLO_s_aq_oqcd(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim)

    if (C2Lim%makecut) then

       FintNNLO_s_aq_oqcd(3) = zero
       kin(3) = zero

    else

       call res_tree_qqb(C2Lim%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       call get_respdf(aq_lumi,1,1,C2Lim,res_lo,respdf)

       EC    = C2Lim%Lim_Ei(1)
       s25   = C2Lim%Lim_sij(2,5)
       z5    = C2Lim%Lim_z(2)

       call fill_sv_logs(C2Lim%muf(1)**2,4*EC**2,sv_logs)

       !-- damp = zero
       intsub = PqgNLO_v
       intsub = intsub - Pqg0*sv_logs
       intsub = xn * intsub

       respdf = intsub * respdf &
              * CF * 2/s25 * Pqg(z5)/(one-z5) &
              * C2Lim%wgt
       respdf = - respdf

       FintNNLO_s_aq_oqcd(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                             Soft Limit 5                              !!
    !!-----------------------------------------------------------------------!!
    S5Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(S5Lim)

    if (S5Lim%makecut) then

       kin(4) = zero
       FintNNLO_s_aq_oqcd(4:6) = zero

    else

       !! ----------------------------- S5Lim ------------------------------ !!
       call res_tree_qqb(S5Lim%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       call get_respdf(aq_lumi,1,1,S5Lim,res_lo,respdf_bak)

       call get_qcd_eik(CF,S5Lim%Lim_etaij,[1,2],5,eik_qcd)
       EC    = S5Lim%Lim_Ei(1)
       E5    = S5Lim%Lim_Ei(2)
       eta51 = S5Lim%Lim_etaij(1,5)
       eta52 = S5Lim%Lim_etaij(2,5)
       damp  = eta52 * (one + eta51)

       call fill_sv_logs(S5Lim%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqgNLO_v + Pqg0 * log(eta51/2) * damp
       intsub = intsub - Pqg0*sv_logs
       intsub = xn * intsub

       respdf_1 = intsub * respdf_bak &
                * eik_qcd/E5**2 &
                * S5Lim%wgt
       respdf_1 = - respdf_1

       FintNNLO_s_aq_oqcd(4) = respdf_1(1)

       !! ---------------------------- S5C1Lim ----------------------------- !!

       !-- use etas from hard process
       eta51 = S5C1Lim%Lim_KinInv(4)

       !-- damp = one
       intsub = PqgNLO_v + Pqg0 * log(eta51/2)
       intsub = intsub - Pqg0*sv_logs
       intsub = xn * intsub

       respdf_2 = intsub * respdf_bak &
                * CF/eta51/E5**2 &
                * S5C1Lim%wgt

       FintNNLO_s_aq_oqcd(5) = respdf_2(1)

       !! ---------------------------- S5C2Lim ----------------------------- !!

       !-- use etas from hard process
       eta52 = S5C2Lim%Lim_KinInv(4)

       !-- damp = zero
       intsub = PqgNLO
       intsub = intsub - Pqg0*sv_logs
       intsub = xn * intsub

       respdf_3 = intsub * respdf_bak &
                * CF/eta52/E5**2 &
                * S5C2Lim%wgt

       FintNNLO_s_aq_oqcd(6) = respdf_3(1)

       !! --------------------------- S5Lim tot ---------------------------- !!

       respdf = respdf_1 + respdf_2 + respdf_3
       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_aq_onlo = FintNNLO_s_aq_oqcd
#endif

  end function xsect_nnlo_s_aq_oqcd

  function xsect_nnlo_s_aq_oewk(yRnd,ff,vegasweight)
    use mod_kinematics_nlo_z
    implicit none
    integer :: xsect_nnlo_s_aq_oewk
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    type(KinConfig) :: HardProc_z,C1Lim_z,C2Lim_z
    !--
    real(dp) :: xx(kNLO_max_full),z
    real(dp) :: kin(6),FintNNLO_s_aq_oewk(6)
    real(dp) :: EC,E5,damp,damp1,damp2,eta51,eta52,s15,s25,z5
    real(dp) :: respdf(ipdf),intsub(ipdf),intsub_el(ipdf),sv_logs(ipdf)
    real(dp) :: Pqq0_R(-1:1),PqqNLO(-1:1),res_nlo(2,2),res_lo(2,2),res_aa

    xsect_nnlo_s_aq_oewk = 0

    ff(1) = zero
    kin = 0

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z = buff+onet*real(yRnd(kNLO_max_full),dp)

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_r)) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !-- Boosted kinematics
    call kinematics_nlo_z_is(xx,one,z,HardProc=HardProc_z,&
      C1Lim=C1Lim_z,C2Lim=C2Lim_z,compute_etas=.false.)

    !-- z=1 kinematics
    call kinematics_nlo_z_is(xx,one,one,HardProc=HardProc,&
      C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.false.)

    !-- Splittings and integrated sub terms
    Pqq0_R = PqqAP_0_R(z)
    PqqNLO = Pqq_NLO(z)

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !!                         FLM[1_a,z.2_q,3,4|5_q]                        !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc_z%ids(1:5) = [0,0,id_el,-id_el,id_q]
    call cut_histo(HardProc_z)

    if (HardProc_z%makecut) then

       FintNNLO_s_aq_oewk(1) = zero
       kin(1) = zero

    else

       call res_tree_a_aq(HardProc_z%AmpMom,res_nlo)

       call get_respdf(aq_lumi,1,1,HardProc_z,res_nlo,respdf)

       EC    = HardProc_z%Lim_Ei(1)
       eta51 = HardProc_z%Lim_etaij(1,5)
       eta52 = HardProc_z%Lim_etaij(2,5)
       damp = eta51 * (one + eta52)

       call fill_sv_logs(HardProc_z%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta52/2) * damp
       intsub = intsub + Pqq0_R(0)*sv_logs
       intsub = CF * intsub

       respdf = intsub * respdf * HardProc_z%wgt

       FintNNLO_s_aq_oewk(1) = respdf(1)
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 51                           !!
    !!-----------------------------------------------------------------------!!
    C1Lim_z%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C1Lim_z)

    if (C1Lim_z%makecut) then

       FintNNLO_s_aq_oewk(2) = zero
       kin(2) = zero

    else

       call res_tree_qqb(C1Lim_z%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       call get_respdf(aq_lumi,1,1,C1Lim_z,res_lo,respdf)

       EC  = C1Lim_z%Lim_Ei(1)
       s15 = C1Lim_z%Lim_sij(1,5)
       z5  = C1Lim_z%Lim_z(1)

       call fill_sv_logs(C1Lim_z%muf(1)**2,4*EC**2,sv_logs)

       !-- damp = zero
       intsub = PqqNLO(0)
       intsub = intsub + Pqq0_R(0)*sv_logs
       intsub = CF * intsub

       respdf = intsub * respdf &
              * xn * 2/s15 * Pgq_spav(z5)/(one-z5) &
              * C1Lim_z%wgt
       respdf = - respdf

       FintNNLO_s_aq_oewk(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim_z%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim_z)

    if (C2Lim_z%makecut) then

       FintNNLO_s_aq_oewk(3) = zero
       kin(3) = zero

    else

       call res_treeAA_aa(C2Lim_z%AmpMom,res_aa)
       res_lo(:,1) = Qdn2 * res_aa
       res_lo(:,2) = Qup2 * res_aa

       call get_respdf(aq_lumi,1,1,C2Lim_z,res_lo,respdf)

       eta52 = C2Lim_z%Lim_KinInv(4)
       EC    = C2Lim_z%Lim_Ei(1)
       s25   = C2Lim_z%Lim_sij(2,5)
       z5    = C2Lim_z%Lim_z(2)

       call fill_sv_logs(C2Lim_z%muf(1)**2,4*EC**2,sv_logs)

       !-- damp = one
       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta52/2)
       intsub = intsub + Pqq0_R(0)*sv_logs
       intsub = CF * intsub

       respdf = intsub * respdf &
              * 2/s25 * Pqq(z5)/(one-z5) &
              * C2Lim_z%wgt
       respdf = - respdf

       FintNNLO_s_aq_oewk(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !!                          FLM[1_a,2_q,3,4|5_q]                         !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    call cut_histo(HardProc)

    if (HardProc%makecut) then

       FintNNLO_s_aq_oewk(4) = zero
       kin(4) = zero

    else

       call res_tree_a_aq(HardProc%AmpMom,res_nlo)

       call get_respdf(aq_lumi,1,1,HardProc,res_nlo,respdf)

       EC    = HardProc%Lim_Ei(1)
       eta51 = HardProc%Lim_etaij(1,5)
       eta52 = HardProc%Lim_etaij(2,5)
       damp = eta51 * (one + eta52)

       call fill_sv_logs(HardProc%muf(1)**2,4*EC**2,sv_logs)

       !-- subtract plus
       intsub = PqqNLO(1) - Pqq0_R(1) * log(eta52/2) * damp
       intsub = intsub + Pqq0_R(1)*sv_logs
       intsub = CF * intsub
       intsub = -intsub

       !-- add elastic component
       call fill_sv_logs(HardProc%muf(1)**2,EC**2,sv_logs)

       E5    = HardProc%Lim_Ei(2)
       damp1 = eta52**2 * (3*eta51+eta52)
       damp2 = eta51**2 * (3*eta52+eta51)
       intsub_el = CF*(calG_CF(EC,E5,eta51,eta52,damp1,damp2) - 1.5_dp*sv_logs)

       respdf = (intsub + intsub_el) * respdf * HardProc%wgt

       FintNNLO_s_aq_oewk(4) = respdf(1)
       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 51                           !!
    !!-----------------------------------------------------------------------!!
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C1Lim)

    if (C1Lim%makecut) then

       FintNNLO_s_aq_oewk(5) = zero
       kin(5) = zero

    else

       call res_tree_qqb(C1Lim%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       call get_respdf(aq_lumi,1,1,C1Lim,res_lo,respdf)

       EC  = C1Lim%Lim_Ei(1)
       s15 = C1Lim%Lim_sij(1,5)
       z5  = C1Lim%Lim_z(1)

       call fill_sv_logs(C1Lim%muf(1)**2,4*EC**2,sv_logs)

       !-- subtract plus, damp = zero
       intsub = PqqNLO(1)
       intsub = intsub + Pqq0_R(1)*sv_logs
       intsub = CF * intsub
       intsub = -intsub

       !-- add elastic component
       call fill_sv_logs(C1Lim%muf(1)**2,EC**2,sv_logs)

       E5    = C1Lim%Lim_Ei(2)
       eta51 = C1Lim%Lim_KinInv(4)
       eta52 = C1Lim%Lim_etaij(2,5)
       intsub_el = CF*(calG_CF(EC,E5,eta51,eta52,1) - 1.5_dp*sv_logs)

       respdf = (intsub + intsub_el) * respdf &
              * xn * 2/s15 * Pgq_spav(z5)/(one-z5) &
              * C1Lim%wgt
       respdf = - respdf

       FintNNLO_s_aq_oewk(5) = respdf(1)
       kin(5) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim)

    if (C2Lim%makecut) then

       FintNNLO_s_aq_oewk(6) = zero
       kin(6) = zero

    else

       call res_treeAA_aa(C2Lim%AmpMom,res_aa)
       res_lo(:,1) = Qdn2 * res_aa
       res_lo(:,2) = Qup2 * res_aa

       call get_respdf(aq_lumi,1,1,C2Lim,res_lo,respdf)

       eta52 = C2Lim%Lim_KinInv(4)
       EC    = C2Lim%Lim_Ei(1)
       s25   = C2Lim%Lim_sij(2,5)
       z5    = C2Lim%Lim_z(2)

       call fill_sv_logs(C2Lim%muf(1)**2,4*EC**2,sv_logs)

       !-- subtract plus, damp = one
       intsub = PqqNLO(1) - Pqq0_R(1) * log(eta52/2)
       intsub = intsub + Pqq0_R(1)*sv_logs
       intsub = CF * intsub
       intsub = - intsub

       !-- add elastic component
       call fill_sv_logs(C2Lim%muf(1)**2,EC**2,sv_logs)

       E5    = C2Lim%Lim_Ei(2)
       eta51 = C2Lim%Lim_etaij(1,5)
       intsub_el = CF*(calG_CF(EC,E5,eta51,eta52,2) - 1.5_dp*sv_logs)

       respdf = (intsub + intsub_el) * respdf &
              * 2/s25 * Pqq(z5)/(one-z5) &
              * C2Lim%wgt
       respdf = - respdf

       FintNNLO_s_aq_oewk(6) = respdf(1)
       kin(6) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_aq_onlo = FintNNLO_s_aq_oewk
#endif

  end function xsect_nnlo_s_aq_oewk

end module mod_xsects_nnlo_s_aq

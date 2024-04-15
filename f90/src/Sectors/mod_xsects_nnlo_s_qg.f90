module mod_xsects_nnlo_s_qg
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
  use mod_splittings_bare
  use mod_eikonals
  use mod_partitions
  use mod_hoppet_tools
  use mod_hoppet_nnlo
  use mod_int_sub_nnlo
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_s_qg(2),FintNNLO_s_qg_onlo(6)
#endif

  public :: xsect_nnlo_s_qg_s12
  public :: xsect_nnlo_s_qg_vewk
  public :: xsect_nnlo_s_qg_oqcd
  public :: xsect_nnlo_s_qg_oewk_is
  public :: xsect_nnlo_s_qg_oewk_fs_53,xsect_nnlo_s_qg_oewk_fs_54

  public :: xsect_nnlo_s_qg_vewknf

  real(dp), parameter :: charges(4,4) = reshape(& 
     [-Qdn,   Qdn,  -Qup,   Qup,   &
       Qdn,  -Qdn,   Qup,  -Qup,   &
       Q_lep, Q_lep, Q_lep, Q_lep, &
      -Q_lep,-Q_lep,-Q_lep,-Q_lep] &
      ,[4,4])

  !-- positions of charged quarks and leptons (ql)
  integer, parameter :: ql_pos(4) = [1,2,3,4]

  private

contains

  !!*************************************************************************!!

  function xsect_nnlo_s_qg_oewk_fs_53(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_s_qg_oewk_fs_53
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nnlo_s_qg_oewk_fs_53 = xsect_nnlo_s_qg_oewk_fs_5i(yRnd,ff,vegasweight,3,4)

  end function xsect_nnlo_s_qg_oewk_fs_53

  function xsect_nnlo_s_qg_oewk_fs_54(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_s_qg_oewk_fs_54
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nnlo_s_qg_oewk_fs_54 = xsect_nnlo_s_qg_oewk_fs_5i(yRnd,ff,vegasweight,4,3)

  end function xsect_nnlo_s_qg_oewk_fs_54


  !!*************************************************************************!!

  function xsect_nnlo_s_qg_s12(yRnd,ff,vegasweight)
    use mod_kinematics_lo
    use mod_kinematics_lo_z
    integer :: xsect_nnlo_s_qg_s12
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc,LOProc_z
    !--
    real(dp) :: xx(kLO_max_full),z
    real(dp) :: respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf),respdf_tmp(ipdf,4)
    real(dp) :: PqgNLO(ipdf),tildePqgNLO(ipdf),sv_logs(ipdf)
    real(dp) :: EC
    real(dp) :: kin(2),FintNNLO_s_qg_s12(2)
    real(dp) :: res_lo(2,2),Pqg0,qg_int_sub(2,2)

    xsect_nnlo_s_qg_s12 = 0

    ff(1) = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    z = buff+onet*real(yRnd(kLO_max_full),dp)

#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()
    !-- z = 1 kinematics
    call kinematics_lo(xx,LOProc)
    !-- boosted kinematics
    call kinematics_lo_z(xx,one,z,LOProc_z)

    Pqg0 = Pgq_spav(z)
    PqgNLO = 2*((one-z)*z + Pqg0*log(one-z))

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(LOProc)

    if (LOProc%makecut) then

       kin(1) = zero
       FintNNLO_s_qg_s12(1) = zero

    else

       call res_tree_qqb(LOProc%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       !-- F[z1.p1,z2.p2]
       call get_respdf_hoppet(xPij_A_2,xPij_A_1,ns_lumi,1,1,LOProc,res_lo,respdf_tmp(:,1),&
          myPDFs1_Lmu=[xPij_A_2_Lmu],myPDFs2_Lmu=[xPij_A_1_Lmu])

       call get_respdf_hoppet(xPij_B_2,xPij_B_1,ns_lumi,1,1,LOProc,res_lo,respdf_tmp(:,2),&
          myPDFs2_Lmu=[xPij_B_1_Lmu])

       call get_respdf_hoppet(xPij_C_2,xPij_C_1,ns_lumi,1,1,LOProc,res_lo,respdf_tmp(:,3),&
          myPDFs1_Lmu=[xPij_C_2_Lmu])

       call get_respdf_hoppet(xPij_D_2,xPij_D_1,ns_lumi,1,1,LOProc,res_lo,respdf_tmp(:,4))

       respdf_1 = sum(respdf_tmp,2)

       !-- F[p1,z.p2]
       call get_respdf_hoppet(PDFs,xPij,ns_lumi,1,1,LOProc,res_lo,respdf_2,&
            myPDFs2_Lmu=[xPij_Lmu,xPij_Lmu2])

       respdf = Tr*(respdf_1 + respdf_2)*LOProc%wgt

       FintNNLO_s_qg_s12(1) = respdf(1)
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !! ----------------------- Soft-photon component ----------------------- !!

    LOProc_z%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(LOProc_z)

    if (LOProc_z%makecut) then

       kin(2) = zero
       FintNNLO_s_qg_s12(2) = zero

    else

       call res_tree_qqb(LOProc_z%AmpMom,res_lo)

       qg_int_sub = reshape(calG_QqQl(LOProc_z,z,[Qdn,-Qdn, Qup,-Qup],Q_lep,1),[2,2])

       res_lo(:,1) = qg_int_sub(:,1) * res_lo(:,1)
       res_lo(:,2) = qg_int_sub(:,2) * res_lo(:,2)

       call get_respdf(qg_lumi,1,1,LOProc_z,res_lo,respdf)

       EC = LOProc_z%Lim_Ei(1)
       call fill_sv_logs(LOProc_z%muf(1)**2,4*EC**2,sv_logs)

       tildePqgNLO = Tr * (PqgNLO - Pqg0 * sv_logs)

       respdf = tildePqgNLO * respdf * LOProc_z%wgt

       FintNNLO_s_qg_s12(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_qg = FintNNLO_s_qg_s12
#endif

  end function xsect_nnlo_s_qg_s12


  function xsect_nnlo_s_qg_vewk(yRnd,ff,vegasweight)
    use mod_kinematics_lo
    integer :: xsect_nnlo_s_qg_vewk
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: respdf(ipdf)
    real(dp) :: res_tree(2,3),res_nlo(2,3),kin(1)

    xsect_nnlo_s_qg_vewk = 0

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

       call res_ewkloop_qqb(LOProc%AmpMom,res_tree,res_nlo)

       !-- F[1,z.2]
       call get_respdf_hoppet(PDFs,xPij,ns_lumi_splitb,1,1,LOProc,res_nlo,respdf,myPDFs2_Lmu=[xPij_Lmu])
       respdf = Tr * respdf * LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_qg(1) = kin(1)
#endif

  end function xsect_nnlo_s_qg_vewk

  function xsect_nnlo_s_qg_oqcd(yRnd,ff,vegasweight)
    use mod_kinematics_nlo_z
    implicit none
    integer :: xsect_nnlo_s_qg_oqcd
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C2Lim
    type(KinConfig) :: HardProc_z,C2Lim_z
    !--
    real(dp) :: xx(kNLO_max_full),z
    real(dp) :: kin(4),FintNNLO_s_qg_oqcd(4)
    real(dp) :: EC,s25,z5
    real(dp) :: respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf),respdf_tmp(ipdf)
    real(dp) :: intsub(ipdf),sv_logs(ipdf)
    real(dp) :: Pqq0_R(-1:1),PqqNLO(-1:1),qg_int_sub(2,2)
    real(dp) :: res_nlo(2,2),res_nlo_1(2,2),res_nlo_2(2,2),res_nlo_tmp(2,2),res_lo(2,2)

    xsect_nnlo_s_qg_oqcd = 0

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
    call kinematics_nlo_z_is(xx,z,one,HardProc=HardProc_z,C2Lim=C2Lim_z,&
      compute_etas=.false.)

    !-- z=1 kinematics
    call kinematics_nlo_z_is(xx,one,one,HardProc=HardProc,C2Lim=C2Lim,&
      compute_etas=.true.)

    !-- Splittings and integrated sub terms
    Pqq0_R = PqqAP_0_R(z)
    PqqNLO = Pqq_NLO(z)

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !!                         FLM[z.1_q,2_g,3,4|5_q]                        !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc_z%ids(1:5) = [0,0,id_el,-id_el,id_q]
    call cut_histo(HardProc_z)

    if (HardProc_z%makecut) then

       FintNNLO_s_qg_oqcd(1) = zero
       kin(1) = zero

    else

       call res_tree_g_qg(HardProc_z%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(qg_lumi,1,1,HardProc_z,res_nlo,respdf)

       EC = HardProc_z%Lim_Ei(1)
       call fill_sv_logs(HardProc_z%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) + Pqq0_R(0)*sv_logs

       respdf = intsub * respdf * HardProc_z%wgt

       FintNNLO_s_qg_oqcd(1) = respdf(1)
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim_z%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim_z)

    if (C2Lim_z%makecut) then

       FintNNLO_s_qg_oqcd(2) = zero
       kin(2) = zero

    else

       call res_tree_qqb(C2Lim_z%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       call get_respdf(qg_lumi,1,1,C2Lim_z,res_lo,respdf)

       EC  = C2Lim_z%Lim_Ei(1)
       s25 = C2Lim_z%Lim_sij(2,5)
       z5  = C2Lim_z%Lim_z(2)

       call fill_sv_logs(C2Lim_z%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) + Pqq0_R(0)*sv_logs

       respdf = intsub * respdf * C2Lim_z%wgt &
              * Tr * 2/s25 * Pgq_spav(z5)/(one-z5)
       respdf = - respdf

       FintNNLO_s_qg_oqcd(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !!                          FLM[1_q,2_g,3,4|5_q]                         !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    call cut_histo(HardProc)

    if (HardProc%makecut) then

       FintNNLO_s_qg_oqcd(3) = zero
       kin(3) = zero

    else

       call res_tree_g_qg(HardProc%AmpMom,res_nlo_tmp)
       res_nlo_1(:,1) = Qdn2 * res_nlo_tmp(:,1)
       res_nlo_1(:,2) = Qup2 * res_nlo_tmp(:,2)

       !-- subtract plus
       EC = HardProc%Lim_Ei(1)
       call fill_sv_logs(HardProc%muf(1)**2,4*EC**2,sv_logs)

       call get_respdf(qg_lumi,1,1,HardProc,res_nlo_1,respdf_tmp)

       intsub = PqqNLO(1) + Pqq0_R(1)*sv_logs
       intsub = -intsub

       respdf_1 = intsub * respdf_tmp 

       !-- add elastic component
       qg_int_sub = reshape(calG_Q2(HardProc,[-Qdn,Qdn, -Qup,Qup],Q_lep,2,1,.false.),[2,2])

       res_nlo_2(:,1) = qg_int_sub(:,1) * res_nlo_tmp(:,1)
       res_nlo_2(:,2) = qg_int_sub(:,2) * res_nlo_tmp(:,2)

       call get_respdf(qg_lumi,1,1,HardProc,res_nlo_2,respdf_2)

       !-- scale variation of the elastic bit
       call fill_sv_logs(HardProc%muf(1)**2,EC**2,sv_logs)
       respdf_2 = respdf_2 - 1.5_dp * respdf_tmp * sv_logs

       respdf = (respdf_1 + respdf_2) * HardProc%wgt

       FintNNLO_s_qg_oqcd(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim)

    if (C2Lim%makecut) then

       FintNNLO_s_qg_oqcd(4) = zero
       kin(4) = zero

    else

       call res_tree_qqb(C2Lim%AmpMom,res_nlo_tmp)
       res_nlo_1(:,1) = Qdn2 * res_nlo_tmp(:,1)
       res_nlo_1(:,2) = Qup2 * res_nlo_tmp(:,2)

       !-- subtract plus
       EC = C2Lim%Lim_Ei(1)
       call fill_sv_logs(C2Lim%muf(1)**2,4*EC**2,sv_logs)

       call get_respdf(qg_lumi,1,1,C2Lim,res_nlo_1,respdf_tmp)

       intsub = PqqNLO(1) + Pqq0_R(1)*sv_logs
       intsub = -intsub

       respdf_1 = intsub * respdf_tmp 

       !-- add elastic component
       C2Lim%Lim_etaij(2,5) = C2Lim%Lim_KinInv(4)
       qg_int_sub = reshape(calG_Q2(C2Lim,[-Qdn,Qdn, -Qup,Qup],Q_lep,2,1,.true.),[2,2])

       res_nlo_2(:,1) = qg_int_sub(:,1) * res_nlo_tmp(:,1)
       res_nlo_2(:,2) = qg_int_sub(:,2) * res_nlo_tmp(:,2)

       call get_respdf(qg_lumi,1,1,C2Lim,res_nlo_2,respdf_2)

       !-- scale variation of the elastic bit
       call fill_sv_logs(C2Lim%muf(1)**2,EC**2,sv_logs)
       respdf_2 = respdf_2 - 1.5_dp * respdf_tmp * sv_logs

       s25 = C2Lim%Lim_sij(2,5)
       z5  = C2Lim%Lim_z(2)

       respdf = (respdf_1 + respdf_2) * C2Lim%wgt &
              * Tr * 2/s25 * Pgq_spav(z5)/(one-z5)
       respdf = - respdf

       FintNNLO_s_qg_oqcd(4) = respdf(1)
       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_qg_onlo(1:4) = FintNNLO_s_qg_oqcd
#endif

  end function xsect_nnlo_s_qg_oqcd


  function xsect_nnlo_s_qg_oewk_is(yRnd,ff,vegasweight)
    use mod_kinematics_nlo_z
    implicit none
    integer :: xsect_nnlo_s_qg_oewk_is
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim,S5Lim,S5C1Lim,S5C2Lim
    !--
    real(dp) :: xx(kNLO_max_full),z
    real(dp) :: kin(4),FintNNLO_s_qg_oewk_is(6)
    real(dp) :: EC,damp,damp_ewk,eta51,eta52,s15,s25,z5
    real(dp) :: eik_qed(4),E5,E5sq
    real(dp) :: respdf(ipdf),intsub(ipdf),sv_logs(ipdf),PqgNLO(ipdf),Pqg0
    real(dp) :: respdf_bak(ipdf),respdf_1(ipdf),respdf_2(ipdf),respdf_3(ipdf)
    real(dp) :: res_nlo(2,2),res_lo(2,2),res_lo_tmp(2,2)

    xsect_nnlo_s_qg_oewk_is = 0

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
    call kinematics_nlo_z_is(xx,one,z,HardProc=HardProc,&
      C1Lim=C1Lim,C2Lim=C2Lim,SLim=S5Lim,SC1Lim=S5C1Lim,SC2Lim=S5C2Lim,&
      compute_etas=.true.)

    Pqg0 = Pgq_spav(z)
    PqgNLO = 2*((one-z)*z + Pqg0*log(one-z))

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    call cut_histo(HardProc)

    if (HardProc%makecut) then

       FintNNLO_s_qg_oewk_is(1) = zero
       kin(1) = zero

    else

       call res_tree_a_qqb(HardProc%AmpMom,res_nlo)
       call get_respdf(qg_lumi,1,1,HardProc,res_nlo,respdf)

       call partition_nnlo_fact(HardProc,damp,iconf_qed=[1,2,3,4,5],i_qed=2)

       EC = HardProc%Lim_Ei(1)
       call fill_sv_logs(HardProc%muf(1)**2,4*EC**2,sv_logs)

       eta52 = HardProc%Lim_etaij(2,5)
       intsub = PqgNLO + Pqg0 * log(eta52/2) * damp
       intsub = Tr * (intsub - Pqg0*sv_logs)

       call partition_nlo_qed(HardProc,damp_ewk,12)
       respdf = intsub * respdf * HardProc%wgt * damp_ewk

       FintNNLO_s_qg_oewk_is(1) = respdf(1)
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 51                           !!
    !!-----------------------------------------------------------------------!!
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C1Lim)

    if (C1Lim%makecut) then

       FintNNLO_s_qg_oewk_is(2) = zero
       kin(2) = zero

    else

       call res_tree_qqb(C1Lim%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       call get_respdf(qg_lumi,1,1,C1Lim,res_lo,respdf)

       s15 = C1Lim%Lim_sij(1,5)
       z5  = C1Lim%Lim_z(1)
       EC  = C1Lim%Lim_Ei(1)
       call fill_sv_logs(C1Lim%muf(1)**2,4*EC**2,sv_logs)

       !-- damp = zero
       intsub = PqgNLO
       intsub = Tr * (intsub - Pqg0*sv_logs)

       !-- damp_ewk = one
       respdf = intsub * respdf * C1Lim%wgt &
              * 2/s15 * Pqg(z5)/(one-z5)
       respdf = - respdf

       FintNNLO_s_qg_oewk_is(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim)

    if (C2Lim%makecut) then

       FintNNLO_s_qg_oewk_is(3) = zero
       kin(3) = zero

    else

       call res_tree_qqb(C2Lim%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       call get_respdf(qg_lumi,1,1,C2Lim,res_lo,respdf)

       !-- use etas from hard process
       eta52 = C2Lim%Lim_KinInv(4)
       EC    = C2Lim%Lim_Ei(1)
       s25   = C2Lim%Lim_sij(2,5)
       z5    = C2Lim%Lim_z(2)

       call fill_sv_logs(C2Lim%muf(1)**2,4*EC**2,sv_logs)

       !-- damp = one
       intsub = PqgNLO + Pqg0 * log(eta52/2)
       intsub = Tr * (intsub - Pqg0*sv_logs)

       !-- damp_ewk = one
       respdf = intsub * respdf &
              * 2/s25 * Pqg(z5)/(one-z5) &
              * C2Lim%wgt
       respdf = - respdf

       FintNNLO_s_qg_oewk_is(3) = respdf(1)
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
       FintNNLO_s_qg_oewk_is(4:6) = zero

    else

       !! ----------------------------- S5Lim ------------------------------ !!
       call res_tree_qqb(S5Lim%AmpMom,res_lo_tmp)

       E5   = S5Lim%Lim_Ei(2)
       E5sq = E5*E5
       call get_qed_eik(charges,S5Lim%Lim_etaij,ql_pos,5,eik_qed)
       eik_qed = eik_qed/E5sq

       res_lo(1:2,1) = res_lo_tmp(1:2,1) * eik_qed(1:2)
       res_lo(1:2,2) = res_lo_tmp(1:2,2) * eik_qed(3:4)

       call get_respdf(qg_lumi,1,1,S5Lim,res_lo,respdf_bak)

       EC    = S5Lim%Lim_Ei(1)
       call fill_sv_logs(S5Lim%muf(1)**2,4*EC**2,sv_logs)
       eta52 = S5Lim%Lim_etaij(2,5)

       call partition_nnlo_fact(S5Lim,damp,iconf_qed=[1,2,3,4,5],i_qed=2)
       intsub = PqgNLO + Pqg0 * log(eta52/2) * damp
       intsub = Tr * (intsub - Pqg0*sv_logs)

       call partition_nlo_qed(S5Lim,damp_ewk,12)
       respdf_1 = intsub * respdf_bak * S5Lim%wgt * damp_ewk
       respdf_1 = - respdf_1

       FintNNLO_s_qg_oewk_is(4) = respdf_1(1)

       !! ---------------------------- S5C1Lim ----------------------------- !!
       eta51 = S5C1Lim%Lim_KinInv(4)

       res_lo(:,1) = Qdn2 * res_lo_tmp(:,1)
       res_lo(:,2) = Qup2 * res_lo_tmp(:,2)

       call get_respdf(qg_lumi,1,1,S5Lim,res_lo,respdf_bak)

       !-- damp = zero
       intsub = PqgNLO
       intsub = Tr * (intsub - Pqg0*sv_logs)

       !-- damp_ewk = one
       respdf_2 = intsub * respdf_bak * S5C1Lim%wgt &
                * one/eta51/E5sq

       FintNNLO_s_qg_oewk_is(5) = respdf_2(1)

       !! ---------------------------- S5C2Lim ----------------------------- !!
       eta52 = S5C2Lim%Lim_KinInv(4)

       !-- damp = one
       intsub = PqgNLO + Pqg0 * log(eta52/2)
       intsub = Tr * (intsub - Pqg0*sv_logs)

       !-- damp_ewk = one
       respdf_3 = intsub * respdf_bak * S5C2Lim%wgt &
                * one/eta52/E5sq

       FintNNLO_s_qg_oewk_is(6) = respdf_3(1)

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
    FintNNLO_s_qg_onlo = FintNNLO_s_qg_oewk_is
#endif

  end function xsect_nnlo_s_qg_oewk_is


  function xsect_nnlo_s_qg_oewk_fs_5i(yRnd,ff,vegasweight,icoll,jother)
    use mod_kinematics_nlo_z
    implicit none
    integer :: xsect_nnlo_s_qg_oewk_fs_5i,icoll,jother
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,CLim,SLim,SCLim
    !--
    real(dp) :: xx(kNLO_max_full),z
    real(dp) :: kin(3),FintNNLO_s_qg_oewk_fs(4)
    real(dp) :: EC,damp,damp_ewk,eta52,eta5i,s5i,z5
    real(dp) :: eik_qed(4),E5,E5sq
    real(dp) :: respdf(ipdf),intsub(ipdf),sv_logs(ipdf),PqgNLO(ipdf),Pqg0
    real(dp) :: respdf_tmp(ipdf,2),res_tmp_vect(2,2,2),respdf_vect(2,ipdf)
    real(dp) :: res_nlo(2,2),res_lo(2,2)

    xsect_nnlo_s_qg_oewk_fs_5i = 0

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
    call kinematics_nlo_z_fs(xx,icoll,jother,one,z,HardProc=HardProc,&
    CLim=CLim,SLim=SLim,SCLim=SCLim)

    Pqg0 = Pgq_spav(z)
    PqgNLO = 2*((one-z)*z + Pqg0*log(one-z))

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    call cut_histo(HardProc)

    if (HardProc%makecut) then

       FintNNLO_s_qg_oewk_fs(1) = zero
       kin(1) = zero

    else

       call res_tree_a_qqb(HardProc%AmpMom,res_nlo)
       call get_respdf(qg_lumi,1,1,HardProc,res_nlo,respdf)

       EC    = HardProc%Lim_Ei(1)
       eta52 = HardProc%Lim_etaij(2,5)
       call partition_nnlo_fact(HardProc,damp,iconf_qed=[1,2,3,4,5],i_qed=2)

       call fill_sv_logs(HardProc%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqgNLO + Pqg0 * log(eta52/2) * damp
       intsub = Tr * (intsub - Pqg0*sv_logs)

       call partition_nlo_qed(HardProc,damp_ewk,icoll)
       respdf = intsub * respdf * HardProc%wgt * damp_ewk

       FintNNLO_s_qg_oewk_fs(1) = respdf(1)
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                            Collinear Limit                            !!
    !!-----------------------------------------------------------------------!!
    CLim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(CLim)

    if (CLim%makecut) then

       FintNNLO_s_qg_oewk_fs(2) = zero
       kin(2) = zero

    else

       call res_tree_qqb(CLim%AmpMom,res_lo)

       call get_respdf(qg_lumi,1,1,CLim,res_lo,respdf)

       EC = CLim%Lim_Ei(1)
       call fill_sv_logs(CLim%muf(1)**2,4*EC**2,sv_logs)

       !-- damp = zero
       intsub = PqgNLO
       intsub = Tr * (intsub - Pqg0*sv_logs)

       !-- damp_ewk = one
       s5i = CLim%Lim_sij(icoll,5)
       z5  = CLim%Lim_z(1)
       respdf = intsub * respdf * CLim%wgt &
              * 2/s5i * Pqg(z5) * Q_lep2
       respdf = - respdf

       FintNNLO_s_qg_oewk_fs(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Soft Limit                               !!
    !!-----------------------------------------------------------------------!!
    SLim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(SLim)

    if (SLim%makecut) then

       kin(3) = zero
       FintNNLO_s_qg_oewk_fs(3:4) = zero

    else

       call res_tree_qqb(SLim%AmpMom,res_lo)

       E5   = SLim%Lim_Ei(2)
       E5sq = E5*E5
       call get_qed_eik(charges,SLim%Lim_etaij,ql_pos,5,eik_qed)
       eik_qed = eik_qed/E5sq

       res_tmp_vect(:,1,1) = res_lo(:,1) * eik_qed(1:2)
       res_tmp_vect(:,2,1) = res_lo(:,2) * eik_qed(3:4)

       !-- prepare PDFs structures, CS
       res_tmp_vect(:,:,2) = res_lo(:,:) * Q_lep2

       !-- get PDFs for all of them
       call get_respdf_vect(qg_lumi,1,1,SLim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))

       !! ----------------------------- SLim ------------------------------- !!
       EC    = SLim%Lim_Ei(1)
       call fill_sv_logs(SLim%muf(1)**2,4*EC**2,sv_logs)
       eta52 = SLim%Lim_etaij(2,5)
       call partition_nnlo_fact(SLim,damp,iconf_qed=[1,2,3,4,5],i_qed=2)
       intsub = PqgNLO + Pqg0 * log(eta52/2) * damp
       intsub = Tr * (intsub - Pqg0*sv_logs)

       call partition_nlo_qed(SLim,damp_ewk,icoll)
       respdf_tmp(:,1) = intsub * respdf_vect(1,:) * SLim%wgt * damp_ewk
       respdf_tmp(:,1) = - respdf_tmp(:,1)

       FintNNLO_s_qg_oewk_fs(3) = respdf_tmp(1,1)

       !! ----------------------------- SCLim ------------------------------ !!
       !-- damp = zero
       intsub = PqgNLO
       intsub = Tr * (intsub - Pqg0*sv_logs)

       !-- damp_ewk = one
       eta5i = SCLim%Lim_etaij(icoll,5)
       respdf_tmp(:,2) = intsub * respdf_vect(2,:) * SCLim%wgt &
                       * one/eta5i/E5sq

       FintNNLO_s_qg_oewk_fs(4) = respdf_tmp(1,2)

       !! --------------------------- SLim tot ---------------------------- !!

       respdf(:) = sum(respdf_tmp(:,1:2),2)
       kin(3) = respdf(1)       

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_qg_onlo(1:4) = FintNNLO_s_qg_oewk_fs
#endif

  end function xsect_nnlo_s_qg_oewk_fs_5i

  function xsect_nnlo_s_qg_vewknf(yRnd,ff,vegasweight)
    use mod_kinematics_lo
    integer :: xsect_nnlo_s_qg_vewknf
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: respdf(ipdf)
    real(dp) :: res_nlo(2,2),kin(1)

    xsect_nnlo_s_qg_vewknf = 0

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

       call res_ewkloop_qqb_nf(LOProc%AmpMom,res_nlo)

       !-- F[1,z.2]
       call get_respdf_hoppet(PDFs,xPij,ns_lumi,1,1,LOProc,res_nlo,respdf,myPDFs2_Lmu=[xPij_Lmu])
       respdf = Tr * respdf * LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_qg(1) = kin(1)
#endif

  end function xsect_nnlo_s_qg_vewknf

end module mod_xsects_nnlo_s_qg

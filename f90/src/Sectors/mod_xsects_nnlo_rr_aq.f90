module mod_xsects_nnlo_rr_aq
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_aux_sectors
  use mod_process
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_splittings_bare
  use mod_eikonals
  use mod_partitions
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_rr_tc_aq(8), FintNNLO_rr_dc_aq(8)
#endif

  !! positions of coloured fermions
  integer, parameter :: q_pos(2) = [2,6]

  public :: xsect_nnlo_rr_5161a_aq,xsect_nnlo_rr_5161c_aq
  public :: xsect_nnlo_rr_5161b_aq,xsect_nnlo_rr_5161d_aq
  public :: xsect_nnlo_rr_5262a_aq,xsect_nnlo_rr_5262c_aq
  public :: xsect_nnlo_rr_5262b_aq,xsect_nnlo_rr_5262d_aq
  public :: xsect_nnlo_rr_5162_aq, xsect_nnlo_rr_5261_aq

  private

  !! -- a(p1) + qb/q(p2) -> l-(p3) l+(p4) g(p5) qb/q(p6)

contains

  function xsect_nnlo_rr_5161b_aq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161b_aq
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5161b_aq = xsect_nnlo_rr_5161bd_aq(&
      yRnd,ff,vegasweight,2)

  end function xsect_nnlo_rr_5161b_aq

  function xsect_nnlo_rr_5161d_aq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161d_aq
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5161d_aq = xsect_nnlo_rr_5161bd_aq(&
      yRnd,ff,vegasweight,4)

  end function xsect_nnlo_rr_5161d_aq

  function xsect_nnlo_rr_5262b_aq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5262b_aq
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5262b_aq = xsect_nnlo_rr_5262bd_aq(&
      yRnd,ff,vegasweight,2)

  end function xsect_nnlo_rr_5262b_aq

  function xsect_nnlo_rr_5262d_aq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5262d_aq
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5262d_aq = xsect_nnlo_rr_5262bd_aq(&
      yRnd,ff,vegasweight,4)

  end function xsect_nnlo_rr_5262d_aq

  !!*************************************************************************!!

  function xsect_nnlo_rr_5161a_aq(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5161a_aq
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: S5Lim, TCLim, TCC6Lim, TCS5Lim, TCC6S5Lim, C6Lim, C6S5Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_aq(8),kin(1:6)
    real(dp)  :: respdf(ipdf),respdf_tmp(ipdf,2),respdf_bak(ipdf)
    real(dp)  :: res_tmp(2,2), damp
    real(dp)  :: eik_qcd,z1,z5,z6,s15,s16,s56,s156,E5,E5sq
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5161a_aq = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_abcd_5i6iac(&
         xx(1:kNNLO_max_full), 1, 2, 1, &
         HardProc, S5Lim=S5Lim, TCLim=TCLim, TCC6Lim=TCC6Lim, &
         TCS5Lim=TCS5Lim, TCC6S5Lim=TCC6S5Lim, C6Lim=C6Lim, C6S5Lim=C6S5Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_g,id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_aq(1) = zero

    else

      call res_tree_ga_aq(HardProc%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,HardProc,res_tmp,respdf)

      call partition_nnlo_qcd(HardProc, damp, iconf, 1,1)

      respdf = respdf * HardProc%wgt * damp

      kin(1) = respdf(1)
      FintNNLO_aq(1) = kin(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 5                                !!
    !!-----------------------------------------------------------------------!!

    S5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    S5Lim%npart = 5

    call cut_histo(S5Lim)

    if (S5Lim%makecut) then

      kin(2) = zero
      FintNNLO_aq(2) = zero

    else

      E5   = S5Lim%Lim_KinInv(1)
      E5sq = E5*E5

      call get_qcd_eik(CF,S5Lim%Lim_etaij,q_pos,5,eik_qcd)
      eik_qcd = eik_qcd/E5sq

      call res_tree_a_aq(S5Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,S5Lim,res_tmp,respdf)

      call partition_nnlo_qcd(S5Lim, damp, iconf, 1,1)

      respdf =  respdf * eik_qcd * S5Lim%wgt * damp
      respdf = -respdf
      FintNNLO_aq(2) = respdf(1)

      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Triple Collinear 5                           !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(3) = zero
      FintNNLO_aq(3:4) = zero

    else

      call res_tree_qqb(TCLim%AmpMom,res_tmp)
      res_tmp(:,1) = res_tmp(:,1) * Qdn2
      res_tmp(:,2) = res_tmp(:,2) * Qup2

      call get_respdf(aq_lumi,1,1,TCLim,res_tmp,respdf_bak)

      !! --------------------------------- TC ------------------------------ !!

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s15 = TCLim%Lim_sij(1,5)
      s16 = TCLim%Lim_sij(1,6)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * xn * CF * (-Pgaq(-s15,s56,-s16,z5,z1,z6)) * TCLim%wgt
      respdf_tmp(:,1) = -respdf_tmp(:,1)

      FintNNLO_aq(3) = respdf_tmp(1,1)

      !! ------------------------------ TC + C6 ---------------------------- !!

      !! (- si5 - si6 + s56) in the x4 -> 0 limit
      s156 = TCC6Lim%Lim_KinInv(1)
      z5  = TCC6Lim%Lim_z(1)
      z6  = TCC6Lim%Lim_z(2)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (-2/s156) * xn * Pqg(z5)  &
                      * (  2/s16) * CF * Pqq(z6) &
                      * TCC6Lim%wgt

      FintNNLO_aq(4) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)
      kin(3) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                      Triple Collinear 5 + Soft 5                      !!
    !!-----------------------------------------------------------------------!!

    TCS5Lim%ids = [0,0,id_el,-id_el,0,0]
    TCS5Lim%npart = 4

    call cut_histo(TCS5Lim)

    if (TCS5Lim%makecut) then

      kin(4) = zero
      FintNNLO_aq(5:6) = zero

    else

      !! ------------------------------- TC + S5 --------------------------- !!

      call res_tree_qqb(TCS5Lim%AmpMom,res_tmp)
      res_tmp(:,1) = res_tmp(:,1) * Qdn2
      res_tmp(:,2) = res_tmp(:,2) * Qup2

      call get_respdf(aq_lumi,1,1,TCS5Lim,res_tmp,respdf_bak)

      z5 = TCS5Lim%Lim_z(1)
      z6 = TCS5Lim%Lim_z(2)
      z1 = TCS5Lim%Lim_z(3)
      s56 = TCS5Lim%Lim_sij(5,6)
      s16 = TCS5Lim%Lim_sij(1,6)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * (4*CF*z6/z5/s56)            &
                      * xn * (2/s16) * Pqg(z1)  &
                      * TCS5Lim%wgt
      FintNNLO_aq(5) = respdf_tmp(1,1)

      !! ----------------------------- TC + C6 + S5 ------------------------ !!

      !! s56 in the x4 -> 0 limit
      s56 = TCC6S5Lim%Lim_sij(5,6)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (4*CF*z6/z5/s56)            &
                      * xn * (2/s16) * Pqg(z1)  &
                      * TCC6S5Lim%wgt
      respdf_tmp(:,2) = -respdf_tmp(:,2)

      FintNNLO_aq(6) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)
      kin(4) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 6                              !!
    !!-----------------------------------------------------------------------!!

    C6Lim%ids = [0,0,id_el,-id_el,id_g,0]
    C6Lim%npart = 5

    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(5) = zero
      FintNNLO_aq(7) = zero

    else

      z6  = C6Lim%Lim_z(2)
      s16 = C6Lim%Lim_sij(1,6)

      call res_tree_g_qqb(C6Lim%AmpMom,res_tmp)
      res_tmp(:,1) = res_tmp(:,1) * Qdn2
      res_tmp(:,2) = res_tmp(:,2) * Qup2
      call get_respdf(aq_lumi,1,1,C6Lim,res_tmp,respdf)

      call partition_nnlo_qcd(C6Lim, damp, iconf, 1,1)

      respdf = (2/s16) * xn * Pqg(z6) * respdf * C6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_aq(7) = respdf(1)
      kin(5) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                         Collinear 6 + Soft 5                          !!
    !!-----------------------------------------------------------------------!!

    C6S5Lim%ids = [0,0,id_el,-id_el,0,0]
    C6S5Lim%npart = 4

    call cut_histo(C6S5Lim)

    if (C6S5Lim%makecut) then

      kin(6) = zero
      FintNNLO_aq(8) = zero

    else

      call res_tree_qqb(C6S5Lim%AmpMom,res_tmp)
      res_tmp(:,1) = res_tmp(:,1) * Qdn2
      res_tmp(:,2) = res_tmp(:,2) * Qup2
      call get_respdf(aq_lumi,1,1,C6S5Lim,res_tmp,respdf)
      
      E5   = C6S5Lim%Lim_KinInv(1)
      E5sq = E5*E5

      call get_qcd_eik(CF,C6S5Lim%Lim_etaij,q_pos,5,eik_qcd)
      eik_qcd = eik_qcd/E5sq

      z6  = C6S5Lim%Lim_z(2)
      s16 = C6S5Lim%Lim_sij(1,6)

      call partition_nnlo_qcd(C6S5Lim, damp, iconf, 1,1)

      respdf = respdf &
             * (2/s16) * xn * Pqg(z6) &
             * eik_qcd * C6S5Lim%wgt * damp

      FintNNLO_aq(8) = respdf(1)
      kin(6) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_aq)

#if(_withchecks == 1)
    FintNNLO_rr_tc_aq = FintNNLO_aq
#endif

  end function xsect_nnlo_rr_5161a_aq

  function xsect_nnlo_rr_5161c_aq(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5161c_aq
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: S5Lim, TCLim, TCS5Lim
    !--
    real(dp) :: xx(kNNLO_max_full)
    real(dp) :: FintNNLO_aq(4),kin(1:4)
    real(dp) :: respdf(ipdf), res_tmp(2,2)
    real(dp) :: damp
    real(dp) :: eik_qcd,z1,z5,z6,s15,s16,s56,E5,E5sq
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5161c_aq = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_abcd_5i6iac(&
         xx(1:kNNLO_max_full), 1, 2, 3, &
         HardProc, S5Lim=S5Lim, TCLim=TCLim, TCS5Lim=TCS5Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_g,id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_aq(1) = zero

    else

      call res_tree_ga_aq(HardProc%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,HardProc,res_tmp,respdf)

      call partition_nnlo_qcd(HardProc, damp, iconf, 1,1)

      respdf = respdf * HardProc%wgt * damp

      kin(1) = respdf(1)
      FintNNLO_aq(1) = kin(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 5                                !!
    !!-----------------------------------------------------------------------!!

    S5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    S5Lim%npart = 5

    call cut_histo(S5Lim)

    if (S5Lim%makecut) then

      kin(2) = zero
      FintNNLO_aq(2) = zero

    else

      E5   = S5Lim%Lim_KinInv(1)
      E5sq = E5*E5

      call get_qcd_eik(CF,S5Lim%Lim_etaij,q_pos,5,eik_qcd)
      eik_qcd = eik_qcd/E5sq

      call res_tree_a_aq(S5Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,S5Lim,res_tmp,respdf)

      call partition_nnlo_qcd(S5Lim, damp, iconf, 1,1)

      respdf =  respdf * eik_qcd * S5Lim%wgt * damp
      respdf = -respdf
      FintNNLO_aq(2) = respdf(1)

      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Triple Collinear 5                           !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(3) = zero
      FintNNLO_aq(3) = zero

    else

      call res_tree_qqb(TCLim%AmpMom,res_tmp)
      res_tmp(:,1) = res_tmp(:,1) * Qdn2
      res_tmp(:,2) = res_tmp(:,2) * Qup2

      call get_respdf(aq_lumi,1,1,TCLim,res_tmp,respdf)

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s15 = TCLim%Lim_sij(1,5)
      s16 = TCLim%Lim_sij(1,6)

      !-- damp = one
      respdf = respdf &
             * xn * CF * (-Pgaq(-s15,s56,-s16,z5,z1,z6)) * TCLim%wgt
      respdf = -respdf

      FintNNLO_aq(3) = respdf(1)
      kin(3) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                      Triple Collinear 5 + Soft 5                      !!
    !!-----------------------------------------------------------------------!!

    TCS5Lim%ids = [0,0,id_el,-id_el,0,0]
    TCS5Lim%npart = 4

    call cut_histo(TCS5Lim)

    if (TCS5Lim%makecut) then

      kin(4) = zero
      FintNNLO_aq(4) = zero

    else

      call res_tree_qqb(TCS5Lim%AmpMom,res_tmp)
      res_tmp(:,1) = res_tmp(:,1) * Qdn2
      res_tmp(:,2) = res_tmp(:,2) * Qup2

      call get_respdf(aq_lumi,1,1,TCS5Lim,res_tmp,respdf)

      z5 = TCS5Lim%Lim_z(1)
      z6 = TCS5Lim%Lim_z(2)
      z1 = TCS5Lim%Lim_z(3)
      s56 = TCS5Lim%Lim_sij(5,6)
      s16 = TCS5Lim%Lim_sij(1,6)

      !-- damp = one
      respdf = respdf &
             * (4*CF*z6/z5/s56)        &
             * xn * (2/s16) * Pqg(z1)  &
             * TCS5Lim%wgt
      FintNNLO_aq(4) = respdf(1)
      kin(4) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_aq)

#if(_withchecks == 1)
    FintNNLO_rr_tc_aq(1:4) = FintNNLO_aq
    FintNNLO_rr_tc_aq(5:8) = zero
#endif

  end function xsect_nnlo_rr_5161c_aq

  function xsect_nnlo_rr_5161bd_aq(yRnd,ff,vegasweight,&
    sec)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5161bd_aq, sec
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc, S5Lim, TCLim, TCS5Lim
    type(KinConfig) :: C6Lim, C6S5Lim, TCC6Lim, TCC6S5Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_aq(8),kin(1:6)
    real(dp)  :: respdf(ipdf),respdf_tmp(ipdf,2),respdf_bak(ipdf)
    real(dp)  :: res_tmp(2,2)
    real(dp)  :: eik_qcd,z1,z5,z6,s15,s16,s56,s156,E5,E5sq,damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5161bd_aq = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_abcd_5i6ibd(&
         xx(1:kNNLO_max_full), 1, 2, sec, &
         HardProc, S5Lim=S5Lim, TCLim=TCLim, TCC6Lim=TCC6Lim, &
         TCS5Lim=TCS5Lim, TCC6S5Lim=TCC6S5Lim, C6Lim=C6Lim, C6S5Lim=C6S5Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_g,id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_aq(1) = zero

    else

      call res_tree_ga_aq(HardProc%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,HardProc,res_tmp,respdf)

      call partition_nnlo_qcd(HardProc, damp, iconf, 1,1)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_aq(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 5                                !!
    !!-----------------------------------------------------------------------!!

    S5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    S5Lim%npart = 5

    call cut_histo(S5Lim)

    if (S5Lim%makecut) then

      kin(2) = zero
      FintNNLO_aq(2) = zero

    else

      E5   = S5Lim%Lim_KinInv(1)
      E5sq = E5*E5

      call get_qcd_eik(CF,S5Lim%Lim_etaij,q_pos,5,eik_qcd)
      eik_qcd = eik_qcd/E5sq

      call res_tree_a_aq(S5Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,S5Lim,res_tmp,respdf)

      call partition_nnlo_qcd(S5Lim, damp, iconf, 1,1)

      respdf =  respdf * eik_qcd * S5Lim%wgt * damp
      respdf = -respdf
      FintNNLO_aq(2) = respdf(1)

      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Triple Collinear                             !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(3) = zero
      FintNNLO_aq(3:4) = zero

    else

      call res_tree_qqb(TCLim%AmpMom,res_tmp)
      res_tmp(:,1) = res_tmp(:,1) * Qdn2
      res_tmp(:,2) = res_tmp(:,2) * Qup2

      call get_respdf(aq_lumi,1,1,TCLim,res_tmp,respdf_bak)

      !! --------------------------------- TC ------------------------------ !!

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s15 = TCLim%Lim_sij(1,5)
      s16 = TCLim%Lim_sij(1,6)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * xn * CF * (-Pgaq(-s15,s56,-s16,z5,z1,z6)) * TCLim%wgt
      respdf_tmp(:,1) = -respdf_tmp(:,1)

      FintNNLO_aq(3) = respdf_tmp(1,1)

      !! ------------------------------ TC + C6 ---------------------------- !!

      s56 = TCC6Lim%Lim_sij(5,6)
      !! (- si5 - si6 + s56) in the x4 -> 0 limit
      s156 = TCC6Lim%Lim_KinInv(1)
      z5   = TCC6Lim%Lim_z(2)
      z6   = TCC6Lim%Lim_z(3)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (-2/s156) * xn * Pqg(z6)  &
                      * ( 2/s56) * CF * Pqq(z5) &
                      * TCC6Lim%wgt

      FintNNLO_aq(4) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)
      kin(3) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif


    !!-----------------------------------------------------------------------!!
    !!                             Collinear 56                              !!
    !!-----------------------------------------------------------------------!!

    C6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C6Lim%npart = 5

    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(4) = zero
      FintNNLO_aq(5) = zero

    else

      z6  = C6Lim%Lim_z(2)
      s56 = C6Lim%Lim_sij(5,6)

      call res_tree_a_aq(C6Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,C6Lim,res_tmp,respdf)

      call partition_nnlo_qcd(C6Lim, damp, iconf, 1,1)

      respdf = (two/s56) * CF * Pqq(z6) * respdf * C6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_aq(5) = respdf(1)
      kin(4) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Soft 5 + C65                             !!
    !!-----------------------------------------------------------------------!!

    C6S5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C6S5Lim%npart = 5

    call cut_histo(C6S5Lim)

    if (C6S5Lim%makecut) then

      kin(5) = zero
      FintNNLO_aq(6) = zero

    else

      z6  = C6S5Lim%Lim_z(2)
      s56 = C6S5Lim%Lim_sij(5,6)

      call res_tree_a_aq(C6S5Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,C6S5Lim,res_tmp,respdf)

      call partition_nnlo_qcd(C6S5Lim, damp, iconf, 1,1)

      respdf = respdf * (4*CF/z6/s56) * C6S5Lim%wgt * damp

      FintNNLO_aq(6)  = respdf(1)
      kin(5) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                      Triple Collinear 5 + Soft 5                      !!
    !!-----------------------------------------------------------------------!!

    TCS5Lim%ids = [0,0,id_el,-id_el,0,0]
    TCS5Lim%npart = 4

    call cut_histo(TCS5Lim)

    if (TCS5Lim%makecut) then

      kin(6) = zero
      FintNNLO_aq(7:8) = zero

    else

      !! ------------------------------- TC + S5 --------------------------- !!

      call res_tree_qqb(TCS5Lim%AmpMom,res_tmp)
      res_tmp(:,1) = res_tmp(:,1) * Qdn2
      res_tmp(:,2) = res_tmp(:,2) * Qup2

      call get_respdf(aq_lumi,1,1,TCS5Lim,res_tmp,respdf_bak)

      z5 = TCS5Lim%Lim_z(1)
      z6 = TCS5Lim%Lim_z(2)
      z1 = TCS5Lim%Lim_z(3)
      s56 = TCS5Lim%Lim_sij(5,6)
      s16 = TCS5Lim%Lim_sij(1,6)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * (4*CF*z6/z5/s56)        &
                      * xn * (2/s16) * Pqg(z1)  &
                      * TCS5Lim%wgt
      FintNNLO_aq(7) = respdf_tmp(1,1)

      !! ----------------------------- TC + C6 + S5 ------------------------ !!

      !! s16 and s56 in the x4 -> 0 limit
      s56 = TCC6S5Lim%Lim_sij(5,6)
      s16 = TCC6S5Lim%Lim_sij(1,6)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (4*CF*z6/z5/s56)        &
                      * xn * (2/s16) * Pqg(z1)  &
                      * TCC6S5Lim%wgt
      respdf_tmp(:,2) = -respdf_tmp(:,2)
      FintNNLO_aq(8) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)
      kin(6) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_aq)

#if(_withchecks == 1)
    FintNNLO_rr_tc_aq = FintNNLO_aq
#endif

  end function xsect_nnlo_rr_5161bd_aq

  function xsect_nnlo_rr_5262a_aq(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5262a_aq
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: S5Lim, TCLim, TCS5Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_aq(4),kin(1:4)
    real(dp)  :: respdf(ipdf), res_tmp(2,2),res_aaee
    real(dp)  :: eik_qcd,z1,z5,z6,s25,s26,s56,E5,E5sq,damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5262a_aq = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_abcd_5i6iac(&
         xx(1:kNNLO_max_full), 2, 1, 1, &
         HardProc, S5Lim=S5Lim, TCLim=TCLim, TCS5Lim=TCS5Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_g,id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_aq(1) = zero

    else

      call res_tree_ga_aq(HardProc%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,HardProc,res_tmp,respdf)

      call partition_nnlo_qcd(HardProc, damp, iconf, 2,2)

      respdf = respdf * HardProc%wgt * damp

      kin(1) = respdf(1)
      FintNNLO_aq(1) = kin(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 5                                !!
    !!-----------------------------------------------------------------------!!

    S5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    S5Lim%npart = 5

    call cut_histo(S5Lim)

    if (S5Lim%makecut) then

      kin(2) = zero
      FintNNLO_aq(2) = zero

    else

      E5   = S5Lim%Lim_KinInv(1)
      E5sq = E5*E5

      call get_qcd_eik(CF,S5Lim%Lim_etaij,q_pos,5,eik_qcd)
      eik_qcd = eik_qcd/E5sq

      call res_tree_a_aq(S5Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,S5Lim,res_tmp,respdf)

      call partition_nnlo_qcd(S5Lim, damp, iconf, 2,2)

      respdf =  respdf * eik_qcd * S5Lim%wgt * damp
      respdf = -respdf

      FintNNLO_aq(2) = respdf(1)
      kin(2) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Triple Collinear 5                           !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(3) = zero
      FintNNLO_aq(3) = zero

    else

      call res_treeAA_aa(TCLim%AmpMom,res_aaee)
      res_tmp(:,1) = res_aaee * Qdn2
      res_tmp(:,2) = res_aaee * Qup2

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 =  TCLim%Lim_sij(5,6)
      s25 = -TCLim%Lim_sij(2,5)
      s26 = -TCLim%Lim_sij(2,6)

      !-- damp = one
      res_tmp = res_tmp &
              * Cf * (-Pgqqb_spav_ab(s25,s56,s26,z5,z1,z6)) &
              * TCLim%wgt

      call get_respdf(aq_lumi,1,1,TCLim,res_tmp,respdf)
      respdf = -respdf

      FintNNLO_aq(3) = respdf(1)
      kin(3) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                      Triple Collinear 5 + Soft 5                      !!
    !!-----------------------------------------------------------------------!!

    TCS5Lim%ids = [0,0,id_el,-id_el,0,0]
    TCS5Lim%npart = 4

    call cut_histo(TCS5Lim)

    if (TCS5Lim%makecut) then

      kin(4) = zero
      FintNNLO_aq(4) = zero

    else

      call res_treeAA_aa(TCS5Lim%AmpMom,res_aaee)
      res_tmp(:,1) = res_aaee * Qdn2
      res_tmp(:,2) = res_aaee * Qup2

      z1 = TCS5Lim%Lim_z(3)
      s25 = TCS5Lim%Lim_sij(2,5)
      s56 = TCS5Lim%Lim_sij(5,6)

      !-- damp = one
      res_tmp = res_tmp &
              * (4*CF/s56) * (2/s25) * Pgq_spav(z1) &
              * TCS5Lim%wgt

      call get_respdf(aq_lumi,1,1,TCS5Lim,res_tmp,respdf)

      FintNNLO_aq(4) = respdf(1)
      kin(4) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_aq)

#if(_withchecks == 1)
    FintNNLO_rr_tc_aq(1:4) = FintNNLO_aq
    FintNNLO_rr_tc_aq(5:8) = zero
#endif

  end function xsect_nnlo_rr_5262a_aq

  function xsect_nnlo_rr_5262c_aq(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5262c_aq
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: S5Lim, TCLim, TCS5Lim, C5Lim, C5S5Lim, TCC5Lim, TCC5S5Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_aq(8),kin(1:6)
    real(dp)  :: respdf(ipdf),respdf_tmp(ipdf,2),respdf_bak(ipdf)
    real(dp)  :: res_tmp(2,2),res_aaee
    real(dp)  :: eik_qcd,z1,z5,z6,s25,s26,s56,s256,E5,E5sq,damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5262c_aq = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_abcd_5i6iac(&
         xx(1:kNNLO_max_full), 2, 1, 3, &
         HardProc, S5Lim=S5Lim, TCLim=TCLim, TCC5Lim=TCC5Lim, &
         TCS5Lim=TCS5Lim, TCC5S5Lim=TCC5S5Lim, C5Lim=C5Lim, C5S5Lim=C5S5Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_g,id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_aq(1) = zero

    else

      call res_tree_ga_aq(HardProc%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,HardProc,res_tmp,respdf)

      call partition_nnlo_qcd(HardProc, damp, iconf, 2,2)

      respdf = respdf * HardProc%wgt * damp

      kin(1) = respdf(1)
      FintNNLO_aq(1) = kin(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 5                                !!
    !!-----------------------------------------------------------------------!!

    S5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    S5Lim%npart = 5

    call cut_histo(S5Lim)

    if (S5Lim%makecut) then

      kin(2) = zero
      FintNNLO_aq(2) = zero

    else

      E5   = S5Lim%Lim_KinInv(1)
      E5sq = E5*E5

      call get_qcd_eik(CF,S5Lim%Lim_etaij,q_pos,5,eik_qcd)
      eik_qcd = eik_qcd/E5sq

      call res_tree_a_aq(S5Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,S5Lim,res_tmp,respdf)

      call partition_nnlo_qcd(S5Lim, damp, iconf, 2,2)

      respdf =  respdf * eik_qcd * S5Lim%wgt * damp
      respdf = -respdf

      FintNNLO_aq(2) = respdf(1)
      kin(2) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Triple Collinear 5                           !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(3) = zero
      FintNNLO_aq(3:4) = zero

    else

    !! --------------------------------- TC -------------------------------- !!

      call res_treeAA_aa(TCLim%AmpMom,res_aaee)
      res_tmp(:,1) = res_aaee * Qdn2
      res_tmp(:,2) = res_aaee * Qup2

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 =  TCLim%Lim_sij(5,6)
      s25 = -TCLim%Lim_sij(2,5)
      s26 = -TCLim%Lim_sij(2,6)

      call get_respdf(aq_lumi,1,1,TCLim,res_tmp,respdf_bak)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * Cf * (-Pgqqb_spav_ab(s25,s56,s26,z5,z1,z6)) &
                      * TCLim%wgt
      respdf_tmp(:,1) = -respdf_tmp(:,1)

      FintNNLO_aq(3) = respdf_tmp(1,1)

    !! ------------------------------- TC + C5 ----------------------------- !!

      z5 =  TCC5Lim%Lim_z(1)
      z6 =  TCC5Lim%Lim_z(2)
      s256 = TCC5Lim%Lim_KinInv(1)

      respdf_tmp(:,2) = respdf_bak &
                      * 2/s256 * CF * Pqg(z5)  &
                      * 2/s25 * Pgq_spav(z6) &
                      * TCC5Lim%wgt

      FintNNLO_aq(4) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)

      kin(3) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                      Triple Collinear 5 + Soft 5                      !!
    !!-----------------------------------------------------------------------!!

    TCS5Lim%ids = [0,0,id_el,-id_el,0,0]
    TCS5Lim%npart = 4

    call cut_histo(TCS5Lim)

    if (TCS5Lim%makecut) then

      kin(4) = zero
      FintNNLO_aq(5:6) = zero

    else

    !! ------------------------------- TC + S5 ----------------------------- !!

      call res_treeAA_aa(TCS5Lim%AmpMom,res_aaee)
      res_tmp(:,1) = res_aaee * Qdn2
      res_tmp(:,2) = res_aaee * Qup2

      z1  = TCS5Lim%Lim_z(3)
      s25 = TCS5Lim%Lim_sij(2,5)
      s56 = TCS5Lim%Lim_sij(5,6)

      call get_respdf(aq_lumi,1,1,TCS5Lim,res_tmp,respdf_bak)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * (4*CF/s56) * (2/s25) * Pgq_spav(z1) &
                      * TCS5Lim%wgt

      FintNNLO_aq(5) = respdf_tmp(1,1)

    !! ---------------------------- TC + S5 + C5 --------------------------- !!

      z5  = TCC5S5Lim%Lim_z(1)
      z6  = TCC5S5Lim%Lim_z(2)
      s25 = TCC5S5Lim%Lim_sij(2,5)
      s26 = TCC5S5Lim%Lim_sij(2,6)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (4*CF/s25) * (2/z5/s26) * Pgq_spav(z6) &
                      * TCC5S5Lim%wgt
      respdf_tmp(:,2) = -respdf_tmp(:,2)

      FintNNLO_aq(6) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)

      kin(4) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 5                              !!
    !!-----------------------------------------------------------------------!!

    C5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C5Lim%npart = 5

    call cut_histo(C5Lim)

    if (C5Lim%makecut) then

      kin(5) = zero
      FintNNLO_aq(7) = zero

    else

      z5  = C5Lim%Lim_z(1)
      s25 = C5Lim%Lim_sij(2,5)

      call res_tree_a_aq(C5Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,C5Lim,res_tmp,respdf)

      call partition_nnlo_qcd(C5Lim, damp, iconf, 2,2)

      respdf =  respdf * Cf * (-2/s25 * Pqq(z5)) * C5Lim%wgt * damp
      respdf = -respdf
      FintNNLO_aq(7) = respdf(1)

      kin(5) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                         Soft 5 + Collinear 5                          !!
    !!-----------------------------------------------------------------------!!

    C5S5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C5S5Lim%npart = 5

    call cut_histo(C5S5Lim)

    if (C5S5Lim%makecut) then

      kin(6) = zero
      FintNNLO_aq(8) = zero

    else

      s25 = C5S5Lim%Lim_sij(2,5)
      z5  = C5S5Lim%Lim_z(1)

      call res_tree_a_aq(C5S5Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,C5S5Lim,res_tmp,respdf)

      call partition_nnlo_qcd(C5S5Lim, damp, iconf, 2,2)

      respdf = respdf * (4*CF/z5/s25) * C5S5Lim%wgt * damp

      FintNNLO_aq(8) = respdf(1)
      kin(6) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_aq)

#if(_withchecks == 1)
    FintNNLO_rr_tc_aq = FintNNLO_aq
#endif

  end function xsect_nnlo_rr_5262c_aq

  function xsect_nnlo_rr_5262bd_aq(yRnd,ff,vegasweight,&
    sec)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5262bd_aq, sec
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc, S5Lim, TCLim, TCS5Lim
    type(KinConfig) :: C6Lim, C6S5Lim, TCC6Lim, TCC6S5Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_aq(8),kin(1:6)
    real(dp)  :: respdf(ipdf),respdf_tmp(ipdf,2),respdf_bak(ipdf)
    real(dp)  :: res_tmp(2,2), res_aaee
    real(dp)  :: eik_qcd,z1,z5,z6,s25,s26,s56,s256,E5,E5sq,damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5262bd_aq = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_abcd_5i6ibd(&
         xx(1:kNNLO_max_full), 2, 1, sec, &
         HardProc, S5Lim=S5Lim, TCLim=TCLim, TCC6Lim=TCC6Lim, &
         TCS5Lim=TCS5Lim, TCC6S5Lim=TCC6S5Lim, C6Lim=C6Lim, C6S5Lim=C6S5Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_g,id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_aq(1) = zero

    else

      call res_tree_ga_aq(HardProc%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,HardProc,res_tmp,respdf)

      call partition_nnlo_qcd(HardProc, damp, iconf, 2,2)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_aq(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 5                                !!
    !!-----------------------------------------------------------------------!!

    S5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    S5Lim%npart = 5

    call cut_histo(S5Lim)

    if (S5Lim%makecut) then

      kin(2) = zero
      FintNNLO_aq(2) = zero

    else

      E5   = S5Lim%Lim_KinInv(1)
      E5sq = E5*E5

      call get_qcd_eik(CF,S5Lim%Lim_etaij,q_pos,5,eik_qcd)
      eik_qcd = eik_qcd/E5sq

      call res_tree_a_aq(S5Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,S5Lim,res_tmp,respdf)

      call partition_nnlo_qcd(S5Lim, damp, iconf, 2,2)

      respdf =  respdf * eik_qcd * S5Lim%wgt * damp
      respdf = -respdf
      FintNNLO_aq(2) = respdf(1)

      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Triple Collinear                             !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(3) = zero
      FintNNLO_aq(3:4) = zero

    else

      call res_treeAA_aa(TCLim%AmpMom,res_aaee)
      res_tmp(:,1) = res_aaee * Qdn2
      res_tmp(:,2) = res_aaee * Qup2

      call get_respdf(aq_lumi,1,1,TCLim,res_tmp,respdf_bak)

      !! --------------------------------- TC ------------------------------ !!

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s25 = TCLim%Lim_sij(2,5)
      s26 = TCLim%Lim_sij(2,6)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * CF * (-Pgqqb_spav_ab(-s25,s56,-s26,z5,z1,z6)) * TCLim%wgt
      respdf_tmp(:,1) = -respdf_tmp(:,1)

      FintNNLO_aq(3) = respdf_tmp(1,1)

      !! ------------------------------ TC + C6 ---------------------------- !!

      s56 = TCC6Lim%Lim_sij(5,6)
      !! (- si5 - si6 + s56) in the x4 -> 0 limit
      s256 = TCC6Lim%Lim_KinInv(1)
      z5   = TCC6Lim%Lim_z(2)
      z6   = TCC6Lim%Lim_z(3)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * 2/s56 * CF * Pqq(z5)  &
                      * (-2/s256) * Pgq_spav(z6) &
                      * TCC6Lim%wgt

      FintNNLO_aq(4) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)
      kin(3) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                             Collinear 56                              !!
    !!-----------------------------------------------------------------------!!

    C6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C6Lim%npart = 5

    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(4) = zero
      FintNNLO_aq(5) = zero

    else

      z6  = C6Lim%Lim_z(2)
      s56 = C6Lim%Lim_sij(5,6)

      call res_tree_a_aq(C6Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,C6Lim,res_tmp,respdf)

      call partition_nnlo_qcd(C6Lim, damp, iconf, 2,2)

      respdf = (two/s56) * CF * Pqq(z6) * respdf * C6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_aq(5) = respdf(1)
      kin(4) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Soft 5 + C65                             !!
    !!-----------------------------------------------------------------------!!

    C6S5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C6S5Lim%npart = 5

    call cut_histo(C6S5Lim)

    if (C6S5Lim%makecut) then

      kin(5) = zero
      FintNNLO_aq(6) = zero

    else

      z6  = C6S5Lim%Lim_z(2)
      s56 = C6S5Lim%Lim_sij(5,6)

      call res_tree_a_aq(C6S5Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,C6S5Lim,res_tmp,respdf)

      call partition_nnlo_qcd(C6S5Lim, damp, iconf, 2,2)

      respdf = respdf * (4*CF/z6/s56) * C6S5Lim%wgt * damp

      FintNNLO_aq(6)  = respdf(1)
      kin(5) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                      Triple Collinear 5 + Soft 5                      !!
    !!-----------------------------------------------------------------------!!

    TCS5Lim%ids = [0,0,id_el,-id_el,0,0]
    TCS5Lim%npart = 4

    call cut_histo(TCS5Lim)

    if (TCS5Lim%makecut) then

      kin(6) = zero
      FintNNLO_aq(7:8) = zero

    else

      !! ------------------------------- TC + S5 --------------------------- !!

      call res_treeAA_aa(TCS5Lim%AmpMom,res_aaee)
      res_tmp(:,1) = res_aaee * Qdn2
      res_tmp(:,2) = res_aaee * Qup2

      call get_respdf(aq_lumi,1,1,TCS5Lim,res_tmp,respdf_bak)

      z5 = TCS5Lim%Lim_z(1)
      z6 = TCS5Lim%Lim_z(2)
      z1 = TCS5Lim%Lim_z(3)
      s56 = TCS5Lim%Lim_sij(5,6)
      s25 = TCS5Lim%Lim_sij(2,5)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * (4*CF/s56)        &
                      * (2/s25) * Pgq_spav(z1)  &
                      * TCS5Lim%wgt
      FintNNLO_aq(7) = respdf_tmp(1,1)

      !! ----------------------------- TC + C6 + S5 ------------------------ !!

      !! s16 and s56 in the x4 -> 0 limit
      s56 = TCC6S5Lim%Lim_sij(5,6)
      s25 = TCC6S5Lim%Lim_sij(2,5)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (4*CF/s56)        &
                      * (2/s25) * Pgq_spav(z1)  &
                      * TCC6S5Lim%wgt
      respdf_tmp(:,2) = -respdf_tmp(:,2)
      FintNNLO_aq(8) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)
      kin(6) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_aq)

#if(_withchecks == 1)
    FintNNLO_rr_tc_aq = FintNNLO_aq
#endif

  end function xsect_nnlo_rr_5262bd_aq

  function xsect_nnlo_rr_5162_aq(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_dc_ii
    implicit none
    integer :: xsect_nnlo_rr_5162_aq
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc, S5Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_aq(2),kin(2)
    real(dp)  :: respdf(ipdf)
    real(dp)  :: res_tmp(2,2)
    real(dp)  :: eik_qcd,E5,E5sq,damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5162_aq = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_ii_5i6j(&
      xx(1:kNNLO_max_full), 1, 2, HardProc, S5Lim=S5Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_g,id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_aq(1) = zero

    else

      call res_tree_ga_aq(HardProc%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,HardProc,res_tmp,respdf)

      call partition_nnlo_qcd(HardProc, damp, iconf, 1,2)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_aq(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 5                                !!
    !!-----------------------------------------------------------------------!!

    S5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    S5Lim%npart = 5

    call cut_histo(S5Lim)

    if (S5Lim%makecut) then

      kin(2) = zero
      FintNNLO_aq(2) = zero

    else

      call res_tree_a_aq(S5Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,S5Lim,res_tmp,respdf)

      E5   = S5Lim%Lim_KinInv(1)
      E5sq = E5*E5
      call get_qcd_eik(CF,S5Lim%Lim_etaij,q_pos,5,eik_qcd)
      eik_qcd = eik_qcd/E5sq

      call partition_nnlo_qcd(S5Lim, damp, iconf, 1,2)

      respdf =  respdf * eik_qcd * S5Lim%wgt * damp
      respdf = -respdf

      FintNNLO_aq(2) = respdf(1)
      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_aq)

#if(_withchecks == 1)
    FintNNLO_rr_dc_aq(1:2) = FintNNLO_aq
    FintNNLO_rr_dc_aq(3:8) = zero
#endif

  end function xsect_nnlo_rr_5162_aq

  function xsect_nnlo_rr_5261_aq(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_dc_ii
    implicit none
    integer :: xsect_nnlo_rr_5261_aq
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: S5Lim, C5Lim, C5S5Lim, C6Lim, C5C6Lim, C6S5Lim, C5C6S5Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_aq(8),kin(1:6)
    real(dp)  :: respdf(ipdf),respdf_tmp(ipdf,2),respdf_bak(ipdf)
    real(dp)  :: res_tmp(2,2)
    real(dp)  :: eik_qcd,z5,z6,s16,s25,E5,E5sq,damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5261_aq = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_ii_5i6j(&
      xx(1:kNNLO_max_full), 2, 1, &
      HardProc, S5Lim=S5Lim, C5Lim=C5Lim, C5S5Lim=C5S5Lim,&
      C6Lim=C6Lim, C5C6Lim=C5C6Lim, C6S5Lim=C6S5Lim, C5C6S5Lim=C5C6S5Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_g,id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_aq(1) = zero

    else

      call res_tree_ga_aq(HardProc%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,HardProc,res_tmp,respdf)

      call partition_nnlo_qcd(HardProc, damp, iconf, 2,1)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_aq(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 5                                !!
    !!-----------------------------------------------------------------------!!

    S5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    S5Lim%npart = 5

    call cut_histo(S5Lim)

    if (S5Lim%makecut) then

      kin(2) = zero
      FintNNLO_aq(2:3) = zero

    else

      call res_tree_a_aq(S5Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,S5Lim,res_tmp,respdf_bak)

    !! --------------------------------- S5 --------------------------------- !!
      E5   = S5Lim%Lim_KinInv(1)
      E5sq = E5*E5
      call get_qcd_eik(CF,S5Lim%Lim_etaij,q_pos,5,eik_qcd)
      eik_qcd = eik_qcd/E5sq

      call partition_nnlo_qcd(S5Lim, damp, iconf, 2,1)

      respdf_tmp(:,1) =  respdf_bak * eik_qcd * S5Lim%wgt * damp
      respdf_tmp(:,1) = -respdf_tmp(:,1)

      FintNNLO_aq(2)  =  respdf_tmp(1,1)

    !! ------------------------------ C5 + S5 ------------------------------ !!

      z5  = C5S5Lim%Lim_z(1)
      s25 = C5S5Lim%Lim_sij(2,5)

      call partition_nnlo_qcd(C5S5Lim, damp, iconf, 2,1)

      respdf_tmp(:,2) = respdf_bak * (4*CF/z5/s25) * C5S5Lim%wgt * damp

      FintNNLO_aq(3)  = respdf_tmp(1,2)

      respdf = respdf_tmp(1,1) + respdf_tmp(1,2)
      kin(2) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 5                              !!
    !!-----------------------------------------------------------------------!!

    C5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C5Lim%npart = 5

    call cut_histo(C5Lim)

    if (C5Lim%makecut) then

      kin(3) = zero
      FintNNLO_aq(4) = zero

    else

      z5  = C5Lim%Lim_z(1)
      s25 = C5Lim%Lim_sij(2,5)

      call res_tree_a_aq(C5Lim%AmpMom,res_tmp)
      call get_respdf(aq_lumi,1,1,C5Lim,res_tmp,respdf)

      call partition_nnlo_qcd(C5Lim, damp, iconf, 2,1)

      respdf = (-two/s25) * CF * Pqq(z5) * respdf * C5Lim%wgt * damp
      respdf = -respdf

      FintNNLO_aq(4) = respdf(1)

      kin(3) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                         Collinear 6 + Soft 5                          !!
    !!-----------------------------------------------------------------------!!

    C6S5Lim%ids = [0,0,id_el,-id_el,0,0]
    C6S5Lim%npart = 4

    call cut_histo(C6S5Lim)

    if (C6S5Lim%makecut.or.C6S5Lim%flag) then

      kin(4) = zero
      FintNNLO_aq(5:6) = zero

    else

      call res_tree_qqb(C6S5Lim%AmpMom,res_tmp)
      res_tmp(:,1) = res_tmp(:,1) * Qdn2
      res_tmp(:,2) = res_tmp(:,2) * Qup2

      call get_respdf(aq_lumi,1,1,C6S5Lim,res_tmp,respdf_bak)

    !! ------------------------------ C6 + S5 ------------------------------ !!
      E5   = C6S5Lim%Lim_KinInv(1)
      E5sq = E5*E5
      call get_qcd_eik(CF,C6S5Lim%Lim_etaij,q_pos,5,eik_qcd)
      eik_qcd = eik_qcd/E5sq

      z6  = C6S5Lim%Lim_z(2)
      s16 = C6S5Lim%Lim_sij(1,6)

      call partition_nnlo_qcd(C6S5Lim, damp, iconf, 2,1)

      respdf_tmp(:,1) = respdf_bak &
                      * eik_qcd &
                      * (two/s16) * xn * Pqg(z6) &
                      * C6S5Lim%wgt * damp

      FintNNLO_aq(5) = respdf_tmp(1,1)

    !! --------------------------- C6 + C5 + S5 ---------------------------- !!

      z5  = C5C6S5Lim%Lim_z(1)
      s25 = C5C6S5Lim%Lim_sij(2,5)

      z6  = C5C6S5Lim%Lim_z(2)
      s16 = C5C6S5Lim%Lim_sij(1,6)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak    &
                      * (4*CF/z5/s25) &
                      * (two/s16) * xn * Pqg(z6) &
                      * C5C6S5Lim%wgt
      respdf_tmp(:,2) = - respdf_tmp(:,2)

      FintNNLO_aq(6) = respdf_tmp(1,2)

      respdf = respdf_tmp(1,1) + respdf_tmp(1,2)
      kin(4) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 6                              !!
    !!-----------------------------------------------------------------------!!

    C6Lim%ids = [0,0,id_el,-id_el,id_g,0]
    C6Lim%npart = 5

    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(5) = zero
      FintNNLO_aq(7) = zero

    else

      z6  = C6Lim%Lim_z(2)
      s16 = C6Lim%Lim_sij(1,6)

      call res_tree_g_qqb(C6Lim%AmpMom,res_tmp)
      res_tmp(:,1) = res_tmp(:,1) * Qdn2
      res_tmp(:,2) = res_tmp(:,2) * Qup2

      call get_respdf(aq_lumi,1,1,C6Lim,res_tmp,respdf)

      call partition_nnlo_qcd(C6Lim, damp, iconf, 2,1)

      respdf = respdf * (two/s16) * xn * Pqg(z6) * C6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_aq(7) = respdf(1)

      kin(5) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                       Collinear 5 + Collinear 6                       !!
    !!-----------------------------------------------------------------------!!

    C5C6Lim%ids = [0,0,id_el,-id_el,0,0]
    C5C6Lim%npart = 4

    call cut_histo(C5C6Lim)

    if (C5C6Lim%makecut) then

      kin(6) = zero
      FintNNLO_aq(8) = zero

    else

      z5  = C5C6Lim%Lim_z(1)
      s25 = C5C6Lim%Lim_sij(2,5)
      z6  = C5C6Lim%Lim_z(2)
      s16 = C5C6Lim%Lim_sij(1,6)

      call res_tree_qqb(C5C6Lim%AmpMom,res_tmp)
      res_tmp(:,1) = res_tmp(:,1) * Qdn2
      res_tmp(:,2) = res_tmp(:,2) * Qup2

      call get_respdf(aq_lumi,1,1,C5C6Lim,res_tmp,respdf)

      !-- damp = one
      respdf = respdf &
             * (-two/s25) * CF * Pqq(z5) &
             * ( two/s16) * xn * Pqg(z6) &
             * C5C6Lim%wgt

      FintNNLO_aq(8) = respdf(1)

      kin(6) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_aq)

#if(_withchecks == 1)
    FintNNLO_rr_dc_aq = FintNNLO_aq
#endif

  end function xsect_nnlo_rr_5261_aq

end module mod_xsects_nnlo_rr_aq

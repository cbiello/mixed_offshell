module mod_xsects_nnlo_rr_qg
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
  use mod_limvals
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_rr_tc_qg(8), FintNNLO_rr_dc_qg(8)
#endif

  !! Electric charges. 
  !! Each column is for the relevant res_jj(i,j)
  !! NB: different assignment wrt to matrix element
  real(dp), parameter :: charges(4,4) = reshape(& 
    [-Qdn,   Qdn,  -Qup,   Qup,   &
      Q_lep, Q_lep, Q_lep, Q_lep, &
     -Q_lep,-Q_lep,-Q_lep,-Q_lep, &
      Qdn,  -Qdn,   Qup,  -Qup  ],&
    [4,4])

  !! positions of charged quarks and leptons (ql)
  integer, parameter :: ql_pos(4) = [1,3,4,5]

  public :: xsect_nnlo_rr_5262a_qg,xsect_nnlo_rr_5262c_qg
  public :: xsect_nnlo_rr_5262b_qg,xsect_nnlo_rr_5262d_qg
  public :: xsect_nnlo_rr_5261_qg
  public :: xsect_nnlo_rr_5263_qg,xsect_nnlo_rr_5264_qg

  private

  !! -- q/qb(p1) + g(p2) -> l-(p3) l+(p4) q/qb(p5) a(p6)

contains

  function xsect_nnlo_rr_5262b_qg(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5262b_qg
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5262b_qg = xsect_nnlo_rr_5262bd_qg(&
      yRnd,ff,vegasweight,2)

  end function xsect_nnlo_rr_5262b_qg

  function xsect_nnlo_rr_5262d_qg(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5262d_qg
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5262d_qg = xsect_nnlo_rr_5262bd_qg(&
      yRnd,ff,vegasweight,4)

  end function xsect_nnlo_rr_5262d_qg

  function xsect_nnlo_rr_5263_qg(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5263_qg
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5263_qg = xsect_nnlo_rr_526k_qg(&
      yRnd,ff,vegasweight,3,4)

  end function xsect_nnlo_rr_5263_qg

  function xsect_nnlo_rr_5264_qg(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5264_qg
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5264_qg = xsect_nnlo_rr_526k_qg(&
      yRnd,ff,vegasweight,4,3)

  end function xsect_nnlo_rr_5264_qg

  !!*************************************************************************!!

  function xsect_nnlo_rr_5262a_qg(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5262a_qg
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: S6Lim, TCLim, TCS6Lim
    !--
    real(dp) :: xx(kNNLO_max_full)
    real(dp) :: FintNNLO_qg(4),kin(4)
    real(dp) :: respdf(ipdf)
    real(dp) :: res_tmp(2,2)
    real(dp) :: res_nnlo(-5:7,-5:7), res_nlo(-5:7,-5:7), res_lo(-5:7,-5:7), res_nlo_eikqed(-5:7,-5:7), res_lo_ischarges(-5:7,-5:7)
    real(dp) :: eik_qed(4),z1,z5,z6,s25,s26,s56,E6,E6sq
    real(dp) :: damp
    logical  :: oldcode

    oldcode = .false.


    xsect_nnlo_rr_5262a_qg = 0

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
         HardProc, S6Lim=S6Lim, TCLim=TCLim, TCS6Lim=TCS6Lim, opt_etas=[3,4])
#if (_Vcharge == 0)
    HardProc%ids(1:6) = [id_q,id_g,id_el,-id_el,id_q,id_a]
    S6Lim%ids(1:5)    = [id_q,id_g,id_el,-id_el,id_q]
    TCLim%ids(1:4)    = [id_q,-id_qp,id_el,-id_el]
    TCS6Lim%ids(1:4)  = [id_q,-id_qp,id_el,-id_el]
#elif  (_Vcharge == -1)  
    HardProc%ids(1:6) = [id_q,id_g,id_el,-id_nue,id_q,id_a]
    S6Lim%ids(1:5)    = [id_q,id_g,id_el,-id_nue,id_q]
    TCLim%ids(1:4)    = [id_q,-id_qp,id_el,-id_nue]
    TCS6Lim%ids(1:4)  = [id_q,-id_qp,id_el,-id_nue]
#elif  (_Vcharge == +1)  
    HardProc%ids(1:6) = [id_q,id_g,id_nue,-id_el,id_q,id_a]
    S6Lim%ids(1:5)    = [id_q,id_g,id_nue,-id_el,id_q]
    TCLim%ids(1:4)    = [id_q,-id_qp,id_nue,-id_el]
    TCS6Lim%ids(1:4)  = [id_q,-id_qp,id_nue,-id_el]
#endif
    HardProc%npart = 6
    S6Lim%npart    = 5
    TCLim%npart    = 4
    TCS6Lim%npart  = 4

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_qg(1) = zero

    else

       
       if (oldcode) then
          call res_tree_ga_qg(HardProc%AmpMom,res_tmp)
          call get_respdf(qg_lumi,1,1,HardProc,res_tmp,respdf)          
       else
          call res_tree_ga_qg_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,HardProc,res_nnlo,respdf)
       endif

       call partition_nnlo_fact(HardProc,damp,iconf_qed=[1,5,3,4,6],i_qed=2)
          

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_qg(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 6                                !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(S6Lim)

    if (S6Lim%makecut) then

      kin(2) = zero
      FintNNLO_qg(2) = zero

    else
       E6   = S6Lim%Lim_KinInv(1)
       E6sq = E6*E6

       if (oldcode) then
          call res_tree_g_qg(S6Lim%AmpMom,res_tmp)
          call get_qed_eik(charges,S6Lim%Lim_etaij,ql_pos,6,eik_qed)
          eik_qed = eik_qed/E6sq
          res_tmp(1:2,1) = res_tmp(1:2,1) * eik_qed(1:2)
          res_tmp(1:2,2) = res_tmp(1:2,2) * eik_qed(3:4)
          call get_respdf(qg_lumi,1,1,S6Lim,res_tmp,respdf)
       else
          call res_tree_g_qg_gen(S6Lim%AmpMom,res_nlo)
          call get_qed_eik_gen(res_nlo,S6Lim%Lim_etaij,[1,3,4,5],6,res_nlo_eikqed)
          call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf)
          respdf = respdf / E6sq
       endif

      call partition_nnlo_fact(S6Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=2)

      respdf =  respdf * S6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_qg(2) = respdf(1)
      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Triple Collinear 5                           !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(3) = zero
      FintNNLO_qg(3) = zero

    else

       if (oldcode) then
          call res_tree_qqb(TCLim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2
          
          call get_respdf(qg_lumi,1,1,TCLim,res_tmp,respdf)
       else
          call res_tree_qqb_gen(TCLim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,2)
          res_lo_ischarges = transition('none', 'g -> q', res_lo_ischarges)
          call get_respdf_gen(1,1,TCLim,res_lo_ischarges,respdf)
       endif

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s25 = TCLim%Lim_sij(2,5)
      s26 = TCLim%Lim_sij(2,6)

      !-- damp = one
      respdf = respdf &
             * Tr * (-Pgaq(-s26,s56,-s25,z6,z1,z5)) * TCLim%wgt
      respdf = -respdf

      FintNNLO_qg(3) = respdf(1)
      kin(3) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                      Triple Collinear 5 + Soft 6                      !!
    !!-----------------------------------------------------------------------!!


    call cut_histo(TCS6Lim)

    if (TCS6Lim%makecut) then

      kin(4) = zero
      FintNNLO_qg(4) = zero

    else

       if (oldcode) then
          call res_tree_qqb(TCS6Lim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2

          call get_respdf(qg_lumi,1,1,TCS6Lim,res_tmp,respdf)
       else
          call res_tree_qqb_gen(TCS6Lim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,2)
          res_lo_ischarges = transition('none', 'g -> q', res_lo_ischarges)
          call get_respdf_gen(1,1,TCS6Lim,res_lo_ischarges,respdf)
       endif

      z5  = TCS6Lim%Lim_z(1)
      z1  = TCS6Lim%Lim_z(3)
      s56 = TCS6Lim%Lim_sij(5,6)
      s25 = TCS6Lim%Lim_sij(2,5)

      !-- damp = one
      respdf = respdf * (4*Tr/z5/s56)  &
             * (2/s25) * Pqg(z1)  &
             * TCS6Lim%wgt

      FintNNLO_qg(4) = respdf(1)
      kin(4) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_qg)

#if(_withchecks == 1)
    FintNNLO_rr_tc_qg(1:4) = FintNNLO_qg
    FintNNLO_rr_tc_qg(5:8) = zero
#endif

  end function xsect_nnlo_rr_5262a_qg

  function xsect_nnlo_rr_5262c_qg(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5262c_qg
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: S6Lim, TCLim, TCS6Lim, C5Lim, C5S6Lim, TCC5Lim, TCC5S6Lim
    !--
    real(dp) :: xx(kNNLO_max_full)
    real(dp) :: FintNNLO_qg(8),kin(6)
    real(dp) :: respdf(ipdf),respdf_tmp(ipdf,2),respdf_bak(ipdf)
    real(dp) :: res_tmp(2,2)
    real(dp) :: res_nnlo(-5:7,-5:7), res_nlo(-5:7,-5:7), res_lo(-5:7,-5:7), res_nlo_eikqed(-5:7,-5:7),res_lo_eikqed(-5:7,-5:7), res_lo_ischarges(-5:7,-5:7)
    real(dp) :: eik_qed(4),z1,z5,z6,s25,s26,s56,s256,E6,E6sq
    real(dp) :: damp
    logical  :: oldcode

    oldcode = .false.

    xsect_nnlo_rr_5262c_qg = 0

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
         HardProc, S6Lim=S6Lim, TCLim=TCLim, TCC5Lim=TCC5Lim, &
         TCS6Lim=TCS6Lim, TCC5S6Lim=TCC5S6Lim, C5Lim=C5Lim, C5S6Lim=C5S6Lim, &
         opt_etas=[3,4])
    
#if (_Vcharge == 0)
    HardProc%ids(1:6) = [id_q,id_g,id_el,-id_el,id_q,id_a]
    S6Lim%ids(1:5)    = [id_q,id_g,id_el,-id_el,id_q]
    TCLim%ids(1:4)    = [id_q,-id_qp,id_el,-id_el]
    C5Lim%ids(1:5)    = [id_q,-id_qp,id_el,-id_el,id_a]
    C5S6Lim%ids(1:4)  = [id_q,-id_qp,id_el,-id_el]
    TCS6Lim%ids(1:4)  = [id_q,-id_qp,id_el,-id_el]    
#elif  (_Vcharge == -1)
    HardProc%ids(1:6) = [id_q,id_g,  id_el,-id_nue,id_q,id_a]
    S6Lim%ids(1:5)    = [id_q,id_g,  id_el,-id_nue,id_q]
    TCLim%ids(1:4)    = [id_q,-id_qp,id_el,-id_nue]
    C5Lim%ids(1:5)    = [id_q,-id_qp,id_el,-id_nue,id_a]
    C5S6Lim%ids(1:4)  = [id_q,-id_qp,id_el,-id_nue]
    TCS6Lim%ids(1:4)  = [id_q,-id_qp,id_el,-id_nue]    
#elif  (_Vcharge == +1)
    HardProc%ids(1:6) = [id_q,id_g,  id_nue,-id_el,id_q,id_a]
    S6Lim%ids(1:5)    = [id_q,id_g,  id_nue,-id_el,id_q]
    TCLim%ids(1:4)    = [id_q,-id_qp,id_nue,-id_el]
    C5Lim%ids(1:5)    = [id_q,-id_qp,id_nue,-id_el,id_a]
    C5S6Lim%ids(1:4)  = [id_q,-id_qp,id_nue,-id_el]
    TCS6Lim%ids(1:4)  = [id_q,-id_qp,id_nue,-id_el]   
#endif
    HardProc%npart = 6
    S6Lim%npart    = 5
    TCLim%npart    = 4
    C5Lim%npart    = 5
    C5S6Lim%npart  = 4
    TCS6Lim%npart  = 4





    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_qg(1) = zero

    else

       if (oldcode) then
          call res_tree_ga_qg(HardProc%AmpMom,res_tmp)
          call get_respdf(qg_lumi,1,1,HardProc,res_tmp,respdf)
       else
          call res_tree_ga_qg_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,HardProc,res_nnlo,respdf)
       endif
       
       call partition_nnlo_fact(HardProc,damp,iconf_qed=[1,5,3,4,6],i_qed=2)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_qg(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 6                                !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(S6Lim)

    if (S6Lim%makecut) then
       
       kin(2) = zero
       FintNNLO_qg(2) = zero
       
    else
       E6   = S6Lim%Lim_KinInv(1)
       E6sq = E6*E6
       if (oldcode) then
          call res_tree_g_qg(S6Lim%AmpMom,res_tmp)         
          call get_qed_eik(charges,S6Lim%Lim_etaij,ql_pos,6,eik_qed)
          eik_qed = eik_qed/E6sq
          res_tmp(1:2,1) = res_tmp(1:2,1) * eik_qed(1:2)
          res_tmp(1:2,2) = res_tmp(1:2,2) * eik_qed(3:4)
          call get_respdf(qg_lumi,1,1,S6Lim,res_tmp,respdf)
       else
          call res_tree_g_qg_gen(S6Lim%AmpMom,res_nlo)         
          call get_qed_eik_gen(res_nlo,S6Lim%Lim_etaij,[1,3,4,5],6,res_nlo_eikqed)
          call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf)
          respdf = respdf / E6sq
       endif
          

      call partition_nnlo_fact(S6Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=2)

      respdf =  respdf * S6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_qg(2) = respdf(1)
      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Triple Collinear 5                           !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(3) = zero
      FintNNLO_qg(3:4) = zero

    else

       if (oldcode) then
          call res_tree_qqb(TCLim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2
          
          call get_respdf(qg_lumi,1,1,TCLim,res_tmp,respdf_bak)
       else
          call res_tree_qqb_gen(TCLim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,2)
          res_lo_ischarges = transition('none', 'g -> q', res_lo_ischarges)
          call get_respdf_gen(1,1,TCLim,res_lo_ischarges,respdf_bak)
       endif
          

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s25 = TCLim%Lim_sij(2,5)
      s26 = TCLim%Lim_sij(2,6)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                 * Tr * (-Pgaq(-s26,s56,-s25,z6,z1,z5)) * TCLim%wgt
      respdf_tmp(:,1) = -respdf_tmp(:,1)

      FintNNLO_qg(3) = respdf_tmp(1,1)

      !! ----------------------------- TC + C5 ----------------------------- !!

      !! (- si5 - si6 + s56) in the x4 -> 0 limit
      s256 = TCC5Lim%Lim_KinInv(1)
      z5  = TCC5Lim%Lim_z(1)
      z6  = TCC5Lim%Lim_z(2)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (2/s256) * Tr * Pqg(z6)  &
                      * (-2/s25) * Pqq(z5) &
                      * TCC5Lim%wgt

      FintNNLO_qg(4) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)
      kin(3) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                      Triple Collinear 5 + Soft 6                      !!
    !!-----------------------------------------------------------------------!!


    call cut_histo(TCS6Lim)

    if (TCS6Lim%makecut) then

      kin(4) = zero
      FintNNLO_qg(5:6) = zero

    else

       if (oldcode) then
          call res_tree_qqb(TCS6Lim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2

          call get_respdf(qg_lumi,1,1,TCS6Lim,res_tmp,respdf_bak)
       else
          call res_tree_qqb_gen(TCS6Lim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,2)
          res_lo_ischarges = transition('none', 'g -> q', res_lo_ischarges)
          call get_respdf_gen(1,1,TCS6Lim,res_lo_ischarges,respdf_bak)
       endif
          

      z5  = TCS6Lim%Lim_z(1)
      z1  = TCS6Lim%Lim_z(3)
      s56 = TCS6Lim%Lim_sij(5,6)
      s25 = TCS6Lim%Lim_sij(2,5)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak * (4*Tr/z5/s56)   &
             * (2/s25) * Pqg(z1)  &
             * TCS6Lim%wgt

      FintNNLO_qg(5) = respdf_tmp(1,1)

      !! --------------------------- TC + C5 + S6 -------------------------- !!

      !! s56 in the x4 -> 0 limit
      s56 = TCC5S6Lim%Lim_sij(5,6)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (4*Tr/z5/s56)          &
                      * (2/s25) * Pqg(z1)  &
                      * TCC5S6Lim%wgt
      respdf_tmp(:,2) = -respdf_tmp(:,2)

      FintNNLO_qg(6) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) +  respdf_tmp(:,2)
      kin(4) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 5                              !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(C5Lim)

    if (C5Lim%makecut) then

      kin(5) = zero
      FintNNLO_qg(7) = zero

    else

      z5  = C5Lim%Lim_z(1)
      s25 = C5Lim%Lim_sij(2,5)

      if (oldcode) then
         call res_tree_a_qqb(C5Lim%AmpMom,res_tmp)
         call get_respdf(qg_lumi,1,1,C5Lim,res_tmp,respdf)
      else
         call res_tree_a_qqb_gen(C5Lim%AmpMom,res_nlo)
         res_nlo = transition('none', 'g -> q',res_nlo)
         call get_respdf_gen(1,1,C5Lim,res_nlo,respdf)
      endif

      call partition_nnlo_fact(C5Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=2)

      respdf = (two/s25) * Tr * Pqg(z5) * respdf * C5Lim%wgt * damp
      respdf = -respdf

      FintNNLO_qg(7) = respdf(1)
      kin(5) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                         Collinear 5 + Soft 6                          !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(C5S6Lim)

    if (C5S6Lim%makecut) then

      kin(6) = zero
      FintNNLO_qg(8) = zero

    else

       E6   = C5S6Lim%Lim_KinInv(1)
       E6sq = E6*E6

       if (oldcode) then
          call res_tree_qqb(C5S6Lim%AmpMom,res_tmp)
          call get_qed_eik(charges,C5S6Lim%Lim_etaij,ql_pos,6,eik_qed)
          eik_qed = eik_qed/E6sq          
          res_tmp(1:2,1) = res_tmp(1:2,1) * eik_qed(1:2)
          res_tmp(1:2,2) = res_tmp(1:2,2) * eik_qed(3:4)         
          call get_respdf(qg_lumi,1,1,C5S6Lim,res_tmp,respdf)
       else
          call res_tree_qqb_gen(C5S6Lim%AmpMom,res_lo)
          res_lo = transition('none','g -> q', res_lo)
          call get_qed_eik_gen(res_lo,C5S6Lim%Lim_etaij,[1,3,4,5],6,res_lo_eikqed)
          call get_respdf_gen(1,1,C5S6Lim,res_lo_eikqed,respdf)
          respdf = respdf / E6sq
       endif

      z5  = C5S6Lim%Lim_z(1)
      s25 = C5S6Lim%Lim_sij(2,5)

      call partition_nnlo_fact(C5S6Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=2)

      respdf = respdf * (two/s25) * Tr * Pqg(z5) * C5S6Lim%wgt * damp

      FintNNLO_qg(8) = respdf(1)
      kin(6) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
      
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_qg)

#if(_withchecks == 1)
    FintNNLO_rr_tc_qg = FintNNLO_qg
#endif

  end function xsect_nnlo_rr_5262c_qg

  function xsect_nnlo_rr_5262bd_qg(yRnd,ff,vegasweight,&
    sec)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5262bd_qg, sec
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: S6Lim, TCLim, TCS6Lim, C6Lim, C6S6Lim, TCC6Lim, TCC6S6Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_qg(8),kin(6)
    real(dp)  :: respdf(ipdf),respdf_tmp(ipdf,2),respdf_bak(ipdf)
    real(dp) :: res_nnlo(-5:7,-5:7), res_nlo(-5:7,-5:7), res_lo(-5:7,-5:7), res_nlo_eikqed(-5:7,-5:7),res_lo_eikqed(-5:7,-5:7), res_lo_ischarges(-5:7,-5:7)
    real(dp)  :: res_tmp(2,2)
    real(dp)  :: eik_qed(4),z1,z5,z6,s25,s26,s56,s256,E6,E6sq,damp,Qsq_FS(2)
    logical   :: oldcode

    oldcode = .false.

    xsect_nnlo_rr_5262bd_qg = 0

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
         HardProc, S6Lim=S6Lim, TCLim=TCLim, TCC6Lim=TCC6Lim, &
         TCS6Lim=TCS6Lim, TCC6S6Lim=TCC6S6Lim, C6Lim=C6Lim, C6S6Lim=C6S6Lim, &
         opt_etas=[3,4])

#if (_Vcharge == 0)
    HardProc%ids(1:6) = [id_q,id_g,  id_el,-id_el,id_q,id_a]
    S6Lim%ids(1:5)    = [id_q,id_g,  id_el,-id_el,id_q]
    TCLim%ids(1:4)    = [id_q,-id_qp,id_el,-id_el]
    TCS6Lim%ids(1:4)  = [id_q,-id_qp,id_el,-id_el]
    C6Lim%ids(1:5)    = [id_q,id_g,  id_el,-id_el,id_q]
    C6S6Lim%ids(1:5)  = [id_q,id_g,  id_el,-id_el,id_q]
#elif  (_Vcharge == -1)
    HardProc%ids(1:6) = [id_q,id_g,  id_el,-id_nue,id_q,id_a]
    S6Lim%ids(1:5)    = [id_q,id_g,  id_el,-id_nue,id_q]
    TCLim%ids(1:4)    = [id_q,-id_qp,id_el,-id_nue]
    TCS6Lim%ids(1:4)  = [id_q,-id_qp,id_el,-id_nue]
    C6Lim%ids(1:5)    = [id_q,id_g,  id_el,-id_nue,id_q]
    C6S6Lim%ids(1:5)  = [id_q,id_g,  id_el,-id_nue,id_q]
#elif  (_Vcharge == +1)
    HardProc%ids(1:6) = [id_q,id_g,  id_nue,-id_el,id_q,id_a]
    S6Lim%ids(1:5)    = [id_q,id_g,  id_nue,-id_el,id_q]
    TCLim%ids(1:4)    = [id_q,-id_qp,id_nue,-id_el]
    TCS6Lim%ids(1:4)  = [id_q,-id_qp,id_nue,-id_el]
    C6Lim%ids(1:5)    = [id_q,id_g,  id_nue,-id_el,id_q]
    C6S6Lim%ids(1:5)  = [id_q,id_g,  id_nue,-id_el,id_q]
#endif
    HardProc%npart = 6
    S6Lim%npart    = 5
    TCLim%npart    = 4
    TCS6Lim%npart  = 4
    C6Lim%npart    = 5
    C6S6Lim%npart  = 5

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_qg(1) = zero

    else

       if (oldcode) then
          call res_tree_ga_qg(HardProc%AmpMom,res_tmp)
          call get_respdf(qg_lumi,1,1,HardProc,res_tmp,respdf)
       else
          call res_tree_ga_qg_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,HardProc,res_nnlo,respdf)
       endif

      call partition_nnlo_fact(HardProc,damp,iconf_qed=[1,5,3,4,6],i_qed=2)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_qg(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 6                                !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(S6Lim)

    if (S6Lim%makecut) then

      kin(2) = zero
      FintNNLO_qg(2) = zero

    else

      E6   = S6Lim%Lim_KinInv(1)
      E6sq = E6*E6
      
      if (oldcode) then
         call res_tree_g_qg(S6Lim%AmpMom,res_tmp)
         call get_qed_eik(charges,S6Lim%Lim_etaij,ql_pos,6,eik_qed)
         eik_qed = eik_qed/E6sq
         res_tmp(1:2,1) = res_tmp(1:2,1) * eik_qed(1:2)
         res_tmp(1:2,2) = res_tmp(1:2,2) * eik_qed(3:4)
         
         call get_respdf(qg_lumi,1,1,S6Lim,res_tmp,respdf)
      else
         call res_tree_g_qg_gen(S6Lim%AmpMom,res_nlo)         
         call get_qed_eik_gen(res_nlo,S6Lim%Lim_etaij,[1,3,4,5],6,res_nlo_eikqed)
         call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf)
         respdf = respdf / E6sq
      endif

      call partition_nnlo_fact(S6Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=2)

      respdf =  respdf * S6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_qg(2) = respdf(1)
      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Triple Collinear                             !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(3) = zero
      FintNNLO_qg(3:4) = zero

    else

       if (oldcode) then
          call res_tree_qqb(TCLim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2
          call get_respdf(qg_lumi,1,1,TCLim,res_tmp,respdf_bak)
       else
          call res_tree_qqb_gen(TCLim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,2)
          res_lo_ischarges = transition('none', 'g -> q', res_lo_ischarges)
          call get_respdf_gen(1,1,TCLim,res_lo_ischarges,respdf_bak)
       endif

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s25 = TCLim%Lim_sij(2,5)
      s26 = TCLim%Lim_sij(2,6)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * Tr * (-Pgaq(-s26,s56,-s25,z6,z1,z5)) * TCLim%wgt
      respdf_tmp(:,1) = -respdf_tmp(:,1)

      FintNNLO_qg(3) = respdf_tmp(1,1)

      !! ----------------------------- TC + C6 ----------------------------- !!

      !! (- si5 - si6 + s56) in the x4 -> 0 limit
      s256 = -TCC6Lim%Lim_KinInv(1)
      s56  =  TCC6Lim%Lim_sij(5,6)
      z5   =  TCC6Lim%Lim_z(1)
      z6   =  TCC6Lim%Lim_z(3)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (2/s56) * Tr * Pqq(z5)  &
                      * (2/s256) * Pqg(z6) &
                      * TCC6Lim%wgt

      FintNNLO_qg(4) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) +  respdf_tmp(:,2)
      kin(3) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 6                              !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(4) = zero
      FintNNLO_qg(5) = zero

    else

      z5  = C6Lim%Lim_z(1)
      s56 = C6Lim%Lim_sij(5,6)
      if (oldcode) then
         call res_tree_g_qg(C6Lim%AmpMom,res_tmp)
         res_tmp(:,1) = res_tmp(:,1) * Qdn2
         res_tmp(:,2) = res_tmp(:,2) * Qup2
         call get_respdf(qg_lumi,1,1,C6Lim,res_tmp,respdf)
      else
         call res_tree_g_qg_gen(C6Lim%AmpMom,res_nlo)
         res_nlo = multiply_FSQ_charges_sq(res_nlo)
         call get_respdf_gen(1,1,C6Lim,res_nlo,respdf)
      endif

      !-- damp = one
      respdf = (two/s56) * Pqq(z5) * respdf * C6Lim%wgt
      respdf = -respdf

      FintNNLO_qg(5) = respdf(1)
      kin(4) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                         Collinear 6 + Soft 6                          !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(C6S6Lim)

    if (C6S6Lim%makecut) then

      kin(5) = zero
      FintNNLO_qg(6) = zero

    else

       if (oldcode) then
          call res_tree_g_qg(C6S6Lim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2                    
          call get_respdf(qg_lumi,1,1,C6S6Lim,res_tmp,respdf)
       else
          call res_tree_g_qg_gen(C6S6Lim%AmpMom,res_nlo)
          res_nlo = multiply_FSQ_charges_sq(res_nlo)
          call get_respdf_gen(1,1,C6S6Lim,res_nlo,respdf)
      endif

      
      z5  = C6S6Lim%Lim_z(1)
      s56 = C6S6Lim%Lim_sij(5,6)

      !-- damp = one
      respdf = respdf * (4/z5/s56) * C6S6Lim%wgt

      FintNNLO_qg(6) = respdf(1)
      kin(5) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                       Triple Collinear + Soft 6                       !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(TCS6Lim)

    if (TCS6Lim%makecut) then

      kin(6) = zero
      FintNNLO_qg(7:8) = zero

    else

       if (oldcode) then
          call res_tree_qqb(TCS6Lim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2
          call get_respdf(qg_lumi,1,1,TCS6Lim,res_tmp,respdf_bak)
       else
          call res_tree_qqb_gen(TCS6Lim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,2)
          res_lo_ischarges = transition('none', 'g -> q', res_lo_ischarges)
          call get_respdf_gen(1,1,TCS6Lim,res_lo_ischarges,respdf_bak)
       endif
       
      z5  = TCS6Lim%Lim_z(1)
      z6  = TCS6Lim%Lim_z(2)
      z1  = TCS6Lim%Lim_z(3)
      s56 = TCS6Lim%Lim_sij(5,6)
      s25 = TCS6Lim%Lim_sij(2,5)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * (4*Tr*z5/z6/s56)  &
                      * (2/s25) * Pqg(z1) &
                      * TCS6Lim%wgt

      FintNNLO_qg(7) = respdf_tmp(1,1)

      !! --------------------------- TC + C6 + S6 -------------------------- !!

      z1 = TCC6S6Lim%Lim_z(3)
      s26 = TCC6S6Lim%Lim_sij(2,6)
      s56 = TCC6S6Lim%Lim_sij(5,6)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (4*Tr/s26)        &
                      * (2/s56) * Pqq(z1) &
                      * TCC6S6Lim%wgt
      respdf_tmp(:,2) = -respdf_tmp(:,2)

      FintNNLO_qg(8) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)
      kin(6) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
 !   if (any(FintNNLO_qg .ne. zero)) then
    !    if (FintNNLO_qg(1) .ne. zero) then
    if (product(FintNNLO_qg) .ne. zero) then
       print *, "FintNNLO_qg",FintNNLO_qg
       pause
    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_qg)

#if(_withchecks == 1)
    FintNNLO_rr_tc_qg = FintNNLO_qg
#endif

  end function xsect_nnlo_rr_5262bd_qg

  function xsect_nnlo_rr_5261_qg(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_dc_ii
    implicit none
    integer :: xsect_nnlo_rr_5261_qg
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: C5Lim, S6Lim, C6Lim, C6S6Lim, C5C6Lim, C5S6Lim, C5C6S6Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_qg(8),kin(1:6)
    real(dp)  :: respdf(ipdf),respdf_tmp(ipdf,2)
    real(dp)  :: res_nlo_old(2,2),res_tmp(2,2)
    real(dp) :: res_nnlo(-5:7,-5:7), res_nlo(-5:7,-5:7), res_lo(-5:7,-5:7), res_nlo_eikqed(-5:7,-5:7),res_lo_eikqed(-5:7,-5:7), res_lo_ischarges(-5:7,-5:7), res_nlo_ischarges(-5:7,-5:7)
    real(dp)  :: eik_qed(4)
    real(dp)  :: z5,z6,s25,s16,E6,E6sq,damp
    logical   :: oldcode

    oldcode = .false.

    xsect_nnlo_rr_5261_qg = 0

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
         HardProc, C5Lim=C5Lim, S6Lim=S6Lim, C6Lim=C6Lim, & 
         C6S6Lim=C6S6Lim, C5C6Lim=C5C6Lim, C5S6Lim=C5S6Lim, C5C6S6Lim=C5C6S6Lim, &
         opt_etas=[3,4])

#if (_Vcharge == 0)
    HardProc%ids(1:6) = [id_q,id_g,  id_el,-id_el,id_q,id_a]
    S6Lim%ids(1:5)    = [id_q,id_g,  id_el,-id_el,id_q]
    C6Lim%ids(1:5)    = [id_q,id_g,  id_el,-id_el,id_q]
    C5Lim%ids(1:5)    = [id_q,-id_qp,id_el,-id_el,id_a]
    C5S6Lim%ids(1:4)  = [id_q,-id_qp,id_el,-id_el]
    C5C6Lim%ids(1:4)  = [id_q,-id_qp,id_el,-id_el]
#elif  (_Vcharge == -1)
    HardProc%ids(1:6) = [id_q,id_g,  id_el,-id_nue,id_q,id_a]
    S6Lim%ids(1:5)    = [id_q,id_g,  id_el,-id_nue,id_q]
    C6Lim%ids(1:5)    = [id_q,id_g,  id_el,-id_nue,id_q]
    C5Lim%ids(1:5)    = [id_q,-id_qp,id_el,-id_nue,id_a]
    C5S6Lim%ids(1:4)  = [id_q,-id_qp,id_el,-id_nue]
    C5C6Lim%ids(1:4)  = [id_q,-id_qp,id_el,-id_nue]
#elif  (_Vcharge == +1)
    HardProc%ids(1:6) = [id_q,id_g,  id_nue,-id_el,id_q,id_a]
    S6Lim%ids(1:5)    = [id_q,id_g,  id_nue,-id_el,id_q]
    C6Lim%ids(1:5)    = [id_q,id_g,  id_nue,-id_el,id_q]
    C5Lim%ids(1:5)    = [id_q,-id_qp,id_nue,-id_el,id_a]
    C5S6Lim%ids(1:4)  = [id_q,-id_qp,id_nue,-id_el]
    C5C6Lim%ids(1:4)  = [id_q,-id_qp,id_nue,-id_el]
#endif
    HardProc%npart = 6
    S6Lim%npart    = 5
    C5Lim%npart    = 5
    C6Lim%npart    = 5
    C5S6Lim%npart  = 4
    C5C6Lim%npart  = 4


    !!-----------------------------------------------------------------------!!


    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_qg(1) = zero

    else

       if (oldcode) then
          call res_tree_ga_qg(HardProc%AmpMom,res_tmp)
          call get_respdf(qg_lumi,1,1,HardProc,res_tmp,respdf)
       else
          call res_tree_ga_qg_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,HardProc,res_nnlo,respdf)
       endif

      call partition_nnlo_fact(HardProc,damp,iconf_qed=[1,5,3,4,6],i_qed=1)

      respdf = respdf * HardProc%wgt * damp

      kin(1) = respdf(1)
      FintNNLO_qg(1) = kin(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 6                                !!
    !!-----------------------------------------------------------------------!!


    call cut_histo(S6Lim)

    if (S6Lim%makecut) then

      kin(2) = zero
      FintNNLO_qg(2:3) = zero

   else
      
       E6   = S6Lim%Lim_KinInv(1)
       E6sq = E6*E6
       if (oldcode) then
          call res_tree_g_qg(S6Lim%AmpMom,res_tmp)          
          !! --------------------------------- S6 --------------------------------- !!
          call get_qed_eik(charges,S6Lim%Lim_etaij,ql_pos,6,eik_qed)
          eik_qed = eik_qed/E6sq          
          res_nlo_old(:,1) = res_tmp(:,1) * eik_qed(1:2)
          res_nlo_old(:,2) = res_tmp(:,2) * eik_qed(3:4)          
          call get_respdf(qg_lumi,1,1,S6Lim,res_nlo_old,respdf_tmp(:,1))
       else
          call res_tree_g_qg_gen(S6Lim%AmpMom,res_nlo)         
          call get_qed_eik_gen(res_nlo,S6Lim%Lim_etaij,[1,3,4,5],6,res_nlo_eikqed)
          call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf_tmp(:,1))
          respdf_tmp = respdf_tmp / E6sq
       endif

      call partition_nnlo_fact(S6Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=1)

      respdf_tmp(:,1) =  respdf_tmp(:,1) * S6Lim%wgt * damp
      respdf_tmp(:,1) = -respdf_tmp(:,1)
      FintNNLO_qg(2)  =  respdf_tmp(1,1)

    !! ------------------------------ C6 + S6 ------------------------------ !!
      C6S6Lim%mur = S6Lim%mur
      C6S6Lim%muf = S6Lim%muf

      z6  = C6S6Lim%Lim_z(2)
      s16 = C6S6Lim%Lim_sij(1,6)

      if (oldcode) then
         res_nlo_old(:,1) = res_tmp(:,1) * Qdn2
         res_nlo_old(:,2) = res_tmp(:,2) * Qup2
         call get_respdf(qg_lumi,1,1,C6S6Lim,res_nlo_old,respdf_tmp(:,2))
      else
         res_nlo_ischarges = multiply_IS_charges_sq(res_nlo,1)
         call get_respdf_gen(1,1,C6S6Lim,res_nlo_ischarges,respdf_tmp(:,2))
      endif
      
      call partition_nnlo_fact(C6S6Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=1)

      respdf_tmp(:,2) = respdf_tmp(:,2) * (4/z6/s16) * C6S6Lim%wgt * damp
      FintNNLO_qg(3)  = respdf_tmp(1,2)

      respdf = respdf_tmp(1,1) + respdf_tmp(1,2)
      kin(2) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 6                              !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(3) = zero
      FintNNLO_qg(4) = zero

    else

      z6  = C6Lim%Lim_z(2)
      s16 = C6Lim%Lim_sij(1,6)

      if (oldcode) then
         call res_tree_g_qg(C6Lim%AmpMom,res_tmp)
         res_tmp(:,1) = res_tmp(:,1) * Qdn2
         res_tmp(:,2) = res_tmp(:,2) * Qup2         
         call get_respdf(qg_lumi,1,1,C6Lim,res_tmp,respdf)
      else
         call res_tree_g_qg_gen(C6Lim%AmpMom,res_nlo)
         res_nlo_ischarges = multiply_IS_charges_sq(res_nlo,1)
         call get_respdf_gen(1,1,C6Lim,res_nlo_ischarges,respdf)
      endif

      call partition_nnlo_fact(C6Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=1)

      respdf = (-two/s16) * Pqq(z6) * respdf * C6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_qg(4) = respdf(1)

      kin(3) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                         Collinear 5 + Soft 6                          !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(C5S6Lim)

    if (C5S6Lim%makecut) then

      kin(4) = zero
      FintNNLO_qg(5:6) = zero

    else

       E6   = C5S6Lim%Lim_KinInv(1)
       E6sq = E6*E6

       if (oldcode) then
          call res_tree_qqb(C5S6Lim%AmpMom,res_tmp)
          !! ------------------------------ C5 + S6 ------------------------------ !!
          call get_qed_eik(charges,C5S6Lim%Lim_etaij,ql_pos,6,eik_qed)
          eik_qed = eik_qed/E6sq
          res_nlo_old(:,1) = res_tmp(:,1) * eik_qed(1:2)
          res_nlo_old(:,2) = res_tmp(:,2) * eik_qed(3:4)
          call get_respdf(qg_lumi,1,1,C5S6Lim,res_nlo_old,respdf_tmp(:,1))
       else
          call res_tree_qqb_gen(C5S6Lim%AmpMom,res_lo)
          res_lo = transition('none','g -> q', res_lo)
          call get_qed_eik_gen(res_lo,C5S6Lim%Lim_etaij,[1,3,4,5],6,res_lo_eikqed)
          call get_respdf_gen(1,1,C5S6Lim,res_lo_eikqed,respdf_tmp(:,1))
          respdf_tmp(:,1)  = respdf_tmp(:,1)/E6sq
       endif

      z5  = C5S6Lim%Lim_z(1)
      s25 = C5S6Lim%Lim_sij(2,5)



      call partition_nnlo_fact(C5S6Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=1)

      respdf_tmp(:,1) = respdf_tmp(:,1) &
                      * (two/s25) * Tr * Pqg(z5) * C5S6Lim%wgt * damp

      FintNNLO_qg(5) = respdf_tmp(1,1)

    !! --------------------------- C5 + C6 + S6 ---------------------------- !!
      C5C6S6Lim%mur = C5S6Lim%mur
      C5C6S6Lim%muf = C5S6Lim%muf

      z5  = C5C6S6Lim%Lim_z(1)
      s25 = C5C6S6Lim%Lim_sij(2,5)
      z6  = C5C6S6Lim%Lim_z(2)
      s16 = C5C6S6Lim%Lim_sij(1,6)

      if (oldcode) then
         res_nlo_old(:,1) = res_tmp(:,1) * Qdn2
         res_nlo_old(:,2) = res_tmp(:,2) * Qup2
         call get_respdf(qg_lumi,1,1,C5C6S6Lim,res_nlo_old,respdf_tmp(:,2))
      else
         res_lo_ischarges = multiply_IS_charges_sq(res_lo,1)
         call get_respdf_gen(1,1,C5C6S6Lim,res_lo_ischarges,respdf_tmp(:,2))
      endif

      !-- damp = one
      respdf_tmp(:,2) = respdf_tmp(:,2) &
                      * (4/z6/s16)      &
                      * (two/s25) * Tr * Pqg(z5) &
                      * C5C6S6Lim%wgt
      respdf_tmp(:,2) = - respdf_tmp(:,2)

      FintNNLO_qg(6) = respdf_tmp(1,2)

      respdf = respdf_tmp(1,1) + respdf_tmp(1,2)
      kin(4) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 5                              !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(C5Lim)

    if (C5Lim%makecut) then

      kin(5) = zero
      FintNNLO_qg(7) = zero

    else

      z5  = C5Lim%Lim_z(1)
      s25 = C5Lim%Lim_sij(2,5)

      if (oldcode) then
         call res_tree_a_qqb(C5Lim%AmpMom,res_tmp)
         call get_respdf(qg_lumi,1,1,C5Lim,res_tmp,respdf)
      else
         call res_tree_a_qqb_gen(C5Lim%AmpMom,res_nlo)
         res_nlo = transition('none','g -> q', res_nlo)
         call get_respdf_gen(1,1,C5Lim,res_nlo,respdf)
      endif
      
      call partition_nnlo_fact(C5Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=1)

      respdf = (two/s25) * Tr * Pqg(z5) * respdf * C5Lim%wgt * damp
      respdf = -respdf

      FintNNLO_qg(7) = respdf(1)

      kin(5) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                       Collinear 5 + Collinear 6                       !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(C5C6Lim)

    if (C5C6Lim%makecut.or.C5C6Lim%flag) then

      kin(6) = zero
      FintNNLO_qg(8) = zero

    else

      z5  = C5C6Lim%Lim_z(1)
      s25 = C5C6Lim%Lim_sij(2,5)
      z6  = C5C6Lim%Lim_z(2)
      s16 = C5C6Lim%Lim_sij(1,6)

      if (oldcode) then
         call res_tree_qqb(C5C6Lim%AmpMom,res_tmp)
         res_tmp(:,1) = res_tmp(:,1) * Qdn2
         res_tmp(:,2) = res_tmp(:,2) * Qup2         
         call get_respdf(qg_lumi,1,1,C5C6Lim,res_tmp,respdf)
      else
         call res_tree_qqb_gen(C5C6Lim%AmpMom,res_lo)
         res_lo = transition('none','g -> q', res_lo)
         res_lo_ischarges = multiply_IS_charges_sq(res_lo,1)
         call get_respdf_gen(1,1,C5C6Lim,res_lo_ischarges,respdf)
      endif

      !-- damp = one
      respdf = respdf &
             * ( two/s25) * Tr * Pqg(z5) &
             * (-two/s16) * Pqq(z6) &
             * C5C6Lim%wgt

      FintNNLO_qg(8) = respdf(1)

      kin(6) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_qg)

#if(_withchecks == 1)
    FintNNLO_rr_dc_qg = FintNNLO_qg
#endif

  end function xsect_nnlo_rr_5261_qg

  function xsect_nnlo_rr_526k_qg(yRnd,ff,vegasweight,k,l)
    use mod_kinematics_nnlo_dc_if
    implicit none
    integer :: xsect_nnlo_rr_526k_qg, k, l
    real(dp15) :: yRnd(30), ff(1), vegasweight
    !--
    integer, parameter :: imax_ipdf = 2, imax_ilim = 2
    type(KinConfig) :: HardProc
    type(KinConfig) :: C6Lim, S6Lim, C6S6Lim, C5Lim, C5C6Lim, C5S6Lim, C5C6S6Lim
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_qg(8),kin(1:6)
    real(dp)  :: respdf(ipdf),respdf_vect(imax_ipdf,ipdf)
    real(dp)  :: respdf_vect_part(imax_ilim,ipdf),respdf_vect1(ipdf), respdf_vect2(ipdf)
    real(dp)  :: res_tmp(2,2),res_tmp_vect(2,2,imax_ipdf)
    real(dp) :: res_nnlo(-5:7,-5:7), res_nlo(-5:7,-5:7), res_lo(-5:7,-5:7), res_nlo_eikqed(-5:7,-5:7),  res_lo_eikqed(-5:7,-5:7), res_lo_ischarges(-5:7,-5:7), res_nlo_ischarges(-5:7,-5:7)
    real(dp)  :: eik_qed(4),Qsq_FS(2)
    real(dp)  :: z5,z6,s25,sk6,E6,E6sq,damp
    logical   :: oldcode

    oldcode = .false.

    xsect_nnlo_rr_526k_qg = 0

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
    call kinematics_nnlo_if_5i6k(&
    xx(1:kNNLO_max_full), 2, 1, k, l, &
    HardProc, S6Lim=S6Lim, C6Lim=C6Lim, C6S6Lim=C6S6Lim,  &
    C5Lim=C5Lim, C5S6Lim=C5S6Lim, C5C6Lim=C5C6Lim, C5C6S6Lim=C5C6S6Lim,  &
    opt_etas=[3,4])

    Qsq_FS = [Q3**2, Q4**2]

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_q,id_a]
    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_qg(1) = zero

    else

       if (oldcode) then
          call res_tree_ga_qg(HardProc%AmpMom,res_tmp)
          call get_respdf(qg_lumi,1,1,HardProc,res_tmp,respdf)
       else
          call res_tree_ga_qg_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,HardProc,res_nnlo,respdf)
       endif

      call partition_nnlo_fact(HardProc,damp,iconf_qed=[1,5,3,4,6],i_qed=k)

      respdf = respdf * HardProc%wgt * damp

      kin(1) = respdf(1)
      FintNNLO_qg(1) = kin(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 6                                !!
    !!-----------------------------------------------------------------------!!
    S6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    call cut_histo(S6Lim)

    if (S6Lim%makecut) then

      kin(2) = zero
      FintNNLO_qg(2:3) = zero

    else

      !-- prepare PDFs structures, S6
      E6   = S6Lim%Lim_KinInv(1)
      E6sq = E6*E6

      if (oldcode) then
         call res_tree_g_qg(S6Lim%AmpMom,res_tmp)         
         call get_qed_eik(charges,S6Lim%Lim_etaij,ql_pos,6,eik_qed)
         eik_qed = eik_qed/E6sq
         res_tmp_vect(:,1,1) = res_tmp(:,1) * eik_qed(1:2)
         res_tmp_vect(:,2,1) = res_tmp(:,2) * eik_qed(3:4)
         !-- prepare PDFs structures, C6S6
         res_tmp_vect(:,:,2) = res_tmp(:,:) * Q_lep2
         call get_respdf_vect(qg_lumi,1,1,S6Lim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
      else
         call res_tree_g_qg_gen(S6Lim%AmpMom,res_nlo)
         call get_qed_eik_gen(res_nlo,S6Lim%Lim_etaij,[1,3,4,5],6,res_nlo_eikqed)
         call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf_vect1)
         respdf_vect1 = respdf_vect1 / E6sq                    
         res_nlo = Qsq_Fs(k-2)**2 * res_nlo
         call get_respdf_gen(1,1,S6Lim,res_nlo,respdf_vect2)
         respdf_vect(1,:) = respdf_vect1(:)
         respdf_vect(2,:) = respdf_vect2(:)
      endif

    !! --------------------------------- S6 -------------------------------- !!

      call partition_nnlo_fact(S6Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=k)

      respdf_vect_part(1,:) =  - respdf_vect(1,:) * S6Lim%wgt * damp
      FintNNLO_qg(2) = respdf_vect_part(1,1)

    !! ------------------------------ C6 + S6 ------------------------------ !!
      z6  = C6S6Lim%Lim_z(2)
      sk6 = C6S6Lim%Lim_sij(k,6)

      call partition_nnlo_fact(C6S6Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=k)

      if (oldcode) then
         respdf_vect_part(2,:) = respdf_vect(2,:) * (4*Q_lep2/z6/sk6) * C6S6Lim%wgt * damp
      else
         respdf_vect_part(2,:) = respdf_vect(2,:) * (4/z6/sk6) * C6S6Lim%wgt * damp
      endif
         
      FintNNLO_qg(3)  = respdf_vect_part(2,1)

      respdf = sum(respdf_vect_part(1:2,:),1)

      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !-----------------------------------------------------------------------!!
    !                              Collinear 6                              !!
    !-----------------------------------------------------------------------!!

    C6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(3) = zero
      FintNNLO_qg(4) = zero

    else

      z6  = C6Lim%Lim_z(2)
      sk6 = C6Lim%Lim_sij(k,6)

      if (oldcode) then
         call res_tree_g_qg(C6Lim%AmpMom,res_tmp)
         call get_respdf(qg_lumi,1,1,C6Lim,res_tmp,respdf)
      else
         call res_tree_g_qg_gen(C6Lim%AmpMom,res_nlo)
         res_nlo = Qsq_Fs(k-2)**2 * res_nlo
         call get_respdf_gen(1,1,C6Lim,res_nlo,respdf)
         
      endif
      call partition_nnlo_fact(C6Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=k)

      respdf = respdf * 2/sk6 * Q_lep2 * Pqg(z6) * C6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_qg(4) = respdf(1)

      kin(3) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 5                              !!
    !!-----------------------------------------------------------------------!!

    C5Lim%ids = [0,0,id_el,-id_el,id_a,0]
    call cut_histo(C5Lim)

    if (C5Lim%makecut) then

      kin(4) = zero
      FintNNLO_qg(5) = zero

    else

      z5  = C5Lim%Lim_z(1)
      s25 = C5Lim%Lim_sij(2,5)

      if (oldcode) then
         call res_tree_a_qqb(C5Lim%AmpMom,res_tmp)
         call get_respdf(qg_lumi,1,1,C5Lim,res_tmp,respdf)
      else
         call res_tree_a_qqb_gen(C5Lim%AmpMom,res_nlo)
         res_nlo = transition('none','g -> q', res_nlo)
         call get_respdf_gen(1,1,C5Lim,res_nlo,respdf)
      endif

      call partition_nnlo_fact(C5Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=k)

      respdf = respdf * (2/s25) * Tr * Pqq(z5) * C5Lim%wgt * damp
      respdf = -respdf

      FintNNLO_qg(5) = respdf(1)

      kin(4) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                         Collinear 5 + Soft 6                          !!
    !!-----------------------------------------------------------------------!!

    C5S6Lim%ids = [0,0,id_el,-id_el,0,0]
    call cut_histo(C5S6Lim)

    if (C5S6Lim%makecut) then

      kin(5) = zero
      FintNNLO_qg(6:7) = zero

    else

      E6   = C5S6Lim%Lim_KinInv(1)
      E6sq = E6*E6

      if (oldcode) then
         call res_tree_qqb(C5S6Lim%AmpMom,res_tmp)
         
         !-- prepare PDFs structures, C5S6
         call get_qed_eik(charges,C5S6Lim%Lim_etaij,ql_pos,6,eik_qed)
         eik_qed = eik_qed/E6sq
         res_tmp_vect(:,1,1) = res_tmp(:,1) * eik_qed(1:2)
         res_tmp_vect(:,2,1) = res_tmp(:,2) * eik_qed(3:4)
         
         !-- prepare PDFs structures, C5C6S6
         res_tmp_vect(:,:,2) = res_tmp(:,:) * Q_lep2
         
         call get_respdf_vect(qg_lumi,1,1,C5S6Lim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
      else
         call res_tree_qqb_gen(C5S6Lim%AmpMom,res_lo)
         res_lo = transition('none','g -> q', res_lo)
         call get_qed_eik_gen(res_lo,C5S6Lim%Lim_etaij,[1,3,4,5],6,res_lo_eikqed)
         call get_respdf_gen(1,1,C5S6Lim,res_lo_eikqed,respdf_vect1)
         respdf_vect1 = respdf_vect1 / E6sq                    
         res_lo = Qsq_Fs(k-2)**2 * res_lo
         call get_respdf_gen(1,1,C5S6Lim,res_lo,respdf_vect2)
         respdf_vect(1,:) = respdf_vect1(:)
         respdf_vect(2,:) = respdf_vect2(:)
      endif

         

    !! ------------------------------ C5 + S6 ------------------------------ !!

      z5  = C5S6Lim%Lim_z(1)
      s25 = C5S6Lim%Lim_sij(2,5)

      call partition_nnlo_fact(C5S6Lim,damp,iconf_qed=[1,5,3,4,6],i_qed=k)

      respdf_vect_part(1,:) = respdf_vect(1,:) &
                            * (2/s25) * Tr * Pqq(z5) * C5S6Lim%wgt * damp

      FintNNLO_qg(6) = respdf_vect_part(1,1)

    !! --------------------------- C5 + C6 + S6 ---------------------------- !!

      z5  = C5C6S6Lim%Lim_z(1)
      s25 = C5C6S6Lim%Lim_sij(2,5)
      z6  = C5C6S6Lim%Lim_z(2)
      sk6 = C5C6S6Lim%Lim_sij(k,6)

      !-- damp = one
      if (oldcode) then
         respdf_vect_part(2,:) = - respdf_vect(2,:) &
              * (4*Q_lep2/z6/sk6) &
              * (two/s25) * Tr * Pqq(z5) &
              * C5C6S6Lim%wgt
      else
         respdf_vect_part(2,:) = - respdf_vect(2,:) &
              * (4/z6/sk6) &
              * (two/s25) * Tr * Pqq(z5) &
              * C5C6S6Lim%wgt
      endif

      FintNNLO_qg(7) = respdf_vect_part(2,1)

      respdf = sum(respdf_vect_part(1:2,:),1)
      kin(5) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                       Collinear 5 + Collinear 6                       !!
    !!-----------------------------------------------------------------------!!

    C5C6Lim%ids = [0,0,id_el,-id_el,0,0]
    call cut_histo(C5C6Lim)

    if (C5C6Lim%makecut) then

      kin(6) = zero
      FintNNLO_qg(8) = zero

    else

      z5  = C5C6Lim%Lim_z(1)
      s25 = C5C6Lim%Lim_sij(2,5)
      z6  = C5C6Lim%Lim_z(2)
      sk6 = C5C6Lim%Lim_sij(k,6)

      if (oldcode) then
         call res_tree_qqb(C5C6Lim%AmpMom,res_tmp)
         call get_respdf(qg_lumi,1,1,C5C6Lim,res_tmp,respdf)
      else
         call res_tree_qqb_gen(C5C6Lim%AmpMom,res_lo)
         res_lo = transition('none','g -> q', res_lo)
         res_lo = Qsq_Fs(k-2)**2 * res_lo
         call get_respdf_gen(1,1,C5C6Lim,res_lo,respdf)
      endif
         
         

      !-- damp = one
      if (oldcode) then
         respdf = respdf &
              * (2/s25) * Tr * Pqq(z5)  &
              * (2/sk6) * Q_lep2 * Pqg(z6) &
              * C5C6Lim%wgt
      else
         respdf = respdf &
              * (2/s25) * Tr * Pqq(z5)  &
              * (2/sk6) * Pqg(z6) &
              * C5C6Lim%wgt
      endif
      
      FintNNLO_qg(8) = respdf(1)

      kin(6) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    if (product(FintNNLO_qg) .ne. zero) then
!    if (any(FintNNLO_qg .ne. zero)) then
       print *, "FintNNLO_qg",FintNNLO_qg
       pause
    endif


    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_qg)

#if(_withchecks == 1)
    FintNNLO_rr_dc_qg = FintNNLO_qg
#endif

  end function xsect_nnlo_rr_526k_qg

end module mod_xsects_nnlo_rr_qg

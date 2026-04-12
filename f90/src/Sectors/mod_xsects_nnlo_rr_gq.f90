module mod_xsects_nnlo_rr_gq
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
  real(dp), public, save :: FintNNLO_rr_tc_gq(8), FintNNLO_rr_dc_gq(8)
#endif

  !! Electric charges. 
  !! Each column is for the relevant res_jj(i,j)
  !! NB: different assignment wrt to matrix element
  real(dp), parameter :: charges(4,4) = reshape(& 
    [Qdn,  -Qdn,   Qup,  -Qup,   &
     Q_lep, Q_lep, Q_lep, Q_lep, &
    -Q_lep,-Q_lep,-Q_lep,-Q_lep, &
    -Qdn,   Qdn,  -Qup,   Qup  ],&
    [4,4])

  !! positions of charged quarks and leptons (ql)
  integer, parameter :: ql_pos(4) = [2,3,4,5]

  public :: xsect_nnlo_rr_5161a_gq,xsect_nnlo_rr_5161c_gq
  public :: xsect_nnlo_rr_5161b_gq,xsect_nnlo_rr_5161d_gq
  public :: xsect_nnlo_rr_5162_gq
  public :: xsect_nnlo_rr_5163_gq,xsect_nnlo_rr_5164_gq

  private

  !! -- g(p1) + qb/q(p2) -> l-(p3) l+(p4) qb/q(p5) a(p6)

contains

  function xsect_nnlo_rr_5161b_gq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161b_gq
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5161b_gq = xsect_nnlo_rr_5161bd_gq(&
      yRnd,ff,vegasweight,2)

  end function xsect_nnlo_rr_5161b_gq

  function xsect_nnlo_rr_5161d_gq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161d_gq
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5161d_gq = xsect_nnlo_rr_5161bd_gq(&
      yRnd,ff,vegasweight,4)

  end function xsect_nnlo_rr_5161d_gq

  function xsect_nnlo_rr_5163_gq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5163_gq
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5163_gq = xsect_nnlo_rr_516k_gq(&
      yRnd,ff,vegasweight,3,4)

  end function xsect_nnlo_rr_5163_gq

  function xsect_nnlo_rr_5164_gq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5164_gq
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5164_gq = xsect_nnlo_rr_516k_gq(&
      yRnd,ff,vegasweight,4,3)

  end function xsect_nnlo_rr_5164_gq

  !!*************************************************************************!!

  function xsect_nnlo_rr_5161a_gq(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5161a_gq
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: S6Lim, TCLim, TCS6Lim
    !--
    real(dp) :: xx(kNNLO_max_full)
    real(dp) :: FintNNLO_gq(4),kin(4)
    real(dp) :: respdf(ipdf)
    real(dp) :: res_tmp(2,2)
    real(dp) :: res_nnlo(-5:7,-5:7), res_nlo(-5:7,-5:7), res_lo(-5:7,-5:7), res_nlo_eikqed(-5:7,-5:7), res_lo_ischarges(-5:7,-5:7)
    real(dp) :: eik_qed(4),z1,z5,z6,s15,s16,s56,E6,E6sq
    real(dp) :: damp
    logical  :: oldcode


    oldcode = .false.
    

    xsect_nnlo_rr_5161a_gq = 0

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
         HardProc, S6Lim=S6Lim, TCLim=TCLim, TCS6Lim=TCS6Lim, opt_etas=[3,4])

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_q,id_a]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_gq(1) = zero

    else

       if (oldcode) then
          call res_tree_ga_gq(HardProc%AmpMom,res_tmp)
          call get_respdf(gq_lumi,1,1,HardProc,res_tmp,respdf)
          
       else
          call res_tree_ga_gq_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,HardProc,res_nnlo,respdf)
       endif
       
      call partition_nnlo_fact(HardProc,damp,iconf_qed=[5,2,3,4,6],i_qed=1)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_gq(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 6                                !!
    !!-----------------------------------------------------------------------!!

    S6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    S6Lim%npart = 5

    call cut_histo(S6Lim)

    if (S6Lim%makecut) then

      kin(2) = zero
      FintNNLO_gq(2) = zero

    else

       E6   = S6Lim%Lim_KinInv(1)
       E6sq = E6*E6
       
       if (oldcode) then
          call res_tree_g_gq(S6Lim%AmpMom,res_tmp)
          
          call get_qed_eik(charges,S6Lim%Lim_etaij,ql_pos,6,eik_qed)
          eik_qed = eik_qed/E6sq
          res_tmp(1:2,1) = res_tmp(1:2,1) * eik_qed(1:2)
          res_tmp(1:2,2) = res_tmp(1:2,2) * eik_qed(3:4)
          call get_respdf(gq_lumi,1,1,S6Lim,res_tmp,respdf)

       else
          call res_tree_g_gq_gen(S6Lim%AmpMom,res_nlo)
          call get_qed_eik_gen(res_nlo,S6Lim%Lim_etaij,[2,3,4,5],6,res_nlo_eikqed)
          call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf)
          respdf = respdf / E6sq
       endif
          


      call partition_nnlo_fact(S6Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=1)

      respdf =  respdf * S6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_gq(2) = respdf(1)
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
      FintNNLO_gq(3) = zero

    else

       if (oldcode) then
          call res_tree_qqb(TCLim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2
          call get_respdf(gq_lumi,1,1,TCLim,res_tmp,respdf)
       else
          call res_tree_qqb_gen(TCLim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,1)
          res_lo_ischarges = transition('g -> q', 'none', res_lo_ischarges)
          call get_respdf_gen(1,1,TCLim,res_lo_ischarges,respdf)
       endif
          
      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s15 = TCLim%Lim_sij(1,5)
      s16 = TCLim%Lim_sij(1,6)

      !-- damp = one
      respdf = respdf &
             * Tr * (-Pgaq(-s16,s56,-s15,z6,z1,z5)) * TCLim%wgt
      respdf = -respdf

      FintNNLO_gq(3) = respdf(1)
      kin(3) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                      Triple Collinear 5 + Soft 6                      !!
    !!-----------------------------------------------------------------------!!

    TCS6Lim%ids = [0,0,id_el,-id_el,0,0]
    TCS6Lim%npart = 4

    call cut_histo(TCS6Lim)

    if (TCS6Lim%makecut) then

      kin(4) = zero
      FintNNLO_gq(4) = zero

    else

       if (oldcode) then
          call res_tree_qqb(TCS6Lim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2
          call get_respdf(gq_lumi,1,1,TCS6Lim,res_tmp,respdf)
       else
          call res_tree_qqb_gen(TCS6Lim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,1)
          res_lo_ischarges = transition('g -> q', 'none', res_lo_ischarges)
          call get_respdf_gen(1,1,TCS6Lim,res_lo_ischarges,respdf)
       endif
       
      z5  = TCS6Lim%Lim_z(1)
      z1  = TCS6Lim%Lim_z(3)
      s56 = TCS6Lim%Lim_sij(5,6)
      s15 = TCS6Lim%Lim_sij(1,5)

      !-- damp = one
      respdf = respdf * (4*Tr/z5/s56)  &
             * (2/s15) * Pqg(z1)  &
             * TCS6Lim%wgt

      FintNNLO_gq(4) = respdf(1)
      kin(4) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_gq)

#if(_withchecks == 1)
    limval_nnlo = zero
    limval_nnlo(1) = FintNNLO_gq(1)    ! hard
    limval_nnlo(3) = FintNNLO_gq(2)    ! S6
    limval_nnlo(4) = FintNNLO_gq(3)    ! TC
    limval_nnlo(9) = FintNNLO_gq(4)    ! TC+S6
    
    FintNNLO_rr_tc_gq(1:4) = FintNNLO_gq
    FintNNLO_rr_tc_gq(5:8) = zero
#endif

  end function xsect_nnlo_rr_5161a_gq

  function xsect_nnlo_rr_5161c_gq(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5161c_gq
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: S6Lim, TCLim, TCS6Lim, C5Lim, C5S6Lim, TCC5Lim, TCC5S6Lim
    !--
    real(dp) :: xx(kNNLO_max_full)
    real(dp) :: FintNNLO_gq(8),kin(6)
    real(dp) :: respdf(ipdf),respdf_tmp(ipdf,2),respdf_bak(ipdf)
    real(dp) :: res_nnlo(-5:7,-5:7), res_nlo(-5:7,-5:7), res_lo(-5:7,-5:7), res_nlo_eikqed(-5:7,-5:7), res_lo_eikqed(-5:7,-5:7), res_lo_ischarges(-5:7,-5:7)
    real(dp) :: res_tmp(2,2)
    real(dp) :: eik_qed(4),z1,z5,z6,s15,s16,s56,s156,E6,E6sq
    real(dp) :: damp
    logical  :: oldcode


    oldcode = .false.

    xsect_nnlo_rr_5161c_gq = 0

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
         HardProc, S6Lim=S6Lim, TCLim=TCLim, TCC5Lim=TCC5Lim, &
         TCS6Lim=TCS6Lim, TCC5S6Lim=TCC5S6Lim, C5Lim=C5Lim, C5S6Lim=C5S6Lim, &
         opt_etas=[3,4])

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_q,id_a]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_gq(1) = zero

    else

       if (oldcode) then
          call res_tree_ga_gq(HardProc%AmpMom,res_tmp)
          call get_respdf(gq_lumi,1,1,HardProc,res_tmp,respdf)
       else
          call res_tree_ga_gq_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,HardProc,res_nnlo,respdf)
       endif
       call partition_nnlo_fact(HardProc,damp,iconf_qed=[5,2,3,4,6],i_qed=1)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_gq(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 6                                !!
    !!-----------------------------------------------------------------------!!

    S6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    S6Lim%npart = 5

    call cut_histo(S6Lim)

    if (S6Lim%makecut) then

      kin(2) = zero
      FintNNLO_gq(2) = zero

   else

      E6   = S6Lim%Lim_KinInv(1)
      E6sq = E6*E6
      
      if (oldcode) then
         call res_tree_g_gq(S6Lim%AmpMom,res_tmp)
         call get_qed_eik(charges,S6Lim%Lim_etaij,ql_pos,6,eik_qed)
         eik_qed = eik_qed/E6sq
         res_tmp(1:2,1) = res_tmp(1:2,1) * eik_qed(1:2) ! g dx and g d
         res_tmp(1:2,2) = res_tmp(1:2,2) * eik_qed(3:4) ! g ux and g u
         call get_respdf(gq_lumi,1,1,S6Lim,res_tmp,respdf)
      else
         call res_tree_g_gq_gen(S6Lim%AmpMom,res_nlo)
         call get_qed_eik_gen(res_nlo,S6Lim%Lim_etaij,[2,3,4,5],6,res_nlo_eikqed)
         call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf)
         respdf = respdf / E6sq
      endif
         
      call partition_nnlo_fact(S6Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=1)

      respdf =  respdf * S6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_gq(2) = respdf(1)
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
      FintNNLO_gq(3:4) = zero

    else

       if (oldcode) then
          call res_tree_qqb(TCLim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2
          call get_respdf(gq_lumi,1,1,TCLim,res_tmp,respdf_bak)
       else
          call res_tree_qqb_gen(TCLim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,1)
          res_lo_ischarges = transition('g -> q', 'none', res_lo_ischarges)
          call get_respdf_gen(1,1,TCLim,res_lo_ischarges,respdf_bak)
       endif

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s15 = TCLim%Lim_sij(1,5)
      s16 = TCLim%Lim_sij(1,6)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                 * Tr * (-Pgaq(-s16,s56,-s15,z6,z1,z5)) * TCLim%wgt
      respdf_tmp(:,1) = -respdf_tmp(:,1)

      FintNNLO_gq(3) = respdf_tmp(1,1)

      !! ----------------------------- TC + C5 ----------------------------- !!

      !! (- si5 - si6 + s56) in the x4 -> 0 limit
      s156 = TCC5Lim%Lim_KinInv(1)
      z5  = TCC5Lim%Lim_z(1)
      z6  = TCC5Lim%Lim_z(2)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (2/s156) * Tr * Pqg(z6)  &
                      * (-2/s15) * Pqq(z5) &
                      * TCC5Lim%wgt

      FintNNLO_gq(4) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)
      kin(3) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                      Triple Collinear 5 + Soft 6                      !!
    !!-----------------------------------------------------------------------!!

    TCS6Lim%ids = [0,0,id_el,-id_el,0,0]
    TCS6Lim%npart = 4

    call cut_histo(TCS6Lim)

    if (TCS6Lim%makecut) then

      kin(4) = zero
      FintNNLO_gq(5:6) = zero

    else

       if (oldcode) then
          call res_tree_qqb(TCS6Lim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2
          call get_respdf(gq_lumi,1,1,TCS6Lim,res_tmp,respdf_bak)
       else
          call res_tree_qqb_gen(TCS6Lim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,1)
          res_lo_ischarges = transition('g -> q', 'none', res_lo_ischarges)
          call get_respdf_gen(1,1,TCS6Lim,res_lo_ischarges,respdf_bak)
       endif


      z5  = TCS6Lim%Lim_z(1)
      z1  = TCS6Lim%Lim_z(3)
      s56 = TCS6Lim%Lim_sij(5,6)
      s15 = TCS6Lim%Lim_sij(1,5)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak * (4*Tr/z5/s56)   &
             * (2/s15) * Pqg(z1)  &
             * TCS6Lim%wgt

      FintNNLO_gq(5) = respdf_tmp(1,1)

      !! --------------------------- TC + C5 + S6 -------------------------- !!

      !! s56 in the x4 -> 0 limit
      s56 = TCC5S6Lim%Lim_sij(5,6)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (4*Tr/z5/s56)          &
                      * (2/s15) * Pqg(z1)  &
                      * TCC5S6Lim%wgt
      respdf_tmp(:,2) = -respdf_tmp(:,2)

      FintNNLO_gq(6) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) +  respdf_tmp(:,2)
      kin(4) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 5                              !!
    !!-----------------------------------------------------------------------!!

    C5Lim%ids = [0,0,id_el,-id_el,id_a,0]
    C5Lim%npart = 5

    call cut_histo(C5Lim)

    if (C5Lim%makecut) then

      kin(5) = zero
      FintNNLO_gq(7) = zero

    else

      z5  = C5Lim%Lim_z(1)
      s15 = C5Lim%Lim_sij(1,5)

      if (oldcode) then
         call res_tree_a_qqb(C5Lim%AmpMom,res_tmp)
         call get_respdf(gq_lumi,1,1,C5Lim,res_tmp,respdf)
      else
         call res_tree_a_qqb_gen(C5Lim%AmpMom,res_nlo)
         res_nlo = transition('g -> q', 'none', res_nlo)
         call get_respdf_gen(1,1,C5Lim,res_nlo,respdf)
      endif
      call partition_nnlo_fact(C5Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=1)

      respdf = (two/s15) * Tr * Pqg(z5) * respdf * C5Lim%wgt * damp
      respdf = -respdf

      FintNNLO_gq(7) = respdf(1)
      kin(5) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                         Collinear 5 + Soft 6                          !!
    !!-----------------------------------------------------------------------!!

    C5S6Lim%ids = [0,0,id_el,-id_el,0,0]
    C5S6Lim%npart = 4

    call cut_histo(C5S6Lim)

    if (C5S6Lim%makecut) then

      kin(6) = zero
      FintNNLO_gq(8) = zero

    else

       E6   = C5S6Lim%Lim_KinInv(1)
       E6sq = E6*E6
       
       if (oldcode) then
          call res_tree_qqb(C5S6Lim%AmpMom,res_tmp)
          call get_qed_eik(charges,C5S6Lim%Lim_etaij,ql_pos,6,eik_qed)
          eik_qed = eik_qed/E6sq
          res_tmp(1:2,1) = res_tmp(1:2,1) * eik_qed(1:2) ! g dx and g d
          res_tmp(1:2,2) = res_tmp(1:2,2) * eik_qed(3:4) ! g ux and g u
          call get_respdf(gq_lumi,1,1,C5S6Lim,res_tmp,respdf)
       else
          call res_tree_qqb_gen(C5S6Lim%AmpMom,res_lo)
          res_lo = transition('g -> q', 'none', res_lo)
          call get_qed_eik_gen(res_lo,C5S6Lim%Lim_etaij,[2,3,4,5],6,res_lo_eikqed)
!          res_lo_eikqed = transition('g -> q', 'none', res_lo_eikqed)
          call get_respdf_gen(1,1,C5S6Lim,res_lo_eikqed,respdf)
          respdf = respdf / E6sq
       endif

      z5  = C5S6Lim%Lim_z(1)
      s15 = C5S6Lim%Lim_sij(1,5)

      call partition_nnlo_fact(C5S6Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=1)

      respdf = respdf * (two/s15) * Tr * Pqg(z5) * C5S6Lim%wgt * damp

      FintNNLO_gq(8) = respdf(1)
      kin(6) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_gq)

#if(_withchecks == 1)
    FintNNLO_rr_tc_gq = FintNNLO_gq
    limval_nnlo = zero
    limval_nnlo(1) = FintNNLO_gq(1)    ! hard
    limval_nnlo(3) = FintNNLO_gq(2)    ! S6
    limval_nnlo(4) = FintNNLO_gq(3)    ! TC
    limval_nnlo(11) = FintNNLO_gq(4)    ! TC+C5
    limval_nnlo(9) = FintNNLO_gq(5)    ! TC+S6
    limval_nnlo(15) = FintNNLO_gq(6)    ! TC+C5+S6
    limval_nnlo(5) = FintNNLO_gq(7)    ! C5
    limval_nnlo(10) = FintNNLO_gq(8)    ! C5+S6

#endif

  end function xsect_nnlo_rr_5161c_gq

  function xsect_nnlo_rr_5161bd_gq(yRnd,ff,vegasweight,&
    sec)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5161bd_gq, sec
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: S6Lim, TCLim, TCS6Lim, C6Lim, C6S6Lim, TCC6Lim, TCC6S6Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_gq(8),kin(6)
    real(dp) :: res_nnlo(-5:7,-5:7), res_nlo(-5:7,-5:7), res_lo(-5:7,-5:7), res_nlo_eikqed(-5:7,-5:7), res_lo_ischarges(-5:7,-5:7)
    real(dp)  :: respdf(ipdf),respdf_tmp(ipdf,2),respdf_bak(ipdf)
    real(dp)  :: res_tmp(2,2)
    real(dp)  :: eik_qed(4),z1,z5,z6,s15,s16,s56,s156,E6,E6sq,damp
    logical   :: oldcode


    oldcode = .true.
    
    xsect_nnlo_rr_5161bd_gq = 0

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
         HardProc, S6Lim=S6Lim, TCLim=TCLim, TCC6Lim=TCC6Lim, &
         TCS6Lim=TCS6Lim, TCC6S6Lim=TCC6S6Lim, C6Lim=C6Lim, C6S6Lim=C6S6Lim, &
         opt_etas=[3,4])

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_q,id_a]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_gq(1) = zero

    else

       if (oldcode) then
          call res_tree_ga_gq(HardProc%AmpMom,res_tmp)
          call get_respdf(gq_lumi,1,1,HardProc,res_tmp,respdf)
       else
          call res_tree_ga_gq_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,HardProc,res_nnlo,respdf)
       endif
       call partition_nnlo_fact(HardProc,damp,iconf_qed=[5,2,3,4,6],i_qed=1)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_gq(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 6                                !!
    !!-----------------------------------------------------------------------!!

    S6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    S6Lim%npart = 5

    call cut_histo(S6Lim)

    if (S6Lim%makecut) then

      kin(2) = zero
      FintNNLO_gq(2) = zero

    else

      E6   = S6Lim%Lim_KinInv(1)
      E6sq = E6*E6

      if (oldcode) then
         call res_tree_g_gq(S6Lim%AmpMom,res_tmp)
         call get_qed_eik(charges,S6Lim%Lim_etaij,ql_pos,6,eik_qed)
         eik_qed = eik_qed/E6sq
         
         res_tmp(1:2,1) = res_tmp(1:2,1) * eik_qed(1:2)
         res_tmp(1:2,2) = res_tmp(1:2,2) * eik_qed(3:4)
         
         call get_respdf(gq_lumi,1,1,S6Lim,res_tmp,respdf)
      else
         call res_tree_g_gq_gen(S6Lim%AmpMom,res_nlo)
         call get_qed_eik_gen(res_nlo,S6Lim%Lim_etaij,[2,3,4,5],6,res_nlo_eikqed)
         call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf)
         respdf = respdf / E6sq
       endif

      call partition_nnlo_fact(S6Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=1)

      respdf =  respdf * S6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_gq(2) = respdf(1)
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
      FintNNLO_gq(3:4) = zero

    else

       if (oldcode) then
          call res_tree_qqb(TCLim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2
          
          call get_respdf(gq_lumi,1,1,TCLim,res_tmp,respdf_bak)
       else
          call res_tree_qqb_gen(TCLim%AmpMom,res_lo)
          res_lo = transition('g -> q', 'none', res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,1)
          call get_respdf_gen(1,1,TCLim,res_lo_ischarges,respdf_bak)
       endif
       
      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s15 = TCLim%Lim_sij(1,5)
      s16 = TCLim%Lim_sij(1,6)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * Tr * (-Pgaq(-s16,s56,-s15,z6,z1,z5)) * TCLim%wgt
      respdf_tmp(:,1) = -respdf_tmp(:,1)

      FintNNLO_gq(3) = respdf_tmp(1,1)

      !! ----------------------------- TC + C6 ----------------------------- !!

      !! (- si5 - si6 + s56) in the x4 -> 0 limit
      s156 = -TCC6Lim%Lim_KinInv(1)
      s56  =  TCC6Lim%Lim_sij(5,6)
      z5   =  TCC6Lim%Lim_z(1)
      z6   =  TCC6Lim%Lim_z(3)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (2/s56) * Tr * Pqq(z5)  &
                      * (2/s156) * Pqg(z6) &
                      * TCC6Lim%wgt

      FintNNLO_gq(4) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) +  respdf_tmp(:,2)
      kin(3) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 6                              !!
    !!-----------------------------------------------------------------------!!

    C6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C6Lim%npart = 5

    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(4) = zero
      FintNNLO_gq(5) = zero

    else

      z5  = C6Lim%Lim_z(1)
      s56 = C6Lim%Lim_sij(5,6)

      if (oldcode) then
         call res_tree_g_gq(C6Lim%AmpMom,res_tmp)
         res_tmp(:,1) = res_tmp(:,1) * Qdn2
         res_tmp(:,2) = res_tmp(:,2) * Qup2

         call get_respdf(gq_lumi,1,1,C6Lim,res_tmp,respdf)
      else
         call res_tree_g_gq_gen(C6Lim%AmpMom,res_lo)
         res_lo = multiply_FSQ_charges_sq(res_lo)
         call get_respdf_gen(1,1,C6Lim,res_lo,respdf)
      endif
         

      !-- damp = one
      respdf = (two/s56) * Pqq(z5) * respdf * C6Lim%wgt
      respdf = -respdf

      FintNNLO_gq(5) = respdf(1)
      kin(4) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                         Collinear 6 + Soft 6                          !!
    !!-----------------------------------------------------------------------!!

    C6S6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C6S6Lim%npart = 5

    call cut_histo(C6S6Lim)

    if (C6S6Lim%makecut) then

      kin(5) = zero
      FintNNLO_gq(6) = zero

    else

       if (oldcode) then
          call res_tree_g_gq(C6S6Lim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2
          call get_respdf(gq_lumi,1,1,C6S6Lim,res_tmp,respdf)
       else
          call res_tree_g_gq_gen(C6S6Lim%AmpMom,res_lo)
         res_lo = multiply_FSQ_charges_sq(res_lo)
         call get_respdf_gen(1,1,C6S6Lim,res_lo,respdf)
      endif

      z5  = C6S6Lim%Lim_z(1)
      s56 = C6S6Lim%Lim_sij(5,6)

      !-- damp = one
      respdf = respdf * (4/z5/s56) * C6S6Lim%wgt

      FintNNLO_gq(6) = respdf(1)
      kin(5) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                       Triple Collinear + Soft 6                       !!
    !!-----------------------------------------------------------------------!!

    TCS6Lim%ids = [0,0,id_el,-id_el,0,0]
    TCS6Lim%npart = 4

    call cut_histo(TCS6Lim)

    if (TCS6Lim%makecut) then

      kin(6) = zero
      FintNNLO_gq(7:8) = zero

    else

       if (oldcode) then
          call res_tree_qqb(TCS6Lim%AmpMom,res_tmp)
          res_tmp(:,1) = res_tmp(:,1) * Qdn2
          res_tmp(:,2) = res_tmp(:,2) * Qup2
          
          call get_respdf(gq_lumi,1,1,TCS6Lim,res_tmp,respdf_bak)
       else
          call res_tree_qqb(TCS6Lim%AmpMom,res_lo)
          res_lo = transition('g -> q', 'none', res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,1)
          call get_respdf_gen(1,1,TCS6Lim,res_lo,respdf_bak)
       endif

      z5  = TCS6Lim%Lim_z(1)
      z6  = TCS6Lim%Lim_z(2)
      z1  = TCS6Lim%Lim_z(3)
      s56 = TCS6Lim%Lim_sij(5,6)
      s15 = TCS6Lim%Lim_sij(1,5)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * (4*Tr*z5/z6/s56)     &
                      * (2/s15) * Pqg(z1) &
                      * TCS6Lim%wgt

      FintNNLO_gq(7) = respdf_tmp(1,1)

      !! --------------------------- TC + C6 + S6 -------------------------- !!

      z1 = TCC6S6Lim%Lim_z(3)
      s16 = TCC6S6Lim%Lim_sij(1,6)
      s56 = TCC6S6Lim%Lim_sij(5,6)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
                      * (4*Tr/s16)        &
                      * (2/s56) * Pqq(z1) &
                      * TCC6S6Lim%wgt
      respdf_tmp(:,2) = -respdf_tmp(:,2)

      FintNNLO_gq(8) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)
      kin(6) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_gq)

#if(_withchecks == 1)
    FintNNLO_rr_tc_gq = FintNNLO_gq
    limval_nnlo = zero
    limval_nnlo(1) = FintNNLO_gq(1)    ! hard
    limval_nnlo(3) = FintNNLO_gq(2)    ! S6
    limval_nnlo(4) = FintNNLO_gq(3)    ! TC
    limval_nnlo(11) = FintNNLO_gq(4)    ! TC+C5
    limval_nnlo(5) = FintNNLO_gq(5)    ! Coll
    limval_nnlo(10) = FintNNLO_gq(6)    ! C5+S6
    limval_nnlo(9) = FintNNLO_gq(7)    ! TC+S6
    limval_nnlo(15) = FintNNLO_gq(8)    ! TC+C6+S6

#endif

  end function xsect_nnlo_rr_5161bd_gq

  function xsect_nnlo_rr_5162_gq(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_dc_ii
    implicit none
    integer :: xsect_nnlo_rr_5162_gq
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: C5Lim, S6Lim, C6Lim, C6S6Lim, C5C6Lim, C5S6Lim, C5C6S6Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_gq(8),kin(1:6)
    real(dp)  :: respdf(ipdf),respdf_tmp(ipdf,2)
    real(dp) :: res_nnlo(-5:7,-5:7), res_nlo(-5:7,-5:7), res_lo(-5:7,-5:7), res_nlo_eikqed(-5:7,-5:7),  res_lo_eikqed(-5:7,-5:7), res_lo_ischarges(-5:7,-5:7), res_nlo_ischarges(-5:7,-5:7)
    real(dp)  :: res_nlo_old(2,2),res_tmp(2,2)
    real(dp)  :: eik_qed(4)
    real(dp)  :: z5,z6,s15,s26,E6,E6sq,damp
    logical   :: oldcode

    oldcode = .false.

    xsect_nnlo_rr_5162_gq = 0

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
         xx(1:kNNLO_max_full), 1, 2, &
         HardProc, C5Lim=C5Lim, S6Lim=S6Lim, C6Lim=C6Lim, & 
         C6S6Lim=C6S6Lim, C5C6Lim=C5C6Lim, C5S6Lim=C5S6Lim, C5C6S6Lim=C5C6S6Lim, &
         opt_etas=[3,4])

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_q,id_a]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_gq(1) = zero

    else

       if (oldcode) then
          call res_tree_ga_gq(HardProc%AmpMom,res_tmp)
          call get_respdf(gq_lumi,1,1,HardProc,res_tmp,respdf)
       else
          call res_tree_ga_gq_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,HardProc,res_nnlo,respdf)
       endif

      call partition_nnlo_fact(HardProc,damp,iconf_qed=[5,2,3,4,6],i_qed=2)

      respdf = respdf * HardProc%wgt * damp

      kin(1) = respdf(1)
      FintNNLO_gq(1) = kin(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 6                                !!
    !!-----------------------------------------------------------------------!!

    S6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    S6Lim%npart = 5

    call cut_histo(S6Lim)

    if (S6Lim%makecut) then

      kin(2) = zero
      FintNNLO_gq(2:3) = zero

    else



    !! --------------------------------- S6 --------------------------------- !!
      E6   = S6Lim%Lim_KinInv(1)
      E6sq = E6*E6
      if (oldcode) then
         call res_tree_g_gq(S6Lim%AmpMom,res_tmp)
         call get_qed_eik(charges,S6Lim%Lim_etaij,ql_pos,6,eik_qed)
         eik_qed = eik_qed/E6sq         
         res_nlo_old(:,1) = res_tmp(:,1) * eik_qed(1:2)
         res_nlo_old(:,2) = res_tmp(:,2) * eik_qed(3:4)         
         call get_respdf(gq_lumi,1,1,S6Lim,res_nlo_old,respdf_tmp(:,1))
      else
         call res_tree_g_gq_gen(S6Lim%AmpMom,res_nlo)
         call get_qed_eik_gen(res_nlo,S6Lim%Lim_etaij,[2,3,4,5],6,res_nlo_eikqed)
         call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf_tmp(:,1))
         respdf_tmp(:,1) = respdf_tmp(:,1) / E6sq
      endif
      
      call partition_nnlo_fact(S6Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=2)

      respdf_tmp(:,1) =  respdf_tmp(:,1) * S6Lim%wgt * damp
      respdf_tmp(:,1) = -respdf_tmp(:,1)
      FintNNLO_gq(2)  =  respdf_tmp(1,1)

    !! ------------------------------ C6 + S6 ------------------------------ !!
      C6S6Lim%mur = S6Lim%mur
      C6S6Lim%muf = S6Lim%muf

      z6  = C6S6Lim%Lim_z(2)
      s26 = C6S6Lim%Lim_sij(2,6)


      if (oldcode) then
         res_nlo_old(:,1) = res_tmp(:,1) * Qdn2
         res_nlo_old(:,2) = res_tmp(:,2) * Qup2
         
         call get_respdf(gq_lumi,1,1,C6S6Lim,res_nlo_old,respdf_tmp(:,2))
      else
         res_nlo_ischarges = multiply_IS_charges_sq(res_nlo,2)
         call get_respdf_gen(1,1,C6S6Lim,res_nlo_ischarges,respdf_tmp(:,2))
      endif

      call partition_nnlo_fact(C6S6Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=2)

      respdf_tmp(:,2) = respdf_tmp(:,2) * (4/z6/s26) * C6S6Lim%wgt * damp
      FintNNLO_gq(3)  = respdf_tmp(1,2)

      respdf = respdf_tmp(1,1) + respdf_tmp(1,2)
      kin(2) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 6                              !!
    !!-----------------------------------------------------------------------!!

    C6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C6Lim%npart = 5

    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(3) = zero
      FintNNLO_gq(4) = zero

    else

      z6  = C6Lim%Lim_z(2)
      s26 = C6Lim%Lim_sij(2,6)

      if (oldcode) then
         call res_tree_g_gq(C6Lim%AmpMom,res_tmp)
         res_tmp(:,1) = res_tmp(:,1) * Qdn2
         res_tmp(:,2) = res_tmp(:,2) * Qup2

         call get_respdf(gq_lumi,1,1,C6Lim,res_tmp,respdf)
      else
         call res_tree_g_gq_gen(C6Lim%AmpMom,res_nlo)
         res_nlo_ischarges = multiply_IS_charges_sq(res_nlo,2)
         call get_respdf_gen(1,1,C6Lim,res_nlo_ischarges,respdf)
      endif

      call partition_nnlo_fact(C6Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=2)

      respdf = (-two/s26) * Pqq(z6) * respdf * C6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_gq(4) = respdf(1)

      kin(3) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                         Collinear 5 + Soft 6                          !!
    !!-----------------------------------------------------------------------!!

    C5S6Lim%ids = [0,0,id_el,-id_el,0,0]
    C5S6Lim%npart = 4

    call cut_histo(C5S6Lim)

    if (C5S6Lim%makecut) then

      kin(4) = zero
      FintNNLO_gq(5:6) = zero

    else



    !! ------------------------------ C5 + S6 ------------------------------ !!
      E6   = C5S6Lim%Lim_KinInv(1)
      E6sq = E6*E6
      if (oldcode) then
         call res_tree_qqb(C5S6Lim%AmpMom,res_tmp)
         call get_qed_eik(charges,C5S6Lim%Lim_etaij,ql_pos,6,eik_qed)
         eik_qed = eik_qed/E6sq
         res_nlo_old(:,1) = res_tmp(:,1) * eik_qed(1:2)
         res_nlo_old(:,2) = res_tmp(:,2) * eik_qed(3:4)
         call get_respdf(gq_lumi,1,1,C5S6Lim,res_nlo_old,respdf_tmp(:,1))
      else
         call res_tree_qqb_gen(C5S6Lim%AmpMom,res_lo)
         res_lo = transition('g -> q', 'none', res_lo)
         call get_qed_eik_gen(res_lo,C5S6Lim%Lim_etaij,[2,3,4,5],6,res_lo_eikqed)
         call get_respdf_gen(1,1,C5S6Lim,res_lo_eikqed,respdf_tmp(:,1))
         respdf_tmp(:,1)  = respdf_tmp(:,1)/E6sq
      endif
         
         
         

      z5  = C5S6Lim%Lim_z(1)
      s15 = C5S6Lim%Lim_sij(1,5)



      call partition_nnlo_fact(C5S6Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=2)

      respdf_tmp(:,1) = respdf_tmp(:,1) &
                      * (two/s15) * Tr * Pqg(z5) * C5S6Lim%wgt * damp

      FintNNLO_gq(5) = respdf_tmp(1,1)

    !! --------------------------- C5 + C6 + S6 ---------------------------- !!
      C5C6S6Lim%mur = C5S6Lim%mur
      C5C6S6Lim%muf = C5S6Lim%muf

      z5  = C5C6S6Lim%Lim_z(1)
      s15 = C5C6S6Lim%Lim_sij(1,5)
      z6  = C5C6S6Lim%Lim_z(2)
      s26 = C5C6S6Lim%Lim_sij(2,6)

      if (oldcode) then
         res_nlo_old(:,1) = res_tmp(:,1) * Qdn2
         res_nlo_old(:,2) = res_tmp(:,2) * Qup2
         
         call get_respdf(gq_lumi,1,1,C5C6S6Lim,res_nlo_old,respdf_tmp(:,2))
      else
         res_lo_ischarges = multiply_IS_charges_sq(res_lo,2)
         call get_respdf_gen(1,1,C5C6S6Lim,res_lo_ischarges,respdf_tmp(:,2))
      endif
         

      !-- damp = one
      respdf_tmp(:,2) = respdf_tmp(:,2) &
                      * (4/z6/s26)      &
                      * (two/s15) * Tr * Pqg(z5) &
                      * C5C6S6Lim%wgt
      respdf_tmp(:,2) = - respdf_tmp(:,2)

      FintNNLO_gq(6) = respdf_tmp(1,2)

      respdf = respdf_tmp(1,1) + respdf_tmp(1,2)
      kin(4) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 5                              !!
    !!-----------------------------------------------------------------------!!

    C5Lim%ids = [0,0,id_el,-id_el,id_a,0]
    C5Lim%npart = 5

    call cut_histo(C5Lim)

    if (C5Lim%makecut) then

      kin(5) = zero
      FintNNLO_gq(7) = zero

    else

      z5  = C5Lim%Lim_z(1)
      s15 = C5Lim%Lim_sij(1,5)

      if (oldcode) then
         call res_tree_a_qqb(C5Lim%AmpMom,res_tmp)
         call get_respdf(gq_lumi,1,1,C5Lim,res_tmp,respdf)
      else
         call res_tree_a_qqb_gen(C5Lim%AmpMom,res_nlo)
         res_nlo = transition('g -> q', 'none', res_nlo)
         call get_respdf_gen(1,1,C5Lim,res_nlo,respdf)
      endif
      
      call partition_nnlo_fact(C5Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=2)

      respdf = (two/s15) * Tr * Pqg(z5) * respdf * C5Lim%wgt * damp
      respdf = -respdf

      FintNNLO_gq(7) = respdf(1)

      kin(5) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                       Collinear 5 + Collinear 6                       !!
    !!-----------------------------------------------------------------------!!

    C5C6Lim%ids = [0,0,id_el,-id_el,0,0]
    C5C6Lim%npart = 4

    call cut_histo(C5C6Lim)

    if (C5C6Lim%makecut.or.C5C6Lim%flag) then

      kin(6) = zero
      FintNNLO_gq(8) = zero

    else

      z5  = C5C6Lim%Lim_z(1)
      s15 = C5C6Lim%Lim_sij(1,5)
      z6  = C5C6Lim%Lim_z(2)
      s26 = C5C6Lim%Lim_sij(2,6)

      if (oldcode) then
         call res_tree_qqb(C5C6Lim%AmpMom,res_tmp)
         res_tmp(:,1) = res_tmp(:,1) * Qdn2
         res_tmp(:,2) = res_tmp(:,2) * Qup2

         call get_respdf(gq_lumi,1,1,C5C6Lim,res_tmp,respdf)
      else
         call res_tree_qqb_gen(C5C6Lim%AmpMom,res_lo)
         res_lo = transition('g -> q', 'none', res_lo)
         res_lo_ischarges = multiply_IS_charges_sq(res_lo,2)
         call get_respdf_gen(1,1,C5C6Lim,res_lo_ischarges,respdf)
      endif
         

      !-- damp = one
      respdf = respdf &
             * ( two/s15) * Tr * Pqg(z5)  &
             * (-two/s26) * Pqq(z6) &
             * C5C6Lim%wgt

      FintNNLO_gq(8) = respdf(1)

      kin(6) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_gq)

#if(_withchecks == 1)
    FintNNLO_rr_dc_gq = FintNNLO_gq
    limval_nnlo = zero
    limval_nnlo(1) = FintNNLO_gq(1)    ! hard
    limval_nnlo(3) = FintNNLO_gq(2)    ! S6
    limval_nnlo(10) = FintNNLO_gq(3)    ! C6+S6
    limval_nnlo(5) = FintNNLO_gq(4)    ! C6
    limval_nnlo(9) = FintNNLO_gq(5)    ! C5+S6
    limval_nnlo(15) = FintNNLO_gq(6)    ! C5+C6+S6
    limval_nnlo(4) = FintNNLO_gq(7)    ! C5
    limval_nnlo(11) = FintNNLO_gq(8)    ! C5 + C6
#endif

  end function xsect_nnlo_rr_5162_gq

  function xsect_nnlo_rr_516k_gq(yRnd,ff,vegasweight,k,l)
    use mod_kinematics_nnlo_dc_if
    implicit none
    integer :: xsect_nnlo_rr_516k_gq, k, l
    real(dp15) :: yRnd(30), ff(1), vegasweight
    !--
    integer, parameter :: imax_ipdf = 2, imax_ilim = 2
    type(KinConfig) :: HardProc
    type(KinConfig) :: C6Lim, S6Lim, C6S6Lim, C5Lim, C5C6Lim, C5S6Lim, C5C6S6Lim
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_gq(8),kin(1:6)
    real(dp)  :: respdf(ipdf),respdf_tmp(imax_ipdf,ipdf)
    real(dp)  :: respdf_vect(imax_ilim,ipdf)
    real(dp) :: res_nnlo(-5:7,-5:7), res_nlo(-5:7,-5:7), res_lo(-5:7,-5:7), res_nlo_eikqed(-5:7,-5:7),  res_lo_eikqed(-5:7,-5:7), res_lo_ischarges(-5:7,-5:7), res_nlo_ischarges(-5:7,-5:7)
    real(dp)  :: res_tmp(2,2),res_tmp_vect(2,2,imax_ipdf)
    real(dp)  :: eik_qed(4), Qsq_FS(2)
    real(dp)  :: z5,z6,s15,sk6,E6,E6sq,damp
    logical   :: oldcode


    oldcode = .false.

    xsect_nnlo_rr_516k_gq = 0

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
    xx(1:kNNLO_max_full), 1, 2, k, l, &
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
      FintNNLO_gq(1) = zero

    else

       if (oldcode) then
          call res_tree_ga_gq(HardProc%AmpMom,res_tmp)
          call get_respdf(gq_lumi,1,1,HardProc,res_tmp,respdf)
       else
          call res_tree_ga_gq_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,HardProc,res_nnlo,respdf)
       endif
      call partition_nnlo_fact(HardProc,damp,iconf_qed=[5,2,3,4,6],i_qed=k)

      respdf = respdf * HardProc%wgt * damp

      kin(1) = respdf(1)
      FintNNLO_gq(1) = kin(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                 Soft 6                                !!
    !!-----------------------------------------------------------------------!!
    S6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    call cut_histo(S6Lim)

    if (S6Lim%makecut) then

      kin(2) = zero
      FintNNLO_gq(2:3) = zero

    else

       E6   = S6Lim%Lim_KinInv(1)
       E6sq = E6*E6

       if (oldcode) then
          call res_tree_g_gq(S6Lim%AmpMom,res_tmp)

      !-- prepare PDFs structures, S6
          call get_qed_eik(charges,S6Lim%Lim_etaij,ql_pos,6,eik_qed)
          eik_qed = eik_qed/E6sq
          res_tmp_vect(:,1,1) = res_tmp(:,1) * eik_qed(1:2)
          res_tmp_vect(:,2,1) = res_tmp(:,2) * eik_qed(3:4)
          
          !-- prepare PDFs structures, C6S6
          res_tmp_vect(:,:,2) = res_tmp(:,:) * Q_lep2
          
          call get_respdf_vect(gq_lumi,1,1,S6Lim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
       else
          call res_tree_g_gq_gen(S6Lim%AmpMom,res_nlo)
          call get_qed_eik_gen(res_nlo,S6Lim%Lim_etaij,[2,3,4,5],6,res_nlo_eikqed)
          call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf_vect(1,:))
          respdf_vect(1,:) = respdf_vect(1,:) / E6sq                    
          res_nlo = Qsq_Fs(k-2)**2 * res_nlo
          call get_respdf_gen(1,1,S6Lim,res_nlo,respdf_vect(2,:))
       endif
          
    !! --------------------------------- S6 -------------------------------- !!

      call partition_nnlo_fact(S6Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=k)

      respdf_vect_part(1,:) = - respdf_vect(1,:) * S6Lim%wgt * damp
      FintNNLO_gq(2) = respdf_vect_part(1,1)

    !! ------------------------------ C6 + S6 ------------------------------ !!
      z6  = C6S6Lim%Lim_z(2)
      sk6 = C6S6Lim%Lim_sij(k,6)

      call partition_nnlo_fact(C6S6Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=k)

      if (oldcode) then
         respdf_vect_part(2,:) = respdf_vect(2,:) * (4*Q_lep2/z6/sk6) * C6S6Lim%wgt * damp
      else
         respdf_vect_part(2,:) = respdf_vect(2,:) * (four/z6/sk6) * C6S6Lim%wgt * damp
      endif
      
      FintNNLO_gq(3)  = respdf_vect_part(2,1)

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
      FintNNLO_gq(4) = zero

    else

      z6  = C6Lim%Lim_z(2)
      sk6 = C6Lim%Lim_sij(k,6)

      if (oldcode) then
         call res_tree_g_gq(C6Lim%AmpMom,res_tmp)
         call get_respdf(gq_lumi,1,1,C6Lim,res_tmp,respdf)
      else
         call res_tree_g_gq_gen(C6Lim%AmpMom,res_nlo)
         res_nlo = Qsq_Fs(k-2)**2 * res_nlo
         call get_respdf_gen(1,1,C6Lim,res_nlo,respdf)
      endif
      call partition_nnlo_fact(C6Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=k)

      if (oldcode) then
         respdf = respdf * 2/sk6 * Q_lep2 * Pqg(z6) * C6Lim%wgt * damp
      else
         respdf = respdf * 2/sk6 * Pqg(z6) * C6Lim%wgt * damp
      endif

      respdf = -respdf

      FintNNLO_gq(4) = respdf(1)

      kin(3) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 5                              !!
    !!-----------------------------------------------------------------------!!

    C5Lim%ids = [0,0,id_el,-id_el,id_a,0]
    C5Lim%npart = 5
    
    call cut_histo(C5Lim)

    if (C5Lim%makecut) then

      kin(4) = zero
      FintNNLO_gq(5) = zero

    else

      z5  = C5Lim%Lim_z(1)
      s15 = C5Lim%Lim_sij(1,5)
      
      if (oldcode) then
         call res_tree_a_qqb(C5Lim%AmpMom,res_tmp)
         call get_respdf(gq_lumi,1,1,C5Lim,res_tmp,respdf)
      else
         call res_tree_a_qqb_gen(C5Lim%AmpMom,res_nlo)
         res_nlo = transition('g -> q', 'none', res_nlo)
         call get_respdf_gen(1,1,C5Lim,res_nlo,respdf)
      endif

      call partition_nnlo_fact(C5Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=k)

      respdf = respdf * (2/s15) * Tr * Pqq(z5) * C5Lim%wgt * damp
      respdf = -respdf

      FintNNLO_gq(5) = respdf(1)

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
      FintNNLO_gq(6:7) = zero

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

          call get_respdf_vect(gq_lumi,1,1,C5S6Lim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
       else
          call res_tree_qqb_gen(C5S6Lim%AmpMom,res_lo)
          res_lo = transition('g -> q', 'none', res_lo)
          call get_qed_eik_gen(res_lo,C5S6Lim%Lim_etaij,[2,3,4,5],6,res_lo_eikqed)
          call get_respdf_gen(1,1,C5S6Lim,res_lo_eikqed,respdf_vect(1,:))
          respdf_vect(1,:) = respdf_vect(1,:) / E6sq
          res_lo = Qsq_Fs(k-2)**2 * res_lo
          call get_respdf_gen(1,1,S6Lim,res_lo,respdf_vect(2,:))
       endif

    !! ------------------------------ C5 + S6 ------------------------------ !!

      z5  = C5S6Lim%Lim_z(1)
      s15 = C5S6Lim%Lim_sij(1,5)

      call partition_nnlo_fact(C5S6Lim,damp,iconf_qed=[5,2,3,4,6],i_qed=k)

      respdf_vect_part(1,:) = respdf_vect(1,:) &
                            * (2/s15) * Tr * Pqq(z5) * C5S6Lim%wgt * damp

      FintNNLO_gq(6) = respdf_vect_part(1,1)

    !! --------------------------- C5 + C6 + S6 ---------------------------- !!

      z5  = C5C6S6Lim%Lim_z(1)
      s15 = C5C6S6Lim%Lim_sij(1,5)
      z6  = C5C6S6Lim%Lim_z(2)
      sk6 = C5C6S6Lim%Lim_sij(k,6)

      !-- damp = one
      if (oldcode) then
         respdf_vect_part(2,:) = - respdf_vect(2,:) &
              * (4*Q_lep2/z6/sk6) &
              * (two/s15) * Tr * Pqq(z5) &
              * C5C6S6Lim%wgt
      else
         respdf_vect_part(2,:) = - respdf_vect(2,:) &
              * (four/z6/sk6) &
              * (two/s15) * Tr * Pqq(z5) &
              * C5C6S6Lim%wgt
      endif
      
      FintNNLO_gq(7) = respdf_vect_part(2,1)

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
      FintNNLO_gq(8) = zero

    else

      z5  = C5C6Lim%Lim_z(1)
      s15 = C5C6Lim%Lim_sij(1,5)
      z6  = C5C6Lim%Lim_z(2)
      sk6 = C5C6Lim%Lim_sij(k,6)

      if (oldcode) then
         call res_tree_qqb(C5C6Lim%AmpMom,res_tmp)
         call get_respdf(gq_lumi,1,1,C5C6Lim,res_tmp,respdf)
      else
         call res_tree_qqb_gen(C5C6Lim%AmpMom,res_lo)
         res_lo = transition('g -> q', 'none', res_lo)
         res_lo = Qsq_Fs(k-2)**2 * res_lo
         call get_respdf_gen(1,1,C5C6Lim,res_lo,respdf)
      endif


      !-- damp = one
      if (oldcode) then
         respdf = respdf &
              * (2/s15) * Tr * Pqq(z5)  &
              * (2/sk6) * Q_lep2 * Pqg(z6) &
              * C5C6Lim%wgt
      else
         respdf = respdf &
              * (2/s15) * Tr * Pqq(z5)  &
              * (2/sk6) *  Pqg(z6) &
              * C5C6Lim%wgt
      endif
      
      FintNNLO_gq(8) = respdf(1)

      kin(6) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_gq)

#if(_withchecks == 1)
    FintNNLO_rr_dc_gq = FintNNLO_gq
    limval_nnlo = zero
    limval_nnlo(1) = FintNNLO_gq(1)    ! hard
    limval_nnlo(3) = FintNNLO_gq(2)    ! S6
    limval_nnlo(10) = FintNNLO_gq(3)    ! C6+S6
    limval_nnlo(5) = FintNNLO_gq(4)    ! C6
    limval_nnlo(4) = FintNNLO_gq(5)    ! C5
    limval_nnlo(9) = FintNNLO_gq(6)    ! C5+S6
    limval_nnlo(15) = FintNNLO_gq(7)    ! C5+C6+S6
    limval_nnlo(11) = FintNNLO_gq(8)    ! C5 + C6
#endif

  end function xsect_nnlo_rr_516k_gq

end module mod_xsects_nnlo_rr_gq

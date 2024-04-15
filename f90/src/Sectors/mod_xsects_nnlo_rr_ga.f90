module mod_xsects_nnlo_rr_ga
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
  real(dp), public, save :: FintNNLO_rr_tc_ga(4), FintNNLO_rr_dc_ga(4)
#endif

  public :: xsect_nnlo_rr_5161a_ga,xsect_nnlo_rr_5161c_ga
  public :: xsect_nnlo_rr_5161b_ga,xsect_nnlo_rr_5161d_ga
  public :: xsect_nnlo_rr_5162_ga,xsect_nnlo_rr_5261_ga
  public :: xsect_nnlo_rr_5262_ga

  private

contains

  function xsect_nnlo_rr_5161b_ga(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161b_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5161b_ga = xsect_nnlo_rr_5161bd_ga(yRnd,ff,vegasweight,2)

  end function xsect_nnlo_rr_5161b_ga

  function xsect_nnlo_rr_5161d_ga(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161d_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5161d_ga = xsect_nnlo_rr_5161bd_ga(yRnd,ff,vegasweight,4)

  end function xsect_nnlo_rr_5161d_ga

  !!*************************************************************************!!

  function xsect_nnlo_rr_5161a_ga(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5161a_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc, C6Lim, TCLim, TCC6Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_ga(4),kin(3)
    real(dp)  :: respdf(ipdf),respdf_bak(ipdf),respdf_tmp(ipdf,2)
    real(dp)  :: ampres(1,1),res_tmp(2,2)
    real(dp)  :: z1,z5,z6,s15,s16,s56,s156,damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5161a_ga = 0

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
         xx(1:kNNLO_max_full),1,2,1, &
         HardProc, TCLim=TCLim, TCC6Lim=TCC6Lim, C6Lim=C6Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_q,id_q]
    HardProc%npart = 6
    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_ga(1) = zero

    else

      call res_tree_ga_ga(HardProc%AmpMom,ampres(1,1))
      call get_respdf(ga_lumi,1,1,HardProc,ampres,respdf)

      call partition_nnlo_qcd(HardProc, damp, iconf, 1,1)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_ga(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Triple Collinear                             !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids(1:6) = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4
    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(2) = zero
      FintNNLO_ga(2:3) = zero

    else

      call res_treeAA_aa(TCLim%AmpMom,ampres(1,1))
      ampres(1,1) = ampres(1,1) * (ndn*Qdn2 + nup*Qup2)

      call get_respdf(ga_lumi,1,1,TCLim,ampres,respdf_bak)

    !! --------------------------------- TC -------------------------------- !!

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s15 = TCLim%Lim_sij(1,5)
      s16 = TCLim%Lim_sij(1,6)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * Tr * Pgqqb_spav_ab(-s15,-s16,s56,z1,z5,z6) &
                      * TCLim%wgt
      respdf_tmp(:,1) = -respdf_tmp(:,1)

      FintNNLO_ga(2) = respdf_tmp(1,1)

    !! ------------------------------- TC + C6 ----------------------------- !!

      z5   = TCC6Lim%Lim_z(1)
      z6   = TCC6Lim%Lim_z(2)
      s156 = TCC6Lim%Lim_KinInv(1)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
              * (2/s156 * Pgq_spav(z5)) &
              * (2/s16  * Tr * Pqg(one-z6))  &
              * TCC6Lim%wgt

      FintNNLO_ga(3) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)

      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 6                              !!
    !!-----------------------------------------------------------------------!!

    C6Lim%ids(1:6) = [0,0,id_el,-id_el,id_q,0]
    C6Lim%npart = 5
    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(3) = zero
      FintNNLO_ga(4) = zero

    else

      z6  = C6Lim%Lim_z(2)
      s16 = C6Lim%Lim_sij(1,6)

      call res_tree_a_qa(C6Lim%AmpMom,res_tmp)

      !! take only quarks, i.e. res_tmp(1,:)
      ampres(1,1) = res_tmp(1,1)*ndn + res_tmp(1,2)*nup
      call get_respdf(ga_lumi,1,1,C6Lim,ampres,respdf)

      call partition_nnlo_qcd(C6Lim, damp, iconf, 1,1)

      respdf = (two/s16) * Tr * Pqg(z6) * respdf * C6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_ga(4) = respdf(1)
      kin(3) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ga)

#if(_withchecks == 1)
    FintNNLO_rr_tc_ga = FintNNLO_ga
#endif

  end function xsect_nnlo_rr_5161a_ga

  function xsect_nnlo_rr_5161c_ga(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5161c_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc, C5Lim, TCLim, TCC5Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_ga(4),kin(3)
    real(dp)  :: respdf(ipdf),respdf_bak(ipdf),respdf_tmp(ipdf,2)
    real(dp)  :: ampres(1,1),res_tmp(2,2)
    real(dp)  :: z1,z5,z6,s15,s16,s56,s156,damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5161c_ga = 0

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
         xx(1:kNNLO_max_full),1,2,3, &
         HardProc, TCLim=TCLim, TCC5Lim=TCC5Lim, C5Lim=C5Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_q,id_q]
    HardProc%npart = 6
    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_ga(1) = zero

    else

      call res_tree_ga_ga(HardProc%AmpMom,ampres(1,1))
      call get_respdf(ga_lumi,1,1,HardProc,ampres,respdf)

      call partition_nnlo_qcd(HardProc, damp, iconf, 1,1)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_ga(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Triple Collinear                             !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids(1:6) = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4
    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(2) = zero
      FintNNLO_ga(2:3) = zero

    else

      call res_treeAA_aa(TCLim%AmpMom,ampres(1,1))
      ampres(1,1) = ampres(1,1) * (ndn*Qdn2 + nup*Qup2)

      call get_respdf(ga_lumi,1,1,TCLim,ampres,respdf_bak)

    !! --------------------------------- TC -------------------------------- !!

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s15 = TCLim%Lim_sij(1,5)
      s16 = TCLim%Lim_sij(1,6)

      !-- damp = one
      respdf_tmp(:,1) = respdf_bak &
                      * Tr * Pgqqb_spav_ab(-s15,-s16,s56,z1,z5,z6) &
                      * TCLim%wgt
      respdf_tmp(:,1) = -respdf_tmp(:,1)

      FintNNLO_ga(2) = respdf_tmp(1,1)

    !! ------------------------------- TC + C5 ----------------------------- !!

      z5   = TCC5Lim%Lim_z(1)
      z6   = TCC5Lim%Lim_z(2)
      s156 = TCC5Lim%Lim_KinInv(1)

      !-- damp = one
      respdf_tmp(:,2) = respdf_bak &
              * (2/s156 * Pgq_spav(z6)) &
              * (2/s15  * Tr * Pqq(z5)) & !-- Pqg(1-z5)
              * TCC5Lim%wgt

      FintNNLO_ga(3) = respdf_tmp(1,2)

      respdf = respdf_tmp(:,1) + respdf_tmp(:,2)

      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 5                              !!
    !!-----------------------------------------------------------------------!!

    C5Lim%ids(1:6) = [0,0,id_el,-id_el,id_q,0]
    C5Lim%npart = 5
    call cut_histo(C5Lim)

    if (C5Lim%makecut) then

      kin(3) = zero
      FintNNLO_ga(4) = zero

    else

      z5  = C5Lim%Lim_z(1)
      s15 = C5Lim%Lim_sij(1,5)

      call res_tree_a_qa(C5Lim%AmpMom,res_tmp)

      !-- take only anti-quarks, i.e. res_tmp(2,:)
      ampres(1,1) = res_tmp(2,1)*ndn + res_tmp(2,2)*nup
      call get_respdf(ga_lumi,1,1,C5Lim,ampres,respdf)

      call partition_nnlo_qcd(C5Lim, damp, iconf, 1,1)

      respdf = (two/s15) * Tr * Pqg(z5) * respdf * C5Lim%wgt * damp
      respdf = -respdf

      FintNNLO_ga(4) = respdf(1)
      kin(3) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ga)

#if(_withchecks == 1)
    FintNNLO_rr_tc_ga = FintNNLO_ga
#endif

  end function xsect_nnlo_rr_5161c_ga


  function xsect_nnlo_rr_5161bd_ga(yRnd,ff,vegasweight,sec)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5161bd_ga
    integer, intent(in) :: sec
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc, TCLim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_ga(2),kin(2)
    real(dp)  :: ampres(1,1), respdf(ipdf)
    real(dp)  :: z1,z5,z6,s15,s16,s56,damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5161bd_ga = 0

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
         xx(1:kNNLO_max_full),1,2,sec,HardProc,TCLim=TCLim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_q,id_q]
    HardProc%npart = 6
    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_ga(1) = zero

    else

      call res_tree_ga_ga(HardProc%AmpMom,ampres(1,1))
      call get_respdf(ga_lumi,1,1,HardProc,ampres,respdf)

      call partition_nnlo_qcd(HardProc, damp, iconf, 1,1)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_ga(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Triple Collinear                             !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids(1:6) = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4
    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(2) = zero
      FintNNLO_ga(2) = zero

    else

      call res_treeAA_aa(TCLim%AmpMom,ampres(1,1))
      ampres(1,1) = ampres(1,1) * (ndn*Qdn2 + nup*Qup2)

      call get_respdf(ga_lumi,1,1,TCLim,ampres,respdf)

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s15 = TCLim%Lim_sij(1,5)
      s16 = TCLim%Lim_sij(1,6)

      !-- damp = one
      respdf = respdf &
             * Tr * Pgqqb_spav_ab(-s15,-s16,s56,z1,z5,z6) &
             * TCLim%wgt
      respdf = -respdf

      FintNNLO_ga(2) = respdf(1)
      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ga)

#if(_withchecks == 1)
    FintNNLO_rr_tc_ga(1:2) = FintNNLO_ga
    FintNNLO_rr_tc_ga(3:4) = zero
#endif

  end function xsect_nnlo_rr_5161bd_ga

  function xsect_nnlo_rr_5262_ga(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_dc_ii
    implicit none
    integer :: xsect_nnlo_rr_5262_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: C5Lim,C6Lim
    !--
    real(dp) :: xx(kNNLO_max_full)
    real(dp) :: FintNNLO_ga(3),kin(3)
    real(dp) :: respdf(ipdf)
    real(dp) :: res_tmp1(1,1),res_tmp2(2,2)
    real(dp) :: z5,z6,s25,s26
    real(dp) :: damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5262_ga = 0

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
    !! Use this even if no soft singularities
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_ii_5i6i(xx(1:kNNLO_max_full), 2, 1, &
         HardProc=HardProc,C5Lim=C5Lim,C6Lim=C6Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_q,-id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

       kin(1) = zero
       FintNNLO_ga(1) = zero
       
    else
       
       call res_tree_ga_ga(HardProc%AmpMom,res_tmp1(1,1))
       call get_respdf(ga_lumi,1,1,HardProc,res_tmp1,respdf)
       
       call partition_nnlo_qcd(HardProc,damp,iconf,2,2)
       
       respdf = respdf * HardProc%wgt * damp
       
       FintNNLO_ga(1) = respdf(1)
       kin(1) = respdf(1)
       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                               C5                                      !!
    !!-----------------------------------------------------------------------!!
    
    C5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C5Lim%npart = 5
    
    call cut_histo(C5Lim)
    
    if (C5Lim%makecut) then
       
       kin(2) = zero
       FintNNLO_ga(2) = zero
    
    else
       
       call res_tree_g_gq(C5Lim%AmpMom,res_tmp2)
       res_tmp1(1,1) = res_tmp2(1,1)*ndn*Qdn2 + res_tmp2(1,2)*nup*Qup2
       call get_respdf(ga_lumi,1,1,C5Lim,res_tmp1,respdf)
       
       call partition_nnlo_qcd(C5Lim,damp,iconf,2,2)
       
       s25 = C5Lim%Lim_sij(2,5)
       z5  = C5Lim%Lim_z(1)

       respdf = -respdf*xn*(two/s25*Pqg(z5))*C5Lim%wgt*damp
       
       FintNNLO_ga(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    
    !!-----------------------------------------------------------------------!!
    !!                               C6                                      !!
    !!-----------------------------------------------------------------------!!
    
    C6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C6Lim%npart = 5
    
    call cut_histo(C6Lim)
    
    if (C6Lim%makecut) then

       kin(3) = zero
       FintNNLO_ga(3) = zero
       
    else

       call res_tree_g_gq(C6Lim%AmpMom,res_tmp2)
       res_tmp1(1,1) = res_tmp2(2,1)*ndn*Qdn2 + res_tmp2(2,2)*nup*Qup2
       call get_respdf(ga_lumi,1,1,C6Lim,res_tmp1,respdf)
       
       call partition_nnlo_qcd(C6Lim,damp,iconf,2,2)
       
       s26 = C6Lim%Lim_sij(2,6)
       z6  = C6Lim%Lim_z(2)
       
       respdf = -respdf*xn*(two/s26*Pqg(z6))*C6Lim%wgt*damp
       
       FintNNLO_ga(3) = respdf(1)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif
    
    ff(1) = sum(kin)
    call close_histo()
    
    call check_ff(ff,xx,FintNNLO_ga)

#if(_withchecks == 1)
    FintNNLO_rr_tc_ga(1:3) = FintNNLO_ga
#endif

  end function xsect_nnlo_rr_5262_ga

  function xsect_nnlo_rr_5162_ga(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_dc_ii
    implicit none
    integer :: xsect_nnlo_rr_5162_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc, C5Lim, C6Lim, C5C6Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_ga(4),kin(4)
    real(dp)  :: respdf(ipdf)
    real(dp)  :: ampres(1,1),res_tmp(2,2)
    real(dp)  :: z5,z6,s15,s26,damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5162_ga = 0

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
         HardProc, C5Lim=C5Lim, C6Lim=C6Lim, C5C6Lim=C5C6Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_q,id_q]
    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_ga(1) = zero

    else

      call res_tree_ga_ga(HardProc%AmpMom,ampres(1,1))
      call get_respdf(ga_lumi,1,1,HardProc,ampres,respdf)

      call partition_nnlo_qcd(HardProc, damp, iconf, 1,2)

      respdf = respdf * HardProc%wgt * damp

      kin(1) = respdf(1)
      FintNNLO_ga(1) = kin(1)

      call fill_histo(respdf,vegasweight)

    endif


    !!-----------------------------------------------------------------------!!
    !!                              Collinear 5                              !!
    !!-----------------------------------------------------------------------!!

    C5Lim%ids(1:6) = [0,0,id_el,-id_el,id_q,0]
    call cut_histo(C5Lim)

    if (C5Lim%makecut) then

      kin(2) = zero
      FintNNLO_ga(2) = zero

    else

      z5  = C5Lim%Lim_z(1)
      s15 = C5Lim%Lim_sij(1,5)

      call res_tree_a_qa(C5Lim%AmpMom,res_tmp)

      !-- take only quarks, i.e. res_tmp(1,:)
      ampres(1,1) = res_tmp(2,1)*ndn + res_tmp(2,2)*nup
      call get_respdf(ga_lumi,1,1,C5Lim,ampres,respdf)

      call partition_nnlo_qcd(C5Lim, damp, iconf, 1,2)

      respdf = (two/s15) * Tr * Pqg(z5) * respdf * C5Lim%wgt * damp
      respdf = -respdf

      FintNNLO_ga(2) = respdf(1)
      kin(2) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 6                              !!
    !!-----------------------------------------------------------------------!!

    C6Lim%ids(1:6) = [0,0,id_el,-id_el,id_q,0]
    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(3) = zero
      FintNNLO_ga(3) = zero

    else

      z6  = C6Lim%Lim_z(2)
      s26 = C6Lim%Lim_sij(2,6)

      call res_tree_g_gq(C6Lim%AmpMom,res_tmp)

      !-- take only quarks, i.e. res_tmp(1,:)
      ampres(1,1) = res_tmp(2,1)*ndn*Qdn2 + res_tmp(2,2)*nup*Qup2
      call get_respdf(ga_lumi,1,1,C6Lim,ampres,respdf)

      call partition_nnlo_qcd(C6Lim, damp, iconf, 1,2)

      respdf = (two/s26) * xn * Pqg(z6) * respdf * C6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_ga(3) = respdf(1)
      kin(3) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                       Collinear 5 + Collinear 6                       !!
    !!-----------------------------------------------------------------------!!

    C5C6Lim%ids(1:6) = [0,0,id_el,-id_el,0,0]
    call cut_histo(C5C6Lim)

    if (C5C6Lim%makecut) then

      kin(4) = zero
      FintNNLO_ga(4) = zero

    else

      z5  = C5C6Lim%Lim_z(1)
      z6  = C5C6Lim%Lim_z(2)
      s15 = C5C6Lim%Lim_sij(1,5)
      s26 = C5C6Lim%Lim_sij(2,6)

      call res_tree_qqb(C5C6Lim%AmpMom,res_tmp)
      ampres = (res_tmp(2,1)*ndn*Qdn2 + res_tmp(2,2)*nup*Qup2)
      call get_respdf(ga_lumi,1,1,C5C6Lim,ampres,respdf)

      !-- damp = one
      respdf = (two/s15) * Tr * Pqg(z5) &
             * (two/s26) * xn * Pqg(z6) &
             * respdf * C5C6Lim%wgt

      FintNNLO_ga(4) = respdf(1)
      kin(4) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ga)

#if(_withchecks == 1)
    FintNNLO_rr_dc_ga = FintNNLO_ga
#endif

  end function xsect_nnlo_rr_5162_ga

  function xsect_nnlo_rr_5261_ga(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_dc_ii
    implicit none
    integer :: xsect_nnlo_rr_5261_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc, C5Lim, C6Lim, C5C6Lim
    !--
    real(dp)  :: xx(kNNLO_max_full)
    real(dp)  :: FintNNLO_ga(4),kin(4)
    real(dp)  :: respdf(ipdf)
    real(dp)  :: ampres(1,1),res_tmp(2,2)
    real(dp)  :: z5,z6,s25,s16,damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5261_ga = 0

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
         HardProc, C5Lim=C5Lim, C6Lim=C6Lim, C5C6Lim=C5C6Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_q,id_q]
    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_ga(1) = zero

    else

      call res_tree_ga_ga(HardProc%AmpMom,ampres(1,1))
      call get_respdf(ga_lumi,1,1,HardProc,ampres,respdf)

      call partition_nnlo_qcd(HardProc, damp, iconf, 2,1)

      respdf = respdf * HardProc%wgt * damp

      kin(1) = respdf(1)
      FintNNLO_ga(1) = kin(1)

      call fill_histo(respdf,vegasweight)

    endif


    !!-----------------------------------------------------------------------!!
    !!                              Collinear 5                              !!
    !!-----------------------------------------------------------------------!!

    C5Lim%ids(1:6) = [0,0,id_el,-id_el,id_q,0]
    call cut_histo(C5Lim)

    if (C5Lim%makecut) then

      kin(2) = zero
      FintNNLO_ga(2) = zero

    else

      z5  = C5Lim%Lim_z(1)
      s25 = C5Lim%Lim_sij(2,5)

      call res_tree_g_gq(C5Lim%AmpMom,res_tmp)

      !-- take only quarks, i.e. res_tmp(1,:)
      ampres(1,1) = res_tmp(1,1)*ndn*Qdn2 + res_tmp(1,2)*nup*Qup2
      call get_respdf(ga_lumi,1,1,C5Lim,ampres,respdf)

      call partition_nnlo_qcd(C5Lim, damp, iconf, 2,1)

      respdf = (two/s25) * xn * Pqg(z5) * respdf * C5Lim%wgt * damp
      respdf = -respdf

      FintNNLO_ga(2) = respdf(1)
      kin(2) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              Collinear 6                              !!
    !!-----------------------------------------------------------------------!!

    C6Lim%ids(1:6) = [0,0,id_el,-id_el,id_q,0]
    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(3) = zero
      FintNNLO_ga(3) = zero

    else

      z6  = C6Lim%Lim_z(2)
      s16 = C6Lim%Lim_sij(1,6)

      call res_tree_a_qa(C6Lim%AmpMom,res_tmp)

      !-- take only anti-quarks, i.e. res_tmp(1,:)
      ampres(1,1) = res_tmp(1,1)*ndn + res_tmp(1,2)*nup
      call get_respdf(ga_lumi,1,1,C6Lim,ampres,respdf)

      call partition_nnlo_qcd(C6Lim, damp, iconf, 2,1)

      respdf = (two/s16) * Tr * Pqg(z6) * respdf * C6Lim%wgt * damp
      respdf = -respdf

      FintNNLO_ga(3) = respdf(1)
      kin(3) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                       Collinear 5 + Collinear 6                       !!
    !!-----------------------------------------------------------------------!!

    C5C6Lim%ids(1:6) = [0,0,id_el,-id_el,0,0]
    call cut_histo(C5C6Lim)

    if (C5C6Lim%makecut) then

      kin(4) = zero
      FintNNLO_ga(4) = zero

    else

      z5  = C5C6Lim%Lim_z(1)
      z6  = C5C6Lim%Lim_z(2)
      s25 = C5C6Lim%Lim_sij(2,5)
      s16 = C5C6Lim%Lim_sij(1,6)

      call res_tree_qqb(C5C6Lim%AmpMom,res_tmp)
      ampres = (res_tmp(1,1)*ndn*Qdn2 + res_tmp(1,2)*nup*Qup2)
      call get_respdf(ga_lumi,1,1,C5C6Lim,ampres,respdf)

      !-- damp = one
      respdf = (two/s25) * xn * Pqg(z5) &
             * (two/s16) * Tr * Pqg(z6) &
             * respdf * C5C6Lim%wgt

      FintNNLO_ga(4) = respdf(1)
      kin(4) = respdf(1)

      call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ga)

#if(_withchecks == 1)
    FintNNLO_rr_dc_ga = FintNNLO_ga
#endif

  end function xsect_nnlo_rr_5261_ga

end module mod_xsects_nnlo_rr_ga

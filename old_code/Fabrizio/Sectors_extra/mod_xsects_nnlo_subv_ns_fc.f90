module mod_xsects_nnlo_subv_ns_fc
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_kinematics_lo
  use mod_process
  use mod_aux_sectors
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_amplitudes_loop_ppll
  use mod_hoppet_tools
  use mod_hoppet_nlo
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_subv(1)
#endif

  private

  public :: xsect_nnlo_subvqcd_ns_fc
  public :: xsect_nnlo_subvewk_ns_fc

contains
  
  function xsect_nnlo_subvqcd_ns_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_subvqcd_ns_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf),respdf_3(ipdf)
    real(dp)    :: res_tree(2,2),res_loop(2,2),res_loop_tmp(2,2),bit1,bit2,eta13,eta14

    xsect_nnlo_subvqcd_ns_fc = 0

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
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero
       
    else    

       call res_qcdloop_qqb(LOProc%AmpMom,res_tree,res_loop)
       res_loop_tmp(:,1) = res_loop(:,1) * Qdn2
       res_loop_tmp(:,2) = res_loop(:,2) * Qup2

       !-- z-dependent bit
       call get_respdf_hoppet(xPij,PDFs,ns_lumi,1,1,LOProc,res_loop_tmp,respdf_1,myPDFs1_Lmu=[xPij_Lmu])
       respdf_1 = respdf_1*LOProc%wgt

       call get_respdf_hoppet(PDFs,xPij,ns_lumi,1,1,LOProc,res_loop_tmp,respdf_2,myPDFs2_Lmu=[xPij_Lmu])
       respdf_2 = respdf_2*LOProc%wgt

       !-- FLM[1,2] bit, assuming Emax = sqrt(q2)/2
       eta13 = half*scr(LOProc%AmpMom(:,1),LOProc%AmpMom(:,3))/LoProc%AmpMom(1,1)/LOProc%AmpMom(1,3)
       eta14 = half*scr(LOProc%AmpMom(:,1),LOProc%AmpMom(:,4))/LoProc%AmpMom(1,1)/LOProc%AmpMom(1,4)
       bit1 = 13._dp - two/three*pisq
       bit2 = log(eta13/eta14)*(three/two)+real(dilog2(one-eta13),kind=dp)-real(dilog2(one-eta14),kind=dp)

       res_loop_tmp(1,:) = (Q_lep2 * bit1 + four*Q_lep*[Qdn,Qup] * bit2)*res_loop(1,:)
       res_loop_tmp(2,:) = (Q_lep2 * bit1 - four*Q_lep*[Qdn,Qup] * bit2)*res_loop(2,:)

       call get_respdf(ns_lumi,1,1,LOProc,res_loop_tmp,respdf_3)
       respdf_3 = respdf_3*LOproc%wgt
       
       respdf = respdf_1 + respdf_2 + respdf_3
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNNLO_subv = kin
#endif

  end function xsect_nnlo_subvqcd_ns_fc

  function xsect_nnlo_subvewk_ns_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_subvewk_ns_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf)
    real(dp)    :: res_tree(2,3),res_loop(2,3)

    xsect_nnlo_subvewk_ns_fc = 0

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
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero
       
    else    

       call res_ewkloop_qqb(LOProc%AmpMom,res_tree,res_loop)

       call get_respdf_hoppet(xPij,PDFs,ns_lumi_splitb,1,1,LOProc,res_loop,respdf_1,myPDFs1_Lmu=[xPij_Lmu])
       respdf_1 = Cf*respdf_1*LOProc%wgt

       call get_respdf_hoppet(PDFs,xPij,ns_lumi_splitb,1,1,LOProc,res_loop,respdf_2,myPDFs2_Lmu=[xPij_Lmu])
       respdf_2 = Cf*respdf_2*LOProc%wgt

       respdf = respdf_1 + respdf_2
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNNLO_subv = kin
#endif

  end function xsect_nnlo_subvewk_ns_fc
  
end module mod_xsects_nnlo_subv_ns_fc
  

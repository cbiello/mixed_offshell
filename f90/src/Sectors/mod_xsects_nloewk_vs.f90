module mod_xsects_nloewk_vs
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
  real(dp), public, save :: FintNLOEWK_vs(1)
#endif

  private

  public :: xsect_nloewk_v_ns,xsect_nloewk_v_aa
  public :: xsect_nloewk_s_ns,xsect_nloewk_s_aa
  public :: xsect_nloewk_s_aq,xsect_nloewk_s_qa

contains

  function xsect_nloewk_v_ns(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_v_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_tree(2,3),res_loop(2,3)

    xsect_nloewk_v_ns = 0

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
       call get_respdf(ns_lumi_splitb,0,1,LOProc,res_loop,respdf)

       respdf = respdf*LOProc%wgt
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOEWK_vs = kin
#endif

  end function xsect_nloewk_v_ns

  function xsect_nloewk_v_aa(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_v_aa
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_tree,res_loop(1,1)

    xsect_nloewk_v_aa = 0

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

       call res_ewkloopAA_aa(LOProc%AmpMom,res_tree,res_loop(1,1))
       call get_respdf(aa_lumi,0,1,LOProc,res_loop,respdf)

       respdf = respdf*LOProc%wgt
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOEWK_vs = kin
#endif

  end function xsect_nloewk_v_aa

  !-- subtractions below
  
  function xsect_nloewk_s_ns(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_s_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf),respdf_3(ipdf)
    real(dp)    :: res_lo(2,2),res_lo_tmp(2,2),bit1,bit2,eta13,eta14

    xsect_nloewk_s_ns = 0

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

       call res_tree_qqb(LOProc%AmpMom,res_lo)
       res_lo_tmp(:,1) = res_lo(:,1) * Qdn2
       res_lo_tmp(:,2) = res_lo(:,2) * Qup2

       !-- z-dependent bit
       call get_respdf_hoppet(xPij,PDFs,ns_lumi,0,1,LOProc,res_lo_tmp,respdf_1,myPDFs1_Lmu=[xPij_Lmu])
       respdf_1 = respdf_1*LOProc%wgt

       call get_respdf_hoppet(PDFs,xPij,ns_lumi,0,1,LOProc,res_lo_tmp,respdf_2,myPDFs2_Lmu=[xPij_Lmu])
       respdf_2 = respdf_2*LOProc%wgt

       !-- FLM[1,2] bit, assuming Emax = sqrt(q2)/2
       eta13 = half*scr(LOProc%AmpMom(:,1),LOProc%AmpMom(:,3))/LoProc%AmpMom(1,1)/LOProc%AmpMom(1,3)
       eta14 = half*scr(LOProc%AmpMom(:,1),LOProc%AmpMom(:,4))/LoProc%AmpMom(1,1)/LOProc%AmpMom(1,4)
       bit1 = 13._dp - two/three*pisq
       bit2 = log(eta13/eta14)*(three/two)+real(dilog2(one-eta13),kind=dp)-real(dilog2(one-eta14),kind=dp)

       res_lo_tmp(1,:) = (Q_lep2 * bit1 + four*Q_lep*[Qdn,Qup] * bit2)*res_lo(1,:)
       res_lo_tmp(2,:) = (Q_lep2 * bit1 - four*Q_lep*[Qdn,Qup] * bit2)*res_lo(2,:)

       call get_respdf(ns_lumi,0,1,LOProc,res_lo_tmp,respdf_3)
       respdf_3 = respdf_3*LOproc%wgt
       
       respdf = respdf_1 + respdf_2 + respdf_3
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOEWK_vs = kin
#endif

  end function xsect_nloewk_s_ns

  !--

  function xsect_nloewk_s_aa(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_s_aa
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf),aa_int_sub(ipdf),sv_logs(ipdf)
    real(dp)    :: res_lo(1,1),eta13,eta14,EC

    xsect_nloewk_s_aa = 0

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

       call res_treeAA_aa(LOProc%AmpMom,res_lo(1,1))
       call get_respdf(aa_lumi,0,1,LOProc,res_lo,respdf)

       EC = LOProc%Lim_Ei(1)
       call fill_sv_logs(LOProc%muf(1)**2,4*EC**2,sv_logs)
       eta13 = LOProc%Lim_etaij(1,3)
       eta14 = LOProc%Lim_etaij(1,4)
       aa_int_sub = (13._dp - two/three*pisq)*Q_lep**2 - 2*sv_logs * gamma_a

       respdf = aa_int_sub * respdf * LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOEWK_vs = kin
#endif

  end function xsect_nloewk_s_aa
  
  !--

  function xsect_nloewk_s_aq(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_s_aq
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf_vect(1:2,ipdf),respdf(ipdf)
    real(dp)    :: res_lo(2,2)

    xsect_nloewk_s_aq = 0

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

       call res_tree_qqb(LOProc%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       !-- Pqa*sigma_qq, ns_lumi because the convolution is done in xPij
       call get_respdf_hoppet(xPij,PDFs,ns_lumi,0,1,LOProc,res_lo,respdf_vect(1,:),myPDFs1_Lmu=[xPij_Lmu])
       respdf_vect(1,:) = xn*respdf_vect(1,:)*LOProc%wgt

       !--

       !-- Paq*sigma_aa, quark charges already in the luminosity
       call res_treeAA_aa(LOProc%AmpMom,res_lo(1,1))
       call get_respdf_hoppet(PDFs,xPij_2,aa_lumi,0,1,LOProc,res_lo,respdf_vect(2,:),myPDFs2_Lmu=[xPij_2_Lmu])
       respdf_vect(2,:) = respdf_vect(2,:)*LOProc%wgt

       respdf = sum(respdf_vect(1:2,:),1)
       
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOEWK_vs = kin
#endif

  end function xsect_nloewk_s_aq

  function xsect_nloewk_s_qa(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_s_qa
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf_vect(2,ipdf),respdf(ipdf)
    real(dp)    :: res_lo(2,2)

    xsect_nloewk_s_qa = 0

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

       call res_tree_qqb(LOProc%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       !-- Pqa*sigma_qq, ns_lumi because the convolution is done in xPij
       call get_respdf_hoppet(PDFs,xPij,ns_lumi,0,1,LOProc,res_lo,respdf_vect(1,:),myPDFs2_Lmu=[xPij_Lmu])
       respdf_vect(1,:) = xn*respdf_vect(1,:)*LOProc%wgt

       !--

       !-- Paq*sigma_aa, quark charges already in the luminosity
       call res_treeAA_aa(LOProc%AmpMom,res_lo(1,1))
       call get_respdf_hoppet(xPij_2,PDFs,aa_lumi,0,1,LOProc,res_lo,respdf_vect(2,:),myPDFs1_Lmu=[xPij_2_Lmu])
       respdf_vect(2,:) = respdf_vect(2,:)*LOProc%wgt

       respdf = sum(respdf_vect(1:2,:),1)

       kin(1) = respdf(1)
      
       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOEWK_vs = kin
#endif

  end function xsect_nloewk_s_qa
  
end module mod_xsects_nloewk_vs
  

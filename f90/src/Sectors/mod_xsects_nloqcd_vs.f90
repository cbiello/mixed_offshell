module mod_xsects_nloqcd_vs
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
  use mod_amplitudes_tree_ppnul
  use mod_amplitudes_loop_ppll
  use mod_amplitudes_loop_pplnu
  use mod_hoppet_tools
  use mod_hoppet_nlo
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNLOQCD_vs(1)
#endif

  private

  public :: xsect_nloqcd_v_ns
  public :: xsect_nloqcd_s_ns
  public :: xsect_nloqcd_s_gq,xsect_nloqcd_s_qg

  public :: xsect_nloqcd_v_ns_wp ! deprecated -> remove
  public :: xsect_nloqcd_s_ns_wp ! deprecated -> remove

contains

  function xsect_nloqcd_v_ns(yRnd,ff,vegasweight)
    integer :: xsect_nloqcd_v_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_tree(-5:7,-5:7),res_loop(-5:7,-5:7)
    !--
    real(dp)    :: res_tree_old(2,2),res_loop_old(2,2)
    logical :: oldcode

    oldcode = .false.

    xsect_nloqcd_v_ns = 0

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

    ! define process specific partons
#if (_Vcharge == 0)
    LOProc%part(1:4) = [id_q,-id_q,id_el,-id_el]
#elif  (_Vcharge == -1)
    LOProc%part(1:4) = [id_q,-id_qp,id_el,-id_nue]
#elif  (_Vcharge == +1)
    LOProc%part(1:4) = [id_q,-id_qp,id_nue,-id_el]
#endif

    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero
       
    else    

       if(oldcode) then
          call res_qcdloop_qqb(LOProc%AmpMom,res_tree_old,res_loop_old)
          call get_respdf(ns_lumi,1,0,LOProc,res_loop_old,respdf)
       else 
          call res_qcdloop_qqb_gen(LOProc%AmpMom,res_tree,res_loop)
          call get_respdf_gen(1,0,LOProc,res_loop,respdf)
       endif

       respdf = respdf*LOProc%wgt
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOQCD_vs = kin
#endif

  end function xsect_nloqcd_v_ns



  function xsect_nloqcd_v_ns_wp(yRnd,ff,vegasweight)
    integer :: xsect_nloqcd_v_ns_wp
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_tree(2,2),res_loop(1,2)

    xsect_nloqcd_v_ns_wp = 0

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

       call res_qcdloop_qqb_wp(LOProc%AmpMom,res_tree,res_loop)
       call get_respdf(qQpb_lumi_wp,1,0,LOProc,res_loop,respdf)

       respdf = respdf*LOProc%wgt
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOQCD_vs = kin
#endif

  end function xsect_nloqcd_v_ns_wp


  
  function xsect_nloqcd_s_ns(yRnd,ff,vegasweight)
    integer :: xsect_nloqcd_s_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf)
    real(dp)    :: res_lo(2,2)

    xsect_nloqcd_s_ns = 0

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

       call get_respdf_hoppet(xPij,PDFs,ns_lumi,1,0,LOProc,res_lo,respdf_1,myPDFs1_Lmu=[xPij_Lmu])
       respdf_1 = Cf*respdf_1*LOProc%wgt

       call get_respdf_hoppet(PDFs,xPij,ns_lumi,1,0,LOProc,res_lo,respdf_2,myPDFs2_Lmu=[xPij_Lmu])
       respdf_2 = Cf*respdf_2*LOProc%wgt

       respdf = respdf_1 + respdf_2
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOQCD_vs = kin
#endif

  end function xsect_nloqcd_s_ns


  
  function xsect_nloqcd_s_ns_wp(yRnd,ff,vegasweight)
    integer :: xsect_nloqcd_s_ns_wp
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf)
    real(dp)    :: res_lo(1,2)

    xsect_nloqcd_s_ns_wp = 0

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

       call res_tree_qqb_wp(LOProc%AmpMom,res_lo)

       call get_respdf_hoppet(xPij,PDFs,qQpb_lumi_wp,1,0,LOProc,res_lo,respdf_1,myPDFs1_Lmu=[xPij_Lmu])
       respdf_1 = Cf*respdf_1*LOProc%wgt

       call get_respdf_hoppet(PDFs,xPij,qQpb_lumi_wp,1,0,LOProc,res_lo,respdf_2,myPDFs2_Lmu=[xPij_Lmu])
       respdf_2 = Cf*respdf_2*LOProc%wgt

       respdf = respdf_1 + respdf_2
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOQCD_vs = kin
#endif

  end function xsect_nloqcd_s_ns_wp
  






  function xsect_nloqcd_s_gq(yRnd,ff,vegasweight)
    integer :: xsect_nloqcd_s_gq
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_lo(2,2)

    xsect_nloqcd_s_gq = 0

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

       !-- ns_lumi because the convolution is done in xPij
       call get_respdf_hoppet(xPij,PDFs,ns_lumi,1,0,LOProc,res_lo,respdf,myPDFs1_Lmu=[xPij_Lmu])
       respdf = tr*respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOQCD_vs = kin
#endif

  end function xsect_nloqcd_s_gq

  function xsect_nloqcd_s_qg(yRnd,ff,vegasweight)
    integer :: xsect_nloqcd_s_qg
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_lo(2,2)

    xsect_nloqcd_s_qg = 0

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

       !-- ns_lumi because the convolution is done in xPij
       call get_respdf_hoppet(PDFs,xPij,ns_lumi,1,0,LOProc,res_lo,respdf,myPDFs2_Lmu=[xPij_Lmu])
       respdf = tr*respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOQCD_vs = kin
#endif

  end function xsect_nloqcd_s_qg

end module mod_xsects_nloqcd_vs
  

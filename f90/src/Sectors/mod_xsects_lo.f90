module mod_xsects_lo
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_kinematics_gen
  use mod_kinematics_lo
  use mod_process
  use mod_auxfunctions
  use mod_aux_sectors
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppnul
  use mod_amplitudes_tree_ppll
  
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintLO(1)
#endif

  private

  public :: xsect_lo_ns,xsect_lo_aa
  public :: xsect_lo_ns_w

contains

  function xsect_lo_ns(yRnd,ff,vegasweight)
    integer :: xsect_lo_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_lo(2,2)

    xsect_lo_ns = 0

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
       call get_respdf(ns_lumi,0,0,LOProc,res_lo,respdf)

       respdf = respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintLO = kin
#endif

  end function xsect_lo_ns

  function xsect_lo_aa(yRnd,ff,vegasweight)
    integer :: xsect_lo_aa
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_lo(1,1)

    xsect_lo_aa = 0

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
       call get_respdf(aa_lumi,0,0,LOProc,res_lo,respdf)

       respdf = respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintLO = kin
#endif

  end function xsect_lo_aa

  
  function xsect_lo_ns_w(yRnd,ff,vegasweight)
    integer :: xsect_lo_ns_w
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_lo(2,2)

    xsect_lo_ns_w = 0

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

       !TEST MADGRAPH
       !LOProc%AmpMom(:,1) = (/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03/)
       !LOProc%AmpMom(:,2) = (/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  -0.5000000E+03/)
       !LOProc%AmpMom(:,3) = (/0.5000000E+03,  0.1109243E+03,  0.4448308E+03, -0.1995529E+03/)
       !LOProc%AmpMom(:,4) = (/0.5000000E+03,  -0.1109243E+03,  -0.4448308E+03, 0.1995529E+03/)

       !call res_tree_qqb_w(LOProc%AmpMom,res_lo)
       
       !print*, 'LO amp MG udx_veep  ', (0.094835522759998875_dp)**2*res_lo(1,1)/eesq2
       !print*, 'MG udx_veep            ', '1.3875944806836576E-003' 
       !print*, 'LO amp MG dux_vexem ', (0.094835522759998875_dp)**2*res_lo(2,1)/eesq2
       !print*, 'MG dux_vexem           ', '7.5225966654526560E-003'
       !print*, 'LO amp MG uxd_vexem ', (0.094835522759998875_dp)**2*res_lo(1,2)/eesq2
       !print*, 'MG uxd_vexem           ', '1.3875944806836580E-003'
       !print*, 'LO amp MG dxu_veep  ', (0.094835522759998875_dp)**2*res_lo(2,2)/eesq2
       !print*, 'MG dxu_veep            ', '7.5225966654526560E-003'
       !stop

       call res_tree_qqb_w(LOProc%AmpMom,res_lo)
       call get_respdf(qQpb_lumi_w,0,0,LOProc,res_lo,respdf)

       respdf = respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintLO = kin
#endif

  end function xsect_lo_ns_w

end module mod_xsects_lo

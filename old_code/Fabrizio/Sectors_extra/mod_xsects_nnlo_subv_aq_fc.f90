module mod_xsects_nnlo_subv_aq_fc
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
  real(dp), public, save :: Fintnnlo_subv_aq(1)
#endif

  private

  public :: xsect_nnlo_subvqcd_aq_fc

contains

  function xsect_nnlo_subvqcd_aq_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_subvqcd_aq_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_tree(2,2),res_loop(2,2)

    xsect_nnlo_subvqcd_aq_fc = 0

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
       res_loop(:,1) = Qdn2 * res_loop(:,1)
       res_loop(:,2) = Qup2 * res_loop(:,2)

       !-- Pqa*sigma_qq, ns_lumi because the convolution is done in xPij
       call get_respdf_hoppet(xPij,PDFs,ns_lumi,1,1,LOProc,res_loop,respdf,myPDFs1_Lmu=[xPij_Lmu])
       respdf = xn*respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    Fintnnlo_subv_aq = kin
#endif

  end function xsect_nnlo_subvqcd_aq_fc

end module mod_xsects_nnlo_subv_aq_fc
  

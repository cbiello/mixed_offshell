module mod_run_proc_s_qg
  use mod_types
  use mod_consts_dp
#if ( _MPIbuild )
  use mod_vegas_mpi
#else
  use mod_my_vegas
#endif
  use mod_proc_parms
  implicit none
  private

  public :: run_proc_s_qg,help_run_proc_s_qg

contains

  subroutine run_proc_s_qg(vg_result,vg_error,vg_chi2)
    use mod_kinematics_lo, only: kLO_max_full
    use mod_kinematics_nlo_z, only: kNLO_max_full
    use mod_hoppet_nnlo
    use mod_xsects_nnlo_s_qg
    use mod_amplitudes_loop_ppll, only: ew_renorm_nf_1L
    implicit none
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2

    call init_xPij_nnlo()
    if     (sec.eq.'s_12')         then; call vegas_integrate(kLO_max_full,xsect_nnlo_s_qg_s12,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_vewk')       then; call vegas_integrate(kLO_max_full,xsect_nnlo_s_qg_vewk,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_oqcd')       then; call vegas_integrate(kNLO_max_full,xsect_nnlo_s_qg_oqcd,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_oewk_is')    then; call vegas_integrate(kNLO_max_full,xsect_nnlo_s_qg_oewk_is,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_oewk_fs_53') then; call vegas_integrate(kNLO_max_full,xsect_nnlo_s_qg_oewk_fs_53,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_oewk_fs_54') then; call vegas_integrate(kNLO_max_full,xsect_nnlo_s_qg_oewk_fs_54,vg_result,vg_error,vg_chi2)
    !-- nf
    elseif (sec.eq.'s_vewknf') then
      call ew_renorm_nf_1L()
      call vegas_integrate(kLO_max_full,xsect_nnlo_s_qg_vewknf,vg_result,vg_error,vg_chi2)
    else
      print *, 'wrong channel for -ch qg -sec s, type -h for help'
      stop
    endif

  end subroutine run_proc_s_qg

  subroutine help_run_proc_s_qg(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -corr nnlo   -ch qg -sec s_12, s_vewk, s_oqcd, s_oewk_is, s_oewk_fs_5[3,4] '
    write(idev,*) ''

  end subroutine help_run_proc_s_qg

end module mod_run_proc_s_qg

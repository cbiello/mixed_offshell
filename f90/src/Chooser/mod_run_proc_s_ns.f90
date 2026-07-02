module mod_run_proc_s_ns
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

  public :: run_proc_s_ns,help_run_proc_s_ns

contains

  subroutine run_proc_s_ns(vg_result,vg_error,vg_chi2)
    use mod_kinematics_lo,  only: kLO_max_full
    use mod_kinematics_nlo, only: kNLO_max_full
    use mod_hoppet_nnlo
    use mod_xsects_nnlo_s_ns
    use mod_amplitudes_loop_ppll, only: ew_renorm_nf_1L
    implicit none
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2

    call init_xPij_nnlo()
    print *, sec
    if     (sec.eq.'s_12')   then; call vegas_integrate(kLO_max_full+1,xsect_nnlo_s_ns_s12,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_vqcd') then; call vegas_integrate(kLO_max_full,xsect_nnlo_s_ns_vqcd,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_vewk') then; call vegas_integrate(kLO_max_full,xsect_nnlo_s_ns_vewk,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_qqb')  then; call vegas_integrate(kLO_max_full,xsect_nnlo_s_ns_qqb,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_qq')   then; call vegas_integrate(kLO_max_full,xsect_nnlo_s_ns_qq,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_oqcd') then; call vegas_integrate(kNLO_max_full,xsect_nnlo_s_ns_oqcd,vg_result,vg_error,vg_chi2)

    elseif (sec.eq.'s_oewk_is') then; call vegas_integrate(kNLO_max_full,xsect_nnlo_s_ns_oewk_is,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_oewk_is_ra') then; call vegas_integrate(kNLO_max_full,xsect_nnlo_s_ns_oewk_is_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_oewk_fs_53') then; call vegas_integrate(kNLO_max_full,xsect_nnlo_s_ns_oewk_fs_53_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_oewk_fs_54') then; call vegas_integrate(kNLO_max_full,xsect_nnlo_s_ns_oewk_fs_54_raoul,vg_result,vg_error,vg_chi2)
    !-- nf
    elseif (sec.eq.'s_vewknf') then
      call ew_renorm_nf_1L()
      call vegas_integrate(kLO_max_full,xsect_nnlo_s_ns_vewknf,vg_result,vg_error,vg_chi2)
    else
      print *, 'wrong channel for -ch ns/ns_qqb/ns_qq -sec s, type -h for help'
      stop
    endif

  end subroutine run_proc_s_ns

  subroutine help_run_proc_s_ns(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -corr nnlo   -ch ns  -sec s_12, s_vqcd, s_vewk, s_oqcd, s_oewk_is, s_oewk_fs_5[3,4] '
    write(idev,*) ' -corr nnlo   -ch ns  -sec s_qqb'
    write(idev,*) ' -corr nnlo   -ch ns  -sec s_qq '
    write(idev,*) ''
    
  end subroutine help_run_proc_s_ns
  
end module mod_run_proc_s_ns

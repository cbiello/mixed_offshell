module mod_run_proc_s_aq
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

  public :: run_proc_s_aq,help_run_proc_s_aq

contains

  subroutine run_proc_s_aq(vg_result,vg_error,vg_chi2)
    use mod_kinematics_lo,  only: kLO_max_full
    use mod_kinematics_nlo, only: kNLO_max_full
    use mod_hoppet_nnlo
    use mod_xsects_nnlo_s_aq
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2

    call init_xPij_nnlo()
    if     (sec.eq.'s_12')   then; call vegas_integrate(kLO_max_full,xsect_nnlo_s_aq_s12,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_vqcd') then; call vegas_integrate(kLO_max_full,xsect_nnlo_s_aq_vqcd,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_oqcd') then; call vegas_integrate(kNLO_max_full,xsect_nnlo_s_aq_oqcd,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'s_oewk') then; call vegas_integrate(kNLO_max_full,xsect_nnlo_s_aq_oewk,vg_result,vg_error,vg_chi2)
    else
      print *, 'wrong channel for -ch aq -sec s, type -h for help'
      stop
    endif

  end subroutine run_proc_s_aq

  subroutine help_run_proc_s_aq(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -corr nnlo   -ch aq -sec s_12, s_vqcd, s_oqcd, s_oewk'
    
  end subroutine help_run_proc_s_aq
  
end module mod_run_proc_s_aq

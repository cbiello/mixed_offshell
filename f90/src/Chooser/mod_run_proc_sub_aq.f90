module mod_run_proc_sub_aq
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

  public :: run_proc_sub_aq,help_run_proc_sub_aq
  
contains

  subroutine run_proc_sub_aq(vg_result,vg_error,vg_chi2)
    use mod_kinematics_lo
    use mod_kinematics_nlo_z_raoul
    use mod_xsects_nnlo_s_aq_raoul
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2


    if (sec.eq.'sub12') then
       call vegas_integrate(LOdim_vegas+2,xsect_nnlo_sub12_aq_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub_onloqcd') then
       call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_onloqcd_aq_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub_onloewk') then
       call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_onloewk_aq_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'subv_qcd') then
       call vegas_integrate(LOdim_vegas+1,xsect_nnlo_subvqcd_aq_raoul,vg_result,vg_error,vg_chi2)

    else
       print *, 'wrong sector, type -h for help'
       stop
    endif
    
  end subroutine run_proc_sub_aq

  subroutine help_run_proc_sub_aq(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -ch aq -sec sub_onloqcd sub_onloewk sub12 subv_qcd'
    write(idev,*) ''
    
  end subroutine help_run_proc_sub_aq
  
end module mod_run_proc_sub_aq


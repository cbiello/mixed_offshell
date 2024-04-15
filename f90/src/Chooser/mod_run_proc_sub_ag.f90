module mod_run_proc_sub_ag
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

  public :: run_proc_sub_ag,help_run_proc_sub_ag
  
contains

  subroutine run_proc_sub_ag(vg_result,vg_error,vg_chi2)
    use mod_kinematics_lo
    use mod_kinematics_nlo_z_raoul
    use mod_xsects_nnlo_s_ag_raoul
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2


    if (sec.eq.'sub12') then
       call vegas_integrate(LOdim_vegas+2,xsect_nnlo_sub12_ag_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub_onloqcd') then
       call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_onloqcd_ag_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub_onloewk') then
       call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_onloewk_ag_raoul,vg_result,vg_error,vg_chi2)
    else
       print *, 'wrong sector, type -h for help'
       stop
    endif
    
  end subroutine run_proc_sub_ag

  subroutine help_run_proc_sub_ag(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -ch ag -sec sub_onloqcd sub_onloewk sub12'
    write(idev,*) ''
    
  end subroutine help_run_proc_sub_ag
  
end module mod_run_proc_sub_ag


module mod_run_proc_sub_ns
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

  public :: run_proc_sub_ns,help_run_proc_sub_ns
  
contains

  subroutine run_proc_sub_ns(vg_result,vg_error,vg_chi2)
    use mod_kinematics_lo
    use mod_kinematics_nlo_z_raoul
    use mod_xsects_nnlo_s_qqb_raoul
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2


    if (sec.eq.'sub12') then
       call vegas_integrate(LOdim_vegas+2,xsect_nnlo_sub12_ns_ga_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub_onloqcd') then
       call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_onloqcd_ns_ga_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub_onloewk_is') then
       call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_onloewk_is_ns_ga_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub_onloewk_53') then
       call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_onloewk_fs_53_ns_ga_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub_onloewk_54') then
       call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_onloewk_fs_54_ns_ga_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'subv_ewk') then
       call vegas_integrate(LOdim_vegas+1,xsect_nnlo_subvewk_ns_ga_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'subv_qcd') then
       call vegas_integrate(LOdim_vegas+1,xsect_nnlo_subvqcd_ns_ga_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub12_qqb') then
       call vegas_integrate(LOdim_vegas+1,xsect_nnlo_sub12_ns_qqb_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub12_qq') then
       call vegas_integrate(LOdim_vegas+1,xsect_nnlo_sub12_ns_qq_raoul,vg_result,vg_error,vg_chi2)

    else
       print *, 'wrong sector, type -h for help'
       stop
    endif
    
  end subroutine run_proc_sub_ns

  subroutine help_run_proc_sub_ns(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -ch ns -sec sub12 sub_onloqcd sub_onloewk_is sub_onloewk_53 sub_onloewk_54 subv_ewk subv_qcd sub12_qqb sub12_qq'
    write(idev,*) ''
    
  end subroutine help_run_proc_sub_ns
  
end module mod_run_proc_sub_ns


module mod_run_proc_sub_gq
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

  public :: run_proc_sub_gq,help_run_proc_sub_gq
  
contains

  subroutine run_proc_sub_gq(vg_result,vg_error,vg_chi2)
    use mod_kinematics_lo
    use mod_kinematics_nlo_z_raoul
    use mod_xsects_nnlo_s_gq_raoul
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2


    if (sec.eq.'sub12') then
       call vegas_integrate(LOdim_vegas+2,xsect_nnlo_sub12_gq_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub_onloqcd') then
       call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_onloqcd_gq_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub_onloewk_is') then
       call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_onloewk_is_gq_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub_onloewk_53') then
       call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_onloewk_fs_53_gq_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'sub_onloewk_54') then
       call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_onloewk_fs_54_gq_raoul,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'subv_ewk') then
       call vegas_integrate(LOdim_vegas+1,xsect_nnlo_subvewk_gq_raoul,vg_result,vg_error,vg_chi2)

    else
       print *, 'wrong sector, type -h for help'
       stop
    endif
    
  end subroutine run_proc_sub_gq

  subroutine help_run_proc_sub_gq(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -ch gq -sec sub_onloqcd sub_onloewk_is sub_onloewk_53 sub_onloewk_54 sub12 subv_ewk'
    write(idev,*) ''
    
  end subroutine help_run_proc_sub_gq
  
end module mod_run_proc_sub_gq


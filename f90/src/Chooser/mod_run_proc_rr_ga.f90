module mod_run_proc_rr_ga
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

  public :: run_proc_rr_ga,help_run_proc_rr_ga

contains

  subroutine run_proc_rr_ga(vg_result,vg_error,vg_chi2)
    use mod_kinematics_nnlo_tc_is_eu_ac, only: NNLOdim_vegas
    use mod_xsects_nnlo_rr_ga
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2

    if (sec.eq.'rr_5161a') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161a_ga,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5161b') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161b_ga,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5161c') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161c_ga,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5161d') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161d_ga,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5262') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262_ga,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5162') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5162_ga,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5261') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5261_ga,vg_result,vg_error,vg_chi2)
    else
       print *, 'wrong sector, type -h for help'
       stop
    endif

  end subroutine run_proc_rr_ga

  subroutine help_run_proc_rr_ga(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -corr nnlo   -ch ga -sec rr_5161[a,b,c,d] rr_5262, rr_5162 rr_5261'
    write(idev,*) ''

  end subroutine help_run_proc_rr_ga

end module mod_run_proc_rr_ga

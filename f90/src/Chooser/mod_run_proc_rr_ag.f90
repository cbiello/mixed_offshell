module mod_run_proc_rr_ag
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

  public :: run_proc_rr_ag,help_run_proc_rr_ag

contains

  subroutine run_proc_rr_ag(vg_result,vg_error,vg_chi2)
    use mod_kinematics_nnlo_tc_is_eu_ac, only: NNLOdim_vegas
    use mod_xsects_nnlo_rr_ag
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2

    if (sec.eq.'rr_5161') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161_ag,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5262a') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262a_ag,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5262b') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262b_ag,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5262c') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262c_ag,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5262d') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262d_ag,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5261') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5261_ag,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5162') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5162_ag,vg_result,vg_error,vg_chi2)
    else
       print *, 'wrong sector, type -h for help'
       stop
    endif

  end subroutine run_proc_rr_ag

  subroutine help_run_proc_rr_ag(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -corr nnlo   -ch ag -sec rr_5161, rr_5262[a,b,c,d] rr_5162 rr_5261'

  end subroutine help_run_proc_rr_ag

end module mod_run_proc_rr_ag

module mod_run_proc_rr_qa
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

  public :: run_proc_rr_qa,help_run_proc_rr_qa

contains

  subroutine run_proc_rr_qa(vg_result,vg_error,vg_chi2)
    use mod_kinematics_nnlo_tc_is_eu_abcd, only: NNLOdim_vegas
    use mod_xsects_nnlo_rr_qa
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2

    if (sec.eq.'rr_5161a') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161a_qa,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5161b') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161b_qa,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5161c') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161c_qa,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5161d') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161d_qa,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5262a') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262a_qa,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5262b') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262b_qa,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5262c') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262c_qa,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5262d') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262d_qa,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5162') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5162_qa,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5261') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5261_qa,vg_result,vg_error,vg_chi2)
    else
       print *, 'wrong sector, type -h for help'
       stop
    endif

  end subroutine run_proc_rr_qa

  subroutine help_run_proc_rr_qa(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -corr nnlo   -ch qa -sec rr_5161[a,b,c,d] rr_5262[a,b,c,d] rr_5162 rr_5261'
    write(idev,*) ''

  end subroutine help_run_proc_rr_qa

end module mod_run_proc_rr_qa

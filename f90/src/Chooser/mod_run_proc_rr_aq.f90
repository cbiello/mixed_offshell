module mod_run_proc_rr_aq
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

  public :: run_proc_rr_aq,help_run_proc_rr_aq

contains

  subroutine run_proc_rr_aq(vg_result,vg_error,vg_chi2)
    use mod_kinematics_nnlo_tc_is_eu_abcd, only: NNLOdim_vegas
    use mod_xsects_nnlo_rr_aq
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2

    if (sec.eq.'rr_5161a') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161a_aq,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5161b') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161b_aq,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5161c') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161c_aq,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5161d') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161d_aq,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5262a') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262a_aq,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5262b') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262b_aq,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5262c') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262c_aq,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5262d') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262d_aq,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5162') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5162_aq,vg_result,vg_error,vg_chi2)
    else if (sec.eq.'rr_5261') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5261_aq,vg_result,vg_error,vg_chi2)
    else
       print *, 'wrong sector, type -h for help'
       stop
    endif

  end subroutine run_proc_rr_aq

  subroutine help_run_proc_rr_aq(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -corr nnlo   -ch aq -sec rr_5161[a,b,c,d] rr_5262[a,b,c,d] rr_5162 rr_5261'

  end subroutine help_run_proc_rr_aq

end module mod_run_proc_rr_aq

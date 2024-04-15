module mod_run_proc_rr_qg
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

  public :: run_proc_rr_qg,help_run_proc_rr_qg

contains

  subroutine run_proc_rr_qg(vg_result,vg_error,vg_chi2)
    use mod_kinematics_nnlo_tc_is_eu_abcd, only: NNLOdim_vegas
    use mod_xsects_nnlo_rr_qg
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2

    if (sec.eq.'rr_5262a') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262a_qg,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'rr_5262b') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262b_qg,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'rr_5262c') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262c_qg,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'rr_5262d') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262d_qg,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'rr_5261') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5261_qg,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'rr_5263') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5263_qg,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'rr_5264') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5264_qg,vg_result,vg_error,vg_chi2)
    else
       print *, 'wrong sector, type -h for help'
       stop
    endif

  end subroutine run_proc_rr_qg

  subroutine help_run_proc_rr_qg(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -corr nnlo   -ch qg -sec rr_5262[a,b,c,d] rr_526[1,3,4]'
    write(idev,*) ''

  end subroutine help_run_proc_rr_qg

end module mod_run_proc_rr_qg

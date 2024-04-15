module mod_run_proc_rr_gq
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

  public :: run_proc_rr_gq,help_run_proc_rr_gq

contains

  subroutine run_proc_rr_gq(vg_result,vg_error,vg_chi2)
    use mod_kinematics_nnlo_tc_is_eu_abcd, only: NNLOdim_vegas
    use mod_xsects_nnlo_rr_gq
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2

    if (sec.eq.'rr_5161a') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161a_gq,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'rr_5161b') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161b_gq,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'rr_5161c') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161c_gq,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'rr_5161d') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161d_gq,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'rr_5162') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5162_gq,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'rr_5163') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5163_gq,vg_result,vg_error,vg_chi2)
    elseif (sec.eq.'rr_5164') then
       call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5164_gq,vg_result,vg_error,vg_chi2)
    else
       print *, 'wrong sector, type -h for help'
       stop
    endif

  end subroutine run_proc_rr_gq

  subroutine help_run_proc_rr_gq(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -corr nnlo   -ch gq -sec rr_5161[a,b,c,d] rr_516[2,3,4]'

  end subroutine help_run_proc_rr_gq

end module mod_run_proc_rr_gq

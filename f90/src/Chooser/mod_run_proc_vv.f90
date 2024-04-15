module mod_run_proc_vv
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

  public :: run_proc_vv,help_run_proc_vv
  
contains

  subroutine run_proc_vv(vg_result,vg_error,vg_chi2)
    use mod_kinematics_lo
    use mod_xsects_nnlo_vv
    use mod_amplitudes_loop_ppll, only: ew_renorm_nf_1L
    use mod_amplitudes_twol_nf_ppll, only: ew_renorm_nf_2L
    implicit none
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2

    if (ch.eq.'ns') then
       if     (sec.eq.'vvnf') then
         call ew_renorm_nf_2L()
         call vegas_integrate(LOdim_vegas,xsect_nnlo_vv_ns_nf,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'vvfc') then; call vegas_integrate(LOdim_vegas,xsect_nnlo_vv_ns_fc,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'vvxf') then; call vegas_integrate(LOdim_vegas,xsect_nnlo_vv_ns_xf,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'vvnf_full') then
         call ew_renorm_nf_1L()
         call ew_renorm_nf_2L()
         call vegas_integrate(LOdim_vegas,xsect_nnlo_vv_ns_nf_full,vg_result,vg_error,vg_chi2)
       else
         print *, 'wrong sector for ns, type -h for help'
         stop
       endif
    endif

  end subroutine run_proc_vv

  subroutine help_run_proc_vv(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -corr nnlo   -ch ns -sec vvnf vvfc vvxf'
    write(idev,*) ''
    
  end subroutine help_run_proc_vv
  
end module mod_run_proc_vv

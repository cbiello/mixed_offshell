module mod_run_proc_rv
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

  public :: run_proc_rv,help_run_proc_rv
  
contains

  subroutine run_proc_rv(vg_result,vg_error,vg_chi2)
    use mod_kinematics_nlo, only: NLOdim_vegas
    use mod_xsects_lo
    use mod_xsects_nnlo_rvewk
    use mod_xsects_nnlo_rvqcd
    use mod_xsects_checklim_gen
    use mod_amplitudes_loop_ppll, only: ew_renorm_nf_1L
    implicit none
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2

    if (checklim) then
       call vegas_integrate(NLOdim_vegas,xsect_checklim,vg_result,vg_error,vg_chi2)
    endif
        
    if (corr.eq.'nnlo' .and. ch.eq.'ns') then
       
       if (sec.eq.'rvewk_is') then
          call vegas_integrate(NLOdim_vegas,xsect_nnlo_rvewk_is_ns,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rvqcd_is') then
          call vegas_integrate(NLOdim_vegas,xsect_nnlo_rvqcd_is_ns,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rvqcd_fs_53') then
          call vegas_integrate(NLOdim_vegas,xsect_nnlo_rvqcd_fs_53_ns,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rvqcd_fs_54') then
          call vegas_integrate(NLOdim_vegas,xsect_nnlo_rvqcd_fs_54_ns,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rvewknf_is') then
          call ew_renorm_nf_1L()
          call vegas_integrate(NLOdim_vegas,xsect_nnlo_rvewknf_is_ns,vg_result,vg_error,vg_chi2)
       endif

    elseif (corr.eq.'nnlo' .and. ch.eq.'gq' .and. sec.eq.'rvewk_is') then
       call vegas_integrate(NLOdim_vegas,xsect_nnlo_rvewk_is_gq,vg_result,vg_error,vg_chi2)
    elseif (corr.eq.'nnlo' .and. ch.eq.'qg' .and. sec.eq.'rvewk_is') then
       call vegas_integrate(NLOdim_vegas,xsect_nnlo_rvewk_is_qg,vg_result,vg_error,vg_chi2)

    elseif (corr.eq.'nnlo' .and. ch.eq.'gq' .and. sec.eq.'rvewknf_is') then
       call ew_renorm_nf_1L()
       call vegas_integrate(NLOdim_vegas,xsect_nnlo_rvewknf_is_gq,vg_result,vg_error,vg_chi2)
    elseif (corr.eq.'nnlo' .and. ch.eq.'qg' .and. sec.eq.'rvewknf_is') then
       call ew_renorm_nf_1L()
       call vegas_integrate(NLOdim_vegas,xsect_nnlo_rvewknf_is_qg,vg_result,vg_error,vg_chi2)

    elseif (corr.eq.'nnlo' .and. ch.eq.'aq' .and. sec.eq.'rvqcd_is') then
       call vegas_integrate(NLOdim_vegas,xsect_nnlo_rvqcd_is_aq,vg_result,vg_error,vg_chi2)
    elseif (corr.eq.'nnlo' .and. ch.eq.'qa' .and. sec.eq.'rvqcd_is') then
       call vegas_integrate(NLOdim_vegas,xsect_nnlo_rvqcd_is_qa,vg_result,vg_error,vg_chi2)

    else
          print *, 'wrong sector, type -h for help'
    endif
       
  end subroutine run_proc_rv

  subroutine help_run_proc_rv(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -corr nnlo   -ch ns -sec rvewk_is'
    write(idev,*) ' -corr nnlo   -ch ns -sec rvqcd_is rvqcd_fs_5[3,4]'
    write(idev,*) ''
    write(idev,*) ' -corr nnlo   -ch gq -sec rvewk_is'
    write(idev,*) ' -corr nnlo   -ch qg -sec rvewk_is'
    write(idev,*) ''
    write(idev,*) ' -corr nnlo   -ch aq -sec rvqcd_is'
    write(idev,*) ' -corr nnlo   -ch qa -sec rvqcd_is'
    write(idev,*) ''
    
  end subroutine help_run_proc_rv
  
end module mod_run_proc_rv


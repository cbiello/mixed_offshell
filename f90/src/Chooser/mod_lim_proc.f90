module mod_lim_proc
  use mod_types
  use mod_consts_dp
  use mod_my_vegas
  use mod_proc_parms
  use mod_lim_proc_nlo
  use mod_check_lim
  
  implicit none
  private

  public :: lim_proc
  
contains

  subroutine lim_proc()
    use mod_kinematics_nlo, only: NLOdim_vegas
    use mod_kinematics_nnlo_tc_is_eu_ac, only: NNLOdim_vegas
    use mod_xsects_checklim_gen
    real(dp15)  :: vg_result,vg_error,vg_chi2

!    vg_result = 0
!    vg_error = 0
!    vg_chi2 = -1

!!    if (checklim) then
!!       if(corr(1:3).eq.'nlo')                    call vegas_integrate(NLOdim_vegas,xsect_checklim,vg_result,vg_error,vg_chi2)
!!       if(corr.eq.'nnlo' .and. sec(1:2).eq.'rr') call vegas_integrate(NNLOdim_vegas,xsect_checklim,vg_result,vg_error,vg_chi2)
!!    endif

    
    if (corr(1:3).eq.'nlo') then
       call run_check_lim_nlo(vg_result,vg_error,vg_chi2)
    elseif (corr.eq.'nnlo') then
       call vegas_integrate(NNLOdim_vegas,check_lim_nnlo,vg_result,vg_error,vg_chi2)
!       if (ch(1:2).eq.'ns' .and. sec(1:2) .eq. 'rr') call run_proc_rr_ns(vg_result,vg_error,vg_chi2)
!       if (ch.eq.'gq' .and. sec(1:2) .eq. 'rr') call run_proc_rr_gq(vg_result,vg_error,vg_chi2)
!       if (ch.eq.'qg' .and. sec(1:2) .eq. 'rr') call run_proc_rr_qg(vg_result,vg_error,vg_chi2)
!       if (ch.eq.'aq' .and. sec(1:2) .eq. 'rr') call run_proc_rr_aq(vg_result,vg_error,vg_chi2)
!       if (ch.eq.'qa' .and. sec(1:2) .eq. 'rr') call run_proc_rr_qa(vg_result,vg_error,vg_chi2)
!       if (ch.eq.'ag' .and. sec(1:2) .eq. 'rr') call run_proc_rr_ag(vg_result,vg_error,vg_chi2)
!       if (ch.eq.'ga' .and. sec(1:2) .eq. 'rr') call run_proc_rr_ga(vg_result,vg_error,vg_chi2)
!       if (sec(1:2) .eq. 'rv') call run_proc_rv(vg_result,vg_error,vg_chi2)
!       if (sec(1:2) .eq. 'vv') call run_proc_vv(vg_result,vg_error,vg_chi2)
!       !--
!       if (ch.eq.'ns' .and. sec(1:2) .eq. 's_')   call run_proc_s_ns(vg_result,vg_error,vg_chi2)
!       if (ch.eq.'gq' .and. sec(1:2) .eq. 's_')   call run_proc_s_gq(vg_result,vg_error,vg_chi2)
!       if (ch.eq.'qg' .and. sec(1:2) .eq. 's_')   call run_proc_s_qg(vg_result,vg_error,vg_chi2)
!       if (ch.eq.'aq' .and. sec(1:2) .eq. 's_')   call run_proc_s_aq(vg_result,vg_error,vg_chi2)
!       if (ch.eq.'qa' .and. sec(1:2) .eq. 's_')   call run_proc_s_qa(vg_result,vg_error,vg_chi2)
!       if (ch.eq.'ag' .and. sec(1:2) .eq. 's_')   call run_proc_s_ag(vg_result,vg_error,vg_chi2)
!       if (ch.eq.'ga' .and. sec(1:2) .eq. 's_')   call run_proc_s_ga(vg_result,vg_error,vg_chi2)
    else
       print *, 'wrong sector, type -h for help'
       stop
    endif
    
  end subroutine lim_proc

end module mod_lim_proc

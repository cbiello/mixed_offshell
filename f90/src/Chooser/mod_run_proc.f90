module mod_run_proc
  use mod_types
  use mod_consts_dp
  use mod_my_vegas
  use mod_proc_parms
  use mod_run_proc_lonlo
  use mod_run_proc_rr_ns
  use mod_run_proc_rr_gq
  use mod_run_proc_rr_qg
  use mod_run_proc_rr_aq
  use mod_run_proc_rr_qa
  use mod_run_proc_rr_ag
  use mod_run_proc_rr_ga
  use mod_run_proc_rv
  use mod_run_proc_vv
  !--
  use mod_run_proc_s_ns
  use mod_run_proc_s_gq
  use mod_run_proc_s_qg
  use mod_run_proc_s_aq
  use mod_run_proc_s_qa
  use mod_run_proc_s_ag
  use mod_run_proc_s_ga
  !--
  implicit none
  private

  public :: run_proc,help_run_proc
  
contains

  subroutine run_proc(vg_result,vg_error,vg_chi2)
    use mod_kinematics_nlo, only: NLOdim_vegas
    use mod_kinematics_nnlo_tc_is_eu_ac, only: NNLOdim_vegas
    use mod_xsects_checklim_gen
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2

    vg_result = 0
    vg_error = 0
    vg_chi2 = -1

    if (checklim) then
       if(corr(1:3).eq.'nlo')                    call vegas_integrate(NLOdim_vegas,xsect_checklim,vg_result,vg_error,vg_chi2)
       if(corr.eq.'nnlo' .and. sec(1:2).eq.'rr') call vegas_integrate(NNLOdim_vegas,xsect_checklim,vg_result,vg_error,vg_chi2)
    endif

    
    if (corr.eq.'lo') then
       call run_proc_lonlo(vg_result,vg_error,vg_chi2)
    elseif (corr(1:3).eq.'nlo') then
       call run_proc_lonlo(vg_result,vg_error,vg_chi2)
    elseif (corr.eq.'nnlo') then
       if (ch(1:2).eq.'ns' .and. sec(1:2) .eq. 'rr') call run_proc_rr_ns(vg_result,vg_error,vg_chi2)
       if (ch.eq.'gq' .and. sec(1:2) .eq. 'rr') call run_proc_rr_gq(vg_result,vg_error,vg_chi2)
       if (ch.eq.'qg' .and. sec(1:2) .eq. 'rr') call run_proc_rr_qg(vg_result,vg_error,vg_chi2)
       if (ch.eq.'aq' .and. sec(1:2) .eq. 'rr') call run_proc_rr_aq(vg_result,vg_error,vg_chi2)
       if (ch.eq.'qa' .and. sec(1:2) .eq. 'rr') call run_proc_rr_qa(vg_result,vg_error,vg_chi2)
       if (ch.eq.'ag' .and. sec(1:2) .eq. 'rr') call run_proc_rr_ag(vg_result,vg_error,vg_chi2)
       if (ch.eq.'ga' .and. sec(1:2) .eq. 'rr') call run_proc_rr_ga(vg_result,vg_error,vg_chi2)
       if (sec(1:2) .eq. 'rv') call run_proc_rv(vg_result,vg_error,vg_chi2)
       if (sec(1:2) .eq. 'vv') call run_proc_vv(vg_result,vg_error,vg_chi2)
       !--
       if (ch.eq.'ns' .and. sec(1:2) .eq. 's_')   call run_proc_s_ns(vg_result,vg_error,vg_chi2)
       if (ch.eq.'gq' .and. sec(1:2) .eq. 's_')   call run_proc_s_gq(vg_result,vg_error,vg_chi2)
       if (ch.eq.'qg' .and. sec(1:2) .eq. 's_')   call run_proc_s_qg(vg_result,vg_error,vg_chi2)
       if (ch.eq.'aq' .and. sec(1:2) .eq. 's_')   call run_proc_s_aq(vg_result,vg_error,vg_chi2)
       if (ch.eq.'qa' .and. sec(1:2) .eq. 's_')   call run_proc_s_qa(vg_result,vg_error,vg_chi2)
       if (ch.eq.'ag' .and. sec(1:2) .eq. 's_')   call run_proc_s_ag(vg_result,vg_error,vg_chi2)
       if (ch.eq.'ga' .and. sec(1:2) .eq. 's_')   call run_proc_s_ga(vg_result,vg_error,vg_chi2)
    else
       print *, 'wrong sector, type -h for help'
       stop
    endif
    
  end subroutine run_proc

  subroutine help_run_proc(idev)
    integer, intent(in) :: idev

    write(idev,*) ''
    write(idev,*) ' -corr lo/nlo/nnlo'
    write(idev,*) ''
    call help_run_proc_lonlo(idev)
    call help_run_proc_rr_ns(idev)
    call help_run_proc_rr_gq(idev)
    call help_run_proc_rr_qg(idev)
    call help_run_proc_rr_aq(idev)
    call help_run_proc_rr_qa(idev)
    call help_run_proc_rr_ag(idev)
    call help_run_proc_rr_ga(idev)
    call help_run_proc_rv(idev)
    call help_run_proc_vv(idev)
    !--
    call help_run_proc_s_ns(idev)
    call help_run_proc_s_gq(idev)
    call help_run_proc_s_qg(idev)
    call help_run_proc_s_aq(idev)
    call help_run_proc_s_qa(idev)
    call help_run_proc_s_ag(idev)
    call help_run_proc_s_ga(idev)
  end subroutine help_run_proc
  
end module mod_run_proc

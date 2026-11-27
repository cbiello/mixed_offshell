module mod_tc_proc
  use mod_types; use common_def; use mod_consts_dp
  use mod_vegas
  use mod_tripleclimits_q_gg_qstar
  use mod_tripleclimits_q_qqb_qstar
  use mod_tripleclimits_g_qg_qstar
  
  implicit none
  private

  public :: run_tc_proc

contains
  
  subroutine run_tc_proc(corr,sec,ch,gridfile,vg_result,vg_error,vg_chi2)
    character, intent(in) :: corr*(8),sec*(10),ch*(10),gridfile*(100)
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2
    integer :: ndim

    ndim=4
    if (corr .eq. 'tcsplit') then
       if (ch .eq. "q_gg_q") then
          call go_vegas_gen(ndim,seed,tripcoll_subtr_q_gg_qstar,vg_result,vg_error,vg_chi2,gridfile)
       elseif (ch .eq. "g_qg_q") then
          call go_vegas_gen(ndim,seed,tripcoll_subtr_g_qg_qstar,vg_result,vg_error,vg_chi2,gridfile)
       elseif (ch .eq. "q_qq_q") then
          call go_vegas_gen(ndim,seed,tripcoll_subtr_q_qqb_qstar,vg_result,vg_error,vg_chi2,gridfile)

       endif
    endif

  end subroutine run_tc_proc

end module mod_tc_proc

  

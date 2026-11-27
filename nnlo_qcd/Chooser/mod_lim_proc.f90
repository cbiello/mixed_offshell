module mod_lim_proc
  use mod_types; use common_def; use mod_consts_dp
  use mod_consts_debug
  use mod_vegas
  !--
  use mod_check_lim_nlo
  use mod_check_lim_rv
  use mod_check_lim_rr_wp_ns
  use mod_check_lim_rr_wp_qg
  use mod_check_lim_rr_wp_gq
  use mod_check_lim_rr_wp_s
  use mod_check_lim_rr_wp_gg
  use mod_check_lim_sub
  !--
  implicit none
  private

  public :: run_limits_proc

contains

  subroutine run_limits_proc(corr,sec,ch,gridfile,vg_result,vg_error,vg_chi2)
    character, intent(in) :: corr*(8),sec*(10),ch*(10),gridfile*(100)
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2
    integer :: ndim

    if (corr(1:3) .eq. 'lim') then

       print *, '#---------------------------------------'
       print *, '# limcheck: ', limcheck
       print *, '#---------------------------------------'
       
    endif

    !---

    if (corr .eq. 'limnlo' .and. sec .eq. 'r' ) then

       ndim = LOdim + 3
       call go_vegas_gen(ndim,seed,xsect_nlo_checklim_r,vg_result,vg_error,vg_chi2,gridfile)

    endif

    if (corr .eq. 'limnnlo') then


       ! subtraction ~ FLM[1,2,4]
       if (sec .eq. 'sub124') then
       
          ndim = LOdim+4
          if (ch .eq. 'ns') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_sub124_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'gq') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_sub124_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'qg') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_sub124_wp_qg,vg_result,vg_error,vg_chi2,gridfile)

       !-- RV
       elseif (sec .eq. 'rv') then

          ndim = LOdim + 3
          if (ch .eq. 'ns') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rv_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'qg') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rv_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'gq') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rv_wp_gq,vg_result,vg_error,vg_chi2,gridfile)

       endif

       !-- RR, ns

       if (ch .eq. 'ns') then

          globalsec = sec
          ndim = LOdim + 6

          if (sec .eq. 'rr4151a') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151b') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151c') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151d') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4152') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_dc_wp_ns,vg_result,vg_error,vg_chi2,gridfile)

          if (sec .eq. 'rr4252a') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252b') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252c') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252d') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4251') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_dc_wp_ns,vg_result,vg_error,vg_chi2,gridfile)

       elseif (ch .eq. 'qg') then

          globalsec = sec
          ndim = LOdim + 6

          if (sec .eq. 'rr4252a') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252b') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252c') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252d') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4251') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_dc_wp_qg,vg_result,vg_error,vg_chi2,gridfile)

          if (sec .eq. 'rr4151a') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151b') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151c') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151d') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4152') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_dc_wp_qg,vg_result,vg_error,vg_chi2,gridfile)

       elseif (ch .eq. 'gq') then

          globalsec = sec
          ndim = LOdim + 6

          if (sec .eq. 'rr4252a') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252b') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252c') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252d') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4251') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_dc_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
!
          if (sec .eq. 'rr4151a') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151b') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151c') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151d') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4152') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_dc_wp_gq,vg_result,vg_error,vg_chi2,gridfile)

       elseif (ch .eq. 's') then

          globalsec = sec
          ndim = LOdim + 6

          if (sec .eq. 'rr4252a') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252b') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252c') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252d') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4251') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_dc_wp_s,vg_result,vg_error,vg_chi2,gridfile)

          if (sec .eq. 'rr4151a') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151b') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151c') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151d') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_tc_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4152') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_dc_wp_s,vg_result,vg_error,vg_chi2,gridfile)

       elseif (ch .eq. 'gg') then

          globalsec = sec
          ndim = LOdim + 6

          if (sec .eq. 'rr') call go_vegas_gen(ndim,seed,xsect_nnlo_checklim_rr_wp_gg,vg_result,vg_error,vg_chi2,gridfile)

       endif

    endif
       
  end subroutine run_limits_proc

end module mod_lim_proc

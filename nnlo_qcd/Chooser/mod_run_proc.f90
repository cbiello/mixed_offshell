module mod_run_proc
  use mod_types; use common_def; use mod_consts_dp
  use mod_vegas
  !--
  use mod_xsects_lo
  use mod_xsects_nlo
  !--
  use mod_xsects_nnlo_vv_wp_ns
  !--
  use mod_xsects_nnlo_rv_wp_ns
  use mod_xsects_nnlo_rv_wp_gq
  use mod_xsects_nnlo_rv_wp_qg
  !--
  use mod_xsects_nnlo_rr_wp_ns
  use mod_xsects_nnlo_rr_wp_qg
  use mod_xsects_nnlo_rr_wp_gq
  use mod_xsects_nnlo_rr_wp_s
  use mod_xsects_nnlo_rr_wp_gg
  !--
  use mod_xsects_nnlo_sub_wp_ns
  use mod_xsects_nnlo_sub_wp_gq
  use mod_xsects_nnlo_sub_wp_qg
  use mod_xsects_nnlo_sub_wp_gg
  implicit none
  private

  public :: run_sector_proc

contains

  subroutine run_sector_proc(corr,sec,ch,gridfile,vg_result,vg_error,vg_chi2)
    character, intent(in) :: corr*(8),sec*(10),ch*(10),gridfile*(100)
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2
    integer :: ndim

    if (corr .eq. 'lo') then

       ndim = LOdim
       call go_vegas_gen(ndim,seed,xsect_lo,vg_result,vg_error,vg_chi2,gridfile)
       
    endif

    if (corr .eq. 'nlo') then

       if (sec.eq.'r') then
          ndim = LOdim + 3
          ndim = ndim + 1
          call go_vegas_gen(ndim,seed,xsect_nlo_r,vg_result,vg_error,vg_chi2,gridfile)
       endif

       if (sec.eq.'v') then
          ndim = LOdim
          call go_vegas_gen(ndim,seed,xsect_nlo_v,vg_result,vg_error,vg_chi2,gridfile)
       endif

       if (sec.eq.'s') then
          ndim = LOdim + 1
          call go_vegas_gen(ndim,seed,xsect_nlo_s,vg_result,vg_error,vg_chi2,gridfile)
       endif
       
    endif

    !-- NNLO

    if (corr .eq. 'nnlo') then

       !-- VV
       if (sec .eq. 'vv') then
          
          ndim = LOdim

          if (ch .eq. 'ns') call go_vegas_gen(ndim,seed,xsect_nnlo_vv_wp_ns,vg_result,vg_error,vg_chi2,gridfile)

       endif

       !-- RV
       if (sec .eq. 'rv') then

          ndim = LOdim + 3
          if (ch .eq. 'ns') call go_vegas_gen(ndim,seed,xsect_nnlo_rv_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'qg') call go_vegas_gen(ndim,seed,xsect_nnlo_rv_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'gq') call go_vegas_gen(ndim,seed,xsect_nnlo_rv_wp_gq,vg_result,vg_error,vg_chi2,gridfile)

       ! -- VV
       elseif (sec .eq. 'vv') then
!          if (ch .eq. 'gg') call go_vegas_gen(ndim,seed,xsect_nnlo_vv_h,vg_result,vg_error,vg_chi2,gridfile)          

       ! -- subtraction  ~ FLM[1,2]
       elseif (sec .eq. 'sub12') then
          ndim = LOdim + 2    ! z,zb
          ndim = ndim + 3     ! for triple collinear integration
          ndim = ndim + 1     ! max
          if (ch .eq. 'ns') call go_vegas_gen(ndim,seed,xsect_nnlo_sub12_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'gq') call go_vegas_gen(ndim,seed,xsect_nnlo_sub12_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'qg') call go_vegas_gen(ndim,seed,xsect_nnlo_sub12_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'gg') call go_vegas_gen(ndim,seed,xsect_nnlo_sub12_wp_gg,vg_result,vg_error,vg_chi2,gridfile)

       ! --  subtraction ~ FLVfin[1,2]
       elseif (sec .eq. 'subv') then
          ndim = LOdim + 1    
          if (ch .eq. 'ns') call go_vegas_gen(ndim,seed,xsect_nnlo_subv_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'gq') call go_vegas_gen(ndim,seed,xsect_nnlo_subv_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'qg') call go_vegas_gen(ndim,seed,xsect_nnlo_subv_wp_qg,vg_result,vg_error,vg_chi2,gridfile)

       ! --  subtraction ~ FLM[1,2,4]
       elseif (sec .eq. 'sub124') then
          ndim = LOdim + 3 + 1
          if (ch .eq. 'ns') call go_vegas_gen(ndim,seed,xsect_nnlo_sub124_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'gq') call go_vegas_gen(ndim,seed,xsect_nnlo_sub124_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'qg') call go_vegas_gen(ndim,seed,xsect_nnlo_sub124_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (ch .eq. 'gg') call go_vegas_gen(ndim,seed,xsect_nnlo_sub124_wp_gg,vg_result,vg_error,vg_chi2,gridfile)

       endif

       !-- RR, gg

       if (ch .eq. 'ns') then

          ndim = LOdim + 6

          if (sec .eq. 'rr4151a') call go_vegas_gen(ndim,seed,xsect_nnlo_4151a_rr_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151b') call go_vegas_gen(ndim,seed,xsect_nnlo_4151b_rr_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151c') call go_vegas_gen(ndim,seed,xsect_nnlo_4151c_rr_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151d') call go_vegas_gen(ndim,seed,xsect_nnlo_4151d_rr_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4152') call go_vegas_gen(ndim,seed,xsect_nnlo_4152_rr_wp_ns,vg_result,vg_error,vg_chi2,gridfile)

          if (sec .eq. 'rr4252a') call go_vegas_gen(ndim,seed,xsect_nnlo_4252a_rr_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252b') call go_vegas_gen(ndim,seed,xsect_nnlo_4252b_rr_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252c') call go_vegas_gen(ndim,seed,xsect_nnlo_4252c_rr_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252d') call go_vegas_gen(ndim,seed,xsect_nnlo_4252d_rr_wp_ns,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4251') call go_vegas_gen(ndim,seed,xsect_nnlo_4251_rr_wp_ns,vg_result,vg_error,vg_chi2,gridfile)

       elseif (ch .eq. 'qg') then

          ndim = LOdim + 6

          if (sec .eq. 'rr4252a') call go_vegas_gen(ndim,seed,xsect_nnlo_4252a_rr_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252b') call go_vegas_gen(ndim,seed,xsect_nnlo_4252b_rr_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252c') call go_vegas_gen(ndim,seed,xsect_nnlo_4252c_rr_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252d') call go_vegas_gen(ndim,seed,xsect_nnlo_4252d_rr_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4251') call go_vegas_gen(ndim,seed,xsect_nnlo_4251_rr_wp_qg,vg_result,vg_error,vg_chi2,gridfile)

          if (sec .eq. 'rr4151a') call go_vegas_gen(ndim,seed,xsect_nnlo_4151a_rr_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151b') call go_vegas_gen(ndim,seed,xsect_nnlo_4151b_rr_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151c') call go_vegas_gen(ndim,seed,xsect_nnlo_4151c_rr_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151d') call go_vegas_gen(ndim,seed,xsect_nnlo_4151d_rr_wp_qg,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4152') call go_vegas_gen(ndim,seed,xsect_nnlo_4152_rr_wp_qg,vg_result,vg_error,vg_chi2,gridfile)

       elseif (ch .eq. 'gq') then

          ndim = LOdim + 6

          if (sec .eq. 'rr4252a') call go_vegas_gen(ndim,seed,xsect_nnlo_4252a_rr_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252b') call go_vegas_gen(ndim,seed,xsect_nnlo_4252b_rr_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252c') call go_vegas_gen(ndim,seed,xsect_nnlo_4252c_rr_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252d') call go_vegas_gen(ndim,seed,xsect_nnlo_4252d_rr_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4251') call go_vegas_gen(ndim,seed,xsect_nnlo_4251_rr_wp_gq,vg_result,vg_error,vg_chi2,gridfile)

          if (sec .eq. 'rr4151a') call go_vegas_gen(ndim,seed,xsect_nnlo_4151a_rr_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151b') call go_vegas_gen(ndim,seed,xsect_nnlo_4151b_rr_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151c') call go_vegas_gen(ndim,seed,xsect_nnlo_4151c_rr_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151d') call go_vegas_gen(ndim,seed,xsect_nnlo_4151d_rr_wp_gq,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4152') call go_vegas_gen(ndim,seed,xsect_nnlo_4152_rr_wp_gq,vg_result,vg_error,vg_chi2,gridfile)

       elseif (ch .eq. 's') then

          ndim = LOdim + 6

          if (sec .eq. 'rr4252a') call go_vegas_gen(ndim,seed,xsect_nnlo_4252a_rr_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252b') call go_vegas_gen(ndim,seed,xsect_nnlo_4252b_rr_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252c') call go_vegas_gen(ndim,seed,xsect_nnlo_4252c_rr_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4252d') call go_vegas_gen(ndim,seed,xsect_nnlo_4252d_rr_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4251') call go_vegas_gen(ndim,seed,xsect_nnlo_4251_rr_wp_s,vg_result,vg_error,vg_chi2,gridfile)

          if (sec .eq. 'rr4151a') call go_vegas_gen(ndim,seed,xsect_nnlo_4151a_rr_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151b') call go_vegas_gen(ndim,seed,xsect_nnlo_4151b_rr_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151c') call go_vegas_gen(ndim,seed,xsect_nnlo_4151c_rr_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4151d') call go_vegas_gen(ndim,seed,xsect_nnlo_4151d_rr_wp_s,vg_result,vg_error,vg_chi2,gridfile)
          if (sec .eq. 'rr4152') call go_vegas_gen(ndim,seed,xsect_nnlo_4152_rr_wp_s,vg_result,vg_error,vg_chi2,gridfile)

       elseif (ch .eq. 'gg') then

          ndim = LOdim + 6

          if (sec .eq. 'rr') call go_vegas_gen(ndim,seed,xsect_nnlo_rr_wp_gg,vg_result,vg_error,vg_chi2,gridfile)

       endif

    endif

  end subroutine run_sector_proc


end module mod_run_proc

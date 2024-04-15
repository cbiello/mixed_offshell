module mod_run_proc_sub_fc
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

  public :: run_proc_sub_fc,help_run_proc_sub_fc
  
contains

  subroutine run_proc_sub_fc(vg_result,vg_error,vg_chi2)
    use mod_kinematics_lo_z_fc
    use mod_kinematics_nlo_z_fc
    use mod_kinematics_nnlo_tc_is_eu_abcd
    use mod_xsects_nnlo_sub12_ns_fc
    use mod_xsects_nnlo_subonloqcd_ns_fc
    use mod_xsects_nnlo_subonloewk_ns_fc
    use mod_xsects_nnlo_subv_ns_fc
    !
    use mod_xsects_nnlo_subonloqcd_aq_fc
    use mod_xsects_nnlo_subonloewk_aq_fc
    use mod_xsects_nnlo_sub12_aq_fc
    use mod_xsects_nnlo_subv_aq_fc
    !
    use mod_xsects_nnlo_subonloqcd_ag_fc
    use mod_xsects_nnlo_sub12_ag_fc
    use mod_xsects_nnlo_rr_ag_fc
    use mod_xsects_nnlo_rr_ga_fc
    !
    use mod_xsects_nnlo_subonloqcd_gq_fc
    !
    use mod_hoppet_nlo
    use mod_hoppet_nnlo_fc
    use mod_ol_interface
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2

    if (ch.eq.'ns') then
       if (sec.eq.'fcsub12') then
          call vegas_integrate(LOdim_vegas+2,xsect_nnlo_sub12_ns_ga_fc,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'fcsub12test') then
          call vegas_integrate(LOdim_vegas+2,xsect_nnlo_sub12_test,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'fcsub_onloqcd') then
          call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_subonloqcd_ns_ga_fc,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'fcsub_onloewk_is') then
          call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_subonloewk_is_ns_ga_fc,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'fcsub_onloewk_53') then
          call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_subonloewk_fs_53_ns_ga_fc,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'fcsub_onloewk_54') then
          call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_subonloewk_fs_54_ns_ga_fc,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'fcsubv_qcd') then
          call init_xPij_nlo()
          call vegas_integrate(LOdim_vegas,xsect_nnlo_subvqcd_ns_fc,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'fcsubv_ewk') then
          sec = 'subv_ewk'
          call initialise_ol(0)
          sec = 'fcsubv_ewk'
          call init_xPij_nlo()
          call vegas_integrate(LOdim_vegas,xsect_nnlo_subvewk_ns_fc,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'fcsub12_qqb') then
          call init_xPij_nnlo_fc()
          call vegas_integrate(LOdim_vegas,xsect_nnlo_sub12_ns_qqbORqq_fc,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'fcsub12_qq') then
          call init_xPij_nnlo_fc()
          call vegas_integrate(LOdim_vegas,xsect_nnlo_sub12_ns_qqbORqq_fc,vg_result,vg_error,vg_chi2)       
       else
          print *, 'wrong sector, type -h for help'
          stop
       endif
    endif

    if (ch.eq.'aq') then
       if (sec.eq.'fcsub_onloqcd') then
          call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_subonloqcd_aq_fc,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'fcsub_onloewk') then
          call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_subonloewk_is_aq_fc,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'fcsub12') then
          call vegas_integrate(LOdim_vegas+2,xsect_nnlo_sub12_aq_fc,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'fcsubv_qcd') then
          call init_xPij_nlo()
          call vegas_integrate(LOdim_vegas,xsect_nnlo_subvqcd_aq_fc,vg_result,vg_error,vg_chi2)
      else
          print *, 'wrong sector, type -h for help'
          stop
       endif
    endif

    if (ch.eq.'ag') then
       if (sec.eq.'fcsub_onloqcd') then
          !call init_xPij_nnlo_fc()
          !call vegas_integrate(NLOdim_vegas,xsect_nnlo_subonloqcd_ag_fc,vg_result,vg_error,vg_chi2) !-- hoppet
          call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_subonloqcd_ag_fc_nohoppet,vg_result,vg_error,vg_chi2) !-- no hoppet
       elseif (sec.eq.'fcsub_onloqcd_split') then
          call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_subonloqcd_ag_fc_nohoppet_split,vg_result,vg_error,vg_chi2) !-- no hoppet
       elseif (sec.eq.'fcsub12') then
          call vegas_integrate(LOdim_vegas+2,xsect_nnlo_sub12_ag_fc,vg_result,vg_error,vg_chi2)
       else if (sec.eq.'fcrr_5262a') then
          call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262a_ag_fc,vg_result,vg_error,vg_chi2)
       else if (sec.eq.'fcrr_5262b') then
          call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262b_ag_fc,vg_result,vg_error,vg_chi2)
       else if (sec.eq.'fcrr_5262c') then
          call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262c_ag_fc,vg_result,vg_error,vg_chi2)
       else if (sec.eq.'fcrr_5262d') then
          call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262d_ag_fc,vg_result,vg_error,vg_chi2)
       else if (sec.eq.'fcrr_5161') then !-- as in the paper, no splitting
          call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161_ag_fc,vg_result,vg_error,vg_chi2)
       else
          print *, 'wrong sector, type -h for help'
          stop
       endif
    endif

    if (ch.eq.'ga') then
       if (sec.eq.'fcrr_5161a') then
          call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161a_ga_fc,vg_result,vg_error,vg_chi2)
       else if (sec.eq.'fcrr_5161b') then
          call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161b_ga_fc,vg_result,vg_error,vg_chi2)
       else if (sec.eq.'fcrr_5161c') then
          call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161c_ga_fc,vg_result,vg_error,vg_chi2)
       else if (sec.eq.'fcrr_5161d') then
          call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161d_ga_fc,vg_result,vg_error,vg_chi2)
       else
          print *, 'wrong sector, type -h for help'
          stop
       endif
    endif

    if (ch.eq.'gq') then
       if (sec.eq.'fcsub_onloqcd') then
          call vegas_integrate(NLOdim_vegas+1,xsect_nnlo_subonloqcd_gq_fc,vg_result,vg_error,vg_chi2)
       else
          print *, 'wrong sector, type -h for help'
          stop
       endif
    endif

  end subroutine run_proc_sub_fc

  subroutine help_run_proc_sub_fc(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -ch ns -sec fcsub12 fcsub_onloqcd fcsub_onloewk_is fcsub_onloewk_53 fcsub_onloewk_54 fcsubv_ewk fcsubv_qcd fcsub12_qqb fcsub12_qq'
    write(idev,*) ''
    
  end subroutine help_run_proc_sub_fc
  
end module mod_run_proc_sub_fc


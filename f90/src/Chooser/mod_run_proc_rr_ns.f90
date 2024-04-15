module mod_run_proc_rr_ns
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

  public :: run_proc_rr_ns,help_run_proc_rr_ns
  
contains

  subroutine run_proc_rr_ns(vg_result,vg_error,vg_chi2)
    use mod_kinematics_nnlo_tc_is_eu_ac, only: NNLOdim_vegas
    use mod_xsects_nnlo_rr_ns
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2


    if (ch.eq.'ns_ga') then
       if     (sec.eq.'rr_5161a') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161a_ns_ga,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5161c') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161c_ns_ga,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5262a') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262a_ns_ga,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5262c') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262c_ns_ga,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5162') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5162_ns_ga,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5261') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5261_ns_ga,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5163') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5163_ns_ga,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5263') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5263_ns_ga,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5164') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5164_ns_ga,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5264') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5264_ns_ga,vg_result,vg_error,vg_chi2)
       else
          print *, 'wrong sector for ns_ga, type -h for help'
          stop
       endif
    endif

    if (ch.eq.'ns_qqb') then
       if     (sec.eq.'rr_5161a') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161a_ns_qqb,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5161c') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161c_ns_qqb,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5262a') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262a_ns_qqb,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5262c') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262c_ns_qqb,vg_result,vg_error,vg_chi2)
       else
          print *, 'wrong sector for ns_4q_qqb, type -h for help'
          stop
       endif
    endif

    if (ch.eq.'ns_qq') then
       if     (sec.eq.'rr_5161a') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161a_ns_qq,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5161c') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5161c_ns_qq,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5262a') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262a_ns_qq,vg_result,vg_error,vg_chi2)
       elseif (sec.eq.'rr_5262c') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_5262c_ns_qq,vg_result,vg_error,vg_chi2)
       else
          print *, 'wrong sector for ns_4q_qq, type -h for help'
          stop
       endif
    endif

    if (ch.eq.'ns_qqb_w') then
       if     (sec.eq.'rr') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_ns_qqb_w,vg_result,vg_error,vg_chi2)
       else
          print *, 'wrong sector for ns_4q_qqb_W, type -h for help'
          stop
       endif
    endif

    if (ch.eq.'ns_qqp_w') then
       if     (sec.eq.'rr') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_ns_qqp_w,vg_result,vg_error,vg_chi2)
       else
          print *, 'wrong sector for ns_4q_qqp_W, type -h for help'
          stop
       endif
    endif

    if (ch.eq.'ns_qqpb_w') then
       if     (sec.eq.'rr') then; call vegas_integrate(NNLOdim_vegas,xsect_nnlo_rr_ns_qqpb_w,vg_result,vg_error,vg_chi2)
       else
          print *, 'wrong sector for ns_4q_qqp_W, type -h for help'
          stop
       endif
    endif

  end subroutine run_proc_rr_ns

  subroutine help_run_proc_rr_ns(idev)
    integer, intent(in) :: idev

    write(idev,*) ' -corr nnlo   -ch ns_ga     -sec rr_5161[a,c] rr_5262[a,c] rr_516[2,3,4] rr_526[1,3,4] #qqb / qbq -> ga channels'
    write(idev,*) ' -corr nnlo   -ch ns_qqb    -sec rr_5161[a,c] rr_5262[a,c]                             #qqb -> qqb / qb q  -> q qb  channels'
    write(idev,*) ' -corr nnlo   -ch ns_qq     -sec rr_5161[a,c] rr_5262[a,c]                             #qq  -> qq  / qb qb -> qb qb channels'
    write(idev,*) ' -corr nnlo   -ch ns_qqb_w  -sec rr'
    write(idev,*) ' -corr nnlo   -ch ns_qqp_w  -sec rr'
    write(idev,*) ' -corr nnlo   -ch ns_qqpb_w -sec rr'
    write(idev,*) ''
    
  end subroutine help_run_proc_rr_ns
  
end module mod_run_proc_rr_ns


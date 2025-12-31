module mod_lim_proc_nlo
  use mod_types
  use mod_consts_dp
#if ( _MPIbuild )
  use mod_vegas_mpi
#else
  use mod_my_vegas
#endif
  use mod_proc_parms
  use mod_check_lim
  implicit none
  private

  public :: run_check_lim_nlo!,help_run_proc_lonlo
  
contains

  subroutine run_check_lim_nlo(vg_result,vg_error,vg_chi2)
    use mod_kinematics_lo
    use mod_kinematics_nlo
    use mod_xsects_lo
    use mod_xsects_nloqcd_r
    use mod_xsects_nloqcd_vs
    use mod_hoppet_nlo
    use mod_xsects_nloewk_r
    use mod_xsects_nloewk_vs
    use mod_xsects_checklim_gen
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2
!    real(dp15)  :: xx(30) 

    if (corr.eq.'nloqcd') then

       if (ch.eq.'ns') then
!          call random_number( xx(1:NLOdim_vegas) )
          if     (sec.eq.'r_is') then
             call vegas_integrate(NLOdim_vegas,check_lim_nlo,vg_result,vg_error,vg_chi2)
!             call vegas_integrate(NLOdim_vegas,xsect_nloqcd_r_is_ns,vg_result,vg_error,vg_chi2)
!          elseif     (sec.eq.'r_is_wp') then
!             call vegas_integrate(NLOdim_vegas,xsect_nloqcd_r_is_ns_wp,vg_result,vg_error,vg_chi2)
!          elseif     (sec.eq.'r_is_wm') then
!             call vegas_integrate(NLOdim_vegas,xsect_nloqcd_r_is_ns_wm,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'s')    then
!             call init_xPij_nlo()
!             call vegas_integrate(LOdim_vegas,xsect_nloqcd_s_ns,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'v')    then
!             call vegas_integrate(LOdim_vegas,xsect_nloqcd_v_ns,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'s_wp')    then
!             call init_xPij_nlo()
!             call vegas_integrate(LOdim_vegas,xsect_nloqcd_s_ns_wp,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'v_wp')    then
!             call vegas_integrate(LOdim_vegas,xsect_nloqcd_v_ns_wp,vg_result,vg_error,vg_chi2)
          else
             print *, 'wrong sector for nloqcd ns -> -sec r_is,s,v'
             stop
          endif
!       elseif (ch.eq.'gq') then
!          if     (sec.eq.'r_is') then
!             call vegas_integrate(NLOdim_vegas,xsect_nloqcd_r_is_gq,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'s')    then
!             call init_xPij_nlo()
!             call vegas_integrate(LOdim_vegas,xsect_nloqcd_s_gq,vg_result,vg_error,vg_chi2)
!          else
!             print *, 'wrong sector for nloqcd gq -> -sec r_is,s'
!             stop
!          endif
!       elseif (ch.eq.'qg') then
!          if     (sec.eq.'r_is') then
!             call vegas_integrate(NLOdim_vegas,xsect_nloqcd_r_is_qg,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'s')    then
!             call init_xPij_nlo()
!             call vegas_integrate(LOdim_vegas,xsect_nloqcd_s_qg,vg_result,vg_error,vg_chi2)
!          else
!             print *, 'wrong sector for nloqcd qg -> -sec r_is,s'
!             stop
!          endif
       else
          print *, 'wrong channel for nloqcd -> -ch ns gq qg'
          stop
       endif
       
    elseif (corr.eq.'nloewk') then

       if (ch.eq.'ns') then
          if     (sec.eq.'r_is') then
             call vegas_integrate(NLOdim_vegas,check_lim_nlo_is,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'r_is_wp') then
!             call vegas_integrate(NLOdim_vegas,xsect_nloewk_r_is_ns_wp,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'r_is_wm') then
!             call vegas_integrate(NLOdim_vegas,xsect_nloewk_r_is_ns_wm,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'r_fs_53') then
!             call vegas_integrate(NLOdim_vegas,xsect_nloewk_r_fs_53_ns,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'r_fs_54') then
!             call vegas_integrate(NLOdim_vegas,xsect_nloewk_r_fs_54_ns,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'s')    then
!             call init_xPij_nlo()
!             call vegas_integrate(LOdim_vegas,xsect_nloewk_s_ns,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'v')    then
!             call vegas_integrate(LOdim_vegas,xsect_nloewk_v_ns,vg_result,vg_error,vg_chi2)
          else
             print *, 'wrong sector for nloewk ns -> -sec r_is,r_fs_53,r_fs_54,s,v'
             stop
          endif
!       elseif (ch.eq.'aq') then
!          if     (sec.eq.'r_is') then
!             call vegas_integrate(NLOdim_vegas,xsect_nloewk_r_is_aq,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'s')    then
!             call init_xPij_nlo()
!             call vegas_integrate(LOdim_vegas,xsect_nloewk_s_aq,vg_result,vg_error,vg_chi2)
!          else
!             print *, 'wrong sector for nloewk aq -> -sec r_is,s'
!             stop
!          endif
!       elseif (ch.eq.'qa') then
!          if     (sec.eq.'r_is') then
!             call vegas_integrate(NLOdim_vegas,xsect_nloewk_r_is_qa,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'s')    then
!             call init_xPij_nlo()
!             call vegas_integrate(LOdim_vegas,xsect_nloewk_s_qa,vg_result,vg_error,vg_chi2)
!          else
!             print *, 'wrong sector for nloewk qa -> -sec r_is,s'
!             stop
!          endif
!       elseif (ch.eq.'aa') then
!          if (sec.eq.'v') then
!             call vegas_integrate(LOdim_vegas,xsect_nloewk_v_aa,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'r_fs_53') then
!             call vegas_integrate(NLOdim_vegas,xsect_nloewk_r_fs_53_aa,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'r_fs_54') then
!             call vegas_integrate(NLOdim_vegas,xsect_nloewk_r_fs_54_aa,vg_result,vg_error,vg_chi2)
!          elseif (sec.eq.'s')    then
!             call vegas_integrate(LOdim_vegas,xsect_nloewk_s_aa,vg_result,vg_error,vg_chi2)
!          endif
       else
          print *, 'wrong channel for nloqewk -> -ch ns aq qa aa'
          stop
       endif

    else
       print *, 'wrong sector for lo/nlo, type -h for help'
       stop
    endif
    
  end subroutine run_check_lim_nlo

!!  subroutine help_run_proc_lonlo(idev)
!!    integer, intent(in) :: idev
!!
!!    write(idev,*) 'Z: -corr lo -ch ns/aa -sec na'
!!    write(idev,*) ''
!!    write(idev,*) 'W+: -corr lo -ch ns -sec wp'
!!    write(idev,*) ''
!!    write(idev,*) 'Z: -corr nloqcd -ch ns -sec r_is s v'
!!    write(idev,*) '   -corr nloqcd -ch gq -sec r_is s'
!!    write(idev,*) '   -corr nloqcd -ch qg -sec r_is s'
!!    write(idev,*) ''
!!    write(idev,*) 'W+: -corr nloqcd -ch ns -sec r_is_wp s_wp v_wp'
!!    write(idev,*) ''
!!    write(idev,*) 'Z: -corr nloewk -ch ns -sec r_is r_fs_5[3,4] s v'
!!    write(idev,*) '   -corr nloewk -ch aq -sec r_is s'
!!    write(idev,*) '   -corr nloewk -ch qa -sec r_is s'
!!    write(idev,*) ''
!!    write(idev,*) 'Z: -corr nloewk -ch aa -sec r_fs_5[3,4] s v'
!!    write(idev,*) ''
!!    
!!  end subroutine help_run_proc_lonlo
  
end module mod_lim_proc_nlo


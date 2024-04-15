module mod_xsects_checklim_gen
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_auxfunctions
  use mod_xsects_checklim_nlo
  use mod_xsects_checklim_nnlo_ns
  implicit none

  private

  !-- nlor_g_ns    --> H,C1,C2
  !-- nlor_a_is_ns --> H,C1,C2,S,SC1,SC2
  !-- nlor_a_fs_ns --> H,C,S,SC
  !-- nlor_gq      --> H,C1
  !-- nlor_aq      --> H,C1,C2
  !-- nlor_qg      --> H,C2
  
  public :: xsect_checklim

contains

  function xsect_checklim(yRnd,ff,vegasweight)
    integer :: xsect_checklim
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_checklim = 0
    
    print*, 'limit depth: ',lim
    print *, ''

    !-- NLO
    if (corr.eq.'nloqcd') then
       if     (ch.eq.'ns' .and. sec.eq.'r_is') then; xsect_checklim = xsect_checklim_HC1C2(yRnd,ff,vegasweight)
       elseif (ch.eq.'gq' .and. sec.eq.'r_is') then; xsect_checklim = xsect_checklim_HC1(yRnd,ff,vegasweight)
       elseif (ch.eq.'qg' .and. sec.eq.'r_is') then; xsect_checklim = xsect_checklim_HC2(yRnd,ff,vegasweight)
       else
          call error('nloqcd')
       endif

    elseif (corr.eq.'nloewk') then
       if     (ch.eq.'ns' .and. sec.eq.'r_is')      then; xsect_checklim = xsect_checklim_nlo6(yRnd,ff,vegasweight)
       elseif (ch.eq.'ns' .and. sec(1:4).eq.'r_fs') then; xsect_checklim = xsect_checklim_nlo4(yRnd,ff,vegasweight)
       elseif (ch.eq.'aa' .and. sec(1:4).eq.'r_fs') then; xsect_checklim = xsect_checklim_nlo4(yRnd,ff,vegasweight)
       elseif (ch.eq.'aq' .and. sec.eq.'r_is')      then; xsect_checklim = xsect_checklim_HC1C2(yRnd,ff,vegasweight)
       elseif (ch.eq.'qa' .and. sec.eq.'r_is')      then; xsect_checklim = xsect_checklim_HC1C2(yRnd,ff,vegasweight)
       else
          call error('nloewk')
       endif
       
    endif

    !-- NNLO, RR
    if (corr.eq.'nnlo' .and. sec(1:2).eq.'rr') then
       if (ch.eq.'ns_qqb'.or.ch.eq.'ns_qq') then; xsect_checklim = xsect_checklim_nnlo4q(yRnd,ff,vegasweight)
       endif
    endif

    !-- NNLO, RV
    if (corr.eq.'nnlo' .and. sec(1:2).eq.'rv') then
       if     (sec.eq.'rvewk_is'.and.ch.eq.'ns') then; xsect_checklim = xsect_checklim_HC1C2(yRnd,ff,vegasweight)
       elseif (sec.eq.'rvewk_is'.and.ch.eq.'gq') then; xsect_checklim = xsect_checklim_HC1(yRnd,ff,vegasweight)
       elseif (sec.eq.'rvewk_is'.and.ch.eq.'qg') then; xsect_checklim = xsect_checklim_HC2(yRnd,ff,vegasweight)
       elseif (sec.eq.'rvqcd_is'.and.ch.eq.'ns') then; xsect_checklim = xsect_checklim_nlo6(yRnd,ff,vegasweight)
       elseif (sec(1:8).eq.'rvqcd_fs'.and.ch.eq.'ns') then; xsect_checklim = xsect_checklim_nlo4(yRnd,ff,vegasweight)
       elseif (sec.eq.'rvqcd_is'.and.ch.eq.'aq') then; xsect_checklim = xsect_checklim_HC1C2(yRnd,ff,vegasweight)          
       elseif (sec.eq.'rvqcd_is'.and.ch.eq.'qa') then; xsect_checklim = xsect_checklim_HC1C2(yRnd,ff,vegasweight)
       endif
    endif
    
    ff(1) = zero

  contains

    subroutine error(mycorr)
      character(*), intent(in) :: mycorr

      print *, 'error in the limit routine, for corr ',mycorr
      print *, 'sec = ',sec,'corr = ',corr, 'ch = ',ch
      print *, 'type -h for help'
      stop
    end subroutine error
    
  end function xsect_checklim
  
end module mod_xsects_checklim_gen
  


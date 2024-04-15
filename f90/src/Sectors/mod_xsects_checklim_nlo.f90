module mod_xsects_checklim_nlo
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_auxfunctions
  use mod_kinematics_nlo
  use mod_xsects_nloqcd_r
  use mod_xsects_nloewk_r
  use mod_xsects_nnlo_rvewk
  use mod_xsects_nnlo_rvqcd
  implicit none

  private

  public :: xsect_checklim_HC1C2
  public :: xsect_checklim_HC1,xsect_checklim_HC2
  public :: xsect_checklim_nlo6 !-- H,C1,C2,S,SC1,SC2
  public :: xsect_checklim_nlo4 !-- H,C,S,SC
  
contains

  function xsect_checklim_HC1C2(yRnd,ff,vegasweight)
    integer :: xsect_checklim_HC1C2
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)
    real(dp)   :: myFint(16)

    xsect_checklim_HC1C2 = 0
    ff(1) = zero

    buff   = 0
    buff_r = 0
    call set_buff()

    print *, 'coll1'
    yRnd_tmp = yRnd
    yRnd_tmp(xrho) = yRnd(xrho)*lim
    xsect_checklim_HC1C2 = get_xsec(yRnd_tmp)
    print *, 'Hard  ', myFint(1)
    print *, 'Coll1 ',-myFint(2),-myFint(2)-myFint(1),(-myFint(2)-myFint(1))/(myFint(1))
    print *, ''

    print *, 'coll2'
    yRnd_tmp = yRnd
    yRnd_tmp(xrho) = one-yRnd(xrho)*lim
    xsect_checklim_HC1C2 = get_xsec(yRnd_tmp)
    print *, 'Hard  ', myFint(1)
    print *, 'Coll2 ',-myFint(3),-myFint(3)-myFint(1),(-myFint(3)-myFint(1))/(myFint(1))
    print *, ''
    
    call mypause()

  contains

    function get_xsec(yRnd)
      real(dp15), intent(in) :: yRnd(:)
      integer :: get_xsec

      get_xsec = 0
      
      myFint = zero

#if(_withchecks == 1)
      if (corr.eq.'nloqcd'.and.sec.eq.'r_is' .and. ch.eq.'ns') then
         get_xsec = xsect_nloqcd_r_is_ns(yRnd,ff,vegasweight)
         myFint(1:size(FintNLO_HC1C2)) = FintNLO_HC1C2
      elseif (corr.eq.'nloewk'.and.sec.eq.'r_is' .and. ch.eq.'aq') then
         get_xsec = xsect_nloewk_r_is_aq(yRnd,ff,vegasweight)
         myFint(1:3) = FintNLO_ew(1:3)
      elseif (corr.eq.'nloewk'.and.sec.eq.'r_is' .and. ch.eq.'qa') then
         get_xsec = xsect_nloewk_r_is_qa(yRnd,ff,vegasweight)
         myFint(1:3) = FintNLO_ew(1:3)
      elseif (sec.eq.'rvewk_is' .and. ch.eq.'ns') then
         get_xsec = xsect_nnlo_rvewk_is_ns(yRnd,ff,vegasweight)
         myFint(1:3) = FintRVEWK_ns(1:3)
      elseif (sec.eq.'rvqcd_is' .and. ch.eq.'aq') then
         get_xsec = xsect_nnlo_rvqcd_is_aq(yRnd,ff,vegasweight)
         myFint(1:3) = FintRVQCD_aq(1:3)
      elseif (sec.eq.'rvqcd_is' .and. ch.eq.'qa') then
         get_xsec = xsect_nnlo_rvqcd_is_qa(yRnd,ff,vegasweight)
         myFint(1:3) = FintRVQCD_qa(1:3)
      endif
#else
      print *, 'yRnd',yRnd
      print *, 'vegasweight',vegasweight
#endif
      
    end function get_xsec
    
  end function xsect_checklim_HC1C2

  function xsect_checklim_HC1(yRnd,ff,vegasweight)
    integer :: xsect_checklim_HC1
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)
    real(dp)   :: myFint(16)

    xsect_checklim_HC1 = 0
    ff(1) = zero

    buff   = 0
    buff_r = 0
    call set_buff()

    print *, 'coll1'
    yRnd_tmp = yRnd
    yRnd_tmp(xrho) = yRnd(xrho)*lim
    xsect_checklim_HC1 = get_xsec(yRnd_tmp)
    print *, 'Hard  ', myFint(1)
    print *, 'Coll1 ',-myFint(2),-myFint(2)-myFint(1),(-myFint(2)-myFint(1))/(myFint(1))
    print *, ''
    
    call mypause()

  contains

    function get_xsec(yRnd)
      real(dp15), intent(in) :: yRnd(:)
      integer :: get_xsec

      get_xsec = 0
      
      myFint = zero

#if(_withchecks == 1)      
      if (sec.eq.'r_is' .and. ch.eq.'gq') then
         get_xsec = xsect_nloqcd_r_is_gq(yRnd,ff,vegasweight)
         myFint(1:2) = FintNLO_HC
      elseif (sec.eq.'rvewk_is' .and. ch.eq.'gq') then
         get_xsec = xsect_nnlo_rvewk_is_gq(yRnd,ff,vegasweight)
         myFint(1:2) = FintRVEWK_gq
      endif
#else
      print *, 'yRnd',yRnd
      print *, 'vegasweight',vegasweight
#endif
      
    end function get_xsec
    
  end function xsect_checklim_HC1

  function xsect_checklim_HC2(yRnd,ff,vegasweight)
    integer :: xsect_checklim_HC2
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)
    real(dp)   :: myFint(16)

    xsect_checklim_HC2 = 0
    ff(1) = zero

    buff   = 0
    buff_r = 0
    call set_buff()

    print *, 'coll2'
    yRnd_tmp = yRnd
    yRnd_tmp(xrho) = one-yRnd(xrho)*lim
    xsect_checklim_HC2 = get_xsec(yRnd_tmp)
    print *, 'Hard  ', myFint(1)
    print *, 'Coll2 ',-myFint(2),-myFint(2)-myFint(1),(-myFint(2)-myFint(1))/(myFint(1))
    print *, ''
    
    call mypause()

  contains

    function get_xsec(yRnd)
      real(dp15), intent(in) :: yRnd(:)
      integer :: get_xsec

      get_xsec = 0
      
      myFint = zero

#if(_withchecks == 1)      
      if (sec.eq.'r_is' .and. ch.eq.'qg') then
         get_xsec = xsect_nloqcd_r_is_qg(yRnd,ff,vegasweight)
         myFint(1:size(FintNLO_HC)) = FintNLO_HC
      elseif (sec.eq.'rvewk_is' .and. ch.eq.'qg') then
         get_xsec = xsect_nnlo_rvewk_is_qg(yRnd,ff,vegasweight)
         myFint(1:2) = FintRVEWK_qg
      endif
#else
      print *, 'yRnd',yRnd
      print *, 'vegasweight',vegasweight
#endif      
      
    end function get_xsec
    
  end function xsect_checklim_HC2

  function xsect_checklim_nlo6(yRnd,ff,vegasweight)
    integer :: xsect_checklim_nlo6
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)
    real(dp)   :: myFint(16)

    xsect_checklim_nlo6 = 0
    ff(1) = zero

    buff   = 0
    buff_r = 0
    call set_buff()

    print *, 'coll1'
    yRnd_tmp = yRnd
    yRnd_tmp(xrho) = yRnd(xrho)*lim
    xsect_checklim_nlo6 = get_xsec(yRnd_tmp)
    print *, 'Hard  ', myFint(1)
    print *, 'Coll1 ',-myFint(2),-myFint(2)-myFint(1),(-myFint(2)-myFint(1))/(myFint(1))
    print *, 'S SC1 ',-myFint(4),-myFint(4)-myFint(5),(-myFint(4)-myFint(5))/(myFint(4))
    print *, ''

    print *, 'coll2'
    yRnd_tmp = yRnd
    yRnd_tmp(xrho) = one-yRnd(xrho)*lim
    xsect_checklim_nlo6 = get_xsec(yRnd_tmp)
    print *, 'Hard  ', myFint(1)
    print *, 'Coll2 ',-myFint(3),-myFint(3)-myFint(1),(-myFint(3)-myFint(1))/(myFint(1))
    print *, 'S SC2 ',-myFint(4),-myFint(4)-myFint(6),(-myFint(4)-myFint(6))/(myFint(4))
    print *, ''

    print *, 'soft'
    yRnd_tmp = yRnd
    yRnd_tmp(xE) = yRnd(xE)*lim
    xsect_checklim_nlo6 = get_xsec(yRnd_tmp)
    print *, 'Hard  ', myFint(1)
    print *, 'Soft  ',-myFint(4),-myFint(4)-myFint(1),(-myFint(4)-myFint(1))/(myFint(1))
    print *, 'C1 SC1',-myFint(2),-myFint(2)-myFint(5),(-myFint(2)-myFint(5))/(myFint(2))
    print *, 'C2 SC2',-myFint(3),-myFint(3)-myFint(6),(-myFint(3)-myFint(6))/(myFint(3))
    print *, ''

    print *, 'soft coll1'
    yRnd_tmp = yRnd
    yRnd_tmp(xE) = yRnd(xE)*lim
    yRnd_tmp(xrho) = yRnd(xrho)*lim
    xsect_checklim_nlo6 = get_xsec(yRnd_tmp)
    print *, 'Hard  ', myFint(1)
    print *, 'SC1   ', myFint(5), myFint(5)-myFint(1),( myFint(5)-myFint(1))/(myFint(1))
    print *, 'C1 SC1',-myFint(2),-myFint(2)-myFint(5),(-myFint(2)-myFint(5))/(myFint(2))
    print *, ''

    print *, 'soft coll2'
    yRnd_tmp = yRnd
    yRnd_tmp(xE) = yRnd(xE)*lim
    yRnd_tmp(xrho) = one-yRnd(xrho)*lim
    xsect_checklim_nlo6 = get_xsec(yRnd_tmp)
    print *, 'Hard  ', myFint(1)
    print *, 'SC2   ', myFint(6), myFint(6)-myFint(1),( myFint(6)-myFint(1))/(myFint(1))
    print *, 'C2 SC2',-myFint(3),-myFint(3)-myFint(6),(-myFint(3)-myFint(6))/(myFint(3))
    print *, ''

    call mypause()

  contains

    function get_xsec(yRnd)
      real(dp15), intent(in) :: yRnd(:)
      integer :: get_xsec

      get_xsec = 0
      
      myFint = zero

#if(_withchecks == 1)      
      if (sec.eq.'r_is' .and. ch.eq.'ns') then
         get_xsec = xsect_nloewk_r_is_ns(yRnd,ff,vegasweight)
         myFint(1:6) = FintNLO_ew(1:6)
      elseif (sec.eq.'rvqcd_is' .and. ch.eq.'ns') then
         get_xsec = xsect_nnlo_rvqcd_is_ns(yRnd,ff,vegasweight)
         myFint(1:6) = FintRVQCD_ns(1:6)
      endif
#else
      print *, 'yRnd',yRnd
      print *, 'vegasweight',vegasweight
#endif      
      
    end function get_xsec
    
  end function xsect_checklim_nlo6

  function xsect_checklim_nlo4(yRnd,ff,vegasweight)
    integer :: xsect_checklim_nlo4
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)
    real(dp)   :: myFint(16)

    xsect_checklim_nlo4 = 0
    ff(1) = zero

    buff   = 0
    buff_r = 0
    call set_buff()

    print *, 'coll1'
    yRnd_tmp = yRnd
    yRnd_tmp(xrho) = yRnd(xrho)*lim
    xsect_checklim_nlo4 = get_xsec(yRnd_tmp)
    print *, 'Hard  ', myFint(1)
    print *, 'Coll1 ',-myFint(2),-myFint(2)-myFint(1),(-myFint(2)-myFint(1))/(myFint(1))
    print *, 'S SC1 ',-myFint(3),-myFint(3)-myFint(4),(-myFint(3)-myFint(4))/(myFint(3))
    print *, ''

    print *, 'soft'
    yRnd_tmp = yRnd
    yRnd_tmp(xE) = yRnd(xE)*lim
    xsect_checklim_nlo4 = get_xsec(yRnd_tmp)
    print *, 'Hard  ', myFint(1)
    print *, 'Soft  ',-myFint(3),-myFint(3)-myFint(1),(-myFint(3)-myFint(1))/(myFint(1))
    print *, 'C1 SC1',-myFint(2),-myFint(2)-myFint(4),(-myFint(2)-myFint(4))/(myFint(2))
    print *, ''

    print *, 'soft coll'
    yRnd_tmp = yRnd
    yRnd_tmp(xE) = yRnd(xE)*lim
    yRnd_tmp(xrho) = yRnd(xrho)*lim
    xsect_checklim_nlo4 = get_xsec(yRnd_tmp)
    print *, 'Hard  ', myFint(1)
    print *, 'SC1   ', myFint(4), myFint(4)-myFint(1),( myFint(4)-myFint(1))/(myFint(1))
    print *, 'C1 SC1',-myFint(2),-myFint(2)-myFint(4),(-myFint(2)-myFint(4))/(myFint(2))
    print *, ''

    call mypause()

  contains

    function get_xsec(yRnd)
      real(dp15), intent(in) :: yRnd(:)
      integer :: get_xsec

      get_xsec = 0
      
      myFint = zero

#if(_withchecks == 1)      
      if (sec.eq.'r_fs_53' .and. ch.eq.'ns') then
         get_xsec = xsect_nloewk_r_fs_53_ns(yRnd,ff,vegasweight)
         myFint(1:4) = FintNLO_ew(1:4)
      elseif (sec.eq.'r_fs_54' .and. ch.eq.'ns') then
         get_xsec = xsect_nloewk_r_fs_54_ns(yRnd,ff,vegasweight)
         myFint(1:4) = FintNLO_ew(1:4)

      elseif (sec.eq.'r_fs_53' .and. ch.eq.'aa') then
         get_xsec = xsect_nloewk_r_fs_53_aa(yRnd,ff,vegasweight)
         myFint(1:4) = FintNLO_ew(1:4)
      elseif (sec.eq.'r_fs_54' .and. ch.eq.'aa') then
         get_xsec = xsect_nloewk_r_fs_54_aa(yRnd,ff,vegasweight)
         myFint(1:4) = FintNLO_ew(1:4)

      elseif (sec.eq.'rvqcd_fs_53' .and. ch.eq.'ns') then
         get_xsec = xsect_nnlo_rvqcd_fs_53_ns(yRnd,ff,vegasweight)
         myFint(1:4) = FintRVQCD_ns(1:4)
      elseif (sec.eq.'rvqcd_fs_54' .and. ch.eq.'ns') then
         get_xsec = xsect_nnlo_rvqcd_fs_54_ns(yRnd,ff,vegasweight)
         myFint(1:4) = FintRVQCD_ns(1:4)
        
      endif
#else
      print *, 'yRnd',yRnd
      print *, 'vegasweight',vegasweight
#endif
      
    end function get_xsec
    
  end function xsect_checklim_nlo4

end module mod_xsects_checklim_nlo
  


module mod_xsects_checklim_nnlo_ns
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_auxfunctions
  use mod_kinematics_nnlo_tc_is_eu_ac
  use mod_xsects_nnlo_rr_ns
  implicit none

  private

  public :: xsect_checklim_nnlo4q
  
contains

  function xsect_checklim_nnlo4q(yRnd,ff,vegasweight)
    integer :: xsect_checklim_nnlo4q
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)
    real(dp)   :: myFint(16)

    xsect_checklim_nnlo4q = 0
    ff(1) = zero

    buff    = zero
    buff_rr = zero
    call set_buff()

    print *, 'coll1'
    yRnd_tmp = yRnd
    yRnd_tmp(xrho5) = yRnd(xrho5)*lim
    xsect_checklim_nnlo4q = get_xsec(yRnd_tmp)
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
      if (corr.eq.'nnlo'.and.sec.eq.'rr_5161a' .and. ch.eq.'ns_qqb') then
         get_xsec = xsect_nnlo_rr_5161a_ns_qqb(yRnd,ff,vegasweight)
         myFint(1:2) = FintNNLO_rr_tc_ns(1:2)
      elseif (corr.eq.'nnlo'.and.sec.eq.'rr_5262a' .and. ch.eq.'ns_qqb') then
         get_xsec = xsect_nnlo_rr_5262a_ns_qqb(yRnd,ff,vegasweight)
         myFint(1:2) = FintNNLO_rr_tc_ns(1:2)
      elseif (corr.eq.'nnlo'.and.sec.eq.'rr_5161c' .and. ch.eq.'ns_qqb') then
         get_xsec = xsect_nnlo_rr_5161c_ns_qqb(yRnd,ff,vegasweight)
         myFint(1:2) = FintNNLO_rr_tc_ns(1:2)
      elseif (corr.eq.'nnlo'.and.sec.eq.'rr_5262c' .and. ch.eq.'ns_qqb') then
         get_xsec = xsect_nnlo_rr_5262c_ns_qqb(yRnd,ff,vegasweight)
         myFint(1:2) = FintNNLO_rr_tc_ns(1:2)
      !--
      elseif (corr.eq.'nnlo'.and.sec.eq.'rr_5161a' .and. ch.eq.'ns_qq') then
         get_xsec = xsect_nnlo_rr_5161a_ns_qq(yRnd,ff,vegasweight)
         myFint(1:2) = FintNNLO_rr_tc_ns(1:2)
      elseif (corr.eq.'nnlo'.and.sec.eq.'rr_5262a' .and. ch.eq.'ns_qq') then
         get_xsec = xsect_nnlo_rr_5262a_ns_qq(yRnd,ff,vegasweight)
         myFint(1:2) = FintNNLO_rr_tc_ns(1:2)
      elseif (corr.eq.'nnlo'.and.sec.eq.'rr_5161c' .and. ch.eq.'ns_qq') then
         get_xsec = xsect_nnlo_rr_5161c_ns_qq(yRnd,ff,vegasweight)
         myFint(1:2) = FintNNLO_rr_tc_ns(1:2)
      elseif (corr.eq.'nnlo'.and.sec.eq.'rr_5262c' .and. ch.eq.'ns_qq') then
         get_xsec = xsect_nnlo_rr_5262c_ns_qq(yRnd,ff,vegasweight)
         myFint(1:2) = FintNNLO_rr_tc_ns(1:2)
      endif
#else
      print *, 'yRnd',yRnd
      print *, 'vegasweight',vegasweight      
#endif
      
    end function get_xsec
    
  end function xsect_checklim_nnlo4q

end module mod_xsects_checklim_nnlo_ns
  


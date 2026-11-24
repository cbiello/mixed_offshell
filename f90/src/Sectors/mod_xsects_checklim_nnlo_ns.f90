module mod_xsects_checklim_nnlo_ns
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_auxfunctions
  use mod_kinematics_nnlo_tc_is_eu_ac
  use mod_xsects_nnlo_rr_ns
  use mod_xsects_nnlo_rr_gq
  use mod_xsects_nnlo_rr_qg
  use mod_limvals
  implicit none

  private

  public :: xsect_checklim_nnlo4q
  public :: xsect_checklim_rr_gq       !-- RR gq, qg: H and all counterterms at one point, limits 1e-1..1e-10

contains

  !-- One Vegas point: x2 = E6 (S6), x3 (X = TC or C5), x4 (Y = C or C6) set to 1e-1..1e-10 for each limit, printing
  !-- H and the counterterms in the slots of limval_nnlo with sign (-1)^(number of limits): the terms that must cancel agree.
  !-- The sectors draw their last two variables themselves: fixed per point and passed with override.
  function xsect_checklim_rr_gq(yRnd,ff,vegasweight) result(r)
    integer :: r
    real(dp15) :: yRnd(30),ff(1),vegasweight
    real(dp15) :: y0(30),y(30)
    real(dp)   :: F(8)
    integer    :: il,k,m
    integer, parameter :: slots(8) = [1,3,4,9,5,10,11,15]     !-- H, S6, X, X+S6, Y, Y+S6, X+Y, X+Y+S6
    integer, parameter :: sgn(8) = [1,-1,-1,1,-1,1,1,-1]
    integer, parameter :: vars(3) = [xE6,xRHO5,xRHO6]          !-- bits 0, 1, 2 of the limit m
    character(2) :: nx,ny
    character(7) :: names(8)

    r = 0
    buff    = zero
    buff_rr = zero
    call set_buff()

    y0 = yRnd
    call random_number(y0(kNNLO_max_full-1:kNNLO_max_full))
    override = .true.
    call get_xsec(y0)
    if (F(1) /= zero) then
       nx = merge('TC','C5',scan(sec,'abcd') > 0)
       ny = merge('C ','C6',scan(sec,'abcd') > 0)
       names = [character(7) :: 'H','S6',nx,trim(nx)//'+S6',ny,trim(ny)//'+S6',trim(nx)//'+'//trim(ny),trim(nx)//'+'//trim(ny)//'+S6']
       write(*,'(a,3es14.6)') ' point: x2 (E6), x3, x4 =', y0(vars)
       do m = 1,7
          if (btest(m,2) .and. scan(sec,'a') > 0) cycle    !-- sector a: no Y
          write(*,'(/,a,/,a,*(a16))') ' limit '//trim(names(m+1)), '            x2            x3            x4', adjustr(names)
          do il = 1,10
             y = y0
             do k = 1,3
                if (btest(m,k-1)) y(vars(k)) = 10._dp15**(-il)
             enddo
             call get_xsec(y)
             write(*,'(3es14.6,*(es16.8))') y(vars), sgn*F
          enddo
       enddo
       call mypause()
    endif
    override = .false.
    ff(1) = zero

  contains

    subroutine get_xsec(yy)
      real(dp15), intent(in) :: yy(:)
      integer :: g
      limval_nnlo = zero
#if(_withchecks == 1)
      select case (trim(ch)//' '//trim(sec))
      case ('gq rr_5161a'); g = xsect_nnlo_rr_5161a_gq(yy,ff,vegasweight)
      case ('gq rr_5161b'); g = xsect_nnlo_rr_5161b_gq(yy,ff,vegasweight)
      case ('gq rr_5161c'); g = xsect_nnlo_rr_5161c_gq(yy,ff,vegasweight)
      case ('gq rr_5161d'); g = xsect_nnlo_rr_5161d_gq(yy,ff,vegasweight)
      case ('gq rr_5162');  g = xsect_nnlo_rr_5162_gq(yy,ff,vegasweight)
      case ('gq rr_5163');  g = xsect_nnlo_rr_5163_gq(yy,ff,vegasweight)
      case ('gq rr_5164');  g = xsect_nnlo_rr_5164_gq(yy,ff,vegasweight)
      case ('qg rr_5262a'); g = xsect_nnlo_rr_5262a_qg(yy,ff,vegasweight)
      case ('qg rr_5262b'); g = xsect_nnlo_rr_5262b_qg(yy,ff,vegasweight)
      case ('qg rr_5262c'); g = xsect_nnlo_rr_5262c_qg(yy,ff,vegasweight)
      case ('qg rr_5262d'); g = xsect_nnlo_rr_5262d_qg(yy,ff,vegasweight)
      case ('qg rr_5261');  g = xsect_nnlo_rr_5261_qg(yy,ff,vegasweight)
      case ('qg rr_5263');  g = xsect_nnlo_rr_5263_qg(yy,ff,vegasweight)
      case ('qg rr_5264');  g = xsect_nnlo_rr_5264_qg(yy,ff,vegasweight)
      case default; print *, 'xsect_checklim_rr_gq: unknown sector ',ch,sec; stop
      end select
#else
      stop 'xsect_checklim_rr_gq: build with _withchecks = 1'
#endif
      F = limval_nnlo(slots)
    end subroutine get_xsec

  end function xsect_checklim_rr_gq

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
  


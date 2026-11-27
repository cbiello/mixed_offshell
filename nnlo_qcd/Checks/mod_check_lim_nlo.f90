module mod_check_lim_nlo
  use mod_types; use common_def; use mod_consts_dp
  use mod_consts_debug
  use mod_histograms
  use mod_xsects_nlo
  implicit none
  
contains

  function xsect_nlo_checklim_r(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nlo_checklim_r,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)

    ! !-- soft limit
    yRnd_tmp = yRnd

    yRnd_tmp(xE_nlo) = (yRnd_tmp(xE_nlo) * limcheck - buff)/onet
    xsect_nlo_checklim_r = xsect_nlo_r(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    print *, 'hard ', FintNLO_CheckLim(1)
    print *, 'soft ', FintNLO_CheckLim(2), one + FintNLO_CheckLim(2)/FintNLO_CheckLim(1)

    print *, ''

    !-- c41 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xRHO_nlo) = (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nlo_checklim_r = xsect_nlo_r(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    print *, 'hard ', FintNLO_CheckLim(1)
    print *, 'c41  ', FintNLO_CheckLim(3), one + FintNLO_CheckLim(3)/FintNLO_CheckLim(1)

    print *, ''

    !-- c42 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xRHO_nlo) = one - (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nlo_checklim_r = xsect_nlo_r(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    print *, 'hard ', FintNLO_CheckLim(1)
    print *, 'c42  ', FintNLO_CheckLim(4), one + FintNLO_CheckLim(4)/FintNLO_CheckLim(1)

    print *, ''
    
    print *, 'press a key to continue'
    read(*,*)

    return

  end function xsect_nlo_checklim_r
    
end module mod_check_lim_nlo


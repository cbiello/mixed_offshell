module mod_check_lim_rv
  use mod_types; use common_def; use mod_consts_dp
  use mod_consts_debug
  use mod_histograms
  implicit none

contains

  function xsect_nnlo_checklim_rv_wp_ns(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    use mod_xsects_nnlo_rv_wp_ns
    integer :: xsect_nnlo_checklim_rv_wp_ns,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)

    !-- soft limit
    yRnd_tmp = yRnd

    yRnd_tmp(xE_nlo) = (yRnd_tmp(xE_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_rv_wp_ns = xsect_nnlo_rv_wp_ns(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    print *, 'hard ', Fint_fullONLO_rv_CheckLim(1)
    print *, 'soft ', Fint_fullONLO_rv_CheckLim(2), one + Fint_fullONLO_rv_CheckLim(2)/Fint_fullONLO_rv_CheckLim(1)

    print *, ''

    !-- c41s limit
    yRnd_tmp = yRnd

    yRnd_tmp(xE_nlo) = (yRnd_tmp(xE_nlo) * limcheck - buff)/onet
    yRnd_tmp(xRHO_nlo) = (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_rv_wp_ns = xsect_nnlo_rv_wp_ns(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    print *, 'hard ', Fint_fullONLO_rv_CheckLim(1)
    print *, 'c41s ', Fint_fullONLO_rv_CheckLim(3), one - Fint_fullONLO_rv_CheckLim(3)/Fint_fullONLO_rv_CheckLim(1)

    print *, ''

    !-- c42s limit
    yRnd_tmp = yRnd

    yRnd_tmp(xE_nlo) = (yRnd_tmp(xE_nlo) * limcheck - buff)/onet
    yRnd_tmp(xRHO_nlo) = one - (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_rv_wp_ns = xsect_nnlo_rv_wp_ns(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    print *, 'hard ', Fint_fullONLO_rv_CheckLim(1)
    print *, 'c42s ', Fint_fullONLO_rv_CheckLim(4), one - Fint_fullONLO_rv_CheckLim(4)/Fint_fullONLO_rv_CheckLim(1)

    print *, ''

    !-- c41 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xRHO_nlo) = (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_rv_wp_ns = xsect_nnlo_rv_wp_ns(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    print *, 'hard ', Fint_fullONLO_rv_CheckLim(1)
    print *, 'c41  ', Fint_fullONLO_rv_CheckLim(5), one + Fint_fullONLO_rv_CheckLim(5)/Fint_fullONLO_rv_CheckLim(1)

    print *, ''

    !-- c42 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xRHO_nlo) = one - (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_rv_wp_ns = xsect_nnlo_rv_wp_ns(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    print *, 'hard ', Fint_fullONLO_rv_CheckLim(1)
    print *, 'c42  ', Fint_fullONLO_rv_CheckLim(6), one + Fint_fullONLO_rv_CheckLim(6)/Fint_fullONLO_rv_CheckLim(1)

    print *, ''
    
    print *, 'press a key to continue'
    read(*,*)

    return

  end function xsect_nnlo_checklim_rv_wp_ns

  function xsect_nnlo_checklim_rv_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    use mod_xsects_nnlo_rv_wp_gq
    integer :: xsect_nnlo_checklim_rv_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)

    !-- c41 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xRHO_nlo) = (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_rv_wp_gq = xsect_nnlo_rv_wp_gq(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    print *, 'hard ', Fint_c41only_rv_CheckLim(1)
    print *, 'c41  ', Fint_c41only_rv_CheckLim(2), one + Fint_c41only_rv_CheckLim(2)/Fint_c41only_rv_CheckLim(1)

    print *, ''
    
    print *, 'press a key to continue'
    read(*,*)

    return

  end function xsect_nnlo_checklim_rv_wp_gq

  function xsect_nnlo_checklim_rv_wp_qg(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    use mod_xsects_nnlo_rv_wp_qg
    integer :: xsect_nnlo_checklim_rv_wp_qg,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)

    !-- c42 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xRHO_nlo) = one-(yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_rv_wp_qg = xsect_nnlo_rv_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    print *, 'hard ', Fint_c42only_rv_CheckLim(1)
    print *, 'c42  ', Fint_c42only_rv_CheckLim(2), one + Fint_c42only_rv_CheckLim(2)/Fint_c42only_rv_CheckLim(1)

    print *, ''
    
    print *, 'press a key to continue'
    read(*,*)

    return

  end function xsect_nnlo_checklim_rv_wp_qg

end module mod_check_lim_rv


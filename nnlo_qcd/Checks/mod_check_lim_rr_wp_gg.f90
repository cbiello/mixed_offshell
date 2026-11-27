module mod_check_lim_rr_wp_gg
  use mod_types; use common_def; use mod_consts_dp
  use mod_consts_debug
  use mod_histograms
  use mod_xsects_nnlo_rr_wp_gg
  implicit none

contains

  function xsect_nnlo_checklim_rr_wp_gg(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_checklim_rr_wp_gg,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)
#if (_DebugFeatures == 0)
    character :: globalsec*(10) !-- to compile without global variables defined
#endif

    buffrr = zero

    !-- c41
    yRnd_tmp = yRnd
    yRnd_tmp(xRHO1) = (yRnd_tmp(xRHO1) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_wp_gg = xsect_nnlo_rr_wp_gg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    print *, 'hard    ', FintNNLO_CheckLim_wp_gg(1)
    print *, 'c41     ', FintNNLO_CheckLim_wp_gg(2), one + FintNNLO_CheckLim_wp_gg(2)/FintNNLO_CheckLim_wp_gg(1)
    !print *, ''

    !-- c51
    yRnd_tmp = yRnd
    yRnd_tmp(xRHO2) = (yRnd_tmp(xRHO2) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_wp_gg = xsect_nnlo_rr_wp_gg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_gg(1)
    print *, 'c51     ', FintNNLO_CheckLim_wp_gg(3), one + FintNNLO_CheckLim_wp_gg(3)/FintNNLO_CheckLim_wp_gg(1)
    !print *, ''

    !-- c42
    yRnd_tmp = yRnd
    yRnd_tmp(xRHO1) = one-(yRnd_tmp(xRHO1) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_wp_gg = xsect_nnlo_rr_wp_gg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_gg(1)
    print *, 'c42     ', FintNNLO_CheckLim_wp_gg(4), one + FintNNLO_CheckLim_wp_gg(4)/FintNNLO_CheckLim_wp_gg(1)
    !print *, ''

    !-- c52
    yRnd_tmp = yRnd
    yRnd_tmp(xRHO2) = one-(yRnd_tmp(xRHO2) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_wp_gg = xsect_nnlo_rr_wp_gg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_gg(1)
    print *, 'c52     ', FintNNLO_CheckLim_wp_gg(5), one + FintNNLO_CheckLim_wp_gg(5)/FintNNLO_CheckLim_wp_gg(1)
    !print *, ''

    !-- c41c52
    yRnd_tmp = yRnd
    yRnd_tmp(xRHO1) = (yRnd_tmp(xRHO1) * limcheck - buff)/onet
    yRnd_tmp(xRHO2) = one-(yRnd_tmp(xRHO2) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_wp_gg = xsect_nnlo_rr_wp_gg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_gg(1)
    print *, 'c41c52  ', FintNNLO_CheckLim_wp_gg(6), one - FintNNLO_CheckLim_wp_gg(6)/FintNNLO_CheckLim_wp_gg(1)
    !print *, ''

    !-- c42c51
    yRnd_tmp = yRnd
    yRnd_tmp(xRHO1) = one-(yRnd_tmp(xRHO1) * limcheck - buff)/onet
    yRnd_tmp(xRHO2) = (yRnd_tmp(xRHO2) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_wp_gg = xsect_nnlo_rr_wp_gg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_gg(1)
    print *, 'c42c51  ', FintNNLO_CheckLim_wp_gg(7), one - FintNNLO_CheckLim_wp_gg(7)/FintNNLO_CheckLim_wp_gg(1)
    !print *, ''

    print *, 'press a key to continue'
    read(*,*)

    return

  end function xsect_nnlo_checklim_rr_wp_gg
        
end module mod_check_lim_rr_wp_gg


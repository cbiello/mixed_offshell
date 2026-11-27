module mod_check_lim_rr_wp_qg
  use mod_types; use common_def; use mod_consts_dp
  use mod_consts_debug
  use mod_histograms
  use mod_xsects_nnlo_rr_wp_qg
  implicit none

contains

  function xsect_nnlo_checklim_rr_tc_wp_qg(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_checklim_rr_tc_wp_qg,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)
#if (_DebugFeatures == 0)
    character :: globalsec*(10) !-- to compile without global variables defined
#endif

    buffrr = zero

    !-- s5
    yRnd_tmp = yRnd
    yRnd_tmp(xE5) = (yRnd_tmp(xE5) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_tc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 's5       ', FintNNLO_CheckLim_wp_qg(2), one + FintNNLO_CheckLim_wp_qg(2)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    !-- s5c5
    yRnd_tmp = yRnd
    yRnd_tmp(xE5) = (yRnd_tmp(xE5) * limcheck - buff)/onet
    yRnd_tmp(xRHO2) = (yRnd_tmp(xRHO2) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_tc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 's5c5     ', FintNNLO_CheckLim_wp_qg(3), one - FintNNLO_CheckLim_wp_qg(3)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    !-- c5
    yRnd_tmp = yRnd
    yRnd_tmp(xRHO2) = (yRnd_tmp(xRHO2) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_tc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 'c5       ', FintNNLO_CheckLim_wp_qg(4), one + FintNNLO_CheckLim_wp_qg(4)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    !-- tc
    yRnd_tmp = yRnd
    yRnd_tmp(xRHO1) = (yRnd_tmp(xRHO1) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_tc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 'tc       ', FintNNLO_CheckLim_wp_qg(5), one + FintNNLO_CheckLim_wp_qg(5)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    !-- tcc5
    yRnd_tmp = yRnd
    yRnd_tmp(xRHO1) = (yRnd_tmp(xRHO1) * limcheck - buff)/onet
    yRnd_tmp(xRHO2) = (yRnd_tmp(xRHO2) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_tc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 'tcc5     ', FintNNLO_CheckLim_wp_qg(6), one - FintNNLO_CheckLim_wp_qg(6)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    !-- tcs5
    yRnd_tmp = yRnd
    yRnd_tmp(xE5) = (yRnd_tmp(xE5) * limcheck - buff)/onet
    yRnd_tmp(xRHO1) = (yRnd_tmp(xRHO1) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_tc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 'tcs5     ', FintNNLO_CheckLim_wp_qg(7), one - FintNNLO_CheckLim_wp_qg(7)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    !-- tcc5s5
    yRnd_tmp = yRnd
    yRnd_tmp(xE5) = (yRnd_tmp(xE5) * limcheck - buff)/onet
    yRnd_tmp(xRHO1) = (yRnd_tmp(xRHO1) * limcheck - buff)/onet
    yRnd_tmp(xRHO2) = (yRnd_tmp(xRHO2) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_tc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 'tcc5s5   ', FintNNLO_CheckLim_wp_qg(8), one + FintNNLO_CheckLim_wp_qg(8)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    print *, 'press a key to continue'
    read(*,*)

    return

  end function xsect_nnlo_checklim_rr_tc_wp_qg

  function xsect_nnlo_checklim_rr_dc_wp_qg(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_checklim_rr_dc_wp_qg,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)
#if (_DebugFeatures == 0)
    character :: globalsec*(10) !-- to compile without global variables defined
#endif

    buffrr = zero

    !-- s5
    yRnd_tmp = yRnd
    yRnd_tmp(xE5) = (yRnd_tmp(xE5) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_dc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 's5       ', FintNNLO_CheckLim_wp_qg(2), one + FintNNLO_CheckLim_wp_qg(2)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    !-- s5c5
    yRnd_tmp = yRnd
    yRnd_tmp(xE5) = (yRnd_tmp(xE5) * limcheck - buff)/onet
    yRnd_tmp(xRHO2) = (yRnd_tmp(xRHO2) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_dc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 's5c5     ', FintNNLO_CheckLim_wp_qg(3), one - FintNNLO_CheckLim_wp_qg(3)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    !-- c4
    yRnd_tmp = yRnd
    yRnd_tmp(xRHO1) = (yRnd_tmp(xRHO1) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_dc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 'c4       ', FintNNLO_CheckLim_wp_qg(4), one + FintNNLO_CheckLim_wp_qg(4)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    !-- c5
    yRnd_tmp = yRnd
    yRnd_tmp(xRHO2) = (yRnd_tmp(xRHO2) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_dc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 'c5       ', FintNNLO_CheckLim_wp_qg(5), one + FintNNLO_CheckLim_wp_qg(5)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    !-- c4s5
    yRnd_tmp = yRnd
    yRnd_tmp(xE5) = (yRnd_tmp(xE5) * limcheck - buff)/onet
    yRnd_tmp(xRHO1) = (yRnd_tmp(xRHO1) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_dc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 'c4s5     ', FintNNLO_CheckLim_wp_qg(6), one - FintNNLO_CheckLim_wp_qg(6)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    !-- c4c5s5
    yRnd_tmp = yRnd
    yRnd_tmp(xE5) = (yRnd_tmp(xE5) * limcheck - buff)/onet
    yRnd_tmp(xRHO1) = (yRnd_tmp(xRHO1) * limcheck - buff)/onet
    yRnd_tmp(xRHO2) = (yRnd_tmp(xRHO2) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_dc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 'c4c5s5   ', FintNNLO_CheckLim_wp_qg(7), one + FintNNLO_CheckLim_wp_qg(7)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    !-- c4c5
    yRnd_tmp = yRnd
    yRnd_tmp(xRHO1) = (yRnd_tmp(xRHO1) * limcheck - buff)/onet
    yRnd_tmp(xRHO2) = (yRnd_tmp(xRHO2) * limcheck - buff)/onet
    xsect_nnlo_checklim_rr_dc_wp_qg = getfunc_rr_wp_qg(globalsec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    !
    !print *, 'hard     ', FintNNLO_CheckLim_wp_qg(1)
    print *, 'c4c5     ', FintNNLO_CheckLim_wp_qg(8), one - FintNNLO_CheckLim_wp_qg(8)/FintNNLO_CheckLim_wp_qg(1)
    !print *, ''

    print *, 'press a key to continue'
    read(*,*)

    return

  end function xsect_nnlo_checklim_rr_dc_wp_qg

  function getfunc_rr_wp_qg(sec,ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    character :: sec*(10)
    integer :: getfunc_rr_wp_qg,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: ff(1),vegasweight,yRnd_tmp(30)

    if (sec .eq. 'rr4252a') getfunc_rr_wp_qg = xsect_nnlo_4252a_rr_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    if (sec .eq. 'rr4252b') getfunc_rr_wp_qg = xsect_nnlo_4252b_rr_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    if (sec .eq. 'rr4252c') getfunc_rr_wp_qg = xsect_nnlo_4252c_rr_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    if (sec .eq. 'rr4252d') getfunc_rr_wp_qg = xsect_nnlo_4252d_rr_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    if (sec .eq. 'rr4251') getfunc_rr_wp_qg = xsect_nnlo_4251_rr_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    if (sec .eq. 'rr4151a') getfunc_rr_wp_qg = xsect_nnlo_4151a_rr_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    if (sec .eq. 'rr4151b') getfunc_rr_wp_qg = xsect_nnlo_4151b_rr_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    if (sec .eq. 'rr4151c') getfunc_rr_wp_qg = xsect_nnlo_4151c_rr_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    if (sec .eq. 'rr4151d') getfunc_rr_wp_qg = xsect_nnlo_4151d_rr_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    if (sec .eq. 'rr4152') getfunc_rr_wp_qg = xsect_nnlo_4152_rr_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

  end function getfunc_rr_wp_qg
        
end module mod_check_lim_rr_wp_qg


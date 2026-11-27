module mod_check_lim_sub
  use mod_types; use common_def; use mod_consts_dp
  use mod_consts_debug
  use mod_histograms
  use mod_xsects_nnlo_sub_wp_ns
  use mod_xsects_nnlo_sub_wp_gq
  use mod_xsects_nnlo_sub_wp_qg
  implicit none
  
contains


  function xsect_nnlo_checklim_sub124_wp_ns(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_checklim_sub124_wp_ns,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)
    integer :: i
    logical :: check_qqbg_z124,check_qqbg_1z24,check_qqbg_124,check_gqq_z124,check_qgq_1z24

    bufftc = zero
    buffrr = zero

    check_qqbg_z124= .false.
    check_qqbg_1z24= .false.
    check_qqbg_124 = .false.
    check_gqq_z124= .false.
    check_qgq_1z24= .true.
    

    print *, 'careful, if sping is included the soft and soft/collinear individually should not work properly'
    print *, 'since we assume cancellation at the sping level so we do not include them'
    print *, ''

!    do i =1,10
    !-- soft limit
    yRnd_tmp = yRnd

    yRnd_tmp(xE_nlo) = (yRnd_tmp(xE_nlo) * limcheck - buff)/onet

    xsect_nnlo_checklim_sub124_wp_ns = xsect_nnlo_sub124_wp_ns(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

       if (check_qqbg_z124) then
          print *, 'hard z124', Fint_qqb_ns_sub_QQBG_z124_CheckLim(1)
          print *, 'soft z124', Fint_qqb_ns_sub_QQBG_z124_CheckLim(2), one + Fint_qqb_ns_sub_QQBG_z124_CheckLim(2)/Fint_qqb_ns_sub_QQBG_z124_CheckLim(1)
       endif
       if (check_qqbg_1z24) then
          print *, 'hard 1z24', Fint_qqb_ns_sub_QQBG_1z24_CheckLim(1)
          print *, 'soft 1z24', Fint_qqb_ns_sub_QQBG_1z24_CheckLim(2), one + Fint_qqb_ns_sub_QQBG_1z24_CheckLim(2)/Fint_qqb_ns_sub_QQBG_1z24_CheckLim(1)
        endif
       if (check_qqbg_124) then
          print *, 'hard 124', Fint_qqb_ns_sub_QQBG_124_CheckLim(1)
          print *, 'soft 124', Fint_qqb_ns_sub_QQBG_124_CheckLim(2), one + Fint_qqb_ns_sub_QQBG_124_CheckLim(2)/Fint_qqb_ns_sub_QQBG_124_CheckLim(1)
       endif
       if (check_gqq_z124) then
          print *, 'hard gqq z124', Fint_qqb_ns_sub_GQQ_z124_CheckLim(1)
       endif
       if (check_qgq_1z24) then
          print *, 'hard qgq 1z24', Fint_qqb_ns_sub_Qgq_1z24_CheckLim(1)
       endif
       print *, ''

    !-- c41 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xRHO_nlo) = (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet

    xsect_nnlo_checklim_sub124_wp_ns = xsect_nnlo_sub124_wp_ns(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

       if (check_qqbg_z124) then
          print *, 'hard z124', Fint_qqb_ns_sub_QQBG_z124_CheckLim(1)
          print *, 'c41  z124', Fint_qqb_ns_sub_QQBG_z124_CheckLim(3), one + Fint_qqb_ns_sub_QQBG_z124_CheckLim(3)/Fint_qqb_ns_sub_QQBG_z124_CheckLim(1)
       endif
       if (check_qqbg_1z24) then
          print *, 'hard 1z24', Fint_qqb_ns_sub_QQBG_1z24_CheckLim(1)
          print *, 'c41  1z24', Fint_qqb_ns_sub_QQBG_1z24_CheckLim(3), one + Fint_qqb_ns_sub_QQBG_1z24_CheckLim(3)/Fint_qqb_ns_sub_QQBG_1z24_CheckLim(1)
       endif
       if (check_qqbg_124) then
          print *, 'hard 124', Fint_qqb_ns_sub_QQBG_124_CheckLim(1)
          print *, 'c41  124', Fint_qqb_ns_sub_QQBG_124_CheckLim(3), one + Fint_qqb_ns_sub_QQBG_124_CheckLim(3)/Fint_qqb_ns_sub_QQBG_124_CheckLim(1)
       endif
       if (check_gqq_z124) then
          print *, 'hard gqq z124', Fint_qqb_ns_sub_GQQ_z124_CheckLim(1)
          print *, 'c41  gqq z124', Fint_qqb_ns_sub_GQQ_z124_CheckLim(2), one + Fint_qqb_ns_sub_GQQ_z124_CheckLim(2)/Fint_qqb_ns_sub_GQQ_z124_CheckLim(1)
        endif
        if (check_qgq_1z24) then
           print *, 'hard qgq 1z24', Fint_qqb_ns_sub_Qgq_1z24_CheckLim(1)
        endif
           
    print *, ''


    !-- c42 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xRHO_nlo) = one - (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_sub124_wp_ns = xsect_nnlo_sub124_wp_ns(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

       if (check_qqbg_z124) then
          print *, 'hard z124', Fint_qqb_ns_sub_QQBG_z124_CheckLim(1)
          print *, 'c42  z114', Fint_qqb_ns_sub_QQBG_z124_CheckLim(4), one + Fint_qqb_ns_sub_QQBG_z124_CheckLim(4)/Fint_qqb_ns_sub_QQBG_z124_CheckLim(1)
       endif
       if (check_qqbg_1z24) then
          print *, 'hard 1z24', Fint_qqb_ns_sub_QQBG_1z24_CheckLim(1)
          print *, 'c42  1z24', Fint_qqb_ns_sub_QQBG_1z24_CheckLim(4), one + Fint_qqb_ns_sub_QQBG_1z24_CheckLim(4)/Fint_qqb_ns_sub_QQBG_1z24_CheckLim(1)
       endif
       if (check_qqbg_124) then
          print *, 'hard 124', Fint_qqb_ns_sub_QQBG_124_CheckLim(1)
          print *, 'c42  124', Fint_qqb_ns_sub_QQBG_124_CheckLim(4), one + Fint_qqb_ns_sub_QQBG_124_CheckLim(4)/Fint_qqb_ns_sub_QQBG_124_CheckLim(1)
       endif
       if (check_gqq_z124) then 
          print *, 'hard gqq z124', Fint_qqb_ns_sub_GQQ_z124_CheckLim(1)
       endif
        if (check_qgq_1z24) then
           print *, 'hard qgq 1z24', Fint_qqb_ns_sub_Qgq_1z24_CheckLim(1)
           print *, 'c42  qgq 1z24', Fint_qqb_ns_sub_Qgq_1z24_CheckLim(2), one + Fint_qqb_ns_sub_Qgq_1z24_CheckLim(2)/Fint_qqb_ns_sub_Qgq_1z24_CheckLim(1)
        endif
    print *, ''


    !-- sc41 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xE_nlo) = (yRnd_tmp(xE_nlo) * limcheck - buff)/onet
    yRnd_tmp(xRHO_nlo) = (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_sub124_wp_ns = xsect_nnlo_sub124_wp_ns(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

       if (check_qqbg_z124) then
          print *, 'hard z124', Fint_qqb_ns_sub_QQBG_z124_CheckLim(1)
          print *, 'sc41 z124', Fint_qqb_ns_sub_QQBG_z124_CheckLim(5), one - Fint_qqb_ns_sub_QQBG_z124_CheckLim(5)/Fint_qqb_ns_sub_QQBG_z124_CheckLim(1)
       endif
       if (check_qqbg_1z24) then
          print *, 'hard 1z24', Fint_qqb_ns_sub_QQBG_1z24_CheckLim(1)
          print *, 'sc41 1z24', Fint_qqb_ns_sub_QQBG_1z24_CheckLim(5), one - Fint_qqb_ns_sub_QQBG_1z24_CheckLim(5)/Fint_qqb_ns_sub_QQBG_1z24_CheckLim(1)
       endif
       if (check_qqbg_124) then
          print *, 'hard 124', Fint_qqb_ns_sub_QQBG_124_CheckLim(1)
          print *, 'sc41 124', Fint_qqb_ns_sub_QQBG_124_CheckLim(5), one - Fint_qqb_ns_sub_QQBG_124_CheckLim(5)/Fint_qqb_ns_sub_QQBG_124_CheckLim(1)
        endif
        if (check_gqq_z124) then
           print *, 'hard gqq z124', Fint_qqb_ns_sub_GQQ_z124_CheckLim(1)
           print *, 'sc41 gqq z124', Fint_qqb_ns_sub_GQQ_z124_CheckLim(2), one + Fint_qqb_ns_sub_GQQ_z124_CheckLim(2)/Fint_qqb_ns_sub_GQQ_z124_CheckLim(1)
        endif
        if (check_qgq_1z24) then
          print *, 'hard qgq 1z24', Fint_qqb_ns_sub_Qgq_1z24_CheckLim(1)
           print *, ''
        endif

    !-- sc42 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xE_nlo) = (yRnd_tmp(xE_nlo) * limcheck - buff)/onet
    yRnd_tmp(xRHO_nlo) = one - (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_sub124_wp_ns = xsect_nnlo_sub124_wp_ns(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

       if (check_qqbg_z124) then
          print *, 'hard z124', Fint_qqb_ns_sub_QQBG_z124_CheckLim(1)
          print *, 'sc42 z124', Fint_qqb_ns_sub_QQBG_z124_CheckLim(6), one - Fint_qqb_ns_sub_QQBG_z124_CheckLim(6)/Fint_qqb_ns_sub_QQBG_z124_CheckLim(1)
       endif
       if (check_qqbg_1z24) then
          print *, 'hard 1z24', Fint_qqb_ns_sub_QQBG_1z24_CheckLim(1)
          print *, 'sc42 1z24', Fint_qqb_ns_sub_QQBG_1z24_CheckLim(6), one - Fint_qqb_ns_sub_QQBG_1z24_CheckLim(6)/Fint_qqb_ns_sub_QQBG_1z24_CheckLim(1)
        endif
       if (check_qqbg_124) then
          print *, 'hard 124', Fint_qqb_ns_sub_QQBG_124_CheckLim(1)
          print *, 'sc42 124', Fint_qqb_ns_sub_QQBG_124_CheckLim(6), one - Fint_qqb_ns_sub_QQBG_124_CheckLim(6)/Fint_qqb_ns_sub_QQBG_124_CheckLim(1)
       endif
       if (check_gqq_z124) then
          print *, 'hard gqq z124', Fint_qqb_ns_sub_GQQ_z124_CheckLim(1)
       endif
       if (check_qgq_1z24) then    
          print *, 'hard qgq 1z24', Fint_qqb_ns_sub_Qgq_1z24_CheckLim(1)
          print *, 'sc42  qgq 1z24', Fint_qqb_ns_sub_Qgq_1z24_CheckLim(2), one + Fint_qqb_ns_sub_Qgq_1z24_CheckLim(2)/Fint_qqb_ns_sub_Qgq_1z24_CheckLim(1)
       endif
    print *, ''

    
    print *, 'press a key to continue'
    read(*,*)


    return
    
  end function xsect_nnlo_checklim_sub124_wp_ns
  




  function xsect_nnlo_checklim_sub124_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_checklim_sub124_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)
    logical :: check_gqq_z124,check_gqq_1z24,check_gqq_124,check_qqbg_z124

    bufftc = zero
    buffrr = zero

    check_gqq_z124= .false.
    check_gqq_1z24= .false.
    check_gqq_124 = .false.
    check_qqbg_z124= .true.


    print *, "Checking limits for  sub124 terms in gq->qg channel"
    print *, " "

    print *, 'careful, if sping is included the soft and soft/collinear individually should not work properly'
    print *, 'since we assume cancellation at the sping level so we do not include them'
    print *, ''

    !-- soft limit
    yRnd_tmp = yRnd

    yRnd_tmp(xE_nlo) = (yRnd_tmp(xE_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_sub124_wp_gq = xsect_nnlo_sub124_wp_gq(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

       if (check_gqq_z124) then
          print *, 'hard gqq z124', Fint_gq_qg_sub_GQQ_z124_CheckLim(1)
       endif
       if (check_gqq_1z24) then
          print *, 'hard gqq 1z24', Fint_gq_qg_sub_GQQ_1z24_CheckLim(1)
       endif
       if (check_gqq_124) then
          print *, 'hard gqq 124', Fint_gq_qg_sub_GQQ_124_CheckLim(1)
       endif
       if (check_qqbg_z124) then
          print *, 'hard qqbg z124', Fint_gq_qg_sub_QQBG_z124_CheckLim(1)
          print *, 'soft qqbg z124', Fint_gq_qg_sub_QQBG_z124_CheckLim(2), one + Fint_gq_qg_sub_QQBG_z124_CheckLim(2)/Fint_gq_qg_sub_QQBG_z124_CheckLim(1)
       endif

    print *, ''

    !-- c41 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xRHO_nlo) = (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_sub124_wp_gq = xsect_nnlo_sub124_wp_gq(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    if (check_gqq_z124) then
       print *, 'hard gqq z124', Fint_gq_qg_sub_GQQ_z124_CheckLim(1)
       print *, 'c41  gqq z124', Fint_gq_qg_sub_GQQ_z124_CheckLim(2), one + Fint_gq_qg_sub_GQQ_z124_CheckLim(2)/Fint_gq_qg_sub_GQQ_z124_CheckLim(1)
    endif
    if (check_gqq_1z24) then
       print *, 'hard gqq 1z24', Fint_gq_qg_sub_GQQ_1z24_CheckLim(1)
       print *, 'c41  gqq 1z24', Fint_gq_qg_sub_GQQ_1z24_CheckLim(2), one + Fint_gq_qg_sub_GQQ_1z24_CheckLim(2)/Fint_gq_qg_sub_GQQ_1z24_CheckLim(1)
    endif
    if (check_gqq_124) then
       print *, 'hard gqq 124', Fint_gq_qg_sub_GQQ_124_CheckLim(1)
       print *, 'c41  gqq 124', Fint_gq_qg_sub_GQQ_124_CheckLim(2), one + Fint_gq_qg_sub_GQQ_124_CheckLim(2)/Fint_gq_qg_sub_GQQ_124_CheckLim(1)
    endif
    if (check_qqbg_z124) then
       print *, 'hard qqbg z124', Fint_gq_qg_sub_QQBG_z124_CheckLim(1)
       print *, 'c41  qqbg z124', Fint_gq_qg_sub_QQBG_z124_CheckLim(3), one + Fint_gq_qg_sub_QQBG_z124_CheckLim(3)/Fint_gq_qg_sub_QQBG_z124_CheckLim(1)
    endif
    print *, ' '


    !-- c42 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xRHO_nlo) = one - (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_sub124_wp_gq = xsect_nnlo_sub124_wp_gq(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    if (check_gqq_z124) then
       print *, 'hard gqq z124', Fint_gq_qg_sub_GQQ_z124_CheckLim(1)
    endif
    if (check_gqq_1z24) then
       print *, 'hard gqq 1z24', Fint_gq_qg_sub_GQQ_1z24_CheckLim(1)
    endif
    if (check_gqq_124) then
       print *, 'hard gqq 124', Fint_gq_qg_sub_GQQ_124_CheckLim(1)
    endif
    if (check_qqbg_z124) then
       print *, 'hard qqbg z124', Fint_gq_qg_sub_QQBG_z124_CheckLim(1)
       print *, 'col2 qqbg z124', Fint_gq_qg_sub_QQBG_z124_CheckLim(4), one + Fint_gq_qg_sub_QQBG_z124_CheckLim(4)/Fint_gq_qg_sub_QQBG_z124_CheckLim(1)
    endif

    print *, ''


    !-- sc41 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xE_nlo) = (yRnd_tmp(xE_nlo) * limcheck - buff)/onet
    yRnd_tmp(xRHO_nlo) = (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_sub124_wp_gq = xsect_nnlo_sub124_wp_gq(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    if (check_gqq_z124) then
       print *, 'hard gqq z124', Fint_gq_qg_sub_GQQ_z124_CheckLim(1)
       print *, 'c41  gqq z124', Fint_gq_qg_sub_GQQ_z124_CheckLim(2), one + Fint_gq_qg_sub_GQQ_z124_CheckLim(2)/Fint_gq_qg_sub_GQQ_z124_CheckLim(1)
    endif
    if (check_gqq_1z24) then
       print *, 'hard gqq 1z24', Fint_gq_qg_sub_GQQ_1z24_CheckLim(1)
       print *, 'c41  gqq 1z24', Fint_gq_qg_sub_GQQ_1z24_CheckLim(2), one + Fint_gq_qg_sub_GQQ_1z24_CheckLim(2)/Fint_gq_qg_sub_GQQ_1z24_CheckLim(1)
    endif
    if (check_gqq_124) then
       print *, 'hard gqq 124', Fint_gq_qg_sub_GQQ_124_CheckLim(1)
       print *, 'c41  gqq 124', Fint_gq_qg_sub_GQQ_124_CheckLim(2), one + Fint_gq_qg_sub_GQQ_124_CheckLim(2)/Fint_gq_qg_sub_GQQ_124_CheckLim(1)
    endif
    if (check_qqbg_z124) then
       print *, 'hard qqbg z124', Fint_gq_qg_sub_QQBG_z124_CheckLim(1)
       print *, 'scol41 qqbg z124', Fint_gq_qg_sub_QQBG_z124_CheckLim(5), one - Fint_gq_qg_sub_QQBG_z124_CheckLim(5)/Fint_gq_qg_sub_QQBG_z124_CheckLim(1)
    endif

    print *, ''


    !-- sc42 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xE_nlo) = (yRnd_tmp(xE_nlo) * limcheck - buff)/onet
    yRnd_tmp(xRHO_nlo) = one - (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_sub124_wp_gq = xsect_nnlo_sub124_wp_gq(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    if (check_gqq_z124) then
       print *, 'hard gqq z124', Fint_gq_qg_sub_GQQ_z124_CheckLim(1)
    endif
    if (check_gqq_1z24) then
       print *, 'hard gqq 1z24', Fint_gq_qg_sub_GQQ_1z24_CheckLim(1)
    endif
    if (check_gqq_124) then
       print *, 'hard gqq 124', Fint_gq_qg_sub_GQQ_124_CheckLim(1)
    endif
    if (check_qqbg_z124) then
       print *, 'hard qqbg z124', Fint_gq_qg_sub_QQBG_z124_CheckLim(1)
       print *, 'scol42 qqbg z124', Fint_gq_qg_sub_QQBG_z124_CheckLim(6), one - Fint_gq_qg_sub_QQBG_z124_CheckLim(6)/Fint_gq_qg_sub_QQBG_z124_CheckLim(1)
    endif
    print *, ''

    
    print *, 'press a key to continue'
    read(*,*)

    return

  end function xsect_nnlo_checklim_sub124_wp_gq




  function xsect_nnlo_checklim_sub124_wp_qg(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_checklim_sub124_wp_qg,ndim,ncomp,userdata,iternumber,binarray(nobs)
    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)
    logical :: check_qgq_z124,check_qgq_1z24,check_qgq_124,check_qqbg_1z24

    bufftc = zero
    buffrr = zero

    check_qgq_z124= .false.
    check_qgq_1z24= .false.
    check_qgq_124 = .false.
    check_qqbg_1z24= .true.


    print *, "Checking limits for  sub124 terms in qg->qg channel"
    print *, " "

    print *, 'careful, if sping is included the soft and soft/collinear individually should not work properly'
    print *, 'since we assume cancellation at the sping level so we do not include them'
    print *, ''

    !-- soft limit
    yRnd_tmp = yRnd

    yRnd_tmp(xE_nlo) = (yRnd_tmp(xE_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_sub124_wp_qg = xsect_nnlo_sub124_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

       if (check_qgq_z124) then
          print *, 'hard gqq z124', Fint_qg_qg_sub_QGQ_z124_CheckLim(1)
       endif
       if (check_qgq_1z24) then
          print *, 'hard gqq 1z24', Fint_qg_qg_sub_QGQ_1z24_CheckLim(1)
       endif
       if (check_qgq_124) then
          print *, 'hard gqq 124', Fint_qg_qg_sub_QGQ_124_CheckLim(1)
       endif
       if (check_qqbg_1z24) then
          print *, 'hard qqbg z124', Fint_qg_qg_sub_Qqbg_1z24_CheckLim(1)
          print *, 'soft qqbg z124', Fint_qg_qg_sub_Qqbg_1z24_CheckLim(2), one + Fint_qg_qg_sub_Qqbg_1z24_CheckLim(2)/Fint_qg_qg_sub_Qqbg_1z24_CheckLim(1)
       endif

    print *, ''

    !-- c41 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xRHO_nlo) = (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_sub124_wp_qg = xsect_nnlo_sub124_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
    if (check_qgq_z124) then
       print *, 'hard gqq z124', Fint_qg_qg_sub_QGQ_z124_CheckLim(1)
    endif
    if (check_qgq_1z24) then
       print *, 'hard gqq z124', Fint_qg_qg_sub_QGQ_z124_CheckLim(1)
    endif
    if (check_qgq_124) then
       print *, 'hard gqq 124', Fint_qg_qg_sub_QGQ_124_CheckLim(1)
    endif
    if (check_qqbg_1z24) then
       print *, 'hard qqbg z124', Fint_qg_qg_sub_Qqbg_1z24_CheckLim(1)
       print *, 'c41  qqbg z124', Fint_qg_qg_sub_Qqbg_1z24_CheckLim(3), one + Fint_qg_qg_sub_Qqbg_1z24_CheckLim(3)/Fint_qg_qg_sub_Qqbg_1z24_CheckLim(1)
    endif
    print *, ' '


    !-- c42 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xRHO_nlo) = one - (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_sub124_wp_qg = xsect_nnlo_sub124_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    if (check_qgq_z124) then
       print *, 'hard gqq z124', Fint_qg_qg_sub_QGQ_z124_CheckLim(1)
       print *, 'c42  gqq z124', Fint_qg_qg_sub_QGQ_z124_CheckLim(2), one + Fint_qg_qg_sub_QGQ_z124_CheckLim(2)/Fint_qg_qg_sub_QGQ_z124_CheckLim(1)
    endif
    if (check_qgq_1z24) then
       print *, 'hard gqq 1z24', Fint_qg_qg_sub_QGQ_1z24_CheckLim(1)
       print *, 'c42  gqq 1z24', Fint_qg_qg_sub_QGQ_1z24_CheckLim(2), one + Fint_qg_qg_sub_QGQ_1z24_CheckLim(2)/Fint_qg_qg_sub_QGQ_1z24_CheckLim(1)
    endif
    if (check_qgq_124) then
       print *, 'hard gqq 124', Fint_qg_qg_sub_QGQ_124_CheckLim(1)
       print *, 'c42  gqq 124', Fint_qg_qg_sub_QGQ_124_CheckLim(2), one + Fint_qg_qg_sub_QGQ_124_CheckLim(2)/Fint_qg_qg_sub_QGQ_124_CheckLim(1)
    endif
    if (check_qqbg_1z24) then
       print *, 'hard qqbg z124', Fint_qg_qg_sub_Qqbg_1z24_CheckLim(1)
       print *, 'col2 qqbg z124', Fint_qg_qg_sub_Qqbg_1z24_CheckLim(4), one + Fint_qg_qg_sub_Qqbg_1z24_CheckLim(4)/Fint_qg_qg_sub_Qqbg_1z24_CheckLim(1)
    endif

    print *, ''


    !-- sc41 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xE_nlo) = (yRnd_tmp(xE_nlo) * limcheck - buff)/onet
    yRnd_tmp(xRHO_nlo) = (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_sub124_wp_qg = xsect_nnlo_sub124_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    if (check_qgq_z124) then
       print *, 'hard gqq z124', Fint_qg_qg_sub_QGQ_z124_CheckLim(1)
    endif
    if (check_qgq_1z24) then
       print *, 'hard gqq 1z24', Fint_qg_qg_sub_QGQ_1z24_CheckLim(1)
    endif
    if (check_qgq_124) then
       print *, 'hard gqq 124', Fint_qg_qg_sub_QGQ_124_CheckLim(1)
    endif
    if (check_qqbg_1z24) then
       print *, 'hard qqbg z124', Fint_qg_qg_sub_Qqbg_1z24_CheckLim(1)
       print *, 'scol41 qqbg z124', Fint_qg_qg_sub_Qqbg_1z24_CheckLim(5), one - Fint_qg_qg_sub_Qqbg_1z24_CheckLim(5)/Fint_qg_qg_sub_Qqbg_1z24_CheckLim(1)
    endif

    print *, ''


    !-- sc42 limit
    yRnd_tmp = yRnd

    yRnd_tmp(xE_nlo) = (yRnd_tmp(xE_nlo) * limcheck - buff)/onet
    yRnd_tmp(xRHO_nlo) = one - (yRnd_tmp(xRHO_nlo) * limcheck - buff)/onet
    xsect_nnlo_checklim_sub124_wp_qg = xsect_nnlo_sub124_wp_qg(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)

    if (check_qgq_z124) then
       print *, 'hard gqq z124', Fint_qg_qg_sub_QGQ_z124_CheckLim(1)
       print *, 'sc42  gqq z124', Fint_qg_qg_sub_QGQ_z124_CheckLim(2), one + Fint_qg_qg_sub_QGQ_z124_CheckLim(2)/Fint_qg_qg_sub_QGQ_z124_CheckLim(1)
    endif
    if (check_qgq_1z24) then
       print *, 'hard gqq 1z24', Fint_qg_qg_sub_QGQ_1z24_CheckLim(1)
       print *, 'sc42  gqq 1z24', Fint_qg_qg_sub_QGQ_1z24_CheckLim(2), one + Fint_qg_qg_sub_QGQ_1z24_CheckLim(2)/Fint_qg_qg_sub_QGQ_1z24_CheckLim(1)
    endif
    if (check_qgq_124) then
       print *, 'hard gqq 124', Fint_qg_qg_sub_QGQ_124_CheckLim(1)
       print *, 'sc42  gqq 124', Fint_qg_qg_sub_QGQ_124_CheckLim(2), one + Fint_qg_qg_sub_QGQ_124_CheckLim(2)/Fint_qg_qg_sub_QGQ_124_CheckLim(1)
    endif
    if (check_qqbg_1z24) then
       print *, 'hard qqbg z124', Fint_qg_qg_sub_Qqbg_1z24_CheckLim(1)
       print *, 'scol42 qqbg z124', Fint_qg_qg_sub_Qqbg_1z24_CheckLim(6), one - Fint_qg_qg_sub_Qqbg_1z24_CheckLim(6)/Fint_qg_qg_sub_Qqbg_1z24_CheckLim(1)
    endif
    print *, ''

    
    print *, 'press a key to continue'
    read(*,*)

    return

  end function xsect_nnlo_checklim_sub124_wp_qg

!!  function xsect_sub_checklim_2q(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
!!    integer :: xsect_sub_checklim_2q,ndim,ncomp,userdata,iternumber,binarray(nobs)
!!    real(dp15) :: yRnd(30),ff(1),vegasweight,yRnd_tmp(30)
!!
!!    bufftc = zero
!!    buffrr = zero
!!
!!    print *, 'careful, if sping is included the soft and soft/collinear individually should not work properly'
!!    print *, 'since we assume cancellation at the sping level so we do not include them'
!!    print *, ''
!!
!!    !-- soft limit
!!    yRnd_tmp = yRnd
!!
!!    yRnd_tmp(3) = (yRnd_tmp(3) * limcheck - buff)/onet
!!    xsect_sub_checklim_2q = xsect_2q_nnlo_sub_124(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
!!
!!    print *, 'hard 124', Fint_2q_sub_124_CheckLim(1)
!!    print *, 'soft 124', Fint_2q_sub_124_CheckLim(2), one + Fint_2q_sub_124_CheckLim(2)/Fint_2q_sub_124_CheckLim(1)
!!    print *, ''
!!
!!    !-- c41 limit
!!    yRnd_tmp = yRnd
!!
!!    yRnd_tmp(4) = (yRnd_tmp(4) * limcheck - buff)/onet
!!    xsect_sub_checklim_2q = xsect_2q_nnlo_sub_124(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
!!
!!    print *, 'hard 124', Fint_2q_sub_124_CheckLim(1)
!!    print *, 'c41  124', Fint_2q_sub_124_CheckLim(3), one + Fint_2q_sub_124_CheckLim(3)/Fint_2q_sub_124_CheckLim(1)
!!    print *, ''
!!    print *, 'hard z124', Fint_2q_sub_z124_CheckLim(1)
!!    print *, 'c41  z124', Fint_2q_sub_z124_CheckLim(2), one + Fint_2q_sub_z124_CheckLim(2)/Fint_2q_sub_z124_CheckLim(1)
!!    print *, ''
!!
!!    !-- c42 limit
!!    yRnd_tmp = yRnd
!!
!!    yRnd_tmp(4) = one - (yRnd_tmp(4) * limcheck - buff)/onet
!!    xsect_sub_checklim_2q = xsect_2q_nnlo_sub_124(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
!!
!!    print *, 'hard 124', Fint_2q_sub_124_CheckLim(1)
!!    print *, 'c42  124', Fint_2q_sub_124_CheckLim(4), one + Fint_2q_sub_124_CheckLim(4)/Fint_2q_sub_124_CheckLim(1)
!!    print *, ''
!!    print *, 'hard 1z4', Fint_2q_sub_1z24_CheckLim(1)
!!    print *, 'c42  1z4', Fint_2q_sub_1z24_CheckLim(2), one + Fint_2q_sub_1z24_CheckLim(2)/Fint_2q_sub_1z24_CheckLim(1)
!!    print *, ''
!!
!!    !-- sc41 limit
!!    yRnd_tmp = yRnd
!!
!!    yRnd_tmp(3) = (yRnd_tmp(3) * limcheck - buff)/onet
!!    yRnd_tmp(4) = (yRnd_tmp(4) * limcheck - buff)/onet
!!    xsect_sub_checklim_2q = xsect_2q_nnlo_sub_124(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
!!
!!    print *, 'hard 124', Fint_2q_sub_124_CheckLim(1)
!!    print *, 'sc41 124', Fint_2q_sub_124_CheckLim(5), one - Fint_2q_sub_124_CheckLim(5)/Fint_2q_sub_124_CheckLim(1)
!!    print *, ''
!!
!!    !-- sc42 limit
!!    yRnd_tmp = yRnd
!!
!!    yRnd_tmp(3) = (yRnd_tmp(3) * limcheck - buff)/onet
!!    yRnd_tmp(4) = one - (yRnd_tmp(4) * limcheck - buff)/onet
!!    xsect_sub_checklim_2q = xsect_2q_nnlo_sub_124(ndim,yRnd_tmp,ncomp,ff,userdata,vegasweight,iternumber)
!!
!!    print *, 'hard 124', Fint_2q_sub_124_CheckLim(1)
!!    print *, 'sc42 124', Fint_2q_sub_124_CheckLim(6), one - Fint_2q_sub_124_CheckLim(6)/Fint_2q_sub_124_CheckLim(1)
!!    print *, ''
!!
!!    
!!    print *, 'press a key to continue'
!!    read(*,*)
!!
!!    return
!!
!!  end function xsect_sub_checklim_2q
    
end module mod_check_lim_sub


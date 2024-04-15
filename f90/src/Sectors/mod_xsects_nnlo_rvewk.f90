!-- extra gluon emission, ewk loop
module mod_xsects_nnlo_rvewk
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_kinematics_nlo
  use mod_process
  use mod_auxfunctions
  use mod_aux_sectors
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_loop_ppll
  use mod_limits_rv
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintRVEWK_ns(3)
  real(dp), public, save :: FintRVEWK_gq(2)
  real(dp), public, save :: FintRVEWK_qg(2)
#endif

  private

  public :: xsect_nnlo_rvewk_is_ns
  public :: xsect_nnlo_rvewk_is_gq,xsect_nnlo_rvewk_is_qg

  !-- nf sectors
  public :: xsect_nnlo_rvewknf_is_ns
  public :: xsect_nnlo_rvewknf_is_gq,xsect_nnlo_rvewknf_is_qg

contains

  function xsect_nnlo_rvewk_is_ns(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rvewk_is_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintRV_ns(3),kin(3)
    real(dp)    :: res0(2,3),res1(2,3) !-- dn,qup and b
    real(dp)    :: res0red(2,3),res1red(2,3),restmp(2,3)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: limcol(3)

    xsect_nnlo_rvewk_is_ns = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    if ((xx(xE)*xx(xRHO).lt.buff_rv) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_rv)) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.false.)
    !-- these must go before the cut_histo routine
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_g]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintRV_ns(1) = zero

    else    

       call res_ewkloop_g_qqb(HardProc%AmpMom,res0,res1)
       call get_respdf(ns_lumi_splitb,1,1,HardProc,res1,respdf)
       
       respdf = respdf*HardProc%wgt
       
       kin(1) = respdf(1)
       FintRV_ns(1) = kin(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintRV_ns(2) = zero

    else

       call res_ewkloop_qqb(C1Lim%AmpMom,res0red,res1red)

       limcol = Cf*col_rvewk_is_ns(1,C1Lim)

       restmp(1,:) = res0red(1,:) * (limcol(1)*[Qdn2,Qup2,Qdn2] + limcol(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red(1,:) * limcol(3)
       !
       restmp(2,:) = res0red(2,:) * (limcol(1)*[Qdn2,Qup2,Qdn2] - limcol(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red(2,:) * limcol(3)

       call get_respdf(ns_lumi_splitb,1,1,C1Lim,restmp,respdf)
       respdf = -respdf*C1Lim%wgt

       FintRV_ns(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(3) = zero
       FintRV_ns(3) = zero

    else

       call res_ewkloop_qqb(C2Lim%AmpMom,res0red,res1red)

       limcol = Cf*col_rvewk_is_ns(2,C2Lim)

       restmp(1,:) = res0red(1,:) * (limcol(1)*[Qdn2,Qup2,Qdn2] - limcol(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red(1,:) * limcol(3)
       !
       restmp(2,:) = res0red(2,:) * (limcol(1)*[Qdn2,Qup2,Qdn2] + limcol(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red(2,:) * limcol(3)
       
       call get_respdf(ns_lumi_splitb,1,1,C2Lim,restmp,respdf)
       respdf = -respdf*C2Lim%wgt

       FintRV_ns(3) = respdf(1)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintRV_ns)
    
#if(_withchecks == 1)
    FintRVEWK_ns = FintRV_ns
#endif

  end function xsect_nnlo_rvewk_is_ns

  !--

  function xsect_nnlo_rvewk_is_gq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rvewk_is_gq
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintRV_gq(2),kin(2)
    real(dp)    :: res0(2,3),res1(2,3) !-- dn,qup and b
    real(dp)    :: res0red(2,3),res1red(2,3),restmp(2,3)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: limcol(3)

    xsect_nnlo_rvewk_is_gq = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    if (xx(xE)*xx(xRHO).lt.buff_r) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C1Lim=C1Lim,compute_etas=.false.)
    !-- these must go before the cut_histo routine
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    
    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintRV_gq(1) = zero

    else    

       call res_ewkloop_g_gq(HardProc%AmpMom(:,1:5),res0,res1)
       call get_respdf(gq_lumi_splitb,1,1,HardProc,res1,respdf)
       
       respdf = respdf*HardProc%wgt
       
       kin(1) = respdf(1)
       FintRV_gq(1) = kin(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintRV_gq(2) = zero

    else

       call res_ewkloop_qqb(C1Lim%AmpMom(:,1:4),res0red,res1red)

       limcol = tr*col_rvewk_is_qg(1,C1Lim)

       restmp(1,:) = res0red(1,:) * (limcol(1)*[Qdn2,Qup2,Qdn2] + limcol(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red(1,:) * limcol(3)
       !
       restmp(2,:) = res0red(2,:) * (limcol(1)*[Qdn2,Qup2,Qdn2] - limcol(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red(2,:) * limcol(3)
       
       call get_respdf(gq_lumi_splitb,1,1,C1Lim,restmp,respdf)
       respdf = -respdf*C1Lim%wgt

       FintRV_gq(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    ff(1) = sum(kin)
    call close_histo()
    
    call check_ff(ff,xx,FintRV_gq)
    
#if(_withchecks == 1)
    FintRVEWK_gq = FintRV_gq
#endif

  end function xsect_nnlo_rvewk_is_gq

  function xsect_nnlo_rvewk_is_qg(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rvewk_is_qg
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintRV_qg(2),kin(2)
    real(dp)    :: res0(2,3),res1(2,3) !-- dn,qup and b
    real(dp)    :: res0red(2,3),res1red(2,3),restmp(2,3)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: limcol(3)

    xsect_nnlo_rvewk_is_qg = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    if (xx(xE)*xx(xRHO).lt.buff_r) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C2Lim=C2Lim,compute_etas=.false.)
    !-- these must go before the cut_histo routine
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
    
    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintRV_qg(1) = zero

    else    

       call res_ewkloop_g_qg(HardProc%AmpMom(:,1:5),res0,res1)
       call get_respdf(qg_lumi_splitb,1,1,HardProc,res1,respdf)
       
       respdf = respdf*HardProc%wgt
       
       kin(1) = respdf(1)
       FintRV_qg(1) = kin(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(2) = zero
       FintRV_qg(2) = zero

    else

       call res_ewkloop_qqb(C2Lim%AmpMom(:,1:4),res0red,res1red)

       limcol = tr*col_rvewk_is_qg(2,C2Lim)

       restmp(1,:) = res0red(1,:) * (limcol(1)*[Qdn2,Qup2,Qdn2] - limcol(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red(1,:) * limcol(3)
       !
       restmp(2,:) = res0red(2,:) * (limcol(1)*[Qdn2,Qup2,Qdn2] + limcol(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red(2,:) * limcol(3)
       
       call get_respdf(qg_lumi_splitb,1,1,C2Lim,restmp,respdf)
       respdf = -respdf*C2Lim%wgt

       FintRV_qg(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    ff(1) = sum(kin)
    call close_histo()
    
    call check_ff(ff,xx,FintRV_qg)
    
#if(_withchecks == 1)
    FintRVEWK_qg = FintRV_qg
#endif

  end function xsect_nnlo_rvewk_is_qg


  !---------------------------------------------------------------------------!
  ! --------------              RV for nf sectors              -------------- !
  !---------------------------------------------------------------------------!
  function xsect_nnlo_rvewknf_is_ns(yRnd,ff,vegasweight)
    use mod_splittings_bare
    implicit none
    integer :: xsect_nnlo_rvewknf_is_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintRV_ns(3),kin(3)
    real(dp)    :: res_lo(2,2),res_nlo(2,2),respdf(ipdf)
    real(dp)    :: z,s5i

    xsect_nnlo_rvewknf_is_ns = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if ((xx(xE)*xx(xRHO).lt.buff_rv) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_rv)) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.false.)
    !-- these must go before the cut_histo routine
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_g]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintRV_ns(1) = zero

    else    

       call res_ewkloop_g_qqb_nf(HardProc%AmpMom,res_nlo)
       call get_respdf(ns_lumi,1,1,HardProc,res_nlo,respdf)
       
       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintRV_ns(1) = kin(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintRV_ns(2) = zero

    else

       call res_ewkloop_qqb_nf(C1Lim%AmpMom,res_lo)
       call get_respdf(ns_lumi,1,1,C1Lim,res_lo,respdf)

       z   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z)/(one-z))&
              * C1Lim%wgt

       FintRV_ns(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(3) = zero
       FintRV_ns(3) = zero

    else

       call res_ewkloop_qqb_nf(C2Lim%AmpMom,res_lo)
       call get_respdf(ns_lumi,1,1,C2Lim,res_lo,respdf)

       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z)/(one-z))&
              * C2Lim%wgt

       FintRV_ns(3) = respdf(1)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintRV_ns)

#if(_withchecks == 1)
    FintRVEWK_ns = FintRV_ns
#endif

  end function xsect_nnlo_rvewknf_is_ns

  function xsect_nnlo_rvewknf_is_gq(yRnd,ff,vegasweight)
    use mod_splittings_bare
    implicit none
    integer :: xsect_nnlo_rvewknf_is_gq
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintRV_gq(2),kin(2),z,s5i
    real(dp)    :: res_lo(2,2),res_nlo(2,2),respdf(ipdf)

    xsect_nnlo_rvewknf_is_gq = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    if (xx(xE)*xx(xRHO).lt.buff_r) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C1Lim=C1Lim,compute_etas=.false.)
    !-- these must go before the cut_histo routine
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    
    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintRV_gq(1) = zero

    else    

       call res_ewkloop_g_gq_nf(HardProc%AmpMom,res_nlo)
       call get_respdf(gq_lumi,1,1,HardProc,res_nlo,respdf)
       
       respdf = respdf*HardProc%wgt
       
       kin(1) = respdf(1)
       FintRV_gq(1) = kin(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintRV_gq(2) = zero

    else

       call res_ewkloop_qqb_nf(C1Lim%AmpMom,res_lo)
       call get_respdf(gq_lumi,1,1,C1Lim,res_lo,respdf)

       z   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pgq so that in the soft limit we do not have 1/[1-(1-x)]
       !-- Tr = Cf * aveqg/aveqq
       respdf = (-one)*respdf*(two/s5i)*(tr*Pgq_spav(z)/(one-z)) & 
              * C1Lim%wgt

       FintRV_gq(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintRV_gq)

#if(_withchecks == 1)
    FintRVEWK_gq = FintRV_gq
#endif

  end function xsect_nnlo_rvewknf_is_gq

  function xsect_nnlo_rvewknf_is_qg(yRnd,ff,vegasweight)
    use mod_splittings_bare
    implicit none
    integer :: xsect_nnlo_rvewknf_is_qg
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintRV_qg(2),kin(2),z,s5i
    real(dp)    :: res_lo(2,2),res_nlo(2,2),respdf(ipdf)

    xsect_nnlo_rvewknf_is_qg = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    if (xx(xE)*xx(xRHO).lt.buff_r) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C2Lim=C2Lim,compute_etas=.false.)
    !-- these must go before the cut_histo routine
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
    
    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintRV_qg(1) = zero

    else    

       call res_ewkloop_g_qg_nf(HardProc%AmpMom,res_nlo)
       call get_respdf(qg_lumi,1,1,HardProc,res_nlo,respdf)
       
       respdf = respdf*HardProc%wgt
       
       kin(1) = respdf(1)
       FintRV_qg(1) = kin(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(2) = zero
       FintRV_qg(2) = zero

    else

       call res_ewkloop_qqb_nf(C2Lim%AmpMom,res_lo)
       call get_respdf(qg_lumi,1,1,C2Lim,res_lo,respdf)

       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Tr*Pgq_spav(z)/(one-z))&
              * C2Lim%wgt

       FintRV_qg(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    ff(1) = sum(kin)
    call close_histo()
    
    call check_ff(ff,xx,FintRV_qg)

#if(_withchecks == 1)
    FintRVEWK_qg = FintRV_qg
#endif

  end function xsect_nnlo_rvewknf_is_qg

end module mod_xsects_nnlo_rvewk


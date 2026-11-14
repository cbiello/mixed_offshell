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
    real(dp)    :: res0_old(2,3),res1_old(2,3) !-- dn,qup and b
    real(dp)    :: res0red_old(2,3),res1red_old(2,3),restmp_old(2,3)
    real(dp)    :: res0(2,3),res1(-5:7,-5:7)
    real(dp)    :: res0red(-5:7,-5:7),res1red(-5:7,-5:7),res_tmp(-5:7,-5:7)
    real(dp)    :: Qi(-5:7,-5:7) !-- charge of the emitting leg
    real(dp)    :: respdf(ipdf)
    real(dp)    :: limcol_old(3),limcol(4)

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
    HardProc%ids(1:5) = [0,0,lep_ids(1),lep_ids(2),id_g]
    C1Lim%ids(1:4) = [0,0,lep_ids(1),lep_ids(2)]
    C2Lim%ids(1:4) = [0,0,lep_ids(1),lep_ids(2)]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintRV_ns(1) = zero

    else    

#if (_oldcode == 1)
       call res_ewkloop_g_qqb(HardProc%AmpMom,res0_old,res1_old)
       call get_respdf(ns_lumi_splitb,1,1,HardProc,res1_old,respdf)
#else
       call res_ewkloop_g_qqb_gen(HardProc%AmpMom,res0,res1)
       call get_respdf_gen(1,1,HardProc,res1,respdf)
#endif
       
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

#if (_oldcode == 1)
       call res_ewkloop_qqb(C1Lim%AmpMom,res0red_old,res1red_old)

       limcol_old = Cf*col_rvewk_is_ns(1,C1Lim)

       restmp_old(1,:) = res0red_old(1,:) * (limcol_old(1)*[Qdn2,Qup2,Qdn2] + limcol_old(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red_old(1,:) * limcol_old(3)
       !
       restmp_old(2,:) = res0red_old(2,:) * (limcol_old(1)*[Qdn2,Qup2,Qdn2] - limcol_old(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red_old(2,:) * limcol_old(3)

       call get_respdf(ns_lumi_splitb,1,1,C1Lim,restmp_old,respdf)
#else
       call res_ewkloop_qqb_gen(C1Lim%AmpMom,res0red,res1red)

       limcol = Cf*col_rvewk_is_ns_gen(1,C1Lim)

       Qi = spread(Q_IS,2,13) !-- Q_IS of leg 1
       res_tmp = res0red * (limcol(1)*Qi**2 + (limcol(2)*Q3 + limcol(3)*Q4)*Qi) + &
            res1red * limcol(4)

       call get_respdf_gen(1,1,C1Lim,res_tmp,respdf)
#endif
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

#if (_oldcode == 1)
       call res_ewkloop_qqb(C2Lim%AmpMom,res0red_old,res1red_old)

       limcol_old = Cf*col_rvewk_is_ns(2,C2Lim)

       restmp_old(1,:) = res0red_old(1,:) * (limcol_old(1)*[Qdn2,Qup2,Qdn2] - limcol_old(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red_old(1,:) * limcol_old(3)
       !
       restmp_old(2,:) = res0red_old(2,:) * (limcol_old(1)*[Qdn2,Qup2,Qdn2] + limcol_old(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red_old(2,:) * limcol_old(3)

       call get_respdf(ns_lumi_splitb,1,1,C2Lim,restmp_old,respdf)
#else
       call res_ewkloop_qqb_gen(C2Lim%AmpMom,res0red,res1red)

       limcol = Cf*col_rvewk_is_ns_gen(2,C2Lim)

       Qi = spread(Q_IS,1,13) !-- Q_IS of leg 2
       res_tmp = res0red * (limcol(1)*Qi**2 + (limcol(2)*Q3 + limcol(3)*Q4)*Qi) + &
            res1red * limcol(4)

       call get_respdf_gen(1,1,C2Lim,res_tmp,respdf)
#endif
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
    real(dp)    :: res1gen(-5:7,-5:7),res0gen(-5:7,-5:7),res_tmp(-5:7,-5:7),Qi(-5:7,-5:7),limgen(4)

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
    HardProc%ids(1:5) = [0,0,lep_ids(1),lep_ids(2),id_q]
    C1Lim%ids(1:4) = [0,0,lep_ids(1),lep_ids(2)]
    
    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintRV_gq(1) = zero

    else    

       
#if (_oldcode == 1)
       call res_ewkloop_g_gq(HardProc%AmpMom(:,1:5),res0,res1)
       call get_respdf(gq_lumi_splitb,1,1,HardProc,res1,respdf)
#else
       call res_ewkloop_g_gq_gen(HardProc%AmpMom(:,1:5),res0,res1gen)
       call get_respdf_gen(1,1,HardProc,res1gen,respdf)
#endif
       
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

#if (_oldcode == 1)
       call res_ewkloop_qqb(C1Lim%AmpMom(:,1:4),res0red,res1red)

       limcol = tr*col_rvewk_is_qg(1,C1Lim)

       restmp(1,:) = res0red(1,:) * (limcol(1)*[Qdn2,Qup2,Qdn2] + limcol(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red(1,:) * limcol(3)
       !
       restmp(2,:) = res0red(2,:) * (limcol(1)*[Qdn2,Qup2,Qdn2] - limcol(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red(2,:) * limcol(3)
       
       call get_respdf(gq_lumi_splitb,1,1,C1Lim,restmp,respdf)
#else
       call res_ewkloop_qqb_gen(C1Lim%AmpMom(:,1:4),res0gen,res1gen)
       limgen = tr*col_rvewk_is_qg_gen(1,C1Lim)
       Qi = spread(Q_IS,2,13) !-- Q_IS of leg 1
       res_tmp = res0gen * (limgen(1)*Qi**2 + (limgen(2)*Q3 + limgen(3)*Q4)*Qi) + &
            res1gen * limgen(4)
       res_tmp = transition('g -> q','none',res_tmp)
       call get_respdf_gen(1,1,C1Lim,res_tmp,respdf)
#endif
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
    real(dp)    :: res1gen(-5:7,-5:7),res0gen(-5:7,-5:7),res_tmp(-5:7,-5:7),Qi(-5:7,-5:7),limgen(4)

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
    HardProc%ids(1:5) = [0,0,lep_ids(1),lep_ids(2),id_q]
    C2Lim%ids(1:4) = [0,0,lep_ids(1),lep_ids(2)]
    
    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintRV_qg(1) = zero

    else    

       
#if (_oldcode == 1)
       call res_ewkloop_g_qg(HardProc%AmpMom(:,1:5),res0,res1)
       call get_respdf(qg_lumi_splitb,1,1,HardProc,res1,respdf)
#else
       call res_ewkloop_g_qg_gen(HardProc%AmpMom(:,1:5),res0,res1gen)
       call get_respdf_gen(1,1,HardProc,res1gen,respdf)
#endif
       
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

#if (_oldcode == 1)
       call res_ewkloop_qqb(C2Lim%AmpMom(:,1:4),res0red,res1red)

       limcol = tr*col_rvewk_is_qg(2,C2Lim)

       restmp(1,:) = res0red(1,:) * (limcol(1)*[Qdn2,Qup2,Qdn2] - limcol(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red(1,:) * limcol(3)
       !
       restmp(2,:) = res0red(2,:) * (limcol(1)*[Qdn2,Qup2,Qdn2] + limcol(2)*Q_lep*[Qdn,Qup,Qdn]) + &
            res1red(2,:) * limcol(3)
       
       call get_respdf(qg_lumi_splitb,1,1,C2Lim,restmp,respdf)
#else
       call res_ewkloop_qqb_gen(C2Lim%AmpMom(:,1:4),res0gen,res1gen)
       limgen = tr*col_rvewk_is_qg_gen(2,C2Lim)
       Qi = spread(Q_IS,1,13) !-- Q_IS of leg 2
       res_tmp = res0gen * (limgen(1)*Qi**2 + (limgen(2)*Q3 + limgen(3)*Q4)*Qi) + &
            res1gen * limgen(4)
       res_tmp = transition('none','g -> q',res_tmp)
       call get_respdf_gen(1,1,C2Lim,res_tmp,respdf)
#endif
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
    if (Vcharge /= 0) stop 'xsect_nnlo_rvewknf_is_ns: not implemented for W'

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
    if (Vcharge /= 0) stop 'xsect_nnlo_rvewknf_is_gq: not implemented for W'

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
    if (Vcharge /= 0) stop 'xsect_nnlo_rvewknf_is_qg: not implemented for W'

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


module mod_cut_histo_V_lla_nloqcd
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_histo
  use mod_auxfunctions
  use mod_kinfunctions
  use mod_cut_histo_aux
  use mod_process
  implicit none
  private

  public :: init_user_histo,cut_histo

contains

  !-- l l' + resolved photon, Z, W+ and W- (benchmark against POWHEG NSC-lnuaEW)
  !-- photon: pT > 50, |y| < 2.5, Delta R(photon, charged lepton) > R_phot, Delta R(photon, quark) > 0.4
  !-- W: W_ATLAS cuts on the bare lepton;  Z: mll in (qmin,qmax), ptlep_cut, ylep_cut, pT product
  !-- both: pT(l l') > ptll_cut (bare lepton pair)

  subroutine init_user_histo()

    call new_histo('rate',0._dp,1._dp,0.5_dp)       !-- 1
    call new_histo('pta', 50._dp,1000._dp,25._dp)   !-- 2
    call new_histo('mV',200._dp,1000._dp,20._dp)    !-- 3: mT (W) or mll (Z)

  end subroutine init_user_histo

  subroutine cut_histo(event)
    type(KinConfig), intent(inout) :: event
    real(dp) :: rap_boost
    real(dp) :: rec_mom(4,4)
    integer  :: nreco,ids(4),n_phot,i,nquarks,il,inu
    real(dp) :: p_gamma(4),p_quark(4),pt1,pt2,mV
    logical  :: leptons_recombined
    real(dp), parameter :: ptll_cut = 30._dp
    real(dp), parameter :: ptlep_geom = 35._dp, R_iso = 0.4_dp, pta_cut = 50._dp, ya_cut = 2.5_dp
    real(dp), parameter :: mT_min = 200._dp, mT_max = 5000._dp, ptl_cut_W = 65._dp, ptmiss_cut = 85._dp, yl_cut_W = 2.4_dp

    event%makecut = .true.
    leptons_recombined = .false.

    if (event%flag) return

    rap_boost = half*log(event%PartFrac(1)/event%PartFrac(2))

    !-- a photon within R_phot of a charged lepton is not resolved: reject
    call recombine_photons(event,rec_mom,nreco,ids,leptons_recombined,n_phot=n_phot)
    if (leptons_recombined) return
    if (n_phot < 1) return

    nquarks = 0
    do i = 5, event%npart
       if (event%ids(i) == id_a) then
          p_gamma = event%AmpMom(:,i)
       elseif (event%ids(i) == id_q) then
          p_quark = event%AmpMom(:,i)
          nquarks = nquarks + 1
       endif
    enddo

    if (nquarks == 1) then
       if (get_r(p_quark,p_gamma) < R_iso) return
    endif

    if (abs(get_y(p_gamma) + rap_boost) >= ya_cut) return
    if (get_pt(p_gamma) <= pta_cut) return

    if (get_pt(rec_mom(:,1)+rec_mom(:,2)) <= ptll_cut) return

    if (Vcharge == 0) then
       mV = sqrt(scr(rec_mom(:,1)+rec_mom(:,2),rec_mom(:,1)+rec_mom(:,2)))
       if (mV <= qmin .or. mV >= qmax) return
       pt1 = get_pt(rec_mom(:,1))
       pt2 = get_pt(rec_mom(:,2))
       if (pt1 <= ptlep_cut .or. pt2 <= ptlep_cut) return
       if (pt1*pt2 < ptlep_geom**2) return
       if (abs(get_y(rec_mom(:,1)) + rap_boost) >= ylep_cut) return
       if (abs(get_y(rec_mom(:,2)) + rap_boost) >= ylep_cut) return
    else
       if (event%ids(3) <= 0) stop 'V_lla_nloqcd: lepton ids not set'
       inu = merge(1,2,mod(event%ids(3),2) == 0)
       il  = 3 - inu
       pt1 = get_pt(rec_mom(:,il))
       pt2 = get_pt(rec_mom(:,inu))
       mV = sqrt(2*pt1*pt2*(one-cos(get_dphi(rec_mom(:,1),rec_mom(:,2)))))
       if (mV <= mT_min .or. mV >= mT_max) return
       if (pt1 <= ptl_cut_W) return
       if (pt2 <= ptmiss_cut) return
       if (pt1*pt2 < ptlep_geom**2) return
       if (abs(get_y(rec_mom(:,1)) + rap_boost) >= yl_cut_W) return
       if (abs(get_y(rec_mom(:,2)) + rap_boost) >= yl_cut_W) return
    endif

    if (dynscale) then
       event%mu = mV/two
    else
       event%mu = mu
    endif
    event%mur = event%mu * xmuR * (/1._dp/)
    event%muf = event%mu * xmuF * (/1._dp/)

    event%makecut = .false.

    obs(1) = 0.1_dp
    obs(2) = get_pt(p_gamma)
    obs(3) = mV

  end subroutine cut_histo

end module mod_cut_histo_V_lla_nloqcd

module mod_cut_histo_lla_nloqcd
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

  subroutine init_user_histo()

    !-- name, lower bin, upper bin, bin size

    !-- histograms
    call new_histo('rate',0._dp,1._dp,0.5_dp)       !-- 1
    call new_histo('mll',200._dp,1000._dp,20._dp)   !-- 2
    call new_histo('yll' ,-4._dp,4._dp,0.2_dp)      !-- 3
    call new_histo('ptll',50._dp,1000._dp,25._dp)   !-- 4
    call new_histo('pta', 50._dp,1000._dp,25._dp)   !-- 5
    call new_histo('ptlm',30._dp,980._dp,25._dp)    !-- 6
    call new_histo('ptlp',30._dp,980._dp,25._dp)    !-- 7
    call new_histo('ylm',-2.5_dp,2.5_dp,0.25_dp)    !-- 8
    call new_histo('ylp',-2.5_dp,2.5_dp,0.25_dp)    !-- 9
    call new_histo('ptj',0._dp,1000._dp,25._dp)     !-- 10

  end subroutine init_user_histo

  !-- order is f + f -> e- e+ [g a/q qb]
  subroutine cut_histo(event)
    type(KinConfig), intent(inout) :: event
    real(dp) :: rap_boost
    real(dp) :: rec_mom(4,4),pv(4)
    real(dp) :: mll2,mll,yll,ptll,ptlm,ptlp,ylp,ylm
    integer  :: nreco,ids(4), n_phot, i, nquarks
    real(dp) :: R_qa, pta, ya, p_gamma(4), p_quark(4)
    logical  :: leptons_recombined
    real(dp) :: ptlep_geom = 35._dp, R_iso = 0.4, pta_cut = 50._dp, ya_cut = 2.5_dp
    real(dp) :: ptj, yj, yj_cut = 2.5_dp, ptj_cut = 50._dp, pjet(4)
    logical  :: cut_on_jet = .true.

    !----------------------------------------------------------
    !-- do not touch this part
    !----------------------------------------------------------

    event%makecut = .true.
    leptons_recombined = .false.

    if (event%flag) return
    
    rap_boost = half*log(event%PartFrac(1)/event%PartFrac(2))

    !----------------------------------------------------------
    !-- user-defined cuts below
    !----------------------------------------------------------
    
    !-- dress leptons
    call recombine_photons(event,rec_mom,nreco,ids,leptons_recombined,n_phot=n_phot)

    !-- Recombination of leptons: reject event
    if(leptons_recombined) return

    !-- lepton-photon recombination: reject event
    if(n_phot < 1) return

    !-- Photon and quark momenta
    nquarks = 0
    do i= 5, event%npart
      if (event%ids(i) .eq. id_a) then
        p_gamma = event%AmpMom(:,i)
      elseif (event%ids(i) .eq. id_q) then
        p_quark = event%AmpMom(:,i)
        nquarks = nquarks + 1
      endif
    enddo

    !-- Photon isolation
    if(nquarks .eq. 1) then
      R_qa = get_r(p_quark,p_gamma)
      if(R_qa < R_iso) return 
    endif

    pv(:) = rec_mom(:,1) + rec_mom(:,2)

    !-- compute lepton observables and cut on them
    mll2 = scr(pv,pv)
    mll  = sqrt(mll2)
    if (mll.le.qmin .or. mll.ge.qmax) return

    ptlm  = get_pt(rec_mom(:,1))
    ptlp  = get_pt(rec_mom(:,2))

    !-- Cuts on lepton pTs + product cuts
    if (ptlm .le. ptlep_cut) return
    if (ptlp .le. ptlep_cut) return
    if (ptlp*ptlm < ptlep_geom**2) return

    !-- Cuts on lepton rapidities
    ylm = get_y(rec_mom(:,1)) + rap_boost
    if (abs(ylm).ge.ylep_cut) return

    ylp = get_y(rec_mom(:,2)) + rap_boost
    if (abs(ylp).ge.ylep_cut) return

    !-- Cuts on photon
    ya = get_y(p_gamma) + rap_boost
    if (abs(ya).ge.ya_cut) return

    pta = get_pt(p_gamma)
    if (pta.le.pta_cut) return

    !-- compute scale
    if (dynscale) then
       event%mu = mll/two
    else
       event%mu = mu
    endif
    event%mur = event%mu * xmuR * (/1._dp/)!,0.5_dp,2._dp/)
    event%muf = event%mu * xmuF * (/1._dp/)!,0.5_dp,2._dp/)

    yll = get_y(pv) + rap_boost
    ptll = get_pt(pv)

    pjet = event%AmpMom(:,1)+event%AmpMom(:,2)-pv-p_gamma
    ptj = get_pt(pjet)
    yj  = get_y(pjet) + rap_boost
    if(cut_on_jet) then
        if (ptj.le.ptj_cut) return
        if (abs(yj).ge.yj_cut) return
    endif

    !-- if I get here, all cuts are passed
    event%makecut = .false.

    !-- fill histograms
    obs(1)  = 0.1_dp !-- rate
    obs(2)  = mll
    obs(3)  = yll
    obs(4)  = ptll
    obs(5)  = pta
    obs(6)  = ptlm
    obs(7)  = ptlp
    obs(8)  = ylm
    obs(9)  = ylp
    obs(10) = ptj

    return

  end subroutine cut_histo
  
end module mod_cut_histo_lla_nloqcd

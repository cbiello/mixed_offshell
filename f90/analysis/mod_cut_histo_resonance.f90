module mod_cut_histo_resonance
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

    !-- Inclusive histograms
    call new_histo('rate',0._dp,1._dp,0.5_dp)      !-- 1
    call new_histo('mll',66._dp,116._dp,2._dp)     !-- 2
    call new_histo('yll' ,-5._dp,5._dp,0.1_dp)     !-- 3
    call new_histo('ptll',0._dp,100._dp,2._dp)     !-- 4
    call new_histo('ptlm',27._dp, 60._dp,1._dp)    !-- 5
    call new_histo('ptlp',27._dp, 60._dp,1._dp)    !-- 6

  end subroutine init_user_histo

  !-- order is f + f -> e- e+ [g a/q qb]
  subroutine cut_histo(event)
    type(KinConfig), intent(inout) :: event
    real(dp) :: rap_boost
    real(dp) :: rec_mom(4,4),rec_mom_lab(4,4),pv(4)
    real(dp) :: mll,yll,ptlm,ptlp,ylp,ylm,ptl1,ptl2,dYll,costh_star
    integer :: nreco,ids(4)
    logical :: leptons_recombined
    real(dp) :: ptlep_geom = 27._dp

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
    call recombine_photons(event,rec_mom,nreco,ids,leptons_recombined)
    !-- Recombination of leptons: reject event
    if(leptons_recombined) return

    pv(:) = rec_mom(:,1) + rec_mom(:,2)

    !-- compute lepton observables and cut on them
    mll = sqrt(scr(pv,pv))
    if (mll.le.qmin .or. mll.ge.qmax) return

    yll = get_y(pv) + rap_boost

    ptlm  = get_pt(rec_mom(:,1))
    ptlp  = get_pt(rec_mom(:,2))

    !-- Cuts on leptons pT. Product cuts
    if (ptlm .le. ptlep_cut) return
    if (ptlp .le. ptlep_cut) return

    if (ptlm .le. ptlm_min) return
    if (ptlm .ge. ptlm_max) return

    !-- Harder and softer leptons
    if(ptlm < ptlp) then
      ptl1 = ptlp
      ptl2 = ptlm
    else
      ptl1 = ptlm
      ptl2 = ptlp
    endif
    
    ylm = get_y(rec_mom(:,1)) + rap_boost
    if (abs(ylm).ge.ylep_cut) return

    ylp = get_y(rec_mom(:,2)) + rap_boost
    if (abs(ylp).ge.ylep_cut) return

    !-- compute scale
    if (dynscale) then
       event%mu = sqrt(mll**2 + get_pt(pv)**2)
    else
       event%mu = mu
    endif
    event%mur = event%mu * xmuR * (/1._dp/)!,0.5_dp,2._dp/)
    event%muf = event%mu * xmuF * (/1._dp/)!,0.5_dp,2._dp/)

    !-- if I get here, all cuts are passed
    event%makecut = .false.

    dYll = get_dY(rec_mom(:,1),rec_mom(:,2))

    !-- fill histograms
    obs(1) = 0.1_dp !-- rate
    obs(2) = mll
    obs(3) = yll
    obs(4) = get_pt(pv)
    obs(5) = ptlm
    obs(6) = ptlp

    return

  end subroutine cut_histo
  
end module mod_cut_histo_resonance

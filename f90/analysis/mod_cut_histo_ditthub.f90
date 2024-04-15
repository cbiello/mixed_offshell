module mod_cut_histo_ditthub
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
    call new_histo('rate',0._dp,1._dp,0.5_dp)       !-- 1
    call new_histo('mll',100._dp,1000._dp,10._dp)   !-- 2
    call new_histo('yll' ,-5._dp,5._dp,0.1_dp)      !-- 3

  end subroutine init_user_histo

  !-- order is f + f -> e- e+ [g a/q qb]
  subroutine cut_histo(event)
    type(KinConfig), intent(inout) :: event
    real(dp) :: rap_boost
    real(dp) :: rec_mom(4,4),rec_mom_lab(4,4),pv(4)
    real(dp) :: mll,yll,ptlm,ptlp,ylp,ylm
    integer  :: nreco,ids(4)
    logical  :: leptons_recombined

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
    call recombine_photons(event,rec_mom,nreco,ids,leptons_recombined,boost=rap_boost)
    !-- Recombination of leptons: reject event
    if(leptons_recombined) return

    pv(:) = rec_mom(:,1) + rec_mom(:,2)

    !-- compute lepton observables and cut on them
    mll = sqrt(scr(pv,pv))
    if (mll.le.qmin .or. mll.ge.qmax) return

    yll = get_y(pv) + rap_boost

    ptlm  = get_pt(rec_mom(:,1))
    if (ptlm .le. ptlep_cut) return
    
    ptlp  = get_pt(rec_mom(:,2))
    if (ptlp .le. ptlep_cut) return

    ylm = get_y(rec_mom(:,1)) + rap_boost
    if (abs(ylm).ge.ylep_cut) return

    ylp = get_y(rec_mom(:,2)) + rap_boost
    if (abs(ylp).ge.ylep_cut) return

    !-- compute scale
    if (dynscale) then
       event%mu = mll/two
    else
       event%mu = mu
    endif
    event%mur = event%mu * xmuR * (/1._dp/)!,0.5_dp,2._dp/)
    event%muf = event%mu * xmuF * (/1._dp/)!,0.5_dp,2._dp/)

    !-- if I get here, all cuts are passed
    event%makecut = .false.

    !-- fill histograms
    obs(1) = 0.1_dp !-- rate
    obs(2) = mll
    obs(3) = yll

    return

  end subroutine cut_histo

end module mod_cut_histo_ditthub

module mod_cut_histo_nf_200902229
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
    call new_histo('rate',0._dp,1._dp,0.5_dp)       !-- 1
    call new_histo('mll',100._dp,3000._dp,20._dp)   !-- 2
    call new_histo('yll' ,-5._dp,5._dp,0.1_dp)      !-- 3
    call new_histo('ptll',0._dp,3000._dp,10._dp)    !-- 4
    call new_histo('ptlm',10._dp, 3000._dp,10._dp)  !-- 5
    call new_histo('ptlp',10._dp, 3000._dp,10._dp)  !-- 6
    call new_histo('ptl1',10._dp, 3000._dp,10._dp)  !-- 7
    call new_histo('ptl2',10._dp, 3000._dp,10._dp)  !-- 8
    call new_histo('ylm',-5._dp, 5._dp,0.1_dp)      !-- 9
    call new_histo('ylp',-5._dp, 5._dp,0.1_dp)      !-- 10
    call new_histo('dRll',0._dp, 7._dp,0.1_dp)      !-- 11
    call new_histo('dYll',-5._dp, 5._dp,0.1_dp)     !-- 12
    call new_histo('dphill',0._dp,1._dp,0.05_dp)    !-- 13
    call new_histo('costhCS',-one,one,0.05_dp)      !-- 14
    call new_histo('sigma_b',0._dp,1._dp,0.5_dp)    !-- 15: sigma backward
    call new_histo('sigma_f',0._dp,1._dp,0.5_dp)    !-- 16: sigma forward
    call new_histo('mllB',200._dp,3000._dp,10._dp)  !-- 17
    call new_histo('mllF',200._dp,3000._dp,10._dp)  !-- 18

  end subroutine init_user_histo

  !-- order is f + f -> e- e+ [g a/q qb]
  subroutine cut_histo(event)
    type(KinConfig), intent(inout) :: event
    real(dp) :: rap_boost
    real(dp) :: rec_mom(4,4),rec_mom_lab(4,4),pv(4)
    real(dp) :: mll,yll,ptlm,ptlp,ylp,ylm,ptl1,ptl2,dYll,costh_star
    integer :: nreco,ids(4)
    logical :: leptons_recombined

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
       event%mu = mll/two
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
    obs(7) = ptl1
    obs(8) = ptl2
    obs(9) = ylm
    obs(10) = ylp
    obs(11) = get_R(rec_mom(:,1),rec_mom(:,2))
    obs(12) = dYll
    obs(13) = get_dphi(rec_mom(:,1),rec_mom(:,2))/pi    !-- TODO: check that it works for massive particles

    !-- need boosted momenta
    call boostz(event%PartFrac(1),event%PartFrac(2),2,rec_mom(:,1:2),rec_mom_lab(:,1:2))
    costh_star = get_costh_star(rec_mom_lab(:,1),rec_mom_lab(:,2))  !-- TODO: check that it works for massive particles
    obs(14) = costh_star

    !-- Needed for A_fb
    if(    (-1._dp.le.costh_star) .and. (costh_star.le. 0._dp)) then

      obs(15) = 0.1_dp
      obs(16) = -1._dp !-- reject

      obs(17) = mll
      obs(18) = -1._dp !-- reject
      
    elseif(( 0._dp.le.costh_star) .and. (costh_star.le. 1._dp)) then

      obs(15) = -1._dp !-- reject
      obs(16) = 0.1_dp

      obs(17) = -1._dp !-- reject
      obs(18) = mll

    endif

    return

  end subroutine cut_histo
  
end module mod_cut_histo_nf_200902229

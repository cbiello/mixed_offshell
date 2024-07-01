module mod_cut_histo_w_dflt
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
    call new_histo('rate',   0._dp,  1._dp, 0.5_dp)   !-- 1
    call new_histo('mT',    60._dp,100._dp, 0.5_dp)   !-- 2
    call new_histo('yln',   -5._dp,  5._dp, 0.1_dp)   !-- 3
    call new_histo('ptln',   0._dp,100._dp, 1._dp)    !-- 4
    call new_histo('ptl',   30._dp, 50._dp, 0.5_dp)   !-- 5
    call new_histo('ptn',   30._dp, 50._dp, 0.5_dp)   !-- 6
    call new_histo('yl',    -2.5_dp,2.5_dp, 0.1_dp)   !-- 7
    call new_histo('dphiln', 0._dp, 1._dp, 0.05_dp)   !-- 8

  end subroutine init_user_histo

  !-- order is f + f -> l nu [g a/q qb]
  subroutine cut_histo(event)
    type(KinConfig), intent(inout) :: event
    real(dp) :: rap_boost
    real(dp) :: rec_mom(4,4),pv(4)
    real(dp) :: massT,ptl,yl,ptmiss,dphiln
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

    !-- lepton pT
    ptl  = get_pt(rec_mom(:,1))
    if (ptl .le. ptlep_cut) return

    !-- pT-miss
    ptmiss = get_pT(rec_mom(:,2))
    if (ptmiss .le. ptmiss_cut) return

    !-- transverse mass
    dphiln = get_dphi(rec_mom(:,1),rec_mom(:,2))
    massT = sqrt(2*ptl*ptmiss*(one-cos(dphiln)))
    if (massT.le.mT_min) return

    yl = get_y(rec_mom(:,1)) + rap_boost
    if (abs(yl).ge.ylep_cut) return

    !-- compute scale
    if (dynscale) then
       event%mu = massT/two
    else
       event%mu = mu
    endif
    event%mur = event%mu * xmuR * (/1._dp/)!,0.5_dp,2._dp/)
    event%muf = event%mu * xmuF * (/1._dp/)!,0.5_dp,2._dp/)

    !-- if I get here, all cuts are passed
    event%makecut = .false.

    pv(:) = rec_mom(:,1) + rec_mom(:,2)

    !-- fill histograms
    obs(1) = 0.1_dp !-- rate
    obs(2) = massT
    obs(3) = get_y(pv) + rap_boost
    obs(4) = get_pt(pv)
    obs(5) = ptl
    obs(6) = ptmiss
    obs(7) = yl
    obs(8) = dphiln/pi

  end subroutine cut_histo
  
end module mod_cut_histo_w_dflt

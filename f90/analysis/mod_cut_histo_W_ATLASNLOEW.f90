module mod_cut_histo_W_ATLASNLOEW
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
    call new_histo('rate',0._dp,1._dp,0.8_dp)        !-- 1
    call new_histo('massT',200._dp,3000._dp,10._dp)  !-- 2
    call new_histo('pT-lep',5._dp, 1500._dp,5._dp)   !-- 3
    call new_histo('pT-nu',5._dp, 1500._dp,5._dp)    !-- 4
    call new_histo('mass',200._dp,3000._dp,10._dp)   !-- 5
    call new_histo('massbare',0._dp,3000._dp,10._dp) !-- 6
    call new_histo('y-lep',-4.9_dp, 4.9_dp,0.2_dp)   !-- 7

  end subroutine init_user_histo

  !-- order is f + f -> e- e+ [g a/q qb]
  subroutine cut_histo(event)
    type(KinConfig), intent(inout) :: event
    real(dp) :: rap_boost
    real(dp) :: rec_mom(4,4),rec_mom_lab(4,4),pv(4),pvbare(4)
    real(dp) :: mll,yll,ptlm,ptlp,ylp,ylm,ptl1,ptl2,dYll,costh_star
    integer :: nreco,ids(4)
    logical :: leptons_recombined
    real(dp) :: ptlep_geom = 35._dp, dYll_cut = 3.5_dp
    real(dp) :: ptl, ptmiss, dphiln, massT
    real(dp) :: mT_min = 200._dp, mT_max = 5000._dp
    real(dp) :: ptlep_cut = 65._dp, ptmiss_cut = 85._dp 
    real(dp) :: ylep_cut = 2.4_dp
    real(dp) :: invmass, invmassbare
    integer :: ihep, il, inu, iphot,n
    real(dp) :: ptphot, drlepphot
    real(dp) :: ptll
    real(dp) :: dr12, dr23, dr13
    integer :: icluster

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

    !-- dress leptons à la CB
    !first collect the particles
    n = event%npart
    il = 0
    inu = 0
    iphot = 0
#if (_Vcharge == 0)
    do ihep=3,n
       if(event%ids(ihep).eq.11) then
          il=ihep
       endif
       if(event%ids(ihep).eq.-11) then
          inu=ihep !nu in the Z case is e+
       endif
       if(event%ids(ihep).eq.22) then
          iphot=ihep
       endif
    enddo
#else
    do ihep=3,n
       if(abs(event%ids(ihep)).eq.11) then
          il=ihep
       endif
       if(abs(event%ids(ihep)).eq.12) then
          inu=ihep
       endif
       if(event%ids(ihep).eq.22) then
          iphot=ihep
       endif
    enddo
#endif

    !define the bare leptons
    rec_mom(:,1) = event%AmpMom(:,il)
    rec_mom(:,2) = event%AmpMom(:,inu)

    !-- invariant mass
    pvbare(:) = rec_mom(:,1) + rec_mom(:,2)
    invmassbare = sqrt(scr(pvbare,pvbare))
    
    dr12=get_r(event%AmpMom(:,il),event%AmpMom(:,inu))
    if(dr12.lt.0.3) return    

    if(iphot.gt.0) then
        dr13=get_r(event%AmpMom(:,il),event%AmpMom(:,iphot))
        dr23=get_r(event%AmpMom(:,inu),event%AmpMom(:,iphot))
        icluster=0
        if(dr13.lt.0.3 .and. dr23.ge.0.3) then
           icluster=1
        elseif(dr13.ge.0.3 .and. dr23.lt.0.3) then
           icluster=2
        elseif(dr13.lt.0.3 .and. dr23.lt.0.3) then
           if(dr13.lt.dr23) then
              icluster=1
           else
              icluster=2
           endif
        endif
        if(icluster.eq.1) rec_mom(:,1) = rec_mom(:,1) + event%AmpMom(:,iphot)
        if(icluster.eq.2) rec_mom(:,2) = rec_mom(:,2) + event%AmpMom(:,iphot)
    endif
    
    !-- invariant mass
    pv(:) = rec_mom(:,1) + rec_mom(:,2)
    invmass = sqrt(scr(pv,pv))
    if(invmass .lt. 220d0) return
    
    !-- compute the transverse mass from bare leptons
    ptl  = get_pt(rec_mom(:,1))
    ptmiss = get_pT(rec_mom(:,2))
    dphiln = get_dphi(rec_mom(:,1),rec_mom(:,2))
    massT = sqrt(2*ptl*ptmiss*(one-cos(dphiln)))
    if (massT.le.mT_min .or. massT.ge.mT_max ) return

    !-- Cuts on leptons pT. Product cuts
    if (ptl .le. ptlep_cut) return
    if (ptmiss .le. ptmiss_cut) return
    if (ptl*ptmiss < ptlep_geom**2) return

    !-- Cuts on rapidities
    ylm = get_y(rec_mom(:,1)) + rap_boost
    if (abs(ylm).ge.ylep_cut) return

    ylp = get_y(rec_mom(:,2)) + rap_boost
    if (abs(ylp).ge.ylep_cut) return

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

    !-- fill histograms
    obs(1) = 0.1_dp !-- rate
    obs(2) = massT
    obs(3) = ptl
    obs(4) = ptmiss
    obs(5) = invmass
    obs(6) = invmassbare
#if (_Vcharge == 1)
    obs(7) = ylp
#else
    obs(7) = ylm
#endif
    
    !-- need boosted momenta
    call boostz(event%PartFrac(1),event%PartFrac(2),2,rec_mom(:,1:2),rec_mom_lab(:,1:2))
    costh_star = get_costh_star(rec_mom_lab(:,1),rec_mom_lab(:,2))  !-- TODO: check that it works for massive particles

    return

  end subroutine cut_histo
  
end module mod_cut_histo_W_ATLASNLOEW 

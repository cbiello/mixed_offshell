module mod_cut_histo_W_ATLAS
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
    call new_histo('y-lep',-4.9_dp, 4.9_dp,0.2_dp)   !-- 6
    
    !call new_histo('yll' ,-5._dp,5._dp,0.1_dp)      !-- 3
    !call new_histo('ptll',0._dp,3000._dp,10._dp)    !-- 4
    !call new_histo('ptlp',10._dp, 3000._dp,10._dp)  !-- 6
    !call new_histo('ptl1',10._dp, 3000._dp,10._dp)  !-- 7
    !call new_histo('ptl2',10._dp, 3000._dp,10._dp)  !-- 8
    !call new_histo('ylm',-5._dp, 5._dp,0.1_dp)      !-- 9
    !call new_histo('ylp',-5._dp, 5._dp,0.1_dp)      !-- 10
    !call new_histo('dRll',0._dp, 7._dp,0.1_dp)      !-- 11
    !call new_histo('dYll',-5._dp, 5._dp,0.1_dp)     !-- 12
    !call new_histo('dphill',0._dp,1._dp,0.05_dp)    !-- 13
    !call new_histo('costhCS',-one,one,0.05_dp)      !-- 14
    !call new_histo('sigma_b',0._dp,1._dp,0.5_dp)    !-- 15: sigma backward
    !call new_histo('sigma_f',0._dp,1._dp,0.5_dp)    !-- 16: sigma forward
    !call new_histo('mllB',200._dp,5000._dp,10._dp)  !-- 17
    !call new_histo('mllF',200._dp,5000._dp,10._dp)  !-- 18
    
    !-- |\Delta_yll| cut histograms
    !call new_histo('dY_rate',0._dp,1._dp,0.5_dp)       !-- 19
    !call new_histo('dY_mll',200._dp,3000._dp,10._dp)   !-- 20
    !call new_histo('dY_yll' ,-5._dp,5._dp,0.1_dp)      !-- 21
    !call new_histo('dY_ptll',0._dp,3000._dp,10._dp)    !-- 22
    !call new_histo('dY_ptlm',10._dp, 3000._dp,10._dp)  !-- 23
    !call new_histo('dY_ptlp',10._dp, 3000._dp,10._dp)  !-- 24
    !call new_histo('dY_ptl1',10._dp, 3000._dp,10._dp)  !-- 25
    !call new_histo('dY_ptl2',10._dp, 3000._dp,10._dp)  !-- 26
    !call new_histo('dY_ylm',-5._dp, 5._dp,0.1_dp)      !-- 27
    !call new_histo('dY_ylp',-5._dp, 5._dp,0.1_dp)      !-- 28
    !call new_histo('dY_dRll',0._dp, 7._dp,0.1_dp)      !-- 29
    !call new_histo('dY_dYll',-5._dp, 5._dp,0.1_dp)     !-- 30
    !call new_histo('dY_dphill',0._dp,1._dp,0.05_dp)    !-- 31
    !call new_histo('dY_costhCS',-one,one,0.05_dp)      !-- 32
    !call new_histo('dY_sigmaf',0._dp,1._dp,0.5_dp)     !-- 33
    !call new_histo('dY_sigmab',0._dp,1._dp,0.5_dp)     !-- 34
    !call new_histo('mllB',200._dp,5000._dp,10._dp)     !-- 35
    !call new_histo('mllF',200._dp,5000._dp,10._dp)     !-- 36


  end subroutine init_user_histo



  !-- we need 
  !-- 1) |eta| < 2.4
  !-- 2) p_T^l > 65 GeV
  !-- 3) p_T^v > 85 GeV
  !-- 4) 200 < m_T^W < 5000 GeV [single differential measurements]
  !-- 5) 200 < m_T^W < 2000 GeV


  !-- order is f + f -> e- e+ [g a/q qb]
  subroutine cut_histo(event)
    type(KinConfig), intent(inout) :: event
    real(dp) :: rap_boost
    real(dp) :: rec_mom(4,4),rec_mom_lab(4,4),pv(4)
    real(dp) :: mll,yll,ptlm,ptlp,ylp,ylm,ptl1,ptl2,dYll,costh_star
    integer :: nreco,ids(4)
    logical :: leptons_recombined
    real(dp) :: ptlep_geom = 35._dp, dYll_cut = 3.5_dp
    real(dp) :: ptl, ptmiss, dphiln, massT
    real(dp) :: mT_min = 200._dp, mT_max = 5000._dp
    real(dp) :: ptlep_cut = 65._dp, ptmiss_cut = 85._dp 
    real(dp) :: ylep_cut = 2.4_dp
    real(dp) :: invmass

    !----------------------------------------------------------
    !-- consistency test
    !----------------------------------------------------------

#if (_Vcharge == 0)
    !print*, '******************************************************'
    !print*, '******* This is NOT the correct analys to run! *******'
    !print*, '******************************************************'
    !print*, ''
    !stop
#endif

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
    !call recombine_photons(event,rec_mom,nreco,ids,leptons_recombined)
    !-- Recombination of leptons: reject event
    !if(leptons_recombined) return

    rec_mom(:,1) = event%AmpMom(:,3)
    rec_mom(:,2) = event%AmpMom(:,4)
    !print*, 'ids= ', event%part(3), event%part(4)

    !invariant mass
    pv(:) = rec_mom(:,1) + rec_mom(:,2)
    invmass = sqrt(scr(pv,pv))

    !-- q qb -> l lx    
    !-- lepton pT
    if ( (event%part(3) .ne. 0) .and. (mod(event%part(3), 2) .eq. 0)) then
       !print*, 'ptl is second'
       ptl  = get_pt(rec_mom(:,2))
       ptmiss = get_pT(rec_mom(:,1))
    elseif ((event%part(3) .ne. 0) .and. (mod(event%part(3), 2) .ne. 0) ) then
       !print*, 'ptl is first'
       ptl  = get_pt(rec_mom(:,1))
       ptmiss = get_pT(rec_mom(:,2))
    endif
    
#if (_Vcharge == 0)
    !in case of Z
    !print*, 'ids= ', event%part(3), event%part(4)
    if( event%part(3) .eq. 11) then
       !print*, 'ptl is first'
       ptl = get_pt(rec_mom(:,1))
       ptmiss = get_pT(rec_mom(:,2))
    else
       !print*, 'ptl is second'
       ptl = get_pT(rec_mom(:,2))
       ptmiss = get_pT(rec_mom(:,1))
    endif
#endif
    
    !-- transverse mass
    dphiln = get_dphi(rec_mom(:,1),rec_mom(:,2))
    massT = sqrt(2*ptl*ptmiss*(one-cos(dphiln)))
    if (massT.le.mT_min .or. massT.ge.mT_max ) return

    !-- compute lepton observables and cut on them
    !mll = sqrt(scr(pv,pv))
    !if (mll.le.qmin .or. mll.ge.qmax) return
    
    !-- Momentum of the combined final state
    pv(:) = rec_mom(:,1) + rec_mom(:,2)

    !rapidity of the lepton-neutrino pair
    yll = get_y(pv) + rap_boost

    !-- Cuts on leptons pT. Product cuts
    if (ptl .le. ptlep_cut) return
    if (ptmiss .le. ptmiss_cut) return
    if (ptl*ptmiss < ptlep_geom**2) return

    ! !-- Harder and softer leptons
    ! if(ptlm < ptlp) then
    !   ptl1 = ptlp
    !   ptl2 = ptlm
    ! else
    !   ptl1 = ptlm
    !   ptl2 = ptlp
    ! endif
    
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

    !dYll = get_dY(rec_mom(:,1),rec_mom(:,2))

    !-- fill histograms
    obs(1) = 0.1_dp !-- rate
    obs(2) = massT
    obs(3) = ptl
    obs(4) = ptmiss
    obs(5) = invmass
#if (_Vcharge == 1)
    obs(6) = ylp
#else
    obs(6) = ylm
#endif
    
    !-- need boosted momenta
    call boostz(event%PartFrac(1),event%PartFrac(2),2,rec_mom(:,1:2),rec_mom_lab(:,1:2))
    costh_star = get_costh_star(rec_mom_lab(:,1),rec_mom_lab(:,2))  !-- TODO: check that it works for massive particles
    !obs(14) = costh_star

    !-- Needed for A_fb
    !if(    (-1._dp.le.costh_star) .and. (costh_star.le. 0._dp)) then
    !
    !  obs(15) = 0.1_dp
    !  obs(16) = -1._dp !-- reject
    ! 
    !  obs(17) = mll
    !  obs(18) = -1._dp !-- reject
      
    !elseif(( 0._dp.le.costh_star) .and. (costh_star.le. 1._dp)) then
    !
    !  obs(15) = -1._dp !-- reject
    !  obs(16) = 0.1_dp
    !
    !  obs(17) = -1._dp !-- reject
    !  obs(18) = mll
    !
    !endif

    !if(dYll < dYll_cut) then
    !  obs(19:36) = obs(1:18)
    !else
    !  obs(19:36) = 1E15_dp
    !endif

    return

  end subroutine cut_histo
  
end module mod_cut_histo_W_ATLAS 

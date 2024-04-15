module mod_cut_histo_fc_atlas
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
  !-- these cuts overrule the command-line ones
  real(dp) :: pt_lead  = 40._dp
  real(dp) :: pt_subl  = 30._dp
  real(dp) :: yl_cut   = 2.5_dp
  real(dp) :: dyll_cut = 3.5_dp
  private

  public :: init_user_histo,cut_histo

contains

  subroutine init_user_histo()

    !-- name, lower bin, upper bin, bin size
    !-- full
    call new_histo('rateFidNoDY',0._dp,1._dp,0.5_dp)   
    call new_histo('mllFidNoDY',200._dp,3000._dp,10._dp)
    call new_histo('yllFidNoDY',-5._dp,5._dp,0.1_dp)    
    call new_histo('ptllFidNoDY',0._dp,3000._dp,10._dp)  
    call new_histo('ptlmFidNoDY',10._dp, 3000._dp,10._dp)
    call new_histo('ptlpFidNoDY',10._dp, 3000._dp,10._dp)
    call new_histo('ylmFidNoDY',-5._dp, 5._dp,0.1_dp)    
    call new_histo('ylpFidNoDY',-5._dp, 5._dp,0.1_dp)    
    call new_histo('dRllFidNoDY',0._dp, 7._dp,0.1_dp)    
    call new_histo('dYllFidNoDY',-5._dp, 5._dp,0.1_dp)   
    call new_histo('dphillFidNoDY',0._dp,1._dp,0.05_dp)  
    call new_histo('costhCSFidNoDY',-one,one,0.05_dp)      
    !
    call new_histo('rateFidDY',0._dp,1._dp,0.5_dp)   
    call new_histo('mllFidDY',200._dp,3000._dp,10._dp)
    call new_histo('yllFidDY',-5._dp,5._dp,0.1_dp)    
    call new_histo('ptllFidDY',0._dp,3000._dp,10._dp)  
    call new_histo('ptlmFidDY',10._dp, 3000._dp,10._dp)
    call new_histo('ptlpFidDY',10._dp, 3000._dp,10._dp)
    call new_histo('ylmFidDY',-5._dp, 5._dp,0.1_dp)    
    call new_histo('ylpFidDY',-5._dp, 5._dp,0.1_dp)    
    call new_histo('dRllFidDY',0._dp, 7._dp,0.1_dp)    
    call new_histo('dYllFidDY',-5._dp, 5._dp,0.1_dp)   
    call new_histo('dphillFidDY',0._dp,1._dp,0.05_dp)  
    call new_histo('costhCSFidDY',-one,one,0.05_dp)      

    !-- 200-300
    call new_histo('rateFidNoDYw1',0._dp,1._dp,0.5_dp)   
    call new_histo('mllFidNoDYw1',200._dp,3000._dp,10._dp)
    call new_histo('yllFidNoDYw1',-5._dp,5._dp,0.1_dp)    
    call new_histo('ptllFidNoDYw1',0._dp,3000._dp,10._dp)  
    call new_histo('ptlmFidNoDYw1',10._dp, 3000._dp,10._dp)
    call new_histo('ptlpFidNoDYw1',10._dp, 3000._dp,10._dp)
    call new_histo('ylmFidNoDYw1',-5._dp, 5._dp,0.1_dp)    
    call new_histo('ylpFidNoDYw1',-5._dp, 5._dp,0.1_dp)    
    call new_histo('dRllFidNoDYw1',0._dp, 7._dp,0.1_dp)    
    call new_histo('dYllFidNoDYw1',-5._dp, 5._dp,0.1_dp)   
    call new_histo('dphillFidNoDYw1',0._dp,1._dp,0.05_dp)  
    call new_histo('costhCSFidNoDYw1',-one,one,0.05_dp)      
    !
    call new_histo('rateFidDYw1',0._dp,1._dp,0.5_dp)   
    call new_histo('mllFidDYw1',200._dp,3000._dp,10._dp)
    call new_histo('yllFidDYw1',-5._dp,5._dp,0.1_dp)    
    call new_histo('ptllFidDYw1',0._dp,3000._dp,10._dp)  
    call new_histo('ptlmFidDYw1',10._dp, 3000._dp,10._dp)
    call new_histo('ptlpFidDYw1',10._dp, 3000._dp,10._dp)
    call new_histo('ylmFidDYw1',-5._dp, 5._dp,0.1_dp)    
    call new_histo('ylpFidDYw1',-5._dp, 5._dp,0.1_dp)    
    call new_histo('dRllFidDYw1',0._dp, 7._dp,0.1_dp)    
    call new_histo('dYllFidDYw1',-5._dp, 5._dp,0.1_dp)   
    call new_histo('dphillFidDYw1',0._dp,1._dp,0.05_dp)  
    call new_histo('costhCSFidDYw1',-one,one,0.05_dp)      

    !-- 300-500
    call new_histo('rateFidNoDYw2',0._dp,1._dp,0.5_dp)   
    call new_histo('mllFidNoDYw2',200._dp,3000._dp,10._dp)
    call new_histo('yllFidNoDYw2',-5._dp,5._dp,0.1_dp)    
    call new_histo('ptllFidNoDYw2',0._dp,3000._dp,10._dp)  
    call new_histo('ptlmFidNoDYw2',10._dp, 3000._dp,10._dp)
    call new_histo('ptlpFidNoDYw2',10._dp, 3000._dp,10._dp)
    call new_histo('ylmFidNoDYw2',-5._dp, 5._dp,0.1_dp)    
    call new_histo('ylpFidNoDYw2',-5._dp, 5._dp,0.1_dp)    
    call new_histo('dRllFidNoDYw2',0._dp, 7._dp,0.1_dp)    
    call new_histo('dYllFidNoDYw2',-5._dp, 5._dp,0.1_dp)   
    call new_histo('dphillFidNoDYw2',0._dp,1._dp,0.05_dp)  
    call new_histo('costhCSFidNoDYw2',-one,one,0.05_dp)      
    !
    call new_histo('rateFidDYw2',0._dp,1._dp,0.5_dp)   
    call new_histo('mllFidDYw2',200._dp,3000._dp,10._dp)
    call new_histo('yllFidDYw2',-5._dp,5._dp,0.1_dp)    
    call new_histo('ptllFidDYw2',0._dp,3000._dp,10._dp)  
    call new_histo('ptlmFidDYw2',10._dp, 3000._dp,10._dp)
    call new_histo('ptlpFidDYw2',10._dp, 3000._dp,10._dp)
    call new_histo('ylmFidDYw2',-5._dp, 5._dp,0.1_dp)    
    call new_histo('ylpFidDYw2',-5._dp, 5._dp,0.1_dp)    
    call new_histo('dRllFidDYw2',0._dp, 7._dp,0.1_dp)    
    call new_histo('dYllFidDYw2',-5._dp, 5._dp,0.1_dp)   
    call new_histo('dphillFidDYw2',0._dp,1._dp,0.05_dp)  
    call new_histo('costhCSFidDYw2',-one,one,0.05_dp)      

    !-- 500-1500
    call new_histo('rateFidNoDYw3',0._dp,1._dp,0.5_dp)   
    call new_histo('mllFidNoDYw3',200._dp,3000._dp,10._dp)
    call new_histo('yllFidNoDYw3',-5._dp,5._dp,0.1_dp)    
    call new_histo('ptllFidNoDYw3',0._dp,3000._dp,10._dp)  
    call new_histo('ptlmFidNoDYw3',10._dp, 3000._dp,10._dp)
    call new_histo('ptlpFidNoDYw3',10._dp, 3000._dp,10._dp)
    call new_histo('ylmFidNoDYw3',-5._dp, 5._dp,0.1_dp)    
    call new_histo('ylpFidNoDYw3',-5._dp, 5._dp,0.1_dp)    
    call new_histo('dRllFidNoDYw3',0._dp, 7._dp,0.1_dp)    
    call new_histo('dYllFidNoDYw3',-5._dp, 5._dp,0.1_dp)   
    call new_histo('dphillFidNoDYw3',0._dp,1._dp,0.05_dp)  
    call new_histo('costhCSFidNoDYw3',-one,one,0.05_dp)      
    !
    call new_histo('rateFidDYw3',0._dp,1._dp,0.5_dp)   
    call new_histo('mllFidDYw3',200._dp,3000._dp,10._dp)
    call new_histo('yllFidDYw3',-5._dp,5._dp,0.1_dp)    
    call new_histo('ptllFidDYw3',0._dp,3000._dp,10._dp)  
    call new_histo('ptlmFidDYw3',10._dp, 3000._dp,10._dp)
    call new_histo('ptlpFidDYw3',10._dp, 3000._dp,10._dp)
    call new_histo('ylmFidDYw3',-5._dp, 5._dp,0.1_dp)    
    call new_histo('ylpFidDYw3',-5._dp, 5._dp,0.1_dp)    
    call new_histo('dRllFidDYw3',0._dp, 7._dp,0.1_dp)    
    call new_histo('dYllFidDYw3',-5._dp, 5._dp,0.1_dp)   
    call new_histo('dphillFidDYw3',0._dp,1._dp,0.05_dp)  
    call new_histo('costhCSFidDYw3',-one,one,0.05_dp)      

    !-- 1500 - inf
    call new_histo('rateFidNoDYw4',0._dp,1._dp,0.5_dp)   
    call new_histo('mllFidNoDYw4',200._dp,3000._dp,10._dp)
    call new_histo('yllFidNoDYw4',-5._dp,5._dp,0.1_dp)    
    call new_histo('ptllFidNoDYw4',0._dp,3000._dp,10._dp)  
    call new_histo('ptlmFidNoDYw4',10._dp, 3000._dp,10._dp)
    call new_histo('ptlpFidNoDYw4',10._dp, 3000._dp,10._dp)
    call new_histo('ylmFidNoDYw4',-5._dp, 5._dp,0.1_dp)    
    call new_histo('ylpFidNoDYw4',-5._dp, 5._dp,0.1_dp)    
    call new_histo('dRllFidNoDYw4',0._dp, 7._dp,0.1_dp)    
    call new_histo('dYllFidNoDYw4',-5._dp, 5._dp,0.1_dp)   
    call new_histo('dphillFidNoDYw4',0._dp,1._dp,0.05_dp)  
    call new_histo('costhCSFidNoDYw4',-one,one,0.05_dp)      
    !
    call new_histo('rateFidDYw4',0._dp,1._dp,0.5_dp)   
    call new_histo('mllFidDYw4',200._dp,3000._dp,10._dp)
    call new_histo('yllFidDYw4',-5._dp,5._dp,0.1_dp)    
    call new_histo('ptllFidDYw4',0._dp,3000._dp,10._dp)  
    call new_histo('ptlmFidDYw4',10._dp, 3000._dp,10._dp)
    call new_histo('ptlpFidDYw4',10._dp, 3000._dp,10._dp)
    call new_histo('ylmFidDYw4',-5._dp, 5._dp,0.1_dp)    
    call new_histo('ylpFidDYw4',-5._dp, 5._dp,0.1_dp)    
    call new_histo('dRllFidDYw4',0._dp, 7._dp,0.1_dp)    
    call new_histo('dYllFidDYw4',-5._dp, 5._dp,0.1_dp)   
    call new_histo('dphillFidDYw4',0._dp,1._dp,0.05_dp)  
    call new_histo('costhCSFidDYw4',-one,one,0.05_dp)      
    
    !call new_histo('mll',80._dp,120._dp,1._dp,cml=.true.) cml not yet implemented

  end subroutine init_user_histo

  !-- order is f + f -> e- e+ [g a/q qb]
  subroutine cut_histo(event)
    type(KinConfig), intent(inout) :: event
    real(dp) :: rap_boost
    real(dp) :: rec_mom(4,4),rec_mom_lab(4,4),pv(4)
    real(dp) :: mll,yll,ptlm,ptlp,ylp,ylm
    integer :: nreco,ids(4)
    logical :: leptons_recombined
    real(dp) :: dYll,dRll,ptv,costhCS

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

    if (max(ptlp,ptlm) .lt. pt_lead) return
    if (min(ptlp,ptlm) .lt. pt_subl) return
    
    ylm = get_y(rec_mom(:,1)) + rap_boost
    ylp = get_y(rec_mom(:,2)) + rap_boost

    if (abs(ylp).ge.yl_cut .or. abs(ylm).ge.yl_cut) return

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

    !-- prepare observables to plot
    ptv  = get_pt(pv)
    dYll = get_dY(rec_mom(:,1),rec_mom(:,2))
    dRll = get_R( rec_mom(:,1),rec_mom(:,2))

    call boostz(event%PartFrac(1),event%PartFrac(2),2,rec_mom(:,1:2),rec_mom_lab(:,1:2))
    costhCS = get_costh_star(rec_mom_lab(:,1),rec_mom_lab(:,2))  !-- TODO: check that it works for massive particles

    !-- fill histograms
    obs(1) = 0.1_dp !-- rate
    obs(2) = mll
    obs(3) = yll
    obs(4) = ptv
    obs(5) = ptlm
    obs(6) = ptlp
    obs(7) = ylm
    obs(8) = ylp
    obs(9) = dRll
    obs(10) = dYll
    obs(11) = get_dphi(rec_mom(:,1),rec_mom(:,2))/pi    !-- TODO: check that it works for massive particles
    obs(12) = costhCS

    if (abs(dYll).lt.dyll_cut) then
       obs(13:24) = obs(1:12)
    else
       obs(13:24) = -99._dp
    endif

    if (mll > 200._dp .and. mll < 300._dp) then
       obs(25:48) = obs(1:24)
    else
       obs(25:48) = -99._dp
    endif

    if (mll > 300._dp .and. mll < 500._dp) then
       obs(49:72) = obs(1:24)
    else
       obs(49:72) = -99._dp
    endif

    if (mll > 500._dp .and. mll < 1500._dp) then
       obs(73:96) = obs(1:24)
    else
       obs(73:96) = -99._dp
    endif

    if (mll > 1500._dp) then
       obs(97:120) = obs(1:24)
    else
       obs(97:120) = -99._dp
    endif
    
    return

  end subroutine cut_histo
  
end module mod_cut_histo_fc_atlas

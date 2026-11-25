module mod_cut_histo_W_incl
  !-- SETUP.md setup A: W -> l nu, dressed lepton, mT > 200 GeV, no other cut
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
    call new_histo('rate',0._dp,1._dp,0.5_dp)              !-- 1
    call new_histo('rate_window',0.5_dp,4.5_dp,1._dp)      !-- 2: mT windows 200-300, 300-500, 500-1500, >1500
    call new_histo('mT',200._dp,1000._dp,20._dp)           !-- 3
    call new_histo('mT_high',1000._dp,3000._dp,100._dp)    !-- 4
    call new_histo('mlnu',200._dp,1000._dp,20._dp)         !-- 5
    call new_histo('mlnu_high',1000._dp,3000._dp,100._dp)  !-- 6
    call new_histo('ptl',0._dp,1000._dp,20._dp)            !-- 7
    call new_histo('yl',-5._dp,5._dp,0.2_dp)               !-- 8
    call new_histo('ptlnu',0._dp,500._dp,10._dp)           !-- 9
    call new_histo('costhCS',-one,one,0.05_dp)             !-- 10: Collins-Soper angle of the charged lepton
  end subroutine init_user_histo

  subroutine cut_histo(event)
    type(KinConfig), intent(inout) :: event
    real(dp) :: rec_mom(4,4),pl(4),pnu(4),pl_lab(4),pnu_lab(4),mT
    integer  :: nreco,ids(4)
    logical  :: lept_recomb

    if (Vcharge == 0) stop 'W_incl: W only'
    event%makecut = .true.
    lept_recomb = .false.
    if (event%flag) return

    call recombine_photons(event,rec_mom,nreco,ids,lept_recomb)   !-- dresses the charged lepton
    if (mod(event%ids(3),2) == 0) then                            !-- leg 3 is the neutrino
       pnu = rec_mom(:,1); pl = rec_mom(:,2)
    else
       pl = rec_mom(:,1); pnu = rec_mom(:,2)
    endif

    mT = sqrt(2*get_pt(pl)*get_pt(pnu)*(one-cos(get_dphi(pl,pnu))))
    if (mT <= 200._dp) return

    event%mu = mu
    if (dynscale) event%mu = mT/two
    event%mur = event%mu * xmuR * (/1._dp/)
    event%muf = event%mu * xmuF * (/1._dp/)
    event%makecut = .false.

    obs(1) = 0.1_dp
    obs(2) = 1
    if (mT > 300._dp)  obs(2) = 2
    if (mT > 500._dp)  obs(2) = 3
    if (mT > 1500._dp) obs(2) = 4
    obs(3:4) = mT
    obs(5:6) = sqrt(scr(pl+pnu,pl+pnu))
    obs(7) = get_pt(pl)
    obs(8) = get_y(pl) + half*log(event%PartFrac(1)/event%PartFrac(2))   !-- lab frame
    obs(9) = get_pt(pl+pnu)
    call boostz(event%PartFrac(1),event%PartFrac(2),1,pl,pl_lab)      !-- lab frame: the sign of p_z(l nu)
    call boostz(event%PartFrac(1),event%PartFrac(2),1,pnu,pnu_lab)    !-- orients the angle
    obs(10) = get_costh_star(pl_lab,pnu_lab)

  end subroutine cut_histo

end module mod_cut_histo_W_incl

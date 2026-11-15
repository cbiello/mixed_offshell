module mod_cut_histo_V_lljj_4q
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

  !-- l l' + two jets, Z, W+ and W- (four-quark channels against POWHEG NSC-lnujjMIX)
  !-- jets: the two final partons, each with pT > ptj_cut, Delta R(j,j) > Rjj_cut
  !-- W: W_ATLAS cuts on the lepton;  Z: mll in (qmin,qmax), ptlep_cut, ylep_cut, pT product

  subroutine init_user_histo()

    call new_histo('rate',0._dp,1._dp,0.5_dp)       !-- 1
    call new_histo('ptj1',30._dp,1000._dp,25._dp)   !-- 2: hardest jet
    call new_histo('mjj',0._dp,2000._dp,50._dp)     !-- 3

  end subroutine init_user_histo

  subroutine cut_histo(event)
    type(KinConfig), intent(inout) :: event
    real(dp) :: rap_boost,p(4,6),pt1,pt2,mV
    integer  :: il,inu
    real(dp), parameter :: ptj_cut = 30._dp, Rjj_cut = 0.4_dp, ptlep_geom = 35._dp
    real(dp), parameter :: mT_min = 200._dp, mT_max = 5000._dp, ptl_cut_W = 65._dp, ptmiss_cut = 85._dp, yl_cut_W = 2.4_dp

    event%makecut = .true.

    if (event%flag) return
    if (event%npart /= 6) return

    rap_boost = half*log(event%PartFrac(1)/event%PartFrac(2))

    p = event%AmpMom(:,1:6)

    if (min(get_pt(p(:,5)),get_pt(p(:,6))) <= ptj_cut) return
    if (get_r(p(:,5),p(:,6)) <= Rjj_cut) return

    if (Vcharge /= 0) then
       if (event%ids(3) <= 0) stop 'V_lljj_4q: lepton ids not set'
       inu = merge(3,4,mod(event%ids(3),2) == 0)
       il  = 7 - inu
    endif

    if (Vcharge == 0) then
       mV = sqrt(scr(p(:,3)+p(:,4),p(:,3)+p(:,4)))
       if (mV <= qmin .or. mV >= qmax) return
       pt1 = get_pt(p(:,3))
       pt2 = get_pt(p(:,4))
       if (pt1 <= ptlep_cut .or. pt2 <= ptlep_cut) return
       if (pt1*pt2 < ptlep_geom**2) return
       if (abs(get_y(p(:,3)) + rap_boost) >= ylep_cut) return
       if (abs(get_y(p(:,4)) + rap_boost) >= ylep_cut) return
    else
       pt1 = get_pt(p(:,il))
       pt2 = get_pt(p(:,inu))
       mV = sqrt(2*pt1*pt2*(one-cos(get_dphi(p(:,3),p(:,4)))))
       if (mV <= mT_min .or. mV >= mT_max) return
       if (pt1 <= ptl_cut_W) return
       if (pt2 <= ptmiss_cut) return
       if (pt1*pt2 < ptlep_geom**2) return
       if (abs(get_y(p(:,3)) + rap_boost) >= yl_cut_W) return
       if (abs(get_y(p(:,4)) + rap_boost) >= yl_cut_W) return
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
    obs(2) = max(get_pt(p(:,5)),get_pt(p(:,6)))
    obs(3) = sqrt(scr(p(:,5)+p(:,6),p(:,5)+p(:,6)))

  end subroutine cut_histo

end module mod_cut_histo_V_lljj_4q

module mod_cut_histo_V_llj_nloew
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

  !-- l l' + jet, Z, W+ and W- (benchmark against POWHEG NSC-lnujQCD)
  !-- photon: recombined with the nearest charged lepton or final quark if Delta R < R_phot; never with a gluon
  !-- jet: a final quark (dressed) or gluon with pT > ptj_cut; a photon is not a jet
  !-- W: W_ATLAS cuts on the dressed lepton;  Z: mll in (qmin,qmax), ptlep_cut, ylep_cut, pT product

  subroutine init_user_histo()

    call new_histo('rate',0._dp,1._dp,0.5_dp)       !-- 1
    call new_histo('ptj', 30._dp,1000._dp,25._dp)   !-- 2: hardest jet
    call new_histo('mV',200._dp,1000._dp,20._dp)    !-- 3: mT (W) or mll (Z)

  end subroutine init_user_histo

  subroutine cut_histo(event)
    type(KinConfig), intent(inout) :: event
    real(dp) :: rap_boost,p(4,6),r,rmin,ptj,pt1,pt2,mV
    integer  :: i,ia,imin,il,inu
    logical  :: charged(6)
    real(dp), parameter :: ptj_cut = 30._dp, ptlep_geom = 35._dp
    real(dp), parameter :: mT_min = 200._dp, mT_max = 5000._dp, ptl_cut_W = 65._dp, ptmiss_cut = 85._dp, yl_cut_W = 2.4_dp

    event%makecut = .true.

    if (event%flag) return

    rap_boost = half*log(event%PartFrac(1)/event%PartFrac(2))

    p(:,1:event%npart) = event%AmpMom(:,1:event%npart)

    !-- particles a photon can be recombined with
    charged = .false.
    charged(3:4) = .true.
    if (Vcharge /= 0) then
       if (event%ids(3) <= 0) stop 'V_llj_nloew: lepton ids not set'
       inu = merge(3,4,mod(event%ids(3),2) == 0)
       il  = 7 - inu
       charged(inu) = .false.
    endif
    ia = 0
    do i = 5, event%npart
       if (event%ids(i) == id_q) charged(i) = .true.
       if (event%ids(i) == id_a) ia = i
    enddo

    if (ia > 0) then
       rmin = R_phot
       imin = 0
       do i = 3, event%npart
          if (.not. charged(i)) cycle
          r = get_r(p(:,i),p(:,ia))
          if (r < rmin) then
             rmin = r
             imin = i
          endif
       enddo
       if (imin > 0) p(:,imin) = p(:,imin) + p(:,ia)
    endif

    ptj = zero
    do i = 5, event%npart
       if (event%ids(i) == id_q .or. event%ids(i) == id_g) ptj = max(ptj,get_pt(p(:,i)))
    enddo
    if (ptj <= ptj_cut) return

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
    obs(2) = ptj
    obs(3) = mV

  end subroutine cut_histo

end module mod_cut_histo_V_llj_nloew

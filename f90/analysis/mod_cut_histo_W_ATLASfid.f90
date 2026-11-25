module mod_cut_histo_W_ATLASfid
  !-- SETUP.md setup B, ATLAS fiducial region (arXiv:2502.21088): dressed lepton pT > 65 GeV, |eta_l| < 2.4,
  !-- pT_nu > 85 GeV (no rapidity cut on the neutrino), 200 < mT < 5000 GeV
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
    call new_histo('rate',0._dp,1._dp,0.5_dp)               !-- 1
    call new_histo('mT',200._dp,5000._dp,10._dp)            !-- 2: regroup into the 12 ATLAS bins (edges from HEPData)
    call new_histo('etal_mT200-300',0._dp,2.4_dp,0.2_dp)    !-- 3-7: |eta_l| per mT bin, ATLAS Table 1
    call new_histo('etal_mT300-425',0._dp,2.4_dp,0.2_dp)    !-- 4
    call new_histo('etal_mT425-600',0._dp,2.4_dp,0.4_dp)    !-- 5
    call new_histo('etal_mT600-900',0._dp,2.4_dp,0.4_dp)    !-- 6
    call new_histo('etal_mT900-2000',0._dp,2.4_dp,0.6_dp)   !-- 7
    call new_histo('ptl',0._dp,2500._dp,10._dp)             !-- 8
    call new_histo('ptnu',0._dp,2500._dp,10._dp)            !-- 9
    call new_histo('etal',-2.4_dp,2.4_dp,0.1_dp)            !-- 10
    call new_histo_required()                               !-- 11: pta or ptj with -require photon/jet
  end subroutine init_user_histo

  subroutine cut_histo(event)
    type(KinConfig), intent(inout) :: event
    real(dp), parameter :: mTbin(6) = [200._dp, 300._dp, 425._dp, 600._dp, 900._dp, 2000._dp]
    real(dp) :: pl(4),pnu(4),pl_lab(4),mT,ptl,ptnu,etal,ptobj
    integer  :: i
    logical  :: pass

    if (Vcharge == 0) stop 'W_ATLASfid: W only'
    event%makecut = .true.
    if (event%flag) return

    call W_leptons(event,pl,pnu,pass,ptobj)                       !-- dressed charged lepton, neutrino, -require
    if (.not. pass) return
    call boostz(event%PartFrac(1),event%PartFrac(2),1,pl,pl_lab)  !-- eta of the massive dressed lepton: lab frame

    ptl  = get_pt(pl)
    ptnu = get_pt(pnu)
    etal = get_eta(pl_lab)
    mT   = sqrt(2*ptl*ptnu*(one-cos(get_dphi(pl,pnu))))
    if (ptl <= 65._dp .or. abs(etal) >= 2.4_dp .or. ptnu <= 85._dp) return
    if (mT <= 200._dp .or. mT >= 5000._dp) return

    event%mu = mu
    if (dynscale) event%mu = mT/two
    event%mur = event%mu * xmuR * (/1._dp/)
    event%muf = event%mu * xmuF * (/1._dp/)
    event%makecut = .false.

    obs(1) = 0.1_dp
    obs(2) = mT
    do i = 1,5                                                   !-- outside its mT bin: -1, the underflow
       obs(2+i) = -1._dp
       if (mT > mTbin(i) .and. mT <= mTbin(i+1)) obs(2+i) = abs(etal)
    enddo
    obs(8)  = ptl
    obs(9)  = ptnu
    obs(10) = etal
    if (require /= 'none') obs(11) = ptobj

  end subroutine cut_histo

end module mod_cut_histo_W_ATLASfid

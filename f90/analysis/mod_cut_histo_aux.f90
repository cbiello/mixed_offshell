module mod_cut_histo_aux
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_auxfunctions
  use mod_kinfunctions
  use mod_process
  use mod_histo
  implicit none
  private

  public :: recombine_photons,W_leptons,new_histo_required

contains

  !-- W -> l nu: charged lepton and neutrino after the photon treatment chosen by -require
  !-- none:   a photon within R_phot is added to the nearest of charged lepton and final quark (collinear safe: a photon
  !--         collinear to a quark goes to the quark, PLAN team question 10)
  !-- photon: resolved photon (POWHEG NSC-lnuaEW): pT > 50, |y| < 2.5, R(photon,l) > R_phot, R(photon,quark) > 0.4,
  !--         pT(l nu) > 30; ptobj = pT of the photon
  !-- jet:    (POWHEG NSC-lnujQCD) a photon within R_phot is added to the nearest of charged lepton and final quark;
  !--         jet = final quark or gluon with pT > 30; ptobj = pT of the hardest one
  subroutine W_leptons(event,pl,pnu,pass,ptobj)
    type(KinConfig), intent(in) :: event
    real(dp), intent(out) :: pl(4),pnu(4),ptobj
    logical,  intent(out) :: pass
    real(dp) :: p(4,6),r,rmin,rap_boost
    integer  :: i,il,inu,ia,imin
    logical  :: charged(6)
    real(dp), parameter :: pta_cut = 50._dp, ya_cut = 2.5_dp, R_iso = 0.4_dp, ptlnu_cut = 30._dp, ptj_cut = 30._dp

    if (event%ids(3) <= 0) stop 'W_leptons: lepton ids not set'
    inu = merge(3,4,mod(event%ids(3),2) == 0)
    il  = 7 - inu
    p(:,1:event%npart) = event%AmpMom(:,1:event%npart)
    rap_boost = half*log(event%PartFrac(1)/event%PartFrac(2))
    pass  = .false.
    ptobj = zero
    pl    = p(:,il)
    pnu   = p(:,inu)

    ia = 0
    charged = .false.
    charged(il) = .true.
    do i = 5,event%npart
       if (event%ids(i) == id_a) ia = i
       if (event%ids(i) == id_q .and. require /= 'photon') charged(i) = .true.
    enddo

    imin = 0                                  !-- nearest charged particle within R_phot of the photon
    if (ia > 0) then
       rmin = R_phot
       do i = 3,event%npart
          if (.not. charged(i)) cycle
          r = get_r(p(:,i),p(:,ia))
          if (r < rmin) then
             rmin = r
             imin = i
          endif
       enddo
       if (rec_phot_beam .and. abs(get_y(p(:,ia)) + rap_boost) > ya_beam) imin = 0
    endif

    if (require == 'photon') then
       if (ia == 0 .or. imin > 0) return
       if (get_pt(p(:,ia)) <= pta_cut .or. abs(get_y(p(:,ia)) + rap_boost) >= ya_cut) return
       do i = 5,event%npart
          if (event%ids(i) == id_q) then
             if (get_r(p(:,i),p(:,ia)) < R_iso) return
          endif
       enddo
       if (get_pt(pl+pnu) <= ptlnu_cut) return
       ptobj = get_pt(p(:,ia))
    else
       if (imin > 0) p(:,imin) = p(:,imin) + p(:,ia)
       pl = p(:,il)
       if (require == 'jet') then
          do i = 5,event%npart
             if (event%ids(i) == id_q .or. event%ids(i) == id_g) ptobj = max(ptobj,get_pt(p(:,i)))
          enddo
          if (ptobj <= ptj_cut) return
       endif
    endif
    pass = .true.

  end subroutine W_leptons

  !-- histogram of the object asked by -require; booked last, filled with ptobj of W_leptons
  subroutine new_histo_required()
    if (require == 'photon') call new_histo('pta',50._dp,1000._dp,25._dp)
    if (require == 'jet')    call new_histo('ptj',30._dp,1030._dp,25._dp)
  end subroutine new_histo_required

  !-- recombine photons, return list of final state momenta and the number of
  !-- total reconstructed particles
  !-- only works with up to 1 extra photon
  subroutine recombine_photons(event,recomb_mom,nreco,ids,lept_recomb,boost,n_phot)
    type(KinConfig), intent(in) :: event
    real(dp), intent(out) :: recomb_mom(4,4)
    integer, intent(out) :: nreco,ids(4)
    logical, intent(inout) :: lept_recomb
    integer, optional, intent(out) :: n_phot
    real(dp), optional, intent(in) :: boost
    real(dp) :: dressed_lep(4,2), lboost
    integer :: i,n,k,mygamma,nphot

    n = event%npart

    k = 3 !-- first 2 are always l, lb
    nphot   = 0
    mygamma = 0
    lboost = zero !-- longitudinal boost

    if(present(boost)) lboost = boost

    do i = 5,n !-- run over non-lepton momenta
       if (event%ids(i) .ne. id_a) then
          recomb_mom(:,k) = event%AmpMom(:,i)
          ids(k) = event%ids(i)
          k = k + 1
       else
          mygamma = i
          if (event%ids(3) == -event%ids(4)) then
             call recombine_lla(event%AmpMom(:,3),event%AmpMom(:,4),event%AmpMom(:,i),&
              dressed_lep,nphot,lept_recomb,lboost)
          else
             call recombine_la([event%ids(3),event%ids(4)],event%AmpMom(:,3),event%AmpMom(:,4),event%AmpMom(:,i),&
              dressed_lep,nphot,lept_recomb,lboost) 
          endif 
          !-- Recombination of leptons: reject event
          if(lept_recomb) return
          recomb_mom(:,1) = dressed_lep(:,1)
          recomb_mom(:,2) = dressed_lep(:,2)
          if (nphot == 1) then
             recomb_mom(:,k) = event%AmpMom(:,i) !-- if I have resolved photon, load it
             ids(k) = event%ids(i)
             k = k+1
          endif
       endif
    enddo

    !-- if mygamma == 0, need to load lepton momenta
    if (mygamma == 0) then
       recomb_mom(:,1) = event%AmpMom(:,3)
       recomb_mom(:,2) = event%AmpMom(:,4)
    endif

    ids(1) =  id_el
    ids(2) = -id_el

    nreco = k - 1

    if(present(n_phot)) n_phot = nphot

  end subroutine recombine_photons

  subroutine recombine_lla(pl1,pl2,pa,dress_lep,nphot,lept_recomb,lboost)
    implicit none
    real(dp), intent(in) :: pl1(4),pl2(4),pa(4),lboost
    real(dp), intent(out) :: dress_lep(4,2)
    logical,  intent(inout) :: lept_recomb
    integer,  intent(out) :: nphot
    real(dp) :: r_vec(3)
    integer  :: irmin(1)
    logical :: phot_in_beam

    phot_in_beam = .false.
    dress_lep(:,1) = pl1
    dress_lep(:,2) = pl2

    r_vec  = [get_r(pl1,pa),get_r(pl2,pa),get_r(pl1,pl2)]
    irmin = minloc(r_vec)

    !-- Recombination of leptons: reject event
    if((irmin(1) == 3) .and. (r_vec(3) .lt. R_phot)) then
      lept_recomb = .true.
      return
    endif

    if(rec_phot_beam) phot_in_beam = abs((get_y(pa) + lboost)) > ya_beam

    if((r_vec(irmin(1)).lt.R_phot) .and. (.not.phot_in_beam)) then
       dress_lep(:,irmin(1)) = dress_lep(:,irmin(1)) + pa
       nphot = 0
    else
       nphot = 1
    endif

    return
  end subroutine recombine_lla


  subroutine recombine_la(idx,pl1,pl2,pa,dress_lep,nphot,lept_recomb,lboost)
  implicit none
  integer,  intent(in) :: idx(2)
  real(dp), intent(in) :: pl1(4),pl2(4),pa(4),lboost
  real(dp), intent(out) :: dress_lep(4,2)
  logical,  intent(inout) :: lept_recomb
  integer,  intent(out) :: nphot

  real(dp) :: r_vec(2)
  integer  :: irmin(1)
  logical :: phot_in_beam
  logical :: nu1, nu2

  phot_in_beam = .false.
  dress_lep(:,1) = pl1
  dress_lep(:,2) = pl2

  !-- identify neutrino (your rule)
  nu1 = (idx(1) .gt. 0) .and. (mod(idx(1),2) .eq. 0)
  nu2 = (idx(2) .gt. 0) .and. (mod(idx(2),2) .eq. 0)

  !--------------------------------------------------
  ! Only one particle can recombine (the non-neutrino)
  !--------------------------------------------------
  if (nu1) then
     !-- pl1 is neutrino → only pl2 can recombine
     r_vec = [1d10, get_r(pl2,pa)]
  else
     !-- pl2 is neutrino → only pl1 can recombine
     r_vec = [get_r(pl1,pa), 1d10]
  endif

  irmin = minloc(r_vec)

  if(rec_phot_beam) phot_in_beam = abs((get_y(pa) + lboost)) > ya_beam

  if((r_vec(irmin(1)).lt.R_phot) .and. (.not.phot_in_beam)) then
     dress_lep(:,irmin(1)) = dress_lep(:,irmin(1)) + pa
     nphot = 0
  else
     nphot = 1
  endif

  return
end subroutine recombine_la


  
end module mod_cut_histo_aux

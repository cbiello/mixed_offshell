!-- observables and cumulative observables. Format must be !-- xxx = number, with these spaces
!-- Nobs = 17
!-- Nnoncum = 16
!-- cutname = 'fid1' !-- maximum 50 characters

subroutine init_histo_w()
  integer   :: nhisto
  character :: c1*(1),c2*(1)
  integer   :: i,j,ntot

  nhisto = 1

  !#1
  Histo(nhisto)%obs_code="rate"
  Histo(nhisto)%LowLimit=zero
  Histo(nhisto)%UpLimit=one
  Histo(nhisto)%Step=one

  !#2
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="ptl"
  Histo(nhisto)%LowLimit=zero
  Histo(nhisto)%UpLimit=60.0_dp
  Histo(nhisto)%Step=0.5_dp

  !#3
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="ptmiss"
  Histo(nhisto)%LowLimit=zero
  Histo(nhisto)%UpLimit=100.0_dp
  Histo(nhisto)%Step=1.0_dp

  !#4
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="eta(l)"
  Histo(nhisto)%LowLimit=-3.0_dp
  Histo(nhisto)%UpLimit=3.0_dp
  Histo(nhisto)%Step=0.1_dp

  !#5
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="pt(W)"
  Histo(nhisto)%LowLimit=zero
  Histo(nhisto)%UpLimit=80.0_dp
  Histo(nhisto)%Step=1.0_dp

  !#6
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="pt(gamma)"
  Histo(nhisto)%LowLimit=zero
  Histo(nhisto)%UpLimit=100.0_dp
  Histo(nhisto)%Step=1.0_dp

  !#7
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="y(gamma)"
  Histo(nhisto)%LowLimit=-5.0_dp
  Histo(nhisto)%UpLimit=5.0_dp
  Histo(nhisto)%Step=0.2_dp

  !#8
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="pt(j)"
  Histo(nhisto)%LowLimit=0_dp
  Histo(nhisto)%UpLimit=100.0_dp
  Histo(nhisto)%Step=1.0_dp

  !#9
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="y(j)"
  Histo(nhisto)%LowLimit=-5.0_dp
  Histo(nhisto)%UpLimit=5.0_dp
  Histo(nhisto)%Step=0.2_dp

  !#10
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="mll"
  Histo(nhisto)%LowLimit=50.0_dp
  Histo(nhisto)%UpLimit=115.0_dp
  Histo(nhisto)%Step=1.0_dp

  !#11
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="mT(W)"
  Histo(nhisto)%LowLimit=zero
  Histo(nhisto)%UpLimit=100.0_dp
  Histo(nhisto)%Step=1.0_dp

  !#12
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="ptl moment"
  Histo(nhisto)%LowLimit=zero
  Histo(nhisto)%UpLimit=one
  Histo(nhisto)%Step=one

  !#13
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="ptl^2 moment"
  Histo(nhisto)%LowLimit=zero
  Histo(nhisto)%UpLimit=one
  Histo(nhisto)%Step=one

  !#14
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="ptl^3 moment"
  Histo(nhisto)%LowLimit=zero
  Histo(nhisto)%UpLimit=one
  Histo(nhisto)%Step=one

  !#15
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="<ptl>/dptl"
  Histo(nhisto)%LowLimit=zero
  Histo(nhisto)%UpLimit=100.0_dp
  Histo(nhisto)%Step=0.5_dp

  !#16
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="<ptl^2>/dptl"
  Histo(nhisto)%LowLimit=zero
  Histo(nhisto)%UpLimit=100.0_dp
  Histo(nhisto)%Step=0.5_dp

  !-- cumulative below

  !#17
  nhisto = nhisto + 1
  Histo(nhisto)%obs_code="c_pT(W)"
  Histo(nhisto)%LowLimit=zero
  Histo(nhisto)%UpLimit=80.0_dp
  Histo(nhisto)%Step=1.0_dp

  if (nhisto.ne.nobs) then
     print *, 'uninitialized histograms'
     print *, nobs, nhisto
     stop
  endif

  Histo(:)%Nbins = (Histo(:)%UpLimit - Histo(:)%LowLimit)/Histo(:)%Step

  return

end subroutine init_histo_w

subroutine cut_histo_w(xa,xb,nprod,pprod,makecut,binarray,weightarray,mu,partid)
! fiducial = 0: no cuts
! fiducial = 1: cuts for paper 2009.xxxxx
! fiducial = 2: cuts per MG
! fiducial = 3: cuts per MCFM
! fiducial = 4/5/6: cuts for Wmass determination, see ATLAS, EPJC (2018) 78:110
  integer, intent(in)            :: nprod
  real(dp), intent(in)           :: xa,xb,pprod(4,nprod)
  integer, intent(in), optional  :: partid(:)
  real(dp), intent(out)          :: weightarray(nobs),mu
  integer, intent(out)           :: binarray(nobs)
  logical, intent(out)           :: makecut
  integer, parameter :: fiducial = 0
  real(dp), parameter :: delta0=0.1_dp
  real(dp), parameter :: ptjetmin = 0.0_dp
  real(dp), parameter :: ptphotonmin=0.0_dp
  real(dp), parameter :: yjetmax = 99.0_dp
  real(dp), parameter :: etaphotonmax = 99.0_dp
  integer, parameter :: njetmin = 0
  real(dp) :: ptlmin,ptmissmin,etalmax,mtwmin
  ! end cuts
  real(dp)  :: pprod_lab(4,nprod),pinp(4,nprod-2),pjetp(4,3)
  real(dp)  :: pdec_lab(4,3),pind(4,3),pjetd(4,3)
  integer   :: njetp,njetd
  integer   :: charges(3),ch_jet(3)
  real(dp)  :: Obs(nobs)
  !-- process-dependent part below
  real(dp)  :: pll(4),lpos(4),lneg(4),lnu(4),lel(4),l1(4),l2(4),lgamma(4), lpartons(1:4,1:2),ll(1:4,1:2)
  real(dp)  :: mll,yll,ptll,mtW
  real(dp)  :: ptpos,ptneg,ptjem3,ptl1,ptl2
  real(dp)  :: ypos,yneg,yjem3,yl1,yl2
  real(dp)  :: costh,dphi,Rgammaj,ptjet,yjet,etaphoton,ptphoton,ptl,etal,ptlnu,ptmiss,pt4,pt5
  integer   :: i,nzero, npartons,hardjet,softjet,ngamma
  logical   :: dodecay, frix_isol
  logical   :: flagrecomb
  logical, save   :: first = .true.

  makecut = .true.
  pprod_lab(:,:) = zero
  pdec_lab(:,:)  = zero
  binarray(:)    = 0
  weightarray(:) = one

  if (fiducial .eq. 1) then
     ptlmin = 15.0_dp
     ptmissmin = 15.0_dp
     etalmax = 2.4_dp
  elseif (fiducial .eq. 4) then
     ptlmin = 30.0_dp
     ptmissmin = 30.0_dp
     etalmax = 2.4_dp
     mtwmin = 60.0_dp
  elseif (fiducial .eq. 5) then
     ptlmin = 20.0_dp
     ptmissmin = 20.0_dp
     etalmax = 2.4_dp
     mtwmin = zero
  elseif (fiducial .eq. 6) then
     ptlmin = 25.44_dp
     ptmissmin = 25.44_dp
     etalmax = 2.4_dp
     mtwmin = zero
  endif
  call boostz(xa,xb,nprod,pprod,pprod_lab)

  pinp(1:4,1:nprod - 4) = pprod_lab(1:4,3:nprod-2)
  lel = zero
  lnu = zero
  lgamma = zero
  lpartons = zero
  npartons = 0
  hardjet = 0
  softjet = 0

  ptjet =zero
  yjet = 99.0_dp


  if (first) then
     print *, "Doing cuts with option fiducial = ", fiducial
     print *, "pt(l) >",ptlmin
     print *, "pt(miss) >",ptmissmin
     print *, "|eta(l) | < ", etalmax
     print *, "mT(W) >",mtwmin
     print *, "pt(jet) >", ptjetmin
     print *, "|y(jet)| >", yjetmax
     print *, "pt(photon) >", ptphotonmin
     print *, "|eta(photon)| >", etaphotonmax
     print *, "Number of resolved jets", njetmin
     if (fiducial .eq. 3) print*, "Frixione isolation with delta=",delta0
     first = .false.
  endif





!  if (.not. present(partid)) then
!     print *, "Error: partid not passed to histogram routine"
!     stop
!  endif
!  call identify_particles(pprod_lab(1:4,3:nprod),nprod-2,partid,lel,lnu,lgamma,lpartons,ngamma,npartons)

! since this is only QCD corrections, there is never ambiguity about the final state particles
  lgamma=zero
  ngamma=0
  npartons=nprod-4
  lel(1:4)=pprod_lab(1:4,nprod)
  lnu(1:4)=pprod_lab(1:4,nprod-1)
  if (npartons .eq. 0) then
     lpartons = zero
  elseif (npartons .eq. 1) then
     lpartons(1:4,1) = pprod_lab(1:4,nprod-2)
  elseif (npartons .eq. 2) then
     lpartons(1:4,1) = pprod_lab(1:4,nprod-2)
     lpartons(1:4,2) = pprod_lab(1:4,nprod-3)
  endif


!  print *, "ngamma",ngamma

  if (fiducial .eq. 2) then    ! Treat photon as a jet
     ll = zero
     if (npartons .gt. 0) ll(1:4,1:npartons) = lpartons(1:4,1:npartons)
     if (ngamma .gt. 0) ll(1:4,npartons+ngamma) = lgamma(1:4)
     call ktjet_w(ll,npartons+ngamma,pjetp,njetp)
     do i=1,njetp
        ptjet=get_pt(pjetp(1:4,i))
        yjet=get_y(pjetp(1:4,i))
        if (ptjet .gt. ptjetmin .and. abs(yjet).lt. yjetmax) then
           hardjet=hardjet+1
        endif
     enddo
     if (hardjet .lt. njetmin) return

  elseif (fiducial .eq. 3) then      ! Frixione clustering (a la MCFM) with jet

     if (ngamma .eq. 1) then
        call frixione_isolation_mcfm(lgamma,lpartons,npartons,delta0,frix_isol)
     else
        frix_isol = .true.
        return
     endif
     if (frix_isol .eqv. .false.) then
!        print *, "failed frix isol"
        return
     endif

     call ktjet_w(lpartons,npartons,pjetp,njetp)
     do i=1,njetp
        Rgammaj=get_r(pjetp(1:4,i),lgamma)
!        if (ngamma .eq. 0 .or. Rgammaj .gt. delta0) then
        if (Rgammaj .gt. delta0) then
           ptjet=get_pt(pjetp(1:4,i))
           yjet=get_y(pjetp(1:4,i))
           if (ptjet .gt. ptjetmin .and. abs(yjet).lt. yjetmax) then
              hardjet=hardjet+1
           endif
        else
           softjet = softjet+1
        endif
     enddo
     if (hardjet .lt. njetmin) return
     ptphoton = get_pt(lgamma)
     etaphoton = get_eta(lgamma)
     if (ptphoton .lt. ptphotonmin .or. abs(etaphoton) .gt. etaphotonmax) return

  else
     call ktjet_w(lpartons,npartons,pjetp,njetp)
     do i=1,njetp
        ptjet=get_pt(pjetp(1:4,i))
        yjet=get_y(pjetp(1:4,i))
        if (ptjet .gt. ptjetmin .and. abs(yjet).lt. yjetmax) then
           hardjet=hardjet+1
        endif
     enddo
  endif

  ptphoton = get_pt(lgamma)
  etaphoton = get_eta(lgamma)

  mu = sqrt(musq_in)

  pll = lnu + lel

  ptl = get_pt(lel)
  ptmiss=get_pt(lnu)
  etal = get_eta(lel)
  ptlnu=get_pt(lel+lnu)
  dphi = get_dphi(lel,lnu)
  mll = sqrt(scr(pll,pll))
  mtw = sqrt(two*ptl*ptmiss*(one-cos(dphi)))

  !  cuts
  if (fiducial .eq. 1) then
     if (ptl .lt. ptlmin) return
     if (ptmiss .lt. ptmissmin) return
     if (abs(etal) .gt. etalmax) return
  endif
  if (fiducial .eq. 4 .or. fiducial .eq. 5 .or. fiducial .eq. 6) then
     if (ptl .lt. ptlmin) return
     if (ptmiss .lt. ptmissmin) return
     if (abs(etal) .gt. etalmax) return
     if (mtw .lt. mtwmin) return
  endif

  !-----------------------------------------------------------------
  !-- fill observables

  !-- standard
  Obs(1)  = 0.5_dp !Actually the "1" moment
  Obs(2)  = ptl
  Obs(3)  = ptmiss
  Obs(4)  = etal
  Obs(5)  = ptlnu
  Obs(6)  = ptphoton
  Obs(7)  = etaphoton
  Obs(8)  = ptjet
  Obs(9)  = yjet
  Obs(10) = mll
  Obs(11) = mtW
  !-- moments
  Obs(12)  = 0.5_dp
  Obs(13)  = 0.5_dp
  Obs(14)  = 0.5_dp
  Obs(15)  = ptl
  Obs(16)  = ptl
  !-- cumulative
  Obs(17) = ptll ! new since phton tag rv check

  call find_bin(Obs,binarray,weightarray)

  !Weight for the moments
!  weightarray(12) = weightarray(12)*(ptl/40.0_dp)
!  weightarray(13) = weightarray(13)*(ptl/40.0_dp)**2
!  weightarray(14) = weightarray(14)*(ptl/40.0_dp)**3

  weightarray(12) = weightarray(12)*(ptl)
  weightarray(13) = weightarray(13)*(ptl)**2
  weightarray(14) = weightarray(14)*(ptl)**3
  weightarray(15) = weightarray(15)*(ptl)
  weightarray(16) = weightarray(16)*(ptl)**2

  makecut = .false.

  return

end subroutine cut_histo_w


!subroutine identify_particles(pin,nin,partid,plep,pneut,pgamma,ppartons,ngamma,npartons)
!  real(dp), intent(in) :: pin(1:4,1:nin)
!  integer, intent(in)  :: nin
!  integer, intent(in)  :: partid(:)
!  real(dp), intent(out):: plep(4), pneut(4),pgamma(4), ppartons(:,:)
!  integer, intent(out) :: ngamma,npartons
!  integer              :: i,nlep,nneut
!
!  !! sanity checks
!  if (size(pin,dim=2) .ne. nin .or. size(partid,dim=1) .ne. nin) then
!     print *, "Error in particle identification"
!     print *, "Number of input particles: ", nin
!     print *, "Size of momentum array: ", size(pin,dim=2)
!     print *, "Size of particle ID array", size(partid,dim=1)
!     stop
!  endif
!
!
!  nlep     = 0
!  nneut    = 0
!  ngamma   = 0
!  npartons = 0
!
!  do i=1,nin
!     if (abs(partid(i)) .eq. elec) then
!        plep(:) = pin(:,i)
!        nlep = nlep+1
!     elseif (abs(partid(i)) .eq. neu) then
!        pneut(:) = pin(:,i)
!        nneut = nneut+1
!     elseif (partid(i) .eq. pho) then
!        pgamma(:) = pin(:,i)
!        ngamma = ngamma+1
!     elseif (partid(i) .eq. glu .or. partid(i) .eq. qua) then
!        npartons = npartons+1
!        ppartons(:,npartons) = pin(:,i)
!     else
!        print *, "Error in particle identification"
!        print *, "Unidentified particle with ID ",partid(i)
!        stop
!     endif
!  enddo
!
!  !! more sanity checks
!  if (nlep .ne. 1 .or. nneut .ne. 1) then
!     print *, "Error in particle identification"
!     print *, "nlep: ", nlep
!     print *, "nneut: ", nneut
!     stop
!  endif
!
!  if (ngamma .gt. 1) then
!     print *, "Error in particle identification"
!     print *, "ngamma: ", ngamma
!     stop
!  endif
!
!  if (ngamma + npartons + nlep + nneut .ne. nin) then
!     print *, "Error in particle identification"
!     print *, "nlep: ", nlep
!     print *, "nneut: ", nneut
!     print *, "ngamma: ", ngamma
!     print *, "npartons: ", npartons
!     print *, "nparticles", nin
!     stop
!  endif
!
!end subroutine identify_particles

subroutine ktjet_w(pin,  nin, pjet, njet)
  integer, intent(in) :: nin
  real(dp), intent(in) :: pin(4,nin)
  integer, intent(out) :: njet
  real(dp), intent(out) :: pjet(4,3)
  real(dp) :: pti,dbi,ptj,dbj
  real(dp) :: yi,yj,dy,d_cosphi,dphi,r12sq,dij
  real(dp) :: dmin,rsq
  real(dp) :: prehard(4,nin),prejet(4,nin),etai
  integer  :: prehard_ch(nin),prejet_ch(nin)
  integer :: nh,nj,ijmin(2)
  integer :: i,j,le
  real(dp), parameter :: ktclass = -one   ! antikt
!  real(dp), parameter :: ptjetcut = 20.0_dp
!  real(dp), parameter :: etajetcut = 99_dp
  real(dp), parameter :: Rjet = 0.4_dp

  rsq = Rjet**2

  le = 0

  njet = 0
  pjet = zero

  do i=1,nin
     if (pin(1,i).ne.zero.and.(abs(pin(4,i)).ne.abs(pin(1,i)))) then
        le = le+1
        prehard(:,le)  = pin(:,i)
     endif
  enddo

  nj = 0
  nh = le
  prejet = zero

  !   recombine

  do while (nh.gt.0) ! recombination loop

     pti = sqrt(prehard(2,nh)**2+prehard(3,nh)**2)
     dbi = pti**(2*ktclass)
     dmin = dbi
     ijmin = (/nh,0/)

     do i = 1, nh-1 ! compute db, dij

        pti = sqrt(prehard(2,i)**2+prehard(3,i)**2)
        dbi = pti**(2*ktclass)

        if (dbi.lt.dmin) then
           dmin = dbi
           ijmin = (/i,0/)
        endif

        do j = i+1,nh
           ptj = sqrt(prehard(2,j)**2+prehard(3,j)**2)
           dbj = ptj**(2*ktclass)

           yi = half * log((prehard(1,i)+prehard(4,i))/(prehard(1,i)-prehard(4,i)))
           yj = half * log((prehard(1,j)+prehard(4,j))/(prehard(1,j)-prehard(4,j)))

           dy = yi-yj

           ! cos(a-b) = cos(a)cos(b)+sin(a)sin(b)
           d_cosphi = (prehard(2,i)*prehard(2,j)+prehard(3,i)*prehard(3,j))/pti/ptj
           if (d_cosphi.gt.one) then
              dphi = zero
           elseif (d_cosphi.lt.-one) then
              dphi = pi
           else
              dphi = acos(d_cosphi) ! 0<dphi<pi
           endif

           r12sq = dy**2 + dphi**2

           dij = min(dbi,dbj) * r12sq/rsq

           if (dij.lt.dmin) then
              dmin = dij
              ijmin = (/min(i,j),max(i,j)/)
           endif

        enddo

     enddo ! compute db, dij

     if (ijmin(2).eq.0) then ! no recombination, 1 new jet candidate
        nj = nj+1
        prejet(:,nj)  = prehard(:,ijmin(1))

        do i = ijmin(1),nh-1
           prehard(:,i)  = prehard(:,i+1)
        enddo

     else ! recombine, no new jet candidates

        prehard(:,ijmin(1))  = prehard(:,ijmin(1))+prehard(:,ijmin(2))

        do i = ijmin(2),nh-1
           prehard(:,i)  = prehard(:,i+1)
        enddo

     endif

     nh = nh - 1

  enddo ! recombination loop

  do i = 1, nj
     pti  = sqrt(prejet(2,i)**2+prejet(3,i)**2)
     etai = get_eta(prejet(:,i))
!     if (pti.gt.ptjetcut .and. abs(etai) .lt. etajetcut) then
        njet = njet + 1
        pjet(:,njet) = prejet(:,i)
!     endif
  enddo

  return

end subroutine ktjet_w


subroutine frixione_isolation_mcfm(p_gamma,ppartons,npartons,delta0,passed)
  ! Frixione photon isolation a la MCFM, for direct comparison
  ! see MCFM/src/User/frix.f -- adapted from these routines
  real(dp), intent(in)      :: p_gamma(4), ppartons(:,:),delta0
  integer, intent(in)       :: npartons
  logical, intent(out)      :: passed
  integer                   :: i
  real(dp)                  :: ETgamma,xi,Rigamma,Ejet
  real(dp),parameter        :: eps_gamma = 1.0_dp
!  real(dp),parameter        :: delta0 = 0.4_dp
  integer, parameter        :: n_pow  = 1


  ETgamma = ret_ET(p_gamma)
  passed = .true.

  do i =1,npartons
     Rigamma=get_r(p_gamma,ppartons(1:4,i))
     if (Rigamma .lt. delta0) then

        Ejet = ret_ET(ppartons(:,i))
        xi = ETgamma*eps_gamma*( (one-cos(Rigamma))/(one-cos(delta0)) )**n_pow

        if (Ejet  .lt. xi) then
           passed =.true.
        else
           passed =.false.
           return
        endif
     endif
  end do

end subroutine frixione_isolation_mcfm


function ret_ET(p)
  real(dp) :: p(4)
  real(dp) :: ret_ET
  real(dp) :: ptsq

  ptsq = p(2)**2 + p(3)**2
  ret_ET = p(1)*sqrt(ptsq)/(sqrt(ptsq+p(4)**2))

end function ret_ET


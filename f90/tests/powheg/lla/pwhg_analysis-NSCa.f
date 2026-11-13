c     Fixed-order analysis for NSC-lnuaEW (pp -> l l' gamma at NLO QCD),
c     twin of the mixed_ll analyses W_lla_nloqcd (W) and lla_nloqcd (Z, no jet cut).
c     Needs detailedNLO 1 (flavour of the radiated parton).
c     Photon: pT > 50, |y| < 2.5, Delta R(gamma, charged lepton) > 0.1,
c             Delta R(gamma, final quark) > 0.4 (no condition for a gluon).
c     W: bare lepton, mT in (200,5000), pT(l) > 65, pT(nu) > 85,
c        pT(l) pT(nu) >= 35^2, |y(l)|, |y(nu)| < 2.4.
c     Z: m(ll) in (200,13000), pT(l) > 20, pT(l+) pT(l-) >= 35^2, |y(l)| < 2.5.
c     Both: pT(l l') > ptllcut (bare lepton pair; powheg.input key, default 30 GeV).
c     Histogram rate_ptll50: the same with pT(l l') > max(ptllcut,50).

      subroutine init_hist
      implicit none
      call inihists
      call bookupeqbins('total',1d0,0d0,1d0)
      call bookupeqbins('rate',1d0,0d0,1d0)
      call bookupeqbins('rate_ptll50',1d0,0d0,1d0)
      call bookupeqbins('rate_g',1d0,0d0,1d0)
      call bookupeqbins('rate_q',1d0,0d0,1d0)
      call bookupeqbins('pta',25d0,50d0,1000d0)
      call bookupeqbins('mV',20d0,200d0,1000d0)
      end

      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0
      include 'hepevt.h'
      include 'pwhg_math.h'
      include 'LesHouches.h'
      include 'pwhg_weights.h'
      include 'pwhg_rwl.h'
      real * 8, allocatable, save :: dsig(:)
      logical, save :: ini = .true.
      integer ihep,il(2),nl,inu,ia,ip
      real * 8 pt1,pt2,mV,pta,pv(4),ptll
      real * 8, save :: ptllcut
      real * 8, external :: powheginput
      real * 8, external :: nsc_pt,nsc_y,nsc_r,nsc_dphi
      real * 8, parameter :: rphot=0.1d0,riso=0.4d0,ptacut=50d0,yacut=2.5d0

      if(ini) then
         call setupmulti(max(1,rwl_num_weights))
         allocate(dsig(max(10,rwl_num_weights+1)))
         ptllcut = powheginput('#ptllcut')
         if(ptllcut.lt.0d0) ptllcut = 30d0
         ini = .false.
      endif

      dsig = 0d0
      if(rwl_num_weights.eq.0) then
         dsig(1) = dsig0
      else
         dsig(1:rwl_num_weights) = rwl_weights(1:rwl_num_weights)
      endif
      if(sum(abs(dsig)).eq.0d0) return

      call filld('total',0.5d0,dsig)

      nl = 0
      inu = 0
      ia = 0
      ip = 0
      do ihep=3,nhep
         if(isthep(ihep).ne.1) cycle
         if(abs(idhep(ihep)).eq.11) then
            nl = nl + 1
            il(nl) = ihep
         elseif(abs(idhep(ihep)).eq.12) then
            inu = ihep
         elseif(idhep(ihep).eq.22) then
            ia = ihep
         else
            ip = ihep
         endif
      enddo
      if(ia.eq.0) return

c     photon
      pta = nsc_pt(phep(1,ia))
      if(pta.le.ptacut) return
      if(abs(nsc_y(phep(1,ia))).ge.yacut) return
      if(ip.gt.0) then
         if(abs(idhep(ip)).ge.1.and.abs(idhep(ip)).le.5) then
            if(nsc_r(phep(1,ip),phep(1,ia)).lt.riso) return
         endif
      endif

      if(inu.gt.0) then
c     W
         pv = phep(1:4,il(1))+phep(1:4,inu)
         if(nsc_r(phep(1,il(1)),phep(1,ia)).lt.rphot) return
         pt1 = nsc_pt(phep(1,il(1)))
         pt2 = nsc_pt(phep(1,inu))
         mV = sqrt(2d0*pt1*pt2*(1d0-cos(nsc_dphi(phep(1,il(1)),phep(1,inu)))))
         if(mV.le.200d0 .or. mV.ge.5000d0) return
         if(pt1.le.65d0) return
         if(pt2.le.85d0) return
         if(pt1*pt2.lt.35d0**2) return
         if(abs(nsc_y(phep(1,il(1)))).ge.2.4d0) return
         if(abs(nsc_y(phep(1,inu))).ge.2.4d0) return
      else
c     Z
         if(min(nsc_r(phep(1,il(1)),phep(1,ia)),nsc_r(phep(1,il(2)),phep(1,ia)),
     &        nsc_r(phep(1,il(1)),phep(1,il(2)))).lt.rphot) return
         pv = phep(1:4,il(1))+phep(1:4,il(2))
         mV = sqrt(pv(4)**2-pv(1)**2-pv(2)**2-pv(3)**2)
         if(mV.le.200d0 .or. mV.ge.13000d0) return
         pt1 = nsc_pt(phep(1,il(1)))
         pt2 = nsc_pt(phep(1,il(2)))
         if(pt1.le.20d0) return
         if(pt2.le.20d0) return
         if(pt1*pt2.lt.35d0**2) return
         if(abs(nsc_y(phep(1,il(1)))).ge.2.5d0) return
         if(abs(nsc_y(phep(1,il(2)))).ge.2.5d0) return
      endif

      ptll = nsc_pt(pv)
      if(ptll.le.ptllcut) return

      call filld('rate',0.5d0,dsig)
      if(ptll.gt.50d0) call filld('rate_ptll50',0.5d0,dsig)
      call filld('pta',pta,dsig)
      call filld('mV',mV,dsig)
      if(ip.gt.0) then
         if(abs(idhep(ip)).ge.1.and.abs(idhep(ip)).le.5) then
            call filld('rate_q',0.5d0,dsig)
         else
            call filld('rate_g',0.5d0,dsig)
         endif
      endif
      end

      real * 8 function nsc_pt(p)
      implicit none
      real * 8 p(4)
      nsc_pt = sqrt(p(1)**2+p(2)**2)
      end

      real * 8 function nsc_y(p)
      implicit none
      real * 8 p(4)
      nsc_y = 0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      end

      real * 8 function nsc_dphi(p1,p2)
      implicit none
      real * 8 p1(4),p2(4),c
      c = (p1(1)*p2(1)+p1(2)*p2(2))
     &     /sqrt((p1(1)**2+p1(2)**2)*(p2(1)**2+p2(2)**2))
      nsc_dphi = acos(max(-1d0,min(1d0,c)))
      end

      real * 8 function nsc_r(p1,p2)
      implicit none
      real * 8 p1(4),p2(4),nsc_y,nsc_dphi
      external nsc_y,nsc_dphi
      nsc_r = sqrt((nsc_y(p1)-nsc_y(p2))**2+nsc_dphi(p1,p2)**2)
      end

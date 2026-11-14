c     Fixed-order analysis for NSC-lnujQCD (pp -> l l' + jet at NLO EW),
c     twin of the mixed_ll analysis V_llj_nloew.  Needs detailedNLO 1.
c     Photon: recombined with the nearest of {charged leptons, final-state quark}
c             if Delta R < 0.1 (four-momentum sum); never with a gluon or the neutrino.
c     Jet: the final QCD parton (dressed quark or gluon), pT > ptjcut (input key, default 30),
c          no rapidity cut.  A photon is never a jet.
c     W: dressed lepton, mT in (200,5000), pT(l) > 65, pT(nu) > 85,
c        pT(l) pT(nu) >= 35^2, |y(l)|, |y(nu)| < 2.4.
c     Z: m(ll) in (200,13000), pT(l) > 20, pT(l+) pT(l-) >= 35^2, |y(l)| < 2.5.
c     Histograms: rate, rate_ptj50 (jet pT > max(ptjcut,50)), rate_qq / rate_g (no / one
c     incoming gluon), rate_a (events with a photon outside all cones), ptj, mV.

      subroutine init_hist
      implicit none
      call inihists
      call bookupeqbins('total',1d0,0d0,1d0)
      call bookupeqbins('rate',1d0,0d0,1d0)
      call bookupeqbins('rate_ptj50',1d0,0d0,1d0)
      call bookupeqbins('rate_qq',1d0,0d0,1d0)
      call bookupeqbins('rate_g',1d0,0d0,1d0)
      call bookupeqbins('rate_a',1d0,0d0,1d0)
      call bookupeqbins('ptj',38.8d0,30d0,1000d0)
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
      integer ihep,il(2),nl,inu,ia,ip,k,kmin
      real * 8 pt1,pt2,mV,ptj,pv(4),pl(4,2),pj(4),r,rmin
      real * 8, save :: ptjcut
      logical isq,loose
      real * 8, external :: powheginput
      real * 8, external :: nsc_pt,nsc_y,nsc_r,nsc_dphi
      real * 8, parameter :: rphot=0.1d0

      if(ini) then
         call setupmulti(max(1,rwl_num_weights))
         allocate(dsig(max(10,rwl_num_weights+1)))
         ptjcut = powheginput('#ptjcut')
         if(ptjcut.lt.0d0) ptjcut = 30d0
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
      if(ip.eq.0) return
      isq = abs(idhep(ip)).ge.1.and.abs(idhep(ip)).le.5
      do k=1,nl
         pl(:,k) = phep(1:4,il(k))
      enddo
      pj = phep(1:4,ip)

c     photon recombination
      loose = .false.
      if(ia.gt.0) then
         if(phep(4,ia).gt.0d0) then
            rmin = rphot
            kmin = 0
            do k=1,nl
               r = nsc_r(pl(1,k),phep(1,ia))
               if(r.lt.rmin) then
                  rmin = r
                  kmin = k
               endif
            enddo
            if(isq) then
               r = nsc_r(pj,phep(1,ia))
               if(r.lt.rmin) then
                  rmin = r
                  kmin = 3
               endif
            endif
            if(kmin.eq.3) then
               pj = pj + phep(1:4,ia)
            elseif(kmin.gt.0) then
               pl(:,kmin) = pl(:,kmin) + phep(1:4,ia)
            else
               loose = .true.
            endif
         endif
      endif

      ptj = nsc_pt(pj)
      if(ptj.le.ptjcut) return

      if(inu.gt.0) then
c     W
         pt1 = nsc_pt(pl(1,1))
         pt2 = nsc_pt(phep(1,inu))
         mV = sqrt(2d0*pt1*pt2*(1d0-cos(nsc_dphi(pl(1,1),phep(1,inu)))))
         if(mV.le.200d0 .or. mV.ge.5000d0) return
         if(pt1.le.65d0) return
         if(pt2.le.85d0) return
         if(pt1*pt2.lt.35d0**2) return
         if(abs(nsc_y(pl(1,1))).ge.2.4d0) return
         if(abs(nsc_y(phep(1,inu))).ge.2.4d0) return
      else
c     Z
         pv = pl(:,1)+pl(:,2)
         mV = sqrt(pv(4)**2-pv(1)**2-pv(2)**2-pv(3)**2)
         if(mV.le.200d0 .or. mV.ge.13000d0) return
         pt1 = nsc_pt(pl(1,1))
         pt2 = nsc_pt(pl(1,2))
         if(pt1.le.20d0) return
         if(pt2.le.20d0) return
         if(pt1*pt2.lt.35d0**2) return
         if(abs(nsc_y(pl(1,1))).ge.2.5d0) return
         if(abs(nsc_y(pl(1,2))).ge.2.5d0) return
      endif

      call filld('rate',0.5d0,dsig)
      if(ptj.gt.50d0) call filld('rate_ptj50',0.5d0,dsig)
      if(idhep(1).eq.21.or.idhep(1).eq.0.or.
     &   idhep(2).eq.21.or.idhep(2).eq.0) then
         call filld('rate_g',0.5d0,dsig)
      else
         call filld('rate_qq',0.5d0,dsig)
      endif
      if(loose) call filld('rate_a',0.5d0,dsig)
      call filld('ptj',ptj,dsig)
      call filld('mV',mV,dsig)
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

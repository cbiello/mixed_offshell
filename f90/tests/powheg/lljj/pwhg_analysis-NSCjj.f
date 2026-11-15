c     Fixed-order analysis for NSC-lnujjMIX (pp -> l l' + 2 quarks, tree interference of
c     O(alpha_s alpha^3)).  Needs detailedNLO 1.
c     Jets: the two final partons, each pT > ptjcut (input key, default 30), no rapidity cut,
c           Delta R(j,j) > rjjcut (input key, default 0.4).
c     Lepton cuts as in pwhg_analysis-NSCj.f (no photon here, bare leptons):
c     W: mT in (200,5000), pT(l) > 65, pT(nu) > 85, pT(l) pT(nu) >= 35^2, |y(l)|, |y(nu)| < 2.4.
c     Z: m(ll) in (200,13000), pT(l) > 20, pT(l+) pT(l-) >= 35^2, |y(l)| < 2.5.
c     Histograms: total, rate, rate_os / rate_ss (incoming flavour ids of opposite / same sign),
c     ptj1 (harder jet), mjj.

      subroutine init_hist
      implicit none
      call inihists
      call bookupeqbins('total',1d0,0d0,1d0)
      call bookupeqbins('rate',1d0,0d0,1d0)
      call bookupeqbins('rate_os',1d0,0d0,1d0)
      call bookupeqbins('rate_ss',1d0,0d0,1d0)
      call bookupeqbins('ptj1',38.8d0,30d0,1000d0)
      call bookupeqbins('mjj',50d0,0d0,2000d0)
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
      integer ihep,il(2),nl,inu,ip(2),np,k
      real * 8 pt1,pt2,mV,pv(4),pl(4,2),pjj(4),mjj
      real * 8, save :: ptjcut,rjjcut
      real * 8, external :: powheginput
      real * 8, external :: nsc_pt,nsc_y,nsc_r,nsc_dphi

      if(ini) then
         call setupmulti(max(1,rwl_num_weights))
         allocate(dsig(max(10,rwl_num_weights+1)))
         ptjcut = powheginput('#ptjcut')
         if(ptjcut.lt.0d0) ptjcut = 30d0
         rjjcut = powheginput('#rjjcut')
         if(rjjcut.lt.0d0) rjjcut = 0.4d0
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
      np = 0
      do ihep=3,nhep
         if(isthep(ihep).ne.1) cycle
         if(abs(idhep(ihep)).eq.11) then
            nl = nl + 1
            il(nl) = ihep
         elseif(abs(idhep(ihep)).eq.12) then
            inu = ihep
         else
            np = np + 1
            if(np.le.2) ip(np) = ihep
         endif
      enddo
      if(np.ne.2) then
         write(*,*) ' analysis: expected two final partons, found ',np
         call pwhg_exit(-1)
      endif
      do k=1,nl
         pl(:,k) = phep(1:4,il(k))
      enddo

      if(nsc_pt(phep(1,ip(1))).le.ptjcut) return
      if(nsc_pt(phep(1,ip(2))).le.ptjcut) return
      if(nsc_r(phep(1,ip(1)),phep(1,ip(2))).le.rjjcut) return

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
      if(idhep(1)*idhep(2).lt.0) then
         call filld('rate_os',0.5d0,dsig)
      else
         call filld('rate_ss',0.5d0,dsig)
      endif
      pjj = phep(1:4,ip(1))+phep(1:4,ip(2))
      mjj = sqrt(max(0d0,pjj(4)**2-pjj(1)**2-pjj(2)**2-pjj(3)**2))
      call filld('ptj1',max(nsc_pt(phep(1,ip(1))),nsc_pt(phep(1,ip(2)))),dsig)
      call filld('mjj',mjj,dsig)
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

c     Fixed-order analysis for the paper benchmarks (SETUP.md section 4): twin of the mixed_ll
c     analyses W_incl (setup A) and W_ATLASfid (setup B) with -require none/photon/jet.
c     Codes: NSC-lnuQCD, NSC-lnuEW (nscrequire 0), NSC-lnuaEW (1), NSC-lnujQCD (2).  Needs detailedNLO 1.
c     powheg.input keys:
c     nscsetup   1: mT > 200 (default);  2: pT(l) > 65, |eta(l)| < 2.4, pT(nu) > 85, 200 < mT < 5000
c     nscrequire 0: a photon within 0.1 of the charged lepton is added to it (default)
c                1: photon pT > 50, |y| < 2.5, R(photon,l) > 0.1, R(photon,final quark) > 0.4, pT(l nu) > 30
c                2: a photon within 0.1 is added to the nearest of charged lepton and final quark;
c                   jet = final quark or gluon with pT > 30

      subroutine init_hist
      implicit none
      integer, external :: nsc_key
      call inihists
      call bookupeqbins('total',1d0,0d0,1d0)
      call bookupeqbins('rate',1d0,0d0,1d0)
      if(nsc_key('#nscsetup',1).eq.1) then
         call bookupeqbins('rate_window',1d0,0.5d0,4.5d0)
         call bookupeqbins('mT',20d0,200d0,1000d0)
         call bookupeqbins('mT_high',100d0,1000d0,3000d0)
         call bookupeqbins('mlnu',20d0,200d0,1000d0)
         call bookupeqbins('mlnu_high',100d0,1000d0,3000d0)
         call bookupeqbins('ptl',20d0,0d0,1000d0)
         call bookupeqbins('yl',0.2d0,-5d0,5d0)
         call bookupeqbins('ptlnu',10d0,0d0,500d0)
         call bookupeqbins('costhCS',0.05d0,-1d0,1d0)
      else
         call bookupeqbins('mT',10d0,200d0,5000d0)
         call bookupeqbins('etal_mT200-300',0.2d0,0d0,2.4d0)
         call bookupeqbins('etal_mT300-425',0.2d0,0d0,2.4d0)
         call bookupeqbins('etal_mT425-600',0.4d0,0d0,2.4d0)
         call bookupeqbins('etal_mT600-900',0.4d0,0d0,2.4d0)
         call bookupeqbins('etal_mT900-2000',0.6d0,0d0,2.4d0)
         call bookupeqbins('ptl',10d0,0d0,2500d0)
         call bookupeqbins('ptnu',10d0,0d0,2500d0)
         call bookupeqbins('etal',0.1d0,-2.4d0,2.4d0)
      endif
      if(nsc_key('#nscrequire',0).eq.1) call bookupeqbins('pta',25d0,50d0,1000d0)
      if(nsc_key('#nscrequire',0).eq.2) call bookupeqbins('ptj',25d0,30d0,1030d0)
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
      integer, save :: setup,require
      integer ihep,il,inu,ia,ip,kmin,i
      real * 8 pl(4),pnu(4),pa(4),pj(4),pv(4),r,rmin,ptobj,ptl,ptnu,mT,etal,mlnu
      logical isq
      integer, external :: nsc_key
      real * 8, external :: nsc_pt,nsc_y,nsc_eta,nsc_r,nsc_dphi,nsc_costh
      real * 8, parameter :: rphot=0.1d0,riso=0.4d0,ptacut=50d0,yacut=2.5d0,ptlnucut=30d0,ptjcut=30d0
      real * 8, parameter :: mTbin(6) = (/200d0,300d0,425d0,600d0,900d0,2000d0/)
      character * 15, parameter :: etaname(5) = (/'etal_mT200-300 ','etal_mT300-425 ',
     &     'etal_mT425-600 ','etal_mT600-900 ','etal_mT900-2000'/)

      if(ini) then
         call setupmulti(max(1,rwl_num_weights))
         allocate(dsig(max(10,rwl_num_weights+1)))
         setup   = nsc_key('#nscsetup',1)
         require = nsc_key('#nscrequire',0)
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

      il = 0
      inu = 0
      ia = 0
      ip = 0
      do ihep=3,nhep
         if(isthep(ihep).ne.1 .or. phep(4,ihep).le.0d0) cycle
         if(abs(idhep(ihep)).eq.11) then
            il = ihep
         elseif(abs(idhep(ihep)).eq.12) then
            inu = ihep
         elseif(idhep(ihep).eq.22) then
            ia = ihep
         else
            ip = ihep
         endif
      enddo
      pl  = phep(1:4,il)
      pnu = phep(1:4,inu)
      isq = .false.
      if(ip.gt.0) then
         pj = phep(1:4,ip)
         isq = abs(idhep(ip)).ge.1 .and. abs(idhep(ip)).le.5
      endif

c     nearest charged particle within rphot of the photon: 1 lepton, 2 final quark (require 2 only)
      kmin = 0
      if(ia.gt.0) then
         pa = phep(1:4,ia)
         rmin = rphot
         r = nsc_r(pl,pa)
         if(r.lt.rmin) then
            rmin = r
            kmin = 1
         endif
         if(require.eq.2 .and. isq) then
            if(nsc_r(pj,pa).lt.rmin) kmin = 2
         endif
      endif

      ptobj = 0d0
      if(require.eq.1) then
         if(ia.eq.0 .or. kmin.gt.0) return
         ptobj = nsc_pt(pa)
         if(ptobj.le.ptacut .or. abs(nsc_y(pa)).ge.yacut) return
         if(isq) then
            if(nsc_r(pj,pa).lt.riso) return
         endif
         if(nsc_pt(pl+pnu).le.ptlnucut) return
      else
         if(kmin.eq.1) pl = pl + pa
         if(kmin.eq.2) pj = pj + pa
         if(require.eq.2) then
            if(ip.eq.0) return
            ptobj = nsc_pt(pj)
            if(ptobj.le.ptjcut) return
         endif
      endif

      ptl  = nsc_pt(pl)
      ptnu = nsc_pt(pnu)
      mT   = sqrt(2d0*ptl*ptnu*(1d0-cos(nsc_dphi(pl,pnu))))
      etal = nsc_eta(pl)
      if(mT.le.200d0) return
      if(setup.eq.2) then
         if(ptl.le.65d0 .or. abs(etal).ge.2.4d0 .or. ptnu.le.85d0 .or. mT.ge.5000d0) return
      endif

      call filld('rate',0.5d0,dsig)
      call filld('mT',mT,dsig)
      call filld('ptl',ptl,dsig)
      if(setup.eq.1) then
         i = 1
         if(mT.gt.300d0)  i = 2
         if(mT.gt.500d0)  i = 3
         if(mT.gt.1500d0) i = 4
         pv = pl + pnu
         mlnu = sqrt(pv(4)**2-pv(1)**2-pv(2)**2-pv(3)**2)
         call filld('rate_window',dble(i),dsig)
         call filld('mT_high',mT,dsig)
         call filld('mlnu',mlnu,dsig)
         call filld('mlnu_high',mlnu,dsig)
         call filld('yl',nsc_y(pl),dsig)
         call filld('ptlnu',nsc_pt(pv),dsig)
         call filld('costhCS',nsc_costh(pl,pnu),dsig)
      else
         do i=1,5
            if(mT.gt.mTbin(i) .and. mT.le.mTbin(i+1)) call filld(trim(etaname(i)),abs(etal),dsig)
         enddo
         call filld('ptnu',ptnu,dsig)
         call filld('etal',etal,dsig)
      endif
      if(require.eq.1) call filld('pta',ptobj,dsig)
      if(require.eq.2) call filld('ptj',ptobj,dsig)
      end

      integer function nsc_key(key,default)
      implicit none
      character *(*) key
      integer default
      real * 8, external :: powheginput
      nsc_key = default
      if(powheginput(key).ge.0d0) nsc_key = nint(powheginput(key))
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

      real * 8 function nsc_eta(p)
      implicit none
      real * 8 p(4),pmod
      pmod = sqrt(p(1)**2+p(2)**2+p(3)**2)
      nsc_eta = 0.5d0*log((pmod+p(3))/(pmod-p(3)))
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
      real * 8 p1(4),p2(4)
      real * 8, external :: nsc_y,nsc_dphi
      nsc_r = sqrt((nsc_y(p1)-nsc_y(p2))**2+nsc_dphi(p1,p2)**2)
      end

c     Collins-Soper angle of p1 in the rest frame of p1+p2 (mixed_ll get_costh_star)
      real * 8 function nsc_costh(p1,p2)
      implicit none
      real * 8 p1(4),p2(4),pv(4),m2,pt2
      pv = p1 + p2
      m2 = pv(4)**2-pv(1)**2-pv(2)**2-pv(3)**2
      pt2 = pv(1)**2+pv(2)**2
      nsc_costh = ((p1(4)+p1(3))*(p2(4)-p2(3))-(p1(4)-p1(3))*(p2(4)+p2(3)))
     &     /sqrt(m2*(m2+pt2))*sign(1d0,pv(3))
      end

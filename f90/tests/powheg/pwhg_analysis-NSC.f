c     Fixed-order analysis for NSC-lnuQCD / NSC-lnuEW that mirrors
c     f90/analysis/mod_cut_histo_W_ATLAS.f90 of mixed_ll:
c     the charged lepton is dressed with a photon within R_phot (0.1),
c     then mT in (200,5000), pT(l) > 65, pT(nu) > 85,
c     pT(l) pT(nu) >= 35^2, |y(l)| < 2.4, |y(nu)| < 2.4.

      subroutine init_hist
      implicit none
      call inihists
      call bookupeqbins('total',1d0,0d0,1d0)
      call bookupeqbins('rate',1d0,0d0,1d0)
      call bookupeqbins('massT',10d0,200d0,3000d0)
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
      real * 8, save :: rphot
      integer ihep,il,inu,iphot
      real * 8 pl(4),pnu(4),ptl,ptnu,yl,ynu,dphi,massT
      real * 8, external :: powheginput,nsc_pt,nsc_y,nsc_dphi

      if(ini) then
         call setupmulti(max(1,rwl_num_weights))
         allocate(dsig(max(10,rwl_num_weights+1)))
         rphot = powheginput('#R_phot')
         if(rphot.lt.0d0) rphot = 0.1d0
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
      iphot = 0
      do ihep=3,nhep
         if(abs(idhep(ihep)).eq.11) il = ihep
         if(abs(idhep(ihep)).eq.12) inu = ihep
         if(idhep(ihep).eq.22) iphot = ihep
      enddo

      pl  = phep(1:4,il)
      pnu = phep(1:4,inu)
c     dress the charged lepton
      if(iphot.gt.0) then
         if(sqrt((nsc_y(pl)-nsc_y(phep(1,iphot)))**2
     &        + nsc_dphi(pl,phep(1,iphot))**2).lt.rphot)
     &        pl = pl + phep(1:4,iphot)
      endif

      ptl  = nsc_pt(pl)
      ptnu = nsc_pt(pnu)
      dphi = nsc_dphi(pl,pnu)
      massT = sqrt(2d0*ptl*ptnu*(1d0-cos(dphi)))
      yl  = nsc_y(pl)
      ynu = nsc_y(pnu)

      if(massT.le.200d0 .or. massT.ge.5000d0) return
      if(ptl.le.65d0) return
      if(ptnu.le.85d0) return
      if(ptl*ptnu.lt.35d0**2) return
      if(abs(yl).ge.2.4d0) return
      if(abs(ynu).ge.2.4d0) return

      call filld('rate',0.5d0,dsig)
      call filld('massT',massT,dsig)
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

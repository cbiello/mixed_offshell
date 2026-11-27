module mod_xsects_nnlo_rr_wp_gg
  use mod_types; use common_def; use mod_consts_dp
#if (_DebugFeatures == 1)
  use mod_consts_debug
#endif
  use mod_kinematics_gen; use mod_kinematicsNNLO_EunordAunord
  use mod_process
  use mod_auxfunctions; use mod_procfunctions
  use mod_histograms
  use mod_lumi
  use mod_auxlumi_wp
  use mod_amplitudes_tree_wp
  use mod_splitting
  implicit none
  integer, parameter :: namp_lo   = 2
  integer, parameter :: namp_nlo  = 6
  integer, parameter :: namp_nnlo = 1

  private

  public ::  xsect_nnlo_rr_wp_gg

contains

  function xsect_nnlo_rr_wp_gg(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_rr_wp_gg,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,&
         C41Lim,C51Lim,C42Lim,C52Lim,&
         C4152Lim,C4251Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,mudum,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b
    complex(dp) :: respdf_sping(-1:1,-1:1),respdf_sping2(-1:1,-1:1,-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping2(namp_lo,-1:1,-1:1,-1:1,-1:1),res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    complex(dp) :: hell_mat4(-1:1,-1:1),hell_mat5(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z4,z5,s14,s15,s24,s25,s45,s4_15,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_gg(7),kin(1:7)
    logical  :: makecut

    callcount = callcount + 1

    xsect_nnlo_rr_wp_gg = 0
    ff(1) = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    C41Lim=>MyKinSet%AllKinConfig(2)
    C51Lim=>MyKinSet%AllKinConfig(3)
    C42Lim=>MyKinSet%AllKinConfig(4)
    C52Lim=>MyKinSet%AllKinConfig(5)
    C4152Lim=>MyKinSet%AllKinConfig(6)
    C4251Lim=>MyKinSet%AllKinConfig(7)

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if ( xx(xRHO1)*xx(xRHO2).lt.buffrr .or. &
         xx(xRHO1)*(one-xx(xRHO2)).lt.buffrr .or. &
         (one-xx(xRHO1))*xx(xRHO2).lt.buffrr .or. &
         (one-xx(xRHO1))*(one-xx(xRHO2)).lt.buffrr ) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,&
         C41Lim,C51Lim,C42Lim,C52Lim,&
         C4152Lim,C4251Lim)
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_gg_qqb(HardProc%AmpMom,res_jj)
       respdf = gg_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !C41Lim=>MyKinSet%AllKinConfig(2)

    xi1 = C41Lim%PartFrac(1) 
    xi2 = C41Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C41Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C41Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(2) = zero       
#endif

       kin(2) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       !call restree_j_sping(1,C41Lim%AmpMom,res_sping_j,res_j)
       call restree_j(C41Lim%AmpMom,res_sping_j,res_j)
       respdf = gg_lumi_qxg(res_j,f1,f2)

       s14 = C41Lim%KinInv(1) 
       z1 = C41Lim%KinInv(4)

       kin(2) = -(-two/s14)*respdf*Pqg_flip(z1)*C41Lim%wgt*jac*pref

!        nperp=C41Lim%LimMom(2:4,1)
!        ncoll=C41Lim%LimMom(2:4,2)

!        hell_mat = comp_hell_mat_pgq_flip(z1,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
!        call zipp(respdf_sping,hell_mat,respdf)

!        kin(2) = -(-two/s14)*respdf*C41Lim%wgt*jac*pref*half !-- identical particles
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(2) = kin(2)
#endif

       call HistoFill(5,binarray,weightarray,kin(2),vegasweight,iternumber)

    endif
 
    !-- 

    !C51Lim=>MyKinSet%AllKinConfig(3)

    xi1 = C51Lim%PartFrac(1) 
    xi2 = C51Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C51Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C51Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(3) = zero
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       !call restree_j_sping(1,C51Lim%AmpMom,res_sping_j,res_j)
       call restree_j(C51Lim%AmpMom,res_sping_j,res_j)
       respdf = gg_lumi_qg(res_j,f1,f2)
       !respdf_sping(:,:) = res_sping_j(3,:,:) * get_lumi_wph_gg(f1,f2)

       s15 = C51Lim%KinInv(1) 
       z1 = C51Lim%KinInv(4)

       kin(3) = -(-two/s15)*respdf*Pqg_flip(z1)*C51Lim%wgt*jac*pref

!        nperp=C51Lim%LimMom(2:4,1)
!        ncoll=C51Lim%LimMom(2:4,2)

!        hell_mat = comp_hell_mat_pgq_flip(z1,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
!        call zipp(respdf_sping,hell_mat,respdf)

!        kin(3) = -(-two/s15)*respdf*C51Lim%wgt*jac*pref*half 

#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(3) = kin(3)
#endif

       call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber)

    endif

    !--

    !C42Lim=>MyKinSet%AllKinConfig(4)

    xi1 = C42Lim%PartFrac(1) 
    xi2 = C42Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C42Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C42Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(4) = zero
#endif

       kin(4) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       !call restree_j_sping(2,C42Lim%AmpMom,res_sping_j,res_j)
       call restree_j(C42Lim%AmpMom,res_sping_j,res_j)
       respdf = gg_lumi_gqx(res_j,f1,f2)

       s24 = C42Lim%KinInv(1) 
       z2 = C42Lim%KinInv(4)

       kin(4) = -(-two/s24)*respdf*Pqg_flip(z2)*C42Lim%wgt*jac*pref

!        nperp=C42Lim%LimMom(2:4,1)
!        ncoll=C42Lim%LimMom(2:4,2)

!        hell_mat = comp_hell_mat_pgq_flip(z2,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
!        call zipp(respdf_sping,hell_mat,respdf)
       
!        kin(4) = -(-two/s24)*respdf*C42Lim%wgt*jac*pref*half !-- the different channel average cancels out

#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(4) = kin(4)
#endif

       call HistoFill(5,binarray,weightarray,kin(4),vegasweight,iternumber)

    endif
 
    !-- 

    !C52Lim=>MyKinSet%AllKinConfig(5)

    xi1 = C52Lim%PartFrac(1) 
    xi2 = C52Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C52Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C52Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(5) = zero
#endif

       kin(5) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       !call restree_j_sping(2,C52Lim%AmpMom,res_sping_j,res_j)
       call restree_j(C52Lim%AmpMom,res_sping_j,res_j)
       respdf = gg_lumi_gq(res_j,f1,f2)

       s25 = C52Lim%KinInv(1) 
       z2 = C52Lim%KinInv(4)

       kin(5) = -(-two/s25)*respdf*Pqg_flip(z2)*C52Lim%wgt*jac*pref

!        nperp=C52Lim%LimMom(2:4,1)
!        ncoll=C52Lim%LimMom(2:4,2)

!        hell_mat = comp_hell_mat_pgq_flip(z2,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
!        call zipp(respdf_sping,hell_mat,respdf)

!        kin(5) = -(-two/s25)*respdf*C52Lim%wgt*jac*pref*half !-- the different channel average cancels out

#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(5) = kin(5)
#endif

       call HistoFill(5,binarray,weightarray,kin(5),vegasweight,iternumber)

    endif


    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!
    
    !C4152Lim=>MyKinSet%AllKinConfig(6)

    xi1 = C4152Lim%PartFrac(1) 
    xi2 = C4152Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,C4152Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C4152Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(6) = zero
#endif

       kin(6) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       !call restree_sping2(C4152Lim%AmpMom,res_sping2)
       call restree(C4152Lim%AmpMom,res)
       respdf = gg_lumi_qxq(res,f1,f2)

       s14 = C4152Lim%KinInv(1) 
       s25 = C4152Lim%KinInv(2) 
       z1 =  C4152Lim%KinInv(4) 
       z2 =  C4152Lim%KinInv(5) 
       !
       
       kin(6) = (two/s14)*(two/s25)*Pqg_flip(z1)*Pqg_flip(z2)*respdf*C4152Lim%wgt*jac*pref

!        nperp=C4152Lim%LimMom(2:4,1)
!        ncoll=C4152Lim%LimMom(2:4,2)
!        hell_mat4 = comp_hell_mat_pgq_flip(z1,nperp,ncoll) ! n1: correlation vector, n2: collinear direction

!        nperp=C4152Lim%LimMom(2:4,3)
!        ncoll=C4152Lim%LimMom(2:4,4)
!        hell_mat5 = comp_hell_mat_pgq_flip(z2,nperp,ncoll) ! n1: correlation vector, n2: collinear direction

!        call zipp2(hell_mat4,respdf_sping2,hell_mat5,respdf)

!        kin(6) = (two/s14)*(two/s25)*respdf*C4152Lim%wgt*jac*pref*half !-- the different channel average cancels out
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(6) = kin(6)
#endif

       call HistoFill(6,binarray,weightarray,kin(6),vegasweight,iternumber)
       
    endif

    !--

    !C4251Lim=>MyKinSet%AllKinConfig(7)

    xi1 = C4251Lim%PartFrac(1) 
    xi2 = C4251Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,C4251Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C4251Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(7) = zero
#endif

       kin(7) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       !call restree_sping2(C4251Lim%AmpMom,res_sping2)
       call restree(C4251Lim%AmpMom,res)
       respdf = gg_lumi_qqx(res,f1,f2)

       s24 = C4251Lim%KinInv(1) 
       s15 = C4251Lim%KinInv(2) 
       z2 =  C4251Lim%KinInv(4) 
       z1 =  C4251Lim%KinInv(5) 
       !

       kin(7) = (two/s24)*(two/s15)*Pqg_flip(z2)*Pqg_flip(z1)*respdf*C4251Lim%wgt*jac*pref

!        nperp=C4251Lim%LimMom(2:4,1)
!        ncoll=C4251Lim%LimMom(2:4,2)
!        hell_mat4 = comp_hell_mat_pgq_flip(z2,nperp,ncoll) ! n1: correlation vector, n2: collinear direction

!        nperp=C4251Lim%LimMom(2:4,3)
!        ncoll=C4251Lim%LimMom(2:4,4)
!        hell_mat5 = comp_hell_mat_pgq_flip(z1,nperp,ncoll) ! n1: correlation vector, n2: collinear direction

!        call zipp2(hell_mat5,respdf_sping2,hell_mat4,respdf)

!        kin(7) = (two/s24)*(two/s15)*respdf*C4251Lim%wgt*jac*pref*half !-- the different channel average cancels out

#if (_DebugFeatures == 1)
       FintNNLO_wp_gg(7) = kin(7)
#endif
       
    endif
    call HistoFill(10,binarray,weightarray,kin(7),vegasweight,iternumber)

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,7
          print *, i, FintNNLO_wp_gg(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_gg = FintNNLO_wp_gg
#endif

    ! do i = 1,7
    !    print *, i,FintNNLO_wp_gg(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause
       
    return

  end function xsect_nnlo_rr_wp_gg

end module mod_xsects_nnlo_rr_wp_gg

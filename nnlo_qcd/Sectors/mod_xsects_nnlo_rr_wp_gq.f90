module mod_xsects_nnlo_rr_wp_gq
  use mod_types; use common_def; use mod_consts_dp
#if (_DebugFeatures == 1)
  use mod_consts_debug
#endif
  use mod_kinematics_gen; use mod_kinematicsNNLO_EunordAord
  use mod_process
  use mod_auxfunctions; use mod_procfunctions
  use mod_histograms
  use mod_lumi
  use mod_auxlumi_wp
  use mod_amplitudes_tree_wp
  use mod_limits_dy_gq
  use mod_splitting
  implicit none
  integer, parameter :: namp_lo   = 2
  integer, parameter :: namp_nlo  = 6
  integer, parameter :: namp_nnlo = 2

  private

  public :: xsect_nnlo_4151a_rr_wp_gq, xsect_nnlo_4151b_rr_wp_gq, xsect_nnlo_4151c_rr_wp_gq, xsect_nnlo_4151d_rr_wp_gq
  public :: xsect_nnlo_4251_rr_wp_gq

  public :: xsect_nnlo_4252a_rr_wp_gq, xsect_nnlo_4252b_rr_wp_gq, xsect_nnlo_4252c_rr_wp_gq, xsect_nnlo_4252d_rr_wp_gq
  public :: xsect_nnlo_4152_rr_wp_gq

contains

  function xsect_nnlo_4151a_rr_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4151a_rr_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,S5Lim,TCLim,TCS5Lim,C5Lim,C5S5Lim,TCC5Lim,TCC5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b,mylumi
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z4,z5,s14,s15,s45,s4_15,kt1kt2sq,reseik
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_gq(8),kin(1:5)
    logical  :: makecut

    callcount = callcount + 1
    xsect_nnlo_4151a_rr_wp_gq = 0
    ff(1) = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    S5Lim=>MyKinSet%AllKinConfig(2)
    TCLim=>MyKinSet%AllKinConfig(3)
    TCS5Lim=>MyKinSet%AllKinConfig(4)
    C5Lim=>MyKinSet%AllKinConfig(5)
    C5S5Lim=>MyKinSet%AllKinConfig(6)
    TCC5Lim=>MyKinSet%AllKinConfig(7)
    TCC5S5Lim=>MyKinSet%AllKinConfig(8)

    xx = zero    
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4151a(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S5Lim,TCLim,TCS5Lim,&
         C5Lim,C5S5Lim,TCC5Lim,TCC5S5Lim)
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    
    if (makecut.or.HardProc%flag) then


#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_gq_qg(HardProc%AmpMom,res_jj)
       respdf = gq_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5
    !-- 3 -> S5C5

    xi1 = S5Lim%PartFrac(1) 
    xi2 = S5Lim%PartFrac(2) 
    call cut_histo_w(xi1,xi2,5,S5Lim%AmpMom,makecut,binarray,weightarray,mu) 

    if (makecut.or.S5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(2) = zero
       FintNNLO_wp_gq(3) = zero
#endif       

       kin(2) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(S5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       reseik = Eik_s5(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3),S5Lim%LimMom(:,4)) 
       
       FintNNLO_wp_gq(2)  = (-one)*respdf*reseik*S5Lim%wgt*jac*pref

       s15 = C5S5Lim%KinInv(2)
       z5 = C5S5Lim%KinInv(5)
       FintNNLO_wp_gq(3) =  -four*Ca/z5/s15 * respdf * C5S5Lim%wgt*jac*pref

       kin(2) = FintNNLO_wp_gq(2) + FintNNLO_wp_gq(3)
       call HistoFill(5,binarray,weightarray,kin(2),vegasweight,iternumber)

    endif

    !! -----------------------------!!
    !!     calling NLO ME  w sp corr!!
    !! -----------------------------!!

    !-- 4 -> C5

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4)  = zero
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C5Lim%AmpMom,res_sping_j,res_j)
       respdf_sping = gq_lumi_qg_sping(res_sping_j,f1,f2)

       s15 = C5Lim%KinInv(2) 
       z2 = C5Lim%KinInv(5)

       nperp=C5Lim%LimMom(2:4,1)
       ncoll=C5Lim%LimMom(2:4,2)

       hell_mat = comp_hell_mat_pgg(z2,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
       call zipp(respdf_sping,hell_mat,respdf)

       kin(3) = (-two/s15)*respdf*C5Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4)  = kin(3)
#endif

      call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber)
    endif

    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!
    
    !-- 5 -> TC
    !-- 6 -> TCC5


    xi1 = TCLim%PartFrac(1) 
    xi2 = TCLim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCLim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCLim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(5:6)  = zero
#endif
       
       kin(4) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCLim%AmpMom,res)
       respdf = gq_lumi_lo(res,f1,f2)

       s14 = TCLim%KinInv(1)
       s15 = TCLim%KinInv(2)
       s45 = TCLim%KinInv(3)
       !
       z4 = TCLim%KinInv(4)
       z5 = TCLim%KinInv(5)
       z2 = one-z4-z5

       FintNNLO_wp_gq(5) = (-one)*respdf*(Pggq_flip(s15,s45,s14,z5,z2,z4))*TCLim%wgt*jac*pref 

       s15 = TCC5Lim%KinInv(1)
       s4_15 = TCC5Lim%KinInv(2)
       z1 = TCC5Lim%KinInv(4)
       z2 = TCC5Lim%KinInv(5)
       kt1kt2sq = TCC5Lim%KinInv(6)

       FintNNLO_wp_gq(6) =  Ptcso_gg_q_flip(s15,s4_15,z1,z2,kt1kt2sq)
       FintNNLO_wp_gq(6) =  FintNNLO_wp_gq(6) *respdf * TCC5Lim%wgt*jac*pref

       
       kin(4) = sum(FintNNLO_wp_gq(5:6))

       call HistoFill(5,binarray,weightarray,kin(4),vegasweight,iternumber)
    endif


   !-- 7 -> TCS5
   !-- 8 -> TCC5S5

    xi1 = TCS5Lim%PartFrac(1) 
    xi2 = TCS5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCS5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCS5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(7:8)  = zero
#endif
       
       kin(5) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCS5Lim%AmpMom,res)
       respdf = gq_lumi_lo(res,f1,f2)
       
       FintNNLO_wp_gq(7) = Pggq_s5(TCS5Lim%KinInv)
       FintNNLO_wp_gq(7) = FintNNLO_wp_gq(7)*respdf*TCS5Lim%wgt*jac*pref

       s14 = TCC5S5Lim%KinInv(1)
       s15 = TCC5S5Lim%KinInv(2)
       z4 = TCC5S5Lim%KinInv(4)
       z5 = TCC5S5Lim%KinInv(5)
       FintNNLO_wp_gq(8) = (-one)*(four*Ca/s15/z5*two/s14*Pqg_i(one-z4)) *respdf * TCC5S5Lim%wgt*jac*pref

       kin(5) = sum(FintNNLO_wp_gq(7:8))

    endif
    call HistoFill(10,binarray,weightarray,kin(5),vegasweight,iternumber)

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1) .or. abs(ff(1)) .gt. 1d200) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,8
          print *, i, FintNNLO_wp_gq(i)
       enddo
       ff(1) = zero
       stop
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_gq = FintNNLO_wp_gq
#endif

    ! do i = 1,8
    !    print *, i,FintNNLO_wp_gq(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return

  end function xsect_nnlo_4151a_rr_wp_gq

  function xsect_nnlo_4151b_rr_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4151b_rr_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,S5Lim,TCLim,TCS5Lim,&
         C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z3,z4,z5,z45,s12,s14,s24,s15,s45,s1_45,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_gq(8),kin(1:5)
    logical :: makecut

    callcount = callcount + 1

    xsect_nnlo_4151b_rr_wp_gq = 0
    ff(1) = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    S5Lim=>MyKinSet%AllKinConfig(2)
    TCLim=>MyKinSet%AllKinConfig(3)
    TCS5Lim=>MyKinSet%AllKinConfig(4)
    C54Lim=>MyKinSet%AllKinConfig(5)
    C54S5Lim=>MyKinSet%AllKinConfig(6)
    TCC54Lim=>MyKinSet%AllKinConfig(7)
    TCC54S5Lim=>MyKinSet%AllKinConfig(8)

    xx = zero
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*(xx(xRHO2)**buffrr_pow).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4151b(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S5Lim,TCLim,TCS5Lim,&
         C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim)
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_gq_qg(HardProc%AmpMom,res_jj)
       respdf = gq_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not


    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5
    !-- 3 -> C54S5

    xi1 = S5Lim%PartFrac(1)
    xi2 = S5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,S5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.S5Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(2:3) = zero       
#endif

       kin(2) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(S5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       reseik = Eik_s5(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3),S5Lim%LimMom(:,4)) 
       FintNNLO_wp_gq(2)  = (-one)*respdf*reseik*S5Lim%wgt*jac*pref

       s45 = C54S5Lim%KinInv(3) 
       z5 =  C54S5Lim%KinInv(5)
       !
       FintNNLO_wp_gq(3) =  four*Cf/z5/s45 * respdf * C54S5Lim%wgt*jac*pref

       kin(2) = sum(FintNNLO_wp_gq(2:3))
       call HistoFill(5,binarray,weightarray,kin(2),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    endif

    !-- 4 -> C5

    xi1 = C54Lim%PartFrac(1)
    xi2 = C54Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C54Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C54Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4) = zero       
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C54Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       s45 = C54Lim%KinInv(3) 
       z5  = C54Lim%KinInv(5) 
       !
       kin(3) = (-two)/s45*Pqq(z5)*respdf*C54Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4) = kin(3)
#endif

       call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    endif


    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!

    !-- 5 -> TC
    !-- 6 -> TCC54

    xi1 = TCLim%PartFrac(1)
    xi2 = TCLim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCLim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCLim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(5:6) = zero       
#endif

       kin(4) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCLim%AmpMom,res)
       respdf = gq_lumi_lo(res,f1,f2)

       s14 = TCLim%KinInv(1)
       s15 = TCLim%KinInv(2)
       s45 = TCLim%KinInv(3)
       !
       z4 = TCLim%KinInv(4)
       z5 = TCLim%KinInv(5)
       z2 = one-z4-z5

       FintNNLO_wp_gq(5) = (-one)*respdf*(Pggq_flip(-s15,s45,-s14,z5,z2,z4))*TCLim%wgt*jac*pref 

       s45 = TCC54Lim%KinInv(1)
       s1_45 = TCC54Lim%KinInv(3)
       z1 = TCC54Lim%KinInv(4)
       z2 = TCC54Lim%KinInv(5)

       FintNNLO_wp_gq(6) =  four/s45/s1_45 * Pqg(z1) * Pqg_flip(one-z2)
       FintNNLO_wp_gq(6) =  FintNNLO_wp_gq(6) *respdf * TCC54Lim%wgt*jac*pref

       kin(4) = sum(FintNNLO_wp_gq(5:6))
       call HistoFill(5,binarray,weightarray,kin(4),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    endif

    !-- 7 -> TCS5
    !-- 8 -> TCC54S5

    xi1 = TCS5Lim%PartFrac(1)
    xi2 = TCS5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCS5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCS5Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(7:8) = zero       
#endif

       kin(5) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCS5Lim%AmpMom,res)
       respdf = gq_lumi_lo(res,f1,f2)

       FintNNLO_wp_gq(7) = Pggq_s5(TCS5Lim%KinInv)
       FintNNLO_wp_gq(7) = FintNNLO_wp_gq(7)*respdf*TCS5Lim%wgt*jac*pref

       s14 = TCC54S5Lim%KinInv(1)
       s15 = TCC54S5Lim%KinInv(2)
       s45 = TCC54S5Lim%KinInv(3)
       z4 = TCC54S5Lim%KinInv(4)
       z5 = TCC54S5Lim%KinInv(5)
       FintNNLO_wp_gq(8) = -(four*Cf/s45/z5*two/s14*Pqg_flip(one-z4)) *respdf * TCC54S5Lim%wgt*jac*pref

       kin(5) = sum(FintNNLO_wp_gq(7:8))

    endif
    call HistoFill(10,binarray,weightarray,kin(5),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,8
          print *, i, FintNNLO_wp_gq(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_gq = FintNNLO_wp_gq
#endif

    ! do i = 1,8
    !    print *, i,FintNNLO_wp_gq(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return

  end function xsect_nnlo_4151b_rr_wp_gq

  function xsect_nnlo_4151c_rr_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4151c_rr_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,S5Lim,TCLim,TCS5Lim,C4Lim,C4S5Lim,TCC4Lim,TCC4S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b,mylumi
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s24,s25,s45,s5_14,s1_45,kt1kt2sq,reseik
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_gq(8),kin(1:5)
    logical  :: makecut

    callcount = callcount + 1

    xsect_nnlo_4151c_rr_wp_gq = 0
    ff(1) = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    S5Lim=>MyKinSet%AllKinConfig(2)
    TCLim=>MyKinSet%AllKinConfig(3)
    TCS5Lim=>MyKinSet%AllKinConfig(4)
    C4Lim=>MyKinSet%AllKinConfig(5)
    C4S5Lim=>MyKinSet%AllKinConfig(6)
    TCC4Lim=>MyKinSet%AllKinConfig(7)
    TCC4S5Lim=>MyKinSet%AllKinConfig(8)

    xx = zero
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4151c(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S5Lim,TCLim,TCS5Lim,&
         C4Lim,C4S5Lim,TCC4Lim,TCC4S5Lim)
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_gq_qg(HardProc%AmpMom,res_jj)
       respdf = gq_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    !-- 2 -> C4

    xi1 = C4Lim%PartFrac(1)
    xi2 = C4Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C4Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C4Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(2) = zero       
#endif

       kin(2) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C4Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_ns(res_j,f1,f2)

       s14 = C4Lim%KinInv(2) 
       z4 = C4Lim%KinInv(5)

       kin(2) = -(two/s14)*respdf*Pqg_flip(one-z4)*C4Lim%wgt*jac*pref

!        call restree_j_sping(2,C4Lim%AmpMom,res_sping_j,res_j)
!        respdf_sping(:,:) = res_sping_j(1,:,:) * get_lumi_wp_qg(f1,f2)

!        s24 = C4Lim%KinInv(2) 
!        z4 = C4Lim%KinInv(5)

!        nperp=C4Lim%LimMom(2:4,1)
!        ncoll=C4Lim%LimMom(2:4,2)

!        hell_mat = comp_hell_mat_pgq_flip(z4,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
!        call zipp(respdf_sping,hell_mat,respdf)

!        kin(2) = (two/s24)*respdf*C4Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(2) = kin(2)
#endif

       call HistoFill(5,binarray,weightarray,kin(2),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    endif

    !-- 3 -> S5

    xi1 = S5Lim%PartFrac(1)
    xi2 = S5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,S5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.S5Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(3) = zero       
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(S5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       reseik = Eik_s5(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3),S5Lim%LimMom(:,4)) 

       kin(3) = (-one)*respdf*reseik*S5Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(3) = kin(3)
#endif

       call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    endif
    
    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!
    
    !-- 4 -> TC
    !-- 5 -> TCC4

    xi1 = TCLim%PartFrac(1)
    xi2 = TCLim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCLim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCLim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4:5) = zero       
#endif

       kin(4) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCLim%AmpMom,res)
       respdf = gq_lumi_lo(res,f1,f2)

       s14 = TCLim%KinInv(1)
       s15 = TCLim%KinInv(2)
       s45 = TCLim%KinInv(3)
       !
       z4 = TCLim%KinInv(4)
       z5 = TCLim%KinInv(5)
       z2 = one-z4-z5

       FintNNLO_wp_gq(4) = (-one)*respdf*(Pggq_flip(s15,s45,s14,z5,z2,z4))*TCLim%wgt*jac*pref 

       s14 = TCC4Lim%KinInv(1) 
       s5_14 = TCC4Lim%KinInv(2) 
       z1  = TCC4Lim%KinInv(4) 
       z2 = TCC4Lim%KinInv(5)
       kt1kt2sq = TCC4Lim%KinInv(6)

       FintNNLO_wp_gq(5) = four/s14/s5_14*Pqg_flip(one-z1)*Pqq(z2) * respdf * TCC4Lim%wgt*jac*pref 

       kin(4) = sum(FintNNLO_wp_gq(4:5))
       call HistoFill(5,binarray,weightarray,kin(4),vegasweight,iternumber) 

    endif
    
    !-- 6 -> TCS5
    !-- 7 -> C4S5
    !-- 8 -> TCC4S5

    xi1 = TCS5Lim%PartFrac(1)
    xi2 = TCS5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCS5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCS5Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(6:8) = zero       
#endif

       kin(5) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCS5Lim%AmpMom,res)
       respdf = gq_lumi_lo(res,f1,f2)

       FintNNLO_wp_gq(6) = Pggq_s5(TCS5Lim%KinInv)
       FintNNLO_wp_gq(6) = FintNNLO_wp_gq(6)*respdf*TCS5Lim%wgt*jac*pref

       s12 = C4S5Lim%EikInv(1)
       s15 = C4S5Lim%EikInv(2)
       s25 = C4S5Lim%EikInv(3)
       s24 = C4S5Lim%KinInv(1)
       z4 =  C4S5Lim%KinInv(4)


!       s12 = C4S5Lim%KinInv(1) 
!       s15 = C4S5Lim%KinInv(2) 
!       s25 = C4S5Lim%KinInv(3) 
!       s24 = C4S5Lim%KinInv(4) 
!       z4 =  C4S5Lim%KinInv(5) 
       !
       FintNNLO_wp_gq(7) =  -8.0_dp*Cf * s12/s15/s25 * one/s24*Pqg_flip(one-z4) *  respdf * C4S5Lim%wgt*jac*pref

       s14 = TCC4S5Lim%KinInv(1) 
       s15 = TCC4S5Lim%KinInv(2) 
       s45 = TCC4S5Lim%KinInv(3) 
       z4 =  TCC4S5Lim%KinInv(4) 
       z5 =  TCC4S5Lim%KinInv(5) 
       !
       FintNNLO_wp_gq(8) =  (one)*Cf*8.0_dp/s15/z5/s14*Pqg_flip(one-z4) *respdf *TCC4S5Lim%wgt*jac*pref

       kin(5) = sum(FintNNLO_wp_gq(6:8))

    endif
    call HistoFill(10,binarray,weightarray,kin(5),vegasweight,iternumber) 

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,8
          print *, i, FintNNLO_wp_gq(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_gq(1) = FintNNLO_wp_gq(1)
    FintNNLO_CheckLim_wp_gq(2) = FintNNLO_wp_gq(3) !-- S5
    FintNNLO_CheckLim_wp_gq(3) = FintNNLO_wp_gq(7) !-- C4S5
    FintNNLO_CheckLim_wp_gq(4) = FintNNLO_wp_gq(2) !-- C4
    FintNNLO_CheckLim_wp_gq(5) = FintNNLO_wp_gq(4) !-- TC
    FintNNLO_CheckLim_wp_gq(6) = FintNNLO_wp_gq(5) !-- TCC4
    FintNNLO_CheckLim_wp_gq(7) = FintNNLO_wp_gq(6) !-- TCS5
    FintNNLO_CheckLim_wp_gq(8) = FintNNLO_wp_gq(8) !-- TCC4S5
#endif

    ! do i = 1,8
    !    print *, i,FintNNLO_wp_gq(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return

  end function xsect_nnlo_4151c_rr_wp_gq

  function xsect_nnlo_4151d_rr_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4151d_rr_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,S5Lim,TCLim,TCS5Lim,&
         C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z3,z4,z5,z45,s12,s14,s24,s15,s45,s1_45,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_gq(8),kin(1:6)
    logical :: makecut

    callcount = callcount + 1

    xsect_nnlo_4151d_rr_wp_gq = 0
    ff(1) = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    S5Lim=>MyKinSet%AllKinConfig(2)
    TCLim=>MyKinSet%AllKinConfig(3)
    TCS5Lim=>MyKinSet%AllKinConfig(4)
    C54Lim=>MyKinSet%AllKinConfig(5)
    C54S5Lim=>MyKinSet%AllKinConfig(6)
    TCC54Lim=>MyKinSet%AllKinConfig(7)
    TCC54S5Lim=>MyKinSet%AllKinConfig(8)

    xx = zero
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*(xx(xRHO2)**buffrr_pow).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4151d(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S5Lim,TCLim,TCS5Lim,&
         C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim)
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_gq_qg(HardProc%AmpMom,res_jj)
       respdf = gq_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not


    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5

    xi1 = S5Lim%PartFrac(1)
    xi2 = S5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,S5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.S5Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(2) = zero
#endif

       kin(2) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(S5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       reseik = Eik_s5(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3),S5Lim%LimMom(:,4)) 
       kin(2)  = (-one)*respdf*reseik*S5Lim%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(2) = kin(2)
#endif

       call HistoFill(5,binarray,weightarray,kin(2),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    endif

    !-- 3 -> C54S5

    xi1 = C54S5Lim%PartFrac(1)
    xi2 = C54S5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C54S5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C54S5Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(3) = zero
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C54S5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       s45 = C54S5Lim%KinInv(3) 
       z5 =  C54S5Lim%KinInv(5)
       !
       kin(3) =  four*Cf/z5/s45 * respdf * C54S5Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(3) = kin(3)
#endif

       call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    endif


    !-- 4 -> C5

    xi1 = C54Lim%PartFrac(1)
    xi2 = C54Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C54Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C54Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4) = zero       
#endif

       kin(4) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C54Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       s45 = C54Lim%KinInv(3) 
       z5  = C54Lim%KinInv(5) 
       !
       kin(4) = (-two)/s45*Pqq(z5)*respdf*C54Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4) = kin(4)
#endif

       call HistoFill(5,binarray,weightarray,kin(4),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    endif


    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!

    !-- 5 -> TC
    !-- 6 -> TCC54

    xi1 = TCLim%PartFrac(1)
    xi2 = TCLim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCLim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCLim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(5:6) = zero       
#endif

       kin(5) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCLim%AmpMom,res)
       respdf = gq_lumi_lo(res,f1,f2)

       s14 = TCLim%KinInv(1)
       s15 = TCLim%KinInv(2)
       s45 = TCLim%KinInv(3)
       !
       z4 = TCLim%KinInv(4)
       z5 = TCLim%KinInv(5)
       z2 = one-z4-z5

       FintNNLO_wp_gq(5) = (-one)*respdf*(Pggq_flip(s15,s45,s14,z5,z2,z4))*TCLim%wgt*jac*pref 

       s1_45 = TCC54Lim%KinInv(1)
       s45 = TCC54Lim%KinInv(3)

       z1 = TCC54Lim%KinInv(4)
       z2 = TCC54Lim%KinInv(5)

       FintNNLO_wp_gq(6) = four/s45/s1_45 * Pqg(z1) * Pqg_flip(one-z2)
       FintNNLO_wp_gq(6) =  FintNNLO_wp_gq(6) *respdf * TCC54Lim%wgt*jac*pref

       kin(5) = sum(FintNNLO_wp_gq(5:6))
       call HistoFill(5,binarray,weightarray,kin(5),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    endif

    !-- 7 -> TCS5
    !-- 8 -> TCC54S5

    xi1 = TCS5Lim%PartFrac(1)
    xi2 = TCS5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCS5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCS5Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(7:8) = zero       
#endif

       kin(6) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCS5Lim%AmpMom,res)
       respdf = gq_lumi_lo(res,f1,f2)

       FintNNLO_wp_gq(7) = Pggq_s5(TCS5Lim%KinInv)
       FintNNLO_wp_gq(7) = FintNNLO_wp_gq(7)*respdf*TCS5Lim%wgt*jac*pref

       s14 = TCC54S5Lim%KinInv(1)
       s15 = TCC54S5Lim%KinInv(2)
       s45 = TCC54S5Lim%KinInv(3)
       z4 = TCC54S5Lim%KinInv(4)
       z5 = TCC54S5Lim%KinInv(5)
       FintNNLO_wp_gq(8) = -(four*Cf/s45/z5*two/s14*Pqg_flip(one-z4)) *respdf * TCC54S5Lim%wgt*jac*pref

       kin(6) = sum(FintNNLO_wp_gq(7:8))

    endif
    call HistoFill(10,binarray,weightarray,kin(6),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,8
          print *, i, FintNNLO_wp_gq(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_gq = FintNNLO_wp_gq
#endif

    ! do i = 1,8
    !    print *, i,FintNNLO_wp_gq(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return
    
  end function xsect_nnlo_4151d_rr_wp_gq
  
  function xsect_nnlo_4152_rr_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4152_rr_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,S5Lim,C4Lim,C4S5Lim,&
         C5Lim,C5S5Lim,C4C5Lim,C4C5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b,mylumi
    complex(dp) :: respdf_sping(-1:1,-1:1),respdf_sping2(-1:1,-1:1,-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping2(namp_lo,-1:1,-1:1,-1:1,-1:1),res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    complex(dp) :: hell_mat42(-1:1,-1:1),hell_mat51(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s24,s25,s45,kt1kt2sq,reseik
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_gq(8),kin(1:6)
    logical  :: makecut

    callcount = callcount + 1

    xsect_nnlo_4152_rr_wp_gq = 0
    ff(1) = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    S5Lim=>MyKinSet%AllKinConfig(2)
    C5Lim=>MyKinSet%AllKinConfig(3)
    C5S5Lim=>MyKinSet%AllKinConfig(4)
    C4Lim=>MyKinSet%AllKinConfig(5)
    C4S5Lim=>MyKinSet%AllKinConfig(6)
    C4C5Lim=>MyKinSet%AllKinConfig(7)
    C4C5S5Lim=>MyKinSet%AllKinConfig(8)

    xx = zero
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4152(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S5Lim,C4Lim,C4S5Lim,&
         C5Lim,C5S5Lim,C4C5Lim,C4C5S5Lim)  
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_gq_qg(HardProc%AmpMom,res_jj)
       respdf = gq_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    
    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5
    !-- 3 -> S5C5

    xi1 = S5Lim%PartFrac(1) 
    xi2 = S5Lim%PartFrac(2) 
    call cut_histo_w(xi1,xi2,5,S5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.S5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(2) = zero
       FintNNLO_wp_gq(3) = zero
#endif       

       kin(2) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(S5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       reseik = Eik_s5(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3),S5Lim%LimMom(:,4)) 
       FintNNLO_wp_gq(2)  = (-one)*respdf*reseik*S5Lim%wgt*jac*pref

       s25 = C5S5Lim%KinInv(2)
       z5 = C5S5Lim%KinInv(5)
       FintNNLO_wp_gq(3) =  four*Cf/z5/s25 * respdf * C5S5Lim%wgt*jac*pref

       kin(2) = FintNNLO_wp_gq(2) + FintNNLO_wp_gq(3)
       call HistoFill(5,binarray,weightarray,kin(2),vegasweight,iternumber)

    endif

    !-- 4 -> C4

    xi1 = C4Lim%PartFrac(1) 
    xi2 = C4Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C4Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C4Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4)  = zero
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C4Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_ns(res_j,f1,f2)

       s24 = C4Lim%KinInv(1) 
       z4 = C4Lim%KinInv(4)

       kin(3) = (two/s24)*respdf*Pqg_flip(one-z4)*C4Lim%wgt*jac*pref

!        call restree_j_sping(2,C4Lim%AmpMom,res_sping_j,res_j)
!        respdf_sping(:,:) = res_sping_j(1,:,:) * get_lumi_wp_qg(f1,f2)

!        s24 = C4Lim%KinInv(1) 
!        z4 = C4Lim%KinInv(2)

!        nperp=C4Lim%LimMom(2:4,1)
!        ncoll=C4Lim%LimMom(2:4,2)

!        hell_mat = comp_hell_mat_pgq_flip(z4,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
!        call zipp(respdf_sping,hell_mat,respdf)

!        kin(3) = two/s24*respdf*C4Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4)  = kin(3)
#endif

      call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber)
    endif

    !-- 5 -> C5

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(5)  = zero
#endif

       kin(4) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       s15 = C5Lim%KinInv(2) 
       z1  = C5Lim%KinInv(5)

       kin(4) = -two/s15*Pqq(z1)*respdf*C5Lim%wgt*jac*pref

!        call restree_j_sping(1,C5Lim%AmpMom,res_sping_j,res_j)
!        respdf_sping(:,:) = res_sping_j(3,:,:) * get_lumi_wp_qg(f1,f2)

!        s15 = C5Lim%KinInv(1) 
!        z1  = C5Lim%KinInv(2)

!        nperp=C5Lim%LimMom(2:4,1)
!        ncoll=C5Lim%LimMom(2:4,2)

!        hell_mat = comp_hell_mat_pgg(z1,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
!        call zipp(respdf_sping,hell_mat,respdf)

!        kin(4) = two/s15*respdf*C5Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(5)  = kin(4)
#endif

      call HistoFill(5,binarray,weightarray,kin(4),vegasweight,iternumber)
    endif


    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!

    !-- 6 -> C4S5
    !-- 7 -> C4C5S5

    xi1 = C4S5Lim%PartFrac(1) 
    xi2 = C4S5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,C4S5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C4S5Lim%flag) then

       kin(5) = zero

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(6:7)  = zero
#endif

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(C4S5Lim%AmpMom,res)
       respdf = gq_lumi_lo(res,f1,f2)

       s12 = C4S5Lim%EikInv(1) 
       s25 = C4S5Lim%EikInv(2) 
       s15 = C4S5Lim%EikInv(3) 
       s14 = C4S5Lim%KinInv(1) 
       z4  = C4S5Lim%KinInv(4) 
       !
       FintNNLO_wp_gq(6) =  -8.0_dp*Cf/s14 * s12/s15/s25 * Pqg_flip(one-z4) * respdf*C4S5Lim%wgt*jac*pref

       s14 = C4C5S5Lim%KinInv(1) 
       s25 = C4C5S5Lim%KinInv(2) 
       z4 =  C4C5S5Lim%KinInv(4) 
       z5 =  C4C5S5Lim%KinInv(5) 
       !
       FintNNLO_wp_gq(7) =  8.0_dp*Cf/s14/s25/z5  * Pqg_flip(one-z4)  *respdf * C4C5S5Lim%wgt*jac*pref

       kin(5) = sum(FintNNLO_wp_gq(6:7))

       call HistoFill(5,binarray,weightarray,kin(5),vegasweight,iternumber)

    endif


    !-- 8 -> C4C5

    xi1 = C4C5Lim%PartFrac(1) 
    xi2 = C4C5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,C4C5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C4C5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(8)  = zero
#endif

       kin(6) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(C4C5Lim%AmpMom,res)
       respdf = gq_lumi_lo(res,f1,f2)
       
       s24 = C4C5Lim%KinInv(1) 
       s15 = C4C5Lim%KinInv(2) 
       z4 =  C4C5Lim%KinInv(4) 
       z5 =  C4C5Lim%KinInv(5) 

       kin(6) =  four/s24/s15 * Pqq(z5)  * Pqg_flip(one-z4)  *respdf * C4C5Lim%wgt*jac*pref

!        call restree_sping2(C4C5Lim%AmpMom,res_sping2)
!        respdf_sping2(:,:,:,:) = res_sping2(1,:,:,:,:) * get_lumi_wp_gq(f1,f2)

!        s24 = C4C5Lim%KinInv(1) 
!        s15 = C4C5Lim%KinInv(2) 
!        z4 =  C4C5Lim%KinInv(4) 
!        z5 =  C4C5Lim%KinInv(5) 
!        !

!        nperp=C4C5Lim%LimMom(2:4,1)
!        ncoll=C4C5Lim%LimMom(2:4,2)
!        hell_mat42 = comp_hell_mat_pgq_flip(z4,nperp,ncoll) ! n1: correlation vector, n2: collinear direction

!        nperp=C4C5Lim%LimMom(2:4,3)
!        ncoll=C4C5Lim%LimMom(2:4,4)
!        hell_mat51 = comp_hell_mat_pgg(z5,nperp,ncoll) ! n1: correlation vector, n2: collinear direction

!        call zipp2(hell_mat51,respdf_sping2,hell_mat42,respdf)

!        kin(6) =  +four/s24/s15 * respdf * C4C5Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(8)  = kin(6)
#endif

    endif
    call HistoFill(10,binarray,weightarray,kin(6),vegasweight,iternumber)

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,8
          print *, i, FintNNLO_wp_gq(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_gq = FintNNLO_wp_gq
#endif

!   if (abs(ff(1)) .gt. 1d10) then
!      print *, "xX",xx
!      do i = 1,8
!         print *, i,FintNNLO_wp_gq(i)
!      enddo
!      print *, 'ff(1)', ff(1)
!      pause
!   endif

    return

  end function xsect_nnlo_4152_rr_wp_gq

  !-- 

  function xsect_nnlo_4252a_rr_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4252a_rr_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,S5Lim,TCLim,TCS5Lim,C5Lim,C5S5Lim,TCC5Lim,TCC5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b,mylumi
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z4,z5,s14,s15,s24,s25,s45,s4_15,kt1kt2sq,reseik
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_gq(8),kin(1:3)
    logical  :: makecut

    callcount = callcount + 1

    xsect_nnlo_4252a_rr_wp_gq = 0
    ff(1) = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    S5Lim=>MyKinSet%AllKinConfig(2)
    TCLim=>MyKinSet%AllKinConfig(3)
    TCS5Lim=>MyKinSet%AllKinConfig(4)
    C5Lim=>MyKinSet%AllKinConfig(5)
    C5S5Lim=>MyKinSet%AllKinConfig(6)
    TCC5Lim=>MyKinSet%AllKinConfig(7)
    TCC5S5Lim=>MyKinSet%AllKinConfig(8)

    xx = zero
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4252a(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S5Lim,TCLim,TCS5Lim,&
         C5Lim,C5S5Lim,TCC5Lim,TCC5S5Lim)
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_gq_qg(HardProc%AmpMom,res_jj)
       respdf = gq_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5
    !-- 3 -> S5C5

    xi1 = S5Lim%PartFrac(1) 
    xi2 = S5Lim%PartFrac(2) 
    call cut_histo_w(xi1,xi2,5,S5Lim%AmpMom,makecut,binarray,weightarray,mu) 

    if (makecut.or.S5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(2) = zero
       FintNNLO_wp_gq(3) = zero
#endif       

       kin(2) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(S5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       reseik = Eik_s5(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3),S5Lim%LimMom(:,4)) 
       FintNNLO_wp_gq(2)  = (-one)*respdf*reseik*S5Lim%wgt*jac*pref

       s15 = -C5S5Lim%KinInv(2)
       z5 = C5S5Lim%KinInv(5)
       FintNNLO_wp_gq(3) =  four*Cf/z5/s15 * respdf * C5S5Lim%wgt*jac*pref

       kin(2) = FintNNLO_wp_gq(2) + FintNNLO_wp_gq(3)
       call HistoFill(5,binarray,weightarray,kin(2),vegasweight,iternumber)

    endif

    !! -----------------------------!!
    !!     calling NLO ME  w sp corr!!
    !! -----------------------------!!

    !-- 4 -> C5

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4)  = zero
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       s25 = -C5Lim%KinInv(2) 
       z1  = C5Lim%KinInv(5)

       kin(3) = (-one)*(-two/s25)*Pqg(z1)*respdf*C5Lim%wgt*jac*pref

!        call restree_j_sping(1,C5Lim%AmpMom,res_sping_j,res_j)
!        respdf_sping = res_sping_j(3,:,:) * get_lumi_wp_gq(f1,f2)

!        s15 = C5Lim%KinInv(2) 
!        z1  = C5Lim%KinInv(5)

!        nperp=C5Lim%LimMom(2:4,1)
!        ncoll=C5Lim%LimMom(2:4,2)

!        hell_mat = comp_hell_mat_pgg(z1,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
!        call zipp(respdf_sping,hell_mat,respdf)

!        kin(3) = (-one)*(-two/s15)*respdf*C5Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4)  = kin(3)
       FintNNLO_wp_gq(5:8) = zero
#endif

    endif
    call HistoFill(10,binarray,weightarray,kin(3),vegasweight,iternumber)

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,8
          print *, i, FintNNLO_wp_gq(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_gq = FintNNLO_wp_gq
#endif

    ! do i = 1,8
    !    print *, i,FintNNLO_wp_gq(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return
    
  end function xsect_nnlo_4252a_rr_wp_gq

  function xsect_nnlo_4252b_rr_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4252b_rr_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,S5Lim,TCLim,TCS5Lim,&
         C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z3,z4,z5,z45,s12,s14,s24,s25,s45,s2_45,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_gq(8),kin(1:3)
    logical :: makecut

    callcount = callcount + 1

    xsect_nnlo_4252b_rr_wp_gq = 0
    ff(1) = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    S5Lim=>MyKinSet%AllKinConfig(2)
    TCLim=>MyKinSet%AllKinConfig(3)
    TCS5Lim=>MyKinSet%AllKinConfig(4)
    C54Lim=>MyKinSet%AllKinConfig(5)
    C54S5Lim=>MyKinSet%AllKinConfig(6)
    TCC54Lim=>MyKinSet%AllKinConfig(7)
    TCC54S5Lim=>MyKinSet%AllKinConfig(8)

    xx = zero
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*(xx(xRHO2)**buffrr_pow).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4252b(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S5Lim,TCLim,TCS5Lim,&
         C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim)
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_gq_qg(HardProc%AmpMom,res_jj)
       respdf = gq_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not


    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5
    !-- 3 -> C54S5

    xi1 = S5Lim%PartFrac(1)
    xi2 = S5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,S5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.S5Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(2:3) = zero       
#endif

       kin(2) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(S5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       reseik = Eik_s5(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3),S5Lim%LimMom(:,4)) 
       FintNNLO_wp_gq(2)  = (-one)*respdf*reseik*S5Lim%wgt*jac*pref

       s45 = C54S5Lim%KinInv(3) 
       z5 =  C54S5Lim%KinInv(5)
       !
       FintNNLO_wp_gq(3) =  four*Cf/z5/s45 * respdf * C54S5Lim%wgt*jac*pref

       kin(2) = sum(FintNNLO_wp_gq(2:3))
       call HistoFill(5,binarray,weightarray,kin(2),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    endif

    !-- 4 -> C5

    xi1 = C54Lim%PartFrac(1)
    xi2 = C54Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C54Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C54Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4) = zero       
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C54Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       s45 = C54Lim%KinInv(3) 
       z5  = C54Lim%KinInv(5) 
       !
       kin(3) = (-two)/s45*Pqq(z5)*respdf*C54Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4) = kin(3)
       FintNNLO_wp_gq(5:8) = zero
#endif

    endif
    call HistoFill(10,binarray,weightarray,kin(3),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,8
          print *, i, FintNNLO_wp_gq(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_gq = FintNNLO_wp_gq
#endif

    ! do i = 1,8
    !    print *, i,FintNNLO_wp_gq(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return

  end function xsect_nnlo_4252b_rr_wp_gq

  function xsect_nnlo_4252c_rr_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4252c_rr_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,S5Lim,TCLim,TCS5Lim,C4Lim,C4S5Lim,TCC4Lim,TCC4S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b,mylumi
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s24,s25,s45,s5_24,kt1kt2sq,reseik
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_gq(8),kin(1:3)
    logical  :: makecut

    callcount = callcount + 1

    xsect_nnlo_4252c_rr_wp_gq = 0
    ff(1) = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    S5Lim=>MyKinSet%AllKinConfig(2)
    TCLim=>MyKinSet%AllKinConfig(3)
    TCS5Lim=>MyKinSet%AllKinConfig(4)
    C4Lim=>MyKinSet%AllKinConfig(5)
    C4S5Lim=>MyKinSet%AllKinConfig(6)
    TCC4Lim=>MyKinSet%AllKinConfig(7)
    TCC4S5Lim=>MyKinSet%AllKinConfig(8)

    xx = zero
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4252c(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S5Lim,TCLim,TCS5Lim,&
         C4Lim,C4S5Lim,TCC4Lim,TCC4S5Lim)
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_gq_qg(HardProc%AmpMom,res_jj)
       respdf = gq_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    !-- 2 -> C4 --> here for H, not here for WP
    kin(2) = zero
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(2) = zero       
#endif

!     xi1 = C4Lim%PartFrac(1)
!     xi2 = C4Lim%PartFrac(2)
!     call cut_histo_w(xi1,xi2,5,C4Lim%AmpMom,makecut,binarray,weightarray,mu)

!     if (makecut.or.C4Lim%flag) then
       
! #if (_DebugFeatures == 1)
!        FintNNLO_wp_gq(2) = zero       
! #endif

!        kin(2) = zero

!     else

!        call setscale(mu,muF,muR)
!        call get_pdf(xi1,xi2,muF,f1,f2)
!        call get_prefactor(muR,2,pref)

!        call restree_j(C4Lim%AmpMom,res_j)
!        respdf = res_j(4) * get_lumi_wp_gq(f1,f2)

!        s14 = C4Lim%KinInv(2) 
!        z4 = C4Lim%KinInv(5)

!        kin(2) = (two/s14)*Pqg_flip(one-z4)*respdf*C4Lim%wgt*jac*pref

! #if (_DebugFeatures == 1)
!        FintNNLO_wp_gq(2) = kin(2)
! #endif

!        call HistoFill(5,binarray,weightarray,kin(2),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

!     endif

    !-- 3 -> S5

    xi1 = S5Lim%PartFrac(1)
    xi2 = S5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,S5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.S5Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(3) = zero       
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(S5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       reseik = Eik_s5(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3),S5Lim%LimMom(:,4)) 

       kin(3) = (-one)*respdf*reseik*S5Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(3) = kin(3)
       FintNNLO_wp_gq(4:8) = zero
#endif

    endif
    call HistoFill(10,binarray,weightarray,kin(3),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,8
          print *, i, FintNNLO_wp_gq(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_gq(1) = FintNNLO_wp_gq(1)
    FintNNLO_CheckLim_wp_gq(2) = FintNNLO_wp_gq(3) !-- S5
    FintNNLO_CheckLim_wp_gq(3) = FintNNLO_wp_gq(7) !-- C4S5
    FintNNLO_CheckLim_wp_gq(4) = FintNNLO_wp_gq(2) !-- C4
    FintNNLO_CheckLim_wp_gq(5) = FintNNLO_wp_gq(4) !-- TC
    FintNNLO_CheckLim_wp_gq(6) = FintNNLO_wp_gq(5) !-- TCC4
    FintNNLO_CheckLim_wp_gq(7) = FintNNLO_wp_gq(6) !-- TCS5
    FintNNLO_CheckLim_wp_gq(8) = FintNNLO_wp_gq(8) !-- TCC4S5
#endif

    ! do i = 1,8
    !    print *, i,FintNNLO_wp_gq(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return

  end function xsect_nnlo_4252c_rr_wp_gq

  function xsect_nnlo_4252d_rr_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4252d_rr_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,S5Lim,TCLim,TCS5Lim,&
         C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z3,z4,z5,z45,s12,s14,s24,s25,s45,s2_45,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_gq(8),kin(1:4)
    logical :: makecut

    callcount = callcount + 1

    xsect_nnlo_4252d_rr_wp_gq = 0
    ff(1) = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    S5Lim=>MyKinSet%AllKinConfig(2)
    TCLim=>MyKinSet%AllKinConfig(3)
    TCS5Lim=>MyKinSet%AllKinConfig(4)
    C54Lim=>MyKinSet%AllKinConfig(5)
    C54S5Lim=>MyKinSet%AllKinConfig(6)
    TCC54Lim=>MyKinSet%AllKinConfig(7)
    TCC54S5Lim=>MyKinSet%AllKinConfig(8)

    xx = zero
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*(xx(xRHO2)**buffrr_pow).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4252d(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S5Lim,TCLim,TCS5Lim,&
         C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim)
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_gq_qg(HardProc%AmpMom,res_jj)
       respdf = gq_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not


    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5

    xi1 = S5Lim%PartFrac(1)
    xi2 = S5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,S5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.S5Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(2) = zero
#endif

       kin(2) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(S5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       reseik = Eik_s5(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3),S5Lim%LimMom(:,4)) 
       kin(2)  = (-one)*respdf*reseik*S5Lim%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(2) = kin(2)
#endif

       call HistoFill(5,binarray,weightarray,kin(2),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    endif

    !-- 3 -> C54S5

    xi1 = C54S5Lim%PartFrac(1)
    xi2 = C54S5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C54S5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C54S5Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(3) = zero
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C54S5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       s45 = C54S5Lim%KinInv(3) 
       z5 =  C54S5Lim%KinInv(5)
       !
       kin(3) =  four*Cf/z5/s45 * respdf * C54S5Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(3) = kin(3)
#endif

       call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    endif


    !-- 4 -> C5

    xi1 = C54Lim%PartFrac(1)
    xi2 = C54Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C54Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C54Lim%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4) = zero       
#endif

       kin(4) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C54Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       s45 = C54Lim%KinInv(3) 
       z5  = C54Lim%KinInv(5) 
       !
       kin(4) = (-two)/s45*Pqq(z5)*respdf*C54Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4) = kin(4)
       FintNNLO_wp_gq(5:8) = zero
#endif

    endif
    call HistoFill(10,binarray,weightarray,kin(4),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,8
          print *, i, FintNNLO_wp_gq(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_gq = FintNNLO_wp_gq
#endif

    ! do i = 1,8
    !    print *, i,FintNNLO_wp_gq(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return

  end function xsect_nnlo_4252d_rr_wp_gq



  function xsect_nnlo_4251_rr_wp_gq(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4251_rr_wp_gq,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,S5Lim,C4Lim,C4S5Lim,&
         C5Lim,C5S5Lim,C4C5Lim,C4C5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b,mylumi
    complex(dp) :: respdf_sping(-1:1,-1:1),respdf_sping2(-1:1,-1:1,-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping2(namp_lo,-1:1,-1:1,-1:1,-1:1),res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    complex(dp) :: hell_mat42(-1:1,-1:1),hell_mat51(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s24,s25,s45,kt1kt2sq,reseik
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_gq(8),kin(1:4)
    logical  :: makecut

    callcount = callcount + 1

    xsect_nnlo_4251_rr_wp_gq = 0
    ff(1) = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    S5Lim=>MyKinSet%AllKinConfig(2)
    C5Lim=>MyKinSet%AllKinConfig(3)
    C5S5Lim=>MyKinSet%AllKinConfig(4)
    C4Lim=>MyKinSet%AllKinConfig(5)
    C4S5Lim=>MyKinSet%AllKinConfig(6)
    C4C5Lim=>MyKinSet%AllKinConfig(7)
    C4C5S5Lim=>MyKinSet%AllKinConfig(8)

    xx = zero
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4251(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S5Lim,C4Lim,C4S5Lim,&
         C5Lim,C5S5Lim,C4C5Lim,C4C5S5Lim)  
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_gq_qg(HardProc%AmpMom,res_jj)
       respdf = gq_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    
    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5
    !-- 3 -> S5C5

    xi1 = S5Lim%PartFrac(1) 
    xi2 = S5Lim%PartFrac(2) 
    call cut_histo_w(xi1,xi2,5,S5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.S5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(2) = zero
       FintNNLO_wp_gq(3) = zero
#endif       

       kin(2) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(S5Lim%AmpMom,res_sping_j,res_j)
       respdf = gq_lumi_nlo(res_j,f1,f2)

       reseik = Eik_s5(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3),S5Lim%LimMom(:,4)) 
       FintNNLO_wp_gq(2)  = (-one)*respdf*reseik*S5Lim%wgt*jac*pref

       s15 = C5S5Lim%KinInv(2)
       z5 = C5S5Lim%KinInv(5)
       FintNNLO_wp_gq(3) =  four*Ca/z5/s15 * respdf * C5S5Lim%wgt*jac*pref

       kin(2) = FintNNLO_wp_gq(2) + FintNNLO_wp_gq(3)
       call HistoFill(5,binarray,weightarray,kin(2),vegasweight,iternumber)

    endif

    !-- 4 -> C4 -- Here for Higgs, not here for WP
    kin(3) = zero
#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(4)  = zero
#endif

!     xi1 = C4Lim%PartFrac(1) 
!     xi2 = C4Lim%PartFrac(2)
!     call cut_histo_w(xi1,xi2,5,C4Lim%AmpMom,makecut,binarray,weightarray,mu)

!     if (makecut.or.C4Lim%flag) then

! #if (_DebugFeatures == 1)
!        FintNNLO_wp_gq(4)  = zero
! #endif

!        kin(3) = zero

!     else

!        call setscale(mu,muF,muR)
!        call get_pdf(xi1,xi2,muF,f1,f2)
!        call get_prefactor(muR,2,pref)

!        call restree_j(C4Lim%AmpMom,res_j)
!        respdf = res_j(4) * get_lumi_wp_qg(f1,f2)

!        s14 = C4Lim%KinInv(1) 
!        z4 = C4Lim%KinInv(2)

!        kin(3) = two/s14*respdf*Pqg_flip(one-z4)*C4Lim%wgt*jac*pref

! #if (_DebugFeatures == 1)
!        FintNNLO_wp_gq(4)  = kin(3)
! #endif

!       call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber)
!     endif

    !-- 5 -> C5

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(5)  = zero
#endif

       kin(4) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

!---
       call restree_j(C5Lim%AmpMom,res_sping_j,res_j)
       respdf_sping = gq_lumi_qg_sping(res_sping_j,f1,f2)

       s15 = C5Lim%KinInv(1) 
       z2 = C5Lim%KinInv(5)

       nperp=C5Lim%LimMom(2:4,1)
       ncoll=C5Lim%LimMom(2:4,2)

       hell_mat = comp_hell_mat_pgg(z2,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
       call zipp(respdf_sping,hell_mat,respdf)

       kin(4) = (-one)*(-two/s15)*respdf*C5Lim%wgt*jac*pref
!---

!        call restree_j(C5Lim%AmpMom,res_j)
!        respdf = res_j(3) * get_lumi_wp_gq(f1,f2)

!        s25 = C5Lim%KinInv(1) 
!        z2  = C5Lim%KinInv(2)

!        kin(4) = two/s25*respdf*Pqq(z2)*C5Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_gq(5)  = kin(4)
       FintNNLO_wp_gq(6:8) = zero
#endif

    endif
    call HistoFill(10,binarray,weightarray,kin(4),vegasweight,iternumber)

    ff(1) = sum(kin)
    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,8
          print *, i, FintNNLO_wp_gq(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_gq = FintNNLO_wp_gq
#endif

    ! do i = 1,8
    !    print *, i,FintNNLO_wp_qg(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return
    
  end function xsect_nnlo_4251_rr_wp_gq

end module mod_xsects_nnlo_rr_wp_gq

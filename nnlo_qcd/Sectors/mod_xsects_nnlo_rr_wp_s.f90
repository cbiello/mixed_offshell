module mod_xsects_nnlo_rr_wp_s
  use mod_types; use common_def; use mod_consts_dp
#if (_DebugFeatures == 1)
  use mod_consts_debug
#endif
  use mod_kinematics_gen; use mod_kinematicsNNLO_EordAord
  use mod_process
  use mod_auxfunctions; use mod_procfunctions
  use mod_histograms
  use mod_lumi
  use mod_auxlumi_wp
  use mod_amplitudes_tree_wp
  use mod_limits_dy_ns
  use mod_splitting
  implicit none
  integer, parameter :: namp_lo   = 2
  integer, parameter :: namp_nlo  = 6
  integer, parameter :: namp_nnlo = 26

  private

  public :: xsect_nnlo_4151a_rr_wp_s, xsect_nnlo_4151b_rr_wp_s, xsect_nnlo_4151c_rr_wp_s, xsect_nnlo_4151d_rr_wp_s
  public :: xsect_nnlo_4152_rr_wp_s
  !
  public :: xsect_nnlo_4252a_rr_wp_s, xsect_nnlo_4252b_rr_wp_s, xsect_nnlo_4252c_rr_wp_s, xsect_nnlo_4252d_rr_wp_s
  public :: xsect_nnlo_4251_rr_wp_s


contains

  function xsect_nnlo_4151a_rr_wp_s(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4151a_rr_wp_s,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C5Lim,C5DSLim,C5DSS5Lim,C5S5Lim,TCC5Lim,TCC5DSLim,TCC5DSS5Lim,TCC5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z4,z5,s14,s15,s45,s4_15,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_s(16),kin(1:6)
    logical  :: makecut

    callcount = callcount + 1

    xsect_nnlo_4151a_rr_wp_s = 0
    ff(1) = zero
    xx = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    DSLim=>MyKinSet%AllKinConfig(2)
    S5Lim=>MyKinSet%AllKinConfig(3)
    DSS5Lim=>MyKinSet%AllKinConfig(4)
    TCLim=>MyKinSet%AllKinConfig(5)
    TCDSLim=>MyKinSet%AllKinConfig(6)
    TCS5Lim=>MyKinSet%AllKinConfig(7)
    TCDSS5Lim=>MyKinSet%AllKinConfig(8)
    C5Lim=>MyKinSet%AllKinConfig(9)
    C5DSLim=>MyKinSet%AllKinConfig(10)
    C5S5Lim=>MyKinSet%AllKinConfig(11)
    C5DSS5Lim=>MyKinSet%AllKinConfig(12)
    TCC5Lim=>MyKinSet%AllKinConfig(13)
    TCC5DSLim=>MyKinSet%AllKinConfig(14)
    TCC5S5Lim=>MyKinSet%AllKinConfig(15)
    TCC5DSS5Lim=>MyKinSet%AllKinConfig(16)

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4151a(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
         C5Lim,C5DSLim,C5S5Lim,C5DSS5Lim,TCC5Lim,TCC5DSLim,TCC5S5Lim,TCC5DSS5Lim)
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!

    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)
       call restree_s_jj(HardProc%AmpMom,res_jj)
       respdf = s_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    if (vn_tagc+vn_tagce .ne. 0) then    

    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5
    !-- 3 -> S5C5

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(2) = zero
    FintNNLO_wp_s(3) = zero
#endif       

    kin(2) = zero

    !! -----------------------------!!
    !!     calling NLO ME  w sp corr!!
    !! -----------------------------!!
    
    !-- 4 -> C5

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(4)  = zero
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       !-- first argument: which parton is correlated (1,2 --> initial, 3 --> final)
       call restree_j(C5Lim%AmpMom,res_sping_j,res_j)
       respdf_sping = s_lumi_gq_sping(res_sping_j,f1,f2)

       s15 = C5Lim%KinInv(2) 
       z1  = C5Lim%KinInv(5)

       nperp=C5Lim%LimMom(2:4,1)
       ncoll=C5Lim%LimMom(2:4,2)

       hell_mat = comp_hell_mat_pgq_flip(z1,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
       call zipp(respdf_sping,hell_mat,respdf)

       kin(3) = (-one)*(-two/s15)*(respdf*tag_s_s) * C5Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(4)  = kin(3)
#endif

      call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber)
    endif


    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!
    
    !-- 5 -> DS
    !-- 6 -> DSS5
    !-- 7 -> DSTC
    !-- 8 -> DSTCS5
    !-- 9 -> DSC5
    !-- 10 -> DSS5C5
    !-- 11 -> DSTCC5
    !-- 12 -> DSTCC5S5

    kin(4) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(5:12)  = zero
#endif

    !-- 13 -> TC
    !-- 14 -> TCC5

    xi1 = TCLim%PartFrac(1) 
    xi2 = TCLim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCLim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCLim%flag) then

       kin(5) = zero

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(13:14)  = zero
#endif

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       !-- no spin correlation for Hgg amplitude
       call restree(TCLim%AmpMom,res)
       respdf = s_lumi_xj(res,f1,f2)
       respdf_b = s_lumi_xj_id(res,f1,f2)

       s14 = TCLim%KinInv(1)
       s15 = TCLim%KinInv(2)
       s45 = TCLim%KinInv(3)
       !
       z4 = TCLim%KinInv(4)
       z5 = TCLim%KinInv(5)
       z1 = one-z4-z5

       FintNNLO_wp_s(13) = (-one)*( &
            (PqbqR(-s14,s45,-s15,z4,z1,z5)+PqbqR(-s15,s45,-s14,z5,z1,z4))*respdf*tag_s_s+ & !-- 4 <-> 5
             Pqbqq_id_symm(-s14,-s15,s45,z1,z4,z5)*respdf_b*tag_s_v &
            )*TCLim%wgt*jac*pref 

       s15 = TCC5Lim%KinInv(1)
       s4_15 = TCC5Lim%KinInv(2)
       z1 = TCC5Lim%KinInv(3)
       z2 = TCC5Lim%KinInv(4)
       kt1kt2sq = TCC5Lim%KinInv(5)

       FintNNLO_wp_s(14) = Ptcso_qbq_r(s15,s4_15,z1,z2,kt1kt2sq)*tag_s_s*respdf * TCC5Lim%wgt*jac*pref 

       kin(5) = sum(FintNNLO_wp_s(13:14))

    endif

    else
       kin(2:5) = zero
    endif

    call HistoFill(10,binarray,weightarray,kin(5),vegasweight,iternumber)

    !-- 15 -> TCS5
    !-- 16 -> TCC5S5

    kin(6) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(15:16)  = zero
#endif

     ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,16
          print *, i, FintNNLO_wp_s(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_ns = FintNNLO_wp_s
#endif

!     do i = 1,16
!        print *, i,FintNNLO_wp_s(i)
!     enddo
!     print *, 'ff(1)', ff(1)
!     pause

    return

  end function xsect_nnlo_4151a_rr_wp_s

  function xsect_nnlo_4151b_rr_wp_s(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4151b_rr_wp_s,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
          C54Lim,C54DSLim,C54DSS5Lim,C54S5Lim,TCC54Lim,TCC54DSLim,TCC54DSS5Lim,TCC54S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z3,z4,z5,z45,s12,s14,s24,s15,s45,s1_45,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_s(16),kin(1:6)
    logical :: makecut

    callcount = callcount + 1

    xsect_nnlo_4151b_rr_wp_s = 0
    ff(1) = zero
    xx = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    DSLim=>MyKinSet%AllKinConfig(2)
    S5Lim=>MyKinSet%AllKinConfig(3)
    DSS5Lim=>MyKinSet%AllKinConfig(4)
    TCLim=>MyKinSet%AllKinConfig(5)
    TCDSLim=>MyKinSet%AllKinConfig(6)
    TCS5Lim=>MyKinSet%AllKinConfig(7)
    TCDSS5Lim=>MyKinSet%AllKinConfig(8)
    C54Lim=>MyKinSet%AllKinConfig(9)
    C54DSLim=>MyKinSet%AllKinConfig(10)
    C54S5Lim=>MyKinSet%AllKinConfig(11)
    C54DSS5Lim=>MyKinSet%AllKinConfig(12)
    TCC54Lim=>MyKinSet%AllKinConfig(13)
    TCC54DSLim=>MyKinSet%AllKinConfig(14)
    TCC54S5Lim=>MyKinSet%AllKinConfig(15)
    TCC54DSS5Lim=>MyKinSet%AllKinConfig(16)

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*(xx(xRHO2)**buffrr_pow).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4151b(M2,ylab,xx(kNNLO_min:kNNLO_max_full),&
         HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
         C54Lim,C54DSLim,C54S5Lim,C54DSS5Lim,TCC54Lim,TCC54DSLim,TCC54S5Lim,TCC54DSS5Lim) 

    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_s_jj(HardProc%AmpMom,res_jj)
       respdf = s_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    if (vn_tagc+vn_tagce .ne. 0) then    

    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5
    !-- 3 -> S5C5

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(2) = zero
    FintNNLO_wp_s(3) = zero
#endif       

    kin(2) = zero


#if (_DebugFeatures == 1)
    FintNNLO_wp_s(4)  = zero
#endif

    kin(3) = zero


    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!
    
    !-- 5 -> DS
    !-- 6 -> DSS5
    !-- 7 -> DSTC
    !-- 8 -> DSTCS5
    !-- 9 -> DSC5
    !-- 10 -> DSS5C5
    !-- 11 -> DSTCC5
    !-- 12 -> DSTCC5S5

    kin(4) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(5:12)  = zero
#endif

    !-- 13 -> TC
    !-- 14 -> TCC5

    xi1 = TCLim%PartFrac(1) 
    xi2 = TCLim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCLim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCLim%flag) then

       kin(5) = zero

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(13:14)  = zero
#endif

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCLim%AmpMom,res)
       respdf = s_lumi_xj(res,f1,f2)
       respdf_b = s_lumi_xj_id(res,f1,f2)

       s14 = TCLim%KinInv(1)
       s15 = TCLim%KinInv(2)
       s45 = TCLim%KinInv(3)
       !
       z4 = TCLim%KinInv(4)
       z5 = TCLim%KinInv(5)
       z1 = one-z4-z5

       FintNNLO_wp_s(13) = (-one)*( &
            (PqbqR(-s14,s45,-s15,z4,z1,z5)+PqbqR(-s15,s45,-s14,z5,z1,z4))*respdf*tag_s_s+ & !-- 4 <-> 5
             Pqbqq_id_symm(-s14,-s15,s45,z1,z4,z5)*respdf_b*tag_s_v &
            )*TCLim%wgt*jac*pref 

       FintNNLO_wp_s(14) = zero
       
       kin(5) = sum(FintNNLO_wp_s(13:14))

    endif

    else
       kin(2:5) = zero
    endif

    call HistoFill(10,binarray,weightarray,kin(5),vegasweight,iternumber)

    !-- 15 -> TCS5
    !-- 16 -> TCC5S5

    kin(6) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(15:16)  = zero
#endif

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,16
          print *, i, FintNNLO_wp_s(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_ns = FintNNLO_wp_s
#endif

!     do i = 1,16
!        print *, i,FintNNLO_wp_s(i)
!     enddo
!     print *, 'ff(1)', ff(1)
!     pause

    return

  end function xsect_nnlo_4151b_rr_wp_s

  function xsect_nnlo_4151c_rr_wp_s(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4151c_rr_wp_s,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C4Lim,C4DSLim,C4DSS5Lim,C4S5Lim,TCC4Lim,TCC4DSLim,TCC4DSS5Lim,TCC4S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s25,s45,s5_14,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_s(16),kin(1:6)
    logical  :: makecut

    callcount = callcount + 1

    xsect_nnlo_4151c_rr_wp_s = 0
    ff(1) = zero
    xx = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    DSLim=>MyKinSet%AllKinConfig(2)
    S5Lim=>MyKinSet%AllKinConfig(3)
    DSS5Lim=>MyKinSet%AllKinConfig(4)
    TCLim=>MyKinSet%AllKinConfig(5)
    TCDSLim=>MyKinSet%AllKinConfig(6)
    TCS5Lim=>MyKinSet%AllKinConfig(7)
    TCDSS5Lim=>MyKinSet%AllKinConfig(8)
    C4Lim=>MyKinSet%AllKinConfig(9)
    C4DSLim=>MyKinSet%AllKinConfig(10)
    C4S5Lim=>MyKinSet%AllKinConfig(11)
    C4DSS5Lim=>MyKinSet%AllKinConfig(12)
    TCC4Lim=>MyKinSet%AllKinConfig(13)
    TCC4DSLim=>MyKinSet%AllKinConfig(14)
    TCC4S5Lim=>MyKinSet%AllKinConfig(15)
    TCC4DSS5Lim=>MyKinSet%AllKinConfig(16)

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
   call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4151c(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
         C4Lim,C4DSLim,C4S5Lim,C4DSS5Lim,TCC4Lim,TCC4DSLim,TCC4S5Lim,TCC4DSS5Lim)

    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_s_jj(HardProc%AmpMom,res_jj)
       respdf = s_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    if (vn_tagc+vn_tagce .ne. 0) then    

    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(2) = zero
#endif       

    kin(2) = zero

    !! -----------------------------!!
    !!     calling NLO ME  w sp corr!!
    !! -----------------------------!!
    
    !-- 4 -> C4

    xi1 = C4Lim%PartFrac(1) 
    xi2 = C4Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C4Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C4Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(4)  = zero
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C4Lim%AmpMom,res_sping_j,res_j)
       respdf_sping = s_lumi_gq_sping(res_sping_j,f1,f2)

       s14 = C4Lim%KinInv(2) 
       z1  = C4Lim%KinInv(5)

       nperp=C4Lim%LimMom(2:4,1)
       ncoll=C4Lim%LimMom(2:4,2)

       hell_mat = comp_hell_mat_pgq_flip(z1,nperp,ncoll)
       call zipp(respdf_sping,hell_mat,respdf)

       kin(3) = (-one)*(-two/s14)*(respdf*tag_s_s &
            ) * C4Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(4)  = kin(3)
#endif

      call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber)
    endif



    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!
    
    !-- 5 -> DS
    !-- 6 -> DSS5
    !-- 7 -> DSTC
    !-- 8 -> DSTCS5
    !-- 9 -> DSC4
    !-- 10 -> DSS5C4
    !-- 11 -> DSTCC4
    !-- 12 -> DSTCC4S5

    kin(4) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(5:12)  = zero
#endif

    !-- 13 -> TC
    !-- 14 -> TCC4

    xi1 = TCLim%PartFrac(1) 
    xi2 = TCLim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCLim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCLim%flag) then

       kin(5) = zero

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(13:14)  = zero
#endif

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCLim%AmpMom,res)
       respdf = s_lumi_xj(res,f1,f2)
       respdf_b = s_lumi_xj_id(res,f1,f2)

       s14 = TCLim%KinInv(1)
       s15 = TCLim%KinInv(2)
       s45 = TCLim%KinInv(3)
       !
       z4 = TCLim%KinInv(4)
       z5 = TCLim%KinInv(5)
       z1 = one-z4-z5

       FintNNLO_wp_s(13) = (-one)*( &
            (PqbqR(-s14,s45,-s15,z4,z1,z5)+PqbqR(-s15,s45,-s14,z5,z1,z4))*respdf*tag_s_s+ & !-- 4 <-> 5
             Pqbqq_id_symm(-s14,-s15,s45,z1,z4,z5)*respdf_b*tag_s_v &
            )*TCLim%wgt*jac*pref 

       s14 = TCC4Lim%KinInv(1) 
       s5_14 = TCC4Lim%KinInv(2) 
       z1  = TCC4Lim%KinInv(3) 
       z2 = TCC4Lim%KinInv(4)
       kt1kt2sq = TCC4Lim%KinInv(5)

       FintNNLO_wp_s(14) = Ptcso_qbq_r(s14,s5_14,z1,z2,kt1kt2sq)*tag_s_s*respdf * TCC4Lim%wgt*jac*pref 

       kin(5) = sum(FintNNLO_wp_s(13:14))

    endif

    else
       kin(2:5) = zero
    endif

    call HistoFill(10,binarray,weightarray,kin(5),vegasweight,iternumber)

    !-- 3 --> C4S5
    !-- 15 -> TCS5
    !-- 16 -> TCC4S5

    kin(6) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(3)  = zero
    FintNNLO_wp_s(15:16)  = zero
#endif

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,16
          print *, i, FintNNLO_wp_s(i)
       enddo
       ff(1) = zero
    endif


#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_ns = FintNNLO_wp_s
#endif

!     do i = 1,16
!        print *, i,FintNNLO_wp_s(i)
!     enddo
!     print *, 'ff(1)', ff(1)
!     pause

    return
    
  end function xsect_nnlo_4151c_rr_wp_s

  function xsect_nnlo_4151d_rr_wp_s(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4151d_rr_wp_s,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
          C54Lim,C54DSLim,C54DSS5Lim,C54S5Lim,TCC54Lim,TCC54DSLim,TCC54DSS5Lim,TCC54S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z3,z4,z5,z45,s12,s14,s24,s15,s45,s1_45,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_s(16),kin(1:7)
    logical :: makecut

    callcount = callcount + 1

    xsect_nnlo_4151d_rr_wp_s = 0
    ff(1) = zero
    xx = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    DSLim=>MyKinSet%AllKinConfig(2)
    S5Lim=>MyKinSet%AllKinConfig(3)
    DSS5Lim=>MyKinSet%AllKinConfig(4)
    TCLim=>MyKinSet%AllKinConfig(5)
    TCDSLim=>MyKinSet%AllKinConfig(6)
    TCS5Lim=>MyKinSet%AllKinConfig(7)
    TCDSS5Lim=>MyKinSet%AllKinConfig(8)
    C54Lim=>MyKinSet%AllKinConfig(9)
    C54DSLim=>MyKinSet%AllKinConfig(10)
    C54S5Lim=>MyKinSet%AllKinConfig(11)
    C54DSS5Lim=>MyKinSet%AllKinConfig(12)
    TCC54Lim=>MyKinSet%AllKinConfig(13)
    TCC54DSLim=>MyKinSet%AllKinConfig(14)
    TCC54S5Lim=>MyKinSet%AllKinConfig(15)
    TCC54DSS5Lim=>MyKinSet%AllKinConfig(16)

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
   call random_number(xx(kNNLO_max_full))
 
    if (xx(xE4)*xx(xE5)*xx(xRHO1)*(xx(xRHO2)**buffrr_pow).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4151d(M2,ylab,xx(kNNLO_min:kNNLO_max_full),&
         HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
         C54Lim,C54DSLim,C54S5Lim,C54DSS5Lim,TCC54Lim,TCC54DSLim,TCC54S5Lim,TCC54DSS5Lim) 

    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_s_jj(HardProc%AmpMom,res_jj)
       respdf = s_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    if (vn_tagc+vn_tagce .ne. 0) then    

    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(2) = zero
#endif       

    kin(2) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(3) = zero
#endif       

    kin(3) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(4)  = zero
#endif

    kin(4) = zero

    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!
    
    !-- 5 -> DS
    !-- 6 -> DSS5
    !-- 7 -> DSTC
    !-- 8 -> DSTCS5
    !-- 9 -> DSC5
    !-- 10 -> DSS5C5
    !-- 11 -> DSTCC5
    !-- 12 -> DSTCC5S5

    kin(5) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(5:12)  = zero
#endif

    !-- 13 -> TC
    !-- 14 -> TCC5

    xi1 = TCLim%PartFrac(1) 
    xi2 = TCLim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCLim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCLim%flag) then

       kin(6) = zero

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(13:14)  = zero
#endif

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCLim%AmpMom,res)
       respdf = s_lumi_xj(res,f1,f2)
       respdf_b = s_lumi_xj_id(res,f1,f2)

       s14 = TCLim%KinInv(1)
       s15 = TCLim%KinInv(2)
       s45 = TCLim%KinInv(3)
       !
       z4 = TCLim%KinInv(4)
       z5 = TCLim%KinInv(5)
       z1 = one-z4-z5

       FintNNLO_wp_s(13) = (-one)*( &
            (PqbqR(-s14,s45,-s15,z4,z1,z5)+PqbqR(-s15,s45,-s14,z5,z1,z4))*respdf*tag_s_s+ & !-- 4 <-> 5
             Pqbqq_id_symm(-s14,-s15,s45,z1,z4,z5)*respdf_b*tag_s_v &
            )*TCLim%wgt*jac*pref 


       FintNNLO_wp_s(14) = zero

       kin(6) = sum(FintNNLO_wp_s(13:14))

    endif

    else
       kin(2:6) = zero
    endif

    call HistoFill(10,binarray,weightarray,kin(6),vegasweight,iternumber)

    !-- 15 -> TCS5
    !-- 16 -> TCC5S5

    kin(7) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(15:16)  = zero
#endif

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,16
          print *, i, FintNNLO_wp_s(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_ns = FintNNLO_wp_s
#endif

!     do i = 1,16
!        print *, i,FintNNLO_wp_s(i)
!     enddo
!     print *, 'ff(1)', ff(1)
!     pause

    return

  end function xsect_nnlo_4151d_rr_wp_s

  function xsect_nnlo_4152_rr_wp_s(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4152_rr_wp_s,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,DSLim,DSS5Lim,S5Lim,C4Lim,C4DSLim,C4DSS5Lim,C4S5Lim,&
         C5Lim,C5DSLim,C5DSS5Lim,C5S5Lim,C4C5Lim,C4C5DSLim,C4C5DSS5Lim,C4C5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b
    complex(dp) :: respdf_sping(-1:1,-1:1),respdf_sping2(-1:1,-1:1,-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping2(namp_lo,-1:1,-1:1,-1:1,-1:1),res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    complex(dp) :: hell_mat41(-1:1,-1:1),hell_mat52(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s24,s25,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_s(16),kin(1:7)
    logical  :: makecut

    callcount = callcount + 1

    xsect_nnlo_4152_rr_wp_s = 0
    ff(1) = zero

    xx = zero
    
    HardProc=>MyKinSet%AllKinConfig(1)
    DSLim=>MyKinSet%AllKinConfig(2)
    S5Lim=>MyKinSet%AllKinConfig(3)
    DSS5Lim=>MyKinSet%AllKinConfig(4)
    C4Lim=>MyKinSet%AllKinConfig(5)
    C4DSLim=>MyKinSet%AllKinConfig(6)
    C4S5Lim=>MyKinSet%AllKinConfig(7)
    C4DSS5Lim=>MyKinSet%AllKinConfig(8)
    C5Lim=>MyKinSet%AllKinConfig(9)
    C5DSLim=>MyKinSet%AllKinConfig(10)
    C5S5Lim=>MyKinSet%AllKinConfig(11)
    C5DSS5Lim=>MyKinSet%AllKinConfig(12)
    C4C5Lim=>MyKinSet%AllKinConfig(13)
    C4C5DSLim=>MyKinSet%AllKinConfig(14)
    C4C5S5Lim=>MyKinSet%AllKinConfig(15)
    C4C5DSS5Lim=>MyKinSet%AllKinConfig(16)

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
   call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4152(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,DSLim,S5Lim,DSS5Lim,C4Lim,C4DSLim,C4S5Lim,C4DSS5Lim,&
         C5Lim,C5DSLim,C5S5Lim,C5DSS5Lim,C4C5Lim,C4C5DSLim,C4C5S5Lim,C4C5DSS5Lim)  
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_s_jj(HardProc%AmpMom,res_jj)
       respdf = s_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    if (vn_tagc+vn_tagce .ne. 0) then    

    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5
    !-- 3 -> S5C5

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(2) = zero
    FintNNLO_wp_s(3) = zero
#endif       

    kin(2) = zero

    !-- 4 -> C4

    xi1 = C4Lim%PartFrac(1) 
    xi2 = C4Lim%PartFrac(2) 
    call cut_histo_w(xi1,xi2,5,C4Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C4Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(4) = zero
#endif       

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C4Lim%AmpMom,res_sping_j,res_j)
       respdf_sping(:,:) = s_lumi_gq_sping(res_sping_j,f1,f2)

       s14 = C4Lim%KinInv(1) 
       z1  = C4Lim%KinInv(4) 
       !
       nperp=C4Lim%LimMom(2:4,1)
       ncoll=C4Lim%LimMom(2:4,2)

       hell_mat = comp_hell_mat_pgq_flip(z1,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
       call zipp(respdf_sping,hell_mat,respdf)

       kin(3) = (-one)*(-two/s14)*(respdf*tag_s_s &
            ) * C4Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(4) = kin(3)
#endif

       call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber)

    endif

    !-- 5 -> C5

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2) 
    call cut_histo_w(xi1,xi2,5,C5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(5) = zero
#endif       

       kin(4) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C5Lim%AmpMom,res_sping_j,res_j)
       respdf_sping(:,:) = s_lumi_qg_sping(res_sping_j,f1,f2)

       s25 = C5Lim%KinInv(2) 
       z2  = C5Lim%KinInv(5)
       !
       nperp=C5Lim%LimMom(2:4,1)
       ncoll=C5Lim%LimMom(2:4,2)

       hell_mat = comp_hell_mat_pgq_flip(z2,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
       call zipp(respdf_sping,hell_mat,respdf)

       kin(4) = (-one)*(-two/s25)*(respdf*tag_s_s &
            ) * C5Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(5) = kin(4)
#endif

    endif

    else
       kin(2:4) = zero
    endif

    call HistoFill(10,binarray,weightarray,kin(4),vegasweight,iternumber)

    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!
    
    !-- 6 -> DS
    !-- 7 -> DSS5
    !-- 8 -> C4DS
    !-- 9 -> C4DSS5
    !-- 10 -> DSC5
    !-- 11 -> C5DS
    !-- 12 -> C4C5DS
    !-- 13 -> C4C5DSS5

    kin(5) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(6:13)  = zero
#endif

    !-- 14 -> C4S5
    !-- 15 -> C4C5S5

    kin(6) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(14:15)  = zero
#endif

    !-- 16 -> C4C5

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(16) = zero
#endif       

    kin(7) = zero

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,16
          print *, i, FintNNLO_wp_s(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_ns = FintNNLO_wp_s
#endif

!     do i = 1,16
!        print *, i,FintNNLO_wp_s(i)
!     enddo
!     print *, 'ff(1)', ff(1)
!     pause

    return

  end function xsect_nnlo_4152_rr_wp_s


  !------------------------------------------------------------------------------------


  function xsect_nnlo_4252a_rr_wp_s(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4252a_rr_wp_s,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C5Lim,C5DSLim,C5DSS5Lim,C5S5Lim,TCC5Lim,TCC5DSLim,TCC5DSS5Lim,TCC5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z4,z5,s24,s25,s45,s4_25,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_s(16),kin(1:6)
    logical  :: makecut

    callcount = callcount + 1

    xsect_nnlo_4252a_rr_wp_s = 0
    ff(1) = zero
    xx = zero
    
    HardProc=>MyKinSet%AllKinConfig(1)
    DSLim=>MyKinSet%AllKinConfig(2)
    S5Lim=>MyKinSet%AllKinConfig(3)
     DSS5Lim=>MyKinSet%AllKinConfig(4)
    TCLim=>MyKinSet%AllKinConfig(5)
    TCDSLim=>MyKinSet%AllKinConfig(6)
    TCS5Lim=>MyKinSet%AllKinConfig(7)
    TCDSS5Lim=>MyKinSet%AllKinConfig(8)
    C5Lim=>MyKinSet%AllKinConfig(9)
    C5DSLim=>MyKinSet%AllKinConfig(10)
    C5S5Lim=>MyKinSet%AllKinConfig(11)
    C5DSS5Lim=>MyKinSet%AllKinConfig(12)
    TCC5Lim=>MyKinSet%AllKinConfig(13)
    TCC5DSLim=>MyKinSet%AllKinConfig(14)
    TCC5S5Lim=>MyKinSet%AllKinConfig(15)
    TCC5DSS5Lim=>MyKinSet%AllKinConfig(16)

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4252a(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
         C5Lim,C5DSLim,C5S5Lim,C5DSS5Lim,TCC5Lim,TCC5DSLim,TCC5S5Lim,TCC5DSS5Lim)
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_s_jj(HardProc%AmpMom,res_jj)
       respdf = s_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    if (vn_tagc+vn_tagce .ne. 0) then    
    
    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5
    !-- 3 -> S5C5

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(2) = zero
    FintNNLO_wp_s(3) = zero
#endif       

    kin(2) = zero

    !! -----------------------------!!
    !!     calling NLO ME  w sp corr!!
    !! -----------------------------!!
    
    !-- 4 -> C5

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(4)  = zero
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C5Lim%AmpMom,res_sping_j,res_j)
       respdf_sping = s_lumi_qg_sping(res_sping_j,f1,f2)

       s25 = C5Lim%KinInv(2) 
       z2  = C5Lim%KinInv(5)
       
       nperp=C5Lim%LimMom(2:4,1)
       ncoll=C5Lim%LimMom(2:4,2)

       hell_mat = comp_hell_mat_pgq_flip(z2,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
       call zipp(respdf_sping,hell_mat,respdf)

       kin(3) = (-one)*(-two/s25)*(respdf*tag_s_s) * C5Lim%wgt*jac*pref


#if (_DebugFeatures == 1)
       FintNNLO_wp_s(4)  = kin(3)
#endif

      call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber)
    endif


    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!
    
    !-- 5 -> DS
    !-- 6 -> DSS5
    !-- 7 -> DSTC
    !-- 8 -> DSTCS5
    !-- 9 -> DSC5
    !-- 10 -> DSS5C5
    !-- 11 -> DSTCC5
    !-- 12 -> DSTCC5S5

    kin(4) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(5:12)  = zero
#endif

    !-- 13 -> TC
    !-- 14 -> TCC5

    xi1 = TCLim%PartFrac(1) 
    xi2 = TCLim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCLim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCLim%flag) then

       kin(5) = zero

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(13:14)  = zero
#endif

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       !-- no spin correlation for Hgg amplitude
       call restree(TCLim%AmpMom,res)
       respdf = s_lumi_ix(res,f1,f2)
       respdf_b = s_lumi_ix_id(res,f1,f2)

       s24 = TCLim%KinInv(1)
       s25 = TCLim%KinInv(2)
       s45 = TCLim%KinInv(3)
       !
       z4 = TCLim%KinInv(4)
       z5 = TCLim%KinInv(5)
       z2 = one-z4-z5

       FintNNLO_wp_s(13) = (-one)*( &
            (PqbqR(-s24,s45,-s25,z4,z2,z5)+PqbqR(-s25,s45,-s24,z5,z2,z4))*respdf*tag_s_s+ & !-- 4 <-> 5
             Pqbqq_id_symm(-s24,-s25,s45,z2,z4,z5)*respdf_b*tag_s_v &
            )*TCLim%wgt*jac*pref 

       s25 = TCC5Lim%KinInv(1)
       s4_25 = TCC5Lim%KinInv(2)
       z1 = TCC5Lim%KinInv(3)
       z2 = TCC5Lim%KinInv(4)
       kt1kt2sq = TCC5Lim%KinInv(5)

       FintNNLO_wp_s(14) = Ptcso_qbq_r(s25,s4_25,z1,z2,kt1kt2sq)*tag_s_s*respdf * TCC5Lim%wgt*jac*pref 

       kin(5) = sum(FintNNLO_wp_s(13:14))

    endif

    else
       kin(2:5) = zero
    endif

    call HistoFill(10,binarray,weightarray,kin(5),vegasweight,iternumber)

    !-- 15 -> TCS5
    !-- 16 -> TCC5S5

    kin(6) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(15:16)  = zero
#endif

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,16
          print *, i, FintNNLO_wp_s(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_ns = FintNNLO_wp_s
#endif

    ! do i = 1,16
    !    print *, i,FintNNLO_wp_s(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return

  end function xsect_nnlo_4252a_rr_wp_s

  function xsect_nnlo_4252b_rr_wp_s(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4252b_rr_wp_s,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
          C54Lim,C54DSLim,C54DSS5Lim,C54S5Lim,TCC54Lim,TCC54DSLim,TCC54DSS5Lim,TCC54S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z3,z4,z5,z45,s12,s14,s24,s25,s45,s2_45,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_s(16),kin(1:6)
    logical :: makecut

    callcount = callcount + 1

    xsect_nnlo_4252b_rr_wp_s = 0
    ff(1) = zero
    xx = zero 

    HardProc=>MyKinSet%AllKinConfig(1)
    DSLim=>MyKinSet%AllKinConfig(2)
    S5Lim=>MyKinSet%AllKinConfig(3)
    DSS5Lim=>MyKinSet%AllKinConfig(4)
    TCLim=>MyKinSet%AllKinConfig(5)
    TCDSLim=>MyKinSet%AllKinConfig(6)
    TCS5Lim=>MyKinSet%AllKinConfig(7)
    TCDSS5Lim=>MyKinSet%AllKinConfig(8)
    C54Lim=>MyKinSet%AllKinConfig(9)
    C54DSLim=>MyKinSet%AllKinConfig(10)
    C54S5Lim=>MyKinSet%AllKinConfig(11)
    C54DSS5Lim=>MyKinSet%AllKinConfig(12)
    TCC54Lim=>MyKinSet%AllKinConfig(13)
    TCC54DSLim=>MyKinSet%AllKinConfig(14)
    TCC54S5Lim=>MyKinSet%AllKinConfig(15)
    TCC54DSS5Lim=>MyKinSet%AllKinConfig(16)

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*(xx(xRHO2)**buffrr_pow).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4252b(M2,ylab,xx(kNNLO_min:kNNLO_max_full),&
         HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
         C54Lim,C54DSLim,C54S5Lim,C54DSS5Lim,TCC54Lim,TCC54DSLim,TCC54S5Lim,TCC54DSS5Lim) 

    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_s_jj(HardProc%AmpMom,res_jj)
       respdf = s_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    if (vn_tagc+vn_tagce .ne. 0) then    

    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5
    !-- 3 -> S5C5

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(2) = zero
    FintNNLO_wp_s(3) = zero
#endif       

    kin(2) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(4)  = zero
#endif

    kin(3) = zero


    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!
    
    !-- 5 -> DS
    !-- 6 -> DSS5
    !-- 7 -> DSTC
    !-- 8 -> DSTCS5
    !-- 9 -> DSC5
    !-- 10 -> DSS5C5
    !-- 11 -> DSTCC5
    !-- 12 -> DSTCC5S5

    kin(4) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(5:12)  = zero
#endif

    !-- 13 -> TC
    !-- 14 -> TCC5

    xi1 = TCLim%PartFrac(1) 
    xi2 = TCLim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCLim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCLim%flag) then

       kin(5) = zero

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(13:14)  = zero
#endif

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCLim%AmpMom,res)
       respdf = s_lumi_ix(res,f1,f2)
       respdf_b = s_lumi_ix_id(res,f1,f2)

       s24 = TCLim%KinInv(1)
       s25 = TCLim%KinInv(2)
       s45 = TCLim%KinInv(3)
       !
       z4 = TCLim%KinInv(4)
       z5 = TCLim%KinInv(5)
       z2 = one-z4-z5

       FintNNLO_wp_s(13) = (-one)*( &
            (PqbqR(-s24,s45,-s25,z4,z2,z5)+PqbqR(-s25,s45,-s24,z5,z2,z4))*respdf*tag_s_s+ & !-- 4 <-> 5
             Pqbqq_id_symm(-s24,-s25,s45,z2,z4,z5)*respdf_b*tag_s_v &
            )*TCLim%wgt*jac*pref 

       FintNNLO_wp_s(14) = zero

       kin(5) = sum(FintNNLO_wp_s(13:14))

    endif

    else
       kin(2:5) = zero
    endif

    call HistoFill(10,binarray,weightarray,kin(5),vegasweight,iternumber)

    !-- 15 -> TCS5
    !-- 16 -> TCC5S5

    kin(6) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(15:16)  = zero
#endif

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,16
          print *, i, FintNNLO_wp_s(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_ns = FintNNLO_wp_s
#endif

    ! do i = 1,16
    !    print *, i,FintNNLO_wp_s(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return

  end function xsect_nnlo_4252b_rr_wp_s

  function xsect_nnlo_4252c_rr_wp_s(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4252c_rr_wp_s,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C4Lim,C4DSLim,C4DSS5Lim,C4S5Lim,TCC4Lim,TCC4DSLim,TCC4DSS5Lim,TCC4S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z4,z5,s12,s24,s25,s15,s45,s5_24,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_s(16),kin(1:6)
    logical  :: makecut

    callcount = callcount + 1

    xsect_nnlo_4252c_rr_wp_s = 0
    ff(1) = zero
    xx = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    DSLim=>MyKinSet%AllKinConfig(2)
    S5Lim=>MyKinSet%AllKinConfig(3)
    DSS5Lim=>MyKinSet%AllKinConfig(4)
    TCLim=>MyKinSet%AllKinConfig(5)
    TCDSLim=>MyKinSet%AllKinConfig(6)
    TCS5Lim=>MyKinSet%AllKinConfig(7)
    TCDSS5Lim=>MyKinSet%AllKinConfig(8)
    C4Lim=>MyKinSet%AllKinConfig(9)
    C4DSLim=>MyKinSet%AllKinConfig(10)
    C4S5Lim=>MyKinSet%AllKinConfig(11)
    C4DSS5Lim=>MyKinSet%AllKinConfig(12)
    TCC4Lim=>MyKinSet%AllKinConfig(13)
    TCC4DSLim=>MyKinSet%AllKinConfig(14)
    TCC4S5Lim=>MyKinSet%AllKinConfig(15)
    TCC4DSS5Lim=>MyKinSet%AllKinConfig(16)

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buffrr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4252c(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
         C4Lim,C4DSLim,C4S5Lim,C4DSS5Lim,TCC4Lim,TCC4DSLim,TCC4S5Lim,TCC4DSS5Lim)
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_s_jj(HardProc%AmpMom,res_jj)
       respdf = s_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    if (vn_tagc+vn_tagce .ne. 0) then    

    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(2) = zero
#endif       

    kin(2) = zero

    !! -----------------------------!!
    !!     calling NLO ME  w sp corr!!
    !! -----------------------------!!
    
    !-- 4 -> C4

    xi1 = C4Lim%PartFrac(1) 
    xi2 = C4Lim%PartFrac(2)
    call cut_histo_w(xi1,xi2,5,C4Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C4Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(4)  = zero
#endif

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C4Lim%AmpMom,res_sping_j,res_j)
       respdf_sping = s_lumi_qg_sping(res_sping_j,f1,f2)

       s24 = C4Lim%KinInv(2) 
       z2  = C4Lim%KinInv(5)
       
       nperp=C4Lim%LimMom(2:4,1)
       ncoll=C4Lim%LimMom(2:4,2)

       hell_mat = comp_hell_mat_pgq_flip(z2,nperp,ncoll)
       call zipp(respdf_sping,hell_mat,respdf)

       kin(3) = (-one)*(-two/s24)*(respdf*tag_s_s &
            ) * C4Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(4)  = kin(3)
#endif

      call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber)
    endif



    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!
    
    !-- 5 -> DS
    !-- 6 -> DSS5
    !-- 7 -> DSTC
    !-- 8 -> DSTCS5
    !-- 9 -> DSC4
    !-- 10 -> DSS5C4
    !-- 11 -> DSTCC4
    !-- 12 -> DSTCC4S5

    kin(4) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(5:12)  = zero
#endif

    !-- 13 -> TC
    !-- 14 -> TCC4

    xi1 = TCLim%PartFrac(1) 
    xi2 = TCLim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCLim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCLim%flag) then

       kin(5) = zero

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(13:14)  = zero
#endif

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCLim%AmpMom,res)
       respdf = s_lumi_ix(res,f1,f2)
       respdf_b = s_lumi_ix_id(res,f1,f2)

       s24 = TCLim%KinInv(1)
       s25 = TCLim%KinInv(2)
       s45 = TCLim%KinInv(3)
       !
       z4 = TCLim%KinInv(4)
       z5 = TCLim%KinInv(5)
       z2 = one-z4-z5

       FintNNLO_wp_s(13) = (-one)*( &
            (PqbqR(-s24,s45,-s25,z4,z2,z5)+PqbqR(-s25,s45,-s24,z5,z2,z4))*respdf*tag_s_s+ & !-- 4 <-> 5
             Pqbqq_id_symm(-s24,-s25,s45,z2,z4,z5)*respdf_b*tag_s_v &
            )*TCLim%wgt*jac*pref 

       s24 = TCC4Lim%KinInv(1) 
       s5_24 = TCC4Lim%KinInv(2) 
       z1  = TCC4Lim%KinInv(3) 
       z2 = TCC4Lim%KinInv(4)
       kt1kt2sq = TCC4Lim%KinInv(5)

       FintNNLO_wp_s(14) = Ptcso_qbq_r(s24,s5_24,z1,z2,kt1kt2sq)*tag_s_s*respdf * TCC4Lim%wgt*jac*pref 
       
       kin(5) = sum(FintNNLO_wp_s(13:14))

    endif

    else
       kin(2:5) = zero
    endif

    call HistoFill(10,binarray,weightarray,kin(5),vegasweight,iternumber)

    !-- 3 --> C4S5
    !-- 15 -> TCS5
    !-- 16 -> TCC4S5

    kin(6) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(3)  = zero
    FintNNLO_wp_s(15:16)  = zero
#endif

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,16
          print *, i, FintNNLO_wp_s(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_ns = FintNNLO_wp_s
#endif

    ! do i = 1,16
    !    print *, i,FintNNLO_wp_s(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return
    
  end function xsect_nnlo_4252c_rr_wp_s

  function xsect_nnlo_4252d_rr_wp_s(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4252d_rr_wp_s,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
          C54Lim,C54DSLim,C54DSS5Lim,C54S5Lim,TCC54Lim,TCC54DSLim,TCC54DSS5Lim,TCC54S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z3,z4,z5,z45,s12,s24,s14,s25,s45,s2_45,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_s(16),kin(1:7)
    logical :: makecut

    callcount = callcount + 1

    xsect_nnlo_4252d_rr_wp_s = 0
    ff(1) = zero
    xx = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    DSLim=>MyKinSet%AllKinConfig(2)
    S5Lim=>MyKinSet%AllKinConfig(3)
    DSS5Lim=>MyKinSet%AllKinConfig(4)
    TCLim=>MyKinSet%AllKinConfig(5)
    TCDSLim=>MyKinSet%AllKinConfig(6)
    TCS5Lim=>MyKinSet%AllKinConfig(7)
    TCDSS5Lim=>MyKinSet%AllKinConfig(8)
    C54Lim=>MyKinSet%AllKinConfig(9)
    C54DSLim=>MyKinSet%AllKinConfig(10)
    C54S5Lim=>MyKinSet%AllKinConfig(11)
    C54DSS5Lim=>MyKinSet%AllKinConfig(12)
    TCC54Lim=>MyKinSet%AllKinConfig(13)
    TCC54DSLim=>MyKinSet%AllKinConfig(14)
    TCC54S5Lim=>MyKinSet%AllKinConfig(15)
    TCC54DSS5Lim=>MyKinSet%AllKinConfig(16)

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*(xx(xRHO2)**buffrr_pow).lt.buffrr) return 

    call get_tauy(xx(1:2),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4252d(M2,ylab,xx(kNNLO_min:kNNLO_max_full),&
         HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
         C54Lim,C54DSLim,C54S5Lim,C54DSS5Lim,TCC54Lim,TCC54DSLim,TCC54S5Lim,TCC54DSS5Lim) 

    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_s_jj(HardProc%AmpMom,res_jj)
       respdf = s_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    if (vn_tagc+vn_tagce .ne. 0) then    

    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(2) = zero
#endif       
    
    kin(2) = zero

    !-- 3 -> S5C54

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(3) = zero
#endif       

    kin(3) = zero


    !! -----------------------------!!
    !!     calling NLO ME  w sp corr!!
    !! -----------------------------!!
    
    !-- 4 -> C54

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(4)  = zero
#endif

    kin(4) = zero

    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!
    
    !-- 5 -> DS
    !-- 6 -> DSS5
    !-- 7 -> DSTC
    !-- 8 -> DSTCS5
    !-- 9 -> DSC5
    !-- 10 -> DSS5C5
    !-- 11 -> DSTCC5
    !-- 12 -> DSTCC5S5

    kin(5) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(5:12)  = zero
#endif

    !-- 13 -> TC
    !-- 14 -> TCC5

    xi1 = TCLim%PartFrac(1) 
    xi2 = TCLim%PartFrac(2)
    call cut_histo_w(xi1,xi2,4,TCLim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.TCLim%flag) then

       kin(6) = zero

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(13:14)  = zero
#endif

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree(TCLim%AmpMom,res)
       respdf = s_lumi_ix(res,f1,f2)
       respdf_b = s_lumi_ix_id(res,f1,f2)

       s24 = TCLim%KinInv(1)
       s25 = TCLim%KinInv(2)
       s45 = TCLim%KinInv(3)
       !
       z4 = TCLim%KinInv(4)
       z5 = TCLim%KinInv(5)
       z2 = one-z4-z5

       FintNNLO_wp_s(13) = (-one)*( &
            (PqbqR(-s24,s45,-s25,z4,z2,z5)+PqbqR(-s25,s45,-s24,z5,z2,z4))*respdf*tag_s_s+ & !-- 4 <-> 5
             Pqbqq_id_symm(-s24,-s25,s45,z2,z4,z5)*respdf_b*tag_s_v &
            )*TCLim%wgt*jac*pref 

       FintNNLO_wp_s(14) = zero

       kin(6) = sum(FintNNLO_wp_s(13:14))

    endif

    else
       kin(2:6) = zero
    endif

    call HistoFill(10,binarray,weightarray,kin(6),vegasweight,iternumber)

    !-- 15 -> TCS5
    !-- 16 -> TCC5S5

    kin(7) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(15:16)  = zero
#endif

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,16
          print *, i, FintNNLO_wp_s(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_ns = FintNNLO_wp_s
#endif

    ! do i = 1,16
    !    print *, i,FintNNLO_wp_s(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return

  end function xsect_nnlo_4252d_rr_wp_s

  function xsect_nnlo_4251_rr_wp_s(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_4251_rr_wp_s,ndim,ncomp,userdata,iternumber,binarray(nobs),i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig),pointer :: HardProc,DSLim,DSS5Lim,S5Lim,C4Lim,C4DSLim,C4DSS5Lim,C4S5Lim,&
         C5Lim,C5DSLim,C5DSS5Lim,C5S5Lim,C4C5Lim,C4C5DSLim,C4C5DSS5Lim,C4C5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR,weightarray(nobs)
    real(dp)    :: xi1,xi2,f1(-6:6),f2(-6:6),my_lumi,respdf,respdf_b
    complex(dp) :: respdf_sping(-1:1,-1:1),respdf_sping2(-1:1,-1:1,-1:1,-1:1)
    !--
    real(dp)    :: res(namp_lo),res_j(namp_nlo),res_jj(namp_nnlo)
    complex(dp) :: res_sping2(namp_lo,-1:1,-1:1,-1:1,-1:1),res_sping_j(namp_nlo,-1:1,-1:1),hell_mat(-1:1,-1:1)
    complex(dp) :: hell_mat42(-1:1,-1:1),hell_mat51(-1:1,-1:1)
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s24,s25,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_wp_s(16),kin(1:7)
    logical  :: makecut

    callcount = callcount + 1

    xsect_nnlo_4251_rr_wp_s = 0
    ff(1) = zero
    xx = zero

    HardProc=>MyKinSet%AllKinConfig(1)
    DSLim=>MyKinSet%AllKinConfig(2)
    S5Lim=>MyKinSet%AllKinConfig(3)
    DSS5Lim=>MyKinSet%AllKinConfig(4)
    C4Lim=>MyKinSet%AllKinConfig(5)
    C4DSLim=>MyKinSet%AllKinConfig(6)
    C4S5Lim=>MyKinSet%AllKinConfig(7)
    C4DSS5Lim=>MyKinSet%AllKinConfig(8)
    C5Lim=>MyKinSet%AllKinConfig(9)
    C5DSLim=>MyKinSet%AllKinConfig(10)
    C5S5Lim=>MyKinSet%AllKinConfig(11)
    C5DSS5Lim=>MyKinSet%AllKinConfig(12)
    C4C5Lim=>MyKinSet%AllKinConfig(13)
    C4C5DSLim=>MyKinSet%AllKinConfig(14)
    C4C5S5Lim=>MyKinSet%AllKinConfig(15)
    C4C5DSS5Lim=>MyKinSet%AllKinConfig(16)

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full))

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buffrr) return 
    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_4251(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,DSLim,S5Lim,DSS5Lim,C4Lim,C4DSLim,C4S5Lim,C4DSS5Lim,&
         C5Lim,C5DSLim,C5S5Lim,C5DSS5Lim,C4C5Lim,C4C5DSLim,C4C5S5Lim,C4C5DSS5Lim)  
    
    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!
    
    !1-- regular piece
    
    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo_w(xi1,xi2,6,HardProc%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.HardProc%flag) then
       
#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = zero       
#endif

       kin(1) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_s_jj(HardProc%AmpMom,res_jj)
       respdf = s_lumi_nnlo(res_jj,f1,f2)

       kin(1) = respdf*HardProc%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(1) = kin(1)
#endif

    endif
    call HistoFill(0,binarray,weightarray,kin(1),vegasweight,iternumber) !-- initialize histograms, no matter whether we pass cut or not

    if (vn_tagc+vn_tagce .ne. 0) then    

    !! -----------------------------!!
    !!     calling NLO ME           !!
    !! -----------------------------!!

    !-- 2 -> S5
    !-- 3 -> S5C5

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(2) = zero
    FintNNLO_wp_s(3) = zero
#endif       

    kin(2) = zero

    !-- 4 -> C4

    xi1 = C4Lim%PartFrac(1) 
    xi2 = C4Lim%PartFrac(2) 
    call cut_histo_w(xi1,xi2,5,C4Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C4Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(4) = zero
#endif       

       kin(3) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C4Lim%AmpMom,res_sping_j,res_j)
       respdf_sping = s_lumi_qg_sping(res_sping_j,f1,f2)

       s24 = C4Lim%KinInv(1) 
       z2  = C4Lim%KinInv(4) 
       !
       nperp=C4Lim%LimMom(2:4,1)
       ncoll=C4Lim%LimMom(2:4,2)

       hell_mat = comp_hell_mat_pgq_flip(z2,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
       call zipp(respdf_sping,hell_mat,respdf)

       kin(3) = (-one)*(-two/s24)*(respdf*tag_s_s &
            ) * C4Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(4) = kin(3)
#endif

       call HistoFill(5,binarray,weightarray,kin(3),vegasweight,iternumber)

    endif

    !-- 5 -> C5

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2) 
    call cut_histo_w(xi1,xi2,5,C5Lim%AmpMom,makecut,binarray,weightarray,mu)

    if (makecut.or.C5Lim%flag) then

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(5) = zero
#endif       

       kin(4) = zero

    else

       call setscale(mu,muF,muR)
       call get_pdf(xi1,xi2,muF,f1,f2)
       call get_prefactor(muR,2,pref)

       call restree_j(C5Lim%AmpMom,res_sping_j,res_j)
       respdf_sping = s_lumi_gq_sping(res_sping_j,f1,f2)

       s15 = C5Lim%KinInv(2) 
       z1  = C5Lim%KinInv(5)
       !
       nperp=C5Lim%LimMom(2:4,1)
       ncoll=C5Lim%LimMom(2:4,2)

       hell_mat = comp_hell_mat_pgq_flip(z1,nperp,ncoll) ! n1: correlation vector, n2: collinear direction
       call zipp(respdf_sping,hell_mat,respdf)

       kin(4) = (-one)*(-two/s15)*(respdf*tag_s_s &
            ) * C5Lim%wgt*jac*pref

#if (_DebugFeatures == 1)
       FintNNLO_wp_s(5) = kin(4)
#endif

    endif

    else
       kin(2:4) = zero
    endif

    call HistoFill(10,binarray,weightarray,kin(4),vegasweight,iternumber)

    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!
    
    !-- 6 -> DS
    !-- 7 -> DSS5
    !-- 8 -> C4DS
    !-- 9 -> C4DSS5
    !-- 10 -> DSC5
    !-- 11 -> C5DS
    !-- 12 -> C4C5DS
    !-- 13 -> C4C5DSS5

    kin(5) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(6:13)  = zero
#endif

    !-- 14 -> C4S5
    !-- 15 -> C4C5S5

    kin(6) = zero

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(14:15)  = zero
#endif

    !-- 16 -> C4C5

#if (_DebugFeatures == 1)
    FintNNLO_wp_s(16) = zero
#endif       

    kin(7) = zero

    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, 'nan'
       print *, 'xx',xx
       print *, ''
       do i = 1,16
          print *, i, FintNNLO_wp_s(i)
       enddo
       ff(1) = zero
    endif

#if (_DebugFeatures == 1)
    FintNNLO_CheckLim_wp_ns = FintNNLO_wp_s
#endif

    ! do i = 1,16
    !    print *, i,FintNNLO_wp_s(i)
    ! enddo
    ! print *, 'ff(1)', ff(1)
    ! pause

    return

  end function xsect_nnlo_4251_rr_wp_s

end module mod_xsects_nnlo_rr_wp_s

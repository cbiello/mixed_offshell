module mod_xsects_nnlo_rr_z_ga
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_kinematics_gen
  use mod_kinematicsNNLO_EunordAord_ga
  use mod_kinematics_nnlo_is_eu_tc_ac, only: kinematics_nnlo_ac_5i6iac
  use mod_process
  use mod_auxfunctions
  use mod_aux_sectors
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_aux_limits
  use mod_limits_z_ga
  implicit none
  integer, save :: callcount
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO(16)
#endif

  private

  public :: xsect_nnlo_5161a_rr_z_ga,xsect_nnlo_5161c_rr_z_ga !-- checked
  public :: xsect_nnlo_5262a_rr_z_ga,xsect_nnlo_5262c_rr_z_ga

  public :: xsect_nnlo_5162_rr_z_ga,xsect_nnlo_5261_rr_z_ga
  public :: xsect_nnlo_5163_rr_z_ga,xsect_nnlo_5263_rr_z_ga
  public :: xsect_nnlo_5164_rr_z_ga,xsect_nnlo_5264_rr_z_ga

contains

  function xsect_nnlo_5161a_rr_z_ga(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_5161a_rr_z_ga,ndim,ncomp,userdata,iternumber
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,S5Lim,S5S6Lim,S6Lim,TCLim,TCS5Lim,TCS5S6Lim,TCS6Lim, &
         C6Lim,C6S5Lim,C6S5S6Lim,C6S6Lim,TCC6Lim,TCC6S6Lim,TCC6S5S6Lim,TCC6S5Lim
    type(KinConfig) :: HardProc_new,S5Lim_new,S5S6Lim_new,S6Lim_new,TCLim_new,TCS5Lim_new,TCS5S6Lim_new,TCS6Lim_new, &
         C6Lim_new,C6S5Lim_new,C6S5S6Lim_new,C6S6Lim_new,TCC6Lim_new,TCC6S6Lim_new,TCC6S5S6Lim_new,TCC6S5Lim_new, &
         C5Lim_new,C5S5Lim_new,C5S6Lim_new,TCC5Lim_new,TCC5S5Lim_new,TCC5S6Lim_new
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,ylab,jac,m2
    real(dp)    :: xi1,xi2
    real(dp)    :: FintNNLO_dy_ns(16),kin(1:9)
    real(dp)    :: respdf(ipdf),respdf_tmp(ipdf,4)
    real(dp)    :: resnnlo(2,2),res_nlo(2,2),res_a_nlo(2,2),res_g_nlo(2,2)
    real(dp)    :: res_tmp(2,2),res_lo(2,2)
    real(dp)    :: reseik_qcd,reseik_qed,eik(4)
    real(dp)    :: z1,z2,z5,z6,s15,s16,s56,s5_16,split(2,2),ker(2,2)

    xsect_nnlo_5161a_rr_z_ga = 0

    ff(1) = zero

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call open_histo()
    
    call kinematics_nnlo_5161a(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S5Lim,S6Lim,S5S6Lim,TCLim,TCS5Lim,TCS6Lim,TCS5S6Lim,&
         C6Lim,C6S5Lim,C6S6Lim,C6S5S6Lim,TCC6Lim,TCC6S5Lim,TCC6S6Lim,TCC6S5S6Lim)  

    ! call kinematics_nnlo_ac_5i6iac(xx,1,2,1, &
    !      HardProc=HardProc_new,S5Lim=S5Lim_new,S6Lim=S6Lim_new,S5S6Lim=S5S6Lim_new,TCLim=TCLim_new,TCS5Lim=TCS5Lim_new,&
    !      TCS6Lim=TCS6Lim_new,TCS5S6Lim=TCS5S6Lim_new,C6Lim=C6Lim_new,C6S5Lim=C6S5Lim_new,C6S6Lim=C6S6Lim_new,&
    !      C6S5S6Lim=C6S5S6Lim_new,TCC6Lim=TCC6Lim_new,TCC6S5Lim=TCC6S5Lim_new,TCC6S6Lim=TCC6S6Lim_new,TCC6S5S6Lim=TCC6S5S6Lim_new)  

    ! print *, 'wgt'
    ! print *, 'H       ',2*HardProc%wgt*jac,HardProc_new%wgt
    ! print *, 'S5      ',2*S5Lim%wgt*jac,S5Lim_new%wgt
    ! print *, 'S6      ',2*S6Lim%wgt*jac,S6Lim_new%wgt
    ! print *, 'S5S6    ',2*S5S6Lim%wgt*jac,S5S6Lim_new%wgt
    ! print *, 'TC      ',2*TCLim%wgt*jac,TCLim_new%wgt
    ! print *, 'TCS5    ',2*TCS5Lim%wgt*jac,TCS5Lim_new%wgt
    ! print *, 'TCS6    ',2*TCS6Lim%wgt*jac,TCS6Lim_new%wgt
    ! print *, 'TCS5S6  ',2*TCS5S6Lim%wgt*jac,TCS5S6Lim_new%wgt
    ! print *, 'C6      ',2*C6Lim%wgt*jac,C6Lim_new%wgt
    ! print *, 'C6S5    ',2*C6S5Lim%wgt*jac,C6S5Lim_new%wgt
    ! print *, 'C6S6    ',2*C6S6Lim%wgt*jac,C6S6Lim_new%wgt
    ! print *, 'C6S5S6  ',2*C6S5S6Lim%wgt*jac,C6S5S6Lim_new%wgt
    ! print *, 'TCC6    ',2*TCC6Lim%wgt*jac,TCC6Lim_new%wgt
    ! print *, 'TCC6S5  ',2*TCC6S5Lim%wgt*jac,TCC6S5Lim_new%wgt
    ! print *, 'TCC6S6  ',2*TCC6S6Lim%wgt*jac,TCC6S6Lim_new%wgt
    ! print *, 'TCC6S5S6',2*TCC6S5S6Lim%wgt*jac,TCC6S5S6Lim_new%wgt

    ! call mypause()

    !--------------------------------------
    !    NNLO ME :   gluon4,  photon 5
    !--------------------------------------

    ! 1-- regular piece ; quark-quark
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNNLO_dy_ns(1) = zero

    else    

       call res_tree_ga_qqb(HardProc%AmpMom,resnnlo)
       call get_respdf_ns(1,1,HardProc,resnnlo,respdf)

       respdf = respdf*HardProc%wgt*jac

       kin(1) = respdf(1)
       FintNNLO_dy_ns(1) = kin(1)

       call fill_histo([respdf],vegasweight)

    endif

    !! -----------------------------!!
    !!    SOFT   5, gluon           !!
    !! -----------------------------!!

    call cut_histo(S5Lim)

    if (S5Lim%makecut.or.S5Lim%flag) then

       kin(2) = zero
       FintNNLO_dy_ns(2) = zero

    else

       call res_tree_a_qqb(S5Lim%AmpMom,res_a_nlo)
       call get_respdf_ns(1,1,S5Lim,res_a_nlo,respdf)

       reseik_qcd = myEik_g(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3),Cf)

       respdf = (-one)*respdf*reseik_qcd*S5Lim%wgt*jac

       FintNNLO_dy_ns(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo([respdf],vegasweight)

    endif

    !! ---------------------------------------!!
    !!  SOFT 5 (gluon), COLLINEAR 6 (photon)  !!
    !! ---------------------------------------!!

    call cut_histo(C6S5Lim)

    if (C6S5Lim%makecut.or.C6S5Lim%flag) then

       kin(3) = zero
       FintNNLO_dy_ns(3) = zero

    else

       call res_tree_qqb(C6S5Lim%AmpMom,res_tmp)

       reseik_qcd = myEik_g(C6S5Lim%LimMom(:,1),C6S5Lim%LimMom(:,2),C6S5Lim%LimMom(:,3),Cf)

       s16 = C6S5Lim%Lim_KinInv(2) 
       z1  = C6S5Lim%Lim_KinInv(5)
       call AP_nlo_ph(z1,one,split) 
       call product_for_pdf_matrix(res_tmp,split,res_lo)

       call get_respdf_ns(1,1,C6S5Lim,res_lo,respdf)
       respdf = respdf*reseik_qcd*(-two/s16)*C6S5Lim%wgt*jac

       FintNNLO_dy_ns(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo([respdf],vegasweight)

    endif

    !! -------------------------------------!!
    !!     SOFT    6, photon     &   C6 S6  !!
    !! -------------------------------------!!

    call cut_histo(S6Lim)

    if (S6Lim%makecut.or.S6Lim%flag) then

       kin(4) = zero
       FintNNLO_dy_ns(4:5) = zero

    else

       !-- S6
       call res_tree_g_qqb(S6Lim%AmpMom,res_tmp)
       call get_Eik_ph_z_ga(S6Lim%LimMom(:,1),S6Lim%LimMom(:,2),S6Lim%LimMom(:,3),S6Lim%LimMom(:,4),S6Lim%LimMom(:,5),eik)
       call product_for_pdf(res_tmp,eik,res_g_nlo)

       call get_respdf_ns(1,1,S6Lim,res_g_nlo,respdf_tmp(:,1))

       respdf_tmp(:,1) = (-one)*respdf_tmp(:,1)*S6Lim%wgt*jac
       FintNNLO_dy_ns(4) = respdf_tmp(1,1)

       !-- C6S6
       s16 = C6S6Lim%Lim_KinInv(2)
       z6 = C6S6Lim%Lim_KinInv(5)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2) 

       call product_for_pdf_matrix(res_tmp,ker,res_g_nlo)
       call get_respdf_ns(1,1,S6Lim,res_g_nlo,respdf_tmp(:,2))

       respdf_tmp(:,2) = respdf_tmp(:,2)*four/z6/s16 * C6S6Lim%wgt*jac
       FintNNLO_dy_ns(5) = respdf_tmp(1,2)

       respdf(:) = sum(respdf_tmp(:,1:2),2)

       kin(4) = respdf(1)

       call fill_histo([respdf],vegasweight)

    endif

    !! -----------------------------!!
    !!     C6 -- photon             !!
    !! -----------------------------!!

    call cut_histo(C6Lim)

    if (C6Lim%makecut.or.C6Lim%flag) then

       kin(5) = zero
       FintNNLO_dy_ns(6)  = zero
    else

       call res_tree_g_qqb(C6Lim%AmpMom,res_tmp)

       ker(1,1) = Qdn**2 
       ker(1,2) = Qup**2 
       ker(2,1) = ker(1,1) 
       ker(2,2) = ker(1,2) 

       call product_for_pdf_matrix(res_tmp,ker,res_g_nlo)

       call get_respdf_ns(1,1,C6Lim,res_g_nlo,respdf)

       s16 = C6Lim%Lim_KinInv(2) 
       z1  = C6Lim%Lim_KinInv(5)

       respdf = (-one)*(-two/s16)*(respdf*Pqq(z1,one))*C6Lim%wgt*jac

       FintNNLO_dy_ns(6)  = respdf(1)
       kin(5) = respdf(1)

       call fill_histo([respdf],vegasweight)

    endif

    !-- here

    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!

    !-- 7 -> S5S6
    !-- 8 -> TCS5S6
    !-- 9 -> C6S5S6 
    !-- 10 -> TCC6S5S6

    call cut_histo(S5S6Lim)

    if (S5S6Lim%makecut.or.S5S6Lim%flag) then

       kin(6) = zero
       FintNNLO_dy_ns(7:10) = zero

    else

       call res_tree_qqb(S5S6Lim%AmpMom,res_tmp) 

       !-- S5S6
       call get_Eik_ds_z_ga(S5S6Lim%LimMom(:,1),S5S6Lim%LimMom(:,2),S5S6Lim%LimMom(:,3),&
            S5S6Lim%LimMom(:,4),S5S6Lim%LimMom(:,5),S5S6Lim%LimMom(:,6),eik) 
       call product_for_pdf(res_tmp,eik,res_lo)

       call get_respdf_ns(1,1,S5S6Lim,res_lo,respdf_tmp(:,1))
       respdf_tmp(:,1) = (-one)*respdf_tmp(:,1)*S5S6Lim%wgt*jac
       FintNNLO_dy_ns(7) = respdf_tmp(1,1)

       !-- TCS5S6
       s15 = TCS5S6Lim%Lim_KinInv(1)
       s16 = TCS5S6Lim%Lim_KinInv(2)
       s56 = TCS5S6Lim%Lim_KinInv(3)
       !
       z5 = TCS5S6Lim%Lim_KinInv(4)
       z6 = TCS5S6Lim%Lim_KinInv(5)

       call get_Pgaq_ds_z_ga(s56,s15,s16,z5,z6,split)
       call product_for_pdf_matrix(res_tmp,split,res_lo)

       call get_respdf_ns(1,1,S5S6Lim,res_lo,respdf_tmp(:,2))
       respdf_tmp(:,2) = (-one)*respdf_tmp(:,2)*TCS5S6Lim%wgt*jac
       FintNNLO_dy_ns(8) = respdf_tmp(1,2)

       !-- C6S5S6
       s16 = C6S5S6Lim%Lim_KinInv(2)
       z6  = C6S5S6Lim%Lim_KinInv(5)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2)

       call product_for_pdf_matrix(res_tmp,ker,res_lo)
       call get_respdf_ns(1,1,S5S6Lim,res_lo,respdf_tmp(:,3))
       respdf_tmp(:,4) = respdf_tmp(:,3)

       respdf_tmp(:,3) = (respdf_tmp(:,3)) &
            *myEik_g(C6S5S6Lim%LimMom(:,1),C6S5S6Lim%LimMom(:,2),C6S5S6Lim%LimMom(:,3),Cf)*four/z6/s16 &
            *C6S5S6Lim%wgt*jac

       FintNNLO_dy_ns(9) = respdf_tmp(1,3)

       !-- TCC6S5S6
       s15 = TCC6S5S6Lim%Lim_KinInv(1)
       s16 = TCC6S5S6Lim%Lim_KinInv(2)
       !
       z5 = TCC6S5S6Lim%Lim_KinInv(4)
       z6 = TCC6S5S6Lim%Lim_KinInv(5)

       respdf_tmp(:,4) = respdf_tmp(:,4) * &
            four*Cf/(-s15)/z5*four/(-s16)/z6 * TCC6S5S6Lim%wgt*jac

       FintNNLO_dy_ns(10) = respdf_tmp(1,4)

       respdf(:) = sum(respdf_tmp(:,1:4),2)
       kin(6) = respdf(1)

       call fill_histo([respdf],vegasweight)

    endif

    !-- 11 -> TC
    !-- 12 -> TCC6

    call cut_histo(TCLim)

    if (TCLim%makecut.or.TCLim%flag) then

       kin(7) = zero
       FintNNLO_dy_ns(11:12) = zero 

    else

       call res_tree_qqb(TCLim%AmpMom,res_tmp)

       !-- TC
       s15 = TCLim%Lim_KinInv(1)
       s16 = TCLim%Lim_KinInv(2)
       s56 = TCLim%Lim_KinInv(3)
       !
       z5 = TCLim%Lim_KinInv(4)
       z6 = TCLim%Lim_KinInv(5)
       z1 = one-z5-z6

       call AP_nnlo_ga(s56,s15,s16,z5,z6,z1,split)
       call product_for_pdf_matrix(res_tmp,split,res_lo)

       call get_respdf_ns(1,1,TCLim,res_lo,respdf_tmp(:,1))
       respdf_tmp(:,1) = (-one)*respdf_tmp(:,1)*TCLim%wgt*jac
       FintNNLO_dy_ns(11) = respdf_tmp(1,1)

       !-- TCC6
       s16 = TCC6Lim%Lim_KinInv(1)
       s5_16 = TCC6Lim%Lim_KinInv(2)
       z1 = TCC6Lim%Lim_KinInv(3)
       z2 = TCC6Lim%Lim_KinInv(4)

       !print *, 'removed because of esq, careful'
       !call get_Pgaq_C_z_ga(s16,s5_16,z1,z2,esq,split)
       call get_Pgaq_C_z_ga(s16,s5_16,z1,z2,one,split)
       call product_for_pdf_matrix(res_tmp,split,res_lo)

       call get_respdf_ns(1,1,TCLim,res_lo,respdf_tmp(:,2))
       respdf_tmp(:,2) = respdf_tmp(:,2)*TCC6Lim%wgt*jac
       FintNNLO_dy_ns(12) = respdf_tmp(1,2)

       respdf(:) = sum(respdf_tmp(:,1:2),2)
       kin(7) = respdf(1)

       call fill_histo([respdf],vegasweight)

    endif

    !-- 13 -> TCS6
    !-- 14 -> TCC6S6

    call cut_histo(TCS6Lim)

    if (TCS6Lim%makecut.or.TCS6Lim%flag) then

       kin(8) = zero
       FintNNLO_dy_ns(13:14)  = zero

    else

       call res_tree_qqb(TCS6Lim%AmpMom,res_tmp)

       !-- TCS6
       s15 = TCS6Lim%Lim_KinInv(1)
       s16 = TCS6Lim%Lim_KinInv(2)
       z5 = TCS6Lim%Lim_KinInv(4)
       z6 = TCS6Lim%Lim_KinInv(5)

       call get_Pgaq_S_z_ga(TCS6Lim%Lim_KinInv,split)
       call product_for_pdf_matrix(res_tmp,split,res_lo)

       call get_respdf_ns(1,1,TCS6Lim,res_lo,respdf_tmp(:,1))
       respdf_tmp(:,1) = respdf_tmp(:,1)*TCS6Lim%wgt*jac
       FintNNLO_dy_ns(13) = respdf_tmp(1,1)



       !-- TCC6S6
       s15 = TCC6S6Lim%Lim_KinInv(1)
       s16 = TCC6S6Lim%Lim_KinInv(2)
       z5 =  TCC6S6Lim%Lim_KinInv(4)
       z6 =  TCC6S6Lim%Lim_KinInv(5)

       !print *, 'removed because of esq, careful'
       call get_Pgaq_SC_z_ga(s15,s16,z5,z6,one,split)          
       call product_for_pdf_matrix(res_tmp,split,res_lo)

       call get_respdf_ns(1,1,TCS6Lim,res_lo,respdf_tmp(:,2))
       respdf_tmp(:,2) = respdf_tmp(:,2) * TCC6S6Lim%wgt*jac
       FintNNLO_dy_ns(14) = respdf_tmp(1,2)

       respdf(:) = sum(respdf_tmp(:,1:2),2)
       kin(8) = respdf(1)

       call fill_histo([respdf],vegasweight)

    endif

    !-- 15 -> TCS5
    !-- 16 -> TCC6S5
    call cut_histo(TCS5Lim)

    if (TCS5Lim%makecut.or.TCS5Lim%flag) then

       kin(9) = zero
       FintNNLO_dy_ns(15:16)  = zero

    else

       call res_tree_qqb(TCS5Lim%AmpMom,res_tmp)

       ker(1,1) = Qdn**2 
       ker(1,2) = Qup**2 
       ker(2,1) = ker(1,1) 
       ker(2,2) = ker(1,2) 

       call product_for_pdf_matrix(res_tmp,ker,res_lo)
       call get_respdf_ns(1,1,TCS5Lim,res_lo,respdf_tmp(:,1))
       respdf_tmp(:,2) = respdf_tmp(:,1)

       !-- TCS5
       s15 = TCS5Lim%Lim_KinInv(1)
       s16 = TCS5Lim%Lim_KinInv(2)
       z5 = TCS5Lim%Lim_KinInv(4)
       z6 = TCS5Lim%Lim_KinInv(5)

       respdf_tmp(:,1) = (-one)*respdf_tmp(:,1)* four*CF/s15/z5 * two/s16*Pqq(z6,one)*TCS5Lim%wgt*jac
       FintNNLO_dy_ns(15) = respdf_tmp(1,1)

       !-- TCC6S5
       s15 = TCC6S5Lim%Lim_KinInv(1)
       s16 = TCC6S5Lim%Lim_KinInv(2)
       z5 =  TCC6S5Lim%Lim_KinInv(4)
       z6 =  TCC6S5Lim%Lim_KinInv(5)

       respdf_tmp(:,2) = (-four)*CF/s15/z5*two/s16*Pqq(z6,one)*respdf_tmp(:,2)*TCC6S5Lim%wgt*jac
       FintNNLO_dy_ns(16) = respdf_tmp(1,2)

       respdf(:) = sum(respdf_tmp(:,1:2),2)
       kin(9) = respdf(1)

       call fill_histo([respdf],vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_dy_ns)

#if(_withchecks == 1)
    FintNNLO = FintNNLO_dy_ns
#endif

  end function xsect_nnlo_5161a_rr_z_ga


  function xsect_nnlo_5161c_rr_z_ga(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_5161c_rr_z_ga,ndim,ncomp,userdata,iternumber
    real(dp15) :: yRnd(30),ff(1),vegasweight
    real(dp), parameter :: esq = one
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR
    real(dp)    :: xi1,xi2,f1(-6:7),f2(-6:7),my_lumi,respdf(1),respdf_b(1)
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--
    real(dp)    :: respdf_qed,respdf_qcd
    complex(dp) :: respdf_sping_qed(-1:1,-1:1),respdf_sping_qcd(-1:1,-1:1)
    real(dp)    :: respdf_ec,respdf_ce,respdf_b_qed(1),respdf_b_qcd
    real(dp)    :: reseik_qed,reseik_qcd
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s25,s45,s5_14,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_dy_ns(16),kin(1:9), FintNNLO_dy_ns2(1)
    logical  :: makecut
    !--
    real(dp) :: mysymm = half
    real(dp) :: qcdmom(4,3),qedmom(4,3)
    !
    real(dp) :: wmom(4)
    real(dp) :: respdf_ch1(1)
    !
    real(dp) :: resnnlo_2q(32)
    real(dp) :: resnnlo_4q_tc(4),resnnlo_4q_fin(8)
    real(dp) :: res_g_nlo(12)
    real(dp) :: vec_a(4),vec_g(4),vec_ag(4),vec_ga(4)
    real(dp) :: vec_ag_a(4),vec_ag_g(4),vec_ga_a(4),vec_ga_g(4)
    real(dp) :: resvec(4)
    real(dp) :: res_ag(4),res_ga(4)
    real(dp) :: respdf_ag(1),respdf_ga(1)

    real(dp) ::resnnlo(2,2), res_a_nlo(2,2), temp(2,2), res_S6_nlo(2,2), eik(4), res2(2,2), res_S6S5(2,2)
    real(dp) :: split(2,2), respdf2, respdf3, respdf4, respdf1, ker(4), res_TC(2,2), res_TC2(2,2)
    real(dp) :: res_C5S6(2,2), res_C5S62(2,2), res_TCS62(2,2), res_TCS6(2,2), res_TCS5(2,2), res_TCS52(2,2)

    type(KinConfig) :: HardProc,S5Lim,S6Lim,S5S6Lim,TCLim,TCS5Lim,TCS6Lim,TCS5S6Lim,C5Lim, &
         C5S5Lim,C5S6Lim,C5S5S6Lim,TCC5Lim,TCC5S5Lim,TCC5S6Lim,TCC5S5S6Lim
    integer :: i

    xsect_nnlo_5161c_rr_z_ga = 0
    ff(1) = zero
    kin = zero
    xx = zero

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)

    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buff_rr) return 

    call open_histo()

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call kinematics_nnlo_5161c(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S5Lim,S6Lim,S5S6Lim,TCLim,TCS5Lim,TCS6Lim,TCS5S6Lim,&
         C5Lim,C5S5Lim,C5S6Lim,C5S5S6Lim,TCC5Lim,TCC5S5Lim,TCC5S6Lim,TCC5S5S6Lim)

    !! ----------------------------------------!!
    !!    calling NNLO ME; gluon4 photon 5      !!
    !! ----------------------------------------!!

    !1-- regular piece 

    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo(HardProc)

    if (makecut.or.HardProc%flag) then
       kin(1) = zero

       FintNNLO_dy_ns(1) = zero

    else

       call res_tree_ga_qqb(HardProc%AmpMom,resnnlo)

       call get_respdf_ns(1,1,HardProc,resnnlo,respdf)

       kin(1) = respdf(1)*HardProc%wgt*jac

       FintNNLO_dy_ns(1) = kin(1)



       call fill_histo([kin(1)],vegasweight)

    endif


    !! -----------------------------!!
    !!     SOFT   4, gluon          !!
    !! -----------------------------!!
    !-- 2 -> S5       particle ''4'' is a gluon 

    xi1 = S5Lim%PartFrac(1) 
    xi2 = S5Lim%PartFrac(2) 

    call cut_histo(S5Lim)

    if (makecut.or.S5Lim%flag) then
       kin(2) = zero
       FintNNLO_dy_ns(2) = zero

    else

       call res_tree_a_qqb(S5Lim%AmpMom,res_a_nlo)
       call get_respdf_ns(1,1,S5Lim,res_a_nlo,respdf)

       reseik_qcd = myEik_g(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3),Cf)

       FintNNLO_dy_ns(2) = (-one)*respdf(1) *reseik_qcd*S5Lim%wgt*jac
       kin(2) = FintNNLO_dy_ns(2)


       call fill_histo([kin(2)],vegasweight)

    endif

    !-------- S5 C5

    xi1 = C5S5Lim%PartFrac(1) 
    xi2 = C5S5Lim%PartFrac(2) 

    call cut_histo(C5S5Lim)

    !    print*, 'makecut', makecut
    !    print*, 'C5S5Lim%flag', C5S5Lim%flag


    if (makecut.or.C5S5Lim%flag) then
       FintNNLO_dy_ns(3) = zero

       kin(3) = zero

       !       print*, 'I am inside the if'

    else

       call res_tree_a_qqb(C5S5Lim%AmpMom,res_a_nlo)

       call get_respdf_ns(1,1,C5S5Lim,res_a_nlo,respdf)

       s14 = C5S5Lim%Lim_KinInv(2)  
       z4 =  C5S5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(3) = respdf(1) *four*Cf/s14/z4*C5S5Lim%wgt*jac!correct


       !          print*, '***'
       !     print*, 'respdf,Cf,s14,z4,C5S5Lim%wgt,jac,pref', respdf,Cf,s14,z4,C5S5Lim%wgt,jac,pref
       !print*, '***'     
       kin(3) = FintNNLO_dy_ns(3)

       call fill_histo([kin(3)],vegasweight)

    endif


    !! -----------------------------!!
    !!    Soft 5,  photon            !!
    !! -----------------------------!!

    !-- 2 -> S6 ; 5 is a photon

    xi1 = S6Lim%PartFrac(1) 
    xi2 = S6Lim%PartFrac(2) 
    call cut_histo(S6Lim)

    if (makecut.or.S6Lim%flag) then

       kin(4) = zero
       FintNNLO_dy_ns(4) = zero

    else

       call res_tree_g_qqb(S6Lim%AmpMom,res_S6_nlo)
       call get_Eik_ph_z_ga(S6Lim%LimMom(:,1),S6Lim%LimMom(:,2),S6Lim%LimMom(:,3), &
            S6Lim%LimMom(:,4), S6Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S6_nlo,eik,res2)
       call get_respdf_ns(1,1,S6Lim,res2,respdf_ga)

       FintNNLO_dy_ns(4) = (-one)*respdf_ga(1)*S6Lim%wgt*jac

       kin(4)  = FintNNLO_dy_ns(4)

       call fill_histo([kin(4)],vegasweight)

    endif

    !! -----------------------------!!
    !!     C5 -- gluon !!
    !! -----------------------------!!

    !-- 4 -> C5 ; 

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)
    call cut_histo(C5Lim)

    if (makecut.or.C5Lim%flag) then
       kin(5) = zero
       FintNNLO_dy_ns(5) = zero  

    else

       call res_tree_a_qqb(C5Lim%AmpMom,res_a_nlo)
       call get_respdf_ns(1,1,C5Lim,res_a_nlo,respdf_b_qed)

       s14 = C5Lim%Lim_KinInv(2) 
       z1  = C5Lim%Lim_KinInv(5)


       FintNNLO_dy_ns(5)  = (-one)*(-two/s14)* Pqq(z1,Cf ) * respdf_b_qed(1) *C5Lim%wgt*jac


       kin(5) = FintNNLO_dy_ns(5)

       call fill_histo([kin(5)],vegasweight)

    endif


    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!


    !-- 6 -> S5S6
    !-- 7 -> TCS5S6
    !-- 8 -> C5S5S6
    !-- 9 -> TCC5S5S6

    xi1 = S5S6Lim%PartFrac(1) 
    xi2 = S5S6Lim%PartFrac(2)

    call cut_histo(S5S6Lim)

    if (makecut.or.S5S6Lim%flag) then
       kin(6) = zero
       FintNNLO_dy_ns(6:9)  = zero


    else

       call res_tree_qqb(S5S6Lim%AmpMom,res_S6S5) 

       call get_Eik_ds_z_ga(S5S6Lim%LimMom(:,1),S5S6Lim%LimMom(:,2),S5S6Lim%LimMom(:,3),&
            S5S6Lim%LimMom(:,4),S5S6Lim%LimMom(:,5),S5S6Lim%LimMom(:,6),eik)
       call product_for_pdf(res_S6S5,eik,res2)
       call get_respdf_ns(1,1,S5S6Lim,res2,respdf_ag)

       FintNNLO_dy_ns(6) = (-one)*(respdf_ag(1))*S5S6Lim%wgt*jac


       !------- TCS5S6

       s14 = TCS5S6Lim%Lim_KinInv(1)
       s15 = TCS5S6Lim%Lim_KinInv(2)
       s45 = TCS5S6Lim%Lim_KinInv(3)
       !
       z4 = TCS5S6Lim%Lim_KinInv(4)
       z5 = TCS5S6Lim%Lim_KinInv(5)

       call get_Pggq_ds_z_ga(s45,s14,s15,z4,z5,split)
       call product_for_pdf_matrix(res_S6S5,split,res2)
       call get_respdf_ns(1,1,S5S6Lim,res2,respdf)
       !respdf = ns_lumi_z_ga(res2,f1,f2) !correct

       !!print *,'TCS5S6Lim%wgt', TCS5S6Lim%wgt
       !!print *, 'xi1', xi1
       !!print *, 'xi2', xi2
       !!print *, 'respdf', respdf
       !!print *, 'jac', jac
       !!print *, 'pref', pref

       FintNNLO_dy_ns(7) = (-one)*respdf(1)*TCS5S6Lim%wgt*jac


       !----->  C5S5S6

       call get_Eik_ph_z_ga( C5S5S6Lim%LimMom(:,1), C5S5S6Lim%LimMom(:,2), C5S5S6Lim%LimMom(:,3), &
            C5S5S6Lim%LimMom(:,4),  C5S5S6Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S6S5,eik,res2)

       call get_respdf_ns(1,1,S5S6Lim,res2,respdf_ga)
       !respdf_ga = ns_lumi_z_ga(res2,f1,f2)

       FintNNLO_dy_ns(8) = respdf_ga(1)*4.0_dp*Cf*esq/z4/s14*C5S5S6Lim%wgt*jac


       !----->  TCC5S5S6

       s14 = TCC5S5S6Lim%Lim_KinInv(1) 
       s15 = TCC5S5S6Lim%Lim_KinInv(2) 

       z4 =  TCC5S5S6Lim%Lim_KinInv(4)
       z5 =  TCC5S5S6Lim%Lim_KinInv(5)

       ker(1) = Qdn**2
       ker(2) = Qup**2
       ker(3) = ker(1)
       ker(4) = ker(2)

       call product_for_pdf(res_S6S5,ker,res2)
       call get_respdf_ns(1,1,S5S6Lim,res2,respdf_ch1)
       !respdf_ch1 = ns_lumi_z_ga(res2,f1,f2)


       FintNNLO_dy_ns(9) = four*Cf/(-s14)/z4*four/(-s15)/z5 *respdf_ch1(1) * TCC5S5S6Lim%wgt*jac

       kin(6) = sum(FintNNLO_dy_ns(6:9))

       call fill_histo([kin(6)],vegasweight)

    endif

    !-- 10 -> TC
    !-- 11 -> TCC5

    xi1 = TCLim%PartFrac(1) 
    xi2 = TCLim%PartFrac(2)
    call cut_histo(TCLim)

    if (makecut.or.TCLim%flag) then

       kin(7) = zero
       FintNNLO_dy_ns(10:11)  = zero


    else

       s14 = TCLim%Lim_KinInv(1)
       s15 = TCLim%Lim_KinInv(2)
       s45 = TCLim%Lim_KinInv(3)
       !
       z4 = TCLim%Lim_KinInv(4)
       z5 = TCLim%Lim_KinInv(5)
       z1 = one-z4-z5

       call res_tree_qqb(TCLim%AmpMom,res_TC)
       call AP_nnlo_ga(s45,s14,s15,z4,z5,z1,split)
       call product_for_pdf_matrix(res_TC,split,res_TC2)
       call get_respdf_ns(1,1,TCLim,res_TC2,respdf)

       FintNNLO_dy_ns(10) = (-one)*respdf(1)*TCLim%wgt*jac


       s14 = TCC5Lim%Lim_KinInv(1) 
       s5_14 = TCC5Lim%Lim_KinInv(2) 
       z1  = TCC5Lim%Lim_KinInv(3) 
       z2 = TCC5Lim%Lim_KinInv(4)

       call get_Pggq_C_z_ga(s14,s5_14,z1,z2,esq,split)
       call product_for_pdf_matrix(res_TC,split,res_TC2)
       call get_respdf_ns(1,1,TCLim,res_TC2,respdf_ch1)
       !respdf_ch1 = ns_lumi_z_ga(res_TC2,f1,f2)

       FintNNLO_dy_ns(11) = respdf_ch1(1)*TCC5Lim%wgt*jac 

       kin(7) = sum(FintNNLO_dy_ns(10:11))

       call fill_histo([kin(7)],vegasweight)

    endif



    !-- 12 --> C5S6
    !-- 13 -> TCS6
    !-- 14 -> TCC5S6

    xi1 = TCS6Lim%PartFrac(1) 
    xi2 = TCS6Lim%PartFrac(2)
    call cut_histo(TCS6Lim)

    if (makecut.or.TCS6Lim%flag) then

       kin(8) = zero
       FintNNLO_dy_ns(12:14)  = zero

    else

       call res_tree_qqb(C5S6Lim%AmpMom,res_C5S6)
       call get_Eik_ph_z_ga(C5S6Lim%LimMom(:,1),C5S6Lim%LimMom(:,2),C5S6Lim%LimMom(:,3), &
            C5S6Lim%LimMom(:,4),C5S6Lim%LimMom(:,5),eik)
       call product_for_pdf(res_C5S6,eik,res_C5S62)
       call get_respdf_ns(1,1,TCS6Lim,res_C5S62,respdf_ga)

       s12 = C5S6Lim%Lim_KinInv(1) 
       s15 = C5S6Lim%Lim_KinInv(2) 
       s25 = C5S6Lim%Lim_KinInv(3) 
       s14 = C5S6Lim%Lim_KinInv(4) 
       z4 =  C5S6Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(12) = (-one)*(respdf_ga(1))*two/s14*Pqq(z4,Cf)* C5S6Lim%wgt*jac

       call res_tree_qqb(TCS6Lim%AmpMom,res_TCS6)

       s14 = TCS6Lim%Lim_KinInv(1)
       s15 = TCS6Lim%Lim_KinInv(2)
       z4 = TCS6Lim%Lim_KinInv(4)
       z5 = TCS6Lim%Lim_KinInv(5)

       call get_Pggq_S_z_ga(TCS6Lim%Lim_KinInv,split)
       call product_for_pdf_matrix(res_TCS6,split,res_TCS62)

       call get_respdf_ns(1,1,TCS6Lim,res_TCS62,respdf_ch1)
       !respdf_ch1 = ns_lumi_z_ga(res_TCS62,f1,f2) 
       FintNNLO_dy_ns(13) = respdf_ch1(1)*TCS6Lim%wgt*jac


       split(1,1) = Qdn**2
       split(1,2) = Qup**2
       split(2,1) = split(1,1)
       split(2,2) = split(1,2)


       call product_for_pdf_matrix(res_TCS6,split,res_TCS62)
       call get_respdf_ns(1,1,TCS6Lim,res_TCS62,respdf_ch1)
       !respdf_ch1 = ns_lumi_z_ga(res_TCS62,f1,f2) 

       s14 = TCC5S6Lim%Lim_KinInv(1) 
       s15 = TCC5S6Lim%Lim_KinInv(2) 
       s45 = TCC5S6Lim%Lim_KinInv(3) 
       z4 =  TCC5S6Lim%Lim_KinInv(4) 
       z5 =  TCC5S6Lim%Lim_KinInv(5) 
       !

       FintNNLO_dy_ns(14) = (-one)*8.0_dp/s14/s15/z5*Pqq(z4,Cf )/z4 *respdf_ch1(1) *TCC5S6Lim%wgt*jac


       kin(8) = sum(FintNNLO_dy_ns(12:14))

       call fill_histo([kin(8)],vegasweight)

    endif




    !-- 15 ->  TCS5
    !-- 16 -> TCC5S5

    xi1 = TCS5Lim%PartFrac(1) 
    xi2 = TCS5Lim%PartFrac(2)
    call cut_histo(TCS5Lim)

    if (makecut.or.TCS5Lim%flag) then

       kin(9) = zero
       FintNNLO_dy_ns(15:16)  = zero


    else

       split(1,1) = Qdn**2
       split(1,2) = Qup**2
       split(2,1) = split(1,1)
       split(2,2) = split(1,2)


       call res_tree_qqb(TCS5Lim%AmpMom,res_TCS5)

       ! print*, 'res_TCS5', res_TCS5(1,1)
       ! print*, 'res_TCS5', res_TCS5(1,2)
       ! print*, 'res_TCS5', res_TCS5(2,1)
       ! print*, 'res_TCS5', res_TCS5(2,2)


       call product_for_pdf_matrix(res_TCS5,split,res_TCS52)

       !print*, 'res_TCS5', res_TCS5(1,1)
       !print*, 'res_TCS5', res_TCS5(1,2)
       !print*, 'res_TCS5', res_TCS5(2,1)
       !print*, 'res_TCS5', res_TCS5(2,2)


       call get_respdf_ns(1,1,TCS5Lim,res_TCS52,respdf_ch1)

       s14 = TCS5Lim%Lim_KinInv(1) 
       s15 = TCS5Lim%Lim_KinInv(2) 
       s45 = TCS5Lim%Lim_KinInv(3) 
       z4 =  TCS5Lim%Lim_KinInv(4) 
       z5 =  TCS5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(15) = (-one)*respdf_ch1(1) *  four /s15 * Pqq(z5,one) * two*Cf/s14/z4 * TCS5Lim%wgt*jac

       !print*, 'respdf_ch1', respdf_ch1,  four /s15 * Pqq(z5,one) * two*Cf/s14/z4 * TCS5Lim%wgt*jac


       s14 = TCC5S5Lim%Lim_KinInv(1) 
       s15 = TCC5S5Lim%Lim_KinInv(2) 
       s45 = TCC5S5Lim%Lim_KinInv(3) 
       z4 =  TCC5S5Lim%Lim_KinInv(4) 
       z5 =  TCC5S5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(16) = respdf_ch1(1)* four /s15 * Pqq(z5,one) * two*Cf/s14/z4 *TCC5S5Lim%wgt*jac

       kin(9) = sum(FintNNLO_dy_ns(15:16))
       call fill_histo([kin(9)],vegasweight)

    endif

    print *, FintNNLO_dy_ns(1)
#if(_withchecks == 1)
    FintNNLO = FintNNLO_dy_ns
#endif

    ff(1) = sum(kin)
    call close_histo()

    if (ff(1) .ne. ff(1)) then
       print  *, 'nan'
       print  *, 'xx',xx
       print  *, ''
       do i = 1,16
          print  *, i, FintNNLO_dy_ns(i)
       enddo
       ff(1)=zero
    endif

    return

  end function xsect_nnlo_5161c_rr_z_ga

  !--

  function xsect_nnlo_5262a_rr_z_ga(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer ::  xsect_nnlo_5262a_rr_z_ga,ndim,ncomp,userdata,iternumber,i

    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,S4Lim,S4S5Lim,S5Lim,TCLim,TCS4Lim,TCS4S5Lim,TCS5Lim, &
         C5Lim,C5S4Lim,C5S4S5Lim,C5S5Lim,TCC5Lim,TCC5S5Lim,TCC5S4S5Lim,TCC5S4Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR
    real(dp)    :: xi1,xi2,f1(-6:7),f2(-6:7),my_lumi,respdf(1),respdf_b(1)
    !--


    !--



    real(dp)    :: respdf_qed(1),respdf_qcd(1)
    real(dp)    :: reseik_qed,reseik_qcd
    !--
    real(dp) :: z1,z2,z4,z5,s24,s25,s45,s4_25,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_dy_ns(16),kin(1:9),skin(1:6)
    logical  :: makecut
    !--
    real(dp) :: qcdmom(4,3),qedmom(4,3)
    real(dp) :: mysymm = half
    !
    real(dp) :: eikvec_s5(4),respdf_s5
    real(dp) :: splitvec_c5(4),respdf_c5
    real(dp) :: eikvec_c5s5(4),respdf_c5s5
    real(dp) :: eikvec_ds(4),respdf_ds
    real(dp) :: eikvec_dss5(4),respdf_dss5
    real(dp) :: eikvec_c5ds(4),respdf_c5ds
    real(dp) :: eikvec_c5dss5(4),respdf_c5dss5
    real(dp) :: eikvec_tcds(4),respdf_tcds
    real(dp) :: eikvec_tcdss5(4),respdf_tcdss5
    real(dp) :: eikvec_tcc5ds(4),respdf_tcc5ds
    real(dp) :: eikvec_tcc5dss5(4),respdf_tcc5dss5
    real(dp) :: eikvec_tcs5(4),respdf_tcs5
    real(dp) :: eikvec_tcc5s5(4),respdf_tcc5s5
    !
    real(dp) :: wmom(4)
    real(dp) :: respdf_ch1(1)
    !
    real(dp) :: resnnlo_2q(32)

    real(dp) :: resnnlo(2,2)     

    real(dp) :: resnnlo_4q_tc(4),resnnlo_4q_fin(8)
    real(dp) :: res_a_nlo(2,2), res_g_nlo(32)
    real(dp) :: vec_a(4),vec_g(4),vec_ag(4),vec_ga(4)
    real(dp) :: vec_ag_a(4),vec_ag_g(4),vec_ga_a(4),vec_ga_g(4)
    real(dp) :: resvec(4)
    real(dp) :: res_ag(4),res_ga(4)
    real(dp) :: respdf_ag(1),respdf_ga(1)
    real(dp) :: myz
    real(dp) :: myres(4),myrespdf(1),res_gg(4),res_qqb(4)
    !
    real(dp) :: myx1,myx2
    real(dp) :: myza1,myza2,myzb1,myzb2
    real(dp) :: mys145,dbsing

    real(dp) :: tag_ns_wga

    real(dp) :: eik(4), res2(2,2), split(2,2), res_s4c5_2(2,2), res_s4c5(2,2), res_S5_nlo(2,2), ker(2,2)
    real(dp) :: res_S5S4(2,2), res_TC(2,2), res_TC2(2,2), res_TCS5(2,2),res_TCS52(2,2), res_TCC5S52(2,2)
    real(dp) :: res_TCS4(2,2),res_TCS42(2,2)
    real(dp) :: res_C5(2,2), res_C52(2,2), test(2,2), temp(2,2)

#if(_identification_check == 1)
    logical    :: check, check_temp
    !the debug arrays contains hard and single unresolved contirbutions
    ! ( qqb -> ga  , qqb -> ag  , qqb->qqb , a+S5g , g+S5a , a+C5g , g+C5a )
    integer :: debug_binarray(7,nobs)
    real(dp) :: debug_weightarray(7,nobs)
    real(dp) :: debug_kin(7)
    debug_binarray    = zero
    debug_weightarray = zero
    debug_kin         = zero
    check = .true.
#endif    

    tag_ns_wga = one

    callcount = callcount + 1

    xsect_nnlo_5262a_rr_z_ga = 0
    ff(1) = zero
    xx = zero


    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buff_rr) return 


    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call open_histo()    
    
    call kinematics_nnlo_5262a(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S4Lim,S5Lim,S4S5Lim,TCLim,TCS4Lim,TCS5Lim,TCS4S5Lim,&
         C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,TCC5Lim,TCC5S4Lim,TCC5S5Lim,TCC5S4S5Lim)  

    !--------------------------------------
    !    NNLO ME :   gluon4,  photon 5
    !--------------------------------------

    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)

    call cut_histo(HardProc)

    if (makecut.or.HardProc%flag) then
       kin(1) = zero
       FintNNLO_dy_ns(1) = zero 

    else    

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_ga_qqb(HardProc%AmpMom,resnnlo)

       !call get_respdf(ns_lumi,1,1,HardProc,resnnlo,respdf)
       call get_respdf_ns(1,1,HardProc,resnnlo,respdf)

       kin(1) = respdf(1)*HardProc%wgt*jac
       FintNNLO_dy_ns(1) = kin(1)
       call fill_histo([skin(1)],vegasweight)

    endif

    !! -----------------------------!!
    !!     SOFT   4, gluon           !!
    !! -----------------------------!!


    !-- 2 -> S4       particle ''4'' is a gluon 
    xi1 = S4Lim%PartFrac(1) 
    xi2 = S4Lim%PartFrac(2)

    call cut_histo(S4Lim)

    if (makecut.or.S4Lim%flag) then
       kin(2) = zero
       FintNNLO_dy_ns(2) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)
       call res_tree_a_qqb(S4Lim%AmpMom,res_a_nlo)

       !call get_respdf(ns_lumi,1,1,S4Lim,res_a_nlo,respdf)
       call get_respdf_ns(1,1,S4Lim,res_a_nlo,respdf)

       reseik_qcd = myEik_g(S4Lim%LimMom(:,1),S4Lim%LimMom(:,2),S4Lim%LimMom(:,3),Cf)
       print *, S4Lim%muf,S4Lim%mur

       FintNNLO_dy_ns(2) = (-one)*respdf(1) *reseik_qcd*S4Lim%wgt*jac 
       kin(2) = FintNNLO_dy_ns(2)

       call fill_histo([FintNNLO_dy_ns(2)],vegasweight)

    endif



    !! -----------------------------
    !!      SOFT   4 (gluon)  COLLINEAR 5 (photon)
    !!
    !! -----------------------------

    !print *, ''
    !print *, 'S4 C5'
    !print *, ''


    xi1 = C5S4Lim%PartFrac(1) 
    xi2 = C5S4Lim%PartFrac(2)

    call cut_histo(C5S4Lim)

    if (makecut.or.C5S4Lim%flag) then
       kin(3) = zero
       FintNNLO_dy_ns(3) = zero 
    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C5S4Lim%AmpMom,res_s4c5)          
       reseik_qcd = myEik_g(C5S4Lim%LimMom(:,1),C5S4Lim%LimMom(:,2),C5S4Lim%LimMom(:,3),Cf)

       s25 = C5S4Lim%Lim_KinInv(2) 
       z1  = C5S4Lim%Lim_KinInv(5)

       call AP_nlo_ph(z1,one,split)
       call product_for_pdf_matrix(res_s4c5,split,res_s4c5_2)

       !call get_respdf(ns_lumi,1,1,C5S4Lim,res_s4c5_2,respdf)
       call get_respdf_ns(1,1,C5S4Lim,res_s4c5_2,respdf)
       respdf(1) = respdf(1) * reseik_qcd*(-two/s25)

       FintNNLO_dy_ns(3) = respdf(1)*C5S4Lim%wgt*jac
       kin(3) = FintNNLO_dy_ns(3)

       call fill_histo([FintNNLO_dy_ns(2)],vegasweight)

    endif



    !! --------------------------------------!!
    !!     SOFT    5, photon     &   C5 S5  !!
    !! -------------------------------------!!

    !print *, ''
    !print *, 'S5'
    !print *, ''


    xi1 = S5Lim%PartFrac(1) 
    xi2 = S5Lim%PartFrac(2)

    call cut_histo(S5Lim)

    if (makecut.or.S5Lim%flag) then
       kin(4) = zero
       FintNNLO_dy_ns(4) = zero
       FintNNLO_dy_ns(5) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_g_qqb(S5Lim%AmpMom,res_S5_nlo)   
       call get_Eik_ph_z_ga(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3), S5Lim%LimMom(:,4), S5Lim%LimMom(:,5),eik) 
       call product_for_pdf(res_S5_nlo,eik,res2)

       !call get_respdf(ns_lumi,1,1,S5Lim,res2,respdf_ga)
       call get_respdf_ns(1,1,S5Lim,res2,respdf_ga)

       FintNNLO_dy_ns(4) = (-one)*respdf_ga(1)*S5Lim%wgt*jac


       !-------- C5S5


       !print *, ''
       !print *, 'C5S5'
       !print *, ''

       s25 = C5S5Lim%Lim_KinInv(2)
       z5 = C5S5Lim%Lim_KinInv(5)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2) 

       call product_for_pdf_matrix(res_S5_nlo,ker,res2)          
       !call get_respdf(ns_lumi,1,1,S5Lim,res2,respdf_ga)
       call get_respdf_ns(1,1,S5Lim,res2,respdf_ga)

       FintNNLO_dy_ns(5) =  respdf_ga(1)*four/z5/s25 * C5S5Lim%wgt*jac

       kin(4) = FintNNLO_dy_ns(4) + FintNNLO_dy_ns(5)          
       call fill_histo([FintNNLO_dy_ns(2) + FintNNLO_dy_ns(3)],vegasweight)

    endif


    !! -----------------------------!!
    !!     C5 -- photon             !!
    !! -----------------------------!!

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)

    call cut_histo(C5Lim)

    if (makecut.or.C5Lim%flag) then
       kin(5) = zero
       FintNNLO_dy_ns(6)  = zero
    else


       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2) 

       call res_tree_g_qqb(C5Lim%AmpMom,res_C5)
       call product_for_pdf_matrix(res_C5,ker,res_C52)

       !call get_respdf(ns_lumi,1,1,C5Lim,res_C52,respdf_ch1)
       call get_respdf_ns(1,1,C5Lim,res_C52,respdf_ch1)

       s25 = C5Lim%Lim_KinInv(2) 
       z1  = C5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(6)  = (-one)*(-two/s25)*( respdf_ch1(1) * Pqq(z1,one)) * C5Lim%wgt*jac

       kin(5) = FintNNLO_dy_ns(6) 
       call fill_histo([kin(5)],vegasweight)

    endif


    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!

    !-- 7 -> S4S5
    !-- 8 -> TCS4S5
    !-- 9 -> C5S4S5 
    !-- 10 -> TCC5S4S5

    !print *, ''
    !print *, 'S4S5'
    !print *, ''

    xi1 = S4S5Lim%PartFrac(1) 
    xi2 = S4S5Lim%PartFrac(2)

    call cut_histo(S4S5Lim) 

    if (makecut.or.S4S5Lim%flag) then

       kin(6) = zero
       FintNNLO_dy_ns(7:10) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)
       call res_tree_qqb(S4S5Lim%AmpMom,res_S5S4) 

       call get_Eik_ds_z_ga(S4S5Lim%LimMom(:,1),S4S5Lim%LimMom(:,2),S4S5Lim%LimMom(:,3),&
            S4S5Lim%LimMom(:,4),S4S5Lim%LimMom(:,5),S4S5Lim%LimMom(:,6),eik)
       call product_for_pdf(res_S5S4,eik,res2)

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ag)
       call get_respdf_ns(1,1,S4S5Lim,res2,respdf_ag)     
       FintNNLO_dy_ns(7)  = (-one)*respdf_ag(1)*S4S5Lim%wgt*jac


       !print *, ''
       !print *, 'TCS4S5'
       !print *, ''

       s24 = TCS4S5Lim%Lim_KinInv(1)
       s25 = TCS4S5Lim%Lim_KinInv(2)
       s45 = TCS4S5Lim%Lim_KinInv(3)
       !
       z4 = TCS4S5Lim%Lim_KinInv(4)
       z5 = TCS4S5Lim%Lim_KinInv(5)

       call get_Pggq_ds_z_ga(s45,s24,s25,z4,z5,split)
       call product_for_pdf_matrix(res_S5S4,split,res2)

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf) ! to be checked
       call get_respdf_ns(1,1,S4S5Lim,res2,respdf) ! to be checked
       FintNNLO_dy_ns(8) = (-one)*respdf(1)*TCS4S5Lim%wgt*jac ! to be checked

       !print *, ''
       !print *, 'S4S5C5'
       !print *, ''

       s25 = C5S4S5Lim%Lim_KinInv(2)
       z5  = C5S4S5Lim%Lim_KinInv(5)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2) 

       call product_for_pdf_matrix(res_S5S4,ker,res2)

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf)
       call get_respdf_ns(1,1,S4S5Lim,res2,respdf)

       FintNNLO_dy_ns(9) = (respdf(1))* &
            myEik_g(C5S4S5Lim%LimMom(:,1),C5S4S5Lim%LimMom(:,2),C5S4S5Lim%LimMom(:,3),Cf)*four/z5/s25 &
            *C5S4S5Lim%wgt*jac

       !print *, ''
       !print *, 'TCC5S4S5'
       !print *, ''

       s24 = TCC5S4S5Lim%Lim_KinInv(1)
       s25 = TCC5S4S5Lim%Lim_KinInv(2)
       !
       z4 = TCC5S4S5Lim%Lim_KinInv(4)
       z5 = TCC5S4S5Lim%Lim_KinInv(5)


       FintNNLO_dy_ns(10) = (respdf(1))*four*Cf/(-s24)/z4*four/(-s25)/z5 * TCC5S4S5Lim%wgt*jac

       kin(6) = sum(FintNNLO_dy_ns(7:10))

       call fill_histo([kin(6)],vegasweight)

    endif

    !-- 11 -> TC
    !-- 12 -> TCC5

    !print *, ''
    !print *, 'TC'
    !print *, ''

    xi1 = TCLim%PartFrac(1) 
    xi2 = TCLim%PartFrac(2)
    call cut_histo(TCLim)

    if (makecut.or.TCLim%flag) then

       kin(7) = zero
       FintNNLO_dy_ns(11:12) = zero 

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       s24 = TCLim%Lim_KinInv(1)
       s25 = TCLim%Lim_KinInv(2)
       s45 = TCLim%Lim_KinInv(3)
       !
       z4 = TCLim%Lim_KinInv(4)
       z5 = TCLim%Lim_KinInv(5)
       z1 = one-z4-z5

       call res_tree_qqb(TCLim%AmpMom,res_TC)
       call AP_nnlo_ga(s45,s24,s25,z4,z5,z1,split)
       call product_for_pdf_matrix(res_TC,split,res_TC2)

       !call get_respdf(ns_lumi,1,1,TCLim,res_TC2,respdf)
       call get_respdf_ns(1,1,TCLim,res_TC2,respdf)

       FintNNLO_dy_ns(11) = (-one)*respdf(1)*TCLim%wgt*jac


       !print *, ''
       !print *, 'TC C5'
       !print *, ''

       s25 = TCC5Lim%Lim_KinInv(1)
       s4_25 = TCC5Lim%Lim_KinInv(2)
       z1 = TCC5Lim%Lim_KinInv(3)
       z2 = TCC5Lim%Lim_KinInv(4)

       call get_Pggq_C_z_ga(s25,s4_25,z1,z2,esq,split)
       call product_for_pdf_matrix(res_TC,split,res_TC2)
       !call get_respdf(ns_lumi,1,1,TCLim,res_TC2,respdf_ch1)
       call get_respdf_ns(1,1,TCLim,res_TC2,respdf_ch1)

       FintNNLO_dy_ns(12) = respdf_ch1(1)*TCC5Lim%wgt*jac

       kin(7) = sum(FintNNLO_dy_ns(11:12))

       call fill_histo([kin(7)],vegasweight)



    endif


    !-- 13 -> TCS5
    !-- 14 -> TCC5S5

    !print *, ''
    !print *, 'TC S5'
    !print *, ''

    xi1 = TCS5Lim%PartFrac(1) 
    xi2 = TCS5Lim%PartFrac(2)
    call cut_histo(TCS5Lim)

    if (makecut.or.TCS5Lim%flag) then

       kin(8) = zero
       FintNNLO_dy_ns(13:14)  = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(TCS5Lim%AmpMom,res_TCS5)

       s24 = TCS5Lim%Lim_KinInv(1)
       s25 = TCS5Lim%Lim_KinInv(2)
       z4 = TCS5Lim%Lim_KinInv(4)
       z5 = TCS5Lim%Lim_KinInv(5)

       call get_Pggq_S_z_ga(TCS5Lim%Lim_KinInv,split)
       call product_for_pdf_matrix(res_TCS5,split,res_TCS52)

       !call get_respdf(ns_lumi,1,1,TCS5Lim,res_TCS52,respdf_ch1 )
       call get_respdf_ns(1,1,TCS5Lim,res_TCS52,respdf_ch1) 
       FintNNLO_dy_ns(13) = respdf_ch1(1)*TCS5Lim%wgt*jac

       !print *, ''

       s24 = TCC5S5Lim%Lim_KinInv(1)
       s25 = TCC5S5Lim%Lim_KinInv(2)
       z4 =  TCC5S5Lim%Lim_KinInv(4)
       z5 =  TCC5S5Lim%Lim_KinInv(5)


       call get_Pggq_SC_z_ga(s24,s25,z4,z5,esq,split)          
       call product_for_pdf_matrix(res_TCS5,split,res_TCC5S52)

       !call get_respdf(ns_lumi,1,1,TCS5Lim,res_TCC5S52,respdf_ch1)
       call get_respdf_ns(1,1,TCS5Lim,res_TCC5S52,respdf_ch1)
       FintNNLO_dy_ns(14) = respdf_ch1(1) * TCC5S5Lim%wgt*jac

       kin(8) = sum(FintNNLO_dy_ns(13:14))

       call fill_histo([kin(8)],vegasweight)

    endif

    !-- 15 -> TCS4
    !-- 16 -> TCC5S4

    !print *, ''
    !print *, 'TC S4'
    !print *, ''

    xi1 = TCS4Lim%PartFrac(1) 
    xi2 = TCS4Lim%PartFrac(2)
    call cut_histo(TCS4Lim)

    if (makecut.or.TCS4Lim%flag) then

       kin(9) = zero
       FintNNLO_dy_ns(15:16)  = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2) 

       call res_tree_qqb(TCS4Lim%AmpMom,res_TCS4)
       call product_for_pdf_matrix(res_TCS4,ker,res_TCS42)          
       !call get_respdf(ns_lumi,1,1,TCS4Lim,res_TCS42,respdf_ch1)
       call get_respdf_ns(1,1,TCS4Lim,res_TCS42,respdf_ch1)

       s24 = TCS4Lim%Lim_KinInv(1)
       s25 = TCS4Lim%Lim_KinInv(2)
       z4 = TCS4Lim%Lim_KinInv(4)
       z5 = TCS4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(15) = (-one)*respdf_ch1(1)* four*CF/s24/z4 * two/s25*Pqq(z5,one)*TCS4Lim%wgt*jac

       s24 = TCC5S4Lim%Lim_KinInv(1)
       s25 = TCC5S4Lim%Lim_KinInv(2)
       z4 =  TCC5S4Lim%Lim_KinInv(4)
       z5 =  TCC5S4Lim%Lim_KinInv(5)


       FintNNLO_dy_ns(16) = (-four)*CF/s24/z4*two/s25*Pqq(z5,one) *respdf_ch1(1) * TCC5S4Lim%wgt*jac


       kin(9) = sum(FintNNLO_dy_ns(15:16))

    endif

    call fill_histo([kin(6)],vegasweight)

    ff(1) = sum(kin)
    call close_histo()

    !     print*, 'xx=', xx

    !    print  *, 1, FintNNLO_dy_ns(1), 'hard'
    !    print  *, '------------------------------------------------------'
    !    print  *, 2, FintNNLO_dy_ns(2), FintNNLO_dy_ns(2)/FintNNLO_dy_ns(1), 'x1=0'
    !    print  *, 3, FintNNLO_dy_ns(4), FintNNLO_dy_ns(4)/FintNNLO_dy_ns(1), 'x2=0'
    !    print  *, 4, FintNNLO_dy_ns(11), FintNNLO_dy_ns(11)/FintNNLO_dy_ns(1), 'x3=0'
    !    print  *, 5, FintNNLO_dy_ns(6), FintNNLO_dy_ns(6)/FintNNLO_dy_ns(1), 'x4=0'
    !    print  *, '------------------------------------------------------'
    !    print  *, 6, FintNNLO_dy_ns(7), FintNNLO_dy_ns(7)/FintNNLO_dy_ns(1),    'x1=0,x2=0'
    !    print  *, 7, FintNNLO_dy_ns(15), FintNNLO_dy_ns(15)/FintNNLO_dy_ns(1), 'x1=0,x3=0'
    !    print  *, 8, FintNNLO_dy_ns(3), FintNNLO_dy_ns(3)/FintNNLO_dy_ns(1), 'x1=0,x4=0'
    !    print  *, 9, FintNNLO_dy_ns(13), FintNNLO_dy_ns(13)/FintNNLO_dy_ns(1), 'x2=0, x3=0'
    !    print  *, 10, FintNNLO_dy_ns(5), FintNNLO_dy_ns(5)/FintNNLO_dy_ns(1), 'x2=0,x4=0'
    !    print  *, 11, FintNNLO_dy_ns(12), FintNNLO_dy_ns(12)/FintNNLO_dy_ns(1), 'x3=0,x4=0'
    !    print  *, '------------------------------------------------------'     
    !    print  *, 12, FintNNLO_dy_ns(8), FintNNLO_dy_ns(8)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0'
    !    print  *, 13, FintNNLO_dy_ns(9), FintNNLO_dy_ns(9)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x4=0'
    !    print  *, 14, FintNNLO_dy_ns(16), FintNNLO_dy_ns(16)/FintNNLO_dy_ns(1), 'x1=0,x3=0,x4=0'
    !    print  *, 15, FintNNLO_dy_ns(14), FintNNLO_dy_ns(14)/FintNNLO_dy_ns(1), 'x2=0,x3=0,x4=0'
    !    print  *, '------------------------------------------------------'     
    !    print  *, 16, FintNNLO_dy_ns(10), FintNNLO_dy_ns(10)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0,x4=0'

    ! !!$    !         do i = 1,16
    ! !!$    !            !print  *, i,FintNNLO_dy_ns(i)
    ! !!$    !         enddo

    ! print  *, 'ff(1)', ff(1)
    !          pause
    ! !!$


    !      pause

    if (ff(1) .ne. ff(1)) then
       print  *, 'nan'
       print  *, 'xx',xx
       print  *, ''
       do i = 1,16
          print  *, i, FintNNLO_dy_ns(i)
       enddo
       ff(1) = zero
    endif

#if(_withchecks == 1)
    FintNNLO = FintNNLO_dy_ns
#endif

    return

  end function xsect_nnlo_5262a_rr_z_ga



  function xsect_nnlo_5262c_rr_z_ga(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer ::  xsect_nnlo_5262c_rr_z_ga,ndim,ncomp,userdata,iternumber,i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,S4Lim,S4S5Lim,S5Lim,TCLim,TCS4Lim,TCS4S5Lim,TCS5Lim,&
         C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,TCC4Lim,TCC4S4Lim,TCC4S4S5Lim,TCC4S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR
    real(dp)    :: xi1,xi2,f1(-6:7),f2(-6:7),my_lumi,respdf(1),respdf_b(1)
    complex(dp) :: respdf_sping(-1:1,-1:1)
    !--


    !--



    real(dp)    :: respdf_qed(1),respdf_qcd(1)
    complex(dp) :: respdf_sping_qed(-1:1,-1:1),respdf_sping_qcd(-1:1,-1:1)
    real(dp)    :: respdf_ec,respdf_ce,respdf_b_qed(1),respdf_b_qcd(1)
    real(dp)    :: reseik_qed,reseik_qcd
    !--
    real(dp) ::z1,z2,z4,z5,s24,s25,s45,s5_24,s12,s15,kt1kt2sq,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_dy_ns(16),kin(1:9), FintNNLO_dy_ns2(1)
    logical  :: makecut
    !--
    real(dp) :: mysymm = half
    real(dp) :: qcdmom(4,3),qedmom(4,3)
    !
    real(dp) :: wmom(4)
    real(dp) :: respdf_ch1(1)
    !
    real(dp) :: resnnlo_2q(32)
    real(dp) :: resnnlo_4q_tc(4),resnnlo_4q_fin(8)
    real(dp) :: res_g_nlo(12)
    real(dp) :: vec_a(4),vec_g(4),vec_ag(4),vec_ga(4)
    real(dp) :: vec_ag_a(4),vec_ag_g(4),vec_ga_a(4),vec_ga_g(4)
    real(dp) :: resvec(4)
    real(dp) :: res_ag(4),res_ga(4)
    real(dp) :: respdf_ag(1),respdf_ga(1)

    real(dp) :: tag_ns_wga

    real(dp) ::resnnlo(2,2), res_a_nlo(2,2), temp(2,2), res_S5_nlo(2,2), eik(4), res2(2,2), res_S5S4(2,2)
    real(dp) :: split(2,2), ker(2,2), res_TC(2,2), res_TC2(2,2)
    real(dp) :: res_C4S5(2,2), res_C4S52(2,2), res_TCS52(2,2), res_TCS5(2,2), res_TCS4(2,2), res_TCS42(2,2)

#if(_identification_check == 1)
    logical    :: check, check_temp
    !the debug arrays contains hard and single unresolved contirbutions
    ! ( qqb -> ga  , qqb -> ag  , qqb->qqb , a+S5g , g+S5a a+C4g , g+C4a )
    integer :: debug_binarray(7,nobs)
    real(dp) :: debug_weightarray(7,nobs)
    real(dp) :: debug_kin(7)
    debug_binarray    = zero
    debug_weightarray = zero
    debug_kin         = zero
    check = .true.
#endif

    callcount = callcount + 1


    tag_ns_wga = one


    xsect_nnlo_5262c_rr_z_ga = 0
    ff(1) = zero
    kin = zero
    xx = zero


    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buff_rr) return 

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call open_histo()

    call kinematics_nnlo_5262c(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S4Lim,S5Lim,S4S5Lim,TCLim,TCS4Lim,TCS5Lim,TCS4S5Lim,&
         C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,TCC4Lim,TCC4S4Lim,TCC4S5Lim,TCC4S4S5Lim)

    !! ----------------------------------------!!
    !!    calling NNLO ME; gluon4 photon 5      !!
    !! ----------------------------------------!!

    !1-- regular piece 

    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo(HardProc)

    if (makecut.or.HardProc%flag) then
       kin(1) = zero

       FintNNLO_dy_ns(1) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_ga_qqb(HardProc%AmpMom,resnnlo)

       !call get_respdf(ns_lumi,1,1,HardProc,resnnlo,respdf)
       call get_respdf_ns(1,1,HardProc,resnnlo,respdf)

       kin(1) = respdf(1)*HardProc%wgt*jac
       FintNNLO_dy_ns(1) = kin(1)

       call fill_histo([kin(1)],vegasweight)

    endif


    !! -----------------------------!!
    !!     SOFT   4, gluon          !!
    !! -----------------------------!!
    !-- 2 -> S4       particle ''4'' is a gluon 

    xi1 = S4Lim%PartFrac(1) 
    xi2 = S4Lim%PartFrac(2) 

    call cut_histo(S4Lim)

    if (makecut.or.S4Lim%flag) then
       kin(2) = zero
       FintNNLO_dy_ns(2) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_a_qqb(S4Lim%AmpMom,res_a_nlo)
       !call get_respdf(ns_lumi,1,1,S4Lim,res_a_nlo,respdf)
       call get_respdf_ns(1,1,S4Lim,res_a_nlo,respdf)

       reseik_qcd = myEik_g(S4Lim%LimMom(:,1),S4Lim%LimMom(:,2),S4Lim%LimMom(:,3),Cf)

       FintNNLO_dy_ns(2) = (-one)*respdf(1) *reseik_qcd*S4Lim%wgt*jac
       kin(2) = FintNNLO_dy_ns(2)

       call fill_histo([kin(2)],vegasweight)

    endif



    !-------- S4 C4

    xi1 = C4S4Lim%PartFrac(1) 
    xi2 = C4S4Lim%PartFrac(2) 

    call cut_histo(C4S4Lim)

    if (makecut.or.C4S4Lim%flag) then
       FintNNLO_dy_ns(3) = zero

       kin(3) = zero 
    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_a_qqb(C4S4Lim%AmpMom,res_a_nlo)

       !call get_respdf(ns_lumi,1,1,C4S4Lim,res_a_nlo,respdf)
       call get_respdf_ns(1,1,C4S4Lim,res_a_nlo,respdf)

       s24 = C4S4Lim%Lim_KinInv(2)  
       z4 =  C4S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(3) = respdf(1) *four*Cf/s24/z4*C4S4Lim%wgt*jac
       kin(3) = FintNNLO_dy_ns(3)


       call fill_histo([kin(3)],vegasweight)

    endif


    !! -----------------------------!!
    !!    Soft 5,  photon            !!
    !! -----------------------------!!

    !-- 2 -> S5 ; 5 is a photon

    xi1 = S5Lim%PartFrac(1) 
    xi2 = S5Lim%PartFrac(2) 
    call cut_histo(S5Lim)

    if (makecut.or.S5Lim%flag) then

       kin(4) = zero
       FintNNLO_dy_ns(4) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_g_qqb(S5Lim%AmpMom,res_S5_nlo)
       call get_Eik_ph_z_ga(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3), &
            S5Lim%LimMom(:,4), S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5_nlo,eik,res2)
       !call get_respdf(ns_lumi,1,1,S5Lim,res2,respdf_ga)
       call get_respdf_ns(1,1,S5Lim,res2,respdf_ga)

       FintNNLO_dy_ns(4) = (-one)*respdf_ga(1)*S5Lim%wgt*jac*tag_ns_wga

       kin(4)  = FintNNLO_dy_ns(4) 
       call fill_histo([kin(4)],vegasweight)
    endif

    !! -----------------------------!!
    !!     C4 -- gluon !!
    !! -----------------------------!!

    !-- 4 -> C4 ; 

    xi1 = C4Lim%PartFrac(1) 
    xi2 = C4Lim%PartFrac(2)
    call cut_histo(C4Lim)


    if (makecut.or.C4Lim%flag) then
       kin(5) = zero
       FintNNLO_dy_ns(5) = zero  

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_a_qqb(C4Lim%AmpMom,res_a_nlo)
       !call get_respdf(ns_lumi,1,1,C4Lim,res_a_nlo,respdf_b_qed)
       call get_respdf_ns(1,1,C4Lim,res_a_nlo,respdf_b_qed)

       s24 = C4Lim%Lim_KinInv(2) 
       z1  = C4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(5)  = (-one)*(-two/s24)* Pqq(z1,Cf) * respdf_b_qed(1) *C4Lim%wgt*jac


       kin(5) = FintNNLO_dy_ns(5) 
       call fill_histo([kin(4)],vegasweight)

    endif


    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!


    !-- 6 -> S4S5
    !-- 7 -> TCS4S5
    !-- 8 -> C4S4S5
    !-- 9 -> TCC4S4S5

    xi1 = S4S5Lim%PartFrac(1) 
    xi2 = S4S5Lim%PartFrac(2)

    call cut_histo(S4S5Lim)

    if (makecut.or.S4S5Lim%flag) then
       kin(6) = zero
       FintNNLO_dy_ns(6:9)  = zero


    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(S4S5Lim%AmpMom,res_S5S4) 

       call get_Eik_ds_z_ga(S4S5Lim%LimMom(:,1),S4S5Lim%LimMom(:,2),S4S5Lim%LimMom(:,3),&
            S4S5Lim%LimMom(:,4),S4S5Lim%LimMom(:,5),S4S5Lim%LimMom(:,6),eik)
       call product_for_pdf(res_S5S4,eik,res2)
       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ag)
       call get_respdf_ns(1,1,S4S5Lim,res2,respdf_ag)

       FintNNLO_dy_ns(6) = (-one)*(respdf_ag(1))*S4S5Lim%wgt*jac


       !------- TCS4S5

       s24 = TCS4S5Lim%Lim_KinInv(1)
       s25 = TCS4S5Lim%Lim_KinInv(2)
       s45 = TCS4S5Lim%Lim_KinInv(3)
       !
       z4 = TCS4S5Lim%Lim_KinInv(4)
       z5 = TCS4S5Lim%Lim_KinInv(5)

       call get_Pggq_ds_z_ga(s45,s24,s25,z4,z5,split)
       call product_for_pdf_matrix(res_S5S4,split,res2)

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf)
       call get_respdf_ns(1,1,S4S5Lim,res2,respdf)

       FintNNLO_dy_ns(7) = (-one)*respdf(1)*TCS4S5Lim%wgt*jac 


       !----->  C4S4S5

       call get_Eik_ph_z_ga(C4S4S5Lim%LimMom(:,1),C4S4S5Lim%LimMom(:,2),C4S4S5Lim%LimMom(:,3), &
            C4S4S5Lim%LimMom(:,4), C4S4S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5S4,eik,res2)

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ga)
       call get_respdf_ns(1,1,S4S5Lim,res2,respdf_ga)

       FintNNLO_dy_ns(8) = respdf_ga(1)*4.0_dp*Cf*esq/z4/s24*C4S4S5Lim%wgt*jac

       !----->  TCC4S4S5

       s12 = TCC4S4S5Lim%Lim_KinInv(1) 
       s24 = TCC4S4S5Lim%Lim_KinInv(2) 
       s25 = TCC4S4S5Lim%Lim_KinInv(3) 
       s15 = TCC4S4S5Lim%Lim_KinInv(4) 
       z4 =  TCC4S4S5Lim%Lim_KinInv(5)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2)

       call product_for_pdf_matrix(res_S5S4,ker,res2)
       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ch1)
       call get_respdf_ns(1,1,S4S5Lim,res2,respdf_ch1)


       FintNNLO_dy_ns(9) =  (-one)*  16.0_dp*Cf/s24/z4 * s12/s25/s15  *respdf_ch1(1) * TCC4S4S5Lim%wgt*jac


       kin(6) = sum(FintNNLO_dy_ns(6:9))

       call fill_histo([kin(6)],vegasweight)

    endif

    !-- 10 -> TC
    !-- 11 -> TCC4

    xi1 = TCLim%PartFrac(1) 
    xi2 = TCLim%PartFrac(2)
    call cut_histo(TCLim)

    if (makecut.or.TCLim%flag) then

       kin(7) = zero
       FintNNLO_dy_ns(10:11)  = zero


    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       s24 = TCLim%Lim_KinInv(1)
       s25 = TCLim%Lim_KinInv(2)
       s45 = TCLim%Lim_KinInv(3)
       !
       z4 = TCLim%Lim_KinInv(4)
       z5 = TCLim%Lim_KinInv(5)
       z1 = one-z4-z5

       call res_tree_qqb(TCLim%AmpMom,res_TC)
       call AP_nnlo_ga(s45,s24,s25,z4,z5,z1,split)
       call product_for_pdf_matrix(res_TC,split,res_TC2)
       !call get_respdf(ns_lumi,1,1,TCLim,res_TC2,respdf)
       call get_respdf_ns(1,1,TCLim,res_TC2,respdf)


       FintNNLO_dy_ns(10) = (-one)*respdf(1)*TCLim%wgt*jac


       s24 = TCC4Lim%Lim_KinInv(1) 
       s5_24 = TCC4Lim%Lim_KinInv(2) 
       z1  = TCC4Lim%Lim_KinInv(3) 
       z2 = TCC4Lim%Lim_KinInv(4)

       call get_Pggq_C_z_ga(s24,s5_24,z1,z2,esq,split)
       call product_for_pdf_matrix(res_TC,split,res_TC2)
       !call get_respdf(ns_lumi,1,1,TCLim,res_TC2,respdf_ch1)
       call get_respdf_ns(1,1,TCLim,res_TC2,respdf_ch1)

       FintNNLO_dy_ns(11) = respdf_ch1(1)*TCC4Lim%wgt*jac 

       kin(7) = sum(FintNNLO_dy_ns(10:11))

       call fill_histo([kin(5)],vegasweight)
    endif



    !-- 12 --> C4S5
    !-- 13 -> TCS5
    !-- 14 -> TCC4S5

    xi1 = TCS5Lim%PartFrac(1) 
    xi2 = TCS5Lim%PartFrac(2)
    call cut_histo(TCS5Lim)

    if (makecut.or.TCS5Lim%flag) then

       kin(8) = zero
       FintNNLO_dy_ns(12:14)  = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C4S5Lim%AmpMom,res_C4S5)
       call get_Eik_ph_z_ga(C4S5Lim%LimMom(:,1),C4S5Lim%LimMom(:,2),C4S5Lim%LimMom(:,3), &
            C4S5Lim%LimMom(:,4),C4S5Lim%LimMom(:,5),eik)

       call product_for_pdf(res_C4S5,eik,res_C4S52)
       !call get_respdf(ns_lumi,1,1,TCS5Lim,res_C4S52,respdf_ga)
       call get_respdf_ns(1,1,TCS5Lim,res_C4S52,respdf_ga)

       s12 = C4S5Lim%Lim_KinInv(1) 
       s25 = C4S5Lim%Lim_KinInv(2) 
       s15 = C4S5Lim%Lim_KinInv(3) 
       s24 = C4S5Lim%Lim_KinInv(4) 
       z4 =  C4S5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(12) = (- respdf_ga(1))*two/s24*Pqq(z4,Cf)* C4S5Lim%wgt*jac


       call res_tree_qqb(TCS5Lim%AmpMom,res_TCS5)

       s24 = TCS5Lim%Lim_KinInv(1)
       s25 = TCS5Lim%Lim_KinInv(2)
       z4 = TCS5Lim%Lim_KinInv(4)
       z5 = TCS5Lim%Lim_KinInv(5)

       call get_Pggq_S_z_ga(TCS5Lim%Lim_KinInv,split)
       call product_for_pdf_matrix(res_TCS5,split,res_TCS52)

       !call get_respdf(ns_lumi,1,1,TCS5Lim,res_TCS52,respdf_ch1)
       call get_respdf_ns(1,1,TCS5Lim,res_TCS52,respdf_ch1) 
       FintNNLO_dy_ns(13) = respdf_ch1(1)*TCS5Lim%wgt*jac


       call product_for_pdf_matrix(res_TCS5,ker,res_TCS52)
       !call get_respdf(ns_lumi,1,1,TCS5Lim,res_TCS52,respdf_ch1)
       call get_respdf_ns(1,1,TCS5Lim,res_TCS52,respdf_ch1) 

       s24 = TCC4S5Lim%Lim_KinInv(1) 
       s25 = TCC4S5Lim%Lim_KinInv(2) 
       s45 = TCC4S5Lim%Lim_KinInv(3) 
       z4 =  TCC4S5Lim%Lim_KinInv(4) 
       z5 =  TCC4S5Lim%Lim_KinInv(5) 
       !

       FintNNLO_dy_ns(14) = (-one)*8.0_dp/s24/s25/z5*Pqq(z4,Cf)/z4 *respdf_ch1(1) *TCC4S5Lim%wgt*jac


       kin(8) = sum(FintNNLO_dy_ns(12:14))

    endif




    !-- 15 ->  TCS4
    !-- 16 -> TCC4S4

    xi1 = TCS4Lim%PartFrac(1) 
    xi2 = TCS4Lim%PartFrac(2)
    call cut_histo(TCS4Lim)

    if (makecut.or.TCS4Lim%flag) then

       kin(9) = zero
       FintNNLO_dy_ns(15:16)  = zero


    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(TCS4Lim%AmpMom,res_TCS4)
       call product_for_pdf_matrix(res_TCS4,ker,res_TCS42)


       !call get_respdf(ns_lumi,1,1,TCS4Lim,res_TCS42,respdf_ch1)
       call get_respdf_ns(1,1,TCS4Lim,res_TCS42,respdf_ch1)

       s24 = TCS4Lim%Lim_KinInv(1) 
       s25 = TCS4Lim%Lim_KinInv(2) 
       s45 = TCS4Lim%Lim_KinInv(3) 
       z4 =  TCS4Lim%Lim_KinInv(4) 
       z5 =  TCS4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(15) = (-one)*respdf_ch1(1) *  four /s25 * Pqq(z5,one) * two*Cf/s24/z4 * TCS4Lim%wgt*jac


       s24 = TCC4S4Lim%Lim_KinInv(1) 
       s25 = TCC4S4Lim%Lim_KinInv(2) 
       s45 = TCC4S4Lim%Lim_KinInv(3) 
       z4 =  TCC4S4Lim%Lim_KinInv(4) 
       z5 =  TCC4S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(16) = respdf_ch1(1)* four /s25 * Pqq(z5,one) * two*Cf/s24/z4 *TCC4S4Lim%wgt*jac*tag_ns_wga

       kin(9) = sum(FintNNLO_dy_ns(15:16))

    endif


    call fill_histo([kin(9)],vegasweight)

    ff(1) = sum(kin)
    call close_histo()

    ! print*, 'xx=', xx
    !      print  *, 1, FintNNLO_dy_ns(1), 'hard'
    !      print  *, '------------------------------------------------------'
    !      print  *, 2, FintNNLO_dy_ns(2)/FintNNLO_dy_ns(1), 'x1=0'
    !      print  *, 3, FintNNLO_dy_ns(4)/FintNNLO_dy_ns(1), 'x2=0'
    !      print  *, 4, FintNNLO_dy_ns(10)/FintNNLO_dy_ns(1), 'x3=0'
    !      print  *, 5, FintNNLO_dy_ns(5)/FintNNLO_dy_ns(1), 'x4=0'
    !      print  *, '------------------------------------------------------'
    !      print  *, 6, FintNNLO_dy_ns(6)/FintNNLO_dy_ns(1),    'x1=0,x2=0'
    !      print  *, 7, FintNNLO_dy_ns(15)/FintNNLO_dy_ns(1), 'x1=0,x3=0'     
    !      print  *, 8, FintNNLO_dy_ns(3)/FintNNLO_dy_ns(1), 'x1=0,x4=0'
    !      print  *, 9, FintNNLO_dy_ns(13)/FintNNLO_dy_ns(1), 'x2=0,x3=0'
    !      print  *, 10, FintNNLO_dy_ns(12)/FintNNLO_dy_ns(1), 'x2=0,x4=0'
    !      print  *, 11, FintNNLO_dy_ns(11)/FintNNLO_dy_ns(1), 'x3=0,x4=0'
    !      print  *, '------------------------------------------------------'     
    !      print  *, 12, FintNNLO_dy_ns(7)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0'
    !      print  *, 13, FintNNLO_dy_ns(8)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x4=0'
    !      print  *, 14, FintNNLO_dy_ns(16)/FintNNLO_dy_ns(1), 'x1=0,x3=0,x4=0'  
    !      print  *, 15, FintNNLO_dy_ns(14)/FintNNLO_dy_ns(1), 'x2=0,x3=0,x4=0'
    !      print  *, '------------------------------------------------------'     
    !      print  *, 16, FintNNLO_dy_ns(9)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0,x4=0'


    !    print*, 'ff(1)', ff(1)

    !   pause


    !#if(_identification_check == 1)  
    !    if( .not. check ) stop
    !#endif
    !
    !    if (ff(1) .ne. ff(1)) then
    !       !!print  *, 'nan'
    !       !!print  *, 'xx',xx
    !       !!print  *, ''
    !       do i = 1,16
    !          !!print  *, i, FintNNLO_dy_ns(i)
    !       enddo
    !       ff(1) = zero
    !    endif

    !#if (_DebugFeatures == 1)
    !    FintNNLO_CheckLim_dy_ns = FintNNLO_dy_ns
    !#endif

    !     do i = 1,16
    !        !!print  *, i,FintNNLO_dy_ns(i)
    !     enddo
    !     !!print  *, 'ff(1)', ff(1)
    !     pause

    !   pause


    !do i =1,9

    !print*, i, kin(i)

    !enddo

    !print *, ''

    if (ff(1) .ne. ff(1)) then
       print  *, 'nan'
       print  *, 'xx',xx
       print  *, ''
       do i = 1,16
          print  *, i, FintNNLO_dy_ns(i)
       enddo
       ff(1) = zero
    endif

#if(_withchecks == 1)
    FintNNLO = FintNNLO_dy_ns
#endif
    
    return
    
  end function xsect_nnlo_5262c_rr_z_ga



  function xsect_nnlo_5162_rr_z_ga(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_5162_rr_z_ga,ndim,ncomp,userdata,iternumber,i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,S4Lim,S4S5Lim,S5Lim,C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,&
         C5Lim,C5S4Lim,C5S4S5Lim,C5S5Lim,C4C5Lim,C4C5S4Lim,C4C5S4S5Lim,C4C5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR
    real(dp)    :: xi1,xi2,f1(-6:7),f2(-6:7),my_lumi,respdf(1),respdf_b(1)
    complex(dp) :: respdf_sping(-1:1,-1:1),respdf_sping2(-1:1,-1:1,-1:1,-1:1)
    !--


    complex(dp) :: hell_mat42(-1:1,-1:1),hell_mat51(-1:1,-1:1)
    !--



    real(dp)    :: respdf_qed(1),respdf_qcd(1)
    complex(dp) :: respdf_sping_qed(-1:1,-1:1),respdf_sping_qcd(-1:1,-1:1)
    real(dp)    :: respdf_ec,respdf_ce,respdf_b_qed(1),respdf_b_qcd(1)
    real(dp)    :: reseik_qed,reseik_qcd
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s24,s25,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_dy_ns(16),kin(1:9),skin(1:7)
    logical  :: makecut
    !--
    real(dp) :: mysymm = half
    real(dp) :: qcdmom(4,3),qedmom(4,3)
    !
    real(dp) :: wmom(4)
    real(dp) :: respdf_ch1(1),respdf_ch2
    !
    real(dp) :: resnnlo(2,2), res_S5S4(2,2), res2(2,2), eik(4), ker(2,2), res_S5_nlo(2,2)
    real(dp) :: res_a_nlo(2,2), res_C4S5(2,2), res_C4S52(2,2), res_C5(2,2), res_C52(2,2)
    real(dp) :: res_C4C5(2,2)
    real(dp) :: resnnlo_4q_tc(4),resnnlo_4q_fin(8)
    real(dp) :: res_g_nlo(12)
    real(dp) :: vec_a(4),vec_g(4),vec_ag(4),vec_ga(4)
    real(dp) :: vec_ag_a(4),vec_ag_g(4),vec_ga_a(4),vec_ga_g(4)
    real(dp) :: resvec(4)
    real(dp) :: res_ag(4),res_ga(4)
    real(dp) :: respdf_a(1),respdf_g(1)
    real(dp) :: respdf_ag(1),respdf_ga(1), tag_ns_wga


#if(_identification_check == 1)
    logical    :: check, check_temp
    !the debug arrays contains hard and single unresolved contirbutions
    ! ( qqb -> ga  , qqb -> ag  , qqb->qqb , a+S5g , g+S5a a+C5g , g+C5a )
    integer :: debug_binarray(9,nobs)
    real(dp) :: debug_weightarray(9,nobs)
    real(dp) :: debug_kin(9)
    debug_binarray    = zero
    debug_weightarray = zero
    debug_kin         = zero
    check = .true.
#endif

    tag_ns_wga = one 

    callcount = callcount + 1

    xsect_nnlo_5162_rr_z_ga = 0
    ff(1) = zero
    xx = zero

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)

    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buff_rr) return


    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call open_histo()

    call kinematics_nnlo_5162(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S4Lim,S5Lim,S4S5Lim,C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,&
         C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim)  

    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!

    !1-- regular piece ; gluon-photon

    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo(HardProc)

    if (makecut.or.HardProc%flag) then
       kin(1) = zero

       FintNNLO_dy_ns(1) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_ga_qqb(HardProc%AmpMom,resnnlo)

       !call get_respdf(ns_lumi,1,1,HardProc,resnnlo,respdf)
       call get_respdf_ns(1,1,HardProc,resnnlo,respdf)

       kin(1) = respdf(1)*HardProc%wgt*jac
       FintNNLO_dy_ns(1) = kin(1)


       call fill_histo([kin(1)],vegasweight)


    endif



    !! -----------------------------!!
    !!     S5 + S5/C5          !!
    !! -----------------------------!!

    !-- 2 -> S5   ; photon
    !-- 3 -> S5C5 ; photon

    xi1 = S5Lim%PartFrac(1) 
    xi2 = S5Lim%PartFrac(2) 
    print *, 'here'
call cut_histo(S5Lim)

    if (makecut.or.S5Lim%flag) then

       kin(2) = zero
       FintNNLO_dy_ns(2:3) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_g_qqb(S5Lim%AmpMom,res_S5_nlo)
       call get_Eik_ph_z_ga(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3), &
            S5Lim%LimMom(:,4), S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5_nlo,eik,res2)
       !call get_respdf(ns_lumi,1,1,S5Lim,res2,respdf_ga)
       call get_respdf_ns(1,1,S5Lim,res2,respdf_ga)

       FintNNLO_dy_ns(2) = (-one)*respdf_ga(1)*S5Lim%wgt*jac



       s15 = C5S5Lim%Lim_KinInv(2)
       z5 = C5S5Lim%Lim_KinInv(5)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2) 

       call product_for_pdf_matrix(res_S5_nlo,ker,res2)          
       !call get_respdf(ns_lumi,1,1,S5Lim,res2,respdf_qcd)
       call get_respdf_ns(1,1,S5Lim,res2,respdf_qcd)

       FintNNLO_dy_ns(3) = four/z5/s15 * respdf_qcd(1) * C5S5Lim%wgt*jac

       kin(2) = FintNNLO_dy_ns(2) + FintNNLO_dy_ns(3)


       call fill_histo([kin(2)],vegasweight)

    endif



    !-- 4 -> C4 ; 

    xi1 = C4Lim%PartFrac(1) 
    xi2 = C4Lim%PartFrac(2)
    print *, 'here'
call cut_histo(C4Lim)

    if (makecut.or.C4Lim%flag) then
       kin(3) = zero
       FintNNLO_dy_ns(4) = zero  

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)
       call res_tree_a_qqb(C4Lim%AmpMom,res_a_nlo)
       !call get_respdf(ns_lumi,1,1,C4Lim,res_a_nlo,respdf_b_qed)
       call get_respdf_ns(1,1,C4Lim,res_a_nlo,respdf_b_qed)

       s14 = C4Lim%Lim_KinInv(1) 
       z1  = C4Lim%Lim_KinInv(4)

       FintNNLO_dy_ns(4)  = (-one)*(-two/s14)* Pqq(z1,Cf) * respdf_b_qed(1) *C4Lim%wgt*jac


       kin(3) = FintNNLO_dy_ns(4) 
       call fill_histo([kin(3)],vegasweight)

    endif

    !-- 5 -> C5 ; photon

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)

    print *, 'here'
call cut_histo(C5Lim)

    if (makecut.or.C5Lim%flag) then
       kin(4) = zero
       FintNNLO_dy_ns(5)  = zero
    else


       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2) 

       call res_tree_g_qqb(C5Lim%AmpMom,res_C5)
       call product_for_pdf_matrix(res_C5,ker,res_C52)

       !call get_respdf(ns_lumi,1,1,C5Lim,res_C52,respdf_ch1)
       call get_respdf_ns(1,1,C5Lim,res_C52,respdf_ch1)

       s25 = C5Lim%Lim_KinInv(2) 
       z2  = C5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(5)  = (-one)*(-two/s25)* Pqq(z2,one)*( respdf_ch1(1) ) * C5Lim%wgt*jac

       kin(4) = FintNNLO_dy_ns(5) 
       call fill_histo([kin(5)],vegasweight)

    endif



    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!

    !-- 6 -> S4S5
    !-- 7 -> C4DS
    !-- 8 -> C5S4S5
    !-- 9 -> C4C5S4S5

    xi1 = S4S5Lim%PartFrac(1) 
    xi2 = S4S5Lim%PartFrac(2)

    print *, 'here'
call cut_histo(S4S5Lim)

    if (makecut.or.S4S5Lim%flag) then
       kin(6) = zero
       FintNNLO_dy_ns(6:9)  = zero


    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(S4S5Lim%AmpMom,res_S5S4) 

       call get_Eik_ds_z_ga(S4S5Lim%LimMom(:,1),S4S5Lim%LimMom(:,2),S4S5Lim%LimMom(:,3),&
            S4S5Lim%LimMom(:,4),S4S5Lim%LimMom(:,5),S4S5Lim%LimMom(:,6),eik)
       call product_for_pdf(res_S5S4,eik,res2)
       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ag)
       call get_respdf_ns(1,1,S4S5Lim,res2,respdf_ag)

       FintNNLO_dy_ns(6) = (-one)*(respdf_ag(1))*S4S5Lim%wgt*jac



       s12 = C4S4S5Lim%Lim_KinInv(1) 
       s25 = C4S4S5Lim%Lim_KinInv(2) 
       s15 = C4S4S5Lim%Lim_KinInv(3) 
       s14 = C4S4S5Lim%Lim_KinInv(4) 
       z4  = C4S4S5Lim%Lim_KinInv(5) 

       call get_Eik_ph_z_ga(C4S4S5Lim%LimMom(:,1),C4S4S5Lim%LimMom(:,2),C4S4S5Lim%LimMom(:,3), &
            C4S4S5Lim%LimMom(:,4), C4S4S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5S4,eik,res2)

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ga)


       FintNNLO_dy_ns(7) =(-one)*(respdf_ga(1))*four*Cf/s14/z4*C4S4S5Lim%wgt*jac



       !C5S4S5
       s12 = C5S4S5Lim%Lim_KinInv(1) 
       s24 = C5S4S5Lim%Lim_KinInv(2) 
       s14 = C5S4S5Lim%Lim_KinInv(3) 
       s25 = C5S4S5Lim%Lim_KinInv(4) 
       z5  = C5S4S5Lim%Lim_KinInv(5)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2)

       call product_for_pdf_matrix(res_S5S4,ker,res2)
       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf)

       FintNNLO_dy_ns(8) = (-one)*respdf(1) * four/s25/z5 * Cf*four*s12/s14/s24*C5S4S5Lim%wgt*jac

       !C4C5S4S5

       s14 = C4C5S4S5Lim%Lim_KinInv(1) 
       s25 = C4C5S4S5Lim%Lim_KinInv(2) 
       z5 =  C4C5S4S5Lim%Lim_KinInv(4) 
       z4 =  C4C5S4S5Lim%Lim_KinInv(5) 

       FintNNLO_dy_ns(9) =  (one)*16.0_dp*Cf/s14/z4/s25/z5 &
            * respdf(1)* C4C5S4S5Lim%wgt*jac


       kin(5) = sum(FintNNLO_dy_ns(6:9))


       call fill_histo([kin(5)],vegasweight)

    endif



    !-- 10 -> C4S5
    !-- 11 -> C4C5S5

    xi1 = C4S5Lim%PartFrac(1) 
    xi2 = C4S5Lim%PartFrac(2)

    print *, 'here'
call cut_histo(C4S5Lim)

    if (makecut.or.C4S5Lim%flag) then

       kin(6) = zero
       FintNNLO_dy_ns(10:11)  = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C4S5Lim%AmpMom,res_C4S5)         
       call get_Eik_ph_z_ga(C4S5Lim%LimMom(:,1),C4S5Lim%LimMom(:,2),C4S5Lim%LimMom(:,3), &
            C4S5Lim%LimMom(:,4),C4S5Lim%LimMom(:,5),eik)

       call product_for_pdf(res_C4S5,eik,res_C4S52)


       !call get_respdf(ns_lumi,1,1,C4S5Lim,res_C4S52,respdf)          

       s14 = C4S5Lim%Lim_KinInv(4) 
       z4 =  C4S5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(10) = (- respdf(1))*two/s14*Pqq(z4,Cf)* C4S5Lim%wgt*jac



       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2) 

       call product_for_pdf_matrix(res_C4S5,ker,res2)          

       !call get_respdf(ns_lumi,1,1,C4S5Lim,res2,respdf_ch1)

       s14 = C4C5S5Lim%Lim_KinInv(1) 
       s25 = C4C5S5Lim%Lim_KinInv(2) 
       z4 =  C4C5S5Lim%Lim_KinInv(4) 
       z5 =  C4C5S5Lim%Lim_KinInv(5) 
       !

       FintNNLO_dy_ns(11) = 8.0_dp/s14/s25/z5 * Pqq(z4,Cf )  *respdf_ch1(1) * C4C5S5Lim%wgt*jac


       kin(6) = sum(FintNNLO_dy_ns(10:11))

       call fill_histo([kin(6)],vegasweight)

    endif


    !-- 12 -> C4C5

    xi1 = C4C5Lim%PartFrac(1) 
    xi2 = C4C5Lim%PartFrac(2) 
    print *, 'here'
call cut_histo(C4C5Lim)
    if (makecut.or.C4C5Lim%flag) then

       FintNNLO_dy_ns(12) = zero
       kin(7) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C4C5Lim%AmpMom,res_C4C5)  

       s14 = C4C5Lim%Lim_KinInv(1) 
       s25 = C4C5Lim%Lim_KinInv(2) 
       z4 =  C4C5Lim%Lim_KinInv(4) 
       z5 =  C4C5Lim%Lim_KinInv(5)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2)

       call product_for_pdf_matrix(res_C4C5,ker,res2)
       !call get_respdf(ns_lumi,1,1,C4C5Lim,res2,respdf)


       FintNNLO_dy_ns(12) = four/s25/s14 * Pqq(z5,esq) * Pqq(z4,Cf )  *respdf(1) * C4C5Lim%wgt*jac


       kin(7) = FintNNLO_dy_ns(12)

       call fill_histo([kin(7)],vegasweight)

    endif



    !! -----------------------------!!
    !!     S4 + S4/C4         !!
    !! -----------------------------!!

    !--  -> S4   ; gluon
    !--  -> S4C4 ; gluon

    xi1 = S4Lim%PartFrac(1) 
    xi2 = S4Lim%PartFrac(2) 
    print *, 'here'
call cut_histo(S4Lim)


    if (makecut.or.S4Lim%flag) then
       kin(8) = zero
       FintNNLO_dy_ns(13) = zero
       FintNNLO_dy_ns(14) = zero
    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)
       call res_tree_a_qqb(S4Lim%AmpMom,res_a_nlo)

       !call get_respdf(ns_lumi,1,1,S4Lim,res_a_nlo,respdf)

       reseik_qcd = myEik_g(S4Lim%LimMom(:,1),S4Lim%LimMom(:,2),S4Lim%LimMom(:,3),Cf)

       FintNNLO_dy_ns(13) = (-one)*respdf(1) *reseik_qcd*S4Lim%wgt*jac

       s14 = C4S4Lim%Lim_KinInv(2)
       z4 = C4S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(14) = four*Cf/z4/s14 * respdf(1) * C4S4Lim%wgt*jac


       kin(8) = FintNNLO_dy_ns(13) + FintNNLO_dy_ns(14)

       call fill_histo([kin(8)],vegasweight)


    endif


    ! - - - - - - C5S4 and  C4C5S4 ; photon

    xi1 = C5S4Lim%PartFrac(1) 
    xi2 = C5S4Lim%PartFrac(2) 
    print *, 'here'
call cut_histo(C5S4Lim)

    if (makecut.or.C5S4Lim%flag) then

       kin(9) = zero
       FintNNLO_dy_ns(15) = zero
       FintNNLO_dy_ns(16) = zero 

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2) 

       call res_tree_qqb(C5S4Lim%AmpMom,res_C5)
       call product_for_pdf_matrix(res_C5,ker,res_C52)

       !call get_respdf(ns_lumi,1,1,C5S4Lim,res_C52,respdf_ch1)


       s12 = C5S4Lim%Lim_KinInv(1)
       s14 = C5S4Lim%Lim_KinInv(2)
       s24 = C5S4Lim%Lim_KinInv(3) 
       s25 = C5S4Lim%Lim_KinInv(4) 
       z5  = C5S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(15)  = (-one)*respdf_ch1(1)*two*Cf*s12/s14/s24*four/s25*Pqq(z5,one) * C5S4Lim%wgt*jac


       z4 = C4C5S4Lim%Lim_KinInv(1)
       s14 = C4C5S4Lim%Lim_KinInv(2)
       s24 = C4C5S4Lim%Lim_KinInv(3) 
       s15 = C4C5S4Lim%Lim_KinInv(4) 
       z1  = C4C5S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(16)  = respdf_ch1(1)*four*Cf/s24/z4*two/s15*Pqq(z1,one) * C4C5S4Lim%wgt*jac*tag_ns_wga


       kin(9) = FintNNLO_dy_ns(15)+FintNNLO_dy_ns(16)

    endif

    call fill_histo([kin(9)],vegasweight)

    ff(1) = sum(kin)
    call close_histo()

    !  print*, 'xx=', xx    
    !     print  *, 1, FintNNLO_dy_ns(1), 'hard'
    !     print  *, '------------------------------------------------------'
    !     print  *, 2, FintNNLO_dy_ns(13)/FintNNLO_dy_ns(1), 'x1=0'
    !     print  *, 3, FintNNLO_dy_ns(2)/FintNNLO_dy_ns(1), 'x2=0'
    !     print  *, 4, FintNNLO_dy_ns(4)/FintNNLO_dy_ns(1), 'x3=0'
    !     print  *, 5, FintNNLO_dy_ns(5)/FintNNLO_dy_ns(1), 'x4=0'
    !     print  *, '------------------------------------------------------'
    !     print  *, 6, FintNNLO_dy_ns(6)/FintNNLO_dy_ns(1),    'x1=0,x2=0'
    !     print  *, 7, FintNNLO_dy_ns(14)/FintNNLO_dy_ns(1), 'x1=0,x3=0'     
    !     print  *, 8, FintNNLO_dy_ns(15)/FintNNLO_dy_ns(1), 'x1=0,x4=0'
    !     print  *, 9, FintNNLO_dy_ns(10)/FintNNLO_dy_ns(1), 'x2=0, x3=0'
    !     print  *, 10, FintNNLO_dy_ns(3)/FintNNLO_dy_ns(1), 'x2=0,x4=0'
    !     print  *, 11, FintNNLO_dy_ns(12)/FintNNLO_dy_ns(1), 'x3=0,x4=0'
    !     print  *, '------------------------------------------------------'     
    !     print  *, 12, FintNNLO_dy_ns(7)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0'
    !     print  *, 13, FintNNLO_dy_ns(8)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x4=0'
    !     print  *, 14, FintNNLO_dy_ns(16)/FintNNLO_dy_ns(1), 'x1=0,x3=0,x4=0'  
    !     print  *, 15, FintNNLO_dy_ns(11)/FintNNLO_dy_ns(1), 'x2=0,x3=0,x4=0'
    !     print  *, '------------------------------------------------------'     
    !     print  *, 16, FintNNLO_dy_ns(9)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0,x4=0'

    ! !!$    print *, ''
    ! !!$    print*, 'kin(1)', kin(1)
    ! !!$    print*, 'kin(2)', kin(2)
    ! !!$    print*, 'kin(3)', kin(3)
    ! !!$    print*, 'kin(4)', kin(4)
    ! !!$    print*, 'kin(5)', kin(5)
    ! !!$    print*, 'kin(6)', kin(6)
    ! !!$    print*, 'kin(7)', kin(7)
    ! !!$    print*, 'kin(8)', kin(8)
    ! !!$    print*, 'kin(9)', kin(9)
    ! !!$    print*

    !    print*, 'ff(1)', ff(1)
    !        pause
    ! !!$
!!$
!!$
!!$
!!$
!!$    
!!$    if (ff(1) .ne. ff(1)) then
!!$       !print  *, 'nan'
!!$       !print  *, 'xx',xx
!!$       !print  *, ''
!!$       do i = 1,16
!!$          !print  *, i, FintNNLO_dy_ns(i)
!!$       enddo
!!$       ff(1) = zero
!!$    endif
!!$#if(_identification_check == 1)  
!!$    if( .not. check ) stop
!!$#endif
!!$
!!$#if (_DebugFeatures == 1)
!!$    FintNNLO_CheckLim_dy_ns(1) = FintNNLO_dy_ns(1)
!!$#endif
!!$#if (_DebugFeatures == 1)
!!$    FintNNLO_CheckLim_dy_ns(4:16) = FintNNLO_dy_ns(4:16)
!!$#endif
!!$
!!$    ! do i = 1,16
!!$    !    !print  *, i,FintNNLO_dy_ns(i)
!!$    ! enddo
!!$    ! !print  *, 'ff(1)', ff(1)
!!$    ! pause


    if (ff(1) .ne. ff(1)) then
       print  *, 'nan'
       print  *, 'xx',xx
       print  *, ''
       do i = 1,16
          print  *, i, FintNNLO_dy_ns(i)
       enddo
       ff(1)=zero
    endif

#if(_withchecks == 1)
    FintNNLO = FintNNLO_dy_ns
#endif

    return

  end function xsect_nnlo_5162_rr_z_ga






  function xsect_nnlo_5261_rr_z_ga(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_5261_rr_z_ga,ndim,ncomp,userdata,iternumber,i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig)  :: HardProc,S4Lim,S4S5Lim,S5Lim,C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,&
         C5Lim,C5S4Lim,C5S4S5Lim,C5S5Lim,C4C5Lim,C4C5S4Lim,C4C5S4S5Lim,C4C5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR
    real(dp)    :: xi1,xi2,f1(-6:7),f2(-6:7),my_lumi,respdf(1),respdf_b(1)
    complex(dp) :: respdf_sping(-1:1,-1:1),respdf_sping2(-1:1,-1:1,-1:1,-1:1)
    !--


    complex(dp) :: hell_mat42(-1:1,-1:1),hell_mat51(-1:1,-1:1)
    !--



    real(dp)    :: respdf_qed(1),respdf_qcd(1)
    complex(dp) :: respdf_sping_qed(-1:1,-1:1),respdf_sping_qcd(-1:1,-1:1)
    real(dp)    :: respdf_ec,respdf_ce,respdf_b_qed(1),respdf_b_qcd(1)
    real(dp)    :: reseik_qed,reseik_qcd
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s24,s25,reseik,reseik_ds,reseik_dss5
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_dy_ns(16),kin(1:9),skin(1:7)
    logical  :: makecut
    !--
    real(dp) :: mysymm = half
    real(dp) :: qcdmom(4,3),qedmom(4,3)
    !
    real(dp) :: wmom(4)
    real(dp) :: respdf_ch1(1),respdf_ch2
    !
    real(dp) :: resnnlo(2,2), res_S5S4(2,2), res2(2,2), eik(4), ker(2,2), res_S5_nlo(2,2)
    real(dp) :: res_a_nlo(2,2), res_C4S5(2,2), res_C4S52(2,2), res_C5(2,2), res_C52(2,2)
    real(dp) :: res_C4C5(2,2)
    real(dp) :: resnnlo_4q_tc(4),resnnlo_4q_fin(8)
    real(dp) :: res_g_nlo(12)
    real(dp) :: vec_a(4),vec_g(4),vec_ag(4),vec_ga(4)
    real(dp) :: vec_ag_a(4),vec_ag_g(4),vec_ga_a(4),vec_ga_g(4)
    real(dp) :: resvec(4)
    real(dp) :: res_ag(4),res_ga(4)
    real(dp) :: respdf_a(1),respdf_g(1)
    real(dp) :: respdf_ag(1),respdf_ga(1), tag_ns_wga

#if(_identification_check == 1)
    logical    :: check, check_temp
    !the debug arrays contains hard and single unresolved contirbutions
    ! ( qqb -> ga  , qqb -> ag  , qqb->qqb , a+S5g , g+S5a a+C5g , g+C5a )
    integer :: debug_binarray(9,nobs)
    real(dp) :: debug_weightarray(9,nobs)
    real(dp) :: debug_kin(9)
    debug_binarray    = zero
    debug_weightarray = zero
    debug_kin         = zero
    check = .true.
#endif

    tag_ns_wga = one 

    callcount = callcount + 1

    xsect_nnlo_5261_rr_z_ga = 0
    ff(1) = zero
    xx = zero


    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buff_rr) return 
    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call open_histo()

    call kinematics_nnlo_5261(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S4Lim,S5Lim,S4S5Lim,C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,&
         C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim)  

    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!

    !1-- regular piece ; gluon-photon

    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)
    call cut_histo(HardProc)

    if (makecut.or.HardProc%flag) then
       kin(1) = zero

       FintNNLO_dy_ns(1) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_ga_qqb(HardProc%AmpMom,resnnlo)

       !call get_respdf(ns_lumi,1,1,HardProc,resnnlo,respdf)

       kin(1) = respdf(1)*HardProc%wgt*jac
       FintNNLO_dy_ns(1) = kin(1)


       call fill_histo([kin(1)],vegasweight)


    endif

    !! -----------------------------!!
    !!     S5 + S5/C5          !!
    !! -----------------------------!!

    !-- 2 -> S5   ; photon
    !-- 3 -> S5C5 ; photon

    xi1 = S5Lim%PartFrac(1) 
    xi2 = S5Lim%PartFrac(2) 
    call cut_histo(S5Lim)

    if (makecut.or.S5Lim%flag) then

       kin(2) = zero
       FintNNLO_dy_ns(2:3) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_g_qqb(S5Lim%AmpMom,res_S5_nlo)
       call get_Eik_ph_z_ga(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3), &
            S5Lim%LimMom(:,4), S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5_nlo,eik,res2)
       !call get_respdf(ns_lumi,1,1,S5Lim,res2,respdf_ga)

       FintNNLO_dy_ns(2) = (-one)*respdf_ga(1)*S5Lim%wgt*jac


       s15 = C5S5Lim%Lim_KinInv(2)
       z5 = C5S5Lim%Lim_KinInv(5)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2) 

       call product_for_pdf_matrix(res_S5_nlo,ker,res2)          
       !call get_respdf(ns_lumi,1,1,S5Lim,res2,respdf_qcd)

       FintNNLO_dy_ns(3) = four/z5/s15 * respdf_qcd(1) * C5S5Lim%wgt*jac

       kin(2) = FintNNLO_dy_ns(2) + FintNNLO_dy_ns(3)

       call fill_histo([kin(2)],vegasweight)

    endif


    !-- 4 -> C4 ; 

    xi1 = C4Lim%PartFrac(1) 
    xi2 = C4Lim%PartFrac(2)
    call cut_histo(C4Lim)

    if (makecut.or.C4Lim%flag) then
       kin(3) = zero
       FintNNLO_dy_ns(4) = zero  

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_a_qqb(C4Lim%AmpMom,res_a_nlo)
       !call get_respdf(ns_lumi,1,1,C4Lim,res_a_nlo,respdf_b_qed)

       s24 = C4Lim%Lim_KinInv(1) 
       z2  = C4Lim%Lim_KinInv(4)

       FintNNLO_dy_ns(4)  = (-one)*(-two/s24)* Pqq(z2,Cf) * respdf_b_qed(1) *C4Lim%wgt*jac


       kin(3) = FintNNLO_dy_ns(4) 
       call fill_histo([kin(3)],vegasweight)

    endif


    !-- 5 -> C5 ; photon

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)

    call cut_histo(C5Lim)

    if (makecut.or.C5Lim%flag) then
       kin(4) = zero
       FintNNLO_dy_ns(5)  = zero
    else


       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2) 

       call res_tree_g_qqb(C5Lim%AmpMom,res_C5)
       call product_for_pdf_matrix(res_C5,ker,res_C52)

       !call get_respdf(ns_lumi,1,1,C5Lim,res_C52,respdf_ch1)

       s15 = C5Lim%Lim_KinInv(2) 
       z1  = C5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(5)  = (-one)*(-two/s15)*( respdf_ch1(1) * Pqq(z1,one)) * C5Lim%wgt*jac

       kin(4) = FintNNLO_dy_ns(5) 
       call fill_histo([kin(5)],vegasweight)

    endif



    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!

    !-- 6 -> S4S5
    !-- 7 -> C4DS
    !-- 8 -> C5S4S5
    !-- 9 -> C4C5S4S5

    xi1 = S4S5Lim%PartFrac(1) 
    xi2 = S4S5Lim%PartFrac(2)

    call cut_histo(S4S5Lim)

    if (makecut.or.S4S5Lim%flag) then
       kin(5) = zero
       FintNNLO_dy_ns(6:9)  = zero


    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(S4S5Lim%AmpMom,res_S5S4) 

       call get_Eik_ds_z_ga(S4S5Lim%LimMom(:,1),S4S5Lim%LimMom(:,2),S4S5Lim%LimMom(:,3),&
            S4S5Lim%LimMom(:,4),S4S5Lim%LimMom(:,5),S4S5Lim%LimMom(:,6),eik)
       call product_for_pdf(res_S5S4,eik,res2)
       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ag)

       FintNNLO_dy_ns(6) = (-one)*(respdf_ag(1))*S4S5Lim%wgt*jac

       s12 = C4S4S5Lim%Lim_KinInv(1) 
       s25 = C4S4S5Lim%Lim_KinInv(2) 
       s15 = C4S4S5Lim%Lim_KinInv(3) 
       s14 = C4S4S5Lim%Lim_KinInv(4) 
       z4  = C4S4S5Lim%Lim_KinInv(5) 

       call get_Eik_ph_z_ga(C4S4S5Lim%LimMom(:,1),C4S4S5Lim%LimMom(:,2),C4S4S5Lim%LimMom(:,3), &
            C4S4S5Lim%LimMom(:,4), C4S4S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5S4,eik,res2)

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ga)


       FintNNLO_dy_ns(7) =(-one)*(respdf_ga(1))*four*Cf/s14/z4*C4S4S5Lim%wgt*jac*tag_ns_wga
       !C5S4S5
       s12 = C5S4S5Lim%Lim_KinInv(1) 
       s24 = C5S4S5Lim%Lim_KinInv(2) 
       s14 = C5S4S5Lim%Lim_KinInv(3) 
       s25 = C5S4S5Lim%Lim_KinInv(4) 
       z5  = C5S4S5Lim%Lim_KinInv(5)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2)

       call product_for_pdf_matrix(res_S5S4,ker,res2)
       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf)

       FintNNLO_dy_ns(8) = (-one)*respdf(1) * four/s25/z5 * Cf*four*s12/s14/s24*C5S4S5Lim%wgt*jac


       !C4C5S4S5

       s14 = C4C5S4S5Lim%Lim_KinInv(1) 
       s25 = C4C5S4S5Lim%Lim_KinInv(2) 
       z4 =  C4C5S4S5Lim%Lim_KinInv(4) 
       z5 =  C4C5S4S5Lim%Lim_KinInv(5) 

       FintNNLO_dy_ns(9) =  (one)*16.0_dp*Cf/s14/z4/s25/z5 &
            * respdf(1)* C4C5S4S5Lim%wgt*jac


       kin(5) = sum(FintNNLO_dy_ns(6:9))

       call fill_histo([kin(5)],vegasweight)

    endif



    !-- 10 -> C4S5
    !-- 11 -> C4C5S5

    xi1 = C4S5Lim%PartFrac(1) 
    xi2 = C4S5Lim%PartFrac(2)

    call cut_histo(C4S5Lim)

    if (makecut.or.C4S5Lim%flag) then

       kin(6) = zero
       FintNNLO_dy_ns(10:11)  = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C4S5Lim%AmpMom,res_C4S5)         
       call get_Eik_ph_z_ga(C4S5Lim%LimMom(:,1),C4S5Lim%LimMom(:,2),C4S5Lim%LimMom(:,3), &
            C4S5Lim%LimMom(:,4),C4S5Lim%LimMom(:,5),eik)

       call product_for_pdf(res_C4S5,eik,res_C4S52)


       !call get_respdf(ns_lumi,1,1,C4S5Lim,res_C4S52,respdf)          

       s24 = C4S5Lim%Lim_KinInv(4) 
       z4 =  C4S5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(10) = (- respdf(1))*two/s24*Pqq(z4,Cf)* C4S5Lim%wgt*jac


       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2) 

       call product_for_pdf_matrix(res_C4S5,ker,res2)          

       !call get_respdf(ns_lumi,1,1,C4S5Lim,res2,respdf_ch1)

       s24 = C4C5S5Lim%Lim_KinInv(1) 
       s15 = C4C5S5Lim%Lim_KinInv(2) 
       z4 =  C4C5S5Lim%Lim_KinInv(4) 
       z5 =  C4C5S5Lim%Lim_KinInv(5) 
       !

       FintNNLO_dy_ns(11) = 8.0_dp/s15/s24/z5 * Pqq(z4,Cf )  *respdf_ch1(1) * C4C5S5Lim%wgt*jac

       kin(6) = sum(FintNNLO_dy_ns(10:11))

       call fill_histo([kin(6)],vegasweight)

    endif


    !-- 12 -> C4C5

    xi1 = C4C5Lim%PartFrac(1) 
    xi2 = C4C5Lim%PartFrac(2) 
    call cut_histo(C4C5Lim)
    if (makecut.or.C4C5Lim%flag) then

       FintNNLO_dy_ns(12) = zero
       kin(7) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C4C5Lim%AmpMom,res_C4C5)  

       s24 = C4C5Lim%Lim_KinInv(1) 
       s15 = C4C5Lim%Lim_KinInv(2) 
       z4 =  C4C5Lim%Lim_KinInv(4) 
       z5 =  C4C5Lim%Lim_KinInv(5)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2)

       call product_for_pdf_matrix(res_C4C5,ker,res2)
       !call get_respdf(ns_lumi,1,1,C4C5Lim,res2,respdf)


       FintNNLO_dy_ns(12) = four/s24/s15 * Pqq(z5,esq) * Pqq(z4,Cf )  *respdf(1) * C4C5Lim%wgt*jac


       kin(7) = FintNNLO_dy_ns(12)

       call fill_histo([kin(7)],vegasweight)

    endif



    !! -----------------------------!!
    !!     S4 + S4/C4         !!
    !! -----------------------------!!

    !here here here here

    !--  -> S4   ; gluon
    !--  -> S4C4 ; gluon

    xi1 = S4Lim%PartFrac(1) 
    xi2 = S4Lim%PartFrac(2) 
    call cut_histo(S4Lim)


    if (makecut.or.S4Lim%flag) then
       kin(8) = zero
       FintNNLO_dy_ns(13) = zero
       FintNNLO_dy_ns(14) = zero
    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)
       call res_tree_a_qqb(S4Lim%AmpMom,res_a_nlo)

       !call get_respdf(ns_lumi,1,1,S4Lim,res_a_nlo,respdf)

       reseik_qcd = myEik_g(S4Lim%LimMom(:,1),S4Lim%LimMom(:,2),S4Lim%LimMom(:,3),Cf)

       FintNNLO_dy_ns(13) = (-one)*respdf(1) *reseik_qcd*S4Lim%wgt*jac

       s24 = C4S4Lim%Lim_KinInv(2)
       z4 = C4S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(14) = four*Cf/z4/s24 * respdf(1) * C4S4Lim%wgt*jac

       kin(8) = FintNNLO_dy_ns(13) + FintNNLO_dy_ns(14)

       call fill_histo([kin(8)],vegasweight)


    endif


    ! - - - - - - C5S4 and  C4C5S4 ; photon

    xi1 = C5S4Lim%PartFrac(1) 
    xi2 = C5S4Lim%PartFrac(2) 
    call cut_histo(C5S4Lim)

    if (makecut.or.C5S4Lim%flag) then

       kin(9) = zero
       FintNNLO_dy_ns(15) = zero
       FintNNLO_dy_ns(16) = zero 

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       ker(1,1) = Qdn**2
       ker(1,2) = Qup**2
       ker(2,1) = ker(1,1)
       ker(2,2) = ker(1,2) 

       call res_tree_qqb(C5S4Lim%AmpMom,res_C5)
       call product_for_pdf_matrix(res_C5,ker,res_C52)

       !call get_respdf(ns_lumi,1,1,C5S4Lim,res_C52,respdf_ch1)


       s12 = C5S4Lim%Lim_KinInv(1)
       s14 = C5S4Lim%Lim_KinInv(2)
       s24 = C5S4Lim%Lim_KinInv(3) 
       s15 = C5S4Lim%Lim_KinInv(4) 
       z1  = C5S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(15)  = (-one)*respdf_ch1(1)*four*Cf*s12/s14/s24*two/s15*Pqq(z1,one) * C5S4Lim%wgt*jac

       z4 = C4C5S4Lim%Lim_KinInv(1)
       s14 = C4C5S4Lim%Lim_KinInv(2)
       s24 = C4C5S4Lim%Lim_KinInv(3) 
       s15 = C4C5S4Lim%Lim_KinInv(4) 
       z1  = C4C5S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(16)  = respdf_ch1(1)*four*Cf/s24/z4*two/s15*Pqq(z1,one) * C4C5S4Lim%wgt*jac*tag_ns_wga

       kin(9) = FintNNLO_dy_ns(15)+FintNNLO_dy_ns(16)


    endif

    call fill_histo([kin(9)],vegasweight)

    ff(1) = sum(kin)
    call close_histo()
    
    !  print*, 'xx=', xx
    !     print  *, 1, FintNNLO_dy_ns(1), 'hard'
    !     print  *, '------------------------------------------------------'
    !     print  *, 2, FintNNLO_dy_ns(13)/FintNNLO_dy_ns(1), 'x1=0'
    !     print  *, 3, FintNNLO_dy_ns(2)/FintNNLO_dy_ns(1), 'x2=0'
    !     print  *, 4, FintNNLO_dy_ns(4)/FintNNLO_dy_ns(1), 'x3=0'
    !     print  *, 5, FintNNLO_dy_ns(5)/FintNNLO_dy_ns(1), 'x4=0'
    !     print  *, '------------------------------------------------------'
    !     print  *, 6, FintNNLO_dy_ns(6)/FintNNLO_dy_ns(1),    'x1=0,x2=0'
    !     print  *, 7, FintNNLO_dy_ns(14)/FintNNLO_dy_ns(1), 'x1=0,x3=0'     
    !     print  *, 8, FintNNLO_dy_ns(15)/FintNNLO_dy_ns(1), 'x1=0,x4=0'
    !     print  *, 9, FintNNLO_dy_ns(10)/FintNNLO_dy_ns(1), 'x2=0, x3=0'
    !     print  *, 10, FintNNLO_dy_ns(3)/FintNNLO_dy_ns(1), 'x2=0,x4=0'
    !     print  *, 11, FintNNLO_dy_ns(12)/FintNNLO_dy_ns(1), 'x3=0,x4=0'
    !     print  *, '------------------------------------------------------'     
    !     print  *, 12, FintNNLO_dy_ns(7)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0'
    !     print  *, 13, FintNNLO_dy_ns(8)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x4=0'
    !     print  *, 14, FintNNLO_dy_ns(16)/FintNNLO_dy_ns(1), 'x1=0,x3=0,x4=0'  
    !     print  *, 15, FintNNLO_dy_ns(11)/FintNNLO_dy_ns(1), 'x2=0,x3=0,x4=0'
    !     print  *, '------------------------------------------------------'     
    !     print  *, 16, FintNNLO_dy_ns(9)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0,x4=0'



    ! !!$    print *, ''
    ! !!$    print*, 'kin(1)', kin(1)
    ! !!$    print*, 'kin(2)', kin(2)
    ! !!$    print*, 'kin(3)', kin(3)
    ! !!$    print*, 'kin(4)', kin(4)
    ! !!$    print*, 'kin(5)', kin(5)
    ! !!$    print*, 'kin(6)', kin(6)
    ! !!$    print*, 'kin(7)', kin(7)
    ! !!$    print*, 'kin(8)', kin(8)
    ! !!$    print*, 'kin(9)', kin(9)
    ! !!$    print*
    ! !!$    
    !     print*, 'ff(1)', ff(1)
    !     pause
    ! !!$
!!$
!!$
!!$
!!$
!!$    
!!$    if (ff(1) .ne. ff(1)) then
!!$       !print  *, 'nan'
!!$       !print  *, 'xx',xx
!!$       !print  *, ''
!!$       do i = 1,16
!!$          !print  *, i, FintNNLO_dy_ns(i)
!!$       enddo
!!$       ff(1) = zero
!!$    endif
!!$#if(_identification_check == 1)  
!!$    if( .not. check ) stop
!!$#endif
!!$
!!$#if (_DebugFeatures == 1)
!!$    FintNNLO_CheckLim_dy_ns(1) = FintNNLO_dy_ns(1)
!!$#endif
!!$#if (_DebugFeatures == 1)
!!$    FintNNLO_CheckLim_dy_ns(4:16) = FintNNLO_dy_ns(4:16)
!!$#endif
!!$
!!$    ! do i = 1,16
!!$    !    !print  *, i,FintNNLO_dy_ns(i)
!!$    ! enddo
!!$    ! !print  *, 'ff(1)', ff(1)
!!$    ! pause


    if (ff(1) .ne. ff(1)) then
       print  *, 'nan'
       print  *, 'xx',xx
       print  *, ''
       do i = 1,16
          print  *, i, FintNNLO_dy_ns(i)
       enddo
       ff(1)=zero
    endif
    
#if(_withchecks == 1)
    FintNNLO = FintNNLO_dy_ns
#endif


    return

  end function xsect_nnlo_5261_rr_z_ga







  function xsect_nnlo_5163_rr_z_ga(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_5163_rr_z_ga,ndim,ncomp,userdata,iternumber,i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,S4Lim,S4S5Lim,S5Lim,C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,&
         C5Lim,C5S4Lim,C5S4S5Lim,C5S5Lim,C4C5Lim,C4C5S4Lim,C4C5S4S5Lim,C4C5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR
    real(dp)    :: xi1,xi2,f1(-6:7),f2(-6:7),my_lumi,respdf(ipdf),respdf_b(1)
    complex(dp) :: respdf_sping(-1:1,-1:1),respdf_sping2(-1:1,-1:1,-1:1,-1:1)
    !--


    complex(dp) :: hell_mat42(-1:1,-1:1),hell_mat51(-1:1,-1:1)
    !--



    real(dp)    :: respdf_qed(1),respdf_qcd(1)
    complex(dp) :: respdf_sping_qed(-1:1,-1:1),respdf_sping_qcd(-1:1,-1:1)
    real(dp)    :: respdf_ec,respdf_ce,respdf_b_qed(1),respdf_b_qcd(1)
    real(dp)    :: reseik_qed,reseik_qcd
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s24,s25,reseik,reseik_ds,reseik_dss5, s65
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_dy_ns(16),kin(1:9),skin(1:7)
    logical  :: makecut
    !--
    real(dp) :: mysymm = half
    real(dp) :: qcdmom(4,3),qedmom(4,3)
    !
    real(dp) :: wmom(4)
    real(dp) :: respdf_ch1(1),respdf_ch2
    !
    real(dp) :: resnnlo(2,2), res_S5S4(2,2), res2(2,2), eik(4), ker(2,2), res_S5_nlo(2,2)
    real(dp) :: res_a_nlo(2,2), res_C4S5(2,2), res_C4S52(2,2), res_C5(2,2), res_C52(2,2)
    real(dp) :: res_C4C5(2,2)
    real(dp) :: resnnlo_4q_tc(4),resnnlo_4q_fin(8)
    real(dp) :: res_g_nlo(12)
    real(dp) :: vec_a(4),vec_g(4),vec_ag(4),vec_ga(4)
    real(dp) :: vec_ag_a(4),vec_ag_g(4),vec_ga_a(4),vec_ga_g(4)
    real(dp) :: resvec(4)
    real(dp) :: res_ag(4),res_ga(4)
    real(dp) :: respdf_a(1),respdf_g(1)
    real(dp) :: respdf_ag(1),respdf_ga(1), tag_ns_wga

    real(dp) :: pdec(4,3)

#if(_identification_check == 1)
    logical    :: check, check_temp
    !the debug arrays contains hard and single unresolved contirbutions
    ! ( qqb -> ga  , qqb -> ag  , qqb->qqb , a+S5g , g+S5a a+C5g , g+C5a )
    integer :: debug_binarray(9,nobs)
    real(dp) :: debug_weightarray(9,nobs)
    real(dp) :: debug_kin(9)
    debug_binarray    = zero
    debug_weightarray = zero
    debug_kin         = zero
    check = .true.
#endif


    tag_ns_wga = one 

    callcount = callcount + 1

    xsect_nnlo_5163_rr_z_ga= 0
    ff(1) = zero
    xx = zero


    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)

    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif



    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buff_rr) return



    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call open_histo()

    call kinematics_nnlo_5163(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S4Lim,S5Lim,S4S5Lim,C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,&
         C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim)  

    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!

    !1-- regular piece ; gluon-photon

    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)

    pdec(:,1) = HardProc%AmpMom(:,3)
    pdec(:,2) = HardProc%AmpMom(:,4)
    pdec(:,3) = HardProc%AmpMom(:,6)

    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then
       kin(1) = zero
       FintNNLO_dy_ns(1) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)


       call res_tree_ga_qqb(HardProc%AmpMom,resnnlo)

       !call get_respdf(ns_lumi,1,1,HardProc,resnnlo,respdf)

       kin(1) = respdf(1)*HardProc%wgt*jac

       FintNNLO_dy_ns(1) = kin(1)

       call fill_histo([kin(1)],vegasweight)

    endif


    !! -----------------------------!!
    !!     S5 + S5/C5          !!
    !! -----------------------------!!

    !-- 2 -> S5   ; photon

    xi1 = S5Lim%PartFrac(1) 
    xi2 = S5Lim%PartFrac(2)


    !print*, 'S5/S5C5'
    call cut_histo(S5Lim)

    if (S5Lim%makecut.or.S5Lim%flag) then

       kin(2) = zero
       FintNNLO_dy_ns(2) = zero
       FintNNLO_dy_ns(3) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_g_qqb(S5Lim%AmpMom,res_S5_nlo)
       call get_Eik_ph_z_ga(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3), &
            S5Lim%LimMom(:,4), S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5_nlo,eik,res2)
       !call get_respdf(ns_lumi,1,1,S5Lim,res2,respdf_ga)

       FintNNLO_dy_ns(2) = (-one)*respdf_ga(1)*S5Lim%wgt*jac


       !-- 3 -> S5C5 ; photon

       s65 = C5S5Lim%Lim_KinInv(2)
       z5 = C5S5Lim%Lim_KinInv(5)

       !call get_respdf(ns_lumi,1,1,S5Lim,res_S5_nlo,respdf_qcd)

       FintNNLO_dy_ns(3) = four/z5/s65 * respdf_qcd(1) * C5S5Lim%wgt*jac

       kin(2) = FintNNLO_dy_ns(2) + FintNNLO_dy_ns(3)

       call fill_histo([kin(2)],vegasweight)

    endif


    !-- 4 -> C4 ; 

    xi1 = C4Lim%PartFrac(1) 
    xi2 = C4Lim%PartFrac(2)


    call cut_histo(C4Lim)

    if (C4Lim%makecut.or.C4Lim%flag) then
       kin(3) = zero
       FintNNLO_dy_ns(4) = zero  

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_a_qqb(C4Lim%AmpMom,res_a_nlo)
       !call get_respdf(ns_lumi,1,1,C4Lim,res_a_nlo,respdf_b_qed)

       s14 = C4Lim%Lim_KinInv(1) 
       z4  = C4Lim%Lim_KinInv(4)

       FintNNLO_dy_ns(4)  = (-one)*(-two/s14)* Pqq(z4,Cf) * respdf_b_qed(1) *C4Lim%wgt*jac


       kin(3) = FintNNLO_dy_ns(4) 
       call fill_histo([kin(3)],vegasweight)

    endif


    !-- 5 -> C5 ; photon

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)

    !print*, 'C5'
    call cut_histo(C5Lim)

    if (C5Lim%makecut.or.C5Lim%flag) then
       kin(4) = zero
       FintNNLO_dy_ns(5)  = zero
    else


       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_g_qqb(C5Lim%AmpMom,res_C5)

       !call get_respdf(ns_lumi,1,1,C5Lim,res_C5,respdf_ch1)

       s65 = C5Lim%Lim_KinInv(2) 
       z5  = C5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(5)  = (-two/s65)* Pqq(z5,one)*( respdf_ch1(1) ) * C5Lim%wgt*jac

       kin(4) = FintNNLO_dy_ns(5) 
       call fill_histo([kin(5)],vegasweight)

    endif

    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!

    !-- 6 -> S4S5
    !-- 7 -> C4DS
    !-- 8 -> C5S4S5
    !-- 9 -> C4C5S4S5

    xi1 = S4S5Lim%PartFrac(1) 
    xi2 = S4S5Lim%PartFrac(2)

    !print*, 'S4S5'
    call cut_histo(S4S5Lim)
    if (S4S5Lim%makecut.or.S4S5Lim%flag) then
       kin(5) = zero
       FintNNLO_dy_ns(6:9)  = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(S4S5Lim%AmpMom,res_S5S4) 

       call get_Eik_ds_z_ga(S4S5Lim%LimMom(:,1),S4S5Lim%LimMom(:,2),S4S5Lim%LimMom(:,3),&
            S4S5Lim%LimMom(:,4),S4S5Lim%LimMom(:,5),S4S5Lim%LimMom(:,6),eik)
       call product_for_pdf(res_S5S4,eik,res2)
       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ag)

       FintNNLO_dy_ns(6) = (-one)*(respdf_ag(1))*S4S5Lim%wgt*jac

       s14 = C4S4S5Lim%Lim_KinInv(4) 
       z4  = C4S4S5Lim%Lim_KinInv(5) 

       call get_Eik_ph_z_ga(C4S4S5Lim%LimMom(:,1),C4S4S5Lim%LimMom(:,2),C4S4S5Lim%LimMom(:,3), &
            C4S4S5Lim%LimMom(:,4), C4S4S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5S4,eik,res2)  

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ga)


       FintNNLO_dy_ns(7) =(-one)*(respdf_ga(1))*four*Cf/s14/z4*C4S4S5Lim%wgt*jac


       !C5S4S5

       s12 = C5S4S5Lim%Lim_KinInv(1) 
       s14 = C5S4S5Lim%Lim_KinInv(2) 
       s24 = C5S4S5Lim%Lim_KinInv(3) 
       s65 = C5S4S5Lim%Lim_KinInv(4) 
       z5  = C5S4S5Lim%Lim_KinInv(5)

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res_S5S4,respdf)

       FintNNLO_dy_ns(8) = (-one)*respdf(1) * four/s65/z5 * Cf*four*s12/s14/s24*C5S4S5Lim%wgt*jac


       !C4C5S4S5

       s14 = C4C5S4S5Lim%Lim_KinInv(1) 
       s65 = C4C5S4S5Lim%Lim_KinInv(2) 
       z5 =  C4C5S4S5Lim%Lim_KinInv(5) 
       z4 =  C4C5S4S5Lim%Lim_KinInv(4)

       FintNNLO_dy_ns(9) =  16.0_dp*Cf/s14/z4/s65/z5 &
            * respdf(1)* C4C5S4S5Lim%wgt*jac


       kin(5) = sum(FintNNLO_dy_ns(6:9))

       call fill_histo([kin(5)],vegasweight)

    endif


    !-- 10 -> C4S5
    !-- 11 -> C4C5S5

    xi1 = C4S5Lim%PartFrac(1) 
    xi2 = C4S5Lim%PartFrac(2)

    !print*, 'C4S5'
    call cut_histo(C4S5Lim)

    if (C4S5Lim%makecut.or.C4S5Lim%flag) then

       kin(6) = zero
       FintNNLO_dy_ns(10:11)  = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C4S5Lim%AmpMom,res_C4S5)         
       call get_Eik_ph_z_ga(C4S5Lim%LimMom(:,1),C4S5Lim%LimMom(:,2),C4S5Lim%LimMom(:,3), &
            C4S5Lim%LimMom(:,4),C4S5Lim%LimMom(:,5),eik)

       call product_for_pdf(res_C4S5,eik,res_C4S52)


       !call get_respdf(ns_lumi,1,1,C4S5Lim,res_C4S52,respdf)

       s14 = C4S5Lim%Lim_KinInv(4) 
       z4 =  C4S5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(10) = (- respdf(1))*two/s14*Pqq(z4,Cf)* C4S5Lim%wgt*jac

       !call get_respdf(ns_lumi,1,1,C4S5Lim,res_C4S5,respdf_ch1)

       s14 = C4C5S5Lim%Lim_KinInv(1) 
       s65 = C4C5S5Lim%Lim_KinInv(2) 
       z4 =  C4C5S5Lim%Lim_KinInv(4) 
       z5 =  C4C5S5Lim%Lim_KinInv(5) 
       !

       FintNNLO_dy_ns(11) = 8.0_dp/s14/s65/z5 * Pqq(z4,Cf )  *respdf_ch1(1) * C4C5S5Lim%wgt*jac

       kin(6) = sum(FintNNLO_dy_ns(10:11))

       call fill_histo([kin(6)],vegasweight)

    endif


    !-- 12 -> C4C5

    xi1 = C4C5Lim%PartFrac(1) 
    xi2 = C4C5Lim%PartFrac(2)

    !print*, 'C4C5' 
    call cut_histo(C4C5Lim)
    if (C4C5Lim%makecut.or.C4C5Lim%flag) then

       FintNNLO_dy_ns(12) = zero
       kin(7) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C4C5Lim%AmpMom,res_C4C5)

       s14 = C4C5Lim%Lim_KinInv(1) 
       s65 = C4C5Lim%Lim_KinInv(2) 
       z4 =  C4C5Lim%Lim_KinInv(4) 
       z5 =  C4C5Lim%Lim_KinInv(5)

       !call get_respdf(ns_lumi,1,1,C4C5Lim,res_C4C5,respdf)


       FintNNLO_dy_ns(12) = (-one)*four/s65/s14 * Pqq(z5,esq) * Pqq(z4,Cf )  *respdf(1) * C4C5Lim%wgt*jac


       kin(7) = FintNNLO_dy_ns(12)

       call fill_histo([kin(7)],vegasweight)

    endif



    !! -----------------------------!!
    !!     S4 + S4/C4         !!
    !! -----------------------------!!

    !--  -> S4   ; gluon
    !--  -> S4C4 ; gluon

    xi1 = S4Lim%PartFrac(1) 
    xi2 = S4Lim%PartFrac(2)

    !print*, 'S4' 
    call cut_histo(S4Lim)

    if (S4Lim%makecut.or.S4Lim%flag) then
       kin(8) = zero
       FintNNLO_dy_ns(13) = zero
       FintNNLO_dy_ns(14) = zero
    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_a_qqb(S4Lim%AmpMom,res_a_nlo)

       !call get_respdf(ns_lumi,1,1,S4Lim,res_a_nlo,respdf)

       reseik_qcd = myEik_g(S4Lim%LimMom(:,1),S4Lim%LimMom(:,2),S4Lim%LimMom(:,3),Cf)


       FintNNLO_dy_ns(13) = (-one)*respdf(1) *reseik_qcd*S4Lim%wgt*jac


       s14 = C4S4Lim%Lim_KinInv(2)
       z4 = C4S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(14) = four*Cf/z4/s14 * respdf(1) * C4S4Lim%wgt*jac


       kin(8) = FintNNLO_dy_ns(13) + FintNNLO_dy_ns(14)

       call fill_histo([kin(8)],vegasweight)


    endif


    ! - - - - - - C5S4 and  C4C5S4 ; photon

    xi1 = C5S4Lim%PartFrac(1) 
    xi2 = C5S4Lim%PartFrac(2) 

    !print*, 'C5S4'
    call cut_histo(C5S4Lim)

    if (C5S4Lim%makecut.or.C5S4Lim%flag) then



       kin(9) = zero
       FintNNLO_dy_ns(15) = zero
       FintNNLO_dy_ns(16) = zero 

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C5S4Lim%AmpMom,res_C5)      
       !call get_respdf(ns_lumi,1,1,C5S4Lim,res_C5,respdf_ch1)


       s12 = C5S4Lim%Lim_KinInv(1)
       s14 = C5S4Lim%Lim_KinInv(2)
       s24 = C5S4Lim%Lim_KinInv(3) 
       s25 = C5S4Lim%Lim_KinInv(4) 
       z5  = C5S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(15)  = respdf_ch1(1)*two*Cf*s12/s14/s24*four/s25*Pqq(z5,one) * C5S4Lim%wgt*jac



       s14 = C4C5S4Lim%Lim_KinInv(1)
       s65 = C4C5S4Lim%Lim_KinInv(2)


       z4 = C4C5S4Lim%Lim_KinInv(4)
       z5 = C4C5S4Lim%Lim_KinInv(5)


       FintNNLO_dy_ns(16)  = (-one)*respdf_ch1(1)*four*Cf/s65/z4*two/s14*Pqq(z5,one) * C4C5S4Lim%wgt*jac


       kin(9) = FintNNLO_dy_ns(15)+FintNNLO_dy_ns(16)


    endif

    call fill_histo([kin(9)],vegasweight)

    ff(1) = sum(kin)
    call close_histo()

    ! print*, 'xx=', xx
    ! print  *, 1, FintNNLO_dy_ns(1), 'hard'
    ! print  *, '------------------------------------------------------'
    ! print  *, 2, FintNNLO_dy_ns(13)/FintNNLO_dy_ns(1), 'x1=0'
    ! print  *, 3, FintNNLO_dy_ns(2)/FintNNLO_dy_ns(1), 'x2=0'
    ! print  *, 4, FintNNLO_dy_ns(4)/FintNNLO_dy_ns(1), 'x3=0'
    ! print  *, 5, FintNNLO_dy_ns(5)/FintNNLO_dy_ns(1), 'x4=0'
    ! print  *, '------------------------------------------------------'
    ! print  *, 6, FintNNLO_dy_ns(6)/FintNNLO_dy_ns(1),    'x1=0,x2=0'
    ! print  *, 7, FintNNLO_dy_ns(14)/FintNNLO_dy_ns(1), 'x1=0,x3=0'     
    ! print  *, 8, FintNNLO_dy_ns(15)/FintNNLO_dy_ns(1), 'x1=0,x4=0'
    ! print  *, 9, FintNNLO_dy_ns(10)/FintNNLO_dy_ns(1), 'x2=0, x3=0'
    ! print  *, 10, FintNNLO_dy_ns(3)/FintNNLO_dy_ns(1), 'x2=0,x4=0'
    ! print  *, 11, FintNNLO_dy_ns(12)/FintNNLO_dy_ns(1), 'x3=0,x4=0'
    ! print  *, '------------------------------------------------------'     
    ! print  *, 12, FintNNLO_dy_ns(7)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0'
    ! print  *, 13, FintNNLO_dy_ns(8)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x4=0'
    ! print  *, 14, FintNNLO_dy_ns(16)/FintNNLO_dy_ns(1), 'x1=0,x3=0,x4=0'  
    ! print  *, 15, FintNNLO_dy_ns(11)/FintNNLO_dy_ns(1), 'x2=0,x3=0,x4=0'
    ! print  *, '------------------------------------------------------'     
    ! print  *, 16, FintNNLO_dy_ns(9)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0,x4=0'


    ! print*, 'ff(1)', ff(1)

    ! pause

    !print*, '5163 ff(1)', ff(1)
    !pause

!!$
!!$     if (ff(1) .ne. ff(1)) then
!!$       print  *, 'nan'
!!$       print  *, 'xx',xx
!!$       print  *, ''
!!$       do i = 1,16
!!$          print  *, i, FintNNLO_dy_ns(i)
!!$       enddo
!!$
!!$       pause
!!$       
!!$    endif
!!$
    ! pause

!!$    print *, ''
!!$    print*, 'kin(1)', kin(1)
!!$    print*, 'kin(2)', kin(2)
!!$    print*, 'kin(3)', kin(3)
!!$    print*, 'kin(4)', kin(4)
!!$    print*, 'kin(5)', kin(5)
!!$    print*, 'kin(6)', kin(6)
!!$    print*, 'kin(7)', kin(7)
!!$    print*, 'kin(8)', kin(8)
!!$    print*, 'kin(9)', kin(9)
!!$    print*

    !   print*, 'ff(1)', ff(1)
    !   pause
!!$
!!$
!!$
!!$
!!$
!!$    
!!$    if (ff(1) .ne. ff(1)) then
!!$       !print  *, 'nan'
!!$       !print  *, 'xx',xx
!!$       !print  *, ''
!!$       do i = 1,16
!!$          !print  *, i, FintNNLO_dy_ns(i)
!!$       enddo
!!$       ff(1) = zero
!!$    endif
!!$#if(_identification_check == 1)  
!!$    if( .not. check ) stop
!!$#endif
!!$
!!$#if (_DebugFeatures == 1)
!!$    FintNNLO_CheckLim_dy_ns(1) = FintNNLO_dy_ns(1)
!!$#endif
!!$#if (_DebugFeatures == 1)
!!$    FintNNLO_CheckLim_dy_ns(4:16) = FintNNLO_dy_ns(4:16)
!!$#endif
!!$
!!$    ! do i = 1,16
!!$    !    !print  *, i,FintNNLO_dy_ns(i)
!!$    ! enddo
!!$    ! !print  *, 'ff(1)', ff(1)
!!$    ! pause


    !do i=1,9

    !print*, i, kin(i)

    !enddo

    !print*, 'xx', xx

    !print *, ''

    !if (ff(1) .gt. ff1max) then
    !        print*, ff(1)
    !        print*, HardProc%wgt
    !        print*, HardProc%AmpMom
    !endif 

    !if ((HardProc%wgt .gt. wmax) .and. (ff(1) .gt. zero)) then
    !        print*, 'ff1', ff(1)
    !        print*, 'wgt', HardProc%wgt
    !        print*, HardProc%AmpMom
    !endif 

    !pause

    if (ff(1) .ne. ff(1)) then
       print  *, 'nan'
       print  *, 'xx',xx
       print  *, ''
       do i = 1,16
          print  *, i, FintNNLO_dy_ns(i)
       enddo
       ff(1)=zero
    endif

#if(_withchecks == 1)
    FintNNLO = FintNNLO_dy_ns
#endif

    return

  end function xsect_nnlo_5163_rr_z_ga





  function xsect_nnlo_5164_rr_z_ga(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_5164_rr_z_ga,ndim,ncomp,userdata,iternumber,i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,S4Lim,S4S5Lim,S5Lim,C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,&
         C5Lim,C5S4Lim,C5S4S5Lim,C5S5Lim,C4C5Lim,C4C5S4Lim,C4C5S4S5Lim,C4C5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR
    real(dp)    :: xi1,xi2,f1(-6:7),f2(-6:7),my_lumi,respdf(1),respdf_b(1)
    complex(dp) :: respdf_sping(-1:1,-1:1),respdf_sping2(-1:1,-1:1,-1:1,-1:1)
    !--


    complex(dp) :: hell_mat42(-1:1,-1:1),hell_mat51(-1:1,-1:1)
    !--



    real(dp)    :: respdf_qed(1),respdf_qcd(1)
    complex(dp) :: respdf_sping_qed(-1:1,-1:1),respdf_sping_qcd(-1:1,-1:1)
    real(dp)    :: respdf_ec,respdf_ce,respdf_b_qed(1),respdf_b_qcd(1)
    real(dp)    :: reseik_qed,reseik_qcd
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s24,s25,reseik,reseik_ds,reseik_dss5, s75
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_dy_ns(16),kin(1:9),skin(1:7)
    logical  :: makecut
    !--
    real(dp) :: mysymm = half
    real(dp) :: qcdmom(4,3),qedmom(4,3)
    !
    real(dp) :: wmom(4)
    real(dp) :: respdf_ch1(1),respdf_ch2
    !
    real(dp) :: resnnlo(2,2), res_S5S4(2,2), res2(2,2), eik(4), ker(2,2), res_S5_nlo(2,2)
    real(dp) :: res_a_nlo(2,2), res_C4S5(2,2), res_C4S52(2,2), res_C5(2,2), res_C52(2,2)
    real(dp) :: res_C4C5(2,2)
    real(dp) :: resnnlo_4q_tc(4),resnnlo_4q_fin(8)
    real(dp) :: res_g_nlo(12)
    real(dp) :: vec_a(4),vec_g(4),vec_ag(4),vec_ga(4)
    real(dp) :: vec_ag_a(4),vec_ag_g(4),vec_ga_a(4),vec_ga_g(4)
    real(dp) :: resvec(4)
    real(dp) :: res_ag(4),res_ga(4)
    real(dp) :: respdf_a(1),respdf_g(1)
    real(dp) :: respdf_ag(1),respdf_ga(1), tag_ns_wga

    real(dp) :: pdec(4,3)

#if(_identification_check == 1)
    logical    :: check, check_temp
    !the debug arrays contains hard and single unresolved contirbutions
    ! ( qqb -> ga  , qqb -> ag  , qqb->qqb , a+S5g , g+S5a a+C5g , g+C5a )
    integer :: debug_binarray(9,nobs)
    real(dp) :: debug_weightarray(9,nobs)
    real(dp) :: debug_kin(9)
    debug_binarray    = zero
    debug_weightarray = zero
    debug_kin         = zero
    check = .true.
#endif

    tag_ns_wga = one 

    callcount = callcount + 1

    xsect_nnlo_5164_rr_z_ga = 0
    ff(1) = zero
    xx = zero


    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))
    
#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif


    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buff_rr) return

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call open_histo()    

    call kinematics_nnlo_5164(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S4Lim,S5Lim,S4S5Lim,C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,&
         C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim)  

    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!

    !1-- regular piece ; gluon-photon

    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)

    pdec(:,1) = HardProc%AmpMom(:,3)
    pdec(:,2) = HardProc%AmpMom(:,4)
    pdec(:,3) = HardProc%AmpMom(:,6)

    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then
       kin(1) = zero

       FintNNLO_dy_ns(1) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)


       call res_tree_ga_qqb(HardProc%AmpMom,resnnlo)

       !call get_respdf(ns_lumi,1,1,HardProc,resnnlo,respdf)

       FintNNLO_dy_ns(1) = respdf(1)*HardProc%wgt*jac


       kin(1) = FintNNLO_dy_ns(1)


       call fill_histo([kin(1)],vegasweight)


    endif



    !! -----------------------------!!
    !!     S5 + S5/C5          !!
    !! -----------------------------!!

    !-- 2 -> S5   ; photon

    xi1 = S5Lim%PartFrac(1) 
    xi2 = S5Lim%PartFrac(2)

    call cut_histo(S5Lim)

    if (S5Lim%makecut.or.S5Lim%flag) then

       kin(2) = zero
       FintNNLO_dy_ns(2) = zero
       FintNNLO_dy_ns(3) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_g_qqb(S5Lim%AmpMom,res_S5_nlo)
       call get_Eik_ph_z_ga(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3), &
            S5Lim%LimMom(:,4), S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5_nlo,eik,res2)
       !call get_respdf(ns_lumi,1,1,S5Lim,res2,respdf_ga)

       FintNNLO_dy_ns(2) = (-one)*respdf_ga(1)*S5Lim%wgt*jac


       !-- 3 -> S5C5 ; photon

       s75 = C5S5Lim%Lim_KinInv(2)
       z5 = C5S5Lim%Lim_KinInv(5)

       !call get_respdf(ns_lumi,1,1,S5Lim,res_S5_nlo,respdf_qcd)

       FintNNLO_dy_ns(3) = four/z5/s75 * respdf_qcd(1) * C5S5Lim%wgt*jac

       kin(2) = FintNNLO_dy_ns(2) + FintNNLO_dy_ns(3)

       call fill_histo([kin(2)],vegasweight)

    endif



    !-- 4 -> C4 ; 

    xi1 = C4Lim%PartFrac(1) 
    xi2 = C4Lim%PartFrac(2)
    call cut_histo(C4Lim)

    if (C4Lim%makecut.or.C4Lim%flag) then
       kin(3) = zero
       FintNNLO_dy_ns(4) = zero  

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_a_qqb(C4Lim%AmpMom,res_a_nlo)
       !call get_respdf(ns_lumi,1,1,C4Lim,res_a_nlo,respdf_b_qed)

       s14 = C4Lim%Lim_KinInv(1) 
       z4  = C4Lim%Lim_KinInv(4)

       FintNNLO_dy_ns(4)  = (-one)*(-two/s14)* Pqq(z4,Cf) * respdf_b_qed(1) *C4Lim%wgt*jac


       kin(3) = FintNNLO_dy_ns(4) 
       call fill_histo([kin(3)],vegasweight)

    endif

    !-- 5 -> C5 ; photon

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)

    call cut_histo(C5Lim)

    if (C5Lim%makecut.or.C5Lim%flag) then
       kin(4) = zero
       FintNNLO_dy_ns(5)  = zero
    else


       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_g_qqb(C5Lim%AmpMom,res_C5)

       !call get_respdf(ns_lumi,1,1,C5Lim,res_C5,respdf_ch1)

       s75 = C5Lim%Lim_KinInv(2) 
       z5  = C5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(5)  = (-two/s75)* Pqq(z5,one)*( respdf_ch1(1) ) * C5Lim%wgt*jac

       kin(4) = FintNNLO_dy_ns(5) 
       call fill_histo([kin(5)],vegasweight)

    endif

    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!

    !-- 6 -> S4S5
    !-- 7 -> C4DS
    !-- 8 -> C5S4S5
    !-- 9 -> C4C5S4S5

    xi1 = S4S5Lim%PartFrac(1) 
    xi2 = S4S5Lim%PartFrac(2)

    call cut_histo(S4S5Lim)

    if (S4S5Lim%makecut.or.S4S5Lim%flag) then
       kin(5) = zero
       FintNNLO_dy_ns(6:9)  = zero


    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(S4S5Lim%AmpMom,res_S5S4) 

       call get_Eik_ds_z_ga(S4S5Lim%LimMom(:,1),S4S5Lim%LimMom(:,2),S4S5Lim%LimMom(:,3),&
            S4S5Lim%LimMom(:,4),S4S5Lim%LimMom(:,5),S4S5Lim%LimMom(:,6),eik)
       call product_for_pdf(res_S5S4,eik,res2)
       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ag)

       FintNNLO_dy_ns(6) = (-one)*(respdf_ag(1))*S4S5Lim%wgt*jac

       s14 = C4S4S5Lim%Lim_KinInv(4) 
       z4  = C4S4S5Lim%Lim_KinInv(5) 

       call get_Eik_ph_z_ga(C4S4S5Lim%LimMom(:,1),C4S4S5Lim%LimMom(:,2),C4S4S5Lim%LimMom(:,3), &
            C4S4S5Lim%LimMom(:,4), C4S4S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5S4,eik,res2)  

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ga)


       FintNNLO_dy_ns(7) =(-one)*(respdf_ga(1))*four*Cf/s14/z4*C4S4S5Lim%wgt*jac


       !C5S4S5

       s12 = C5S4S5Lim%Lim_KinInv(1) 
       s14 = C5S4S5Lim%Lim_KinInv(2) 
       s24 = C5S4S5Lim%Lim_KinInv(3) 
       s75 = C5S4S5Lim%Lim_KinInv(4) 
       z5  = C5S4S5Lim%Lim_KinInv(5)

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res_S5S4,respdf)

       FintNNLO_dy_ns(8) = (-one)*respdf(1) * four/s75/z5 * Cf*four*s12/s14/s24*C5S4S5Lim%wgt*jac


       !C4C5S4S5

       s14 = C4C5S4S5Lim%Lim_KinInv(1) 
       s75 = C4C5S4S5Lim%Lim_KinInv(2) 
       z5 =  C4C5S4S5Lim%Lim_KinInv(5) 
       z4 =  C4C5S4S5Lim%Lim_KinInv(4)

       FintNNLO_dy_ns(9) =  16.0_dp*Cf/s14/z4/s75/z5 &
            * respdf(1)* C4C5S4S5Lim%wgt*jac


       kin(5) = sum(FintNNLO_dy_ns(6:9))


       call fill_histo([kin(5)],vegasweight)

    endif



    !-- 10 -> C4S5
    !-- 11 -> C4C5S5

    xi1 = C4S5Lim%PartFrac(1) 
    xi2 = C4S5Lim%PartFrac(2)

    call cut_histo(C4S5Lim)

    if (C4S5Lim%makecut.or.C4S5Lim%flag) then

       kin(6) = zero
       FintNNLO_dy_ns(10:11)  = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C4S5Lim%AmpMom,res_C4S5)         
       call get_Eik_ph_z_ga(C4S5Lim%LimMom(:,1),C4S5Lim%LimMom(:,2),C4S5Lim%LimMom(:,3), &
            C4S5Lim%LimMom(:,4),C4S5Lim%LimMom(:,5),eik)

       call product_for_pdf(res_C4S5,eik,res_C4S52)


       !call get_respdf(ns_lumi,1,1,C4S5Lim,res_C4S52,respdf)

       s14 = C4S5Lim%Lim_KinInv(4) 
       z4 =  C4S5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(10) = (- respdf(1))*two/s14*Pqq(z4,Cf)* C4S5Lim%wgt*jac


       !call get_respdf(ns_lumi,1,1,C4S5Lim,res_C4S5,respdf_ch1)

       s14 = C4C5S5Lim%Lim_KinInv(1) 
       s75 = C4C5S5Lim%Lim_KinInv(2) 
       z4 =  C4C5S5Lim%Lim_KinInv(4) 
       z5 =  C4C5S5Lim%Lim_KinInv(5) 
       !

       FintNNLO_dy_ns(11) = 8.0_dp/s14/s75/z5 * Pqq(z4,Cf )  *respdf_ch1(1) * C4C5S5Lim%wgt*jac

       kin(6) = sum(FintNNLO_dy_ns(10:11))

       call fill_histo([kin(6)],vegasweight)

    endif


    !-- 12 -> C4C5

    xi1 = C4C5Lim%PartFrac(1) 
    xi2 = C4C5Lim%PartFrac(2)
    call cut_histo(C4C5Lim)
    if (C4C5Lim%makecut.or.C4C5Lim%flag) then

       FintNNLO_dy_ns(12) = zero
       kin(7) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C4C5Lim%AmpMom,res_C4C5)

       s14 = C4C5Lim%Lim_KinInv(1) 
       s75 = C4C5Lim%Lim_KinInv(2) 
       z4 =  C4C5Lim%Lim_KinInv(4) 
       z5 =  C4C5Lim%Lim_KinInv(5)

       !call get_respdf(ns_lumi,1,1,C4C5Lim,res_C4C5,respdf)


       FintNNLO_dy_ns(12) = (-one)*four/s75/s14 * Pqq(z5,esq) * Pqq(z4,Cf )  *respdf(1) * C4C5Lim%wgt*jac



       kin(7) = FintNNLO_dy_ns(12)

       call fill_histo([kin(7)],vegasweight)

    endif



    !! -----------------------------!!
    !!     S4 + S4/C4         !!
    !! -----------------------------!!

    !--  -> S4   ; gluon
    !--  -> S4C4 ; gluon

    xi1 = S4Lim%PartFrac(1) 
    xi2 = S4Lim%PartFrac(2) 
    call cut_histo(S4Lim)


    if (S4Lim%makecut.or.S4Lim%flag) then
       kin(8) = zero
       FintNNLO_dy_ns(13) = zero
       FintNNLO_dy_ns(14) = zero
    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_a_qqb(S4Lim%AmpMom,res_a_nlo)

       !call get_respdf(ns_lumi,1,1,S4Lim,res_a_nlo,respdf)

       reseik_qcd = myEik_g(S4Lim%LimMom(:,1),S4Lim%LimMom(:,2),S4Lim%LimMom(:,3),Cf)


       FintNNLO_dy_ns(13) = (-one)*respdf(1) *reseik_qcd*S4Lim%wgt*jac


       s14 = C4S4Lim%Lim_KinInv(2)
       z4 = C4S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(14) = four*Cf/z4/s14 * respdf(1) * C4S4Lim%wgt*jac


       kin(8) = FintNNLO_dy_ns(13) + FintNNLO_dy_ns(14)

       call fill_histo([kin(8)],vegasweight)


    endif


    ! - - - - - - C5S4 and  C4C5S4 ; photon

    xi1 = C5S4Lim%PartFrac(1) 
    xi2 = C5S4Lim%PartFrac(2) 
    call cut_histo(C5S4Lim)

    if (C5S4Lim%makecut.or.C5S4Lim%flag) then



       kin(9) = zero
       FintNNLO_dy_ns(15) = zero
       FintNNLO_dy_ns(16) = zero 

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C5S4Lim%AmpMom,res_C5)      
       !call get_respdf(ns_lumi,1,1,C5S4Lim,res_C5,respdf_ch1)


       s12 = C5S4Lim%Lim_KinInv(1)
       s14 = C5S4Lim%Lim_KinInv(2)
       s24 = C5S4Lim%Lim_KinInv(3) 
       s75 = C5S4Lim%Lim_KinInv(4) 
       z5  = C5S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(15)  = respdf_ch1(1)*two*Cf*s12/s14/s24*four/s75*Pqq(z5,one) * C5S4Lim%wgt*jac

       s14 = C4C5S4Lim%Lim_KinInv(1)
       s75 = C4C5S4Lim%Lim_KinInv(2)


       z4 = C4C5S4Lim%Lim_KinInv(4)
       z5 = C4C5S4Lim%Lim_KinInv(5)


       FintNNLO_dy_ns(16)  = (-one)*respdf_ch1(1)*four*Cf/s75/z4*two/s14*Pqq(z5,one) * C4C5S4Lim%wgt*jac


       kin(9) = FintNNLO_dy_ns(15)+FintNNLO_dy_ns(16)


    endif

    call fill_histo([kin(9)],vegasweight)

    ff(1) = sum(kin)
    call close_histo()
    
    !  print*, 'xx=', xx
    !     print  *, 1, FintNNLO_dy_ns(1), 'hard'
    !     print  *, '------------------------------------------------------'
    !     print  *, 2, FintNNLO_dy_ns(13)/FintNNLO_dy_ns(1), 'x1=0'
    !     print  *, 3, FintNNLO_dy_ns(2)/FintNNLO_dy_ns(1), 'x2=0'
    !     print  *, 4, FintNNLO_dy_ns(4)/FintNNLO_dy_ns(1), 'x3=0'
    !     print  *, 5, FintNNLO_dy_ns(5)/FintNNLO_dy_ns(1), 'x4=0'
    !     print  *, '------------------------------------------------------'
    !     print  *, 6, FintNNLO_dy_ns(6)/FintNNLO_dy_ns(1),    'x1=0,x2=0'
    !     print  *, 7, FintNNLO_dy_ns(14)/FintNNLO_dy_ns(1), 'x1=0,x3=0'     
    !     print  *, 8, FintNNLO_dy_ns(15)/FintNNLO_dy_ns(1), 'x1=0,x4=0'
    !     print  *, 9, FintNNLO_dy_ns(10)/FintNNLO_dy_ns(1), 'x2=0, x3=0'
    !     print  *, 10, FintNNLO_dy_ns(3)/FintNNLO_dy_ns(1), 'x2=0,x4=0'
    !     print  *, 11, FintNNLO_dy_ns(12)/FintNNLO_dy_ns(1), 'x3=0,x4=0'
    !     print  *, '------------------------------------------------------'     
    !     print  *, 12, FintNNLO_dy_ns(7)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0'
    !     print  *, 13, FintNNLO_dy_ns(8)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x4=0'
    !     print  *, 14, FintNNLO_dy_ns(16)/FintNNLO_dy_ns(1), 'x1=0,x3=0,x4=0'  
    !     print  *, 15, FintNNLO_dy_ns(11)/FintNNLO_dy_ns(1), 'x2=0,x3=0,x4=0'
    !     print  *, '------------------------------------------------------'     
    !     print  *, 16, FintNNLO_dy_ns(9)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0,x4=0'


    !     print*, 'ff(1)', ff(1)
    !     pause
    ! !!$
!!$     if (ff(1) .ne. ff(1)) then
!!$       print  *, 'nan'
!!$       print  *, 'xx',xx
!!$       print  *, ''
!!$       do i = 1,16
!!$          print  *, i, FintNNLO_dy_ns(i)
!!$       enddo
!!$
!!$       pause
!!$       
!!$    endif
!!$
!!$    pause
!!$    

!!$    pause

!!$    print *, ''
!!$    print*, 'kin(1)', kin(1)
!!$    print*, 'kin(2)', kin(2)
!!$    print*, 'kin(3)', kin(3)
!!$    print*, 'kin(4)', kin(4)
!!$    print*, 'kin(5)', kin(5)
!!$    print*, 'kin(6)', kin(6)
!!$    print*, 'kin(7)', kin(7)
!!$    print*, 'kin(8)', kin(8)
!!$    print*, 'kin(9)', kin(9)
!!$    print*

    !print*, 'pref', pref

    !   print*, 'ff(1)', ff(1)
    !   pause
!!$
!!$
!!$
!!$
!!$
!!$    
!!$    if (ff(1) .ne. ff(1)) then
!!$       !print  *, 'nan'
!!$       !print  *, 'xx',xx
!!$       !print  *, ''
!!$       do i = 1,16
!!$          !print  *, i, FintNNLO_dy_ns(i)
!!$       enddo
!!$       ff(1) = zero
!!$    endif
!!$#if(_identification_check == 1)  
!!$    if( .not. check ) stop
!!$#endif
!!$
!!$#if (_DebugFeatures == 1)
!!$    FintNNLO_CheckLim_dy_ns(1) = FintNNLO_dy_ns(1)
!!$#endif
!!$#if (_DebugFeatures == 1)
!!$    FintNNLO_CheckLim_dy_ns(4:16) = FintNNLO_dy_ns(4:16)
!!$#endif
!!$
!!$    ! do i = 1,16
!!$    !    !print  *, i,FintNNLO_dy_ns(i)
!!$    ! enddo
!!$    ! !print  *, 'ff(1)', ff(1)
!!$    ! pause


    if (ff(1) .ne. ff(1)) then
       print  *, 'nan'
       print  *, 'xx',xx
       print  *, ''
       do i = 1,16
          print  *, i, FintNNLO_dy_ns(i)
       enddo
       ff(1)=zero
    endif

#if(_withchecks == 1)
    FintNNLO = FintNNLO_dy_ns
#endif


    return

  end function xsect_nnlo_5164_rr_z_ga







  function xsect_nnlo_5263_rr_z_ga(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_5263_rr_z_ga,ndim,ncomp,userdata,iternumber,i, j
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,S4Lim,S4S5Lim,S5Lim,C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,&
         C5Lim,C5S4Lim,C5S4S5Lim,C5S5Lim,C4C5Lim,C4C5S4Lim,C4C5S4S5Lim,C4C5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR
    real(dp)    :: xi1,xi2,f1(-6:7),f2(-6:7),my_lumi,respdf(1),respdf_b(1)
    complex(dp) :: respdf_sping(-1:1,-1:1),respdf_sping2(-1:1,-1:1,-1:1,-1:1)
    !--


    complex(dp) :: hell_mat42(-1:1,-1:1),hell_mat51(-1:1,-1:1)
    !--



    real(dp)    :: respdf_qed(1),respdf_qcd(1)
    complex(dp) :: respdf_sping_qed(-1:1,-1:1),respdf_sping_qcd(-1:1,-1:1)
    real(dp)    :: respdf_ec,respdf_ce,respdf_b_qed(1),respdf_b_qcd(1)
    real(dp)    :: reseik_qed,reseik_qcd
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s24,s25,reseik,reseik_ds,reseik_dss5, s75
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_dy_ns(16),kin(1:9),skin(1:7)
    logical  :: makecut
    !--
    real(dp) :: mysymm = half
    real(dp) :: qcdmom(4,3),qedmom(4,3)
    !
    real(dp) :: wmom(4)
    real(dp) :: respdf_ch1(1),respdf_ch2
    !
    real(dp) :: resnnlo(2,2), res_S5S4(2,2), res2(2,2), eik(4), ker(2,2), res_S5_nlo(2,2)
    real(dp) :: res_a_nlo(2,2), res_C4S5(2,2), res_C4S52(2,2), res_C5(2,2), res_C52(2,2)
    real(dp) :: res_C4C5(2,2)
    real(dp) :: resnnlo_4q_tc(4),resnnlo_4q_fin(8)
    real(dp) :: res_g_nlo(12)
    real(dp) :: vec_a(4),vec_g(4),vec_ag(4),vec_ga(4)
    real(dp) :: vec_ag_a(4),vec_ag_g(4),vec_ga_a(4),vec_ga_g(4)
    real(dp) :: resvec(4)
    real(dp) :: res_ag(4),res_ga(4)
    real(dp) :: respdf_a(1),respdf_g(1)
    real(dp) :: respdf_ag(1),respdf_ga(1), tag_ns_wga

    real(dp) :: pdec(4,3)

#if(_identification_check == 1)
    logical    :: check, check_temp
    !the debug arrays contains hard and single unresolved contirbutions
    ! ( qqb -> ga  , qqb -> ag  , qqb->qqb , a+S5g , g+S5a a+C5g , g+C5a )
    integer :: debug_binarray(9,nobs)
    real(dp) :: debug_weightarray(9,nobs)
    real(dp) :: debug_kin(9)
    debug_binarray    = zero
    debug_weightarray = zero
    debug_kin         = zero
    check = .true.
#endif

    tag_ns_wga = one 

    callcount = callcount + 1

    xsect_nnlo_5263_rr_z_ga = 0
    ff(1) = zero
    xx = zero


    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buff_rr) return

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call open_histo()    
    
    call kinematics_nnlo_5263(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S4Lim,S5Lim,S4S5Lim,C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,&
         C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim)  

    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!

    !1-- regular piece ; gluon-photon

    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)

    pdec(:,1) = HardProc%AmpMom(:,3)
    pdec(:,2) = HardProc%AmpMom(:,4)
    pdec(:,3) = HardProc%AmpMom(:,6)

    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then
       kin(1) = zero

       FintNNLO_dy_ns(1) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)


       call res_tree_ga_qqb(HardProc%AmpMom,resnnlo)

       !call get_respdf(ns_lumi,1,1,HardProc,resnnlo,respdf)

       FintNNLO_dy_ns(1) = respdf(1)*HardProc%wgt*jac


       kin(1) = FintNNLO_dy_ns(1)


       call fill_histo([kin(1)],vegasweight)


    endif



    !! -----------------------------!!
    !!     S5 + S5/C5          !!
    !! -----------------------------!!

    !-- 2 -> S5   ; photon

    xi1 = S5Lim%PartFrac(1) 
    xi2 = S5Lim%PartFrac(2)

    call cut_histo(S5Lim)

    if (S5Lim%makecut.or.S5Lim%flag) then

       kin(2) = zero
       FintNNLO_dy_ns(2) = zero
       FintNNLO_dy_ns(3) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_g_qqb(S5Lim%AmpMom,res_S5_nlo)
       call get_Eik_ph_z_ga(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3), &
            S5Lim%LimMom(:,4), S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5_nlo,eik,res2)
       !call get_respdf(ns_lumi,1,1,S5Lim,res2,respdf_ga)

       FintNNLO_dy_ns(2) = (-one)*respdf_ga(1)*S5Lim%wgt*jac


       !-- 3 -> S5C5 ; photon

       s75 = C5S5Lim%Lim_KinInv(2)
       z5 = C5S5Lim%Lim_KinInv(5)

       !call get_respdf(ns_lumi,1,1,S5Lim,res_S5_nlo,respdf_qcd)

       FintNNLO_dy_ns(3) = four/z5/s75 * respdf_qcd(1) * C5S5Lim%wgt*jac

       kin(2) = FintNNLO_dy_ns(2) + FintNNLO_dy_ns(3)

       call fill_histo([kin(2)],vegasweight)

    endif



    !-- 4 -> C4 ; 

    xi1 = C4Lim%PartFrac(1) 
    xi2 = C4Lim%PartFrac(2)
    call cut_histo(C4Lim)

    if (C4Lim%makecut.or.C4Lim%flag) then
       kin(3) = zero
       FintNNLO_dy_ns(4) = zero  

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_a_qqb(C4Lim%AmpMom,res_a_nlo)
       !call get_respdf(ns_lumi,1,1,C4Lim,res_a_nlo,respdf_b_qed)

       s14 = C4Lim%Lim_KinInv(1) 
       z4  = C4Lim%Lim_KinInv(4)

       FintNNLO_dy_ns(4)  = (-one)*(-two/s14)* Pqq(z4,Cf) * respdf_b_qed(1) *C4Lim%wgt*jac


       kin(3) = FintNNLO_dy_ns(4) 
       call fill_histo([kin(3)],vegasweight)

    endif

    !-- 5 -> C5 ; photon

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)

    call cut_histo(C5Lim)

    if (C5Lim%makecut.or.C5Lim%flag) then
       kin(4) = zero
       FintNNLO_dy_ns(5)  = zero
    else


       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_g_qqb(C5Lim%AmpMom,res_C5)

       !call get_respdf(ns_lumi,1,1,C5Lim,res_C5,respdf_ch1)

       s75 = C5Lim%Lim_KinInv(2) 
       z5  = C5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(5)  = (-two/s75)* Pqq(z5,one)*( respdf_ch1(1) ) * C5Lim%wgt*jac

!!$        print *, ''
!!$
!!$       print*, 'C5'
!!$       print*, '(-two/s75)* Pqq(z5,one)*( respdf_ch1(1) )',(-two/s75)* Pqq(z5,one)*( respdf_ch1(1) )
!!$       print *, ''
!!$
!!$       print*, 's75',s75
!!$       print*, ' Pqq(z5,one)',  Pqq(z5,one)
!!$       print*, 'respdf_ch1(1)', respdf_ch1(1) 
!!$       print *, ''
!!$
!!$       
!!$       print*, 'res_C5(1,1)', res_C5(1,1)
!!$       print*, 'res_C5(1,2)', res_C5(1,2)
!!$       print*, 'res_C5(2,1)', res_C5(2,1)
!!$       print*, 'res_C5(2,2)', res_C5(2,2)
!!$       print *, ''
!!$
!!$       print*, 'C5Lim%AmpMom(:,1)',C5Lim%AmpMom(:,1)
!!$       print*, 'C5Lim%AmpMom(:,2)',C5Lim%AmpMom(:,2)
!!$       print*, 'C5Lim%AmpMom(:,3)',C5Lim%AmpMom(:,3)
!!$       print*, 'C5Lim%AmpMom(:,4)',C5Lim%AmpMom(:,4)
!!$       print*, 'C5Lim%AmpMom(:,5)',C5Lim%AmpMom(:,5)
!!$       
!!$       print *, ''
!!$       print*, 'C5Lim%wgt', C5Lim%wgt
!!$       print*, 'FintNNLO_dy_ns(5)', FintNNLO_dy_ns(5)
!!$
!!$       print *, ''



       kin(4) = FintNNLO_dy_ns(5) 
       call fill_histo([kin(5)],vegasweight)

    endif

    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!

    !-- 6 -> S4S5
    !-- 7 -> C4DS
    !-- 8 -> C5S4S5
    !-- 9 -> C4C5S4S5

    xi1 = S4S5Lim%PartFrac(1) 
    xi2 = S4S5Lim%PartFrac(2)

    call cut_histo(S4S5Lim)

    if (S4S5Lim%makecut.or.S4S5Lim%flag) then
       kin(5) = zero
       FintNNLO_dy_ns(6:9)  = zero


    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(S4S5Lim%AmpMom,res_S5S4) 

       call get_Eik_ds_z_ga(S4S5Lim%LimMom(:,1),S4S5Lim%LimMom(:,2),S4S5Lim%LimMom(:,3),&
            S4S5Lim%LimMom(:,4),S4S5Lim%LimMom(:,5),S4S5Lim%LimMom(:,6),eik)
       call product_for_pdf(res_S5S4,eik,res2)
       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ag)

       FintNNLO_dy_ns(6) = (-one)*(respdf_ag(1))*S4S5Lim%wgt*jac

       s14 = C4S4S5Lim%Lim_KinInv(4) 
       z4  = C4S4S5Lim%Lim_KinInv(5) 

       call get_Eik_ph_z_ga(C4S4S5Lim%LimMom(:,1),C4S4S5Lim%LimMom(:,2),C4S4S5Lim%LimMom(:,3), &
            C4S4S5Lim%LimMom(:,4), C4S4S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5S4,eik,res2)  

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ga)


       FintNNLO_dy_ns(7) =(-one)*(respdf_ga(1))*four*Cf/s14/z4*C4S4S5Lim%wgt*jac


       !C5S4S5

       s12 = C5S4S5Lim%Lim_KinInv(1) 
       s14 = C5S4S5Lim%Lim_KinInv(2) 
       s24 = C5S4S5Lim%Lim_KinInv(3) 
       s75 = C5S4S5Lim%Lim_KinInv(4) 
       z5  = C5S4S5Lim%Lim_KinInv(5)

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res_S5S4,respdf)

       FintNNLO_dy_ns(8) = (-one)*respdf(1) * four/s75/z5 * Cf*four*s12/s14/s24*C5S4S5Lim%wgt*jac


       !C4C5S4S5

       s14 = C4C5S4S5Lim%Lim_KinInv(1) 
       s75 = C4C5S4S5Lim%Lim_KinInv(2) 
       z5 =  C4C5S4S5Lim%Lim_KinInv(5) 
       z4 =  C4C5S4S5Lim%Lim_KinInv(4)

       FintNNLO_dy_ns(9) =  16.0_dp*Cf/s14/z4/s75/z5 &
            * respdf(1)* C4C5S4S5Lim%wgt*jac


       kin(5) = sum(FintNNLO_dy_ns(6:9))


       call fill_histo([kin(5)],vegasweight)

    endif



    !-- 10 -> C4S5
    !-- 11 -> C4C5S5

    xi1 = C4S5Lim%PartFrac(1) 
    xi2 = C4S5Lim%PartFrac(2)

    call cut_histo(C4S5Lim)

    if (C4S5Lim%makecut.or.C4S5Lim%flag) then

       kin(6) = zero
       FintNNLO_dy_ns(10:11)  = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C4S5Lim%AmpMom,res_C4S5)         
       call get_Eik_ph_z_ga(C4S5Lim%LimMom(:,1),C4S5Lim%LimMom(:,2),C4S5Lim%LimMom(:,3), &
            C4S5Lim%LimMom(:,4),C4S5Lim%LimMom(:,5),eik)

       call product_for_pdf(res_C4S5,eik,res_C4S52)


       !call get_respdf(ns_lumi,1,1,C4S5Lim,res_C4S52,respdf)

       s14 = C4S5Lim%Lim_KinInv(4) 
       z4 =  C4S5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(10) = (- respdf(1))*two/s14*Pqq(z4,Cf)* C4S5Lim%wgt*jac


       !call get_respdf(ns_lumi,1,1,C4S5Lim,res_C4S5,respdf_ch1)

       s14 = C4C5S5Lim%Lim_KinInv(1) 
       s75 = C4C5S5Lim%Lim_KinInv(2) 
       z4 =  C4C5S5Lim%Lim_KinInv(4) 
       z5 =  C4C5S5Lim%Lim_KinInv(5) 
       !

       FintNNLO_dy_ns(11) = 8.0_dp/s14/s75/z5 * Pqq(z4,Cf )  *respdf_ch1(1) * C4C5S5Lim%wgt*jac

       kin(6) = sum(FintNNLO_dy_ns(10:11))

       call fill_histo([kin(6)],vegasweight)

    endif


    !-- 12 -> C4C5

    xi1 = C4C5Lim%PartFrac(1) 
    xi2 = C4C5Lim%PartFrac(2) 
    call cut_histo(C4C5Lim)
    if (C4C5Lim%makecut.or.C4C5Lim%flag) then

       FintNNLO_dy_ns(12) = zero
       kin(7) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C4C5Lim%AmpMom,res_C4C5)

       s14 = C4C5Lim%Lim_KinInv(1) 
       s75 = C4C5Lim%Lim_KinInv(2) 
       z4 =  C4C5Lim%Lim_KinInv(4) 
       z5 =  C4C5Lim%Lim_KinInv(5)

       !call get_respdf(ns_lumi,1,1,C4C5Lim,res_C4C5,respdf)


       FintNNLO_dy_ns(12) = (-one)*four/s75/s14 * Pqq(z5,esq) * Pqq(z4,Cf )  *respdf(1) * C4C5Lim%wgt*jac



       kin(7) = FintNNLO_dy_ns(12)

       call fill_histo([kin(7)],vegasweight)

    endif



    !! -----------------------------!!
    !!     S4 + S4/C4         !!
    !! -----------------------------!!

    !--  -> S4   ; gluon
    !--  -> S4C4 ; gluon

    xi1 = S4Lim%PartFrac(1) 
    xi2 = S4Lim%PartFrac(2) 
    call cut_histo(S4Lim)


    if (S4Lim%makecut.or.S4Lim%flag) then
       kin(8) = zero
       FintNNLO_dy_ns(13) = zero
       FintNNLO_dy_ns(14) = zero
    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_a_qqb(S4Lim%AmpMom,res_a_nlo)

       !call get_respdf(ns_lumi,1,1,S4Lim,res_a_nlo,respdf)

       reseik_qcd = myEik_g(S4Lim%LimMom(:,1),S4Lim%LimMom(:,2),S4Lim%LimMom(:,3),Cf)


       FintNNLO_dy_ns(13) = (-one)*respdf(1) *reseik_qcd*S4Lim%wgt*jac


       s14 = C4S4Lim%Lim_KinInv(2)
       z4 = C4S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(14) = four*Cf/z4/s14 * respdf(1) * C4S4Lim%wgt*jac


       kin(8) = FintNNLO_dy_ns(13) + FintNNLO_dy_ns(14)

       call fill_histo([kin(8)],vegasweight)


    endif


    ! - - - - - - C5S4 and  C4C5S4 ; photon

    xi1 = C5S4Lim%PartFrac(1) 
    xi2 = C5S4Lim%PartFrac(2) 
    call cut_histo(C5S4Lim)

    if (C5S4Lim%makecut.or.C5S4Lim%flag) then



       kin(9) = zero
       FintNNLO_dy_ns(15) = zero
       FintNNLO_dy_ns(16) = zero 

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C5S4Lim%AmpMom,res_C5)      
       !call get_respdf(ns_lumi,1,1,C5S4Lim,res_C5,respdf_ch1)


       s12 = C5S4Lim%Lim_KinInv(1)
       s14 = C5S4Lim%Lim_KinInv(2)
       s24 = C5S4Lim%Lim_KinInv(3) 
       s75 = C5S4Lim%Lim_KinInv(4) 
       z5  = C5S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(15)  = respdf_ch1(1)*two*Cf*s12/s14/s24*four/s75*Pqq(z5,one) * C5S4Lim%wgt*jac

       s14 = C4C5S4Lim%Lim_KinInv(1)
       s75 = C4C5S4Lim%Lim_KinInv(2)


       z4 = C4C5S4Lim%Lim_KinInv(4)
       z5 = C4C5S4Lim%Lim_KinInv(5)


       FintNNLO_dy_ns(16)  = (-one)*respdf_ch1(1)*four*Cf/s75/z4*two/s14*Pqq(z5,one) * C4C5S4Lim%wgt*jac


       kin(9) = FintNNLO_dy_ns(15)+FintNNLO_dy_ns(16)


    endif

    call fill_histo([kin(9)],vegasweight)

    ff(1) = sum(kin)
    call close_histo()
    
    ! print*, 'xx=', xx
    ! print  *, 1, FintNNLO_dy_ns(1), 'hard'
    ! print  *, '------------------------------------------------------'
    ! print  *, 2, FintNNLO_dy_ns(13)/FintNNLO_dy_ns(1), 'x1=0'
    ! print  *, 3, FintNNLO_dy_ns(2)/FintNNLO_dy_ns(1), 'x2=0'
    ! print  *, 4, FintNNLO_dy_ns(4)/FintNNLO_dy_ns(1), 'x3=0'
    ! print  *, 5, FintNNLO_dy_ns(5)/FintNNLO_dy_ns(1), 'x4=0'
    ! print  *, '------------------------------------------------------'
    ! print  *, 6, FintNNLO_dy_ns(6)/FintNNLO_dy_ns(1),    'x1=0,x2=0'
    ! print  *, 7, FintNNLO_dy_ns(14)/FintNNLO_dy_ns(1), 'x1=0,x3=0'     
    ! print  *, 8, FintNNLO_dy_ns(15)/FintNNLO_dy_ns(1), 'x1=0,x4=0'
    ! print  *, 9, FintNNLO_dy_ns(10)/FintNNLO_dy_ns(1), 'x2=0, x3=0'
    ! print  *, 10, FintNNLO_dy_ns(3)/FintNNLO_dy_ns(1), 'x2=0,x4=0'
    ! print  *, 11, FintNNLO_dy_ns(12)/FintNNLO_dy_ns(1), 'x3=0,x4=0'
    ! print  *, '------------------------------------------------------'     
    ! print  *, 12, FintNNLO_dy_ns(7)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0'
    ! print  *, 13, FintNNLO_dy_ns(8)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x4=0'
    ! print  *, 14, FintNNLO_dy_ns(16)/FintNNLO_dy_ns(1), 'x1=0,x3=0,x4=0'  
    ! print  *, 15, FintNNLO_dy_ns(11)/FintNNLO_dy_ns(1), 'x2=0,x3=0,x4=0'
    ! print  *, '------------------------------------------------------'     
    ! print  *, 16, FintNNLO_dy_ns(9)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0,x4=0'

    ! print*, 'ff(1)', ff(1)

    ! pause

    !print*, pref

    !    print*, '5263 ff(1)', ff(1)

    !    pause
!!$
!!$    do i = 1,16
!!$
!!$       if (FintNNLO_dy_ns(i) .ne. FintNNLO_dy_ns(i)) then
!!$       print  *, 'STOP'
!!$       print  *, 'xx',xx
!!$       do j = 1,16
!!$          print  *, j, FintNNLO_dy_ns(j)
!!$       enddo
!!$    endif
!!$    
!!$
!!$ enddo
!!$ 
!!$
!!$     if (ff(1) .ne. ff(1)) then
!!$       print  *, 'nan'
!!$       print  *, 'xx',xx
!!$       print  *, ''
!!$       do i = 1,16
!!$          print  *, i, FintNNLO_dy_ns(i)
!!$       enddo
!!$
!!$       pause
!!$       
!!$    endif
!!$
!!$    pause


!!$    pause

!!$    print *, ''
!!$    print*, 'kin(1)', kin(1)
!!$    print*, 'kin(2)', kin(2)
!!$    print*, 'kin(3)', kin(3)
!!$    print*, 'kin(4)', kin(4)
!!$    print*, 'kin(5)', kin(5)
!!$    print*, 'kin(6)', kin(6)
!!$    print*, 'kin(7)', kin(7)
!!$    print*, 'kin(8)', kin(8)
!!$    print*, 'kin(9)', kin(9)
!!$    print*

    !   print*, 'ff(1)', ff(1)
    !   pause
!!$
!!$
!!$
!!$
!!$
!!$    
!!$    if (ff(1) .ne. ff(1)) then
!!$       !print  *, 'nan'
!!$       !print  *, 'xx',xx
!!$       !print  *, ''
!!$       do i = 1,16
!!$          !print  *, i, FintNNLO_dy_ns(i)
!!$       enddo
!!$       ff(1) = zero
!!$    endif
!!$#if(_identification_check == 1)  
!!$    if( .not. check ) stop
!!$#endif
!!$
!!$#if (_DebugFeatures == 1)
!!$    FintNNLO_CheckLim_dy_ns(1) = FintNNLO_dy_ns(1)
!!$#endif
!!$#if (_DebugFeatures == 1)
!!$    FintNNLO_CheckLim_dy_ns(4:16) = FintNNLO_dy_ns(4:16)
!!$#endif
!!$
!!$    ! do i = 1,16
!!$    !    !print  *, i,FintNNLO_dy_ns(i)
!!$    ! enddo
!!$    ! !print  *, 'ff(1)', ff(1)
!!$    ! pause


    if (ff(1) .ne. ff(1)) then
       print  *, 'nan'
       print  *, 'xx',xx
       print  *, ''
       do i = 1,16
          print  *, i, FintNNLO_dy_ns(i)
       enddo
       ff(1)=zero
    endif

#if(_withchecks == 1)
    FintNNLO = FintNNLO_dy_ns
#endif

    return

  end function xsect_nnlo_5263_rr_z_ga





  function xsect_nnlo_5264_rr_z_ga(ndim,yRnd,ncomp,ff,userdata,vegasweight,iternumber)
    integer :: xsect_nnlo_5264_rr_z_ga,ndim,ncomp,userdata,iternumber,i
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,S4Lim,S4S5Lim,S5Lim,C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,&
         C5Lim,C5S4Lim,C5S4S5Lim,C5S5Lim,C4C5Lim,C4C5S4Lim,C4C5S4S5Lim,C4C5S5Lim
    !--
    real(dp)    :: xx(kNNLO_max_full)
    real(dp)    :: tau,M2,ylab,jac,pref,mu,muF,muR
    real(dp)    :: xi1,xi2,f1(-6:7),f2(-6:7),my_lumi,respdf(1),respdf_b(1)
    complex(dp) :: respdf_sping(-1:1,-1:1),respdf_sping2(-1:1,-1:1,-1:1,-1:1)
    !--


    complex(dp) :: hell_mat42(-1:1,-1:1),hell_mat51(-1:1,-1:1)
    !--



    real(dp)    :: respdf_qed(1),respdf_qcd(1)
    complex(dp) :: respdf_sping_qed(-1:1,-1:1),respdf_sping_qcd(-1:1,-1:1)
    real(dp)    :: respdf_ec,respdf_ce,respdf_b_qed(1),respdf_b_qcd(1)
    real(dp)    :: reseik_qed,reseik_qcd
    !--
    real(dp) :: z1,z2,z4,z5,s12,s14,s15,s24,s25,reseik,reseik_ds,reseik_dss5, s75
    real(dp) :: ncoll(3),nperp(3)
    !--
    real(dp) :: FintNNLO_dy_ns(16),kin(1:9),skin(1:7)
    logical  :: makecut
    !--
    real(dp) :: mysymm = half
    real(dp) :: qcdmom(4,3),qedmom(4,3)
    !
    real(dp) :: wmom(4)
    real(dp) :: respdf_ch1(1),respdf_ch2
    !
    real(dp) :: resnnlo(2,2), res_S5S4(2,2), res2(2,2), eik(4), ker(2,2), res_S5_nlo(2,2)
    real(dp) :: res_a_nlo(2,2), res_C4S5(2,2), res_C4S52(2,2), res_C5(2,2), res_C52(2,2)
    real(dp) :: res_C4C5(2,2)
    real(dp) :: resnnlo_4q_tc(4),resnnlo_4q_fin(8)
    real(dp) :: res_g_nlo(12)
    real(dp) :: vec_a(4),vec_g(4),vec_ag(4),vec_ga(4)
    real(dp) :: vec_ag_a(4),vec_ag_g(4),vec_ga_a(4),vec_ga_g(4)
    real(dp) :: resvec(4)
    real(dp) :: res_ag(4),res_ga(4)
    real(dp) :: respdf_a(1),respdf_g(1)
    real(dp) :: respdf_ag(1),respdf_ga(1), tag_ns_wga

    real(dp) :: pdec(4,3)

#if(_identification_check == 1)
    logical    :: check, check_temp
    !the debug arrays contains hard and single unresolved contirbutions
    ! ( qqb -> ga  , qqb -> ag  , qqb->qqb , a+S5g , g+S5a a+C5g , g+C5a )
    integer :: debug_binarray(9,nobs)
    real(dp) :: debug_weightarray(9,nobs)
    real(dp) :: debug_kin(9)
    debug_binarray    = zero
    debug_weightarray = zero
    debug_kin         = zero
    check = .true.
#endif

    tag_ns_wga = one 

    callcount = callcount + 1

    xsect_nnlo_5264_rr_z_ga = 0
    ff(1) = zero
    xx = zero


    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    

    if (xx(xE4)*xx(xE5)*xx(xRHO1)*xx(xRHO2).lt.buff_rr) return

    call get_tauy(xx(1:tauy_max),tau,ylab,jac)
    M2 = tau*sh

    call open_histo()    
    
    call kinematics_nnlo_5264(M2,ylab,xx(kNNLO_min:kNNLO_max_full),HardProc,S4Lim,S5Lim,S4S5Lim,C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,&
         C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim)  

    !! -----------------------------!!
    !!    calling NNLO ME           !!
    !! -----------------------------!!

    !1-- regular piece ; gluon-photon

    xi1 = HardProc%PartFrac(1)
    xi2 = HardProc%PartFrac(2)

    pdec(:,1) = HardProc%AmpMom(:,3)
    pdec(:,2) = HardProc%AmpMom(:,4)
    pdec(:,3) = HardProc%AmpMom(:,6)

    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then
       kin(1) = zero

       FintNNLO_dy_ns(1) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)


       call res_tree_ga_qqb(HardProc%AmpMom,resnnlo)

       !call get_respdf(ns_lumi,1,1,HardProc,resnnlo,respdf)
       call get_respdf_ns(1,1,HardProc,resnnlo,respdf)

       FintNNLO_dy_ns(1) = respdf(1)*HardProc%wgt*jac


       kin(1) = FintNNLO_dy_ns(1)


       call fill_histo([kin(1)],vegasweight)


    endif



    !! -----------------------------!!
    !!     S5 + S5/C5          !!
    !! -----------------------------!!

    !-- 2 -> S5   ; photon

    xi1 = S5Lim%PartFrac(1) 
    xi2 = S5Lim%PartFrac(2)

    call cut_histo(S5Lim)
    if (S5Lim%makecut.or.S5Lim%flag) then

       kin(2) = zero
       FintNNLO_dy_ns(2) = zero
       FintNNLO_dy_ns(3) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_g_qqb(S5Lim%AmpMom,res_S5_nlo)
       call get_Eik_ph_z_ga(S5Lim%LimMom(:,1),S5Lim%LimMom(:,2),S5Lim%LimMom(:,3), &
            S5Lim%LimMom(:,4), S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5_nlo,eik,res2)
       !call get_respdf(ns_lumi,1,1,S5Lim,res2,respdf_ga)
       call get_respdf_ns(1,1,S5Lim,res2,respdf_ga)

       FintNNLO_dy_ns(2) = (-one)*respdf_ga(1)*S5Lim%wgt*jac


       !-- 3 -> S5C5 ; photon

       s75 = C5S5Lim%Lim_KinInv(2)
       z5 = C5S5Lim%Lim_KinInv(5)

       !call get_respdf(ns_lumi,1,1,S5Lim,res_S5_nlo,respdf_qcd)
       call get_respdf_ns(1,1,S5Lim,res_S5_nlo,respdf_qcd)

       FintNNLO_dy_ns(3) = four/z5/s75 * respdf_qcd(1) * C5S5Lim%wgt*jac

       kin(2) = FintNNLO_dy_ns(2) + FintNNLO_dy_ns(3)

       call fill_histo([kin(2)],vegasweight)

    endif



    !-- 4 -> C4 ; 

    xi1 = C4Lim%PartFrac(1) 
    xi2 = C4Lim%PartFrac(2)
    call cut_histo(C4Lim)

    if (C4Lim%makecut.or.C4Lim%flag) then
       kin(3) = zero
       FintNNLO_dy_ns(4) = zero  

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_a_qqb(C4Lim%AmpMom,res_a_nlo)
       !call get_respdf(ns_lumi,1,1,C4Lim,res_a_nlo,respdf_b_qed)
       call get_respdf_ns(1,1,C4Lim,res_a_nlo,respdf_b_qed)

       s14 = C4Lim%Lim_KinInv(1) 
       z4  = C4Lim%Lim_KinInv(4)

       FintNNLO_dy_ns(4)  = (-one)*(-two/s14)* Pqq(z4,Cf) * respdf_b_qed(1) *C4Lim%wgt*jac


       kin(3) = FintNNLO_dy_ns(4) 
       call fill_histo([kin(3)],vegasweight)

    endif

    !-- 5 -> C5 ; photon

    xi1 = C5Lim%PartFrac(1) 
    xi2 = C5Lim%PartFrac(2)

    call cut_histo(C5Lim)

    if (C5Lim%makecut.or.C5Lim%flag) then
       kin(4) = zero
       FintNNLO_dy_ns(5)  = zero
    else


       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_g_qqb(C5Lim%AmpMom,res_C5)

       !call get_respdf(ns_lumi,1,1,C5Lim,res_C5,respdf_ch1)
       call get_respdf_ns(1,1,C5Lim,res_C5,respdf_ch1)

       s75 = C5Lim%Lim_KinInv(2) 
       z5  = C5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(5)  = (-two/s75)* Pqq(z5,one)*( respdf_ch1(1) ) * C5Lim%wgt*jac

!!$        print *, ''
!!$
!!$       print*, 'C5'
!!$       print*, '(-two/s75)* Pqq(z5,one)*( respdf_ch1(1) )',(-two/s75)* Pqq(z5,one)*( respdf_ch1(1) )
!!$       print *, ''
!!$
!!$       print*, 's75',s75
!!$       print*, ' Pqq(z5,one)',  Pqq(z5,one)
!!$       print*, 'respdf_ch1(1)', respdf_ch1(1) 
!!$       print *, ''
!!$
!!$       
!!$       print*, 'res_C5(1,1)', res_C5(1,1)
!!$       print*, 'res_C5(1,2)', res_C5(1,2)
!!$       print*, 'res_C5(2,1)', res_C5(2,1)
!!$       print*, 'res_C5(2,2)', res_C5(2,2)
!!$       print *, ''
!!$
!!$       print*, 'C5Lim%AmpMom(:,1)',C5Lim%AmpMom(:,1)
!!$       print*, 'C5Lim%AmpMom(:,2)',C5Lim%AmpMom(:,2)
!!$       print*, 'C5Lim%AmpMom(:,3)',C5Lim%AmpMom(:,3)
!!$       print*, 'C5Lim%AmpMom(:,4)',C5Lim%AmpMom(:,4)
!!$       print*, 'C5Lim%AmpMom(:,5)',C5Lim%AmpMom(:,5)
!!$       
!!$       print *, ''
!!$       print*, 'C5Lim%wgt', C5Lim%wgt
!!$       print*, 'FintNNLO_dy_ns(5)', FintNNLO_dy_ns(5)
!!$
!!$       print *, ''



       kin(4) = FintNNLO_dy_ns(5) 
       call fill_histo([kin(5)],vegasweight)

    endif

    !! -----------------------------!!
    !!     calling LO ME            !!
    !! -----------------------------!!

    !-- 6 -> S4S5
    !-- 7 -> C4DS
    !-- 8 -> C5S4S5
    !-- 9 -> C4C5S4S5

    xi1 = S4S5Lim%PartFrac(1) 
    xi2 = S4S5Lim%PartFrac(2)

    call cut_histo(S4S5Lim)

    if (S4S5Lim%makecut.or.S4S5Lim%flag) then
       kin(5) = zero
       FintNNLO_dy_ns(6:9)  = zero


    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(S4S5Lim%AmpMom,res_S5S4) 

       call get_Eik_ds_z_ga(S4S5Lim%LimMom(:,1),S4S5Lim%LimMom(:,2),S4S5Lim%LimMom(:,3),&
            S4S5Lim%LimMom(:,4),S4S5Lim%LimMom(:,5),S4S5Lim%LimMom(:,6),eik)
       call product_for_pdf(res_S5S4,eik,res2)
       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ag)
       call get_respdf_ns(1,1,S4S5Lim,res2,respdf_ag)

       FintNNLO_dy_ns(6) = (-one)*(respdf_ag(1))*S4S5Lim%wgt*jac

       s14 = C4S4S5Lim%Lim_KinInv(4) 
       z4  = C4S4S5Lim%Lim_KinInv(5) 

       call get_Eik_ph_z_ga(C4S4S5Lim%LimMom(:,1),C4S4S5Lim%LimMom(:,2),C4S4S5Lim%LimMom(:,3), &
            C4S4S5Lim%LimMom(:,4), C4S4S5Lim%LimMom(:,5),eik)
       call product_for_pdf(res_S5S4,eik,res2)  

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res2,respdf_ga)
       call get_respdf_ns(1,1,S4S5Lim,res2,respdf_ga)


       FintNNLO_dy_ns(7) =(-one)*(respdf_ga(1))*four*Cf/s14/z4*C4S4S5Lim%wgt*jac


       !C5S4S5

       s12 = C5S4S5Lim%Lim_KinInv(1) 
       s14 = C5S4S5Lim%Lim_KinInv(2) 
       s24 = C5S4S5Lim%Lim_KinInv(3) 
       s75 = C5S4S5Lim%Lim_KinInv(4) 
       z5  = C5S4S5Lim%Lim_KinInv(5)

       !call get_respdf(ns_lumi,1,1,S4S5Lim,res_S5S4,respdf)
       call get_respdf_ns(1,1,S4S5Lim,res_S5S4,respdf)

       FintNNLO_dy_ns(8) = (-one)*respdf(1) * four/s75/z5 * Cf*four*s12/s14/s24*C5S4S5Lim%wgt*jac


       !C4C5S4S5

       s14 = C4C5S4S5Lim%Lim_KinInv(1) 
       s75 = C4C5S4S5Lim%Lim_KinInv(2) 
       z5 =  C4C5S4S5Lim%Lim_KinInv(5) 
       z4 =  C4C5S4S5Lim%Lim_KinInv(4)

       FintNNLO_dy_ns(9) =  16.0_dp*Cf/s14/z4/s75/z5 &
            * respdf(1)* C4C5S4S5Lim%wgt*jac


       kin(5) = sum(FintNNLO_dy_ns(6:9))


       call fill_histo([kin(5)],vegasweight)

    endif



    !-- 10 -> C4S5
    !-- 11 -> C4C5S5

    xi1 = C4S5Lim%PartFrac(1) 
    xi2 = C4S5Lim%PartFrac(2)

    call cut_histo(C4S5Lim)

    if (C4S5Lim%makecut.or.C4S5Lim%flag) then

       kin(6) = zero
       FintNNLO_dy_ns(10:11)  = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C4S5Lim%AmpMom,res_C4S5)         
       call get_Eik_ph_z_ga(C4S5Lim%LimMom(:,1),C4S5Lim%LimMom(:,2),C4S5Lim%LimMom(:,3), &
            C4S5Lim%LimMom(:,4),C4S5Lim%LimMom(:,5),eik)

       call product_for_pdf(res_C4S5,eik,res_C4S52)


       !call get_respdf(ns_lumi,1,1,C4S5Lim,res_C4S52,respdf)
       call get_respdf_ns(1,1,C4S5Lim,res_C4S52,respdf)

       s14 = C4S5Lim%Lim_KinInv(4) 
       z4 =  C4S5Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(10) = (- respdf(1))*two/s14*Pqq(z4,Cf)* C4S5Lim%wgt*jac


       !call get_respdf(ns_lumi,1,1,C4S5Lim,res_C4S5,respdf_ch1)
       call get_respdf_ns(1,1,C4S5Lim,res_C4S5,respdf_ch1)

       s14 = C4C5S5Lim%Lim_KinInv(1) 
       s75 = C4C5S5Lim%Lim_KinInv(2) 
       z4 =  C4C5S5Lim%Lim_KinInv(4) 
       z5 =  C4C5S5Lim%Lim_KinInv(5) 
       !

       FintNNLO_dy_ns(11) = 8.0_dp/s14/s75/z5 * Pqq(z4,Cf )  *respdf_ch1(1) * C4C5S5Lim%wgt*jac

       kin(6) = sum(FintNNLO_dy_ns(10:11))

       call fill_histo([kin(6)],vegasweight)

    endif


    !-- 12 -> C4C5

    xi1 = C4C5Lim%PartFrac(1) 
    xi2 = C4C5Lim%PartFrac(2) 
    call cut_histo(C4C5Lim)
    if (C4C5Lim%makecut.or.C4C5Lim%flag) then

       FintNNLO_dy_ns(12) = zero
       kin(7) = zero

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C4C5Lim%AmpMom,res_C4C5)

       s14 = C4C5Lim%Lim_KinInv(1) 
       s75 = C4C5Lim%Lim_KinInv(2) 
       z4 =  C4C5Lim%Lim_KinInv(4) 
       z5 =  C4C5Lim%Lim_KinInv(5)

       !call get_respdf(ns_lumi,1,1,C4C5Lim,res_C4C5,respdf)
       call get_respdf_ns(1,1,C4C5Lim,res_C4C5,respdf)


       FintNNLO_dy_ns(12) = (-one)*four/s75/s14 * Pqq(z5,esq) * Pqq(z4,Cf )  *respdf(1) * C4C5Lim%wgt*jac



       kin(7) = FintNNLO_dy_ns(12)

       call fill_histo([kin(7)],vegasweight)

    endif



    !! -----------------------------!!
    !!     S4 + S4/C4         !!
    !! -----------------------------!!

    !--  -> S4   ; gluon
    !--  -> S4C4 ; gluon

    xi1 = S4Lim%PartFrac(1) 
    xi2 = S4Lim%PartFrac(2) 
    call cut_histo(S4Lim)


    if (S4Lim%makecut.or.S4Lim%flag) then
       kin(8) = zero
       FintNNLO_dy_ns(13) = zero
       FintNNLO_dy_ns(14) = zero
    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_a_qqb(S4Lim%AmpMom,res_a_nlo)

       !call get_respdf(ns_lumi,1,1,S4Lim,res_a_nlo,respdf)
       call get_respdf_ns(1,1,S4Lim,res_a_nlo,respdf)

       reseik_qcd = myEik_g(S4Lim%LimMom(:,1),S4Lim%LimMom(:,2),S4Lim%LimMom(:,3),Cf)


       FintNNLO_dy_ns(13) = (-one)*respdf(1) *reseik_qcd*S4Lim%wgt*jac


       s14 = C4S4Lim%Lim_KinInv(2)
       z4 = C4S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(14) = four*Cf/z4/s14 * respdf(1) * C4S4Lim%wgt*jac


       kin(8) = FintNNLO_dy_ns(13) + FintNNLO_dy_ns(14)

       call fill_histo([kin(8)],vegasweight)


    endif


    ! - - - - - - C5S4 and  C4C5S4 ; photon

    xi1 = C5S4Lim%PartFrac(1) 
    xi2 = C5S4Lim%PartFrac(2) 
    call cut_histo(C5S4Lim)

    if (C5S4Lim%makecut.or.C5S4Lim%flag) then



       kin(9) = zero
       FintNNLO_dy_ns(15) = zero
       FintNNLO_dy_ns(16) = zero 

    else

       !call setscale(mu,muF,muR)
       !call get_pdf(xi1,xi2,muF,f1,f2)
       !call get_prefactor(muR,2,pref)

       call res_tree_qqb(C5S4Lim%AmpMom,res_C5)      
       !call get_respdf(ns_lumi,1,1,C5S4Lim,res_C5,respdf_ch1)
       call get_respdf_ns(1,1,C5S4Lim,res_C5,respdf_ch1)


       s12 = C5S4Lim%Lim_KinInv(1)
       s14 = C5S4Lim%Lim_KinInv(2)
       s24 = C5S4Lim%Lim_KinInv(3) 
       s75 = C5S4Lim%Lim_KinInv(4) 
       z5  = C5S4Lim%Lim_KinInv(5)

       FintNNLO_dy_ns(15)  = respdf_ch1(1)*two*Cf*s12/s14/s24*four/s75*Pqq(z5,one) * C5S4Lim%wgt*jac

       s14 = C4C5S4Lim%Lim_KinInv(1)
       s75 = C4C5S4Lim%Lim_KinInv(2)


       z4 = C4C5S4Lim%Lim_KinInv(4)
       z5 = C4C5S4Lim%Lim_KinInv(5)


       FintNNLO_dy_ns(16)  = (-one)*respdf_ch1(1)*four*Cf/s75/z4*two/s14*Pqq(z5,one) * C4C5S4Lim%wgt*jac


       kin(9) = FintNNLO_dy_ns(15)+FintNNLO_dy_ns(16)


    endif

    call fill_histo([kin(9)],vegasweight)

    ff(1) = sum(kin)
    call close_histo()
    
    !  print*, 'xx=', xx

    !     print  *, 1, FintNNLO_dy_ns(1), 'hard'
    !     print  *, '------------------------------------------------------'
    !     print  *, 2, FintNNLO_dy_ns(13)/FintNNLO_dy_ns(1), 'x1=0'
    !     print  *, 3, FintNNLO_dy_ns(2)/FintNNLO_dy_ns(1), 'x2=0'
    !     print  *, 4, FintNNLO_dy_ns(4)/FintNNLO_dy_ns(1), 'x3=0'
    !     print  *, 5, FintNNLO_dy_ns(5)/FintNNLO_dy_ns(1), 'x4=0'
    !     print  *, '------------------------------------------------------'
    !     print  *, 6, FintNNLO_dy_ns(6)/FintNNLO_dy_ns(1),    'x1=0,x2=0'
    !     print  *, 7, FintNNLO_dy_ns(14)/FintNNLO_dy_ns(1), 'x1=0,x3=0'     
    !     print  *, 8, FintNNLO_dy_ns(15)/FintNNLO_dy_ns(1), 'x1=0,x4=0'
    !     print  *, 9, FintNNLO_dy_ns(10)/FintNNLO_dy_ns(1), 'x2=0, x3=0'
    !     print  *, 10, FintNNLO_dy_ns(3)/FintNNLO_dy_ns(1), 'x2=0,x4=0'
    !     print  *, 11, FintNNLO_dy_ns(12)/FintNNLO_dy_ns(1), 'x3=0,x4=0'
    !     print  *, '------------------------------------------------------'     
    !     print  *, 12, FintNNLO_dy_ns(7)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0'
    !     print  *, 13, FintNNLO_dy_ns(8)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x4=0'
    !     print  *, 14, FintNNLO_dy_ns(16)/FintNNLO_dy_ns(1), 'x1=0,x3=0,x4=0'  
    !     print  *, 15, FintNNLO_dy_ns(11)/FintNNLO_dy_ns(1), 'x2=0,x3=0,x4=0'
    !     print  *, '------------------------------------------------------'     
    !     print  *, 16, FintNNLO_dy_ns(9)/FintNNLO_dy_ns(1), 'x1=0,x2=0,x3=0,x4=0'
    !     print*, 'ff(1)', ff(1)

    !     pause
    ! !!$
!!$     if (ff(1) .ne. ff(1)) then
!!$       print  *, 'nan'
!!$       print  *, 'xx',xx
!!$       print  *, ''
!!$       do i = 1,16
!!$          print  *, i, FintNNLO_dy_ns(i)
!!$       enddo
!!$
!!$       pause
!!$       
!!$    endif
!!$
!!$    pause
!!$    

!!$    pause

!!$    print *, ''
!!$    print*, 'kin(1)', kin(1)
!!$    print*, 'kin(2)', kin(2)
!!$    print*, 'kin(3)', kin(3)
!!$    print*, 'kin(4)', kin(4)
!!$    print*, 'kin(5)', kin(5)
!!$    print*, 'kin(6)', kin(6)
!!$    print*, 'kin(7)', kin(7)
!!$    print*, 'kin(8)', kin(8)
!!$    print*, 'kin(9)', kin(9)
!!$    print*

    !print*, pref

    !    print*, '5264 ff(1)', ff(1)

    !    pause


    !   print*, 'ff(1)', ff(1)
    !   pause
!!$
!!$
!!$
!!$
!!$
!!$    
!!$    if (ff(1) .ne. ff(1)) then
!!$       !print  *, 'nan'
!!$       !print  *, 'xx',xx
!!$       !print  *, ''
!!$       do i = 1,16
!!$          !print  *, i, FintNNLO_dy_ns(i)
!!$       enddo
!!$       ff(1) = zero
!!$    endif
!!$#if(_identification_check == 1)  
!!$    endif
!!$#if(_identification_check == 1)  
!!$    if( .not. check ) stop
!!$#endif
!!$
!!$#if (_DebugFeatures == 1)
!!$    FintNNLO_CheckLim_dy_ns(1) = FintNNLO_dy_ns(1)
!!$#endif
!!$#if (_DebugFeatures == 1)
!!$    FintNNLO_CheckLim_dy_ns(4:16) = FintNNLO_dy_ns(4:16)
!!$#endif
!!$
!!$    ! do i = 1,16
!!$    !    !print  *, i,FintNNLO_dy_ns(i)
!!$    ! enddo
!!$    ! !print  *, 'ff(1)', ff(1)
!!$    ! pause


    if (ff(1) .ne. ff(1)) then
       print  *, 'nan'
       print  *, 'xx',xx
       print  *, ''
       do i = 1,16
          print  *, i, FintNNLO_dy_ns(i)
       enddo
       ff(1)=zero
    endif

#if(_withchecks == 1)
    FintNNLO = FintNNLO_dy_ns
#endif

    return

  end function xsect_nnlo_5264_rr_z_ga

end module mod_xsects_nnlo_rr_z_ga



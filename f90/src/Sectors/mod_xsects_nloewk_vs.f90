module mod_xsects_nloewk_vs
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_kinematics_lo
  use mod_process
  use mod_aux_sectors
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_amplitudes_loop_ppll
  use mod_hoppet_tools
  use mod_hoppet_nlo
  use, intrinsic :: ieee_arithmetic
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNLOEWK_vs(1)
#endif

  private

  public :: xsect_nloewk_v_ns ![WORK IN PROGRESS]
  public :: xsect_nloewk_s_ns ![WORK IN PROGRESS]

  public :: xsect_nloewk_v_aa
  public :: xsect_nloewk_s_aa
  public :: xsect_nloewk_s_aq,xsect_nloewk_s_qa

contains

  function xsect_nloewk_v_ns(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_v_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_tree(-5:7,-5:7),res_loop(-5:7,-5:7)
    real(dp)    :: res_tmp(-5:7,-5:7)
    !--
    real(dp)    :: res_tree_old(2,3),res_loop_old(2,3)
    real(dp)    :: res_loop_old_tmp(2,3)
    integer :: i, j 
    logical :: oldcode

    oldcode = .true.

    xsect_nloewk_v_ns = 0
    res_loop(:,:) = 0
    res_tree(:,:) = 0

    ff(1) = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()

    call kinematics_lo(xx,LOProc)

    ! define process specific partons
#if (_Vcharge == 0)
    LOProc%part(1:4) = [id_q,-id_q,id_el,-id_el]
#elif  (_Vcharge == -1)
    LOProc%part(1:4) = [id_q,-id_qp,id_el,-id_nue]
#elif  (_Vcharge == +1)
    LOProc%part(1:4) = [id_q,-id_qp,id_nue,-id_el]
#endif

    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero
       
    else    
    
       if (oldcode) then
          call res_ewkloop_qqb(LOProc%AmpMom,res_tree_old,res_loop_old)
          call get_respdf(ns_lumi_splitb,0,1,LOProc,res_loop_old,respdf)
       else
          call res_ewkloop_qqb_gen(LOProc%AmpMom,res_tree,res_loop)
          call get_respdf_gen(0,1,LOProc,res_loop,respdf)
       endif 
       
       respdf = respdf*LOProc%wgt
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOEWK_vs = kin
#endif

  end function xsect_nloewk_v_ns

  function xsect_nloewk_v_aa(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_v_aa
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_tree,res_loop(1,1)

    xsect_nloewk_v_aa = 0

    ff(1) = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()

    call kinematics_lo(xx,LOProc)
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero

    else

       call res_ewkloopAA_aa(LOProc%AmpMom,res_tree,res_loop(1,1))
       call get_respdf(aa_lumi,0,1,LOProc,res_loop,respdf)

       respdf = respdf*LOProc%wgt
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOEWK_vs = kin
#endif

  end function xsect_nloewk_v_aa

  !-- subtractions below
  
  function xsect_nloewk_s_ns(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_s_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf),respdf_3(ipdf)
    real(dp)    :: res_lo_old(2,2),res_lo_tmp_old(2,2),bit1,bit2,eta13,eta14
    !--
    real(dp)    :: res_lo(-5:7,-5:7), res_lo_tmp(-5:7,-5:7), res_lo_tmp1(-5:7,-5:7), res_lo_tmp2(-5:7,-5:7)
    !--
    integer, parameter :: n = 4 ! number of momenta
    real(dp) :: eta(n,n), PolyLogij(n,n), Lij(n,n) 
    real(dp) :: fin_elasticZ_dn, fin_elasticZ_up, fin_elasticW_du, fin_elasticW_ud
    real(dp) :: fin_elasticZ_ddx, fin_elasticZ_uux, fin_elasticZ_dxd, fin_elasticZ_uxu
    real(dp) :: fin_elasticW_dxu, fin_elasticW_uxd, fin_elasticW_dux

    real(dp) :: Emax
    integer, parameter :: nf = 5
    real*8 :: charge_factor(nf)
    real*8 :: factor
    integer :: a, b, fa, fb
    real*8 :: Qcharge(nf)
    real*8 :: Q2(nf)

    integer :: i, j
    logical :: oldcode

    oldcode = .false.
    
    Qcharge = (/ Qdn, Qup, Qdn, Qup, Qdn /)
    Q2 = Qcharge**2
    
    res_lo_old(:,:) = 0
    res_lo_tmp_old(:,:) = 0
    
    res_lo(:,:) = zero
    res_lo_tmp(:,:) = zero
    res_lo_tmp1(:,:) = zero
    res_lo_tmp2(:,:) = zero 
    
    xsect_nloewk_s_ns = 0

    ff(1) = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

!    xx = (/ 5.7767012882822426D-002, &
!        0.24444789910714418D0, &
!        0.94960547288538688D0, &
!        1.3608257177144900D-002, &
!        3.5008535261739970D-002 /)

    call open_histo()

    call kinematics_lo(xx,LOProc)
       ! define process specific partons
#if (_Vcharge == 0)
    LOProc%part(1:4) = [id_q,-id_q,id_el,-id_el]
#elif  (_Vcharge == -1)
    LOProc%part(1:4) = [id_q,-id_qp,id_el,-id_nue]
#elif  (_Vcharge == +1)
    LOProc%part(1:4) = [id_q,-id_qp,id_nue,-id_el]
#endif

    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero
       
    else    
        
        if (oldcode) then 

             call res_tree_qqb(LOProc%AmpMom,res_lo_old)
             res_lo_tmp_old(:,1) = res_lo_old(:,1) * Qdn2
             res_lo_tmp_old(:,2) = res_lo_old(:,2) * Qup2

             !-- z-dependent bit
             call get_respdf_hoppet(xPij,PDFs,ns_lumi,0,1,LOProc,res_lo_tmp_old,respdf_1,myPDFs1_Lmu=[xPij_Lmu])
             respdf_1 = respdf_1*LOProc%wgt
             
             call get_respdf_hoppet(PDFs,xPij,ns_lumi,0,1,LOProc,res_lo_tmp_old,respdf_2,myPDFs2_Lmu=[xPij_Lmu])
             respdf_2 = respdf_2*LOProc%wgt

             !-- FLM[1,2] bit, assuming Emax = sqrt(q2)/2
             eta13 = half*scr(LOProc%AmpMom(:,1),LOProc%AmpMom(:,3))/LoProc%AmpMom(1,1)/LOProc%AmpMom(1,3)
             eta14 = half*scr(LOProc%AmpMom(:,1),LOProc%AmpMom(:,4))/LoProc%AmpMom(1,1)/LOProc%AmpMom(1,4)
             bit1 = 13._dp - two/three*pisq
             bit2 = log(eta13/eta14)*(three/two)+real(dilog2(one-eta13),kind=dp)-real(dilog2(one-eta14),kind=dp)

             res_lo_tmp_old(1,:) = (Q_lep2 * bit1 + four*Q_lep*[Qdn,Qup] * bit2)*res_lo_old(1,:)
             res_lo_tmp_old(2,:) = (Q_lep2 * bit1 - four*Q_lep*[Qdn,Qup] * bit2)*res_lo_old(2,:)

             call get_respdf(ns_lumi,0,1,LOProc,res_lo_tmp_old,respdf_3)
             
             respdf_3 = respdf_3*LOproc%wgt
       
       else 
             call res_tree_qqb_gen(LOProc%AmpMom,res_lo)

             do a = -nf, nf
               if (a == 0) cycle
               fa = abs(a)
               do b = -nf, nf
                 if (b == 0) cycle
                 fb = abs(b)
                 res_lo_tmp1(a,b) = Q2(fa) * res_lo(a,b)
                 res_lo_tmp2(a,b) = Q2(fb) * res_lo(a,b)
               end do
             end do

             !-- z-dependent bit
             call get_respdf_hoppet_gen(xPij,PDFs,0,1,LOProc,res_lo_tmp1,respdf_1,myPDFs1_Lmu=[xPij_Lmu])
!             print*, 'respdf_1= ', respdf_1
!             print*, 'LOProc%wgt= ', LOProc%wgt  

             respdf_1 = respdf_1*LOProc%wgt  

!             print*, 'respdf_1= ', respdf_1

             call get_respdf_hoppet_gen(PDFs,xPij,0,1,LOProc,res_lo_tmp2,respdf_2,myPDFs2_Lmu=[xPij_Lmu])
             respdf_2 = respdf_2*LOProc%wgt   
             
!             print*, 'respdf_2= ', respdf_2

             !-- FLM[1,2] bit, assuming E1 = E2 = E3 = E4 = Emax
             Emax = LOProc%AmpMom(1,1)
             eta = get_eta(LOProc%AmpMom(:,:), n)    
             do i = 1, n
                do j = 1, n
                if (i .eq. j) cycle
                Lij(i,j) = log(eta(i,j))
                PolyLogij(i,j) = real(dilog2(one-eta(i,j)),kind=dp)
                enddo
             enddo

#if (_Vcharge == 0)
             fin_elasticZ_ddx = get_subtra_elastic_ewk_qqbllb_gen(1,2,3,4,Lij,Emax,mu,PolyLogij,Qdn,-Qdn,Q_lep,-Q_lep)
             fin_elasticZ_uux = get_subtra_elastic_ewk_qqbllb_gen(1,2,3,4,Lij,Emax,mu,PolyLogij,Qup,-Qup,Q_lep,-Q_lep)
             fin_elasticZ_dxd = get_subtra_elastic_ewk_qqbllb_gen(1,2,3,4,Lij,Emax,mu,PolyLogij,-Qdn,Qdn,Q_lep,-Q_lep)
             fin_elasticZ_uxu = get_subtra_elastic_ewk_qqbllb_gen(1,2,3,4,Lij,Emax,mu,PolyLogij,-Qup,Qup,Q_lep,-Q_lep)
             ! down-type (1,3,5)
             do i = 1,5,2
                res_lo(i,-i) = res_lo(i,-i) * fin_elasticZ_ddx
                res_lo(-i,i) = res_lo(-i,i) * fin_elasticZ_dxd
             end do
             ! up-type (2,4)
             do i = 2,4,2
                res_lo(i,-i) = res_lo(i,-i) * fin_elasticZ_uux
                res_lo(-i,i) = res_lo(-i,i) * fin_elasticZ_uxu
             end do

#elif (_Vcharge == +1)
             fin_elasticW_dxu = get_subtra_elastic_ewk_qqbllb_gen(1,2,3,4,Lij,Emax,mu,PolyLogij,-Qdn,Qup,0d0,-Q_lep)
             fin_elasticW_udx = get_subtra_elastic_ewk_qqbllb_gen(1,2,3,4,Lij,Emax,mu,PolyLogij,Qup,-Qdn,0d0,-Q_lep)
             ! u dbar  → fin_elasticW_udx
             do i = 2,4,2              ! up-type: 2,4
                j = i - 1              ! corresponding down-type: 1,3
             res_lo(i,-j) = res_lo(i,-j) * fin_elasticW_udx
             end do
             ! dbar u  → fin_elasticW_dxu
             do j = 1,3,2              ! down-type: 1,3
                i = j + 1              ! corresponding up-type: 2,4
                res_lo(-j,i) = res_lo(-j,i) * fin_elasticW_dxu
             end do

#elif (_Vcharge == -1)
             fin_elasticW_dux = get_subtra_elastic_ewk_qqbllb_gen(1,2,3,4,Lij,Emax,mu,PolyLogij,Qdn,-Qup,Q_lep,0d0)
             fin_elasticW_uxd = get_subtra_elastic_ewk_qqbllb_gen(1,2,3,4,Lij,Emax,mu,PolyLogij,-Qup,Qdn,Q_lep,0d0)
             ! d ubar → fin_elasticW_dux
             do j = 1,3,2           ! down-type: 1,3
                i = j + 1           ! corresponding up-type: 2,4
                res_lo(j,-i) = res_lo(j,-i) * fin_elasticW_dux
             end do
             ! ubar d → fin_elasticW_uxd
             do i = 2,4,2           ! up-type: 2,4
                j = i - 1           ! corresponding down-type: 1,3
                res_lo(-i,j) = res_lo(-i,j) * fin_elasticW_uxd
             end do

#endif
             call get_respdf_gen(0,1,LOProc,res_lo,respdf_3)
             respdf_3 = respdf_3*LOproc%wgt
             
       endif 

       respdf = respdf_1 + respdf_2 + respdf_3
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNLOEWK_vs = kin
#endif


  end function xsect_nloewk_s_ns




  !-----------------------------------------------------------------
  !--- generic routines
  !-----------------------------------------------------------------

  function get_eta(p, n) result(eta)
  implicit none
  integer,  intent(in)  :: n
  real(dp), intent(in)  :: p(:,:)          ! expected shape (dim, n)
  real(dp)              :: eta(n,n)

  integer :: i, j

  do i = 1, n
     do j = i, n
        eta(i,j) = scr(p(:,i), p(:,j)) / ( p(1,i) * p(1,j) * two )
        eta(j,i) = eta(i,j)
     end do
  end do

  end function get_eta


  function get_subtra_elastic_ewk_qqbllb_gen(i1,i2,i3,i4,Lij,Emax,mu,PolyLogij,Qq,Qqbp,Ql,Qlb) result(res)
  implicit none

  integer, intent(in)   :: i1,i2,i3,i4
  real(dp), intent(in)  :: Lij(:,:)
  real(dp), intent(in)  :: Emax, mu
  real(dp), intent(in)  :: PolyLogij(:,:)
  real(dp) :: res

  real(dp) :: Qq, Qqbp, Ql, Qlb
  real(dp) :: logETAbit,logENERGYbit,polylogbit,constantbit

  !---------------------------------------------------
  ! Compute logarithmic pieces
  !---------------------------------------------------
  ! Log[eta_ij]
  logETAbit  = ewk_log_ETA_bit(Qq, Qqbp, Ql, Qlb, i1, i2, i3, i4, Lij)
  
  ! Log[2EC/mu] -> actually this is always proportional to the sum on the charges
  ! which by charge conservation is zero
  !logENERGYbit = ewk_log_ENERGY_bit(Qq, Qqbp, Ql, Qlb, i1, i2, i3, i4, Emax, mu)
  
  ! PolyLog[1-eta_ij]
  polylogbit = ewk_polylog_bit(Qq, Qqbp, Ql, Qlb, i1, i2, i3, i4, PolyLogij)
  ! pi^2 + const
  constantbit = ewk_constant_bit(Qq, Qqbp, Ql, Qlb)

  res = 3.0_dp*logETAbit &
          !+ logENERGYbit 
  + two*polylogbit + constantbit 

  !print*, 'constantbit= ', constantbit
  !print*, 'logENERGYbit= ', logENERGYbit
  !print*, 'polylogbit= ', two*polylogbit
  !print*, 'logETAbit= ', 3.0_dp*logETAbit
contains

  !---------------------------------------------------
  ! Log-eta contribution
  !---------------------------------------------------
  function ewk_log_ETA_bit(Qq, Qqbp, Ql, Qlb, i1, i2, i3, i4, Lij) result(val)
    implicit none
    real(dp), intent(in) :: Qq, Qqbp, Ql, Qlb
    integer, intent(in)  :: i1, i2, i3, i4
    real(dp), intent(in) :: Lij(:,:) ! eta[i,j]
    real(dp)             :: val

    val = - Qq *Qqbp*Lij(i1,i2) + Ql *Qq  *Lij(i1,i3) &
          + Qlb*Qq  *Lij(i1,i4) + Ql *Qqbp*Lij(i2,i3) &
          + Qlb*Qqbp*Lij(i2,i4) - Ql *Qlb *Lij(i3,i4)

  end function ewk_log_ETA_bit

  !---------------------------------------------------
  ! Log-energy contribution
  !---------------------------------------------------
  function ewk_log_ENERGY_bit(Qq, Qqbp, Ql, Qlb, i1, i2, i3, i4, Emax, mu) result(val)
    implicit none
    real(dp), intent(in) :: Qq, Qqbp, Ql, Qlb
    integer, intent(in)  :: i1, i2, i3, i4
    real(dp), intent(in) :: Emax, mu
    real(dp)             :: logtwoEConMu
    real(dp)             :: val

    logtwoEConMu= Log(two*Emax/mu)

    val = + 2.0_dp * (Ql + Qlb - Qq - Qqbp)**2 * logtwoEConMu**2  &
          - 3.0_dp * (Ql + Qlb - Qq - Qqbp)**2 * logtwoEConMu

  end function ewk_log_ENERGY_bit

  !---------------------------------------------------
  ! Constant contribution
  !---------------------------------------------------
  function ewk_constant_bit(Qq, Qqbp, Ql, Qlb) result(val)
    implicit none
    real(dp), intent(in) :: Qq, Qqbp, Ql, Qlb
    real(dp)             :: val

    val = 13.0_dp/2.0_dp * Ql**2 + 13.0_dp/2.0_dp * Qlb**2  &
          - ( 5.0_dp/6.0_dp * Ql**2 + 5.0_dp/6.0_dp * Qlb**2  &
          + 1.0_dp/6.0_dp * Qq**2 + 1.0_dp/6.0_dp * Qqbp**2 &
          + Ql * Qlb + Qq * Qqbp ) * pisq

  ! Important: in the convolutions by hoppet a factor two*zeta2
  !            per each initial state spitting is already included
  !            see line 136 in src/IntSub/mod_hoppet_nlo.f90 
  !            we need to remove it here to avoid double counting

  val = val - (two*zeta2) * Qq**2 - (two*zeta2) * Qqbp**2
  
  end function ewk_constant_bit

  !---------------------------------------------------
  ! PolyLog contribution
  !---------------------------------------------------
  function ewk_polylog_bit(Qq, Qqbp, Ql, Qlb, i1, i2, i3, i4, PolyLogij) result(val)
    implicit none
    real(dp), intent(in) :: Qq, Qqbp, Ql, Qlb
    integer, intent(in)  :: i1, i2, i3, i4
    real(dp), intent(in) :: PolyLogij(:,:) ! PolyLog(1-eta[i,j])
    real(dp)             :: val

    val = - Qq *Qqbp*PolyLogij(i1,i2) + Ql *Qq  *PolyLogij(i1,i3) &
          + Qlb*Qq  *PolyLogij(i1,i4) + Ql *Qqbp*PolyLogij(i2,i3) &
          + Qlb*Qqbp*PolyLogij(i2,i4) - Ql *Qlb *PolyLogij(i3,i4)

  end function ewk_polylog_bit

end function get_subtra_elastic_ewk_qqbllb_gen

















  !--

  function xsect_nloewk_s_aa(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_s_aa
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf),aa_int_sub(ipdf),sv_logs(ipdf)
    real(dp)    :: res_lo(1,1),eta13,eta14,EC

    xsect_nloewk_s_aa = 0

    ff(1) = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()

    call kinematics_lo(xx,LOProc)
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero
       
    else    

       call res_treeAA_aa(LOProc%AmpMom,res_lo(1,1))
       call get_respdf(aa_lumi,0,1,LOProc,res_lo,respdf)

       EC = LOProc%Lim_Ei(1)
       call fill_sv_logs(LOProc%muf(1)**2,4*EC**2,sv_logs)
       eta13 = LOProc%Lim_etaij(1,3)
       eta14 = LOProc%Lim_etaij(1,4)
       aa_int_sub = (13._dp - two/three*pisq)*Q_lep**2 - 2*sv_logs * gamma_a

       respdf = aa_int_sub * respdf * LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOEWK_vs = kin
#endif

  end function xsect_nloewk_s_aa
  
  !--

  function xsect_nloewk_s_aq(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_s_aq
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf_vect(1:2,ipdf),respdf(ipdf)
    real(dp)    :: res_lo(2,2)

    xsect_nloewk_s_aq = 0

    ff(1) = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()

    call kinematics_lo(xx,LOProc)
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero
       
    else    

       call res_tree_qqb(LOProc%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       !-- Pqa*sigma_qq, ns_lumi because the convolution is done in xPij
       call get_respdf_hoppet(xPij,PDFs,ns_lumi,0,1,LOProc,res_lo,respdf_vect(1,:),myPDFs1_Lmu=[xPij_Lmu])
       respdf_vect(1,:) = xn*respdf_vect(1,:)*LOProc%wgt

       !--

       !-- Paq*sigma_aa, quark charges already in the luminosity
       call res_treeAA_aa(LOProc%AmpMom,res_lo(1,1))
       call get_respdf_hoppet(PDFs,xPij_2,aa_lumi,0,1,LOProc,res_lo,respdf_vect(2,:),myPDFs2_Lmu=[xPij_2_Lmu])
       respdf_vect(2,:) = respdf_vect(2,:)*LOProc%wgt

       respdf = sum(respdf_vect(1:2,:),1)
       
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOEWK_vs = kin
#endif

  end function xsect_nloewk_s_aq

  function xsect_nloewk_s_qa(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_s_qa
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf_vect(2,ipdf),respdf(ipdf)
    real(dp)    :: res_lo(2,2)

    xsect_nloewk_s_qa = 0

    ff(1) = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()

    call kinematics_lo(xx,LOProc)
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero
       
    else    

       call res_tree_qqb(LOProc%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       !-- Pqa*sigma_qq, ns_lumi because the convolution is done in xPij
       call get_respdf_hoppet(PDFs,xPij,ns_lumi,0,1,LOProc,res_lo,respdf_vect(1,:),myPDFs2_Lmu=[xPij_Lmu])
       respdf_vect(1,:) = xn*respdf_vect(1,:)*LOProc%wgt

       !--

       !-- Paq*sigma_aa, quark charges already in the luminosity
       call res_treeAA_aa(LOProc%AmpMom,res_lo(1,1))
       call get_respdf_hoppet(xPij_2,PDFs,aa_lumi,0,1,LOProc,res_lo,respdf_vect(2,:),myPDFs1_Lmu=[xPij_2_Lmu])
       respdf_vect(2,:) = respdf_vect(2,:)*LOProc%wgt

       respdf = sum(respdf_vect(1:2,:),1)

       kin(1) = respdf(1)
      
       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNLOEWK_vs = kin
#endif

  end function xsect_nloewk_s_qa
  
end module mod_xsects_nloewk_vs
  

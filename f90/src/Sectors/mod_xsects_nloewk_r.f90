module mod_xsects_nloewk_r
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_limvals
  use mod_proc_parms
  use mod_kinematics_nlo
  use mod_process
  use mod_auxfunctions
  use mod_aux_sectors
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_amplitudes_tree_ppnul
  use mod_partitions
  use mod_splittings_bare
  use mod_eikonals
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNLO_ew(6)
#endif
  real(dp), parameter :: charges_ns(4,4) = reshape(& 
       [-Qdn,   Qdn,  -Qup,   Qup,   &
         Qdn,  -Qdn,   Qup,  -Qup,   &
         Q_lep, Q_lep, Q_lep, Q_lep, &
        -Q_lep,-Q_lep,-Q_lep,-Q_lep]  &
        , [4,4])

  private

  public :: xsect_nloewk_r_is_ns
  public :: xsect_nloewk_r_fs_53_ns, xsect_nloewk_r_fs_54_ns
  public :: xsect_nloewk_r_fs_53_aa, xsect_nloewk_r_fs_54_aa

  public :: xsect_nloewk_r_is_aq, xsect_nloewk_r_is_qa

  !CB: W xsect
  public :: xsect_nloewk_r_is_ns_wp ! deprecated -> delete
  public :: xsect_nloewk_r_is_ns_wm ! deprecated -> delete


contains

  !-------------------------------------------------------------
  !-- initial-state sector
  !-------------------------------------------------------------

  function xsect_nloewk_r_is_ns(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_r_is_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim,SLim,SC1Lim,SC2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(6),kin(4)
    real(dp)    :: respdf(ipdf),respdf_tmp(ipdf,3)
    real(dp)    :: res_nlo_old(2,2),res_lo_old(2,2),res_lo_tmp_old(2,2)
    real(dp)    :: res_nlo(-5:7,-5:7),res_lo(-5:7,-5:7),res_lo_tmp(-5:7,-5:7)
    real(dp)    :: z,s5i,eik(4),e5,eta5i,eik_gen(4,4)
    real(dp)    :: damp
    logical     :: oldcode

    xsect_nloewk_r_is_ns = 0

    oldcode = .false.

    ff(1) = zero
    limval_nlo_is = zero
    
    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_r)) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,&
         C1Lim=C1Lim,C2Lim=C2Lim,SLim=SLim,SC1Lim=SC1Lim,SC2Lim=SC2Lim,compute_etas=.true.)
    
    ! define process specific partons
#if (_Vcharge == 0)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    C1Lim%ids(1:4)    = [0,0,id_el,-id_el]
    C2Lim%ids(1:4)    = [0,0,id_el,-id_el]
    SLim%ids(1:4)     = [0,0,id_el,-id_el]
    SC1Lim%ids(1:4)   = [0,0,id_el,-id_el]
    SC2Lim%ids(1:4)   = [0,0,id_el,-id_el]

    HardProc%part = [id_q,-id_q,id_el,-id_el,id_a]
    C1Lim%part    = [id_q,-id_q,id_el,-id_el]
    C2Lim%part    = [id_q,-id_q,id_el,-id_el]
    SLim%part     = [id_q,-id_q,id_el,-id_el]
    SC1Lim%part   = [id_q,-id_q,id_el,-id_el]
    SC2Lim%part   = [id_q,-id_q,id_el,-id_el]

#elif  (_Vcharge == -1)
    HardProc%ids(1:5) = [0,0,id_el,-id_nue,id_a]
    C1Lim%ids(1:4)    = [0,0,id_el,-id_nue]
    C2Lim%ids(1:4)    = [0,0,id_el,-id_nue]
    SLim%ids(1:4)     = [0,0,id_el,-id_nue]
    SC1Lim%ids(1:4)   = [0,0,id_el,-id_nue]
    SC2Lim%ids(1:4)   = [0,0,id_el,-id_nue]

    HardProc%part = [id_q,-id_qp,id_el,-id_nue,id_a]
    C1Lim%part    = [id_q,-id_qp,id_el,-id_nue]
    C2Lim%part    = [id_q,-id_qp,id_el,-id_nue]
    SLim%part     = [id_q,-id_qp,id_el,-id_nue]
    SC1Lim%part   = [id_q,-id_qp,id_el,-id_nue]
    SC2Lim%part   = [id_q,-id_qp,id_el,-id_nue]

#elif  (_Vcharge == +1)
    HardProc%ids(1:5) = [0,0,id_nue,-id_el,id_a]
    C1Lim%ids(1:4)    = [0,0,id_nue,-id_el]
    C2Lim%ids(1:4)    = [0,0,id_nue,-id_el]
    SLim%ids(1:4)     = [0,0,id_nue,-id_el]
    SC1Lim%ids(1:4)   = [0,0,id_nue,-id_el]
    SC2Lim%ids(1:4)   = [0,0,id_nue,-id_el]

    HardProc%part = [id_q,-id_qp,id_nue,-id_el,id_a]
    C1Lim%part    = [id_q,-id_qp,id_nue,-id_el]
    C2Lim%part    = [id_q,-id_qp,id_nue,-id_el]
    SLim%part     = [id_q,-id_qp,id_nue,-id_el]
    SC1Lim%part   = [id_q,-id_qp,id_nue,-id_el]
    SC2Lim%part   = [id_q,-id_qp,id_nue,-id_el]
#endif

    !-- Hard
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_ns(1) = zero

    else    

       if (oldcode) then
          call res_tree_a_qqb(HardProc%AmpMom,res_nlo_old)       
          call get_respdf(ns_lumi,0,1,HardProc,res_nlo_old,respdf)
       else
          call res_tree_a_qqb_gen(HardProc%AmpMom,res_nlo)
          call get_respdf_gen(0,1,HardProc,res_nlo,respdf)
       endif

       call partition_nlo_qed(HardProc,damp,12)
       respdf = respdf*HardProc%wgt*damp


       kin(1) = respdf(1)

       FintNLO_ns(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintNLO_ns(2) = zero

    else

       if (oldcode) then
          call res_tree_qqb(C1Lim%AmpMom,res_lo_old)
          res_lo_old(1,:) = [Qdn2,Qup2] * res_lo_old(1,:)
          res_lo_old(2,:) = [Qdn2,Qup2] * res_lo_old(2,:)
          call get_respdf(ns_lumi,0,1,C1Lim,res_lo_old,respdf)
       else
          call res_tree_qqb_gen(C1Lim%AmpMom,res_lo)
          res_lo = multiply_IS_charges_sq(res_lo,1)
          call get_respdf_gen(0,1,C1Lim,res_lo,respdf)
       endif


       z   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z)/(one-z))&
            * C1Lim%wgt

       FintNLO_ns(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(3) = zero
       FintNLO_ns(3) = zero

    else

       if (oldcode) then
          call res_tree_qqb(C2Lim%AmpMom,res_lo_old)
          res_lo_old(1,:) = [Qdn2,Qup2] * res_lo_old(1,:)
          res_lo_old(2,:) = [Qdn2,Qup2] * res_lo_old(2,:)
          call get_respdf(ns_lumi,0,1,C2Lim,res_lo_old,respdf)
       else
          call res_tree_qqb_gen(C2Lim%AmpMom,res_lo)
          res_lo = multiply_IS_charges_sq(res_lo,2)
          call get_respdf_gen(0,1,C2Lim,res_lo,respdf)
       endif
       
       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z)/(one-z))&
            * C2Lim%wgt

       FintNLO_ns(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- S and CS
       
    !-- S
    call cut_histo(SLim)
    
    if (SLim%makecut.or.SLim%flag) then

       kin(4) = zero
       FintNLO_ns(4:6) = zero

    else
       e5 = SLim%Lim_KinInv(1)
       if (oldcode) then
          call res_tree_qqb(SLim%AmpMom,res_lo_old)
          call get_qed_eik(charges_ns,SLim%Lim_etaij,[1,2,3,4],5,eik)
          res_lo_tmp_old(:,1) = res_lo_old(:,1) * eik(1:2) !-- dn
          res_lo_tmp_old(:,2) = res_lo_old(:,2) * eik(3:4) !-- up
          call get_respdf(ns_lumi,0,1,SLim,res_lo_tmp_old,respdf_tmp(:,1))
       else
          call res_tree_qqb_gen(SLim%AmpMom,res_lo_tmp)
          call get_qed_eik_gen(res_lo_tmp,SLim%Lim_etaij,[1,2,3,4],5,res_lo)
          call get_respdf_gen(0,1,SLim,res_lo,respdf_tmp(:,1))
       endif
       
       call partition_nlo_qed(SLim,damp,12)
       
       respdf_tmp(:,1) = -respdf_tmp(:,1)/e5**2 * SLim%wgt * damp

       FintNLO_ns(4) = respdf_tmp(1,1)

       !-- SC1
       if (oldcode) then
          res_lo_tmp_old(:,1) = res_lo_old(:,1) * Qdn2
          res_lo_tmp_old(:,2) = res_lo_old(:,2) * Qup2
          call get_respdf(ns_lumi,0,1,SLim,res_lo_tmp_old,respdf_tmp(:,2))
          respdf_tmp(:,3) = respdf_tmp(:,2)
       else
          res_lo = multiply_IS_charges_sq(res_lo_tmp,1)
          call get_respdf_gen(0,1,SLim,res_lo,respdf_tmp(:,2))
          res_lo = multiply_IS_charges_sq(res_lo_tmp,2)
          call get_respdf_gen(0,1,SLim,res_lo,respdf_tmp(:,3))
       endif
       
       
       e5    = SC1Lim%Lim_KinInv(1)
       eta5i = SC1Lim%Lim_KinInv(2)

       respdf_tmp(:,2) = respdf_tmp(:,2)/e5**2/eta5i * SC1Lim%wgt
       FintNLO_ns(5) = respdf_tmp(1,2)

       !-- SC2
       e5    = SC2Lim%Lim_KinInv(1)
       eta5i = SC2Lim%Lim_KinInv(2)

       respdf_tmp(:,3) = respdf_tmp(:,3)/e5**2/eta5i * SC2Lim%wgt
       
       FintNLO_ns(6) = respdf_tmp(1,3)

       respdf(:) = sum(respdf_tmp(:,1:3),2)
       kin(4) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_ns)
    
#if(_withchecks == 1)
    limval_nlo_is = FintNLO_ns
    limval_nlo_is(2) = FintNLO_ns(4)  ! soft
    limval_nlo_is(3:4) = FintNLO_ns(2:3)  ! coll
#endif

  end function xsect_nloewk_r_is_ns

  
  !-------------------------------------------------------------
  !-- final-state sectors below
  !-------------------------------------------------------------
  
  function xsect_nloewk_r_fs_53_ns(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_r_fs_53_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nloewk_r_fs_53_ns = xsect_nlo_5i_a_ns(yRnd,ff,vegasweight,3,4)

  end function xsect_nloewk_r_fs_53_ns

  function xsect_nloewk_r_fs_54_ns(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_r_fs_54_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nloewk_r_fs_54_ns = xsect_nlo_5i_a_ns(yRnd,ff,vegasweight,4,3)

  end function xsect_nloewk_r_fs_54_ns

  !-- same for aa channel
  
  function xsect_nloewk_r_fs_53_aa(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_r_fs_53_aa
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nloewk_r_fs_53_aa = xsect_nlo_5i_a_aa(yRnd,ff,vegasweight,3,4)

  end function xsect_nloewk_r_fs_53_aa

  function xsect_nloewk_r_fs_54_aa(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_r_fs_54_aa
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nloewk_r_fs_54_aa = xsect_nlo_5i_a_aa(yRnd,ff,vegasweight,4,3)

  end function xsect_nloewk_r_fs_54_aa

  !----------------------------------------------------------------------------------------------------
  
  !-- master sector for final-state sctors below
  function xsect_nlo_5i_a_ns(yRnd,ff,vegasweight,icoll,jother)
    integer :: xsect_nlo_5i_a_ns,icoll,jother
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,CLim,CSLim,SLim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(4),kin(3)
    real(dp)    :: respdf(ipdf),respdf_tmp(ipdf,2)
    real(dp)    :: res_nlo(-5:7,-5:7),res_lo(-5:7,-5:7),res_lo_tmp(-5:7,-5:7)
    real(dp)    :: res_tmp_vect(2,2,2),respdf_vect(2,ipdf)    
    real(dp)    :: z,s5i,eik(4),e5,eta5i
    real(dp)    :: damp,Qsq_FS(2)
    !--
    real(dp)    :: res_nlo_old(2,2),res_lo_old(2,2),res_lo_tmp_old(2,2)
    logical     :: oldcode

    oldcode = .false.

    xsect_nlo_5i_a_ns = 0

    limval_nlo = zero

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    if (xx(xE)*xx(xRHO).lt.buff_r) then
       failed_points = failed_points + 1
       return
    endif
   
    call open_histo()

    call kinematics_nlo_fs(xx,icoll,jother,HardProc,CLim,CSLim,SLim)
    
#if (_Vcharge == 0)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    CLim%ids(1:4)     = [0,0,id_el,-id_el]
    SLim%ids(1:4)     = [0,0,id_el,-id_el]

    HardProc%part = [id_q,-id_q,id_el,-id_el,id_a]
    CLim%part     = [id_q,-id_q,id_el,-id_el]
    SLim%part     = [id_q,-id_q,id_el,-id_el]

#elif  (_Vcharge == -1)
    HardProc%ids(1:5) = [0,0,id_el,-id_nue,id_a]
    CLim%ids(1:4)     = [0,0,id_el,-id_nue]
    SLim%ids(1:4)     = [0,0,id_el,-id_nue]

    HardProc%part = [id_q,-id_qp,id_el,-id_nue,id_a]
    CLim%part     = [id_q,-id_qp,id_el,-id_nue]
    SLim%part     = [id_q,-id_qp,id_el,-id_nue]

#elif  (_Vcharge == +1)
    HardProc%ids(1:5) = [0,0,id_nue,-id_el,id_a]
    CLim%ids(1:4)     = [0,0,id_nue,-id_el]
    SLim%ids(1:4)     = [0,0,id_nue,-id_el]

    HardProc%part = [id_q,-id_qp,id_nue,-id_el,id_a]
    CLim%part     = [id_q,-id_qp,id_nue,-id_el]
    SLim%part     = [id_q,-id_qp,id_nue,-id_el]
#endif

    Qsq_FS = [Q3**2, Q4**2]
    
    !-- Hard
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_ns(1) = zero

    else    

       if (oldcode) then
          call res_tree_a_qqb(HardProc%AmpMom,res_nlo_old)
          call get_respdf(ns_lumi,0,1,HardProc,res_nlo_old,respdf)
       else
          call res_tree_a_qqb_gen(HardProc%AmpMom,res_nlo)
          call get_respdf_gen(0,1,HardProc,res_nlo,respdf)
       endif

       call partition_nlo_qed(HardProc,damp,icoll)
       
       respdf = respdf*HardProc%wgt*damp

       kin(1) = respdf(1)
       FintNLO_ns(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C
    call cut_histo(CLim)
    
    if (CLim%makecut.or.CLim%flag) then

       kin(2) = zero
       FintNLO_ns(2) = zero

    else

       if (oldcode) then

          call res_tree_qqb(CLim%AmpMom,res_lo_old)
          res_lo_old = Q_lep2*res_lo_old
          call get_respdf(ns_lumi,0,1,CLim,res_lo_old,respdf)
       else
          
          call res_tree_qqb_gen(CLim%AmpMom,res_lo)
          res_lo = Qsq_Fs(icoll-2)**2 * res_lo
          call get_respdf_gen(0,1,CLim,res_lo,respdf)

       endif

       z   = CLim%Lim_KinInv(1)
       s5i = CLim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z))&
            * CLim%wgt

       FintNLO_ns(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- S and CS
    call cut_histo(SLim)
    
    if (SLim%makecut.or.SLim%flag) then

       kin(3) = zero
       FintNLO_ns(3:4) = zero

    else
       if (oldcode) then
          call res_tree_qqb(SLim%AmpMom,res_lo_old)
          call get_qed_eik(charges_ns,SLim%Lim_etaij,[1,2,3,4],5,eik)
          res_tmp_vect(:,1,1) = res_lo_old(:,1) * eik(1:2) !-- dn
          res_tmp_vect(:,2,1) = res_lo_old(:,2) * eik(3:4) !-- up
          res_tmp_vect(:,:,2) = res_lo_old(:,:) * Q_lep2
          call get_respdf_vect(ns_lumi,0,1,SLim,res_tmp_vect(:,:,1:2), &
                  respdf_vect(1:2,:))

       else
          call res_tree_qqb_gen(SLim%AmpMom,res_lo_tmp)
          call get_qed_eik_gen(res_lo_tmp,SLim%Lim_etaij,[1,2,3,4],5,res_lo)
          call get_respdf_gen(0,1,SLim,res_lo,respdf_tmp(:,1))

       endif 
       
       !-- S
       call partition_nlo_qed(SLim,damp,icoll)
       
       e5 = SLim%Lim_KinInv(1)

       if (oldcode) then
           respdf_tmp(:,1) = -respdf_vect(1,:)/e5**2 * SLim%wgt * damp
       else        
           respdf_tmp(:,1) = -respdf_tmp(:,1)/e5**2 * SLim%wgt * damp
       endif

       FintNLO_ns(3) = respdf_tmp(1,1)

       !-- CS
       e5    = CSLim%Lim_KinInv(1)
       eta5i = CSLim%Lim_KinInv(2)

       if (oldcode) then
               respdf_tmp(:,2) = respdf_vect(2,:)/e5**2/eta5i * CSLim%wgt
       else
               res_lo = res_lo_tmp * Qsq_Fs(icoll-2)**2 /e5**2/eta5i * CSLim%wgt
               call get_respdf_gen(0,1,SLim,res_lo,respdf_tmp(:,2))
       endif
       
       FintNLO_ns(4) = respdf_tmp(1,2)

       respdf(:) = sum(respdf_tmp(:,1:2),2)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_ns)


    
#if(_withchecks == 1)
    FintNLO_ew(1:4) = FintNLO_ns
    limval_nlo = FintNLO_ns
    limval_nlo(2) = FintNLO_ns(3)   ! soft
    limval_nlo(3) = FintNLO_ns(2)   ! soft
#endif

  end function xsect_nlo_5i_a_ns

  !-- master sector for final-state sctors below, for a a -> e-e+ a
  function xsect_nlo_5i_a_aa(yRnd,ff,vegasweight,icoll,jother)
    integer :: xsect_nlo_5i_a_aa,icoll,jother
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,CLim,CSLim,SLim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_aa(4),kin(3)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(1,1),res_lo(1,1),respdf_tmp(ipdf,2)
    real(dp)    :: z,s5i,eik,e5,eta5i
    real(dp)    :: damp

    xsect_nlo_5i_a_aa = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    if (xx(xE)*xx(xRHO).lt.buff_r) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_fs(xx,icoll,jother,HardProc,CLim,CSLim,SLim)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    CLim%ids(1:4)     = [0,0,id_el,-id_el]
    SLim%ids(1:4)     = [0,0,id_el,-id_el]

    HardProc%part = [id_a,id_a,id_el,-id_el,id_a]
    CLim%part     = [id_a,id_a,id_el,-id_el]
    SLim%part     = [id_a,id_a,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_aa(1) = zero

    else    

       call res_treeAA_a_aa(HardProc%AmpMom,res_nlo(1,1))
       call get_respdf(aa_lumi,0,1,HardProc,res_nlo,respdf)

       call partition_nloAA_qed(HardProc,damp,icoll)
       
       respdf = respdf*HardProc%wgt*damp

       kin(1) = respdf(1)
       FintNLO_aa(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C
    call cut_histo(CLim)
    
    if (CLim%makecut.or.CLim%flag) then

       kin(2) = zero
       FintNLO_aa(2) = zero

    else

       call res_treeAA_aa(CLim%AmpMom,res_lo(1,1))
       call get_respdf(aa_lumi,0,1,CLim,res_lo,respdf)

       z   = CLim%Lim_KinInv(1)
       s5i = CLim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision
       respdf = (-one)*respdf*(two/s5i)*(Q_lep2*Pqg(z))&
            * CLim%wgt

       FintNLO_aa(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- S and CS
    call cut_histo(SLim)
    
    if (SLim%makecut.or.SLim%flag) then

       kin(3) = zero
       FintNLO_aa(3:4) = zero

    else

       call res_treeAA_aa(SLim%AmpMom,res_lo(1,1))
       call get_respdf(aa_lumi,0,1,SLim,res_lo,respdf)

       !-- S, use QCD eikonal not to return a matrix
       call get_qcd_eik(Q_lep**2,SLim%Lim_etaij,[3,4],5,eik)

       call partition_nloAA_qed(SLim,damp,icoll)
       
       e5 = SLim%Lim_KinInv(1)
       respdf_tmp(:,1) = -respdf(:)*eik/e5**2 * SLim%wgt * damp
       FintNLO_aa(3) = respdf_tmp(1,1)

       !-- CS
       e5    = CSLim%Lim_KinInv(1)
       eta5i = CSLim%Lim_KinInv(2)
       respdf_tmp(:,2) = respdf*Q_lep2/e5**2/eta5i * CSLim%wgt
       FintNLO_aa(4) = respdf_tmp(1,2)

       respdf(:) = sum(respdf_tmp(:,1:2),2)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_aa)
    
#if(_withchecks == 1)
    FintNLO_ew(1:4) = FintNLO_aa
#endif

  end function xsect_nlo_5i_a_aa

  !-------------------------------------------------------------
  !-- aq channel below
  !-------------------------------------------------------------
  
  function xsect_nloewk_r_is_aq(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_r_is_aq
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(3),kin(3)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo_old(2,2),res_lo_old(2,2),res_tmp_old(2,2),res_loAA_old
    real(dp)    :: res_nlo(-5:7,-5:7),res_lo(-5:7,-5:7),res_tmp(-5:7,-5:7),res_loAA(-5:7,-5:7)
    real(dp)    :: z,s5i
    logical     :: oldcode

    xsect_nloewk_r_is_aq = 0

    limval_nlo_is = zero

    oldcode = .false.

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_r)) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,&
         C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.false.)

#if (_Vcharge == 0)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
#elif  (_Vcharge == -1)
    HardProc%ids(1:5) = [0,0,id_el,-id_nue,-id_q]
    C1Lim%ids(1:4)     = [0,0,id_el,-id_nue]
#elif  (_Vcharge == +1)
    HardProc%ids(1:5) = [0,0,id_nue,-id_el,-id_q]
    C1Lim%ids(1:4)     = [0,0,id_nue,-id_el]
#endif
    
    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_ns(1) = zero

    else

       if (oldcode) then
          call res_tree_a_aq(HardProc%AmpMom,res_nlo_old)
          call get_respdf(aq_lumi,0,1,HardProc,res_nlo_old,respdf)
       else
          call res_tree_a_aq_gen(HardProc%AmpMom,res_nlo)
          call get_respdf_gen(0,1,HardProc,res_nlo,respdf)
       endif
       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintNLO_ns(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintNLO_ns(2) = zero

    else
       if (oldcode) then
          call res_tree_qqb(C1Lim%AmpMom,res_lo_old)
          res_lo_old(1,:) = [Qdn2,Qup2] * res_lo_old(1,:)
          res_lo_old(2,:) = [Qdn2,Qup2] * res_lo_old(2,:)
          
          call get_respdf(aq_lumi,0,1,C1Lim,res_lo_old,respdf)
       else
          call res_tree_qqb_gen(C1Lim%AmpMom,res_lo)
          res_lo = multiply_IS_charges_sq(res_lo,1)
          res_lo = transition('ga -> q', 'none', res_lo)
          call get_respdf_gen(0,1,C1Lim,res_lo,respdf)
          print*, 'respdf= ', respdf
       endif
          
       z   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       !-- xn = aveqa/aveqq
       respdf = (-one)*respdf*(two/s5i)*(xn*Pgq_spav(z)/(one-z))&
            * C1Lim%wgt

       FintNLO_ns(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C2 -> only present in the neutral-change case
#if (_Vcharge == 0)
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(3) = zero
       FintNLO_ns(3) = zero

    else

       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- here we use spin-averages since there are no spin correlations for the aa -> e-e+ amplitude
       if (oldcode) then
          call res_treeAA_aa(C2Lim%AmpMom,res_loAA_old)
          !-- 1=xn*aveqa/aveaa
          res_loAA_old = res_loAA_old * Pqq(z)/(one-z) !-- use Pqq/(1-z) because of definition of z

          res_tmp_old(:,1) = Qdn2 * res_loAA_old
          res_tmp_old(:,2) = Qup2 * res_loAA_old
       
          call get_respdf(aq_lumi,0,1,C2Lim,res_tmp_old,respdf)
       else
          call res_treeAA_aa_gen(C2Lim%AmpMom,res_lo)
          res_lo = transition('none','q -> ga',  res_lo)
          res_lo = multiply_IS_charges_sq(res_lo,2)

          call get_respdf_gen(0,1,C2Lim,res_lo,respdf)
          respdf = respdf * Pqq(z) / (one-z)
       endif
       
       respdf = (-one)*respdf*(two/s5i)&
            * C2Lim%wgt

       FintNLO_ns(3) = respdf(1)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()
#endif
    
    call check_ff(ff,xx,FintNLO_ns)
    
#if(_withchecks == 1)
    FintNLO_ew(1:3) = FintNLO_ns
    limval_nlo_is(1) = FintNLO_ns(1)
    limval_nlo_is(3:4) = FintNLO_ns(2:3)   ! coll
#endif
    
  end function xsect_nloewk_r_is_aq

  !-------------------------------------------------------------
  !-- qa channel below
  !-------------------------------------------------------------
  
  function xsect_nloewk_r_is_qa(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_r_is_qa
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(3),kin(3)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo_old(2,2),res_lo_old(2,2),res_tmp(2,2),res_loAA,res_nlo(-5:7,-5:7), res_lo(-5:7,-5:7)
    real(dp)    :: z,s5i
    logical     :: oldcode


    oldcode = .false.

    xsect_nloewk_r_is_qa = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_r)) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,&
         C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.false.)

#if (_Vcharge == 0)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
#elif  (_Vcharge == -1)
    HardProc%ids(1:5) = [0,0,id_el,-id_nue,id_qp]
    C2Lim%ids(1:4)     = [0,0,id_el,-id_nue]
#elif  (_Vcharge == +1)
    HardProc%ids(1:5) = [0,0,id_nue,-id_el,id_qp]
    C2Lim%ids(1:4)     = [0,0,id_nue,-id_el]
#endif
    
    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_ns(1) = zero

    else    

       if (oldcode) then
          call res_tree_a_qa(HardProc%AmpMom,res_nlo_old)
          call get_respdf(qa_lumi,0,1,HardProc,res_nlo_old,respdf)
       else
          call res_tree_a_qa_gen(HardProc%AmpMom,res_nlo)
          call get_respdf_gen(0,1,HardProc,res_nlo,respdf)
       endif
       
       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintNLO_ns(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1 -> only present in the neutral-change case
#if (_Vcharge == 0)
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintNLO_ns(2) = zero

    else

       z   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       if (oldcode) then
          !-- here we use spin-averages since there are no spin correlations for the aa -> e-e+ amplitude
          call res_treeAA_aa(C1Lim%AmpMom,res_loAA)
          !-- 1=xn*aveqa/aveaa
          res_loAA = res_loAA * Pqq(z)/(one-z) !-- use Pqq/(1-z) because of definition of z
          
          res_tmp(1,:) = [Qdn2,Qup2] * res_loAA
          res_tmp(2,:) = [Qdn2,Qup2] * res_loAA
       
          call get_respdf(qa_lumi,0,1,C1Lim,res_tmp,respdf)
       else
          call res_treeAA_aa_gen(C1Lim%AmpMom,res_lo)
          res_lo = transition('q -> ga','none', res_lo)
          res_lo = multiply_IS_charges_sq(res_lo,1)

          call get_respdf_gen(0,1,C1Lim,res_lo,respdf)
       endif
       
       respdf = respdf * Pqq(z) / (one-z)
          
       respdf = (-one)*respdf*(two/s5i)&
            * C1Lim%wgt

       FintNLO_ns(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif
#endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(3) = zero
       FintNLO_ns(3) = zero

    else

       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       if (oldcode) then
          call res_tree_qqb(C2Lim%AmpMom,res_lo_old)
          res_lo_old(1,:) = [Qdn2,Qup2] * res_lo_old(1,:)
          res_lo_old(2,:) = [Qdn2,Qup2] * res_lo_old(2,:)

          call get_respdf(qa_lumi,0,1,C2Lim,res_lo_old,respdf)
       else
          call res_tree_qqb_gen(C2Lim%AmpMom,res_lo)
          res_lo = multiply_IS_charges_sq(res_lo,1)
          res_lo = transition('none','ga -> q', res_lo)
          call get_respdf_gen(0,1,C2Lim,res_lo,respdf)
       endif


       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       !-- xn = aveqa/aveqq
       respdf = (-one)*respdf*(two/s5i)*(xn*Pgq_spav(z)/(one-z))&
            * C2Lim%wgt

       FintNLO_ns(3) = respdf(1)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_ns)
    
#if(_withchecks == 1)
    FintNLO_ew(1:3) = FintNLO_ns
#endif
    
  end function xsect_nloewk_r_is_qa


  function xsect_nloewk_r_is_ns_wp(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_r_is_ns_wp
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim,SLim,SC1Lim,SC2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(6),kin(4)
    real(dp)    :: respdf(ipdf),respdf_tmp(ipdf,3)
    real(dp)    :: res_nlo(1,2),res_lo(1,2),res_tmp(1,2)
    real(dp)    :: z,s5i,eik(4),e5,eta5i
    real(dp)    :: damp


    print*, 'I AM HERE'
    
    xsect_nloewk_r_is_ns_wp = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_r)) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,&
        C1Lim=C1Lim,C2Lim=C2Lim,SLim=SLim,SC1Lim=SC1Lim,SC2Lim=SC2Lim,compute_etas=.true.)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    C1Lim%ids(1:4)    = [0,0,id_el,-id_el]
    C2Lim%ids(1:4)    = [0,0,id_el,-id_el]
    SLim%ids(1:4)     = [0,0,id_el,-id_el]
    SC1Lim%ids(1:4)   = [0,0,id_el,-id_el]
    SC2Lim%ids(1:4)   = [0,0,id_el,-id_el]

    !!!!!! TO DO:
    !check ids for final state neutrinos -> relevant for analysis 



    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_ns(1) = zero

    else    
       call res_tree_a_qqb_wp(HardProc%AmpMom,res_nlo)

       print*, 'res_nlo= ', res_nlo
       res_nlo(1,2) = zero



      print*, 'res_nlo= ', res_nlo 


       call get_respdf(qQpb_lumi_wp,1,0,HardProc,res_nlo,respdf)
       call partition_nlo_qed(HardProc,damp,12)
       
       respdf = respdf*HardProc%wgt*damp




       print*, 'HardProc%AmpMom(:,1)', HardProc%AmpMom(:,1)
       print*, 'HardProc%AmpMom(:,2)', HardProc%AmpMom(:,2)
       print*, 'HardProc%AmpMom(:,3)', HardProc%AmpMom(:,3)
       print*, 'HardProc%AmpMom(:,4)', HardProc%AmpMom(:,4)
       print*, 'HardProc%AmpMom(:,5)', HardProc%AmpMom(:,5)

       print*, 'HardProc%wgt= ', HardProc%wgt
       print*, 'damp= ', damp
       print*, 'respdf= ', respdf



       kin(1) = respdf(1)
       FintNLO_ns(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintNLO_ns(2) = zero

    else

       call res_tree_qqb_wp(C1Lim%AmpMom,res_lo)
       res_lo(1,1) = Qup2 * res_lo(1,1)
       res_lo(1,2) = Qdn2 * res_lo(1,2)
      


       res_lo(1,2) = zero

       print*, 'res_lo = ', res_lo 
       
       print*, 'C1Lim%AmpMom(:,1)', C1Lim%AmpMom(:,1)
       print*, 'C1Lim%AmpMom(:,2)', C1Lim%AmpMom(:,2)
       print*, 'C1Lim%AmpMom(:,3)', C1Lim%AmpMom(:,3)
       print*, 'C1Lim%AmpMom(:,4)', C1Lim%AmpMom(:,4)





       call get_respdf(qQpb_lumi_wp,1,0,C1Lim,res_lo,respdf)

       z   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z)/(one-z))&
            * C1Lim%wgt

       print*, 'C1Lim%wgt= ', C1Lim%wgt
       print*, 'respdf*(two/s5i)*(Pqg(z)/(one-z)) = ', respdf*(two/s5i)*(Pqg(z)/(one-z))
       print*, 'z = ', z
       print*, 's5i = ', s5i


       FintNLO_ns(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(3) = zero
       FintNLO_ns(3) = zero

    else

       call res_tree_qqb_wp(C2Lim%AmpMom,res_lo)
       res_lo(1,1) = Qdn2 * res_lo(1,1)
       res_lo(1,2) = Qup2 * res_lo(1,2)
       
       call get_respdf(qQpb_lumi_wp,1,0,C2Lim,res_lo,respdf)

       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z)/(one-z))&
            * C2Lim%wgt

       FintNLO_ns(3) = respdf(1)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- S and CS
    !call cut_histo(SLim)
    !
    !if (SLim%makecut.or.SLim%flag) then
    !
    !   kin(4) = zero
    !   FintNLO_ns(4:6) = zero
    !
    !else
    !
    !   call res_tree_qqb_w(SLim%AmpMom,res_lo)
    !
    !   !-- S
    !   call get_qed_eik(charges_ns,SLim%Lim_etaij,[1,2,3,4],5,eik)
    !   res_tmp(:,1) = res_lo(:,1) * eik(1:2) !-- dn
    !   res_tmp(:,2) = res_lo(:,2) * eik(3:4) !-- up
    !   call get_respdf(qQpb_lumi_wp,1,0,SLim,res_tmp,respdf_tmp(:,1))
    !
    !   call partition_nlo_qed(SLim,damp,12)
    !   
    !   e5 = SLim%Lim_KinInv(1)
    !   respdf_tmp(:,1) = -respdf_tmp(:,1)/e5**2 * SLim%wgt * damp
    !   FintNLO_ns(4) = respdf_tmp(1,1)
    !
    !   !-- SC1
    !   res_tmp(:,1) = res_lo(:,1) * Qdn2
    !   res_tmp(:,2) = res_lo(:,2) * Qup2
    !   call get_respdf(qQpb_lumi_wp,1,0,SLim,res_tmp,respdf_tmp(:,2))
    !   respdf_tmp(:,3) = respdf_tmp(:,2)
    !
    !   e5    = SC1Lim%Lim_KinInv(1)
    !   eta5i = SC1Lim%Lim_KinInv(2)
    !   respdf_tmp(:,2) = respdf_tmp(:,2)/e5**2/eta5i * SC1Lim%wgt
    !   FintNLO_ns(5) = respdf_tmp(1,2)
    !
    !   !-- SC2
    !   e5    = SC2Lim%Lim_KinInv(1)
    !   eta5i = SC2Lim%Lim_KinInv(2)
    !   respdf_tmp(:,3) = respdf_tmp(:,3)/e5**2/eta5i * SC2Lim%wgt
    !   FintNLO_ns(6) = respdf_tmp(1,3)
    !
    !   respdf(:) = sum(respdf_tmp(:,1:3),2)
    !   kin(4) = respdf(1)
    !   
    !   call fill_histo(respdf,vegasweight)
    !
    !endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_ns)
    
#if(_withchecks == 1)
    FintNLO_ew(1:6) = FintNLO_ns
#endif




    !!!!!!!!!!!!!debug CS
    if (FintNLO_ns(1) .ne. 0 ) then
    print*, 'FintNLO_ns(1) = ', FintNLO_ns(1)
    print*, 'FintNLO_ns(2) = ', FintNLO_ns(2), '   ', FintNLO_ns(2)/FintNLO_ns(1)


    stop

    
    print*, ''

    endif




  end function xsect_nloewk_r_is_ns_wp


  
  function xsect_nloewk_r_is_ns_wm(yRnd,ff,vegasweight)
    integer :: xsect_nloewk_r_is_ns_wm
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim,SLim,SC1Lim,SC2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(6),kin(4)
    real(dp)    :: respdf(ipdf),respdf_tmp(ipdf,3)
    real(dp)    :: res_nlo(2,2),res_lo(2,2),res_tmp(2,2)
    real(dp)    :: z,s5i,eik(4),e5,eta5i
    real(dp)    :: damp

    xsect_nloewk_r_is_ns_wm = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_r)) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,&
        C1Lim=C1Lim,C2Lim=C2Lim,SLim=SLim,SC1Lim=SC1Lim,SC2Lim=SC2Lim,compute_etas=.true.)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    C1Lim%ids(1:4)    = [0,0,id_el,-id_el]
    C2Lim%ids(1:4)    = [0,0,id_el,-id_el]
    SLim%ids(1:4)     = [0,0,id_el,-id_el]
    SC1Lim%ids(1:4)   = [0,0,id_el,-id_el]
    SC2Lim%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_ns(1) = zero

    else    

       !TEST MADGRAPH
       !MG ps
       !HardProc%AmpMom(:,1) = (/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03/)
       !HardProc%AmpMom(:,2) = (/0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  -0.5000000E+03/)
       !HardProc%AmpMom(:,3) = (/0.4585788E+03,  0.1694532E+03,  0.3796537E+03,  -0.1935025E+03/)
       !HardProc%AmpMom(:,4) = (/0.3640666E+03, -0.1832987E+02, -0.3477043E+03,  0.1063496E+03/)
       !HardProc%AmpMom(:,5) = (/0.1773546E+03, -0.1511234E+03, -0.3194936E+02,  0.8715287E+02/)
       
       call res_tree_a_qqb_wm(HardProc%AmpMom,res_nlo)

       !print*, 'NLO amp udx_epvea  ', (0.094835522759998875_dp)**2*res_nlo(1,1)/eesq2*(four*pi)/132.50700000000001_dp
       !print*, 'NLO amp dxu_epvea  ', (0.094835522759998875_dp)**2*res_nlo(1,2)/eesq2*(four*pi)/132.50700000000001_dp
       !print*, 'NLO amp dux_epvea  ', (0.094835522759998875_dp)**2*res_nlo(2,1)/eesq2*(four*pi)/132.50700000000001_dp
       !print*, 'NLO amp uxd_epvea  ', (0.094835522759998875_dp)**2*res_nlo(2,2)/eesq2*(four*pi)/132.50700000000001_dp
       !stop
       
       call get_respdf(qQpb_lumi_wm,1,0,HardProc,res_nlo,respdf)

       call partition_nlo_qed(HardProc,damp,12)
       
       respdf = respdf*HardProc%wgt*damp

       print*, 'HardProc%wgt = ', HardProc%wgt
       stop

       kin(1) = respdf(1)
       FintNLO_ns(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintNLO_ns(2) = zero

    else

       call res_tree_qqb_wm(C1Lim%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)
       
       call get_respdf(qQpb_lumi_wm,1,0,C1Lim,res_lo,respdf)

       z   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z)/(one-z))&
            * C1Lim%wgt



       print*, 'C1Lim%wgt = ', C1Lim%wgt

       FintNLO_ns(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(3) = zero
       FintNLO_ns(3) = zero

    else

       call res_tree_qqb_wm(C2Lim%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)
       
       call get_respdf(qQpb_lumi_wm,1,0,C2Lim,res_lo,respdf)

       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       respdf = (-one)*respdf*(two/s5i)*(Pqg(z)/(one-z))&
            * C2Lim%wgt

       FintNLO_ns(3) = respdf(1)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- S and CS
    call cut_histo(SLim)
    
    if (SLim%makecut.or.SLim%flag) then

       kin(4) = zero
       FintNLO_ns(4:6) = zero

    else

       call res_tree_qqb_wm(SLim%AmpMom,res_lo)

       !-- S
       call get_qed_eik(charges_ns,SLim%Lim_etaij,[1,2,3,4],5,eik)
       res_tmp(:,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_tmp(:,2) = res_lo(:,2) * eik(3:4) !-- up
       call get_respdf(qQpb_lumi_wm,1,0,SLim,res_tmp,respdf_tmp(:,1))

       call partition_nlo_qed(SLim,damp,12)
       
       e5 = SLim%Lim_KinInv(1)
       respdf_tmp(:,1) = -respdf_tmp(:,1)/e5**2 * SLim%wgt * damp
       FintNLO_ns(4) = respdf_tmp(1,1)

       !-- SC1
       res_tmp(:,1) = res_lo(:,1) * Qdn2
       res_tmp(:,2) = res_lo(:,2) * Qup2
       call get_respdf(qQpb_lumi_wm,1,0,SLim,res_tmp,respdf_tmp(:,2))
       respdf_tmp(:,3) = respdf_tmp(:,2)

       e5    = SC1Lim%Lim_KinInv(1)
       eta5i = SC1Lim%Lim_KinInv(2)
       respdf_tmp(:,2) = respdf_tmp(:,2)/e5**2/eta5i * SC1Lim%wgt
       FintNLO_ns(5) = respdf_tmp(1,2)

       !-- SC2
       e5    = SC2Lim%Lim_KinInv(1)
       eta5i = SC2Lim%Lim_KinInv(2)
       respdf_tmp(:,3) = respdf_tmp(:,3)/e5**2/eta5i * SC2Lim%wgt
       FintNLO_ns(6) = respdf_tmp(1,3)

       respdf(:) = sum(respdf_tmp(:,1:3),2)
       kin(4) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_ns)
    
#if(_withchecks == 1)
    FintNLO_ew(1:6) = FintNLO_ns
#endif

  end function xsect_nloewk_r_is_ns_wm

  
end module mod_xsects_nloewk_r
  

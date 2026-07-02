!-- extra photon emission, qcd loop
module mod_xsects_nnlo_rvqcd
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_kinematics_nlo
  use mod_process
  use mod_auxfunctions
  use mod_aux_sectors
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_amplitudes_loop_ppll
  use mod_limits_rv
  use mod_splittings_bare
  use mod_eikonals
  use mod_partitions
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintRVQCD_ns(6)
  real(dp), public, save :: FintRVQCD_aq(3)
  real(dp), public, save :: FintRVQCD_qa(3)
#endif
  real(dp), parameter :: charges_ns(4,4) = reshape(& 
       [-Qdn,   Qdn,  -Qup,   Qup,   &
         Qdn,  -Qdn,   Qup,  -Qup,   &
         Q_lep, Q_lep, Q_lep, Q_lep, &
        -Q_lep,-Q_lep,-Q_lep,-Q_lep]  &
        , [4,4])
  
  private

  public :: xsect_nnlo_rvqcd_is_ns
  public :: xsect_nnlo_rvqcd_fs_53_ns,xsect_nnlo_rvqcd_fs_54_ns

  public :: xsect_nnlo_rvqcd_is_aq,xsect_nnlo_rvqcd_is_qa
  
contains

  !-------------------------------------------------------------
  !-- initial-state sector, ns
  !-------------------------------------------------------------

  function xsect_nnlo_rvqcd_is_ns(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rvqcd_is_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim,SLim,SC1Lim,SC2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintRV_ns(6),kin(4)
    real(dp)    :: respdf(ipdf),respdf_tmp(ipdf,3)
    real(dp)    :: res0_old(2,2),res1_old(2,2)
    real(dp)    :: res0(-5:7,-5:7),res1(-5:7,-5:7)
    real(dp)    :: res0red_old(2,2),res1red_old(2,2),res_tmp_old(2,2)
    real(dp)    :: res0red(-5:7,-5:7),res1red(-5:7,-5:7),res_tmp(-5:7,-5:7)
    real(dp)    :: limcol(2),eik(4),e5,eta5i
    real(dp)    :: damp
    logical :: oldcode
    integer :: i 

    oldcode = .false.

    xsect_nnlo_rvqcd_is_ns = 0

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

    ! define process specific partons
#if (_Vcharge == 0)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
    SLim%ids(1:4)     = [0,0,id_el,-id_el]
    SC1Lim%ids(1:4)   = [0,0,id_el,-id_el]
    SC2Lim%ids(1:4)   = [0,0,id_el,-id_el]

    HardProc%part = [id_q,-id_q,id_el,-id_el,id_a]
    C1Lim%part = [id_q,-id_q,id_el,-id_el]
    C2Lim%part = [id_q,-id_q,id_el,-id_el]
    SLim%part = [id_q,-id_q,id_el,-id_el]
    SC1Lim%part = [id_q,-id_q,id_el,-id_el]
    SC2Lim%part = [id_q,-id_q,id_el,-id_el]

#elif  (_Vcharge == -1)
    HardProc%ids(1:5) = [0,0,id_el,-id_nue,id_a]
    C1Lim%ids(1:4) = [0,0,id_el,-id_nue]
    C2Lim%ids(1:4) = [0,0,id_el,-id_nue]
    SLim%ids(1:4) = [0,0,id_el,-id_nue]
    SC1Lim%ids(1:4) = [0,0,id_el,-id_nue]
    SC2Lim%ids(1:4) = [0,0,id_el,-id_nue]

    HardProc%part = [id_q,-id_qp,id_el,-id_nue,id_a]
    C1Lim%part = [id_q,-id_qp,id_el,-id_nue]
    C2Lim%part = [id_q,-id_qp,id_el,-id_nue]
    SLim%part = [id_q,-id_qp,id_el,-id_nue]
    SC1Lim%part = [id_q,-id_qp,id_el,-id_nue]
    SC2Lim%part = [id_q,-id_qp,id_el,-id_nue]

#elif  (_Vcharge == +1)
    HardProc%ids(1:5) = [0,0,id_nue,-id_el,id_a]
    C1Lim%ids(1:4) = [0,0,id_nue,-id_el]
    C2Lim%ids(1:4) = [0,0,id_nue,-id_el]
    SLim%ids(1:4) = [0,0,id_nue,-id_el]
    SC1Lim%ids(1:4) = [0,0,id_nue,-id_el]
    SC2Lim%ids(1:4) = [0,0,id_nue,-id_el]

    HardProc%part = [id_q,-id_qp,id_nue,-id_el,id_a]
    C1Lim%part = [id_q,-id_qp,id_nue,-id_el]
    C2Lim%part = [id_q,-id_qp,id_nue,-id_el]
    SLim%part = [id_q,-id_qp,id_nue,-id_el]
    SC1Lim%part = [id_q,-id_qp,id_nue,-id_el]
    SC2Lim%part = [id_q,-id_qp,id_nue,-id_el]
#endif


    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintRV_ns(1) = zero

    else    

       if ( oldcode ) then 
          call res_qcdloop_a_qqb(HardProc%AmpMom,res0_old,res1_old)
          call get_respdf(ns_lumi,1,1,HardProc,res1_old,respdf)
       else 
          call ol_res_qcdloop_a_qqb_gen(HardProc%AmpMom,res0,res1)
          call get_respdf_gen(1,1,HardProc,res1,respdf)
       endif

       call partition_nlo_qed(HardProc,damp,12)
       
       respdf = respdf*HardProc%wgt*damp

       kin(1) = respdf(1)
       FintRV_ns(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintRV_ns(2) = zero

    else

       if ( oldcode ) then
          call res_qcdloop_qqb(C1Lim%AmpMom,res0red_old,res1red_old)
          limcol = col_rvqcd_is_ns(C1Lim)
          res_tmp_old(:,1) = ( res0red_old(:,1)*limcol(1) * Cf + res1red_old(:,1)*limcol(2) )*Qdn2
          res_tmp_old(:,2) = ( res0red_old(:,2)*limcol(1) * Cf + res1red_old(:,2)*limcol(2) )*Qup2
          call get_respdf(ns_lumi,1,1,C1Lim,res_tmp_old,respdf)
       else 
          call res_qcdloop_qqb_gen(C1Lim%AmpMom,res0red,res1red)
          limcol = col_rvqcd_is_ns(C1Lim)
          res_tmp = res0red(:,:)*limcol(1) * Cf + res1red(:,:)*limcol(2)
          res_tmp = multiply_IS_charges_sq(res_tmp,1)
          call get_respdf_gen(1,1,C1Lim,res_tmp,respdf)
       endif 

       respdf = -respdf*C1Lim%wgt

       FintRV_ns(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(3) = zero
       FintRV_ns(3) = zero

    else

       if ( oldcode ) then
          call res_qcdloop_qqb(C2Lim%AmpMom,res0red_old,res1red_old)
          limcol = col_rvqcd_is_ns(C2Lim)
          res_tmp_old(:,1) = ( res0red_old(:,1)*limcol(1) * Cf + res1red_old(:,1)*limcol(2) )*Qdn2
          res_tmp_old(:,2) = ( res0red_old(:,2)*limcol(1) * Cf + res1red_old(:,2)*limcol(2) )*Qup2
          call get_respdf(ns_lumi,1,1,C2Lim,res_tmp_old,respdf)
!          print*, 'res1red_old= ', res1red_old
!          print*, 'limcol= ', limcol
!          print*, 'res_tmp(1,1)= ', res_tmp(1,1)
!          print*, 'respdf= ', respdf

       else 
          call res_qcdloop_qqb_gen(C2Lim%AmpMom,res0red,res1red)
          limcol = col_rvqcd_is_ns(C2Lim)
          res_tmp = res0red(:,:)*limcol(1) * Cf + res1red(:,:)*limcol(2)
          res_tmp = multiply_IS_charges_sq(res_tmp,1)
          call get_respdf_gen(1,1,C2Lim,res_tmp,respdf)

!          print*, 'res1red= ', res1red
!          print*, 'limcol= ', limcol
!          print*, 'res_tmp(-1,1)= ', res_tmp(-1,1)
!          print*, 'respdf= ', respdf
       endif

       respdf = -respdf*C2Lim%wgt

       FintRV_ns(3) = respdf(1)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- S and CS
    call cut_histo(SLim)
    
    if (SLim%makecut.or.SLim%flag) then

       kin(4) = zero
       FintRV_ns(4:6) = zero

    else

       !--S
       if ( oldcode ) then
          call res_qcdloop_qqb(SLim%AmpMom,res0red_old,res1red_old)
          call get_qed_eik(charges_ns,SLim%Lim_etaij,[1,2,3,4],5,eik)
          res_tmp_old(:,1) = res1red_old(:,1) * eik(1:2) !-- dn
          res_tmp_old(:,2) = res1red_old(:,2) * eik(3:4) !-- up
          call get_respdf(ns_lumi,1,1,SLim,res_tmp_old,respdf_tmp(:,1))
       else 
          call res_qcdloop_qqb_gen(SLim%AmpMom,res0red,res_tmp)
          call get_qed_eik_gen(res_tmp,SLim%Lim_etaij,[1,2,3,4],5,res1red)
          call get_respdf_gen(1,1,SLim,res1red,respdf_tmp(:,1))
       endif

       call partition_nlo_qed(SLim,damp,12)
       
       e5 = SLim%Lim_KinInv(1)
       respdf_tmp(:,1) = -respdf_tmp(:,1)/e5**2 * SLim%wgt * damp
       FintRV_ns(4) = respdf_tmp(1,1)

       !-- SC1
       if ( oldcode ) then
          res_tmp_old(:,1) = res1red_old(:,1) * Qdn2
          res_tmp_old(:,2) = res1red_old(:,2) * Qup2
          call get_respdf(ns_lumi,1,1,SLim,res_tmp_old,respdf_tmp(:,2))
          respdf_tmp(:,3) = respdf_tmp(:,2)
       else 
          res_tmp = multiply_IS_charges_sq(res_tmp,1)
          call get_respdf_gen(1,1,SLim,res_tmp,respdf_tmp(:,2))
          respdf_tmp(:,3) = respdf_tmp(:,2)
       endif 

       e5    = SC1Lim%Lim_KinInv(1)
       eta5i = SC1Lim%Lim_KinInv(2)
       respdf_tmp(:,2) = respdf_tmp(:,2)/e5**2/eta5i * SC1Lim%wgt
       FintRV_ns(5) = respdf_tmp(1,2)

       !-- SC2
       e5    = SC2Lim%Lim_KinInv(1)
       eta5i = SC2Lim%Lim_KinInv(2)
       respdf_tmp(:,3) = respdf_tmp(:,3)/e5**2/eta5i * SC2Lim%wgt
       FintRV_ns(6) = respdf_tmp(1,3)

       respdf(:) = sum(respdf_tmp(:,1:3),2)
       kin(4) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif





!if (kin(1)*kin(2)*kin(3)*kin(4) .ne. zero ) then
!    print*, 'xx=', xx
!    do i = 1 , 4
!    print*, 'kin(', i, '=', kin(i)
!    enddo
!    pause
!endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintRV_ns)
  
#if(_withchecks == 1)
    FintRVQCD_ns(1:6) = FintRV_ns
#endif

  end function xsect_nnlo_rvqcd_is_ns

  !-------------------------------------------------------------
  !-- final-state sectors below
  !-------------------------------------------------------------

  function xsect_nnlo_rvqcd_fs_53_ns(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rvqcd_fs_53_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nnlo_rvqcd_fs_53_ns = xsect_nnlo_rvqcd_fs_5i_ns(yRnd,ff,vegasweight,3,4)

  end function xsect_nnlo_rvqcd_fs_53_ns

  function xsect_nnlo_rvqcd_fs_54_ns(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rvqcd_fs_54_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nnlo_rvqcd_fs_54_ns = xsect_nnlo_rvqcd_fs_5i_ns(yRnd,ff,vegasweight,4,3)

  end function xsect_nnlo_rvqcd_fs_54_ns

  !----------------------------------------------------------------------------------------------------
  
  !-- master sector for final-state sctors below
  function xsect_nnlo_rvqcd_fs_5i_ns(yRnd,ff,vegasweight,icoll,jother)
    integer :: xsect_nnlo_rvqcd_fs_5i_ns,icoll,jother
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,CLim,CSLim,SLim
    integer, parameter :: imax_ipdf = 2, imax_ilim = 4
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintRV_ns(4),kin(3)
    real(dp)    :: respdf(ipdf),respdf_tmp(ipdf,2)
    !--
    real(dp)    :: res0_old(2,2),res1_old(2,2)
    real(dp)    :: res0(-5:7,-5:7),res1(-5:7,-5:7), res_eikqed(-5:7,-5:7)
    real(dp)    :: res0red_old(2,2),res1red_old(2,2),res_tmp_old(2,2)
    real(dp)    :: res0red(-5:7,-5:7),res1red(-5:7,-5:7),res_tmp(-5:7,-5:7), res1red_fs_charge(-5:7,-5:7)
    real(dp) :: respdf_vect(imax_ipdf,ipdf),res_tmp_vect(2,2,imax_ipdf),respdf_vect_rev(ipdf,imax_ipdf)
    !--
    real(dp)    :: z,s5i,eik(4),e5,eta5i
    real(dp)    :: damp
    logical :: oldcode
    integer :: i

    oldcode = .false.

    xsect_nnlo_rvqcd_fs_5i_ns = 0

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

    ! for check
    ! z = one-1E-5_dp
    ! for check
    ! xx(xE) = 1E-5_dp
    ! for check
    ! xx(xRHO)=1E-6_dp

    !xx = [ &
    !0.27384198205924926_dp, &
    !0.29177341328902862_dp, &
    !2.6509538461396708E-002_dp, &
    !0.16909047503513841_dp, &
    !0.69571821718859106_dp, &
    !0.43621888932298258_dp, &
    !0.80514737395748393_dp, &
    !0.66869645956556378_dp  &
    !]

    call kinematics_nlo_fs(xx,icoll,jother,HardProc,CLim,CSLim,SLim)

        ! define process specific partons
#if (_Vcharge == 0)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    CLim%ids(1:4) = [0,0,id_el,-id_el]
    SLim%ids(1:4)     = [0,0,id_el,-id_el]

    HardProc%part = [id_q,-id_q,id_el,-id_el,id_a]
    CLim%part = [id_q,-id_q,id_el,-id_el]
    SLim%part = [id_q,-id_q,id_el,-id_el]

#elif  (_Vcharge == -1)
    HardProc%ids(1:5) = [0,0,id_el,-id_nue,id_a]
    CLim%ids(1:4) = [0,0,id_el,-id_nue]
    SLim%ids(1:4) = [0,0,id_el,-id_nue]

    HardProc%part = [id_q,-id_qp,id_el,-id_nue,id_a]
    CLim%part = [id_q,-id_qp,id_el,-id_nue]
    SLim%part = [id_q,-id_qp,id_el,-id_nue]

#elif  (_Vcharge == +1)
    HardProc%ids(1:5) = [0,0,id_nue,-id_el,id_a]
    CLim%ids(1:4) = [0,0,id_nue,-id_el]
    SLim%ids(1:4) = [0,0,id_nue,-id_el]

    HardProc%part = [id_q,-id_qp,id_nue,-id_el,id_a]
    CLim%part = [id_q,-id_qp,id_nue,-id_el]
    SLim%part = [id_q,-id_qp,id_nue,-id_el]
#endif


    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintRV_ns(1) = zero

    else
      
       if ( oldcode ) then
          call res_qcdloop_a_qqb(HardProc%AmpMom,res0_old,res1_old)
          call get_respdf(ns_lumi,1,1,HardProc,res1_old,respdf)
       else 
          call ol_res_qcdloop_a_qqb_gen(HardProc%AmpMom,res0,res1)
          call get_respdf_gen(1,1,HardProc,res1,respdf)
       endif

       call partition_nlo_qed(HardProc,damp,icoll)
       
       respdf = respdf*HardProc%wgt*damp

       kin(1) = respdf(1)
       FintRV_ns(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C
    call cut_histo(CLim)
    
    if (CLim%makecut.or.CLim%flag) then

       kin(2) = zero
       FintRV_ns(2) = zero

    else

       if ( oldcode ) then
          call res_qcdloop_qqb(CLim%AmpMom,res0red_old,res1red_old)
          call get_respdf(ns_lumi,1,1,CLim,res1red_old,respdf)
       else
          call res_qcdloop_qqb_gen(CLim%AmpMom,res0red,res1red)
          call get_respdf_gen(1,1,CLim,res1red,respdf)
       endif

       z   = CLim%Lim_KinInv(1)
       s5i = CLim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision
       respdf = (-one)*respdf*(two/s5i)*(Q_lep2*Pqg(z))&
            * CLim%wgt

       FintRV_ns(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- S and CS
    call cut_histo(SLim)
    
    if (SLim%makecut.or.SLim%flag) then

       kin(3) = zero
       FintRV_ns(3:4) = zero

    else

       if ( oldcode ) then
          call res_qcdloop_qqb(SLim%AmpMom,res0red_old,res1red_old)

          !-- prepare PDFs structures, S
          call get_qed_eik(charges_ns,SLim%Lim_etaij,[1,2,3,4],5,eik)
          res_tmp_vect(:,1,1) = res1red_old(:,1) * eik(1:2) !-- dn
          res_tmp_vect(:,2,1) = res1red_old(:,2) * eik(3:4) !-- dn

          !-- prepare PDFs structures, CS
          res_tmp_vect(:,:,2) = res1red_old(:,:) * Q_lep2

          !-- get PDFs for all of them
          call get_respdf_vect(ns_lumi,1,1,SLim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
       else
          call res_qcdloop_qqb_gen(SLim%AmpMom,res0red,res1red)
          call get_qed_eik_gen(res1red,SLim%Lim_etaij,[1,2,3,4],5,res_eikqed)
          call get_respdf_gen(1,1,SLim,res_eikqed,respdf_vect_rev(:,1))
         
          res1red_fs_charge(:,:) = res1red(:,:) * Q_lep2
          call get_respdf_gen(1,1,SLim,res1red_fs_charge,respdf_vect_rev(:,2))
          
          respdf_vect(1,:) = respdf_vect_rev(:,1)
          respdf_vect(2,:) = respdf_vect_rev(:,2)

       endif


       !-- S
       call partition_nlo_qed(SLim,damp,icoll)
       
       e5 = SLim%Lim_KinInv(1)
       respdf_tmp(:,1) = -respdf_vect(1,:)/e5**2 * SLim%wgt * damp
       FintRV_ns(3) = respdf_tmp(1,1)

       !-- CS
       e5    = CSLim%Lim_KinInv(1)
       eta5i = CSLim%Lim_KinInv(2)
       respdf_tmp(:,2) = respdf_vect(2,:)/e5**2/eta5i * CSLim%wgt
       FintRV_ns(4) = respdf_tmp(1,2)

       respdf(:) = sum(respdf_tmp(:,1:2),2)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintRV_ns)

    !if (product(kin) .ne. zero) then
    !   print*,  "hard         ", FintRV_ns(1) 
    !   print *, "collinear    ", FintRV_ns(2), FintRV_ns(2)/FintRV_ns(1)
    !   print *, "soft         ", FintRV_ns(3), FintRV_ns(3)/FintRV_ns(1)
    !   print *, "soft-coll    ", FintRV_ns(4), FintRV_ns(4)/FintRV_ns(1)
    !   print *, "ff",ff(1)
    !   pause
    !endif

#if(_withchecks == 1)
    FintRVQCD_ns(1:4) = FintRV_ns
#endif

  end function xsect_nnlo_rvqcd_fs_5i_ns
  
  !-------------------------------------------------------------
  !-- initial-state sector, aq and qa
  !-------------------------------------------------------------
 
  function xsect_nnlo_rvqcd_is_aq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rvqcd_is_aq
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintRV_aq(3),kin(3)
    real(dp)    :: res0(2,2),res1(2,2)
    real(dp)    :: res1gen(-5:7,-5:7)
    real(dp)    :: res0red(2,2),res1red(2,2),restmp(2,2),res0red_AA
    real(dp)    :: respdf(ipdf)
    real(dp)    :: limcol(2),z,s5i

    xsect_nnlo_rvqcd_is_aq = 0

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

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.false.)
    !-- these must go before the cut_histo routine
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintRV_aq(1) = zero

    else    


       call ol_res_qcdloop_a_aq_gen(HardProc%AmpMom(:,1:5),res0,res1gen)
       call res_qcdloop_a_aq(HardProc%AmpMom(:,1:5),res0,res1)
       call get_respdf(aq_lumi,1,1,HardProc,res1,respdf)

       respdf = respdf*HardProc%wgt
       
       kin(1) = respdf(1)
       FintRV_aq(1) = kin(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintRV_aq(2) = zero

    else

       call res_qcdloop_qqb(C1Lim%AmpMom(:,1:4),res0red,res1red)

       limcol = xn*col_rvqcd_is_qa(C1Lim)

       restmp(:,1) = (res0red(:,1) * limcol(1) * Cf + res1red(:,1) * limcol(2)) * Qdn2
       restmp(:,2) = (res0red(:,2) * limcol(1) * Cf + res1red(:,2) * limcol(2)) * Qup2
       
       call get_respdf(aq_lumi,1,1,C1Lim,restmp,respdf)
       respdf = -respdf*C1Lim%wgt

       FintRV_aq(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(3) = zero
       FintRV_aq(3) = zero

    else

       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- here we use spin-averages since there are no spin correlations for the aa -> e-e+ amplitude
       call res_treeAA_aa(C2Lim%AmpMom,res0red_AA)
       !-- 1=xn*aveqa/aveaa
       !-- use Pqq/(1-z) because of definition of z
       !-- to get 1L splitting function: recall that finite part of FF is -8*Cf
       res0red_AA = res0red_AA * Q_lep2 * Pqq(z)/(one-z) * (-8._dp*Cf) 

       res0red(:,1) = Qdn2 * res0red_AA
       res0red(:,2) = Qup2 * res0red_AA
       
       call get_respdf(aq_lumi,1,1,C2Lim,res0red,respdf)
       
       respdf = (-one)*respdf*(two/s5i)&
            * C2Lim%wgt

       FintRV_aq(3) = respdf(1)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()
    
    call check_ff(ff,xx,FintRV_aq)
    
#if(_withchecks == 1)
    FintRVQCD_aq = FintRV_aq
#endif

  end function xsect_nnlo_rvqcd_is_aq

  function xsect_nnlo_rvqcd_is_qa(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rvqcd_is_qa
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintRV_qa(3),kin(3)
    real(dp)    :: res0(2,2),res1(2,2)
    real(dp)    :: res1gen(-5:7,-5:7)
    real(dp)    :: res0red(2,2),res1red(2,2),restmp(2,2),res0red_AA
    real(dp)    :: respdf(ipdf)
    real(dp)    :: limcol(2),z,s5i

    xsect_nnlo_rvqcd_is_qa = 0

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

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.false.)
    !-- these must go before the cut_histo routine
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintRV_qa(1) = zero

    else    


       call ol_res_qcdloop_a_qa_gen(HardProc%AmpMom(:,1:5),res0,res1gen)
       call res_qcdloop_a_qa(HardProc%AmpMom(:,1:5),res0,res1)
       call get_respdf(qa_lumi,1,1,HardProc,res1,respdf)

       respdf = respdf*HardProc%wgt
       
       kin(1) = respdf(1)
       FintRV_qa(1) = kin(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintRV_qa(2) = zero

    else

       z   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- here we use spin-averages since there are no spin correlations for the aa -> e-e+ amplitude
       call res_treeAA_aa(C1Lim%AmpMom,res0red_AA)
       !-- 1=xn*aveqa/aveaa
       !-- use Pqq/(1-z) because of definition of z
       !-- to get 1L splitting function: recall that finite part of FF is -8*Cf
       res0red_AA = res0red_AA * Q_lep2 * Pqq(z)/(one-z) * (-8._dp*Cf) 

       res0red(:,1) = Qdn2 * res0red_AA
       res0red(:,2) = Qup2 * res0red_AA
       
       call get_respdf(qa_lumi,1,1,C1Lim,res0red,respdf)
       
       respdf = (-one)*respdf*(two/s5i)&
            * C1Lim%wgt

       FintRV_qa(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(3) = zero
       FintRV_qa(3) = zero

    else

       call res_qcdloop_qqb(C2Lim%AmpMom(:,1:4),res0red,res1red)

       limcol = xn*col_rvqcd_is_qa(C2Lim)

       restmp(:,1) = (res0red(:,1) * limcol(1) * Cf + res1red(:,1) * limcol(2)) * Qdn2
       restmp(:,2) = (res0red(:,2) * limcol(1) * Cf + res1red(:,2) * limcol(2)) * Qup2
       
       call get_respdf(qa_lumi,1,1,C2Lim,restmp,respdf)
       respdf = -respdf*C2Lim%wgt

       FintRV_qa(3) = respdf(1)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()
    
    call check_ff(ff,xx,FintRV_qa)
    
#if(_withchecks == 1)
    FintRVQCD_qa = FintRV_qa
#endif

  end function xsect_nnlo_rvqcd_is_qa

end module mod_xsects_nnlo_rvqcd


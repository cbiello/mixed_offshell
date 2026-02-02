module mod_xsects_nloqcd_r
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_limvals
  use mod_proc_parms
  use mod_kinematics_nlo
  use mod_process
  use mod_aux_sectors
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_amplitudes_tree_ppnul
  use mod_splittings_bare
  implicit none
  integer, parameter :: nFint  = 3 !-- 6 if full, 3 if use S(C1+C2)-S = 0
  integer, parameter :: nkin   = 3 !-- 4 if full, 3 if use S(C1+C2)-S = 0
#if(_withchecks == 1)
  real(dp), public, save :: FintNLO_HC1C2(3)
  real(dp), public, save :: FintNLO_HC(2)
#endif

  private

  public :: xsect_nloqcd_r_is_ns
  public :: xsect_nloqcd_r_is_gq,xsect_nloqcd_r_is_qg

  !!!!! W xsect
  public :: xsect_nloqcd_r_is_ns_wp ! deprecated -> remove
  public :: xsect_nloqcd_r_is_ns_wm ! deprecated -> remove

contains

  function xsect_nloqcd_r_is_ns(yRnd,ff,vegasweight)
    integer :: xsect_nloqcd_r_is_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(nFint),kin(nkin)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo_old(2,2),res_lo_old(2,2)
    real(dp)    :: res_nlo(-5:7,-5:7),res_lo(-5:7,-5:7)
    real(dp)    :: z,s5i
    integer     :: i, j
    !--
    logical :: oldcode

    oldcode = .false.

    !----- initialisation
    xsect_nloqcd_r_is_ns = 0
    res_nlo_old(:,:) = zero
    res_lo_old(:,:) = zero
    res_nlo(:,:) = zero
    res_lo(:,:) = zero
    !----- end

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

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.false.)

    ! define process specific partons
#if (_Vcharge == 0)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_g]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]

    HardProc%part = [id_q,-id_q,id_el,-id_el,id_g]
    C1Lim%part = [id_q,-id_q,id_el,-id_el]
    C2Lim%part = [id_q,-id_q,id_el,-id_el]
#elif  (_Vcharge == -1)
    HardProc%ids(1:5) = [0,0,id_el,-id_nue,id_g]
    C1Lim%ids(1:4) = [0,0,id_el,-id_nue]
    C2Lim%ids(1:4) = [0,0,id_el,-id_nue]

    HardProc%part = [id_q,-id_qp,id_el,-id_nue,id_g]
    C1Lim%part = [id_q,-id_qp,id_el,-id_nue]
    C2Lim%part = [id_q,-id_qp,id_el,-id_nue]

#elif  (_Vcharge == +1)
    HardProc%ids(1:5) = [0,0,id_nue,-id_el,id_g]
    C1Lim%ids(1:4) = [0,0,id_nue,-id_el]
    C2Lim%ids(1:4) = [0,0,id_nue,-id_el]

    HardProc%part = [id_q,-id_qp,id_nue,-id_el,id_g]
    C1Lim%part = [id_q,-id_qp,id_nue,-id_el]
    C2Lim%part = [id_q,-id_qp,id_nue,-id_el]
#endif


    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_ns(1) = zero

    else    

       if(oldcode) then
          call res_tree_g_qqb(HardProc%AmpMom,res_nlo_old)
          call get_respdf(ns_lumi,1,0,HardProc,res_nlo_old,respdf)
       else 
          call res_tree_g_qqb_gen(HardProc%AmpMom,res_nlo)
          
          ! print*, 'res_nlo'
          ! do i = -5, 7
          !    do j = -5, 7
          !        write(*,'(2I4,1X,ES24.16)') i, j, res_nlo(i,j)
          !    end do
          ! end do
       
          call get_respdf_gen(1,0,HardProc,res_nlo,respdf)
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
       if(oldcode) then
          call res_tree_qqb(C1Lim%AmpMom,res_lo_old)
          call get_respdf(ns_lumi,1,0,C1Lim,res_lo_old,respdf)
       else
          call res_tree_qqb_gen(C1Lim%AmpMom,res_lo)
          call get_respdf_gen(1,0,C1Lim,res_lo,respdf)
       endif 

       z   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z)/(one-z))&
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
       if(oldcode) then
          call res_tree_qqb(C2Lim%AmpMom,res_lo_old)
          call get_respdf(ns_lumi,1,0,C2Lim,res_lo_old,respdf)
       else
          call res_tree_qqb_gen(C2Lim%AmpMom,res_lo)
          call get_respdf_gen(1,0,C2Lim,res_lo,respdf)
       endif

       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z)/(one-z))&
            * C2Lim%wgt

       FintNLO_ns(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif


    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_ns)
    
#if(_withchecks == 1)
    FintNLO_HC1C2 = FintNLO_ns
    limval_nlo_is = zero
    limval_nlo_is(1) = FintNLO_ns(1)
    limval_nlo_is(3:4) = FintNLO_ns(2:3)  ! coll
#endif

  end function xsect_nloqcd_r_is_ns

  function xsect_nloqcd_r_is_gq(yRnd,ff,vegasweight)
    integer :: xsect_nloqcd_r_is_gq
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_gq(2),kin(2)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo_old(2,2),res_lo_old(2,2),res_nlo(-5:7,-5:7),res_lo(-5:7,-5:7)
    real(dp)    :: z,s5i
    logical     :: oldcode

    xsect_nloqcd_r_is_gq = 0

    oldcode = .false.
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

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C1Lim=C1Lim,compute_etas=.false.)

    ! define process specific partons
#if (_Vcharge == 0)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]

    HardProc%part = [id_g,id_q,id_el,-id_el,id_q]
    C1Lim%part = [-id_q,id_q,id_el,-id_el]
#elif  (_Vcharge == -1)
    HardProc%ids(1:5) = [0,0,id_el,-id_nue,id_q]
    C1Lim%ids(1:4) = [0,0,id_el,-id_nue]

    HardProc%part = [id_g,id_qp,id_el,-id_nue,id_q]
    C1Lim%part = [-id_q,id_qp,id_el,-id_nue]
#elif  (_Vcharge == +1)
    HardProc%ids(1:5) = [0,0,id_nue,-id_el,id_q]
    C1Lim%ids(1:4) = [0,0,id_nue,-id_el]
    
    HardProc%part = [id_g,id_qp,id_nue,-id_el,id_q]
    C1Lim%part = [-id_q,id_qp,id_nue,-id_el]
#endif


    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_gq(1) = zero
       

    else    
       if (oldcode) then
          call res_tree_g_gq(HardProc%AmpMom,res_nlo_old)
          call get_respdf(gq_lumi,1,0,HardProc,res_nlo_old,respdf)
       else
          call res_tree_g_gq_gen(HardProc%AmpMom,res_nlo)          
          call get_respdf_gen(1,0,HardProc,res_nlo,respdf)
       endif

       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintNLO_gq(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintNLO_gq(2) = zero

    else

       if (oldcode) then
          call res_tree_qqb(C1Lim%AmpMom,res_lo_old)
          call get_respdf(gq_lumi,1,0,C1Lim,res_lo_old,respdf)
       else
          call res_tree_qqb_gen(C1Lim%AmpMom,res_lo)
          res_lo = transition('g -> q', 'none', res_lo)
          call get_respdf_gen(1,0,C1Lim,res_lo,respdf)
       endif 

       z   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pgq so that in the soft limit we do not have 1/[1-(1-x)]
       !-- Tr = Cf * aveqg/aveqq
       respdf = (-one)*respdf*(two/s5i)*(tr*Pgq_spav(z)/(one-z)) & 
            * C1Lim%wgt

       FintNLO_gq(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_gq)
    
#if(_withchecks == 1)
    FintNLO_HC = FintNLO_gq
    limval_nlo(1) = FintNLO_gq(1)
    limval_nlo(3) = FintNLO_gq(2)
#endif

  end function xsect_nloqcd_r_is_gq

  function xsect_nloqcd_r_is_qg(yRnd,ff,vegasweight)
    integer :: xsect_nloqcd_r_is_qg
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_qg(2),kin(2)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(-5:7,-5:7),res_lo(-5:7,-5:7)
    !--
    real(dp)    :: res_nlo_old(2,2),res_lo_old(2,2)
    real(dp)    :: z,s5i
    logical     :: oldcode

    
    oldcode = .false.

    xsect_nloqcd_r_is_qg = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    if (xx(xE)*(one-xx(xRHO)).lt.buff_r) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C2Lim=C2Lim,compute_etas=.false.)

#if (_Vcharge == 0)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]

    HardProc%part = [id_q,id_g,id_el,-id_el,id_q]
    C2Lim%part = [id_q,-id_q,id_el,-id_el]
#elif  (_Vcharge == -1)
    HardProc%ids(1:5) = [0,0,id_el,-id_nue,id_q]
    C2Lim%ids(1:4) = [0,0,id_el,-id_nue]

    HardProc%part = [id_qp,id_g,id_el,-id_nue,id_q]
    C2Lim%part = [id_qp,-id_q,id_el,-id_nue]
#elif  (_Vcharge == +1)
    HardProc%ids(1:5) = [0,0,id_nue,-id_el,id_q]
    C2Lim%ids(1:4) = [0,0,id_nue,-id_el]

    HardProc%part = [id_qp,id_g,id_nue,-id_el,id_q]
    C2Lim%part = [id_qp,-id_q,id_nue,-id_el]
#endif

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_qg(1) = zero

    else    

       if (oldcode) then
          call res_tree_g_qg(HardProc%AmpMom,res_nlo_old)
          call get_respdf(qg_lumi,1,0,HardProc,res_nlo_old,respdf)
       else
          call res_tree_g_qg_gen(HardProc%AmpMom,res_nlo)
          call get_respdf_gen(1,0,HardProc,res_nlo,respdf)
       endif


       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintNLO_qg(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(2) = zero
       FintNLO_qg(2) = zero

    else

       if (oldcode) then
          call res_tree_qqb(C2Lim%AmpMom,res_lo_old)
          call get_respdf(qg_lumi,1,0,C2Lim,res_lo_old,respdf)
       else   
          call res_tree_qqb_gen(C2Lim%AmpMom,res_lo)
          res_lo = transition('none', 'g -> q', res_lo)
          call get_respdf_gen(1,0,C2Lim,res_lo,respdf)
       endif

       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Tr*Pgq_spav(z)/(one-z))&
            * C2Lim%wgt

       FintNLO_qg(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_qg)
    
#if(_withchecks == 1)
    FintNLO_HC = FintNLO_qg
#endif

  end function xsect_nloqcd_r_is_qg


!!!!!!!!! W exchanged process !!!!!!!!!!!

function xsect_nloqcd_r_is_ns_wp(yRnd,ff,vegasweight)
    integer :: xsect_nloqcd_r_is_ns_wp
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(nFint),kin(nkin)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(1,2),res_lo(1,2)
    real(dp)    :: z,s5i
    
    xsect_nloqcd_r_is_ns_wp = 0

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

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.false.)
    HardProc%ids(1:5) = [0,0,-id_el,id_nue,id_g]  
    C1Lim%ids(1:4) = [0,0,-id_el,id_nue] 
    C2Lim%ids(1:4) = [0,0,-id_el,id_nue] 

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

       call res_tree_g_qqb_wp(HardProc%AmpMom,res_nlo)
       call get_respdf(qQpb_lumi_wp,1,0,HardProc,res_nlo,respdf)


       !print*, 'res11= ', res_nlo(1,1)
       !print*, 'res12= ', res_nlo(1,2)
       !stop
              
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

       call res_tree_qqb_wp(C1Lim%AmpMom,res_lo)
       call get_respdf(qQpb_lumi_wp,1,0,C1Lim,res_lo,respdf)

       z   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z)/(one-z))&
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

       call res_tree_qqb_wp(C2Lim%AmpMom,res_lo)
       call get_respdf(qQpb_lumi_wp,1,0,C2Lim,res_lo,respdf)

       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z)/(one-z))&
            * C2Lim%wgt

       FintNLO_ns(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    
    call close_histo()


!    if(FintNLO_ns(1) .ne. zero  ) then
!            print*, 'FintNLO_ns(1) = ', FintNLO_ns(1)
!            print*, 'FintNLO_ns(2) = ', FintNLO_ns(2), '', FintNLO_ns(2)/FintNLO_ns(1)
!            print*, 'FintNLO_ns(3) = ', FintNLO_ns(3), '', FintNLO_ns(3)/FintNLO_ns(1)
!
!            pause
!    endif


    call check_ff(ff,xx,FintNLO_ns)

#if(_withchecks == 1)
    FintNLO_HC1C2 = FintNLO_ns
#endif

  end function xsect_nloqcd_r_is_ns_wp





function xsect_nloqcd_r_is_ns_wm(yRnd,ff,vegasweight)
    integer :: xsect_nloqcd_r_is_ns_wm
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(nFint),kin(nkin)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(1,2),res_lo(1,2)
    real(dp)    :: z,s5i

    xsect_nloqcd_r_is_ns_wm = 0

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

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.false.)
    HardProc%ids(1:5) = [0,0,id_el,-id_nue,id_g]  
    C1Lim%ids(1:4) = [0,0,id_el,-id_nue] 
    C2Lim%ids(1:4) = [0,0,id_el,-id_nue] 

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

       call res_tree_g_qqb_wm(HardProc%AmpMom,res_nlo)
       call get_respdf(qQpb_lumi_wm,1,0,HardProc,res_nlo,respdf)
       
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

       call res_tree_qqb_wm(C1Lim%AmpMom,res_lo)
       call get_respdf(qQpb_lumi_wm,1,0,C1Lim,res_lo,respdf)

       z   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z)/(one-z))&
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

       call res_tree_qqb_wm(C2Lim%AmpMom,res_lo)
       call get_respdf(qQpb_lumi_wm,1,0,C2Lim,res_lo,respdf)

       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(z)/(one-z))&
            * C2Lim%wgt

       FintNLO_ns(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    
    call close_histo()

    call check_ff(ff,xx,FintNLO_ns)


!        if(FintNLO_ns(1) .ne. zero  ) then
!            print*, 'FintNLO_ns(1) = ', FintNLO_ns(1)
!            print*, 'FintNLO_ns(2) = ', FintNLO_ns(2), '', FintNLO_ns(2)/FintNLO_ns(1)
!            print*, 'FintNLO_ns(3) = ', FintNLO_ns(3), '', FintNLO_ns(3)/FintNLO_ns(1)
!
!            pause
!    endif


#if(_withchecks == 1)
    FintNLO_HC1C2 = FintNLO_ns
#endif

  end function xsect_nloqcd_r_is_ns_wm


end module mod_xsects_nloqcd_r
  

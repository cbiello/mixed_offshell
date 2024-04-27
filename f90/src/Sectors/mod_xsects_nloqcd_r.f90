module mod_xsects_nloqcd_r
  use mod_types
  use mod_consts_dp
  use mod_parms
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
  public :: xsect_nloqcd_r_is_ns_wp
  public :: xsect_nloqcd_r_is_ns_wm

contains

  function xsect_nloqcd_r_is_ns(yRnd,ff,vegasweight)
    integer :: xsect_nloqcd_r_is_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(nFint),kin(nkin)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(2,2),res_lo(2,2)
    real(dp)    :: z,s5i

    xsect_nloqcd_r_is_ns = 0

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
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_g]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_ns(1) = zero

    else    

       call res_tree_g_qqb(HardProc%AmpMom,res_nlo)
       call get_respdf(ns_lumi,1,0,HardProc,res_nlo,respdf)

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

       call res_tree_qqb(C1Lim%AmpMom,res_lo)
       call get_respdf(ns_lumi,1,0,C1Lim,res_lo,respdf)

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

       call res_tree_qqb(C2Lim%AmpMom,res_lo)
       call get_respdf(ns_lumi,1,0,C2Lim,res_lo,respdf)

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
    real(dp)    :: res_nlo(2,2),res_lo(2,2)
    real(dp)    :: z,s5i

    xsect_nloqcd_r_is_gq = 0

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
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_gq(1) = zero

    else    

       call res_tree_g_gq(HardProc%AmpMom,res_nlo)
       call get_respdf(gq_lumi,1,0,HardProc,res_nlo,respdf)

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

       call res_tree_qqb(C1Lim%AmpMom,res_lo)
       call get_respdf(gq_lumi,1,0,C1Lim,res_lo,respdf)

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
    real(dp)    :: res_nlo(2,2),res_lo(2,2)
    real(dp)    :: z,s5i

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
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_qg(1) = zero

    else    

       call res_tree_g_qg(HardProc%AmpMom,res_nlo)
       call get_respdf(qg_lumi,1,0,HardProc,res_nlo,respdf)

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

       call res_tree_qqb(C2Lim%AmpMom,res_lo)
       call get_respdf(qg_lumi,1,0,C2Lim,res_lo,respdf)

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

       call res_tree_g_qqb_w(HardProc%AmpMom,res_nlo)
       call get_respdf(qQpb_lumi_wp,1,0,HardProc,res_nlo,respdf)
       
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

       call res_tree_qqb_w(C1Lim%AmpMom,res_lo)
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

       call res_tree_qqb_w(C2Lim%AmpMom,res_lo)
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

       call res_tree_g_qqb_w(HardProc%AmpMom,res_nlo)
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

       call res_tree_qqb_w(C1Lim%AmpMom,res_lo)
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

       call res_tree_qqb_w(C2Lim%AmpMom,res_lo)
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

#if(_withchecks == 1)
    FintNLO_HC1C2 = FintNLO_ns
#endif

  end function xsect_nloqcd_r_is_ns_wm


end module mod_xsects_nloqcd_r
  

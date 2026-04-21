module mod_xsects_nnlo_rr_ns
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_aux_sectors
  use mod_process
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_splittings_bare
  use mod_eikonals
  use mod_partitions
  use mod_ol_interface, only: OL_rr_id,gsq_ol,eesq_ol
  use openloops
  use mod_limvals
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_rr_tc_ns(16)
  real(dp), public, save :: FintNNLO_rr_dc_ns(16)
#endif

  !! Electric charges. 
  !! Each column is for the relevant res_jj(i,j)
  !! NB: different assignment wrt to matrix element
  real(dp), parameter :: charges(4,4) = reshape([ & 
       -Qdn,   Qdn,  -Qup,   Qup,     &
        Qdn,  -Qdn,   Qup,  -Qup,     &
        Q_lep, Q_lep, Q_lep, Q_lep,   &
       -Q_lep,-Q_lep,-Q_lep,-Q_lep ], &
       [4,4])

  public :: xsect_nnlo_rr_5161a_ns_ga,xsect_nnlo_rr_5262a_ns_ga
  public :: xsect_nnlo_rr_5161c_ns_ga,xsect_nnlo_rr_5262c_ns_ga
  public :: xsect_nnlo_rr_5162_ns_ga,xsect_nnlo_rr_5261_ns_ga
  public :: xsect_nnlo_rr_5163_ns_ga,xsect_nnlo_rr_5164_ns_ga, &
       xsect_nnlo_rr_5263_ns_ga,xsect_nnlo_rr_5264_ns_ga

  public :: xsect_nnlo_rr_5161a_ns_qqb,xsect_nnlo_rr_5262a_ns_qqb, &
       xsect_nnlo_rr_5161c_ns_qqb,xsect_nnlo_rr_5262c_ns_qqb

  public :: xsect_nnlo_rr_5161a_ns_qq,xsect_nnlo_rr_5262a_ns_qq, &
       xsect_nnlo_rr_5161c_ns_qq,xsect_nnlo_rr_5262c_ns_qq

  public :: xsect_nnlo_rr_ns_qqb_w,xsect_nnlo_rr_ns_qqp_w,xsect_nnlo_rr_ns_qqpb_w

  private

  !! -- q/qb(p1) + qb/q(p2) -> l-(p3) l+(p4) g(p5) a(p6)
  !! -- q/qb(p1) + qb/q(p2) -> l-(p3) l+(p4) q(p5) qb(p6)
  !! -- q/qb(p1) + q/qb(p2) -> l-(p3) l+(p4) q/qb(p5) q/qb(p6)
  !! -- q/qb(p1) + qb/q(p2) -> l-(p3) l+(p4) Q/Qb(p5) Qb/Q(p6)
  !! -- q/qb(p1) + Q/Qb(p2) -> l-(p3) l+(p4) q/qb(p5) Q/Qb(p6)
  !! -- q/qb(p1) + Qb/Q(p2) -> l-(p3) l+(p4) q/qb(p5) Qb/Q(p6)

contains

  function xsect_nnlo_rr_5161a_ns_ga(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161a_ns_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5161a_ns_ga = xsect_nnlo_rr_5i6ia_ns_ga(yRnd,ff,vegasweight,1,2)
    
  end function xsect_nnlo_rr_5161a_ns_ga

  function xsect_nnlo_rr_5161c_ns_ga(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161c_ns_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5161c_ns_ga = xsect_nnlo_rr_5i6ic_ns_ga(yRnd,ff,vegasweight,1,2)
    
  end function xsect_nnlo_rr_5161c_ns_ga

  function xsect_nnlo_rr_5262a_ns_ga(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5262a_ns_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5262a_ns_ga = xsect_nnlo_rr_5i6ia_ns_ga(yRnd,ff,vegasweight,2,1)
    
  end function xsect_nnlo_rr_5262a_ns_ga

  function xsect_nnlo_rr_5262c_ns_ga(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5262c_ns_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5262c_ns_ga = xsect_nnlo_rr_5i6ic_ns_ga(yRnd,ff,vegasweight,2,1)
    
  end function xsect_nnlo_rr_5262c_ns_ga

  !--

  function xsect_nnlo_rr_5162_ns_ga(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5162_ns_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5162_ns_ga = xsect_nnlo_rr_ii_5i6j_ns_ga(yRnd,ff,vegasweight,1,2)
    
  end function xsect_nnlo_rr_5162_ns_ga

  function xsect_nnlo_rr_5261_ns_ga(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5261_ns_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5261_ns_ga = xsect_nnlo_rr_ii_5i6j_ns_ga(yRnd,ff,vegasweight,2,1)
    
  end function xsect_nnlo_rr_5261_ns_ga

  function xsect_nnlo_rr_5163_ns_ga(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5163_ns_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5163_ns_ga = xsect_nnlo_rr_if_5i6k_ns_ga(yRnd,ff,vegasweight,1,2,3,4)
    
  end function xsect_nnlo_rr_5163_ns_ga

  function xsect_nnlo_rr_5164_ns_ga(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5164_ns_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5164_ns_ga = xsect_nnlo_rr_if_5i6k_ns_ga(yRnd,ff,vegasweight,1,2,4,3)
    
  end function xsect_nnlo_rr_5164_ns_ga

  function xsect_nnlo_rr_5263_ns_ga(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5263_ns_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5263_ns_ga = xsect_nnlo_rr_if_5i6k_ns_ga(yRnd,ff,vegasweight,2,1,3,4)
    
  end function xsect_nnlo_rr_5263_ns_ga

    function xsect_nnlo_rr_5264_ns_ga(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5264_ns_ga
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5264_ns_ga = xsect_nnlo_rr_if_5i6k_ns_ga(yRnd,ff,vegasweight,2,1,4,3)
    
  end function xsect_nnlo_rr_5264_ns_ga

  !-- 4q sectors below
  function xsect_nnlo_rr_5161a_ns_qqb(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161a_ns_qqb
    real(dp15) :: yRnd(30), ff(1), vegasweight
    integer, parameter :: iac = 1

    xsect_nnlo_rr_5161a_ns_qqb = xsect_nnlo_rr_5i6iac_ns_4q_qqb(yRnd,ff,vegasweight,1,2,iac)
    
  end function xsect_nnlo_rr_5161a_ns_qqb

  function xsect_nnlo_rr_5262a_ns_qqb(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5262a_ns_qqb
    real(dp15) :: yRnd(30), ff(1), vegasweight
    integer, parameter :: iac = 1

    xsect_nnlo_rr_5262a_ns_qqb = xsect_nnlo_rr_5i6iac_ns_4q_qqb(yRnd,ff,vegasweight,2,1,iac)
    
  end function xsect_nnlo_rr_5262a_ns_qqb

  function xsect_nnlo_rr_5161c_ns_qqb(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161c_ns_qqb
    real(dp15) :: yRnd(30), ff(1), vegasweight
    integer, parameter :: iac = 3

    xsect_nnlo_rr_5161c_ns_qqb = xsect_nnlo_rr_5i6iac_ns_4q_qqb(yRnd,ff,vegasweight,1,2,iac)
    
  end function xsect_nnlo_rr_5161c_ns_qqb

  function xsect_nnlo_rr_5262c_ns_qqb(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5262c_ns_qqb
    real(dp15) :: yRnd(30), ff(1), vegasweight
    integer, parameter :: iac = 3

    xsect_nnlo_rr_5262c_ns_qqb = xsect_nnlo_rr_5i6iac_ns_4q_qqb(yRnd,ff,vegasweight,2,1,iac)
    
  end function xsect_nnlo_rr_5262c_ns_qqb

  !--
  
  function xsect_nnlo_rr_5161a_ns_qq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161a_ns_qq
    real(dp15) :: yRnd(30), ff(1), vegasweight
    integer, parameter :: iac = 1

    xsect_nnlo_rr_5161a_ns_qq = xsect_nnlo_rr_5i6iac_ns_4q_qq(yRnd,ff,vegasweight,1,2,iac)
    
  end function xsect_nnlo_rr_5161a_ns_qq

  function xsect_nnlo_rr_5262a_ns_qq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5262a_ns_qq
    real(dp15) :: yRnd(30), ff(1), vegasweight
    integer, parameter :: iac = 1

    xsect_nnlo_rr_5262a_ns_qq = xsect_nnlo_rr_5i6iac_ns_4q_qq(yRnd,ff,vegasweight,2,1,iac)
    
  end function xsect_nnlo_rr_5262a_ns_qq

  function xsect_nnlo_rr_5161c_ns_qq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161c_ns_qq
    real(dp15) :: yRnd(30), ff(1), vegasweight
    integer, parameter :: iac = 3

    xsect_nnlo_rr_5161c_ns_qq = xsect_nnlo_rr_5i6iac_ns_4q_qq(yRnd,ff,vegasweight,1,2,iac)
    
  end function xsect_nnlo_rr_5161c_ns_qq

  function xsect_nnlo_rr_5262c_ns_qq(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5262c_ns_qq
    real(dp15) :: yRnd(30), ff(1), vegasweight
    integer, parameter :: iac = 3

    xsect_nnlo_rr_5262c_ns_qq = xsect_nnlo_rr_5i6iac_ns_4q_qq(yRnd,ff,vegasweight,2,1,iac)
    
  end function xsect_nnlo_rr_5262c_ns_qq

  !!*************************************************************************!!

  function xsect_nnlo_rr_5i6ia_ns_ga(yRnd,ff,vegasweight,i,j)
    use mod_kinematics_nnlo_tc_is_eu_ac
    integer :: xsect_nnlo_rr_5i6ia_ns_ga,i,j
    real(dp15) :: yRnd(30),ff(1),vegasweight
    real(dp)   :: xx(kNNLO_max_full)
    type(KinConfig) :: HardProc,S5Lim,S6Lim,C6S6Lim, &
         TCLim,TCC6Lim,TCS5Lim,TCC6S5Lim,C6S5Lim, &
         TCS6Lim,TCC6S6Lim,S5S6Lim,C6S5S6Lim,TCS5S6Lim,TCC6S5S6Lim,C6Lim
    integer, parameter :: imax_ipdf = 2, imax_ilim = 4
    real(dp) :: FintNNLO_ns(16),kin(8)
    real(dp) :: res_lo(-5:7,-5:7),res_nlo(-5:7,-5:7),res_nnlo(-5:7,-5:7),res_nlo_tmp(-5:7,-5:7),res_lo_tmp(-5:7,-5:7)
    real(dp) :: res_lo_old(2,2),res_nlo_old(2,2),res_nnlo_old(2,2),res_tmp_old(2,2)
    real(dp) :: respdf(ipdf), res_nlo_eikqed(-5:7,-5:7), res_nlo_ischarges(-5:7,-5:7)
    real(dp) :: res_lo_eikqed(-5:7,-5:7),res_lo_ischarges(-5:7,-5:7)
    real(dp) :: respdf_tmp_part(ipdf,1:4)
    real(dp) :: respdf_vect(imax_ipdf,ipdf),res_tmp_vect(2,2,imax_ipdf), res_tmp_vect_old(2,2,imax_ipdf)
    real(dp) :: respdf_vect_part(ipdf,imax_ilim),respdf_tmp(ipdf,3)
    real(dp) :: e5,e6,eta5i,eta6i,eik_qcd,eik_qed(4),si5,si6,s56
    real(dp) :: z,z5,z6,zi
    real(dp) :: damp
    logical  :: oldcode

    oldcode = .false.

    xsect_nnlo_rr_5i6ia_ns_ga = 0

    ff(1) = zero

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

    limval_nnlo = zero

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()
    
    !! Generate kinematics, isec = 1
    call kinematics_nnlo_ac_5i6iac(&
         xx(1:kNNLO_max_full),i,j,1, &
         HardProc, &
         S5Lim = S5Lim, S6Lim = S6Lim, C6S6Lim = C6S6Lim, &
         TCLim = TCLim, TCC6Lim = TCC6Lim, TCS5Lim = TCS5Lim, &
         TCC6S5Lim = TCC6S5Lim, C6S5Lim = C6S5Lim, TCS6Lim = TCS6Lim, &
         TCC6S6Lim = TCC6S6Lim, S5S6Lim = S5S6Lim, C6S5S6Lim = C6S5S6Lim, &
         TCS5S6Lim = TCS5S6Lim, TCC6S5S6Lim = TCC6S5S6Lim, C6Lim = C6Lim, &
         opt_etas=[3,4])

#if (_Vcharge == 0)
    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a]
    S5Lim%ids(1:5)    = [0,0,id_el,-id_el,id_a]
    S6Lim%ids(1:5)    = [0,0,id_el,-id_el,id_g]
    C6Lim%ids(1:5)    = [0,0,id_el,-id_el,id_g]
    S5S6Lim%ids(1:4)  = [0,0,id_el,-id_el]
    TCLim%ids(1:4)    = [0,0,id_el,-id_el]
    TCS6Lim%ids(1:4)  = [0,0,id_el,-id_el]
    TCS5Lim%ids(1:4)  = [0,0,id_el,-id_el]

    HardProc%part = [id_q,-id_q,id_el,-id_el,id_g,id_a]
    S5Lim%part    = [id_q,-id_q,id_el,-id_el,id_a]
    S6Lim%part    = [id_q,-id_q,id_el,-id_el,id_g]
    C6Lim%part    = [id_q,-id_q,id_el,-id_el,id_g]
    S5S6Lim%part  = [id_q,-id_q,id_el,-id_el]
    TCLim%part    = [id_q,-id_q,id_el,-id_el]
    TCS6Lim%part  = [id_q,-id_q,id_el,-id_el]
    TCS5Lim%part  = [id_q,-id_q,id_el,-id_el]

#elif  (_Vcharge == -1)
    HardProc%ids(1:6) = [0,0,id_el,-id_nue,id_g,id_a]
    S5Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_a]
    S6Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_g]
    C6Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_g]
    S5S6Lim%ids(1:4)  = [0,0,id_el,-id_nue]
    TCLim%ids(1:4)    = [0,0,id_el,-id_nue]
    TCS6Lim%ids(1:4)  = [0,0,id_el,-id_nue]
    TCS5Lim%ids(1:4)  = [0,0,id_el,-id_nue]

    HardProc%part = [id_q,-id_qp,id_el,-id_nue,id_g,id_a]
    S5Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_a]
    S6Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_g]
    C6Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_g]
    S5S6Lim%part  = [id_q,-id_qp,id_el,-id_nue]
    TCLim%part    = [id_q,-id_qp,id_el,-id_nue]
    TCS6Lim%part  = [id_q,-id_qp,id_el,-id_nue]
    TCS5Lim%part  = [id_q,-id_qp,id_el,-id_nue]

#elif  (_Vcharge == +1)
    HardProc%ids(1:6) = [0,0,id_nue,-id_el,id_g,id_a]
    S5Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_a]
    S6Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_g]
    C6Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_g]
    S5S6Lim%ids(1:4)  = [0,0,id_nue,-id_el]
    TCLim%ids(1:4)    = [0,0,id_nue,-id_el]
    TCS6Lim%ids(1:4)  = [0,0,id_nue,-id_el]
    TCS5Lim%ids(1:4)  = [0,0,id_nue,-id_el]

    HardProc%part = [id_q,-id_qp,id_nue,-id_el,id_g,id_a]
    S5Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_a]
    S6Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_g]
    C6Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_g]
    S5S6Lim%part  = [id_q,-id_qp,id_nue,-id_el]
    TCLim%part    = [id_q,-id_qp,id_nue,-id_el]
    TCS6Lim%part  = [id_q,-id_qp,id_nue,-id_el]
    TCS5Lim%part  = [id_q,-id_qp,id_nue,-id_el]

#endif



    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

!    HardProc%ids(1:6) = [0,0,id_el,id_elbar,id_g,id_a]
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNNLO_ns(1) = zero
       
    else

       call partition_nnlo_fact(HardProc,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=i)
   
       if (oldcode) then
          call res_tree_ga_qqb(HardProc%AmpMom,res_nnlo_old)
          call get_respdf(ns_lumi,1,1,Hardproc,res_nnlo_old,respdf)
       else
          call res_tree_ga_qqb_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,Hardproc,res_nnlo,respdf)
       endif 

       respdf = respdf*HardProc%wgt*damp

       kin(1) = respdf(1)
       FintNNLO_ns(1)= respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                             S5                                        !!
    !!-----------------------------------------------------------------------!!

!    S5Lim%ids(1:5) = [0,0,id_el,id_elbar,id_a]
    call cut_histo(S5Lim)

    if (S5Lim%makecut.or.S5Lim%flag) then

       kin(2) = zero
       FintNNLO_ns(2) = zero
       
    else

       call partition_nnlo_fact(S5Lim,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=i)
       
       if (oldcode) then
          call res_tree_a_qqb(S5Lim%AmpMom,res_nlo_old)
          call get_respdf(ns_lumi,1,1,S5Lim,res_nlo_old,respdf)
       else
          call res_tree_a_qqb_gen(S5Lim%AmpMom,res_nlo)
          call get_respdf_gen(1,1,S5Lim,res_nlo,respdf)
       endif

       call get_qcd_eik(Cf,S5Lim%Lim_etaij,[1,2],5,eik_qcd)
       e5 = S5Lim%Lim_KinInv(1)
       respdf = respdf*eik_qcd/e5**2
       
       respdf = -respdf*S5Lim%wgt*damp

       kin(2) = respdf(1)
       FintNNLO_ns(2)= respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                             S6 + C6S6                                 !!
    !!-----------------------------------------------------------------------!!

!    S6Lim%ids(1:5) = [0,0,id_el,id_elbar,id_g]
    call cut_histo(S6Lim)

    if (S6Lim%makecut.or.S6Lim%flag) then

       kin(3) = zero
       FintNNLO_ns(3:4) = zero
       
    else

       if ( oldcode ) then
          !-- prepare PDFs structures, S6
          call get_qed_eik(charges,S6Lim%Lim_etaij,[1,2,3,4],6,eik_qed)
          res_tmp_vect_old(1:2,1,1) = res_nlo_old(1:2,1) * eik_qed(1:2)
          res_tmp_vect_old(1:2,2,1) = res_nlo_old(1:2,2) * eik_qed(3:4)
          !-- prepare PDFs structures, CS6
          res_tmp_vect_old(:,1,2) = res_nlo_old(:,1) * Qdn2
          res_tmp_vect_old(:,2,2) = res_nlo_old(:,2) * Qup2
          !-- get PDFs for all of them
          call get_respdf_vect(ns_lumi,1,1,S6Lim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
       
       else
          call res_tree_g_qqb_gen(S6Lim%AmpMom,res_nlo_tmp)
          call get_qed_eik_gen(res_nlo_tmp,S6Lim%Lim_etaij,[1,2,3,4],6,res_nlo_eikqed)
          call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf_tmp(:,1))
          res_nlo_ischarges = multiply_IS_charges_sq(res_nlo_tmp,i)
          call get_respdf_gen(1,1,S6Lim,res_nlo_ischarges,respdf_tmp(:,2))
                 
       endif 

       !-- S6
       call partition_nnlo_fact(S6Lim,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=i)

       e6 = S6Lim%Lim_KinInv(1)
       respdf_tmp(:,1) = -respdf_tmp(:,1)/e6**2*S6Lim%wgt*damp
       FintNNLO_ns(3) = respdf_tmp(1,1)

       !-- C6S6
       call partition_nnlo_fact(S6Lim,damp,iconf_qcd=[1,2,5],&
            i_qcd=i)

       e6    = C6S6Lim%Lim_KinInv(1)
       eta6i = C6S6Lim%Lim_etaij(i,6)

       respdf_tmp(:,2) = respdf_tmp(:,2)/e6**2/eta6i*C6S6Lim%wgt*damp
       FintNNLO_ns(4) = respdf_tmp(1,2)

       respdf = sum(respdf_tmp(:,1:2),2)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                             C6                                        !!
    !!-----------------------------------------------------------------------!!

!    C6Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]
    call cut_histo(C6Lim)

    if (C6Lim%makecut.or.C6Lim%flag) then

       kin(4) = zero
       FintNNLO_ns(5) = zero
       
    else
       
       call partition_nnlo_fact(C6Lim,damp,iconf_qcd=[1,2,5],&
            i_qcd=i)

       if ( oldcode ) then 
          call res_tree_g_qqb(C6Lim%AmpMom,res_nlo_old)
          res_nlo_old(:,1) = res_nlo_old(:,1) * Qdn2
          res_nlo_old(:,2) = res_nlo_old(:,2) * Qup2
          call get_respdf(ns_lumi,1,1,C6Lim,res_nlo_old,respdf)
       else 
          call res_tree_g_qqb_gen(C6Lim%AmpMom,res_nlo_tmp)
          res_nlo_ischarges = multiply_IS_charges_sq(res_nlo_tmp,i)
          call get_respdf_gen(1,1,C6Lim,res_nlo_ischarges,respdf)
       endif

       z   = C6Lim%Lim_z(2)
       si6 = C6Lim%Lim_sij(i,6)
       
       respdf = -respdf*(-two*Pqg(z)/si6)*C6Lim%wgt*damp

       FintNNLO_ns(5) = respdf(1)
       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                      S5S6 + TCS5S6 + C6S5S6 + TCC6S5S6                !!
    !!-----------------------------------------------------------------------!!

!    S5S6Lim%ids(1:4) = [0,0,id_el,id_elbar]
    call cut_histo(S5S6Lim)

    if (S5S6Lim%makecut.or.S5S6Lim%flag) then

       kin(5) = zero
       FintNNLO_ns(6:9) = zero
       
    else

       if ( oldcode ) then
          call res_tree_qqb(S5S6Lim%AmpMom,res_lo_old)
          !-- Prepare for PDFs, S5S6
          call get_qcd_eik(Cf,S5S6Lim%Lim_etaij,[1,2],5,eik_qcd)
          call get_qed_eik(charges,S5S6Lim%Lim_etaij,[1,2,3,4],6,eik_qed)
          res_tmp_vect_old(1:2,1,1) = res_lo_old(1:2,1) * eik_qed(1:2)
          res_tmp_vect_old(1:2,2,1) = res_lo_old(1:2,2) * eik_qed(3:4)
          !-- Prepare for PDFs, TCS5S6
          res_tmp_vect_old(:,1,2) = res_lo_old(:,1) * Qdn2
          res_tmp_vect_old(:,2,2) = res_lo_old(:,2) * Qup2
          call get_respdf_vect(ns_lumi,1,1,S5S6Lim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
       else 
          call res_tree_qqb_gen(S5S6Lim%AmpMom,res_lo_tmp)
          !-- Prepare for PDFs, S5S6
          call get_qcd_eik(Cf,S5S6Lim%Lim_etaij,[1,2],5,eik_qcd) 
          call get_qed_eik_gen(res_lo_tmp,S5S6Lim%Lim_etaij,[1,2,3,4],6,res_lo_eikqed)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo_tmp,i)
          call get_respdf_gen(1,1,S5S6Lim,res_lo_eikqed,respdf_tmp(:,1))
          call get_respdf_gen(1,1,S5S6Lim,res_lo_ischarges,respdf_tmp(:,2))
       endif 
       
       !-- S5S6
       call partition_nnlo_fact(S5S6Lim,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=i) 
       e5 = S5S6Lim%Lim_KinInv(1)
       e6 = S5S6Lim%Lim_KinInv(2)

       respdf_tmp_part(:,1) = respdf_tmp(:,1)*eik_qcd/e5**2/e6**2*S5S6Lim%wgt*damp
       FintNNLO_ns(6) = respdf_tmp_part(1,1)

       !-- TCS5S6
       si5 = TCS5S6Lim%Lim_sij(i,5)
       si6 = TCS5S6Lim%Lim_sij(i,6)
       !
       z5 = TCS5S6Lim%Lim_z(1)
       z6 = TCS5S6Lim%Lim_z(2)

       respdf_tmp_part(:,2) = -respdf_tmp(:,2)*Cf*Pgaq_ds(si5,si6,z5,z6)*TCS5S6Lim%wgt
       FintNNLO_ns(7) = respdf_tmp_part(1,2)

       !-- C6S5S6
       call partition_nnlo_fact(S5S6Lim,damp,iconf_qcd=[1,2,5],i_qcd=i)
       e5    = C6S5S6Lim%Lim_KinInv(1)
       e6    = C6S5S6Lim%Lim_KinInv(2)
       eta6i = C6S5S6Lim%Lim_etaij(i,6)

       respdf_tmp_part(:,3) = -respdf_tmp(:,2)*eik_qcd/e5**2/eta6i/e6**2*C6S5S6Lim%wgt*damp
       FintNNLO_ns(8) = respdf_tmp_part(1,3)

       !-- TCC6S5S6
       e5    = TCC6S5S6Lim%Lim_KinInv(1)
       e6    = TCC6S5S6Lim%Lim_KinInv(2)
       eta5i = TCC6S5S6Lim%Lim_etaij(i,5)
       eta6i = TCC6S5S6Lim%Lim_etaij(i,6)

       respdf_tmp_part(:,4) = respdf_tmp(:,2)*Cf/eta5i/e5**2/eta6i/e6**2*TCC6S5S6Lim%wgt
       FintNNLO_ns(9) = respdf_tmp_part(1,4)

       respdf = sum(respdf_tmp_part(:,1:4),2)
       kin(5) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                              TC + TCC6                                !!
    !!-----------------------------------------------------------------------!!

!    TCLim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(TCLim)

    if (TCLim%makecut.or.TCLim%flag) then

       kin(6) = zero
       FintNNLO_ns(10:11) = zero
       
    else

       if ( oldcode ) then 
          call res_tree_qqb(TCLim%AmpMom,res_lo_old)
          res_tmp_old(:,1) = res_lo_old(:,1) * Qdn2
          res_tmp_old(:,2) = res_lo_old(:,2) * Qup2
          call get_respdf(ns_lumi,1,1,TCLim,res_tmp_old,respdf)    
       else 
          call res_tree_qqb_gen(TCLim%AmpMom,res_lo_tmp)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo_tmp,i)
          call get_respdf_gen(1,1,TCLim,res_lo_ischarges,respdf_tmp(:,1))
       endif 


       !-- TC
       si5 = TCLim%Lim_sij(i,5)
       si6 = TCLim%Lim_sij(i,6)
       s56 = TCLim%Lim_sij(5,6)

       z5 = TCLim%Lim_z(1)
       z6 = TCLim%Lim_z(2)
       zi = TCLim%Lim_z(3)

       respdf_vect_part(:,1) = -respdf_tmp(:,1)*Cf*Pgaq(s56,-si5,-si6,z5,z6,zi)*TCLim%wgt
       FintNNLO_ns(10) = respdf_vect_part(1,1)

       !-- TCC6
       z6 = TCC6Lim%Lim_z(2)
       z5 = TCC6Lim%Lim_z(1)
       si5 = TCC6Lim%Lim_KinInv(1)
       
       respdf_vect_part(:,2) = respdf_tmp(:,1)*(Cf*two*Pqg(z5)/si5)*(two*Pqg(z6)/si6)*TCC6Lim%wgt
       FintNNLO_ns(11) = respdf_vect_part(1,2)
       
       respdf = sum(respdf_vect_part(:,1:2),2)
       kin(6) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                            TCS6 + TCC6S6                              !!
    !!-----------------------------------------------------------------------!!

!    TCS6Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(TCS6Lim)

    if (TCS6Lim%makecut.or.TCS6Lim%flag) then

       kin(7) = zero
       FintNNLO_ns(12:13) = zero
       
    else

       if ( oldcode ) then
          call res_tree_qqb(TCS6Lim%AmpMom,res_lo_old)
          res_tmp_old(:,1) = res_lo_old(:,1) * Qdn2
          res_tmp_old(:,2) = res_lo_old(:,2) * Qup2
          call get_respdf(ns_lumi,1,1,TCS6Lim,res_tmp_old,respdf)
       else
          call res_tree_qqb_gen(TCS6Lim%AmpMom,res_lo_tmp)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo_tmp,i)
          call get_respdf_gen(1,1,TCS6Lim,res_lo_ischarges,respdf_tmp(:,1))
       endif

       !-- TCS6
       z5  = TCS6Lim%Lim_z(1)
       z6  = TCS6Lim%Lim_z(2)
       
       si5 = TCS6Lim%Lim_sij(i,5)
       si6 = TCS6Lim%Lim_sij(i,6)

       respdf_tmp_part(:,1) = respdf_tmp(:,1)*Pgaq_s2(si5,si6,z5,z6)*Cf*TCS6Lim%wgt
       FintNNLO_ns(12) = respdf_tmp_part(1,1)

       !-- TCC6S6
       e6 = TCC6S6Lim%Lim_KinInv(1)
       eta6i = TCC6S6Lim%Lim_etaij(i,6)
       
       respdf_tmp_part(:,2) = respdf_tmp(:,1)*(Cf*two*Pqg(z5)/si5)/e6**2/eta6i*TCC6S6Lim%wgt
       FintNNLO_ns(13) = respdf_tmp_part(1,2)

       respdf = sum(respdf_tmp_part(:,1:2),2)
       kin(7) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                     TCS5 + TCC6S5 + C6S5                              !!
    !!-----------------------------------------------------------------------!!

!    TCS5Lim%ids(1:4) = [0,0,id_el,id_elbar]
    call cut_histo(TCS5Lim)

    if (TCS5Lim%makecut.or.TCS5Lim%flag) then

       kin(8) = zero
       FintNNLO_ns(14:16) = zero
       
    else

       if ( oldcode ) then
          call res_tree_qqb(TCS5Lim%AmpMom,res_lo_old)
          res_tmp_old(:,1) = res_lo_old(:,1) * Qdn2
          res_tmp_old(:,2) = res_lo_old(:,2) * Qup2
          call get_respdf(ns_lumi,1,1,TCS5Lim,res_tmp_old,respdf)
       else 
          call res_tree_qqb_gen(TCS5Lim%AmpMom,res_lo_tmp)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo_tmp,i)
          call get_respdf_gen(1,1,TCS5Lim,res_lo_ischarges,respdf_tmp(:,1))
       endif 

       !-- TCS5
       z5  = TCS5Lim%Lim_z(1)
       z6  = TCS5Lim%Lim_z(2)
       
       si5 = TCS5Lim%Lim_sij(i,5)
       si6 = TCS5Lim%Lim_sij(i,6)

       respdf_tmp_part(:,1) = respdf_tmp(:,1)*Pgaq_s2(si6,si5,z6,z5)*Cf*TCS5Lim%wgt
       FintNNLO_ns(14) = respdf_tmp_part(1,1)

       !-- TCC6S5
       e5 = TCC6S5Lim%Lim_KinInv(1)
       eta5i = TCC6S5Lim%Lim_etaij(i,5)
       
       respdf_tmp_part(:,2) = respdf_tmp(:,1)*(Cf/e5**2/eta5i)*(two*Pqg(z6)/si6)*TCC6S5Lim%wgt
       FintNNLO_ns(15) = respdf_tmp_part(1,2)

       !-- C6S5
       call partition_nnlo_fact(C6S5Lim,damp,iconf_qcd=[1,2,5],i_qcd=i) 

       call get_qcd_eik(Cf,C6S5Lim%Lim_etaij,[1,2],5,eik_qcd)
       e5 = C6S5Lim%Lim_KinInv(1)
       
       z6  = C6S5Lim%Lim_z(2)
       si6 = C6S5Lim%Lim_sij(i,6)

       respdf_tmp_part(:,3) = -respdf_tmp(:,1)*(eik_qcd/e5**2)*(two*Pqg(z6)/si6)*C6S5Lim%wgt*damp
       FintNNLO_ns(16) = respdf_tmp_part(1,3)

       !       respdf = sum(respdf_vect_part(:,1:3),2)
       respdf = sum(respdf_tmp_part(:,1:3),2)
       kin(8) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ns)

 !!   if (FintNNLO_ns(1) .ne. zero) then
 !!      print *, "FintNNLO_ns",FintNNLO_ns
 !!      print *, "S5",FintNNLO_ns(2),one + FintNNLO_ns(1) / FintNNLO_ns(2)
 !!      print *, "S6",FintNNLO_ns(3),one + FintNNLO_ns(1) / FintNNLO_ns(3)
 !!      print *, "S6C6",FintNNLO_ns(4),one - FintNNLO_ns(1) / FintNNLO_ns(4)
 !!      print *, "C6",FintNNLO_ns(5),one + FintNNLO_ns(1) / FintNNLO_ns(5)
 !!      print *, "S5S6",FintNNLO_ns(6),one - FintNNLO_ns(1) / FintNNLO_ns(6)
 !!      print *, "TC S5 S6",FintNNLO_ns(7),one + FintNNLO_ns(1) / FintNNLO_ns(7)
 !!      print *, "C6 S5 S6",FintNNLO_ns(8),one + FintNNLO_ns(1) / FintNNLO_ns(8)
 !!      print *, "TC C6 S5 S6",FintNNLO_ns(9),one - FintNNLO_ns(1) / FintNNLO_ns(9)
 !!      print *, "TC",FintNNLO_ns(10),one + FintNNLO_ns(1) / FintNNLO_ns(10)
 !!      print *, "TC C6",FintNNLO_ns(11),one - FintNNLO_ns(1) / FintNNLO_ns(11)
 !!      print *, "TC S6",FintNNLO_ns(12),one - FintNNLO_ns(1) / FintNNLO_ns(12)
 !!      print *, "TC C6 S6",FintNNLO_ns(13),one + FintNNLO_ns(1) / FintNNLO_ns(13)
 !!      print *, "TC S5",FintNNLO_ns(14),one - FintNNLO_ns(1) / FintNNLO_ns(14)
 !!      print *, "TC C6 S5",FintNNLO_ns(15),one + FintNNLO_ns(1) / FintNNLO_ns(15)
 !!      print *, "C6 S5",FintNNLO_ns(16),one - FintNNLO_ns(1) / FintNNLO_ns(16)
 !!      print *, "sum", sum(FintNNLO_ns)/FintNNLO_ns(1)
!!!       stop
 !!   endif

#if(_withchecks == 1)
    FintNNLO_rr_tc_ns = FintNNLO_ns
    limval_nnlo = FintNNLO_ns          ! HARD, S5, S6 --> 1,2,3
    limval_nnlo(10) = FintNNLO_ns(4)    ! S6 C6
    limval_nnlo(5) = FintNNLO_ns(5)    ! C6
    limval_nnlo(6) = FintNNLO_ns(6)    ! S5 S6
    limval_nnlo(12) = FintNNLO_ns(7)    ! S5 S6 TC
    limval_nnlo(13) = FintNNLO_ns(8)    ! S5 S6 C6 
    limval_nnlo(16) = FintNNLO_ns(9)    ! S5 S6 TC C6 
    limval_nnlo(4) = FintNNLO_ns(10)    ! TC
    limval_nnlo(11) = FintNNLO_ns(11)    ! TC C6
    limval_nnlo(9) = FintNNLO_ns(12)    ! S6 TC
    limval_nnlo(15) = FintNNLO_ns(13)    ! S6 C6 TC
    limval_nnlo(7) = FintNNLO_ns(14)    ! S5 TC
    limval_nnlo(14) = FintNNLO_ns(15)    ! S5 TC C6
    limval_nnlo(8) = FintNNLO_ns(16)    ! S5 C6
#endif

  end function xsect_nnlo_rr_5i6ia_ns_ga

  function xsect_nnlo_rr_5i6ic_ns_ga(yRnd,ff,vegasweight,i,j)
    use mod_kinematics_nnlo_tc_is_eu_ac
    integer :: xsect_nnlo_rr_5i6ic_ns_ga,i,j
    real(dp15) :: yRnd(30),ff(1),vegasweight
    real(dp)   :: xx(kNNLO_max_full)
    type(KinConfig) :: HardProc,S5Lim,S6Lim,&
       TCLim,TCC5Lim,TCS5Lim,TCC5S5Lim,TCS6Lim,TCC5S6Lim,C5S6Lim, &
       S5S6Lim,C5S5S6Lim,TCS5S6Lim,TCC5S5S6Lim,C5Lim,C5S5Lim
    integer, parameter :: imax_ipdf = 3, imax_ilim = 4
    real(dp) :: FintNNLO_ns(16),kin(9)
    real(dp) :: res_lo_old(2,2),res_nlo_old(2,2),res_nnlo_old(2,2),res_tmp_old(2,2)
    real(dp) :: res_nnlo(-5:7,-5:7), res_nlo(-5:7,-5:7), res_lo(-5:7,-5:7), res_nlo_eikqed(-5:7,-5:7), res_nlo_ischarges(-5:7,-5:7), res_lo_eikqed(-5:7,-5:7), res_lo_ischarges(-5:7,-5:7)
    real(dp) :: respdf(ipdf), respdf_tmp(ipdf,3), respdf_vect1(ipdf), respdf_vect2(ipdf), respdf_vect3(ipdf)
    real(dp) :: respdf_vect(imax_ipdf,ipdf),res_tmp_vect_old(2,2,imax_ipdf)
    real(dp) :: respdf_vect_part(imax_ilim,ipdf)
    real(dp) :: e5,e6,eta5i,eta6i,eik_qcd,eik_qed(4),si5,si6,s56
    real(dp) :: z,z5,z6,zi
    real(dp) :: damp
    logical  :: oldcode

    xsect_nnlo_rr_5i6ic_ns_ga = 0

    ff(1) = zero

    oldcode = .false.

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))
    
#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()
    
    !! Generate kinematics, isec = 3
    call kinematics_nnlo_ac_5i6iac(&
         xx(1:kNNLO_max_full),i,j,3, &
         HardProc, &
         S5Lim = S5Lim, S6Lim = S6Lim, C5S6Lim = C5S6Lim, &
         TCLim = TCLim, TCC5Lim = TCC5Lim, TCS5Lim = TCS5Lim, &
         TCC5S5Lim = TCC5S5Lim, C5S5Lim = C5S5Lim, TCS6Lim = TCS6Lim, &
         TCC5S6Lim = TCC5S6Lim, S5S6Lim = S5S6Lim, C5S5S6Lim = C5S5S6Lim, &
         TCS5S6Lim = TCS5S6Lim, TCC5S5S6Lim = TCC5S5S6Lim, C5Lim = C5Lim, &
         opt_etas=[3,4])

#if (_Vcharge == 0)
    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a]
    S5Lim%ids(1:5)    = [0,0,id_el,-id_el,id_a]
    C5S5Lim%ids(1:5)  = [0,0,id_el,-id_el,id_a]
    S6Lim%ids(1:5)    = [0,0,id_el,-id_el,id_g]
    C5Lim%ids(1:5)    = [0,0,id_el,-id_el,id_a]
    S5S6Lim%ids(1:4)  = [0,0,id_el,-id_el]    
    TCLim%ids(1:4)    = [0,0,id_el,-id_el]
    TCS6Lim%ids(1:4)  = [0,0,id_el,-id_el]
    TCS5Lim%ids(1:4)  = [0,0,id_el,-id_el]

    HardProc%part = [id_q,-id_q,id_el,-id_el,id_g,id_a]
    S5Lim%part    = [id_q,-id_q,id_el,-id_el,id_a]
    C5S5Lim%part  = [id_q,-id_q,id_el,-id_el,id_a]
    S6Lim%part    = [id_q,-id_q,id_el,-id_el,id_g]
    C5Lim%part    = [id_q,-id_q,id_el,-id_el,id_a]
    S5S6Lim%part  = [id_q,-id_q,id_el,-id_el]
    TCLim%part    = [id_q,-id_q,id_el,-id_el]
    TCS6Lim%part  = [id_q,-id_q,id_el,-id_el]
    TCS5Lim%part  = [id_q,-id_q,id_el,-id_el]

#elif  (_Vcharge == -1)
    HardProc%ids(1:6) = [0,0,id_el,-id_nue,id_g,id_a]
    S5Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_a]
    C5S5Lim%ids(1:5)  = [0,0,id_el,-id_nue,id_a]
    S6Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_g]
    C5Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_a]
    S5S6Lim%ids(1:4)  = [0,0,id_el,-id_nue]
    TCLim%ids(1:4)    = [0,0,id_el,-id_nue]
    TCS6Lim%ids(1:4)  = [0,0,id_el,-id_nue]
    TCS5Lim%ids(1:4)  = [0,0,id_el,-id_nue]
    
    HardProc%part = [id_q,-id_qp,id_el,-id_nue,id_g,id_a]
    S5Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_a]
    C5S5Lim%part  = [id_q,-id_qp,id_el,-id_nue,id_a]
    S6Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_g]
    C5Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_a]
    S5S6Lim%part  = [id_q,-id_qp,id_el,-id_nue]
    TCLim%part    = [id_q,-id_qp,id_el,-id_nue]
    TCS6Lim%part  = [id_q,-id_qp,id_el,-id_nue]
    TCS5Lim%part  = [id_q,-id_qp,id_el,-id_nue]

#elif  (_Vcharge == +1)
    HardProc%ids(1:6) = [0,0,id_nue,-id_el,id_g,id_a]
    S5Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_a]
    C5S5Lim%ids(1:5)  = [0,0,id_nue,-id_el,id_a]
    S6Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_g]
    C5Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_a]
    S5S6Lim%ids(1:4)  = [0,0,id_nue,-id_el]
    TCLim%ids(1:4)    = [0,0,id_nue,-id_el]
    TCS6Lim%ids(1:4)  = [0,0,id_nue,-id_el]
    TCS5Lim%ids(1:4)  = [0,0,id_nue,-id_el]
    
    HardProc%part = [id_q,-id_qp,id_nue,-id_el,id_g,id_a]
    S5Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_a]
    C5S5Lim%part  = [id_q,-id_qp,id_nue,-id_el,id_a]
    S6Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_g]
    C5Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_a]
    S5S6Lim%part  = [id_q,-id_qp,id_nue,-id_el]
    TCLim%part    = [id_q,-id_qp,id_nue,-id_el]
    TCS6Lim%part  = [id_q,-id_qp,id_nue,-id_el]
    TCS5Lim%part  = [id_q,-id_qp,id_nue,-id_el]

#endif

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    ! HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a]
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNNLO_ns(1) = zero
       
    else

       call partition_nnlo_fact(HardProc,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=i)

       if (oldcode) then
          call res_tree_ga_qqb(HardProc%AmpMom,res_nnlo_old)
          call get_respdf(ns_lumi,1,1,Hardproc,res_nnlo_old,respdf)
       else
          call res_tree_ga_qqb_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,Hardproc,res_nnlo,respdf)
       endif

       respdf = respdf*HardProc%wgt*damp

       kin(1) = respdf(1)
       FintNNLO_ns(1)= respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                             S5                                        !!
    !!-----------------------------------------------------------------------!!

    ! S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]
    call cut_histo(S5Lim)

    if (S5Lim%makecut.or.S5Lim%flag) then

       kin(2) = zero
       FintNNLO_ns(2) = zero
       
    else

       call partition_nnlo_fact(S5Lim,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=i)

       if (oldcode) then
          call res_tree_a_qqb(S5Lim%AmpMom,res_nlo_old)
          call get_respdf(ns_lumi,1,1,S5Lim,res_nlo_old,respdf)
       else
          call res_tree_a_qqb_gen(S5Lim%AmpMom,res_nlo)
          call get_respdf_gen(1,1,S5Lim,res_nlo,respdf)
       endif

       call get_qcd_eik(Cf,S5Lim%Lim_etaij,[1,2],5,eik_qcd)
       e5 = S5Lim%Lim_KinInv(1)
       respdf = respdf*eik_qcd/e5**2
       
       respdf = -respdf*S5Lim%wgt*damp

       kin(2) = respdf(1)
       FintNNLO_ns(2)= respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif
    
    !!-----------------------------------------------------------------------!!
    !!                             C5S5                                      !!
    !!-----------------------------------------------------------------------!!

    C5S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]
    call cut_histo(C5S5Lim)

    if (C5S5Lim%makecut.or.C5S5Lim%flag) then

       kin(3) = zero
       FintNNLO_ns(3) = zero
       
    else

       call partition_nnlo_fact(C5S5Lim,damp,iconf_qed=[1,2,3,4,6],&
            i_qed=i)

       if (oldcode) then
          call res_tree_a_qqb(C5S5Lim%AmpMom,res_nlo_old)
          call get_respdf(ns_lumi,1,1,C5S5Lim,res_nlo_old,respdf)
       else
          call res_tree_a_qqb_gen(C5S5Lim%AmpMom,res_nlo)
          call get_respdf_gen(1,1,C5S5Lim,res_nlo,respdf)
       endif

       e5    = C5S5Lim%Lim_KinInv(1)
       eta5i = C5S5Lim%Lim_etaij(i,5)
       
       respdf = respdf*Cf/e5**2/eta5i*C5S5Lim%wgt*damp

       kin(3) = respdf(1)
       FintNNLO_ns(3)= respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                                 S6                                    !!
    !!-----------------------------------------------------------------------!!

    ! S6Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]
    call cut_histo(S6Lim)

    if (S6Lim%makecut.or.S6Lim%flag) then

       kin(4) = zero
       FintNNLO_ns(4) = zero
       
    else
       if (oldcode) then
          call res_tree_g_qqb(S6Lim%AmpMom,res_nlo_old)
       
          !-- prepare PDFs structures, S6
          call get_qed_eik(charges,S6Lim%Lim_etaij,[1,2,3,4],6,eik_qed)
          res_tmp_old(1:2,1) = res_nlo_old(1:2,1) * eik_qed(1:2)
          res_tmp_old(1:2,2) = res_nlo_old(1:2,2) * eik_qed(3:4)
          call get_respdf(ns_lumi,1,1,S6Lim,res_tmp_old,respdf)
       else
          call res_tree_g_qqb_gen(S6Lim%AmpMom,res_nlo)
          call get_qed_eik_gen(res_nlo,S6Lim%Lim_etaij,[1,2,3,4],6,res_nlo_eikqed)
          call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf)
       endif

       !-- S6
       call partition_nnlo_fact(S6Lim,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=i)

       e6 = S6Lim%Lim_KinInv(1)
       respdf = -respdf/e6**2*S6Lim%wgt*damp

       FintNNLO_ns(4) = respdf(1)
       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                             C5                                        !!
    !!-----------------------------------------------------------------------!!

    ! C5Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]
    call cut_histo(C5Lim)

    if (C5Lim%makecut.or.C5Lim%flag) then

       kin(5) = zero
       FintNNLO_ns(5) = zero
       
    else
       
       call partition_nnlo_fact(C5Lim,damp,iconf_qed=[1,2,3,4,6],&
            i_qed=i)

       if (oldcode) then
          call res_tree_a_qqb(C5Lim%AmpMom,res_nlo_old)
          call get_respdf(ns_lumi,1,1,C5Lim,res_nlo_old,respdf)
       else
          call res_tree_a_qqb_gen(C5Lim%AmpMom,res_lo)
          call get_respdf_gen(1,1,C5Lim,res_lo,respdf)
       endif

       z   = C5Lim%Lim_z(1)
       si5 = C5Lim%Lim_sij(i,5)
       
       respdf = -respdf*(-two*Cf*Pqg(z)/si5)*C5Lim%wgt*damp

       FintNNLO_ns(5) = respdf(1)
       kin(5) = respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                      S5S6 + TCS5S6 + C5S5S6 + TCC5S5S6                !!
    !!-----------------------------------------------------------------------!!

    ! S5S6Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(S5S6Lim)

    if (S5S6Lim%makecut.or.S5S6Lim%flag) then

       kin(6) = zero
       FintNNLO_ns(6:9) = zero
       
    else

       if (oldcode) then
          call res_tree_qqb(S5S6Lim%AmpMom,res_lo_old)

          !-- Prepare for PDFs, S5S6
          call get_qcd_eik(Cf,S5S6Lim%Lim_etaij,[1,2],5,eik_qcd) 
          call get_qed_eik(charges,S5S6Lim%Lim_etaij,[1,2,3,4],6,eik_qed)
          res_tmp_vect_old(1:2,1,1) = res_lo_old(1:2,1) * eik_qed(1:2) 
          res_tmp_vect_old(1:2,2,1) = res_lo_old(1:2,2) * eik_qed(3:4) 

          !-- Prepare for PDFs, TCS5S6
          res_tmp_vect_old(:,1,2) = res_lo_old(:,1) * Qdn2
          res_tmp_vect_old(:,2,2) = res_lo_old(:,2) * Qup2

          !-- Prepare for PDFs, C5S5S
          call get_qed_eik(charges,C5S5S6Lim%Lim_etaij,[1,2,3,4],6,eik_qed)
          
          res_tmp_vect_old(1:2,1,3) = res_lo_old(1:2,1) * eik_qed(1:2) 
          res_tmp_vect_old(1:2,2,3) = res_lo_old(1:2,2) * eik_qed(3:4) 
              
          call get_respdf_vect(ns_lumi,1,1,S5S6Lim,res_tmp_vect_old(:,:,1:3),respdf_vect(1:3,:))
       else
          
          call res_tree_qqb_gen(S5S6Lim%AmpMom,res_lo)

          !-- Prepare for PDFs, S5S6 
          call get_qcd_eik(Cf,S5S6Lim%Lim_etaij,[1,2],5,eik_qcd)
          call get_qed_eik_gen(res_lo,S5S6Lim%Lim_etaij,[1,2,3,4],6,res_lo_eikqed)
          !-- Prepare for PDFs, TCS5S6
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,i)
          
          call get_respdf_gen(1,1,S5S6Lim,res_lo_eikqed,respdf_vect1(:))      ! S5 S6
          call get_respdf_gen(1,1,S5S6Lim,res_lo_ischarges,respdf_vect2(:))   ! TC S5 S6

          ! Prepare for PDFS, C5 S5S6
          call get_qed_eik_gen(res_lo,C5S5S6Lim%Lim_etaij,[1,2,3,4],6,res_lo_eikqed)
          call get_respdf_gen(1,1,S5S6Lim,res_lo_eikqed,respdf_vect3(:))      ! C5 S5 S6
          respdf_vect(1,:) = respdf_vect1(:)
          respdf_vect(2,:) = respdf_vect2(:)
          respdf_vect(3,:) = respdf_vect3(:)

       endif
       
       !-- S5S6
       call partition_nnlo_fact(S5S6Lim,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=i) 
       e5 = S5S6Lim%Lim_KinInv(1)
       e6 = S5S6Lim%Lim_KinInv(2)

       respdf_vect_part(1,:) = respdf_vect(1,:)*eik_qcd/e5**2/e6**2*S5S6Lim%wgt*damp
       FintNNLO_ns(6) = respdf_vect_part(1,1)

       !-- TCS5S6
       si5 = TCS5S6Lim%Lim_sij(i,5)
       si6 = TCS5S6Lim%Lim_sij(i,6)
       !
       z5 = TCS5S6Lim%Lim_z(1)
       z6 = TCS5S6Lim%Lim_z(2)

       respdf_vect_part(2,:) = -respdf_vect(2,:)*Cf*Pgaq_ds(si5,si6,z5,z6)*TCS5S6Lim%wgt
       FintNNLO_ns(7) = respdf_vect_part(2,1)

       !-- C5S5S6
       call partition_nnlo_fact(C5S5S6Lim,damp,iconf_qed=[1,2,3,4,6],i_qed=i)
       e5    = C5S5S6Lim%Lim_KinInv(1)
       e6    = C5S5S6Lim%Lim_KinInv(2)
       eta5i = C5S5S6Lim%Lim_KinInv(3)
       
       respdf_vect_part(3,:) = -respdf_vect(3,:)/e6**2*Cf/eta5i/e5**2*C5S5S6Lim%wgt*damp
       FintNNLO_ns(8) = respdf_vect_part(3,1)

       !-- TCC5S5S6
       e5    = TCC5S5S6Lim%Lim_KinInv(1)
       e6    = TCC5S5S6Lim%Lim_KinInv(2)
       eta5i = TCC5S5S6Lim%Lim_etaij(i,5)
       eta6i = TCC5S5S6Lim%Lim_etaij(i,6)

       respdf_vect_part(4,:) = respdf_vect(2,:)*Cf/eta5i/e5**2/eta6i/e6**2*TCC5S5S6Lim%wgt
       FintNNLO_ns(9) = respdf_vect_part(4,1)

       respdf = sum(respdf_vect_part(1:4,:),1)
       kin(6) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif


    !!-----------------------------------------------------------------------!!
    !!                              TC + TCC5                                !!
    !!-----------------------------------------------------------------------!!

    ! TCLim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(TCLim)

    if (TCLim%makecut.or.TCLim%flag) then

       kin(7) = zero
       FintNNLO_ns(10:11) = zero
       
    else

       if (oldcode) then
          call res_tree_qqb(TCLim%AmpMom,res_lo_old)

          res_tmp_old(:,1) = res_lo_old(:,1) * Qdn2
          res_tmp_old(:,2) = res_lo_old(:,2) * Qup2
          call get_respdf(ns_lumi,1,1,TCLim,res_tmp_old,respdf)
       else
          
          call res_tree_qqb_gen(TCLim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,i)
          call get_respdf_gen(1,1,TCLim,res_lo_ischarges,respdf)
       endif


       !-- TC
       si5 = TCLim%Lim_sij(i,5)
       si6 = TCLim%Lim_sij(i,6)
       s56 = TCLim%Lim_sij(5,6)

       z5 = TCLim%Lim_z(1)
       z6 = TCLim%Lim_z(2)
       zi = TCLim%Lim_z(3)

       respdf_vect_part(1,:) = -respdf*Cf*Pgaq(s56,-si5,-si6,z5,z6,zi)*TCLim%wgt
       FintNNLO_ns(10) = respdf_vect_part(1,1)

       !-- TCC5
       z6 = TCC5Lim%Lim_z(2)
       z5 = TCC5Lim%Lim_z(1)
       si6 = TCC5Lim%Lim_KinInv(1)
       
       respdf_vect_part(2,:) = respdf*(Cf*two*Pqg(z5)/si5)*(two*Pqg(z6)/si6)*TCC5Lim%wgt
       FintNNLO_ns(11) = respdf_vect_part(2,1)
       
       respdf = sum(respdf_vect_part(1:2,:),1)
       kin(7) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                         TCS6 + TCC5S6 + C5S6                          !!
    !!-----------------------------------------------------------------------!!

    ! TCS6Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(TCS6Lim)

    if (TCS6Lim%makecut.or.TCS6Lim%flag) then

       kin(8) = zero
       FintNNLO_ns(12:14) = zero
       
    else

       if (oldcode) then
          call res_tree_qqb(TCS6Lim%AmpMom,res_lo_old)

          !-- prepare PDFs structures, TCS6 + tCC5S6
          res_tmp_vect_old(:,1,1) = res_lo_old(:,1) * Qdn2
          res_tmp_vect_old(:,2,1) = res_lo_old(:,2) * Qup2

          !-- prepare PDFs structures, C5S6
          call get_qed_eik(charges,C5S6Lim%Lim_etaij,[1,2,3,4],6,eik_qed)
          res_tmp_vect_old(1:2,1,2) = res_lo_old(1:2,1) * eik_qed(1:2) 
          res_tmp_vect_old(1:2,2,2) = res_lo_old(1:2,2) * eik_qed(3:4) 

          call get_respdf_vect(ns_lumi,1,1,TCS6Lim,res_tmp_vect_old(:,:,1:2),respdf_vect(1:2,:))
       else
          call res_tree_qqb_gen(TCS6Lim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,i)
          call get_respdf_gen(1,1,TCS6Lim,res_lo_ischarges,respdf_vect1(:))
          call get_qed_eik_gen(res_lo,C5S6Lim%Lim_etaij,[1,2,3,4],6,res_lo_eikqed)
          call get_respdf_gen(1,1,TCS6Lim,res_lo_eikqed,respdf_vect2(:))
          respdf_vect(1,:) = respdf_vect1(:)
          respdf_vect(2,:) = respdf_vect2(:)
       endif

       !-- TCS6
       z5  = TCS6Lim%Lim_z(1)
       z6  = TCS6Lim%Lim_z(2)
       
       si5 = TCS6Lim%Lim_sij(i,5)
       si6 = TCS6Lim%Lim_sij(i,6)

       respdf_vect_part(1,:) = respdf_vect(1,:)*Pgaq_s2(si5,si6,z5,z6)*Cf*TCS6Lim%wgt
       FintNNLO_ns(12) = respdf_vect_part(1,1)

       !-- TCC5S6
       e6 = TCC5S6Lim%Lim_KinInv(1)
       eta6i = TCC5S6Lim%Lim_etaij(i,6)
       
       respdf_vect_part(2,:) = respdf_vect(1,:)*(Cf*two*Pqg(z5)/si5)/e6**2/eta6i*TCC5S6Lim%wgt
       FintNNLO_ns(13) = respdf_vect_part(2,1)

       !-- C5S6
       call partition_nnlo_fact(C5S6Lim,damp,iconf_qed=[1,2,3,4,6],&
            i_qed=i) 

       respdf_vect_part(3,:) = respdf_vect(2,:)*(-Cf*two*Pqg(z5)/si5)/e6**2*C5S6Lim%wgt*damp
       FintNNLO_ns(14) = respdf_vect_part(3,1)

       respdf = sum(respdf_vect_part(1:3,:),1)
       kin(8) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                           TCS5 + TCC5S5                               !!
    !!-----------------------------------------------------------------------!!

    ! TCS5Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(TCS5Lim)

    if (TCS5Lim%makecut.or.TCS5Lim%flag) then

       kin(9) = zero
       FintNNLO_ns(15:16) = zero
       
    else

       if (oldcode) then
          call res_tree_qqb(TCS5Lim%AmpMom,res_lo_old)

          res_tmp_old(:,1) = res_lo_old(:,1) * Qdn2
          res_tmp_old(:,2) = res_lo_old(:,2) * Qup2
          call get_respdf(ns_lumi,1,1,TCS5Lim,res_tmp_old,respdf)
       else
          call res_tree_qqb_gen(TCS5Lim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,i)
          call get_respdf_gen(1,1,TCS5Lim,res_lo_ischarges,respdf)
       endif
       
          
          
       
       !-- TCS5
       z5  = TCS5Lim%Lim_z(1)
       z6  = TCS5Lim%Lim_z(2)
       
       si5 = TCS5Lim%Lim_sij(i,5)
       si6 = TCS5Lim%Lim_sij(i,6)

       !-- Pgaq is symmetric in 1<->2
       respdf_vect_part(1,:) = respdf*Pgaq_s2(si6,si5,z6,z5)*Cf*TCS5Lim%wgt
       FintNNLO_ns(15) = respdf_vect_part(1,1)

       !-- TCC5S5
       e5 = TCC5S5Lim%Lim_KinInv(1)
       eta5i = TCC5S5Lim%Lim_etaij(i,5)
       
       respdf_vect_part(2,:) = respdf*(Cf/e5**2/eta5i)*(two*Pqg(z6)/si6)*TCC5S5Lim%wgt
       FintNNLO_ns(16) = respdf_vect_part(2,1)

       respdf = sum(respdf_vect_part(1:2,:),1)
       kin(9) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif
   

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ns)

#if(_withchecks == 1)
    FintNNLO_rr_tc_ns = FintNNLO_ns
    limval_nnlo(1:2) = FintNNLO_ns(1:2)     ! hard, s5
    limval_nnlo(8) = FintNNLO_ns(3)  ! s5 c5
    limval_nnlo(3) = FintNNLO_ns(4)  ! s6
    limval_nnlo(5) = FintNNLO_ns(5)  ! c5
    limval_nnlo(6) = FintNNLO_ns(6)  ! s5 s6
    limval_nnlo(12) = FintNNLO_ns(7)  ! tc s5 s6
    limval_nnlo(13) = FintNNLO_ns(8)  ! c5 s5 s6
    limval_nnlo(16) = FintNNLO_ns(9)  ! tc c5 s5 s6
    limval_nnlo(4) = FintNNLO_ns(10)  ! tc
    limval_nnlo(11) = FintNNLO_ns(11)  ! tc c5 
    limval_nnlo(9) = FintNNLO_ns(12)  ! tc  s6
    limval_nnlo(15) = FintNNLO_ns(13)  ! tc c5 s6
    limval_nnlo(10) = FintNNLO_ns(14)  ! c5 s6
    limval_nnlo(7) = FintNNLO_ns(15)  ! tc s5
    limval_nnlo(14) = FintNNLO_ns(16)  ! s5 tc c5 
    
#endif
    
  end function xsect_nnlo_rr_5i6ic_ns_ga

  !-- 
  
  function xsect_nnlo_rr_ii_5i6j_ns_ga(yRnd,ff,vegasweight,i,j)
    use mod_kinematics_nnlo_dc_ii
    integer :: xsect_nnlo_rr_ii_5i6j_ns_ga,i,j
    real(dp15) :: yRnd(30),ff(1),vegasweight
    real(dp)   :: xx(kNNLO_max_full)
    type(KinConfig) :: HardProc,S5Lim,C5S5Lim,S6Lim,C6S6Lim,C5Lim, &
         C6Lim,C5C6Lim,C5S6Lim,C5C6S6Lim,C6S5Lim,C5C6S5Lim,S5S6Lim, &
         C5S5S6Lim,C6S5S6Lim,C5C6S5S6Lim 
    integer, parameter :: imax_ipdf = 2, imax_ilim = 4
    real(dp) :: FintNNLO_ns(16),kin(9)
    real(dp) :: res_lo_old(2,2),res_nlo_old(2,2),res_nnlo_old(2,2)
    real(dp) :: res_lo(-5:7,-5:7),res_nlo(-5:7,-5:7),res_nnlo(-5:7,-5:7), res_nlo_eikqed(-5:7,-5:7), res_nlo_ischarges(-5:7,-5:7), res_lo_eikqed(-5:7,-5:7), res_lo_ischarges(-5:7,-5:7)
    real(dp) :: respdf(ipdf)
    real(dp) :: respdf_vect(imax_ipdf,ipdf),res_tmp_vect(2,2,imax_ipdf),respdf_vect_rev(ipdf,imax_ipdf)
    real(dp) :: respdf_vect_part(imax_ilim,ipdf)
    real(dp) :: e5,e6,eik_qcd,eik_qed(4),si5,sj6
    real(dp) :: z5,z6
    real(dp) :: damp
    logical  :: oldcode


    xsect_nnlo_rr_ii_5i6j_ns_ga = 0

    ff(1) = zero

    oldcode = .false.

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))
    
#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()
    
    call kinematics_nnlo_ii_5i6j( &
         xx(1:kNNLO_max_full),i,j, &
         HardProc, &
         S5Lim = S5Lim, C5S5Lim = C5S5Lim, &
         S6Lim = S6Lim, C6S6Lim = C6S6Lim, &
         C5Lim = C5Lim, C6Lim = C6Lim, C5C6Lim = C5C6Lim, &
         C5S6Lim = C5S6Lim, C5C6S6Lim = C5C6S6Lim, &
         C6S5Lim = C6S5Lim, C5C6S5Lim = C5C6S5Lim, &
         S5S6Lim = S5S6Lim, C5S5S6Lim = C5S5S6Lim, &
         C6S5S6Lim = C6S5S6Lim, C5C6S5S6Lim = C5C6S5S6Lim, &
         opt_etas=[3,4])

#if (_Vcharge == 0)
    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a]
    S5Lim%ids(1:5)    = [0,0,id_el,-id_el,id_a]
    S6Lim%ids(1:5)    = [0,0,id_el,-id_el,id_g]
    C5Lim%ids(1:5)    = [0,0,id_el,-id_el,id_a]
    C6Lim%ids(1:5)    = [0,0,id_el,-id_el,id_g]
    C5C6Lim%ids(1:4)  = [0,0,id_el,-id_el]
    C5S6Lim%ids(1:4)  = [0,0,id_el,-id_el]
    C6S5Lim%ids(1:4)  = [0,0,id_el,-id_el]
    S5S6Lim%ids(1:4)  = [0,0,id_el,-id_el]

    HardProc%part = [id_q,-id_q,id_el,-id_el,id_g,id_a]
    S5Lim%part    = [id_q,-id_q,id_el,-id_el,id_a]
    S6Lim%part    = [id_q,-id_q,id_el,-id_el,id_g]
    C5Lim%part    = [id_q,-id_q,id_el,-id_el,id_a]
    C6Lim%part    = [id_q,-id_q,id_el,-id_el,id_g]
    C5C6Lim%part  = [id_q,-id_q,id_el,-id_el]
    C5S6Lim%part  = [id_q,-id_q,id_el,-id_el]
    C6S5Lim%part  = [id_q,-id_q,id_el,-id_el]
    S5S6Lim%part  = [id_q,-id_q,id_el,-id_el]

#elif  (_Vcharge == -1)
    HardProc%ids(1:6) = [0,0,id_el,-id_nue,id_g,id_a]
    S5Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_a]
    S6Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_g]
    C5Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_a]
    C6Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_g]
    C5C6Lim%ids(1:4)  = [0,0,id_el,-id_nue]
    C5S6Lim%ids(1:4)  = [0,0,id_el,-id_nue]
    C6S5Lim%ids(1:4)  = [0,0,id_el,-id_nue]
    S5S6Lim%ids(1:4)  = [0,0,id_el,-id_nue]

    HardProc%part = [id_q,-id_qp,id_el,-id_nue,id_g,id_a]
    S5Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_a]
    S6Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_g]
    C5Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_a]
    C6Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_g]
    C5C6Lim%part  = [id_q,-id_qp,id_el,-id_nue]
    C5S6Lim%part  = [id_q,-id_qp,id_el,-id_nue]
    C6S5Lim%part  = [id_q,-id_qp,id_el,-id_nue]
    S5S6Lim%part  = [id_q,-id_qp,id_el,-id_nue]

#elif  (_Vcharge == +1)
    HardProc%ids(1:6) = [0,0,id_nue,-id_el,id_g,id_a]
    S5Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_a]
    S6Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_g]
    C5Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_a]
    C6Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_g]
    C5C6Lim%ids(1:4)  = [0,0,id_nue,-id_el]
    C5S6Lim%ids(1:4)  = [0,0,id_nue,-id_el]
    C6S5Lim%ids(1:4)  = [0,0,id_nue,-id_el]
    S5S6Lim%ids(1:4)  = [0,0,id_nue,-id_el]

    HardProc%part = [id_q,-id_qp,id_nue,-id_el,id_g,id_a]
    S5Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_a]
    S6Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_g]
    C5Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_a]
    C6Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_g]
    C5C6Lim%part  = [id_q,-id_qp,id_nue,-id_el]
    C5S6Lim%part  = [id_q,-id_qp,id_nue,-id_el]
    C6S5Lim%part  = [id_q,-id_qp,id_nue,-id_el]
    S5S6Lim%part  = [id_q,-id_qp,id_nue,-id_el]

#endif


    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNNLO_ns(1) = zero
       
    else

       call partition_nnlo_fact(HardProc,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=j)

       if (oldcode) then
          call res_tree_ga_qqb(HardProc%AmpMom,res_nnlo_old)
          call get_respdf(ns_lumi,1,1,Hardproc,res_nnlo_old,respdf)
       else
          call res_tree_ga_qqb_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,Hardproc,res_nnlo,respdf)
       endif
       
       respdf = respdf*HardProc%wgt*damp

       kin(1) = respdf(1)
       FintNNLO_ns(1)= respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                             S5 + C5S5                                 !!
    !!-----------------------------------------------------------------------!!

!    S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]
    call cut_histo(S5Lim)

    if (S5Lim%makecut.or.S5Lim%flag) then

       kin(2) = zero
       FintNNLO_ns(2:3) = zero
       
    else

       if (oldcode) then
          call res_tree_a_qqb(S5Lim%AmpMom,res_nlo_old)
          call get_respdf(ns_lumi,1,1,S5Lim,res_nlo_old,respdf)
       else
          call res_tree_a_qqb_gen(S5Lim%AmpMom,res_nlo)
          call get_respdf_gen(1,1,S5Lim,res_nlo,respdf)
       endif

       !-- S5
       call partition_nnlo_fact(S5Lim,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=j)
       
       call get_qcd_eik(Cf,S5Lim%Lim_etaij,[1,2],5,eik_qcd)
       e5 = S5Lim%Lim_KinInv(1)

       respdf_vect_part(1,:) = -respdf*(eik_qcd/e5**2)*S5Lim%wgt*damp
       FintNNLO_ns(2)= respdf_vect_part(1,1)

       !-- C5S5
       call partition_nnlo_fact(S5Lim,damp,iconf_qed=[1,2,3,4,6],&
            i_qed=j)
       
       z5  = C5S5Lim%Lim_z(1)
       si5 = C5S5Lim%Lim_sij(i,5)
       
       respdf_vect_part(2,:) = respdf*(four*Cf/si5/z5)*C5S5Lim%wgt*damp
       FintNNLO_ns(3)= respdf_vect_part(2,1)

       respdf = sum(respdf_vect_part(1:2,:),1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                             S6 + C6S6                                 !!
    !!-----------------------------------------------------------------------!!

 !   S6Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]
    call cut_histo(S6Lim)

    if (S6Lim%makecut.or.S6Lim%flag) then

       kin(3) = zero
       FintNNLO_ns(4:5) = zero
       
    else

       if (oldcode) then
          call res_tree_g_qqb(S6Lim%AmpMom,res_nlo_old)

          !-- prepare PDFs structures, S6
          call get_qed_eik(charges,S6Lim%Lim_etaij,[1,2,3,4],6,eik_qed)
          res_tmp_vect(1:2,1,1) = res_nlo_old(1:2,1) * eik_qed(1:2)
          res_tmp_vect(1:2,2,1) = res_nlo_old(1:2,2) * eik_qed(3:4)

          !-- prepare PDFs structures, C6S6
          res_tmp_vect(:,1,2) = res_nlo_old(:,1) * Qdn2
          res_tmp_vect(:,2,2) = res_nlo_old(:,2) * Qup2

          !-- get PDFs for all of them
          call get_respdf_vect(ns_lumi,1,1,S6Lim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
       else
          call res_tree_g_qqb_gen(S6Lim%AmpMom,res_nlo)
          call get_qed_eik_gen(res_nlo,S6Lim%Lim_etaij,[1,2,3,4],6,res_nlo_eikqed)
          call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf_vect_rev(:,1))
          
          res_nlo_ischarges = multiply_IS_charges_sq(res_nlo,j)
          call get_respdf_gen(1,1,S6Lim,res_nlo_ischarges,respdf_vect_rev(:,2))
          respdf_vect(1,:) = respdf_vect_rev(:,1)
          respdf_vect(2,:) = respdf_vect_rev(:,2)
       endif
          

       !-- S6
       call partition_nnlo_fact(S6Lim,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=j)

       e6 = S6Lim%Lim_KinInv(1)
       respdf_vect_part(1,:) = -respdf_vect(1,:)/e6**2*S6Lim%wgt*damp
       FintNNLO_ns(4) = respdf_vect_part(1,1)

       !-- C6S6
       call partition_nnlo_fact(S6Lim,damp,iconf_qcd=[1,2,5],&
            i_qcd=i)

       z6  = C6S6Lim%Lim_z(2)
       sj6 = C6S6Lim%Lim_sij(j,6)

       respdf_vect_part(2,:) = respdf_vect(2,:)*(four/z6/sj6)*C6S6Lim%wgt*damp
       FintNNLO_ns(5) = respdf_vect_part(2,1)

       respdf = sum(respdf_vect_part(1:2,:),1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                                 C5                                    !!
    !!-----------------------------------------------------------------------!!

!    C5Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]
    call cut_histo(C5Lim)

    if (C5Lim%makecut.or.C5Lim%flag) then

       kin(4) = zero
       FintNNLO_ns(6) = zero
       
    else

       call partition_nnlo_fact(C5Lim,damp,iconf_qed=[1,2,3,4,6],&
            i_qed=j)

       if (oldcode) then
          call res_tree_a_qqb(C5Lim%AmpMom,res_nlo_old)
          call get_respdf(ns_lumi,1,1,C5Lim,res_nlo_old,respdf)
       else
          call res_tree_a_qqb_gen(C5Lim%AmpMom,res_nlo)
          call get_respdf_gen(1,1,C5Lim,res_nlo,respdf)
       endif
       
       z5  = one - C5Lim%Lim_z(1) !!TODO: change this in kinematics
       si5 = C5Lim%Lim_sij(i,5)
       
       respdf = -respdf*(-two*Cf*Pqg(z5)/si5)*C5Lim%wgt*damp
       FintNNLO_ns(6)= respdf(1)

       kin(4) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                                 C6                                    !!
    !!-----------------------------------------------------------------------!!

!    C6Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]
    call cut_histo(C6Lim)

    if (C6Lim%makecut.or.C6Lim%flag) then

       kin(5) = zero
       FintNNLO_ns(7) = zero
       
    else

       call partition_nnlo_fact(C6Lim,damp,iconf_qcd=[1,2,5],&
            i_qcd=i)

       if (oldcode) then
          call res_tree_g_qqb(C6Lim%AmpMom,res_nlo_old)
          res_nlo_old(:,1) = res_nlo_old(:,1) * Qdn2
          res_nlo_old(:,2) = res_nlo_old(:,2) * Qup2
          call get_respdf(ns_lumi,1,1,C6Lim,res_nlo_old,respdf)
       else
          call res_tree_g_qqb_gen(C6Lim%AmpMom,res_nlo)
          res_nlo_ischarges = multiply_IS_charges_sq(res_nlo,j)
          call get_respdf_gen(1,1,C6Lim,res_nlo_ischarges,respdf)
       endif

       z6  = one - C6Lim%Lim_z(2) !!TODO: change this in kinematics
       sj6 = C6Lim%Lim_sij(j,6)

       respdf = -respdf*(-two*Pqg(z6)/sj6)*C6Lim%wgt*damp
       FintNNLO_ns(7)= respdf(1)

       kin(5) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                               C5C6                                    !!
    !!-----------------------------------------------------------------------!!

!    C5C6Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C5C6Lim)

    if (C5C6Lim%makecut.or.C5C6Lim%flag) then

       kin(6) = zero
       FintNNLO_ns(8) = zero
       
    else

       if (oldcode) then
          call res_tree_qqb(C5C6Lim%AmpMom,res_lo_old)
          res_lo_old(:,1) = res_lo_old(:,1) * Qdn2
          res_lo_old(:,2) = res_lo_old(:,2) * Qup2          
          call get_respdf(ns_lumi,1,1,C5C6Lim,res_lo_old,respdf) 
       else
          call res_tree_qqb_gen(C5C6Lim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,j)
          call get_respdf_gen(1,1,C5C6Lim,res_lo_ischarges,respdf)
       endif
       
       z5  = one - C5C6Lim%Lim_z(1) !!TODO: change this in kinematics
       si5 = C5C6Lim%Lim_sij(i,5)
       
       z6  = one - C5C6Lim%Lim_z(2) !!TODO: change this in kinematics
       sj6 = C5C6Lim%Lim_sij(j,6)

       respdf = respdf*(two*Cf*Pqg(z5)/si5)*(two*Pqg(z6)/sj6)*C5C6Lim%wgt
       ! raoul debug -- not sure of this
!       respdf = -respdf
       
       FintNNLO_ns(8)= respdf(1)

       kin(6) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                          C5S6 + C5C6S6                                !!
    !!-----------------------------------------------------------------------!!

!    C5S6Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C5S6Lim)

    if (C5S6Lim%makecut.or.C5S6Lim%flag) then

       kin(7) = zero
       FintNNLO_ns(9:10) = zero
       
    else

       if (oldcode) then
          call res_tree_qqb(C5S6Lim%AmpMom,res_lo_old)

       !-- prepare PDFs structures, S6
       call get_qed_eik(charges,C5S6Lim%Lim_etaij,[1,2,3,4],6,eik_qed)
       res_tmp_vect(1:2,1,1) = res_lo_old(1:2,1) * eik_qed(1:2)
       res_tmp_vect(1:2,2,1) = res_lo_old(1:2,2) * eik_qed(3:4)
       
       !-- prepare PDFs structures, C6S6
       res_tmp_vect(:,1,2) = res_lo_old(:,1) * Qdn2
       res_tmp_vect(:,2,2) = res_lo_old(:,2) * Qup2

       !-- get PDFs for all of them
       call get_respdf_vect(ns_lumi,1,1,C5S6Lim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
    else
       call res_tree_qqb_gen(C5S6Lim%AmpMom,res_lo)
       call get_qed_eik_gen(res_lo,C5S6Lim%Lim_etaij,[1,2,3,4],6,res_lo_eikqed)       
       res_lo_ischarges = multiply_IS_charges_sq(res_lo,j)
       call get_respdf_gen(1,1,C5S6Lim,res_lo_eikqed,respdf_vect_rev(:,1))
       call get_respdf_gen(1,1,C5S6Lim,res_lo_ischarges,respdf_vect_rev(:,2))
       respdf_vect(1,:) = respdf_vect_rev(:,1)
       respdf_vect(2,:) = respdf_vect_rev(:,2)
    endif
    
       !-- C5S6
       call partition_nnlo_fact(C5S6Lim,damp,iconf_qed=[1,2,3,4,6],&
            i_qed=j)
       
       z5  = one - C5S6Lim%Lim_z(1) !!TODO: change in kinematics
       si5 = C5S6Lim%Lim_sij(i,5)
       e6  = C5S6Lim%Lim_KinInv(1)
       
       respdf_vect_part(1,:) = -respdf_vect(1,:)*(two*Cf*Pqg(z5)/si5)/e6**2*C5S6Lim%wgt*damp
       FintNNLO_ns(9) = respdf_vect_part(1,1)

       !-- C5C6S6
       z6  = C5C6S6Lim%Lim_z(2)
       sj6 = C5C6S6Lim%Lim_sij(j,6)

       respdf_vect_part(2,:) = respdf_vect(2,:)*(two*Cf*Pqg(z5)/si5)*(four/z6/sj6)*C5C6S6Lim%wgt
       FintNNLO_ns(10) = respdf_vect_part(2,1)

       respdf = sum(respdf_vect_part(1:2,:),1)
       kin(7) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                          C6S5 + C5C6S5                                !!
    !!-----------------------------------------------------------------------!!

  !  C6S5Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C6S5Lim)

    if (C6S5Lim%makecut.or.C6S5Lim%flag) then

       kin(8) = zero
       FintNNLO_ns(11:12) = zero
       
    else

       if (oldcode) then
          call res_tree_qqb(C6S5Lim%AmpMom,res_lo_old)
          res_lo_old(:,1) = res_lo_old(:,1) * Qdn2
          res_lo_old(:,2) = res_lo_old(:,2) * Qup2

          call get_respdf(ns_lumi,1,1,C6S5Lim,res_lo_old,respdf)
       else
          call res_tree_qqb_gen(C6S5Lim%AmpMom,res_lo)
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,j)
          call get_respdf_gen(1,1,C6S5Lim,res_lo_ischarges,respdf)
       endif


       !-- C6S5
       call partition_nnlo_fact(C6S5Lim,damp,iconf_qcd=[1,2,5],&
            i_qcd=i)

       call get_qcd_eik(Cf,C6S5Lim%Lim_etaij,[1,2],5,eik_qcd)
       e5 = C6S5Lim%Lim_KinInv(1)
       
       z6  = one - C6S5Lim%Lim_z(2) !!TODO: change in kinematics
       sj6 = C6S5Lim%Lim_sij(j,6)

       respdf_vect_part(1,:) = -respdf*(two*Pqg(z6)/sj6)*eik_qcd/e5**2*C6S5Lim%wgt*damp
       FintNNLO_ns(11) = respdf_vect_part(1,1)

       !-- C5C6S5
       z5  = C5C6S5Lim%Lim_z(1)
       si5 = C5C6S5Lim%Lim_sij(i,5)

       respdf_vect_part(2,:) = respdf*(Cf*four/z5/si5)*(two*Pqg(z6)/sj6)*C5C6S5Lim%wgt
       FintNNLO_ns(12) = respdf_vect_part(2,1)

       respdf = sum(respdf_vect_part(1:2,:),1)
       kin(8) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                    S5S6 + C5S5S6 + C6S5S6 + C5C6S5S6                  !!
    !!-----------------------------------------------------------------------!!

!    S5S6Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(S5S6Lim)

    if (S5S6Lim%makecut.or.S5S6Lim%flag) then

       kin(9) = zero
       FintNNLO_ns(13:16) = zero
       
    else

       if (oldcode) then
          call res_tree_qqb(S5S6Lim%AmpMom,res_lo_old)

       !-- prepare PDFs structures, S6
          call get_qed_eik(charges,S5S6Lim%Lim_etaij,[1,2,3,4],6,eik_qed)
          res_tmp_vect(1:2,1,1) = res_lo_old(1:2,1) * eik_qed(1:2)
          res_tmp_vect(1:2,2,1) = res_lo_old(1:2,2) * eik_qed(3:4)
       
          !-- prepare PDFs structures, C6S6
          res_tmp_vect(:,1,2) = res_lo_old(:,1) * Qdn2
          res_tmp_vect(:,2,2) = res_lo_old(:,2) * Qup2

          !-- get PDFs for all of them
          call get_respdf_vect(ns_lumi,1,1,S5S6Lim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))

       else
          call res_tree_qqb_gen(S5S6Lim%AmpMom,res_lo)
          call get_qed_eik_gen(res_lo,S5S6Lim%Lim_etaij,[1,2,3,4],6,res_lo_eikqed)       
          res_lo_ischarges = multiply_IS_charges_sq(res_lo,j)
          call get_respdf_gen(1,1,S5S6Lim,res_lo_eikqed,respdf_vect_rev(:,1))
          call get_respdf_gen(1,1,S5S6Lim,res_lo_ischarges,respdf_vect_rev(:,2))
          respdf_vect(1,:) = respdf_vect_rev(:,1)
          respdf_vect(2,:) = respdf_vect_rev(:,2)
       endif
       
       !-- S5S6
       call partition_nnlo_fact(S5S6Lim,damp,iconf_qed=[1,2,3,4,6],iconf_qcd=[1,2,5],&
            i_qcd=i,i_qed=j)

       call get_qcd_eik(Cf,S5S6Lim%Lim_etaij,[1,2],5,eik_qcd)

       e5  = S5S6Lim%Lim_KinInv(1)
       e6  = S5S6Lim%Lim_KinInv(2)
       
       respdf_vect_part(1,:) = respdf_vect(1,:)*eik_qcd/e5**2/e6**2*S5S6Lim%wgt*damp
       FintNNLO_ns(13) = respdf_vect_part(1,1)

       !-- C5S5S6
       call partition_nnlo_fact(C5S5S6Lim,damp,iconf_qed=[1,2,3,4,6],&
            i_qed=j)

       z5  = C5S5S6Lim%Lim_z(1)
       si5 = C5S5S6Lim%Lim_sij(i,5)

       respdf_vect_part(2,:) = -respdf_vect(1,:)*(four*Cf/z5/si5)/e6**2*C5S5S6Lim%wgt*damp
       FintNNLO_ns(14) = respdf_vect_part(2,1)

       !-- C6S5S6
       call partition_nnlo_fact(C6S5S6Lim,damp,iconf_qcd=[1,2,5],&
            i_qcd=i)

       z6  = C6S5S6Lim%Lim_z(2)
       sj6 = C6S5S6Lim%Lim_sij(j,6)

       respdf_vect_part(3,:) = -respdf_vect(2,:)*(eik_qcd/e5**2)*(four/z6/sj6)*C6S5S6Lim%wgt*damp
       FintNNLO_ns(15) = respdf_vect_part(3,1)

       !-- C5C6S5S6
       respdf_vect_part(4,:) = respdf_vect(2,:)*(four*Cf/z5/si5)*(four/z6/sj6)*C5C6S5S6Lim%wgt
       FintNNLO_ns(16) = respdf_vect_part(4,1)

       respdf = sum(respdf_vect_part(1:4,:),1)
       kin(9) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ns)


#if(_withchecks == 1)
    FintNNLO_rr_dc_ns = FintNNLO_ns
    limval_nnlo(1:2) = FintNNLO_ns(1:2)     ! hard, s5
    limval_nnlo(7) = FintNNLO_ns(3)  ! s5 c5
    limval_nnlo(3) = FintNNLO_ns(4)  ! s6
    limval_nnlo(10) = FintNNLO_ns(5)  ! C6 S6
    limval_nnlo(4) = FintNNLO_ns(6)  ! C5
    limval_nnlo(5) = FintNNLO_ns(7)  ! C6
    limval_nnlo(11) = FintNNLO_ns(8)  ! c5 c6
    limval_nnlo(9) = FintNNLO_ns(9)  ! C5 S6
    limval_nnlo(15) = FintNNLO_ns(10)  ! C5 C6 S6
    limval_nnlo(8) = FintNNLO_ns(11)  ! C6 S5
    limval_nnlo(14) = FintNNLO_ns(12)  ! S5 C5 C6
    limval_nnlo(6) = FintNNLO_ns(13)  ! S5 S6
    limval_nnlo(12) = FintNNLO_ns(14)  ! S5 S6 C5
    limval_nnlo(13) = FintNNLO_ns(15)  ! S5 S6 C6
    limval_nnlo(16) = FintNNLO_ns(16)  ! S5 S6 C5 C6
#endif

  end function xsect_nnlo_rr_ii_5i6j_ns_ga

  function xsect_nnlo_rr_if_5i6k_ns_ga(yRnd,ff,vegasweight,i,j,k,l)
    use mod_kinematics_nnlo_dc_if
    integer :: xsect_nnlo_rr_if_5i6k_ns_ga,i,j,k,l
    real(dp15) :: yRnd(30),ff(1),vegasweight
    real(dp)   :: xx(kNNLO_max_full)
    type(KinConfig) :: HardProc,S5Lim,C5S5Lim,S6Lim,C6S6Lim,C5Lim, &
         C6Lim,C5C6Lim,C5S6Lim,C5C6S6Lim,C6S5Lim,C5C6S5Lim,S5S6Lim, &
         C5S5S6Lim,C6S5S6Lim,C5C6S5S6Lim 
    integer, parameter :: imax_ipdf = 2, imax_ilim = 4
    real(dp) :: FintNNLO_ns(16),kin(9)
    real(dp) :: res_lo_old(2,2),res_nlo_old(2,2),res_nnlo_old(2,2)
    real(dp) :: res_lo(-5:7,-5:7),res_nlo(-5:7,-5:7),res_nnlo(-5:7,-5:7), res_nlo_eikqed(-5:7,-5:7), res_nlo_fscharges(-5:7,-5:7), res_lo_eikqed(-5:7,-5:7), res_lo_fscharges(-5:7,-5:7)
    real(dp) :: respdf(ipdf)
    real(dp) :: respdf_vect(imax_ipdf,ipdf),res_tmp_vect(2,2,imax_ipdf),respdf_vect_rev(ipdf,imax_ipdf)
    real(dp) :: respdf_vect_part(imax_ilim,ipdf)
    real(dp) :: e5,e6,eik_qcd,eik_qed(4),si5,sk6,Qsq_FS(2)
    real(dp) :: z5,z6
    real(dp) :: damp
    logical  :: oldcode
    
    xsect_nnlo_rr_if_5i6k_ns_ga = 0

    ff(1) = zero

    oldcode = .false.

    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))
    
#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()
    
    call kinematics_nnlo_if_5i6k( &
         xx(1:kNNLO_max_full),i,j,k,l, &
         HardProc, &
         S5Lim = S5Lim, C5S5Lim = C5S5Lim, &
         S6Lim = S6Lim, C6S6Lim = C6S6Lim, &
         C5Lim = C5Lim, C6Lim = C6Lim, C5C6Lim = C5C6Lim, &
         C5S6Lim = C5S6Lim, C5C6S6Lim = C5C6S6Lim, &
         C6S5Lim = C6S5Lim, C5C6S5Lim = C5C6S5Lim, &
         S5S6Lim = S5S6Lim, C5S5S6Lim = C5S5S6Lim, &
         C6S5S6Lim = C6S5S6Lim, C5C6S5S6Lim = C5C6S5S6Lim, &
         opt_etas=[3,4])

    Qsq_FS = [Q3**2, Q4**2]

#if (_Vcharge == 0)
    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a]
    S5Lim%ids(1:5)    = [0,0,id_el,-id_el,id_a]
    S6Lim%ids(1:5)    = [0,0,id_el,-id_el,id_g]
    C5Lim%ids(1:5)    = [0,0,id_el,-id_el,id_a]
    C6Lim%ids(1:5)    = [0,0,id_el,-id_el,id_g]
    C5C6Lim%ids(1:4)  = [0,0,id_el,-id_el]
    C5S6Lim%ids(1:4)  = [0,0,id_el,-id_el]
    C6S5Lim%ids(1:4)  = [0,0,id_el,-id_el]
    S5S6Lim%ids(1:4)  = [0,0,id_el,-id_el]

    HardProc%part = [id_q,-id_q,id_el,-id_el,id_g,id_a]
    S5Lim%part    = [id_q,-id_q,id_el,-id_el,id_a]
    S6Lim%part    = [id_q,-id_q,id_el,-id_el,id_g]
    C5Lim%part    = [id_q,-id_q,id_el,-id_el,id_a]
    C6Lim%part    = [id_q,-id_q,id_el,-id_el,id_g]
    C5C6Lim%part  = [id_q,-id_q,id_el,-id_el]
    C5S6Lim%part  = [id_q,-id_q,id_el,-id_el]
    C6S5Lim%part  = [id_q,-id_q,id_el,-id_el]
    S5S6Lim%part  = [id_q,-id_q,id_el,-id_el]

#elif  (_Vcharge == -1)
    HardProc%ids(1:6) = [0,0,id_el,-id_nue,id_g,id_a]
    S5Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_a]
    S6Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_g]
    C5Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_a]
    C6Lim%ids(1:5)    = [0,0,id_el,-id_nue,id_g]
    C5C6Lim%ids(1:4)  = [0,0,id_el,-id_nue]
    C5S6Lim%ids(1:4)  = [0,0,id_el,-id_nue]
    C6S5Lim%ids(1:4)  = [0,0,id_el,-id_nue]
    S5S6Lim%ids(1:4)  = [0,0,id_el,-id_nue]

    HardProc%part = [id_q,-id_qp,id_el,-id_nue,id_g,id_a]
    S5Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_a]
    S6Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_g]
    C5Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_a]
    C6Lim%part    = [id_q,-id_qp,id_el,-id_nue,id_g]
    C5C6Lim%part  = [id_q,-id_qp,id_el,-id_nue]
    C5S6Lim%part  = [id_q,-id_qp,id_el,-id_nue]
    C6S5Lim%part  = [id_q,-id_qp,id_el,-id_nue]
    S5S6Lim%part  = [id_q,-id_qp,id_el,-id_nue]
#elif  (_Vcharge == +1)

    HardProc%ids(1:6) = [0,0,id_nue,-id_el,id_g,id_a]
    S5Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_a]
    S6Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_g]
    C5Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_a]
    C6Lim%ids(1:5)    = [0,0,id_nue,-id_el,id_g]
    C5C6Lim%ids(1:4)  = [0,0,id_nue,-id_el]
    C5S6Lim%ids(1:4)  = [0,0,id_nue,-id_el]
    C6S5Lim%ids(1:4)  = [0,0,id_nue,-id_el]
    S5S6Lim%ids(1:4)  = [0,0,id_nue,-id_el]

    HardProc%part = [id_q,-id_qp,id_nue,-id_el,id_g,id_a]
    S5Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_a]
    S6Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_g]
    C5Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_a]
    C6Lim%part    = [id_q,-id_qp,id_nue,-id_el,id_g]
    C5C6Lim%part  = [id_q,-id_qp,id_nue,-id_el]
    C5S6Lim%part  = [id_q,-id_qp,id_nue,-id_el]
    C6S5Lim%part  = [id_q,-id_qp,id_nue,-id_el]
    S5S6Lim%part  = [id_q,-id_qp,id_nue,-id_el]

#endif

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    !HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a]
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNNLO_ns(1) = zero
       
    else

       call partition_nnlo_fact(HardProc,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=k)

       if (oldcode) then
          call res_tree_ga_qqb(HardProc%AmpMom,res_nnlo_old)
          call get_respdf(ns_lumi,1,1,Hardproc,res_nnlo_old,respdf)
       else
          call res_tree_ga_qqb_gen(HardProc%AmpMom,res_nnlo)
          call get_respdf_gen(1,1,Hardproc,res_nnlo,respdf)
       endif

       respdf = respdf*HardProc%wgt*damp

       kin(1) = respdf(1)
       FintNNLO_ns(1)= respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                             S5 + C5S5                                 !!
    !!-----------------------------------------------------------------------!!

    !S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]
    call cut_histo(S5Lim)

    if (S5Lim%makecut.or.S5Lim%flag) then

       kin(2) = zero
       FintNNLO_ns(2:3) = zero
       
    else

       if (oldcode) then
          call res_tree_a_qqb(S5Lim%AmpMom,res_nlo_old)
          call get_respdf(ns_lumi,1,1,S5Lim,res_nlo_old,respdf)
       else
          call res_tree_a_qqb_gen(S5Lim%AmpMom,res_nlo)
          call get_respdf_gen(1,1,S5Lim,res_nlo,respdf)
       endif

       !-- S5
       call partition_nnlo_fact(S5Lim,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=k)
       
       call get_qcd_eik(Cf,S5Lim%Lim_etaij,[1,2],5,eik_qcd)
       e5 = S5Lim%Lim_KinInv(1)

       respdf_vect_part(1,:) = -respdf*(eik_qcd/e5**2)*S5Lim%wgt*damp
       FintNNLO_ns(2)= respdf_vect_part(1,1)

       !-- C5S5
       call partition_nnlo_fact(S5Lim,damp,iconf_qed=[1,2,3,4,6],&
            i_qed=k)
       
       z5  = C5S5Lim%Lim_z(1)
       si5 = C5S5Lim%Lim_sij(i,5)
       
       respdf_vect_part(2,:) = respdf*(four*Cf/si5/z5)*C5S5Lim%wgt*damp
       FintNNLO_ns(3)= respdf_vect_part(2,1)

       respdf = sum(respdf_vect_part(1:2,:),1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                             S6 + C6S6                                 !!
    !!-----------------------------------------------------------------------!!

    !S6Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]
    call cut_histo(S6Lim)

    if (S6Lim%makecut.or.S6Lim%flag) then

       kin(3) = zero
       FintNNLO_ns(4:5) = zero
       
    else

       if (oldcode) then
          call res_tree_g_qqb(S6Lim%AmpMom,res_nlo_old)

          !-- prepare PDFs structures, S6
          call get_qed_eik(charges,S6Lim%Lim_etaij,[1,2,3,4],6,eik_qed)
          res_tmp_vect(1:2,1,1) = res_nlo_old(1:2,1) * eik_qed(1:2)
          res_tmp_vect(1:2,2,1) = res_nlo_old(1:2,2) * eik_qed(3:4)
          
          !-- prepare PDFs structures, C6S6
          res_tmp_vect(:,:,2) = res_nlo_old(:,:) * Q_lep2

          !-- get PDFs for all of them
          call get_respdf_vect(ns_lumi,1,1,S6Lim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
       else
          call res_tree_g_qqb_gen(S6Lim%AmpMom,res_nlo)
          call get_qed_eik_gen(res_nlo,S6Lim%Lim_etaij,[1,2,3,4],6,res_nlo_eikqed)
          call get_respdf_gen(1,1,S6Lim,res_nlo_eikqed,respdf_vect_rev(:,1))
          res_nlo_fscharges = Qsq_Fs(k-2)**2 * res_nlo
          call get_respdf_gen(1,1,S6Lim,res_nlo_fscharges,respdf_vect_rev(:,2))
          respdf_vect(1,:) = respdf_vect_rev(:,1)
          respdf_vect(2,:) = respdf_vect_rev(:,2)
          
       endif
          
       !-- S6
       call partition_nnlo_fact(S6Lim,damp,iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],&
            i_qcd=i,i_qed=k)

       e6 = S6Lim%Lim_KinInv(1)
       respdf_vect_part(1,:) = -respdf_vect(1,:)/e6**2*S6Lim%wgt*damp
       FintNNLO_ns(4) = respdf_vect_part(1,1)

       !-- C6S6
       call partition_nnlo_fact(S6Lim,damp,iconf_qcd=[1,2,5],&
            i_qcd=i)

       z6  = C6S6Lim%Lim_z(2)
       sk6 = C6S6Lim%Lim_sij(k,6)

       respdf_vect_part(2,:) = respdf_vect(2,:)*(four/z6/sk6)*C6S6Lim%wgt*damp
       FintNNLO_ns(5) = respdf_vect_part(2,1)

       respdf = sum(respdf_vect_part(1:2,:),1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                                 C5                                    !!
    !!-----------------------------------------------------------------------!!

    !C5Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]
    call cut_histo(C5Lim)

    if (C5Lim%makecut.or.C5Lim%flag) then

       kin(4) = zero
       FintNNLO_ns(6) = zero
       
    else

       call partition_nnlo_fact(C5Lim,damp,iconf_qed=[1,2,3,4,6],&
            i_qed=k)
       if (oldcode) then
          call res_tree_a_qqb(C5Lim%AmpMom,res_nlo_old)
          call get_respdf(ns_lumi,1,1,C5Lim,res_nlo_old,respdf)
       else
          call res_tree_a_qqb_gen(C5Lim%AmpMom,res_nlo)
          call get_respdf_gen(1,1,C5Lim,res_nlo,respdf)
       endif

       z5  = C5Lim%Lim_z(1)
       si5 = C5Lim%Lim_sij(i,5)
       
       respdf = -respdf*(-two*Cf*Pqg(z5)/si5)*C5Lim%wgt*damp
       FintNNLO_ns(6)= respdf(1)

       kin(4) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                                 C6                                    !!
    !!-----------------------------------------------------------------------!!

    !C6Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]
    call cut_histo(C6Lim)

    if (C6Lim%makecut.or.C6Lim%flag) then

       kin(5) = zero
       FintNNLO_ns(7) = zero
       
    else

       call partition_nnlo_fact(C6Lim,damp,iconf_qcd=[1,2,5],&
            i_qcd=i)

       if (oldcode) then
          call res_tree_g_qqb(C6Lim%AmpMom,res_nlo_old)
          res_nlo_old = res_nlo_old * Q_lep2
          call get_respdf(ns_lumi,1,1,C6Lim,res_nlo_old,respdf)
       else
          call res_tree_g_qqb_gen(C6Lim%AmpMom,res_nlo)
          res_nlo_fscharges = Qsq_Fs(k-2)**2 * res_nlo
          call get_respdf_gen(1,1,C6Lim,res_nlo_fscharges,respdf)
       endif

          

       z6  = C6Lim%Lim_z(2)
       sk6 = C6Lim%Lim_sij(k,6)

       respdf = -respdf*(two*Pqg(z6)/sk6)*C6Lim%wgt*damp
       FintNNLO_ns(7)= respdf(1)

       kin(5) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                               C5C6                                    !!
    !!-----------------------------------------------------------------------!!

    !C5C6Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C5C6Lim)

    if (C5C6Lim%makecut.or.C5C6Lim%flag) then

       kin(6) = zero
       FintNNLO_ns(8) = zero
       
    else
       if (oldcode) then
          call res_tree_qqb(C5C6Lim%AmpMom,res_lo_old)
          res_lo_old = res_lo_old * Q_lep2

          call get_respdf(ns_lumi,1,1,C5C6Lim,res_lo_old,respdf)
       else
          call res_tree_qqb_gen(C5C6Lim%AmpMom,res_lo)
          res_lo_fscharges = Qsq_Fs(k-2)**2 * res_lo
          call get_respdf_gen(1,1,C5C6Lim,res_lo_fscharges,respdf)
       endif
       
       z5  = C5C6Lim%Lim_z(1)
       si5 = C5C6Lim%Lim_sij(i,5)
       
       z6  = C5C6Lim%Lim_z(2)
       sk6 = C5C6Lim%Lim_sij(k,6)

       respdf = -respdf*(two*Cf*Pqg(z5)/si5)*(two*Pqg(z6)/sk6)*C5C6Lim%wgt
       FintNNLO_ns(8)= respdf(1)

       kin(6) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                          C5S6 + C5C6S6                                !!
    !!-----------------------------------------------------------------------!!

    !C5S6Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C5S6Lim)

    if (C5S6Lim%makecut.or.C5S6Lim%flag) then

       kin(7) = zero
       FintNNLO_ns(9:10) = zero
       
    else

       if (oldcode) then
          call res_tree_qqb(C5S6Lim%AmpMom,res_lo_old)

       !-- prepare PDFs structures, S6
          call get_qed_eik(charges,C5S6Lim%Lim_etaij,[1,2,3,4],6,eik_qed)
          res_tmp_vect(1:2,1,1) = res_lo_old(1:2,1) * eik_qed(1:2)
          res_tmp_vect(1:2,2,1) = res_lo_old(1:2,2) * eik_qed(3:4)
          
          !-- prepare PDFs structures, C6S6
          res_tmp_vect(:,:,2) = res_lo_old(:,:) * Q_lep2

          !-- get PDFs for all of them
          call get_respdf_vect(ns_lumi,1,1,C5S6Lim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
       else
          call res_tree_qqb_gen(C5S6Lim%AmpMom,res_lo)
          call get_qed_eik_gen(res_lo,C5S6Lim%Lim_etaij,[1,2,3,4],6,res_lo_eikqed)
          res_lo_fscharges = Qsq_Fs(k-2)**2 * res_lo
          call get_respdf_gen(1,1,C5S6Lim,res_lo_eikqed,respdf_vect_rev(:,1))
          call get_respdf_gen(1,1,C5S6Lim,res_lo_fscharges,respdf_vect_rev(:,2))
          respdf_vect(1,:) = respdf_vect_rev(:,1)
          respdf_vect(2,:) = respdf_vect_rev(:,2)
       endif
       
       !-- C5S6
       call partition_nnlo_fact(C5S6Lim,damp,iconf_qed=[1,2,3,4,6],&
            i_qed=k)
       
       z5  = C5S6Lim%Lim_z(1)
       si5 = C5S6Lim%Lim_sij(i,5)
       e6  = C5S6Lim%Lim_KinInv(1)
       
       respdf_vect_part(1,:) = -respdf_vect(1,:)*(two*Cf*Pqg(z5)/si5)/e6**2*C5S6Lim%wgt*damp
       FintNNLO_ns(9) = respdf_vect_part(1,1)

       !-- C5C6S6
       z6  = C5C6S6Lim%Lim_z(2)
       sk6 = C5C6S6Lim%Lim_sij(k,6)

       respdf_vect_part(2,:) = respdf_vect(2,:)*(two*Cf*Pqg(z5)/si5)*(four/z6/sk6)*C5C6S6Lim%wgt
       FintNNLO_ns(10) = respdf_vect_part(2,1)

       respdf = sum(respdf_vect_part(1:2,:),1)
       kin(7) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                          C6S5 + C5C6S5                                !!
    !!-----------------------------------------------------------------------!!

    !C6S5Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C6S5Lim)

    if (C6S5Lim%makecut.or.C6S5Lim%flag) then

       kin(8) = zero
       FintNNLO_ns(11:12) = zero
       
    else

       if (oldcode) then
          call res_tree_qqb(C6S5Lim%AmpMom,res_lo_old)
          res_lo_old = res_lo_old * Q_lep2
          call get_respdf(ns_lumi,1,1,C6S5Lim,res_lo_old,respdf)
       else
          call res_tree_qqb_gen(C6S5Lim%AmpMom,res_lo)
          res_lo_fscharges = Qsq_Fs(k-2)**2 * res_lo
          call get_respdf_gen(1,1,C6S5Lim,res_lo_fscharges,respdf)
       endif
       


       !-- C6S5
       call partition_nnlo_fact(C6S5Lim,damp,iconf_qcd=[1,2,5],&
            i_qcd=i)

       call get_qcd_eik(Cf,C6S5Lim%Lim_etaij,[1,2],5,eik_qcd)
       e5 = C6S5Lim%Lim_KinInv(1)
       
       z6  = C6S5Lim%Lim_z(2)
       sk6 = C6S5Lim%Lim_sij(k,6)

       respdf_vect_part(1,:) = respdf*(two*Pqg(z6)/sk6)*eik_qcd/e5**2*C6S5Lim%wgt*damp
       FintNNLO_ns(11) = respdf_vect_part(1,1)

       !-- C5C6S5
       z5  = C5C6S5Lim%Lim_z(1)
       si5 = C5C6S5Lim%Lim_sij(i,5)

       respdf_vect_part(2,:) = -respdf*(Cf*four/z5/si5)*(two*Pqg(z6)/sk6)*C5C6S5Lim%wgt
       FintNNLO_ns(12) = respdf_vect_part(2,1)

       respdf = sum(respdf_vect_part(1:2,:),1)
       kin(8) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                    S5S6 + C5S5S6 + C6S5S6 + C5C6S5S6                  !!
    !!-----------------------------------------------------------------------!!

    !S5S6Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(S5S6Lim)

    if (S5S6Lim%makecut.or.S5S6Lim%flag) then

       kin(9) = zero
       FintNNLO_ns(13:16) = zero
       
    else

       if (oldcode) then
          call res_tree_qqb(S5S6Lim%AmpMom,res_lo_old)

          !-- prepare PDFs structures, S5S6
          call get_qed_eik(charges,S5S6Lim%Lim_etaij,[1,2,3,4],6,eik_qed)
          res_tmp_vect(1:2,1,1) = res_lo_old(1:2,1) * eik_qed(1:2)
          res_tmp_vect(1:2,2,1) = res_lo_old(1:2,2) * eik_qed(3:4)
          
          !-- prepare PDFs structures, S5C6S6
          res_tmp_vect(:,:,2) = res_lo_old(:,:) * Q_lep2

          !-- get PDFs for all of them
          call get_respdf_vect(ns_lumi,1,1,S5S6Lim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))

       else
          call res_tree_qqb_gen(S5S6Lim%AmpMom,res_lo)
          call get_qed_eik_gen(res_lo,S5S6Lim%Lim_etaij,[1,2,3,4],6,res_lo_eikqed)
          call get_respdf_gen(1,1,S5S6Lim,res_lo_eikqed,respdf_vect_rev(:,1))
          res_lo_fscharges = Qsq_Fs(k-2)**2 * res_lo
          call get_respdf_gen(1,1,S5S6Lim,res_lo_fscharges,respdf_vect_rev(:,2))
          respdf_vect(1,:) = respdf_vect_rev(:,1)
          respdf_vect(2,:) = respdf_vect_rev(:,2)
       endif
       
       !-- S5S6
       call partition_nnlo_fact(S5S6Lim,damp,iconf_qed=[1,2,3,4,6],iconf_qcd=[1,2,5],&
            i_qcd=i,i_qed=k)

       call get_qcd_eik(Cf,S5S6Lim%Lim_etaij,[1,2],5,eik_qcd)

       e5  = S5S6Lim%Lim_KinInv(1)
       e6  = S5S6Lim%Lim_KinInv(2)
       
       respdf_vect_part(1,:) = respdf_vect(1,:)*eik_qcd/e5**2/e6**2*S5S6Lim%wgt*damp
       FintNNLO_ns(13) = respdf_vect_part(1,1)

       !-- C5S5S6
       call partition_nnlo_fact(C5S5S6Lim,damp,iconf_qed=[1,2,3,4,6],&
            i_qed=k)

       z5  = C5S5S6Lim%Lim_z(1)
       si5 = C5S5S6Lim%Lim_sij(i,5)

       respdf_vect_part(2,:) = -respdf_vect(1,:)*(four*Cf/z5/si5)/e6**2*C5S5S6Lim%wgt*damp
       FintNNLO_ns(14) = respdf_vect_part(2,1)

       !-- C6S5S6
       call partition_nnlo_fact(C6S5S6Lim,damp,iconf_qcd=[1,2,5],&
            i_qcd=i)

       z6  = C6S5S6Lim%Lim_z(2)
       sk6 = C6S5S6Lim%Lim_sij(k,6)

       respdf_vect_part(3,:) = -respdf_vect(2,:)*(eik_qcd/e5**2)*(four/z6/sk6)*C6S5S6Lim%wgt*damp
       FintNNLO_ns(15) = respdf_vect_part(3,1)

       !-- C5C6S5S6
       respdf_vect_part(4,:) = respdf_vect(2,:)*(four*Cf/z5/si5)*(four/z6/sk6)*C5C6S5S6Lim%wgt
       FintNNLO_ns(16) = respdf_vect_part(4,1)

       respdf = sum(respdf_vect_part(1:4,:),1)
       kin(9) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ns)

#if(_withchecks == 1)
    FintNNLO_rr_dc_ns = FintNNLO_ns
    limval_nnlo(1:2) = FintNNLO_ns(1:2)     ! hard, s5
    limval_nnlo(7) = FintNNLO_ns(3)  ! s5 c5
    limval_nnlo(3) = FintNNLO_ns(4)  ! s6
    limval_nnlo(10) = FintNNLO_ns(5)  ! C6 S6
    limval_nnlo(4) = FintNNLO_ns(6)  ! C5
    limval_nnlo(5) = FintNNLO_ns(7)  ! C6
    limval_nnlo(11) = FintNNLO_ns(8)  ! c5 c6
    limval_nnlo(9) = FintNNLO_ns(9)  ! C5 S6
    limval_nnlo(15) = FintNNLO_ns(10)  ! C5 C6 S6
    limval_nnlo(8) = FintNNLO_ns(11)  ! C6 S5
    limval_nnlo(14) = FintNNLO_ns(12)  ! S5 C5 C6
    limval_nnlo(6) = FintNNLO_ns(13)  ! S5 S6
    limval_nnlo(12) = FintNNLO_ns(14)  ! S5 S6 C5
    limval_nnlo(13) = FintNNLO_ns(15)  ! S5 S6 C6
    limval_nnlo(16) = FintNNLO_ns(16)  ! S5 S6 C5 C6
#endif

  end function xsect_nnlo_rr_if_5i6k_ns_ga

  !-- 4q sectors below

  function xsect_nnlo_rr_5i6iac_ns_4q_qqb(yRnd,ff,vegasweight,i,j,ac)
    use mod_kinematics_nnlo_tc_is_eu_ac
    integer :: xsect_nnlo_rr_5i6iac_ns_4q_qqb,i,j,ac
    real(dp15) :: yRnd(30),ff(1),vegasweight
    real(dp)   :: xx(kNNLO_max_full)
    type(KinConfig) :: HardProc,TCLim
    real(dp) :: FintNNLO_ns(2),kin(2)
    real(dp) :: res_lo(2,2),res_nnlo(2,2)
    real(dp) :: respdf(ipdf)
    real(dp) :: si5,si6,s56,z5,z6,zi
    real(dp) :: damp
    
    xsect_nnlo_rr_5i6iac_ns_4q_qqb = 0

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

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()
    
    call kinematics_nnlo_ac_5i6iac(&
         xx(1:kNNLO_max_full),i,j,ac, &
         HardProc,TCLim = TCLim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_q,id_q]
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNNLO_ns(1) = zero
       
    else

       call partition_nnlo_fact(HardProc,damp,iconf_qcd=[1,2,5],i_qcd=i)
       
       call res_tree_qqb_qqb(HardProc%AmpMom,res_nnlo)
       call get_respdf(ns_lumi,1,1,Hardproc,res_nnlo,respdf)

       respdf = respdf*HardProc%wgt*damp

       kin(1) = respdf(1)
       FintNNLO_ns(1)= respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                                   TC                                  !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(TCLim)

    if (TCLim%makecut.or.TCLim%flag) then

       kin(2) = zero
       FintNNLO_ns(2) = zero
       
    else

       call res_tree_qqb(TCLim%AmpMom,res_lo)

       res_lo(:,1) = res_lo(:,1) * Qdn2
       res_lo(:,2) = res_lo(:,2) * Qup2

       !-- TC
       si5 = TCLim%Lim_sij(i,5)
       si6 = TCLim%Lim_sij(i,6)
       s56 = TCLim%Lim_sij(5,6)

       z5 = TCLim%Lim_z(1)
       z6 = TCLim%Lim_z(2)
       zi = TCLim%Lim_z(3)
              
       res_lo(i,:) = res_lo(i,:) * two*Cf*Pqbqq_id_symm(s56,-si5,-si6,z5,z6,zi) !-- qqb -> qqb
       res_lo(j,:) = res_lo(j,:) * two*Cf*Pqbqq_id_symm(s56,-si6,-si5,z6,z5,zi) !-- qbq -> qqb
       
       call get_respdf(ns_lumi,1,1,TCLim,res_lo,respdf)

       respdf = -respdf*TCLim%wgt
       FintNNLO_ns(2) = respdf(1)

       kin(2) = respdf(1)
       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ns)

#if(_withchecks == 1)
    FintNNLO_rr_tc_ns(1:2) = FintNNLO_ns
#endif

  end function xsect_nnlo_rr_5i6iac_ns_4q_qqb

  function xsect_nnlo_rr_5i6iac_ns_4q_qq(yRnd,ff,vegasweight,i,j,ac)
    use mod_kinematics_nnlo_tc_is_eu_ac
    integer :: xsect_nnlo_rr_5i6iac_ns_4q_qq,i,j,ac
    real(dp15) :: yRnd(30),ff(1),vegasweight
    real(dp)   :: xx(kNNLO_max_full)
    type(KinConfig) :: HardProc,TCLim
    real(dp) :: FintNNLO_ns(2),kin(2)
    real(dp) :: res_lo(2,2),res_nnlo(2,2),res_tmp(2,2),split
    real(dp) :: respdf(ipdf)
    real(dp) :: si5,si6,s56,z5,z6,zi
    real(dp) :: damp
    
    xsect_nnlo_rr_5i6iac_ns_4q_qq = 0

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

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()
    
    call kinematics_nnlo_ac_5i6iac(&
         xx(1:kNNLO_max_full),i,j,ac, &
         HardProc,TCLim = TCLim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_q,id_q]
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNNLO_ns(1) = zero
       
    else

       call partition_nnlo_fact(HardProc,damp,iconf_qcd=[1,2,5],i_qcd=i)
       
       call res_tree_qq_qq(HardProc%AmpMom,res_nnlo)
       call get_respdf(qq_lumi,1,1,Hardproc,res_nnlo,respdf)

       respdf = respdf*HardProc%wgt*damp

       kin(1) = respdf(1)
       FintNNLO_ns(1)= respdf(1)

       call fill_histo(respdf,vegasweight)
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                                   TC                                  !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(TCLim)

    if (TCLim%makecut.or.TCLim%flag) then

       kin(2) = zero
       FintNNLO_ns(2) = zero
       
    else

       call res_tree_qqb(TCLim%AmpMom,res_lo)

       res_lo(:,1) = res_lo(:,1) * Qdn2
       res_lo(:,2) = res_lo(:,2) * Qup2

       !-- TC
       si5 = TCLim%Lim_sij(i,5)
       si6 = TCLim%Lim_sij(i,6)
       s56 = TCLim%Lim_sij(5,6)

       z5 = TCLim%Lim_z(1)
       z6 = TCLim%Lim_z(2)
       zi = TCLim%Lim_z(3)

       split = Cf*Pqbqq_id_symm(-si5,-si6,s56,zi,z5,z6) !-- extra factor 1/2 w.r.t qqb case from 1/2!

       if (i.eq.1) then
          res_tmp(i,:) = res_lo(j,:) * split !-- qq channel needs qb q reduced me2
          res_tmp(j,:) = res_lo(i,:) * split
       elseif (i.eq.2) then
          res_tmp = res_lo * split !-- qq channel need q qb reduced me2
       endif
       
       call get_respdf(qq_lumi,1,1,TCLim,res_tmp,respdf)

       respdf = -respdf*TCLim%wgt
       FintNNLO_ns(2) = respdf(1)

       kin(2) = respdf(1)
       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ns)

#if(_withchecks == 1)
    FintNNLO_rr_tc_ns(1:2) = FintNNLO_ns
#endif

  end function xsect_nnlo_rr_5i6iac_ns_4q_qq

  function xsect_nnlo_rr_ns_qqb_w(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_dc_ii
    integer :: xsect_nnlo_rr_ns_qqb_w
    real(dp15) :: yRnd(30),ff(1),vegasweight
    real(dp)   :: xx(kNNLO_max_full)
    type(KinConfig) :: HardProc
    real(dp) :: FintNNLO_ns(1),kin(1),res_nnlo(2,3)
    real(dp15) :: p_ol(4,6), res_ol(2,2)
    real(dp) :: respdf(ipdf)

    xsect_nnlo_rr_ns_qqb_w = 0

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

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! choose a double-collinear type of parametrisation
    call kinematics_nnlo_ii_5i6j(xx(1:kNNLO_max_full),1,2,HardProc)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_q,id_q]
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNNLO_ns(1) = zero
       
    else

       p_ol = HardProc%AmpMom

       call evaluate_tree(OL_rr_id(1),p_ol,res_ol(1,1)) !-- q  qx, q = d,s
       call evaluate_tree(OL_rr_id(2),p_ol,res_ol(2,1)) !-- qx q
       call evaluate_tree(OL_rr_id(3),p_ol,res_ol(1,2)) !-- Q  Qx, Q = u,c
       call evaluate_tree(OL_rr_id(4),p_ol,res_ol(2,2)) !-- Qx Q
       res_ol = res_ol/gsq_ol/eesq_ol
       res_nnlo(1:2,1:2) = res_ol(1:2,1:2)
       res_nnlo(:,3) = zero

       call get_respdf(ns_lumi_splitb,1,1,Hardproc,res_nnlo,respdf)

       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintNNLO_ns(1)= respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ns)

#if(_withchecks == 1)
    FintNNLO_rr_dc_ns(1) = FintNNLO_ns(1)
#endif

  end function xsect_nnlo_rr_ns_qqb_w

  function xsect_nnlo_rr_ns_qqp_w(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_dc_ii
    integer :: xsect_nnlo_rr_ns_qqp_w
    real(dp15) :: yRnd(30),ff(1),vegasweight
    real(dp)   :: xx(kNNLO_max_full)
    type(KinConfig) :: HardProc
    real(dp) :: FintNNLO_ns(1),kin(1),res_nnlo(2,2)
    real(dp15) :: p_ol(4,6), res_ol(2,2)
    real(dp) :: respdf(ipdf)

    xsect_nnlo_rr_ns_qqp_w = 0

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

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! choose a double-collinear type of parametrisation
    call kinematics_nnlo_ii_5i6j(xx(1:kNNLO_max_full),1,2,HardProc)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_q,id_q]
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNNLO_ns(1) = zero
       
    else

       p_ol = HardProc%AmpMom

       call evaluate_tree(OL_rr_id(5),p_ol,res_ol(1,1)) !-- d  u
       call evaluate_tree(OL_rr_id(6),p_ol,res_ol(2,1)) !-- dx ux
       call evaluate_tree(OL_rr_id(7),p_ol,res_ol(1,2)) !-- u  d
       call evaluate_tree(OL_rr_id(8),p_ol,res_ol(2,2)) !-- ux dx
       res_ol = res_ol/gsq_ol/eesq_ol
       res_nnlo = res_ol

       call get_respdf(qQp_lumi,1,1,Hardproc,res_nnlo,respdf)

       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintNNLO_ns(1)= respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ns)

#if(_withchecks == 1)
    FintNNLO_rr_dc_ns(1) = FintNNLO_ns(1)
#endif

  end function xsect_nnlo_rr_ns_qqp_w

  function xsect_nnlo_rr_ns_qqpb_w(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_dc_ii
    integer :: xsect_nnlo_rr_ns_qqpb_w
    real(dp15) :: yRnd(30),ff(1),vegasweight
    real(dp)   :: xx(kNNLO_max_full)
    type(KinConfig) :: HardProc
    real(dp) :: FintNNLO_ns(1),kin(1),res_nnlo(2,2)
    real(dp15) :: p_ol(4,6), res_ol(2,2)
    real(dp) :: respdf(ipdf)

    xsect_nnlo_rr_ns_qqpb_w = 0

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

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    call kinematics_nnlo_ii_5i6j(xx(1:kNNLO_max_full),1,2,HardProc)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_q,id_q]
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNNLO_ns(1) = zero
       
    else

       p_ol = HardProc%AmpMom

       call evaluate_tree(OL_rr_id(9), p_ol,res_ol(1,1)) !-- d  ux
       call evaluate_tree(OL_rr_id(10),p_ol,res_ol(2,1)) !-- dx u
       call evaluate_tree(OL_rr_id(11),p_ol,res_ol(1,2)) !-- u  dx
       call evaluate_tree(OL_rr_id(12),p_ol,res_ol(2,2)) !-- ux d
       res_ol = res_ol/gsq_ol/eesq_ol
       res_nnlo = res_ol

       call get_respdf(qQpb_lumi,1,1,Hardproc,res_nnlo,respdf)

       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintNNLO_ns(1)= respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ns)

#if(_withchecks == 1)
    FintNNLO_rr_dc_ns(1) = FintNNLO_ns(1)
#endif

  end function xsect_nnlo_rr_ns_qqpb_w

end module mod_xsects_nnlo_rr_ns

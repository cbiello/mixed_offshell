module mod_xsects_lo
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_kinematics_gen
  use mod_kinematics_lo
  use mod_process
  use mod_auxfunctions
  use mod_aux_sectors
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppnul
  use mod_amplitudes_tree_ppll
  
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintLO(1)
#endif

  private

  public :: xsect_lo_ns,xsect_lo_aa
  public :: xsect_lo_ns_wp
  public :: xsect_lo_ns_wm

contains

  function xsect_lo_ns(yRnd,ff,vegasweight)
    integer :: xsect_lo_ns
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_lo(-5:7,-5:7)
    real(dp)    :: res_lo_old(2,2)
    logical :: oldcode

    oldcode = .false.
    
    xsect_lo_ns = 0

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

       if(oldcode) then
          call res_tree_qqb(LOProc%AmpMom,res_lo_old)
          call get_respdf(ns_lumi,0,0,LOProc,res_lo_old,respdf)
       else
          call res_tree_qqb_gen(LOProc%AmpMom,res_lo)
          call get_respdf_gen(0,0,LOProc,res_lo,respdf)
       endif

       respdf = respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintLO = kin
#endif

  end function xsect_lo_ns

  function xsect_lo_aa(yRnd,ff,vegasweight)
    integer :: xsect_lo_aa
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_lo(1,1)

    xsect_lo_aa = 0

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

    LOProc%part(1:4) = [22,22,id_el,-id_el]
    
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero

    else

       call res_treeAA_aa(LOProc%AmpMom,res_lo(1,1))
       call get_respdf(aa_lumi,0,0,LOProc,res_lo,respdf)

       respdf = respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintLO = kin
#endif

  end function xsect_lo_aa


  !!!!!!!! W !!!!!!!!!!!!!
  
  function xsect_lo_ns_wp(yRnd,ff,vegasweight)
    integer :: xsect_lo_ns_wp
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_lo(1,2)



    print*, 'I am in xsect_lo_ns_wp'


    stop
    xsect_lo_ns_wp = 0

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

       call res_tree_qqb_wp(LOProc%AmpMom,res_lo)
       
       call get_respdf(qQpb_lumi_wp,0,0,LOProc,res_lo,respdf)

       respdf = respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintLO = kin
#endif

  end function xsect_lo_ns_wp
  
  
  
  function xsect_lo_ns_wm(yRnd,ff,vegasweight)
    integer :: xsect_lo_ns_wm
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_lo(1,2)

    xsect_lo_ns_wm = 0

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

       call res_tree_qqb_wm(LOProc%AmpMom,res_lo)
       call get_respdf(qQpb_lumi_wm,0,0,LOProc,res_lo,respdf)

       respdf = respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintLO = kin
#endif

  end function xsect_lo_ns_wm


end module mod_xsects_lo

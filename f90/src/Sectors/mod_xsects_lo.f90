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
  use mod_amplitudes_tree_ppll
  
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintLO(1)
#endif

  private

  public :: xsect_lo_ns,xsect_lo_aa

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
 if (Vcharge == 0) then
    LOProc%part(1:4) = [id_q,-id_q,id_el,-id_el]
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
 elseif (Vcharge == -1) then
    LOProc%part(1:4) = [id_q,-id_qp,id_el,-id_nue]
    LOProc%ids(1:4) = [0,0,id_el,-id_nue]
 elseif (Vcharge == 1) then
    LOProc%part(1:4) = [id_q,-id_qp,id_nue,-id_el]
    LOProc%ids(1:4) = [0,0,id_nue,-id_el]
 endif

    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero

    else

#if (_oldcode == 1)
          call res_tree_qqb(LOProc%AmpMom,res_lo_old)
          call get_respdf(ns_lumi,0,0,LOProc,res_lo_old,respdf)
#else
          call res_tree_qqb_gen(LOProc%AmpMom,res_lo)
          call get_respdf_gen(0,0,LOProc,res_lo,respdf)
#endif

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
    if (Vcharge /= 0) stop 'xsect_lo_aa: not implemented for W'

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
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    
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

end module mod_xsects_lo

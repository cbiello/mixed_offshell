module mod_xsects_nnlo_vv
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
  use mod_amplitudes_loop_ppll
  use mod_amplitudes_twol_nf_ppll  
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_vv_ns(1)
#endif

  private

  public :: xsect_nnlo_vv_ns_nf !-- closed fermionic loops. 1PI two-loop diagrams only
  public :: xsect_nnlo_vv_ns_fc !-- factorisable contribution
  public :: xsect_nnlo_vv_ns_xf !-- non-factorisable contribution

  public :: xsect_nnlo_vv_ns_nf_full !-- complete closed fermionic loops

  !-- coefficients for scale variation of non-factorisable contribution
  real(dp), parameter :: TClog0 = -one/8 + 29*zeta2 - 30*zeta3 - 44*zeta2**2/5
  real(dp), parameter :: TClog1 = -three/2 + 12*zeta2 - 24*zeta3

contains

  function xsect_nnlo_vv_ns_nf(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_vv_ns_nf
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res0(2,2),res2(2,2)

    xsect_nnlo_vv_ns_nf = 0

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

       call res_twol_nf_qqb(LOProc%AmpMom,res0,res2)
       call get_respdf(ns_lumi,1,1,LOProc,res2,respdf)

       respdf = respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_vv_ns = kin
#endif

  end function xsect_nnlo_vv_ns_nf


  function xsect_nnlo_vv_ns_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_vv_ns_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res_tree(2,3),res_loop(2,3)

    xsect_nnlo_vv_ns_fc = 0

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

    call kinematics_lo(xx,HardProc)
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero

    else

       call res_ewkloop_qqb(HardProc%AmpMom,res_tree,res_loop)
       call get_respdf(ns_lumi_splitb,1,1,HardProc,res_loop,respdf)

       !-- (-8*CF) * FLVfinEWK
       respdf = (-8*CF)*respdf*HardProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_vv_ns = kin
#endif

  end function xsect_nnlo_vv_ns_fc


  function xsect_nnlo_vv_ns_xf(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_vv_ns_xf
    real(dp15) :: yRnd(30),ff(1),vegasweight
#if (_with2L == 1)
    type(KinConfig) :: LOproc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: kin(1),respdf(ipdf),respdf_nnlo(ipdf),respdf_sv(ipdf),sv_logs(ipdf)
    real(dp) :: res_nnlo(2,3),res_lo(2,2)
    real(dp) :: s,t,u
#endif

    xsect_nnlo_vv_ns_xf = 0

    ff(1) = zero

#if (_with2L == 0)
    print *, 'Two-loop amplitudes not linked, please change the config file if you want them. Aborting.'
    vegasweight = 0.d0
    yRnd        = 0.d0
    stop

#elif (_with2L == 1)

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()

    call kinematics_lo(xx,LOproc)
    call cut_histo(LOproc)

    if (LOproc%makecut.or.LOproc%flag) then

       kin(1) = zero

    else

       s =  2*scr(LOproc%AmpMom(:,1),LOproc%AmpMom(:,2))
       t = -2*scr(LOproc%AmpMom(:,1),LOproc%AmpMom(:,3))
       u = -s -t

       call public_interference_2re_h11_h00(s,t,real(mw,kind=dp),real(mz,kind=dp),res_nnlo(1,2),res_nnlo(1,1)) !-- q qb
       call public_interference_2re_h11_h00(s,u,real(mw,kind=dp),real(mz,kind=dp),res_nnlo(2,2),res_nnlo(2,1)) !-- qb q

       res_nnlo(:,3) = 0 !-- b bx not available

       !-- include LO couplings
       res_nnlo = eesq2 * res_nnlo

       !-- extra factor of 1/4 coming from interference
       res_nnlo = CF * res_nnlo/4

       call get_respdf(ns_lumi_splitb,1,1,LOproc,res_nnlo,respdf_nnlo)

       respdf_nnlo = respdf_nnlo * LOproc%wgt

       !-- scale variation
       call res_tree_qqb(LOproc%AmpMom,res_lo)
       res_lo(:,1) = (CF * Qdn2) * res_lo(:,1)
       res_lo(:,2) = (CF * Qup2) * res_lo(:,2)

       call get_respdf(ns_lumi,1,1,LOproc,res_lo,respdf_sv)

       call fill_sv_logs(LOproc%muf(1)**2,s,sv_logs)
       respdf_sv = (TClog0 + TClog1*sv_logs) * respdf_sv * LOproc%wgt

       !-- total
       respdf = respdf_nnlo + respdf_sv

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_vv_ns = kin
#endif

#endif

  end function xsect_nnlo_vv_ns_xf

  function xsect_nnlo_vv_ns_nf_full(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_vv_ns_nf_full
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf)
    real(dp)    :: res0(2,2),res1(2,2),res2(2,2)

    xsect_nnlo_vv_ns_nf_full = 0

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

       call res_ewkloop_qqb_nf(LOProc%AmpMom,res1)
       call res_twol_nf_qqb(LOProc%AmpMom,res0,res2)

       res2 = res2 + (-8*CF)*res1

       call get_respdf(ns_lumi,1,1,LOProc,res2,respdf)

       respdf = respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_vv_ns = kin
#endif

  end function xsect_nnlo_vv_ns_nf_full

end module mod_xsects_nnlo_vv

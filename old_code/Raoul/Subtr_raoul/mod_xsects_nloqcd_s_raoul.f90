module mod_xsects_nloqcd_s_raoul
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_kinematics_lo
  use mod_kinematics_nlo
  use mod_process
  use mod_auxfunctions
  use mod_aux_sectors
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  implicit none
  integer, parameter :: nFint  = 3 !-- 6 if full, 3 if use S(C1+C2)-S = 0
  integer, parameter :: nkin   = 3 !-- 4 if full, 3 if use S(C1+C2)-S = 0

#if(_withchecks == 1)
  real(dp), public, save :: FintNLO_qcd(nFint)
  real(dp), public, save :: FintNLO_qcd_vs(3)
#endif

  private

  public :: xsect_nloqcd_s_ns_raoul
  public :: xsect_nloqcd_s_qg_raoul
  public :: xsect_nloqcd_s_gq_raoul

contains


  !--------------------------------------------------
  !-- subtractions
  !--------------------------------------------------


  function xsect_nloqcd_s_ns_raoul(yRnd,ff,vegasweight)
    use mod_subtrfn_nlo_z
    integer :: xsect_nloqcd_s_ns_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintNLO_s(3),kin(3,nintsub)
    real(dp) :: res_lo(2,2), res_1l(2,2)
    real(dp) :: respdf_1l(ipdf), respdf_qqb(ipdf)
    real(dp) :: z,xis(2),q2,intsub_qqb(-1:1,nintsub)

    xsect_nloqcd_s_ns_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    z=buff+onet*real(yRnd(kLO_max_full),dp)
    
#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif
        
    call open_histo()

    call kinematics_lo(xx,LOProc)
    
    LOProc%ids(1:4) = [0,0,id_el,-id_el]

    ! save the unboosted xi1 and xi2
    xis(1:2) = LOProc%PartFrac(1:2)
    q2 = two*scr(LOProc%AmpMom(1:4,1),LOProc%AmpMom(1:4,2))

    call cut_histo(LOProc)
    if (LOProc%makecut.or.LOProc%flag) then

       kin = zero
       FintNLO_s = zero

    else

       call res_tree_qqb(LOProc%AmpMom,res_lo)

       ! get the relevant int. subtr. functions
       intsub_qqb = sub_g_qqb_qqb_z(q2,LOProc%muf**2,z)           
       
       !-- [1,2]
       call get_respdf(ns_lumi,1,0,LOProc,res_1l,respdf_1l)
       call get_respdf(ns_lumi,1,0,LOProc,res_lo,respdf_qqb)

       ! -- int. subtraction term
       kin(1,:) = respdf_qqb(1)*LOProc%wgt * (intsub_qqb(-1,:) - intsub_qqb(1,:))*two     ! delta-pls, two legs


       !-- [z,2]
       
       ! change xi1 and recompute pdfs
       LOProc%PartFrac(1) = xis(1)/z
       LOProc%PartFrac(2) = xis(2)

       call get_respdf(ns_lumi,1,0,LOProc,res_lo,respdf_qqb)
       kin(2,:) = respdf_qqb(1)*LOProc%wgt/z * intsub_qqb(0,:) 
            
       !-- [1,z]

       ! change xi1 and recompute pdfs
       LOProc%PartFrac(1) = xis(1)
       LOProc%PartFrac(2) = xis(2)/z

       call get_respdf(ns_lumi,1,0,LOProc,res_lo,respdf_qqb)
       kin(3,:) = respdf_qqb(1)*LOProc%wgt/z * intsub_qqb(0,:) 

       call fill_histo(sum(kin(1:3,:),dim=1),vegasweight)


    endif

    ff(1) = sum(kin(1:3,1),dim=1)

!    if (ff(1) .ne. zero) then
!       print *, "xx",xx,z
!       print *, "kin", kin(1:3,1)
!       print *, "ff",ff(1)
!       stop
!    endif



    call close_histo()

!    call check_ff(ff,kLO_max_full+1,3,(/xx,z/),FintNLO_s)
    
#if(_withchecks == 1)
    FintNLO_qcd_vs = FintNLO_s
#endif
    
  end function xsect_nloqcd_s_ns_raoul





  function xsect_nloqcd_s_qg_raoul(yRnd,ff,vegasweight)
    use mod_subtrfn_nlo_z
    integer :: xsect_nloqcd_s_qg_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintNLO_s(3),kin(1,nintsub)
    real(dp) :: res_lo(2,2)
    real(dp) :: respdf_qg(ipdf)
    real(dp) :: z,xis(2),q2,intsub_qg(-1:1,nintsub)

    xsect_nloqcd_s_qg_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    z=buff+onet*real(yRnd(kLO_max_full),dp)
    
#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif
        
    call open_histo()

    call kinematics_lo(xx,LOProc)
    
 !   LOProc%ids(1:4) = [0,0,id_el,-id_el]

    ! save the unboosted xi1 and xi2
    xis(1:2) = LOProc%PartFrac(1:2)
    q2 = two*scr(LOProc%AmpMom(1:4,1),LOProc%AmpMom(1:4,2))

    call cut_histo(LOProc)
    if (LOProc%makecut.or.LOProc%flag) then

       kin = zero
       FintNLO_s = zero

    else

       call res_tree_qqb(LOProc%AmpMom,res_lo)

       intsub_qg = sub_g_gq_qqb_z(q2,LOProc%muf**2,z)

       !-- [1,z]

       ! change xi1 and recompute pdfs
       LOProc%PartFrac(1) = xis(1)
       LOProc%PartFrac(2) = xis(2)/z

       call get_respdf(qg_lumi,1,0,LOProc,res_lo,respdf_qg)           
       kin(1,:) = respdf_qg(1)*LOProc%wgt/z * intsub_qg(0,:)        ! reg+pls


       call fill_histo(kin(1,:),vegasweight)


    endif

    ff(1) = kin(1,1)

!    if (ff(1) .ne. zero) then
!       print *, "xx",xx,z
!       print *, "kin", kin(1,1)
!       print *, "ff",ff(1)
!       stop
!    endif



    call close_histo()

!    call check_ff(ff,kLO_max_full+1,3,(/xx,z/),FintNLO_s)
    
#if(_withchecks == 1)
    FintNLO_qcd_vs = FintNLO_s
#endif
    
  end function xsect_nloqcd_s_qg_raoul




  function xsect_nloqcd_s_gq_raoul(yRnd,ff,vegasweight)
    use mod_subtrfn_nlo_z
    integer :: xsect_nloqcd_s_gq_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintNLO_s(3),kin(1,nintsub)
    real(dp) :: res_lo(2,2)
    real(dp) :: respdf_gq(ipdf)
    real(dp) :: z,xis(2),q2,intsub_gq(-1:1,nintsub)

    xsect_nloqcd_s_gq_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    z=buff+onet*real(yRnd(kLO_max_full),dp)
    
#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif
        
    call open_histo()

    call kinematics_lo(xx,LOProc)
    
 !   LOProc%ids(1:4) = [0,0,id_el,-id_el]

    ! save the unboosted xi1 and xi2
    xis(1:2) = LOProc%PartFrac(1:2)
    q2 = two*scr(LOProc%AmpMom(1:4,1),LOProc%AmpMom(1:4,2))

    call cut_histo(LOProc)
    if (LOProc%makecut.or.LOProc%flag) then

       kin = zero
       FintNLO_s = zero

    else

       call res_tree_qqb(LOProc%AmpMom,res_lo)

       intsub_gq = sub_g_gq_qqb_z(q2,LOProc%muf**2,z)


       !-- [z,2]
       
       ! change xi1 and recompute pdfs
       LOProc%PartFrac(1) = xis(1)/z
       LOProc%PartFrac(2) = xis(2)

       call get_respdf(gq_lumi,1,0,LOProc,res_lo,respdf_gq)           
       kin(1,:) = respdf_gq(1)*LOProc%wgt/z * intsub_gq(0,:)        ! reg+pls
            

       call fill_histo(kin(1,:),vegasweight)


    endif

    ff(1) = kin(1,1)

!    if (ff(1) .ne. zero) then
!       print *, "xx",xx,z
!       print *, "kin", kin(1:3,1)
!       print *, "ff",ff(1)
!       stop
!    endif



    call close_histo()

!    call check_ff(ff,kLO_max_full+1,3,(/xx,z/),FintNLO_s)
    
#if(_withchecks == 1)
    FintNLO_qcd_vs = FintNLO_s
#endif
    
  end function xsect_nloqcd_s_gq_raoul

  
  
end module mod_xsects_nloqcd_s_raoul
  

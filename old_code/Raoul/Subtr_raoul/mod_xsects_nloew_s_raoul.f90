module mod_xsects_nloew_s_raoul
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

  private

  public :: xsect_nloew_s_ns_raoul
  public :: xsect_nloew_s_qa_raoul
  public :: xsect_nloew_s_aq_raoul

contains


  !--------------------------------------------------
  !-- subtractions
  !--------------------------------------------------


  function xsect_nloew_s_ns_raoul(yRnd,ff,vegasweight)
    use mod_subtrfn_nlo_z
    integer :: xsect_nloew_s_ns_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintNLO_s(3),kin(3,nintsub)
    real(dp) :: res_lo(2,2), res_tmp(2,2)
    real(dp) :: respdf_qqb(ipdf)
    real(dp) :: z,xis(2),q2,intsub_u(-1:1,nintsub),intsub_d(-1:1,nintsub),intsub_ub(-1:1,nintsub),intsub_db(-1:1,nintsub)
    real(dp) :: eta31, eta41, eta32, eta42

    xsect_nloew_s_ns_raoul = 0

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

    ! values of etaij for int.  subtraction functions
    eta31 =  scr(LOProc%AmpMom(1:4,1),LOProc%AmpMom(1:4,3))/two/LOProc%AmpMom(1,1)/LOProc%AmpMom(1,3)
    eta41 =  scr(LOProc%AmpMom(1:4,1),LOProc%AmpMom(1:4,4))/two/LOProc%AmpMom(1,1)/LOProc%AmpMom(1,4)
    eta32 =  scr(LOProc%AmpMom(1:4,2),LOProc%AmpMom(1:4,3))/two/LOProc%AmpMom(1,2)/LOProc%AmpMom(1,3)
    eta42 =  scr(LOProc%AmpMom(1:4,2),LOProc%AmpMom(1:4,4))/two/LOProc%AmpMom(1,2)/LOProc%AmpMom(1,4)


    
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


       ! get the relevant int. subtr. functions for different quark charges
       intsub_u = sub_a_qqb_qqb_z(q2,LOProc%muf**2,z,eta31,eta32,eta41,eta42,Qup,Q_lep)           
       intsub_d = sub_a_qqb_qqb_z(q2,LOProc%muf**2,z,eta31,eta32,eta41,eta42,Qdn,Q_lep)           
       intsub_ub = sub_a_qqb_qqb_z(q2,LOProc%muf**2,z,eta31,eta32,eta41,eta42,-Qup,Q_lep)           
       intsub_db = sub_a_qqb_qqb_z(q2,LOProc%muf**2,z,eta31,eta32,eta41,eta42,-Qdn,Q_lep)           


       !-- [1,2]
       res_tmp(1,1)  = res_lo(1,1)*(intsub_d(-1,1)  - intsub_d(+1,1))*two ! delta-pls, two legs
       res_tmp(1,2)  = res_lo(1,2)*(intsub_u(-1,1)  - intsub_u(+1,1))*two ! delta-pls, two legs

       res_tmp(2,1)  = res_lo(2,1)*(intsub_db(-1,1)  - intsub_db(+1,1))*two ! delta-pls, two legs
       res_tmp(2,2)  = res_lo(2,2)*(intsub_ub(-1,1)  - intsub_ub(+1,1))*two ! delta-pls, two legs
       
       call get_respdf(ns_lumi,0,1,LOProc,res_tmp,respdf_qqb)

       kin(1,:) = respdf_qqb(1)*LOProc%wgt


       !-- [z,2]       
       ! change xi1 and recompute pdfs
       LOProc%PartFrac(1) = xis(1)/z
       LOProc%PartFrac(2) = xis(2)

       res_tmp(1,1)  = res_lo(1,1)*intsub_d(0,1)      ! reg
       res_tmp(1,2)  = res_lo(1,2)*intsub_u(0,1)      ! reg

       res_tmp(2,1)  = res_lo(2,1)*intsub_db(0,1)      ! reg
       res_tmp(2,2)  = res_lo(2,2)*intsub_ub(0,1)      ! reg

       call get_respdf(ns_lumi,0,1,LOProc,res_tmp,respdf_qqb)
       kin(2,:) = respdf_qqb(1)*LOProc%wgt/z 
            
       !-- [1,z]

       ! change xi1 and recompute pdfs
       LOProc%PartFrac(1) = xis(1)
       LOProc%PartFrac(2) = xis(2)/z

       res_tmp(1,1)  = res_lo(1,1)*intsub_d(0,1)      ! reg
       res_tmp(1,2)  = res_lo(1,2)*intsub_u(0,1)      ! reg

       res_tmp(2,1)  = res_lo(2,1)*intsub_db(0,1)      ! reg
       res_tmp(2,2)  = res_lo(2,2)*intsub_ub(0,1)      ! reg

       call get_respdf(ns_lumi,0,1,LOProc,res_tmp,respdf_qqb)
       kin(3,:) = respdf_qqb(1)*LOProc%wgt/z 


       call fill_histo(sum(kin(1:3,:),dim=1),vegasweight)


    endif

    ff(1) = sum(kin(1:3,1),dim=1)

!    if (ff(1) .ne. zero) then
!       print *, "xx",xx,z
!       print *, "kin", kin(1:3)
!       print *, "ff",ff(1)
!       stop
!    endif



    call close_histo()

!    call check_ff(ff,kLO_max_full+1,3,(/xx,z/),FintNLO_s)
    
    
  end function xsect_nloew_s_ns_raoul





  function xsect_nloew_s_qa_raoul(yRnd,ff,vegasweight)
    use mod_subtrfn_nlo_z
    integer :: xsect_nloew_s_qa_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintNLO_s(3),kin(1,nintsub)
    real(dp) :: res_lo(2,2)
    real(dp) :: respdf_qa(ipdf)
    real(dp) :: z,xis(2),q2,intsub_qa(-1:1,nintsub)


    xsect_nloew_s_qa_raoul = 0

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

       intsub_qa = sub_a_aq_qqb_z(q2,LOProc%muf**2,z)

       !-- [1,z]

       ! change xi1 and recompute pdfs
       LOProc%PartFrac(1) = xis(1)
       LOProc%PartFrac(2) = xis(2)/z

       call get_respdf(qa_lumi,0,1,LOProc,res_lo,respdf_qa)           
       kin(1,:) = respdf_qa(1)*LOProc%wgt/z * intsub_qa(0,:)        ! reg+pls


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
    
    
  end function xsect_nloew_s_qa_raoul




  function xsect_nloew_s_aq_raoul(yRnd,ff,vegasweight)
    use mod_subtrfn_nlo_z
    integer :: xsect_nloew_s_aq_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintNLO_s(3),kin(2,nintsub)
    real(dp) :: res_lo(2,2), res_aa,res_tmp(2,2)
    real(dp) :: respdf_aq(ipdf),respdf_aa(ipdf)
    real(dp) :: z,xis(2),q2,intsub_qqb(-1:1,nintsub),intsub_aa(-1:1,nintsub)

    xsect_nloew_s_aq_raoul = 0

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
       call res_treeAA_aa(LOProc%AmpMom,res_aa)

       intsub_qqb = sub_a_aq_qqb_z(q2,LOProc%muf**2,z)
       intsub_aa = sub_a_aq_aa_z(q2,LOProc%muf**2,z)
!       intsub_aa_dn = sub_a_aq_aa_z(q2,LOProc%muf**2,z,Qdn)


       !-- [z,2]
       
       ! change xi1 and recompute pdfs
       LOProc%PartFrac(1) = xis(1)/z
       LOProc%PartFrac(2) = xis(2)

       res_lo(:,1) = Qdn**2*res_lo(:,1)
       res_lo(:,2) = Qup**2*res_lo(:,2)


       call get_respdf(aq_lumi,0,1,LOProc,res_lo,respdf_aq)           
       kin(1,:) = respdf_aq(1)*LOProc%wgt/z * intsub_qqb(0,:)

       !-- [1,z]
       ! change xi1 and recompute pdfs
       LOProc%PartFrac(1) = xis(1)
       LOProc%PartFrac(2) = xis(2)/z

       res_tmp(:,1) = res_aa*Qdn**2                                    ! include charges sq
       res_tmp(:,2) = res_aa*Qup**2
       call get_respdf(aq_lumi,0,1,LOProc,res_tmp,respdf_aa)           

       kin(2,:) = respdf_aa(1)*LOProc%wgt/z * intsub_aa(0,:)

       call fill_histo(sum(kin(1:2,:),dim=1),vegasweight)


    endif

    ff(1) = sum(kin(1:2,1),dim=1)

!    if (ff(1) .ne. zero) then
!       print *, "xx",xx,z
!       print *, "kin", kin(1:3,1)
!       print *, "ff",ff(1)
!       stop
!    endif



    call close_histo()

!    call check_ff(ff,kLO_max_full+1,3,(/xx,z/),FintNLO_s)
    
  end function xsect_nloew_s_aq_raoul

  
  
end module mod_xsects_nloew_s_raoul
  

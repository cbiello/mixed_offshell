module mod_xsects_nnlo_s_ag_raoul
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_kinematics_lo
  use mod_kinematics_lo_z_raoul
  use mod_kinematics_nlo_z_raoul
!  use mod_kinematics_nlo
  use mod_process
  use mod_auxfunctions
  use mod_aux_sectors
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_amplitudes_loop_ppll
  use mod_eikonals
  use mod_partitions
  use mod_splittings_bare
  use mod_subtrfn_nnlo_z_ag

  implicit none
  integer, parameter :: nFint  = 18 !-- 6 if full, 3 if use S(C1+C2)-S = 0
  integer, parameter :: nkin   = 18 !-- 4 if full, 3 if use S(C1+C2)-S = 0
  logical, parameter :: fncheck = .false.     !-- print out values for int. sub. functions and their arguments as a check

  real(dp), parameter :: charges_ns(4,4) = reshape(& 
       [-Qdn,   Qdn,  -Qup,   Qup,   &
         Qdn,  -Qdn,   Qup,  -Qup,   &
         Q_lep, Q_lep, Q_lep, Q_lep, &
        -Q_lep,-Q_lep,-Q_lep,-Q_lep]  &
        , [4,4])


#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_onloqcd_ag(nFint)
  real(dp), public, save :: FintNNLO_onloewk_ag(nFint)
#endif

  private

  public :: xsect_nnlo_sub12_ag_raoul
  public :: xsect_nnlo_onloewk_ag_raoul
  public :: xsect_nnlo_onloqcd_ag_raoul


contains


  !--------------------------------------------------
  !-- subtraction counterterms
  !--------------------------------------------------

 function xsect_nnlo_onloewk_ag_raoul(yRnd,ff,vegasweight)
   ! The coefficients of ONLO^gamma FLM[1_gamma,z2_q,3,4,|5_q]  
    integer :: xsect_nnlo_onloewk_ag_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc_1z2,C1Lim_1z2,C2Lim_1z2
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: kin(3)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(2,2),res_lo(2,2),res_tmp(2,2),res_tmp_aq(1,1),res_loAA
    real(dp)    :: z,z5i,s5i,eta52,Ec,eta52lim
    real(dp)    :: subint_1z24_hard,subint_1z24_c1lim,subint_1z24_c2lim

    xsect_nnlo_onloewk_ag_raoul = 0

    ff(1) = zero
    kin = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z=buff+onet*real(yRnd(kNLO_max_full),dp)

 
#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if ( (xx(xRHO).lt.buff_r) .or. ((one-xx(xRHO)).lt.buff_r) .or. ((one-z) .lt. buff_z) ) then
       failed_points = failed_points + 1
       return
    endif

!    xx(xrho) = 1E-10_dp
!    z = one - 1E-8_dp



    call open_histo()

    !!!!------------------------!!!!!
    !!!! FLM[1_a,z.2_q,3,4|5_q] !!!!!
    !!!!------------------------!!!!!

    
    call kinematics_nlo_z_is(yr=xx,z1=one,z2=z,HardProc=HardProc_1z2,&
         C1Lim=C1Lim_1z2,C2Lim=C2Lim_1z2,compute_etas=.true.)
    HardProc_1z2%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim_1z2%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim_1z2%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard                                                                                                                                                                                                 
    call cut_histo(HardProc_1z2)
    if (HardProc_1z2%makecut.or.HardProc_1z2%flag) then

       kin(1) = zero
       FintNNLO_onloewk_ag(1) = zero
    else

       call res_tree_a_aq(HardProc_1z2%AmpMom,res_nlo)

       res_tmp_aq(1,1) = ndn*(res_nlo(1,1)+res_nlo(2,1))+nup*(res_nlo(1,2)+res_nlo(2,2)) !  sum matrix elements weighted by number of flavors

       call get_respdf(ag_lumi,1,1,HardProc_1z2,res_tmp_aq,respdf)

       Ec = HardProc_1z2%Lim_Ei(1)                      ! Ec = E1 = E2

       eta52 = HardProc_1z2%Lim_etaij(2,5)
       eta52lim = eta52

       subint_1z24_hard = ag_qqb_sub_aqq_zi(Ec,HardProc_1z2%muf(1),z,eta52,eta52lim)       
       if (fncheck) then
          print *, "EC->",EC,",mu->",HardProc_1z2%muf(1),",eta52-",eta52,",z->",z
          print *, "1z24", subint_1z24_hard
          pause
       endif

       respdf = respdf*HardProc_1z2%wgt * subint_1z24_hard


       kin(1) = respdf(1)
       FintNNLO_onloewk_ag(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C1                                                                                                                                                                                                   
    call cut_histo(C1Lim_1z2)

    if (C1Lim_1z2%makecut.or.C1Lim_1z2%flag) then

       kin(2) = zero
       FintNNLO_onloewk_ag(2) = zero

    else

       call res_tree_qqb(C1Lim_1z2%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       res_tmp_aq(1,1) = ndn*(res_lo(1,1)+res_lo(2,1))+nup*(res_lo(1,2)+res_lo(2,2)) !  sum matrix elements weighted by number of flavors

       call get_respdf(ag_lumi,1,1,C1Lim_1z2,res_tmp_aq,respdf)


       z5i   = C1Lim_1z2%Lim_KinInv(1)
       s5i = C1Lim_1z2%Lim_KinInv(2)

       Ec = C1Lim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2

       eta52 = one             
       eta52lim = eta52

       subint_1z24_c1lim = ag_qqb_sub_aqq_zi(Ec,C1Lim_1z2%muf(1),z,eta52,eta52lim)       


       respdf = (-one)*respdf*(two/s5i)*(xn*Pgq_spav(z5i)/(one-z5i))&
            * C1Lim_1z2%wgt * subint_1z24_c1lim

       FintNNLO_onloewk_ag(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C2                                                                                                                                                                                                   
    call cut_histo(C2Lim_1z2)

    if (C2Lim_1z2%makecut.or.C2Lim_1z2%flag) then

       kin(3) = zero
       FintNNLO_onloewk_ag(3) = zero

    else

       z5i = C2Lim_1z2%Lim_KinInv(1)
       s5i = C2Lim_1z2%Lim_KinInv(2)

       Ec = C2Lim_1z2%Lim_Ei(1)                      ! Ec = E1 = E2
       eta52 = zero
       eta52lim = C2Lim_1z2%Lim_etaij(2,5)

       subint_1z24_c2lim = ag_qqb_sub_aqq_zi(Ec,C2Lim_1z2%muf(1),z,eta52,eta52lim)       


      !-- here we use spin-averages since there are no spin correlations for the aa -> e-e+ amplitude
       call res_treeAA_aa(C2Lim_1z2%AmpMom,res_loAA)
       !-- 1=xn*aveqa/aveaa                                                                                                      
       res_loAA = res_loAA * Q_lep2 * Pqq(z5i)/(one-z5i) !-- use Pqq/(1-z) because of definition of z  

       res_tmp(:,1) = Qdn2 * res_loAA
       res_tmp(:,2) = Qup2 * res_loAA

       res_tmp_aq(1,1) = ndn*(res_tmp(1,1)+res_tmp(2,1))+nup*(res_tmp(1,2)+res_tmp(2,2)) !  sum matrix elements weighted by number of flavors

       call get_respdf(ag_lumi,1,1,C2Lim_1z2,res_tmp_aq,respdf)

       respdf = (-one)*respdf*(two/s5i)&
            * C2Lim_1z2%wgt * subint_1z24_c2lim

       FintNNLO_onloewk_ag(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif


    


    ff(1) = sum(kin)

    if (ff(1) .ne. ff(1)) then
       print *, "NaN"
       print *, "xx",xx
       print *, "z",z
       print *, "kin",kin
       stop
    endif
    call close_histo()

!    if (kin(1)*kin(2)*kin(3) .ne. zero) then
!
!       print *, "kin in onlo z12", kin(1:3)
!       print *, "ff",ff(1),ff(1)/kin(1)
!       pause
!   endif

    call check_ff(ff,xx,FintNNLO_onloewk_ag)

#if(_withchecks == 1)
!!    FintNLO_ew(1:3) = FintNNLO_onloewk_ag
#endif

  end function xsect_nnlo_onloewk_ag_raoul




 function xsect_nnlo_onloqcd_ag_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_onloqcd_ag_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc_z12,C2Lim_z12
    !--                                      
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: kin(2)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(2,2),res_lo(2,2),res_tmp_qg(1,1)
    real(dp)    :: z,z5i,s5i,Ec
    real(dp)    :: subint_z124_hard_u,subint_z124_c1lim_u,subint_z124_hard_d,subint_z124_c1lim_d

    xsect_nnlo_onloqcd_ag_raoul = 0

    ff(1) = zero
    kin = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z=buff+onet*real(yRnd(kNLO_max_full),dp)


#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if ( ((one-xx(xRHO)).lt.buff_r) .or. ((one-z) .lt. buff_z) )  then
       failed_points = failed_points + 1
       return
    endif

!    xx(xrho) = one-1E-10_dp
!    z = one - 1E-8_dp

    call open_histo()

    !!!!------------------------!!!!!                                                                                              
    !!!! FLM[z.1_q,2_g,3,4|5_q] !!!!!                                                                                              
    !!!!------------------------!!!!!

    call kinematics_nlo_z_is(yr=xx,z1=z,z2=one,HardProc=HardProc_z12,&
         C2Lim=C2Lim_z12,compute_etas=.true.)
    HardProc_z12%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C2Lim_z12%ids(1:4) = [0,0,id_el,-id_el]

!-- Hard
    call cut_histo(HardProc_z12)
    if (HardProc_z12%makecut.or.HardProc_z12%flag) then

       kin(1) = zero
       FintNNLO_onloqcd_ag(1) = zero

    else

       call res_tree_g_qg(HardProc_z12%AmpMom,res_nlo)

       Ec = HardProc_z12%Lim_Ei(1)                      ! Ec = E1 = E2                                                              
       subint_z124_hard_u = ag_qqb_sub_qgq_zi(Ec,HardProc_z12%muf(1),z,Qup)
       subint_z124_hard_d = ag_qqb_sub_qgq_zi(Ec,HardProc_z12%muf(1),z,Qdn)
       if (fncheck) then
          print *, "EC->",EC,",mu->",HardProc_z12%muf(1),",z->",z
          print *, "z124 u", subint_z124_hard_u
          print *, "z124 d", subint_z124_hard_d
          pause
       endif


       res_nlo(:,1) = res_nlo(:,1)*subint_z124_hard_d
       res_nlo(:,2) = res_nlo(:,2)*subint_z124_hard_u

       res_tmp_qg(1,1) = ndn*(res_nlo(1,1)+res_nlo(2,1))+nup*(res_nlo(1,2)+res_nlo(2,2)) !  sum matrix elements weighted by number of flavors

       call get_respdf(ag_lumi,1,1,HardProc_z12,res_tmp_qg,respdf)

       respdf = respdf*HardProc_z12%wgt

       kin(1) = respdf(1)
       FintNNLO_onloqcd_ag(1) = kin(1)

       call fill_histo(respdf,vegasweight)
endif

   !-- C2
    call cut_histo(C2Lim_z12)

    if (C2Lim_z12%makecut.or.C2Lim_z12%flag) then

       kin(2) = zero
       FintNNLO_onloqcd_ag(2) = zero

    else

       call res_tree_qqb(C2Lim_z12%AmpMom,res_lo)

       z5i   = C2Lim_z12%Lim_KinInv(1)
       s5i = C2Lim_z12%Lim_KinInv(2)

       Ec = C2Lim_z12%Lim_Ei(1)                      ! Ec = E1 = E2

       subint_z124_c1lim_u = ag_qqb_sub_qgq_zi(Ec,C2Lim_z12%muf(1),z,Qup)
       subint_z124_c1lim_d = ag_qqb_sub_qgq_zi(Ec,C2Lim_z12%muf(1),z,Qdn)

       res_lo(:,1) = res_lo(:,1)*subint_z124_c1lim_d
       res_lo(:,2) = res_lo(:,2)*subint_z124_c1lim_u

       res_tmp_qg(1,1) = ndn*(res_lo(1,1)+res_lo(2,1))+nup*(res_lo(1,2)+res_lo(2,2)) !  sum matrix elements weighted by number of flavors

       call get_respdf(ag_lumi,1,1,C2Lim_z12,res_tmp_qg,respdf)

       !-- use Pgq so that in the soft limit we do not have 1/[1-(1-x)]                                                            
       !-- Tr = Cf * aveqg/aveqq                                       
                                                            
       respdf = (-one)*respdf*(two/s5i)*(tr*Pgq_spav(z5i)/(one-z5i)) &
            * C2Lim_z12%wgt

       FintNNLO_onloqcd_ag(2) = respdf(1)
       kin(2) = respdf(1)
       call fill_histo(respdf,vegasweight)

    endif


    ff(1) = sum(kin)
    call close_histo()

!!    if (kin(1)*kin(2) .ne. zero) then
!!       print *, "kin in onlo z12", kin(1:2)
!!       print *, "ff",ff(1)
!!       pause
!!
!!    endif

  end function xsect_nnlo_onloqcd_ag_raoul

  function xsect_nnlo_sub12_ag_raoul(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_sub12_ag_raoul
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: FintLO_s(2),kin(2)
    real(dp) :: res_lo(2,2), res_tmp(2,2),respdf(ipdf), res_aa
    real(dp) :: Ec,mu,z,zbar,xis(1:2)
    real(dp) :: intsub_qqb_u, intsub_qqb_d,intsub_aa_u, intsub_aa_d

    xsect_nnlo_sub12_ag_raoul = 0

    ff(1) = zero
    kin  = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    z   =buff+onet*real(yRnd(kLO_max_full),dp)
    zbar=buff+onet*real(yRnd(kLO_max_full+1),dp)

    if ( ((one-z) .lt. buff_z) .or. ((one-zbar) .lt. buff_z) ) then
       failed_points = failed_points + 1
       return
    endif


#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()

    ! no eta-dependence in the int. subtr. functions
    ! so we can use the usual trick of putting the boosts in the pdfs
    call kinematics_lo(xx,LOProc)

    Ec = LOProc%AmpMom(1,1)

    xis(1:2) = LOProc%PartFrac(1:2)
    LOProc%ids(1:4) = [0,0,id_el,-id_el]

    call cut_histo(LOProc)
    mu = LOProc%muf(1)

    if (LOProc%makecut.or.LOProc%flag) then

       kin = zero
       FintLO_s = zero

    else

       call res_tree_qqb(LOProc%AmpMom,res_lo)
       call res_treeAA_aa(LOProc%AmpMom,res_aa)


       intsub_qqb_u = ag_qqb_sub_qqbz1zbar2(mu,z,zbar,Ec,Qup)
       intsub_qqb_d = ag_qqb_sub_qqbz1zbar2(mu,z,zbar,Ec,Qdn)

       intsub_aa_u = ag_qqb_sub_aa1z2(mu,zbar,Ec,Qup)    
       intsub_aa_d = ag_qqb_sub_aa1z2(mu,zbar,Ec,Qdn)

       if (fncheck) then
          print *, "EC->",EC,",mu->",mu,",z->",z,",zb->",zbar
          print *, "z1zb2 u", intsub_qqb_u
          print *, "z1zb2 d", intsub_qqb_d

          print *, "1zb2 u", intsub_aa_u
          print *, "1zb2 d", intsub_aa_d
          pause
       endif

!!       print *,"mu->",mu,",EC->",Ec,",z->",z,",zb->",zbar
!!
!!       print *, "zzb u",intsub_qqb_u
!!       print *, "zzb d",intsub_qqb_d
!!       print *, "aa u",intsub_aa_u
!!       print *, "aa d",intsub_aa_d
!!       pause


       ! FLM[z1_q,z2_qb]
       res_tmp(1,1) = (res_lo(1,1)+res_lo(2,1))*ndn*intsub_qqb_d + (res_lo(1,2)+res_lo(2,2))*nup*intsub_qqb_u

       LOProc%PartFrac(1) = xis(1)/z
       LOProc%PartFrac(2) = xis(2)/zbar

       call get_respdf(ag_lumi,1,1,LOProc,res_tmp,respdf)
       respdf = respdf*LOProc%wgt/z/zbar

       kin(1) = respdf(1)
       call fill_histo(respdf,vegasweight)


       ! FLM[1_a,zbar 2_a]
       res_tmp = zero
       res_tmp(1,1) = res_aa * (intsub_aa_d *ndn +  intsub_aa_u *nup)
       

       LOProc%PartFrac(1) = xis(1)
       LOProc%PartFrac(2) = xis(2)/zbar

       call get_respdf(ag_lumi,1,1,LOProc,res_tmp,respdf)
       respdf = respdf*LOProc%wgt/zbar

       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin(1:2))

    call close_histo()


  end function xsect_nnlo_sub12_ag_raoul







  
end module mod_xsects_nnlo_s_ag_raoul
  

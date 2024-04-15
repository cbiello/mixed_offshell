module mod_xsects_nnlo_subonloewk_aq_fc
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_kinematics_nlo_z_fc
  use mod_process
  use mod_auxfunctions
  use mod_aux_sectors_fc
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_partitions
  use mod_splittings_bare
  use mod_eikonals
  use mod_intsub_fc
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_subonloewk_aq(6)
#endif
  private

  public :: xsect_nnlo_subonloewk_is_aq_fc

contains

  !-------------------------------------------------------------
  !-- aq channel below
  !-------------------------------------------------------------
  
  function xsect_nnlo_subonloewk_is_aq_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_subonloewk_is_aq_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(6),kin(6)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(2,2),res_lo(2,2),res_tmp(2,2),res_loAA
    real(dp)    :: zz,s5i
    real(dp)    :: z,Lmu2(ipdf),mu2(ipdf),w22,w56_11,w56_22,int_sub(ipdf)

    xsect_nnlo_subonloewk_is_aq_fc = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z = buff+onet*real(yRnd(kNLO_max+1),dp)

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       z = yRnd(kNLO_max+1)
       print *, 'overriding input'
    endif
#endif

    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_r)) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    !-- FLM[1,zb,3,4,5]
    
    call kinematics_nlo_is_z(yr=xx,z1=one,z2=z,HardProc=HardProc,&
         C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.false.)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_ns(1) = zero

    else    

       call get_Lmu2(HardProc,Lmu2)
       
       call res_tree_a_aq(HardProc%AmpMom,res_nlo)
       call get_respdf(aq_lumi,1,1,HardProc,res_nlo,respdf)

       call get_w22(HardProc,w22)
       int_sub = intsub_aq_onloewk_1zb(z,Lmu2,w22,HardProc%Lim_KinInv,0)
       respdf  = respdf*int_sub

       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintNLO_ns(1) = kin(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 1',FintNLO_ns(1)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintNLO_ns(2) = zero

    else

       call get_Lmu2(C1Lim,Lmu2)

       call res_tree_qqb(C1Lim%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       call get_respdf(aq_lumi,1,1,C1Lim,res_lo,respdf)

       w22 = zero
       int_sub = intsub_aq_onloewk_1zb(z,Lmu2,w22,C1Lim%Lim_KinInv,0)
       respdf  = respdf*int_sub
       
       zz   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       !-- xn = aveqa/aveqq
       respdf = (-one)*respdf*(two/s5i)*(xn*Pgq_spav(zz)/(one-zz))&
            * C1Lim%wgt

       FintNLO_ns(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 2',FintNLO_ns(2)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(3) = zero
       FintNLO_ns(3) = zero

    else

       call get_Lmu2(C2Lim,Lmu2)
       
       zz   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- here we use spin-averages since there are no spin correlations for the aa -> e-e+ amplitude
       call res_treeAA_aa(C2Lim%AmpMom,res_loAA)
       !-- 1=xn*aveqa/aveaa
       res_loAA = res_loAA * Pqq(zz)/(one-zz) !-- use Pqq/(1-z) because of definition of z

       res_tmp(:,1) = Qdn2 * res_loAA
       res_tmp(:,2) = Qup2 * res_loAA
       
       call get_respdf(aq_lumi,1,1,C2Lim,res_tmp,respdf)

       w22 = one
       int_sub = intsub_aq_onloewk_1zb(z,Lmu2,w22,C2Lim%Lim_KinInv,0)
       respdf  = respdf*int_sub
       
       respdf = (-one)*respdf*(two/s5i)&
            * C2Lim%wgt

       FintNLO_ns(3) = respdf(1)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 3',FintNLO_ns(3)
       !print *, '[]'
       
    endif

    !---

    !-- FLM[1,2,3,4,5]
    
    call kinematics_nlo_is_z(yr=xx,z1=one,z2=one,HardProc=HardProc,&
         C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.false.)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(4) = zero
       FintNLO_ns(4) = zero

    else    

       call get_Lmu2(HardProc,Lmu2,mu2)
       
       call res_tree_a_aq(HardProc%AmpMom,res_nlo)
       call get_respdf(aq_lumi,1,1,HardProc,res_nlo,respdf)

       call get_w22(HardProc,w22)
       call get_w56_ii(HardProc,w56_11,w56_22)
       int_sub = intsub_aq_onloewk_1zb(z,Lmu2,w22,HardProc%Lim_KinInv,1) + &
            intsub_aq_onloewk_12(mu2,w56_11,w56_22,HardProc)
       respdf  = respdf*int_sub

       respdf = respdf*HardProc%wgt

       kin(4) = respdf(1)
       FintNLO_ns(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 4',FintNLO_ns(4)
       
    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(5) = zero
       FintNLO_ns(5) = zero

    else

       call get_Lmu2(C1Lim,Lmu2,mu2)

       call res_tree_qqb(C1Lim%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)

       call get_respdf(aq_lumi,1,1,C1Lim,res_lo,respdf)

       w22 = zero
       w56_11 = one
       w56_22 = zero
       int_sub = intsub_aq_onloewk_1zb(z,Lmu2,w22,C1Lim%Lim_KinInv,1) + &
            intsub_aq_onloewk_12(mu2,w56_11,w56_22,C1Lim)

       respdf  = respdf*int_sub
       
       zz   = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       !-- xn = aveqa/aveqq
       respdf = (-one)*respdf*(two/s5i)*(xn*Pgq_spav(zz)/(one-zz))&
            * C1Lim%wgt

       FintNLO_ns(5) = respdf(1)
       kin(5) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 5',FintNLO_ns(5)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(6) = zero
       FintNLO_ns(6) = zero

    else

       call get_Lmu2(C2Lim,Lmu2,mu2)

       zz   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- here we use spin-averages since there are no spin correlations for the aa -> e-e+ amplitude
       call res_treeAA_aa(C2Lim%AmpMom,res_loAA)
       !-- 1=xn*aveqa/aveaa
       res_loAA = res_loAA * Pqq(zz)/(one-zz) !-- use Pqq/(1-z) because of definition of z

       res_tmp(:,1) = Qdn2 * res_loAA
       res_tmp(:,2) = Qup2 * res_loAA
       
       call get_respdf(aq_lumi,1,1,C2Lim,res_tmp,respdf)

       w22 = one
       w56_11 = zero
       w56_22 = one
       int_sub = intsub_aq_onloewk_1zb(z,Lmu2,w22,C2Lim%Lim_KinInv,1) + &
            intsub_aq_onloewk_12(mu2,w56_11,w56_22,C2Lim)
       respdf  = respdf*int_sub
       
       respdf = (-one)*respdf*(two/s5i)&
            * C2Lim%wgt

       FintNLO_ns(6) = respdf(1)
       kin(6) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 6',FintNLO_ns(6)

    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_ns)
    
#if(_withchecks == 1)
    FintNNLO_subonloewk_aq = FintNLO_ns
#endif
    
  end function xsect_nnlo_subonloewk_is_aq_fc

  !-- TC damping factor with 6||2
  subroutine get_w22(proc,damp)
    type(KinConfig), intent(in) :: proc
    real(dp), intent(out) :: damp
    real(dp) :: eta15,eta25

    eta15 = proc%Lim_etaij(1,5)
    eta25 = proc%Lim_etaij(2,5)

    damp = eta15*(one+eta25)

  end subroutine get_w22

  subroutine get_w56_ii(proc,w56_11,w56_22)
    type(KinConfig), intent(in) :: proc
    real(dp), intent(out) :: w56_11,w56_22
    real(dp) :: eta15,eta25

    eta15 = proc%Lim_etaij(1,5)
    eta25 = proc%Lim_etaij(2,5)

    w56_11 = eta25**2*(1 + 2*eta15)
    w56_22 = eta15**2*(1 + 2*eta25)

  end subroutine get_w56_ii
  
end module mod_xsects_nnlo_subonloewk_aq_fc
  

module mod_xsects_nnlo_subonloqcd_gq_fc
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_kinematics_nlo_z_fc
  use mod_process
  use mod_aux_sectors_fc
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_splittings_bare
  use mod_intsub_fc
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_subonloqcd_gq(4)
#endif

  private

  public :: xsect_nnlo_subonloqcd_gq_fc

contains

  function xsect_nnlo_subonloqcd_gq_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_subonloqcd_gq_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_gq(4),kin(4)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(2,2),res_lo(2,2)
    real(dp)    :: zz,s5i
    real(dp)    :: z,Lmu2(ipdf),my_mu2(ipdf),int_sub(2,2,ipdf)

    xsect_nnlo_subonloqcd_gq_fc = 0

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
    
    if (xx(xE)*xx(xRHO).lt.buff_r) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    !---
    
    !-- FLM[1,zb,3,4,5]
    
    call kinematics_nlo_is_z(yr=xx,z1=one,z2=z,HardProc=HardProc,C1Lim=C1Lim,compute_etas=.true.)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_gq(1) = zero

    else    

       call get_Lmu2(HardProc,Lmu2)
       
       call res_tree_g_gq(HardProc%AmpMom,res_nlo)
       int_sub = intsub_gq_onloqcd_1zb(z,Lmu2,0)

       call get_respdf(gq_lumi,1,1,HardProc,res_nlo,respdf,int_sub)

       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintNLO_gq(1) = kin(1)

       call fill_histo(respdf,vegasweight)

       !print *, 'FLM[1,zb] 1',FintNLO_gq(1)
       
    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintNLO_gq(2) = zero

    else

       call get_Lmu2(C1Lim,Lmu2)
       
       call res_tree_qqb(C1Lim%AmpMom,res_lo)
       int_sub = intsub_gq_onloqcd_1zb(z,Lmu2,0)
       
       call get_respdf(gq_lumi,1,1,C1Lim,res_lo,respdf,int_sub)

       zz  = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pgq so that in the soft limit we do not have 1/[1-(1-x)]
       !-- Tr = Cf * aveqg/aveqq
       respdf = (-one)*respdf*(two/s5i)*(tr*Pgq_spav(zz)/(one-zz)) & 
            * C1Lim%wgt

       FintNLO_gq(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, 'FLM[1,zb] 2',FintNLO_gq(2)

    endif

    !---
    
    !-- FLM[1,2,3,4,5]
    
    call kinematics_nlo_is_z(yr=xx,z1=one,z2=one,HardProc=HardProc,C1Lim=C1Lim,compute_etas=.true.)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(3) = zero
       FintNLO_gq(3) = zero

    else    

       call get_Lmu2(HardProc,Lmu2,my_mu2)
       
       call res_tree_g_gq(HardProc%AmpMom,res_nlo)
       int_sub = intsub_gq_onloqcd_1zb(z,Lmu2,1) + &
            intsub_gq_onloqcd_12(my_mu2,HardProc)

       call get_respdf(gq_lumi,1,1,HardProc,res_nlo,respdf,int_sub)

       respdf = respdf*HardProc%wgt

       kin(3) = respdf(1)
       FintNLO_gq(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, 'FLM[1,2 ] 3',FintNLO_gq(3)
       
    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(4) = zero
       FintNLO_gq(4) = zero

    else

       call get_Lmu2(C1Lim,Lmu2,my_mu2)
       
       call res_tree_qqb(C1Lim%AmpMom,res_lo)
       int_sub = intsub_gq_onloqcd_1zb(z,Lmu2,1) + &
            intsub_gq_onloqcd_12(my_mu2,C1Lim)
       
       call get_respdf(gq_lumi,1,1,C1Lim,res_lo,respdf,int_sub)

       zz  = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pgq so that in the soft limit we do not have 1/[1-(1-x)]
       !-- Tr = Cf * aveqg/aveqq
       respdf = (-one)*respdf*(two/s5i)*(tr*Pgq_spav(zz)/(one-zz)) & 
            * C1Lim%wgt

       FintNLO_gq(4) = respdf(1)
       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, 'FLM[1,2 ] 4',FintNLO_gq(4)
       !pause

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_gq)
    
#if(_withchecks == 1)
    FintNNLO_subonloqcd_gq = FintNLO_gq
#endif

  end function xsect_nnlo_subonloqcd_gq_fc
  
end module mod_xsects_nnlo_subonloqcd_gq_fc
  

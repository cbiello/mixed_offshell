module mod_xsects_nnlo_subonloqcd_aq_fc
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
  use mod_eikonals
  use mod_partitions
  implicit none
  integer, parameter :: nFint  = 6 !-- 6 if full, don't care about redundancies
  integer, parameter :: nkin   = 4 
  integer, parameter :: iscale = 1  !-- which scale to load into FintNNLO_subonloqcd
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_subonloqcd_aq(6)
#endif

  private

  public :: xsect_nnlo_subonloqcd_aq_fc

contains

  function xsect_nnlo_subonloqcd_aq_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_subonloqcd_aq_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim,SLim,SC1Lim,SC2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(nFint),kin(nkin)
    real(dp)    :: respdf(ipdf),respdf_tmp(ipdf,3)
    real(dp)    :: res_nlo(2,2),res_lo(2,2),int_sub(2,2,ipdf),res_tmp(2,2)
    real(dp)    :: z,zz,s5i,eta5i,e5,eik
    real(dp)    :: Lmu2(ipdf),w11

    !integer,save :: icount = 0 !-- FC debug
    !icount = icount+1
    
    xsect_nnlo_subonloqcd_aq_fc = 0

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

    !!-- debug
    !nohistos = .false.; print*, 'debug'
    
    call open_histo()

    !---
    
    !-- FLM[z,2,3,4,5]
    
    call kinematics_nlo_is_z(yr=xx,z1=z,z2=one,HardProc=HardProc,C1Lim=C1Lim,C2Lim=C2Lim,&
         SLim=SLim,SC1Lim=SC1Lim,SC2Lim=SC2Lim,compute_etas=.false.)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_g]
    C1Lim%ids(1:4)    = [0,0,id_el,-id_el]
    C2Lim%ids(1:4)    = [0,0,id_el,-id_el]
    SLim%ids(1:4)     = [0,0,id_el,-id_el]
    SC1Lim%ids(1:4)   = [0,0,id_el,-id_el]
    SC2Lim%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_ns(1) = zero

    else    

       call get_Lmu2(HardProc,Lmu2)
       
       call res_tree_g_qqb(HardProc%AmpMom,res_nlo)
       call get_w11(HardProc,w11)
       int_sub = intsub_aq_onloqcd_z2(z,Lmu2,w11,HardProc%Lim_KinInv)
       
       call get_respdf(aq_lumi,1,1,HardProc,res_nlo,respdf,int_sub)

       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintNLO_ns(1) = respdf(iscale)

       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 1',FintNLO_ns(1)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(2) = zero
       FintNLO_ns(2) = zero

    else

       call get_Lmu2(C1Lim,Lmu2)
       
       call res_tree_qqb(C1Lim%AmpMom,res_lo)
       call get_w11(C1Lim,w11)
       int_sub = intsub_aq_onloqcd_z2(z,Lmu2,w11,C1Lim%Lim_KinInv)
       
       call get_respdf(aq_lumi,1,1,C1Lim,res_lo,respdf,int_sub)

       zz  = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(zz)/(one-zz))&
            * C1Lim%wgt

       FintNLO_ns(2) = respdf(iscale)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 2',FintNLO_ns(2)
       
    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(3) = zero
       FintNLO_ns(3) = zero

    else

       call get_Lmu2(C2Lim,Lmu2)

       call res_tree_qqb(C2Lim%AmpMom,res_lo)
       call get_w11(C2Lim,w11)
       int_sub = intsub_aq_onloqcd_z2(z,Lmu2,w11,C2Lim%Lim_KinInv)

       call get_respdf(aq_lumi,1,1,C2Lim,res_lo,respdf,int_sub)

       zz   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Cf*Pqg(zz)/(one-zz))&
            * C2Lim%wgt

       FintNLO_ns(3) = respdf(iscale)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 3',FintNLO_ns(3)
       
    endif

    !-- S and CS
    call cut_histo(SLim)
    
    if (SLim%makecut.or.SLim%flag) then

       kin(4) = zero
       FintNLO_ns(4:6) = zero

    else

       call get_Lmu2(SLim,Lmu2)

       call res_tree_qqb(SLim%AmpMom,res_lo)

       !-- S
       call get_qcd_eik(Cf,SLim%Lim_etaij,[1,2],5,eik)
       res_tmp = res_lo * eik
       call get_w11(SLim,w11)
       int_sub = intsub_aq_onloqcd_z2(z,Lmu2,w11,SLim%Lim_KinInv)

       call get_respdf(aq_lumi,1,1,SLim,res_tmp,respdf_tmp(:,1),int_sub)

       e5 = SLim%Lim_KinInv(1)
       respdf_tmp(:,1) = -respdf_tmp(:,1)/e5**2 * SLim%wgt
       FintNLO_ns(4) = respdf_tmp(iscale,1)

       !-- SC1
       res_tmp = res_lo
       call get_w11(SC1Lim,w11)
       int_sub = intsub_aq_onloqcd_z2(z,Lmu2,w11,SC1Lim%Lim_KinInv)
       call get_respdf(aq_lumi,1,1,SLim,res_tmp,respdf_tmp(:,2),int_sub)

       e5    = SC1Lim%Lim_KinInv(1)
       eta5i = SC1Lim%Lim_KinInv(2)
       respdf_tmp(:,2) = Cf*respdf_tmp(:,2)/e5**2/eta5i * SC1Lim%wgt
       FintNLO_ns(5) = respdf_tmp(iscale,2)

       !-- SC2
       res_tmp = res_lo
       call get_w11(SC2Lim,w11)
       int_sub = intsub_aq_onloqcd_z2(z,Lmu2,w11,SC2Lim%Lim_KinInv)
       call get_respdf(aq_lumi,1,1,SLim,res_tmp,respdf_tmp(:,3),int_sub)

       e5    = SC2Lim%Lim_KinInv(1)
       eta5i = SC2Lim%Lim_KinInv(2)
       respdf_tmp(:,3) = Cf*respdf_tmp(:,3)/e5**2/eta5i * SC2Lim%wgt
       FintNLO_ns(6) = respdf_tmp(iscale,3)

       respdf(:) = sum(respdf_tmp(:,1:3),2)
       kin(4) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 4',FintNLO_ns(4)
       !print *, '[z,2] 5',FintNLO_ns(5)
       !print *, '[z,2] 6',FintNLO_ns(6)
       !print *, '[]'
       !pause
       
    endif

    !--
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_ns)
    
#if(_withchecks == 1)
    FintNNLO_subonloqcd_aq = FintNLO_ns
#endif

  end function xsect_nnlo_subonloqcd_aq_fc

  !-- TC damping factor with 6||1
  subroutine get_w11(proc,damp)
   type(KinConfig), intent(in) :: proc
   real(dp), intent(out) :: damp
   real(dp) :: eta15,eta25
  
   eta15 = proc%Lim_etaij(1,5)
   eta25 = proc%Lim_etaij(2,5)
  
   damp = eta25*(one+eta15)
  
  end subroutine get_w11
  
end module mod_xsects_nnlo_subonloqcd_aq_fc
  

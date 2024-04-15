module mod_xsects_nnlo_subonloewk_ns_fc
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
  integer, parameter :: iscale = 2  !-- which scale to load into FintNNLO_subonloqcd
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_subonloewk(18)
#endif
  real(dp), parameter :: charges_ns(4,4) = reshape(& 
       [-Qdn,   Qdn,  -Qup,   Qup,   &
         Qdn,  -Qdn,   Qup,  -Qup,   &
         Q_lep, Q_lep, Q_lep, Q_lep, &
        -Q_lep,-Q_lep,-Q_lep,-Q_lep]  &
        , [4,4])

  private

  public :: xsect_nnlo_subonloewk_is_ns_ga_fc
  public :: xsect_nnlo_subonloewk_fs_53_ns_ga_fc, xsect_nnlo_subonloewk_fs_54_ns_ga_fc

contains

  !-------------------------------------------------------------
  !-- initial-state sector
  !-------------------------------------------------------------

  function xsect_nnlo_subonloewk_is_ns_ga_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_subonloewk_is_ns_ga_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim,SLim,SC1Lim,SC2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(18),kin(12)
    real(dp)    :: respdf(ipdf),respdf_tmp(ipdf,3)
    real(dp)    :: res_nlo(2,2),res_lo(2,2),res_tmp(2,2)
    real(dp)    :: zz,s5i,eik(4),e5,eta5i
    real(dp)    :: damp
    real(dp)    :: z,Lmu2(ipdf),w11,w22,int_sub(ipdf)
    
    !integer,save :: icount = 0 !-- FC debug
    !icount = icount+1
    !nohistos = .false.; print *, 'debug'

    xsect_nnlo_subonloewk_is_ns_ga_fc = 0

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

    !---

    !-- FLM[z,2,3,4]
    
    call kinematics_nlo_is_z(yr=xx,z1=z,z2=one,HardProc=HardProc,&
         C1Lim=C1Lim,C2Lim=C2Lim,SLim=SLim,SC1Lim=SC1Lim,SC2Lim=SC2Lim,compute_etas=.true.)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
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
       
       call res_tree_a_qqb(HardProc%AmpMom,res_nlo)
       call get_respdf(ns_lumi,1,1,HardProc,res_nlo,respdf)

       call partition_nlo_qed(HardProc,damp,12)

       call partition_nnlo_fact(HardProc,w11,iconf_qed=[1,2,3,4,5],i_qed=1)
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,HardProc%Lim_KinInv,0)
       respdf  = respdf*int_sub
       
       respdf = respdf*HardProc%wgt*damp

       kin(1) = respdf(1)
       FintNLO_ns(1) = respdf(iscale)

       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 1',icount, FintNLO_ns(1)

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
       
       call get_respdf(ns_lumi,1,1,C1Lim,res_lo,respdf)

       w11 = one
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,C1Lim%Lim_KinInv,0)
       respdf  = respdf*int_sub

       zz  = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       respdf = (-one)*respdf*(two/s5i)*(Pqg(zz)/(one-zz))&
            * C1Lim%wgt

       FintNLO_ns(2) = respdf(iscale)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 2',icount, FintNLO_ns(2)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(3) = zero
       FintNLO_ns(3) = zero

    else

       call get_Lmu2(C2Lim,Lmu2)
       
       call res_tree_qqb(C2Lim%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)
       
       call get_respdf(ns_lumi,1,1,C2Lim,res_lo,respdf)

       w11 = zero
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,C2Lim%Lim_KinInv,0)
       respdf  = respdf*int_sub
      
       zz  = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       respdf = (-one)*respdf*(two/s5i)*(Pqg(zz)/(one-zz))&
            * C2Lim%wgt

       FintNLO_ns(3) = respdf(iscale)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 3',icount, FintNLO_ns(3)

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
       call get_qed_eik(charges_ns,SLim%Lim_etaij,[1,2,3,4],5,eik)
       res_tmp(:,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_tmp(:,2) = res_lo(:,2) * eik(3:4) !-- up
       call get_respdf(ns_lumi,1,1,SLim,res_tmp,respdf_tmp(:,1))

       call partition_nnlo_fact(SLim,w11,iconf_qed=[1,2,3,4,5],i_qed=1)
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,SLim%Lim_KinInv,0)
       respdf_tmp(:,1)  = respdf_tmp(:,1)*int_sub
       
       call partition_nlo_qed(SLim,damp,12)
       
       e5 = SLim%Lim_KinInv(1)
       respdf_tmp(:,1) = -respdf_tmp(:,1)/e5**2 * SLim%wgt * damp
       FintNLO_ns(4) = respdf_tmp(iscale,1)

       !-- SC1
       res_tmp(:,1) = res_lo(:,1) * Qdn2
       res_tmp(:,2) = res_lo(:,2) * Qup2
       call get_respdf(ns_lumi,1,1,SLim,res_tmp,respdf_tmp(:,2))
       respdf_tmp(:,3) = respdf_tmp(:,2)
       
       w11 = one
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,SC1Lim%Lim_KinInv,0)
       respdf_tmp(:,2)  = respdf_tmp(:,2)*int_sub

       e5    = SC1Lim%Lim_KinInv(1)
       eta5i = SC1Lim%Lim_KinInv(2)
       respdf_tmp(:,2) = respdf_tmp(:,2)/e5**2/eta5i * SC1Lim%wgt
       FintNLO_ns(5) = respdf_tmp(iscale,2)

       !-- SC2
       w11 = zero
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,SC2Lim%Lim_KinInv,0)
       respdf_tmp(:,3)  = respdf_tmp(:,3)*int_sub

       e5    = SC2Lim%Lim_KinInv(1)
       eta5i = SC2Lim%Lim_KinInv(2)
       respdf_tmp(:,3) = respdf_tmp(:,3)/e5**2/eta5i * SC2Lim%wgt
       FintNLO_ns(6) = respdf_tmp(iscale,3)

       respdf(:) = sum(respdf_tmp(:,1:3),2)
       kin(4) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 4',icount, FintNLO_ns(4)
       !print *, '[z,2] 5',icount, FintNLO_ns(5)
       !print *, '[z,2] 6',icount, FintNLO_ns(6)
       !print *, '[]'

    endif

    !--
    
    !-- FLM[1,z,3,4]
    
    call kinematics_nlo_is_z(yr=xx,z1=one,z2=z,HardProc=HardProc,&
         C1Lim=C1Lim,C2Lim=C2Lim,SLim=SLim,SC1Lim=SC1Lim,SC2Lim=SC2Lim,compute_etas=.true.)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    C1Lim%ids(1:4)    = [0,0,id_el,-id_el]
    C2Lim%ids(1:4)    = [0,0,id_el,-id_el]
    SLim%ids(1:4)     = [0,0,id_el,-id_el]
    SC1Lim%ids(1:4)   = [0,0,id_el,-id_el]
    SC2Lim%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(5) = zero
       FintNLO_ns(7) = zero

    else    

       call get_Lmu2(HardProc,Lmu2)

       call res_tree_a_qqb(HardProc%AmpMom,res_nlo)
       call get_respdf(ns_lumi,1,1,HardProc,res_nlo,respdf)

       call partition_nlo_qed(HardProc,damp,12)

       call partition_nnlo_fact(HardProc,w22,iconf_qed=[1,2,3,4,5],i_qed=2)
       int_sub = intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,HardProc%Lim_KinInv,0)
       respdf  = respdf*int_sub
       
       respdf = respdf*HardProc%wgt*damp

       kin(5) = respdf(1)
       FintNLO_ns(7) = respdf(iscale)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 7',icount, FintNLO_ns(7)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(6) = zero
       FintNLO_ns(8) = zero

    else

       call get_Lmu2(C1Lim,Lmu2)

       call res_tree_qqb(C1Lim%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)
       
       call get_respdf(ns_lumi,1,1,C1Lim,res_lo,respdf)

       w22 = zero
       int_sub = intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,C1Lim%Lim_KinInv,0)
       respdf  = respdf*int_sub

       zz  = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       respdf = (-one)*respdf*(two/s5i)*(Pqg(zz)/(one-zz))&
            * C1Lim%wgt

       FintNLO_ns(8) = respdf(iscale)
       kin(6) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 8',icount, FintNLO_ns(8)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(7) = zero
       FintNLO_ns(9) = zero

    else

       call get_Lmu2(C2Lim,Lmu2)
       
       call res_tree_qqb(C2Lim%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)
       
       call get_respdf(ns_lumi,1,1,C2Lim,res_lo,respdf)

       w22 = one
       int_sub = intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,C2Lim%Lim_KinInv,0)
       respdf  = respdf*int_sub
      
       zz  = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       respdf = (-one)*respdf*(two/s5i)*(Pqg(zz)/(one-zz))&
            * C2Lim%wgt

       FintNLO_ns(9) = respdf(iscale)
       kin(7) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 9',icount, FintNLO_ns(9)
       
    endif

    !-- S and CS
    call cut_histo(SLim)
    
    if (SLim%makecut.or.SLim%flag) then

       kin(8) = zero
       FintNLO_ns(10:12) = zero

    else

       call get_Lmu2(SLim,Lmu2)
       
       call res_tree_qqb(SLim%AmpMom,res_lo)

       !-- S
       call get_qed_eik(charges_ns,SLim%Lim_etaij,[1,2,3,4],5,eik)
       res_tmp(:,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_tmp(:,2) = res_lo(:,2) * eik(3:4) !-- up
       call get_respdf(ns_lumi,1,1,SLim,res_tmp,respdf_tmp(:,1))

       call partition_nnlo_fact(SLim,w22,iconf_qed=[1,2,3,4,5],i_qed=2)
       int_sub = intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,SLim%Lim_KinInv,0)
       respdf_tmp(:,1)  = respdf_tmp(:,1)*int_sub
       
       call partition_nlo_qed(SLim,damp,12)
       
       e5 = SLim%Lim_KinInv(1)
       respdf_tmp(:,1) = -respdf_tmp(:,1)/e5**2 * SLim%wgt * damp
       FintNLO_ns(10) = respdf_tmp(iscale,1)

       !-- SC1
       res_tmp(:,1) = res_lo(:,1) * Qdn2
       res_tmp(:,2) = res_lo(:,2) * Qup2
       call get_respdf(ns_lumi,1,1,SLim,res_tmp,respdf_tmp(:,2))
       respdf_tmp(:,3) = respdf_tmp(:,2)
       
       w22 = zero
       int_sub = intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,SC1Lim%Lim_KinInv,0)
       respdf_tmp(:,2)  = respdf_tmp(:,2)*int_sub

       e5    = SC1Lim%Lim_KinInv(1)
       eta5i = SC1Lim%Lim_KinInv(2)
       respdf_tmp(:,2) = respdf_tmp(:,2)/e5**2/eta5i * SC1Lim%wgt
       FintNLO_ns(11) = respdf_tmp(iscale,2)

       !-- SC2
       w22 = one
       int_sub = intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,SC2Lim%Lim_KinInv,0)
       respdf_tmp(:,3)  = respdf_tmp(:,3)*int_sub

       e5    = SC2Lim%Lim_KinInv(1)
       eta5i = SC2Lim%Lim_KinInv(2)
       respdf_tmp(:,3) = respdf_tmp(:,3)/e5**2/eta5i * SC2Lim%wgt
       FintNLO_ns(12) = respdf_tmp(iscale,3)

       respdf(:) = sum(respdf_tmp(:,1:3),2)
       kin(8) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 10',icount, FintNLO_ns(10)
       !print *, '[1,zb] 11',icount, FintNLO_ns(11)
       !print *, '[1,zb] 12',icount, FintNLO_ns(12)
       !print *, '[]'

    endif

    !---
    
    !-- FLM[1,2,3,4]
    
    call kinematics_nlo_is_z(yr=xx,z1=one,z2=one,HardProc=HardProc,&
         C1Lim=C1Lim,C2Lim=C2Lim,SLim=SLim,SC1Lim=SC1Lim,SC2Lim=SC2Lim,compute_etas=.true.)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    C1Lim%ids(1:4)    = [0,0,id_el,-id_el]
    C2Lim%ids(1:4)    = [0,0,id_el,-id_el]
    SLim%ids(1:4)     = [0,0,id_el,-id_el]
    SC1Lim%ids(1:4)   = [0,0,id_el,-id_el]
    SC2Lim%ids(1:4)   = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(9) = zero
       FintNLO_ns(13) = zero

    else    

       call get_Lmu2(HardProc,Lmu2)

       call res_tree_a_qqb(HardProc%AmpMom,res_nlo)
       call get_respdf(ns_lumi,1,1,HardProc,res_nlo,respdf)

       call partition_nlo_qed(HardProc,damp,12)

       call partition_nnlo_fact(HardProc,w11,iconf_qed=[1,2,3,4,5],i_qed=1)
       call partition_nnlo_fact(HardProc,w22,iconf_qed=[1,2,3,4,5],i_qed=2)
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,HardProc%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,HardProc%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_12(Lmu2)
       respdf  = respdf*int_sub
       
       respdf = respdf*HardProc%wgt*damp

       kin(9) = respdf(1)
       FintNLO_ns(13) = respdf(iscale)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 13',icount, FintNLO_ns(13)

    endif

    !-- C1
    call cut_histo(C1Lim)
    
    if (C1Lim%makecut.or.C1Lim%flag) then

       kin(10) = zero
       FintNLO_ns(14) = zero

    else

       call get_Lmu2(C1Lim,Lmu2)

       call res_tree_qqb(C1Lim%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)
       
       call get_respdf(ns_lumi,1,1,C1Lim,res_lo,respdf)

       w11 = one
       w22 = zero
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,C1Lim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,C1Lim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_12(Lmu2)
       respdf  = respdf*int_sub

       zz  = C1Lim%Lim_KinInv(1)
       s5i = C1Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       respdf = (-one)*respdf*(two/s5i)*(Pqg(zz)/(one-zz))&
            * C1Lim%wgt

       FintNLO_ns(14) = respdf(iscale)
       kin(10) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 14',icount, FintNLO_ns(14)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(11) = zero
       FintNLO_ns(15) = zero

    else

       call get_Lmu2(C2Lim,Lmu2)
       
       call res_tree_qqb(C2Lim%AmpMom,res_lo)
       res_lo(1,:) = [Qdn2,Qup2] * res_lo(1,:)
       res_lo(2,:) = [Qdn2,Qup2] * res_lo(2,:)
       
       call get_respdf(ns_lumi,1,1,C2Lim,res_lo,respdf)

       w11 = zero
       w22 = one
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,C2Lim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,C2Lim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_12(Lmu2)
       respdf  = respdf*int_sub
      
       zz  = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision -> flux factor is 1-z and not z
       respdf = (-one)*respdf*(two/s5i)*(Pqg(zz)/(one-zz))&
            * C2Lim%wgt

       FintNLO_ns(15) = respdf(iscale)
       kin(11) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 15',icount, FintNLO_ns(15)

    endif

    !-- S and CS
    call cut_histo(SLim)
    
    if (SLim%makecut.or.SLim%flag) then

       kin(12) = zero
       FintNLO_ns(16:18) = zero

    else

       call get_Lmu2(SLim,Lmu2)
       
       call res_tree_qqb(SLim%AmpMom,res_lo)

       !-- S
       call get_qed_eik(charges_ns,SLim%Lim_etaij,[1,2,3,4],5,eik)
       res_tmp(:,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_tmp(:,2) = res_lo(:,2) * eik(3:4) !-- up
       call get_respdf(ns_lumi,1,1,SLim,res_tmp,respdf_tmp(:,1))

       call partition_nnlo_fact(SLim,w11,iconf_qed=[1,2,3,4,5],i_qed=1)
       call partition_nnlo_fact(SLim,w22,iconf_qed=[1,2,3,4,5],i_qed=2)
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,SLim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,SLim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_12(Lmu2)
       respdf_tmp(:,1)  = respdf_tmp(:,1)*int_sub
       
       call partition_nlo_qed(SLim,damp,12)
       
       e5 = SLim%Lim_KinInv(1)
       respdf_tmp(:,1) = -respdf_tmp(:,1)/e5**2 * SLim%wgt * damp
       FintNLO_ns(16) = respdf_tmp(iscale,1)

       !-- SC1
       res_tmp(:,1) = res_lo(:,1) * Qdn2
       res_tmp(:,2) = res_lo(:,2) * Qup2
       call get_respdf(ns_lumi,1,1,SLim,res_tmp,respdf_tmp(:,2))
       respdf_tmp(:,3) = respdf_tmp(:,2)
       
       w11 = one
       w22 = zero
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,SC1Lim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,SC1Lim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_12(Lmu2)
       respdf_tmp(:,2)  = respdf_tmp(:,2)*int_sub

       e5    = SC1Lim%Lim_KinInv(1)
       eta5i = SC1Lim%Lim_KinInv(2)
       respdf_tmp(:,2) = respdf_tmp(:,2)/e5**2/eta5i * SC1Lim%wgt
       FintNLO_ns(17) = respdf_tmp(iscale,2)

       !-- SC2
       w11 = zero
       w22 = one
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,SC2Lim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,SC2Lim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_12(Lmu2)
       respdf_tmp(:,3)  = respdf_tmp(:,3)*int_sub

       e5    = SC2Lim%Lim_KinInv(1)
       eta5i = SC2Lim%Lim_KinInv(2)
       respdf_tmp(:,3) = respdf_tmp(:,3)/e5**2/eta5i * SC2Lim%wgt
       FintNLO_ns(18) = respdf_tmp(iscale,3)

       respdf(:) = sum(respdf_tmp(:,1:3),2)
       kin(12) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 16',icount, FintNLO_ns(16)
       !print *, '[1,2] 17',icount, FintNLO_ns(17)
       !print *, '[1,2] 18',icount, FintNLO_ns(18)
       !print *, '[]'
       !pause

    endif
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_ns)
    
#if(_withchecks == 1)
    FintNNLO_subonloewk = FintNLO_ns
#endif

  end function xsect_nnlo_subonloewk_is_ns_ga_fc

  
  !-------------------------------------------------------------
  !-- final-state sectors below
  !-------------------------------------------------------------
  
  function xsect_nnlo_subonloewk_fs_53_ns_ga_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_subonloewk_fs_53_ns_ga_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nnlo_subonloewk_fs_53_ns_ga_fc = xsect_subonloewk_5i_ns_ga_fc(yRnd,ff,vegasweight,3,4)

  end function xsect_nnlo_subonloewk_fs_53_ns_ga_fc

  function xsect_nnlo_subonloewk_fs_54_ns_ga_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_subonloewk_fs_54_ns_ga_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nnlo_subonloewk_fs_54_ns_ga_fc = xsect_subonloewk_5i_ns_ga_fc(yRnd,ff,vegasweight,4,3)

  end function xsect_nnlo_subonloewk_fs_54_ns_ga_fc

  !----------------------------------------------------------------------------------------------------
  
  !-- master sector for final-state sctors below
  function xsect_subonloewk_5i_ns_ga_fc(yRnd,ff,vegasweight,icoll,jother)
    integer :: xsect_subonloewk_5i_ns_ga_fc,icoll,jother
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,CLim,CSLim,SLim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_ns(12),kin(9)
    real(dp)    :: respdf(ipdf),respdf_tmp(ipdf,2)
    real(dp)    :: res_nlo(2,2),res_lo(2,2)
    real(dp)    :: res_tmp_vect(2,2,2),respdf_vect(2,ipdf)    
    real(dp)    :: zz,s5i,eik(4),e5,eta5i
    real(dp)    :: damp
    real(dp)    :: z,Lmu2(ipdf),w11,w22,int_sub(ipdf)

    !integer,save :: icount = 0 !-- FC debug
    !icount = icount+1
    !nohistos = .false.; print *, 'debug'
    
    xsect_subonloewk_5i_ns_ga_fc = 0

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

    !-- FLM[z,2,3,4]
    
    call kinematics_nlo_fs_z(xx,z,one,icoll,jother,HardProc,CLim,CSLim,SLim)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    CLim%ids(1:4)     = [0,0,id_el,-id_el]
    SLim%ids(1:4)     = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_ns(1) = zero

    else    

       call get_Lmu2(HardProc,Lmu2)

       call res_tree_a_qqb(HardProc%AmpMom,res_nlo)
       call get_respdf(ns_lumi,1,1,HardProc,res_nlo,respdf)

       call partition_nlo_qed(HardProc,damp,icoll)

       call partition_nnlo_fact(HardProc,w11,iconf_qed=[1,2,3,4,5],i_qed=1)
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,HardProc%Lim_KinInv,0)
       respdf  = respdf*int_sub
       
       respdf = respdf*HardProc%wgt*damp

       kin(1) = respdf(1)
       FintNLO_ns(1) = respdf(iscale)

       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 1',icount, FintNLO_ns(1)

    endif

    !-- C
    call cut_histo(CLim)
    
    if (CLim%makecut.or.CLim%flag) then

       kin(2) = zero
       FintNLO_ns(2) = zero

    else

       call get_Lmu2(CLim,Lmu2)

       call res_tree_qqb(CLim%AmpMom,res_lo)

       call get_respdf(ns_lumi,1,1,CLim,res_lo,respdf)

       w11 = zero
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,CLim%Lim_KinInv,0)
       respdf = respdf*int_sub

       zz  = CLim%Lim_KinInv(1)
       s5i = CLim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision
       respdf = (-one)*respdf*(two/s5i)*(Q_lep2*Pqg(zz))&
            * CLim%wgt

       FintNLO_ns(2) = respdf(iscale)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 2',icount, FintNLO_ns(2)

    endif

    !-- S and CS
    call cut_histo(SLim)
    
    if (SLim%makecut.or.SLim%flag) then

       kin(3) = zero
       FintNLO_ns(3:4) = zero

    else

       call get_Lmu2(SLim,Lmu2)
       
       call res_tree_qqb(SLim%AmpMom,res_lo)

       !-- prepare PDFs structures, S
       call get_qed_eik(charges_ns,SLim%Lim_etaij,[1,2,3,4],5,eik)
       res_tmp_vect(:,1,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_tmp_vect(:,2,1) = res_lo(:,2) * eik(3:4) !-- dn

       !-- prepare PDFs structures, CS
       res_tmp_vect(:,:,2) = res_lo(:,:) * Q_lep2

       !-- get PDFs for all of them
       call get_respdf_vect(ns_lumi,1,1,SLim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
       
       !-- S
       call partition_nlo_qed(SLim,damp,icoll)

       call partition_nnlo_fact(SLim,w11,iconf_qed=[1,2,3,4,5],i_qed=1)
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,SLim%Lim_KinInv,0)
       respdf_vect(1,:) = respdf_vect(1,:) * int_sub
       
       e5 = SLim%Lim_KinInv(1)
       respdf_tmp(:,1) = -respdf_vect(1,:)/e5**2 * SLim%wgt * damp
       FintNLO_ns(3) = respdf_tmp(iscale,1)

       !-- CS
       e5    = CSLim%Lim_KinInv(1)
       eta5i = CSLim%Lim_KinInv(2)

       w11 = zero
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,CSLim%Lim_KinInv,0)
       respdf_vect(2,:) = respdf_vect(2,:) * int_sub
      
       respdf_tmp(:,2) = respdf_vect(2,:)/e5**2/eta5i * CSLim%wgt
       FintNLO_ns(4) = respdf_tmp(iscale,2)

       respdf(:) = sum(respdf_tmp(:,1:2),2)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[z,2] 3',icount, FintNLO_ns(3)
       !print *, '[z,2] 4',icount, FintNLO_ns(4)
       !print *, '[]'

    endif

    !---

    !-- FLM[1,z,3,4]
    
    call kinematics_nlo_fs_z(xx,one,z,icoll,jother,HardProc,CLim,CSLim,SLim)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    CLim%ids(1:4)     = [0,0,id_el,-id_el]
    SLim%ids(1:4)     = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(4) = zero
       FintNLO_ns(5) = zero

    else    

       call get_Lmu2(HardProc,Lmu2)

       call res_tree_a_qqb(HardProc%AmpMom,res_nlo)
       
       call get_respdf(ns_lumi,1,1,HardProc,res_nlo,respdf)
       
       call partition_nnlo_fact(HardProc,w22,iconf_qed=[1,2,3,4,5],i_qed=2)
       int_sub = intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,HardProc%Lim_KinInv,0)
       respdf = respdf*int_sub

       call partition_nlo_qed(HardProc,damp,icoll)
       
       respdf = respdf*HardProc%wgt*damp

       kin(4) = respdf(1)
       FintNLO_ns(5) = respdf(iscale)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 5',icount, FintNLO_ns(5)

    endif

    !-- C
    call cut_histo(CLim)
    
    if (CLim%makecut.or.CLim%flag) then

       kin(5) = zero
       FintNLO_ns(6) = zero

    else

       call get_Lmu2(CLim,Lmu2)

       call res_tree_qqb(CLim%AmpMom,res_lo)

       call get_respdf(ns_lumi,1,1,CLim,res_lo,respdf)

       w22 = zero
       int_sub = intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,CLim%Lim_KinInv,0)
       respdf = respdf*int_sub

       zz  = CLim%Lim_KinInv(1)
       s5i = CLim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision
       respdf = (-one)*respdf*(two/s5i)*(Q_lep2*Pqg(zz))&
            * CLim%wgt

       FintNLO_ns(6) = respdf(iscale)
       kin(5) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 6',icount, FintNLO_ns(6)

    endif

    !-- S and CS
    call cut_histo(SLim)
    
    if (SLim%makecut.or.SLim%flag) then

       kin(6) = zero
       FintNLO_ns(7:8) = zero

    else

       call get_Lmu2(SLim,Lmu2)
       
       call res_tree_qqb(SLim%AmpMom,res_lo)

       !-- prepare PDFs structures, S
       call get_qed_eik(charges_ns,SLim%Lim_etaij,[1,2,3,4],5,eik)
       res_tmp_vect(:,1,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_tmp_vect(:,2,1) = res_lo(:,2) * eik(3:4) !-- dn

       !-- prepare PDFs structures, CS
       res_tmp_vect(:,:,2) = res_lo(:,:) * Q_lep2

       !-- get PDFs for all of them
       call get_respdf_vect(ns_lumi,1,1,SLim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
       
       !-- S
       call partition_nlo_qed(SLim,damp,icoll)

       call partition_nnlo_fact(SLim,w22,iconf_qed=[1,2,3,4,5],i_qed=2)
       int_sub = intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,SLim%Lim_KinInv,0)
       respdf_vect(1,:) = respdf_vect(1,:) * int_sub
       
       e5 = SLim%Lim_KinInv(1)
       respdf_tmp(:,1) = -respdf_vect(1,:)/e5**2 * SLim%wgt * damp
       FintNLO_ns(7) = respdf_tmp(iscale,1)

       !-- CS
       e5    = CSLim%Lim_KinInv(1)
       eta5i = CSLim%Lim_KinInv(2)

       w22 = zero
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w22,CSLim%Lim_KinInv,0)
       respdf_vect(2,:) = respdf_vect(2,:) * int_sub
      
       respdf_tmp(:,2) = respdf_vect(2,:)/e5**2/eta5i * CSLim%wgt
       FintNLO_ns(8) = respdf_tmp(iscale,2)

       respdf(:) = sum(respdf_tmp(:,1:2),2)
       kin(6) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 7',icount, FintNLO_ns(7)
       !print *, '[1,zb] 8',icount, FintNLO_ns(8)
       !print *, '[]'

    endif

    !---

    !-- FLM[1,2,3,4]
    
    call kinematics_nlo_fs_z(xx,one,one,icoll,jother,HardProc,CLim,CSLim,SLim)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    CLim%ids(1:4)     = [0,0,id_el,-id_el]
    SLim%ids(1:4)     = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(7) = zero
       FintNLO_ns(9) = zero

    else    

       call get_Lmu2(HardProc,Lmu2)

       call res_tree_a_qqb(HardProc%AmpMom,res_nlo)
       
       call get_respdf(ns_lumi,1,1,HardProc,res_nlo,respdf)

       call partition_nnlo_fact(HardProc,w11,iconf_qed=[1,2,3,4,5],i_qed=1)
       call partition_nnlo_fact(HardProc,w22,iconf_qed=[1,2,3,4,5],i_qed=2)
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,HardProc%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,HardProc%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_12(Lmu2)
       respdf = respdf*int_sub

       call partition_nlo_qed(HardProc,damp,icoll)

       respdf = respdf*HardProc%wgt*damp

       kin(7) = respdf(1)
       FintNLO_ns(9) = respdf(iscale)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 9',icount, FintNLO_ns(9)

    endif

    !-- C
    call cut_histo(CLim)
    
    if (CLim%makecut.or.CLim%flag) then

       kin(8) = zero
       FintNLO_ns(10) = zero

    else

       call get_Lmu2(CLim,Lmu2)

       call res_tree_qqb(CLim%AmpMom,res_lo)

       call get_respdf(ns_lumi,1,1,CLim,res_lo,respdf)

       w11 = zero
       w22 = zero
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,CLim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,CLim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_12(Lmu2)
       
       respdf = respdf*int_sub

       zz  = CLim%Lim_KinInv(1)
       s5i = CLim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       !-- which loses precision
       respdf = (-one)*respdf*(two/s5i)*(Q_lep2*Pqg(zz))&
            * CLim%wgt

       FintNLO_ns(10) = respdf(iscale)
       kin(8) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 10',icount, FintNLO_ns(10)

    endif

    !-- S and CS
    call cut_histo(SLim)
    
    if (SLim%makecut.or.SLim%flag) then

       kin(9) = zero
       FintNLO_ns(11:12) = zero

    else

       call get_Lmu2(SLim,Lmu2)
       
       call res_tree_qqb(SLim%AmpMom,res_lo)

       !-- prepare PDFs structures, S
       call get_qed_eik(charges_ns,SLim%Lim_etaij,[1,2,3,4],5,eik)
       res_tmp_vect(:,1,1) = res_lo(:,1) * eik(1:2) !-- dn
       res_tmp_vect(:,2,1) = res_lo(:,2) * eik(3:4) !-- dn

       !-- prepare PDFs structures, CS
       res_tmp_vect(:,:,2) = res_lo(:,:) * Q_lep2

       !-- get PDFs for all of them
       call get_respdf_vect(ns_lumi,1,1,SLim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))
       
       !-- S
       call partition_nlo_qed(SLim,damp,icoll)

       call partition_nnlo_fact(SLim,w11,iconf_qed=[1,2,3,4,5],i_qed=1)
       call partition_nnlo_fact(SLim,w22,iconf_qed=[1,2,3,4,5],i_qed=2)
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,SLim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,SLim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_12(Lmu2)
       
       respdf_vect(1,:) = respdf_vect(1,:) * int_sub
       
       e5 = SLim%Lim_KinInv(1)
       respdf_tmp(:,1) = -respdf_vect(1,:)/e5**2 * SLim%wgt * damp
       FintNLO_ns(11) = respdf_tmp(iscale,1)

       !-- CS
       e5    = CSLim%Lim_KinInv(1)
       eta5i = CSLim%Lim_KinInv(2)

       w11 = zero
       w22 = zero
       int_sub = intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,CSLim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_1zb(z,Lmu2,w22,CSLim%Lim_KinInv,1) + &
            intsub_qqb_ga_onloewk_12(Lmu2)

       respdf_vect(2,:) = respdf_vect(2,:) * int_sub
      
       respdf_tmp(:,2) = respdf_vect(2,:)/e5**2/eta5i * CSLim%wgt
       FintNLO_ns(12) = respdf_tmp(iscale,2)

       respdf(:) = sum(respdf_tmp(:,1:2),2)
       kin(9) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '[1,2] 11',icount, FintNLO_ns(11)
       !print *, '[1,2] 12',icount, FintNLO_ns(12)
       !print *, '[]'
       !pause

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNLO_ns)
    
#if(_withchecks == 1)
    FintNNLO_subonloewk(1:12) = FintNLO_ns
#endif

  end function xsect_subonloewk_5i_ns_ga_fc

end module mod_xsects_nnlo_subonloewk_ns_fc
  

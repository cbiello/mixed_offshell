module mod_xsects_nnlo_rr_ga_fc
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_aux_sectors
  use mod_process
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_splittings_bare
  use mod_eikonals
  use mod_partitions
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_rr_tc_ga_fc(4)
#endif

  public :: xsect_nnlo_rr_5161a_ga_fc,xsect_nnlo_rr_5161c_ga_fc
  public :: xsect_nnlo_rr_5161b_ga_fc,xsect_nnlo_rr_5161d_ga_fc

  private

  !! -- g(p1) a(p2) -> e-(p3) e+(p4) q(p5) qb(p6) 

contains

  function xsect_nnlo_rr_5161b_ga_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161b_ga_fc
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5161b_ga_fc = xsect_nnlo_rr_5161bd_ga_fc(&
      yRnd,ff,vegasweight,2)

  end function xsect_nnlo_rr_5161b_ga_fc

  function xsect_nnlo_rr_5161d_ga_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5161d_ga_fc
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5161d_ga_fc = xsect_nnlo_rr_5161bd_ga_fc(&
      yRnd,ff,vegasweight,4)

  end function xsect_nnlo_rr_5161d_ga_fc

  !!*************************************************************************!!

  function xsect_nnlo_rr_5161a_ga_fc(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5161a_ga_fc
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: TCLim, TCC6Lim,C6Lim
    !--
    real(dp) :: xx(kNNLO_max_full)
    real(dp) :: FintNNLO_ga_fc(4),kin(3)
    real(dp) :: respdf(ipdf),respdf_vect(2,ipdf)
    real(dp) :: res_tmp1(1,1),res_tmp2(2,2)
    real(dp) :: z1,z5,z6,s15,s16,s56,si1
    real(dp) :: damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5161a_ga_fc = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_abcd_5i6iac(&
         xx(1:kNNLO_max_full), 1, 2, 1, &
         HardProc,TCLim=TCLim,TCC6Lim=TCC6Lim,C6Lim=C6Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_q,-id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_ga_fc(1) = zero

    else

      call res_tree_ga_ga(HardProc%AmpMom,res_tmp1(1,1))
      call get_respdf(ga_lumi,1,1,HardProc,res_tmp1,respdf)

      call partition_nnlo_qcd(HardProc,damp,iconf,1,1)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_ga_fc(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

      !print *, '1 ',FintNNLO_ga_fc(1)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          TC + TCC6                                    !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(2) = zero
      FintNNLO_ga_fc(2:3) = zero

    else
       
      call res_treeAA_aa(TCLim%AmpMom,res_tmp1(1,1))
      res_tmp1 = res_tmp1 * (Qup2*nup + Qdn2*ndn) !-- qqb/qbq taken care in Pgqqb

      call get_respdf(ga_lumi,1,1,TCLim,res_tmp1,respdf)

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s15 = TCLim%Lim_sij(1,5)
      s16 = TCLim%Lim_sij(1,6)

      !-- damp = one
      respdf_vect(1,:) = -Tr*(Pgqqb_spav_ab(-s15,-s16,s56,z1,z5,z6))*TCLim%wgt*respdf
      FintNNLO_ga_fc(2) = respdf_vect(1,1)

      !-- TC + C6
      !respdf_vect_part(2,:) = respdf*(Cf*two*Pqg(z5)/si5)*(two*Pqg(z6)/si6)*TCC6Lim%wgt
      z5 = TCC6Lim%Lim_z(1)
      z6 = TCC6Lim%Lim_z(2)
      si1 = TCC6Lim%Lim_KinInv(1)
      respdf_vect(2,:) = +Tr*(Pgq_spav(z5)*two/si1)*(Pqq(z6)*two/s16)*TCC6Lim%wgt*respdf
      FintNNLO_ga_fc(3) = respdf_vect(2,1)

      respdf = sum(respdf_vect(1:2,:),1)
      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

      !print *, '2 ',FintNNLO_ga_fc(2)
      !print *, '3 ',FintNNLO_ga_fc(3)

    endif

    !!-----------------------------------------------------------------------!!
    !!                               C6                                      !!
    !!-----------------------------------------------------------------------!!

    C6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C6Lim%npart = 5

    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(3) = zero
      FintNNLO_ga_fc(4) = zero

    else
       
      call res_tree_a_qa(C6Lim%AmpMom,res_tmp2)
      res_tmp1(1,1) = res_tmp2(1,1)*ndn + res_tmp2(1,2)*nup
      call get_respdf(ga_lumi,1,1,C6Lim,res_tmp1,respdf)

      call partition_nnlo_qcd(C6Lim,damp,iconf,1,1)

      s16 = C6Lim%Lim_sij(1,6)
      z6  = C6Lim%Lim_z(2)
      
      respdf = -respdf*Tr*(two/s16*Pqg(z6))*C6Lim%wgt*damp
            
      FintNNLO_ga_fc(4) = respdf(1)
      kin(3) = respdf(1)

      call fill_histo(respdf,vegasweight)
      
      !print *, '4 ',FintNNLO_ga_fc(4)
      
    endif

    !!-----------------------------------------------------------------------!!
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ga_fc)

#if(_withchecks == 1)
    FintNNLO_rr_tc_ga_fc = FintNNLO_ga_fc
#endif

    !if (ff(1).ne.zero) then
    !   print *, 'ff', ff(1)
    !   pause
    !endif

  end function xsect_nnlo_rr_5161a_ga_fc

  function xsect_nnlo_rr_5161c_ga_fc(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5161c_ga_fc
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: TCLim, TCC5Lim,C5Lim
    !--
    real(dp) :: xx(kNNLO_max_full)
    real(dp) :: FintNNLO_ga_fc(4),kin(3)
    real(dp) :: respdf(ipdf),respdf_vect(2,ipdf)
    real(dp) :: res_tmp1(1,1),res_tmp2(2,2)
    real(dp) :: z1,z5,z6,s15,s16,s56,si1
    real(dp) :: damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5161c_ga_fc = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_abcd_5i6iac(&
         xx(1:kNNLO_max_full), 1, 2, 3, &
         HardProc,TCLim=TCLim,TCC5Lim=TCC5Lim,C5Lim=C5Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_q,-id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_ga_fc(1) = zero

    else

      call res_tree_ga_ga(HardProc%AmpMom,res_tmp1(1,1))
      call get_respdf(ga_lumi,1,1,HardProc,res_tmp1,respdf)

      call partition_nnlo_qcd(HardProc,damp,iconf,1,1)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_ga_fc(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

      !print *, '1 ',FintNNLO_ga_fc(1)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          TC + TCC5                                    !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(2) = zero
      FintNNLO_ga_fc(2:3) = zero

    else
       
      call res_treeAA_aa(TCLim%AmpMom,res_tmp1(1,1))
      res_tmp1 = res_tmp1 * (Qup2*nup + Qdn2*ndn) !-- qqb/qbq taken care in Pgqqb

      call get_respdf(ga_lumi,1,1,TCLim,res_tmp1,respdf)

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s15 = TCLim%Lim_sij(1,5)
      s16 = TCLim%Lim_sij(1,6)

      !-- damp = one
      respdf_vect(1,:) = -Tr*(Pgqqb_spav_ab(-s15,-s16,s56,z1,z5,z6))*TCLim%wgt*respdf
      FintNNLO_ga_fc(2) = respdf_vect(1,1)

      !-- TC + C5
      !respdf_vect_part(2,:) = respdf*(Cf*two*Pqg(z5)/si5)*(two*Pqg(z6)/si6)*TCC6Lim%wgt
      z5 = TCC5Lim%Lim_z(1)
      z6 = TCC5Lim%Lim_z(2)
      si1 = TCC5Lim%Lim_KinInv(1)
      respdf_vect(2,:) = +Tr*(Pgq_spav(z6)*two/si1)*(Pqq(z5)*two/s15)*TCC5Lim%wgt*respdf
      FintNNLO_ga_fc(3) = respdf_vect(2,1)

      respdf = sum(respdf_vect(1:2,:),1)
      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

      !print *, '2 ',FintNNLO_ga_fc(2)
      !print *, '3 ',FintNNLO_ga_fc(3)

    endif

    !!-----------------------------------------------------------------------!!
    !!                               C5                                      !!
    !!-----------------------------------------------------------------------!!

    C5Lim%ids = [0,0,id_el,-id_el,-id_q,0]
    C5Lim%npart = 5

    call cut_histo(C5Lim)

    if (C5Lim%makecut) then

      kin(3) = zero
      FintNNLO_ga_fc(4) = zero

    else
       
      call res_tree_a_qa(C5Lim%AmpMom,res_tmp2)
      res_tmp1(1,1) = res_tmp2(2,1)*ndn + res_tmp2(2,2)*nup
      call get_respdf(ga_lumi,1,1,C5Lim,res_tmp1,respdf)

      call partition_nnlo_qcd(C5Lim,damp,iconf,1,1)

      s15 = C5Lim%Lim_sij(1,5)
      z5  = C5Lim%Lim_z(1)
      
      respdf = -respdf*Tr*(two/s15*Pqg(z5))*C5Lim%wgt*damp
            
      FintNNLO_ga_fc(4) = respdf(1)
      kin(3) = respdf(1)

      call fill_histo(respdf,vegasweight)
      
      !print *, '4 ',FintNNLO_ga_fc(4)
      
    endif

    !!-----------------------------------------------------------------------!!
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ga_fc)

#if(_withchecks == 1)
    FintNNLO_rr_tc_ga_fc = FintNNLO_ga_fc
#endif

    !if (ff(1).ne.zero) then
    !   print *, 'ff', ff(1)
    !   pause
    !endif

  end function xsect_nnlo_rr_5161c_ga_fc

  function xsect_nnlo_rr_5161bd_ga_fc(yRnd,ff,vegasweight,sec)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5161bd_ga_fc,sec
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: TCLim
    !--
    real(dp) :: xx(kNNLO_max_full)
    real(dp) :: FintNNLO_ga_fc(2),kin(2)
    real(dp) :: respdf(ipdf)
    real(dp) :: res_tmp1(1,1)
    real(dp) :: z1,z5,z6,s15,s16,s56
    real(dp) :: damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5161bd_ga_fc = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_abcd_5i6ibd(&
         xx(1:kNNLO_max_full), 1, 2, sec, &
         HardProc,TCLim=TCLim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_q,-id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_ga_fc(1) = zero

    else

      call res_tree_ga_ga(HardProc%AmpMom,res_tmp1(1,1))
      call get_respdf(ga_lumi,1,1,HardProc,res_tmp1,respdf)

      call partition_nnlo_qcd(HardProc,damp,iconf,1,1)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_ga_fc(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

      !print *, '1 ',FintNNLO_ga_fc(1)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                TC                                     !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(2) = zero
      FintNNLO_ga_fc(2) = zero

    else
       
      call res_treeAA_aa(TCLim%AmpMom,res_tmp1(1,1))
      res_tmp1 = res_tmp1 * (Qup2*nup + Qdn2*ndn) !-- qqb/qbq taken care in Pgqqb

      call get_respdf(ga_lumi,1,1,TCLim,res_tmp1,respdf)

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s15 = TCLim%Lim_sij(1,5)
      s16 = TCLim%Lim_sij(1,6)

      !-- damp = one
      respdf = -Tr*(Pgqqb_spav_ab(-s15,-s16,s56,z1,z5,z6))*TCLim%wgt*respdf

      FintNNLO_ga_fc(2) = respdf(1)
      kin(2) = respdf(1)
      
      call fill_histo(respdf,vegasweight)

      !print *, '2 ',FintNNLO_ga_fc(2)
      
    endif

    !!-----------------------------------------------------------------------!!
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ga_fc)

#if(_withchecks == 1)
    FintNNLO_rr_tc_ga_fc(1:2) = FintNNLO_ga_fc
    FintNNLO_rr_tc_ga_fc(3:4) = zero
#endif

    !if (ff(1).ne.zero) then
    !   print *, 'ff', ff(1)
    !   pause
    !endif

  end function xsect_nnlo_rr_5161bd_ga_fc  

end module mod_xsects_nnlo_rr_ga_fc

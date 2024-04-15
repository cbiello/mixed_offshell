module mod_xsects_nnlo_s_ag
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_aux_sectors
  use mod_process
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_hoppet_tools
  use mod_hoppet_nnlo
  use mod_splittings_bare
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_s_ag(1), FintNNLO_s_ag_onlo(3)
#endif

  public :: xsect_nnlo_s_ag_s12
  public :: xsect_nnlo_s_ag_oewk,xsect_nnlo_s_ag_oqcd

  private

contains

  !!*************************************************************************!!

  function xsect_nnlo_s_ag_s12(yRnd,ff,vegasweight)
    use mod_kinematics_lo
    integer :: xsect_nnlo_s_ag_s12
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: kin(1)
    real(dp) :: respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf),respdf_tmp(ipdf,4)
    real(dp) :: res_lo_1(1,1),res_lo_2(2,2)

    xsect_nnlo_s_ag_s12 = 0

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

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    LOProc%npart = 4

    call cut_histo(LOProc)
    if (LOProc%makecut) then

       kin(1) = zero

    else

       !-- single boost
       call res_treeAA_aa(LOProc%AmpMom,res_lo_1(1,1))
       call get_respdf_hoppet(PDFs,xPij,aa_lumi,1,1,LOProc,res_lo_1,respdf_1,myPDFs2_Lmu=[xPij_Lmu,xPij_Lmu2])
       respdf_1 = Tr * (Qup2*nup + Qdn2*ndn) * respdf_1 * LOProc%wgt

       !-- F[z1.p1,z2.p2]
       call res_tree_qqb(LOProc%AmpMom,res_lo_2)
       res_lo_2(:,1) = Qdn2 * res_lo_2(:,1)
       res_lo_2(:,2) = Qup2 * res_lo_2(:,2)

       call get_respdf_hoppet(xPij_A_1,xPij_A_2,ns_lumi,1,1,LOProc,res_lo_2,respdf_tmp(:,1),&
          myPDFs1_Lmu=[xPij_A_1_Lmu],myPDFs2_Lmu=[xPij_A_2_Lmu])

       call get_respdf_hoppet(xPij_B_1,xPij_B_2,ns_lumi,1,1,LOProc,res_lo_2,respdf_tmp(:,2),&
          myPDFs1_Lmu=[xPij_B_1_Lmu])

       call get_respdf_hoppet(xPij_C_1,xPij_C_2,ns_lumi,1,1,LOProc,res_lo_2,respdf_tmp(:,3),&
          myPDFs2_Lmu=[xPij_C_2_Lmu])

       call get_respdf_hoppet(xPij_D_1,xPij_D_2,ns_lumi,1,1,LOProc,res_lo_2,respdf_tmp(:,4))

       respdf_2 = sum(respdf_tmp,2)
       respdf_2 = Tr * xn * respdf_2 * LOProc%wgt

       respdf = respdf_1 + respdf_2

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_ag = kin
#endif

  end function xsect_nnlo_s_ag_s12

  function xsect_nnlo_s_ag_oewk(yRnd,ff,vegasweight)
    use mod_kinematics_nlo_z
    implicit none
    integer :: xsect_nnlo_s_ag_oewk
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    !--
    real(dp) :: xx(kNLO_max_full),z
    real(dp) :: kin(3),FintNNLO_s_ag_oewk(3)
    real(dp) :: EC,damp,eta51,eta52,s15,s25,z5
    real(dp) :: respdf(ipdf),intsub(ipdf),sv_logs(ipdf),PqgNLO(ipdf)
    real(dp) :: Pqg0,res_nlo(2,2),res_lo(2,2),res_aa,res_ag(1,1)

    xsect_nnlo_s_ag_oewk = 0

    ff(1) = zero
    kin = 0

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z = buff+onet*real(yRnd(kNLO_max_full),dp)

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    if ((xx(xE)*xx(xRHO).lt.buff_r) .or. (xx(xE)*(one-xx(xRHO)).lt.buff_r)) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()
    call kinematics_nlo_z_is(xx,one,z,HardProc=HardProc,&
      C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.false.)

    Pqg0 = Pgq_spav(z)
    PqgNLO = 2*((one-z)*z + Pqg0*log(one-z))

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    call cut_histo(HardProc)

    if (HardProc%makecut) then

       FintNNLO_s_ag_oewk(1) = zero
       kin(1) = zero

    else

       call res_tree_a_aq(HardProc%AmpMom,res_nlo)
       res_ag = ndn*(res_nlo(1,1)+res_nlo(2,1))+&
                nup*(res_nlo(1,2)+res_nlo(2,2))

       call get_respdf(ag_lumi,1,1,HardProc,res_ag,respdf)

       EC    = HardProc%Lim_Ei(1)
       eta51 = HardProc%Lim_etaij(1,5)
       eta52 = HardProc%Lim_etaij(2,5)
       damp = eta51 * (one + eta52)

       call fill_sv_logs(HardProc%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqgNLO + Pqg0 * log(eta52/2) * damp
       intsub = intsub - Pqg0*sv_logs
       intsub = Tr * intsub

       respdf = intsub * respdf * HardProc%wgt

       FintNNLO_s_ag_oewk(1) = respdf(1)
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 51                           !!
    !!-----------------------------------------------------------------------!!
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C1Lim)

    if (C1Lim%makecut) then

       FintNNLO_s_ag_oewk(2) = zero
       kin(2) = zero

    else

       call res_tree_qqb(C1Lim%AmpMom,res_lo)
       res_ag = Qdn2*ndn*(res_lo(1,1)+res_lo(2,1))+&
                Qup2*nup*(res_lo(1,2)+res_lo(2,2))

       call get_respdf(ag_lumi,1,1,C1Lim,res_ag,respdf)

       EC    = C1Lim%Lim_Ei(1)
       s15   = C1Lim%Lim_sij(1,5)
       z5    = C1Lim%Lim_z(1)

       call fill_sv_logs(C1Lim%muf(1)**2,4*EC**2,sv_logs)

       !-- damp = zero
       intsub = PqgNLO
       intsub = intsub - Pqg0*sv_logs
       intsub = Tr * intsub

       respdf = intsub * respdf &
              * xn * 2/s15 * Pgq_spav(z5)/(one-z5) &
              * C1Lim%wgt
       respdf = - respdf

       FintNNLO_s_ag_oewk(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim)

    if (C2Lim%makecut) then

       FintNNLO_s_ag_oewk(3) = zero
       kin(3) = zero

    else

       call res_treeAA_aa(C2Lim%AmpMom,res_aa)
       !-- factors 2 from q-qb and qb-q
       res_ag = 2 * (Qdn2*ndn + Qup2*nup) * res_aa

       call get_respdf(ag_lumi,1,1,C2Lim,res_ag,respdf)

       eta52 = C2Lim%Lim_KinInv(4)
       EC    = C2Lim%Lim_Ei(1)
       s25   = C2Lim%Lim_sij(2,5)
       z5    = C2Lim%Lim_z(2)

       call fill_sv_logs(C2Lim%muf(1)**2,4*EC**2,sv_logs)

       !-- damp = one
       intsub = PqgNLO + Pqg0 * log(eta52/2)
       intsub = intsub - Pqg0*sv_logs
       intsub = Tr * intsub

       respdf = intsub * respdf &
              * 2/s25 * Pqq(z5)/(one-z5) &
              * C2Lim%wgt
       respdf = - respdf

       FintNNLO_s_ag_oewk(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_ag_onlo = FintNNLO_s_ag_oewk
#endif

  end function xsect_nnlo_s_ag_oewk

  function xsect_nnlo_s_ag_oqcd(yRnd,ff,vegasweight)
    use mod_kinematics_nlo_z
    implicit none
    integer :: xsect_nnlo_s_ag_oqcd
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNNLO_s_ag_oqcd(2),kin(2)
    real(dp)    :: respdf(ipdf),PqgNLO(ipdf),tildePqgNLO(ipdf),sv_logs(ipdf)
    real(dp)    :: res_nlo(2,2),res_nlo_ag(1,1)
    real(dp)    :: z,zz,s25,Pqg0,EC

    xsect_nnlo_s_ag_oqcd = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z = buff+onet*real(yRnd(kNLO_max_full),dp)

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    if (xx(xE)*xx(xRHO).lt.buff_r) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_z_is(xx,z,one,HardProc=HardProc,C2Lim=C2Lim,compute_etas=.false.)

    Pqg0 = Pgq_spav(z)
    PqgNLO = 2*((one-z)*z + Pqg0*log(one-z))

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNNLO_s_ag_oqcd(1) = zero

    else    

       call res_tree_g_qg(HardProc%AmpMom,res_nlo)

       res_nlo_ag(1,1) = ndn * Qdn2 * sum(res_nlo(:,1)) &
                       + nup * Qup2 * sum(res_nlo(:,2))

       call get_respdf(ag_lumi,1,1,HardProc,res_nlo_ag,respdf)

       EC = HardProc%Lim_Ei(1)
       call fill_sv_logs(HardProc%muf(1)**2,4*EC**2,sv_logs)

       tildePqgNLO = xn * (PqgNLO - Pqg0 * sv_logs)

       respdf = tildePqgNLO * respdf * HardProc%wgt

       kin(1) = respdf(1)
       FintNNLO_s_ag_oqcd(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(2) = zero
       FintNNLO_s_ag_oqcd(2) = zero

    else

       call res_tree_qqb(C2Lim%AmpMom,res_nlo)

       res_nlo_ag(1,1) = ndn * Qdn2 * sum(res_nlo(:,1)) &
                       + nup * Qup2 * sum(res_nlo(:,2))

       call get_respdf(ag_lumi,1,1,C2Lim,res_nlo_ag,respdf)

       EC = C2Lim%Lim_Ei(1)
       call fill_sv_logs(C2Lim%muf(1)**2,4*EC**2,sv_logs)

       tildePqgNLO = xn * (PqgNLO - Pqg0 * sv_logs)

       zz  = C2Lim%Lim_z(2)
       s25 = C2Lim%Lim_sij(2,5)

       !-- use Pgq so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = tildePqgNLO * respdf * C2Lim%wgt &
              * (2/s25)*(Tr*Pgq_spav(zz)/(one-zz))
       respdf = - respdf

       FintNNLO_s_ag_oqcd(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_s_ag_oqcd)
    
#if(_withchecks == 1)
    FintNNLO_s_ag_onlo(1:2) = FintNNLO_s_ag_oqcd
#endif

  end function xsect_nnlo_s_ag_oqcd

  !-- Hoppet version
  !---------------------------------------------------------------------------!
  function xsect_nnlo_s_ag_oqcd_hppt(yRnd,ff,vegasweight)
    use mod_kinematics_nlo
    implicit none
    integer :: xsect_nnlo_s_ag_oqcd_hppt
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNNLO_s_ag_oqcd(2),kin(2)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(2,2),res_lo(2,2)
    real(dp)    :: z,s5i

    xsect_nnlo_s_ag_oqcd_hppt = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    if (xx(xE)*xx(xRHO).lt.buff_r) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C2Lim=C2Lim,compute_etas=.false.)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    call cut_histo(HardProc)

    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNNLO_s_ag_oqcd(1) = zero

    else    

       call res_tree_g_qg(HardProc%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf_hoppet(xPij,PDFs,qg_lumi,1,1,HardProc,res_nlo,respdf,myPDFs1_Lmu=[xPij_Lmu])

       respdf = xn * respdf * HardProc%wgt

       kin(1) = respdf(1)
       FintNNLO_s_ag_oqcd(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(2) = zero
       FintNNLO_s_ag_oqcd(2) = zero

    else

       call res_tree_qqb(C2Lim%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       call get_respdf_hoppet(xPij,PDFs,qg_lumi,1,1,C2Lim,res_lo,respdf,myPDFs1_Lmu=[xPij_Lmu])

       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pgq so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = xn * respdf * C2Lim%wgt &
              * (2/s5i)*(Tr*Pgq_spav(z)/(one-z))
       respdf = - respdf

       FintNNLO_s_ag_oqcd(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_s_ag_oqcd)
    
#if(_withchecks == 1)
    FintNNLO_s_ag_onlo(1:2) = FintNNLO_s_ag_oqcd
#endif

  end function xsect_nnlo_s_ag_oqcd_hppt

end module mod_xsects_nnlo_s_ag

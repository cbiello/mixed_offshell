module mod_xsects_nnlo_s_ns
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_aux_sectors
  use mod_process
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_amplitudes_loop_ppll
  use mod_splittings_bare
  use mod_eikonals
  use mod_partitions
  use mod_hoppet_tools
  use mod_hoppet_nnlo
  use mod_int_sub_nnlo
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_s_ns_qqx(3),FintNNLO_s_ns_onlo(18)
#endif

  public :: xsect_nnlo_s_ns_s12
  public :: xsect_nnlo_s_ns_vqcd,xsect_nnlo_s_ns_vewk
  public :: xsect_nnlo_s_ns_oqcd
  public :: xsect_nnlo_s_ns_oewk_is,xsect_nnlo_s_ns_oewk_fs_53,xsect_nnlo_s_ns_oewk_fs_54
  public :: xsect_nnlo_s_ns_qqb,xsect_nnlo_s_ns_qq

  public :: xsect_nnlo_s_ns_vewknf

  private

  real(dp), parameter :: charges(4,4) = reshape([ & 
       -Qdn,   Qdn,  -Qup,   Qup,     &
        Qdn,  -Qdn,   Qup,  -Qup,     &
        Q_lep, Q_lep, Q_lep, Q_lep,   &
       -Q_lep,-Q_lep,-Q_lep,-Q_lep ], &
       [4,4])

contains

  !!*************************************************************************!!

  function xsect_nnlo_s_ns_oewk_fs_53(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_s_ns_oewk_fs_53
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nnlo_s_ns_oewk_fs_53 = xsect_nnlo_s_ns_oewk_fs_5i(yRnd,ff,vegasweight,3,4)

  end function xsect_nnlo_s_ns_oewk_fs_53

  function xsect_nnlo_s_ns_oewk_fs_54(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_s_ns_oewk_fs_54
    real(dp15) :: yRnd(30),ff(1),vegasweight

    xsect_nnlo_s_ns_oewk_fs_54 = xsect_nnlo_s_ns_oewk_fs_5i(yRnd,ff,vegasweight,4,3)

  end function xsect_nnlo_s_ns_oewk_fs_54


  !!*************************************************************************!!

  function xsect_nnlo_s_ns_s12(yRnd,ff,vegasweight)
    use mod_kinematics_lo_z
    implicit none
    integer :: xsect_nnlo_s_ns_s12
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc,LOProc_z1,LOProc_z2
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: kin(3),FintNNLO_s_ns_s12(3)
    real(dp) :: respdf(ipdf),respdf_bak(ipdf)
    real(dp) :: respdf_1(ipdf),respdf_2(ipdf),respdf_3(ipdf),respdf_4(ipdf)
    real(dp) :: respdf_tmp(ipdf,4)
    real(dp) :: z1,z2,t,u,EC
    real(dp) :: ns_int_sub(2,2),sv_logs(ipdf),SplitF_1(ipdf),SplitF_2(ipdf)
    real(dp) :: Pqq0_R_1(-1:1),PqqNLO_1(-1:1),Pqq0_R_2(-1:1),PqqNLO_2(-1:1)
    real(dp) :: res_lo(2,2),res_tmp(2,2)

    xsect_nnlo_s_ns_s12 = 0

    ff(1) = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    z1 = buff+onet*real(yRnd(kLO_max_full),dp)
    z2 = buff+onet*real(yRnd(kLO_max_full+1),dp)

#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()
    !-- z = 1 kinematics
    call kinematics_lo_z(xx,one,one,LOProc)
    !-- boosted kinematics
    call kinematics_lo_z(xx,z1,one,LOProc_z1)
    call kinematics_lo_z(xx,one,z2,LOProc_z2)

    !-- Splittings and integrated sub terms
    Pqq0_R_1 = PqqAP_0_R(z1)
    PqqNLO_1 = Pqq_NLO(z1)
    Pqq0_R_2 = PqqAP_0_R(z2)
    PqqNLO_2 = Pqq_NLO(z2)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    LOProc%npart = 4

    call cut_histo(LOProc)
    if (LOProc%makecut) then

       kin(1) = zero
       FintNNLO_s_ns_s12(1) = zero

    else

       call res_tree_qqb(LOProc%AmpMom,res_tmp)
       res_lo(:,1) = Qdn2 * res_tmp(:,1)
       res_lo(:,2) = Qup2 * res_tmp(:,2)

       !-- F[z1.p1,z2.p2]
       call get_respdf_hoppet(xPij_A_1,xPij_A_2,ns_lumi,1,1,LOProc,res_lo,respdf_tmp(:,1),&
          myPDFs1_Lmu=[xPij_A_1_Lmu],myPDFs2_Lmu=[xPij_A_2_Lmu])

       call get_respdf_hoppet(xPij_B_1,xPij_B_2,ns_lumi,1,1,LOProc,res_lo,respdf_tmp(:,2),&
          myPDFs1_Lmu=[xPij_B_1_Lmu])

       call get_respdf_hoppet(xPij_C_1,xPij_C_2,ns_lumi,1,1,LOProc,res_lo,respdf_tmp(:,3),&
          myPDFs2_Lmu=[xPij_C_2_Lmu])

       call get_respdf_hoppet(xPij_D_1,xPij_D_2,ns_lumi,1,1,LOProc,res_lo,respdf_tmp(:,4))

       respdf_1 = 2 * sum(respdf_tmp,2)

       !-- F[z1.p1,p2]: Qq^2 part
       call get_respdf_hoppet(xPij,PDFs,ns_lumi,1,1,LOProc,res_lo,respdf_tmp(:,1),&
            myPDFs1_Lmu=[xPij_Lmu,xPij_Lmu2])

       !-- F[p1,z2.p2]: Qq^2 part
       call get_respdf_hoppet(PDFs,xPij,ns_lumi,1,1,LOProc,res_lo,respdf_tmp(:,2),&
            myPDFs2_Lmu=[xPij_Lmu,xPij_Lmu2])

       !-- F[z1.p1,p2]/F[p1,z2.p2]: Qq*Qe part
       res_lo(1,1) =  Qdn * res_tmp(1,1); res_lo(2,1) = -Qdn * res_tmp(2,1)
       res_lo(1,2) =  Qup * res_tmp(1,2); res_lo(2,2) = -Qup * res_tmp(2,2)
       res_lo = Q_lep * res_lo

       call get_respdf_hoppet(xPij_B_1,PDFs,ns_lumi,1,1,LOProc,res_lo,respdf_tmp(:,3),&
            myPDFs1_Lmu=[xPij_B_1_Lmu])

       call get_respdf_hoppet(PDFs,xPij_C_2,ns_lumi,1,1,LOProc,res_lo,respdf_tmp(:,4),&
            myPDFs2_Lmu=[xPij_C_2_Lmu])

       t = -2*scr(LOProc%AmpMom(:,1),LOProc%AmpMom(:,3))
       u = -2*scr(LOProc%AmpMom(:,1),LOProc%AmpMom(:,4))

       respdf_tmp(:,3:4) = respdf_tmp(:,3:4) * 2*log(u/t)

       respdf_2 = sum(respdf_tmp,2)

       !-- subtract plus
       ns_int_sub = reshape(calG_QqQl(LOProc,one,[Qdn,-Qdn, Qup,-Qup],Q_lep),[2,2])
       res_lo(:,1) = ns_int_sub(:,1) * res_tmp(:,1)
       res_lo(:,2) = ns_int_sub(:,2) * res_tmp(:,2)

       call get_respdf_hoppet(PDFs,PDFs,ns_lumi,1,1,LOProc,res_lo,respdf_bak)

       EC = LOProc%Lim_Ei(1)
       call fill_sv_logs(LOProc%muf(1)**2,4*EC**2,sv_logs)
       SplitF_1 = - (PqqNLO_1(1) + Pqq0_R_1(1)*sv_logs)
       SplitF_2 = - (PqqNLO_2(1) + Pqq0_R_2(1)*sv_logs)

       respdf_3 = (SplitF_1 + SplitF_2) * respdf_bak

       !-- Elastic contribution
       respdf_4 = 4 * zeta2 * respdf_bak

       ns_int_sub = reshape(calG_QqQl_L(LOProc,[Qdn,-Qdn, Qup,-Qup],Q_lep),[2,2])
       res_lo(:,1) = ns_int_sub(:,1) * res_tmp(:,1)
       res_lo(:,2) = ns_int_sub(:,2) * res_tmp(:,2)

       call get_respdf_hoppet(PDFs,PDFs,ns_lumi,1,1,LOProc,res_lo,respdf_bak)
       respdf_4 = respdf_4 + respdf_bak*sv_logs

       !-- Total
       respdf = CF * (respdf_1 + respdf_2 + respdf_3 + respdf_4) * LOProc%wgt

       FintNNLO_s_ns_s12(1) = respdf(1)
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    !-- F[z1.p1,p2]
    LOProc_z1%ids(1:4) = [0,0,id_el,-id_el]
    LOProc_z1%npart = 4

    call cut_histo(LOProc_z1)
    if (LOProc_z1%makecut) then

       kin(2) = zero
       FintNNLO_s_ns_s12(2) = zero

    else

       call res_tree_qqb(LOProc_z1%AmpMom,res_lo)

       ns_int_sub = reshape(calG_QqQl(LOProc_z1,one,[Qdn,-Qdn, Qup,-Qup],Q_lep),[2,2])
       res_lo(:,1) = ns_int_sub(:,1) * res_lo(:,1)
       res_lo(:,2) = ns_int_sub(:,2) * res_lo(:,2)

       call get_respdf_hoppet(PDFs,PDFs,ns_lumi,1,1,LOProc_z1,res_lo,respdf)

       EC = LOProc_z1%Lim_Ei(1)
       call fill_sv_logs(LOProc_z1%muf(1)**2,4*EC**2,sv_logs)
       SplitF_1 = PqqNLO_1(0) + Pqq0_R_1(0)*sv_logs

       respdf = CF * SplitF_1 * respdf * LOProc_z1%wgt

       FintNNLO_s_ns_s12(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    !-- F[p1,z2.p2]
    LOProc_z2%ids(1:4) = [0,0,id_el,-id_el]
    LOProc_z2%npart = 4

    call cut_histo(LOProc_z2)
    if (LOProc_z2%makecut) then

       kin(3) = zero
       FintNNLO_s_ns_s12(3) = zero

    else

       call res_tree_qqb(LOProc_z2%AmpMom,res_lo)

       ns_int_sub = reshape(calG_QqQl(LOProc_z2,one,[Qdn,-Qdn, Qup,-Qup],Q_lep),[2,2])
       res_lo(:,1) = ns_int_sub(:,1) * res_lo(:,1)
       res_lo(:,2) = ns_int_sub(:,2) * res_lo(:,2)

       call get_respdf_hoppet(PDFs,PDFs,ns_lumi,1,1,LOProc_z2,res_lo,respdf)

       EC = LOProc_z2%Lim_Ei(2)
       call fill_sv_logs(LOProc_z2%muf(1)**2,4*EC**2,sv_logs)
       SplitF_2 = PqqNLO_2(0) + Pqq0_R_2(0)*sv_logs

       respdf = CF * SplitF_2 * respdf * LOProc_z2%wgt

       FintNNLO_s_ns_s12(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_ns_qqx = FintNNLO_s_ns_s12
#endif

  end function xsect_nnlo_s_ns_s12

  function xsect_nnlo_s_ns_vqcd(yRnd,ff,vegasweight)
    use mod_kinematics_lo
    implicit none
    integer :: xsect_nnlo_s_ns_vqcd
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf),respdf_3(ipdf)
    real(dp) :: etab13,etab14,IntSubs(2,2)
    real(dp) :: res_Isub(2,2),res_tree(2,2),res_nlo(2,2),kin(1)

    xsect_nnlo_s_ns_vqcd = 0

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

       call res_qcdloop_qqb(LOProc%AmpMom,res_tree,res_nlo)

       !-- Elastic contribution
       !-- etabij = one - etaij
       etab13 = half*(one + LOProc%AmpMom(4,3)/LOProc%AmpMom(1,3))
       etab14 = half*(one + LOProc%AmpMom(4,4)/LOProc%AmpMom(1,4))

       IntSubs = reshape(fin_ns_vqcd([etab13,etab14],[Qdn,-Qdn, Qup,-Qup],Qel),[2,2])
       res_Isub(:,1) = IntSubs(:,1) * res_nlo(:,1) !-- d,dx
       res_Isub(:,2) = IntSubs(:,2) * res_nlo(:,2) !-- u,ux

       call get_respdf_hoppet(PDFs,PDFs,ns_lumi,1,1,LOProc,res_Isub,respdf_1)

       !-- Add electric charges
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       !-- F[z1,2]
       call get_respdf_hoppet(xPij,PDFs,ns_lumi,1,1,LOProc,res_nlo,respdf_2,myPDFs1_Lmu=[xPij_Lmu])

       !-- F[1,z2]
       call get_respdf_hoppet(PDFs,xPij,ns_lumi,1,1,LOProc,res_nlo,respdf_3,myPDFs2_Lmu=[xPij_Lmu])

       !-- Total
       respdf = (respdf_1 + respdf_2 + respdf_3) * LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_ns_qqx(1) = kin(1)
#endif

  end function xsect_nnlo_s_ns_vqcd

  function xsect_nnlo_s_ns_vewk(yRnd,ff,vegasweight)
    use mod_kinematics_lo
    implicit none
    integer :: xsect_nnlo_s_ns_vewk
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf)
    real(dp) :: res_tree(2,3),res_nlo(2,3),kin(1)

    xsect_nnlo_s_ns_vewk = 0

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

       call res_ewkloop_qqb(LOProc%AmpMom,res_tree,res_nlo)

       !-- F[z1,2]
       call get_respdf_hoppet(xPij,PDFs,ns_lumi_splitb,1,1,LOProc,res_nlo,respdf_1,myPDFs1_Lmu=[xPij_Lmu])

       !-- F[1,z2]
       call get_respdf_hoppet(PDFs,xPij,ns_lumi_splitb,1,1,LOProc,res_nlo,respdf_2,myPDFs2_Lmu=[xPij_Lmu])

       !-- Total
       respdf = CF * (respdf_1 + respdf_2) * LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_ns_qqx(1) = kin(1)
#endif

  end function xsect_nnlo_s_ns_vewk

  function xsect_nnlo_s_ns_qqb(yRnd,ff,vegasweight)
    use mod_kinematics_lo
    implicit none
    integer :: xsect_nnlo_s_ns_qqb
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf)
    real(dp)    :: res_lo(2,2)

    xsect_nnlo_s_ns_qqb = 0

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

       call res_tree_qqb(LOProc%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       call get_respdf_hoppet(xPij,PDFs,ns_lumi,1,1,LOProc,res_lo,respdf_1,myPDFs1_Lmu=[xPij_Lmu])

       call get_respdf_hoppet(PDFs,xPij,ns_lumi,1,1,LOProc,res_lo,respdf_2,myPDFs2_Lmu=[xPij_Lmu])

       respdf = CF * (respdf_1 + respdf_2) * LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_ns_qqx(1) = kin(1)
#endif

  end function xsect_nnlo_s_ns_qqb

  function xsect_nnlo_s_ns_qq(yRnd,ff,vegasweight)
    use mod_kinematics_lo
    implicit none
    integer :: xsect_nnlo_s_ns_qq
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf)
    real(dp)    :: res_lo(2,2), res_swp(2,2)
    
    xsect_nnlo_s_ns_qq = 0

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

       call res_tree_qqb(LOProc%AmpMom,res_lo)
       res_lo(:,1) = Qdn2 * res_lo(:,1)
       res_lo(:,2) = Qup2 * res_lo(:,2)

       res_swp(1,1) = res_lo(2,1); res_swp(2,1) = res_lo(1,1); !-- d <--> db
       res_swp(1,2) = res_lo(2,2); res_swp(2,2) = res_lo(1,2); !-- u <--> ub

       call get_respdf_hoppet(xPij,PDFs,ns_lumi,1,1,LOProc,res_swp,respdf_1,myPDFs1_Lmu=[xPij_Lmu])

       call get_respdf_hoppet(PDFs,xPij,ns_lumi,1,1,LOProc,res_lo,respdf_2,myPDFs2_Lmu=[xPij_Lmu])

       respdf = CF * (respdf_1 + respdf_2) * LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_ns_qqx(1) = kin(1)
#endif

  end function xsect_nnlo_s_ns_qq

  function xsect_nnlo_s_ns_oqcd(yRnd,ff,vegasweight)
    use mod_kinematics_nlo_z
    implicit none
    integer :: xsect_nnlo_s_ns_oqcd
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim
    type(KinConfig) :: HardProc_z1,C1Lim_z1,C2Lim_z1
    type(KinConfig) :: HardProc_z2,C1Lim_z2,C2Lim_z2
    !--
    real(dp) :: xx(kNLO_max_full),z
    real(dp) :: kin(9),FintNNLO_s_ns_oqcd(9)
    real(dp) :: damp,EC,eta51,eta52,s15,s25,z5
    real(dp) :: respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf),respdf_bak(ipdf)
    real(dp) :: intsub(ipdf),intsub_1(ipdf),intsub_2(ipdf),sv_logs(ipdf)
    real(dp) :: Pqq0_R(-1:1),PqqNLO(-1:1)
    real(dp) :: ns_int_sub(2,2)
    real(dp) :: res_nlo(2,2),res_tmp(2,2)

    xsect_nnlo_s_ns_oqcd = 0

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

    !-- Boosted kinematics - z1
    call kinematics_nlo_z_is(xx,z,one,HardProc=HardProc_z1,&
      C1Lim=C1Lim_z1,C2Lim=C2Lim_z1,compute_etas=.true.)

    !-- Boosted kinematics - z2
    call kinematics_nlo_z_is(xx,one,z,HardProc=HardProc_z2,&
      C1Lim=C1Lim_z2,C2Lim=C2Lim_z2,compute_etas=.true.)

    !-- z1=z2=1 kinematics
    call kinematics_nlo_z_is(xx,one,one,HardProc=HardProc,  &
      C1Lim=C1Lim,C2Lim=C2Lim,compute_etas=.true.)

    !-- Splittings and integrated sub terms
    Pqq0_R = PqqAP_0_R(z)
    PqqNLO = Pqq_NLO(z)

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !!                        FLM[z.1_q,2_qb,3,4|5_g]                        !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc_z1%ids(1:5) = [0,0,id_el,-id_el,id_g]
    call cut_histo(HardProc_z1)

    if (HardProc_z1%makecut) then

       FintNNLO_s_ns_oqcd(1) = zero
       kin(1) = zero

    else

       call res_tree_g_qqb(HardProc_z1%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(ns_lumi,1,1,HardProc_z1,res_nlo,respdf)

       !-- partition function in the limit: 6 || 1
       HardProc_z1%Lim_etaij(1,6) = zero
       HardProc_z1%Lim_etaij(2,6) = HardProc_z1%Lim_etaij(1,2)
       HardProc_z1%Lim_etaij(3,6) = HardProc_z1%Lim_etaij(1,3)
       HardProc_z1%Lim_etaij(4,6) = HardProc_z1%Lim_etaij(1,4)

       call partition_nnlo_fact(HardProc_z1,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=1,i_qed=1)

       EC = HardProc_z1%Lim_Ei(1)
       eta51 = HardProc_z1%Lim_etaij(1,5)
       call fill_sv_logs(HardProc_z1%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta51) * damp
       intsub = intsub + Pqq0_R(0)*sv_logs

       respdf = intsub * respdf * HardProc_z1%wgt

       FintNNLO_s_ns_oqcd(1) = respdf(1)
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 51                           !!
    !!-----------------------------------------------------------------------!!
    C1Lim_z1%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C1Lim_z1)

    if (C1Lim_z1%makecut) then

       FintNNLO_s_ns_oqcd(2) = zero
       kin(2) = zero

    else

       call res_tree_qqb(C1Lim_z1%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(ns_lumi,1,1,C1Lim_z1,res_nlo,respdf)

       !-- damp = one
       EC = C1Lim_z1%Lim_Ei(1)
       eta51 = C1Lim_z1%Lim_KinInv(4)
       call fill_sv_logs(C1Lim_z1%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta51)
       intsub = intsub + Pqq0_R(0)*sv_logs

       s15 = C1Lim_z1%Lim_sij(1,5)
       z5  = C1Lim_z1%Lim_z(1)
       respdf = intsub * respdf * C1Lim_z1%wgt &
              * CF * 2/s15 * Pqg(z5)/(one-z5)
       respdf = -respdf

       FintNNLO_s_ns_oqcd(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim_z1%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim_z1)

    if (C2Lim_z1%makecut) then

       FintNNLO_s_ns_oqcd(3) = zero
       kin(3) = zero

    else

       call res_tree_qqb(C2Lim_z1%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(ns_lumi,1,1,C2Lim_z1,res_nlo,respdf)

       !-- damp = zero
       EC = C2Lim_z1%Lim_Ei(1)
       eta51 = C2Lim_z1%Lim_KinInv(4)
       call fill_sv_logs(C2Lim_z1%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0)
       intsub = intsub + Pqq0_R(0)*sv_logs

       s25 = C2Lim_z1%Lim_sij(2,5)
       z5  = C2Lim_z1%Lim_z(2)
       respdf = intsub * respdf * C2Lim_z1%wgt &
              * CF * 2/s25 * Pqg(z5)/(one-z5)
       respdf = -respdf

       FintNNLO_s_ns_oqcd(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !!                        FLM[1_q,z.2_qb,3,4|5_g]                        !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc_z2%ids(1:5) = [0,0,id_el,-id_el,id_g]
    call cut_histo(HardProc_z2)

    if (HardProc_z2%makecut) then

       FintNNLO_s_ns_oqcd(4) = zero
       kin(4) = zero

    else

       call res_tree_g_qqb(HardProc_z2%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(ns_lumi,1,1,HardProc_z2,res_nlo,respdf)

       !-- partition function in the limit: 6 || 2
       HardProc_z2%Lim_etaij(1,6) = HardProc_z2%Lim_etaij(1,2)
       HardProc_z2%Lim_etaij(2,6) = zero
       HardProc_z2%Lim_etaij(3,6) = HardProc_z2%Lim_etaij(2,3)
       HardProc_z2%Lim_etaij(4,6) = HardProc_z2%Lim_etaij(2,4)

       call partition_nnlo_fact(HardProc_z2,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=2,i_qed=2)

       EC = HardProc_z2%Lim_Ei(1)
       eta52 = HardProc_z2%Lim_etaij(2,5)
       call fill_sv_logs(HardProc_z2%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta52) * damp
       intsub = intsub + Pqq0_R(0)*sv_logs

       respdf = intsub * respdf * HardProc_z2%wgt

       FintNNLO_s_ns_oqcd(4) = respdf(1)
       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 51                           !!
    !!-----------------------------------------------------------------------!!
    C1Lim_z2%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C1Lim_z2)

    if (C1Lim_z2%makecut) then

       FintNNLO_s_ns_oqcd(5) = zero
       kin(5) = zero

    else

       call res_tree_qqb(C1Lim_z2%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(ns_lumi,1,1,C1Lim_z2,res_nlo,respdf)

       !-- damp = zero
       EC = C1Lim_z2%Lim_Ei(1)
       call fill_sv_logs(C1Lim_z2%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0)
       intsub = intsub + Pqq0_R(0)*sv_logs

       s15 = C1Lim_z2%Lim_sij(1,5)
       z5  = C1Lim_z2%Lim_z(1)
       respdf = intsub * respdf * C1Lim_z2%wgt &
              * CF * 2/s15 * Pqg(z5)/(one-z5)
       respdf = -respdf

       FintNNLO_s_ns_oqcd(5) = respdf(1)
       kin(5) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim_z2%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim_z2)

    if (C2Lim_z2%makecut) then

       FintNNLO_s_ns_oqcd(6) = zero
       kin(6) = zero

    else

       call res_tree_qqb(C2Lim_z2%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(ns_lumi,1,1,C2Lim_z2,res_nlo,respdf)

       !-- damp = one
       EC = C2Lim_z2%Lim_Ei(1)
       call fill_sv_logs(C2Lim_z2%muf(1)**2,4*EC**2,sv_logs)

       eta52 = C2Lim_z2%Lim_KinInv(4)
       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta52)
       intsub = intsub + Pqq0_R(0)*sv_logs

       s25 = C2Lim_z2%Lim_sij(2,5)
       z5  = C2Lim_z2%Lim_z(2)
       respdf = intsub * respdf * C2Lim_z2%wgt &
              * CF * 2/s25 * Pqg(z5)/(one-z5)
       respdf = -respdf

       FintNNLO_s_ns_oqcd(6) = respdf(1)
       kin(6) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !!                         FLM[1_q,2_qb,3,4|5_g]                         !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_g]
    call cut_histo(HardProc)

    if (HardProc%makecut) then

       FintNNLO_s_ns_oqcd(7) = zero
       kin(7) = zero

    else

       call res_tree_g_qqb(HardProc%AmpMom,res_tmp)
       res_nlo(:,1) = Qdn2 * res_tmp(:,1)
       res_nlo(:,2) = Qup2 * res_tmp(:,2)

       call get_respdf(ns_lumi,1,1,HardProc,res_nlo,respdf_bak)

       EC = HardProc%Lim_Ei(1)
       call fill_sv_logs(HardProc%muf(1)**2,4*EC**2,sv_logs)

       !-- subtract plus, leg 1
       HardProc%Lim_etaij(1,6) = zero
       HardProc%Lim_etaij(2,6) = HardProc%Lim_etaij(1,2)
       HardProc%Lim_etaij(3,6) = HardProc%Lim_etaij(1,3)
       HardProc%Lim_etaij(4,6) = HardProc%Lim_etaij(1,4)

       call partition_nnlo_fact(HardProc,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=1,i_qed=1)
       eta51 = HardProc%Lim_etaij(1,5)
       intsub_1 = PqqNLO(1) - Pqq0_R(1) * log(eta51) * damp
       intsub_1 = intsub_1 + Pqq0_R(1)*sv_logs
       intsub_1 = - intsub_1

       !-- subtract plus, leg 2
       HardProc%Lim_etaij(1,6) = HardProc%Lim_etaij(1,2)
       HardProc%Lim_etaij(2,6) = zero
       HardProc%Lim_etaij(3,6) = HardProc%Lim_etaij(2,3)
       HardProc%Lim_etaij(4,6) = HardProc%Lim_etaij(2,4)

       call partition_nnlo_fact(HardProc,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=2,i_qed=2)
       eta52 = HardProc%Lim_etaij(2,5)
       intsub_2 = PqqNLO(1) - Pqq0_R(1) * log(eta52) * damp
       intsub_2 = intsub_2 + Pqq0_R(1)*sv_logs
       intsub_2 = - intsub_2

       respdf_1 = (intsub_1 + intsub_2) * respdf_bak * HardProc%wgt

       !-- elastic component (soft photon)
       ns_int_sub = reshape(calG_QqQl(HardProc,one,[Qdn,-Qdn, Qup,-Qup],Q_lep),[2,2])
       res_nlo(:,1) = ns_int_sub(:,1) * res_tmp(:,1)
       res_nlo(:,2) = ns_int_sub(:,2) * res_tmp(:,2)

       call get_respdf(ns_lumi,1,1,HardProc,res_nlo,respdf_2)
       respdf_2 = (respdf_2 - 3*respdf_bak*sv_logs) * HardProc%wgt

       !-- total
       respdf = respdf_1 + respdf_2

       FintNNLO_s_ns_oqcd(7) = respdf(1)
       kin(7) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 51                           !!
    !!-----------------------------------------------------------------------!!
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C1Lim)

    if (C1Lim%makecut) then

       FintNNLO_s_ns_oqcd(8) = zero
       kin(8) = zero

    else

       call res_tree_qqb(C1Lim%AmpMom,res_tmp)
       res_nlo(:,1) = Qdn2 * res_tmp(:,1)
       res_nlo(:,2) = Qup2 * res_tmp(:,2)

       call get_respdf(ns_lumi,1,1,C1Lim,res_nlo,respdf_bak)

       EC = C1Lim%Lim_Ei(1)
       call fill_sv_logs(C1Lim%muf(1)**2,4*EC**2,sv_logs)

       z5  = C1Lim%Lim_z(1)
       s15 = C1Lim%Lim_sij(1,5)

       !-- subtract plus, leg 1, damp = one
       eta51 = C1Lim%Lim_KinInv(4)
       intsub_1 = PqqNLO(1) - Pqq0_R(1) * log(eta51)
       intsub_1 = intsub_1 + Pqq0_R(1)*sv_logs
       intsub_1 = - intsub_1

       !-- subtract plus, leg 2, damp = zero
       intsub_2 = PqqNLO(1)
       intsub_2 = intsub_2 + Pqq0_R(1)*sv_logs
       intsub_2 = - intsub_2

       respdf_1 = (intsub_1 + intsub_2) * respdf_bak * C1Lim%wgt &
                * CF * 2/s15 * Pqg(z5)/(one-z5)

       !-- elastic component (soft photon)
       ns_int_sub = reshape(calG_QqQl(C1Lim,one,[Qdn,-Qdn, Qup,-Qup],Q_lep),[2,2])
       res_nlo(:,1) = ns_int_sub(:,1) * res_tmp(:,1)
       res_nlo(:,2) = ns_int_sub(:,2) * res_tmp(:,2)

       call get_respdf(ns_lumi,1,1,C1Lim,res_nlo,respdf_2)
       respdf_2 = (respdf_2 - 3*respdf_bak*sv_logs) * C1Lim%wgt &
                * CF * 2/s15 * Pqg(z5)/(one-z5)

       !-- total
       respdf = respdf_1 + respdf_2
       respdf = -respdf

       FintNNLO_s_ns_oqcd(8) = respdf(1)
       kin(8) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim)

    if (C2Lim%makecut) then

       FintNNLO_s_ns_oqcd(9) = zero
       kin(9) = zero

    else

       call res_tree_qqb(C2Lim%AmpMom,res_tmp)
       res_nlo(:,1) = Qdn2 * res_tmp(:,1)
       res_nlo(:,2) = Qup2 * res_tmp(:,2)

       call get_respdf(ns_lumi,1,1,C2Lim,res_nlo,respdf_bak)

       EC = C2Lim%Lim_Ei(1)
       call fill_sv_logs(C2Lim%muf(1)**2,4*EC**2,sv_logs)

       z5  = C2Lim%Lim_z(2)
       s25 = C2Lim%Lim_sij(2,5)

       !-- subtract plus, leg 1, damp = zero
       intsub_1 = PqqNLO(1)
       intsub_1 = intsub_1 + Pqq0_R(1)*sv_logs
       intsub_1 = - intsub_1

       !-- subtract plus, leg 2, damp = one
       eta52 = C2Lim%Lim_KinInv(4)
       intsub_2 = PqqNLO(1) - Pqq0_R(1) * log(eta52)
       intsub_2 = intsub_2 + Pqq0_R(1)*sv_logs
       intsub_2 = - intsub_2

       respdf_1 = (intsub_1 + intsub_2) * respdf_bak * C2Lim%wgt &
                * CF * 2/s25 * Pqg(z5)/(one-z5)

       !-- elastic component (soft photon)
       ns_int_sub = reshape(calG_QqQl(C2Lim,one,[Qdn,-Qdn, Qup,-Qup],Q_lep),[2,2])
       res_nlo(:,1) = ns_int_sub(:,1) * res_tmp(:,1)
       res_nlo(:,2) = ns_int_sub(:,2) * res_tmp(:,2)

       call get_respdf(ns_lumi,1,1,C2Lim,res_nlo,respdf_2)
       respdf_2 = (respdf_2 - 3*respdf_bak*sv_logs) * C2Lim%wgt &
                * CF * 2/s25 * Pqg(z5)/(one-z5)

       !-- total
       respdf = respdf_1 + respdf_2
       respdf = - respdf

       FintNNLO_s_ns_oqcd(9) = respdf(1)
       kin(9) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_ns_onlo(1:9) = FintNNLO_s_ns_oqcd
#endif

  end function xsect_nnlo_s_ns_oqcd


  function xsect_nnlo_s_ns_oewk_is(yRnd,ff,vegasweight)
    use mod_kinematics_nlo_z
    implicit none
    integer :: xsect_nnlo_s_ns_oewk_is
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C1Lim,C2Lim,SLim,SC1Lim,SC2Lim
    type(KinConfig) :: HardProc_z1,C1Lim_z1,C2Lim_z1,SLim_z1,SC1Lim_z1,SC2Lim_z1
    type(KinConfig) :: HardProc_z2,C1Lim_z2,C2Lim_z2,SLim_z2,SC1Lim_z2,SC2Lim_z2
    !--
    real(dp) :: xx(kNLO_max_full),z
    real(dp) :: kin(12),FintNNLO_s_ns_oewk(18)
    real(dp) :: damp,damp_qed,EC,E5,E5sq,eta51,eta52,eta61,eta62,s15,s25,z5
    real(dp) :: respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf),respdf_3(ipdf),respdf_bak(ipdf)
    real(dp) :: intsub(ipdf),intsub_1(ipdf),intsub_2(ipdf),sv_logs(ipdf)
    real(dp) :: intsub_pls(ipdf),intsub_els(ipdf)
    real(dp) :: Pqq0_R(-1:1),PqqNLO(-1:1)
    real(dp) :: res_nlo(2,2),res_tmp(2,2),eik_qed(4) !FIX ME: add "old" in the name 
    !---
    real(dp) :: res_nlo_new(-5:7,-5:7),res_tmp_new(-5:7,-5:7) !FIX ME: remove "new" in the name 
    logical  :: oldcode

    oldcode = .false.

    xsect_nnlo_s_ns_oewk_is = 0

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

    !-- Boosted kinematics - z1
    call kinematics_nlo_z_is(xx,z,one,HardProc=HardProc_z1,&
      C1Lim=C1Lim_z1,C2Lim=C2Lim_z1,SLim=SLim_z1,SC1Lim=SC1Lim_z1,SC2Lim=SC2Lim_z1,&
      compute_etas=.true.)

    !-- Boosted kinematics - z2
    call kinematics_nlo_z_is(xx,one,z,HardProc=HardProc_z2,&
      C1Lim=C1Lim_z2,C2Lim=C2Lim_z2,SLim=SLim_z2,SC1Lim=SC1Lim_z2,SC2Lim=SC2Lim_z2,&
      compute_etas=.true.)

    !-- z1=z2=1 kinematics
    call kinematics_nlo_z_is(xx,one,one,HardProc=HardProc, &
      C1Lim=C1Lim,C2Lim=C2Lim,SLim=SLim,SC1Lim=SC1Lim,SC2Lim=SC2Lim, &
      compute_etas=.true.)

    !-- Splittings and integrated sub terms
    Pqq0_R = PqqAP_0_R(z)
    PqqNLO = Pqq_NLO(z)

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !!                        FLM[z.1_q,2_qb,3,4|5_a]                        !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    ! define process specific partons
#if (_Vcharge == 0)
    HardProc_z1%ids(1:5) = [0,0,id_el,-id_el,id_a]
    C1Lim_z1%ids(1:4)    = [0,0,id_el,-id_el]
    C2Lim_z1%ids(1:4)    = [0,0,id_el,-id_el]
    SLim_z1%ids(1:4)     = [0,0,id_el,-id_el]
    SC1Lim_z1%ids(1:4)   = [0,0,id_el,-id_el]
    SC2Lim_z1%ids(1:4)   = [0,0,id_el,-id_el]

    HardProc_z1%part = [id_q,-id_q,id_el,-id_el,id_a]
    C1Lim_z1%part    = [id_q,-id_q,id_el,-id_el]
    C2Lim_z1%part    = [id_q,-id_q,id_el,-id_el]
    SLim_z1%part     = [id_q,-id_q,id_el,-id_el]
    SC1Lim_z1%part   = [id_q,-id_q,id_el,-id_el]
    SC2Lim_z1%part   = [id_q,-id_q,id_el,-id_el]

#elif  (_Vcharge == -1)
    HardProc_z1%ids(1:5) = [0,0,id_el,-id_nue,id_a]
    C1Lim_z1%ids(1:4)    = [0,0,id_el,-id_nue]
    C2Lim_z1%ids(1:4)    = [0,0,id_el,-id_nue]
    SLim_z1%ids(1:4)     = [0,0,id_el,-id_nue]
    SC1Lim_z1%ids(1:4)   = [0,0,id_el,-id_nue]
    SC2Lim_z1%ids(1:4)   = [0,0,id_el,-id_nue]

    HardProc_z1%part = [id_q,-id_qp,id_el,-id_nue,id_a]
    C1Lim_z1%part    = [id_q,-id_qp,id_el,-id_nue]
    C2Lim_z1%part    = [id_q,-id_qp,id_el,-id_nue]
    SLim_z1%part     = [id_q,-id_qp,id_el,-id_nue]
    SC1Lim_z1%part   = [id_q,-id_qp,id_el,-id_nue]
    SC2Lim_z1%part   = [id_q,-id_qp,id_el,-id_nue]

#elif  (_Vcharge == +1)
    HardProc_z1%ids(1:5) = [0,0,id_nue,-id_el,id_a]
    C1Lim_z1%ids(1:4)    = [0,0,id_nue,-id_el]
    C2Lim_z1%ids(1:4)    = [0,0,id_nue,-id_el]
    SLim_z1%ids(1:4)     = [0,0,id_nue,-id_el]
    SC1Lim_z1%ids(1:4)   = [0,0,id_nue,-id_el]
    SC2Lim_z1%ids(1:4)   = [0,0,id_nue,-id_el]

    HardProc_z1%part = [id_q,-id_qp,id_nue,-id_el,id_a]
    C1Lim_z1%part    = [id_q,-id_qp,id_nue,-id_el]
    C2Lim_z1%part    = [id_q,-id_qp,id_nue,-id_el]
    SLim_z1%part     = [id_q,-id_qp,id_nue,-id_el]
    SC1Lim_z1%part   = [id_q,-id_qp,id_nue,-id_el]
    SC2Lim_z1%part   = [id_q,-id_qp,id_nue,-id_el]
#endif


!HardProc_z1%ids(1:5) = [0,0,id_el,-id_el,id_a]


    call cut_histo(HardProc_z1)

    if (HardProc_z1%makecut) then

       FintNNLO_s_ns_oewk(1) = zero
       kin(1) = zero

    else

       call res_tree_a_qqb(HardProc_z1%AmpMom,res_nlo)

       call get_respdf(ns_lumi,1,1,HardProc_z1,res_nlo,respdf)

       call partition_nlo_qed(HardProc_z1,damp_qed,12)

       !-- partition function in the limit: 5 || 1
       HardProc_z1%Lim_etaij(:,6) = HardProc_z1%Lim_etaij(:,5)
       HardProc_z1%Lim_etaij(1,5) = zero
       HardProc_z1%Lim_etaij(2,5) = HardProc_z1%Lim_etaij(1,2)
       HardProc_z1%Lim_etaij(3,5) = HardProc_z1%Lim_etaij(1,3)
       HardProc_z1%Lim_etaij(4,5) = HardProc_z1%Lim_etaij(1,4)

       call partition_nnlo_fact(HardProc_z1,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=1,i_qed=1)

       EC = HardProc_z1%Lim_Ei(1)
       eta61 = HardProc_z1%Lim_etaij(1,6)
       call fill_sv_logs(HardProc_z1%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta61) * damp
       intsub = intsub + Pqq0_R(0)*sv_logs

       respdf = CF * intsub * respdf * HardProc_z1%wgt * damp_qed

       FintNNLO_s_ns_oewk(1) = respdf(1)
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 51                           !!
    !!-----------------------------------------------------------------------!!
    C1Lim_z1%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C1Lim_z1)

    if (C1Lim_z1%makecut) then

       FintNNLO_s_ns_oewk(2) = zero
       kin(2) = zero

    else

       call res_tree_qqb(C1Lim_z1%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(ns_lumi,1,1,C1Lim_z1,res_nlo,respdf)

       !-- damp = one
       EC = C1Lim_z1%Lim_Ei(1)
       eta61 = C1Lim_z1%Lim_KinInv(4)
       call fill_sv_logs(C1Lim_z1%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta61)
       intsub = intsub + Pqq0_R(0)*sv_logs

       s15 = C1Lim_z1%Lim_sij(1,5)
       z5  = C1Lim_z1%Lim_z(1)

       !-- damp_qed = one
       respdf = CF * intsub * respdf * C1Lim_z1%wgt &
              * 2/s15 * Pqg(z5)/(one-z5)
       respdf = -respdf

       FintNNLO_s_ns_oewk(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim_z1%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim_z1)

    if (C2Lim_z1%makecut) then

       FintNNLO_s_ns_oewk(3) = zero
       kin(3) = zero

    else

       call res_tree_qqb(C2Lim_z1%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(ns_lumi,1,1,C2Lim_z1,res_nlo,respdf)

       !-- damp = zero
       EC = C2Lim_z1%Lim_Ei(1)
       call fill_sv_logs(C2Lim_z1%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0)
       intsub = intsub + Pqq0_R(0)*sv_logs

       s25 = C2Lim_z1%Lim_sij(2,5)
       z5  = C2Lim_z1%Lim_z(2)

       !-- damp_qed = one
       respdf = CF * intsub * respdf * C2Lim_z1%wgt &
              * 2/s25 * Pqg(z5)/(one-z5)
       respdf = -respdf

       FintNNLO_s_ns_oewk(3) = respdf(1)
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                           Soft-Photon Limit                           !!
    !!-----------------------------------------------------------------------!!
    SLim_z1%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(SLim_z1)

    if (SLim_z1%makecut) then

       FintNNLO_s_ns_oewk(4:6) = zero
       kin(4) = zero

    else

       call res_tree_qqb(SLim_z1%AmpMom,res_tmp)

    !! ------------------------------- SLim -------------------------------- !!
       E5 = SLim_z1%Lim_Ei(2)
       E5sq = E5**2
       eta51 = SLim_z1%Lim_etaij(1,5)
       eta52 = SLim_z1%Lim_etaij(2,5)
       call get_qed_eik(charges,SLim_z1%Lim_etaij,[1,2,3,4],5,eik_qed)
       eik_qed = eik_qed/E5sq

       res_nlo(:,1) = eik_qed(1:2) * res_tmp(:,1)
       res_nlo(:,2) = eik_qed(3:4) * res_tmp(:,2)

       call get_respdf(ns_lumi,1,1,SLim_z1,res_nlo,respdf_1)

       call partition_nlo_qed(SLim_z1,damp_qed,12)

       !-- partition function in the limit: 5 || 1
       SLim_z1%Lim_etaij(:,6) = SLim_z1%Lim_etaij(:,5)
       SLim_z1%Lim_etaij(1,5) = zero
       SLim_z1%Lim_etaij(2,5) = SLim_z1%Lim_etaij(1,2)
       SLim_z1%Lim_etaij(3,5) = SLim_z1%Lim_etaij(1,3)
       SLim_z1%Lim_etaij(4,5) = SLim_z1%Lim_etaij(1,4)

       call partition_nnlo_fact(SLim_z1,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=1,i_qed=1)

       EC = SLim_z1%Lim_Ei(1)
       eta61 = SLim_z1%Lim_etaij(1,6)
       call fill_sv_logs(SLim_z1%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta61) * damp
       intsub = intsub + Pqq0_R(0)*sv_logs

       respdf_1 = CF * intsub * respdf_1 * SLim_z1%wgt * damp_qed
       respdf_1 = - respdf_1

       FintNNLO_s_ns_oewk(4) = respdf_1(1)

    !! ------------------------------- SC1Lim ------------------------------ !!
       res_nlo(:,1) = Qdn2 * res_tmp(:,1)
       res_nlo(:,2) = Qup2 * res_tmp(:,2)

       call get_respdf(ns_lumi,1,1,SLim_z1,res_nlo,respdf_bak)

       !-- damp = one
       eta61  = SC1Lim_z1%Lim_KinInv(4)
       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta61)
       intsub = intsub + Pqq0_R(0)*sv_logs

       !-- damp_qed = one
       respdf_2 = CF * intsub * respdf_bak * SC1Lim_z1%wgt &
                * one/eta51/E5sq

       FintNNLO_s_ns_oewk(5) = respdf_2(1)

    !! ------------------------------- SC2Lim ------------------------------ !!

       !-- damp = zero
       intsub = PqqNLO(0)
       intsub = intsub + Pqq0_R(0)*sv_logs

       !-- damp_qed = one
       respdf_3 = CF * intsub * respdf_bak * SC2Lim_z1%wgt &
                * one/eta52/E5sq

       FintNNLO_s_ns_oewk(6) = respdf_3(1)

    !! ------------------------------- Total ------------------------------- !!

       respdf = respdf_1 + respdf_2 + respdf_3
       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !!                        FLM[1_q,z.2_qb,3,4|5_a]                        !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc_z2%ids(1:5) = [0,0,id_el,-id_el,id_a]
    call cut_histo(HardProc_z2)

    if (HardProc_z2%makecut) then

       FintNNLO_s_ns_oewk(7) = zero
       kin(5) = zero

    else

       call res_tree_a_qqb(HardProc_z2%AmpMom,res_nlo)

       call get_respdf(ns_lumi,1,1,HardProc_z2,res_nlo,respdf)

       call partition_nlo_qed(HardProc_z2,damp_qed,12)

       !-- partition function in the limit: 5 || 2
       HardProc_z2%Lim_etaij(:,6) = HardProc_z2%Lim_etaij(:,5)
       HardProc_z2%Lim_etaij(1,5) = HardProc_z2%Lim_etaij(1,2)
       HardProc_z2%Lim_etaij(2,5) = zero
       HardProc_z2%Lim_etaij(3,5) = HardProc_z2%Lim_etaij(2,3)
       HardProc_z2%Lim_etaij(4,5) = HardProc_z2%Lim_etaij(2,4)

       call partition_nnlo_fact(HardProc_z2,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=2,i_qed=2)

       EC = HardProc_z2%Lim_Ei(1)
       eta62 = HardProc_z2%Lim_etaij(2,6)
       call fill_sv_logs(HardProc_z2%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta62) * damp
       intsub = intsub + Pqq0_R(0)*sv_logs

       respdf = CF * intsub * respdf * HardProc_z2%wgt * damp_qed

       FintNNLO_s_ns_oewk(7) = respdf(1)
       kin(5) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 51                           !!
    !!-----------------------------------------------------------------------!!
    C1Lim_z2%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C1Lim_z2)

    if (C1Lim_z2%makecut) then

       FintNNLO_s_ns_oewk(8) = zero
       kin(6) = zero

    else

       call res_tree_qqb(C1Lim_z2%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(ns_lumi,1,1,C1Lim_z2,res_nlo,respdf)

       !-- damp = zero
       EC = C1Lim_z2%Lim_Ei(1)
       call fill_sv_logs(C1Lim_z2%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0)
       intsub = intsub + Pqq0_R(0)*sv_logs

       s15 = C1Lim_z2%Lim_sij(1,5)
       z5  = C1Lim_z2%Lim_z(1)

       !-- damp_qed = one
       respdf = CF * intsub * respdf * C1Lim_z2%wgt &
              * 2/s15 * Pqg(z5)/(one-z5)
       respdf = -respdf

       FintNNLO_s_ns_oewk(8) = respdf(1)
       kin(6) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim_z2%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim_z2)

    if (C2Lim_z2%makecut) then

       FintNNLO_s_ns_oewk(9) = zero
       kin(7) = zero

    else

       call res_tree_qqb(C2Lim_z2%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(ns_lumi,1,1,C2Lim_z2,res_nlo,respdf)

       !-- damp = one
       EC = C2Lim_z2%Lim_Ei(1)
       eta62 = C2Lim_z2%Lim_KinInv(4)
       call fill_sv_logs(C2Lim_z2%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta62)
       intsub = intsub + Pqq0_R(0)*sv_logs

       s25 = C2Lim_z2%Lim_sij(2,5)
       z5  = C2Lim_z2%Lim_z(2)

       !-- damp_qed = one
       respdf = CF * intsub * respdf * C2Lim_z2%wgt &
              * 2/s25 * Pqg(z5)/(one-z5)
       respdf = -respdf

       FintNNLO_s_ns_oewk(9) = respdf(1)
       kin(7) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                           Soft-Photon Limit                           !!
    !!-----------------------------------------------------------------------!!
    SLim_z2%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(SLim_z2)

    if (SLim_z2%makecut) then

       FintNNLO_s_ns_oewk(10:12) = zero
       kin(8) = zero

    else

       call res_tree_qqb(SLim_z2%AmpMom,res_tmp)

    !! ------------------------------- SLim -------------------------------- !!
       E5 = SLim_z2%Lim_Ei(2)
       E5sq = E5**2
       eta51 = SLim_z2%Lim_etaij(1,5)
       eta52 = SLim_z2%Lim_etaij(2,5)
       call get_qed_eik(charges,SLim_z2%Lim_etaij,[1,2,3,4],5,eik_qed)
       eik_qed = eik_qed/E5sq

       res_nlo(:,1) = eik_qed(1:2) * res_tmp(:,1)
       res_nlo(:,2) = eik_qed(3:4) * res_tmp(:,2)

       call get_respdf(ns_lumi,1,1,SLim_z2,res_nlo,respdf_1)

       call partition_nlo_qed(SLim_z2,damp_qed,12)

       !-- partition function in the limit: 5 || 2
       SLim_z2%Lim_etaij(:,6) = SLim_z2%Lim_etaij(:,5)
       SLim_z2%Lim_etaij(1,5) = SLim_z2%Lim_etaij(1,2)
       SLim_z2%Lim_etaij(2,5) = zero
       SLim_z2%Lim_etaij(3,5) = SLim_z2%Lim_etaij(2,3)
       SLim_z2%Lim_etaij(4,5) = SLim_z2%Lim_etaij(2,4)

       call partition_nnlo_fact(SLim_z2,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=2,i_qed=2)

       EC = SLim_z2%Lim_Ei(1)
       eta62 = SLim_z2%Lim_etaij(2,6)
       call fill_sv_logs(SLim_z2%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta62) * damp
       intsub = intsub + Pqq0_R(0)*sv_logs

       respdf_1 = CF * intsub * respdf_1 * SLim_z2%wgt * damp_qed
       respdf_1 = - respdf_1

       FintNNLO_s_ns_oewk(10) = respdf_1(1)

    !! ------------------------------- SC1Lim ------------------------------ !!
       res_nlo(:,1) = Qdn2 * res_tmp(:,1)
       res_nlo(:,2) = Qup2 * res_tmp(:,2)

       call get_respdf(ns_lumi,1,1,SLim_z2,res_nlo,respdf_bak)

       !-- damp = zero
       intsub = PqqNLO(0)
       intsub = intsub + Pqq0_R(0)*sv_logs

       !-- damp_qed = one
       respdf_2 = CF * intsub * respdf_bak * SC1Lim_z2%wgt &
                * one/eta51/E5sq

       FintNNLO_s_ns_oewk(11) = respdf_2(1)

    !! ------------------------------- SC2Lim ------------------------------ !!

       !-- damp = one
       eta62  = SC2Lim_z2%Lim_KinInv(4)
       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta62)
       intsub = intsub + Pqq0_R(0)*sv_logs

       !-- damp_qed = one
       respdf_3 = CF * intsub * respdf_bak * SC2Lim_z2%wgt &
                * one/eta52/E5sq

       FintNNLO_s_ns_oewk(12) = respdf_3(1)

    !! ------------------------------- Total ------------------------------- !!

       respdf = respdf_1 + respdf_2 + respdf_3
       kin(8) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !!                         FLM[1_q,2_qb,3,4|5_a]                         !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    ! define process specific partons
#if (_Vcharge == 0)
    HardProc_z1%ids(1:5) = [0,0,id_el,-id_el,id_a]
    C1Lim_z1%ids(1:4)    = [0,0,id_el,-id_el]
    C2Lim_z1%ids(1:4)    = [0,0,id_el,-id_el]
    SLim_z1%ids(1:4)     = [0,0,id_el,-id_el]
    SC1Lim_z1%ids(1:4)   = [0,0,id_el,-id_el]
    SC2Lim_z1%ids(1:4)   = [0,0,id_el,-id_el]

    HardProc_z1%part = [id_q,-id_q,id_el,-id_el,id_a]
    C1Lim_z1%part    = [id_q,-id_q,id_el,-id_el]
    C2Lim_z1%part    = [id_q,-id_q,id_el,-id_el]
    SLim_z1%part     = [id_q,-id_q,id_el,-id_el]
    SC1Lim_z1%part   = [id_q,-id_q,id_el,-id_el]
    SC2Lim_z1%part   = [id_q,-id_q,id_el,-id_el]

#elif  (_Vcharge == -1)
    HardProc_z1%ids(1:5) = [0,0,id_el,-id_nue,id_a]
    C1Lim_z1%ids(1:4)    = [0,0,id_el,-id_nue]
    C2Lim_z1%ids(1:4)    = [0,0,id_el,-id_nue]
    SLim_z1%ids(1:4)     = [0,0,id_el,-id_nue]
    SC1Lim_z1%ids(1:4)   = [0,0,id_el,-id_nue]
    SC2Lim_z1%ids(1:4)   = [0,0,id_el,-id_nue]

    HardProc_z1%part = [id_q,-id_qp,id_el,-id_nue,id_a]
    C1Lim_z1%part    = [id_q,-id_qp,id_el,-id_nue]
    C2Lim_z1%part    = [id_q,-id_qp,id_el,-id_nue]
    SLim_z1%part     = [id_q,-id_qp,id_el,-id_nue]
    SC1Lim_z1%part   = [id_q,-id_qp,id_el,-id_nue]
    SC2Lim_z1%part   = [id_q,-id_qp,id_el,-id_nue]

#elif  (_Vcharge == +1)
    HardProc_z1%ids(1:5) = [0,0,id_nue,-id_el,id_a]
    C1Lim_z1%ids(1:4)    = [0,0,id_nue,-id_el]
    C2Lim_z1%ids(1:4)    = [0,0,id_nue,-id_el]
    SLim_z1%ids(1:4)     = [0,0,id_nue,-id_el]
    SC1Lim_z1%ids(1:4)   = [0,0,id_nue,-id_el]
    SC2Lim_z1%ids(1:4)   = [0,0,id_nue,-id_el]

    HardProc_z1%part = [id_q,-id_qp,id_nue,-id_el,id_a]
    C1Lim_z1%part    = [id_q,-id_qp,id_nue,-id_el]
    C2Lim_z1%part    = [id_q,-id_qp,id_nue,-id_el]
    SLim_z1%part     = [id_q,-id_qp,id_nue,-id_el]
    SC1Lim_z1%part   = [id_q,-id_qp,id_nue,-id_el]
    SC2Lim_z1%part   = [id_q,-id_qp,id_nue,-id_el]
#endif

      oldcode = .false.

!HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]

    call cut_histo(HardProc)

    if (HardProc%makecut) then

       FintNNLO_s_ns_oewk(13) = zero
       kin(9) = zero

    else

       if (oldcode ) then
          call res_tree_a_qqb(HardProc%AmpMom,res_nlo)

       call get_respdf(ns_lumi,1,1,HardProc,res_nlo,respdf)

       EC = HardProc%Lim_Ei(1)
       call fill_sv_logs(HardProc%muf(1)**2,4*EC**2,sv_logs)

       call partition_nlo_qed(HardProc,damp_qed,12)

       HardProc%Lim_etaij(:,6) = HardProc%Lim_etaij(:,5)

       !-- subtract plus, leg 1
       HardProc%Lim_etaij(1,5) = zero
       HardProc%Lim_etaij(2,5) = HardProc%Lim_etaij(1,2)
       HardProc%Lim_etaij(3,5) = HardProc%Lim_etaij(1,3)
       HardProc%Lim_etaij(4,5) = HardProc%Lim_etaij(1,4)

       call partition_nnlo_fact(HardProc,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=1,i_qed=1)
       eta61 = HardProc%Lim_etaij(1,6)
       intsub_1 = PqqNLO(1) - Pqq0_R(1) * log(eta61) * damp
       intsub_1 = intsub_1 + Pqq0_R(1)*sv_logs
       intsub_1 = - intsub_1

       !-- subtract plus, leg 2
       HardProc%Lim_etaij(1,5) = HardProc%Lim_etaij(1,2)
       HardProc%Lim_etaij(2,5) = zero
       HardProc%Lim_etaij(3,5) = HardProc%Lim_etaij(2,3)
       HardProc%Lim_etaij(4,5) = HardProc%Lim_etaij(2,4)

       call partition_nnlo_fact(HardProc,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=2,i_qed=2)
       eta62 = HardProc%Lim_etaij(2,6)
       intsub_2 = PqqNLO(1) - Pqq0_R(1) * log(eta62) * damp
       intsub_2 = intsub_2 + Pqq0_R(1)*sv_logs
       intsub_2 = - intsub_2

       intsub_pls = intsub_1 + intsub_2

       !-- elastic component
       intsub_els = 4*zeta2 - 3*sv_logs

       !-- total
       respdf = CF * (intsub_pls + intsub_els) * respdf * HardProc%wgt * damp_qed

       FintNNLO_s_ns_oewk(13) = respdf(1)
       kin(9) = respdf(1)

       call fill_histo(respdf,vegasweight)

       else 
          print*, 'not implemented yet'
          call res_tree_a_qqb_gen(HardProc%AmpMom,res_nlo_new)
          stop
       endif 


    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 51                           !!
    !!-----------------------------------------------------------------------!!
    C1Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C1Lim)

    if (C1Lim%makecut) then

       FintNNLO_s_ns_oewk(14) = zero
       kin(10) = zero

    else

       call res_tree_qqb(C1Lim%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(ns_lumi,1,1,C1Lim,res_nlo,respdf)

       EC = C1Lim%Lim_Ei(1)
       call fill_sv_logs(C1Lim%muf(1)**2,4*EC**2,sv_logs)

       z5  = C1Lim%Lim_z(1)
       s15 = C1Lim%Lim_sij(1,5)

       !-- subtract plus, leg 1, damp = one
       eta61 = C1Lim%Lim_KinInv(4)
       intsub_1 = PqqNLO(1) - Pqq0_R(1) * log(eta61)
       intsub_1 = intsub_1 + Pqq0_R(1)*sv_logs
       intsub_1 = - intsub_1

       !-- subtract plus, leg 2, damp = zero
       intsub_2 = PqqNLO(1)
       intsub_2 = intsub_2 + Pqq0_R(1)*sv_logs
       intsub_2 = - intsub_2

       intsub_pls = intsub_1 + intsub_2

       !-- elastic component
       intsub_els = 4*zeta2 - 3*sv_logs

       !-- total, damp_qed = one
       respdf = CF * (intsub_pls + intsub_els) * respdf * C1Lim%wgt &
              * 2/s15 * Pqg(z5)/(one-z5)
       respdf = -respdf

       FintNNLO_s_ns_oewk(14) = respdf(1)
       kin(10) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 52                           !!
    !!-----------------------------------------------------------------------!!
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(C2Lim)

    if (C2Lim%makecut) then

       FintNNLO_s_ns_oewk(15) = zero
       kin(11) = zero

    else

       call res_tree_qqb(C2Lim%AmpMom,res_nlo)
       res_nlo(:,1) = Qdn2 * res_nlo(:,1)
       res_nlo(:,2) = Qup2 * res_nlo(:,2)

       call get_respdf(ns_lumi,1,1,C2Lim,res_nlo,respdf)

       EC = C2Lim%Lim_Ei(1)
       call fill_sv_logs(C2Lim%muf(1)**2,4*EC**2,sv_logs)

       z5  = C2Lim%Lim_z(2)
       s25 = C2Lim%Lim_sij(2,5)

       !-- subtract plus, leg 1, damp = zero
       intsub_1 = PqqNLO(1)
       intsub_1 = intsub_1 + Pqq0_R(1)*sv_logs
       intsub_1 = - intsub_1

       !-- subtract plus, leg 2, damp = one
       eta62 = C2Lim%Lim_KinInv(4)
       intsub_2 = PqqNLO(1) - Pqq0_R(1) * log(eta62)
       intsub_2 = intsub_2 + Pqq0_R(1)*sv_logs
       intsub_2 = - intsub_2

       intsub_pls = intsub_1 + intsub_2

       !-- elastic component
       intsub_els = 4*zeta2 - 3*sv_logs

       !-- total
       respdf = CF * (intsub_pls + intsub_els) * respdf * C2Lim%wgt &
              * 2/s25 * Pqg(z5)/(one-z5)
       respdf = - respdf

       FintNNLO_s_ns_oewk(15) = respdf(1)
       kin(11) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                           Soft-Photon Limit                           !!
    !!-----------------------------------------------------------------------!!
    SLim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(SLim)

    if (SLim%makecut) then

       FintNNLO_s_ns_oewk(16:18) = zero
       kin(12) = zero

    else

       call res_tree_qqb(SLim%AmpMom,res_tmp)

    !! ------------------------------- SLim -------------------------------- !!
       E5 = SLim%Lim_Ei(2)
       E5sq = E5**2
       eta51 = SLim%Lim_etaij(1,5)
       eta52 = SLim%Lim_etaij(2,5)
       call get_qed_eik(charges,SLim%Lim_etaij,[1,2,3,4],5,eik_qed)
       eik_qed = eik_qed/E5sq

       res_nlo(:,1) = eik_qed(1:2) * res_tmp(:,1)
       res_nlo(:,2) = eik_qed(3:4) * res_tmp(:,2)

       call get_respdf(ns_lumi,1,1,SLim,res_nlo,respdf_1)

       call partition_nlo_qed(SLim,damp_qed,12)

       SLim%Lim_etaij(:,6) = SLim%Lim_etaij(:,5)

       EC = SLim%Lim_Ei(1)
       call fill_sv_logs(SLim%muf(1)**2,4*EC**2,sv_logs)

       !-- subtract plus, leg 1
       SLim%Lim_etaij(1,5) = zero
       SLim%Lim_etaij(2,5) = SLim%Lim_etaij(1,2)
       SLim%Lim_etaij(3,5) = SLim%Lim_etaij(1,3)
       SLim%Lim_etaij(4,5) = SLim%Lim_etaij(1,4)

       call partition_nnlo_fact(SLim,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=1,i_qed=1)
       eta61 = SLim%Lim_etaij(1,6)
       intsub_1 = PqqNLO(1) - Pqq0_R(1) * log(eta61) * damp
       intsub_1 = intsub_1 + Pqq0_R(1)*sv_logs
       intsub_1 = - intsub_1

       !-- subtract plus, leg 2
       SLim%Lim_etaij(1,5) = SLim%Lim_etaij(1,2)
       SLim%Lim_etaij(2,5) = zero
       SLim%Lim_etaij(3,5) = SLim%Lim_etaij(2,3)
       SLim%Lim_etaij(4,5) = SLim%Lim_etaij(2,4)

       call partition_nnlo_fact(SLim,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=2,i_qed=2)
       eta62 = SLim%Lim_etaij(2,6)
       intsub_2 = PqqNLO(1) - Pqq0_R(1) * log(eta62) * damp
       intsub_2 = intsub_2 + Pqq0_R(1)*sv_logs
       intsub_2 = - intsub_2

       intsub_pls = intsub_1 + intsub_2

       !-- elastic component
       intsub_els = 4*zeta2 - 3*sv_logs

       !-- total
       respdf_1 = CF * (intsub_pls + intsub_els) * respdf_1 * SLim%wgt * damp_qed
       respdf_1 = - respdf_1

       FintNNLO_s_ns_oewk(16) = respdf_1(1)

    !! ------------------------------- SC1Lim ------------------------------ !!
       res_nlo(:,1) = Qdn2 * res_tmp(:,1)
       res_nlo(:,2) = Qup2 * res_tmp(:,2)

       call get_respdf(ns_lumi,1,1,SLim,res_nlo,respdf_bak)

       !-- subtract plus, leg 1
       !-- damp = one
       eta61  = SC1Lim%Lim_KinInv(4)
       intsub_1 = PqqNLO(1) - Pqq0_R(1) * log(eta61)
       intsub_1 = intsub_1 + Pqq0_R(1)*sv_logs
       intsub_1 = - intsub_1

       !-- subtract plus, leg 2
       !-- damp = zero
       intsub_2 = PqqNLO(1)
       intsub_2 = intsub_2 + Pqq0_R(1)*sv_logs
       intsub_2 = - intsub_2

       intsub_pls = intsub_1 + intsub_2

       !-- damp_qed = one
       respdf_2 = CF * (intsub_pls + intsub_els) * respdf_bak * SC1Lim%wgt &
                * one/eta51/E5sq

       FintNNLO_s_ns_oewk(17) = respdf_2(1)

    !! ------------------------------- SC2Lim ------------------------------ !!

       !-- subtract plus, leg 1
       !-- damp = zero
       intsub_1 = PqqNLO(1)
       intsub_1 = intsub_1 + Pqq0_R(1)*sv_logs
       intsub_1 = - intsub_1

       !-- subtract plus, leg 2
       !-- damp = one
       eta62  = SC2Lim%Lim_KinInv(4)
       intsub_2 = PqqNLO(1) - Pqq0_R(1) * log(eta62)
       intsub_2 = intsub_2 + Pqq0_R(1)*sv_logs
       intsub_2 = - intsub_2

       intsub_pls = intsub_1 + intsub_2

       !-- damp_qed = one
       respdf_3 = CF * (intsub_pls + intsub_els) * respdf_bak * SC2Lim%wgt &
                * one/eta52/E5sq

       FintNNLO_s_ns_oewk(18) = respdf_3(1)

    !! ------------------------------- Total ------------------------------- !!

       respdf = respdf_1 + respdf_2 + respdf_3
       kin(12) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_ns_onlo = FintNNLO_s_ns_oewk
#endif

  end function xsect_nnlo_s_ns_oewk_is

  function xsect_nnlo_s_ns_oewk_fs_5i(yRnd,ff,vegasweight,icoll,jother)
    use mod_kinematics_nlo_z
    implicit none
    integer :: xsect_nnlo_s_ns_oewk_fs_5i,icoll,jother
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,CLim,SLim,SCLim
    type(KinConfig) :: HardProc_z1,CLim_z1,SLim_z1,SCLim_z1
    type(KinConfig) :: HardProc_z2,CLim_z2,SLim_z2,SCLim_z2
    !--
    real(dp) :: xx(kNLO_max_full),z
    real(dp) :: kin(9),FintNNLO_s_ns_oewk(12)
    real(dp) :: damp,damp_qed,EC,E5,E5sq,eta5i,eta61,eta62,si5,z5
    real(dp) :: respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf)
    real(dp) :: intsub(ipdf),intsub_1(ipdf),intsub_2(ipdf),sv_logs(ipdf)
    real(dp) :: intsub_pls(ipdf),intsub_els(ipdf)
    real(dp) :: Pqq0_R(-1:1),PqqNLO(-1:1)
    real(dp) :: res_nlo(2,2),res_lo(2,2),eik_qed(4)
    real(dp) :: res_tmp_vect(2,2,2),respdf_vect(2,ipdf)

    xsect_nnlo_s_ns_oewk_fs_5i = 0

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

    !-- Boosted kinematics - z1
    call kinematics_nlo_z_fs(xx,icoll,jother,z,one,HardProc=HardProc_z1,&
      CLim=CLim_z1,SLim=SLim_z1,SCLim=SCLim_z1)

    !-- Boosted kinematics - z2
    call kinematics_nlo_z_fs(xx,icoll,jother,one,z,HardProc=HardProc_z2,&
      CLim=CLim_z2,SLim=SLim_z2,SCLim=SCLim_z2)

    !-- z1=z2=1 kinematics
    call kinematics_nlo_z_fs(xx,icoll,jother,one,one,HardProc=HardProc, &
      CLim=CLim,SLim=SLim,SCLim=SCLim)

    !-- Splittings and integrated sub terms
    Pqq0_R = PqqAP_0_R(z)
    PqqNLO = Pqq_NLO(z)

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !!                        FLM[z.1_q,2_qb,3,4|5_a]                        !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc_z1%ids(1:5) = [0,0,id_el,-id_el,id_a]
    call cut_histo(HardProc_z1)

    if (HardProc_z1%makecut) then

       FintNNLO_s_ns_oewk(1) = zero
       kin(1) = zero

    else

       call res_tree_a_qqb(HardProc_z1%AmpMom,res_nlo)

       call get_respdf(ns_lumi,1,1,HardProc_z1,res_nlo,respdf)

       call partition_nlo_qed(HardProc_z1,damp_qed,icoll)

       !-- partition function in the limit: 5 || 1
       HardProc_z1%Lim_etaij(:,6) = HardProc_z1%Lim_etaij(:,5)
       HardProc_z1%Lim_etaij(1,5) = zero
       HardProc_z1%Lim_etaij(2,5) = HardProc_z1%Lim_etaij(1,2)
       HardProc_z1%Lim_etaij(3,5) = HardProc_z1%Lim_etaij(1,3)
       HardProc_z1%Lim_etaij(4,5) = HardProc_z1%Lim_etaij(1,4)

       call partition_nnlo_fact(HardProc_z1,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=1,i_qed=1)

       EC = HardProc_z1%Lim_Ei(1)
       eta61 = HardProc_z1%Lim_etaij(1,6)
       call fill_sv_logs(HardProc_z1%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta61) * damp
       intsub = CF * (intsub + Pqq0_R(0)*sv_logs)

       respdf = intsub * respdf * HardProc_z1%wgt * damp_qed

       FintNNLO_s_ns_oewk(1) = respdf(1)
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 5i                           !!
    !!-----------------------------------------------------------------------!!
    CLim_z1%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(CLim_z1)

    if (CLim_z1%makecut) then

       FintNNLO_s_ns_oewk(2) = zero
       kin(2) = zero

    else

       call res_tree_qqb(CLim_z1%AmpMom,res_nlo)
       call get_respdf(ns_lumi,1,1,CLim_z1,res_nlo,respdf)

       !-- damp = zero
       EC = CLim_z1%Lim_Ei(1)
       call fill_sv_logs(CLim_z1%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0)
       intsub = CF * (intsub + Pqq0_R(0)*sv_logs)

       si5 = CLim_z1%Lim_sij(icoll,5)
       z5  = CLim_z1%Lim_z(1)

       !-- damp_qed = one
       respdf = intsub * respdf * CLim_z1%wgt &
              * Q_lep2 * 2/si5 * Pqg(z5)
       respdf = - respdf

       FintNNLO_s_ns_oewk(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                           Soft-Photon Limit                           !!
    !!-----------------------------------------------------------------------!!
    SLim_z1%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(SLim_z1)

    if (SLim_z1%makecut) then

       FintNNLO_s_ns_oewk(3:4) = zero
       kin(3) = zero

    else

       call res_tree_qqb(SLim_z1%AmpMom,res_lo)

       E5 = SLim_z1%Lim_Ei(2)
       E5sq = E5**2
       eta5i = SLim_z1%Lim_etaij(icoll,5)
       call get_qed_eik(charges,SLim_z1%Lim_etaij,[1,2,3,4],5,eik_qed)
       eik_qed = eik_qed/E5sq

       res_tmp_vect(:,1,1) = eik_qed(1:2) * res_lo(:,1)
       res_tmp_vect(:,2,1) = eik_qed(3:4) * res_lo(:,2)

       !-- prepare PDFs structures, CS
       res_tmp_vect(:,:,2) = Q_lep2 * res_lo(:,:)

       call get_respdf_vect(ns_lumi,1,1,SLim_z1,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))

    !! ------------------------------- SLim -------------------------------- !!

       call partition_nlo_qed(SLim_z1,damp_qed,icoll)

       !-- partition function in the limit: 5 || 1
       SLim_z1%Lim_etaij(:,6) = SLim_z1%Lim_etaij(:,5)
       SLim_z1%Lim_etaij(1,5) = zero
       SLim_z1%Lim_etaij(2,5) = SLim_z1%Lim_etaij(1,2)
       SLim_z1%Lim_etaij(3,5) = SLim_z1%Lim_etaij(1,3)
       SLim_z1%Lim_etaij(4,5) = SLim_z1%Lim_etaij(1,4)

       call partition_nnlo_fact(SLim_z1,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=1,i_qed=1)

       EC = SLim_z1%Lim_Ei(1)
       eta61 = SLim_z1%Lim_etaij(1,6)
       call fill_sv_logs(SLim_z1%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta61) * damp
       intsub = CF * (intsub + Pqq0_R(0)*sv_logs)

       respdf_1 = intsub * respdf_vect(1,:) * SLim_z1%wgt * damp_qed
       respdf_1 = - respdf_1

       FintNNLO_s_ns_oewk(3) = respdf_1(1)

    !! ------------------------------- SCLim ------------------------------ !!

       !-- damp = zero
       intsub = PqqNLO(0)
       intsub = CF * (intsub + Pqq0_R(0)*sv_logs)

       !-- damp_qed = one
       respdf_2 = intsub * respdf_vect(2,:) * SCLim_z1%wgt &
                * one/eta5i/E5sq

       FintNNLO_s_ns_oewk(4) = respdf_2(1)

    !! ------------------------------- Total ------------------------------- !!

       respdf = respdf_1 + respdf_2
       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !!                        FLM[1_q,z.2_qb,3,4|5_a]                        !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc_z2%ids(1:5) = [0,0,id_el,-id_el,id_a]
    call cut_histo(HardProc_z2)

    if (HardProc_z2%makecut) then

       FintNNLO_s_ns_oewk(5) = zero
       kin(4) = zero

    else

       call res_tree_a_qqb(HardProc_z2%AmpMom,res_nlo)

       call get_respdf(ns_lumi,1,1,HardProc_z2,res_nlo,respdf)

       call partition_nlo_qed(HardProc_z2,damp_qed,icoll)

       !-- partition function in the limit: 5 || 2
       HardProc_z2%Lim_etaij(:,6) = HardProc_z2%Lim_etaij(:,5)
       HardProc_z2%Lim_etaij(1,5) = HardProc_z2%Lim_etaij(1,2)
       HardProc_z2%Lim_etaij(2,5) = zero
       HardProc_z2%Lim_etaij(3,5) = HardProc_z2%Lim_etaij(2,3)
       HardProc_z2%Lim_etaij(4,5) = HardProc_z2%Lim_etaij(2,4)

       call partition_nnlo_fact(HardProc_z2,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=2,i_qed=2)

       EC = HardProc_z2%Lim_Ei(1)
       eta62 = HardProc_z2%Lim_etaij(2,6)
       call fill_sv_logs(HardProc_z2%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta62) * damp
       intsub = CF * (intsub + Pqq0_R(0)*sv_logs)

       respdf = intsub * respdf * HardProc_z2%wgt * damp_qed

       FintNNLO_s_ns_oewk(5) = respdf(1)
       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 5i                           !!
    !!-----------------------------------------------------------------------!!
    CLim_z2%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(CLim_z2)

    if (CLim_z2%makecut) then

       FintNNLO_s_ns_oewk(6) = zero
       kin(5) = zero

    else

       call res_tree_qqb(CLim_z2%AmpMom,res_nlo)
       call get_respdf(ns_lumi,1,1,CLim_z2,res_nlo,respdf)

       !-- damp = zero
       EC = CLim_z2%Lim_Ei(1)
       call fill_sv_logs(CLim_z2%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0)
       intsub = CF * (intsub + Pqq0_R(0)*sv_logs)

       si5 = CLim_z2%Lim_sij(icoll,5)
       z5  = CLim_z2%Lim_z(1)

       !-- damp_qed = one
       respdf = intsub * respdf * CLim_z2%wgt &
              * Q_lep2 * 2/si5 * Pqg(z5)
       respdf = - respdf

       FintNNLO_s_ns_oewk(6) = respdf(1)
       kin(5) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                           Soft-Photon Limit                           !!
    !!-----------------------------------------------------------------------!!
    SLim_z2%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(SLim_z2)

    if (SLim_z2%makecut) then

       FintNNLO_s_ns_oewk(7:8) = zero
       kin(6) = zero

    else

       call res_tree_qqb(SLim_z2%AmpMom,res_lo)

       E5 = SLim_z2%Lim_Ei(2)
       E5sq = E5**2
       eta5i = SLim_z2%Lim_etaij(icoll,5)
       call get_qed_eik(charges,SLim_z2%Lim_etaij,[1,2,3,4],5,eik_qed)
       eik_qed = eik_qed/E5sq

       res_tmp_vect(:,1,1) = eik_qed(1:2) * res_lo(:,1)
       res_tmp_vect(:,2,1) = eik_qed(3:4) * res_lo(:,2)

       !-- prepare PDFs structures, CS
       res_tmp_vect(:,:,2) = Q_lep2 * res_lo(:,:)

       call get_respdf_vect(ns_lumi,1,1,SLim_z2,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))

    !! ------------------------------- SLim -------------------------------- !!

       call partition_nlo_qed(SLim_z2,damp_qed,icoll)

       !-- partition function in the limit: 5 || 1
       SLim_z2%Lim_etaij(:,6) = SLim_z2%Lim_etaij(:,5)
       SLim_z2%Lim_etaij(1,5) = SLim_z2%Lim_etaij(1,2)
       SLim_z2%Lim_etaij(2,5) = zero
       SLim_z2%Lim_etaij(3,5) = SLim_z2%Lim_etaij(2,3)
       SLim_z2%Lim_etaij(4,5) = SLim_z2%Lim_etaij(2,4)

       call partition_nnlo_fact(SLim_z2,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=2,i_qed=2)

       EC = SLim_z2%Lim_Ei(1)
       eta62 = SLim_z2%Lim_etaij(2,6)
       call fill_sv_logs(SLim_z2%muf(1)**2,4*EC**2,sv_logs)

       intsub = PqqNLO(0) - Pqq0_R(0) * log(eta62) * damp
       intsub = CF * (intsub + Pqq0_R(0)*sv_logs)

       respdf_1 = intsub * respdf_vect(1,:) * SLim_z2%wgt * damp_qed
       respdf_1 = - respdf_1

       FintNNLO_s_ns_oewk(7) = respdf_1(1)

    !! ------------------------------- SCLim ------------------------------ !!

       !-- damp = zero
       intsub = PqqNLO(0)
       intsub = CF * (intsub + Pqq0_R(0)*sv_logs)

       !-- damp_qed = one
       respdf_2 = intsub * respdf_vect(2,:) * SCLim_z2%wgt &
                * one/eta5i/E5sq

       FintNNLO_s_ns_oewk(8) = respdf_2(1)

    !! ------------------------------- Total ------------------------------- !!

       respdf = respdf_1 + respdf_2
       kin(6) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
    !!                         FLM[1_q,2_qb,3,4|5_a]                         !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_a]
    call cut_histo(HardProc)

    if (HardProc%makecut) then

       FintNNLO_s_ns_oewk(9) = zero
       kin(7) = zero

    else

       call res_tree_a_qqb(HardProc%AmpMom,res_nlo)

       call get_respdf(ns_lumi,1,1,HardProc,res_nlo,respdf)

       EC = HardProc%Lim_Ei(1)
       call fill_sv_logs(HardProc%muf(1)**2,4*EC**2,sv_logs)

       call partition_nlo_qed(HardProc,damp_qed,icoll)

       HardProc%Lim_etaij(:,6) = HardProc%Lim_etaij(:,5)

       !-- subtract plus, leg 1
       HardProc%Lim_etaij(1,5) = zero
       HardProc%Lim_etaij(2,5) = HardProc%Lim_etaij(1,2)
       HardProc%Lim_etaij(3,5) = HardProc%Lim_etaij(1,3)
       HardProc%Lim_etaij(4,5) = HardProc%Lim_etaij(1,4)

       call partition_nnlo_fact(HardProc,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=1,i_qed=1)
       eta61 = HardProc%Lim_etaij(1,6)
       intsub_1 = PqqNLO(1) - Pqq0_R(1) * log(eta61) * damp
       intsub_1 = CF * (intsub_1 + Pqq0_R(1)*sv_logs)
       intsub_1 = - intsub_1

       !-- subtract plus, leg 2
       HardProc%Lim_etaij(1,5) = HardProc%Lim_etaij(1,2)
       HardProc%Lim_etaij(2,5) = zero
       HardProc%Lim_etaij(3,5) = HardProc%Lim_etaij(2,3)
       HardProc%Lim_etaij(4,5) = HardProc%Lim_etaij(2,4)

       call partition_nnlo_fact(HardProc,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=2,i_qed=2)
       eta62 = HardProc%Lim_etaij(2,6)
       intsub_2 = PqqNLO(1) - Pqq0_R(1) * log(eta62) * damp
       intsub_2 = CF * (intsub_2 + Pqq0_R(1)*sv_logs)
       intsub_2 = - intsub_2

       intsub_pls = intsub_1 + intsub_2

       !-- elastic component
       intsub_els = CF * (4*zeta2 - 3*sv_logs)

       !-- total
       respdf = (intsub_pls + intsub_els) * respdf * HardProc%wgt * damp_qed

       FintNNLO_s_ns_oewk(9) = respdf(1)
       kin(7) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          Collinear Limit 5i                           !!
    !!-----------------------------------------------------------------------!!
    CLim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(CLim)

    if (CLim%makecut) then

       FintNNLO_s_ns_oewk(10) = zero
       kin(8) = zero

    else

       call res_tree_qqb(CLim%AmpMom,res_nlo)

       call get_respdf(ns_lumi,1,1,CLim,res_nlo,respdf)

       EC = CLim%Lim_Ei(1)
       call fill_sv_logs(CLim%muf(1)**2,4*EC**2,sv_logs)

       si5 = CLim%Lim_sij(icoll,5)
       z5  = CLim%Lim_z(1)

       !-- subtract plus, leg 1, damp = zero
       intsub_1 = PqqNLO(1)
       intsub_1 = CF * (intsub_1 + Pqq0_R(1)*sv_logs)
       intsub_1 = - intsub_1

       !-- subtract plus, leg 2, damp = zero
       intsub_2 = PqqNLO(1)
       intsub_2 = CF * (intsub_2 + Pqq0_R(1)*sv_logs)
       intsub_2 = - intsub_2

       intsub_pls = intsub_1 + intsub_2

       !-- elastic component
       intsub_els = CF * (4*zeta2 - 3*sv_logs)

       !-- total, damp_qed = one
       respdf = (intsub_pls + intsub_els) * respdf * CLim%wgt &
              * Q_lep2 * 2/si5 * Pqg(z5)
       respdf = -respdf

       FintNNLO_s_ns_oewk(10) = respdf(1)
       kin(8) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !!-----------------------------------------------------------------------!!
    !!                           Soft-Photon Limit                           !!
    !!-----------------------------------------------------------------------!!
    SLim%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(SLim)

    if (SLim%makecut) then

       FintNNLO_s_ns_oewk(11:12) = zero
       kin(9) = zero

    else

       call res_tree_qqb(SLim%AmpMom,res_lo)

       SLim%Lim_etaij(:,6) = SLim%Lim_etaij(:,5)

       E5 = SLim%Lim_Ei(2)
       E5sq = E5**2
       eta5i = SLim%Lim_etaij(icoll,5)
       call get_qed_eik(charges,SLim%Lim_etaij,[1,2,3,4],5,eik_qed)
       eik_qed = eik_qed/E5sq

       res_tmp_vect(:,1,1) = eik_qed(1:2) * res_lo(:,1)
       res_tmp_vect(:,2,1) = eik_qed(3:4) * res_lo(:,2)

       !-- prepare PDFs structures, CS
       res_tmp_vect(:,:,2) = Q_lep2 * res_lo(:,:)

       call get_respdf_vect(ns_lumi,1,1,SLim,res_tmp_vect(:,:,1:2),respdf_vect(1:2,:))

       EC = SLim%Lim_Ei(1)
       call fill_sv_logs(SLim%muf(1)**2,4*EC**2,sv_logs)

    !! ------------------------------- SLim -------------------------------- !!

       call partition_nlo_qed(SLim,damp_qed,icoll)

       !-- subtract plus, leg 1
       SLim%Lim_etaij(1,5) = zero
       SLim%Lim_etaij(2,5) = SLim%Lim_etaij(1,2)
       SLim%Lim_etaij(3,5) = SLim%Lim_etaij(1,3)
       SLim%Lim_etaij(4,5) = SLim%Lim_etaij(1,4)
       call partition_nnlo_fact(SLim,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=1,i_qed=1)

       eta61 = SLim%Lim_etaij(1,6)
       intsub_1 = PqqNLO(1) - Pqq0_R(1) * log(eta61) * damp
       intsub_1 = CF * (intsub_1 + Pqq0_R(1)*sv_logs)
       intsub_1 = - intsub_1

       !-- subtract plus, leg 2
       SLim%Lim_etaij(1,5) = SLim%Lim_etaij(1,2)
       SLim%Lim_etaij(2,5) = zero
       SLim%Lim_etaij(3,5) = SLim%Lim_etaij(2,3)
       SLim%Lim_etaij(4,5) = SLim%Lim_etaij(2,4)

       call partition_nnlo_fact(SLim,damp,&
            iconf_qcd=[1,2,5],iconf_qed=[1,2,3,4,6],i_qcd=2,i_qed=2)

       eta62 = SLim%Lim_etaij(2,6)
       intsub_2 = PqqNLO(1) - Pqq0_R(1) * log(eta62) * damp
       intsub_2 = CF * (intsub_2 + Pqq0_R(1)*sv_logs)
       intsub_2 = - intsub_2

       intsub_pls = intsub_1 + intsub_2

       !-- elastic component
       intsub_els = CF * (4*zeta2 - 3*sv_logs)

       !-- total
       respdf_1 = (intsub_pls + intsub_els) * respdf_vect(1,:) * SLim%wgt * damp_qed
       respdf_1 = - respdf_1

       FintNNLO_s_ns_oewk(11) = respdf_1(1)

    !! ------------------------------- SCLim ------------------------------ !!

       !-- damp = zero --> intsub_1 = intsub_2
       intsub = PqqNLO(1)
       intsub = CF * (intsub + Pqq0_R(1)*sv_logs)
       intsub = - intsub

       intsub_pls = 2 * intsub

       !-- damp_qed = one
       respdf_2 = (intsub_pls + intsub_els) * respdf_vect(2,:) * SCLim%wgt &
                * one/eta5i/E5sq

       FintNNLO_s_ns_oewk(12) = respdf_2(1)

    !! ------------------------------- Total ------------------------------- !!

       respdf = respdf_1 + respdf_2
       kin(9) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif


    !!-----------------------------------------------------------------------!!

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_ns_onlo(1:12) = FintNNLO_s_ns_oewk
#endif

  end function xsect_nnlo_s_ns_oewk_fs_5i

  function xsect_nnlo_s_ns_vewknf(yRnd,ff,vegasweight)
    use mod_kinematics_lo
    implicit none
    integer :: xsect_nnlo_s_ns_vewknf
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp) :: xx(kLO_max_full)
    real(dp) :: respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf)
    real(dp) :: res_nlo(2,2),kin(1)

    xsect_nnlo_s_ns_vewknf = 0

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

       call res_ewkloop_qqb_nf(LOProc%AmpMom,res_nlo)

       !-- F[z1,2]
       call get_respdf_hoppet(xPij,PDFs,ns_lumi,1,1,LOProc,res_nlo,respdf_1,myPDFs1_Lmu=[xPij_Lmu])

       !-- F[1,z2]
       call get_respdf_hoppet(PDFs,xPij,ns_lumi,1,1,LOProc,res_nlo,respdf_2,myPDFs2_Lmu=[xPij_Lmu])

       !-- Total
       respdf = CF * (respdf_1 + respdf_2) * LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)

#if(_withchecks == 1)
    FintNNLO_s_ns_qqx(1) = kin(1)
#endif

  end function xsect_nnlo_s_ns_vewknf

end module mod_xsects_nnlo_s_ns

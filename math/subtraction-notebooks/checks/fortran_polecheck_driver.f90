!-- Values of the compiled Fortran pieces transcribed in W_polecheck_fortran_defs.wl; built and run by run_verify_W_polecheck.sh,
!-- read by verify_W_polecheck_compiled_fortran.wls. Kernels: private copies of mod_hoppet_nnlo (kern_copy), divided by x.
program polecheck_driver
  use kern_copy
  use mod_parms, only: Q_IS
  use mod_process, only: KinConfig, initialize_config
  use mod_int_sub_nnlo, only: calG_gq_cc, calG_ONLOQCD_gq, calG_CF
  use mod_splittings_bare, only: Pqq_NLO, PqqAP_0_R, Pgq_spav
  use mod_limits_rv, only: col_rvewk_is_qg_gen, col_rvqcd_is_qa
  implicit none
  integer, parameter :: wp = kind(1d0)
  real(wp), parameter :: EC = 7.1_wp, E3 = 4.7_wp, E4 = 5.3_wp, E5 = 2.3_wp, eta13 = 0.29_wp, eta14 = 0.61_wp, eta35 = 0.43_wp, eta45 = 0.37_wp
  real(wp), parameter :: x2 = 0.27_wp, zz = 0.3_wp
  real(wp) :: xs(5) = [0.11_wp, 0.3_wp, 0.5_wp, 0.7_wp, 0.9_wp], lmu(2) = [0.37_wp, -1.21_wp]
  real(wp) :: amp(-5:7,-5:7), r(-5:7,-5:7,2), q(-1:1), x, y, rv(4), rq(2), s5i
  type(KinConfig) :: proc
  integer :: i, c, al(4), be(4), il(4)
  real(wp) :: ql(4)
  !-- (leg-1, leg-2) flavours of the reduced Born: W+ dbar u, u dbar (ilept 4, Qlept +1); W- d ubar, ubar d (ilept 3, -1)
  al = [-1,2,1,-2]; be = [2,-1,-2,1]; il = [4,4,3,3]; ql = [1,1,-1,-1]

  open(10, file='fortran_polecheck_values.m')
  write(10,'(a)') 'fv = <||>;'
  !-- kernels (regular parts) and splittings
  cc_piece = cc_REAL
  do i = 1, 5
     x = xs(i); y = -log(x)
     call wr('gq',mySub_gq_func(y)/x);          call wr('gq_Lmu',mySub_gq_Lmu_func(y)/x)
     call wr('gq_b',mySub_gq_b_func(y)/x);      call wr('gq_b_Lmu',mySub_gq_b_Lmu_func(y)/x)
     call wr('TCgq',mySub_TCgq_func(y)/x);      call wr('TCgq_Lmu',mySub_TCgq_Lmu_func(y)/x)
     call wr('TCgq_Lmu2',mySub_TCgq_Lmu2_func(y)/x)
     call wr('ns_qq',mySub_ns_qq_func(y)/x);    call wr('ns_qq_Lmu',mySub_ns_qq_Lmu_func(y)/x)
     call wr('ns_qqb',mySub_ns_qqb_func(y)/x);  call wr('ns_qqb_Lmu',mySub_ns_qqb_Lmu_func(y)/x)
     q = Pqq_NLO(x);   call wr('PA0',q(0));  call wr('PA1',q(1))
     q = PqqAP_0_R(x); call wr('PlA0',q(0)); call wr('PlA1',q(1))
     call wr('Pgq_spav',Pgq_spav(x))
  enddo
  write(10,'(a,4(es26.17e3,","),es26.17e3,a)') 'fv["xs"] = {', xs, '};'
  write(10,'(a,1(es26.17e3,","),es26.17e3,a)') 'fv["lmu"] = {', lmu, '};'

  !-- Born configuration in the partonic c.m. frame with the reduced leg z p1 (s_12 boosted, kinematics_lo_z)
  call initialize_config(proc)
  proc%Lim_Ei(1:4) = [EC, EC, E3, E4]
  proc%Lim_etaij(1,2) = 1
  proc%Lim_etaij(1,3) = eta13; proc%Lim_etaij(2,3) = 1-eta13; proc%Lim_etaij(1,4) = eta14; proc%Lim_etaij(2,4) = 1-eta14
  proc%AmpMom(:,1) = zz*EC*[1._wp,0._wp,0._wp, 1._wp]
  proc%AmpMom(:,2) =    EC*[1._wp,0._wp,0._wp,-1._wp]
  proc%AmpMom(:,3) = E3*[1._wp, sqrt(1-(1-2*eta13)**2), 0._wp, 1-2*eta13]
  proc%AmpMom(:,4) = E4*[1._wp,-sqrt(1-(1-2*eta14)**2), 0._wp, 1-2*eta14]
  do c = 1, 4
     amp = 0; amp(al(c),be(c)) = 1
     write(10,'(a,i0,a,es26.17e3,a,es26.17e3,a,i0,a,es26.17e3,a)') 'fv[{"chg",',c,'}] = {', Q_IS(al(c)), ',', Q_IS(be(c)), ',', il(c), ',', ql(c), '};'
     r = calG_gq_cc(proc,amp,lmu,zz,ql(c),il(c),1)
     write(10,'(a,i0,a,es26.17e3,a,es26.17e3,a)') 'fv[{"calG_gq_cc",',c,'}] = {', r(al(c),be(c),1), ',', r(al(c),be(c),2), '};'
  enddo

  !-- hard configuration of s_oqcd (kinematics_nlo_z_is): Lim_Ei(2) = E5, eta_15 = x2, eta_25 = 1 - x2
  proc%Lim_Ei(2) = E5
  proc%Lim_etaij(1,5) = x2; proc%Lim_etaij(2,5) = 1-x2
  proc%Lim_etaij(3,5) = eta35; proc%Lim_etaij(4,5) = eta45
  proc%Lim_etaij(5,:) = proc%Lim_etaij(:,5)
  do c = 1, 4
     amp = 0; amp(0,be(c)) = 1
     r = calG_ONLOQCD_gq(proc,amp,lmu,ql(c),il(c),1,2)
     write(10,'(a,i0,a,es26.17e3,a,es26.17e3,a)') 'fv[{"calG_ONLOQCD_gq",',c,'}] = {', r(0,be(c),1), ',', r(0,be(c),2), '};'
  enddo
  write(10,'(a,7(es26.17e3,","),es26.17e3,a)') 'fv["kin"] = {', EC, E3, E4, E5, eta13, eta14, eta35, eta45, '};'
  write(10,'(a,1(es26.17e3,","),es26.17e3,a)') 'fv["x2z"] = {', x2, zz, '};'

  !-- calG_CF_hard(EC,Ek,etaik,etajk,damp1,damp2) of s_aq_oewk
  write(10,'(a,es26.17e3,a)') 'fv["calG_CF"] = ', calG_CF(EC,E5,0.31_wp,0.69_wp,0.41_wp,0.59_wp), ';'

  !-- collinear limits of the real-virtual sectors: Lim_KinInv(1:3) = [zF, s5i, eta15] (kinematics_nlo_is), reduced momenta
  s5i = 1.7_wp
  proc%Lim_KinInv(1:3) = [0.35_wp, s5i, 0.13_wp]
  rv = col_rvewk_is_qg_gen(1,proc)*s5i/2
  rq = col_rvqcd_is_qa(proc)*s5i/2
  write(10,'(a,3(es26.17e3,","),es26.17e3,a)') 'fv["col_rvewk"] = {', rv, '};'
  write(10,'(a,1(es26.17e3,","),es26.17e3,a)') 'fv["col_rvqcd"] = {', rq, '};'
  write(10,'(a,2(es26.17e3,","),es26.17e3,a)') 'fv["colkin"] = {', proc%Lim_KinInv(1), proc%Lim_KinInv(3), zz, '};'
  close(10)
contains
  subroutine wr(nm, v)
    character(*), intent(in) :: nm
    real(wp), intent(in) :: v
    write(10,'(3a,i0,a,es26.17e3,a)') 'fv[{"',nm,'",',i,'}] = ', v, ';'
  end subroutine wr
end program polecheck_driver

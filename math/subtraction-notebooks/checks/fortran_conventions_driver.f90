!-- Fortran values for verify_fortran_conventions.wls; built and run by run_verify_fortran_conventions.sh
program driver
  use kern_copy
  use mod_parms, only: Q_IS
  use mod_process, only: KinConfig, initialize_config
  use mod_int_sub_nnlo, only: fin_ns_vqcd_cc, genGnloQCD_cc, legfinalcalG_cc, fin_ns_v7_cc
  use mod_splittings_bare, only: Pqq_NLO, PqqAP_0_R
  implicit none
  integer, parameter :: wp = kind(1d0)
  real(wp) :: xs(4) = [0.07_wp, 0.31_wp, 0.62_wp, 0.93_wp], lmu(2) = [0.37_wp, -1.21_wp]
  real(wp) :: amp(-5:7,-5:7), r(-5:7,-5:7,2), q(-1:1), ql(6), y
  type(KinConfig) :: proc
  integer :: i, ip, k, c, pieces(3), al(6), be(6), il(6)
  !-- flavour pairs: W+ u db, db u, c sb (ilept 4, +1); W- d ub, ub d (ilept 3, -1); u ub with Qlept +1 (charges not conserved)
  al = [2,-1,4,1,-2,2]; be = [-1,2,-3,-2,1,-2]; il = [4,4,4,3,3,4]; ql = [1,1,1,-1,-1,1]
  pieces = [cc_REAL, cc_VIRT, cc_DELTA]
  open(10, file='fortran_values.m')
  write(10,'(a)') 'fv = <||>;'
  do i = 1, 4
     y = -log(xs(i))
     do ip = 1, 3
        cc_piece = pieces(ip)
        call wr('ns',mySub_ns_func(y));  call wr('nsL',mySub_ns_Lmu_func(y))
        call wr('P',PqqNLO_func(y));     call wr('Pl',PqqNLO_Lmu_func(y))
        call wr('Pb',PqqNLO_b_func(y));  call wr('Pbl',PqqNLO_b_Lmu_func(y))
        call wr('TC0',mySub_TCqq_func(y)); call wr('TC1',mySub_TCqq_Lmu_func(y)); call wr('TC2',mySub_TCqq_Lmu2_func(y))
     enddo
     q = Pqq_NLO(xs(i))
     do k = -1, 1
        write(10,'(a,i0,a,i0,a,es26.17e3,a)') 'fv[{"PA",',i,',',k,'}] = ', q(k), ';'
     enddo
     q = PqqAP_0_R(xs(i))
     do k = -1, 1
        write(10,'(a,i0,a,i0,a,es26.17e3,a)') 'fv[{"PlA",',i,',',k,'}] = ', q(k), ';'
     enddo
  enddo
  !-- elastic routines: EC = 45.3, E3 = 38.1, E4 = 51.7, angles of the leptons
  call initialize_config(proc)
  proc%Lim_Ei(1) = 45.3_wp; proc%Lim_Ei(2) = 45.3_wp; proc%Lim_Ei(3) = 38.1_wp; proc%Lim_Ei(4) = 51.7_wp
  proc%Lim_etaij(1,2) = 1; proc%Lim_etaij(1,3) = 0.27_wp; proc%Lim_etaij(2,3) = 0.73_wp
  proc%Lim_etaij(1,4) = 0.64_wp; proc%Lim_etaij(2,4) = 0.36_wp; proc%Lim_etaij(3,4) = 0.81_wp
  do c = 1, 6
     amp = 0; amp(al(c),be(c)) = 1
     write(10,'(a,i0,a,es26.17e3,a,es26.17e3,a,i0,a,es26.17e3,a)') 'fv[{"chg",',c,'}] = {', Q_IS(al(c)), ',', Q_IS(be(c)), ',', il(c), ',', ql(c), '};'
     r = fin_ns_vqcd_cc(proc,amp,lmu,ql(c),il(c));        call wre('vqcd',c)
     r = genGnloQCD_cc(proc,amp,lmu,ql(c),il(c));         call wre('genG',c)
     r = legfinalcalG_cc(proc,amp,lmu,ql(c),il(c),1);     call wre('legG1',c)
     r = legfinalcalG_cc(proc,amp,lmu,ql(c),il(c),2);     call wre('legG2',c)
     r = fin_ns_v7_cc(proc,amp,lmu,ql(c),il(c));          call wre('v7',c)
  enddo
  close(10)
contains
  subroutine wr(nm, v)
    character(*), intent(in) :: nm
    real(wp), intent(in) :: v
    write(10,'(3a,i0,a,i0,a,es26.17e3,a)') 'fv[{"',nm,'",',i,',',ip,'}] = ', v, ';'
  end subroutine wr
  subroutine wre(nm, cc)
    character(*), intent(in) :: nm
    integer, intent(in) :: cc
    write(10,'(3a,i0,a,es26.17e3,a,es26.17e3,a)') 'fv[{"',nm,'",',cc,'}] = {', r(al(cc),be(cc),1), ',', r(al(cc),be(cc),2), '};'
  end subroutine wre
end program driver

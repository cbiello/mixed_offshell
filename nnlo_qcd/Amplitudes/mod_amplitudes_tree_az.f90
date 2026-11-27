!-- amplitudes with extra jet: factor out gs**2 per extra emission

module mod_amplitudes_tree_wp
  use mod_types; use mod_consts_dp
  use mod_auxfunctions
  implicit none
  private

  public :: restree, restree_j
  public :: restree_2q_jj, restree_gg_jj, restree_qg_jj, restree_gq_jj, restree_qqb_jj

contains

  !----------------------------------------------------------------------------!
  !-- LO                                                                     --!
  !----------------------------------------------------------------------------!

  !-- f(1) + f(2) --> [V --> l(3) + lb(4)]
  !-- summed / averaged
  subroutine restree(p,res)
    real(dp), intent(in)  :: p(4,4)
    real(dp), intent(out) :: res(2)
    real(dp)              :: propsq_w,prefac

    propsq_w = one/((two*scr(p(:,3),p(:,4))-mwsq)**2+mwsq*gaw**2)
    prefac = gwsq**2*propsq_w*xn

    !--u dx > ve ep, c sx > ve ep
    res(1) = four*scr(p(:,2),p(:,3))*scr(p(:,1),p(:,4))

    !--dx u > ve ep, sx c > ve ep MIRROR
    res(2) = four*scr(p(:,2),p(:,4))*scr(p(:,1),p(:,3))

    res = res  * aveqq * prefac

    return

  end subroutine restree

  !----------------------------------------------------------------------------!
  !-- NLO                                                                    --!
  !----------------------------------------------------------------------------!

  !-- f(1) + f(2) --> f(3) [V --> l(4) + lb(5)]
  !-- summed / averaged
  subroutine restree_j(p,res_sping_j,res_j)
    real(dp), intent(in) :: p(4,5)
    complex(dp), intent(out) :: res_sping_j(6,-1:1,-1:1)
    real(dp), intent(out) :: res_j(6)
    real(dp)              :: prefac,propsq_w
    real(dp) :: sprod(5,5)
    complex(dp) :: za(5,5),zb(5,5)
    integer :: part_conf(3,6),i,ii,jj
    complex(dp) :: amp(-1:1)
    real(dp), parameter :: avecol_qq = aveqq * half * xn * Cf 
    real(dp), parameter :: avecol_qg = aveqg * half * xn * Cf 

    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)

    propsq_w = one/((two*scr(p(:,4),p(:,5))-mwsq)**2+mwsq*gaw**2)
    prefac = propsq_w * gwsq**2 * four


    !-- amplitude is for 0 -> q(1) g(2) qb(3) [V --> l(4) lb(5)]                                                                                                                                             
    part_conf(1:3,1) = (/2,3,1/) !-- q qx -> g                                                                                                                                                               
    part_conf(1:3,2) = (/3,1,2/) !-- g q -> q                                                                                                                                                                
    part_conf(1:3,3) = (/2,1,3/) !-- g qx -> qx                                                                                                                                                              
    !                                                                                                                                                                                                        
    part_conf(1:3,4) = (/1,3,2/) !-- qx q -> g                                                                                                                                                               
    part_conf(1:3,5) = (/3,2,1/) !-- q g -> q                                                                                                                                                                
    part_conf(1:3,6) = (/1,2,3/) !-- qx g -> qx                                                                                                                                                              


    res_j = zero
    res_sping_j = zero

    do i = 1,6

       call master_amp_qgqbV_w(part_conf(1,i),part_conf(2,i),part_conf(3,i),4,5,za,zb,sprod,amp)

       do ii = -1,1,2
       do jj = -1,1,2

          res_sping_j(i,ii,jj) = res_sping_j(i,ii,jj) + amp(ii)*conjg(amp(jj))

       enddo
       enddo

    enddo


    res_sping_j = res_sping_j * prefac

    !--u dx > g ve ep, c sx > g ve ep
    res_sping_j(1,:,:) = res_sping_j(1,:,:) * avecol_qq

    !--g u > d ve ep, g c > s ve ep
    res_sping_j(2,:,:) = res_sping_j(2,:,:) * avecol_qg

    !--g dx > ux ve ep, g sx > cx ve ep
    res_sping_j(3,:,:) = res_sping_j(3,:,:) * avecol_qg

    !--dx u > g ve ep, sx c > g ve ep MIRROR
    res_sping_j(4,:,:) = res_sping_j(4,:,:) * avecol_qq

    !--u g > d ve ep, c g > s ve ep MIRROR
    res_sping_j(5,:,:) = res_sping_j(5,:,:) * avecol_qg

    !--dx g > ux ve ep, sx g > cx ve ep MIRROR
    res_sping_j(6,:,:) = res_sping_j(6,:,:) * avecol_qg

    res_j(:) = res_sping_j(:,+1,+1)+res_sping_j(:,-1,-1)


    return

  end subroutine restree_j

  !-- amplitude for 0 -> q(1) g(2) qb(3) [V --> l(4) lb(5)]
  !-- index: hq,hg,hl,down/up
  subroutine master_amp_qgqbV_wp(coupl,i1,i2,i3,i4,i5,za,zb,sprod,amp)
    integer, intent(in) :: i1,i2,i3,i4,i5
    complex(dp), intent(in) :: za(5,5),zb(5,5),coupl(-1:1,-1:1,1:2)
    real(dp), intent(in) :: sprod(5,5)
    complex(dp), intent(out) :: amp(-1:1,-1:1,-1:1,1:2)
    
    amp(-1,+1,-1,:) = two * za(i1,i4)**2 * zb(i5,i4)/(za(i1,i2)*za(i2,i3)) * coupl(-1,-1,:) / sprod(i4,i5)
    amp(-1,-1,-1,:) = two * za(i4,i5) * zb(i5,i3)**2/(zb(i1,i2)*zb(i2,i3)) * coupl(-1,-1,:) / sprod(i4,i5)
    
    amp(-1,+1,+1,:) = two * za(i1,i5)**2 * zb(i4,i5)/(za(i1,i2)*za(i2,i3)) * coupl(-1,+1,:) / sprod(i4,i5)
    amp(-1,-1,+1,:) = two * za(i5,i4) * zb(i4,i3)**2/(zb(i1,i2)*zb(i2,i3)) * coupl(-1,+1,:) / sprod(i4,i5)
    
    amp(+1,+1,-1,:) = two * za(i3,i4)**2 * zb(i5,i4)/(za(i1,i2)*za(i2,i3)) * coupl(+1,-1,:) / sprod(i4,i5)
    amp(+1,-1,-1,:) = two * za(i4,i5) * zb(i5,i1)**2/(zb(i1,i2)*zb(i2,i3)) * coupl(+1,-1,:) / sprod(i4,i5)
    
    amp(+1,+1,+1,:) = two * za(i3,i5)**2 * zb(i4,i5)/(za(i1,i2)*za(i2,i3)) * coupl(+1,+1,:) / sprod(i4,i5)
    amp(+1,-1,+1,:) = two * za(i5,i4) * zb(i4,i1)**2/(zb(i1,i2)*zb(i2,i3)) * coupl(+1,+1,:) / sprod(i4,i5)

    return

  end subroutine master_amp_qgqbV_wp


  !----------------------------------------------------------------------------!
  !-- NNLO                                                                   --!
  !----------------------------------------------------------------------------!

  !-- f(1) + f(2) -> f(3) + f(4) + [V-->l(5)+lb(6)]
  !-- summed/averaged
  !-- a...d notation: see vN
  !-- all individual partonic channels checked against madgraph
  subroutine restree_2q_jj(p,res_jj)
    real(dp), intent(in) :: p(4,6)
    real(dp), intent(out) :: res_jj(36)
    real(dp) :: sprod(6,6)
    complex(dp) :: za(6,6),zb(6,6)
    complex(dp) :: coupl(-1:1,-1:1,1:2)
    integer :: part_conf(4,6,4) !-- particle ordering, a-f amplitude, (qqb,qbq,qq,qbqb)
    complex(dp) :: amp_a(-1:1,-1:1,-1:1,1:2) !-- hel index, down/up
    complex(dp) :: amp_b(-1:1,-1:1,-1:1,1:2)
    complex(dp) :: amp_c(-1:1,-1:1,-1:1,1:2)
    complex(dp) :: amp_d(-1:1,-1:1,-1:1,1:2)
    complex(dp) :: amp_e(-1:1,-1:1,-1:1,1:2)
    complex(dp) :: amp_f(-1:1,-1:1,-1:1,1:2)
    real(dp) :: res_a(1:2,1:2) !-- down/up, qqb/qbq
    real(dp) :: res_b(1:2,1:2)
    real(dp) :: res_ab(1:4,1:2) !-- down/up/down_in up_out/up_in down_out, qqb/qqb
    real(dp) :: res_ac(1:2,1:2),res_ad(1:2,1:2)
    real(dp) :: res_bc(1:2,1:2),res_bd(1:2,1:2)
    real(dp) :: res_c(1:2,1:4) !-- down/up, qqb/qbq/qq/qbqb
    real(dp) :: res_d(1:2,1:4),res_e(1:2,1:4),res_f(1:2,1:4) 
    real(dp) :: res_cd(1:4,1:4),res_ce(1:2,1:4),res_cf(1:2,1:4)
    real(dp) :: res_de(1:2,1:4),res_df(1:2,1:4)
    real(dp) :: res_ef(1:4,1:4)
    real(dp), parameter :: pref1 = xn*four*Cf*tr*aveqq
    real(dp), parameter :: pref2 = xn*two*Cf*(Ca-two*Cf)*aveqq
    integer :: i1,i3,i5,i,j

    call spinoru(6,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5),p(:,6)/),za,zb,sprod)

    call get_coupl_wp(sprod(5,6),coupl)

    part_conf(1:4,1,1) = (/2,1,3,4/) !-- a,qqb
    part_conf(1:4,2,1) = (/3,4,2,1/) !-- b,qqb
    part_conf(1:4,3,1) = (/3,1,2,4/) !-- c,qqb
    part_conf(1:4,4,1) = (/2,4,3,1/) !-- d,qqb

    part_conf(1:4,1,2) = (/1,2,3,4/) !-- a,qbq
    part_conf(1:4,2,2) = (/3,4,1,2/) !-- b,qbq
    part_conf(1:4,3,2) = (/3,2,1,4/) !-- c,qbq
    part_conf(1:4,4,2) = (/1,4,3,2/) !-- d,qbq

    part_conf(1:4,3,3) = (/3,1,4,2/) !-- c,qq
    part_conf(1:4,4,3) = (/4,2,3,1/) !-- d,qq
    part_conf(1:4,5,3) = (/3,2,4,1/) !-- e,qq
    part_conf(1:4,6,3) = (/4,1,3,2/) !-- f,qq

    part_conf(1:4,3,4) = (/1,3,2,4/) !-- c,qbqb
    part_conf(1:4,4,4) = (/2,4,1,3/) !-- d,qbqb
    part_conf(1:4,5,4) = (/2,3,1,4/) !-- e,qbqb
    part_conf(1:4,6,4) = (/1,4,2,3/) !-- f,qbqb

    res_a = zero
    res_b = zero
    res_c = zero
    res_d = zero
    res_e = zero
    res_f = zero

    res_ab = zero
    res_ac = zero
    res_ad = zero

    res_bc = zero
    res_bd = zero
    
    res_cd = zero
    res_ce = zero
    res_cf = zero
    
    res_de = zero
    res_df = zero

    res_ef = zero

    !-- qqb + qbq

    do i = 1,2
    
       call master_amp_qqbQQBV_wp(coupl,part_conf(1,1,i),part_conf(2,1,i),part_conf(3,1,i),part_conf(4,1,i),5,6,za,zb,sprod,amp_a)
       call master_amp_qqbQQBV_wp(coupl,part_conf(1,2,i),part_conf(2,2,i),part_conf(3,2,i),part_conf(4,2,i),5,6,za,zb,sprod,amp_b)
       call master_amp_qqbQQBV_wp(coupl,part_conf(1,3,i),part_conf(2,3,i),part_conf(3,3,i),part_conf(4,3,i),5,6,za,zb,sprod,amp_c)
       call master_amp_qqbQQBV_wp(coupl,part_conf(1,4,i),part_conf(2,4,i),part_conf(3,4,i),part_conf(4,4,i),5,6,za,zb,sprod,amp_d)
       
       do i5 = -1,1,2
       do i3 = -1,1,2
       do i1 = -1,1,2
                
          res_a(:,i) = res_a(:,i) + abs(amp_a(i1,i3,i5,:))**2
          res_b(:,i) = res_b(:,i) + abs(amp_b(i1,i3,i5,:))**2
          res_c(:,i) = res_c(:,i) + abs(amp_c(i1,i3,i5,:))**2
          res_d(:,i) = res_d(:,i) + abs(amp_d(i1,i3,i5,:))**2
          
          res_ab(1:2,i) = res_ab(1:2,i) + two*real(amp_a(i1,i3,i5,:)*conjg(amp_b(i3,i1,i5,:)),kind=dp)
          res_ab(3,i) = res_ab(3,i) + two*real(amp_a(i1,i3,i5,1)*conjg(amp_b(i3,i1,i5,2)),kind=dp)
          res_ab(4,i) = res_ab(4,i) + two*real(amp_a(i1,i3,i5,2)*conjg(amp_b(i3,i1,i5,1)),kind=dp)
          
          res_cd(1:2,i) = res_cd(1:2,i) + two*real(amp_c(i1,i3,i5,:)*conjg(amp_d(i3,i1,i5,:)),kind=dp)
          res_cd(3,i) = res_cd(3,i) + two*real(amp_c(i1,i3,i5,1)*conjg(amp_d(i3,i1,i5,2)),kind=dp)
          res_cd(4,i) = res_cd(4,i) + two*real(amp_c(i1,i3,i5,2)*conjg(amp_d(i3,i1,i5,1)),kind=dp)
          
       enddo
       
       res_ac(:,i) = res_ac(:,i) + two*real(amp_a(i3,i3,i5,:)*conjg(amp_c(i3,i3,i5,:)),kind=dp)
       res_ad(:,i) = res_ad(:,i) + two*real(amp_a(i3,i3,i5,:)*conjg(amp_d(i3,i3,i5,:)),kind=dp)
       
       res_bc(:,i) = res_bc(:,i) + two*real(amp_b(i3,i3,i5,:)*conjg(amp_c(i3,i3,i5,:)),kind=dp)
       res_bd(:,i) = res_bd(:,i) + two*real(amp_b(i3,i3,i5,:)*conjg(amp_d(i3,i3,i5,:)),kind=dp)
       
       enddo
       enddo
    
    enddo

    !-- qq + qbqb

    do i = 3,4
    
       call master_amp_qqbQQBV_wp(coupl,part_conf(1,3,i),part_conf(2,3,i),part_conf(3,3,i),part_conf(4,3,i),5,6,za,zb,sprod,amp_c)
       call master_amp_qqbQQBV_wp(coupl,part_conf(1,4,i),part_conf(2,4,i),part_conf(3,4,i),part_conf(4,4,i),5,6,za,zb,sprod,amp_d)
       call master_amp_qqbQQBV_wp(coupl,part_conf(1,5,i),part_conf(2,5,i),part_conf(3,5,i),part_conf(4,5,i),5,6,za,zb,sprod,amp_e)
       call master_amp_qqbQQBV_wp(coupl,part_conf(1,6,i),part_conf(2,6,i),part_conf(3,6,i),part_conf(4,6,i),5,6,za,zb,sprod,amp_f)
       
       do i5 = -1,1,2
       do i3 = -1,1,2
       do i1 = -1,1,2
                
          res_c(:,i) = res_c(:,i) + abs(amp_c(i1,i3,i5,:))**2
          res_d(:,i) = res_d(:,i) + abs(amp_d(i1,i3,i5,:))**2
          res_e(:,i) = res_e(:,i) + abs(amp_e(i1,i3,i5,:))**2
          res_f(:,i) = res_f(:,i) + abs(amp_f(i1,i3,i5,:))**2
          
          res_cd(1:2,i) = res_cd(1:2,i) + two*real(amp_c(i1,i3,i5,:)*conjg(amp_d(i3,i1,i5,:)),kind=dp)
          res_cd(3,i) = res_cd(3,i) + two*real(amp_c(i1,i3,i5,1)*conjg(amp_d(i3,i1,i5,2)),kind=dp)
          res_cd(4,i) = res_cd(4,i) + two*real(amp_c(i1,i3,i5,2)*conjg(amp_d(i3,i1,i5,1)),kind=dp)
          
          res_ef(1:2,i) = res_ef(1:2,i) + two*real(amp_e(i1,i3,i5,:)*conjg(amp_f(i3,i1,i5,:)),kind=dp)
          res_ef(3,i) = res_ef(3,i) + two*real(amp_e(i1,i3,i5,1)*conjg(amp_f(i3,i1,i5,2)),kind=dp)
          res_ef(4,i) = res_ef(4,i) + two*real(amp_e(i1,i3,i5,2)*conjg(amp_f(i3,i1,i5,1)),kind=dp)
          
       enddo
       
       res_ce(:,i) = res_ce(:,i) + two*real(amp_c(i3,i3,i5,:)*conjg(amp_e(i3,i3,i5,:)),kind=dp)
       res_cf(:,i) = res_cf(:,i) + two*real(amp_c(i3,i3,i5,:)*conjg(amp_f(i3,i3,i5,:)),kind=dp)
       
       res_de(:,i) = res_de(:,i) + two*real(amp_d(i3,i3,i5,:)*conjg(amp_e(i3,i3,i5,:)),kind=dp)
       res_df(:,i) = res_df(:,i) + two*real(amp_d(i3,i3,i5,:)*conjg(amp_f(i3,i3,i5,:)),kind=dp)
       
       enddo
       enddo
    
    enddo

    !---

    !-- color and symmetry factors
    res_a = res_a * pref1 * vn_taga
    res_b = res_b * pref1 * vn_tagb
    res_c = res_c * pref1 * vn_tagc
    res_d = res_d * pref1 * vn_tagc 
    res_e = res_e * pref1 * vn_tagc
    res_f = res_f * pref1 * vn_tagc
    
    res_ab = res_ab * pref1 * vn_tagab
    res_cd = res_cd * pref1 * vn_tagcd
    res_ef = res_ef * pref1 * vn_tagcd
    
    res_ac = res_ac * pref2 * vn_tagac
    res_ad = res_ad * pref2 * vn_tagac
    res_bc = res_bc * pref2 * vn_tagbc
    res_bd = res_bd * pref2 * vn_tagbc
    res_ce = res_ce * pref2 * vn_tagcf !-- swap of e and f wrt vN
    res_cf = res_cf * pref2 * vn_tagce !-- swap of e and f wrt vN
    res_de = res_de * pref2 * vn_tagce !-- swap of e and f wrt vN
    res_df = res_df * pref2 * vn_tagcf !-- swap of e and f wrt vN

    !-- fill processes

    !--u u > u u em ep, c c > c c em ep
    res_jj(1) = (res_c(2,3)+res_d(2,3)+res_e(2,3)+res_f(2,3)+res_cd(2,3)+res_ce(2,3)+&
         res_cf(2,3)+res_de(2,3)+res_df(2,3)+res_ef(2,3))*half !-- symm factor

    !--u ux > u ux em ep, c cx > c cx em ep
    res_jj(2) = (res_a(2,1)+res_b(2,1)+res_c(2,1)+res_d(2,1)+&
         res_ab(2,1)+res_ac(2,1)+res_ad(2,1)+res_bc(2,1)+res_bd(2,1)+res_cd(2,1))


    !--d d > d d em ep, s s > s s em ep
    res_jj(3) = (res_c(1,3)+res_d(1,3)+res_e(1,3)+res_f(1,3)+res_cd(1,3)+res_ce(1,3)+&
         res_cf(1,3)+res_de(1,3)+res_df(1,3)+res_ef(1,3))*half !-- symm factor

    !--d dx > d dx em ep, s sx > s sx em ep
    res_jj(4) = (res_a(1,1)+res_b(1,1)+res_c(1,1)+res_d(1,1)+&
         res_ab(1,1)+res_ac(1,1)+res_ad(1,1)+res_bc(1,1)+res_bd(1,1)+res_cd(1,1))

    !--ux ux > ux ux em ep, cx cx > cx cx em ep
    res_jj(5) = (res_c(2,4)+res_d(2,4)+res_e(2,4)+res_f(2,4)+res_cd(2,4)+res_ce(2,4)+&
         res_cf(2,4)+res_de(2,4)+res_df(2,4)+res_ef(2,4))*half !-- symm factor

    !--dx dx > dx dx em ep, sx sx > sx sx em ep
    res_jj(6) = (res_c(1,4)+res_d(1,4)+res_e(1,4)+res_f(1,4)+res_cd(1,4)+res_ce(1,4)+&
         res_cf(1,4)+res_de(1,4)+res_df(1,4)+res_ef(1,4))*half !-- symm factor

    !--u c > u c em ep
    res_jj(7) = res_c(2,3)+res_d(2,3)+res_cd(2,3)

    !--u d > u d em ep, u s > u s em ep, c d > c d em ep, c s > c s em ep
    res_jj(8) = res_c(2,3)+res_d(1,3)+res_cd(4,3)

    !--u ux > c cx em ep, c cx > u ux em ep
    res_jj(9) = res_a(2,1)+res_b(2,1)+res_ab(2,1)

    !--u ux > d dx em ep, u ux > s sx em ep, c cx > d dx em ep, c cx > s sx em ep
    res_jj(10) = res_a(2,1)+res_b(1,1)+res_ab(4,1)

    !--u cx > u cx em ep, c ux > c ux em ep
    res_jj(11) = res_c(2,1)+res_d(2,1)+res_cd(2,1)

    !--u dx > u dx em ep, u sx > u sx em ep, c dx > c dx em ep, c sx > c sx em ep
    res_jj(12) = res_c(2,1)+res_d(1,1)+res_cd(4,1)

    !--d s > d s em ep
    res_jj(13) = res_c(1,3)+res_d(1,3)+res_cd(1,3)

    !--d ux > d ux em ep, d cx > d cx em ep, s ux > s ux em ep, s cx > s cx em ep
    res_jj(14) = res_c(1,1)+res_d(2,1)+res_cd(3,1)

    !--d dx > u ux em ep, d dx > c cx em ep, s sx > u ux em ep, s sx > c cx em ep
    res_jj(15) = res_a(1,1)+res_b(2,1)+res_ab(3,1)

    !--d dx > s sx em ep, s sx > d dx em ep
    res_jj(16) = res_a(1,1)+res_b(1,1)+res_ab(1,1)

    !--d sx > d sx em ep, s dx > s dx em ep
    res_jj(17) = res_c(1,1)+res_d(1,1)+res_cd(1,1)

    !--ux cx > ux cx em ep
    res_jj(18) = res_c(2,4)+res_d(2,4)+res_cd(2,4)

    !--ux dx > ux dx em ep, ux sx > ux sx em ep, cx dx > cx dx em ep, cx sx > cx sx em ep
    res_jj(19) = res_c(2,4)+res_d(1,4)+res_cd(4,4)

    !--dx sx > dx sx em ep
    res_jj(20) = res_c(1,4)+res_d(1,4)+res_cd(1,4)

    !--ux u > u ux em ep, cx c > c cx em ep MIRROR
    res_jj(21) = (res_a(2,2)+res_b(2,2)+res_c(2,2)+res_d(2,2)+&
         res_ab(2,2)+res_ac(2,2)+res_ad(2,2)+res_bc(2,2)+res_bd(2,2)+res_cd(2,2))

    !--dx d > d dx em ep, sx s > s sx em ep MIRROR
    res_jj(22) = (res_a(1,2)+res_b(1,2)+res_c(1,2)+res_d(1,2)+&
         res_ab(1,2)+res_ac(1,2)+res_ad(1,2)+res_bc(1,2)+res_bd(1,2)+res_cd(1,2))

    !--c u > u c em ep MIRROR
    res_jj(23) = res_e(2,3)+res_f(2,3)+res_ef(2,3)

    !--d u > u d em ep, s u > u s em ep, d c > c d em ep, s c > c s em ep MIRROR
    res_jj(24) = res_e(2,3)+res_f(1,3)+res_ef(4,3)

    !--ux u > c cx em ep, cx c > u ux em ep MIRROR
    res_jj(25) = res_a(2,2)+res_b(2,2)+res_ab(2,2)

    !--ux u > d dx em ep, ux u > s sx em ep, cx c > d dx em ep, cx c > s sx em ep MIRROR
    res_jj(26) = res_a(2,2)+res_b(1,2)+res_ab(4,2)

    !--cx u > u cx em ep, ux c > c ux em ep MIRROR
    res_jj(27) = res_c(2,2)+res_d(2,2)+res_cd(2,2)

    !--dx u > u dx em ep, sx u > u sx em ep, dx c > c dx em ep, sx c > c sx em ep MIRROR
    res_jj(28) = res_c(2,2)+res_d(1,2)+res_cd(4,2)

    !--s d > d s em ep MIRROR
    res_jj(29) = res_e(1,3)+res_f(1,3)+res_ef(1,3)

    !--ux d > d ux em ep, cx d > d cx em ep, ux s > s ux em ep, cx s > s cx em ep MIRROR
    res_jj(30) = res_c(1,2)+res_d(2,2)+res_cd(3,2)

    !--dx d > u ux em ep, dx d > c cx em ep, sx s > u ux em ep, sx s > c cx em ep MIRROR
    res_jj(31) = res_a(1,2)+res_b(2,2)+res_ab(3,2)

    !--dx d > s sx em ep, sx s > d dx em ep MIRROR
    res_jj(32) = res_a(1,2)+res_b(1,2)+res_ab(1,2)

    !--sx d > d sx em ep, dx s > s dx em ep MIRROR
    res_jj(33) = res_c(1,2)+res_d(1,2)+res_cd(1,2)

    !--cx ux > ux cx em ep MIRROR
    res_jj(34) = res_e(2,4)+res_f(2,4)+res_ef(2,4)

    !--dx ux > ux dx em ep, sx ux > ux sx em ep, dx cx > cx dx em ep, sx cx > cx sx em ep MIRROR
    res_jj(35) = res_e(2,4)+res_f(1,4)+res_ef(4,4)

    !--sx dx > dx sx em ep MIRROR
    res_jj(36) = res_e(1,4)+res_f(1,4)+res_ef(1,4)

#if (_CheckAmps == 1)
    include './CheckAmps/checkmad_wp_2q.f90'
#endif

    return

  end subroutine restree_2q_jj


  !-- g(1) + g(2) -> q(3) + qb(4) + V(5+6)
  subroutine restree_gg_jj(p,res_jj)
    real(dp), intent(in) :: p(4,6)
    real(dp), intent(out) :: res_jj(2)
    real(dp) :: sprod(6,6),res_a(1:2),res_b(1:2)
    complex(dp) :: za(6,6),zb(6,6)
    complex(dp) :: coupl(-1:1,-1:1,1:2)
    complex(dp) :: amp1(-1:1,-1:1,-1:1,-1:1,1:2),amp2(-1:1,-1:1,-1:1,-1:1,1:2)
    integer :: h1,h2,h3,h4

    call spinoru(6,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5),p(:,6)/),za,zb,sprod)

    call get_coupl_wp(sprod(5,6),coupl)

    res_a = zero
    res_b = zero

    call master_amp_qggqbV_wp(coupl,3,1,2,4,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV_wp(coupl,3,2,1,4,5,6,za,zb,sprod,amp2)

    do h4 = -1,1,2
    do h3 = -1,1,2
    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a(1:2) = res_a(1:2) + abs(amp1(h1,h2,h3,h4,1:2))**2 + abs(amp2(h1,h2,h3,h4,1:2))**2
       res_b(1:2) = res_b(1:2) + two * real(amp1(h1,h2,h3,h4,1:2)*conjg(amp2(h1,h3,h2,h4,1:2)),kind=dp)
       
    enddo
    enddo
    enddo
    enddo

    res_a = res_a * 16.0_dp * Cf**2 * xn * tagab!-- colf -> 4*Cf, 4 -> overall factor 2 of amplitudes
    res_b = res_b * 16.0_dp * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--g g > u ux em ep, g g > c cx em ep
    res_jj(1) = (res_a(2)+res_b(2))*avegg

    !--g g > d dx em ep, g g > s sx em ep
    res_jj(2) = (res_a(1)+res_b(1))*avegg

#if (_CheckAmps == 1)
    include './CheckAmps/checkmad_wp_gg.f90'
#endif

    return

  end subroutine restree_gg_jj

  !-- q/qb(1) + g(2) -> q/qb(3) + qb(4) + V(5+6)
  subroutine restree_qg_jj(p,res_jj)
    real(dp), intent(in) :: p(4,6)
    real(dp), intent(out) :: res_jj(4)
    real(dp) :: sprod(6,6),res_a(1:2),res_b(1:2)
    complex(dp) :: za(6,6),zb(6,6)
    complex(dp) :: coupl(-1:1,-1:1,1:2)
    complex(dp) :: amp1(-1:1,-1:1,-1:1,-1:1,1:2),amp2(-1:1,-1:1,-1:1,-1:1,1:2)
    integer :: h1,h2,h3,h4

    call spinoru(6,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5),p(:,6)/),za,zb,sprod)

    call get_coupl_wp(sprod(5,6),coupl)

    !-- q g -> q g

    res_a = zero
    res_b = zero

    call master_amp_qggqbV_wp(coupl,3,4,2,1,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV_wp(coupl,3,2,4,1,5,6,za,zb,sprod,amp2)

    do h4 = -1,1,2
    do h3 = -1,1,2
    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a(1:2) = res_a(1:2) + abs(amp1(h1,h2,h3,h4,1:2))**2 + abs(amp2(h1,h2,h3,h4,1:2))**2
       res_b(1:2) = res_b(1:2) + two * real(amp1(h1,h2,h3,h4,1:2)*conjg(amp2(h1,h3,h2,h4,1:2)),kind=dp)
       
    enddo
    enddo
    enddo
    enddo

    res_a = res_a * 16.0_dp * Cf**2 * xn * tagab!-- colf -> 4*Cf, 4 -> overall factor 2 of amplitudes
    res_b = res_b * 16.0_dp * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--u g > u g em ep, c g > c g em ep 
    res_jj(1) = (res_a(2)+res_b(2))*aveqg

    !--d g > d g em ep, s g > s g em ep 
    res_jj(2) = (res_a(1)+res_b(1))*aveqg

    ! --- 

    !-- qx g -> qx g

    res_a = zero
    res_b = zero

    call master_amp_qggqbV_wp(coupl,1,4,2,3,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV_wp(coupl,1,2,4,3,5,6,za,zb,sprod,amp2)

    do h4 = -1,1,2
    do h3 = -1,1,2
    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a(1:2) = res_a(1:2) + abs(amp1(h1,h2,h3,h4,1:2))**2 + abs(amp2(h1,h2,h3,h4,1:2))**2
       res_b(1:2) = res_b(1:2) + two * real(amp1(h1,h2,h3,h4,1:2)*conjg(amp2(h1,h3,h2,h4,1:2)),kind=dp)
       
    enddo
    enddo
    enddo
    enddo

    res_a = res_a * 16.0_dp * Cf**2 * xn * tagab!-- colf -> 4*Cf, 4 -> overall factor 2 of amplitudes
    res_b = res_b * 16.0_dp * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--ux g > ux g em ep, cx g > g cx em ep 
    res_jj(3) = (res_a(2)+res_b(2))*aveqg

    !--dx g > dx g em ep, sx g > g sx em ep 
    res_jj(4) = (res_a(1)+res_b(1))*aveqg

#if (_CheckAmps == 1)
    include './CheckAmps/checkmad_wp_qg.f90'
#endif

    return

  end subroutine restree_qg_jj

  !-- g(1) + q/qb(2) -> q/qb(3) + qb(4) + V(5+6)
  subroutine restree_gq_jj(p,res_jj)
    real(dp), intent(in) :: p(4,6)
    real(dp), intent(out) :: res_jj(4)
    real(dp) :: sprod(6,6),res_a(1:2),res_b(1:2)
    complex(dp) :: za(6,6),zb(6,6)
    complex(dp) :: coupl(-1:1,-1:1,1:2)
    complex(dp) :: amp1(-1:1,-1:1,-1:1,-1:1,1:2),amp2(-1:1,-1:1,-1:1,-1:1,1:2)
    integer :: h1,h2,h3,h4

    call spinoru(6,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5),p(:,6)/),za,zb,sprod)

    call get_coupl_wp(sprod(5,6),coupl)

    !-- g q -> q g

    res_a = zero
    res_b = zero

    call master_amp_qggqbV_wp(coupl,3,4,1,2,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV_wp(coupl,3,1,4,2,5,6,za,zb,sprod,amp2)

    do h4 = -1,1,2
    do h3 = -1,1,2
    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a(1:2) = res_a(1:2) + abs(amp1(h1,h2,h3,h4,1:2))**2 + abs(amp2(h1,h2,h3,h4,1:2))**2
       res_b(1:2) = res_b(1:2) + two * real(amp1(h1,h2,h3,h4,1:2)*conjg(amp2(h1,h3,h2,h4,1:2)),kind=dp)
       
    enddo
    enddo
    enddo
    enddo

    res_a = res_a * 16.0_dp * Cf**2 * xn * tagab!-- colf -> 4*Cf, 4 -> overall factor 2 of amplitudes
    res_b = res_b * 16.0_dp * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--g u > u g em ep, g c > c g em ep
    res_jj(1) = (res_a(2)+res_b(2))*aveqg

    !--g d > d g em ep, g s > s g em ep
    res_jj(2) = (res_a(1)+res_b(1))*aveqg

    ! --- 

    !-- g qx -> qx g

    res_a = zero
    res_b = zero

    call master_amp_qggqbV_wp(coupl,2,4,1,3,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV_wp(coupl,2,1,4,3,5,6,za,zb,sprod,amp2)

    do h4 = -1,1,2
    do h3 = -1,1,2
    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a(1:2) = res_a(1:2) + abs(amp1(h1,h2,h3,h4,1:2))**2 + abs(amp2(h1,h2,h3,h4,1:2))**2
       res_b(1:2) = res_b(1:2) + two * real(amp1(h1,h2,h3,h4,1:2)*conjg(amp2(h1,h3,h2,h4,1:2)),kind=dp)
       
    enddo
    enddo
    enddo
    enddo

    res_a = res_a * 16.0_dp * Cf**2 * xn * tagab!-- colf -> 4*Cf, 4 -> overall factor 2 of amplitudes
    res_b = res_b * 16.0_dp * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--g ux > g ux em ep, g cx > g cx em ep
    res_jj(3) = (res_a(2)+res_b(2))*aveqg

    !--g dx > g dx em ep, g sx > g sx em ep
    res_jj(4) = (res_a(1)+res_b(1))*aveqg

#if (_CheckAmps == 1)
    include './CheckAmps/checkmad_wp_gq.f90'
#endif

    return

  end subroutine restree_gq_jj


  !-- q/qb(1) + qb/q(2) -> g(3) + g(4) + V(5+6)
  subroutine restree_qqb_jj(p,res_jj)
    real(dp), intent(in) :: p(4,6)
    real(dp), intent(out) :: res_jj(4)
    real(dp) :: sprod(6,6),res_a(1:2),res_b(1:2)
    complex(dp) :: za(6,6),zb(6,6)
    complex(dp) :: coupl(-1:1,-1:1,1:2)
    complex(dp) :: amp1(-1:1,-1:1,-1:1,-1:1,1:2),amp2(-1:1,-1:1,-1:1,-1:1,1:2)
    integer :: h1,h2,h3,h4

    call spinoru(6,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5),p(:,6)/),za,zb,sprod)

    call get_coupl_wp(sprod(5,6),coupl)

    !-- q qx > g g

    res_a = zero
    res_b = zero

    call master_amp_qggqbV_wp(coupl,2,3,4,1,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV_wp(coupl,2,4,3,1,5,6,za,zb,sprod,amp2)

    do h4 = -1,1,2
    do h3 = -1,1,2
    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a(1:2) = res_a(1:2) + abs(amp1(h1,h2,h3,h4,1:2))**2 + abs(amp2(h1,h2,h3,h4,1:2))**2
       res_b(1:2) = res_b(1:2) + two * real(amp1(h1,h2,h3,h4,1:2)*conjg(amp2(h1,h3,h2,h4,1:2)),kind=dp)
       
    enddo
    enddo
    enddo
    enddo

    res_a = res_a * 16.0_dp * Cf**2 * xn * tagab!-- colf -> 4*Cf, 4 -> overall factor 2 of amplitudes
    res_b = res_b * 16.0_dp * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--u ux > g g em ep, c cx > g g em ep
    res_jj(1) = (res_a(2)+res_b(2))*aveqq

    !--d dx > g g em ep, s sx > g g em ep
    res_jj(2) = (res_a(1)+res_b(1))*aveqq

    ! --- 

    !-- qx q > g g 

    res_a = zero
    res_b = zero

    call master_amp_qggqbV_wp(coupl,1,3,4,2,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV_wp(coupl,1,4,3,2,5,6,za,zb,sprod,amp2)

    do h4 = -1,1,2
    do h3 = -1,1,2
    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a(1:2) = res_a(1:2) + abs(amp1(h1,h2,h3,h4,1:2))**2 + abs(amp2(h1,h2,h3,h4,1:2))**2
       res_b(1:2) = res_b(1:2) + two * real(amp1(h1,h2,h3,h4,1:2)*conjg(amp2(h1,h3,h2,h4,1:2)),kind=dp)
       
    enddo
    enddo
    enddo
    enddo

    res_a = res_a * 16.0_dp * Cf**2 * xn * tagab!-- colf -> 4*Cf, 4 -> overall factor 2 of amplitudes
    res_b = res_b * 16.0_dp * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--ux u > g g em ep, cx c > g g em ep MIRROR
    res_jj(3) = (res_a(2)+res_b(2))*aveqq

    !--dx d > g g em ep, sx s > g g em ep MIRROR
    res_jj(4) = (res_a(1)+res_b(1))*aveqq

#if (_CheckAmps == 1)
    include './CheckAmps/checkmad_wp_qqb.f90'
#endif

    return

  end subroutine restree_qqb_jj

  !-- 

  !----------------------------------------------------------------------------!
  !-- MASTER AMPLITUDES                                                      --!
  !----------------------------------------------------------------------------!


  !-- amplitude for 0 -> q(1) qb(2) [g* -> Q(3) QB(4)] qb(4) [V --> l(5) lb(6)]
  !-- V attached to the 12 line
  !-- index: h1,h3,h5,down/up
  subroutine master_amp_qqbQQBV_wp(coupl,i1,i2,i3,i4,i5,i6,za,zb,sprod,amp)
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    complex(dp), intent(in) :: za(6,6),zb(6,6),coupl(-1:1,-1:1,1:2)
    real(dp), intent(in) :: sprod(6,6)
    complex(dp), intent(out) :: amp(-1:1,-1:1,-1:1,1:2)
    complex(dp) :: zab2
    real(dp) :: dena,denb

    zab2(i1,i2,i3,i4) = za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

    dena = (sprod(i1,i3)+sprod(i1,i4)+sprod(i3,i4))*sprod(i3,i4)*sprod(i5,i6) !-- s134*s34*s56
    denb = (sprod(i2,i3)+sprod(i2,i4)+sprod(i3,i4))*sprod(i3,i4)*sprod(i5,i6) !-- s234*s34*s56

    amp(-1,-1,-1,:) =  (za(i1,i3)*zb(i6,i2)*zab2(i5,i1,i3,i4)/dena-za(i1,i5)*zb(i4,i2)*zab2(i3,i2,i4,i6)/denb)*&
         coupl(-1,-1,:)
    amp(-1,-1,+1,:) =  (za(i1,i3)*zb(i5,i2)*zab2(i6,i1,i3,i4)/dena-za(i1,i6)*zb(i4,i2)*zab2(i3,i2,i4,i5)/denb)*&
         coupl(-1,+1,:)
    amp(-1,+1,-1,:) =  (za(i1,i4)*zb(i6,i2)*zab2(i5,i1,i4,i3)/dena-za(i1,i5)*zb(i3,i2)*zab2(i4,i2,i3,i6)/denb)*&
         coupl(-1,-1,:)
    amp(-1,+1,+1,:) =  (za(i1,i4)*zb(i5,i2)*zab2(i6,i1,i4,i3)/dena-za(i1,i6)*zb(i3,i2)*zab2(i4,i2,i3,i5)/denb)*&
         coupl(-1,+1,:)
    !
    amp(+1,-1,-1,:) = (-za(i2,i3)*zb(i6,i1)*zab2(i5,i2,i3,i4)/denb+za(i2,i5)*zb(i4,i1)*zab2(i3,i1,i4,i6)/dena)*&
         coupl(+1,-1,:)
    amp(+1,-1,+1,:) = (-za(i2,i3)*zb(i5,i1)*zab2(i6,i2,i3,i4)/denb+za(i2,i6)*zb(i4,i1)*zab2(i3,i1,i4,i5)/dena)*&
         coupl(+1,+1,:)
    amp(+1,+1,-1,:) = (-za(i2,i4)*zb(i6,i1)*zab2(i5,i2,i4,i3)/denb+za(i2,i5)*zb(i3,i1)*zab2(i4,i1,i3,i6)/dena)*&
         coupl(+1,-1,:)
    amp(+1,+1,+1,:) = (-za(i2,i4)*zb(i5,i1)*zab2(i6,i2,i4,i3)/denb+za(i2,i6)*zb(i3,i1)*zab2(i4,i1,i3,i5)/dena)*&
         coupl(+1,+1,:)

    amp = amp * two

    return

  end subroutine master_amp_qqbQQBV_wp

  !-- 0 -> q(1) g(2) g(3) qb(4) [V --> l(5) lb(6)]
  !-- index: h1,h2,h3,h5
  subroutine master_amp_qggqbV_wp(coupl,i1,i2,i3,i4,i5,i6,za,zb,sprod,ampl)
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    complex(dp), intent(in) :: za(6,6),zb(6,6),coupl(-1:1,-1:1,1:2)
    real(dp), intent(in) :: sprod(6,6)
    complex(dp), intent(out) :: ampl(-1:1,-1:1,-1:1,-1:1,1:2)
    complex(dp) :: zab2,app,amp,apm,amm
    real(dp) :: s3
    
    !-- statement functions

    zab2(i1,i2,i3,i4) = za(i1,i2)*zb(i2,i4) + za(i1,i3)*zb(i3,i4)
    s3(i1,i2,i3) = sprod(i1,i2) + sprod(i1,i3) + sprod(i2,i3)

    !-- 5- 6+
    app(i1,i2,i3,i4,i5,i6) = -za(i1,i5)**2*zb(i6,i5)/(za(i1,i2)*za(i2,i3)*za(i3,i4))
    amp(i1,i2,i3,i4,i5,i6) = (zb(i1,i3)*zb(i5,i6)*zab2(i5,i1,i2,i3)**2/(s3(i1,i2,i3)*zb(i1,i2)*zb(i2,i3)) + &
         za(i2,i4)*za(i5,i6)*zab2(i2,i3,i4,i6)**2/(s3(i2,i3,i4)*za(i2,i3)*za(i3,i4)))/zab2(i4,i2,i3,i1)
    apm(i1,i2,i3,i4,i5,i6) = (zb(i2,i4)**3*za(i1,i5)**2*zb(i5,i6)/(s3(i2,i3,i4)*zb(i2,i3)*zb(i3,i4)) + &
         za(i1,i3)**3*zb(i6,i4)**2*za(i5,i6)/(s3(i1,i2,i3)*za(i1,i2)*za(i2,i3)))/zab2(i1,i2,i3,i4)
    amm(i1,i2,i3,i4,i5,i6) = za(i5,i6)*zb(i6,i4)**2/(zb(i1,i2)*zb(i2,i3)*zb(i3,i4)) !-- up to a sign

    !--

    ampl(-1,+1,+1,-1,:) = app(i1,i2,i3,i4,i5,i6) * coupl(-1,-1,:) / sprod(i5,i6)
    ampl(-1,+1,+1,+1,:) = app(i1,i2,i3,i4,i6,i5) * coupl(-1,+1,:) / sprod(i5,i6)

    ampl(-1,+1,-1,-1,:) = apm(i1,i2,i3,i4,i5,i6) * coupl(-1,-1,:) / sprod(i5,i6)
    ampl(-1,+1,-1,+1,:) = apm(i1,i2,i3,i4,i6,i5) * coupl(-1,+1,:) / sprod(i5,i6)

    ampl(-1,-1,+1,-1,:) = amp(i1,i2,i3,i4,i5,i6) * coupl(-1,-1,:) / sprod(i5,i6)
    ampl(-1,-1,+1,+1,:) = amp(i1,i2,i3,i4,i6,i5) * coupl(-1,+1,:) / sprod(i5,i6)

    ampl(-1,-1,-1,-1,:) = amm(i1,i2,i3,i4,i5,i6) * coupl(-1,-1,:) / sprod(i5,i6)
    ampl(-1,-1,-1,+1,:) = amm(i1,i2,i3,i4,i6,i5) * coupl(-1,+1,:) / sprod(i5,i6)

    !-- 

    ampl(+1,+1,+1,-1,:) = app(i4,i3,i2,i1,i5,i6) * coupl(+1,-1,:) / sprod(i5,i6)
    ampl(+1,+1,+1,+1,:) = app(i4,i3,i2,i1,i6,i5) * coupl(+1,+1,:) / sprod(i5,i6)

    ampl(+1,+1,-1,-1,:) = amp(i4,i3,i2,i1,i5,i6) * coupl(+1,-1,:) / sprod(i5,i6)
    ampl(+1,+1,-1,+1,:) = amp(i4,i3,i2,i1,i6,i5) * coupl(+1,+1,:) / sprod(i5,i6)

    ampl(+1,-1,+1,-1,:) = apm(i4,i3,i2,i1,i5,i6) * coupl(+1,-1,:) / sprod(i5,i6)
    ampl(+1,-1,+1,+1,:) = apm(i4,i3,i2,i1,i6,i5) * coupl(+1,+1,:) / sprod(i5,i6)

    ampl(+1,-1,-1,-1,:) = amm(i4,i3,i2,i1,i5,i6) * coupl(+1,-1,:) / sprod(i5,i6)
    ampl(+1,-1,-1,+1,:) = amm(i4,i3,i2,i1,i6,i5) * coupl(+1,+1,:) / sprod(i5,i6)

    return

  end subroutine master_amp_qggqbV_wp

end module mod_amplitudes_tree_wp

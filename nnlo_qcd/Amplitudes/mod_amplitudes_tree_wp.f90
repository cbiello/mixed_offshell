!-- amplitudes with extra jet: factor out gs**2 per extra emission
module mod_amplitudes_tree_wp
  use mod_types; use mod_consts_dp
  use mod_auxfunctions
  implicit none
  private

  public :: restree, restree_j
  !  public :: restree_gg_qqb, restree_qg_qg, restree_gq_qg
  public :: restree_ns_jj, restree_gg_qqb, restree_qg_qg, restree_gq_qg, restree_s_jj
  public :: master_amp_qgqbV

contains

#include "./inc_coupl_az.f90"

  !----------------------------------------------------------------------------!
  !-- LO                                                                     --!
  !----------------------------------------------------------------------------!

  !-- f(1) + f(2) --> [W --> l(3) + lb(4)] [H --> b(5) bb(6)]
  !-- summed / averaged
  subroutine restree(p,res)
    real(dp), intent(in) :: p(4,4)
    real(dp), intent(out) :: res(2)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4),zb(4,4),coupl(-1:1,-1:1,1:2)
    complex(dp) :: amp
    integer :: i1,i2,i3,i4

!    real(dp) :: res_dec,propsq_w,propsq_vh,propsq_h,prefac

!!    propsq_w = one/((two*scr(p(:,3),p(:,4))-mwsq)**2+mwsq*gaw**2)
!!    prefac = gwsq**2*propsq_w*xn
!!
!!    !--u dx >  ve em b bx, c sx > ve ep b bx
!!    res(1) = scr(p(:,2),p(:,3))*scr(p(:,1),p(:,4)) * prefac
!!
!!    !--dx u > ve ep b bx, sx c > ve ep b bx MIRROR
!!    res(2) = scr(p(:,2),p(:,4))*scr(p(:,1),p(:,3)) * prefac
!!
!!    res = res * aveqq
    
    amp(i1,i2,i3,i4) = two * zb(i2,i4)*za(i3,i1)/sprod(i1,i2)


    call spinoru(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)
    call get_coupl_w(sprod(1,2),coupl)

    ! u(1) db(2) -> W+ nue(3) ep(4)                 
    res(1) = abs(amp(2,1,3,4)*coupl(-1,-1,1))**2

    ! db(1) u(2) -> W+ nue(3) ep(4)                 
    res(2) = abs(amp(1,2,3,4)*coupl(-1,-1,1))**2

    res = res * xn * aveqq

#if (_CheckAmps == 1)
    include './CheckAmps/checkmad_wp_LO.f90'
#endif


    return

  end subroutine restree

  !----------------------------------------------------------------------------!
  !-- NLO                                                                    --!
  !----------------------------------------------------------------------------!

  !-- f(1) + f(2) --> f(3) [W --> l(4) + lb(5)] [H --> b(6) bb(7)]
  !-- summed / averaged
  subroutine restree_j(p,res_sping_j,res_j)
    real(dp), intent(in) :: p(4,5)
    real(dp), intent(out) :: res_j(6)
    complex(dp), intent(out) :: res_sping_j(6,-1:1,-1:1)
    real(dp) :: sprod(5,5),s123
    real(dp) :: propsq_w,propsq_h,propsq_vh,res_dec,prefac
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

       call master_amp_qgqbV(part_conf(1,i),part_conf(2,i),part_conf(3,i),4,5,za,zb,sprod,amp)

       do ii = -1,1,2
       do jj = -1,1,2

          res_sping_j(i,ii,jj) = res_sping_j(i,ii,jj) + amp(ii)*conjg(amp(jj))

       enddo
       enddo

    enddo

    res_sping_j = res_sping_j * prefac

    !--u dx > g ve ep b bx, c sx > g ve ep b bx
    res_sping_j(1,:,:) = res_sping_j(1,:,:) * avecol_qq

    !--g u > d ve ep b bx, g c > s ve ep b bx
    res_sping_j(2,:,:) = res_sping_j(2,:,:) * avecol_qg

    !--g dx > ux ve ep b bx, g sx > cx ve ep b bx
    res_sping_j(3,:,:) = res_sping_j(3,:,:) * avecol_qg

    !--dx u > g ve ep b bx, sx c > g ve ep b bx MIRROR
    res_sping_j(4,:,:) = res_sping_j(4,:,:) * avecol_qq

    !--u g > d ve ep b bx, c g > s ve ep b bx MIRROR
    res_sping_j(5,:,:) = res_sping_j(5,:,:) * avecol_qg

    !--dx g > ux ve ep b bx, sx g > cx ve ep b bx MIRROR
    res_sping_j(6,:,:) = res_sping_j(6,:,:) * avecol_qg

    res_j(:) = res_sping_j(:,+1,+1)+res_sping_j(:,-1,-1)

#if (_CheckAmps == 1)
    include './CheckAmps/checkmad_wp_NLO.f90'
#endif

    return

  end subroutine restree_j

  !----------------------------------------------------------------------------!
  !-- NNLO                                                                   --!
  !----------------------------------------------------------------------------!


  !-- f(1) + f(2) -> f(3) + f(4) + [V-->l(5)+lb(6)]
  !-- summed/averaged
  !-- a...d notation: see vN
  !-- all individual partonic channels checked against madgraph
  subroutine restree_s_jj(p,res_jj)
    real(dp), intent(in) :: p(4,6)
    real(dp), intent(out) :: res_jj(26)
    complex(dp) :: za(6,6),zb(6,6)
    real(dp) :: sprod(6,6),res_dec,s1234
    real(dp) :: propsq_w,propsq_h,propsq_vh,prefac
    integer :: part_conf(4,6,4,2) !-- particle ordering, a-f amplitude, (qqb,qbq,qq,qbqb)
    complex(dp) :: amp_c(-1:1)
    complex(dp) :: amp_d(-1:1)
    complex(dp) :: amp_e(-1:1)
    complex(dp) :: amp_f(-1:1)
    real(dp) :: res_c(1:4) !-- qqb/qbq/qq/qbqb
    real(dp) :: res_d(1:4),res_e(1:4),res_f(1:4) 
    real(dp) :: res_ce(1:4),res_cf(1:4)
    real(dp) :: res_de(1:4),res_df(1:4)
    real(dp), parameter :: pref1 = xn*Cf*tr*aveqq
    real(dp), parameter :: pref2 = xn*half*Cf*(Ca-two*Cf)*aveqq
    integer :: i,j,jj
    !

    call spinoru(6,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5),p(:,6)/),za,zb,sprod)
    propsq_w = one/((two*scr(p(:,5),p(:,6))-mwsq)**2+mwsq*gaw**2)
    prefac = propsq_w * gwsq**2 * four

    !--

    part_conf(1:4,3,1,1) = (/3,1,2,4/) !-- c,qqb
    part_conf(1:4,4,1,1) = (/2,4,3,1/) !-- d,qqb
                     
    part_conf(1:4,3,2,1) = (/3,2,1,4/) !-- c,qbq
    part_conf(1:4,4,2,1) = (/1,4,3,2/) !-- d,qbq
                     
    part_conf(1:4,3,3,1) = (/3,1,4,2/) !-- c,qq
    part_conf(1:4,4,3,1) = (/4,2,3,1/) !-- d,qq
    part_conf(1:4,5,3,1) = (/3,2,4,1/) !-- e,qq
    part_conf(1:4,6,3,1) = (/4,1,3,2/) !-- f,qq
                     
    part_conf(1:4,3,4,1) = (/1,3,2,4/) !-- c,qbqb
    part_conf(1:4,4,4,1) = (/2,4,1,3/) !-- d,qbqb
    part_conf(1:4,5,4,1) = (/2,3,1,4/) !-- e,qbqb
    part_conf(1:4,6,4,1) = (/1,4,2,3/) !-- f,qbqb

    ! symmetrize
    part_conf(1:4,3,1,2) = (/4,1,2,3/) !-- c,qqb
    part_conf(1:4,4,1,2) = (/2,3,4,1/) !-- d,qqb
                      
    part_conf(1:4,3,2,2) = (/4,2,1,3/) !-- c,qbq
    part_conf(1:4,4,2,2) = (/1,3,4,2/) !-- d,qbq
                      
    part_conf(1:4,3,3,2) = (/4,1,3,2/) !-- c,qq
    part_conf(1:4,4,3,2) = (/3,2,4,1/) !-- d,qq
    part_conf(1:4,5,3,2) = (/4,2,3,1/) !-- e,qq
    part_conf(1:4,6,3,2) = (/3,1,4,2/) !-- f,qq
                      
    part_conf(1:4,3,4,2) = (/1,4,2,3/) !-- c,qbqb
    part_conf(1:4,4,4,2) = (/2,3,1,4/) !-- d,qbqb
    part_conf(1:4,5,4,2) = (/2,4,1,3/) !-- e,qbqb
    part_conf(1:4,6,4,2) = (/1,3,2,4/) !-- f,qbqb

    res_c = zero
    res_d = zero
    res_e = zero
    res_f = zero

    res_ce = zero
    res_cf = zero    
    res_de = zero
    res_df = zero


        !-- qqb + qbq
    do j = 1,2
       do i = 1,2
    
          call master_amp_qqbQQBV(part_conf(1,3,i,j),part_conf(2,3,i,j),part_conf(3,3,i,j),part_conf(4,3,i,j),5,6,za,zb,sprod,amp_c)
          call master_amp_qqbQQBV(part_conf(1,4,i,j),part_conf(2,4,i,j),part_conf(3,4,i,j),part_conf(4,4,i,j),5,6,za,zb,sprod,amp_d)
          
          do jj = -1,1,2
             
             res_c(i) = res_c(i) + abs(amp_c(jj))**2
             res_d(i) = res_d(i) + abs(amp_d(jj))**2
             
          enddo
          
       enddo
    enddo

    !-- qq + qbqb
    do j = 1,2
    do i = 3,4
    
       call master_amp_qqbQQBV(part_conf(1,3,i,j),part_conf(2,3,i,j),part_conf(3,3,i,j),part_conf(4,3,i,j),5,6,za,zb,sprod,amp_c)
       call master_amp_qqbQQBV(part_conf(1,4,i,j),part_conf(2,4,i,j),part_conf(3,4,i,j),part_conf(4,4,i,j),5,6,za,zb,sprod,amp_d)
       call master_amp_qqbQQBV(part_conf(1,5,i,j),part_conf(2,5,i,j),part_conf(3,5,i,j),part_conf(4,5,i,j),5,6,za,zb,sprod,amp_e)
       call master_amp_qqbQQBV(part_conf(1,6,i,j),part_conf(2,6,i,j),part_conf(3,6,i,j),part_conf(4,6,i,j),5,6,za,zb,sprod,amp_f)

       do jj = -1,1,2

          res_c(i) = res_c(i) + abs(amp_c(jj))**2
          res_d(i) = res_d(i) + abs(amp_d(jj))**2
          res_e(i) = res_e(i) + abs(amp_e(jj))**2
          res_f(i) = res_f(i) + abs(amp_f(jj))**2
          
       enddo

       res_ce(i) = res_ce(i) + two*real(amp_c(-1)*conjg(amp_e(-1)),kind=dp)
       res_cf(i) = res_cf(i) + two*real(amp_c(-1)*conjg(amp_f(-1)),kind=dp)

       res_de(i) = res_de(i) + two*real(amp_d(-1)*conjg(amp_e(-1)),kind=dp)
       res_df(i) = res_df(i) + two*real(amp_d(-1)*conjg(amp_f(-1)),kind=dp)
    
    enddo
    enddo

    res_c = res_c * pref1 * vn_tagc
    res_d = res_d * pref1 * vn_tagc
    res_e = res_e * pref1 * vn_tagc
    res_f = res_f * pref1 * vn_tagc
    
    res_ce = res_ce * pref2 * vn_tagcf
    res_cf = res_cf * pref2 * vn_tagce

    res_de = res_de * pref2 * vn_tagce
    res_df = res_df * pref2 * vn_tagcf

    !--d u > d d em vex b bx, s c > s s em vex b bx
    res_jj(1) = ( (res_d(3)+res_e(3))*tag_s_s+res_de(3)*tag_s_v)*half

    !--d c > d s em vex b bx, s u > s d em vex b bx, u c > u s em vex b bx
    res_jj(2) = res_d(3)*tag_s_s

    !--d sx > d cx em vex b bx, s dx > s ux em vex b bx, u sx > u cx em vex b bx, c dx > c ux em vex b bx
    ! ?? d dx > d ux , s sx > s cx
    res_jj(3) = res_d(1)*tag_s_s

    !--u u > d u em vex b bx, c c > s c em vex b bx
    res_jj(4) = (res_c(3)+res_e(3))*tag_s_s+res_ce(3)*tag_s_v

    !--u c > d c em vex b bx
    res_jj(5) = res_c(3)*tag_s_s

    !--u sx > d sx em vex b bx, u cx > d cx em vex b bx, c dx > s dx em vex b bx, c ux > s ux em vex b bx
    ! ?? u ux > d ux, c cx > s cx
    res_jj(6) = res_c(1)*tag_s_s

    !--dx dx > dx ux em vex b bx, sx sx > sx cx em vex b bx
    res_jj(7) = (res_d(4)+res_f(4))*tag_s_s+res_df(4)*tag_s_v

    !--dx sx > dx cx em vex b bx 
    res_jj(8) = res_d(4)*tag_s_s
    
    !--dx sx > sx ux em vex b bx, dx cx > cx ux em vex b bx, sx ux > ux cx em vex b bx
    res_jj(9) = res_f(4)*tag_s_s

    !--dx ux > ux ux em vex b bx, sx cx > cx cx em vex b bx
    res_jj(10) = ((res_c(4)+res_f(4))*tag_s_s+res_cf(4)*tag_s_v)*half

    !--u d > d d em vex b bx, c s > s s em vex b bx MIRROR
    res_jj(11) = ((res_c(3)+res_f(3))*tag_s_s +res_cf(3)*tag_s_v)*half

    !--c d > d s em vex b bx, u s > s d em vex b bx, c u > u s em vex b bx MIRROR
    res_jj(12) = res_f(3)*tag_s_s

    !--sx d > d cx em vex b bx, dx s > s ux em vex b bx, sx u > u cx em vex b bx, dx c > c ux em vex b bx MIRROR
    !?? dx d > d ux, sx x s cx
    res_jj(13) = res_d(2)*tag_s_s

    !--c u > d c em vex b bx MIRROR
    res_jj(14) = res_e(3)*tag_s_s

    !--sx u > d sx em vex b bx, cx u > d cx em vex b bx, dx c > s dx em vex b bx, ux c > s ux em vex b bx MIRROR
    !??ux u > d ux, cx c > s cx
    res_jj(15) = res_c(2)*tag_s_s

    !--sx dx > dx cx em vex b bx MIRROR 
    res_jj(16) = res_f(4) *tag_s_s

     !--sx dx > sx ux em vex b bx, cx dx > cx ux em vex b bx, ux sx > ux cx em vex b bx MIRROR
    res_jj(17) = res_d(4)*tag_s_s

    !--ux dx > ux ux em vex b bx, cx sx > cx cx em vex b bx MIRROR
    res_jj(18) = ((res_d(4)+res_e(4))*tag_s_s+res_de(4)*tag_s_v)*half

    if (ndn .eq. 3) then

       !--b u > b d em vex b bx, b c > b s em vex b bx
       res_jj(19) = res_d(3)*tag_s_s
       
       !--u bx > d bx em vex b bx, c bx > s bx em vex b bx
       res_jj(20) = res_c(1)*tag_s_s
       
       !--b dx > b ux em vex b bx, b sx > b cx em vex b bx
       res_jj(21) = res_d(1)*tag_s_s
       
       !--bx dx > bx ux em vex b bx, bx sx > bx x em vex b bx
       res_jj(22) = res_d(4)*tag_s_s


       !- u b > b d em vex b bx, c b > b s em vex b bx
       res_jj(23) = res_f(3)*tag_s_s

       !--bx u > d bx em vex b bx, bx c > s bx em vex b bx
       res_jj(24) = res_c(2)*tag_s_s

       !--dx b > b ux em vex b bx, b sx > b cx em vex b bx
       res_jj(25) = res_d(2)*tag_s_s

       !--dx bx > bx ux em vex b bx , sx bx > bx cx em vex b bx 
       res_jj(26) = res_f(4) *tag_s_s

    endif





    res_jj = res_jj * prefac

#if (_CheckAmps == 1)
    include './CheckAmps/checkmad_wp_2q.f90'
#endif


    return

  end subroutine restree_s_jj


    



  !-- q/qb(1) + g(2) -> q/qb(3) + qb(4) + V(5+6)
  subroutine restree_qg_qg(p,res_jj)
    real(dp), intent(in) :: p(4,6)
    real(dp), intent(out) :: res_jj(2)
    complex(dp) :: za(6,6),zb(6,6)
    real(dp) :: sprod(6,6),res_dec,s1234
    real(dp) :: propsq_w,propsq_h,propsq_vh,prefac
    real(dp) :: res_a,res_b
    integer :: h1,h2
    complex(dp) :: amp1(-1:1,-1:1),amp2(-1:1,-1:1)

    call spinoru(6,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5),p(:,6)/),za,zb,sprod)

    propsq_w = one/((two*scr(p(:,5),p(:,6))-mwsq)**2+mwsq*gaw**2)
    prefac = propsq_w * gwsq**2 * four

    !-- q g -> q g

    res_a = zero
    res_b = zero

    call master_amp_qggqbV(3,4,2,1,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV(3,2,4,1,5,6,za,zb,sprod,amp2)

    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a = res_a + abs(amp1(h1,h2))**2 + abs(amp2(h1,h2))**2
       res_b = res_b + two * real(amp1(h1,h2)*conjg(amp2(h2,h1)),kind=dp)
       
    enddo
    enddo

    res_a = res_a * Cf**2 * xn * tagab
    res_b = res_b * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--u g > d g ve ep b bx, c g > s g ve ep b bx
    res_jj(1) = (res_a+res_b)*aveqg*prefac

    ! --- 

    !-- qx g -> qx g

    res_a = zero
    res_b = zero

    call master_amp_qggqbV(1,4,2,3,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV(1,2,4,3,5,6,za,zb,sprod,amp2)

    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a = res_a + abs(amp1(h1,h2))**2 + abs(amp2(h1,h2))**2
       res_b = res_b + two * real(amp1(h1,h2)*conjg(amp2(h2,h1)),kind=dp)
       
    enddo
    enddo

    res_a = res_a * Cf**2 * xn * tagab
    res_b = res_b * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--dx g > ux g em ep, sx g > g cx em ep 
    res_jj(2) = (res_a+res_b)*aveqg*prefac

#if (_CheckAmps == 1)
    include './CheckAmps/checkmad_wp_qg.f90'
#endif

    return

  end subroutine restree_qg_qg

  !-- g(1) + q/qb(2) -> q/qb(3) + qb(4) + V(5+6)
  subroutine restree_gq_qg(p,res_jj)
    real(dp), intent(in) :: p(4,8)
    real(dp), intent(out) :: res_jj(2)
    complex(dp) :: za(6,6),zb(6,6)
    real(dp) :: sprod(6,6),res_dec,s1234
    real(dp) :: propsq_w,propsq_h,propsq_vh,prefac
    real(dp) :: res_a,res_b
    integer :: h1,h2
    complex(dp) :: amp1(-1:1,-1:1),amp2(-1:1,-1:1)

    call spinoru(6,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5),p(:,6)/),za,zb,sprod)

    propsq_w = one/((two*scr(p(:,5),p(:,6))-mwsq)**2+mwsq*gaw**2)
    prefac = propsq_w * gwsq**2 * four

    !-- g q -> q g

    res_a = zero
    res_b = zero

    call master_amp_qggqbV(3,4,1,2,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV(3,1,4,2,5,6,za,zb,sprod,amp2)

    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a = res_a + abs(amp1(h1,h2))**2 + abs(amp2(h1,h2))**2
       res_b = res_b + two * real(amp1(h1,h2)*conjg(amp2(h2,h1)),kind=dp)
       
    enddo
    enddo

    res_a = res_a * Cf**2 * xn * tagab
    res_b = res_b * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--g u > d g ve ep b bx, g s > c g ve ep b bx
    res_jj(1) = (res_a+res_b)*aveqg*prefac

    ! --- 

    !-- g qx -> qx g

    res_a = zero
    res_b = zero

    call master_amp_qggqbV(2,4,1,3,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV(2,1,4,3,5,6,za,zb,sprod,amp2)

    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a = res_a + abs(amp1(h1,h2))**2 + abs(amp2(h1,h2))**2
       res_b = res_b + two * real(amp1(h1,h2)*conjg(amp2(h2,h1)),kind=dp)
       
    enddo
    enddo

    res_a = res_a * Cf**2 * xn * tagab
    res_b = res_b * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--g dx > ux g ve ep b bx, g cx > sx g ve ep b bx
    res_jj(2) = (res_a+res_b)*aveqg*prefac

#if (_CheckAmps == 1)
    include './CheckAmps/checkmad_wp_gq.f90'
#endif

    return

  end subroutine restree_gq_qg

  !-- g(1) + g(2) -> q(3) + qb(4) + V(5+6)
  subroutine restree_gg_qqb(p,res_jj)
    real(dp), intent(in) :: p(4,8)
    real(dp), intent(out) :: res_jj(1)
    complex(dp) :: za(6,6),zb(6,6)
    real(dp) :: sprod(6,6),res_dec,s1234
    real(dp) :: propsq_w,propsq_h,propsq_vh,prefac
    real(dp) :: res_a,res_b
    integer :: h1,h2
    complex(dp) :: amp1(-1:1,-1:1),amp2(-1:1,-1:1)

    call spinoru(6,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5),p(:,6)/),za,zb,sprod)

    propsq_w = one/((two*scr(p(:,5),p(:,6))-mwsq)**2+mwsq*gaw**2)
    prefac = propsq_w * gwsq**2 * four

    !--

    res_a = zero
    res_b = zero

    call master_amp_qggqbV(3,1,2,4,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV(3,2,1,4,5,6,za,zb,sprod,amp2)

    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a = res_a + abs(amp1(h1,h2))**2 + abs(amp2(h1,h2))**2
       res_b = res_b + two * real(amp1(h1,h2)*conjg(amp2(h2,h1)),kind=dp)
       
    enddo
    enddo

    res_a = res_a * Cf**2 * xn * tagab
    res_b = res_b * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--g g > u dx em vex b bx, g g > c sx em vex b bx
    res_jj(1) = (res_a+res_b)*avegg*prefac

#if (_CheckAmps == 1)
    include './CheckAmps/checkmad_wp_gg.f90'
#endif

    return

  end subroutine restree_gg_qqb
  
  subroutine restree_ns_jj(p,res_jj)
    real(dp), intent(in) :: p(4,6)
    real(dp), intent(out) :: res_jj(12)
    complex(dp) :: za(6,6),zb(6,6)
    real(dp) :: sprod(6,6),res_dec,s1234
    real(dp) :: propsq_w,propsq_h,propsq_vh,prefac
    real(dp) :: res_a_gg,res_b_gg
    integer :: part_conf(4,4,4) !-- particle ordering, a-f amplitude, (qqb,qbq,qq,qbqb)
    complex(dp) :: amp_a(-1:1) !-- hel index
    complex(dp) :: amp_b(-1:1)
    complex(dp) :: amp_c(-1:1)
    complex(dp) :: amp_d(-1:1)
    real(dp) :: res_a(1:2) !-- qqb/qbq
    real(dp) :: res_b(1:2) !-- qqb/qbq    
    real(dp) :: res_bc(1:2),res_ac(1:2),res_ad(1:2),res_bd(1:2)
    real(dp) :: res_c(1:4) !-- qqb/qbq/qq/qbqb
    real(dp) :: res_d(1:4)
    real(dp), parameter :: pref1 = xn*Cf*tr*aveqq
    real(dp), parameter :: pref2 = xn*half*Cf*(Ca-two*Cf)*aveqq
    integer :: h1,h2,i,j
    complex(dp) :: amp1(-1:1,-1:1),amp2(-1:1,-1:1)

    res_jj = zero

    call spinoru(6,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5),p(:,6)/),za,zb,sprod)


    propsq_w = one/((two*scr(p(:,5),p(:,6))-mwsq)**2+mwsq*gaw**2)
    prefac = propsq_w * gwsq**2 * four


    !-- q qx > g g

    res_a_gg = zero
    res_b_gg = zero

    call master_amp_qggqbV(2,3,4,1,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV(2,4,3,1,5,6,za,zb,sprod,amp2)

    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a_gg = res_a_gg + abs(amp1(h1,h2))**2 + abs(amp2(h1,h2))**2
       res_b_gg = res_b_gg + two * real(amp1(h1,h2)*conjg(amp2(h2,h1)),kind=dp)
       
    enddo
    enddo

    res_a_gg = res_a_gg * Cf**2 * xn * tagab
    res_b_gg = res_b_gg * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--u dx > g g ve ep b bx, c sx > g g ve ep b bx
    res_jj(1) = (res_a_gg+res_b_gg)*aveqq*prefac*tag_ns_gg

    ! --- 

    !-- qx q > g g 

    res_a_gg = zero
    res_b_gg = zero

    call master_amp_qggqbV(1,3,4,2,5,6,za,zb,sprod,amp1)
    call master_amp_qggqbV(1,4,3,2,5,6,za,zb,sprod,amp2)

    do h2 = -1,1,2
    do h1 = -1,1,2

       res_a_gg = res_a_gg + abs(amp1(h1,h2))**2 + abs(amp2(h1,h2))**2
       res_b_gg = res_b_gg + two * real(amp1(h1,h2)*conjg(amp2(h2,h1)),kind=dp)
       
    enddo
    enddo

    res_a_gg = res_a_gg * Cf**2 * xn * tagab
    res_b_gg = res_b_gg * Cf*(Cf * tagab - Ca*half*tagnab) * xn

    !--dx u > g g ve ep b bx, sx c > g g ve ep b bx MIRROR
    res_jj(2) = (res_a_gg+res_b_gg)*aveqq*prefac*tag_ns_gg

    ! quark contributions
    part_conf(1:4,1,1) = (/2,1,3,4/) !-- a,qqb
    part_conf(1:4,2,1) = (/3,4,2,1/) !-- b,qqb
    part_conf(1:4,3,1) = (/3,1,2,4/) !-- c,qqb
    part_conf(1:4,4,1) = (/2,4,3,1/) !-- d,qqb

    part_conf(1:4,1,2) = (/1,2,3,4/) !-- a,qbq
    part_conf(1:4,2,2) = (/3,4,1,2/) !-- b,qbq
    part_conf(1:4,3,2) = (/3,2,1,4/) !-- c,qbq
    part_conf(1:4,4,2) = (/1,4,3,2/) !-- d,qbq
    ! symmetrize
    part_conf(1:4,1,3) = (/2,1,4,3/) !-- a,qqb
    part_conf(1:4,2,3) = (/4,3,2,1/) !-- b,qqb
    part_conf(1:4,3,3) = (/4,1,2,3/) !-- c,qqb
    part_conf(1:4,4,3) = (/2,3,4,1/) !-- d,qqb

    part_conf(1:4,1,4) = (/1,2,4,3/) !-- a,qbq
    part_conf(1:4,2,4) = (/4,3,1,2/) !-- b,qbq
    part_conf(1:4,3,4) = (/4,2,1,3/) !-- c,qbq
    part_conf(1:4,4,4) = (/1,3,4,2/) !-- d,qbq

    res_a = zero
    res_b = zero
    res_c = zero
    res_d = zero

    res_ac = zero
    res_ad = zero
    res_bc = zero
    res_bd = zero
    

    !-- qqb + qbq

    do i = 1,2
    
       call master_amp_qqbQQBV(part_conf(1,1,i),part_conf(2,1,i),part_conf(3,1,i),part_conf(4,1,i),5,6,za,zb,sprod,amp_a)
       call master_amp_qqbQQBV(part_conf(1,2,i),part_conf(2,2,i),part_conf(3,2,i),part_conf(4,2,i),5,6,za,zb,sprod,amp_b)
       call master_amp_qqbQQBV(part_conf(1,3,i),part_conf(2,3,i),part_conf(3,3,i),part_conf(4,3,i),5,6,za,zb,sprod,amp_c)
       call master_amp_qqbQQBV(part_conf(1,4,i),part_conf(2,4,i),part_conf(3,4,i),part_conf(4,4,i),5,6,za,zb,sprod,amp_d)

       do j = -1,1,2

          res_a(i) = res_a(i) + abs(amp_a(j))**2
          res_b(i) = res_b(i) + abs(amp_b(j))**2
          res_c(i) = res_c(i) + abs(amp_c(j))**2
          res_d(i) = res_d(i) + abs(amp_d(j))**2
          
       enddo

       res_ac(i) = res_ac(i) + two*real(amp_a(-1)*conjg(amp_c(-1)),kind=dp)
       res_ad(i) = res_ad(i) + two*real(amp_a(-1)*conjg(amp_d(-1)),kind=dp)

       res_bc(i) = res_bc(i) + two*real(amp_b(-1)*conjg(amp_c(-1)),kind=dp)
       res_bd(i) = res_bd(i) + two*real(amp_b(-1)*conjg(amp_d(-1)),kind=dp)

    enddo

    !-- symmetrize

    do i = 1,2
    
       call master_amp_qqbQQBV(part_conf(1,1,i+2),part_conf(2,1,i+2),part_conf(3,1,i+2),part_conf(4,1,i+2),5,6,za,zb,sprod,amp_a)
       call master_amp_qqbQQBV(part_conf(1,2,i+2),part_conf(2,2,i+2),part_conf(3,2,i+2),part_conf(4,2,i+2),5,6,za,zb,sprod,amp_b)
       call master_amp_qqbQQBV(part_conf(1,3,i+2),part_conf(2,3,i+2),part_conf(3,3,i+2),part_conf(4,3,i+2),5,6,za,zb,sprod,amp_c)
       call master_amp_qqbQQBV(part_conf(1,4,i+2),part_conf(2,4,i+2),part_conf(3,4,i+2),part_conf(4,4,i+2),5,6,za,zb,sprod,amp_d)

       do j = -1,1,2

          res_a(i) = res_a(i) + abs(amp_a(j))**2
          res_b(i) = res_b(i) + abs(amp_b(j))**2
          res_c(i) = res_c(i) + abs(amp_c(j))**2
          res_d(i) = res_d(i) + abs(amp_d(j))**2
          
       enddo

       res_ac(i) = res_ac(i) + two*real(amp_a(-1)*conjg(amp_c(-1)),kind=dp)
       res_ad(i) = res_ad(i) + two*real(amp_a(-1)*conjg(amp_d(-1)),kind=dp)

       res_bc(i) = res_bc(i) + two*real(amp_b(-1)*conjg(amp_c(-1)),kind=dp)
       res_bd(i) = res_bd(i) + two*real(amp_b(-1)*conjg(amp_d(-1)),kind=dp)

    enddo


    !-- color and symmetry factors
    res_a = res_a * pref1 * vn_taga * prefac
    res_b = res_b * pref1 * vn_tagb * prefac
    res_c = res_c * pref1 * vn_tagc * prefac
    res_d = res_d * pref1 * vn_tagc * prefac
     
    res_ac = res_ac * pref2 * vn_tagac * prefac
    res_ad = res_ad * pref2 * vn_tagac * prefac
    res_bc = res_bc * pref2 * vn_tagbc * prefac
    res_bd = res_bd * pref2 * vn_tagbc * prefac


    if (do_nnlo_qqb_ns .eq. 1) then

       ! u ux > d ux, u ux > s cx
       res_jj(3) = two*res_b(1)*tag_ns_nf + res_bc(1)*tag_ns_id

       ! d dx
       res_jj(4) = two*res_b(1)*tag_ns_nf + res_bd(1)*tag_ns_id
       
       !--u dx
       res_jj(5) = res_a(1)*tag_ns_nf*four + (res_c(1)+res_d(1))*tag_ns_s + (res_ac(1) + res_ad(1))*tag_ns_id

       ! ux u
       res_jj(6) = two*res_b(2)*tag_ns_nf + res_bc(2)*tag_ns_id

       ! dx d
       res_jj(7) = two*res_b(2)*tag_ns_nf + res_bd(2)*tag_ns_id

       !--dx u
       res_jj(8) = res_a(2)*tag_ns_nf*four + (res_c(2)+res_d(2))*tag_ns_s + (res_ac(2) + res_ad(2))*tag_ns_id

       !! added b-quarks
       
       if (ndn .eq. 3) then
          ! b bx > d ux, b bx > s cx
          res_jj(9) = two*res_b(1)*tag_ns_nf

          ! bx b > d ux, bx b > s cx
          res_jj(10) = two*res_b(2)*tag_ns_nf

          !--u dx > b bx
          res_jj(11) = res_a(1)*tag_ns_nf

          !--dx u > b bx
          res_jj(12) = res_a(2)*tag_ns_nf


       endif
       
    endif


#if (_CheckAmps == 1)
    include './CheckAmps/checkmad_wp_qqb.f90'
#endif

    return

  end subroutine restree_ns_jj

  !--------------------------------

  !-- amplitude for 0 -> q(1) g(2) qb(3) [V --> l(4) lb(5)]
  !-- index: hg
  !-- does not include coupling / propagator
  subroutine master_amp_qgqbV(i1,i2,i3,i4,i5,za,zb,sprod,amp)
    integer, intent(in) :: i1,i2,i3,i4,i5
    complex(dp), intent(in) :: za(5,5),zb(5,5)
    real(dp), intent(in) :: sprod(5,5)
    complex(dp), intent(out) :: amp(-1:1)
    complex(dp) :: zab2
    integer :: j1,j2,j3,j4
    
    zab2(j1,j2,j3,j4) = za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

    amp(+1) = za(i1,i4) * zab2(i1,i2,i3,i5)/(za(i1,i2)*za(i2,i3))
    amp(-1) = -zab2(i4,i1,i2,i3) * zb(i5,i3)/(zb(i1,i2)*zb(i2,i3))

    return

  end subroutine master_amp_qgqbV

  !-- 0 -> q(1) g(2) g(3) qb(4) [V --> l(5) lb(6)
  !-- index: hg1,hg2
  !-- does not include coupling / propagator
  subroutine master_amp_qggqbV(i1,i2,i3,i4,i5,i6,za,zb,sprod,amp)
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    complex(dp), intent(out) :: amp(-1:1,-1:1)
    complex(dp) :: zab2,zab3
    real(dp) :: s3

    zab2(i1,i2,i3,i4) = za(i1,i2)*zb(i2,i4) + za(i1,i3)*zb(i3,i4)
    zab3(i1,i2,i3,i4,i5) = za(i1,i2)*zb(i2,i5) + za(i1,i3)*zb(i3,i5) + za(i1,i4)*zb(i4,i5)
    s3(i1,i2,i3) = sprod(i1,i2) + sprod(i1,i3) + sprod(i2,i3)

    amp(+1,+1) = -za(i1,i5)*zab3(i1,i2,i3,i4,i6)/(za(i1,i2)*za(i2,i3)*za(i3,i4))

    amp(-1,+1) = (-zb(i1,i3)*(zb(i3,i6)*s3(i1,i2,i3)+zab2(i4,i1,i2,i3)*zb(i4,i6))*zab2(i5,i1,i2,i3)/(s3(i1,i2,i3)*zb(i1,i2)*zb(i2,i3)) &
         -za(i2,i4)*zab2(i2,i3,i4,i6)*(s3(i2,i3,i4)*za(i5,i2)+za(i5,i1)*zab2(i2,i3,i4,i1))/(s3(i2,i3,i4)*za(i2,i3)*za(i3,i4)))/zab2(i4,i2,i3,i1)

    amp(+1,-1) = (-zb(i2,i4)**3*za(i1,i5)*zab3(i1,i2,i3,i4,i6)/(s3(i2,i3,i4)*zb(i2,i3)*zb(i3,i4)) &
         -za(i1,i3)**3*zab3(i5,i1,i2,i3,i4)*zb(i6,i4)/(s3(i1,i2,i3)*za(i1,i2)*za(i2,i3)))/zab2(i1,i2,i3,i4)

    amp(-1,-1) = -zab3(i5,i1,i2,i3,i4)*zb(i6,i4)/(zb(i1,i2)*zb(i2,i3)*zb(i3,i4)) !-- up to a sign

    return

  end subroutine master_amp_qggqbV

  !-- amplitude for 0 -> q(1) qb(2) [g* -> Q(3) QB(4)] qb(4) [V --> l(5) lb(6)]
  !-- V attached to the 12 line
  !-- index: h3
  !-- does not include coupling / propagator
  subroutine master_amp_qqbQQBV(i1,i2,i3,i4,i5,i6,za,zb,sprod,amp)
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    complex(dp), intent(out) :: amp(-1:1)
    complex(dp) :: zab2
    real(dp) :: dena,denb

    zab2(i1,i2,i3,i4) = za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

    dena = (sprod(i1,i3)+sprod(i1,i4)+sprod(i3,i4))*sprod(i3,i4) !-- s134*s34
    denb = (sprod(i2,i3)+sprod(i2,i4)+sprod(i3,i4))*sprod(i3,i4) !-- s234*s34

    amp(-1) =  (za(i1,i3)*zb(i6,i2)*zab2(i5,i1,i3,i4)/dena-za(i1,i5)*zb(i4,i2)*zab2(i3,i2,i4,i6)/denb)
    amp(+1) =  (za(i1,i4)*zb(i6,i2)*zab2(i5,i1,i4,i3)/dena-za(i1,i5)*zb(i3,i2)*zab2(i4,i2,i3,i6)/denb)

    return

  end subroutine master_amp_qqbQQBV

end module mod_amplitudes_tree_wp

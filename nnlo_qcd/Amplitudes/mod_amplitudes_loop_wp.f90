!-- amplitudes with extra jet: factor our gs**2 per extra emission
!-- loop amplitudes: factor out asontwopi per loop
!-- finite part according to Catani

!-- if we include the 3rd generation, gg triangles should be added to 1-loop squared

module mod_amplitudes_loop_wp
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_amplitudes_tree_wp
  implicit none
  private

  public :: res_1loop
  public :: res_1loop_sq, res_2loop
  public :: res_1loop_j_qqb, res_1loop_j_qg, res_1loop_j_gq

contains

#include "./inc_coupl_az.f90"


  !-- q(1) + qb(2) -> X
  !-- summed / averaged
  subroutine res_1loop(p,musq,res0,res1)
    real(dp), intent(in) :: p(4,4),musq
    real(dp), intent(out) :: res0(2),res1(2)

    call restree(p,res0)

    res1 = res0 * (-8.0_dp * Cf)

    return

  end subroutine res_1loop

  !-- q(1) + qb(2) -> X
  !-- square of the finite part, coefficient of asontwopi**2
  subroutine res_1loop_sq(p,musq,res0,res1)
    real(dp), intent(in) :: p(4,4),musq
    real(dp), intent(out) :: res0(2),res1(2)

    call restree(p,res0)

    res1 = res0 * 16.0_dp * Cf**2

    return

    return

  end subroutine res_1loop_sq

  !-- q(1) + qb(2) -> V(3+4)                                                                                                                                                                                 
  !-- finite part, coefficient of asontwopi**2                                                                                                                                                               
  subroutine res_2loop(p,musq,res0,res2)
    real(dp), intent(in) :: p(4,4),musq
    real(dp), intent(out) :: res0(2),res2(2)
    real(dp) :: sprod(4,4),q2,Lmu,ff
    complex(dp) :: za(4,4),zb(4,4)

    call restree(p,res0)
    q2 = two*scr(p(:,1),p(:,2))

    Lmu = log(musq/q2)

    ff = Cf**2*(255/16.0_dp + (29*pisq)/12.0_dp - (11*pisq**2)/90.0_dp +  &
         Lmu*(-3/4.0_dp + pisq - 12*zeta3) - 15*zeta3)*tag_ns_gg +  &
         Cf*(nf*tag_ns_nf*(4085/648.0_dp + (7*pisq)/36.0_dp +  &
         Lmu*(119/54.0_dp + pisq/12.0_dp) - zeta3/18.0_dp) +  &
         Ca*tag_ns_gg*(-51157/1296.0_dp - (107*pisq)/72.0_dp +  &
         (31*pisq**2)/240.0_dp + (659*zeta3)/36.0_dp + &
         Lmu*(-1339/108.0_dp - (23*pisq)/24.0_dp + 13*zeta3)))

     res2 = ff * res0

    return

  end subroutine res_2loop


!! old Catani subtr !!  !-- q(1) + qb(2) -> X
!! old Catani subtr !!  !-- finite part, coefficient of asontwopi**2
!! old Catani subtr !!  subroutine res_2loop_qqb(p,musq,res0,res2)
!! old Catani subtr !!    real(dp), intent(in) :: p(4,4),musq
!! old Catani subtr !!    real(dp), intent(out) :: res2(2),res0(2)
!! old Catani subtr !!    real(dp) :: l12,res0(2)
!! old Catani subtr !!
!! old Catani subtr !!
!! old Catani subtr !!    call restree(p,res0)
!! old Catani subtr !!
!! old Catani subtr !!    l12 = log(two*scr(p(:,1),p(:,2))/musq)
!! old Catani subtr !!
!! old Catani subtr !!    if (tagab .eq. zero) then
!! old Catani subtr !!
!! old Catani subtr !!       res2 = two * res0 * ( &
!! old Catani subtr !!            Ca*Cf*(-51157/2592.0_dp - (107*pisq)/144.0_dp + (31*pisq**2)/480.0_dp + (22*l12)/3.0_dp + &
!! old Catani subtr !!            (659*zeta3)/72.0_dp) &
!! old Catani subtr !!            )
!! old Catani subtr !!
!! old Catani subtr !!    elseif (tagnab .eq. zero) then
!! old Catani subtr !!
!! old Catani subtr !!       res2 = two * res0 * ( &
!! old Catani subtr !!            Cf**2*(255/32.0_dp + (29*pisq)/24.0_dp - (11*pisq**2)/180.0_dp - (15*zeta3)/2.0_dp)  &
!! old Catani subtr !!            )
!! old Catani subtr !!
!! old Catani subtr !!    else
!! old Catani subtr !!
!! old Catani subtr !!       res2 = two * res0 * ( &
!! old Catani subtr !!            Cf**2*(255/32.0_dp + (29*pisq)/24.0_dp - (11*pisq**2)/180.0_dp - (15*zeta3)/2.0_dp) + &
!! old Catani subtr !!            Ca*Cf*(-51157/2592.0_dp - (107*pisq)/144.0_dp + (31*pisq**2)/480.0_dp + (22*l12)/3.0_dp + &
!! old Catani subtr !!            (659*zeta3)/72.0_dp) &
!! old Catani subtr !!            )
!! old Catani subtr !!
!! old Catani subtr !!    endif
!! old Catani subtr !!
!! old Catani subtr !!    return
!! old Catani subtr !!
!! old Catani subtr !!  end subroutine res_2loop_qqb
!! old Catani subtr !!
!! old Catani subtr !!  !-- q(1) + qb(2) -> X
!! old Catani subtr !!  subroutine res_2loop_2q(p,res2)
!! old Catani subtr !!    real(dp), intent(in) :: p(4,4)
!! old Catani subtr !!    real(dp), intent(out) :: res2(2)
!! old Catani subtr !!    real(dp) :: Lmu,res0(2),s12
!! old Catani subtr !!    s12 = two*scr(p(:,1),p(:,2))
!! old Catani subtr !!
!! old Catani subtr !!    Lmu = log(musq/s12)
!! old Catani subtr !!    call restree(p,res0)
!! old Catani subtr !!
!! old Catani subtr !!    res2 = res0* two * Cf*(4085/1296.0_dp + (4*Lmu)/3.0_dp + (7*pisq)/72.0_dp - &
!! old Catani subtr !!           zeta3/36.0_dp)
!! old Catani subtr !!
!! old Catani subtr !!    return
!! old Catani subtr !!
!! old Catani subtr !!  end subroutine res_2loop_2q

  !-------------------------------------------------------------
  !-- Amplitudes with extra jet
  !-------------------------------------------------------------

  !-- q(1) + qb(2) -> g(3) + X
  subroutine res_1loop_j_qqb(p,musq,res_loop_j)
    real(dp), intent(in) :: p(4,5),musq
    real(dp), intent(out) :: res_loop_j(2)
    complex(dp) :: za(5,5),zb(5,5),coupl(-1:1,-1:1,1:2)
    real(dp) :: l12,l14,l24
    real(dp) :: sprod(5,5),res0(2,1:2),res1(2,1:2)
    integer :: part_conf(5,2),i
    real(dp) :: propsq_w,s123,prefac
    real(dp), parameter :: avecol_qq = aveqq * 8.0_dp * Cf * xn 


    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)

    propsq_w = one/((two*scr(p(:,4),p(:,5))-mwsq)**2+mwsq*gaw**2)

    call get_coupl_w(sprod(4,5),coupl)

!    prefac = propsq_w * gwsq**2 * four
!    prefac = abs(coupl(-1,-1,1))**2
    s123 = sprod(1,2)+sprod(1,3)+sprod(2,3)

    ! ---

    prefac =  aveqq * two * s123**2 !-- extra terms in MCFM 

!!old    !-- amplitude is for 0 -> qb(1) q(2) [V --> l(3) lb(4)] g(5)
!!old    part_conf(:,1) = (/1,2,4,5,3/) !-- q qx -> g
!!old    part_conf(:,2) = (/2,1,4,5,3/) !-- qx q -> g

    !-- amplitude is for 0 -> qb(1) g(2) q(3) [V --> l(4) lb(5)]
    !-- different notation wrt below because of different hel assignment
    part_conf(1:5,1) = (/1,3,2,4,5/) !-- q qx -> g 
    part_conf(1:5,2) = (/2,3,1,4,5/) !-- qx q -> g

    do i = 1,2
!       call virt5_VH(part_conf(:,i),za,zb,sprod,res0(i),res1(i))  
       call me2_1loop_j(musq,coupl,part_conf(1,i),part_conf(2,i),part_conf(3,i),part_conf(4,i),part_conf(5,i),za,zb,sprod,res0(i,:),res1(i,:))
    enddo
    !-- remember to remove the FDH -> CDR in virt5_VH for this check
    !print *, 'qqb',res1(1)/Ca !-- MCFM check
    !print *, 'qbq',res1(2)/Ca

    !-- check against madloop
#if (_checkmadloop == 1)
    print *, 'res1/res0',res1(:)/res0(:)
!    print *, 'res0_predec',res0(1) * propsq_vh * propsq_w * gwsq**2 * (four*gwsq*mwsq) * aveqq * two * s123**2 * 1.2177157847767197_dp**(2)
#endif


    !-- go to res1_fin normalization
    l12 = log(sprod(1,2)/musq)  !s12
    l14 = log(-sprod(1,3)/musq) !s13
    l24 = log(-sprod(2,3)/musq) !s23

    res1 = res1 + res0 * (half * (Ca-two*Cf)*(pisq+three*l12-l12**2) + &
         Ca/four * (l14*(-three+two*l14)+l24*(-three+two*l24)))


! added b0 term 13/12/2019                                                                                                                                                                               
    res1 = res1*tag_ns_gg - half * (11.0_dp/6.0_dp*Ca*tag_ns_gg - two/three*TR*nf*tag_ns_nf) * res0 * (l14+l24)

    ! !-- to cross-check the asignment against tree-level [should be done with Z coupling in]
    ! print *, 'res0(1)', res0(1) * prefac ! --> res_tree_j(1)
    ! print *, 'res0(4)', res0(2) * prefac ! --> res_tree_j(4)

    !--d ux > g em vex, s cx > g em vex
    res_loop_j(1) = res1(1,1) * avecol_qq


    !--ux d > g em vex, cx s > g em vex MIRROR
    res_loop_j(2) = res1(2,1) * avecol_qq




    return

  end subroutine res_1loop_j_qqb

  !-- q/qb(1) + g(2) -> q/qb(3) + X
  subroutine res_1loop_j_qg(p,musq,res_loop_j)
    real(dp), intent(in) :: p(4,5),musq
    real(dp), intent(out) :: res_loop_j(2)
    complex(dp) :: za(5,5),zb(5,5),coupl(-1:1,-1:1,1:2)
    real(dp) :: sprod(5,5),res0(2,1:2),res1(2,1:2)
    integer :: part_conf(5,2),i
    real(dp) :: l12,l14,l24
    real(dp) :: propsq_w,s123,prefac
    real(dp), parameter :: avecol_qg = aveqg * 8.0_dp * Cf * xn 


    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)

    propsq_w = one/((two*scr(p(:,4),p(:,5))-mwsq)**2+mwsq*gaw**2)

    s123 = sprod(1,2)+sprod(1,3)+sprod(2,3)

    prefac = propsq_w * gwsq**2 * four


    prefac = prefac * aveqg * two * s123**2 !-- extra terms in MCFM

    call get_coupl_w(sprod(4,5),coupl)

    !-- amplitude is for 0 -> qb(1) g(2) q(3) [V --> l(4) lb(5)]
    !-- different notation wrt below because of different hel assignment
    part_conf(1:5,1) = (/1,2,3,4,5/) !-- q g -> q
    part_conf(1:5,2) = (/3,2,1,4,5/) !-- qx g -> qx

!!old    !-- amplitude is for 0 -> qb(1) q(2) [V --> l(3) lb(4)] g(5)
!!old    part_conf(:,1) = (/1,3,4,5,2/) !-- q g -> q
!!old    part_conf(:,2) = (/3,1,4,5,2/) !-- qx g -> qx

    do i = 1,2
!       call virt5_VH(part_conf(:,i),za,zb,sprod,res0(i),res1(i))
       call me2_1loop_j(musq,coupl,part_conf(1,i),part_conf(2,i),part_conf(3,i),part_conf(4,i),part_conf(5,i),za,zb,sprod,res0(i,:),res1(i,:))
    enddo
    !print *, 'qg',res1(1)/Ca !-- check against MCFM
    !print *, 'qbg',res1(2)/Ca


    !-- check against madloop
#if (_checkmadloop == 1)
    print *, 'res1/res0',res1(:)/res0(:)
!    print *, 'res0_predec',res0(1) * propsq_vh * propsq_w * gwsq**2 * (four*gwsq*mwsq) * aveqg * two * s123**2 * 1.2177157847767197_dp**(2)
!    print *, 'res0_predec',res0(2) * propsq_vh * propsq_w * gwsq**2 * (four*gwsq*mwsq) * aveqg * two * s123**2 * 1.2177157847767197_dp**(2)
#endif

   !-- go to res1_fin normalization                                                                                                                                                                         
    l12 = log(+sprod(1,2)/musq) !s13
    l14 = log(-sprod(1,3)/musq)  !s12
    l24 = log(-sprod(2,3)/musq) !s23

    res1 = res1 + res0 * (half * (Ca-two*Cf)*(three*l14-l14**2) + &
         Ca/four * (l12*(-three+two*l12)+l24*(-three+two*l24)-two*pisq))

    ! added b0 term 13/12/2019
    res1 = res1 - half * b0 * res0 * (l12+l24)

    ! !-- to cross-check the asignment against tree-level [should be done with Z coupling in]
    ! print *, 'res0(5)', res0(1) * prefac ! --> res_tree_j(5)
    ! print *, 'res0(6)', res0(2) * prefac ! --> res_tree_j(6)


    
    !--d g > u em vex, s g > c em vex MIRROR
    res_loop_j(1) = res1(1,1) * avecol_qg

    !--ux g > dx em vex, cx g > sx em vex MIRROR
    res_loop_j(2) = res1(2,1) * avecol_qg


    return

  end subroutine res_1loop_j_qg

  !-- g(1) + q/qb(2) -> g(3) + X
  subroutine res_1loop_j_gq(p,musq,res_loop_j)
    real(dp), intent(in) :: p(4,5),musq
    real(dp), intent(out) :: res_loop_j(2)
    complex(dp) :: za(5,5),zb(5,5),coupl(-1:1,-1:1,1:2)
    real(dp) :: sprod(5,5),res0(2,1:2),res1(2,1:2)
    integer :: part_conf(5,2),i
    real(dp) :: l12,l14,l24
    real(dp) :: propsq_w,s123,prefac
    real(dp), parameter :: avecol_gq = aveqg * 8.0_dp * Cf * xn

    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)

    propsq_w = one/((two*scr(p(:,4),p(:,5))-mwsq)**2+mwsq*gaw**2)

    s123 = sprod(1,2)+sprod(1,3)+sprod(2,3)

    prefac = propsq_w * gwsq**2 *  four

    call get_coupl_w(sprod(4,5),coupl)

    ! ---

    prefac = prefac * aveqg * two * s123**2 !-- extra terms in MCFM

!!old    !-- amplitude is for 0 -> qb(1) q(2) [V --> l(3) lb(4)] g(5)
!!old    part_conf(:,1) = (/2,3,4,5,1/) !-- g q -> q
!!old    part_conf(:,2) = (/3,2,4,5,1/) !-- g qx -> qx


    !-- amplitude is for 0 -> qb(1) g(2) q(3) [V --> l(4) lb(5)]
    !-- different notation wrt below because of different hel assignment
    part_conf(1:5,1) = (/2,1,3,4,5/) !-- g q -> q
    part_conf(1:5,2) = (/3,1,2,4,5/) !-- g qx -> qx   

    do i = 1,2
!       call virt5_VH(part_conf(:,i),za,zb,sprod,res0(i),res1(i))
       call me2_1loop_j(musq,coupl,part_conf(1,i),part_conf(2,i),part_conf(3,i),part_conf(4,i),part_conf(5,i),za,zb,sprod,res0(i,:),res1(i,:))
    enddo
    !print *, 'gq',res1(1)/Ca !-- check against MCFM
    !print *, 'gqb',res1(2)/Ca


    !-- check against madloop
#if (_checkmadloop == 1)
    print *, 'res1/res0',res1(:)/res0(:)
!    print *, 'res0_predec',res0(1) * propsq_vh * propsq_w * gwsq**2 * (four*gwsq*mwsq) * aveqg * two * s123**2 * 1.2177157847767197_dp**(2)
#endif

    !-- go to res1_fin normalization
    l24 = log(-sprod(2,3)/musq)  !s12
    l12 = log(+sprod(1,2)/musq) !s13
    l14 = log(-sprod(1,3)/musq) !s23

    res1 = res1 + res0 * (half * (Ca-two*Cf)*(three*l24-l24**2) + &
         Ca/four * (l12*(-three+two*l12)+l14*(-three+two*l14)-two*pisq))

    ! added b0 term 23/12/2019
    res1 = res1 - half * b0 * res0 * (l12+l14)


    ! !-- to cross-check the asignment against tree-level [should be done with Z coupling in]
    ! print *, 'res0(2)', res0(1) * prefac ! --> res_tree_j(2)
    ! print *, 'res0(3)', res0(2) * prefac ! --> res_tree_j(3)

    
    !--g d > u em vex, g s > c em vex
    res_loop_j(1) = res1(1,1) * avecol_gq

    !--g ux > dx em vex, g cx > sx em vex
    res_loop_j(2) = res1(2,1) * avecol_gq


    return

  end subroutine res_1loop_j_gq

  !-- master amplitude routine

  !-- 0 -> q(1) g(2) qb(3) [V --> l(4) lb(5)], MCFM normalization                                                                                                                                            
  subroutine me2_1loop_j(musq,coupl,i1,i2,i3,i4,i5,za,zb,sprod,res0,res1)
    complex(dp), intent(in) :: coupl(-1:1,-1:1,1:2),za(5,5),zb(5,5)
    integer, intent(in) :: i1,i2,i3,i4,i5
    real(dp), intent(in) :: musq,sprod(5,5)
    real(dp), intent(out) :: res0(1:2),res1(1:2)
    complex(dp) :: a0(-1:1,-1:1,-1:1),a1(-1:1,-1:1,-1:1)
    integer :: ii,jj,kk

    res0 = zero
    res1 = zero

    call A5NLO(musq,i1,i2,i3,i4,i5,za,zb,sprod,a0(-1,+1,-1),a1(-1,+1,-1))
    call A5NLO(musq,i3,i2,i1,i5,i4,zb,za,sprod,a0(-1,-1,-1),a1(-1,-1,-1))
    !                                                                                                                                                                                                        
    call A5NLO(musq,i1,i2,i3,i5,i4,za,zb,sprod,a0(-1,+1,+1),a1(-1,+1,+1))
    call A5NLO(musq,i3,i2,i1,i4,i5,zb,za,sprod,a0(-1,-1,+1),a1(-1,-1,+1))
    !                                                                                                                                                                                                        
    call A5NLO(musq,i3,i2,i1,i4,i5,za,zb,sprod,a0(+1,+1,-1),a1(+1,+1,-1))
    call A5NLO(musq,i1,i2,i3,i5,i4,zb,za,sprod,a0(+1,-1,-1),a1(+1,-1,-1))
    !                                                                                                                                                                                                        
    call A5NLO(musq,i3,i2,i1,i5,i4,za,zb,sprod,a0(+1,+1,+1),a1(+1,+1,+1))
    call A5NLO(musq,i1,i2,i3,i4,i5,zb,za,sprod,a0(+1,-1,+1),a1(+1,-1,+1))

    do kk = -1,1,2
    do jj = -1,1,2
    do ii = -1,1,2
     res0(:) = res0(:) + abs(a0(ii,jj,kk)*coupl(ii,kk,:))**2
       res1(:) = res1(:) + real(a0(ii,jj,kk)*coupl(ii,kk,:)*conjg(a1(ii,jj,kk)*coupl(ii,kk,:)),kind=dp)

    enddo
    enddo
    enddo

    !-- MCFM is in FDH --> convert to tHV                                                                                                                                                                    
    res1 = res1 - Cf * res0

    return

  end subroutine me2_1loop_j


  !-- amplitude routines, from https://arxiv.org/pdf/hep-ph/9708239.pdf
  !-- as coded up in MCFM
  !-- 0 -> 1q^+, 2g+, 3qb^-, 4l^-, 5lb^+ [different wrt MCFM]
  !-- 2*e^2*Qq*T*gs*[ALO + (as/4/pi) (4 pi)^ep Gamma(1+ep)*Gamma(1-ep)^2/Gamma(1-2*ep) * ANLO]
  !-- extra Nc absorbed wrt BDK and MCFM
  !-- Ti Tj --> delta_ij                                          
 subroutine A5NLO(musq,j1,j2,j3,j4,j5,za,zb,sprod,A5LOp,A5NLOp)
    integer, intent(in) :: j1,j2,j3,j4,j5
    complex(dp), intent(in) :: za(5,5),zb(5,5)
    real(dp), intent(in) :: musq,sprod(5,5)
    complex(dp), intent(out) :: A5LOp,A5NLOp

    A5LOp = -za(j3,j4)**2/(za(j1,j2)*za(j2,j3)*za(j4,j5))
    !-- note the different ordering of A52                                                
    A5NLOp = Ca*A51(musq,j1,j2,j3,j4,j5,za,zb,sprod) + (Ca-two*Cf)*A52(musq,j1,j3,j2,j4,j5,za,zb,sprod)

    return

  end subroutine A5NLO

  !-- 0 -> 1q^+, 2g+, 3qb^-, 4l^-, 5lb^+ [different l/lb wrt BDK because outgoing leptons] 
  function A51(musq,j1,j2,j3,j4,j5,za,zb,sprod)
    complex(dp) :: A51
    integer, intent(in) :: j1,j2,j3,j4,j5
    complex(dp), intent(in) :: za(5,5),zb(5,5)
    real(dp), intent(in) :: musq,sprod(5,5)
    complex(dp) :: A5lom,Vcc,Fcc,Vsc,Fsc
    complex(dp) :: l12,l23

  !-- -i * A5tree
    A5lom =-za(j3,j4)**2/(za(j1,j2)*za(j2,j3)*za(j4,j5))
    l12=lnrat(musq,-sprod(j1,j2))
    l23=lnrat(musq,-sprod(j2,j3))

    !-- leading N, without poles
    Vcc = -half * (l12**2+l23**2)-two*l23-four

    Fcc=za(j3,j4)**2/(za(j1,j2)*za(j2,j3)*za(j4,j5)) &
         *(Lsm1(-sprod(j1,j2),-sprod(j4,j5),-sprod(j2,j3),-sprod(j4,j5)) &
         -two*za(j3,j1)*zb(j1,j5)*za(j5,j4)/za(j3,j4) &
         *L0(-sprod(j2,j3),-sprod(j4,j5))/sprod(j4,j5))

    Vsc = half*l23 + one

    Fsc =za(j3,j4)*za(j3,j1)*zb(j1,j5)*za(j5,j4) &
         /(za(j1,j2)*za(j2,j3)*za(j4,j5))*L0(-sprod(j2,j3),-sprod(j4,j5))/sprod(j4,j5) &
         +half*(za(j3,j1)*zb(j1,j5))**2*za(j4,j5) &
         /(za(j1,j2)*za(j2,j3))*L1(-sprod(j2,j3),-sprod(j4,j5))/sprod(j4,j5)**2

    A51=(Vcc+Vsc)*A5Lom+Fcc+Fsc

    return

  end function A51

  function A52(musq,j1,j2,j3,j4,j5,za,zb,sprod)
    complex(dp) :: A52
    integer, intent(in) :: j1,j2,j3,j4,j5
    complex(dp), intent(in) :: za(5,5),zb(5,5)
    real(dp), intent(in) :: musq,sprod(5,5)
    complex(dp) :: l12,l45,A5lom
    complex(dp) :: Vcc,Vcctest,Fcc,Vsc,Vsctest,Fsc

    l12=lnrat(musq,-sprod(j1,j2))
    l45=lnrat(musq,-sprod(j4,j5))

    !-- -i * A5tree
    A5lom=za(j2,j4)**2/(za(j2,j3)*za(j3,j1)*za(j4,j5))

    Vcc = -half*l12**2-two*l45-four

    !--subleading N, without poles
    Fcc=-za(j2,j4)**2/(za(j2,j3)*za(j3,j1)*za(j4,j5)) &
         *Lsm1(-sprod(j1,j2),-sprod(j4,j5),-sprod(j1,j3), -sprod(j4,j5)) &
         +za(j2,j4)*(za(j1,j2)*za(j3,j4)-za(j1,j4)*za(j2,j3)) &
         /(za(j2,j3)*za(j1,j3)**2*za(j4,j5)) &
         *Lsm1(-sprod(j1,j2),-sprod(j4,j5),-sprod(j2,j3),-sprod(j4,j5)) &
         +two*zb(j1,j3)*za(j1,j4)*za(j2,j4)/(za(j1,j3)*za(j4,j5)) &
         *L0(-sprod(j2,j3),-sprod(j4,j5))/sprod(j4,j5)

    Vsc = half*l45 + half

    Fsc=za(j1,j4)**2*za(j2,j3)/(za(j1,j3)**3*za(j4,j5)) &
         *Lsm1(-sprod(j1,j2),-sprod(j4,j5),-sprod(j2,j3),-sprod(j4,j5)) &
         -half*(za(j4,j1)*zb(j1,j3))**2*za(j2,j3)/(za(j1,j3)*za(j4,j5)) &
         *L1(-sprod(j4,j5),-sprod(j2,j3))/sprod(j2,j3)**2 &
         +za(j1,j4)**2*za(j2,j3)*zb(j3,j1)/(za(j1,j3)**2*za(j4,j5)) &
         *L0(-sprod(j4,j5),-sprod(j2,j3))/sprod(j2,j3) &
         -za(j2,j1)*zb(j1,j3)*za(j4,j3)*zb(j3,j5)/za(j1,j3) &
         *L1(-sprod(j4,j5),-sprod(j1,j2))/sprod(j1,j2)**2 &
         -za(j2,j1)*zb(j1,j3)*za(j3,j4)*za(j1,j4)/(za(j1,j3)**2*za(j4,j5)) &
         *L0(-sprod(j4,j5),-sprod(j1,j2))/sprod(j1,j2) &
         -half*zb(j3,j5)*(zb(j1,j3)*zb(j2,j5)+zb(j2,j3)*zb(j1,j5)) &
         /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*zb(j4,j5))

    A52=(Vcc+Vsc)*A5lom+Fcc+Fsc

    return

  end function A52





!!old VH code !!  !------------------------------------------------------------
!!old VH code !!
!!old VH code !!  !-- amplitude routines, as coded up in MCFMx
!!old VH code !!  !-- coefficient of asontwopi
!!old VH code !!
!!old VH code !!  !-- 0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_L/R(5)
!!old VH code !!  !-- does not contain couplings and propagators
!!old VH code !!  subroutine virt5_VH(jj,za,zb,sprod,res0,res1)
!!old VH code !!    integer, intent(in) :: jj(5)
!!old VH code !!    complex(dp), intent(in) :: za(5,5),zb(5,5)
!!old VH code !!    real(dp), intent(in) :: sprod(5,5)
!!old VH code !!    real(dp), intent(out) :: res0,res1
!!old VH code !!    complex(dp) :: A5LOm,A5LOp,A5NLOm,A5NLOp
!!old VH code !!
!!old VH code !!    !--  0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_L(5)
!!old VH code !!    call A5NLO_VH(jj(1),jj(2),jj(3),jj(4),jj(5),za,zb,sprod,A5LOm,A5NLOm)
!!old VH code !!
!!old VH code !!    !-- 0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_R(5)
!!old VH code !!    call A5NLO_VH(jj(2),jj(1),jj(4),jj(3),jj(5),zb,za,sprod,A5LOp,A5NLOp)
!!old VH code !!
!!old VH code !!    res0 = abs(A5LOp)**2 + abs(A5LOm)**2
!!old VH code !!    res1 = real(conjg(A5LOp)*A5NLOp,kind=dp) + real(conjg(A5LOm)*A5NLOm,kind=dp)
!!old VH code !!
!!old VH code !!    !-- MCFM is in FDH --> convert to tHV
!!old VH code !!    res1 = res1 - Cf * res0
!!old VH code !!    
!!old VH code !!    return
!!old VH code !!
!!old VH code !!  end subroutine virt5_VH
!!old VH code !!
!!old VH code !!
!!old VH code !!  !--  0--> qb_R(1)+q_L(2)+l_L(3)+lb_R(4)+g_L(5)
!!old VH code !!  subroutine A5NLO_VH(j1,j2,j3,j4,j5,za,zb,sprod,A5LOm,A5NLOm)
!!old VH code !!    integer, intent(in) :: j1,j2,j3,j4,j5
!!old VH code !!    complex(dp), intent(in) :: za(5,5),zb(5,5)
!!old VH code !!    real(dp), intent(in) :: sprod(5,5)
!!old VH code !!    complex(dp), intent(out) :: A5LOm,A5NLOm
!!old VH code !!    complex(dp) :: zab2
!!old VH code !!    real(dp) :: s125
!!old VH code !!
!!old VH code !!    zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
!!old VH code !!
!!old VH code !!    s125=sprod(j1,j2)+sprod(j1,j5)+sprod(j2,j5)
!!old VH code !!
!!old VH code !!    A5LOm=-zb(j1,j4)*zab2(j3,j2,j5,j1)/(zb(j2,j5)*zb(j5,j1)*s125)
!!old VH code !!    A5NLOm=Ca * A51_VH(j2,j5,j1,j4,j3,zb,za,sprod) + (Ca-two*Cf)*A52_VH(j2,j1,j5,j4,j3,zb,za,sprod)
!!old VH code !!
!!old VH code !!    return
!!old VH code !!
!!old VH code !!  end subroutine A5NLO_VH
!!old VH code !!
!!old VH code !!  function A51_VH(j1,j2,j3,j4,j5,za,zb,sprod)
!!old VH code !!    complex(dp) :: A51_VH
!!old VH code !!    integer, intent(in) :: j1,j2,j3,j4,j5
!!old VH code !!    complex(dp), intent(in) :: za(5,5),zb(5,5)
!!old VH code !!    real(dp), intent(in) :: sprod(5,5)
!!old VH code !!    complex(dp) :: zab2,zab3,l12,l23
!!old VH code !!    complex(dp) :: A5lom,Vcc,Fcc,Vsc,Fsc,MPC
!!old VH code !!    real(dp) :: s123
!!old VH code !!
!!old VH code !!    zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
!!old VH code !!    zab3(j1,j2,j3,j4,j5)=za(j1,j2)*zb(j2,j5)+za(j1,j3)*zb(j3,j5)+za(j1,j4)*zb(j4,j5)
!!old VH code !!
!!old VH code !!    s123=sprod(j1,j2)+sprod(j2,j3)+sprod(j3,j1)
!!old VH code !!
!!old VH code !!    A5lom=-za(j3,j4)*zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123)
!!old VH code !!
!!old VH code !!    l12=lnrat(musq,-sprod(j1,j2))
!!old VH code !!    l23=lnrat(musq,-sprod(j2,j3))
!!old VH code !!
!!old VH code !!    Vcc = -half * (l12**2+l23**2) -two*l23-four      
!!old VH code !!
!!old VH code !!    Fcc = zab2(j3,j1,j2,j5)/(za(j1,j2)*za(j2,j3)*s123) &
!!old VH code !!         *(za(j3,j4)*Lsm1(-sprod(j1,j2),-s123,-sprod(j2,j3),-s123) &
!!old VH code !!         +two*za(j3,j1)*zab2(j4,j2,j3,j1) &
!!old VH code !!         *L0(-sprod(j2,j3),-s123)/s123)
!!old VH code !!
!!old VH code !!    Vsc =half*l23+one
!!old VH code !!
!!old VH code !!    Fsc =-za(j3,j4)*za(j3,j1)*zb(j1,j5) &
!!old VH code !!         /(za(j1,j2)*za(j2,j3))*L0(-sprod(j2,j3),-s123)/s123 &
!!old VH code !!         +half*za(j3,j1)**2*zb(j1,j5)*zab2(j4,j2,j3,j1) &
!!old VH code !!         /(za(j1,j2)*za(j2,j3))*L1(-sprod(j2,j3),-s123)/s123**2
!!old VH code !!
!!old VH code !!    MPC=-(3._dp*L0(-sprod(j2,j3),-s123)*za(j1,j3)* &
!!old VH code !!         zab3(j4,j1,j2,j3,j5)*zb(j2,j1))/(two*s123**2*za(j1,j2))
!!old VH code !!
!!old VH code !!    A51_VH=(Vcc+Vsc)*A5Lom+Fcc+Fsc+MPC
!!old VH code !!
!!old VH code !!    return
!!old VH code !!
!!old VH code !!  end function A51_VH
!!old VH code !!  
!!old VH code !!  function A52_VH(j1,j2,j3,j4,j5,za,zb,sprod)
!!old VH code !!    complex(dp) :: A52_VH
!!old VH code !!    integer, intent(in) :: j1,j2,j3,j4,j5
!!old VH code !!    complex(dp), intent(in) :: za(5,5),zb(5,5)
!!old VH code !!    real(dp), intent(in) :: sprod(5,5)
!!old VH code !!    complex(dp) :: zab2,l12,l45
!!old VH code !!    complex(dp) :: A5lom,Vcc,Fcc,Vsc,Fsc,TCor
!!old VH code !!    real(dp) :: s123
!!old VH code !!
!!old VH code !!    zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
!!old VH code !!
!!old VH code !!    s123=sprod(j1,j2)+sprod(j2,j3)+sprod(j3,j1)
!!old VH code !!
!!old VH code !!    l12=lnrat(musq,-sprod(j1,j2))  
!!old VH code !!    l45=lnrat(musq,-s123)  
!!old VH code !!
!!old VH code !!    !    -i * A5tree  
!!old VH code !!    A5lom=za(j2,j4)*zab2(j2,j1,j3,j5)/(za(j2,j3)*za(j3,j1)*s123)
!!old VH code !!
!!old VH code !!    Vcc=-(half*l12**2) - two*l45 - four
!!old VH code !!
!!old VH code !!    Fcc=-za(j2,j4)*zab2(j2,j1,j3,j5)/(za(j2,j3)*za(j3,j1)*s123) &
!!old VH code !!         *Lsm1(-sprod(j1,j2),-s123,-sprod(j1,j3), -s123) &
!!old VH code !!         +zab2(j2,j1,j3,j5)*(za(j1,j2)*za(j3,j4)-za(j1,j4)*za(j2,j3)) &
!!old VH code !!         /(za(j2,j3)*za(j1,j3)**2*s123) &
!!old VH code !!         *Lsm1(-sprod(j1,j2),-s123,-sprod(j2,j3),-s123) &
!!old VH code !!         +two*zb(j1,j3)*za(j1,j4)*zab2(j2,j1,j3,j5)/(za(j1,j3)*s123) &
!!old VH code !!         *L0(-sprod(j2,j3),-s123)/s123
!!old VH code !!
!!old VH code !!    Vsc=half*l45+half
!!old VH code !!    
!!old VH code !!    Fsc=za(j1,j4)*za(j2,j3)*zab2(j1,j2,j3,j5)/(za(j1,j3)**3*s123) &
!!old VH code !!         *Lsm1(-sprod(j1,j2),-s123,-sprod(j2,j3),-s123) &
!!old VH code !!         +half*za(j4,j1)*zb(j1,j3)**2*za(j2,j3)*zab2(j1,j2,j3,j5) &
!!old VH code !!         /(za(j1,j3)*s123)*L1(-s123,-sprod(j2,j3))/sprod(j2,j3)**2 &
!!old VH code !!         +za(j1,j4)*za(j2,j3)*zb(j3,j1)*zab2(j1,j2,j3,j5) &
!!old VH code !!         /(za(j1,j3)**2*s123)*L0(-s123,-sprod(j2,j3))/sprod(j2,j3) &
!!old VH code !!         -za(j2,j1)*zb(j1,j3)*za(j4,j3)*zb(j3,j5)/za(j1,j3) &
!!old VH code !!         *L1(-s123,-sprod(j1,j2))/sprod(j1,j2)**2 &
!!old VH code !!         -za(j2,j1)*zb(j1,j3)*za(j3,j4)*zab2(j1,j2,j3,j5) &
!!old VH code !!         /(za(j1,j3)**2*s123)*L0(-s123,-sprod(j1,j2))/sprod(j1,j2) &
!!old VH code !!         +half*zab2(j4,j1,j2,j3)*(zb(j1,j3)*zb(j2,j5)+zb(j2,j3)*zb(j1,j5)) &
!!old VH code !!         /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*s123) 
!!old VH code !!
!!old VH code !!    TCor=   (Lsm1(-sprod(j1,j2),-s123,-sprod(j2,j3),-s123)*za(j1,j2)*  &
!!old VH code !!         (za(j4,j1)*zb(j1,j5) + za(j4,j2)*zb(j2,j5) +  &
!!old VH code !!         za(j4,j3)*zb(j3,j5)))/(s123*za(j1,j3)**2) +  &
!!old VH code !!         ((-(L0(-s123,-sprod(j1,j2))/sprod(j1,j2)) +  &
!!old VH code !!         L0(-s123,-sprod(j2,j3))/sprod(j2,j3))*za(j1,j2)* &
!!old VH code !!         zb(j3,j1)*(za(j4,j1)*zb(j1,j5) +  &
!!old VH code !!         za(j4,j2)*zb(j2,j5) + za(j4,j3)*zb(j3,j5)))/ &
!!old VH code !!         (s123*za(j1,j3)) -  &
!!old VH code !!         (zb(j3,j1)*(za(j4,j1)*zb(j1,j5) +  &
!!old VH code !!         za(j4,j2)*zb(j2,j5) + za(j4,j3)*zb(j3,j5)))/ &
!!old VH code !!         (s123*za(j1,j3)*zb(j2,j1)) 
!!old VH code !!     
!!old VH code !!    A52_VH=(Vcc+Vsc)*A5lom+Fcc+Fsc+TCor
!!old VH code !!
!!old VH code !!    return
!!old VH code !!
!!old VH code !!  end function A52_VH
  



end module mod_amplitudes_loop_wp

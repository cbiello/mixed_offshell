!-- amplitudes with extra jet: factor our gs**2 per extra emission
!-- loop amplitudes: factor out asontwopi per loop
!-- finite part according to Catani

!-- if we include the 3rd generation, gg triangles should be added to 1-loop squared

module mod_amplitudes_loop_az
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_amplitudes_tree_az
  implicit none
  private

  public :: res_1loop
  public :: res_1loop_sq, res_2loop
  public :: res_1loop_j_ns,res_1loop_j_qg,res_1loop_j_gq

contains

#include "./inc_coupl_az.f90"

  !-- q(1) + qb(2) -> V(3+4)
  !-- summed / averaged
  subroutine res_1loop(p,musq,res0,res1)
    real(dp), intent(in) :: p(4,4),musq
    real(dp), intent(out) :: res0(4),res1(4)

    call restree(p,res0)

    !-- wilson coefficient 
    res1 = -8*Cf*res0

    !-- scale dependent term
    
    return

  end subroutine res_1loop

  !-- q(1) + qb(2) -> V(3+4)
  !-- square of the finite part, coefficient of asontwopi**2
  subroutine res_1loop_sq(p,musq,res0,res1)
    real(dp), intent(in) :: p(4,4),musq
    real(dp), intent(out) :: res0(4),res1(4)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4),zb(4,4)

    call restree(p,res0)

    res1 = res0 * 16.0_dp * Cf**2 * tag_ns_gg

    return

  end subroutine res_1loop_sq

  !-- q(1) + qb(2) -> V(3+4)
  !-- finite part, coefficient of asontwopi**2
  subroutine res_2loop(p,musq,res0,res2)
    real(dp), intent(in) :: p(4,4),musq
    real(dp), intent(out) :: res0(4),res2(4)
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

  !-------------------------------------------------------------
  !-- Amplitudes with extra jet
  !-------------------------------------------------------------

  subroutine res_1loop_j_ns(p,musq,res_j)
    real(dp), intent(in) :: p(4,5),musq
    real(dp), intent(out) :: res_j(4)
    real(dp) :: sprod(5,5)
    complex(dp) :: za(5,5),zb(5,5),coupl(-1:1,-1:1,1:2)
    integer :: part_conf(3,2),i
    complex(dp) :: amp0(-1:1,-1:1,-1:1,1:2),amp1(-1:1,-1:1,-1:1,1:2)
    real(dp) :: res0(2,1:2),res1(2,1:2),l12,l14,l24
    real(dp), parameter :: avecol_qq = aveqq * 8.0_dp * Cf * xn 

    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)

    call get_coupl_az(sprod(4,5),coupl)

    !-- amplitude is for 0 -> qb(1) g(2) q(3) [V --> l(4) lb(5)]
    !-- different notation wrt below because of different hel assignment
    part_conf(1:3,1) = (/1,3,2/) !-- q qx -> g
    part_conf(1:3,2) = (/2,3,1/) !-- qx q -> g

    do i = 1,2
       
       call me2_1loop_j(musq,coupl,part_conf(1,i),part_conf(2,i),part_conf(3,i),4,5,za,zb,sprod,res0(i,:),res1(i,:)) !-- returns the amplitude in MCFM normalization

    enddo

    !-- go to res1_fin normalization
    l12 = log(sprod(1,2)/musq)  !s12
    l14 = log(-sprod(1,3)/musq) !s13
    l24 = log(-sprod(2,3)/musq) !s23

    res1 = res1 + res0 * (half * (Ca-two*Cf)*(pisq+three*l12-l12**2) + &
         Ca/four * (l14*(-three+two*l14)+l24*(-three+two*l24)))
    
!    res1 = res1 - half * b0 * res0 * (l14+l24)
    res1 = res1*tag_ns_gg - half * (11.0_dp/6.0_dp*Ca*tag_ns_gg - two/three*TR*nf*tag_ns_nf) * res0 * (l14+l24)
    

    !-- to cross-check the asignment against tree-level [should be done with Z coupling in]
    !print *, 'res0(1)', res0(1,2) * avecol_qq ! --> res_tree_j(1)
    !print *, 'res0(2)', res0(1,1) * avecol_qq ! --> res_tree_j(2)
    !print *, 'res0(3)', res0(2,2) * avecol_qq ! --> res_tree_j(7)
    !print *, 'res0(4)', res0(2,1) * avecol_qq ! --> res_tree_j(8)

    !--u ux > g em ep, c cx > g em ep
    res_j(1) = res1(1,2) * avecol_qq

    !--d dx > g em ep, s sx > g em ep
    res_j(2) = res1(1,1) * avecol_qq

    !--ux u > g em ep, cx c > g em ep MIRROR
    res_j(3) = res1(2,2) * avecol_qq

    !--dx d > g em ep, sx s > g em ep MIRROR
    res_j(4) = res1(2,1) * avecol_qq

    return

  end subroutine res_1loop_j_ns

  subroutine res_1loop_j_qg(p,musq,res_j)
    real(dp), intent(in) :: p(4,5),musq
    real(dp), intent(out) :: res_j(4)
    real(dp) :: sprod(5,5)
    complex(dp) :: za(5,5),zb(5,5),coupl(-1:1,-1:1,1:2)
    integer :: part_conf(3,2),i
    complex(dp) :: amp0(-1:1,-1:1,-1:1,1:2),amp1(-1:1,-1:1,-1:1,1:2)
    real(dp) :: res0(2,1:2),res1(2,1:2),l12,l14,l24
    real(dp), parameter :: avecol_qg = aveqg * 8.0_dp * Cf * xn 

    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)

    call get_coupl_az(sprod(4,5),coupl)

    !-- amplitude is for 0 -> qb(1) g(2) q(3) [V --> l(4) lb(5)]
    !-- different notation wrt below because of different hel assignment
    part_conf(1:3,1) = (/1,2,3/) !-- q g -> q
    part_conf(1:3,2) = (/3,2,1/) !-- qx g -> qx

    do i = 1,2
       
       call me2_1loop_j(musq,coupl,part_conf(1,i),part_conf(2,i),part_conf(3,i),4,5,za,zb,sprod,res0(i,:),res1(i,:)) !-- returns the amplitude in MCFM normalization

    enddo

    !-- go to res1_fin normalization
    l12 = log(sprod(1,2)/musq)  !s12
    l14 = log(-sprod(1,3)/musq) !s13
    l24 = log(-sprod(2,3)/musq) !s23

    res1 = res1 + res0 * (half * (Ca-two*Cf)*(three*l14-l14**2) + &
         Ca/four * (l12*(-three+two*l12)+l24*(-three+two*l24)-two*pisq))
    res1 = res1 - half * b0 * res0 * (l12+l24)
    
    !-- to cross-check the asignment against tree-level [should be done with Z coupling in]
    ! print *, 'res0(1)', res0(1,2) * avecol_qg ! --> res_tree_j(9)
    ! print *, 'res0(2)', res0(1,1) * avecol_qg ! --> res_tree_j(10)
    ! print *, 'res0(3)', res0(2,2) * avecol_qg ! --> res_tree_j(11)
    ! print *, 'res0(4)', res0(2,1) * avecol_qg ! --> res_tree_j(12)

    !--u g > u em ep, c g > c em ep MIRROR
    res_j(1) = res1(1,2) * avecol_qg

    !--d g > d em ep, s g > s em ep MIRROR
    res_j(2) = res1(1,1) * avecol_qg

    !--ux g > ux em ep, cx g > cx em ep MIRROR
    res_j(3) = res1(2,2) * avecol_qg

    !--dx g > dx em ep, sx g > sx em ep MIRROR
    res_j(4) = res1(2,1) * avecol_qg

    return

  end subroutine res_1loop_j_qg

  subroutine res_1loop_j_gq(p,musq,res_j)
    real(dp), intent(in) :: p(4,5),musq
    real(dp), intent(out) :: res_j(4)
    real(dp) :: sprod(5,5)
    complex(dp) :: za(5,5),zb(5,5),coupl(-1:1,-1:1,1:2)
    integer :: part_conf(3,2),i
    complex(dp) :: amp0(-1:1,-1:1,-1:1,1:2),amp1(-1:1,-1:1,-1:1,1:2)
    real(dp) :: res0(2,1:2),res1(2,1:2),l12,l14,l24
    real(dp), parameter :: avecol_qg = aveqg * 8.0_dp * Cf * xn 

    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)

    call get_coupl_az(sprod(4,5),coupl)

    !-- amplitude is for 0 -> qb(1) g(2) q(3) [V --> l(4) lb(5)]
    !-- different notation wrt below because of different hel assignment
    part_conf(1:3,1) = (/2,1,3/) !-- g q -> q
    part_conf(1:3,2) = (/3,1,2/) !-- g qx -> qx

    do i = 1,2
       
       call me2_1loop_j(musq,coupl,part_conf(1,i),part_conf(2,i),part_conf(3,i),4,5,za,zb,sprod,res0(i,:),res1(i,:)) !-- returns the amplitude in MCFM normalization

    enddo

    !-- go to res1_fin normalization
    l12 = log(sprod(1,2)/musq)  !s12
    l14 = log(-sprod(1,3)/musq) !s13
    l24 = log(-sprod(2,3)/musq) !s23

    res1 = res1 + res0 * (half * (Ca-two*Cf)*(three*l24-l24**2) + &
         Ca/four * (l12*(-three+two*l12)+l14*(-three+two*l14)-two*pisq))
    res1 = res1 - half * b0 * res0 * (l12+l14)
    
    !-- to cross-check the asignment against tree-level [should be done with Z coupling in]
    !print *, 'res0(1)', res0(1,2) * avecol_qg ! --> res_tree_j(3)
    !print *, 'res0(2)', res0(1,1) * avecol_qg ! --> res_tree_j(4)
    !print *, 'res0(3)', res0(2,2) * avecol_qg ! --> res_tree_j(5)
    !print *, 'res0(4)', res0(2,1) * avecol_qg ! --> res_tree_j(6)

    !--g u > u em ep, g c > c em ep
    res_j(1) = res1(1,2) * avecol_qg

    !--g d > d em ep, g s > s em ep
    res_j(2) = res1(1,1) * avecol_qg

    !--g ux > ux em ep, g cx > cx em ep
    res_j(3) = res1(2,2) * avecol_qg

    !--g dx > dx em ep, g sx > sx em ep
    res_j(4) = res1(2,1) * avecol_qg

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

  
  ! !-- q(1) + qb(2) --> g(3) + [V --> l(4)+lb(5)]
  ! subroutine res_1loop_j(p,res1)
  !   real(dp), intent(in) :: p(4,5)
  !   real(dp), intent(out) :: res1
  !   real(dp) :: sprod(5,5)
  !   complex(dp) :: za(5,5),zb(5,5)
  !   real(dp) :: res0,l12,l14,l24

  !   res1 = zero

  !   call spinoru(5,(/-p(:,2),p(:,3),-p(:,1),p(:,4),p(:,5)/),za,zb,sprod)
  !   call me2_1loop_j(za,zb,sprod,res0,res1)

  !   !-- go to res1_fin normalization
  !   l12 = log(sprod(1,3)/musq)  !s12
  !   l14 = log(-sprod(2,3)/musq) !s13
  !   l24 = log(-sprod(1,2)/musq) !s23

  !   res1 = res1 + res0 * (half * (Ca-two*Cf)*(pisq+three*l12-l12**2) + &
  !        Ca/four * (l14*(-three+two*l14)+l24*(-three+two*l24)))
  !   res1 = res1 * aveqq

  ! end subroutine res_1loop_j


  ! !-- q(1) + g(2) --> q(3) + [V --> l(4)+lb(5)]
  ! subroutine res_1loop_j_qg(p,res1)
  !   real(dp), intent(in) :: p(4,5)
  !   real(dp), intent(out) :: res1
  !   real(dp) :: sprod(5,5)
  !   complex(dp) :: za(5,5),zb(5,5)
  !   real(dp) :: res0,l12,l14,l24

  !   res1 = zero

  !   call spinoru(5,(/p(:,3),-p(:,2),-p(:,1),p(:,4),p(:,5)/),za,zb,sprod)
  !   call me2_1loop_j(za,zb,sprod,res0,res1)

  !   !-- go to res1_fin normalization
  !   l12 = log(sprod(2,3)/musq)
  !   l14 = log(-sprod(1,3)/musq)
  !   l24 = log(-sprod(1,2)/musq)

  !   res1 = res1 + (&
  !        -((Ca-two*Cf)*(-3 + l14)*l14)/two + &
  !        (Ca*(-3*l12 + 2*l12**2 - 3*l24 + 2*l24**2 - 2*Pi**2))/four)*res0

  !   res1 = res1 * aveqg

  ! end subroutine res_1loop_j_qg

  ! subroutine res_1loop_j_gq(p,res1)

  !   real(dp), intent(in) :: p(4,5)
  !   real(dp), intent(out) :: res1
  !   real(dp) :: pswap(4,5),res

  !   pswap = p
  !   pswap(:,1) = p(:,2)
  !   pswap(:,2) = p(:,1)
  !   call res_1loop_j_qg(pswap,res)
  !   res1 = res

  ! end subroutine res_1loop_j_gq



end module mod_amplitudes_loop_az

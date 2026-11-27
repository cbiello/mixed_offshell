module mod_subtrfn_nnlo_dy_ns
! integrated subtraction terms for gg->ff channel, where ff=gg or qqb
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_tripleclimits_q_gg_qstar
  use mod_tripleclimits_q_qqb_qstar
  implicit none
  integer, parameter :: colstruc=5                              ! choose the color structure: 1=Ca^2, 2=Ca*Nf, 3=Cf*Nf, 4=Nf^2, 5=total

  integer, parameter :: ab = 1
  integer, parameter :: nab = 2
  integer, parameter :: tot = 4

  private 

  public :: qqb_ns_sub_qqbz124,qqb_ns_sub_qqbz124_pls,qqb_ns_sub_qqb124
  public :: qqb_ns_sub_qqbz1zbar2,qqb_ns_sub_qqbz12,qqb_ns_sub_qqb12
!  public :: qqb_ns_sub_qgz124

  
contains

  
  function qqb_ns_sub_qqbz124(shat,mu2,z,eta41,eta42,eta41log)
    ! Integrated subtraction terms to multiply FLMqqb[z1,2,4] etc
    ! This gives reg+pls
    real(dp), intent(in)   :: shat,mu2,z,eta41,eta42,eta41log
    real(dp)               :: qqb_ns_sub_qqbz124
    real(dp)               :: Lmu,omz,lomz,w51tc

    Lmu = log(mu2/shat)
    omz = one-z
    lomz = log(omz)
    w51tc = eta42*(one+eta41)                          ! w4151tilde51

! scale dependent pieces

    qqb_ns_sub_qqbz124 = Cf*(one + z)*Lmu
    qqb_ns_sub_qqbz124 = qqb_ns_sub_qqbz124 - two*Cf/omz*Lmu

! scale independent pieces
    qqb_ns_sub_qqbz124 = qqb_ns_sub_qqbz124 + &
         Cf*(omz - w51tc*(one+z)*log(eta41log/two) - two*(one+z)*lomz)
    qqb_ns_sub_qqbz124 = qqb_ns_sub_qqbz124 + &
         four*Cf*lomz/omz +two*Cf*w51tc*log(eta41log/two)/omz

    qqb_ns_sub_qqbz124 = qqb_ns_sub_qqbz124 * tag_ns_gg

    


  end function qqb_ns_sub_qqbz124


  
  function qqb_ns_sub_qqbz124_pls(shat,mu2,z,eta41,eta42,eta41log)
    ! Integrated subtraction terms to multiply FLM[z1,2,4],  cf. Eq. 10.10
    ! This gives [pls(z=1)]
    real(dp), intent(in)   :: shat,mu2,z,eta41,eta42,eta41log
    real(dp)               :: qqb_ns_sub_qqbz124_pls
    real(dp)               :: Lmu,omz,lomz,w51tc

    Lmu = log(mu2/shat)
    omz = one-z
    lomz = log(omz)

    w51tc = eta42*(one+eta41)                          ! w4151tilde51 

    qqb_ns_sub_qqbz124_pls = - two*Cf/omz*Lmu
    qqb_ns_sub_qqbz124_pls = qqb_ns_sub_qqbz124_pls + &
         four*Cf*lomz/omz +two*Cf*w51tc*log(eta41log/two)/omz

    qqb_ns_sub_qqbz124_pls = qqb_ns_sub_qqbz124_pls * tag_ns_gg


  end function qqb_ns_sub_qqbz124_pls

    function qqb_ns_sub_qqb124(shat,mu2,eta41,eta42,eta41log,eta42log,E4)
    ! Integrated subtraction terms to multiply FLM[1,2,4],  cf. Eq. 10.10
    ! This gives [pls(z=1)]
    real(dp), intent(in)   :: shat,mu2,eta41,eta42,eta41log,eta42log,E4
    real(dp)               :: qqb_ns_sub_qqb124
    real(dp)               :: Lmu,w4151t51,w4252t52,w4151t45,w4252t45,le4oe1

    Lmu = log(mu2/shat)
    
    w4151t51 = eta42*(one+eta41)
    w4252t52 = eta41*(one+eta42)
    !
    w4151t45 = eta42**2*(one+two*eta41)
    w4252t45 = eta41**2*(one+two*eta42)

    le4oe1 = half*log(E4**2/shat) + log2   ! log(E4/E1)

    ! log term
     qqb_ns_sub_qqb124 = -three*Cf*Lmu*tag_ns_gg

     ! scale independent term
     qqb_ns_sub_qqb124 = qqb_ns_sub_qqb124 + &
          ( tag_ns_gg*(137.0_dp*Ca - 6.0_dp*Ca*le4oe1 + 18.0_dp*Ca*le4oe1**2 + 126.0_dp*Ca*log2 + 36.0_dp*Ca*le4oe1*log2 - 24.0_dp*Ca*pisq +&
          12.0_dp*Cf*pisq  + 18.0_dp*Ca*dilog2(eta41) + 18.0_dp*Ca*dilog2(eta42) -&
          33.0_dp*Ca*w4151t45*Log(eta41log/eta42log)  - 33.0_dp*Ca*w4252t45*Log(eta42log/eta41log) + &
          30.0_dp*Ca*Log(eta41log/two) + 18.0_dp*Ca*le4oe1*Log(eta41log/two) - 36.0_dp*Cf*le4oe1*w4151t51*Log(eta41log/two) + &
          30.0_dp*Ca*Log(eta42log/two) + 18.0_dp*Ca*le4oe1*Log(eta42log/two) - 36.0_dp*Cf*le4oe1*w4252t52*Log(eta42log/two)) + &
          Nf*TR*tag_ns_nf*(-52.0_dp + 12.0_dp*le4oe1 - 36.0_dp*log2+ 12.0_dp*w4151t45*Log(eta41log/eta42log) +&
          12.0_dp*w4252t45*Log(eta42log/eta41log) -  6.0_dp*Log(eta41log/two) -6.0_dp*Log(eta42log/two) ) )/18.0_dp






   end function qqb_ns_sub_qqb124
         


  function qqb_ns_sub_qqbz1zbar2(shat,mu2,z,zbar)
! Integrated subtraction terms to multiply FLM[z1,zbar2] cf. Eq. 10.15.
! log(musq/s) evaluated for {shat,shat/zbar,shat/z,shat/z/zbar} 
! output indices: coensicient of {z1 zbar2,z12, 1z2, 12}
    implicit none
    real(dp), intent(in)  :: shat,mu2,z,zbar
    real(dp)              :: qqb_ns_sub_qqbz1zbar2(1:4)
    real(dp)              :: coeff_noL(1:4),coeff_L(1:4),coeff_Lsq(1:4)
    real(dp)              :: omz,omzbar,lomz,lomzbar,Lmu(1:4)
    integer               :: j

    Lmu(4)  = log(mu2/shat)
    Lmu(3)  = log(mu2/shat) + log(zbar) 
    Lmu(2)  = log(mu2/shat) + log(z)
    Lmu(1)  = log(mu2/shat) + log(z) + log(zbar)

    omz     = one-z
    omzbar  = one-zbar
    lomz    = log(omz)
    lomzbar = log(omzbar)

!! -- musq-independent  terms -- !!

    coeff_noL(1) = Cf**2*(omz**2+two*(one+z**2)*lomz)*(omzbar**2+two*(one+zbar**2)*lomzbar)/omz/omzbar
    coeff_noL(2) = -four*Cf**2*(omz**2+two*(one+z**2)*lomz)*lomzbar/omz/omzbar
    coeff_noL(3) = -four*Cf**2*lomz*(omzbar**2+two*(one+zbar**2)*lomzbar)/omz/omzbar
    coeff_noL(4) = 16.0_dp*Cf**2*lomz*lomzbar/omz/omzbar
    

!! -- log terms  -- !!

    coeff_L(1) = -two*Cf**2*( (one-z+z**2)*(one-zbar+zbar**2) - z*zbar + &
         (lomz+lomzbar)*(one + z**2)*(one + zbar**2) ) /omz/omzbar
    coeff_L(2) = two*Cf**2*(omz**2+two*(lomz+lomzbar)*(one+z**2))/omz/omzbar
    coeff_L(3) = two*Cf**2*(omzbar**2+two*(lomz+lomzbar)*(one+zbar**2))/omz/omzbar
    coeff_L(4) = -8.0_dp*Cf**2*(lomz+lomzbar)/omz/omzbar
    
!! -- log^2 terms  -- !!

    coeff_Lsq(1) = Cf**2*(one+z**2)*(one+zbar**2)/omz/omzbar
    coeff_Lsq(2) = -two*Cf**2*(one+z**2)/omz/omzbar
    coeff_Lsq(3) = -two*Cf**2*(one+zbar**2)/omz/omzbar
    coeff_Lsq(4) = four*Cf**2/omz/omzbar


    do j=1,4
       qqb_ns_sub_qqbz1zbar2(j) = coeff_noL(j) + coeff_L(j) * Lmu(j) + coeff_Lsq(j) * Lmu(j)**2
    enddo

    qqb_ns_sub_qqbz1zbar2 = qqb_ns_sub_qqbz1zbar2 * tag_ns_gg
    qqb_ns_sub_qqbz1zbar2 = qqb_ns_sub_qqbz1zbar2 * tag_z1zb2

    
  end function qqb_ns_sub_qqbz1zbar2


  

  function qqb_ns_sub_qqbz12(shat,mu2,z,yy)
! Integrated subtraction terms to multiply FLM[z1,2] and FLM[1,z2], cf. Eq. 10.17 and 10.18
! log(musq/s) evaluated for s={shat,shat/z}
    ! Output indices: {0=reg+pls,1=pls}
    ! color factor conventions: 1=Ca^2, 2=Ca*Nf, 3=Cf*Nf, 4=tot
    implicit none
    real(dp), intent(in) :: shat,mu2,z,yy(1:4)
    real(dp)             :: qqb_ns_sub_qqbz12(0:1)
    real(dp)             :: omz,opz,lz,lomz,lopz,li2omz,li2z,li2mz,li3omz,li3mz,li3z,Lmu,Lmu_z
    real(dp)             :: pl3noL_Cfsq,pl2noL_CaCf,pl2noL_CfNf,pl1noL_Cfsq,pl1noL_CaCf
    real(dp)             :: pl1noL_CfNf,pl0noL_Cfsq,pl0noL_CaCf,pl0noL_CfNf,regnoL_Cfsq,regnoL_CaCf,regnoL_CfNf,regnoL_Int
    real(dp)             :: pl2L_Cfsq,pl1L_Cfsq,pl1L_CaCf,pl1L_CfNf,pl0L_CaCf,pl0L_CfNf,regL_Cfsq,regL_CaCf,regL_CfNf,regL_Int
    real(dp)             :: pl1Lsq_Cfsq,pl0Lsq_Cfsq,pl0Lsq_CaCf,pl0Lsq_CfNf,regLsq_Cfsq,regLsq_CaCf,regLsq_CfNf
    real(dp)             :: plsnoL_Cfsq,plsnoL_CaCf,plsnoL_CfNf,plsL_Cfsq,plsL_CaCf,plsL_CfNf,plsLsq_Cfsq,plsLsq_CaCf,plsLsq_CfNf
    real(dp)             :: coeff_noL(0:1,1:4),coeff_L(0:1,1:4),coeff_Lsq(0:1,1:4)
    real(dp)             :: coeff_tc_pls_gg(1:3), coeff_tc_pls_qqb(1:4),tripcoll_gg(1:3,1:3),tripcoll_nf(1:4,1:3),coeff_tc_reg_gg(1:3),coeff_tc_reg_qqb(1:4)
    real(dp)             :: myyy(3),myz
    real(dp), parameter :: tcpls_coeff_Cfsq_nonmax=-9.1213_dp
    real(dp), parameter :: tcpls_coeff_CaCf_nonmax=15.920_dp
    real(dp), parameter :: tcpls_coeff_CfNf_nonmax=-9.349_dp*TR


    Lmu  = log(mu2/shat)
    Lmu_z  = log(mu2/shat) + log(z)
    
    omz    = one-z
    opz    = one+z
    lz   = log(z)
    lomz   = log(omz)
    lopz = log(opz)
    li2omz = dilog2(omz)
    li2z = dilog2(z)
    li2mz = dilog2(-z)
    li3omz = trilog(omz)
    li3mz  =  trilog(-z)
    li3z   =  trilog(+z)

    coeff_tc_pls_gg = zero
    coeff_tc_pls_qqb = zero
    
! -- musq-independent terms

    ! coefficients of + distr
    pl3noL_Cfsq = 8.0_dp

    pl2noL_CaCf = -22.0_dp/three
    pl2noL_CfNf = 8.0_dp/three*TR

    pl1noL_Cfsq = four*(one+log2)
    pl1noL_CaCf = 134.0_dp/9.0_dp -two/three*pisq
    pl1noL_CfNf = -40.0_dp/9.0_dp*TR

    pl0noL_Cfsq = four/three*log2*pisq + 16.0_dp*zeta3
    pl0noL_CaCf = -802.0_dp/27.0_dp - log2/three + 11.0_dp*log2**2 + 11.0_dp*pisq/18.0_dp + two*log2*pisq/three + 16.0_dp*zeta3
    pl0noL_CfNf = TR*(260.0_dp/27.0_dp - four*log2**2 - two/9.0_dp*pisq)

    plsnoL_Cfsq =  pl0noL_Cfsq/omz + pl1noL_Cfsq*lomz/omz                           + pl3noL_Cfsq*lomz**3/omz
    plsnoL_CaCf =  pl0noL_CaCf/omz + pl1noL_CaCf*lomz/omz + pl2noL_CaCf*lomz**2/omz
    plsnoL_CfNf =  pl0noL_CfNf/omz + pl1noL_CfNf*lomz/omz + pl2noL_CfNf*lomz**2/omz 

    ! regular function
    regnoL_Cfsq = (-42.0_dp - 60*li3omz - 72*log2 + 42*log2**2 - 24*lomz + 24*log2*lomz -&
         24*lomz**2 - 48*lomz**3 + 30*lz + 12*log2*lz - 21*log2**2*lz + 60*lomz*lz - 18*lomz**2*lz + 12*lz**2 - &
         42*log2*lz**2 - 24*lomz*lz**2 - 4*lz**3 + 20*pisq - 12*log2*pisq +&
         8*lomz*pisq - 4*lz*pisq + 72*z + 144*log2*z - 84*log2**2*z + &
         18*lomz*z - 96*log2*lomz*z + 48*lomz**2*z - &
         96*lz*z - 48*log2*lz*z + 24*lomz*lz*z - 24*lz**2*z - 48*pisq*z - 30*z**2 +&
         36*li3omz*z**2 - 72*log2*z**2 + 42*log2**2*z**2 + 6*lomz*z**2 +&
         72*log2*lomz*z**2 - 24*lomz**2*z**2 + &
         48*lomz**3*z**2 + 66*lz*z**2 + 36*log2*lz*z**2 + 21*log2**2*lz*z**2 + &
         60*lomz*lz*z**2 + 30*lomz**2*lz*z**2 + 12*lz**2*z**2 - 6*log2*lz**2*z**2 - 24*lomz*lz**2*z**2 + &
         4*lz**3*z**2 + 28*pisq*z**2 + 12*log2*pisq*z**2 - 8*lomz*pisq*z**2 + &
         4*lz*pisq*z**2 - 24*li2z*(-1 + z)*(-1 + 3*z + log2*(1 + z) - 2*lomz*(1 + z)) + &
         12*li2omz*((-1 + z)**2 + 3*lomz*(1 + z**2)) - 96*zeta3 + 96*z**2*zeta3)/(12.0_dp*omz)

    regnoL_CaCf = -(-594*log2**2*(-1 + z**2) + 108*li2omz*((-1 + z)**2 + 2*log2*(1 + z**2) + 3*lomz*(1 + z**2)) + &
         18*log2*(-44*lomz*(-1 + z**2) + lz*(4 + 22*z**2) -&
         (-1 + z)*(-3 + 7*z + 2*pisq*(1 + z))) + &
         (-1 + z)*(2048 + 1268*z + 396*lomz**2*(1 + z) - 3*pisq*(17 + 5*z) +&
         6*lomz*(-256 - 295*z + 12*pisq*(1 + z)) - 864*zeta3 - 864*z*zeta3))/(108.0_dp*omz)

    regnoL_CfNf = -TR*(54*log2**2*(-1 + z**2) - (-1 + z)*(36*lomz**2*(1 + z) - 3*pisq*(1 + z) -&
         6*lomz*(17 + 23*z) + 2*(77 + 53*z)) - 9*log2*((-1 + z)**2 - &
         8*lomz*(-1 + z**2) + 4*lz*(1 + z**2)))/(27.0_dp*omz)

    regnoL_Int = -(log2*(8 - 15*z + 7*z**2 + lz*(5 - 2*z**2) + 2*li2omz*(1 + z**2) + lz**2*(1 + z**2)))/omz


    !! triple collinear contribution
    if (withmax) then
       if (tcanalytic) then
          coeff_tc_reg_gg(1) = (  lomz*(6 - 6*z) + (lz**2*z)/two + (2*li3z*(5 + 3*z**2))/(-1 + z) +&
               (lz**3*(19 + 9*z**2))/(12 - 12*z) + log2**2*((7*(-1 + z))/two + (7*lz*(1 + z))/4.0_dp) + &
               li2z*(-4*z - 2*log2*(1 + z) - (2*lz*(1 + z**2))/(-1 + z)) + log2*(6 + 4*lomz*(-1 + z) - &
               6*z + pisq*(1 + z) + lz*(-1 + 3*z) + (lz**2*(7 + z**2))/(2 - 2*z)) + &
               (lz*(3*(-1 + z)*z - pisq*(5 + 3*z**2)))/(three*(-1 + z)) + (33 - 2*(33 + 2*pisq)*z + &
               z**2*(33 + 4*pisq - 36*zeta3) - 60*zeta3)/(6.0_dp*(-1 + z)) )*tag_tcreg
          coeff_tc_reg_gg(2) = (   (990 - 36*li3z - 12*log2 - 396*log2**2 - 608*lomz - 528*log2*lomz - &
               40*lz - 48*log2*lz + 36*lomz*lz + 144*log2*lomz*lz + 144*lomz**2*lz - 18*lz**2 + 18*lomz*lz**2 + 24*pisq - &
               48*log2*pisq - 276*z - 72*log2*z + 120*lomz*z + 114*lz*z - 72*lomz*lz*z + 24*pisq*z - 714*z**2 - &
               36*li3z*z**2 + 84*log2*z**2 + 396*log2**2*z**2 + 488*lomz*z**2 + &
               528*log2*lomz*z**2 + 98*lz*z**2 - 264*log2*lz*z**2 + 36*lomz*lz*z**2 + 144*log2*lomz*lz*z**2 + &
               144*lomz**2*lz*z**2 - 99*lz**2*z**2 + 18*lomz*lz**2*z**2 - 100*pisq*z**2 - &
               48*lomz*pisq*z**2 + 144*li3omz*(1 + z**2) + 12*li2z*(7 + 3*lz - 6*z + 25*z**2 + 3*lz*z**2 + &
               12*log2*(1 + z**2) + 12*lomz*(1 + z**2)) - 216*zeta3 + 288*z**2*zeta3)/(72.0_dp*(-1 + z)) )*tag_tcreg
          
          coeff_tc_reg_qqb(3) =  TR*( (-58 + 52*lomz + 14*lz + 9*lz**2 - 24*lomz*z - 12*lz*z + 58*z**2 - 28*lomz*z**2 - &
               10*lz*z**2 + 9*lz**2*z**2 + 8*pisq*z**2 - 36*log2**2*(-1 + z**2) - 24*li2z*(1 + z**2) + &
               6*log2*(5 - 6*z + z**2 - 8*lomz*(-1 + z**2) + 4*lz*(1 + z**2)))/(18.0_dp*(-1 + z)) )*tag_tcreg 
          coeff_tc_reg_qqb(4) = (  -(-6 - 156*li2z - 96*li3mz - 72*li3omz - 108*li3z + 192*log2 - 48*li2z*log2 + &
               192*lomz - 48*li2z*lomz + 132*lz + 36*li2z*lz + 120*log2*lz - 36*lomz*lz - 48*log2*lomz*lz - &
               48*lomz**2*lz + 75*lz**2 + 24*log2*lz**2 + 6*lomz*lz**2 + 14*lz**3 + 22*pisq + 8*log2*pisq + &
               8*lomz*pisq + 8*lz*pisq - 36*z + 72*li2z*z - 360*log2*z - 360*lomz*z - 114*lz*z + &
               72*lomz*lz*z - 36*lz**2*z - 12*pisq*z + 42*z**2 + 12*li2z*z**2 - 96*li3mz*z**2 - 72*li3omz*z**2 -&
               108*li3z*z**2 + 168*log2*z**2 - 48*li2z*log2*z**2 + 168*lomz*z**2 - &
               48*li2z*lomz*z**2 - 78*lz*z**2 + 36*li2z*lz*z**2 - 48*log2*lz*z**2 - 36*lomz*lz*z**2 - &
               48*log2*lomz*lz*z**2 - 48*lomz**2*lz*z**2 - 12*lz**2*z**2 + 24*log2*lz**2*z**2 + &
               6*lomz*lz**2*z**2 + 14*lz**3*z**2 + 2*pisq*z**2 + 8*log2*pisq*z**2 + 8*lomz*pisq*z**2 + &
               8*lz*pisq*z**2 + 48*li2mz*(-1 + lz + z**2 + lz*z**2) + 36*zeta3 + 36*z**2*zeta3 + &
               48*lz*(-1 + z**2)*Log(1 + z))/(24.0_dp*(-1 + z)) )*tag_tcreg
!          print *, "z",z
!          print *, "reg Cf^2",coeff_tc_reg_gg(1)*Cf**2
!          print *, "reg Cf^2",coeff_tc_reg_gg(2)*Cf*Ca
!          print *, "reg Nf*Cf",coeff_tc_reg_qqb(3)*Cf
!          print *, "reg Cf*(Ca-2Cf)",coeff_tc_reg_qqb(4)*Cf*(Ca-two*Cf)
!          pause
       else
          call get_tripcoll_subtr_q_gg_qstar(z,yy,tripcoll_gg)
          coeff_tc_reg_gg(1:3) = tripcoll_gg(1:3,3) * tag_tcreg
          call get_tripcoll_subtr_q_qqb_qstar(z,yy,tripcoll_nf)
          coeff_tc_reg_qqb(1:4) = tripcoll_nf(1:4,3) * tag_tcreg
       endif
       coeff_tc_pls_gg(1) = -four/three*pisq*log2*tag_tcpls
       coeff_tc_pls_gg(2) = (22.0_dp + 11.0_dp/9.0_dp*pisq -11.0_dp*log2**2-7.0_dp*zeta3 +&
            one/three*log2*(-one-two*pisq))*tag_tcpls
       coeff_tc_pls_qqb(3) = TR*(-four/9.0_dp*(16.0_dp+pisq)+four/three*log2+four*log2**2)*tag_tcpls
!          print *, "pls Cf^2",coeff_tc_pls_gg(1)*Cf**2
!          print *, "pls CaCf",coeff_tc_pls_gg(2)*Cf*Ca
!          print *, "pls Nf*Cf",coeff_tc_pls_qqb(3)*Cf
!          pause
    else
       coeff_tc_pls_gg(1)   = tcpls_coeff_Cfsq_nonmax * tag_tcpls
       coeff_tc_pls_gg(2)   = tcpls_coeff_CaCf_nonmax * tag_tcpls
       coeff_tc_pls_qqb(3) = tcpls_coeff_CfNf_nonmax * tag_tcpls

       call get_tripcoll_subtr_q_gg_qstar(z,yy,tripcoll_gg)
       coeff_tc_reg_gg(1:3) = tripcoll_gg(1:3,3) * tag_tcreg
       call get_tripcoll_subtr_q_qqb_qstar(z,yy,tripcoll_nf)
       coeff_tc_reg_qqb(1:4) = tripcoll_nf(1:4,3) * tag_tcreg
    endif
       

    
    plsnoL_Cfsq = plsnoL_Cfsq*tag_nontc + coeff_tc_pls_gg(1)/omz*tag_tcpls
    plsnoL_CaCf = plsnoL_CaCf*tag_nontc + coeff_tc_pls_gg(2)/omz*tag_tcpls
    plsnoL_CfNf = plsnoL_CfNf*tag_nontc + coeff_tc_pls_qqb(3)/omz*tag_tcpls

    regnoL_Cfsq = regnoL_Cfsq*tag_nontc + coeff_tc_reg_gg(1)
    regnoL_CaCf = regnoL_CaCf*tag_nontc + coeff_tc_reg_gg(2)
    regnoL_CfNf = regnoL_CfNf*tag_nontc + coeff_tc_reg_qqb(3)
    regnoL_Int =  regnoL_Int*tag_nontc +  coeff_tc_reg_qqb(4)



     ! add the different color structures
     coeff_noL(+1,tot)  = Cf**2*plsnoL_Cfsq*tag_ns_gg + Ca*Cf*plsnoL_CaCf*tag_ns_gg + Cf*Nf*plsnoL_CfNf*tag_ns_nf     
     coeff_noL(0,tot)   = Cf**2*regnoL_Cfsq*tag_ns_gg + Ca*Cf*regnoL_CaCf*tag_ns_gg + Cf*Nf*regnoL_CfNf*tag_ns_nf + tag_ns_id*Cf*(Ca-two*Cf)*regnoL_Int

! -- log terms --

    ! plus distribution
    pl2L_Cfsq = -12.0_dp
    
    pl1L_Cfsq = -12.0_dp
    pl1L_CaCf = 22.0_dp/three
    pl1L_CfNf = -8.0_dp/three*TR

    pl0L_CaCf = (-67.0_dp+three*pisq)/9.0_dp
    pl0L_CfNf = 20.0_dp/9.0_dp*TR

    plsL_Cfsq =                  pl1L_Cfsq*lomz/omz + pl2L_Cfsq*lomz**2/omz
    plsL_CaCf =  pl0L_CaCf/omz + pl1L_CaCf*lomz/omz
    plsL_CfNf =  pl0L_CfNf/omz + pl1L_CfNf*lomz/omz 
    
    ! regular function
    regL_Cfsq = -(36 - 48*lomz - 36*lomz**2 + 12*lz - 12*lomz*lz + 3*lz**2 + 2*pisq - 66*z + &
         24*lomz*z + 12*lz*z + 30*z**2 + 24*lomz*z**2 + 36*lomz**2*z**2 + 12*lz*z**2 - 12*lomz*lz*z**2 - &
         3*lz**2*z**2 - 2*pisq*z**2 + 12*li2z*(-1 + z**2) + 12*li2omz*(1 + z**2))/(6.0_dp*omz)
    regL_CaCf = -(-52 + 6*lz + 3*pisq - 39*z + 91*z**2 + 33*lz*z**2 - 3*pisq*z**2 -&
         66*lomz*(-1 + z**2) - 18*li2omz*(1 + z**2))/(18.0_dp*omz)
    regL_CfNf = -two*TR*(two + 6.0_dp*z - 8.0_dp*z**2 + 6.0_dp*lomz*(-1 + z**2) - 3.0_dp*lz*(1 + z**2))/(9.0_dp*omz)
    regL_Int =  -(8 - 15*z + 7*z**2 + lz*(5 - 2*z**2) + 2*li2omz*(1 + z**2) + lz**2*(1 + z**2))/(two*omz)

    plsL_Cfsq = plsL_Cfsq*tag_nontc
    plsL_CaCf = plsL_CaCf*tag_nontc
    plsL_CfNf = plsL_CfNf*tag_nontc
    regL_Cfsq = regL_Cfsq*tag_nontc
    regL_CaCf = regL_CaCf*tag_nontc
    regL_CfNf = regL_CfNf*tag_nontc
    regL_Int = regL_Int*tag_nontc
    
    ! add the different color structures
    coeff_L(+1,tot) = tag_ns_gg*(Cf**2*plsL_Cfsq + Ca*Cf*plsL_CaCf) + tag_ns_nf*Cf*Nf*plsL_CfNf
    coeff_L(0,tot)  = tag_ns_gg*(Cf**2*regL_Cfsq + Ca*Cf*regL_CaCf) + tag_ns_nf*Cf*Nf*regL_CfNf +  tag_ns_id*Cf*(Ca-two*Cf)*regL_Int
    
! -- log^2 terms --

    ! plus distribution
    pl1Lsq_Cfsq = four

    pl0Lsq_Cfsq = 6.0_dp
    pl0Lsq_CaCf = -11.0_dp/6.0_dp
    pl0Lsq_CfNf = two*Tr/three

    plsLsq_Cfsq =  pl0Lsq_Cfsq/omz + pl1Lsq_Cfsq*lomz/omz 
    plsLsq_CaCf =  pl0Lsq_CaCf/omz
    plsLsq_CfNf =  pl0Lsq_CfNf/omz

    ! regular function
    regLsq_Cfsq = -(lz + 3*lz*z**2 - 4*lomz*(-1 + z**2) - 4*(-2 + z + z**2))/(two*omz)
    regLsq_CaCf = 11.0_dp/12.0_dp*(one+Z)
    regLsq_CfNf = -TR*(one + z)/three

    plsLsq_Cfsq = plsLsq_Cfsq*tag_nontc
    plsLsq_CaCf = plsLsq_CaCf*tag_nontc
    plsLsq_CfNf = plsLsq_CfNf*tag_nontc
    regLsq_Cfsq = regLsq_Cfsq*tag_nontc
    regLsq_CaCf = regLsq_CaCf*tag_nontc
    regLsq_CfNf = regLsq_CfNf*tag_nontc
    



! add the different color structures
    coeff_Lsq(+1,tot) = tag_ns_gg*(Cf**2*plsLsq_Cfsq + Ca*Cf*plsLsq_CaCf) + tag_ns_nf*Cf*Nf*plsLsq_CfNf
    coeff_Lsq(0,tot)  = tag_ns_gg*(Cf**2*regLsq_Cfsq + Ca*Cf*regLsq_CaCf) + tag_ns_nf*Cf*Nf*regLsq_CfNf
    ! combine  -- reg+ pls
        qqb_ns_sub_qqbz12(0) = sum(coeff_noL(0:1,tot),dim=1) + sum(coeff_L(0:1,tot),dim=1)*Lmu_z + sum(coeff_Lsq(0:1,tot),dim=1)*Lmu_z**2
! pls  -- evaluated with z=1 in log, i.e. log(musq/shat)
        qqb_ns_sub_qqbz12(1) = coeff_noL(1,tot) + coeff_L(1,tot)*Lmu + coeff_Lsq(1,tot)*Lmu**2


        ! debug
        qqb_ns_sub_qqbz12 = qqb_ns_sub_qqbz12 * tag_z12

      end function qqb_ns_sub_qqbz12



  function qqb_ns_sub_qqb12(shat,mu2)
    ! Integrated subtraction terms to multiply FLM[1,2], cf. Eq. 10.19
    ! color indices: 1=Cf^2, 2=Ca*Cf, 3=Cf*Nf,4=total
    implicit none
    real(dp), intent(in)  :: shat,mu2
    real(dp)              :: qqb_ns_sub_qqb12
    real(dp)              :: Lmu
    real(dp)              :: coeff_noL(1:4),coeff_L(1:4),coeff_Lsq(1:4)
    real(dp)              :: tc_coeff_Cfsq,tc_coeff_CaCf,tc_coeff_CfNf
    real(dp), parameter   :: tc_coeff_Cfsq_nonmax = 1.1184471817543686_dp
    real(dp), parameter   :: tc_coeff_CaCf_nonmax = 8.50406915208553_dp
    real(dp), parameter   :: tc_coeff_CfNf_nonmax = 1.0848_dp*TR

    Lmu  = log(mu2/shat)
    coeff_Lsq = zero
    coeff_L = zero
    coeff_noL = zero
    
! -- musq-independent  terms -- !!
    coeff_noL(1) = 8.0_dp/45.0_dp*pisq**2 + 7.0_dp/three*log2**2*pisq - 16.0_dp*log2*zeta3
    coeff_noL(2) = 3565.0_dp/81.0_dp + (352*log2**3)/9.0_dp - (59*Pi**4)/120.0_dp + (667*pisq)/108.0_dp +&
         (2*log2**2*(70 + 9*pisq))/9.0_dp - (374*zeta3)/9.0_dp - (log2*(445 + 165*pisq + 702*zeta3))/54.0_dp
    coeff_noL(3) = (TR*(-3938.0_dp - 2952.0_dp*log2**2 - 4608.0_dp*log2**3 - 681.0_dp*pisq +&
         12.0_dp*log2*(121.0_dp + 30.0_dp*pisq)  + 4896.0_dp*zeta3))/324.0_dp

    if (withmax) then 
       tc_coeff_Cfsq = (-7*log2**2*pisq)/6.0_dp + 8*log2*zeta3
       tc_coeff_CaCf = (-91080 - 8640*li4half - 84480*log2**3 - 360*log2**4 + 633*Pi**4 - &
            3010*pisq - 120*log2**2*(316 + 33*pisq) + 65340*zeta3 + 40*log2*(913 + 165*pisq + 513*zeta3))/4320.0_dp
       tc_coeff_CfNf = TR*(1018 + 1416*log2**2 + 1536*log2**3 + 97*pisq - 20*log2*(53 + 6*pisq) - 1188*zeta3)/216.0_dp
!          print *, "del Cf^2",tc_coeff_Cfsq*Cf**2
!          print *, "del CaCf",tc_coeff_CaCf*Ca*Cf
!          print *, "del Nf*Cf",tc_coeff_CfNf*Cf
!          pause
    else
       tc_coeff_Cfsq = tc_coeff_Cfsq_nonmax 
       tc_coeff_CaCf = tc_coeff_CaCf_nonmax 
       tc_coeff_CfNf = tc_coeff_CfNf_nonmax
    endif
    coeff_noL(1) = coeff_noL(1)*tag_nontc + tc_coeff_Cfsq * tag_tcdel * two
    coeff_noL(2) = coeff_noL(2)*tag_nontc + tc_coeff_CaCf * tag_tcdel * two
    coeff_noL(3) = coeff_noL(3)*tag_nontc + tc_coeff_CfNf * tag_tcdel * two


    coeff_noL(4) = tag_ns_gg*(Cf**2*coeff_noL(1) + Ca*Cf*coeff_noL(2)) + tag_ns_nf*Cf*Nf*coeff_noL(3)



! --  log terms -- 
    coeff_L(1) = -two*(pisq + 8.0_dp*zeta3)
    coeff_L(2) = -199.0_dp/54.0_dp +23.0_dp/24.0_dp*pisq - 7.0_dp*zeta3
    coeff_L(3) = TR*(68.0_dp - 9.0_dp*pisq)/54.0_dp

    coeff_L(4) = tag_ns_gg*(Cf**2*coeff_L(1) + Ca*Cf*coeff_L(2)) + tag_ns_nf*Cf*Nf*coeff_L(3) 

    coeff_L = coeff_L*tag_nontc
  
    ! --  log^2 terms  -- !!
    coeff_Lsq = zero
    coeff_Lsq(1) = 9.0_dp/two - two/three*pisq
    coeff_Lsq(2) = -11.0_dp/four
    coeff_Lsq(3) = TR

    coeff_Lsq(4) = tag_ns_gg*(Cf**2*coeff_Lsq(1) + Ca*Cf*coeff_Lsq(2)) + tag_ns_nf*Cf*Nf*coeff_Lsq(3) 

    coeff_Lsq = coeff_Lsq*tag_nontc

    qqb_ns_sub_qqb12 = coeff_noL(4) + Lmu * coeff_L(4) + Lmu**2 * coeff_Lsq(4)
 

    ! debug
    qqb_ns_sub_qqb12 = qqb_ns_sub_qqb12 * tag_12

    

  end function qqb_ns_sub_qqb12


  subroutine get_tripcoll_subtr_q_gg_qstar(z,yy,tripcoll)
! input=z, and integration variables yy(1:3)
! output indices: 1st is color structure {ab,nab,tot}, 2nd is {del,pls,reg}
    implicit none
    real(dp), intent(in)  :: z,yy(1:4)
    real(dp), intent(out) :: tripcoll(1:3,1:3)
    real(dp), parameter   :: zmin = zero
    real(dp)              :: x1,x2,x3,x3a,xx2,xx(1:3),xx3
    real(dp)              :: seca(1:3,1:3,0:1),secb(1:3,1:3,0:1),secc(1:3,1:3,0:1),secd(1:3,1:3,0:1)
    real(dp)              :: seca0(1:3,1:3,0:1),secb0(1:3,1:3,0:1),secc0(1:3,1:3,0:1),secd0(1:3,1:3,0:1)
    real(dp)              :: wghtac,wghtbd
    real(dp)              :: res_tc(1:3,1:3,0:1)

    tripcoll = zero

    !--
    seca = zero; seca0 = zero
    secb = zero; secb0 = zero
    secc = zero; secc0 = zero
    secd = zero; secd0 = zero

    x1 =  yy(1)   ! r
    x2 =  yy(2)   ! x4
    x3a = yy(3)  ! lambda
    xx3 = yy(4)  ! x3

    ! remap lambda
    x3 = sin( half * pi * x3a )**2

    ! !-- debug
    ! print*, 'debug'
    ! x1 = 0.123_dp
    ! x2 = 0.456_dp
    ! x3 = 0.678_dp
    ! xx3 = 0.789_dp
    ! !-- end debug

    if (x1*x2*(one-z)*xx3 .lt. bufftc) return

    !--

    xx = (/x1,x2,x3/)
    seca0 = PTC_A_q_gg_qst(z,zmin,xx,zero)
    secc0 = PTC_C_q_gg_qst(z,zmin,xx,zero)

    if (withmax) then 
       seca = PTC_A_q_gg_qst(z,zmin,xx,xx3)
       secc = PTC_C_q_gg_qst(z,zmin,xx,xx3)
    endif

    wghtac = -one/x1/x2

    !--
    
    xx2 = x2**2
    xx = (/x1,xx2,x3/)
    secb0 = PTC_B_q_gg_qst(z,zmin,xx,zero)
    secd0 = PTC_D_q_gg_qst(z,zmin,xx,zero)

    if (withmax) then
       secb = PTC_B_q_gg_qst(z,zmin,xx,xx3)
       secd = PTC_D_q_gg_qst(z,zmin,xx,xx3)
    endif

    wghtbd = -one/x1/x2*two

    if (xx2 .lt. bufftc) then

       res_tc = (seca0 + secc0)*wghtac

       if (withmax) then 
          res_tc(:,:,1) = res_tc(:,:,1) - two*( &
               seca(:,:,0)-seca0(:,:,0) + &
               secc(:,:,0)-secc0(:,:,0))/xx3*wghtac
       endif

    else

       res_tc = (seca0 + secc0)*wghtac + (secb0 + secd0)*wghtbd 

       if (withmax) then
          res_tc(:,:,1) = res_tc(:,:,1) - two*( &
               seca(:,:,0)-seca0(:,:,0) + &
               secc(:,:,0)-secc0(:,:,0))/xx3*wghtac

          res_tc(:,:,1) = res_tc(:,:,1) - two*( &
               secb(:,:,0)-secb0(:,:,0) + &
               secd(:,:,0)-secd0(:,:,0))/xx3*wghtbd
       endif

    endif

    !-- remove color factors for abelian and non-abelian contributions                                                                                                                                            
    tripcoll(1,:) = res_tc(1,:,1)/Cf**2
    tripcoll(2,:) = res_tc(2,:,1)/Cf/Ca
    tripcoll(3,:) = res_tc(1,:,1) + res_tc(2,:,1)


  end subroutine get_tripcoll_subtr_q_gg_qstar

  subroutine get_tripcoll_subtr_q_qqb_qstar(z,yy,tripcoll)
! input=z, and integration variables yy(1:3)
! output indices: 1st is flavor structure: 1=singlet diff, 2=singlet id, 3=non-singlet diff,4=non-singlet id 
!                 2nd is {del,pls,reg}
    implicit none
    real(dp), intent(in)  :: z,yy(1:4)
    real(dp), intent(out) :: tripcoll(1:4,1:3)
    real(dp), parameter   :: zmin = zero
    real(dp)              :: x1,x2,x3,x3a,xx2,xx(1:3),xx3
    real(dp)              :: seca(1:4,1:3,0:2),secb(1:4,1:3,0:2),secc(1:4,1:3,0:2),secd(1:4,1:3,0:2)
    real(dp)              :: seca0(1:4,1:3,0:2),secb0(1:4,1:3,0:2),secc0(1:4,1:3,0:2),secd0(1:4,1:3,0:2)
    real(dp)              :: wghtac,wghtbd
    real(dp)              :: res_tc(1:4,1:3,0:2)

    tripcoll = zero

    !--
    seca = zero; seca0 = zero
    secb = zero; secb0 = zero
    secc = zero; secc0 = zero
    secd = zero; secd0 = zero

    x1 =  yy(1)   ! r
    x2 =  yy(2)   ! x4
    x3a = yy(3)  ! lambda
    xx3 = yy(4)  ! x3

    ! remap lambda
    x3 = sin( half * pi * x3a )**2

    ! !-- debug
    ! print *, 'debug'
    ! x1 = 0.123_dp
    ! x2 = 0.456_dp
    ! x3 = 0.678_dp
    ! xx3 = 0.789_dp
    ! !-- end debug

    if (x1*x2*(one-z)*xx3 .lt. bufftc) return

    !--

    xx = (/x1,x2,x3/)
    seca0 = PTC_A_q_qqb_qst(z,zmin,xx,zero)
    secc0 = PTC_C_q_qqb_qst(z,zmin,xx,zero)

    if (withmax) then 
       seca = PTC_A_q_qqb_qst(z,zmin,xx,xx3) 
       secc = PTC_C_q_qqb_qst(z,zmin,xx,xx3) 
    endif

    wghtac = -one/x1/x2

    !--
    
    xx2 = x2**2
    xx = (/x1,xx2,x3/)
    secb0 = PTC_B_q_qqb_qst(z,zmin,xx,zero)
    secd0 = PTC_D_q_qqb_qst(z,zmin,xx,zero)

    if (withmax) then 
       secb = PTC_B_q_qqb_qst(z,zmin,xx,xx3) 
       secd = PTC_D_q_qqb_qst(z,zmin,xx,xx3) 
    endif

    wghtbd = -one/x1/x2*two

    if (xx2 .lt. bufftc) then

       res_tc = (seca0 + secc0)*wghtac

       if (withmax) then
          res_tc(:,:,1) = res_tc(:,:,1) - two*( &
               seca(:,:,0)-seca0(:,:,0) + &
               secc(:,:,0)-secc0(:,:,0))/xx3*wghtac
       endif

    else

       res_tc = (seca0 + secc0)*wghtac + (secb0 + secd0)*wghtbd 

       if (withmax) then
          
          res_tc(:,:,1) = res_tc(:,:,1) - two*( &
               seca(:,:,0)-seca0(:,:,0) + &
               secc(:,:,0)-secc0(:,:,0))/xx3*wghtac
          
          res_tc(:,:,1) = res_tc(:,:,1) - two*( &
               secb(:,:,0)-secb0(:,:,0) + &
               secd(:,:,0)-secd0(:,:,0))/xx3*wghtbd
          
       endif
          
    endif

    tripcoll(:,:) = res_tc(:,:,1)

    ! remove color factor Cf(Cf-2Ca) for the interference case (myfl=4); otherwisre remove color factor Cf * Tr                                                                          
    tripcoll(1,:) = tripcoll(1,:)/(Cf)
    tripcoll(3,:) = tripcoll(3,:)/(Cf)
    tripcoll(2,:) = tripcoll(2,:)/Cf/(Ca-two*Cf)
    tripcoll(4,:) = tripcoll(4,:)/Cf/(Ca-two*Cf)



  end subroutine get_tripcoll_subtr_q_qqb_qstar

end module mod_subtrfn_nnlo_dy_ns

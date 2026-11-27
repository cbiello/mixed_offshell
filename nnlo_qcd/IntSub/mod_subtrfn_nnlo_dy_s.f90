module mod_subtrfn_nnlo_dy_s
! integrated subtraction terms for gg->ff channel, where ff=gg or qqb
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_tripleclimits_q_qqb_qstar
  implicit none


  private 

  public :: qqb_s_sub_gqz124,qqb_s_sub_qqbz12,qqb_v_sub_qqbz12

  contains 

  function qqb_s_sub_gqz124(shat,mu2,z,eta41,eta42,eta41log)
    ! Integrated subtraction terms to multiply FLMqg[z1,2,4] 
    ! Regular only
    real(dp), intent(in)   :: shat,mu2,z,eta41,eta42,eta41log
    real(dp)               :: qqb_s_sub_gqz124
    real(dp)               :: Lmu,omz,lomz,w51tc

    Lmu = log(mu2/shat)
    omz = one-z
    lomz = log(omz)
    w51tc = eta42*(one+eta41)                          ! w4151tilde51

! scale dependent piece
    qqb_s_sub_gqz124 = -Cf*(two-two*z+z**2)/z * Lmu

! scale independent pieces
    qqb_s_sub_gqz124 = qqb_s_sub_gqz124 + &
         Cf*( w51tc*(two-two*z+z**2)*log(eta41log/2) +&
         z**2 + two*(two-two*z+z**2)*lomz )/z
    

  end function qqb_s_sub_gqz124

  function qqb_s_sub_qqbz12(shat,mu2,z,yy)
    ! Regular only
    implicit none
    real(dp), intent(in) :: shat,mu2,z,yy(1:4)
    real(dp)             :: qqb_s_sub_qqbz12
    real(dp)             :: omz,lz,lomz,li2omz,li2mz,li2z,li3omz,li3mz,li3z,li3omzsq,li3zoopz
    real(dp)             :: Lmu,coeff_noL,coeff_L,coeff_Lsq,tripcoll(1:4,1:3)


    Lmu  = log(mu2/shat) + log(z)
    
    omz    = one-z
    lz   = log(z)
    lomz   = log(omz)

    li2omz = dilog2(omz)
    li2z = dilog2(z)
    li2mz = dilog2(-z)
    li3omz = trilog(omz)
    li3mz  =  trilog(-z)
    li3z   =  trilog(+z)
    li3omzsq = trilog(one-z**2)
    li3zoopz = trilog(z/(one+z))

    coeff_noL  =            (21*log2**2*(-4 - 3*(1 + 2*lz)*z + (3 - 6*lz)*z**2 + 4*z**3) - &
         2*log2*(110 - 177*z - 72*li2z*z + 12*pisq*z + 177*z**2 - 72*li2z*z**2 + 12*pisq*z**2 -&
         110*z**3 + 54*lz**2*z*(1 + z) + 12*lz*(10 + 9*z + 9*z**2 + 2*z**3) - &
         12*lomz*(-4 - 3*z + 3*z**2 + 4*z**3)) - 4*(9 + 4*pisq - 45*z + 72*li2z*z + 72*li3omz*z -&
         9*pisq*z + 39*z**2 + 108*li2z*z**2 + 72*li3omz*z**2 - 21*pisq*z**2 - 3*z**3 + &
         48*li2z*z**3 - 12*pisq*z**3 + 6*lz**3*z*(1 + z) - &
         6*lz*(-1 - (-3 + pisq)*z - (-12 + pisq)*z**2 + z**3) - 6*lz**2*(-2 + 3*z**2 + 2*z**3) + &
         6*lomz**2*(-4 + (-3 + 6*lz)*z + (3 + 6*lz)*z**2 + 4*z**3) + &
         6*lomz*z*(-5 - z + 6*z**2 + 12*li2z*(1 + z) - 2*pisq*(1 + z))))/(36.0_dp*z)
    
    coeff_L = ((-26 - 24*lz + 30*z - 45*lz*z - 9*lz**2*z - 6*pisq*z - 48*z**2 - 63*lz*z**2 -&
         9*lz**2*z**2 - 6*pisq*z**2 + 44*z**3 - 24*lz*z**3 + 36*li2z*z*(1 + z) + &
         6*lomz*(-4 - 3*z + 3*z**2 + 4*z**3)))/(9.0_dp*z)
    
    coeff_Lsq = -((-4 - 3*z + 3*z**2 + 4*z**3 - 6*lz*z*(1 + z)))/(6.0_dp*z)

    
    qqb_s_sub_qqbz12 = Cf*TR*(coeff_noL + coeff_L*Lmu + coeff_Lsq*Lmu**2)*tag_nontc
    
    if (withmax .and. tcanalytic) then
       tripcoll(1,3) =   (788 - 356*z**3 - 432*li3z*z*(1 + z) + 90*lz**3*z*(1 + z) + 24*lomz*(26 - 45*z + 45*z**2 - 26*z**3) -&
            144*li2mz*(2 + 3*z + 3*z**2 + 2*z**3) - 12*pisq*(-8 + 3*z - 3*z**2 + 4*z**3) + &
            9*lz**2*(48 + 33*z + 57*z**2 + 8*z**3) + 63*log2**2*(4 + 3*z - 3*z**2 - 4*z**3 + 6*lz*z*(1 + z)) + &
            6*log2*(122 - 213*z + 213*z**2 - 122*z**3 + 54*lz**2*z*(1 + z) + 12*pisq*z*(1 + z) + &
            12*lomz*(4 + 3*z - 3*z**2 - 4*z**3) + 12*lz*(10 + 9*z + 9*z**2 + 2*z**3)) + 6*z**2*(-41 + 72*zeta3) + &
            6*z*(-31 + 72*zeta3) - 144*li2z*(5 + 3*(1 + log2)*z**2 - z**3 + z*three*log2) +& 
            12*lz*(70 + (57 + 6*pisq)*z + 6*(9 + pisq)*z**2 + 26*z**3 + 6*lomz*(-2 + 3*z - 3*z**2 + 2*z**3) - &
            12*(1 + z)*(2 + z + 2*z**2)*Log(1 + z)))/(108.0_dp*z)        
!       print *, "z",z
!       print *, "S",tripcoll(1,3)
!       pause
    else
       call get_tripcoll_subtr_q_qqb_qstar(z,yy,tripcoll)    
    endif
    
    qqb_s_sub_qqbz12 = qqb_s_sub_qqbz12 + Cf*TR*tripcoll(1,3) *tag_tcreg



  end function qqb_s_sub_qqbz12
    


  function qqb_v_sub_qqbz12(shat,mu2,z,yy)
    ! Regular only
    implicit none
    real(dp), intent(in) :: shat,mu2,z,yy(1:4)
    real(dp)             :: qqb_v_sub_qqbz12
    real(dp)             :: omz,lz,lomz,li2z,li2mz,li3omz,li3mz,li3z,li3omzsq,li3zoopz
    real(dp)             :: Lmu,coeff_noL,coeff_L,coeff_Lsq,tripcoll(1:4,1:3)


    Lmu  = log(mu2/shat) + log(z)
    
    omz    = one-z
    lz   = log(z)
    lomz   = log(omz)
    li2mz = dilog2(-z)
    li2z = dilog2(z)
    li3mz  =  trilog(-z)
    li3z   =  trilog(+z)
    li3omzsq = trilog(one-z**2)
    li3zoopz = trilog(z/(one+z))
    li3omz = trilog(one-z)



    coeff_L = -(-12 - 6*lz - 3*lz**2 + pisq - 12*lz*z + 12*z**2 - 6*lz*z**2 - &
         3*lz**2*z**2 + pisq*z**2 + 12*li2mz*(1 + z**2) + &
         12*lz*(1 + z**2)*Log(1 + z))/(6.0_dp*(1 + z))
    coeff_noL   =  -(log2*(-12 + pisq + 12*z**2 + pisq*z**2 + 12*li2mz*(1 + z**2) - &
         3*lz**2*(1 + z**2) + 6*lz*(-(1 + z)**2 + 2*(1 + z**2)*Log(1 + z))))/(3.0_dp*(1 + z))


    qqb_v_sub_qqbz12 = Cf*(Ca-two*Cf)*(coeff_noL + coeff_L*Lmu) * tag_nontc    

    if (withmax .and. tcanalytic) then


       tripcoll(2,3) =   (-12*lomz*lz*(-1 + z)*(1 + z) + 36*li2mz*(1 + z)**2 - 24*log2*lz*(1 + z)**2 - &
            12*lz**2*(1 + z)*(2 + z) + 12*li2z*(1 + z)*(3 + z) + 108*li3mz*(1 + z**2) - 48*li3omz*(1 + z**2) + &
            24*li3omzsq*(1 + z**2) + 48*li3z*(1 + z**2) + 72*li3zoopz*(1 + z**2) + 48*li2mz*log2*(1 + z**2) + &
            48*li2mz*lomz*(1 + z**2) - 12*li2mz*lz*(1 + z**2) - 12*li2z*lz*(1 + z**2) - 12*log2*lz**2*(1 + z**2) - &
            7*lz**3*(1 + z**2) + pisq*(-3 - 2*z + z**2) + 4*log2*(-12 + pisq + 12*z**2 + pisq*z**2) + 4*lomz*(-12 + pisq + 12*z**2 + pisq*z**2) - &
            lz*(57 + 90*z + 33*z**2 + 4*pisq*(1 + z**2)) - 15*(3 + 2*zeta3 + z**2*(-3 + 2*zeta3)) + &
            36*lz*(1 + z)**2*Log(1 + z) + 48*log2*lz*(1 + z**2)*Log(1 + z) + &
            48*lomz*lz*(1 + z**2)*Log(1 + z) + 6*lz**2*(1 + z**2)*Log(1 + z) + 6*pisq*(1 + z**2)*Log(1 + z) + &
            36*lz*(1 + z**2)*Log(1 + z)**2 - 12*(1 + z**2)*Log(1 + z)**3)/(12.0_dp*(1 + z))
!       print *, "z",z
!       print *, "V",tripcoll(2,3)
!       pause
    else
       call get_tripcoll_subtr_q_qqb_qstar(z,yy,tripcoll)    
    endif
    qqb_v_sub_qqbz12 = qqb_v_sub_qqbz12  + Cf*(Ca-two*Cf)*tripcoll(2,3) *tag_tcreg


  end function qqb_v_sub_qqbz12
      





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
    tripcoll(1,:) = tripcoll(1,:)/(Cf*TR)
    tripcoll(3,:) = tripcoll(3,:)/(Cf*TR)
    tripcoll(2,:) = tripcoll(2,:)/Cf/(Ca-two*Cf)
    tripcoll(4,:) = tripcoll(4,:)/Cf/(Ca-two*Cf)


  end subroutine get_tripcoll_subtr_q_qqb_qstar  


end module mod_subtrfn_nnlo_dy_s

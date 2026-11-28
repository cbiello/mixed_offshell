module mod_subtrfn_nnlo_dy_gq
! integrated subtraction terms for qg->qg channel
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions; use mod_procfunctions
  use mod_tripleclimits_g_qg_qstar
  
  implicit none  

  private 

  public :: gq_qg_sub_gqz124, gq_qg_sub_gq1z24, gq_qg_sub_gq124
  public :: gq_qg_sub_gqz124_pls, gq_qg_sub_gq1z24_pls
  public :: gq_qg_sub_qqbz124
  public :: gq_qg_sub_qqbz1zbar2,gq_qg_sub_qqbz12

contains

   
  function gq_qg_sub_gqz124(shat,mu2,z,eta41,eta42,eta41log)
    ! Integrated subtraction terms to multiply FLMgq[z1,2,4]
    ! This gives reg+pls
    real(dp), intent(in) :: shat,mu2,z,eta41,eta42,eta41log
    real(dp)             :: gq_qg_sub_gqz124      
    real(dp)             :: Lmu,omz,lomz,w51tc

    w51tc = eta42*(one+eta41)                          ! w4151tilde51

    Lmu = log(mu2/shat)
    omz = one - z
    lomz = log(omz)

    ! scale dependent pieces
    gq_qg_sub_gqz124 = two*Ca*(-one/z+two-z+z**2)  * Lmu                    ! reg
    gq_qg_sub_gqz124 = gq_qg_sub_gqz124 -two*Ca/omz * Lmu   ! pls

    ! scale independent pieces
    gq_qg_sub_gqz124 = gq_qg_sub_gqz124 - two*Ca*(-one/z+two-z+z**2)*(w51tc*log(eta41log/two)+2*lomz)          ! reg
    gq_qg_sub_gqz124 = gq_qg_sub_gqz124 + four*Ca*lomz/omz + two*Ca*w51tc*log(eta41log/two)/omz                 ! pls

    
  end function gq_qg_sub_gqz124


  function gq_qg_sub_gqz124_pls(shat,mu2,z,eta41,eta42,eta41log)
    ! Integrated subtraction terms to multiply FLMqg[z1,2,4]
    ! This gives pls(z=1)
    real(dp), intent(in) :: shat,mu2,z,eta41,eta42,eta41log
    real(dp)             :: gq_qg_sub_gqz124_pls     
    real(dp)             :: Lmu,omz,lomz,w51tc

    w51tc = eta42*(one+eta41)                          ! w4151tilde51

    Lmu = log(mu2/shat)
    omz = one - z
    lomz = log(omz)

    ! scale dependent pieces
    gq_qg_sub_gqz124_pls =  -two*Ca/omz * Lmu  

    ! scale independent pieces
    gq_qg_sub_gqz124_pls = gq_qg_sub_gqz124_pls + four*Ca*lomz/omz + two*Ca*w51tc*log(eta41log/two)/omz                 ! pls



    
  end function gq_qg_sub_gqz124_pls


  function gq_qg_sub_gq1z24(shat,mu2,z,eta41,eta42)
    ! Integrated subtraction terms to multiply FLMqg[1,z2,4] -- [reg],[pls]
    real(dp), intent(in) :: shat,mu2,z,eta41,eta42
    real(dp)             :: gq_qg_sub_gq1z24
    real(dp)             :: Lmu,omz,lomz,w52tc

    w52tc = eta41*(one+eta42)                          ! w4252tilde52

    Lmu = log(mu2/shat)
    omz = one - z
    lomz = log(omz)

    ! scale dependent pieces
    gq_qg_sub_gq1z24 = Cf*(one+z)* Lmu                                         ! reg
    gq_qg_sub_gq1z24 = gq_qg_sub_gq1z24 -two*Cf/omz * Lmu                      ! pls

    ! scale independent pieces
    gq_qg_sub_gq1z24 = gq_qg_sub_gq1z24 + Cf*(omz - w52tc*(one+z)*log(eta42/two) - two*(one+z)*lomz)    ! reg
    gq_qg_sub_gq1z24 = gq_qg_sub_gq1z24 + two*Cf*w52tc*log(eta42/two)/omz + four*Cf*lomz/omz            ! pls




  end function gq_qg_sub_gq1z24



  function gq_qg_sub_gq1z24_pls(shat,mu2,z,eta41,eta42)
    ! Integrated subtraction terms to multiply FLMqg[1,z2,4] -- pls(z=1)
    real(dp), intent(in) :: shat,mu2,z,eta41,eta42
    real(dp)             :: gq_qg_sub_gq1z24_pls
    real(dp)             :: Lmu,omz,lomz,w52tc

    w52tc = eta41*(one+eta42)                          ! w4252tilde52

    Lmu = log(mu2/shat)
    omz = one - z
    lomz = log(omz)

    ! scale dependent pieces
    gq_qg_sub_gq1z24_pls = -two*Cf/omz * Lmu                      

    ! scale independent pieces
    gq_qg_sub_gq1z24_pls = gq_qg_sub_gq1z24_pls + two*Cf*w52tc*log(eta42/two)/omz + four*Cf*lomz/omz            ! pls


  end function gq_qg_sub_gq1z24_pls





  function gq_qg_sub_gq124(shat,mu2,E4,eta41,eta42,eta41log)
    ! Integrated subtraction terms to multiply FLMqg[1,2,4]
    real(dp), intent(in) :: shat,mu2,eta41,eta42,E4,eta41log
    real(dp)             :: gq_qg_sub_gq124
    real(dp)             :: tt,Lmu,log2,log2E4osqrts
    real(dp)             :: wt4252C45,wt4151C45,li2ometa41,li2ometa42
    

    wt4151C45 = eta42**2*(one+two*eta41)
    wt4252C45 = eta41**2*(one+two*eta42)

!
    Lmu = log(mu2/shat)
    log2 = log(two)
    log2E4osqrts = log(two*E4/sqrt(shat))

    li2ometa41 = real(dilog2(one-eta41))
    li2ometa42 = real(dilog2(one-eta42))
! scale dependent term
    gq_qg_sub_gq124 = Lmu*(-11.0_dp/6.0_dp * Ca +nf/three -three/two*Cf)
    ! scale independent term
    gq_qg_sub_gq124 = gq_qg_sub_gq124  + &
         Cf*(          (39.0_dp + 12.0_dp *li2ometa42 + 18.0_dp *log2 - 24.0_dp *log2*Log2E4osqrts + &
         6.0_dp *Log2E4osqrts**2 - 6.0_dp *pisq - 9.0_dp *wt4151C45*Log(eta41log) + &
         12.0_dp *Log2E4osqrts*wt4151C45*Log(eta41log) + 9.0_dp *wt4252C45*Log(eta41log) -&
         12.0_dp *Log2E4osqrts*wt4252C45*Log(eta41log) + &
         18.0_dp *Log(eta42) - 12.0_dp *Log2E4osqrts*Log(eta42) + 9.0_dp *wt4151C45*Log(eta42) -&
         12.0_dp *Log2E4osqrts*wt4151C45*Log(eta42) - &
         9.0_dp *wt4252C45*Log(eta42) + 12.0_dp *Log2E4osqrts*wt4252C45*Log(eta42))/6.0_dp  ) + &
         Ca*(        (6.0_dp *li2ometa41 - 6.0_dp *li2ometa42 + Log2E4osqrts + two*pisq + &
         10.0_dp *Log(eta41log) - 6.0_dp *Log2E4osqrts*Log(eta41log) - &
         9.0_dp *Log(eta42) + 6.0_dp *Log2E4osqrts*Log(eta42))/6.0_dp  ) + &
         nf*((-Log2E4osqrts - Log(eta41log))/6.0_dp)




  end function gq_qg_sub_gq124

   function gq_qg_sub_qqbz124(shat,mu2,z,eta41,eta42,eta41log)
     ! Integrated subtraction terms to multiply FLMqqb[z1,2,4] -- regular only
     real(dp), intent(in) :: shat,mu2,z,eta41,eta42,eta41log
     real(dp)             :: gq_qg_sub_qqbz124
     real(dp)             :: Lmu,w51tc,omz,lomz
     
     w51tc = eta42*(one+eta41)                          ! w4151tilde51
     Lmu = log(mu2/shat)
     omz = one-z
     lomz = log(omz)

     gq_qg_sub_qqbz124  = half*(-one + two*z - two*z**2)*Lmu
     gq_qg_sub_qqbz124 =  gq_qg_sub_qqbz124  + &
          half*(two*z*omz + w51tc*(one - two*z + two*z**2)*Log(eta41log/two) + two*(one-two*z+two*z**2)*lomz)


   end function gq_qg_sub_qqbz124



  function gq_qg_sub_qqbz1zbar2(shat,mu2,z,zbar)
! Integrated subtraction terms to multiply FLMqqb[z1,zbar2]
! output indices: 0=reg(z,zbar), 1=reg(z)pls(zbar)
    real(dp), intent(in) :: shat,mu2,z,zbar
    real(dp)             :: gq_qg_sub_qqbz1zbar2(0:1)
    real(dp)             :: Lmu_zzbar,Lmu_z
    real(dp)             :: omz,omzbar,lomz,lomzbar
    integer              :: j

    Lmu_zzbar  = log(mu2/shat) + log(z) + log(zbar)
    Lmu_z  = log(mu2/shat) + log(z) 

    omz     = one-z
    omzbar  = one-zbar
    lomz    = log(omz)
    lomzbar = log(omzbar)

    ! scale-dependent terms
    
    !L^2
    gq_qg_sub_qqbz1zbar2(0) =  Cf*(one-two*z+two*z**2)*(one+zbar**2)/two/omzbar * Lmu_zzbar**2
    gq_qg_sub_qqbz1zbar2(1) =  -Cf*(one-two*z+two*z**2)/omzbar * Lmu_z**2

    !L
    gq_qg_sub_qqbz1zbar2(0) = gq_qg_sub_qqbz1zbar2(0) - &
         Cf*(one-two*zbar + four*z*zbar -four*z**2*zbar + zbar**2 +&
         two*(one-two*z+two*z**2)*(one+zbar**2)*(lomz+lomzbar))/two/omzbar * Lmu_zzbar

    gq_qg_sub_qqbz1zbar2(1) = gq_qg_sub_qqbz1zbar2(1) + &
         two*Cf*(omz*z + (one-two*z+two*z**2)*(lomz+lomzbar))/omzbar * Lmu_z

    ! scale-independent terms

    gq_qg_sub_qqbz1zbar2(0) = gq_qg_sub_qqbz1zbar2(0) + &
         Cf*(omz*z+(one-two*z+two*z**2)*lomz)*(omzbar**2 + two*(one+zbar**2)*lomzbar)/omzbar
    gq_qg_sub_qqbz1zbar2(1) = gq_qg_sub_qqbz1zbar2(1) - &
         four*Cf*(omz*z + (one-two*z+two*z**2)*lomz)*lomzbar/omzbar

    gq_qg_sub_qqbz1zbar2 = gq_qg_sub_qqbz1zbar2 * tag_nontc * tag_z1zb2

  end function gq_qg_sub_qqbz1zbar2
    
    
    
    

   function gq_qg_sub_qqbz12(shat,mu2,z,yy)
! Integrated subtraction terms to multiply FLM[1,z2]
     real(dp), intent(in) :: shat,mu2,z,yy(1:4)
     real(dp)             :: gq_qg_sub_qqbz12
     real(dp)             :: Lmu,li3z,li3mz,li3zoopz,li3omzsq
     real(dp)             :: omz,opz,lz,lomz,lopz,li2z,li2omz,li3omz,li2mz
     real(dp)             :: tripcoll(1:3),coeffL2,coeffL,coeffL0
     

     Lmu = log(mu2/shat) + log(z)

     omz = one-z
     opz = one+z
     lomz = log(omz)
     lopz = log(opz)
     lz = log(z)
     li2z = real(dilog2(z))
     li2mz = real(dilog2(-z))
     li2omz = real(dilog2(omz))
     li3omz = trilog(omz)   
     li3z     = trilog(z)
     li3mz    = trilog(-z)
     li3zoopz = trilog(z/(one+z))
     li3omzsq = trilog(one-z**2)

     !Lmu^2
     coeffL2 =  Ca*(8 + 6*z + 12*lomz*z + 12*lz*z + 48*z**2 - &
          24*lomz*z**2 + 48*lz*z**2 - 62*z**3 + 24*lomz*z**3)/(24.0_dp*z) +&
          Cf*(15*z + 12*lomz*z - 6*lz*z - 24*z**2 - 24*lomz*z**2 + 12*lz*z**2 + &
          36*z**3 + 24*lomz*z**3 - 24*lz*z**3)/(24.0_dp*z)
          
     coeffL = -Ca*(  26 + 24*lomz + 24*lz - 39*z - 36*li2z*z + 18*lomz*z + 27*lomz**2*z + &
          45*lz*z + 9*lz**2*z + 93*z**2 - 144*li2z*z**2 + 144*lomz*z**2 - &
          54*lomz**2*z**2 + 144*lz*z**2 + 54*lz**2*z**2 + 30*pisq*z**2 - &
          71*z**3 - 186*lomz*z**3 + 54*lomz**2*z**3 + 132*lz*z**3 - &
          12*pisq*z**3 - 18*li2mz*z*(1 + 2*z + 2*z**2) - 18*lz*z*Log(1 + z) - &
          36*lz*z**2*Log(1 + z) - 36*lz*z**3*Log(1 + z)  )/(18.0_dp*z) + &
          Cf*(         li2z - 2*li2z*z + (-48 - 4*pisq + 123*z + 8*pisq*z - 90*z**2 - &
          4*pisq*z**2 + lz**2*(-3 + 6*z) + lz*(-9 + 60*z - 60*z**2) - & 
          18*lomz**2*(1 - 2*z + 2*z**2) + &
          6*lomz*(1 - 8*z + 4*z**2 + lz*(2 - 4*z + 4*z**2)))/12.0_dp   )

     coeffL0 = Ca*( 15.0_dp/four + li2mz/two + li2omz/two - 4*li2z - 4*li3omz + (71*log2)/12.0_dp + &
          2*li2mz*log2 + 2*li2omz*log2 + 2*li2z*log2 - (7*log2**2)/8.0_dp + &
          (5*lomz)/three + (3*li2omz*lomz)/two - 4*li2z*lomz - (3*log2*lomz)/two - &
          (7*log2**2*lomz)/four + lomz**2 - (3*log2*lomz**2)/two + (13*lomz**3)/12.0_dp + &
          (3*lz)/four - (7*log2*lz)/two - (7*log2**2*lz)/four - 2*lomz**2*lz - lz**2/four -& 
          (3*log2*lz**2)/two - lz**3/three + (7*pisq)/12.0_dp + (log2*pisq)/6.0_dp + &
          (lomz*pisq)/three - (lz*pisq)/three - 11/(9.0_dp*z) - (61*log2)/(18.0_dp*z) - &
          (7*log2**2)/(6.0_dp*z) - (4*log2*lomz)/(three*z) + (4*lomz**2)/(three*z) - lz/z - &
          (10*log2*lz)/(three*z) - (2*lz**2)/(three*z) - (2*pisq)/(9.0_dp*z) - (17*z)/12.0_dp + &
          li2mz*z - 4*li2z*z - 16*li3omz*z - (85*log2*z)/6.0_dp + 4*li2mz*log2*z - &
          4*li2omz*log2*z + 8*li2z*log2*z - 7*log2**2*z - (83*lomz*z)/12.0_dp - &
          3*li2omz*lomz*z - 16*li2z*lomz*z - 8*log2*lomz*z + (7*log2**2*lomz*z)/two + &
          (15*lomz**2*z)/two + 3*log2*lomz**2*z - (13*lomz**3*z)/6.0_dp + 2*lz*z - &
          17*log2*lz*z - 7*log2**2*lz*z - 8*lomz**2*lz*z - (5*lz**2*z)/two - &
          8*log2*lz**2*z - (4*lz**3*z)/three - (pisq*z)/two - (5*log2*pisq*z)/three + & 
          (10*lomz*pisq*z)/three - (4*lz*pisq*z)/three - (29*z**2)/18.0_dp + li2mz*z**2 - &
          (44*li2z*z**2)/three + (419*log2*z**2)/36.0_dp + 4*li2mz*log2*z**2 + &
          4*li2omz*log2*z**2 + (217*log2**2*z**2)/24.0_dp + 6*lomz*z**2 + &
          3*li2omz*lomz*z**2 + (31*log2*lomz*z**2)/three - (7*log2**2*lomz*z**2)/two - &
          (59*lomz**2*z**2)/6.0_dp - 3*log2*lomz**2*z**2 + (13*lomz**3*z**2)/6.0_dp - &
          (11*lz*z**2)/three - (11*log2*lz*z**2)/three + (11*lz**2*z**2)/three + &
          (25*pisq*z**2)/6.0_dp + log2*pisq*z**2 - (2*lomz*pisq*z**2)/three + &
          ((1 + 4*log2)*lz*(1 + 2*z + 2*z**2)*Log(1 + z))/2.0_dp  )  +  &
          Cf*( -207 - 12*li2omz - 60*li3omz - 54*log2 - 48*li2omz*log2 + 24*li2z*log2 + &
          96*log2**2 + 192*lomz - 36*li2omz*lomz - 48*li2z*lomz + 84*log2*lomz - &
          114*log2**2*lomz - 78*lomz**2 - 108*log2*lomz**2 + 22*lomz**3 + 33*lz + &
          48*log2*lz - 21*log2**2*lz + 60*lomz*lz - 18*lomz**2*lz + 15*lz**2 - &
          42*log2*lz**2 - 24*lomz*lz**2 - 4*lz**3 + 9*pisq + 12*lomz*pisq - &
          4*lz*pisq + 669*z + 120*li3omz*z + 300*log2*z + 96*li2omz*log2*z - &
          48*li2z*log2*z - 213*log2**2*z - 570*lomz*z + 72*li2omz*lomz*z + &
          96*li2z*lomz*z - 264*log2*lomz*z + 228*log2**2*lomz*z + 240*lomz**2*z + &
          216*log2*lomz**2*z - 44*lomz**3*z + 60*lz*z + 24*log2*lz*z + &
          42*log2**2*lz*z - 144*lomz*lz*z + 36*lomz**2*lz*z - 6*lz**2*z +& 
          84*log2*lz**2*z + 48*lomz*lz**2*z + 8*lz**3*z - 14*pisq*z - &
          24*lomz*pisq*z + 8*lz*pisq*z - 588*z**2 - 24*li3omz*z**2 - &
          234*log2*z**2 - 96*li2omz*log2*z**2 + 171*log2**2*z**2 + 504*lomz*z**2 - &
          72*li2omz*lomz*z**2 + 216*log2*lomz*z**2 - 228*log2**2*lomz*z**2 - &
          192*lomz**2*z**2 - 216*log2*lomz**2*z**2 + 44*lomz**3*z**2 - &
          72*lz*z**2 + 144*lomz*lz*z**2 + 12*lomz**2*lz*z**2 + 24*lz**2*z**2 - &
          48*log2*lz**2*z**2 - 48*lomz*lz**2*z**2 + 6*pisq*z**2 + &
          8*log2*pisq*z**2 + 8*lomz*pisq*z**2 + &
          12*(1 - 2*z + 2*z**2)*Log(-one + one/z) - &
          6*(1 - 2*z + 2*z**2)*Log(-one + one/z)**2 + & 
          192*zeta3 - 384*z*zeta3 + 384*z**2*zeta3   )/24.0_dp

     gq_qg_sub_qqbz12 =coeffL2 * Lmu**2 + coeffL * Lmu + coeffL0

     gq_qg_sub_qqbz12 = gq_qg_sub_qqbz12 * tag_nontc


! get the triple collinear splitting functions 
     if (withmax .and. tcanalytic) then
        tripcoll(3) = -Cf*TR*(   lomz + 5*lomz*z + pisq*(-one/three + z - z**2) - (2*lomz*pisq*(1 - 2*z + 2*z**2))/3.0_dp + &
             li3omz*(4 - 8*z + 8*z**2) + li3z*(9 - 18*z + 14*z**2) - (lz**3*(19 - 38*z + 28*z**2))/12.0_dp + lz**2*(9.0_dp/8.0_dp + z + lomz*(half - z + z**2)) + &
             (log2**2*(32 - 71*z + 57*z**2 + 7*lz*(-1 + 2*z) - 38*lomz*(1 - 2*z + 2*z**2)))/4.0_dp +&
             li2z*(-2 + lz*(-1 + 2*z - 2*z**2) + 2*log2*(3 - 6*z + 4*z**2) + lomz*(4 - 8*z + 8*z**2)) + lz*(three/four + lomz - (7*z)/four -&
             (pisq*(5 - 10*z + 8*z**2))/3.0_dp + lomz**2*(4 - 8*z + 8*z**2)) + &
             log2*((-2*pisq*(-1 + z)**2)/3.0_dp + 2*lz*(2 + z) + (-9 + 50*z - 39*z**2)/2.0_dp + lz**2*(-7.0_dp/two + 7*z - 4*z**2) - &
             9*lomz**2*(1 - 2*z + 2*z**2) + lomz*lz*(4 - 8*z + 8*z**2) + lomz*(7 - 22*z + 18*z**2)) + 4*z*(-2 + zeta3) - 2*zeta3 ) - &
             Ca*TR*( 4*li3z*(1 + 8*z) - (lz**3*(5 + 34*z))/6.0_dp + li3omz*(-9 + 2*z - 18*z**2) + 6*li3zoopz*(1 + 2*z + 2*z**2) + &
             li3omzsq*(2 + 4*z + 4*z**2) + li3mz*(5 + 26*z + 10*z**2) + (4*lomz*(-13 + 3*(9 + pisq)*z - 3*(27 + pisq)*z**2 + (67 + 6*pisq)*z**3))/(9.0_dp*z) + &
             li2z*(-1 + 20/(3.0_dp*z) + 32*z - 4*lz*z - 8*log2*(-3 + z)*z - (22*z**2)/3.0_dp + lomz*(-4 + 8*z - 8*z**2)) - &
             (7*log2**2*(4 + 3*z + 24*z**2 - 31*z**3 + 6*lz*z*(1 + 4*z) + 6*lomz*z*(1 - 2*z + 2*z**2)))/(12.0_dp*z) + &
             li2mz*(5 + 8/(3.0_dp*z) + 16*z + (44*z**2)/3.0_dp + lz*(1 - 6*z + 2*z**2) + log2*(4 + 8*z + 8*z**2) + lomz*(4 + 8*z + 8*z**2)) +&
             (-236 + 6*pisq*(-4 + 3*z - 18*z**2 + 11*z**3) + z**3*(398 - 270*zeta3) - 3*z*(7 + 81*zeta3) - 3*z**2*(47 + 144*zeta3))/(27.0_dp*z) + &
             (pisq*(1 + 2*z + 2*z**2)*Log(1 + z))/two - (1 + 2*z + 2*z**2)*Log(1 + z)**3 + &
             lz**2*(-13.0_dp/four - 4/z - 22*z - (11*z**2)/three + (half + z + z**2)*Log(1 + z)) + &
             lz*(-41.0_dp/6.0_dp - 82/(9.0_dp*z) - 31*z - (134*z**2)/9.0_dp - (2*pisq*(1 + 6*z))/three + lomz**2*(-4 + 8*z - 8*z**2) + &
             (5 + 8/(3.0_dp*z) + 16*z + (44*z**2)/3.0_dp)*Log(1 + z) + (3 + 6*z + 6*z**2)*Log(1 + z)**2 + &
             lomz*(-3 + 4/(3.0_dp*z) + 8*z - (22*z**2)/3.0_dp + (4 + 8*z + 8*z**2)*Log(1 + z))) - &
             (log2*(134 - 3*(83 + 6*pisq)*z + 42*(13 + 2*pisq)*z**2 - (431 + 60*pisq)*z**3 + 18*lz**2*z*(3 + 16*z) +&
             54*lomz**2*z*(1 - 2*z + 2*z**2) + 6*lomz*(8 + 9*z + 48*z**2 - 62*z**3) + &
             6*lz*(20 + 21*z + 102*z**2 + 22*z**3 + 3*lomz*z*(4 - 8*z + 8*z**2) - 3*z*(4 + 8*z + 8*z**2)*Log(1 + z)))) /(18.0_dp*z))
!        print *, "z",z
!        print *, "tc" ,tripcoll(3)
!        pause
        

     else
        call get_tripcoll_subtr_g_qg_qstar(z,yy,tripcoll) 
     endif
     gq_qg_sub_qqbz12  = gq_qg_sub_qqbz12*tag_nontc + tripcoll(3)*tag_tcreg

     gq_qg_sub_qqbz12 = gq_qg_sub_qqbz12 * tag_z12
     


     
   end function gq_qg_sub_qqbz12

   subroutine get_tripcoll_subtr_g_qg_qstar(z,yy,tripcoll)
    implicit none
    real(dp), intent(in)  :: z,yy(1:4)
    real(dp), intent(out) :: tripcoll(1:3)
    real(dp), parameter   :: zmin = zero
    integer, parameter :: myep = 1
    integer, parameter :: mycontr = 3 !-- 1 -> delta, 2-> plus, 3-> finite part
    integer, parameter :: mycol = 3   !-- 1 -> abelian, 2-> non-abelian, 3-> full color
    real(dp)              :: x1,x2,x3,xx2,xx(1:3),xx3,x3a
    real(dp)              :: seca(1:3,1:3,0:1),secb(1:3,1:3,0:1),secc(1:3,1:3,0:1),secd(1:3,1:3,0:1)
    real(dp)              :: seca0(1:3,1:3,0:1),secb0(1:3,1:3,0:1),secc0(1:3,1:3,0:1),secd0(1:3,1:3,0:1)
    real(dp)              :: wghtac,wghtbd
    real(dp)              :: res_tc(1:3,1:3,0:1)

    tripcoll = zero

    seca = zero; seca0 = zero
    secb = zero; secb0 = zero
    secc = zero; secc0 = zero
    secd = zero; secd0 = zero

    x1  = yy(1)  ! r
    x2  = yy(2)  ! x4
    x3a = yy(3)  ! lambda
    xx3 = yy(4)  ! x3


    ! remap lambda
    x3 = sin( half * pi * x3a )**2

    if (x1*x2*(one-z)*xx3 .lt. bufftc) return

    xx = (/x1,x2,x3/)


    seca0 = PTC_A_g_qg_qstar(z,zmin,xx,zero)
    secc0 = PTC_C_g_qg_qstar(z,zmin,xx,zero)
    if (withmax) then
       seca = PTC_A_g_qg_qstar(z,zmin,xx,xx3)
       secc = PTC_C_g_qg_qstar(z,zmin,xx,xx3)
    endif


    wghtac = -one/x1/x2

    !-- 

    xx2 = x2**2
    xx = (/x1,xx2,x3/)


    secb0 = PTC_B_g_qg_qstar(z,zmin,xx,zero)
    secd0 = PTC_D_g_qg_qstar(z,zmin,xx,zero)

    if (withmax) then
       secb = PTC_B_g_qg_qstar(z,zmin,xx,xx3)
       secd = PTC_D_g_qg_qstar(z,zmin,xx,xx3)
    endif

    wghtbd = -one/x1/x2*two

    !--

    !-- debug
    !seca = 0
    !secb = 0
    !secc = 0
    !secd = 0
    !-- end debug

    !-- color/[delta,plus,regular]/[finite,O(ep)]


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

    ! remove color factors from Abelian and non-Abelian parts
    res_tc(1,1:3,0:1) = res_tc(1,1:3,0:1)/(Cf**2)
    res_tc(2,1:3,0:1) = res_tc(2,1:3,0:1)/(Ca*Cf)
    res_tc(3,1:3,0:1) = -res_tc(3,1:3,0:1)*tr/Cf

    tripcoll(:) = res_tc(:,3,1)


  end subroutine get_tripcoll_subtr_g_qg_qstar
  
end module mod_subtrfn_nnlo_dy_gq

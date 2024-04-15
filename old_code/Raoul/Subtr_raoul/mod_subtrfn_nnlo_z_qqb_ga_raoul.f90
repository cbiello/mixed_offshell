module mod_subtrfn_nnlo_z_qqb_ga_raoul
! integrated subtraction terms for qqb->g gamma channel
! cf. Notes from Chiara and Federica
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_auxfunctions
  implicit none

  private 

  public :: qqb_ga_sub_qqbg_zi,qqb_ga_sub_qqbg_zi_pls,qqb_ga_sub_qqbg
  public :: qqb_ga_sub_qqba_zi,qqb_ga_sub_qqba_zi_pls, qqb_ga_sub_qqba
  public :: qqb_ga_sub_qqbz1zbar2,qqb_ga_sub_qqbz12,qqb_ga_sub_qqb12

  
contains

  
  function qqb_ga_sub_qqbg_zi(Ec,mu,z,eta_g1,eta_g1log,Qq)
    ! Integrated subtraction terms to multiply FLMqqbg[z.i,3,4|5g ], cf. Eq. 265
    ! This gives reg+pls
    real(dp), intent(in)   :: Ec,mu,z,eta_g1,eta_g1log,Qq
    real(dp)               :: qqb_ga_sub_qqbg_zi
    real(dp)               :: Lmu,omz,lomz,wtgi

    Lmu = log(mu/two/Ec)
    omz = one-z
    lomz = log(omz)
    wtgi = one-eta_g1        ! wtilde^{gamma i, g i}_{gamma || i} = etag2 = 1-etag1, Eq. 132
    

! scale dependent pieces

    qqb_ga_sub_qqbg_zi = two*(one + z)*Lmu                      ! reg
    qqb_ga_sub_qqbg_zi = qqb_ga_sub_qqbg_zi - four/omz*Lmu      ! pls

! scale independent pieces
    qqb_ga_sub_qqbg_zi = qqb_ga_sub_qqbg_zi + &
         omz - two*(one+z)*lomz - log(eta_g1log)*wtgi*(one+z)    ! reg
    qqb_ga_sub_qqbg_zi = qqb_ga_sub_qqbg_zi + &
         four*lomz/omz + two/omz*wtgi*log(eta_g1log)             ! pls


    qqb_ga_sub_qqbg_zi = qqb_ga_sub_qqbg_zi * Qq**2
    


  end function qqb_ga_sub_qqbg_zi


  
  function qqb_ga_sub_qqbg_zi_pls(Ec,mu,z,eta_g1,eta_g1log,Qq)
    ! Integrated subtraction terms to multiply FLMqqbg[z.i,3,4|5g ], cf. Eq. 265
    ! This gives [pls(z=1)]
    real(dp), intent(in)   :: Ec,mu,z,eta_g1,eta_g1log,Qq
    real(dp)               :: qqb_ga_sub_qqbg_zi_pls
    real(dp)               :: Lmu,omz,lomz,wtgi

    Lmu = log(mu/two/Ec)
    omz = one-z
    lomz = log(omz)
    wtgi = one-eta_g1        ! wtilde^{gamma i, g i}_{gamma || i} = etag2 = 1-etag1, Eq. 132
    

! scale dependent pieces
    qqb_ga_sub_qqbg_zi_pls = - four/omz*Lmu      ! pls

! scale independent pieces
    qqb_ga_sub_qqbg_zi_pls = qqb_ga_sub_qqbg_zi_pls + &
         four*lomz/omz + two/omz*wtgi*log(eta_g1log)             ! pls


    qqb_ga_sub_qqbg_zi_pls = qqb_ga_sub_qqbg_zi_pls * Qq**2

    
  end function qqb_ga_sub_qqbg_zi_pls


  function qqb_ga_sub_qqbg(Ec,mu,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qq,Qe)
    ! Integrated subtraction terms to multiply FLMqqbg[1,2,3,4|5g ], cf. Eq. 266
    real(dp), intent(in)   :: mu,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qq,Qe
    real(dp)               :: qqb_ga_sub_qqbg
    real(dp)               :: Lmu,Gesq,Geq,s13,s23,s14,s24,s34

    Lmu = log(mu/two/Ec)

    ! sij
    s13 = four*E1*E3*eta31
    s23 = four*E2*E3*eta32
    s14 = four*E1*E4*eta41
    s24 = four*E2*E4*eta42
    s34 = four*E3*E4*eta34


    ! functions Gesq and Geq
    Gesq = fnGesq(Ec,E3,E4,eta34)
    Geq  = fnGeq(mu,E3,E4,eta31,eta32,eta41,eta42,s13,s14,s23,s24)

    qqb_ga_sub_qqbg = Qq**2*(two/three*pisq - 6.0_dp*Lmu) + Qe*Qq*Geq + Qe**2*Gesq
    

  end function qqb_ga_sub_qqbg
  




  function qqb_ga_sub_qqba_zi(Ec,mu,z,eta51,eta52,eta53,eta54,eta51log)
    ! Integrated subtraction terms to multiply FLMqqba[z.i,3,4|5gamma ], cf. Eq. 265
    ! This gives reg+pls
    real(dp), intent(in)   :: Ec,mu,z,eta51,eta52,eta53,eta54,eta51log
    real(dp)               :: qqb_ga_sub_qqba_zi
    real(dp)               :: Lmu,omz,lomz,wtai
    Lmu = log(mu/two/Ec)
    omz = one-z
    lomz = log(omz)
    wtai = eta52*eta53*eta54/((eta51+eta54)*eta52*eta53+(eta52+eta53)*eta51*eta54)        ! wtilde^{amma i, g i}_{g || i}, Eq. 124

    

! scale dependent pieces

    qqb_ga_sub_qqba_zi = two*(one + z)*Lmu                      ! reg
    qqb_ga_sub_qqba_zi = qqb_ga_sub_qqba_zi - four/omz*Lmu      ! pls

! scale independent pieces
    qqb_ga_sub_qqba_zi = qqb_ga_sub_qqba_zi + &
         omz - two*(one+z)*lomz - log(eta51log)*wtai*(one+z)    ! reg
    qqb_ga_sub_qqba_zi = qqb_ga_sub_qqba_zi + &
         four*lomz/omz + two/omz*wtai*log(eta51log)             ! pls


    qqb_ga_sub_qqba_zi = qqb_ga_sub_qqba_zi * Cf



  end function qqb_ga_sub_qqba_zi


  
  function qqb_ga_sub_qqba_zi_pls(Ec,mu,z,eta51,eta52,eta53,eta54,eta51log)
    ! Integrated subtraction terms to multiply FLMqqba[z.i,3,4|5gamma ], cf. Eq. 265
    ! This gives [pls(z=1)]
    real(dp), intent(in)   :: Ec,mu,z,eta51,eta52,eta53,eta54,eta51log
    real(dp)               :: qqb_ga_sub_qqba_zi_pls
    real(dp)               :: Lmu,omz,lomz,wtai

    Lmu = log(mu/two/Ec)
    omz = one-z
    lomz = log(omz)
    wtai = eta52*eta53*eta54/((eta51+eta54)*eta52*eta53+(eta52+eta53)*eta51*eta54)        ! wtilde^{gamma i, g i}_{g || i}, Eq. 124
    

! scale dependent pieces
    qqb_ga_sub_qqba_zi_pls = - four/omz*Lmu      ! pls

! scale independent pieces
    qqb_ga_sub_qqba_zi_pls = qqb_ga_sub_qqba_zi_pls + &
         four*lomz/omz + two/omz*wtai*log(eta51log)             ! pls


    qqb_ga_sub_qqba_zi_pls = qqb_ga_sub_qqba_zi_pls * Cf

    
  end function qqb_ga_sub_qqba_zi_pls



  function qqb_ga_sub_qqba(Ec,mu)
    ! Integrated subtraction terms to multiply FLMqqba[1,2,3,4|5gamma ], cf. Eq. 266
    real(dp), intent(in)   :: mu,Ec
    real(dp)               :: qqb_ga_sub_qqba
    real(dp)               :: Lmu,Gesq,Ge

    Lmu = log(mu/two/Ec)

    qqb_ga_sub_qqba = Cf*(two/three*pisq - 6.0_dp*Lmu)


  end function qqb_ga_sub_qqba
  



!!  function qqb_ga_sub_qqb124(shat,mu2,eta41,eta42,eta41log,eta42log,E4,myesq)
!!    ! Integrated subtraction terms to multiply FLM[1,2,4],  cf. Eq. 10.10
!!    ! This gives [pls(z=1)]
!!    real(dp), intent(in)   :: shat,mu2,eta41,eta42,eta41log,eta42log,E4
!!    real(dp), intent(in)   :: myesq
!!    real(dp)               :: qqb_ga_sub_qqb124
!!    real(dp)               :: Lmu,w4151t51,w4252t52,w4151t45,w4252t45,le4oe1
!!
!!    Lmu = log(mu2/shat)
!!    
!!    w4151t51 = eta42*(one+eta41)
!!    w4252t52 = eta41*(one+eta42)
!!    !
!!    w4151t45 = eta42**2*(one+two*eta41)
!!    w4252t45 = eta41**2*(one+two*eta42)
!!
!!!    le4oe1 = half*log(E4**2/shat) + log2   ! log(E4/E1)
!!
!!!!    ! log term
!!!!     qqb_ga_sub_qqb124 = -three*myesq*Lmu*tag_ga_gg
!!!!
!!!!     ! scale independent term
!!!!     qqb_ga_sub_qqb124 = qqb_ga_sub_qqb124 + &
!!!!          ( tag_ga_gg*(137.0_dp*Ca - 6.0_dp*Ca*le4oe1 + 18.0_dp*Ca*le4oe1**2 + 126.0_dp*Ca*log2 + 36.0_dp*Ca*le4oe1*log2 - 24.0_dp*Ca*pisq +&
!!!!          12.0_dp*myesq*pisq  + 18.0_dp*Ca*dilog2(eta41) + 18.0_dp*Ca*dilog2(eta42) -&
!!!!          33.0_dp*Ca*w4151t45*Log(eta41log/eta42log)  - 33.0_dp*Ca*w4252t45*Log(eta42log/eta41log) + &
!!!!          30.0_dp*Ca*Log(eta41log/two) + 18.0_dp*Ca*le4oe1*Log(eta41log/two) - 36.0_dp*myesq*le4oe1*w4151t51*Log(eta41log/two) + &
!!!!          30.0_dp*Ca*Log(eta42log/two) + 18.0_dp*Ca*le4oe1*Log(eta42log/two) - 36.0_dp*myesq*le4oe1*w4252t52*Log(eta42log/two)) + &
!!!!          Nf*TR*tag_ga_nf*(-52.0_dp + 12.0_dp*le4oe1 - 36.0_dp*log2+ 12.0_dp*w4151t45*Log(eta41log/eta42log) +&
!!!!          12.0_dp*w4252t45*Log(eta42log/eta41log) -  6.0_dp*Log(eta41log/two) -6.0_dp*Log(eta42log/two) ) )/18.0_dp
!!!!
!!
!!
!!
!!   end function qqb_ga_sub_qqb124
         


  function qqb_ga_sub_qqbz1zbar2(mu,Ec,z,zbar,Qq,Qe)
! Integrated subtraction terms to multiply FLM[z1,zbar2] cf. Eq. 265
! log(musq/s) evaluated for {Ec,Ec/zbar,Ec/z,Ec/z/zbar} 
! output indices: coensicient of {z1 zbar2,z12, 1z2, 12}
    implicit none
    real(dp), intent(in)  :: mu,Ec,z,zbar,Qe,Qq
    real(dp)              :: qqb_ga_sub_qqbz1zbar2(1:4)
    real(dp)              :: coeff_noL(1:4),coeff_L(1:4),coeff_Lsq(1:4)
    real(dp)              :: omz,omzbar,lomz,lomzbar,Lmu(1:4)
    integer               :: j


    Lmu(4)  = log(mu/two/Ec)
    Lmu(3)  = log(mu/two/Ec) + log(zbar)
    Lmu(2)  = log(mu/two/Ec) + log(z)
    Lmu(1)  = log(mu/two/Ec) + log(z) + log(zbar)

    omz     = one-z
    omzbar  = one-zbar
    lomz    = log(omz)
    lomzbar = log(omzbar)

!! -- musq-independent  terms -- !!

    coeff_noL(1) = two*Cf*Qq**2*(omz**2+two*(one+z**2)*lomz)*(omzbar**2+two*(one+zbar**2)*lomzbar)/omz/omzbar
    coeff_noL(2) = -8.0_dp*Cf*Qq**2*(omz**2+two*(one+z**2)*lomz)*lomzbar/omz/omzbar
    coeff_noL(3) = -8.0_dp*Cf*Qq**2*lomz*(omzbar**2+two*(one+zbar**2)*lomzbar)/omz/omzbar
    coeff_noL(4) = 32.0_dp*Cf*Qq**2*lomz*lomzbar/omz/omzbar
    

!! -- log terms  -- !!

    coeff_L(1) = -8.0_dp*Cf*Qq**2*( (one-z+z**2)*(one-zbar+zbar**2) - z*zbar + &
         (lomz+lomzbar)*(one + z**2)*(one + zbar**2) ) /omz/omzbar
    coeff_L(2) = 8.0_dp*Cf*Qq**2*(omz**2+two*(lomz+lomzbar)*(one+z**2))/omz/omzbar
    coeff_L(3) = 8.0_dp*Cf*Qq**2*(omzbar**2+two*(lomz+lomzbar)*(one+zbar**2))/omz/omzbar
    coeff_L(4) = -32.0_dp*Cf*Qq**2*(lomz+lomzbar)/omz/omzbar
    
!! -- log^2 terms  -- !!

    coeff_Lsq(1) = 8.0_dp*Cf*Qq**2*(one+z**2)*(one+zbar**2)/omz/omzbar
    coeff_Lsq(2) = -16.0_dp*Cf*Qq**2*(one+z**2)/omz/omzbar
    coeff_Lsq(3) = -16.0_dp*Cf*Qq**2*(one+zbar**2)/omz/omzbar
    coeff_Lsq(4) = 32.0_dp*Cf*Qq**2/omz/omzbar


    do j=1,4
       qqb_ga_sub_qqbz1zbar2(j) = coeff_noL(j) + coeff_L(j) * Lmu(j) + coeff_Lsq(j) * Lmu(j)**2
    enddo


    
  end function qqb_ga_sub_qqbz1zbar2

  
  

  function qqb_ga_sub_qqbz12(mu,z,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qq,Qe)
    ! Integrated subtraction terms to multiply FLM[z1,2] and FLM[1,z2], cf. Eq. 265
    ! log(musq/s) evaluated for s={shat,shat/z}
    ! Output indices: {0=reg+pls,1=pls}
    implicit none
    real(dp), intent(in) :: mu,z,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,Qq,Qe
    real(dp)             :: qqb_ga_sub_qqbz12(0:1)
    real(dp)             :: omz,lz,lomz,li2omz,li2z,li3omz,li3z,Lmu,Lmu_z
    real(dp)             :: s13,s23,s14,s24,s34,logeta4
    real(dp)             :: Gesq,Geq
    real(dp)             :: pls1_Lsq,pls0_Lsq,pls_Lsq,reg_Lsq
    real(dp)             :: pls2_Lmu,pls1_Lmu,pls0_Lmu,pls_Lmu,reg_Lmu
    real(dp)             :: pls3_noL,pls2_noL,pls1_noL,pls0_noL,pls_noL,reg_noL

    ! sij

    s13 = four*E1*E3*eta31
    s23 = four*E2*E3*eta32
    s14 = four*E1*E4*eta41
    s24 = four*E2*E4*eta42
    s34 = four*E3*E4*eta34

    ! these are the logs of mu that also depend on the initial state energy Ec
    ! the factor of z needs to be taken into correctly
    ! note that these are logs of mu not musq!
    Lmu   = log(mu/two/Ec)
    Lmu_z = log(mu/two/Ec) + log(z)

    ! functions Gesq and Geq
    Gesq = fnGesq(Ec,E3,E4,eta34)
    Geq  = fnGeq(mu,E3,E4,eta31,eta32,eta41,eta42,s13,s14,s23,s24)
    
    omz    = one-z
    lz   = log(z)
    lomz   = log(omz)
    li2omz = dilog2(omz)
    li2z = dilog2(z)
    li3omz = trilog(omz)
    li3z   =  trilog(+z)
    logeta4 = log(eta31*eta42/eta32/eta41)

! Lmu^2
    
    ! plus distributions
    pls1_Lsq = 32.0_dp*Qq**2
    pls0_Lsq = 48.0_dp*Qq**2 + 8.0_dp*logeta4*Qe*Qq

    pls_Lsq = pls1_Lsq * lomz/omz + pls0_Lsq / omz

    ! reg. function
    reg_Lsq= -four*logeta4*Qe*Qq - 32.0_dp*Qq**2 - 16.0_dp*lomz*Qq**2 - (four*lz*Qq**2)/omz - &
         four*logeta4*Qe*Qq*z - 16.0_dp*Qq**2*z - 16.0_dp*lomz*Qq**2*z - (12.0_dp*lz*Qq**2*z**2)/omz

! Lmu
    
    ! plus distributions
    pls2_Lmu = -48.0_dp*Qq**2
    pls1_Lmu = -8.0_dp*logeta4*Qe*Qq - 48.0_dp*Qq**2
    pls0_Lmu = -four*Gesq*Qe**2 - four*Qe*Qq*Geq + 8.0_dp*logeta4*Qe*Qq*ln2 

    pls_Lmu = pls2_Lmu * lomz**2/omz + pls1_Lmu * lomz / omz + pls0_Lmu / omz 

    ! reg. function
    reg_Lmu =          Qe**2*(2*Gesq + 2*Gesq*z) + Qq**2*(-24 + 8*li2z + 32*lomz + 24*lomz**2 - 2*lz**2 - (8*li2omz)/omz - (8*lz)/omz + &
          (8*lomz*lz)/omz - (4*pisq)/3. + 20*z + 8*li2z*z + 16*lomz*z + 24*lomz**2*z - 2*lz**2*z - (8*lz*z)/omz - (4*pisq*z)/3. - &
          (8*li2omz*z**2)/omz - (8*lz*z**2)/omz + (8*lomz*lz*z**2)/omz) + &
       Qe*Qq*(2*Geq - 2*logeta4 - 4*ln2*logeta4 + 4*logeta4*lomz + 2*Geq*z + 2*logeta4*z - 4*ln2*logeta4*z + 4*logeta4*lomz*z - &
          4*lz*Log(s13/s14) + (8*lz*Log(s13/s14))/(1 - z) - 4*lz*z*Log(s13/s14))


! Lmu-independent terms

    ! plus distributions    
    pls3_noL = 16.0_dp*Qq**2
    pls2_noL = zero
    pls1_noL = four*Gesq*Qe**2  +  Qe*Qq*(four*Geq - 8.0_dp*ln2*logeta4)
    pls0_noL = Qq**2* 32.0_dp*zeta3

    pls_noL = pls3_noL * lomz**3/omz + pls2_noL * lomz**2/omz + pls1_noL * lomz / omz + pls0_noL / omz 


    ! reg. function

    reg_noL =   (Qe**2*(6*Gesq - 12*Gesq*lomz - 12*Gesq*z + 6*Gesq*z**2 + 12*Gesq*lomz*z**2))/(6.*omz) +  &
       (Qq**2*(-108 + 12*li2omz - 24*li2z - 60*li3omz - 120*li3z + 72*lomz + 36*li2omz*lomz - 48*li2z*lomz - 24*lomz**2 - 48*lomz**3 + &
            30*lz + 24*li2z*lz + 60*lomz*lz - 18*lomz**2*lz + 12*lz**2 - 24*lomz*lz**2 + 15*lz**3 + 20*pisq + 8*lomz*pisq + 16*lz*pisq + & 
            204*z - 24*li2omz*z + 48*li2z*z - 126*lomz*z + 48*lomz**2*z - 84*lz*z + 24*lomz*lz*z - 18*lz**2*z - 40*pisq*z - 96*z**2 + &
            12*li2omz*z**2 - 24*li2z*z**2 + 36*li3omz*z**2 - 72*li3z*z**2 + 54*lomz*z**2 + 36*li2omz*lomz*z**2 + 48*li2z*lomz*z**2 - &
            24*lomz**2*z**2 + 48*lomz**3*z**2 + 54*lz*z**2 + 24*li2z*lz*z**2 + 60*lomz*lz*z**2 + 30*lomz**2*lz*z**2 + 6*lz**2*z**2 - &
            24*lomz*lz**2*z**2 + 13*lz**3*z**2 + 20*pisq*z**2 - 8*lomz*pisq*z**2 + 16*lz*pisq*z**2 + 24*zeta3 + 168*z**2*zeta3))/(6.*omz) +&
            (Qe*Qq*(6*Geq - 12*ln2*logeta4 - 12*Geq*lomz + 24*ln2*logeta4*lomz - 12*Geq*z + 24*ln2*logeta4*z + 6*Geq*z**2 - &
            12*ln2*logeta4*z**2 + 12*Geq*lomz*z**2 - 24*ln2*logeta4*lomz*z**2 - 12*lz*Log(s13/s14) - 24*lomz*lz*Log(s13/s14) + & 
            24*lz*z*Log(s13/s14) - 12*lz*z**2*Log(s13/s14) - 24*lomz*lz*z**2*Log(s13/s14)))/(6.*omz)

    ! combine regular and plus
    qqb_ga_sub_qqbz12(0) = (pls_Lsq + reg_Lsq) * Lmu_z**2 + (pls_Lmu + reg_Lmu) * Lmu_z +  (pls_noL + reg_noL)    ! plus+reg
    qqb_ga_sub_qqbz12(1) = pls_Lsq*Lmu**2 + pls_Lmu*Lmu + pls_noL                                                 ! plus at z=1



    ! overall Cf
    qqb_ga_sub_qqbz12 = qqb_ga_sub_qqbz12 * Cf  

  end function qqb_ga_sub_qqbz12



  function qqb_ga_sub_qqb12(mu,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,s13,s14,s23,s24,Qq,Qe)
    ! Integrated subtraction terms to multiply FLM[1,2], cf. Eq. 264
    implicit none
    real(dp), intent(in)  :: mu,Ec,E1,E2,E3,E4,eta31,eta32,eta41,eta42,eta34,s13,s14,s23,s24,Qq,Qe
    real(dp)              :: qqb_ga_sub_qqb12
    real(dp)              :: Lmu

    Lmu  = log(two*Ec/mu)     !  note this has mu not musq   !
 
    qqb_ga_sub_qqb12 = Qq**2*(16.0_dp/45.0_dp*pisq + 8.0_dp*(pisq+8.0_dp*zeta3)*Lmu + four/three*(27.0_dp - four*pisq)*Lmu**2)
    qqb_ga_sub_qqb12 = qqb_ga_sub_qqb12 + Qe*Qq*fnGeq(mu,E3,E4,eta31,eta32,eta41,eta42,s13,s14,s23,s24)*(two/three*pisq + 6.0_dp*Lmu)
    qqb_ga_sub_qqb12 = qqb_ga_sub_qqb12 + Qe**2*fnGesq(Ec,E3,E4,eta34)*(two/three*pisq + 6.0_dp*Lmu)

    qqb_ga_sub_qqb12 = Cf*qqb_ga_sub_qqb12                ! overall factor of Cf
    

  end function qqb_ga_sub_qqb12


  function fnGeq(mu,E3,E4,eta31,eta32,eta41,eta42,s13,s14,s23,s24)
    ! Eq 267 in Chiara's notes
    real(dp), intent(in) :: mu,E3,E4,eta31,eta32,eta41,eta42,s13,s14,s23,s24
    real(dp)             :: fnGeq
    
    fnGeq = three*log(eta31*eta42/eta41/eta32) + &
         two*dilog2(one-eta31) - two*dilog2(one-eta41) - two*dilog2(one-eta32) + two*dilog2(one-eta42) + &
         two*log(E3/mu)*log(s23/s13) + two*log(E4/mu)*log(s14/s24)

  end function FnGeq

  function fnGesq(Ec,E3,E4,eta34)
    ! Eq 267 in Chiara's notes
    real(dp), intent(in) :: Ec,E3,E4,eta34
    real(dp)             :: fnGesq
    
    fnGesq = 13.0_dp - two/three*pisq + log(E3/E4)**2 + log(eta34)*(three - two*log(E4*E4/Ec**2)) + two*dilog2(one-eta34)
  
  end function FnGesq


end module mod_subtrfn_nnlo_z_qqb_ga_raoul

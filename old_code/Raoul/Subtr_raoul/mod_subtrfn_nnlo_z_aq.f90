module mod_subtrfn_nnlo_z_aq
! integrated subtraction terms for a qb->qb g  channel
! cf. Notes from Chiara and Federica
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_auxfunctions
  use mod_splittings_bare
  implicit none

  private 

  public :: aqb_qbg_sub_qqbg_zi
  public :: aqb_qbg_sub_aqbq_zi, aqb_qbg_sub_aqbq_zi_pls, aqb_qbg_sub_aqbq
  public :: aqb_qbg_sub_qqbz12,aqb_qbg_sub_qqbz1zbar2,aqb_qbg_sub_aa1z2

  
contains

 
  function aqb_qbg_sub_qqbg_zi(Ec,mu,z,eta51,eta51log,Qq)
    ! Integrated subtraction terms to multiply FLMqqbg[z.1,2,3,4|5g ], cf. Eq. 482
    real(dp), intent(in)   :: Ec,mu,z,eta51,eta51log,Qq
    real(dp)               :: aqb_qbg_sub_qqbg_zi
    real(dp)               :: Lmu,omz,lomz,wt61,Pqg(-1:1),eta52

    Lmu = log(mu/two/Ec)
    omz = one-z
    lomz = log(omz)

    eta52 = one-eta51             ! back-to-back incoming partons
    wt61 = eta52*(one+eta51)      ! coll. limit of partition function
    Pqg = PqgAP(z)

    
! scale dependent pieces
    aqb_qbg_sub_qqbg_zi = -two*Pqg(0)*Lmu

! scale independent pieces
    aqb_qbg_sub_qqbg_zi = aqb_qbg_sub_qqbg_zi + &
         two*(omz*z + lomz*Pqg(0)) +&
         log(eta51log/two)*wt61*Pqg(0)


    aqb_qbg_sub_qqbg_zi = aqb_qbg_sub_qqbg_zi * xn * Qq**2
    


  end function aqb_qbg_sub_qqbg_zi



  function aqb_qbg_sub_aqbq_zi(Ec,mu,z,eta52,eta52log)
    ! Integrated subtraction terms to multiply FLMaqbq[1,z2,3,4|5q ], cf. Eq. 482
    ! pls+reg
    real(dp), intent(in)   :: Ec,mu,z,eta52,eta52log
    real(dp)               :: aqb_qbg_sub_aqbq_zi
    real(dp)               :: Lmu,omz,lomz,wt52

    Lmu = log(mu/two/Ec)
    omz = one-z
    lomz = log(omz)
    wt52 = one-eta52**2      ! coll. limit of partition function

    
! scale dependent pieces
    aqb_qbg_sub_aqbq_zi = two*Cf*(one+z)*Lmu                      ! reg
    aqb_qbg_sub_aqbq_zi = aqb_qbg_sub_aqbq_zi - four*Cf/omz*Lmu   ! pls


! scale independent pieces
    aqb_qbg_sub_aqbq_zi = aqb_qbg_sub_aqbq_zi + &
         Cf*( omz-two*(one+z)*lomz - (one+z)*log(eta52log/two)*wt52 +  & ! reg
         two*log(eta52log/two)*wt52/omz + four*lomz/omz )
    

  end function aqb_qbg_sub_aqbq_zi


  
  function aqb_qbg_sub_aqbq_zi_pls(Ec,mu,z,eta52,eta52log)
    ! Integrated subtraction terms to multiply FLMaqbq[z.i,3,4|5g ], cf. Eq. 482
    ! pls
    real(dp), intent(in)   :: Ec,mu,z,eta52,eta52log
    real(dp)               :: aqb_qbg_sub_aqbq_zi_pls
    real(dp)               :: Lmu,omz,lomz,wt52

    Lmu = log(mu/two/Ec)
    omz = one-z
    lomz = log(omz)
    wt52 = one-eta52**2      ! coll. limit of partition function

    
! scale dependent pieces
    aqb_qbg_sub_aqbq_zi_pls =  -four*Cf/omz*Lmu   ! pls


! scale independent pieces
    aqb_qbg_sub_aqbq_zi_pls = aqb_qbg_sub_aqbq_zi_pls + &
         two*Cf*log(eta52log/two)*wt52/omz + four*Cf*lomz/omz

   
  end function aqb_qbg_sub_aqbq_zi_pls

  function aqb_qbg_sub_aqbq(Ec,mu,E5,eta52,eta51log,eta52log)
    ! Integrated subtraction terms to multiply FLMaqbq[1,2,3,4|5 ], cf. Eq. 482
    ! pls+reg
    real(dp), intent(in)   :: Ec,mu,E5,eta52,eta51log,eta52log
    real(dp)               :: aqb_qbg_sub_aqbq
    real(dp)               :: Lmu,LEcoE5,eta51,wt5161,wt5262

    Lmu = log(mu/Ec)
    LEcoE5 = log(Ec/E5)
    eta51 = one-eta52
    wt5161 = eta52**2*(one+two*eta51)
    wt5262 = eta51**2*(one+two*eta52)

! scale dependent pieces
    aqb_qbg_sub_aqbq = -three*Cf*Lmu
    aqb_qbg_sub_aqbq = aqb_qbg_sub_aqbq + Cf*( &
         13.0_dp/two - pisq + LEcoE5**2 -&
         (three/two+two*LEcoE5)*log(eta51log/eta52log)*wt5161 - &
         (three/two+two*LEcoE5)*log(eta52log/eta51log)*wt5262 + &
         (three+two*LEcoE5)*log(four*eta52log)+2*real(dilog2(one-eta52),kind=dp) )
    
    
  end function aqb_qbg_sub_aqbq

  function aqb_qbg_sub_qqbz1zbar2(mu,Ec,z1,z2,Qq)
    !! Integrated subtraction terms to multiply FLM[z1,zbar2] cf. Eq. 481
    !! log(musq/s) evaluated for {Ec,Ec/z2}
    !! output indices: 0=reg+pls, 1=pls
    implicit none
    real(dp), intent(in)  :: mu,Ec,z1,z2,Qq
    real(dp)              :: aqb_qbg_sub_qqbz1zbar2(0:1)
    real(dp)              :: pls0_Lsq,pls_Lsq,reg_Lsq
    real(dp)              :: pls1_Lmu,pls0_Lmu,pls_Lmu,reg_Lmu
    real(dp)              :: pls1_noL,reg_noL,pls_noL
    real(dp)              :: omz1,omz2,lomz1,lomz2,Lmu_z1,Lmu_z1z2


    omz1 = one-z1
    omz2 = one-z2
    lomz1 = log(one-z1)
    lomz2 = log(one-z2)

    Lmu_z1   = log(mu/two/Ec) + log(z1)/two
    Lmu_z1z2 = log(mu/two/Ec) + log(z1)/two + log(z2)/two

    ! Lmu^2 terms
    pls0_Lsq = 8.0_dp*(one-two*z1+two*z1**2)
    reg_Lsq = -four*(one+z2)*(one-two*z1+two*z1**2)

    pls_Lsq = pls0_Lsq / omz2

    ! Lmu terms
    pls1_Lmu = -8.0_dp*(one-two*z1+two*z1**2)
    pls0_Lmu = -8.0_dp*((one - z1)*z1 + (one - two*z1 + two*z1**2)*lomz1)
    reg_Lmu  = +two*(-one + four*z1 - four*z1**2 + z2 + &
         two*(one - two*z1 + two*z1**2)*(one + z2)*lomz1 + &
         two*(one - two*z1 + two*z1**2)*(one + z2)*lomz2)

    pls_Lmu = pls1_Lmu*lomz2/omz2 + pls0_Lmu/omz2

    ! mu-independent terms
    pls1_noL = 8.0_dp*(-((-one + z1)*z1) + (one - two*z1 + two*z1**2)*lomz1)
    reg_noL = -two*(-((-one + z1)*z1) + &
         (one - two*z1 + two*z1**2)*lomz1)*(-one + z2 + two*(one + z2)*lomz2)

    pls_noL = pls1_noL*lomz2/omz2

    ! combine regular and plus
    aqb_qbg_sub_qqbz1zbar2(0) = (pls_Lsq + reg_Lsq) * Lmu_z1z2**2 + (pls_Lmu + reg_Lmu) * Lmu_z1z2 +  (pls_noL + reg_noL)    ! plus+reg
    aqb_qbg_sub_qqbz1zbar2(1) = pls_Lsq*Lmu_z1**2 + pls_Lmu*Lmu_z1 + pls_noL                                                   ! plus at z2=1

    ! Overall factor Cf * Nc * Qq^2
    aqb_qbg_sub_qqbz1zbar2 = Cf*xn*Qq**2*aqb_qbg_sub_qqbz1zbar2

  end function aqb_qbg_sub_qqbz1zbar2

    
    


  function aqb_qbg_sub_qqbz12(mu,Ec,z,Qq)
!! Integrated subtraction terms to multiply FLM[z1_q,2_q] cf. Eq. 481
    implicit none
    real(dp), intent(in)  :: mu,Ec,z,Qq
    real(dp)              :: aqb_qbg_sub_qqbz12
    real(dp)              :: omz,lz,lomz,li2z,li2omz,li3z,li3omz,Lmu

    omz = one-z
    lomz = log(omz)
    lz = log(z)
    li2omz = real(dilog2(omz),kind=dp)
    li2z = real(dilog2(z),kind=dp)
    li3omz = trilog(omz)
    li3z   =  trilog(+z)

    Lmu = log(mu/two/Ec) + log(z)/two

    ! Lmu^2
    aqb_qbg_sub_qqbz12 = (5.0_dp + four*z*(-two + three*z) + lomz*(four + 8.0_dp*(-one + z)*z) + lz*(-two + four*z - 8.0_dp*z**2))*Lmu**2

    ! Lmu 
    aqb_qbg_sub_qqbz12 = aqb_qbg_sub_qqbz12 + &
         Lmu*( -16.0_dp - lz**2*(one - two*z) + 8.0_dp*li2z*(-one + z)**2 - z*(-41.0_dp + 30.0_dp*z) - lomz*(-two - 8.0_dp*(-two + z)*z) - li2omz*(-four - 8.0_dp*(-one + z)*z) - &
         lomz**2*(6.0_dp + 12.0_dp*(-one + z)*z) - pisq*(two - four*z + (8.0_dp*z**2)/three) - lz*(three + 20.0_dp*(-one + z)*z - 8.0_dp*lomz*(one + two*(-one + z)*z)))

    !mu-independent
    aqb_qbg_sub_qqbz12 = aqb_qbg_sub_qqbz12 - &
         69.0_dp/four + three*li2z + (255.0_dp*z/four - two*lomz*lz*(-one + z)*z - 49.0_dp*z**2 + lomz**2*(-7.0_dp + (21.0_dp - 17.0_dp*z)*z) + &
         (11.0_dp*lomz**3*(one + two*(-one + z)*z))/6.0_dp + lz**3*(5.0_dp/four - (5.0_dp*z)/two + (7.0_dp*z**2)/three) + (pisq*(11.0_dp + two*z*(-13.0_dp + 9.0_dp*z)))/12.0_dp + &
         lz**2*(-three/8.0_dp - z/two + z**2 + lomz*(-5.0_dp/two - 5.0_dp*(-one + z)*z)) + &
         lz*(one + (35.0_dp*z)/four - 8.0_dp*z**2 + lomz**2*(-11.0_dp/two + (11.0_dp - 7.0_dp*z)*z) + li2z*(one + two*(-one + z)*z) + (four*pisq*(one + two*(-one + z)*z))/three + &
         lomz*(6.0_dp + 16.0_dp*(-one + z)*z)) + lomz*(16.0_dp - 8.0_dp*li2z*(-one + z)**2 - (109.0_dp*z)/two + 44.0_dp*z**2 + li2omz*(-three - 6.0_dp*(-one + z)*z) + &
         (pisq*(5.0_dp + two*z*(-5.0_dp + three*z)))/three) + (18.0_dp + 4.0_dp*z*(-9.0_dp + 8.0_dp*z))*zeta3 + (-9.0_dp + two*(9.0_dp - 5.0_dp*z)*z)*li3omz + &
         (-9.0_dp + two*(9.0_dp - 7.0_dp*z)*z)*li3z)

    ! overall factor Cf*Nc*Qq^2
    aqb_qbg_sub_qqbz12 = Cf*xn*Qq**2 * aqb_qbg_sub_qqbz12 

  end function aqb_qbg_sub_qqbz12


  function aqb_qbg_sub_aa1z2(mu,Ec,z,Qq)
!! Integrated subtraction terms to multiply FLM[1_a,z2_a] cf. Eq. 481
    implicit none
    real(dp), intent(in)  :: mu,Ec,z,Qq
    real(dp)              :: aqb_qbg_sub_aa1z2
    real(dp)              :: omz,lz,lomz,li2z,li3z,li3omz,Lmu

    omz = one-z
    lomz = log(omz)
    lz = log(z)
    li2z = real(dilog2(z),kind=dp)
    li3omz = trilog(omz)
    li3z   =  trilog(+z)

    Lmu = log(mu/two/Ec) + log(z)/two

    ! Lmu^2
    aqb_qbg_sub_aa1z2 = (four - two*lz*(-two + z) - z + lomz*(-8.0_dp + 8.0_dp/z + four*z))*Lmu**2

    ! Lmu
    aqb_qbg_sub_aa1z2 = aqb_qbg_sub_aa1z2 + &
         (15.0_dp - 8.0_dp*lz - two*lz**2 - four*pisq - 10.0_dp/z + (8.0_dp*pisq)/(three*z) + four*z -&
         5.0_dp*lz*z + lz**2*z + two*pisq*z + lomz*(-20.0_dp + 12.0_dp/z + 8.0_dp*z) -&
         (6.0_dp*lomz**2*(two - two*z + z**2))/z - four*(-two + z)*li2z)*Lmu

    ! mu-independent
    aqb_qbg_sub_aa1z2 = aqb_qbg_sub_aa1z2 + &
      z*(lz*(5.0_dp/four + lomz**2*(two - four/z) - three/(four*z)) + lz**2*(17.0_dp/8.0_dp + 5.0_dp/(two*z)) - &
      (lz**3*(-two + z))/(12.0_dp*z) + (lomz**2*(-42.0_dp + (58.0_dp - 27.0_dp*z)*z))/(four*z**2) + &
      (13.0_dp*lomz**3*(two + (-two + z)*z))/(6.0_dp*z**2) + (lomz*(108.0_dp + 24.0_dp*li2z*(-two + z)*z -&
      three*z*(46.0_dp + three*z) - two*pisq*(10.0_dp + 7.0_dp*(-two + z)*z)))/(6.0_dp*z**2) + &
      (-324.0_dp + 438.0_dp*z - 87.0_dp*z**2 - 12.0_dp*li2z*z*(8.0_dp + 5.0_dp*z) +&
      pisq*(30.0_dp + z*(-22.0_dp + 27.0_dp*z)) + 192.0_dp*zeta3 + &
      24.0_dp*(-two + z)*z*(two*li3omz + li3z + three*zeta3))/(12.0_dp*z**2))

    ! overall factor Cf*Qq^2
    aqb_qbg_sub_aa1z2 = Cf*Qq**2 * aqb_qbg_sub_aa1z2

  end function aqb_qbg_sub_aa1z2



end module mod_subtrfn_nnlo_z_aq

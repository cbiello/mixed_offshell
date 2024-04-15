module mod_subtrfn_nnlo_z_gq
! integrated subtraction terms for a qb->qb g  channel
! cf. Notes from Chiara and Federica
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_auxfunctions
  use mod_splittings_bare
  implicit none

  private 

  public :: gqb_qba_sub_qqba_zi
  public :: gqb_qba_sub_gqbq_zi,gqb_qba_sub_gqbq_zi_pls,gqb_qba_sub_gqbq
  public :: gqb_qba_sub_qqbz12, gqb_qba_sub_qqbz1zbar2,gqb_qba_sub_qqbz1zbar2_pls
  public :: gqb_qba_sub_qqbz12_vec

  
contains

 
  function gqb_qba_sub_qqba_zi(Ec,mu,z,eta51,eta52,eta53,eta54,eta51log)
    ! Integrated subtraction terms to multiply FLMqqba[z.1,2,3,4|5a ], cf. Eq. 358
    real(dp), intent(in)   :: Ec,mu,z,eta51,eta52,eta53,eta54,eta51log
    real(dp)               :: gqb_qba_sub_qqba_zi
    real(dp)               :: Lmu,omz,lomz,wt51,Pqg(-1:1)

    Lmu = log(mu/two/Ec)
    omz = one-z
    lomz = log(omz)

    wt51 = eta52*eta53*eta54/(eta51*eta53*eta54+eta51*eta52*eta54+eta51*eta52*eta53+eta52*eta53*eta54)
    Pqg = PqgAP(z)

    
! scale dependent pieces
    gqb_qba_sub_qqba_zi = -two*Pqg(0)*Lmu

! scale independent pieces
    gqb_qba_sub_qqba_zi = gqb_qba_sub_qqba_zi + &
         two*(omz*z + lomz*Pqg(0)) +&
         log(eta51log/two)*wt51*Pqg(0)


    gqb_qba_sub_qqba_zi = tr*gqb_qba_sub_qqba_zi
    


  end function gqb_qba_sub_qqba_zi


  function gqb_qba_sub_gqbq_zi(Ec,mu,z,Qq)    
    ! Integrated subtraction terms to multiply FLMgqbq[1,z.2,3,4|5q ], cf. Eq. 358
    ! reg + pls
    real(dp), intent(in)   :: Ec,mu,z,Qq
    real(dp)               :: gqb_qba_sub_gqbq_zi
    real(dp)               :: Lmu,omz,lomz

    Lmu = log(mu/two/Ec)
    omz = one-z
    lomz = log(omz)

    ! scale dependent pieces
     gqb_qba_sub_gqbq_zi =  -four*Lmu/omz  +&         ! pls
          two*(one+z)*Lmu                                ! reg

     ! scale independent pieces
         ! scale dependent pieces
     gqb_qba_sub_gqbq_zi = gqb_qba_sub_gqbq_zi + &
          four*lomz/omz +&                           ! pls
          omz -two*(one+z)*lomz                     ! reg
     
     gqb_qba_sub_gqbq_zi = gqb_qba_sub_gqbq_zi * Qq**2
     
    
   end function gqb_qba_sub_gqbq_zi


  function gqb_qba_sub_gqbq_zi_pls(Ec,mu,z,Qq)    
    ! Integrated subtraction terms to multiply FLMgqbq[1,z.2,3,4|5q ], cf. Eq. 358
    ! pls only
    real(dp), intent(in)   :: Ec,mu,z,Qq
    real(dp)               :: gqb_qba_sub_gqbq_zi_pls
    real(dp)               :: Lmu,omz,lomz

    Lmu = log(mu/two/Ec)
    omz = one-z
    lomz = log(omz)

    ! scale dependent pieces
     gqb_qba_sub_gqbq_zi_pls =  -four*Lmu/omz 


     ! scale dependent pieces
     gqb_qba_sub_gqbq_zi_pls = gqb_qba_sub_gqbq_zi_pls + &
          four*lomz/omz 
     
     gqb_qba_sub_gqbq_zi_pls = gqb_qba_sub_gqbq_zi_pls * Qq**2
     
    
   end function gqb_qba_sub_gqbq_zi_pls




  function gqb_qba_sub_gqbq(Ec,mu,E2,E3,E4,E5,eta32,eta34,eta42,eta51,eta52,eta53,eta45,eta51log,Qq,Qe)    
    ! Integrated subtraction terms to multiply FLMgqbq[1,2,3,4|5q ], cf. Eq. 358
    real(dp), intent(in)   :: Ec,mu,E2,E3,E4,E5,eta32,eta34,eta42,eta51,eta51log,eta52,eta53,eta45,Qq,Qe
    real(dp)               :: gqb_qba_sub_gqbq,s23,s24,s35,s45
    real(dp)               :: Lmu,GQeQq,GQesq,GQqsq

    Lmu = log(mu/two/Ec)

    ! sij
    s23 = four*E2*E3*eta32
    s24 = four*E2*E4*eta42
    s35 = four*E3*E5*eta53
    s45 = four*E4*E5*eta45


!    GQeQq = fnGQeQq(mu,Ec,E3,E4,E5,eta32,eta42,eta45,eta53,s23,s24,s35,s45)
    GQeQq = fnGQeQq(Ec,Ec,E3,E4,E5,eta32,eta42,eta45,eta53,s23,s24,s35,s45)
    GQesq = fnGQesq(Ec,E3,E4,eta34)
    
    GQqsq = fnGQqsq(mu,E5,Ec,eta52)

    !print *, GQqsq-(three/two + two*log(Ec/E5))*log(eta51log/(one-eta51))
    !print *, GqeQq
    !print *, GQesq

    ! scale dependent pieces
     gqb_qba_sub_gqbq = Qe**2*GQesq + Qq**2*GQqsq + Qe*Qq*GQeQq - &
          Qq**2*(three/two + two*log(Ec/E5))*log(eta51log/(one-eta51))
     
          
   end function gqb_qba_sub_gqbq


   function gqb_qba_sub_qqbz1zbar2(mu,z,zb,Ec,Qq)
     ! Integrated subtraction terms to multiply FLM[z1,zbar2], cf. Eq. 357                                                                             ! Gives pls+reg, i.e. coefficient of FLM[z1,zbar2]
     implicit none
     real(dp), intent(in) :: mu,z,zb,Ec,Qq
     real(dp)             :: gqb_qba_sub_qqbz1zbar2
     real(dp)             :: Lmu,omzb,lomzb,lomz
     
     Lmu  = log(mu/two/Ec)! + log(z) + log(zb)

     omzb = one-zb
     lomzb=log(omzb)
     lomz = log(one-z)

     ! Lmu^2 terms
     gqb_qba_sub_qqbz1zbar2 = (four*(one - two*z  + two*z **2)*(one + zb**2))/(omzb)*Lmu**2

     !Lmu terms
     gqb_qba_sub_qqbz1zbar2 = gqb_qba_sub_qqbz1zbar2 + &
          Lmu*(-two*(one - two*zb + four*z*zb - four*z**2*zb + zb**2 + two*lomz*(one - two*z + two*z**2)*(one + zb**2) +  &
          two*lomzb*(one - two*z + two*z**2)*(one + zb**2)))/(omzb)
     
     ! scale independent terms
     gqb_qba_sub_qqbz1zbar2 = gqb_qba_sub_qqbz1zbar2 + &
          (two*(-((-one + z)*z) + lomz*(one - two*z + two*z**2))*((-one + zb)**2 + two*lomzb*(one + zb**2)))/(omzb)

     ! overall factors
     gqb_qba_sub_qqbz1zbar2 = TR * Qq**2 * gqb_qba_sub_qqbz1zbar2
          

   end function gqb_qba_sub_qqbz1zbar2


   function gqb_qba_sub_qqbz1zbar2_pls(mu,z,zb,Ec,Qq)
     ! Integrated subtraction terms to multiply FLM[z1,zbar2], cf. Eq. 357                                                                             ! gives -pls piece, i.e. -FLM[z1,2]
     implicit none
     real(dp), intent(in) :: mu,z,zb,Ec,Qq
     real(dp)             :: gqb_qba_sub_qqbz1zbar2_pls
     real(dp)             :: Lmu,omzb,lomzb,lomz
     
     Lmu  = log(mu/two/Ec)! + log(z)

     omzb = one-zb
     lomzb=log(omzb)
     lomz = log(one-z)

     ! Lmu^2 terms
     gqb_qba_sub_qqbz1zbar2_pls = (-8.0_dp*(one - two*z  + two*z**2))/(omzb)*Lmu**2

     !Lmu terms
     gqb_qba_sub_qqbz1zbar2_pls = gqb_qba_sub_qqbz1zbar2_pls + &
          Lmu*(8.0_dp*(-((-one + z)*z) + lomz*(one - two*z + two*z**2) + lomzb*(one - two*z + two*z**2)))/(omzb)

     gqb_qba_sub_qqbz1zbar2_pls = gqb_qba_sub_qqbz1zbar2_pls + &
          (-8.0_dp*lomzb*(-((-one + z)*z) + lomz*(one - two*z + two*z**2)))/(omzb)

     ! overall factors
     gqb_qba_sub_qqbz1zbar2_pls = TR * Qq**2 * gqb_qba_sub_qqbz1zbar2_pls

          

   end function gqb_qba_sub_qqbz1zbar2_pls

   
!   function gqb_qba_sub_qqbz12(mu,z,Ec,E3,E4,eta13,eta23,eta14,eta24,eta34,eta34b,Qq,Qe)
   function gqb_qba_sub_qqbz12(mu,z,Ec,E3,E4,eta13,eta23,eta14,eta24,eta34,ometa34,Qq,Qe)
     ! Integrated subtraction terms to multiply FLM[z1,2], cf. Eq. 357
     implicit none
     real(dp), intent(in) :: mu,z,Ec,E3,E4,eta13,eta23,eta14,eta24,eta34,ometa34,Qq,Qe
     real(dp)             :: gqb_qba_sub_qqbz12
     real(dp)             :: omz,lz,lomz,li2omz,li2z,li3omz,li3z,Lmu_z
     real(dp)             :: li2eta13,li2eta14,li2eta23,li2eta24,li2eta34


     omz    = one-z
     lz   = log(z)
     lomz   = log(omz)
     li2omz = real(dilog2(omz),kind=dp)
     li2z = real(dilog2(z),kind=dp)
     li3omz = trilog(omz)
     li3z   =  trilog(+z)

     li2eta13 = real(dilog2(eta13),kind=dp)
     li2eta23 = real(dilog2(eta23),kind=dp)
     li2eta14 = real(dilog2(eta14),kind=dp)
     li2eta24 = real(dilog2(eta24),kind=dp)
     li2eta34 = real(dilog2(eta34),kind=dp)

     Lmu_z = log(mu/two/Ec)! + log(z)

     ! Lmu^2
     gqb_qba_sub_qqbz12 = (Qq**2*(5.0_dp - 8.0_dp*z + 12.0_dp*z**2 + lz*(-two + four*z - 8.0_dp*z**2) + lomz*(four - 8.0_dp*z + 8.0_dp*z**2)))/z * Lmu_z**2

     !Lmu
     gqb_qba_sub_qqbz12 = gqb_qba_sub_qqbz12 - &
          ! Qq^2
          Qq**2 * Lmu_z * ((16.0_dp + (four*pisq*(-1 + z)**2)/three + four*li2z*(-one + two*z) - lz**2*(-one + 2*z) + z*(-41.0_dp + 30.0_dp*z) + lomz**2*(6.0_dp + 12.0_dp*(-one + z)*z) -&
          two*lomz*(one - 8.0_dp*z + four*z**2) + lz*(three + 20.0_dp*(-one + z)*z + lomz*(-four - 8.0_dp*(-one + z)*z)))/z) + &
          ! Qq*Qe
          Qq*Qe*Lmu_z*(-four*(one - two*z + two*z**2)*(-li2eta13 + li2eta14 + li2eta23 - li2eta24 - lz*Log(E3/E4) - (lz + Log(one - eta13))*Log(eta13) + &
          (lz + Log(one - eta14))*Log(eta14) + (three/two + Log(EC/E3))*Log(eta13/eta23) + Log(one - eta23)*Log(eta23) - &
          Log(one - eta24)*Log(eta24) + (three/two + Log(EC/E4))*Log(eta24/eta14)))/z + &
          !Qe^2
          Qe**2*Lmu_z*( (two*(one - two*z + two*z**2)*(-39.0_dp +  pisq + 6.0_dp*li2eta34 - &
          three*Log(E3/E4)**2 - 9.0_dp*Log(eta34) + 6.0_dp*Log(E3/EC)*Log(eta34) + &
!          6*Log(E4/EC)*Log(eta34) + 6*Log(eta34b)*Log(eta34)))/(3.*z))
          6.0_dp*Log(E4/EC)*Log(eta34) + 6.0_dp*Log(ometa34)*Log(eta34)))/(three*z))

     !scale independent
     gqb_qba_sub_qqbz12= gqb_qba_sub_qqbz12 + &
          ! Qq^2
          Qq**2 * ( 255.0_dp/four + 18.0_dp*li3omz + 18.0_dp*li3z - (109.0_dp*lomz)/two + 10.0_dp*li2z*lomz + 21.0_dp*lomz**2 - (11.0_dp*lomz**3)/three + (35.0_dp*lz)/four - two*li2z*lz -&
          14.0_dp*lomz*lz + 5.0_dp*lomz**2*lz - lz**2/two + 5.0_dp*lomz*lz**2 - (5.0_dp*lz**3)/two - (13.0_dp*pisq)/6.0_dp - (7.0_dp*lomz*pisq)/three - (8.0_dp*lz*pisq)/three - &
          69.0_dp/(four*z) + (three*li2z)/z - (9.0_dp*li3omz)/z - (9.0_dp*li3z)/z + (16.0_dp*lomz)/z - (5.0_dp*li2z*lomz)/z - (7.0_dp*lomz**2)/z + (11.0_dp*lomz**3)/(6.0_dp*z) + &
          lz/z + (li2z*lz)/z + (6.0_dp*lomz*lz)/z - (5.0_dp*lomz**2*lz)/(two*z) - (three*lz**2)/(8.0_dp*z) - (5.0_dp*lomz*lz**2)/(two*z) + (5.0_dp*lz**3)/(four*z) + &
          (11.0_dp*pisq)/(12.0_dp*z) + (7.0_dp*lomz*pisq)/(6.0_dp*z) + (four*lz*pisq)/(three*z) - 49.0_dp*z - 10.0_dp*li3omz*z - 14.0_dp*li3z*z + 44.0_dp*lomz*z - &
          two*li2z*lomz*z - 17.0_dp*lomz**2*z + (11.0_dp*lomz**3*z)/three - 8.0_dp*lz*z + two*li2z*lz*z + 14.0_dp*lomz*lz*z - lomz**2*lz*z + lz**2*z - &
          5.0_dp*lomz*lz**2*z + (7.0_dp*lz**3*z)/three + (three*pisq*z)/two + lomz*pisq*z + (8.0_dp*lz*pisq*z)/three - 36.0_dp*zeta3 + (18.0_dp*zeta3)/z + 32.0_dp*z*zeta3 ) + &
          ! Qq*Qe
          Qq*Qe * ((four*(lomz + z - two*lomz*z - z**2 + two*lomz*z**2)* &
          (-li2eta13 + li2eta14 + li2eta23 - li2eta24 - lz*Log(E3/E4) - (lz + Log(one - eta13))*Log(eta13) + &
          (lz + Log(one - eta14))*Log(eta14) + (three/two + Log(EC/E3))*Log(eta13/eta23) + Log(one - eta23)*Log(eta23) - &
           Log(one - eta24)*Log(eta24) + (three/two + Log(EC/E4))*Log(eta24/eta14)))/z ) + &
          ! Qe**2
           Qe**2*( (-two*(lomz + z - two*lomz*z - z**2 + two*lomz*z**2)* &
           (-39.0_dp + pisq + 6.0_dp*li2eta34 - three*Log(E3/E4)**2 - 9.0_dp*Log(eta34) + 6.0_dp*Log(E3/EC)*Log(eta34) + 6.0_dp*Log(E4/EC)*Log(eta34) + &
!           6*Log(eta34b)*Log(eta34)))/(3.*z))
           6.0_dp*Log(ometa34)*Log(eta34)))/(three*z))

     ! the above formulas are the coefficients of FLM[z1,2], we want coeffs of FLM[z1,2]/z
     gqb_qba_sub_qqbz12 = z * gqb_qba_sub_qqbz12

     ! overall color factor
     gqb_qba_sub_qqbz12 = Tr * gqb_qba_sub_qqbz12


   end function gqb_qba_sub_qqbz12

   function gqb_qba_sub_qqbz12_vec(mu,z,Ec,E3,E4,eta13,eta23,eta14,eta24,eta34,eta34b,Qq,Qe)
     ! Integrated subtraction terms to multiply FLM[z1,2], cf. Eq. 357
     implicit none
     real(dp), intent(in) :: mu,z,Ec,E3,E4,eta13,eta23,eta14,eta24,eta34,eta34b,Qq(:),Qe
     real(dp)             :: gqb_qba_sub_qqbz12_vec(4)
     real(dp)             :: omz,lz,lomz,li2omz,li2z,li3omz,li3z,Lmu_z
     real(dp)             :: li2eta13,li2eta14,li2eta23,li2eta24,li2eta34
     real(dp)             :: a(7), res(4)
     integer :: i


     omz    = one-z
     lz   = log(z)
     lomz   = log(omz)
     li2omz = real(dilog2(omz),kind=dp)
     li2z = real(dilog2(z),kind=dp)
     li3omz = trilog(omz)
     li3z   =  trilog(+z)

     li2eta13 = real(dilog2(eta13),kind=dp)
     li2eta23 = real(dilog2(eta23),kind=dp)
     li2eta14 = real(dilog2(eta14),kind=dp)
     li2eta24 = real(dilog2(eta24),kind=dp)
     li2eta34 = real(dilog2(eta34),kind=dp)

     Lmu_z = log(mu/two/Ec) + log(z)

     a(1) = (five - 8*z + 12*z**2 + lz*(-two + 4*z - 8*z**2) + lomz*(four - 8*z + 8*z**2))/z
     a(2) = -((16._dp + (4*pisq*omz**2)/3. + 4*li2z*(-one + 2*z) - lz**2*(-one + 2*z) + z*(-41._dp + 30*z) + &
            lomz**2*(six - 12*omz*z) -2*lomz*(one - 8*z + 4*z**2) + lz*(3 - 20*omz*z + lomz*(-four + 8*omz*z)))/z)
     a(3) = (-4*(one - 2*z + 2*z**2)*(-li2eta13 + li2eta14 + li2eta23 - li2eta24 - lz*Log(E3/E4) - (lz + Log(one - eta13))*Log(eta13) + &
            (lz + Log(one - eta14))*Log(eta14) + (1.5_dp + Log(EC/E3))*Log(eta13/eta23) + Log(one - eta23)*Log(eta23) - &
            log(one - eta24)*Log(eta24) + (1.5_dp + Log(EC/E4))*Log(eta24/eta14)))/z
     a(4) = ((2*(one - 2*z + 2*z**2)*(-39._dp + pisq + 6*li2eta34 - 3*Log(E3/E4)**2 - 9*Log(eta34) + &
            6*Log(E3/EC)*Log(eta34) + 6*Log(E4/EC)*Log(eta34) + 6*Log(eta34b)*Log(eta34)))/(3*z))
     a(5) = ( 255.0_dp/four + 18*li3omz + 18*li3z - (109*lomz)/2. + 10*li2z*lomz + 21*lomz**2 - (11*lomz**3)/3. + (35*lz)/4. - 2*li2z*lz -&
            14*lomz*lz + 5*lomz**2*lz - lz**2/2. + 5*lomz*lz**2 - (5*lz**3)/2. - (13*pisq)/6. - (7*lomz*pisq)/3. - (8*lz*pisq)/3. - &
            69/(4.*z) + (3*li2z)/z - (9*li3omz)/z - (9*li3z)/z + (16*lomz)/z - (5*li2z*lomz)/z - (7*lomz**2)/z + (11*lomz**3)/(6.*z) + &
            lz/z + (li2z*lz)/z + (6*lomz*lz)/z - (5*lomz**2*lz)/(2.*z) - (3*lz**2)/(8.*z) - (5*lomz*lz**2)/(2.*z) + (5*lz**3)/(4.*z) + &
            (11*pisq)/(12.*z) + (7*lomz*pisq)/(6.*z) + (4*lz*pisq)/(3.*z) - 49*z - 10*li3omz*z - 14*li3z*z + 44*lomz*z - &
            2*li2z*lomz*z - 17*lomz**2*z + (11*lomz**3*z)/3. - 8*lz*z + 2*li2z*lz*z + 14*lomz*lz*z - lomz**2*lz*z + lz**2*z - &
            5*lomz*lz**2*z + (7*lz**3*z)/3. + (3*pisq*z)/2. + lomz*pisq*z + (8*lz*pisq*z)/3. - 36*zeta3 + (18*zeta3)/z + 32*z*zeta3 )
     a(6) = ((4*(lomz + z - 2*lomz*z - z**2 + 2*lomz*z**2)* &
            (-li2eta13 + li2eta14 + li2eta23 - li2eta24 - lz*Log(E3/E4) - (lz + Log(one - eta13))*Log(eta13) + &
            (lz + Log(one - eta14))*Log(eta14) + (1.5_dp + Log(EC/E3))*Log(eta13/eta23) + Log(one - eta23)*Log(eta23) - &
            Log(one - eta24)*Log(eta24) + (1.5_dp + Log(EC/E4))*Log(eta24/eta14)))/z)
     a(7) = ((-2*(lomz + z - 2*lomz*z - z**2 + 2*lomz*z**2)* &
            (-39._dp + pisq + 6*li2eta34 - 3*Log(E3/E4)**2 - 9*Log(eta34) + 6*Log(E3/EC)*Log(eta34) + 6*Log(E4/EC)*Log(eta34) + &
            6*Log(eta34b)*Log(eta34)))/(3.*z))

     res = zero
     do i = 1, size(Qq)
       res(i) = res(i) + (a(1)*Qq(i)**2)*Lmu_z**2                            !-- Lmu^2
       res(i) = res(i) + (a(2)*Qq(i)**2 + a(3)*Qq(i)*Qe + a(4)*Qe**2)*Lmu_z  !-- Lmu^1
       res(i) = res(i) + (a(5)*Qq(i)**2 + a(6)*Qq(i)*Qe + a(7)*Qe**2)        !-- Lmu^0
     end do

     gqb_qba_sub_qqbz12_vec = z * Tr * res

   end function gqb_qba_sub_qqbz12_vec

  function fnGQeQq(mu,Ec,E3,E4,E5,eta32,eta42,eta45,eta53,s23,s24,s35,s45)
    ! Eq 363 in Chiara's notes
    real(dp), intent(in) :: mu,Ec,E3,E4,E5,eta32,eta42,eta45,eta53,s23,s24,s35,s45
    real(dp)             :: fnGQeQq

    fnGQeQq = three*log(eta53*eta42/eta45/eta32) + &
         two*real(dilog2(one-eta53),kind=dp) - two*real(dilog2(one-eta45),kind=dp) -&
         two*real(dilog2(one-eta32),kind=dp) + two*real(dilog2(one-eta42),kind=dp) + &
         two*log(E3/mu)*log(s23/s35) + two*log(E4/mu)*log(s45/s24) - two*log(E5/Ec)*log(eta53/eta45)
    
  end function FnGQeQq



  function fnGQesq(Ec,E3,E4,eta34)
    ! Eq 267 in Chiara's notes
    real(dp), intent(in) :: Ec,E3,E4,eta34
    real(dp)             :: fnGQesq

    fnGQesq = 13.0_dp - two/three*pisq + log(E3/E4)**2 + log(eta34)*(three - two*log(E3*E4/Ec**2)) + two*real(dilog2(one-eta34),kind=dp)

  end function fnGQesq

  function fnGQqsq(mu,E5,Ec,eta52)
    ! Eq 363 in Chiara's note
    real(dp), intent(in)      :: mu,Ec,E5,eta52
    real(dp)                  :: LE5oEc
    real(dp)                  :: fnGQqsq
    
    LE5oEc = log(E5/Ec)
    
    fnGQqsq = 13.0_dp/two - pisq + two*real(dilog2(one-eta52),kind=dp) + LE5oEc**2 - four*ln2*LE5oEc + &
         three*log(four*Ec/mu) + log(eta52)*(three-two*LE5oEc)

  end function fnGQqsq



end module mod_subtrfn_nnlo_z_gq

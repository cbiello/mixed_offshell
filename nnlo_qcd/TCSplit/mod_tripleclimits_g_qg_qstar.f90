module mod_tripleclimits_g_qg_qstar
  use mod_types; use mod_consts_dp; use common_def
  use mod_kinematics_tc_unord
  use mod_auxfunctions

  implicit none
  private
  
  public :: PTC_A_g_qg_qstar  , PTC_B_g_qg_qstar , PTC_C_g_qg_qstar , PTC_D_g_qg_qstar
  public :: tripcoll_subtr_g_qg_qstar


contains


  function tripcoll_subtr_g_qg_qstar(ndim,yRnd,ncomp,ff,userdata,weight,iternumber)
! output indices: 1st is color {1=ab,2=nab,3=tot},2nd is ep power
    implicit none
    !--
    integer, parameter :: ep = 1
    integer, parameter :: mycol = 3 !-- ab -> 1, nab -> 2, tot -> 3
    integer, parameter :: mycontr = 3 !-- 1 -> delta, 2 -> plus, 3 -> regular
    real(dp), parameter :: myz = 0.123_dp
    real(dp), parameter :: zmin = zero
    !--
    real(dp), parameter :: buff = 1E-10_dp
    real(dp), parameter :: onet = one - two * buff
    real(dp), parameter :: bufftc = 1E-10_dp
    ! -- 
    real(dp15) :: yRnd(30), ff(1), weight
    integer :: ndim, ncomp, userdata, iternumber,tripcoll_subtr_g_qg_qstar
    real(dp)              :: x1,x2,x3,xx2,xx(1:3),x3a,xx3
    real(dp)              :: res_tc(1:3,1:3,0:1)
    real(dp)              :: seca(1:3,1:3,0:1),secb(1:3,1:3,0:1)
    real(dp)              :: secc(1:3,1:3,0:1),secd(1:3,1:3,0:1)
    real(dp)              :: seca0(1:3,1:3,0:1),secb0(1:3,1:3,0:1)
    real(dp)              :: secc0(1:3,1:3,0:1),secd0(1:3,1:3,0:1)
    real(dp)              :: wghtac,wghtbd,z

    tripcoll_subtr_g_qg_qstar = zero
    z = myz
    ff(1) = zero
    res_tc = zero 

    x1=buff+onet*real(yRnd(1),dp)  ! r
    x2=buff+onet*real(yRnd(2),dp)  ! x4
    x3a=buff+onet*real(yRnd(3),dp)            ! lambda
    xx3 = buff+onet*real(yRnd(4),dp)  ! x3

! remap lambda
    x3 = sin( half * pi * x3a )**2
    ! DEBUG -- no buffer dependence yet

    if (x1*x2*(one-z) .lt. bufftc) return


    xx = (/x1,x2,x3/)


    seca0 = PTC_A_g_qg_qstar(z,zmin,xx,zero)
    secc0 = PTC_C_g_qg_qstar(z,zmin,xx,zero)
    if (withmax) then
       seca = PTC_A_g_qg_qstar(z,zmin,xx,xx3)
       secc = PTC_C_g_qg_qstar(z,zmin,xx,xx3)
    endif
       

    wghtac = -one/x1/x2

    xx2 = x2**2
    xx = (/x1,xx2,x3/)

    secb0 = PTC_B_g_qg_qstar(z,zmin,xx,zero)
    secd0 = PTC_D_g_qg_qstar(z,zmin,xx,zero)

    if (withmax) then
       secb = PTC_B_g_qg_qstar(z,zmin,xx,xx3)
       secd = PTC_D_g_qg_qstar(z,zmin,xx,xx3)
    endif

    wghtbd = -one/x1/x2*two

    if (xx2 .lt. bufftc) wghtbd=zero



    ! no buffer dependence yet

    res_tc = (seca0 + secc0 )*wghtac

    res_tc = res_tc + (secb0 + secd0) *wghtbd

    if (withmax) then 
       res_tc(:,:,1) = res_tc(:,:,1) - two*( &
            seca(:,:,0)-seca0(:,:,0) + &
            secc(:,:,0)-secc0(:,:,0))/xx3*wghtac
       
       res_tc(:,:,1) = res_tc(:,:,1) - two*( &
            secb(:,:,0)-secb0(:,:,0) + &
            secd(:,:,0)-secd0(:,:,0))/xx3*wghtbd
    endif

    ! remove color factors from Abelian and non-Abelian parts
    res_tc(1,1:3,0:1) = res_tc(1,1:3,0:1)/(Cf**2)
    res_tc(2,1:3,0:1) = res_tc(2,1:3,0:1)/(Ca*Cf)
    res_tc(3,1:3,0:1) = -res_tc(3,1:3,0:1)*tr/Cf

    ff(1) = res_tc(mycol,mycontr,ep)

    if (ff(1) .ne. ff(1) .or. weight .ne. weight) then
       print *, "xx",xx,x3a
       print *, "wghtbd",wghtbd
       print *, "sec a",seca(mycol,mycontr,ep)
       print *, "sec b",secb(mycol,mycontr,ep)
       print *, "sec c",secc(mycol,mycontr,ep)
       print *, "sec d",secd(mycol,mycontr,ep)
       print *, "ff",ff(1)
       stop
    endif


  end function tripcoll_subtr_g_qg_qstar




!    triple collinear splitting function  Eqs. in sect. C.5

   function TC(kinlist) 
   real(dp) :: TC(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3 

   s14 = kinlist(1) 
   s15 = -kinlist(3) 
   s45 = -kinlist(2) 
   z4 = kinlist(4) 
   z5 = kinlist(5)

   ! convert from the standalone kinematics output to the en. unord. output
   z4 = one-z4-z5

   TC = zero 

   s13 = -s14
   s23 = -s15 
   s12 = s45
   s123 = s12 + s13 + s23 

   z1 = z4 
   z2 = z5 
   z3 = one - z4 - z5

   t123 = two*(z1*s23 - z2*s13)/(z1+z2) + (z1-z2)/(z1+z2)*s12

! -- abelian 


   TC(1,0) = Cf**2*(  s123**2/two/s13/s23*z3*( (one+z3**2)/z1/z2  ) & 
             + s123/s13*(  ( z3*(one-z1) + (one-z2)**3)/z1/z2 ) - s23/s13) 


   TC(1,1) = Cf**2*(  s123**2/two/s13/s23*z3*( -(z1**2+z2**2)/z1/z2 - one  ) & 
             + s123/s13*( -one*(z1**2+z1*z2+z2**2)*(one-z2)/z1/z2 ) +one + two*s23/s13 )



!-- non-abelian 

   TC(2,0)  = Cf*Ca*(   0.25_dp*t123**2/s12**2 + 0.25_dp + s123**2/two/s12/s13*( ( (one-z3)**2 + two*z3 )/z2 & 
     + ( z2**2 + two*(one-z2) )/(one-z3) ) - s123**2/4.0_dp/s13/s23*z3*( ( ( one-z3)**2+two*z3 )/z1/z2  ) & 
  + s123/two/s12*(  ( z1*(two -two*z1 + z1**2) - z2*(6.0_dp-6.0_dp*z2+z2**2 ) )/z2/(one-z3)  )  & 
  + s123/two/s13*( ( (one-z2)**3 + z3**2 - z2)/z2/(one-z3)  & 
  - ( z3*(one-z1)+(one-z2)**3)/z1/z2  ) )


   TC(2,1)  = Cf*Ca*(  -one*( 0.25_dp*t123**2/s12**2 + 0.25_dp) - half  &
     + s123**2/two/s12/s13*( ( -(one-z3)**2 )/z2  + ( -z2**2 )/(one-z3) ) &
   - s123**2/4.0_dp/s13/s23*z3*(    ( -( one-z3)**2)/z1/z2 + one ) & 
  + s123/two/s12*(  -( z1*(two -two*z1 + z1**2) - z2*(6.0_dp-6.0_dp*z2+z2**2 ) )/z2/(one-z3) &
       + two*(z3*(z1-two*z2)-z2)/z2/(one-z3) )  & 
  + s123/two/s13*(    -one*( (one-z2)**3 + z3**2 - z2)/z2/(one-z3)  & 
          -(  two*(one-z2)*(z2-z3)/z2/(one-z3) - z1 + z2 ) + (one-z2)*(z1**2+z2**2)/z1/z2 ) &
   )


   s13 = -s15
   s23 = -s14 
   s12 = s45
   s123 = s12 + s13 + s23 

   z1 = z5 
   z2 = z4 


!   t123 = two*(z1*s23 - z2*s13)/(z1+z2) + (z1-z2)/(z1+z2)*s12   

   t123 =  -t123

   TC(1,0) = TC(1,0) + Cf**2*(  s123**2/two/s13/s23*z3*( (one+z3**2)/z1/z2  ) & 
             + s123/s13*(  ( z3*(one-z1) + (one-z2)**3)/z1/z2 ) - s23/s13) 

   TC(1,1) = TC(1,1) + Cf**2*(  s123**2/two/s13/s23*z3*( -(z1**2+z2**2)/z1/z2 - one  ) & 
             + s123/s13*( -one*(z1**2+z1*z2+z2**2)*(one-z2)/z1/z2 ) +one + two*s23/s13 )


   TC(2,0)  = TC(2,0) + Cf*Ca*(   t123**2/four/s12**2 + 0.25_dp + s123**2/two/s12/s13*( ( (one-z3)**2 + two*z3 )/z2 & 
     + ( z2**2 + two*(one-z2) )/(one-z3) ) - s123**2/4.0_dp/s13/s23*z3*( ( ( one-z3)**2+two*z3 )/z1/z2  ) & 
  + s123/two/s12*(  ( z1*(two -two*z1 + z1**2) - z2*(6.0_dp-6.0_dp*z2+z2**2 ) )/z2/(one-z3)  )  & 
  + s123/two/s13*( ( (one-z2)**3 + z3**2 - z2)/z2/(one-z3)  & 
  - ( z3*(one-z1)+(one-z2)**3)/z1/z2  ) )


   TC(2,1) = TC(2,1) + Cf*Ca*(  -one*( 0.25_dp*t123**2/s12**2 + 0.25_dp) - half  &
     + s123**2/two/s12/s13*( ( -(one-z3)**2 )/z2  + ( -z2**2 )/(one-z3) ) &
   - s123**2/4.0_dp/s13/s23*z3*(    ( -( one-z3)**2)/z1/z2 + one ) & 
  + s123/two/s12*(  -( z1*(two -two*z1 + z1**2) - z2*(6.0_dp-6.0_dp*z2+z2**2 ) )/z2/(one-z3) &
       + two*(z3*(z1-two*z2)-z2)/z2/(one-z3) )  & 
  + s123/two/s13*(    -one*( (one-z2)**3 + z3**2 - z2)/z2/(one-z3)  & 
          -(  two*(one-z2)*(z2-z3)/z2/(one-z3) - z1 + z2 ) + (one-z2)*(z1**2+z2**2)/z1/z2 ) &
   )

 
!-- full 
   
   TC(3,0:1) = TC(1,0:1) + TC(2,0:1) 


  TC = four/s123**2*TC

   end function 



   function TC_C5(lambda,kinlist) 
   real(dp) :: TC_C5(3,0:1) 
   real(dp), intent(in) :: lambda, kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, w, Pqq, Pqq1,  Pgg, Pgg1,  Pggperp, z
   real(dp) :: Pggperp1
   real(dp) :: kperpprod

   Pqq(z) = (one + z**2)/(one-z)
   Pqq1(z) = -(one-z) 

   Pgg(z) = two*(z/(one-z) + (one-z)/z + z*(one-z) )
   Pgg1(z) = -two*z*(1-z) 

   Pggperp(z) =   four*z*(one-z) 
   Pggperp1(z) = -four*z*(one-z) 


   s14 = kinlist(1) 
   s15 = kinlist(3) 
   s45 = kinlist(2) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 
   z1 = kinlist(6)

   ! convert from the standalone kinematics output to the en. unord. output
   z5 = z5/(z1+z5)
   z1 = one/z1

   z = one/z1
   w = z5/(z4+z5)


   TC_C5 = zero 



!    TC_C5(1,0)  =  (Tr/Cf)*four/s45/s15 * z5 /(one - z4) *(   &     
!          (one-z4)*Pqq(z4)*Pgg(z5)*half - Pggkperp(z5)*Cf*z4 * ( one - two*(one-two*lambda)**2)  )   


!   kperpprod = (one-two*lambda)**2
   kperpprod = lambda

   z2 = one - z4


    TC_C5(2,0)  =  -Cf*Ca*four/s45/s15 * z5 /(one - z4) *(   &     
          z2*Pqq(z4)*Pgg(z5)  + half*Pggperp(z5) * ( z2**2 + four*z4*kperpprod - z2*Pqq(z4) ) )


    TC_C5(2,1)  =  -Cf*Ca*four/s45/s15 * z5 /(one - z4) *(   &     
          z2*Pqq1(z4)*Pgg(z5) + z2*Pqq(z4)*Pgg1(z5)   & 
     +  half*Pggperp1(z5) * ( z2**2 + four*z4*kperpprod - z2*Pqq(z4) ) &
     +  half*Pggperp(z5) * ( - z2*Pqq1(z4) )   ) 


   TC_C5(3,0:1) = TC_C5(1,0:1) + TC_C5(2,0:1) 


   end function 




   function TC_C54(kinlist) 
   real(dp) :: TC_C54(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, w, Pqq, Pqq1,  Pgg, Pgg1,  Pggperp, z
   real(dp) :: Pggperp1

   Pqq(z) = (one + z**2)/(one-z)
   Pqq1(z) = -(one-z) 

   

   s14 = kinlist(1) 
   s15 = -kinlist(3) 
   s45 = -kinlist(2) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 
   z1 = kinlist(6)


   ! convert from the standalone kinematics output to the en. unord. output
   z5 = z4 + z5
   z4 = -z4/(z1-one)
   z1 = one/z1

      z = one/z1
   w = z5/(z4+z5)


   TC_C54 = zero 

   TC_C54(1,0) = four*Cf**2 * (one/s15) * Pqq(z4) * one/(s14 - s45) * Pqq(z5)
   TC_C54(1,1) = four*Cf**2 * (one/s15) * ( one/(s14 - s45) ) * ( Pqq1(z4) * Pqq(z5) + Pqq(z4) * Pqq1(z5) )

   TC_C54(3,0:1) = TC_C54(1,0:1) + TC_C54(2,0:1) 


   end function 



   function TC_C4(kinlist) 
   real(dp) :: TC_C4(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, w, Pqq, Pqq1,  Pgg, Pgg1,  Pggperp, z
   real(dp) :: Pggperp1

   Pqq(z) = (one + z**2)/(one-z)
   Pqq1(z) = -(one-z) 

   
   

   s14 = kinlist(1) 
   s15 = kinlist(3) 
   s45 = kinlist(2) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 
   z1 = kinlist(6)

   ! convert from the standalone kinematics output to the en. unord. output
   z5 = z1+z4
   z4 = z4/(z1+z4)
   z1 = one/z1

   z = one/z1
   w = z5/(z4+z5)


   TC_C4 = zero 

   TC_C4(1,0) = four*Cf**2 * (-one/s14) * Pqq(z4) * one/(s15 - s45) * Pqq(z5)
   TC_C4(1,1) = four*Cf**2 * (-one/s14)*(one/(s15-s45) ) * ( Pqq(z4)*Pqq1(z5) + Pqq1(z4)*Pqq(z5) )

   TC_C4(3,0:1) = TC_C4(1,0:1) + TC_C4(2,0:1) 


   end function 




   function TC_S5(kinlist) 
   real(dp) :: TC_S5(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s24,s45,s25,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, z
   real(dp) :: Pqq, Pqq1


   Pqq(z) = (one + z**2)/(one-z)
   Pqq1(z) = -(one-z) 


   s24 = kinlist(1) 
   s45 = kinlist(3) 
   s25 = kinlist(2) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 

   z = kinlist(6)

   ! convert from the standalone kinematics output to the en. unord. output
   z4 = -z4/z
   z5 = -z5/z
   z = one/z
   z1 = (one-z)/(-z)



   TC_S5(1,0) = -Cf**2*( 8.0_dp/s24/s45*z4/z5  ) * Pqq(z1)
   TC_S5(1,1) = -Cf**2*( 8.0_dp/s24/s45*z4/z5  ) * Pqq1(z1)


   TC_S5(2,0) = -Cf*Ca*(four)*( -one/s24/s45*z4/z5  + one/s24/s25/z5 + one/s25/s45 ) *Pqq(z1) 
   TC_S5(2,1) = -Cf*Ca*(four)*( -one/s24/s45*z4/z5  + one/s24/s25/z5 + one/s25/s45 ) *Pqq1(z1) 


    TC_S5(3,0:1) = TC_S5(1,0:1) + TC_S5(2,0:1) 


   end function 


   function TC_S5_C5(kinlist) 
   real(dp) :: TC_S5_C5(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, z, Pqqonemz , Pqq1onemz
   real(dp) :: s24, s25 
   real(dp) :: Pqq, Pqq1


   Pqq(z) = (one + z**2)/(one-z)
   Pqq1(z) = -(one-z) 


   TC_S5_C5 = zero 

   s24 = kinlist(1) 
   s45 = kinlist(3) 
   s25 = kinlist(2) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 

   z = kinlist(6)

   ! convert from the standalone kinematics output to the en. unord. output
   z4 = -z4/z
   z5 = -z5/z
   z = one/z
   z1 = (one-z)/(-z)


   TC_S5_C5(2,0) = -two*Cf*Ca*(four)/s25/s24/z5 *Pqq(z1) 
   TC_S5_C5(2,1) = -two*Cf*Ca*(four)/s25/s24/z5 *Pqq1(z1) 

   TC_S5_C5(3,0:1) = TC_S5_C5(1,0:1) + TC_S5_C5(2,0:1) 


   end function 

!  --  triple collinear, soft Eq.C55 

   function TC_S5_C4(kinlist) 
   real(dp) :: TC_S5_C4(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, z, Pqqonemz, Pqq1onemz 
   real(dp) :: Pqq, Pqq1, s24, s25


   Pqq(z) = (one + z**2)/(one-z)
   Pqq1(z) = -(one-z) 


   TC_S5_C4 = zero 

   s24 = kinlist(1) 
   s45 = kinlist(3) 
   s25 = kinlist(2) 
   z4 = kinlist(4) 
   z5 = kinlist(5)
   z = kinlist(6)

   ! convert from the standalone kinematics output to the en. unord. output
   z4 = -z4/z
   z5 = -z5/z
   z = one/z
   
   z1 = (one-z)/(-z)


    TC_S5_C4(1,0) = -Cf**2*( 8.0_dp/s24/s45*z4/z5  )*Pqq(z1)
    TC_S5_C4(1,1) = -Cf**2*( 8.0_dp/s24/s45*z4/z5 )*Pqq1(z1)

    TC_S5_C4(3,0:1) = TC_S5_C4(1,0:1) + TC_S5_C4(2,0:1) 

   end function 


!  --  triple collinear, soft  from Eq.C55

   function TC_S5_C54(kinlist) 
   real(dp) :: TC_S5_C54(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, z, Pqqonemz, Pqq1onemz 
   real(dp) :: s24, s25 
   real(dp) :: Pqq, Pqq1


   Pqq(z) = (one + z**2)/(one-z)
   Pqq1(z) = -(one-z) 


   TC_S5_C54 = zero 

   s24 = kinlist(1) 
   s45 = kinlist(3) 
   s25 = kinlist(2) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 

   z = kinlist(6)

   ! convert from the standalone kinematics output to the en. unord. output
   z4 = -z4/z
   z5 = -z5/z
   z = one/z
   
   z1 = (one-z)/(-z)



   TC_S5_C54(1,0) = -Cf**2*( 8.0_dp/s24/s45*z4/z5  )*Pqq(z1)
   TC_S5_C54(1,1) = -Cf**2*( 8.0_dp/s24/s45*z4/z5 )*Pqq1(z1)


    TC_S5_C54(3,0:1) = TC_S5_C54(1,0:1) + TC_S5_C54(2,0:1) 


   end function 


! --- The triple collinear Altarelli-Parisi function; sector ``A'' -----------------------


   function PTC_A_g_qg_qstar(z,zmin,yr,x3)
   real(dp) :: PTC_A_g_qg_qstar(1:3,1:3,0:1)
   real(dp), intent(in) :: yr(3), z,zmin,x3
   real(dp) ::  rho51, full(6), fullc51(6),soft5(6), soft5c5(6),dsoft(6), dsord(6), dsoftz(6), dsoftzc5(6)
   real(dp) :: dsoft5(6), dsordc5(6), dsoftc5(6)
   real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4 , PHSP1m1, energy, onemz
   real(dp) :: faux1(1:3,0:1), faux2(1:3,0:1), LogR, LogF, LogF0, faux3(1:3,0:1), faux4(1:3,0:1), faclambda 
   real(dp) :: test1(1:3,0:1), test2(1:3,0:1), test3(1:3,0:1), test4(1:3,0:1),kt1kt2sq 


   r =  yr(1) 
   x4 = yr(2) 
   lambda = yr(3) 

   PTC_A_g_qg_qstar = zero 

   call kinematics_tcl_unord_A(z,(/r,x4,lambda,x3/),rho51,full,fullc51,soft5,soft5c5)
   PHSP1 = (one - r) 
   PHSP2 =  (one-half*x4)/NF1(x3,0.5_dp*x4,lambda)

   faclambda = lambda*(one-lambda) 

   onemz = one-z
   onemz4 = onemz**4


    LogR  = log( r*(one-r))
    LogF  = log( 128.0_dp * PHSP2**2 * faclambda * (one-x3) * (one-x3*x4*half) )
    LogF0 = log( 128.0_dp * faclambda * (one-x3) )

    kt1kt2sq = (two*lambda-one)**2 * (one-x3)

   
!  the regular piece 
    faux1(1:3,0:1) =  rho51 * PHSP2 * PHSP1* onemz**2 * ( r**2 * TC(full)  - TC_S5(soft5)  )

    faux2(1:3,0:1) =                  PHSP1* onemz**2 * ( r**2 * TC_C5(kt1kt2sq,fullc51)  - TC_S5_C5(soft5c5)  )

    faux3(1:3,0:1) =   r*onemz**2 *  rho51 * PHSP2 * TC_S5(soft5)  

    faux4(1:3,0:1) =   r*onemz**2 * TC_S5_C5(soft5c5)  

   PTC_A_g_qg_qstar(1:3,3,0:1) = faux1 - faux2 

   PTC_A_g_qg_qstar(1:3,3,0) = PTC_A_g_qg_qstar(1:3,3,0) + ( faux3(1:3,0) - faux4(1:3,0) ) * ( log(onemz) - one)

 
   PTC_A_g_qg_qstar(1:3,3,1) = PTC_A_g_qg_qstar(1:3,3,1)    -  two*LogR*(faux1(1:3,0) - faux2(1:3,0) ) & 
                                            -  four*log(onemz)*PTC_A_g_qg_qstar(1:3,3,0)  & 
                                     - ( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) )  &
                          + ( faux3(1:3,0) - faux4(1:3,0) ) *( log(onemz)**2 - four + pisq/three ) & 
                                     - ( LogF*faux3(1:3,0) - LogF0*faux4(1:3,0) )*(log(onemz) - one)  &
                                     + ( faux3(1:3,1) - faux4(1:3,1) ) * ( log(onemz) - one)

! Raoul added

   PTC_A_g_qg_qstar(1:3,1:3,1) = PTC_A_g_qg_qstar(1:3,1:3,1) - PTC_A_g_qg_qstar(1:3,1:3,0)*log(x4)
   
! Normalization 

   PTC_A_g_qg_qstar = (one-z)*z * PTC_A_g_qg_qstar 

   end function 

 


   function PTC_B_g_qg_qstar(z,zmin,yr,x3)
   real(dp) :: PTC_B_g_qg_qstar(3,3,0:1)
   real(dp), intent(in) :: yr(3), z,zmin,x3
   real(dp) ::  rho45, full(6), fullc54(6),soft5(6), soft5c5(6),dsoft(6), dsord(6), dsoftz(6), dsoftzc5(6)
   real(dp) :: dsoft5(6), dsordc5(6), dsoftc5(6)
   real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4 , PHSP1m1, energy, onemz
   real(dp) :: faux1(1:3,0:1), faux2(1:3,0:1), LogR, LogF, LogF0, faux3(1:3,0:1), faux4(1:3,0:1), faclambda 
   real(dp) :: test1(1:3,0:1), test2(1:3,0:1), test3(1:3,0:1), test4(1:3,0:1)  

   r =  yr(1) 
   x4 = yr(2) 
   lambda = yr(3) 

   PTC_B_g_qg_qstar = zero 

   call kinematics_tcl_unord_B(z,(/r,x4,lambda,x3/),rho45,full,fullc54,soft5,soft5c5)

   PHSP1 = (one - r) 
   PHSP2 =  one   ! evrything is absorbed into rho45 

   faclambda = lambda*(one-lambda) 

    LogR =  log( r*(one-r))
    LogF =  log( 64.0_dp* (one-half*x4)/NF1(x3,one-half*x4,lambda)**2*faclambda * (one-x3)*(one-x3*(one-x4*half)))
    LogF0 = log( 64.0_dp*one/NF1(x3,one,lambda)**2 * faclambda * (one-x3)**2 )


!  only the regular piece is there 
   onemz = one-z
   
  

    faux1(1:3,0:1) =  rho45 * PHSP2 * PHSP1* onemz**2 * ( r**2 * TC(full)  - TC_S5(soft5)  )
    faux2(1:3,0:1) =                  PHSP1* onemz**2 * ( r**2 * TC_C54(fullc54)  - TC_S5_C54(soft5c5)  )
    faux3(1:3,0:1) =   r * onemz**2 *  rho45 * PHSP2 * TC_S5(soft5)  
    faux4(1:3,0:1) =   r * onemz**2 * TC_S5_C54(soft5c5)  

   PTC_B_g_qg_qstar(1:3,3,0:1) = faux1 - faux2 
   PTC_B_g_qg_qstar(1:3,3,0) = PTC_B_g_qg_qstar(1:3,3,0) + ( faux3(1:3,0) - faux4(1:3,0) ) * ( log(onemz) - one)


    PTC_B_g_qg_qstar(1:3,3,1)  = PTC_B_g_qg_qstar(1:3,3,1)  - two*LogR*(faux1(1:3,0) - faux2(1:3,0) ) &  
                                        - four*log(onemz)*PTC_B_g_qg_qstar(1:3,3,0) & 
                                     - ( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) )  &
                          + ( faux3(1:3,0) - faux4(1:3,0) ) *( log(onemz)**2 - four + pisq/three ) & 
                                     - ( LogF*faux3(1:3,0) - LogF0*faux4(1:3,0) ) * (log(onemz) - one)   &
                                     + ( faux3(1:3,1) - faux4(1:3,1) ) * ( log(onemz) - one)


! Raoul added
   PTC_B_g_qg_qstar(1:3,1:3,1) = PTC_B_g_qg_qstar(1:3,1:3,1) - PTC_B_g_qg_qstar(1:3,1:3,0)*two*log(x4)
   

! Normalization 

   PTC_B_g_qg_qstar = (one-z)*z * PTC_B_g_qg_qstar 

   end function 



   function PTC_C_g_qg_qstar(z,zmin,yr,x3)
   real(dp) :: PTC_C_g_qg_qstar(1:3,1:3,0:1)
   real(dp), intent(in) :: yr(3), z,zmin,x3
   real(dp) ::  rho41, full(6), fullc41(6),soft5(6), soft5c41(6)
   real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4 , PHSP1m1, energy, onemz
   real(dp) :: faux1(1:3,0:1), faux2(1:3,0:1), LogR, LogF, LogF0, faux3(1:3,0:1), faux4(1:3,0:1), faclambda 
   real(dp) :: test1(1:3,0:1), test2(1:3,0:1), test3(1:3,0:1), test4(1:3,0:1)  

   r =  yr(1) 
   x4 = yr(2) 
   lambda = yr(3) 

   PTC_C_g_qg_qstar = zero 

   call kinematics_tcl_unord_C(z,(/r,x4,lambda,x3/),rho41,full,fullc41,soft5,soft5c41)

   PHSP1 = (one - r) 
   PHSP2 =  (one-half*x4)/NF1(x3,0.5_dp*x4,lambda)

   faclambda = lambda*(one-lambda) 

   onemz = one-z
  
    LogR  = log( r*(one-r))
    LogF  = log( 128.0_dp * PHSP2**2 * faclambda * (one-x3) * (one-x3*x4*half) )
    LogF0 = log( 128.0_dp * faclambda * (one-x3)  )

!    kt1kt2sq = (two*lambda-one)**2 * (one-x3)



!  the regular piece 


    faux1(1:3,0:1) =  rho41 * PHSP2 * PHSP1* onemz**2 * ( r**2 * TC(full)  - TC_S5(soft5)  )

    faux2(1:3,0:1) =                  PHSP1* onemz**2 * ( r**2 * TC_C4(fullc41)  - TC_S5_C4(soft5c41)  )

    faux3(1:3,0:1) =   r*onemz**2 *  rho41 * PHSP2 * TC_S5(soft5)  

    faux4(1:3,0:1) =   r*onemz**2 * TC_S5_C4(soft5c41)  

   PTC_C_g_qg_qstar(1:3,3,0:1) = faux1 - faux2 

   PTC_C_g_qg_qstar(1:3,3,0) = PTC_C_g_qg_qstar(1:3,3,0) + ( faux3(1:3,0) - faux4(1:3,0) ) * ( log(onemz) - one)

 
   PTC_C_g_qg_qstar(1:3,3,1) = PTC_C_g_qg_qstar(1:3,3,1)  - two*LogR*(faux1(1:3,0) - faux2(1:3,0) ) &  
                                          - four*log(onemz) *PTC_C_g_qg_qstar(1:3,3,0)  & 
                                     - ( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) )  &
                          + ( faux3(1:3,0) - faux4(1:3,0) ) *( log(onemz)**2 - four + pisq/three) & 
                                     - ( LogF*faux3(1:3,0) - LogF0*faux4(1:3,0) )* ( log(onemz) - one)  &
                                     + ( faux3(1:3,1) - faux4(1:3,1) ) * ( log(onemz) - one)

! Raoul added

   PTC_C_g_qg_qstar(1:3,1:3,1) = PTC_C_g_qg_qstar(1:3,1:3,1) - PTC_C_g_qg_qstar(1:3,1:3,0)*log(x4)
   
! Normalization 

   PTC_C_g_qg_qstar = (one-z)*z * PTC_C_g_qg_qstar 

   end function 


   function PTC_D_g_qg_qstar(z,zmin,yr,x3)
   real(dp) :: PTC_D_g_qg_qstar(1:3,1:3,0:1)
   real(dp), intent(in) :: yr(3), z,zmin,x3
   real(dp) ::  rho45, full(6), fullc54(6),soft5(6), soft5c5(6)
   real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4 , PHSP1m1, energy, onemz
   real(dp) :: faux1(1:3,0:1), faux2(1:3,0:1), LogR, LogF, LogF0, faux3(1:3,0:1), faux4(1:3,0:1), faclambda 
   real(dp) :: test1(1:3,0:1), test2(1:3,0:1), test3(1:3,0:1), test4(1:3,0:1)  

   r =  yr(1) 
   x4 = yr(2) 
   lambda = yr(3) 

   PTC_D_g_qg_qstar = zero 


   call kinematics_tcl_unord_D(z,(/r,x4,lambda,x3/),rho45,full,fullc54,soft5,soft5c5)

   PHSP1 = (one - r) 
   PHSP2 =  one   ! evrything is absorbed into rho45 

   faclambda = lambda*(one-lambda) 

    LogR =  log( r*(one-r))
    LogF =  log( 64.0_dp* (one-half*x4)/NF1(x3,one-0.5_dp*x4,lambda)**2*faclambda * (one-x3)*(one-x3*(one-x4*half)))
    LogF0 = log( 64.0_dp*one/NF1(x3,one,lambda)**2 * faclambda * (one-x3)**2 )



!  only the regular piece is there 
    onemz = one-z
  

    faux1(1:3,0:1) =  rho45 * PHSP2 * PHSP1* onemz**2 * ( r**2 * TC(full)  - TC_S5(soft5)  )
    faux2(1:3,0:1) =                  PHSP1* onemz**2 * ( r**2 * TC_C54(fullc54)  - TC_S5_C54(soft5c5)  )
    faux3(1:3,0:1) =   r * onemz**2 *  rho45 * PHSP2 * TC_S5(soft5)  
    faux4(1:3,0:1) =   r * onemz**2 * TC_S5_C54(soft5c5)  

   PTC_D_g_qg_qstar(1:3,3,0:1) = faux1 - faux2 
   PTC_D_g_qg_qstar(1:3,3,0) = PTC_D_g_qg_qstar(1:3,3,0) + ( faux3(1:3,0) - faux4(1:3,0) ) * ( log(onemz) - one)


    PTC_D_g_qg_qstar(1:3,3,1)  = PTC_D_g_qg_qstar(1:3,3,1)  - two*LogR*(faux1(1:3,0) - faux2(1:3,0) ) &  
                                           -  four*log(onemz)*PTC_D_g_qg_qstar(1:3,3,0) & 
                                     - ( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) )  &
                          + ( faux3(1:3,0) - faux4(1:3,0) ) *( log(onemz)**2 - four + pisq/three ) & 
                                     - ( LogF*faux3(1:3,0) - LogF0*faux4(1:3,0) )* ( log(onemz) - one)  &
                                     + ( faux3(1:3,1) - faux4(1:3,1) ) * ( log(onemz) - one)


! Raoul added
   PTC_D_g_qg_qstar(1:3,1:3,1) = PTC_D_g_qg_qstar(1:3,1:3,1) - two*PTC_D_g_qg_qstar(1:3,1:3,0)*log(x4)
   

! Normalization 

   PTC_D_g_qg_qstar = (one-z)*z * PTC_D_g_qg_qstar 

   end function 



 end module mod_tripleclimits_g_qg_qstar



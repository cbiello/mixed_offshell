module mod_tripleclimits_q_gg_qstar
  use mod_types; use mod_consts_dp; use common_def
  use mod_kinematics_tc_ord
  use mod_auxfunctions

  implicit none

  private

  

  public :: PTC_A_q_gg_qst, PTC_B_q_gg_qst, PTC_C_q_gg_qst, PTC_D_q_gg_qst
  public :: tripcoll_subtr_q_gg_qstar


contains


  function tripcoll_subtr_q_gg_qstar(ndim,yRnd,ncomp,ff,userdata,weight,iternumber)
! output indices: 1st is color {1=ab,2=nab,3=tot},2nd is ep power
    implicit none
   !--
    integer, parameter :: ep = 1
    integer, parameter :: mycol = 1 !-- ab -> 1, nab -> 2, tot -> 3
    integer, parameter :: mycontr =1 !-- 1 -> delta, 2 -> plus, 3 -> regular
    real(dp), parameter :: myz = 0.123_dp
    real(dp), parameter :: zmin = zero
    !--
    real(dp), parameter :: buff = 1E-10_dp
    real(dp), parameter :: onet = one - two * buff
    real(dp), parameter :: bufftc = 1E-10_dp
    !--
    real(dp15) :: yRnd(30), ff(1), weight
    integer :: ndim, ncomp, userdata, iternumber,tripcoll_subtr_q_gg_qstar
    real(dp)              :: x1,x2,x3,xx2,xx(1:3),x3a,xx3
    real(dp)              :: res_tc(1:3,1:3,0:1)
    real(dp)              :: seca(1:3,1:3,0:1),secb(1:3,1:3,0:1),secc(1:3,1:3,0:1),secd(1:3,1:3,0:1),wghtac,wghtbd,z
    real(dp)              :: seca0(1:3,1:3,0:1),secb0(1:3,1:3,0:1),secc0(1:3,1:3,0:1),secd0(1:3,1:3,0:1)

    tripcoll_subtr_q_gg_qstar = zero
    z = myz
    ff(1) = zero
    res_tc = zero

    seca = zero; seca0 = zero
    secb = zero; secb0 = zero
    secc = zero; secc0 = zero
    secd = zero; secd0 = zero

    x1  = buff+onet*real(yRnd(1),dp)  ! r
    x2  = buff+onet*real(yRnd(2),dp)  ! x4
    x3a = buff+onet*real(yRnd(3),dp)  ! lambda
    xx3 = buff +onet*real(yRnd(4),dp) ! x3          

! remap lambda
    x3 = sin( half * pi * x3a )**2

    if (x1*x2*(one-z) .lt. bufftc) return

    xx = (/x1,x2,x3/)
    seca0 = PTC_A_q_gg_qst(z,zmin,xx,zero)
    secc0 = PTC_C_q_gg_qst(z,zmin,xx,zero)

    if (withmax) then
       seca = PTC_A_q_gg_qst(z,zmin,xx,xx3)
       secc = PTC_C_q_gg_qst(z,zmin,xx,xx3)
    endif
    wghtac = -one/x1/x2

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
!    print *, x1,x2,x3a,xx3
!    print *, "seca0",seca0(2,3,1)
!    print *, "seca",seca(2,3,1)
!    print *, "secb0",secb0(2,3,1)
!    print *, "secb",secb(2,3,1)
!    print *, "secc0",secc0(2,3,1)
!    print *, "secc",secc(2,3,1)
!    print *, "secd0",secd0(2,3,1)
!    print *, "secd",secd(2,3,1)
!    print *, "ab 1/ep",res_tc(1,:,0)
!    print *, "ab ep^0",res_tc(1,:,1)
!    print *, "nab 1/ep",res_tc(2,:,0)
!    print *, "nab ep^0",res_tc(2,:,1)
!    pause

! remove color factors
    res_tc(1,:,:) = res_tc(1,:,:)/Cf**2
    res_tc(2,:,:) = res_tc(2,:,:)/Cf/Ca


    ff(1) = res_tc(mycol,mycontr,ep)
    
    if (ff(1) .ne. ff(1) .or. weight .ne. weight) then
       print *, "xx",xx,x3a
       print *, "wghtbd",wghtbd
       print *, "sec a",seca(3,mycontr,ep)
       print *, "sec b",secb(3,mycontr,ep)
       print *, "sec c",secc(3,mycontr,ep)
       print *, "sec d",secd(3,mycontr,ep)
       print *, "ff",ff(1)
       stop
    endif




  end function tripcoll_subtr_q_gg_qstar

!    triple collinear splitting function  Eqs. in sect. C.5

   function TC(kinlist) 
   real(dp) :: TC(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3 

   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 


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


   function TC_C5(kinlist)          ! Eq.(C1.35)
   real(dp) :: TC_C5(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, wz5, wz45, Pqq, Pqq1, z

   Pqq(z) = (one + z**2)/(one-z)
   Pqq1(z) = -(one-z)

   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 
   z1 = kinlist(6)


   TC_C5 = zero 


   TC_C5(1,0) = Cf**2 *four/s15/s14 * z1/(z1+z5) * Pqq(z1/(z1+z5) ) *Pqq(z1+z5)

   TC_C5(1,1) = Cf**2 *four/s15/s14 * z1/(z1+z5) * ( Pqq1(z1/(z1+z5) ) *Pqq(z1+z5)  + Pqq(z1/(z1+z5) ) *Pqq1(z1+z5) )

   TC_C5(2,0) = zero 

   TC_C5(3,0:1) = TC_C5(1,0:1) + TC_C5(2,0:1) 


   end function 


!---  triple collinear, angular ordered 4||1


   function TC_C4(kinlist) 
   real(dp) :: TC_C4(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, wz5, wz45, Pqq, Pqq1, z

   Pqq(z) = (one + z**2)/(one-z)
   Pqq1(z) = -(one-z)

   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 
   z1 = kinlist(6)


   TC_C4 = zero 

!   TC_C4(1) = Cf**2 *four/s15/s14 * z1/(z1+z4) * Pqq(z1/(z1+z4) ) *Pqq(z1+z4)
!( one + (one-z5)**2)/(z5) 

   TC_C4(1,0) = Cf**2 *four/s15/s14 * z1/(z1+z4) * Pqq( z1/(z1+z4) ) * Pqq(one-z5)

   TC_C4(1,1) = Cf**2 *four/s15/s14 * z1/(z1+z4) * ( Pqq( z1/(z1+z4) ) * Pqq1(one-z5) & 
                                                      + Pqq1( z1/(z1+z4) ) * Pqq(one-z5) )

   TC_C4(2,0:1) = zero 

   TC_C4(3,0:1) = TC_C4(1,0:1) + TC_C4(2,0:1) 


   end function 

!--- triple collinear, strong angular, 5||4 Eq.(C.40) 

   function TC_C54(lambda,kinlist) 
   real(dp) :: TC_C54(3,0:1) 
   real(dp), intent(in) :: lambda, kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, w, Pqq, Pqq1,  Pgg, Pgg1,  Pggperp, z
   real(dp) :: Pggperp1

   Pqq(z) = (one + z**2)/(one-z)
   Pqq1(z) = -(one-z) 

   Pgg(z) = two*(z/(one-z) + (one-z)/z + z*(one-z) )
   Pgg1(z) = -two*z*(1-z) 

   Pggperp(z) =   four*z*(one-z) 
   Pggperp1(z) = -four*z*(one-z) 
   

   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 
   z1 = kinlist(6)

   z = one/z1
   w = z5/(z4+z5)

   TC_C54 = zero 

   TC_C54(1,0:1) = zero 

   TC_C54(2,0) = Cf*Ca * four/s45/s14/z * (z4/(z4+z5))**2*(-z1/z5) *w/(one-w)   & 
       *(  ( one-z)*Pqq(z)*Pgg(w) + Pggperp(w)* two*z*(lambda-half)  )


   TC_C54(2,1) = Cf*Ca * four/s45/s14/z * (z4/(z4+z5))**2*(-z1/z5) *w/(one-w)   & 
       *(  ( one-z)*( Pqq1(z)*Pgg(w) + Pqq(z)*Pgg1(w) )  + Pggperp1(w)*two*z*(lambda-half) & 
               + Pggperp(w)*( -two*z + (one+z)**2*half ) )


   TC_C54(3,0:1) = TC_C54(1,0:1) + TC_C54(2,0:1) 


   end function 



!  --  triple collinear, soft Eq.C55

   function TC_S5(kinlist) 
   real(dp) :: TC_S5(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, z, Pqqonemz, Pqq1onemz 

   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 

   z = one/(one-z4)

   Pqqonemz = one + z**2
   Pqq1onemz = -(one -z)**2


   TC_S5(1,0) = Cf**2*( 8.0_dp/s14/s15/z4/z/z5/z )*Pqqonemz

   TC_S5(1,1) = Cf**2*( 8.0_dp/s14/s15/z4/z/z5/z )*Pqq1onemz


   TC_S5(2,0) = Cf*Ca*(-four)*(one/s14/s15/z4/z/z5/z + one/s15/s45/z4/z + one/s14/s45/z5/z )*Pqqonemz
 
   TC_S5(2,1) = Cf*Ca*(-four)*(one/s14/s15/z4/z/z5/z + one/s15/s45/z4/z + one/s14/s45/z5/z )*Pqq1onemz


   TC_S5(3,0:1) = TC_S5(1,0:1) + TC_S5(2,0:1) 


   end function 


!  --  triple collinear, soft Eq.C55 


   function TC_S5_C5(kinlist) 
   real(dp) :: TC_S5_C5(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, z, Pqqonemz , Pqq1onemz

   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 

   z = one/(one-z4)
   Pqqonemz = one + z**2
   Pqq1onemz = -(one-z)**2

   TC_S5_C5 = zero 

   TC_S5_C5(1,0) = Cf**2*( 8.0_dp/s14/s15/z/z4/z5/z )*Pqqonemz

   TC_S5_C5(1,1) = Cf**2*( 8.0_dp/s14/s15/z/z4/z5/z )*Pqq1onemz

   TC_S5_C5(2,0:1) = zero

   TC_S5_C5(3,0:1) = TC_S5_C5(1,0:1) + TC_S5_C5(2,0:1) 

   end function 

!  --  triple collinear, soft Eq.C55 

   function TC_S5_C4(kinlist) 
   real(dp) :: TC_S5_C4(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, z, Pqqonemz, Pqq1onemz 

   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 

   z = one/(one-z4)
   Pqqonemz = one + z**2
   Pqq1onemz = -(one-z)**2

   TC_S5_C4 = zero 

   TC_S5_C4(1,0) = Cf**2*( 8.0_dp/s14/s15/z/z4/z5/z )*Pqqonemz
   TC_S5_C4(1,1) = Cf**2*( 8.0_dp/s14/s15/z/z4/z5/z )*Pqq1onemz

   TC_S5_C4(2,0:1) = zero

   TC_S5_C4(3,0:1) = TC_S5_C4(1,0:1) + TC_S5_C4(2,0:1) 


   end function 


!  --  triple collinear, soft  from Eq.C55

   function TC_S5_C54(kinlist) 
   real(dp) :: TC_S5_C54(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, z, Pqqonemz, Pqq1onemz 

   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 

   z = one/(one-z4)
   Pqqonemz = one + z**2
   Pqq1onemz = -(one-z)**2

   TC_S5_C54(1,0:1) = zero 

   TC_S5_C54(2,0) = Ca*Cf*(8.0_dp)/s45/s14 *(z4-one)/z5 * Pqqonemz 
   TC_S5_C54(2,1) = Ca*Cf*(8.0_dp)/s45/s14 *(z4-one)/z5 * Pqq1onemz 

   TC_S5_C54(3,0:1) = TC_S5_C54(1,0:1) + TC_S5_C54(2,0:1) 


   end function 



! -- double soft, triple collinear splitting function Eqs C.45-C.46

  function DS_TC(kinlist)
   real(dp) :: DS_TC(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5, z45, s145, t145 

   s14 = -kinlist(1) 
   s15 = -kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 
   
   DS_TC = zero 

   z45 = z4 + z5 
   s145 = s14 + s15
   t145 = two*(z4*s15 - z5*s14)/z45

!  Abelian piece

   DS_TC(1,0) = Cf**2*16.0_dp/s14/s15/z4/z5

!  non-Abelian piece ( O(ep^0) ) term 

   DS_TC(2,0) = Cf*Ca*four*( t145**2/two/s45**2/s145**2 + one/s45*(one/z5/s14 + one/z4/s15)  & 
      +  one/s45/z45*( one/s14 + one/s15 )  - one/s14/s15/z4/z5 + one/s145/s45/z45*( z4/z5 + z5/z4 - six)  & 
      + one/s145/z45*(one/z5/s14 + one/z4/s15 ) - one/s145/z4/z5*(one/s14 + one/s15 ) )


   DS_TC(2,1) = -Cf*Ca*four* t145**2/two/s45**2/s145**2


!  full 

   DS_TC(3,0:1) = DS_TC(1,0:1) + DS_TC(2,0:1)


  end function


! double soft, triple collinear, strongly ordered  Eq.(C.59)

  function DS_SO_TC(kinlist)
   real(dp) :: DS_SO_TC(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5

   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 
   
   DS_SO_TC = zero 

!  Abelian piece

   DS_SO_TC(1,0)  = Cf**2*16.0_dp/s14/s15/z4/z5

!  non-Abelian piece  term 

   DS_SO_TC(2,0) = Cf*Ca*(-8.0_dp)*( one/s14/s15/z4/z5 + one/s15/s45/z4 + one/s14/s45/z5 )

!  full 

   DS_SO_TC(3,0:1) = DS_SO_TC(1,0:1) + DS_SO_TC(2,0:1)

  end function

!---- double soft, triple collinear, strong angular Eq.(C51)

  function DS_TC_C5(kinlist)
  real(dp) :: DS_TC_C5(3,0:1) 
  real(dp), intent(in) :: kinlist(6)  
  real(dp) :: s14,s15,s45,z4,z5
 
   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 =  kinlist(4) 
   z5 =  kinlist(5) 
 
   DS_TC_C5 = zero 

   DS_TC_C5(1,0) = Cf**2*16.0_dp/s14/s15/z4/z5
   DS_TC_C5(2,0) = zero 
   DS_TC_C5(3,0:1) = DS_TC_C5(1,0:1) + DS_TC_C5(2,0:1) 
    
  end function 


!---- double soft, triple collinear, strong angular Eq.(C50)

  function DS_TC_C4(kinlist)
  real(dp) :: DS_TC_C4(3,0:1) 
  real(dp), intent(in) :: kinlist(6)  
  real(dp) :: s14,s15,s45,z4,z5
 
   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 
 
   DS_TC_C4 = zero  

   DS_TC_C4(1,0) = Cf**2*16.0_dp/s14/s15/z4/z5
   DS_TC_C4(2,0) = zero 

   DS_TC_C4(3,0:1) = DS_TC_C4(1,0:1) + DS_TC_C4(2,0:1) 
    
  end function 


!---- double soft, triple collinear, strong angular Eq.(C54)

  function DS_TC_C54(lambda,kinlist)
  real(dp) :: DS_TC_C54(3,0:1) 
  real(dp), intent(in) :: kinlist(6), lambda  
  real(dp) :: s14,s15,s45,z4,z5,z, Pgg, Pggperp, Pgg1, Pggperp1
 
  Pgg(z) = two*(z/(one-z) + (one-z)/z + z*(one-z) )
  Pgg1(z) = -two*z*(one-z)

  Pggperp(z) = four*z*(one-z) 
  Pggperp1(z) = -four*z*(one-z)


   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 
 
   z = z5/(z4+z5)

   DS_TC_C54 = zero 

   DS_TC_C54(1,0:1) = zero

   DS_TC_C54(2,0) = Cf*Ca*(8.0_dp)/s45/s14*(-z4)/(z4+z5)**2 & 
        *( Pgg(z) + (lambda - half)*Pggperp(z) )

   DS_TC_C54(2,1) = Cf*Ca*(8.0_dp)/s45/s14*(-z4)/(z4+z5)**2 & 
        *( Pgg1(z) + Pggperp1(z)*(lambda - half) )


   DS_TC_C54(3,0:1) = DS_TC_C54(1,0:1) + DS_TC_C54(2,0:1) 
    
  end function 


! double soft, triple collinear, strongly ordered  Eq.(C.59)

   function DS_SO_TC_C5(kinlist)
   real(dp) :: DS_SO_TC_C5(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5

   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 
   

   DS_SO_TC_C5 = zero 

!  Abelian piece

   DS_SO_TC_C5(1,0)  = Cf**2*16.0_dp/s14/s15/z4/z5

!  non-Abelian piece  term 

   DS_SO_TC_C5(2,0:1) = zero

!  full 

   DS_SO_TC_C5(3,0:1) = DS_SO_TC_C5(1,0:1) + DS_SO_TC_C5(2,0:1)


  end function


  function DS_SO_TC_C4(kinlist)
   real(dp) :: DS_SO_TC_C4(3,0:1) 
   real(dp), intent(in) :: kinlist(6)
   real(dp) :: s14,s15,s45,z4,z5

   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 
   

   DS_SO_TC_C4 = zero 

!  Abelian piece

   DS_SO_TC_C4(1,0)  = Cf**2*16.0_dp/s14/s15/z4/z5

!  non-Abelian piece  term 

   DS_SO_TC_C4(2,0:1) = zero

!  full 

   DS_SO_TC_C4(3,0:1) = DS_SO_TC_C4(1,0:1) + DS_SO_TC_C4(2,0:1)


  end function



!---- double soft, triple collinear, strong angular Eq.(C54)

  function DS_SO_TC_C54(kinlist)
  real(dp) :: DS_SO_TC_C54(3,0:1) 
  real(dp), intent(in) :: kinlist(6)
  real(dp) :: s14,s15,s45,z4,z5
 
   s14 = kinlist(1) 
   s15 = kinlist(2) 
   s45 = kinlist(3) 
   z4 = kinlist(4) 
   z5 = kinlist(5) 
 
   DS_SO_TC_C54 = zero 

   DS_SO_TC_C54(1,0:1) = zero

   DS_SO_TC_C54(2,0) = Cf*Ca*16/s45/s14*(-one/z5)

   DS_SO_TC_C54(3,0:1) = DS_SO_TC_C54(1,0:1) + DS_SO_TC_C54(2,0:1) 
    
  end function 


! --- The triple collinear Altarelli-Parisi function; sector ``A'' -----------------------


   function PTC_A_q_gg_qst(z,zmin,yr,x3)
   real(dp) :: PTC_A_q_gg_qst(1:3,1:3,0:1)
   real(dp), intent(in) :: yr(3), z,zmin,x3
   real(dp) ::  rho51, full(6), fullc51(6),soft5(6), soft5c5(6),dsoft(6), dsord(6), dsoftz(6), dsoftzc5(6)
   real(dp) :: dsoft5(6), dsordc5(6), dsoftc5(6)
   real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4 , PHSP1m1, energy, onemz
   real(dp) :: faux1(1:3,0:1), faux2(1:3,0:1), LogR, LogF, LogF0, faux3(1:3,0:1), faux4(1:3,0:1), faclambda 

   r =  yr(1) 
   x4 = yr(2) 
   lambda = yr(3) 

   PTC_A_q_gg_qst = zero 

   call kinematics_tcl_A(z,(/r,x4,lambda,x3/),rho51,full,fullc51,soft5,soft5c5,dsoft,dsord,dsoftc5,dsordc5,dsoftz,dsoftzc5)

   PHSP1 = (one - half*r) 
   PHSP2 =  (one-half*x4)/NF1(x3,0.5_dp*x4,lambda)

   faclambda = lambda*(one-lambda) 
!   faclambda=one


    LogR  = log( half*r*(one-half*r))
    LogF  = log( 128.0_dp * PHSP2**2 * faclambda * (one-x3) * (one-x3*x4*half) )
    LogF0 = log( 128.0_dp * faclambda * (one-x3) )


! -- We will start with the computation of the +plus distribution 


!  the delta-function coefficient 


   faux1(1:3,0:1) = PHSP2*rho51*( r**2 * DS_TC(dsoftz) - DS_SO_TC(dsord) )
   faux2(1:3,0:1) = ( r**2 * DS_TC_C5(dsoftzc5) - DS_SO_TC_C5(dsordc5) )
   faux3(1:3,0:1) = PHSP2*rho51*( r**2 * DS_TC(dsoftz) )
   faux4(1:3,0:1) = r**2 * DS_TC_C5(dsoftzc5)                              

   PTC_A_q_gg_qst(1:3,1,0:1) = -log(one-zmin)*(faux1 - faux2)  - log(one+r)* (faux3 - faux4)  

   PTC_A_q_gg_qst(1:3,1,1) = PTC_A_q_gg_qst(1:3,1,1) + two*log(one-zmin)**2*( faux1(1:3,0) - faux2(1:3,0) ) & 
                   -two*Log(r)*PTC_A_q_gg_qst(1:3,1,0) - two*log(one+r)**2 * (faux3(1:3,0) - faux4(1:3,0) ) & 
                   +log(one-zmin)*( LogF *faux1(1:3,0) - LogF0*faux2(1:3,0) ) &
                   +log(one+r)*( LogF*faux3(1:3,0) - LogF0*faux4(1:3,0) )          

!  the plus distribution coefficient

     faux1(1:3,0:1) = rho51*PHSP2*(  PHSP1**2*( half*r )**2 * DS_TC(dsoft) - DS_SO_TC(dsord) )
     faux2(1:3,0:1) = PHSP1**2* (half*r )**2 * DS_TC_C5(dsoftc5) - DS_SO_TC_C5(dsordc5) 

! Raoul modified to include ph sp     
!     faux3(1:3,0:1) = r*(rho51*PHSP2*DS_SO_TC(dsord)  - DS_SO_TC_C5(dsordc5))   ! I multiply with r to make integration 
                                                                                ! over r  give 1 when integrated with 1/r dr measure 
     faux3(1:3,0:1) = r*PHSP1*(rho51*PHSP2*DS_SO_TC(dsord)  - DS_SO_TC_C5(dsordc5))   ! I multiply with r to make integration 
                                                                                ! over r  give 1 when integrated with 1/r dr measure 

     PTC_A_q_gg_qst(1:3,2,0:1) =  faux1 - faux2 

     PTC_A_q_gg_qst(1:3,2,1) = PTC_A_q_gg_qst(1:3,2,1) - two * LogR*PTC_A_q_gg_qst(1:3,2,0) - ( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) ) & 
                + pisq/three*faux3(1:3,0) 


!  the regular piece 

   onemz = one-z
   onemz4 = onemz**4

    faux1(1:3,0:1) =  rho51*PHSP2*(  onemz4*(  z*PHSP1**2*(half*r)**2 * TC(full)     -  TC_S5(soft5) )   &        !! 
                           - PHSP1**2*(half*r)**2 * DS_TC(dsoft) +  DS_SO_TC(dsord)  )

    faux2(1:3,0:1) =  (onemz4*(   z*PHSP1**2*(half*r)**2 * TC_C5(fullc51)  -  TC_S5_C5(soft5c5) )   &
                     - PHSP1**2*(half*r)**2 * DS_TC_C5(dsoftc5) +  DS_SO_TC_C5(dsordc5)  )

!    faux3(1:3,0:1) = r* ( rho51* PHSP2*( onemz4* TC_S5(soft5)  - DS_SO_TC(dsord) ) & 
!                         - ( onemz4*TC_S5_C5(soft5c5) - DS_SO_TC_C5(dsordc5) ) )      ! again, multiply with r to get convegent 
                                                                                       ! integration 
    faux3(1:3,0:1) = r*PHSP1* ( rho51* PHSP2*( onemz4* TC_S5(soft5)  - DS_SO_TC(dsord) ) & 
                         - ( onemz4*TC_S5_C5(soft5c5) - DS_SO_TC_C5(dsordc5) ) )      ! again, multiply with r to get convegent 
                                                                                       ! integration 
   PTC_A_q_gg_qst(1:3,3,0:1) = faux1 - faux2 

   PTC_A_q_gg_qst(1:3,3,1) = PTC_A_q_gg_qst(1:3,3,1)   -( two*LogR + four*log(onemz) ) *PTC_A_q_gg_qst(1:3,3,0) & 
                                     -( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) ) &
                                     + pisq/three*faux3(1:3,0) 

   PTC_A_q_gg_qst(1:3,3,0:1) = PTC_A_q_gg_qst(1:3,3,0:1)/onemz

   PTC_A_q_gg_qst(1:3,2:3,0:1) = one/PHSP1*PTC_A_q_gg_qst(1:3,2:3,0:1)

! Raoul added
   PTC_A_q_gg_qst(1:3,1:3,1) = PTC_A_q_gg_qst(1:3,1:3,1) - PTC_A_q_gg_qst(1:3,1:3,0)*log(x4)
   

   end function 
   

! --- The triple collinear Altarelli-Parisi function; sector ``B '' -----------------------


   function PTC_B_q_gg_qst(z,zmin,yr,x3)
   real(dp) :: PTC_B_q_gg_qst(3,3,0:1)
   real(dp), intent(in) :: yr(3), z,zmin,x3
   real(dp) ::  rho45, full(6), fullc54(6),soft5(6), soft5c5(6),dsoft(6), dsord(6)
   real(dp) ::  dsoftc5(6), dsordc5(6), dsoftz(6), dsoftzc5(6)
   real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4, onemz 
   real(dp) :: faux1(1:3,0:1), faux2(1:3,0:1), LogR, LogF, LogF0, faux3(1:3,0:1), faux4(1:3,0:1),faclambda,kt1kt2sq


   r =  yr(1) 
   x4 = yr(2) 
   lambda = yr(3) 

   PTC_B_q_gg_qst = zero 

   call kinematics_tcl_B(z,(/r,x4,lambda,x3/),rho45,full,fullc54,soft5,soft5c5,dsoft,dsord,dsoftc5,dsordc5,dsoftz,dsoftzc5)

   PHSP1 = (one - half*r) 
   PHSP2 = one   ! everything is absorbed into rho45 

   faclambda = lambda * (one - lambda )
!   faclambda=one


    LogR =  log( half*r*(one-half*r))
    LogF =  log( 64.0_dp* (one-half*x4)/NF1(x3,one-0.5_dp*x4,lambda)**2*faclambda * (one-x3)*(one-x3*(one-x4*half)))
    LogF0 = log( 64.0_dp*one/NF1(x3,one,lambda)**2 * faclambda * (one-x3)**2 )
    kt1kt2sq = lambda * (one-x3)

!  the delta-function  

    faux1(1:3,0:1)= PHSP2*rho45*(r**2*DS_TC(dsoftz) - DS_SO_TC(dsord) )
    faux2(1:3,0:1)= r**2*DS_TC_C54(kt1kt2sq,dsoftzc5) - DS_SO_TC_C54(dsordc5)
    faux3(1:3,0:1)= PHSP2*rho45*( r**2 * DS_TC(dsoftz) )
    faux4(1:3,0:1)= r**2 * DS_TC_C54(kt1kt2sq,dsoftzc5)


   PTC_B_q_gg_qst(1:3,1,0:1) = -log(one-zmin)*(faux1 - faux2)  - log(one+r)* (faux3 - faux4)  

   PTC_B_q_gg_qst(1:3,1,1) = PTC_B_q_gg_qst(1:3,1,1) + two*log(one-zmin)**2*( faux1(1:3,0) - faux2(1:3,0) ) & 
                   -two*Log(r)*PTC_B_q_gg_qst(1:3,1,0) - two*log(one+r)**2 * (faux3(1:3,0) - faux4(1:3,0) ) & 
                   +log(one-zmin)*( LogF *faux1(1:3,0) - LogF0*faux2(1:3,0) ) &
                   +log(one+r)*( LogF*faux3(1:3,0) - LogF0*faux4(1:3,0) )     



!  the plus distribution coefficient

   faux1(1:3,0:1) = rho45*PHSP2*( PHSP1**2*(half*r)**2*DS_TC(dsoft) - DS_SO_TC(dsord) ) 
   faux2(1:3,0:1) = PHSP1**2*(half*r)**2*DS_TC_C54(kt1kt2sq,dsoftc5) - DS_SO_TC_C54(dsordc5) 
! Raoul modified to include ph sp
   faux3(1:3,0:1) = r*PHSP1*(rho45*PHSP2*DS_SO_TC(dsord) - DS_SO_TC_C54(dsordc5) )


    PTC_B_q_gg_qst(1:3,2,0:1) =  faux1 - faux2 

    PTC_B_q_gg_qst(1:3,2,1) = PTC_B_q_gg_qst(1:3,2,1) - two * LogR*PTC_B_q_gg_qst(1:3,2,0) & 
                -( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) ) & 
                + pisq/three*faux3(1:3,0) 

!  the regular function 

   onemz = one -z 
   onemz4 = onemz**4

   faux1(1:3,0:1)=rho45*PHSP2*(   onemz4 * ( PHSP1**2*(half*r)**2 *z*TC(full) - TC_S5(soft5) ) &
                             - PHSP1**2*(half*r)**2*DS_TC(dsoft) + DS_SO_TC(dsord)  )

   faux2(1:3,0:1) = onemz4 *( PHSP1**2*(half*r)**2*z*TC_C54(kt1kt2sq,fullc54) -TC_S5_C54(soft5c5) )  &
                - PHSP1**2*(half*r)**2*DS_TC_C54(kt1kt2sq,dsoftc5) + DS_SO_TC_C54(dsordc5) 
! Raoul modified to incl ph sp
   faux3(1:3,0:1) = r *PHSP1*( rho45*PHSP2* ( onemz4 * TC_S5(soft5) -  DS_SO_TC(dsord) ) & 
                         -  ( onemz4 * TC_S5_C54(soft5c5)  -  DS_SO_TC_C54(dsordc5)  ) )

   PTC_B_q_gg_qst(1:3,3,0:1) = faux1 - faux2 

   PTC_B_q_gg_qst(1:3,3,1) = PTC_B_q_gg_qst(1:3,3,1) - ( two*LogR + four*log(onemz) ) *PTC_B_q_gg_qst(1:3,3,0) &
                                   - ( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) ) &
                                   + pisq/three*faux3(1:3,0) 
     
   PTC_B_q_gg_qst(1:3,3,0:1) = one/onemz*PTC_B_q_gg_qst(1:3,3,0:1)

   PTC_B_q_gg_qst(1:3,2:3,0:1) = one/PHSP1*PTC_B_q_gg_qst(1:3,2:3,0:1)

! Raoul added                                                                                                                                                                                              
   PTC_B_q_gg_qst(1:3,1:3,1) = PTC_B_q_gg_qst(1:3,1:3,1) - two*PTC_B_q_gg_qst(1:3,1:3,0)*log(x4)


   end function 



! --- The triple collinear Altarelli-Parisi function; sector ``C'' -----------------------


   function PTC_C_q_gg_qst(z,zmin,yr,x3)
   real(dp) :: PTC_C_q_gg_qst(3,3,0:1)
   real(dp), intent(in) :: yr(3), z,zmin,x3
   real(dp) ::  rho41, full(6), fullc41(6),soft5(6), soft5c41(6),dsoft(6), dsord(6)
   real(dp)  :: dsoftc41(6), dsordc41(6), dsoftz(6), dsoftzc41(6)
   real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4, onemz, LogR, LogF, LogF0
   real(dp) ::  faux1(1:3,0:1), faux2(1:3,0:1), faux3(1:3,0:1), faux4(1:3,0:1), faclambda  


   r =  yr(1) 
   x4 = yr(2) 
   lambda = yr(3) 


   call kinematics_tcl_C(z,(/r,x4,lambda,x3/),rho41,full,fullc41,soft5,soft5c41,dsoft,dsord,dsoftc41,dsordc41,dsoftz,dsoftzc41)

   PTC_C_q_gg_qst = zero 

   PHSP1 = (one - half*r) 
   PHSP2 =  (one-half*x4)/NF1(x3,0.5_dp*x4,lambda)

   faclambda = lambda * (one - lambda )
!   faclambda=one


    LogR  = log( half*r*(one-half*r))
    LogF  = log( 128.0_dp * PHSP2**2 * faclambda * (one-x3) * (one-x3*x4*half) )
    LogF0 = log( 128.0_dp * faclambda * (one-x3)  )



!  the delta-function coefficient 

   faux1(1:3,0:1) =  PHSP2*rho41*( r**2 * DS_TC(dsoftz) - DS_SO_TC(dsord)  )
   faux2(1:3,0:1) =  r**2 * DS_TC_C4(dsoftzc41) - DS_SO_TC_C4(dsordc41) 
   faux3(1:3,0:1) = PHSP2*rho41* r**2 * DS_TC(dsoftz)
   faux4(1:3,0:1) = r**2 * DS_TC_C4(dsoftzc41) 

   PTC_C_q_gg_qst(1:3,1,0:1) = -log(one-zmin)*(faux1 - faux2)  - log(one+r)* (faux3 - faux4)  

   PTC_C_q_gg_qst(1:3,1,1) = PTC_C_q_gg_qst(1:3,1,1) + two*log(one-zmin)**2*( faux1(1:3,0) - faux2(1:3,0) ) & 
                   -two*Log(r)*PTC_C_q_gg_qst(1:3,1,0) - two*log(one+r)**2 * (faux3(1:3,0) - faux4(1:3,0) ) & 
                   +log(one-zmin)*( LogF *faux1(1:3,0) - LogF0*faux2(1:3,0) ) &
                   +log(one+r)*( LogF*faux3(1:3,0) - LogF0*faux4(1:3,0) )          


!  the plus distribution coefficient

   faux1(1:3,0:1) = rho41*PHSP2*( PHSP1**2*( half*r)**2 * DS_TC(dsoft) - DS_SO_TC(dsord) ) 
   faux2(1:3,0:1) = PHSP1**2*(half*r)**2 * DS_TC_C4(dsoftc41) - DS_SO_TC_C4(dsordc41)
! Raoul modified to include ph sp
   faux3(1:3,0:1) = r*PHSP1*(rho41*PHSP2*DS_SO_TC(dsord) - DS_SO_TC_C4(dsordc41) )

   PTC_C_q_gg_qst(1:3,2,0:1) = faux1 - faux2 

   PTC_C_q_gg_qst(1:3,2,1) = PTC_C_q_gg_qst(1:3,2,1) - two * LogR*PTC_C_q_gg_qst(1:3,2,0) & 
                                  -( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) ) & 
                                  + pisq/three*faux3(1:3,0) 


!  the regular function 

   onemz = one-z
   onemz4 = onemz**4

   faux1(1:3,0:1) = rho41*PHSP2*( onemz4*( z*PHSP1**2*(half*r)**2*TC(full)    -  TC_S5(soft5)  )   &
               - PHSP1**2*(half*r)**2*DS_TC(dsoft) + DS_SO_TC(dsord)  )

   faux2(1:3,0:1) = onemz4*( z*PHSP1**2*(half * r)**2*TC_C4(fullc41)  - TC_S5_C4(soft5c41) )  &
                       -PHSP1**2*(half*r)**2*DS_TC_C4(dsoftc41) +  DS_SO_TC_C4(dsordc41) 
! Raoul modified to incl ph sp
   faux3(1:3,0:1) = r*PHSP1*(  rho41*PHSP2*( onemz4*TC_S5(soft5) - DS_SO_TC(dsord)  ) & 
                    - ( onemz4*TC_S5_C4(soft5c41) - DS_SO_TC_C4(dsordc41) ) )


   PTC_C_q_gg_qst(1:3,3,0:1) = faux1 - faux2 


   PTC_C_q_gg_qst(1:3,3,1) =      PTC_C_q_gg_qst(1:3,3,1) - ( two*LogR + four*log(onemz) ) *PTC_C_q_gg_qst(1:3,3,0) &
                       - ( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) ) &
                       + pisq/three*faux3(1:3,0) 

   PTC_C_q_gg_qst(1:3,3,0:1) = PTC_C_q_gg_qst(1:3,3,0:1)/onemz

   PTC_C_q_gg_qst(1:3,2:3,0:1) = one/PHSP1*PTC_C_q_gg_qst(1:3,2:3,0:1)

! Raoul added                                                                                                                                                                                              
   PTC_C_q_gg_qst(1:3,1:3,1) = PTC_C_q_gg_qst(1:3,1:3,1) - PTC_C_q_gg_qst(1:3,1:3,0)*log(x4)


   end function 



! --- The triple collinear Altarelli-Parisi function; sector ``D '' -----------------------


   function PTC_D_q_gg_qst(z,zmin,yr,x3)
   real(dp) :: PTC_D_q_gg_qst(3,3,0:1)
   real(dp), intent(in) :: yr(3), z,zmin,x3
   real(dp) :: rho45, full(6), fullc54(6),soft5(6), soft5c5(6),dsoft(6), dsord(6)
   real(dp)  :: dsoftc5(6), dsordc5(6),dsoftz(6), dsoftzc5(6)
   real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4, onemz 
   real(dp) :: faux1(1:3,0:1), faux2(1:3,0:1), LogR, LogF, LogF0, faux3(1:3,0:1), faux4(1:3,0:1), faclambda,kt1kt2sq


   r =  yr(1) 
   x4 = yr(2) 
   lambda = yr(3) 

   PTC_D_q_gg_qst = zero 

   call kinematics_tcl_D(z,(/r,x4,lambda,x3/),rho45,full,fullc54,soft5,soft5c5,dsoft,dsord,dsoftc5,dsordc5,dsoftz,dsoftzc5)

   PHSP1 = (one - half*r) 
   PHSP2 = one   ! everything is absorbed into rho45 

   faclambda = lambda *( one -  lambda ) 
!   faclambda=one


    LogR =  log( half*r*(one-half*r))
    LogF =  log( 64.0_dp* (one-half*x4)/NF1(x3,one-0.5_dp*x4,lambda)**2*faclambda * (one-x3)*(one-x3*(one-x4*half)))
    LogF0 = log( 64.0_dp*one/NF1(x3,one,lambda)**2 * faclambda * (one-x3)**2 )
    kt1kt2sq = lambda * (one-x3)

   

!  the delta-function  

    faux1(1:3,0:1)= PHSP2*rho45*(r**2*DS_TC(dsoftz) - DS_SO_TC(dsord) )
    faux2(1:3,0:1)= r**2*DS_TC_C54(kt1kt2sq,dsoftzc5 )- DS_SO_TC_C54(dsordc5)
    faux3(1:3,0:1)= PHSP2*rho45*( r**2 * DS_TC(dsoftz) )
    faux4(1:3,0:1)= r**2 * DS_TC_C54(kt1kt2sq,dsoftzc5)


   PTC_D_q_gg_qst(1:3,1,0:1) = -log(one-zmin)*(faux1 - faux2)  - log(one+r)* (faux3 - faux4)  

   PTC_D_q_gg_qst(1:3,1,1) = PTC_D_q_gg_qst(1:3,1,1) + two*log(one-zmin)**2*( faux1(1:3,0) - faux2(1:3,0) ) & 
                   -two*Log(r)*PTC_D_q_gg_qst(1:3,1,0) - two*log(one+r)**2 * (faux3(1:3,0) - faux4(1:3,0) ) & 
                   +log(one-zmin)*( LogF *faux1(1:3,0) - LogF0*faux2(1:3,0) ) &
                   +log(one+r)*( LogF*faux3(1:3,0) - LogF0*faux4(1:3,0) )     



!  the plus distribution coefficient

   faux1(1:3,0:1) = rho45*PHSP2*( PHSP1**2*(half*r)**2*DS_TC(dsoft) - DS_SO_TC(dsord) ) 
   faux2(1:3,0:1) = PHSP1**2*(half*r)**2*DS_TC_C54(kt1kt2sq,dsoftc5) - DS_SO_TC_C54(dsordc5) 
! Raoul modified to include ph sp
   faux3(1:3,0:1) = r *PHSP1* (rho45*PHSP2*DS_SO_TC(dsord) - DS_SO_TC_C54(dsordc5) )


    PTC_D_q_gg_qst(1:3,2,0:1) =  faux1 - faux2 

    PTC_D_q_gg_qst(1:3,2,1) = PTC_D_q_gg_qst(1:3,2,1) - two * LogR*PTC_D_q_gg_qst(1:3,2,0) &
                 -( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) ) & 
                 + pisq/three*faux3(1:3,0) 

!  the regular function 

   onemz = one -z 
   onemz4 = onemz**4

   faux1(1:3,0:1)=rho45*PHSP2*(   onemz4 * ( PHSP1**2*(half*r)**2 *z*TC(full) - TC_S5(soft5) ) &
                             - PHSP1**2*(half*r)**2*DS_TC(dsoft) + DS_SO_TC(dsord)  )
 
    faux2(1:3,0:1) = onemz4 *( PHSP1**2*(half*r)**2*z*TC_C54(kt1kt2sq,fullc54) -TC_S5_C54(soft5c5) )  &
                - PHSP1**2*(half*r)**2*DS_TC_C54(kt1kt2sq,dsoftc5) + DS_SO_TC_C54(dsordc5) 

!Raoul modified to incl ph sp
   faux3(1:3,0:1) = r *PHSP1*( rho45*PHSP2* ( onemz4 * TC_S5(soft5) -  DS_SO_TC(dsord) ) & 
                         -  ( onemz4 * TC_S5_C54(soft5c5)  -  DS_SO_TC_C54(dsordc5)  ) )

   PTC_D_q_gg_qst(1:3,3,0:1) = faux1 - faux2 

   PTC_D_q_gg_qst(1:3,3,1) = PTC_D_q_gg_qst(1:3,3,1) - ( two*LogR + four*log(onemz) ) *PTC_D_q_gg_qst(1:3,3,0) &
                        -( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) ) &
                        + pisq/three*faux3(1:3,0) 
     
   PTC_D_q_gg_qst(1:3,3,0:1) = one/onemz*PTC_D_q_gg_qst(1:3,3,0:1)

   PTC_D_q_gg_qst(1:3,2:3,0:1) = one/PHSP1*PTC_D_q_gg_qst(1:3,2:3,0:1)

! Raoul added                                                                                                                                                                                              
   PTC_D_q_gg_qst(1:3,1:3,1) = PTC_D_q_gg_qst(1:3,1:3,1) - two*PTC_D_q_gg_qst(1:3,1:3,0)*log(x4)


    end function 







  end module mod_tripleclimits_q_gg_qstar

module mod_tripleclimits_q_qg_gstar
  use mod_types; use mod_consts_dp; use common_def
  use mod_kinematics_tc_unord
  use mod_auxfunctions
  implicit none
  real(dp), parameter :: flipcol = -Cf/tr !-- when going from initial to final state
  private
  
  public :: PTC_A_gqqb  , PTC_B_gqqb , PTC_C_gqqb , PTC_D_gqqb
  public :: tripcoll_subtr_gqqb


contains
  
  function tripcoll_subtr_gqqb(ndim,yRnd,ncomp,ff,userdata,weight,iternumber)
    !--
    integer, parameter :: ep = 1
    integer, parameter :: mycol = 3 !-- ab -> 1, nab -> 2, tot -> 3
    real(dp), parameter :: zmin = zero
    real(dp), parameter :: myz = 0.123_dp
    !--
    real(dp), parameter :: buff = 1E-10_dp
    real(dp), parameter :: onet = one - two * buff
    real(dp), parameter :: bufftc = 1E-10_dp
    !--
    real(dp15) :: yRnd(30), ff(1), weight
    integer :: ndim,ncomp,userdata,iternumber,tripcoll_subtr_gqqb
    real(dp)              :: x1,x2,x3,xx2,xx(1:3),x3a,xx3
    real(dp)              :: res_tc(1:3,1:3,0:1)
    real(dp)              :: seca(1:3,1:3,0:1),secb(1:3,1:3,0:1)
    real(dp)              :: secc(1:3,1:3,0:1),secd(1:3,1:3,0:1)
    real(dp)              :: seca0(1:3,1:3,0:1),secb0(1:3,1:3,0:1)
    real(dp)              :: secc0(1:3,1:3,0:1),secd0(1:3,1:3,0:1)
    real(dp)              :: wghtac,wghtbd,z

    tripcoll_subtr_gqqb = zero
    z = myz
    ff(1) = zero
    res_tc = zero 

    seca = zero; seca0 = zero
    secb = zero; secb0 = zero
    secc = zero; secc0 = zero
    secd = zero; secd0 = zero

    x1 = buff+onet*real(yRnd(1),dp)  ! r
    x2 = buff+onet*real(yRnd(2),dp)  ! x4
    x3a = buff+onet*real(yRnd(3),dp) ! lambda
    xx3 = buff+onet*real(yRnd(4),dp)  ! x3

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

    xx = (/x1,x2,x3/)

    seca0 = PTC_A_gqqb(z,xx,zero)
    secc0 = PTC_C_gqqb(z,xx,zero)

    if (withmax) then 
       seca = PTC_A_gqqb(z,xx,xx3)
       secc = PTC_C_gqqb(z,xx,xx3)
    endif

    wghtac = -one/x1/x2

    !-- 

    xx2 = x2**2
    xx = (/x1,xx2,x3/)

    secb0 = PTC_B_gqqb(z,xx,zero)
    secd0 = PTC_D_gqqb(z,xx,zero)

    if (withmax) then
       secb = PTC_B_gqqb(z,xx,xx3)
       secd = PTC_D_gqqb(z,xx,xx3)
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

    !-- 
    
    ff(1) = res_tc(mycol,3,ep)

    if (ff(1) .ne. ff(1) .or. weight .ne. weight) then
       print *, "xx",xx,x3a
       print *, "wghtac",wghtac
       print *, "sec a",seca(mycol,3,ep)
       print *, "sec c",secc(mycol,3,ep)
       print *, "wghtbd",wghtac
       print *, "sec b",secb(mycol,3,ep)
       print *, "sec d",secd(mycol,3,ep)

       print *, "ff",ff(1)
       stop
    endif

    ! print *, 'buff',buff,bufftc
    ! print *, 'seca', seca(3,3,:)
    ! print *, 'secb', secb(3,3,:)
    ! print *, 'secc', secc(3,3,:)
    ! print *, 'secd', secd(3,3,:)
    
    ! print *, 'ff(1)',ff(1)
    ! stop

    return

  end function tripcoll_subtr_gqqb

  !--

  function PTC_A_gqqb(z,yr,x3)
    real(dp) :: PTC_A_gqqb(1:3,1:3,0:1)
    real(dp), intent(in) :: yr(3),z,x3
    real(dp) ::  rho51, full(6), fullc51(6),soft5(6), soft5c5(6),dsoft(6), dsord(6), dsoftz(6), dsoftzc5(6)
    real(dp) :: dsoft5(6), dsordc5(6), dsoftc5(6)
    real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4 , PHSP1m1, energy, onemz
    real(dp) :: faux1(1:3,0:1), faux2(1:3,0:1), LogR, LogF, LogF0, faux3(1:3,0:1), faux4(1:3,0:1), faclambda 
    real(dp) :: test1(1:3,0:1), test2(1:3,0:1), test3(1:3,0:1), test4(1:3,0:1)  
    real(dp) :: kt1kt2sq


    r =  yr(1)
    x4 = yr(2)
    lambda = yr(3) 

    PTC_A_gqqb = zero 

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

    PTC_A_gqqb(1:3,3,0:1) = faux1 - faux2 

    PTC_A_gqqb(1:3,3,0) = PTC_A_gqqb(1:3,3,0) + ( faux3(1:3,0) - faux4(1:3,0) ) * ( log(onemz) - one)


    PTC_A_gqqb(1:3,3,1) = PTC_A_gqqb(1:3,3,1)   -  two*LogR*(faux1(1:3,0)-faux2(1:3,0)) - four*log(onemz)  *PTC_A_gqqb(1:3,3,0)  & 
         - ( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) )  &
         + ( faux3(1:3,0) - faux4(1:3,0) ) *( log(onemz)**2 - four + pisq/three) & 
         - ( LogF*faux3(1:3,0) - LogF0*faux4(1:3,0) )*(log(onemz)-one)  &
         + ( faux3(1:3,1) - faux4(1:3,1) ) * ( log(onemz) - one)

    PTC_A_gqqb(1:3,1:3,1) = PTC_A_gqqb(1:3,1:3,1) - PTC_A_gqqb(1:3,1:3,0)*log(x4)

    ! Normalization 
    PTC_A_gqqb = (one-z)*z * PTC_A_gqqb 

  end function PTC_A_gqqb




  function PTC_B_gqqb(z,yr,x3)
    real(dp) :: PTC_B_gqqb(3,3,0:1)
    real(dp), intent(in) :: yr(3),z,x3
    real(dp) ::  rho45, full(6), fullc54(6),soft5(6), soft5c5(6),dsoft(6), dsord(6), dsoftz(6), dsoftzc5(6)
    real(dp) :: dsoft5(6), dsordc5(6), dsoftc5(6)
    real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4 , PHSP1m1, energy, onemz
    real(dp) :: faux1(1:3,0:1), faux2(1:3,0:1), LogR, LogF, LogF0, faux3(1:3,0:1), faux4(1:3,0:1), faclambda 
    real(dp) :: test1(1:3,0:1), test2(1:3,0:1), test3(1:3,0:1), test4(1:3,0:1)  

    r =  yr(1) 
    x4 = yr(2) 
    lambda = yr(3) 

    PTC_B_gqqb = zero 

    !-- debug
    !r = 0.74792768547143218_dp
    !x4 = 0.36739089737475572_dp
    !lambda = 0.48063689875473148_dp
    !print *, 'fixing r'
    !-- end debug

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


    PTC_B_gqqb(1:3,3,0:1) = faux1 - faux2 
    PTC_B_gqqb(1:3,3,0) = PTC_B_gqqb(1:3,3,0) + ( faux3(1:3,0) - faux4(1:3,0) ) * ( log(onemz) - one)

    PTC_B_gqqb(1:3,3,1) = PTC_B_gqqb(1:3,3,1)   -  two*LogR*(faux1(1:3,0)-faux2(1:3,0)) - four*log(onemz)  *PTC_B_gqqb(1:3,3,0)  & 
         - ( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) )  &
         + ( faux3(1:3,0) - faux4(1:3,0) ) *( log(onemz)**2 - four + pisq/three) & 
         - ( LogF*faux3(1:3,0) - LogF0*faux4(1:3,0) )*(log(onemz)-one)  &
         + ( faux3(1:3,1) - faux4(1:3,1) ) * ( log(onemz) - one)
    
    PTC_B_gqqb(1:3,1:3,1) = PTC_B_gqqb(1:3,1:3,1) - two*PTC_B_gqqb(1:3,1:3,0)*log(x4)


    ! Normalization 

    PTC_B_gqqb = (one-z)*z * PTC_B_gqqb 

  end function PTC_B_gqqb



  function PTC_C_gqqb(z,yr,x3)
    real(dp) :: PTC_C_gqqb(1:3,1:3,0:1)
    real(dp), intent(in) :: yr(3),z,x3
    real(dp) ::  rho41, full(6), fullc41(6),soft5(6), soft5c41(6)
    real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4 , PHSP1m1, energy, onemz
    real(dp) :: faux1(1:3,0:1), faux2(1:3,0:1), LogR, LogF, LogF0, faux3(1:3,0:1), faux4(1:3,0:1), faclambda 
    real(dp) :: test1(1:3,0:1), test2(1:3,0:1), test3(1:3,0:1), test4(1:3,0:1) 
    real(dp) :: kt1kt2sq

    r =  yr(1) 
    x4 = yr(2) 
    lambda = yr(3) 

    PTC_C_gqqb = zero 

    call kinematics_tcl_unord_C(z,(/r,x4,lambda,x3/),rho41,full,fullc41,soft5,soft5c41)

    PHSP1 = (one - r) 
    PHSP2 =  (one-half*x4)/NF1(x3,0.5_dp*x4,lambda)

    faclambda = lambda*(one-lambda) 

    onemz = one-z


    LogR  = log( r*(one-r))
    LogF  = log( 128.0_dp * PHSP2**2 * faclambda * (one-x3) * (one-x3*x4*half) )
    LogF0 = log( 128.0_dp * faclambda * (one-x3)  )

    kt1kt2sq = (two*lambda-one)**2 * (one-x3)

    !  the regular piece 


    faux1(1:3,0:1) =  rho41 * PHSP2 * PHSP1* onemz**2 * ( r**2 * TC(full)  - TC_S5(soft5)  )

    faux2(1:3,0:1) =                  PHSP1* onemz**2 * ( r**2 * TC_C4(kt1kt2sq,fullc41)  - TC_S5_C4(soft5c41)  )

    faux3(1:3,0:1) =   r*onemz**2 *  rho41 * PHSP2 * TC_S5(soft5)  

    faux4(1:3,0:1) =   r*onemz**2 * TC_S5_C4(soft5c41)  

    PTC_C_gqqb(1:3,3,0:1) = faux1 - faux2 

    PTC_C_gqqb(1:3,3,0) = PTC_C_gqqb(1:3,3,0) + ( faux3(1:3,0) - faux4(1:3,0) ) * ( log(onemz) - one)


    PTC_C_gqqb(1:3,3,1) = PTC_C_gqqb(1:3,3,1)   -  two*LogR*(faux1(1:3,0)-faux2(1:3,0)) - four*log(onemz)  *PTC_C_gqqb(1:3,3,0)  & 
         - ( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) )  &
         + ( faux3(1:3,0) - faux4(1:3,0) ) *( log(onemz)**2 - four + pisq/three) & 
         - ( LogF*faux3(1:3,0) - LogF0*faux4(1:3,0) )*(log(onemz)-one)  &
         + ( faux3(1:3,1) - faux4(1:3,1) ) * ( log(onemz) - one)

    PTC_C_gqqb(1:3,1:3,1) = PTC_C_gqqb(1:3,1:3,1) - PTC_C_gqqb(1:3,1:3,0)*log(x4)

    ! Normalization 

    PTC_C_gqqb = (one-z)*z * PTC_C_gqqb 

  end function PTC_C_gqqb


  function PTC_D_gqqb(z,yr,x3)
    real(dp) :: PTC_D_gqqb(1:3,1:3,0:1)
    real(dp), intent(in) :: yr(3),z,x3
    real(dp) ::  rho45, full(6), fullc54(6),soft5(6), soft5c5(6)
    real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4 , PHSP1m1, energy, onemz
    real(dp) :: faux1(1:3,0:1), faux2(1:3,0:1), LogR, LogF, LogF0, faux3(1:3,0:1), faux4(1:3,0:1), faclambda 
    real(dp) :: test1(1:3,0:1), test2(1:3,0:1), test3(1:3,0:1), test4(1:3,0:1)  

    r =  yr(1) 
    x4 = yr(2) 
    lambda = yr(3) 

    PTC_D_gqqb = zero 

    !-- debug
    !r = 0.74792768547143218_dp
    !x4 = 0.36739089737475572_dp
    !lambda = 0.48063689875473148_dp
    !print *, 'fixing r'
    !-- end debug

    call kinematics_tcl_unord_D(z,(/r,x4,lambda,x3/),rho45,full,fullc54,soft5,soft5c5)

    PHSP1 = (one - r) 
    PHSP2 =  one   ! evrything is absorbed into rho45 

    faclambda = lambda*(one-lambda) 

    LogR =  log( r*(one-r))
    LogF =  log( 64.0_dp* (one-half*x4)/NF1(x3,one-0.5_dp*x4,lambda)**2*faclambda * (one-x3)*(one-x3*(one-x4*half)))
    LogF0 = log( 64.0_dp*one/NF1(x3,one,lambda)**2 * faclambda * (one-x3)**2 )

    ! print *, 'logR,logF',logR,logF,logF0

    !  only the regular piece is there 
    onemz = one-z


    faux1(1:3,0:1) =  rho45 * PHSP2 * PHSP1* onemz**2 * ( r**2 * TC(full)  - TC_S5(soft5)  )
    faux2(1:3,0:1) =                  PHSP1* onemz**2 * ( r**2 * TC_C54(fullc54)  - TC_S5_C54(soft5c5)  )
    faux3(1:3,0:1) =   r * onemz**2 *  rho45 * PHSP2 * TC_S5(soft5)  
    faux4(1:3,0:1) =   r * onemz**2 * TC_S5_C54(soft5c5)  

    ! print *, 'faux1(3,:)',faux1(3,:)
    ! print *, 'faux2(3,:)',faux2(3,:)
    ! print *, 'faux3(3,:)',faux3(3,:)
    ! print *, 'faux4(3,:)',faux4(3,:)

    PTC_D_gqqb(1:3,3,0:1) = faux1 - faux2 
    PTC_D_gqqb(1:3,3,0) = PTC_D_gqqb(1:3,3,0) + ( faux3(1:3,0) - faux4(1:3,0) ) * ( log(onemz) - one)



    PTC_D_gqqb(1:3,3,1) = PTC_D_gqqb(1:3,3,1)   -  two*LogR*(faux1(1:3,0)-faux2(1:3,0)) - four*log(onemz)  *PTC_D_gqqb(1:3,3,0)  & 
         - ( LogF*faux1(1:3,0) - LogF0*faux2(1:3,0) )  &
         + ( faux3(1:3,0) - faux4(1:3,0) ) *( log(onemz)**2 - four + pisq/three) & 
         - ( LogF*faux3(1:3,0) - LogF0*faux4(1:3,0) )*(log(onemz)-one)  &
         + ( faux3(1:3,1) - faux4(1:3,1) ) * ( log(onemz) - one)


    PTC_D_gqqb(1:3,1:3,1) = PTC_D_gqqb(1:3,1:3,1) - two*PTC_D_gqqb(1:3,1:3,0)*log(x4)


    ! Normalization 
    PTC_D_gqqb = (one-z)*z * PTC_D_gqqb 


  end function PTC_D_gqqb

  !--

  function TC(kinlist) 
    real(dp) :: TC(3,0:1) 
    real(dp), intent(in) :: kinlist(6)
    real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3 
    real(dp) :: split(0:2)

    s14 = kinlist(1) 
    s15 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 

    TC = zero 

    s13 = -s15
    s23 = -s14 
    s12 = s45
    s123 = s12 + s13 + s23 

    z1 = z5 
    z2 = z4 
    z3 = one - z4 - z5

    ! -- abelian 

    split = Pgqqb_ab_spav(s12,s13,s23,s123,z1,z2,z3)
    TC(1,:) = Cf*Tr*split(0:1)

    !-- non-abelian 

    split = Pgqqb_nab_spav_symm(s12,s13,s23,s123,z1,z2,z3)
    TC(2,:) = Ca*Tr*split(0:1)

    !-- full 

    TC(3,0:1) = TC(1,0:1) + TC(2,0:1) 

    TC = four/s123**2*TC * flipcol

  end function TC

  function TC_C5(kt1kt2sq,kinlist) 
    real(dp) :: TC_C5(3,0:1) 
    real(dp), intent(in) :: kt1kt2sq,kinlist(6)
    real(dp) :: s14,s15,s45,z4,z5,z1,zz1,zz2,sab,sab_c
    real(dp) :: Pqq,Pqq1,Pgq_spav,Pgq_spav1,z
    real(dp) :: kperpprod

    Pqq(z) = Cf * (one + z**2)/(one-z)
    Pqq1(z) = Cf * (-(one-z))

    Pgq_spav(z) = tr * (one-two*z*(one-z))
    Pgq_spav1(z) = -two * tr * z * (one-z)

    s14 = kinlist(1) 
    s15 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 
    z1 = kinlist(6)

    zz1 = z1/(z1+z5)
    zz2 = (z1+z5)/(z1+z4+z5)

    sab = s15
    sab_c = s14 * (z1+z5)/z1

    TC_C5 = zero 

    TC_C5(1,0) = Pgq_spav(zz2)*Pqq(zz1)*two/sab*two/sab_c
    TC_C5(1,1) = (Pgq_spav1(zz2)*Pqq(zz1)+Pgq_spav(zz2)*Pqq1(zz1))*two/sab*two/sab_c

    TC_C5(3,0:1) = (TC_C5(1,0:1) + TC_C5(2,0:1)) * flipcol

  end function TC_C5

  function TC_C54(kinlist) 
    real(dp) :: TC_C54(3,0:1) 
    real(dp), intent(in) :: kinlist(6)
    real(dp) :: s14,s15,s45,z4,z5,z1,zz1,zz2,sab,sab_c
    real(dp) :: Pqq,Pqq1,Pgq_spav,Pgq_spav1,z
    real(dp) :: kperpprod

    Pqq(z) = Cf * (one + z**2)/(one-z)
    Pqq1(z) = Cf * (-(one-z))

    Pgq_spav(z) = tr * (one-two*z*(one-z))
    Pgq_spav1(z) = -two * tr * z * (one-z)

    s14 = kinlist(1) 
    s15 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 
    z1 = kinlist(6)

    zz1 = z4/(z4+z5)
    zz2 = (z4+z5)/(z1+z4+z5)

    sab = s45
    sab_c = -s14 * (z4+z5)/z4

    TC_C54 = zero 

    TC_C54(1,0) = Pgq_spav(zz2)*Pqq(zz1)*two/sab*two/sab_c
    TC_C54(1,1) = (Pgq_spav1(zz2)*Pqq(zz1)+Pgq_spav(zz2)*Pqq1(zz1))*two/sab*two/sab_c

    TC_C54(3,0:1) = (TC_C54(1,0:1) + TC_C54(2,0:1)) * flipcol

  end function TC_C54

  function TC_C4(kt1kt2sq,kinlist) 
    real(dp) :: TC_C4(3,0:1) 
    real(dp), intent(in) :: kt1kt2sq,kinlist(6)
    real(dp) :: s14,s15,s45,z4,z5,z1,zz1,zz2,sab,sab_c
    real(dp) :: Pgg_spav,Pgq_spav,Pgq_spav1,z
    real(dp) :: kperpprod

    Pgg_spav(z) = two * Ca * ((one-z)/z + z/(one-z) + z*(one-z))

    Pgq_spav(z) = tr * (one-two*z*(one-z))
    Pgq_spav1(z) = -two * tr * z * (one-z)

    s14 = kinlist(1) 
    s15 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 
    z1 = kinlist(6)

    zz1 = z4/(z1+z4)
    zz2 = (z1+z4)/(z1+z4+z5)

    sab = -s14
    sab_c = s45 * (z1+z4)/z4

    TC_C4 = zero 

    TC_C4(2,0) = ( Pgg_spav(zz2)*Pgq_spav(zz1) &
         - 8.0_dp * Ca * tr * zz1*(one-zz1)*(one-zz2)/zz2 * (kt1kt2sq-half) ) * two/sab*two/sab_c 

    TC_C4(2,1) = ( Pgg_spav(zz2)*Pgq_spav1(zz1) &
         - 8.0_dp * Ca * tr * zz1*(one-zz1)*(one-zz2)/zz2 * (-half) ) * two/sab*two/sab_c 

    TC_C4(3,0:1) = (TC_C4(1,0:1) + TC_C4(2,0:1)) * flipcol

  end function TC_C4

  !--
  
  function TC_S5(kinlist) 
    real(dp) :: TC_S5(3,0:1) 
    real(dp), intent(in) :: kinlist(6)
    real(dp) :: s24,s25,s45,z4,z5,z1
    real(dp) :: Pgq_spav,Pgq_spav1,z
    real(dp) :: zz1,zz2,zzi,p1p2,pip2,pip1

    Pgq_spav(z) = tr * (one-two*z*(one-z))
    Pgq_spav1(z) = -two * tr * z * (one-z)

    s24 = kinlist(1) 
    s25 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 
    z1 = kinlist(6)

    !--

    z = z1
    zz1 = z1
    zz2 = z5
    zzi = z4

    p1p2 = -half*s25
    pip2 = half*s45
    pip1 = -half*s24

    TC_S5 = 0

    ! -- ((two*Cf-Ca)/p1p2/pip2+Ca*zz1/zz2/p1p2/pip1+Ca*zzi/zz2/pip1/pip2)*Pgq_spav(z)

    TC_S5(1,0) = two*Cf/p1p2/pip2*Pgq_spav(z)
    TC_S5(1,1) = two*Cf/p1p2/pip2*Pgq_spav1(z)

    TC_S5(2,0) = Ca*(-one/p1p2/pip2+zz1/zz2/p1p2/pip1+zzi/zz2/pip1/pip2)*Pgq_spav(z)
    TC_S5(2,1) = Ca*(-one/p1p2/pip2+zz1/zz2/p1p2/pip1+zzi/zz2/pip1/pip2)*Pgq_spav1(z)

    TC_S5(3,0:1) = (TC_S5(1,0:1) + TC_S5(2,0:1))*flipcol

  end function TC_S5

  function TC_S5_C5(kinlist) 
    real(dp) :: TC_S5_C5(3,0:1) 
    real(dp), intent(in) :: kinlist(6)
    real(dp) :: s24,s25,s45,z4,z5,z1
    real(dp) :: Pgq_spav,Pgq_spav1,z
    real(dp) :: zz1,zz2,zzi,p1p2,pip2,pip1

    Pgq_spav(z) = tr * (one-two*z*(one-z))
    Pgq_spav1(z) = -two * tr * z * (one-z)

    s24 = kinlist(1) 
    s25 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 
    z1 = kinlist(6)

    !--

    z = z1
    zz1 = z1
    zz2 = z5
    zzi = z4

    p1p2 = -half*s25
    pip2 = half*s45
    pip1 = -half*s24

    TC_S5_C5 = 0

    TC_S5_C5(1,0) = two*Cf/p1p2/pip2*Pgq_spav(z)
    TC_S5_C5(1,1) = two*Cf/p1p2/pip2*Pgq_spav1(z)

    TC_S5_C5(2,0) = Ca*(-one/p1p2/pip2+zz1/zz2/p1p2/pip1)*Pgq_spav(z)
    TC_S5_C5(2,1) = Ca*(-one/p1p2/pip2+zz1/zz2/p1p2/pip1)*Pgq_spav1(z)

    TC_S5_C5(3,0:1) = (TC_S5_C5(1,0:1) + TC_S5_C5(2,0:1))*flipcol

  end function TC_S5_C5

  function TC_S5_C54(kinlist) 
    real(dp) :: TC_S5_C54(3,0:1) 
    real(dp), intent(in) :: kinlist(6)
    real(dp) :: s24,s25,s45,z4,z5,z1
    real(dp) :: Pgq_spav,Pgq_spav1,z
    real(dp) :: zz1,zz2,zzi,p1p2,pip2,pip1

    Pgq_spav(z) = tr * (one-two*z*(one-z))
    Pgq_spav1(z) = -two * tr * z * (one-z)

    s24 = kinlist(1) 
    s25 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 
    z1 = kinlist(6)

    !--

    z = z1
    zz1 = z1
    zz2 = z5
    zzi = z4

    p1p2 = -half*s25
    pip2 = half*s45
    pip1 = -half*s24

    TC_S5_C54 = 0

    ! -- ((two*Cf-Ca)/p1p2/pip2+Ca*zz1/zz2/p1p2/pip1+Ca*zzi/zz2/pip1/pip2)*Pgq_spav(z)

    TC_S5_C54(1,0) = two*Cf/p1p2/pip2*Pgq_spav(z)
    TC_S5_C54(1,1) = two*Cf/p1p2/pip2*Pgq_spav1(z)

    TC_S5_C54(2,0) = Ca*(-one/p1p2/pip2+zzi/zz2/pip1/pip2)*Pgq_spav(z)
    TC_S5_C54(2,1) = Ca*(-one/p1p2/pip2+zzi/zz2/pip1/pip2)*Pgq_spav1(z)

    TC_S5_C54(3,0:1) = (TC_S5_C54(1,0:1) + TC_S5_C54(2,0:1))*flipcol

  end function TC_S5_C54

  function TC_S5_C4(kinlist) 
    real(dp) :: TC_S5_C4(3,0:1) 
    real(dp), intent(in) :: kinlist(6)
    real(dp) :: s24,s25,s45,z4,z5,z1
    real(dp) :: Pgq_spav,Pgq_spav1,z
    real(dp) :: zz1,zz2,zzi,p1p2,pip2,pip1

    Pgq_spav(z) = tr * (one-two*z*(one-z))
    Pgq_spav1(z) = -two * tr * z * (one-z)

    s24 = kinlist(1) 
    s25 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 
    z1 = kinlist(6)

    !--

    z = z1
    zz1 = z1
    zz2 = z5
    zzi = z4

    p1p2 = -half*s25
    pip2 = half*s45
    pip1 = -half*s24

    TC_S5_C4 = 0

    TC_S5_C4(2,0) = Ca*(zz1/zz2/p1p2/pip1+zzi/zz2/pip1/pip2)*Pgq_spav(z)
    TC_S5_C4(2,1) = Ca*(zz1/zz2/p1p2/pip1+zzi/zz2/pip1/pip2)*Pgq_spav1(z)

    TC_S5_C4(3,0:1) = (TC_S5_C4(1,0:1) + TC_S5_C4(2,0:1))*flipcol

  end function TC_S5_C4

  !--

    !-- different quarks
  function PqbqR(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: PqbqR(0:2)
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3
    real(dp) :: t123,bit

    t123 = tijk(s12,s13,s23,z1,z2,z3)

    bit = z1+z2-s12/s123

    PqbqR(0) = -t123**2/s12/s123 + (four*z3+(z1-z2)**2)/(z1+z2) + bit
    PqbqR(1) = -two*bit
    PqbqR(2) = 0
    
    PqbqR(:) = PqbqR(:) * half * Cf * tr * s123/s12

    return

  end function PqbqR

  !-- identical bit, not symmetrized
  function Pqbqq_id(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pqbqq_id(0:2)
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3

    Pqbqq_id(0) = two*s23/s12 &
         +s123/s12*( (one+z1**2)/(one-z2)-two*z2/(one-z3) ) &
         -s123**2/s12/s13*z1/two*( (one+z1**2)/(one-z2)/(one-z3) )

    Pqbqq_id(1) = -one -two*s23/s12 &
         -s123/s12*( (one-z3)**2/(one-z2) + one+z1 - two*z2/(one-z3) ) &
         +s123**2/s12/s13*z1/two*( one+two*(one-z2)/(one-z3))

    Pqbqq_id(2) = one &
         -s123/s12 * (one-z3) &
         +s123**2/s12/s13*z1/two

    Pqbqq_id(:) = Pqbqq_id(:) * Cf * (Cf-half*Ca)

    return

  end function Pqbqq_id

  function Pqbqq_id_symm(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pqbqq_id_symm(0:2)
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3
    real(dp) :: p1(0:2),p2(0:2)

    p1 = Pqbqq_id(s12,s13,s23,s123,z1,z2,z3)
    p2 = Pqbqq_id(s13,s12,s23,s123,z1,z3,z2)

    Pqbqq_id_symm = p1+p2
    
    return

  end function Pqbqq_id_symm


  function Pggq_ab(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pggq_ab(0:2)
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3

    Pggq_ab(0) = s123**2/two/s13/s23*z3*(one+z3**2)/z1/z2 &
         +s123/s13* (z3*(one-z1)+(one-z2)**3)/z1/z2 &
         - s23/s13

    Pggq_ab(1) = -s123**2/two/s13/s23*z3*( (z1**2+z2**2)/z1/z2 + one ) &
         -s123/s13 * (z1**2 + z1*z2 + z2**2)*(one-z2)/z1/z2 &
         +one + two * s23/s13

    Pggq_ab(2) = -s123**2/two/s13/s23*z3 &
         +s123/s13 * (one+z3) &
         -one - s23/s13

    return

  end function Pggq_ab

  function Pggq_nab(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pggq_nab(0:2)
    real(dp) :: s12,s13,s23,s123,z1,z2,z3,t123
    real(dp) :: bit

    t123 = tijk(s12,s13,s23,z1,z2,z3)

    bit = (z1*(two-two*z1+z1**2)-z2*(six-six*z2+z2**2))/z2/(one-z3)

    Pggq_nab(0) = t123**2/four/s12**2 + one/four + &
         s123**2/two/s12/s13 * ( (one+z3**2)/z2 &
         +(z2**2+two*(one-z2))/(one-z3) ) &
         -s123**2/four/s13/s23*z3 * ( (one+z3**2)/z1/z2 ) &
         +s123/two/s12 * ( bit ) &
         +s123/two/s13 * ( ((one-z2)**3+z3**2-z2)/z2/(one-z3) &
         -(z3*(one-z1)+(one-z2)**3)/z1/z2 )

    Pggq_nab(1) = -t123**2/four/s12**2 - one/four - one/two &
         -s123**2/two/s12/s13 * ( (one-z3)**2/z2 + z2**2/(one-z3) ) &
         +s123**2/four/s13/s23*z3 * ( (one-z3)**2/z1/z2 - one ) &
         +s123/two/s12 * ( -bit + two * (z3*(z1-two*z2) - z2)/z2/(one-z3) ) &
         +s123/two/s13 * ( -((one-z2)**3+z3**2-z2)/z2/(one-z3) &
         -(two*(one-z2)*(z2-z3)/z2/(one-z3) - z1 + z2) &
         + (one-z2) * ((z1**2+z2**2)/z1/z2) )

    Pggq_nab(2) = +half + s123**2/four/s13/s23*z3 - s123/two/s13 * (one-z2)

    return

  end function Pggq_nab

  function Pgqqb_ab_spav(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pgqqb_ab_spav(0:2)
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3
    real(dp) :: overomep

    !-- piece proportional to 1/(1-ep)
    overomep = -two * s123**2/s12/s13 * (z1+two*z2*z3) &
         + s123/s12*two*(z1+z2) + s123/s13*two*(z1+z3)

    Pgqqb_ab_spav(0) = -two - s23*(one/s12+one/s13) + two*s123**2/s12/s13*(one+z1**2) &
         -s123/s12*(one+two*z1) -s123/s13*(one+two*z1) + overomep

    Pgqqb_ab_spav(1) = s23*(one/s12+one/s13) -s123/s12 - s123/s13 + overomep

    Pgqqb_ab_spav(2) = overomep
    
    return

  end function Pgqqb_ab_spav

  function Pgqqb_nab_spav(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pgqqb_nab_spav(0:2)
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3
    real(dp) :: t231,overomep

    t231 = tijk(s23,s12,s13,z2,z3,z1)

    !-- piece proportional to 1/(1-ep)
    overomep = -s123**2/s13/s23*z3**2*(one-z3-two*z1*z2)/z1/(one-z1) &
         -s123/s13*(one-z2)**2*z2/z1/(one-z1) &
         +s123/two/s23*(z1*(z3-z2)**2-two*z2*z3*(one+z1))/z1/(one-z1) &
         +s123**2/two/s12/s13*(z1+two*z2*z3)

    Pgqqb_nab_spav(0) = -t231**2/four/s23**2 + s123**2/two/s13/s23*z3*( &
         ((one-z1)**3-z1**3)/z1/(one-z1) ) &
         +s123/two/s13*(one-z2)*(one+one/z1/(one-z1)) &
         +s123/two/s23*(one+z1**3)/z1/(one-z1) &
          -one/four - s123**2/two/s12/s13*(one+z1**2) + overomep

    Pgqqb_nab_spav(1) = half + overomep

    Pgqqb_nab_spav(2) = overomep

    return

  end function Pgqqb_nab_spav

  function Pgqqb_nab_spav_symm(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pgqqb_nab_spav_symm(0:2)
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3

    Pgqqb_nab_spav_symm = Pgqqb_nab_spav(s12,s13,s23,s123,z1,z2,z3) + &
         Pgqqb_nab_spav(s13,s12,s23,s123,z1,z3,z2)

    return

  end function Pgqqb_nab_spav_symm

  function Pggg_spav(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pggg_spav(0:2)
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3
    real(dp) :: t123

    t123 = tijk(s12,s13,s23,z1,z2,z3)

    Pggg_spav(0) = t123**2/four/s12**2 + three/four + &
         s123/s12 * (four*(z1*z2-one)/(one-z3) + (z1*z2-two)/z3 + three/two + five/two*z3 &
         +(one-z3*(one-z3))**2/z3/z1/(one-z1) ) &
         +s123**2/s12/s13 * ( z1*z2*(one-z2)*(one-two*z3)/z3/(one-z3) + z2*z3 - two &
         + z1*(one+two*z1)/two + (one+two*z1*(one+z1))/two/(one-z2)/(one-z3) + &
         (one-two*z1*(one-z1))/two/z2/z3 )

    Pggg_spav(1) = -t123**2/four/s12**2 - three/four

    Pggg_spav(2) = zero

    Pggg_spav(:) = Pggg_spav(:) * Ca**2

    return

  end function Pggg_spav

  !---

  function tijk(sij,sik,sjk,zi,zj,zk)
    real(dp) :: tijk
    real(dp), intent(in) :: sij,sik,sjk,zi,zj,zk

    tijk = two * (zi*sjk-zj*sik)/(zi+zj) + (zi-zj)/(zi+zj)*sij

    return

  end function tijk
   
 end module mod_tripleclimits_q_qg_gstar



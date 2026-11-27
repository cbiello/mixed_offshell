module mod_tripleclimits_q_qqb_qstar
  use mod_types; use mod_consts_dp; use common_def
  use mod_kinematics_tc_ord
  use mod_auxfunctions

  implicit none



  private
  
  public :: PTC_A_q_qqb_qst, PTC_B_q_qqb_qst, PTC_C_q_qqb_qst, PTC_D_q_qqb_qst
  public :: tripcoll_subtr_q_qqb_qstar


contains


  function tripcoll_subtr_q_qqb_qstar(ndim,yRnd,ncomp,ff,userdata,weight,iternumber)
    implicit none
   !--
    integer, parameter :: ep = 1
    integer, parameter :: myfl = 3 !1=singlet diff, 2=singlet id, 3=nonsinglet diff, 4=nonsinglet id
    integer, parameter :: mycontr = 3 !-- 1 -> delta, 2 -> plus, 3 -> regular
    real(dp), parameter :: myz = 0.123_dp
    real(dp), parameter :: zmin = zero
    !--
    real(dp), parameter :: buff = 1E-10_dp
    real(dp), parameter :: onet = one - two * buff
    real(dp), parameter :: bufftc = 1E-10_dp
    !--
    real(dp15) :: yRnd(30), ff(1), weight
    integer :: ndim, ncomp, userdata, iternumber,tripcoll_subtr_q_qqb_qstar
    real(dp)              :: x1,x2,x3,xx2,xx(1:3),x3a,xx3
    real(dp)              :: res_tc(1:4,1:3,0:2)
    real(dp)              :: seca(1:4,1:3,0:2),secb(1:4,1:3,0:2),secc(1:4,1:3,0:2),secd(1:4,1:3,0:2),wghtac,wghtbd,z
    real(dp)              :: seca0(1:4,1:3,0:2),secb0(1:4,1:3,0:2),secc0(1:4,1:3,0:2),secd0(1:4,1:3,0:2)
    logical, save :: first=.true.


    tripcoll_subtr_q_qqb_qstar = zero
    z = myz
    ff(1) = zero
    res_tc = zero

    seca = zero; seca0 = zero
    secb = zero; secb0 = zero
    secc = zero; secc0 = zero
    secd = zero; secd0 = zero


    x1  = buff+onet*real(yRnd(1),dp)  ! r
    x2  = buff+onet*real(yRnd(2),dp)  ! x4
    x3a = buff+onet*real(yRnd(3),dp)            ! lambda
    xx3 = buff+onet*real(yRnd(4),dp)            ! x3
! remap lambda
    x3 = sin( half * pi * x3a )**2

    if (x1*x2*(one-z) .lt. bufftc) return


    xx = (/x1,x2,x3/)
    seca0 = PTC_A_q_qqb_qst(z,zmin,xx,zero)
    secc0 = PTC_C_q_qqb_qst(z,zmin,xx,zero)

    if (withmax) then
       seca = PTC_A_q_qqb_qst(z,zmin,xx,x3)
       secc = PTC_C_q_qqb_qst(z,zmin,xx,x3)
    endif

    wghtac = -one/x1/x2

    xx2 = x2**2
    xx = (/x1,xx2,x3/)
    secb0 = PTC_B_q_qqb_qst(z,zmin,xx,zero)
    secd0 = PTC_D_q_qqb_qst(z,zmin,xx,zero)
    if (withmax) then
       secb = PTC_B_q_qqb_qst(z,zmin,xx,x3)
       secd = PTC_D_q_qqb_qst(z,zmin,xx,x3)
    endif

    ! debug
!    seca=zero
!    secc=zero
!    secb=zero
!    secd=zero
    
    wghtbd = -one/x1/x2*two
    
    if (xx2 .lt. bufftc) then
       res_tc = (seca0+secc0)*wghtac
       if (withmax) then 
          res_tc(:,:,1) = res_tc(:,:,1) - two*( &
               seca(:,:,0)-seca0(:,:,0) + &
               secc(:,:,0)-secc0(:,:,0))/xx3*wghtac
       endif
    else
       res_tc = (seca0+secc0)*wghtac + (secb0 + secd0)*wghtbd
       if (withmax) then
          res_tc(:,:,1) = res_tc(:,:,1) - two*( &
               seca(:,:,0)-seca0(:,:,0) + &
               secc(:,:,0)-secc0(:,:,0))/xx3*wghtac
          
          res_tc(:,:,1) = res_tc(:,:,1) - two*( &
               secb(:,:,0)-secb0(:,:,0) + &
               secd(:,:,0)-secd0(:,:,0))/xx3*wghtbd
       endif
    endif

    ! remove color factor Cf(Cf-2Ca) for the interference case (myfl=4); otherwisre remove color factor Cf * Tr
    if (myfl .eq. 4 .or. myfl .eq. 2) then
       res_tc = res_tc/Cf/(Ca-two*Cf)
    else
       res_tc = res_tc/Cf/tr
    endif


! from remapping
    
    ff(1) = res_tc(myfl,mycontr,ep) 
    
    if (ff(1) .ne. ff(1) .or. weight .ne. weight) then
       print *, "xx",xx,x3a
       print *, "wghtbd",wghtbd
       print *, "sec a",seca(myfl,mycontr,ep)
       print *, "sec b",secb(myfl,mycontr,ep)
       print *, "sec c",secc(myfl,mycontr,ep)
       print *, "sec d",secd(myfl,mycontr,ep)
       print *, "ff",ff(1)
       stop
    endif




  end function tripcoll_subtr_q_qqb_qstar


  function P_qbqq(s12,s13,s23,z1,z2)
    ! indices: {ep^0,ep^1,ep^2}
    real(dp) :: s12,s13,s23,z1,z2,P_qbqq(0:2)
    real(dp) :: t123, s123, z3

    P_qbqq = zero 

    z3 = one - z1 - z2
    s123 = s12 + s13 + s23
    t123 = (two * (z1*s23-z2*s13)+(z1-z2)*s12)/(z1+z2)

    P_qbqq(0) = half * cf * tr * s123/s12*(-t123**2/s12/s123  &
         + (four*z3+(z1-z2)**2)/(z1+z2)  &
         + z1 + z2 - s12/s123  )

    P_qbqq(1) = half * Cf * tr * s123/s12*(-two*(z1+z2-s12/s123) ) 

    return

  end function P_qbqq

  function P_qbqq_id_bit_unsymm(s12,s13,s23,z1,z2)   ! unsymmetrized relative to Catani
    ! indices: {ep^0,ep^1,ep^2}
    real(dp) :: s12,s13,s23,z1,z2,P_qbqq_id_bit_unsymm(0:2) 
    real(dp) :: z3, s123 

    z3 = one - z1 - z2
    s123 = s12 + s13 + s23

    P_qbqq_id_bit_unsymm(0)  = Cf* (Cf - half*Ca) *(  two*s23/s12            & 
         + s123/s12*( ( one + z1**2)/(one-z2) - two*z2/(one-z3)  ) & 
         - s123**2/s12/s13*z1*half*(one+z1**2)/(one-z2)/(one-z3) )

    P_qbqq_id_bit_unsymm(1)  = Cf* (Cf - half*Ca) *(  -two*s23/s12         & 
         - one  -s123/s12*( (one-z3)**2/(one-z2) + one+z1-two*z2/(one-z3) ) & 
         +s123**2/s12/s13*z1*half*( one+ two*(one-z2)/(one-z3) ) )

    P_qbqq_id_bit_unsymm(2) = Cf*(Cf-half*Ca) *(   one  -s123/s12*(one-z3)  +s123**2/s12/s13*z1*half )

    return

  end function P_qbqq_id_bit_unsymm


  function P_qq_r_stro(z1,z2,kt1_kt2_sq) 
    ! indices: {ep^0,ep^1,ep^2}
    real(dp) :: z1,z2,kt1_kt2_sq,P_qq_r_stro(0:2) 

    P_qq_r_stro = zero 

    P_qq_r_stro(0) = Cf*Tr*(  ( one + (one-z2)**2)/z2 - two*z2*z1*(one-z1) & 
         - 8.0_dp*z1*(one-z1)*(one - z2)/z2*(kt1_kt2_sq) )

    P_qq_r_stro(1) = Cf*Tr*(-z2)
    return 
  end function P_qq_r_stro


  function P_qbqq_ds(s12,s13,s23,z1,z2)
    ! indices: {ep^0,ep^1,ep^2}
    real(dp) :: s12,s13,s23,z1,z2,P_qbqq_ds(0:2)
    real(dp) :: t123, s123, z3

    P_qbqq_ds = zero 

    z3 = one 
    s123 = s13 + s23
    t123 = two * (z1*s23-z2*s13)/(z1+z2)

    P_qbqq_ds(0) = half * cf * tr * s123/s12*(-t123**2/s12/s123 + four/(z1+z2) )

    return

  end function P_qbqq_ds


  function P_qbqq_ds_c54(s12,s13,s23,z1,z2,lambda)
    ! indices: {ep^0,ep^1,ep^2}
    real(dp) :: s12,s13,s23,z1,z2,P_qbqq_ds_c54(0:2)
    real(dp) :: t123, s123, z3, aux, kt1_kt2_sq,lambda

    P_qbqq_ds_c54 = zero 

    kt1_kt2_sq = lambda

    s123 = s13 + s23

    P_qbqq_ds_c54(0) = -cf* tr*s123/s12*( two - 8.0_dp*z1*z2*kt1_kt2_sq )


    return

  end function P_qbqq_ds_c54




  !------------------------------------------------------------
  !   I am considering the following process: 
  !
  !    quark(p1) -> quark(virtual) + quark(p4) + a/quark(p5)   
  !    
  !--------------------------------------------------------------


  !    triple collinear splitting function  Eq.(57) Catani

  function TC(kinlist) 
    ! indices: {ep^0,ep^1,ep^2}
    real(dp) :: TC(4,0:2) 
    real(dp), intent(in) :: kinlist(6)
    real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123,  z1, z2, z3, zbase 

    s14 = kinlist(1) 
    s15 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 
    zbase = one - z4 - z5


    TC = zero 

    !-- singlet_diff

    TC(1,0:2) = (P_qbqq(-s14,-s15,s45,zbase,z4) &
         +  P_qbqq(-s15,-s14,s45,zbase,z5))!*vn_tagc

    !-- singlet_id

    TC(2,0:2) = (P_qbqq_id_bit_unsymm(-s14,-s15,s45,zbase,z4) &
         + P_qbqq_id_bit_unsymm(-s15,-s14,s45,zbase,z5))!*vn_tagce

    !-- non_singlet_diff

    ! add the swap 4<-->5 to include qbar q final states
    TC(3,0:2) = (P_qbqq(s45,-s15,-s14,z5,z4)  + P_qbqq(s45,-s14,-s15,z4,z5))!*vn_taga

    !-- non_singlet_id

    TC(4,0:2) = (P_qbqq_id_bit_unsymm(-s14,s45,-s15,z4,zbase) & 
         + P_qbqq_id_bit_unsymm(s45,-s14,-s15,z4,z5) &
         + P_qbqq_id_bit_unsymm(-s15,s45,-s14,z5,zbase) & 
         + P_qbqq_id_bit_unsymm(s45,-s15,-s14,z5,z4))!*vn_tagac

    ! -- final normalization for everyting 

    TC = four/(s14+s15-s45)**2*TC


  end function TC



  function TC_C41(lambda,kinlist) 
    real(dp) :: TC_C41(4,0:2),lambda 
    real(dp), intent(in) ::  kinlist(6)
    real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, w, Pqq, Pqq1,  Pgg, Pgg1,  Pggperp, z
    real(dp) :: Pggperp1
    real(dp) :: kperpprod, kt1_kt2_sq,ztilde5,ztilde4


    s14 = kinlist(1) 
    s15 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5)
    z1 = kinlist(6)
    ztilde4 = z4/(z1+z4) !=E4/(-E1+E4), i.e. what was called z4 in mod_triplclimits_2q
    ztilde5 = z1+z4 !=(-E1+E4)/(-E1+E4+E5), i.e. what was called z5 in mod_triplclimits_2q
    
    !    kt1_kt2 = kinlist(6)
    kt1_kt2_sq = (one-two*lambda)**2

    TC_C41 = zero 

    !-- singlet_diff

    TC_C41(1,0:2) = -four/s14/(s45-s15)*P_qq_r_stro(ztilde4,ztilde5,kt1_kt2_sq) !* vn_tagc


  end function TC_C41



    function TC_C51(lambda,kinlist) 
    real(dp) :: TC_C51(4,0:2),lambda 
    real(dp), intent(in) ::  kinlist(6)
    real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, w, Pqq, Pqq1,  Pgg, Pgg1,  Pggperp, z
    real(dp) :: Pggperp1
    real(dp) :: kperpprod, kt1_kt2_sq,ztilde5,ztilde4


    s14 = kinlist(1) 
    s15 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5)
    z1 = kinlist(6)
    ztilde4 = z5/(z1+z5) !=E4/(-E1+E4), i.e. what was called z4 in mod_triplclimits_2q
    ztilde5 = z1+z5 !=(-E1+E4)/(-E1+E4+E5), i.e. what was called z5 in mod_triplclimits_2q
    
    !    kt1_kt2 = kinlist(6)
    kt1_kt2_sq = (one-two*lambda)**2

    TC_C51 = zero 

    !-- singlet_diff

    TC_C51(1,0:2) = -four/s15/(s45-s14)*P_qq_r_stro(ztilde4,ztilde5,kt1_kt2_sq) !* vn_tagc


  end function TC_C51



  function TC_C54(lambda,kinlist) 
    real(dp) :: TC_C54(4,0:2) 
    real(dp), intent(in) :: kinlist(6),lambda
    real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, w, Pqq, Pqq1,  Pgg, Pgg1,  Pggperp, z,ztilde5,ztilde4
    real(dp) :: kt1_kt2_sq


! note the convention: kinlist=(s14,s15,s45,z4,z5,z1)
    s14 = kinlist(1) 
    s15 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 

    z45 = z4+z5
    ztilde5 = z5/z45      !=E5/(E4+E5), i.e. what was called z5 in the same fn in mod_triplclimits_2q
    ztilde4 = z4/z45      !=E4/(E4+E5)
    kt1_kt2_sq = lambda


    TC_C54 = zero 

    !-- non_singlet_diff

    ! add the swap 4<-->5 to include qbar q final states
    TC_C54(3,0:2) = -four/s45/(s14+s15)*&
         (P_qq_r_stro(ztilde5,z45,kt1_kt2_sq) +P_qq_r_stro(ztilde4,z45,kt1_kt2_sq)) !*vn_taga

    return 

  end function TC_C54


  function TC_DS(kinlist) 
    real(dp) :: TC_DS(4,0:2) 
    real(dp), intent(in) :: kinlist(6)
    real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, w, Pqq, Pqq1,  Pgg, Pgg1,  Pggperp, z

    s14 = kinlist(1) 
    s15 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 

    TC_DS = zero 

    ! add the swap 4<-->5 to include qbar q final states
    TC_DS(3,0:2) = (P_qbqq_ds(s45,-s15,-s14,z5,z4)+P_qbqq_ds(s45,-s14,-s15,z4,z5))!*vn_taga

    !-- overal normalization 

    TC_DS = four/(s14+s15)**2*TC_DS

    return 
  end function TC_DS


  function TC_C54_DS(lambda,kinlist) 
    real(dp) :: TC_C54_DS(4,0:2) 
    real(dp), intent(in) :: kinlist(6),lambda
    real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3
    real(dp) :: kt1_kt2

    s14 = kinlist(1) 
    s15 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 
!    t123 = kinlist(6)
    t123 = lambda

    TC_C54_DS = zero 

    !-- non_singlet_diff

    TC_C54_DS(3,0:2) = (P_qbqq_ds_c54(s45,-s15,-s14,z5,z4,lambda)+P_qbqq_ds_c54(s45,-s14,-s15,z4,z5,lambda))!*vn_taga

    ! -- overal normalization 

    TC_C54_DS = four/(s14+s15)**2*TC_C54_DS

    return 

  end function TC_C54_DS

  ! these single soft functions are all 0

  function TC_DS_SO(kinlist)
    real(dp) :: TC_DS_SO(4,0:2) 
    real(dp), intent(in) :: kinlist(6)

    TC_DS_SO  = zero

  end function TC_DS_SO

  function TC_DS_SO_C54(kinlist)
    real(dp) :: TC_DS_SO_C54(4,0:2) 
    real(dp), intent(in) :: kinlist(6)

    TC_DS_SO_C54  = zero

  end function TC_DS_SO_C54


    function TC_DS_SO_C41(kinlist)
      real(dp) :: TC_DS_SO_C41(4,0:2) 
      real(dp), intent(in) :: kinlist(6)
      
      TC_DS_SO_C41  = zero
      
    end function TC_DS_SO_C41

    function TC_S5_C41(kinlist)   !????
      real(dp) :: TC_S5_C41(4,0:2) 
      real(dp), intent(in) :: kinlist(6)
      
      TC_S5_C41  = zero
      
    end function TC_S5_C41


    function TC_S5_C54(kinlist)   
      real(dp) :: TC_S5_C54(4,0:2) 
      real(dp), intent(in) :: kinlist(6)
      
      TC_S5_C54  = zero
      
    end function TC_S5_C54

    
    function TC_S5(kinlist)   
      real(dp) :: TC_S5(4,0:2) 
      real(dp), intent(in) :: kinlist(6)
      
      TC_S5  = zero
      
    end function TC_S5
  


    function TC_DS_C41(kinlist)   !????
      real(dp) :: TC_DS_C41(4,0:2) 
      real(dp), intent(in) :: kinlist(6)
      
      TC_DS_C41  = zero
      
    end function TC_DS_C41
  

  
! --- The triple collinear Altarelli-Parisi function; sector ``A'' -----------------------


   function PTC_A_q_qqb_qst(z,zmin,yr,x3)
   real(dp) :: PTC_A(1:4,1:3,0:2),PTC_A_q_qqb_qst(1:4,1:3,0:2)
   real(dp), intent(in) :: yr(3), z,zmin,x3
   real(dp) ::  rho51, full(6), fullc51(6),soft5(6), soft5c5(6),dsoft(6), dsord(6), dsoftz(6), dsoftzc5(6)
   real(dp) :: dsoft5(6), dsordc5(6), dsoftc5(6)
   real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4 , PHSP1m1, energy, onemz
   real(dp) :: faux1(1:4,0:2), faux2(1:4,0:2), LogR, LogF, LogF0, faux3(1:4,0:2), faux4(1:4,0:2), faclambda 


   r =  yr(1) 
   x4 = yr(2) 
   lambda = yr(3)
   PTC_A = zero 

   call kinematics_tcl_A(z,(/r,x4,lambda,x3/),rho51,full,fullc51,soft5,soft5c5,dsoft,dsord,dsoftc5,dsordc5,dsoftz,dsoftzc5)

   PHSP1 = (one - half*r) 
   PHSP2 =  (one-half*x4)/NF1(zero,0.5_dp*x4,lambda)

   faclambda = lambda*(one-lambda) 

   LogR  = log( half*r*(one-half*r))
   LogF  = log( 128.0_dp * PHSP2**2 * faclambda )
   LogF0 = log( 128.0_dp * faclambda   )

!  the delta-function coefficient 


   
   faux1(:,0:2) = PHSP2*rho51*( r**2 * TC_DS(dsoftz) - TC_DS_SO(dsord) ) 
   faux2(:,0:2) = ( r**2 * TC_DS_C41(dsoftzc5) - TC_DS_SO_C41(dsordc5) )   
   faux3(:,0:2) = PHSP2*rho51*( r**2 * TC_DS(dsoftz) )                   
   faux4(:,0:2) = r**2 * TC_DS_C41(dsoftzc5)                              

   PTC_A(1:4,1,0:2) = -log(one-zmin)*(faux1 - faux2)  - log(one+r)* (faux3 - faux4)  

   PTC_A(1:4,1,1) = PTC_A(1:4,1,1) + two*log(one-zmin)**2*( faux1(1:4,0) - faux2(1:4,0) ) & 
                   -two*Log(r)*PTC_A(1:4,1,0) - two*log(one+r)**2 * (faux3(1:4,0) - faux4(1:4,0) ) & 
                   +log(one-zmin)*( LogF *faux1(1:4,0) - LogF0*faux2(1:4,0) ) &
                   +log(one+r)*( LogF*faux3(1:4,0) - LogF0*faux4(1:4,0) )          

!  the plus distribution coefficient

     faux1(:,0:2) = rho51*PHSP2*(  PHSP1**2*( half*r )**2 * TC_DS(dsoft) - TC_DS_SO(dsord) )
     faux2(:,0:2) = PHSP1**2* (half*r )**2 * TC_DS_C41(dsoftc5) - TC_DS_SO_C41(dsordc5)
     faux3(:,0:2) = r*PHSP1*(rho51*PHSP2*TC_DS_SO(dsord)  - TC_DS_SO_C41(dsordc5))   ! I multiply with r to make integration 

     PTC_A(1:4,2,0:2) =  faux1 - faux2

     PTC_A(1:4,2,1) = PTC_A(1:4,2,1) - two * LogR*PTC_A(1:4,2,0) - ( LogF*faux1(1:4,0) - LogF0*faux2(1:4,0) ) & 
                + pisq/three*faux3(1:4,0) 


!  the regular piece 

   onemz = one-z
   onemz4 = onemz**4

   faux1(:,0:2) =   rho51*PHSP2*(  onemz4*(  z*PHSP1**2*(half*r)**2 * TC(full)     -  TC_S5(soft5) )   &        !! 
                           - PHSP1**2*(half*r)**2 * TC_DS(dsoft) +  TC_DS_SO(dsord)  )
   faux2(:,0:2) =  (onemz4*(   z*PHSP1**2*(half*r)**2 * TC_C51(lambda,fullc51)  -  TC_S5_C41(soft5c5) )   &
                     - PHSP1**2*(half*r)**2 * TC_DS_C41(dsoftc5) +  TC_DS_SO_C41(dsordc5)  )
   faux3(:,0:2) =  r*PHSP1* ( rho51* PHSP2*( onemz4* TC_S5(soft5)  - TC_DS_SO(dsord) ) & 
        - ( onemz4*TC_S5_C41(soft5c5) - TC_DS_SO_C41(dsordc5) ) )
    
   PTC_A(1:4,3,0:2) = faux1 - faux2 

   PTC_A(1:4,3,1) = PTC_A(1:4,3,1)   -( two*LogR + four*log(onemz) ) *PTC_A(1:4,3,0) & 
                                     -( LogF*faux1(1:4,0) - LogF0*faux2(1:4,0) ) &
                                     + pisq/three*faux3(1:4,0) 

   PTC_A(1:4,3,0:2) = PTC_A(1:4,3,0:2)/onemz

   PTC_A(1:4,2:3,0:2) = one/PHSP1*PTC_A(1:4,2:3,0:2)

   PTC_A(1:4,1:3,1) = PTC_A(1:4,1:3,1) - PTC_A(1:4,1:3,0)*log(x4)

   PTC_A_q_qqb_qst = PTC_A
   
   
   end function 


! --- The triple collinear Altarelli-Parisi function; sector ``B '' -----------------------


   function PTC_B_q_qqb_qst(z,zmin,yr,x3)
   real(dp) :: PTC_B(1:4,1:3,0:2), PTC_B_q_qqb_qst(1:4,1:3,0:2)
   real(dp), intent(in) :: yr(3), z,zmin,x3
   real(dp) ::  rho45, full(6), fullc54(6),soft5(6), soft5c5(6),dsoft(6), dsord(6)
   real(dp) ::  dsoftc5(6), dsordc5(6), dsoftz(6), dsoftzc5(6)
   real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4, onemz 
   real(dp) :: faux1(1:4,0:2), faux2(1:4,0:2), LogR, LogF, LogF0, faux3(1:4,0:2), faux4(1:4,0:2),faclambda,tmp


   r =  yr(1) 
   x4 = yr(2) 
   lambda = yr(3) 

   PTC_B = zero 

   call kinematics_tcl_B(z,(/r,x4,lambda,x3/),rho45,full,fullc54,soft5,soft5c5,dsoft,dsord,dsoftc5,dsordc5,dsoftz,dsoftzc5)

   PHSP1 = (one - half*r) 
   PHSP2 = one   ! everything is absorbed into rho45 

   faclambda = lambda * (one - lambda )

   LogR =  log( half*r*(one-half*r))
   LogF =  log( 64.0_dp* (one-half*x4)/NF1(zero,one-0.5_dp*x4,lambda)**2*faclambda)
   LogF0 = log( 64.0_dp*one/NF1(zero,one,lambda)**2 * faclambda )

    faux1(:,0:2)= PHSP2*rho45*(r**2*TC_DS(dsoftz) - TC_DS_SO(dsord)  )
    faux2(:,0:2)= r**2*TC_C54_DS(lambda,dsoftzc5) - TC_DS_SO_C54(dsordc5)
    faux3(:,0:2)= PHSP2*rho45*( r**2 * TC_DS(dsoftz) )

! RR: this is a pretty glorious hack to get the right subtraction for the 5||4 double soft limit in TC_C54_DS
! Ideally we want to put this quantity in the kinlist for dsoftzc5 but this may break something somewhere else,
! so leave it as is for now...
    tmp=1.0_dp/8.0_dp*(two/r-2/r/(1+r)+8.0_dp*lambda/(one+r)**3)
    faux4(:,0:2)= r**2 * TC_C54_DS(tmp,dsoftzc5)

! delta coefficient

   PTC_B(1:4,1,0:2) = -log(one-zmin)*(faux1(1:4,0:2) - faux2(1:4,0:2))  - log(one+r)* (faux3(1:4,0:2) - faux4(1:4,0:2))  

   PTC_B(1:4,1,1) = PTC_B(1:4,1,1) + two*log(one-zmin)**2*( faux1(1:4,0) - faux2(1:4,0) ) & 
                   -two*Log(r)*PTC_B(1:4,1,0) - two*log(one+r)**2 * (faux3(1:4,0) - faux4(1:4,0) ) & 
                   +log(one-zmin)*( LogF *faux1(1:4,0) - LogF0*faux2(1:4,0) ) &
                   +log(one+r)*( LogF*faux3(1:4,0) - LogF0*faux4(1:4,0) )     


!  the plus distribution coefficient

   faux1(1:4,0:2) = rho45*PHSP2*( PHSP1**2*(half*r)**2*TC_DS(dsoft) - TC_DS_SO(dsord)  ) 
   faux2(1:4,0:2) = PHSP1**2*(half*r)**2*TC_C54_DS(lambda,dsoftc5)  - TC_DS_SO_C54(dsordc5) 
   faux3(1:4,0:2) = r*PHSP1*(rho45*PHSP2*TC_DS_SO(dsord) - TC_DS_SO_C54(dsordc5) )

    PTC_B(1:4,2,0:2) =  faux1(1:4,0:2) - faux2(1:4,0:2) 

    PTC_B(1:4,2,1) = PTC_B(1:4,2,1) - two * LogR*PTC_B(1:4,2,0) & 
                -( LogF*faux1(1:4,0) - LogF0*faux2(1:4,0) ) & 
                + pisq/three*faux3(1:4,0) 

!  the regular function 

   onemz = one -z 
   onemz4 = onemz**4

   faux1(1:4,0:2)=rho45*PHSP2*(   onemz4 * ( PHSP1**2*(half*r)**2 *z*TC(full)  - TC_S5(soft5)) &
                             - PHSP1**2*(half*r)**2*TC_DS(dsoft) + TC_DS_SO(dsord) )

   faux2(1:4,0:2) = onemz4 *( PHSP1**2*(half*r)**2*z*TC_C54(lambda,fullc54) - TC_S5_C54(soft5c5) )  &
                - PHSP1**2*(half*r)**2*TC_C54_DS(lambda,dsoftc5)  + TC_DS_SO_C54(dsordc5) 
   faux3(1:4,0:2)= r *PHSP1*( rho45*PHSP2* ( onemz4 * TC_S5(soft5) -  TC_DS_SO(dsord) ) & 
                         -  ( onemz4 * TC_S5_C54(soft5c5)  -  TC_DS_SO_C54(dsordc5)  ) )

   PTC_B(1:4,3,0:2) = faux1(1:4,0:2) - faux2(1:4,0:2)

   PTC_B(1:4,3,1) = PTC_B(1:4,3,1) - ( two*LogR + four*log(onemz) ) *PTC_B(1:4,3,0) &
                                   - ( LogF*faux1(1:4,0) - LogF0*faux2(1:4,0) ) &
                                   + pisq/three*faux3(1:4,0) 
     
   PTC_B(1:4,3,0:2) = one/onemz*PTC_B(1:4,3,0:2)

   PTC_B(1:4,2:3,0:2) = one/PHSP1*PTC_B(1:4,2:3,0:2)

   PTC_B(1:4,1:3,1) = PTC_B(1:4,1:3,1) - two*PTC_B(1:4,1:3,0)*log(x4)

   PTC_B_q_qqb_qst = PTC_B

 end function PTC_B_q_qqb_qst



! --- The triple collinear Altarelli-Parisi function; sector ``C'' -----------------------


   function PTC_C_q_qqb_qst(z,zmin,yr,x3)
   real(dp) :: PTC_C(1:4,1:3,0:2),PTC_C_q_qqb_qst(1:4,1:3,0:2)
   real(dp), intent(in) :: yr(3), z,zmin,x3
   real(dp) ::  rho41, full(6), fullc41(6),soft5(6), soft5c41(6),dsoft(6), dsord(6)
   real(dp)  :: dsoftc41(6), dsordc41(6), dsoftz(6), dsoftzc41(6)
   real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4, onemz, LogR, LogF, LogF0
   real(dp) ::  faux1(1:4,0:2), faux2(1:4,0:2), faux3(1:4,0:2), faux4(1:4,0:2), faclambda  


   r =  yr(1) 
   x4 = yr(2) 
   lambda = yr(3)

   call kinematics_tcl_C(z,(/r,x4,lambda,x3/),rho41,full,fullc41,soft5,soft5c41,dsoft,dsord,dsoftc41,dsordc41,dsoftz,dsoftzc41)

   PTC_C = zero 

   PHSP1 = (one - half*r) 
   PHSP2 =  (one-half*x4)/NF1(zero,0.5_dp*x4,lambda)

   faclambda = lambda * (one - lambda )

   LogR  = log( half*r*(one-half*r) )
   LogF  = log( 128.0_dp * PHSP2**2 * faclambda )
   LogF0 = log( 128.0_dp * faclambda )


!  the delta-function coefficient 

   faux1(1:4,0:2) =  PHSP2*rho41*( r**2 * TC_DS(dsoftz) - TC_DS_SO(dsord)   )
   faux2(1:4,0:2) = r**2 * TC_DS_C41(dsoftzc41) - TC_DS_SO_C41(dsordc41) 
   faux3(1:4,0:2) = PHSP2*rho41* r**2 * TC_DS(dsoftz)
   faux4(1:4,0:2) = r**2 * TC_DS_C41(dsoftzc41) 

   PTC_C(1:4,1,0:2) = -log(one-zmin)*(faux1(1:4,0:2) - faux2(1:4,0:2))  - log(one+r)* (faux3(1:4,0:2) - faux4(1:4,0:2))  

   PTC_C(1:4,1,1) = PTC_C(1:4,1,1) + two*log(one-zmin)**2*( faux1(1:4,0) - faux2(1:4,0) ) & 
                   -two*Log(r)*PTC_C(1:4,1,0) - two*log(one+r)**2 * (faux3(1:4,0) - faux4(1:4,0) ) & 
                   +log(one-zmin)*( LogF *faux1(1:4,0) - LogF0*faux2(1:4,0) ) &
                   +log(one+r)*( LogF*faux3(1:4,0) - LogF0*faux4(1:4,0) )          


!  the plus distribution coefficient

   faux1(1:4,0:2) = rho41*PHSP2*( PHSP1**2*( half*r)**2 * TC_DS(dsoft) - TC_DS_SO(dsord) ) 
   faux2(1:4,0:2) = PHSP1**2*(half*r)**2 * TC_DS_C41(dsoftc41) - TC_DS_SO_C41(dsordc41)
   faux3(1:4,0:2) = r*PHSP1*(rho41*PHSP2*TC_DS_SO(dsord) - TC_DS_SO_C41(dsordc41) )


   PTC_C(1:4,2,0:2) = faux1(1:4,0:2) - faux2(1:4,0:2)

   PTC_C(1:4,2,1) = PTC_C(1:4,2,1) - two * LogR*PTC_C(1:4,2,0) & 
                                  -( LogF*faux1(1:4,0) - LogF0*faux2(1:4,0) ) & 
                                  + pisq/three*faux3(1:4,0) 


   onemz = one-z
   onemz4 = onemz**4

   faux1(1:4,0:2) = rho41*PHSP2*( onemz4*( z*PHSP1**2*(half*r)**2*TC(full) -  TC_S5(soft5)   )   &
               - PHSP1**2*(half*r)**2*TC_DS(dsoft)  + TC_DS_SO(dsord) )

   faux2(1:4,0:2) = onemz4*( z*PHSP1**2*(half * r)**2*TC_C41(lambda,fullc41)  - TC_S5_C41(soft5c41) )  &
                       -PHSP1**2*(half*r)**2*TC_DS_C41(dsoftc41) +  TC_DS_SO_C41(dsordc41) 
   faux3(1:4,0:2) = r*PHSP1*(  rho41*PHSP2*( onemz4*TC_S5(soft5) - TC_DS_SO(dsord)  ) & 
                    - ( onemz4*TC_S5_C41(soft5c41) - TC_DS_SO_C41(dsordc41) ) )

   PTC_C(1:4,3,0:2) = faux1(1:4,0:2) - faux2(1:4,0:2)


   PTC_C(1:4,3,1) =      PTC_C(1:4,3,1) - ( two*LogR + four*log(onemz) ) *PTC_C(1:4,3,0) &
                       - ( LogF*faux1(1:4,0) - LogF0*faux2(1:4,0) ) &
                       + pisq/three*faux3(1:4,0) 

   PTC_C(1:4,3,0:1) = PTC_C(1:4,3,0:1)/onemz

   PTC_C(1:4,2:3,0:1) = one/PHSP1*PTC_C(1:4,2:3,0:1)

   PTC_C(1:4,1:3,1) = PTC_C(1:4,1:3,1) - PTC_C(1:4,1:3,0)*log(x4)

   PTC_C_q_qqb_qst = PTC_C

 end function PTC_C_q_qqb_qst



! --- The triple collinear Altarelli-Parisi function; sector ``D '' -----------------------


   function PTC_D_q_qqb_qst(z,zmin,yr,x3)
   real(dp) :: PTC_D(1:4,3,0:2), PTC_D_q_qqb_qst(1:4,3,0:2)
   real(dp), intent(in) :: yr(3), z,zmin,x3
   real(dp) :: rho45, full(6), fullc54(6),soft5(6), soft5c5(6),dsoft(6), dsord(6)
   real(dp)  :: dsoftc5(6), dsordc5(6),dsoftz(6), dsoftzc5(6)
   real(dp) :: x4, r, lambda, PHSP1, PHSP2, onemz4, onemz 
   real(dp) :: faux1(1:4,0:2), faux2(1:4,0:2), LogR, LogF, LogF0, faux3(1:4,0:2), faux4(1:4,0:2), faclambda,tmp


   r =  yr(1) 
   x4 = yr(2) 
   lambda = yr(3) 

   PTC_D = zero 

   call kinematics_tcl_D(z,(/r,x4,lambda,x3/),rho45,full,fullc54,soft5,soft5c5,dsoft,dsord,dsoftc5,dsordc5,dsoftz,dsoftzc5)

   PHSP1 = (one - half*r) 
   PHSP2 = one   ! everything is absorbed into rho45 

   faclambda = lambda *( one -  lambda ) 

   LogR = log( half*r*(one-half*r))
   LogF = Log( 64.0_dp* (one-half*x4)/NF1(zero,one-half*x4,lambda)**2 * faclambda )
   LogF0 = Log( 64.0_dp*one/NF1(zero,one,lambda)**2 * faclambda )
   

!  the delta-function  

    faux1(1:4,0:2)= PHSP2*rho45*(r**2*TC_DS(dsoftz) - TC_DS_SO(dsord) )
    faux2(1:4,0:2)= r**2*TC_C54_DS(lambda,dsoftzc5) - TC_DS_SO_C54(dsordc5) 
    faux3(1:4,0:2)= PHSP2*rho45*( r**2 * TC_DS(dsoftz) )
! see comment in PTC_B about this quantity tmp
    tmp=1.0_dp/8.0_dp*(two/r-2/r/(1+r)+8.0_dp*lambda/(one+r)**3)
    faux4(1:4,0:2)= r**2 * TC_C54_DS(tmp,dsoftzc5)   


   PTC_D(1:4,1,0:2) = -log(one-zmin)*(faux1(1:4,0:2) - faux2(1:4,0:2))  - log(one+r)* (faux3(1:4,0:2) - faux4(1:4,0:2))  

   PTC_D(1:4,1,1) = PTC_D(1:4,1,1) + two*log(one-zmin)**2*( faux1(1:4,0) - faux2(1:4,0) ) & 
                   -two*Log(r)*PTC_D(1:4,1,0) - two*log(one+r)**2 * (faux3(1:4,0) - faux4(1:4,0) ) & 
                   +log(one-zmin)*( LogF *faux1(1:4,0) - LogF0*faux2(1:4,0) ) &
                   +log(one+r)*( LogF*faux3(1:4,0) - LogF0*faux4(1:4,0) )     



!  the plus distribution coefficient

   faux1(1:4,0:2) = rho45*PHSP2*( PHSP1**2*(half*r)**2*TC_DS(dsoft)  - TC_DS_SO(dsord) ) 
   faux2(1:4,0:2) = PHSP1**2*(half*r)**2*TC_C54_DS(lambda,dsoftc5)  - TC_DS_SO_C54(dsordc5)
   faux3(1:4,0:2) = r *PHSP1* (rho45*PHSP2*TC_DS_SO(dsord) - TC_DS_SO_C54(dsordc5) )


   PTC_D(1:4,2,0:2) =  faux1(1:4,0:2) - faux2(1:4,0:2)

    PTC_D(1:4,2,1) = PTC_D(1:4,2,1) - two * LogR*PTC_D(1:4,2,0) &
                 -( LogF*faux1(1:4,0) - LogF0*faux2(1:4,0) ) & 
                 + pisq/three*faux3(1:4,0) 

!  the regular function 

   onemz = one -z 
   onemz4 = onemz**4

   faux1(1:4,0:2)=rho45*PHSP2*(   onemz4 * ( PHSP1**2*(half*r)**2 *z*TC(full)  - TC_S5(soft5)  ) &
                             - PHSP1**2*(half*r)**2*TC_DS(dsoft) +TC_DS_SO(dsord))
 
    faux2(1:4,0:2) = onemz4 *( PHSP1**2*(half*r)**2*z*TC_C54(lambda,fullc54) - TC_S5_C54(soft5c5))  &
         - PHSP1**2*(half*r)**2*TC_C54_DS(lambda,dsoftc5) + TC_DS_SO_C54(dsordc5) 

   faux3(1:4,0:2) = r *PHSP1*( rho45*PHSP2* ( onemz4 * TC_S5(soft5) -  TC_DS_SO(dsord) ) & 
                         -  ( onemz4 * TC_S5_C54(soft5c5)  -  TC_DS_SO_C54(dsordc5)  ) )
    

   PTC_D(1:4,3,0:2) = faux1(1:4,0:2) - faux2(1:4,0:2) 

   PTC_D(1:4,3,1) = PTC_D(1:4,3,1) - ( two*LogR + four*log(onemz) ) *PTC_D(1:4,3,0) &
                        -( LogF*faux1(1:4,0) - LogF0*faux2(1:4,0) ) &
                        + pisq/three*faux3(1:4,0) 
     
   PTC_D(1:4,3,0:2) = one/onemz*PTC_D(1:4,3,0:2)

   PTC_D(1:4,2:3,0:2) = one/PHSP1*PTC_D(1:4,2:3,0:2)

   PTC_D(1:4,1:3,1) = PTC_D(1:4,1:3,1) - two*PTC_D(1:4,1:3,0)*log(x4)

   PTC_D_q_qqb_qst = PTC_D

 end function PTC_D_q_qqb_qst







end module mod_tripleclimits_q_qqb_qstar

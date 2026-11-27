module mod_kinematicsNNLO_EunordAord
  use mod_types; use mod_consts_dp; use common_def
  use mod_process
  use mod_auxfunctions; use mod_procfunctions
  use mod_kinematics_gen; use mod_kinematics_lept
  implicit none
  private

  !-- THESE ARE REDUNDANT, BUT FOR NOW WE DON'T CARE
  !-- LOOK AT ORIGINAL DY CODE FOR MINIMAL SET
  public ::  kinematics_nnlo_4252a, kinematics_nnlo_4252b, kinematics_nnlo_4252c, kinematics_nnlo_4252d 
  public ::  kinematics_nnlo_4251 

  public ::  kinematics_nnlo_4151a, kinematics_nnlo_4151b, kinematics_nnlo_4151c, kinematics_nnlo_4151d 
  public ::  kinematics_nnlo_4152

contains 

  !-- we factor out an extra 1/(8*pisq)**2 from the phase space
  !-- which multiplied by gs**4 in the ME reconstructs asontwopi**2
  subroutine kinematics_nnlo_4252a(MV2,Y,yr,&
       HardProc,S5Lim,TCLim,TCS5Lim,&
       C5Lim,C5S5Lim,TCC5Lim,TCC5S5Lim) 
    real(dp), intent(in) ::  MV2,Y, yr(9)
    type(KinConfig)      ::  HardProc,S5Lim,TCLim,TCS5Lim,&
         C5Lim,C5S5Lim,TCC5Lim,TCC5S5Lim
    !---
    real(dp) :: p1(4), p2(4), p4(4), p5(4)
    real(dp) :: MV, kallenF, sin5, cos5, pVV(4)
    real(dp) :: xx3,x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w1415, s24, s25, s45, z4, z5, ztot, weight4252, z45, ncoll(1:4), nperp(1:4)
    real(dp) :: spart,sqrts, yy, Deltas, acosphi45, asinphi45
    real(dp) :: s4_25,zz1,zz2,kt1kt2sq

    !    statement functions
    weight4252(eta41,eta42,eta51,eta52,eta54)  = eta41*eta51*( one + eta42/(eta54 + eta41 + eta52) & 
         + eta52/(eta54 + eta42 + eta51) )

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S5Lim)          !2
    call initialize_config(TCLim)          !3
    call initialize_config(TCS5Lim)        !4
    call initialize_config(C5Lim)          !5
    call initialize_config(C5S5Lim)        !6
    call initialize_config(TCC5Lim)        !7
    call initialize_config(TCC5S5Lim)      !8

    MV =  sqrt(MV2)

    x1 =   yr(1) 
    x2 =   yr(2) 
    x3 =   yr(3) 
    x4 =   yr(5) 

    if (withmax) then
       xx3 = x3
    else
       xx3 = zero
    endif

    lambda = sin( half * pi * yr(6) )**2



    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = -one + two*x3 
    sin4 = sqrt(one-cos4**2)

    cos5 = -one + x3*x4 
    sin5 = sqrt(one-cos5**2) 


    if (x4.eq.zero) then 
       if(one-two*lambda .ge. zero) then
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif
    elseif(2.0_dp*(one - half*x4)**2/NF1(x3,half*x4,lambda)-two*(one+half*x4*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif




    eta41 = one-x3
    eta51 = one-half*x3*x4
    eta42 = x3 
    eta52 = half*x3*x4
    eta54 = x3*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (HardProc%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) HardProc%Flag = .true.
    else
       HardProc%flag = .true.
    endif


    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif

    ! --------- SOFT 5 ---------

    eta42 = x3
    eta52 = half*x3*x4
    eta41 = (one-x3)
    eta51 = (one - half*x3*x4)
    eta54 = x3*(one-half*x4)**2/NF1(x3,half*x4,lambda)


    spart = MV2/(one - x1)
    if (spart.gt.sh) S5Lim%flag = .true.
    sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41)
    Deltad = one - x1*(one - eta42)

    if ( (S5Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) S5Lim%Flag = .true.
    else 
       S5Lim%flag = .true.
    endif

    if (S5Lim%flag .eqv. .false.) then 

       S5Lim%PartFrac = (/xi1,xi2/)


       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)

       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p4 

       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = p4
       S5Lim%LimMom(:,4) = p5 

       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2   


       ! kin inv and weight for C5 +  S5

       eta42 = x3
       eta52 = zero
       eta41 = (one-x3)
       eta51 = one
       eta54 = x3/NF1(x3,zero,lambda)

       s25 = spart*x2*(half*x3*x4)
       z5 = x2
       C5S5Lim%KinInv = (/zero,-s25,zero,zero,z5,zero/)

       C5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS 
    !-------------------------------------------------------------------------------------------------------------

    ! --------- TRIPLE COLLINEAR -------

    eta42 = zero
    eta52 = zero
    eta41 = one
    eta51 = one
    eta54 = zero 


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( TCLim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCLim%Flag = .true.
    else 
       TCLim%flag = .true.
    endif

    if (TCLim%flag .eqv. .false.) then 

       TCLim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, -one  /)
       p5 = half*sqrts*x2*(/one,zero,zero, -one /)

       s24 = spart *x1*x3 
       s25 = spart *x2*(half*x3*x4)
       s45 = spart *x1*x2*x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

       ztot = x1 + x2 - one 
       z4 = x1/ztot 
       z5 = x2/ztot

       TCLim%KinInv = (/-s24,-s25,s45,z4,z5,zero/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)

       TCLim%AmpMom(:,1) = p1
       TCLim%AmpMom(:,2) = p2 - p4 - p5

       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2


       ! --------- TRIPLE COLLINEAR + C5  -------

       s24 = spart *x1*x3 
       s25 = spart *x2*(half*x3*x4)
       s4_25 = s24 * (one - x2)

       zz1 = x2/(x2-one)
       zz2 = (one-x2)/(one-x1-x2)
       
       kt1kt2sq = (-one+two*lambda)**2 * (one-xx3)

       TCC5Lim%KinInv = (/-s4_25,-s25,zero,zz1,zz2,kt1kt2sq/)   ! the signs here are to compline with the original definition 
       ! q1 -> q g4 g5

       TCC5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    ! -------- TRIPLE COLLINEAR + SOFT 5


    spart = MV2/(one - x1) 
    if (spart.gt.sh) TCS5Lim%flag = .true. 

    sqrts = sqrt(spart)
    Deltan = one - x1*(one - eta41) 
    Deltad = one - x1*(one - eta42)

    if ( ( TCS5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCS5Lim%Flag = .true.
    else 
       TCS5Lim%flag = .true.
    endif

    if (TCS5Lim%flag .eqv. .false.) then 

       TCS5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, -one  /)
       p5 = half*sqrts*x2*(/one,zero,zero,-one /)

       z4 = x1/(x1-one)
       z5 = x2/(x1-one)

       s24 = spart *x1*x3 
       s25 = spart *x2*(half*x3*x4)
       s45 = spart *x1 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

       TCS5Lim%KinInv = (/-s24,-s25,s45,z4,z5,zero/)

       pVV = p1 + p2 - p4 

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1 
       TCS5Lim%AmpMom(:,2) = p2 - p4

       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            *spart/MV2


       ! -------   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT5 ---------
       z4 = x1
       z5 = x2

       s24 = spart *x1*x3 
       s25 = spart *x2*(half*x3*x4)
       s45 = spart *x1*x2*x3/NF1(xx3,zero,lambda)

       TCC5S5Lim%KinInv = (/-s24,-s25,s45,z4,z5,zero/)

       TCC5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C5 FLM 
    !-------------------------------------------------------------------------------------------------------------


    eta42 = x3
    eta52 = zero
    eta41 = (one-x3)
    eta51 = one
    eta54 = x3/NF1(x3,zero,lambda)

    phi4 = twopi*yr(4) 
    cos4 = -one + two*x3 
    sin4 = sqrt(one-cos4**2)


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C5Lim%flag = .true.
    if (C5Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C5Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) -  x2*(one - eta51)
    Deltad = one - x1*(one - eta42) -  x2*(one - eta52)  

    if ( ( C5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C5Lim%Flag = .true.
    else 
       C5Lim%flag = .true.
    endif

    if (C5Lim%flag .eqv. .false.) then 

       C5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, zero, zero, -one /)

       s25 = spart*x2*(half*x3*x4)
!       z5 = one/(one - x2)
       z5 = x2/(x2-one)

       C5Lim%KinInv = (/zero,-s25,zero,zero,z5,zero/)

       ncoll = (/ zero,  zero, zero , -one  /)

       if (lambda.gt.half) then 
          acosphi45 = -(one-two*lambda) 
          asinphi45 = two*sqrt(lambda*(one-lambda))
       else 
          acosphi45 = -(one-two*lambda) 
          asinphi45 = two*sqrt(lambda*(one-lambda))
       endif

       nperp  = (/ zero, cos(phi4)*acosphi45 - sin(phi4)*asinphi45 , & 
            sin(phi4)*acosphi45 + cos(phi4)*asinphi45, zero /)

       C5Lim%LimMom(:,1) = nperp
       C5Lim%LimMom(:,2) = ncoll

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C5Lim)

       C5Lim%AmpMom(:,1) = p1 
       C5Lim%AmpMom(:,2) = p2 - p5
       C5Lim%AmpMom(:,3) = p4 


       C5Lim%LimMom(:,1) = nperp
       C5Lim%LimMom(:,2) = ncoll


       C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif


  end subroutine kinematics_nnlo_4252a

  !-- we factor out an extra 1/(8*pisq)**2 from the phase space
  !-- which multiplied by gs**4 in the ME reconstructs asontwopi**2
  subroutine kinematics_nnlo_4252b(MV2,Y,yr,&
       HardProc,S5Lim,TCLim,TCS5Lim,&
       C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim) 
    real(dp), intent(in) ::  yr(9), MV2, Y
    type(KinConfig)      ::  HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C54Lim,C54DSLim,C54DSS5Lim,C54S5Lim,TCC54Lim,TCC54DSLim,TCC54DSS5Lim,TCC54S5Lim
    !---
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: xx3,x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, Deltas, w1415, s14, s15, s45, z4, z5, ztot, weight4252
    real(dp) :: spart,sqrts, x4p, nperp(4), s12, s24, z45, ncoll(1:4), yy, s25 
    real(dp) :: s2_45,zz1,zz2

    !    statement functions
    weight4252(eta41,eta42,eta51,eta52,eta54)  = eta41*eta51*( one + eta42/(eta54 + eta41 + eta52) & 
         + eta52/(eta54 + eta42 + eta51) )
    
    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S5Lim)          !3
    call initialize_config(TCLim)          !5
    call initialize_config(TCS5Lim)        !7
    call initialize_config(C54Lim)          !9
    call initialize_config(C54S5Lim)        !11
    call initialize_config(TCC54Lim)        !13
    call initialize_config(TCC54S5Lim)      !15

    MV =  sqrt(MV2)

    x1 = yr(1) 
    x2 = yr(2) 
    x3 = yr(3) 
    x4 = yr(5) 

    if (withmax) then
       xx3 = x3
    else
       xx3 = zero
    endif

    lambda = sin( half * pi * yr(6) )**2

    x4p = one - half*x4 

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = -one + two*x3 
    sin4 = sqrt(one-cos4**2)

    cos5 = -one + two*x3*x4p
    sin5 = sqrt(one-cos5**2) 


    if (x4p.eq.zero) then 
       if(one-two*lambda .ge. zero) then
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif
    elseif(2.0_dp*(one - x4p)**2/NF1(x3,x4p,lambda)-two*(one+x4p*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    eta42 = x3
    eta52 = x3*x4p
    eta41 = (one-eta42) 
    eta51 = (one - eta52)
    eta54 = x3*x4**2/four/NF1(x3,x4p,lambda)


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (HardProc%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) HardProc%Flag = .true.
    else
       HardProc%flag = .true.
    endif


    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,x4p,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart  & 
            * spart/MV2

    endif



    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 FLM + C54 S5 FLM   ! angles the same as in hard, do not regenerate; weight is the same 
    !-------------------------------------------------------------------------------------------------------------

    ! ---------- SOFT 5------------
    eta42 = x3
    eta52 = x3*x4p
    eta41 = (one-eta42) 
    eta51 = (one - eta52)
    eta54 = x3*x4**2/four/NF1(x3,x4p,lambda)
    spart = MV2/(one - x1)

    if (spart.gt.sh) S5Lim%flag  = .true.
    sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41)
    Deltad = one - x1*(one - eta42)

    if ( (S5Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) S5Lim%Flag = .true.
    else 
       S5Lim%flag = .true.
    endif

    if (S5Lim%flag .eqv. .false.) then 

       S5Lim%PartFrac= (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p4 

       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = p4 
       S5Lim%LimMom(:,4) = p5 

       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,x4p,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ------  C54 + S5 ---------

       eta42 = x3
       eta52 = x3
       eta41 = (one - eta42) 
       eta51 = (one - eta52)
       eta54 = zero

       s45 = spart *x1*x2 * x3 * x4**2/four/NF1(x3,one,lambda)
       z5 = x2/x1

       C54S5Lim%KinInv = (/zero,zero,s45,zero,z5,zero/)

       C54S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif



    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    ! --------- TRIPLE COLLINEAR -------

    eta42 = zero
    eta52 = zero
    eta41 = one
    eta51 = one
    eta54 = zero 
    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( TCLim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCLim%Flag = .true.
    else 
       TCLim%flag  = .true.
    endif

    if (TCLim%flag  .eqv. .false.) then 

       TCLim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, -one  /)
       p5 = half*sqrts*x2 *(/ one, zero, zero, -one  /)


       s24 = spart *x1*x3 
       s25 = spart *x2*x3*x4p
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       ztot = x1+x2 - one 
       z4 = x1/ztot 
       z5 = x2/ztot

       TCLim%KinInv= (/-s24,-s25,s45,z4,z5,zero/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)

       TCLim%AmpMom(:,1) = p1 
       TCLim%AmpMom(:,2) = p2- p4 - p5

       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

       ! -------   TC  + C54 -------------

       s24 = spart *x1*x3 
       s25 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
       s2_45 = s24 * (x1+x2)/x1
       zz1 = x2/(x1+x2)
       zz2 = (x1+x2)/(x1+x2-one)

       TCC54Lim%KinInv = (/-s2_45,zero,s45,zz1,zz2,zero/)

       TCC54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    ! -------- TRIPLE COLLINEAR + SOFT 5

    spart = MV2/(one - x1) 
    if (spart.gt.sh) TCS5Lim%flag = .true. 

    sqrts = sqrt(spart)
    Deltan = one - x1*(one - eta41) 
    Deltad = one - x1*(one - eta42)

    if ( ( TCS5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCS5Lim%Flag = .true.
    else 
       TCS5Lim%flag = .true.
    endif

    if (TCS5Lim%flag .eqv. .false.) then 

       TCS5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, -one  /)


       z4 = x1/(x1-one)
       z5 = x2/(x1-one)

       s24 = spart *x1*x3 
       s25 = spart *x2*x3*x4p
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       TCS5Lim%KinInv = (/-s24,-s25,s45,z4,z5,zero/)
       pVV = p1 + p2 - p4 
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1 
       TCS5Lim%AmpMom(:,2) = p2- p4 


       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ------   TRIPLE COLLINEAR KINEMATICS  + C54 + SOFT5 -------


       z4 = x1/(x1-one)
       z5 = x2/x1

       s24 = spart *x1*x3 
       s25 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)

       TCC54S5Lim%KinInv = (/-s24,-s25,s45,z4,z5,zero/)

       TCC54S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C54 FLM 
    !-------------------------------------------------------------------------------------------------------------

    cos4 = -one + two*x3 
    sin4 = sqrt(one-cos4**2)

    eta42 = x3
    eta52 = x3
    eta41 = (one - eta42) 
    eta51 = (one - eta52)
    eta54 = zero

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C54Lim%flag = .true.
    if (C54Lim%flag.eqv. .false.) spart = MV2/Deltas
    if (C54Lim%flag.eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( C54Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C54Lim%Flag = .true.
    else 
       C54Lim%flag = .true.
    endif

    if (C54Lim%flag .eqv. .false.) then 

       C54Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)  
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)


       pVV = p1 + p2 - p4 - p5

       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)

       z5 = x1/(x1 + x2) 

       C54Lim%KinInv = (/zero,zero,s45,zero,z5,zero/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C54Lim)


       C54Lim%AmpMom(:,1) = p1
       C54Lim%AmpMom(:,2) = p2
       C54Lim%AmpMom(:,3) = p4 + p5 


       C54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif

  end subroutine kinematics_nnlo_4252b

  !-- we factor out an extra 1/(8*pisq)**2 from the phase space
  !-- which multiplied by gs**4 in the ME reconstructs asontwopi**2
  subroutine kinematics_nnlo_4252c(MV2,Y,yr, &
       HardProc,S5Lim,TCLim,TCS5Lim,  &
       C4Lim,C4S5Lim,TCC4Lim,TCC4S5Lim)
    real(dp), intent(in) ::  MV2, Y, yr(9)
    type(KinConfig) ::  HardProc,S5Lim,TCLim,TCS5Lim, C4Lim,C4S5Lim,TCC4Lim,TCC4S5Lim
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4) 
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5,  pVV(4)
    real(dp) :: xx3,x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w1415, s14, s15, s45, z4, z5, ztot, weight4252
    real(dp) :: spart,sqrts, s12, s25, z45, yy, Deltas, s24  
    real(dp) :: s5_24,zz1,zz2,kt1kt2sq, ncoll(4), nperp(4)


    !    statement functions
    weight4252(eta41,eta42,eta51,eta52,eta54)  = eta41*eta51*( one + eta42/(eta54 + eta41 + eta52) & 
         + eta52/(eta54 + eta42 + eta51) )

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S5Lim)          !3
    call initialize_config(TCLim)          !5
    call initialize_config(TCS5Lim)        !7
    call initialize_config(C4Lim)          !9
    call initialize_config(C4S5Lim)        !11
    call initialize_config(TCC4Lim)        !13
    call initialize_config(TCC4S5Lim)      !15

    MV =  sqrt(MV2)

    x1 =   yr(1) 
    x2 =   yr(2) 
    x3 = yr(3) 
    x4 = yr(5) 

    if (withmax) then
       xx3 = x3
    else
       xx3 = zero
    endif

    lambda = sin( half * pi * yr(6) )**2


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = -one + x3*x4 
    sin4 = sqrt(one-cos4**2)

    cos5 = -one + two*x3 
    sin5 = sqrt(one-cos5**2) 


    if (x4.eq.zero) then 
       if(one-two*lambda .ge. zero) then
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif
    elseif(2.0_dp*(one - half*x4)**2/NF1(x3,half*x4,lambda)-two*(one+half*x4*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    eta42 = half*x3*x4
    eta52 = x3
    eta41 = (one-eta42) 
    eta51 = (one - eta52)
    eta54 = x3*(one-half*x4)**2/NF1(x3,half*x4,lambda)


    !  incoming momenta 

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (HardProc%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) HardProc%Flag = .true.
    else
       HardProc%flag = .true.
    endif


    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2


    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 FLM   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------

    eta42 = half*x3*x4
    eta52 = x3
    eta41 = (one-eta42) 
    eta51 = (one - eta52)
    eta54 = x3*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    spart = MV2/(one - x1)
    if (spart.gt.sh) S5Lim%flag = .true.
    sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41)
    Deltad = one - x1*(one - eta42)

    if ( (S5Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) S5Lim%Flag = .true.
    else 
       S5Lim%flag = .true.
    endif

    if (S5Lim%flag .eqv. .false.) then 

       S5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p4

       S5Lim%LimMom(:,1) =  p1
       S5Lim%LimMom(:,2) =  p2
       S5Lim%LimMom(:,3) =  p4
       S5Lim%LimMom(:,4) =  p5 


       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif



    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    ! --------- TRIPLE COLLINEAR -------

    eta42 = zero
    eta52 = zero
    eta41 = one
    eta51 = one
    eta54 = zero 


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( TCLim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCLim%Flag = .true.
    else 
       TCLim%flag = .true.
    endif

    if (TCLim%flag .eqv. .false.) then 

       TCLim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, -one  /)
       p5 = half*sqrts*x2*(/one,zero,zero,-one /)

       s24 = spart *x1*(half*x3*x4) 
       s25 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

       ztot = x1 +x2 - one 
       z4 = x1/ztot 
       z5 = x2/ztot

       TCLim%KinInv = (/-s24,-s25,s45,z4,z5,zero/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)

       TCLim%AmpMom(:,1) = p1 
       TCLim%AmpMom(:,2) = p2 - p4 - p5 

       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            * spart/MV2


       ! ---------   TC  + C42 ----------------

       s24 = spart * x1 *(half * x3 * x4)  
       s25 = spart * x2* x3
       s5_24 = s25 * (one-x1)
       
       zz1 = x1/(x1-one)
       zz2 = (one-x1)/(one-x1-x2)

       kt1kt2sq = (-one+two*lambda)**2 * (one-xx3)

       TCC4Lim%KinInv = (/-s24,-s5_24,zero,zz1,zz2,kt1kt2sq/)

       TCC4Lim%wgt  = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    ! -------------- TRIPLE COLLINEAR KINEMATICS + SOFT 5  -------------     

    spart = MV2/(one - x1) 
    if (spart.gt.sh) TCS5Lim%flag = .true.

    sqrts = sqrt(spart)
    Deltan = one - x1*(one - eta41) 
    Deltad = one - x1*(one - eta42)

    if ( ( TCS5Lim%flag  .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCS5Lim%Flag = .true.
    else 
       TCS5Lim%flag  = .true.
    endif

    if (TCS5Lim%flag  .eqv. .false.) then 

       TCS5Lim%PartFrac = (/xi1,xi2/)

       s24 = spart *x1*(half*x3*x4) 
       s25 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)


       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, -one  /)
       p5 = half*sqrts*x2*(/one,zero,zero,-one /)

       z4 = x1/(x1-one)
       z5 = x2/(x1-one)

       TCS5Lim%KinInv = (/-s24,-s25,s45,z4,z5,zero/)

       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1  
       TCS5Lim%AmpMom(:,2) = p2 - p4

       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            * spart/MV2


       ! -------   TRIPLE COLLINEAR KINEMATICS  + C4 + SOFT5 -----------

       s24 = spart *x1*(half * x3 * x4)  
       s25 = spart *x2 * x3
       s45 = spart *x1*x2*x3                 ! /NF1(xx3,zero,lambda) ->1

       z4 = x1/(x1 - one) 
       z5 = x2

       TCC4S5Lim%KinInv = (/-s24,-s25,s45,z4,z5,zero/)


       TCC4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1 * x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2



       ! ----------  C42 + S5 

       eta42 = zero
       eta52 = x3
       eta41 = one 
       eta51 = (one - eta52)
       eta54 = x3/NF1(x3,zero,lambda)

       z4 =  x1/(x1 - one) 
       s12 = two*scr(p1,p2) 




       s24 = spart *x1*(half*x3*x4) 
       s25 = spart *x2*x3
       s15 = spart *x2*(one-x3)

       C4S5Lim%KinInv = (/s24,-s25,zero,z4,zero,zero/)
       C4S5Lim%EikInv = (/s12,s15,s25/)

       ! check: do we need this?

       MV2aux = scr(pVV,pVV) 
       MV =  sqrt(MV2aux )


       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4S5Lim)


       C4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart & 
            * spart/MV2

    endif




    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C41 FLM 
    !-------------------------------------------------------------------------------------------------------------

    cos5 = -one + two*x3 
    sin5 = sqrt(one-cos5**2) 

    if(one-two*lambda .ge. zero) then
       sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif


    eta42 = zero
    eta52 = x3
    eta41 = one 
    eta51 = (one - eta52)
    eta54 = x3/NF1(x3,zero,lambda)


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C4Lim%flag  = .true.
    if (C4Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( C4Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C4Lim%Flag = .true.
    else
       C4Lim%flag = .true.
    endif

    if ( C4Lim%flag .eqv. .false.) then 

       C4Lim%PartFrac = (/xi1,xi2/)
       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1*(/one, zero, zero, -one  /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       nperp = (/zero,cos(phi4),sin(phi4),zero/)
       ncoll = (/zero,zero,zero,-one/)

       C4Lim%LimMom(:,1) = nperp
       C4Lim%LimMom(:,2) = ncoll

       s24 = spart*x1*(half*x3*x4)
       z4 = x1/( x1-one ) 

       C4Lim%KinInv = (/zero,-s24,zero,zero,z4,zero/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4Lim)

       C4Lim%AmpMom(:,1) = p1 
       C4Lim%AmpMom(:,2) = p2-p4
       C4Lim%AmpMom(:,3) = p5 


       C4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif

  end subroutine kinematics_nnlo_4252c

  !-- we factor out an extra 1/(8*pisq)**2 from the phase space
  !-- which multiplied by gs**4 in the ME reconstructs asontwopi**2
  subroutine kinematics_nnlo_4252d(MV2,Y,yr,&
       HardProc,S5Lim,TCLim,TCS5Lim,C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim)
    real(dp), intent(in) ::  yr(9), MV2, Y
    type(KinConfig)      ::  HardProc,S5Lim,TCLim,TCS5Lim,&
         C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim
    !---
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV,  p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: xx3,x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w1415, s14, s15, s45, z4, z5, ztot, weight4252
    real(dp) :: spart,sqrts, x4p, nperp(4), s12, s24, s25, z45, ncoll(1:4), yy, Deltas 
    real(dp) :: s2_45,zz1,zz2

    !    statement functions
    weight4252(eta41,eta42,eta51,eta52,eta54)  = eta41*eta51*( one + eta42/(eta54 + eta41 + eta52) & 
         + eta52/(eta54 + eta42 + eta51) )

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S5Lim)          !3
    call initialize_config(TCLim)          !5
    call initialize_config(TCS5Lim)        !7
    call initialize_config(C54Lim)          !9
    call initialize_config(C54S5Lim)        !11
    call initialize_config(TCC54Lim)        !13
    call initialize_config(TCC54S5Lim)      !15

    MV =  sqrt(MV2)

    x1 =   yr(1) 
    x2 =   yr(2) 
    x3 = yr(3) 
    x4 = yr(5) 

    if (withmax) then
       xx3 = x3
    else
       xx3 = zero
    endif

    lambda = sin( half * pi * yr(6) )**2

    x4p = one - half*x4 

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = -one + two*x3*x4p 
    sin4 = sqrt(one-cos4**2)

    cos5 = -one + two*x3
    sin5 = sqrt(one-cos5**2) 


    if (x4p.eq.zero) then 
       if(one-two*lambda .ge. zero) then
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif
    elseif(2.0_dp*(one - x4p)**2/NF1(x3,x4p,lambda)-two*(one+x4p*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    eta42 = x3*x4p
    eta52 = x3
    eta41 = (one-eta42) 
    eta51 = (one - eta52)
    eta54 = x3*x4**2/four/NF1(x3,x4p,lambda)

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (HardProc%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) HardProc%Flag = .true.
    else
       HardProc%flag = .true.
    endif


    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,x4p,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 FLM   ! angles the same as in hard, do not regenerate; weight is the same 
    !-------------------------------------------------------------------------------------------------------------

    eta42 = x3*x4p
    eta52 = x3
    eta41 = (one-eta42) 
    eta51 = (one - eta52)
    eta54 = x3*x4**2/four/NF1(x3,x4p,lambda)

    spart = MV2/(one - x1)
    if (spart.gt.sh) S5Lim%flag = .true.
    sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41)
    Deltad = one - x1*(one - eta42)

    if ( (S5Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) S5Lim%Flag = .true.
    else 
       S5Lim%flag = .true.
    endif

    if (S5Lim%flag .eqv. .false.) then 

       S5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,S5Lim)
       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p4 

       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = p4 
       S5Lim%LimMom(:,4) = p5 

       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,x4p,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2 
    endif


    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C54 + S5 
    !-------------------------------------------------------------------------------------------------------------


    eta42 = x3
    eta52 = x3
    eta41 = (one - eta42) 
    eta51 = (one - eta52)
    eta54 = zero 

    phi4 = twopi*yr(4) 
    cos4 = -one + two*x3 
    sin4 = sqrt(one-cos4**2)

    spart = MV2/(one - x1)
    if (spart.gt.sh) C54S5Lim%flag = .true.

    sqrts = sqrt(spart)
    Deltan = one - x1*(one - eta41)
    Deltad = one - x1*(one - eta42)

    if ( ( C54S5Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then
       yy = Y - half*log(Deltan/Deltad)
       xi1 = sqrt(spart/sh)*exp(yy)
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C54S5Lim%Flag = .true.
    else
       C54S5Lim%flag = .true.
    endif

    if (C54S5Lim%flag .eqv. .false.) then

       C54S5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/)
       p2 =half*sqrts*(/one,zero,zero,-one/)
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)


       z5 = x2/x1
       s45 = spart *x1  * x2 * x3 * x4**2/four/NF1(x3,one,lambda)

       C54S5Lim%KinInv = (/zero,zero,s45,zero,z5,zero/)

       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C54S5Lim)
       C54S5Lim%AmpMom(:,1) = p1
       C54S5Lim%AmpMom(:,2) = p2
       C54S5Lim%AmpMom(:,3) = p4



       C54S5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2
    endif



    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    eta42 = zero
    eta52 = zero
    eta41 = one
    eta51 = one
    eta54 = zero 

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( TCLim%flag .eqv. .false.) .and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCLim%Flag = .true.
    else 
       TCLim%flag = .true.
    endif

    if (TCLim%flag .eqv. .false.) then 

       TCLim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, -one  /)
       p5 = half*sqrts*x2*(/one,zero,zero, -one /)

       s24 = spart *x1*x3*x4p 
       s25 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       ztot = x1 + x2 - one 
       z4 = x1/ztot 
       z5 = x2/ztot

       TCLim%KinInv = (/-s24,-s25,s45,z4,z5,zero/)


       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)


       TCLim%AmpMom(:,1) = p1 
       TCLim%AmpMom(:,2) = p2- p4 - p5


       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

       ! --------------  TC  + C54 ----------

       s24 = spart *x1*x3 
       s25 = spart *x2*x3
       s45 = spart *x1*x2*x3*x4**2/four/NF1(xx3,one,lambda)
       s2_45 = s24 * (x1+x2)/x1

       zz1 = x2/(x1+x2)
       zz2 = (x1+x2)/(x1+x2-one)

       TCC54Lim%KinInv = (/-s2_45,zero,s45,zz1,zz2,zero/)

       TCC54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    spart = MV2/(one - x1) 
    if (spart.gt.sh) TCS5Lim%flag = .true.

    sqrts = sqrt(spart)
    Deltan = one - x1*(one - eta41) 
    Deltad = one - x1*(one - eta42)

    if ( ( TCS5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCS5Lim%Flag = .true.
    else 
       TCS5Lim%flag = .true.
    endif

    if (TCS5Lim%flag .eqv. .false.) then 

       TCS5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, -one  /)

       z4 = x1/(x1-one)
       z5 = x2/(x1-one)

       s24 = spart *x1*x3*x4p 
       s25 = spart * x2*x3
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       TCS5Lim%KinInv = (/-s24,-s25,s45,z4,z5,zero/)

       pVV = p1 + p2 - p4 
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1 
       TCS5Lim%AmpMom(:,2) = p2- p4 

       TCS5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2



       ! --------   TRIPLE COLLINEAR KINEMATICS  + C54 + SOFT5 ------------

       s24 = spart *x1*x3 
       s25 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)

       z4 = x1/(x1-one)
       z5 = x2/x1

       TCC54S5Lim%KinInv = (/-s24,-s25,s45,z4,z5,zero/)

       TCC54S5Lim%wgt = one/8.0_dp/pi * kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C54 FLM 
    !-------------------------------------------------------------------------------------------------------------


    phi4 = twopi*yr(4) 
    cos4 = -one + two*x3 
    sin4 = sqrt(one-cos4**2)

    eta42 = x3
    eta52 = x3
    eta41 = (one - eta42) 
    eta51 = (one - eta52)
    eta54 = zero 

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C54Lim%flag = .true.
    if (C54Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C54Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( C54Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C54Lim%Flag = .true.
    else 
       C54Lim%flag = .true.
    endif

    if (C54Lim%flag .eqv. .false.) then 

       C54Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)


       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)
       z5 = x1/(x1 + x2) 

       C54Lim%KinInv = (/zero,zero,s45,zero,z5,zero/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C54Lim)

       C54Lim%AmpMom(:,1) = p1
       C54Lim%AmpMom(:,2) = p2
       C54Lim%AmpMom(:,3) = p4 + p5 

       C54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart & 
            * spart/MV2

    endif

  end subroutine kinematics_nnlo_4252d

  !-- we factor out an extra 1/(8*pisq)**2 from the phase space
  !-- which multiplied by gs**4 in the ME reconstructs asontwopi**2
  subroutine kinematics_nnlo_4251(MV2,Y,yr,&
       HardProc,S5Lim,C4Lim,C4S5Lim,&
       C5Lim,C5S5Lim,C4C5Lim,C4C5S5Lim) 
    real(dp), intent(in) ::  yr(9), MV2, Y
    type(KinConfig)      :: HardProc,S5Lim,C4Lim,C4S5Lim, C5Lim, C5S5Lim, C4C5Lim, C4C5S5Lim
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: xx3,x1, x2,x3, x4, lambda,  cos4, sin4, phi4, phi5, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w1415, s14, s15, s45, z4, z5, ztot, weight4251
    real(dp) :: spart,sqrts, s25, s12, s24, yy, Deltas

    !    statement functions
    weight4251(eta41,eta42,eta51,eta52,eta54)  = eta41*eta52*eta54/(eta54 + eta51 + eta42)

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)         !1
    call initialize_config(S5Lim)            !2
    call initialize_config(C5Lim)            !3 
    call initialize_config(C5S5Lim)          !4 
    call initialize_config(C4Lim)            !5
    call initialize_config(C4S5Lim)          !6
    call initialize_config(C4C5Lim)          !7
    call initialize_config(C4C5S5Lim)        !8

    !    MV2 and Y are generated

    MV =  sqrt(MV2)

    x1 =   yr(1) 
    x2 =   yr(2) 
    x3 =   yr(3) 
    x4 =   yr(5) 

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = -one +  two*x3 
    sin4 = sqrt(one-cos4**2)

    phi5 = twopi*yr(6)
    cos5 = one - two*x4
    sin5 = sqrt(one-cos5**2) 

    eta41 = (one - x3)
    eta51 = x4
    eta42 = x3 
    eta52 = (one - x4) 
    eta54 =   half*( one - sin4*sin5*cos(phi4-phi5) -cos4*cos5)

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (HardProc%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) HardProc%Flag = .true.
    else
       HardProc%flag = .true.
    endif


    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 



       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4251(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3: S5 FLM   ! angles the same as in hard // do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! --------------- SOFT 5 ----------


    spart = MV2/(one - x1)

    if (spart.gt.sh) S5Lim%flag = .true.
    sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41)
    Deltad = one - x1*(one - eta42)

    if ( (S5Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) S5Lim%Flag = .true.
    else 
       S5Lim%flag = .true.
    endif

    if (S5Lim%flag .eqv. .false.) then 

       S5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)


       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p4 

       ! --      note momentum swap below; this is because of the definition of the Eikonal_s5

       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = p4
       S5Lim%LimMom(:,4) = p5 

       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4251(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

       ! ---------  C51 + S5 --------------

       phi4 = twopi*yr(4) 
       cos4 = -one + two*x3 
       sin4 = sqrt(one-cos4**2)

       eta42 = x3
       eta52 = one
       eta41 = (one-x3) 
       eta51 = zero
       eta54 = (one-x3) 


       s15 = spart*x2*x4
       z5 = x2

       C5S5Lim%KinInv = (/zero,s15,zero,zero,z5,zero/)

       C5S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4251(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3:   C42  
    !-------------------------------------------------------------------------------------------------------------



    eta41 = one
    eta51 = x4 
    eta42 = zero  
    eta52 = (one - x4) 
    eta54 =   eta52


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C4Lim%flag  = .true.
    if (C4Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4Lim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( C4Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C4Lim%Flag = .true.
    else 
       C4Lim%flag = .true.
    endif

    if (C4Lim%flag .eqv. .false.) then 

       C4Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  zero, zero, -one  /)
       p5 = half*sqrts*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)


       s24 = spart *x1*x3 
       z4 = x1/(x1-one)

       C4Lim%KinInv = (/s24,zero,zero,z4,zero,zero/)

       C4Lim%LimMom(:,1) = (/zero,cos(phi4),sin(phi4),zero/) !-- nperp
       C4Lim%LimMom(:,2) = (/zero,zero,zero,-one/) !-- ndir

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4Lim)

       C4Lim%AmpMom(:,1) = p1 
       C4Lim%AmpMom(:,2) = p2 - p4
       C4Lim%AmpMom(:,3) = p5

       C4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4251(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2 

    endif



    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C42  + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------


    spart = MV2/(one - x1) 
    if (spart.gt.sh) C4S5Lim%flag = .true.

    sqrts = sqrt(spart)
    Deltan = one - x1*(one - eta41) 
    Deltad = one - x1*(one - eta42)

    if ( ( C4S5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C4S5Lim%Flag = .true.
    else 
       C4S5Lim%flag = .true.
    endif

    if (C4S5Lim%flag .eqv. .false.) then 

       C4S5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, -one  /)

       s24 = spart *x1*x3 
       s25 = spart*x2*(one-x4) 
       s15 = spart*x2*x4
       s12 = two*scr(p1,p2) 

       z4 = x1/(x1 - one)  

       C4S5Lim%KinInv = (/s24,s15,zero,z4,zero,zero/)
       C4S5Lim%EikInv = (/s12,s15,s25/)

       pVV = p1 + p2 - p4 

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4S5Lim)

       C4S5Lim%AmpMom(:,1) = p1 
       C4S5Lim%AmpMom(:,2) = p2 - p4



       C4S5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4251(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ---------   C42   + C51 + SOFT5 ------------


       eta42 = zero
       eta52 = one
       eta41 = one
       eta51 = zero
       eta54 = one           

       z4 = x1/(x1-one)
       z5 = x2

       s24 = spart *x1*x3 
       s15 = spart *x2*x4

       C4C5S5Lim%KinInv = (/-s24,-s15,zero,z4,z5,zero/)

       C4C5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4251(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: C51 FLM 
    !-------------------------------------------------------------------------------------------------------------


    phi4 = twopi*yr(4) 
    cos4 = -one + two*x3 
    sin4 = sqrt(one-cos4**2)

    eta42 = x3
    eta52 = one
    eta41 = (one-x3) 
    eta51 = zero
    eta54 = (one-x3) 


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C5Lim%flag = .true.
    if (C5Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C5Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( C5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C5Lim%Flag = .true.
    else 
       C5Lim%flag = .true.
    endif

    if (C5Lim%flag .eqv. .false.) then 

       C5Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1*(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, zero, zero, one /)

       s15 = spart*x2*x4
       z5 = one/(one - x2) 

       C5Lim%KinInv = (/s15,zero,zero,zero,z5,zero/)

       C5Lim%LimMom(:,1) = (/zero,cos(phi5),sin(phi5),zero/) !-- nperp
       C5Lim%LimMom(:,2) = (/zero,zero,zero,one/) !-- ndir

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C5Lim)

       C5Lim%AmpMom(:,1) = p1 - p5 
       C5Lim%AmpMom(:,2) = p2 
       C5Lim%AmpMom(:,3) = p4 

       C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4251(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif




    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C42  + C51
    !-------------------------------------------------------------------------------------------------------------

    eta42 = zero
    eta52 = one
    eta41 = one
    eta51 = zero
    eta54 = one  

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C4C5Lim%flag = .true.
    if (C4C5Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4C5Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( C4C5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C4C5Lim%Flag = .true.
    else 
       C4C5Lim%flag = .true.
    endif

    if (C4C5Lim%flag .eqv. .false.) then 

       C4C5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, -one  /)
       p5 = half*sqrts*x2*(/one,zero,zero,one /)

       s24 = spart *x1*x3 
       s15 = spart *x2*x4
       z4 = x1/(x1-one)
       z5 = one/(one - x2) 

       C4C5Lim%KinInv = (/-s24,-s15,zero,z4,z5,zero/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4C5Lim)

       C4C5Lim%AmpMom(:,1) = p1 - p5 
       C4C5Lim%AmpMom(:,2) = p2 - p4

       C4C5Lim%LimMom(:,1) = (/zero,cos(phi4),sin(phi4),zero/) !-- nperp
       C4C5Lim%LimMom(:,2) = (/zero,zero,zero,-one/) !-- ndir
       C4C5Lim%LimMom(:,3) = (/zero,cos(phi5),sin(phi5),zero/) !-- nperp
       C4C5Lim%LimMom(:,4) = (/zero,zero,zero,one/) !-- ndir

       C4C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4251(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif

  end subroutine kinematics_nnlo_4251

  !-- symmetric sectors
  
  !-- we factor out an extra 1/(8*pisq)**2 from the phase space
  !-- which multiplied by gs**4 in the ME reconstructs asontwopi**2
  subroutine kinematics_nnlo_4151a(MV2,Y,yr,&
       HardProc,S5Lim,TCLim,TCS5Lim,&
       C5Lim,C5S5Lim,TCC5Lim,TCC5S5Lim)       
    real(dp), intent(in) ::  MV2,Y, yr(9)
    type(KinConfig)      ::  HardProc,S5Lim,TCLim,TCS5Lim,&
         C5Lim,C5S5Lim,TCC5Lim,TCC5S5Lim
    !---
    real(dp) :: p1(4), p2(4), p4(4), p5(4)
    real(dp) :: MV, kallenF, sin5, cos5, pVV(4)
    real(dp) :: xx3,x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4151, s24, s25, s45, z4, z5, ztot, weight4151, z45, ncoll(1:4), nperp(1:4)
    real(dp) :: spart,sqrts, yy, Deltas, acosphi45, asinphi45, s15, s14,s4_15,zz1,zz2,kt1kt2sq

    !    statement functions
    weight4151(eta41,eta42,eta51,eta52,eta54)  =  eta42*eta52*( one + eta41/(eta54 + eta42 + eta51) & 
         + eta51/(eta54 + eta41 + eta52) ) 

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S5Lim)          !2
    call initialize_config(TCLim)          !3
    call initialize_config(TCS5Lim)        !4
    call initialize_config(C5Lim)          !5
    call initialize_config(C5S5Lim)        !6
    call initialize_config(TCC5Lim)        !7
    call initialize_config(TCC5S5Lim)      !8

    !    MV2 and Y are generated

    MV =  sqrt(MV2)

    x1 =   yr(1) 
    x2 =   yr(2) 
    x3 =   yr(3) 
    x4 =   yr(5) 

    if (withmax) then
       xx3 = x3
    else
       xx3 = zero
    endif

    lambda = sin( half * pi * yr(6) )**2


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = one - two*x3 
    sin4 = sqrt(one-cos4**2)

    cos5 = one - x3*x4 
    sin5 = sqrt(one-cos5**2) 


    if (x4.eq.zero) then 
       if(one-two*lambda .ge. zero) then
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif
    elseif(2.0_dp*(one - half*x4)**2/NF1(x3,half*x4,lambda)-two*(one+half*x4*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    eta42 = one-x3
    eta52 = one-half*x3*x4
    eta41 = x3 
    eta51 = half*x3*x4
    eta54 = x3*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (HardProc%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) HardProc%Flag = .true.
    else
       HardProc%flag = .true.
    endif


    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif

    ! --------- SOFT 5 ---------

    eta41 = x3
    eta51 = half*x3*x4
    eta42 = (one-x3)
    eta52 = (one - half*x3*x4)
    eta54 = x3*(one-half*x4)**2/NF1(x3,half*x4,lambda)


    spart = MV2/(one - x1)
    if (spart.gt.sh) S5Lim%flag = .true.
    sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41)
    Deltad = one - x1*(one - eta42)

    if ( (S5Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) S5Lim%Flag = .true.
    else 
       S5Lim%flag = .true.
    endif

    if (S5Lim%flag .eqv. .false.) then 

       S5Lim%PartFrac = (/xi1,xi2/)


       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)

       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p4 

       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = p4
       S5Lim%LimMom(:,4) = p5 

       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2   


       ! kin inv and weight for C5 +  S5

       eta41 = x3
       eta51 = zero
       eta42 = (one-x3)
       eta52 = one
       eta54 = x3/NF1(x3,zero,lambda)

       s15 = spart*x2*(half*x3*x4)
       z5 = x2
       C5S5Lim%KinInv = (/zero,-s15,zero,zero,z5,zero/)

       C5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS 
    !-------------------------------------------------------------------------------------------------------------

    ! --------- TRIPLE COLLINEAR -------

    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero 


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( TCLim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCLim%Flag = .true.
    else 
       TCLim%flag = .true.
    endif

    if (TCLim%flag .eqv. .false.) then 

       TCLim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)
       p5 = half*sqrts*x2*(/one,zero,zero, one /)

       s14 = spart *x1*x3 
       s15 = spart *x2*(half*x3*x4)
       s45 = spart *x1*x2*x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

       ztot = x1 + x2 - one 
       z4 = x1/ztot 
       z5 = x2/ztot

       TCLim%KinInv = (/-s14,-s15,s45,z4,z5,zero/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)

       TCLim%AmpMom(:,1) = p1- p4 - p5
       TCLim%AmpMom(:,2) = p2 

       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2


       ! --------- TRIPLE COLLINEAR + C5  -------

       s14 = spart *x1*x3 
       s15 = spart *x2*(half*x3*x4)
!       s45 = spart *x1*x2*x3/NF1(xx3,zero,lambda)
       s4_15 = s14 * (one - x2)

!       z45 = x1/(x1 + x2 - one)
!       z5 = x2/(x2 - one)
       zz1 = x2/(x2-one)
       zz2 = (one-x2)/(one-x1-x2)
       kt1kt2sq = (-one+two*lambda)**2 * (one-xx3)


       TCC5Lim%KinInv = (/-s4_15,-s15,zero,zz1,zz2,kt1kt2sq/)   ! the signs here are to compline with the original definition 
       ! q1 -> q g4 g5

       TCC5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    ! -------- TRIPLE COLLINEAR + SOFT 5


    spart = MV2/(one - x1) 
    if (spart.gt.sh) TCS5Lim%flag = .true. 

    sqrts = sqrt(spart)
    Deltan = one - x1*(one - eta41) 
    Deltad = one - x1*(one - eta42)

    if ( ( TCS5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCS5Lim%Flag = .true.
    else 
       TCS5Lim%flag = .true.
    endif

    if (TCS5Lim%flag .eqv. .false.) then 

       TCS5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)
       p5 = half*sqrts*x2*(/one,zero,zero,one /)


!       z4 = one/x1
!       z5 = x2/x1
       z4 = x1/(x1-one)
       z5 = x2/(x1-one)


       s14 = spart *x1*x3 
       s15 = spart *x2*(half*x3*x4)
       s45 = spart *x1 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

       TCS5Lim%KinInv = (/-s14,-s15,s45,z4,z5,zero/)

       pVV = p1 + p2 - p4 

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1 - p4
       TCS5Lim%AmpMom(:,2) = p2 

       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            *spart/MV2


       ! -------   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT5 ---------
!       z4 = x1/(x1-one) 
       z4 = x1
       z5 = x2

       s14 = spart *x1*x3 
       s15 = spart *x2*(half*x3*x4)
       s45 = spart *x1*x2*x3/NF1(xx3,zero,lambda)

       TCC5S5Lim%KinInv = (/-s14,-s15,s45,z4,z5,zero/)

       TCC5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C5 FLM 
    !-------------------------------------------------------------------------------------------------------------


    eta41 = x3
    eta51 = zero
    eta42 = (one-x3)
    eta52 = one
    eta54 = x3/NF1(x3,zero,lambda)

    phi4 = twopi*yr(4) 
    cos4 = one - two*x3 
    sin4 = sqrt(one-cos4**2)


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C5Lim%flag = .true.
    if (C5Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C5Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) -  x2*(one - eta51)
    Deltad = one - x1*(one - eta42) -  x2*(one - eta52)  

    if ( ( C5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C5Lim%Flag = .true.
    else 
       C5Lim%flag = .true.
    endif

    if (C5Lim%flag .eqv. .false.) then 

       C5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, zero, zero, one /)

       s15 = spart*x2*(half*x3*x4)
       z5 = x2/(x2 - one) 

       C5Lim%KinInv = (/zero,-s15,zero,zero,z5,zero/)


       ncoll = (/ zero,  zero, zero , one  /)

       if (lambda.gt.half) then 
          acosphi45 = -(one-two*lambda) 
          asinphi45 = two*sqrt(lambda*(one-lambda))
       else 
          acosphi45 = -(one-two*lambda) 
          asinphi45 = two*sqrt(lambda*(one-lambda))
       endif


       nperp  = (/ zero, cos(phi4)*acosphi45 - sin(phi4)*asinphi45 , & 
            sin(phi4)*acosphi45 + cos(phi4)*asinphi45, zero /)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C5Lim)

       C5Lim%AmpMom(:,1) = p1 - p5
       C5Lim%AmpMom(:,2) = p2 
       C5Lim%AmpMom(:,3) = p4 


       C5Lim%LimMom(:,1) = nperp
       C5Lim%LimMom(:,2) = ncoll


       C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif


  end subroutine kinematics_nnlo_4151a


  !-- we factor out an extra 1/(8*pisq)**2 from the phase space
  !-- which multiplied by gs**4 in the ME reconstructs asontwopi**2
  subroutine kinematics_nnlo_4151b(MV2,Y,yr,&
       HardProc,S5Lim,TCLim,TCS5Lim,&
       C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim) 
    real(dp), intent(in) ::  yr(9), MV2, Y
    type(KinConfig)      ::  HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C54Lim,C54DSLim,C54DSS5Lim,C54S5Lim,TCC54Lim,TCC54DSLim,TCC54DSS5Lim,TCC54S5Lim
    !---
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: xx3,x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, Deltas, w4151, s14, s15, s45, z4, z5, ztot, weight4151
    real(dp) :: spart,sqrts, x4p, nperp(4), s12, s24, z45, ncoll(1:4), yy, s25,s1_45,zz1,zz2

    !    statement functions
    weight4151(eta41,eta42,eta51,eta52,eta54)  =  eta42*eta52*( one + eta41/(eta54 + eta42 + eta51) & 
         + eta51/(eta54 + eta41 + eta52) ) 

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S5Lim)          !3
    call initialize_config(TCLim)          !5
    call initialize_config(TCS5Lim)        !7
    call initialize_config(C54Lim)          !9
    call initialize_config(C54S5Lim)        !11
    call initialize_config(TCC54Lim)        !13
    call initialize_config(TCC54S5Lim)      !15

    !    MV2 and Y are generated

    MV =  sqrt(MV2)

    x1 = yr(1) 
    x2 = yr(2) 
    x3 = yr(3) 
    x4 = yr(5) 

    if (withmax) then
       xx3 = x3
    else
       xx3 = zero
    endif

    lambda = sin( half * pi * yr(6) )**2

    x4p = one - half*x4 

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = one - two*x3 
    sin4 = sqrt(one-cos4**2)

    cos5 = one - two*x3*x4p
    sin5 = sqrt(one-cos5**2) 


    if (x4p.eq.zero) then 
       if(one-two*lambda .ge. zero) then
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif
    elseif(2.0_dp*(one - x4p)**2/NF1(x3,x4p,lambda)-two*(one+x4p*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    eta41 = x3
    eta51 = x3*x4p
    eta42 = (one-eta41) 
    eta52 = (one - eta51)
    eta54 = x3*x4**2/four/NF1(x3,x4p,lambda)


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (HardProc%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) HardProc%Flag = .true.
    else
       HardProc%flag = .true.
    endif


    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,x4p,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart  & 
            * spart/MV2

    endif



    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 FLM + C54 S5 FLM   ! angles the same as in hard, do not regenerate; weight is the same 
    !-------------------------------------------------------------------------------------------------------------

    ! ---------- SOFT 5------------
    eta41 = x3
    eta51 = x3*x4p
    eta42 = (one-eta41) 
    eta52 = (one - eta51)
    eta54 = x3*x4**2/four/NF1(x3,x4p,lambda)
    spart = MV2/(one - x1)

    if (spart.gt.sh) S5Lim%flag  = .true.
    sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41)
    Deltad = one - x1*(one - eta42)

    if ( (S5Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) S5Lim%Flag = .true.
    else 
       S5Lim%flag = .true.
    endif

    if (S5Lim%flag .eqv. .false.) then 

       S5Lim%PartFrac= (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p4 

       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = p4 
       S5Lim%LimMom(:,4) = p5 

       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,x4p,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ------  C54 + S5 ---------

       eta41 = x3
       eta51 = x3
       eta42 = (one - eta41) 
       eta52 = (one - eta51)
       eta54 = zero

       s45 = spart *x1*x2 * x3 * x4**2/four/NF1(x3,one,lambda)
       z5 = x2/x1

       C54S5Lim%KinInv = (/zero,zero,s45,zero,z5,zero/)

       C54S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif



    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    ! --------- TRIPLE COLLINEAR -------

    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero 
    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( TCLim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCLim%Flag = .true.
    else 
       TCLim%flag  = .true.
    endif

    if (TCLim%flag  .eqv. .false.) then 

       TCLim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)
       p5 = half*sqrts*x2 *(/ one, zero, zero, one  /)


       s14 = spart *x1*x3 
       s15 = spart *x2*x3*x4p
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       ztot = x1+x2 - one 
       z4 = x1/ztot 
       z5 = x2/ztot

       TCLim%KinInv= (/s14,s15,s45,z4,z5,zero/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)

       TCLim%AmpMom(:,1) = p1 - p4 - p5 
       TCLim%AmpMom(:,2) = p2

       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

       ! -------   TC  + C54 -------------

       s14 = spart *x1*x3 
       s15 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
       s1_45 = s14 * (x1+x2)/x1
       zz1 = x2/(x1+x2)
       zz2 = (x1+x2)/(x1+x2-one)


       TCC54Lim%KinInv = (/-s1_45,zero, s45, zz1,zz2,zero/)

       TCC54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    ! -------- TRIPLE COLLINEAR + SOFT 5

    spart = MV2/(one - x1) 
    if (spart.gt.sh) TCS5Lim%flag = .true. 

    sqrts = sqrt(spart)
    Deltan = one - x1*(one - eta41) 
    Deltad = one - x1*(one - eta42)

    if ( ( TCS5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCS5Lim%Flag = .true.
    else 
       TCS5Lim%flag = .true.
    endif

    if (TCS5Lim%flag .eqv. .false.) then 

       TCS5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)


!       z4 = one/x1
       !       z5 = x2/x1
       z4 = x1/(x1-one)
       z5 = x2/(x1-one)

       s14 = spart *x1*x3 
       s15 = spart *x2*x3*x4p
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       TCS5Lim%KinInv = (/-s14,-s15,s45,z4,z5,zero/)
       pVV = p1 + p2 - p4 
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1 -p4
       TCS5Lim%AmpMom(:,2) = p2 


       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ------   TRIPLE COLLINEAR KINEMATICS  + C54 + SOFT5 -------


       z45 = x1/(x1-one)
       z5=x2/x1

       s14 = spart *x1*x3 
       s15 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)

       TCC54S5Lim%KinInv = (/-s14,-s25, s45, z45,z5, zero/)

       TCC54S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C54 FLM 
    !-------------------------------------------------------------------------------------------------------------

    cos4 = one - two*x3 
    sin4 = sqrt(one-cos4**2)

    eta41 = x3
    eta51 = x3
    eta42 = (one - eta41) 
    eta52 = (one - eta51)
    eta54 = zero

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C54Lim%flag = .true.
    if (C54Lim%flag.eqv. .false.) spart = MV2/Deltas
    if (C54Lim%flag.eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( C54Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C54Lim%Flag = .true.
    else 
       C54Lim%flag = .true.
    endif

    if (C54Lim%flag .eqv. .false.) then 

       C54Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)  
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)


       pVV = p1 + p2 - p4 - p5

       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)

       z5 = x1/(x1 + x2) 

       C54Lim%KinInv = (/zero,zero,s45,zero,z5,zero/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C54Lim)


       C54Lim%AmpMom(:,1) = p1
       C54Lim%AmpMom(:,2) = p2
       C54Lim%AmpMom(:,3) = p4 + p5 


       C54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3 *x4/NF1(x3,one,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


    endif



  end subroutine kinematics_nnlo_4151b



  !-- we factor out an extra 1/(8*pisq)**2 from the phase space
  !-- which multiplied by gs**4 in the ME reconstructs asontwopi**2
  subroutine kinematics_nnlo_4151c(MV2,Y,yr, &
       HardProc,S5Lim,TCLim,TCS5Lim,  &
       C4Lim,C4S5Lim,TCC4Lim,TCC4S5Lim)
    real(dp), intent(in) ::  MV2, Y, yr(9)
    type(KinConfig) ::  HardProc,S5Lim,TCLim,TCS5Lim, C4Lim,C4S5Lim,TCC4Lim,TCC4S5Lim
    !---
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4) 
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5,  pVV(4)
    real(dp) :: xx3,x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4151, s14, s15, s45, z4, z5, ztot, weight4151,s5_14,zz1,zz2,kt1kt2sq
    real(dp) :: spart,sqrts, s12, s25, z45, yy, Deltas, s24   

    !    statement functions
    weight4151(eta41,eta42,eta51,eta52,eta54)  =  eta42*eta52*( one + eta41/(eta54 + eta42 + eta51) & 
         + eta51/(eta54 + eta41 + eta52) ) 

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S5Lim)          !3
    call initialize_config(TCLim)          !5
    call initialize_config(TCS5Lim)        !7
    call initialize_config(C4Lim)          !9
    call initialize_config(C4S5Lim)        !11
    call initialize_config(TCC4Lim)        !13
    call initialize_config(TCC4S5Lim)      !15

    !    MV2 and Y are generated

    MV =  sqrt(MV2)

    x1 =   yr(1) 
    x2 =   yr(2) 
    x3 = yr(3) 
    x4 = yr(5) 

    if (withmax) then
       xx3 = x3
    else
       xx3 = zero
    endif

    lambda = sin( half * pi * yr(6) )**2


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = one - x3*x4 
    sin4 = sqrt(one-cos4**2)

    cos5 = one - two*x3 
    sin5 = sqrt(one-cos5**2) 


    if (x4.eq.zero) then 
       if(one-two*lambda .ge. zero) then
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif
    elseif(2.0_dp*(one - half*x4)**2/NF1(x3,half*x4,lambda)-two*(one+half*x4*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(x3,half*x4,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    eta41 = half*x3*x4
    eta51 = x3
    eta42 = (one-eta41) 
    eta52 = (one - eta51)
    eta54 = x3*(one-half*x4)**2/NF1(x3,half*x4,lambda)


    !  incoming momenta 

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (HardProc%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) HardProc%Flag = .true.
    else
       HardProc%flag = .true.
    endif


    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2


    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 FLM   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------

    eta41 = half*x3*x4
    eta51 = x3
    eta42 = (one-eta41) 
    eta52 = (one - eta51)
    eta54 = x3*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    spart = MV2/(one - x1)
    if (spart.gt.sh) S5Lim%flag = .true.
    sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41)
    Deltad = one - x1*(one - eta42)

    if ( (S5Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) S5Lim%Flag = .true.
    else 
       S5Lim%flag = .true.
    endif

    if (S5Lim%flag .eqv. .false.) then 

       S5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p4

       S5Lim%LimMom(:,1) =  p1
       S5Lim%LimMom(:,2) =  p2
       S5Lim%LimMom(:,3) =  p4
       S5Lim%LimMom(:,4) =  p5 


       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    ! --------- TRIPLE COLLINEAR -------

    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero 


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( TCLim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCLim%Flag = .true.
    else 
       TCLim%flag = .true.
    endif

    if (TCLim%flag .eqv. .false.) then 

       TCLim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)
       p5 = half*sqrts*x2*(/one,zero,zero,one /)

       s14 = spart *x1*(half*x3*x4) 
       s15 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

       ztot = x1 + x2 - one 
       z4 = x1/ztot 
       z5 = x2/ztot

       TCLim%KinInv = (/-s14,-s15,s45,z4,z5,zero/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)

       TCLim%AmpMom(:,1) = p1 - p4 - p5 
       TCLim%AmpMom(:,2) = p2 

       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            * spart/MV2


       ! ---------   TC  + C42 ----------------

       s14 = spart * x1 *(half * x3 * x4)  
       s15 = spart * x2* x3
       s5_14 = s15 * (one-x1)
!       s45 = spart * x1*x2*x3    ! /NF1(xx3,zero,lambda)   used to divide, but it is 1

!       z4 = x1/( x1 - one) 
!       z45 = (x1 - one)/( x1 + x2 - one) 
       zz1 = x1/(x1-one)
       zz2 = (one-x1)/(one-x1-x2)
       kt1kt2sq = (-one+two*lambda)**2 * (one-xx3)


       TCC4Lim%KinInv = (/-s14,-s5_14,zero,zz1,zz2,kt1kt2sq/)

       TCC4Lim%wgt  = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    ! -------------- TRIPLE COLLINEAR KINEMATICS + SOFT 5  -------------     

    spart = MV2/(one - x1) 
    if (spart.gt.sh) TCS5Lim%flag = .true.

    sqrts = sqrt(spart)
    Deltan = one - x1*(one - eta41) 
    Deltad = one - x1*(one - eta42)

    if ( ( TCS5Lim%flag  .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCS5Lim%Flag = .true.
    else 
       TCS5Lim%flag  = .true.
    endif

    if (TCS5Lim%flag  .eqv. .false.) then 

       TCS5Lim%PartFrac = (/xi1,xi2/)

       s14 = spart *x1*(half*x3*x4) 
       s15 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)


       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)
       p5 = half*sqrts*x2*(/one,zero,zero,one /)

       z4 = x1/(x1-one) 
       z5 = x2/(x1-one)

       TCS5Lim%KinInv = (/-s14,-s15,s45,z4,z5,zero/)

       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1 - p4 
       TCS5Lim%AmpMom(:,2) = p2 

       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            * spart/MV2


       ! -------   TRIPLE COLLINEAR KINEMATICS  + C4 + SOFT5 -----------

       s14 = spart *x1*(half * x3 * x4)  
       s15 = spart *x2 * x3
       s45 = spart *x1*x2*x3                 ! /NF1(xx3,zero,lambda) ->1

       z4 = x1/(x1 - one) 
       z5 = x2

       TCC4S5Lim%KinInv = (/-s14,-s15,s45,z4,z5,zero/)


       TCC4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1 * x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2



       ! ----------  C42 + S5 

       eta41 = zero
       eta51 = x3
       eta42 = one 
       eta52 = (one - eta51)
       eta54 = x3/NF1(x3,zero,lambda)

       z4 =  x1/(x1 - one) 
       s12 = two*scr(p1,p2) 

       s14 = spart *x1*(half*x3*x4) 
       s15 = spart *x2*x3
       s25 = spart *x2*(one-x3)

       C4S5Lim%KinInv = (/s14,-s15,zero,z4,zero,zero/)
       C4S5Lim%EikInv = (/s12,s15,s25/)
! not needed for 4151 in qg channel, check for gq channel and in 4252 for qg

       ! check: do we need this?

       MV2aux = scr(pVV,pVV) 
       MV =  sqrt(MV2aux )


       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4S5Lim)


       C4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart & 
            * spart/MV2

    endif



    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C41 FLM 
    !-------------------------------------------------------------------------------------------------------------

    cos5 = one - two*x3 
    sin5 = sqrt(one-cos5**2) 

    if(one-two*lambda .ge. zero) then
       sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif


    eta41 = zero
    eta51 = x3
    eta42 = one 
    eta52 = (one - eta51)
    eta54 = x3/NF1(x3,zero,lambda)


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C4Lim%flag  = .true.
    if (C4Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( C4Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C4Lim%Flag = .true.
    else
       C4Lim%flag = .true.
    endif

    if ( C4Lim%flag .eqv. .false.) then 

       C4Lim%PartFrac = (/xi1,xi2/)
       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1*(/one, zero, zero, one  /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       C4Lim%LimMom(:,1) = (/zero,cos(phi4),sin(phi4),zero/)
       C4Lim%LimMom(:,2) = (/zero,zero,zero,one/)

       s14 = spart*x1*(half*x3*x4)
       z4 = x1/( x1-one ) 

!       C4Lim%KinInv = (/zero,-s14,zero,z4,zero,zero/)
       C4Lim%KinInv = (/zero,-s14,zero,zero,z4,zero/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4Lim)

       C4Lim%AmpMom(:,1) = p1-p4 
       C4Lim%AmpMom(:,2) = p2
       C4Lim%AmpMom(:,3) = p5 


       C4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif



  end subroutine kinematics_nnlo_4151c




  !-- we factor out an extra 1/(8*pisq)**2 from the phase space
  !-- which multiplied by gs**4 in the ME reconstructs asontwopi**2
  subroutine kinematics_nnlo_4151d(MV2,Y,yr,&
       HardProc,S5Lim,TCLim,TCS5Lim,C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim)
    real(dp), intent(in) ::  yr(9), MV2, Y
    type(KinConfig)      ::  HardProc,S5Lim,TCLim,TCS5Lim,&
         C54Lim,C54S5Lim,TCC54Lim,TCC54S5Lim
    !---
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV,  p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: xx3,x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4151, s14, s15, s45, z4, z5, ztot, weight4151
    real(dp) :: spart,sqrts, x4p, nperp(4), s12, s24, s25, z45, ncoll(1:4), yy, Deltas, s1_45,zz1,zz2

    !    statement functions
    weight4151(eta41,eta42,eta51,eta52,eta54)  =  eta42*eta52*( one + eta41/(eta54 + eta42 + eta51) & 
         + eta51/(eta54 + eta41 + eta52) ) 


    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S5Lim)          !3
    call initialize_config(TCLim)          !5
    call initialize_config(TCS5Lim)        !7
    call initialize_config(C54Lim)          !9
    call initialize_config(C54S5Lim)        !11
    call initialize_config(TCC54Lim)        !13
    call initialize_config(TCC54S5Lim)      !15

    !    MV2 and Y are generated

    MV =  sqrt(MV2)

    x1 =   yr(1) 
    x2 =   yr(2) 
    x3 = yr(3) 
    x4 = yr(5) 

    if (withmax) then
       xx3 = x3
    else
       xx3 = zero
    endif

    lambda = sin( half * pi * yr(6) )**2

    x4p = one - half*x4 

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = one - two*x3*x4p 
    sin4 = sqrt(one-cos4**2)

    cos5 = one - two*x3
    sin5 = sqrt(one-cos5**2) 


    if (x4p.eq.zero) then 
       if(one-two*lambda .ge. zero) then
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif
    elseif(2.0_dp*(one - x4p)**2/NF1(x3,x4p,lambda)-two*(one+x4p*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4p)/NF1(x3,x4p,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    eta41 = x3*x4p
    eta51 = x3
    eta42 = (one-eta41) 
    eta52 = (one - eta51)
    eta54 = x3*x4**2/four/NF1(x3,x4p,lambda)

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (HardProc%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) HardProc%Flag = .true.
    else
       HardProc%flag = .true.
    endif


    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,x4p,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 FLM   ! angles the same as in hard, do not regenerate; weight is the same 
    !-------------------------------------------------------------------------------------------------------------

    eta41 = x3*x4p
    eta51 = x3
    eta42 = (one - eta41) 
    eta52 = (one - eta51)
    eta54 = x3*x4**2/four/NF1(x3,x4p,lambda)

    spart = MV2/(one - x1)
    if (spart.gt.sh) S5Lim%flag = .true.
    sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41)
    Deltad = one - x1*(one - eta42)

    if ( (S5Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) S5Lim%Flag = .true.
    else 
       S5Lim%flag = .true.
    endif

    if (S5Lim%flag .eqv. .false.) then 

       S5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,S5Lim)
       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p4 

       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = p4 
       S5Lim%LimMom(:,4) = p5 

       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,x4p,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2 
    endif


    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C54 + S5 
    !-------------------------------------------------------------------------------------------------------------


    eta41 = x3
    eta51 = x3
    eta42 = (one - eta41) 
    eta52 = (one - eta51)
    eta54 = zero 

    phi4 = twopi*yr(4) 
    cos4 = one - two*x3 
    sin4 = sqrt(one-cos4**2)

    spart = MV2/(one - x1)
    if (spart.gt.sh) C54S5Lim%flag = .true.

    sqrts = sqrt(spart)
    Deltan = one - x1*(one - eta41)
    Deltad = one - x1*(one - eta42)

    if ( ( C54S5Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then
       yy = Y - half*log(Deltan/Deltad)
       xi1 = sqrt(spart/sh)*exp(yy)
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C54S5Lim%Flag = .true.
    else
       C54S5Lim%flag = .true.
    endif

    if (C54S5Lim%flag .eqv. .false.) then

       C54S5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/)
       p2 =half*sqrts*(/one,zero,zero,-one/)
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)


       z5 = x2/x1
       s45 = spart *x1  * x2 * x3 * x4**2/four/NF1(x3,one,lambda)

       C54S5Lim%KinInv = (/zero,zero,s45,zero,z5,zero/)

       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C54S5Lim)
       C54S5Lim%AmpMom(:,1) = p1
       C54S5Lim%AmpMom(:,2) = p2
       C54S5Lim%AmpMom(:,3) = p4



       C54S5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2
    endif



    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero 

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( TCLim%flag .eqv. .false.) .and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCLim%Flag = .true.
    else 
       TCLim%flag = .true.
    endif

    if (TCLim%flag .eqv. .false.) then 

       TCLim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)
       p5 = half*sqrts*x2*(/one,zero,zero, one /)

       s14 = spart *x1*x3*x4p 
       s15 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       ztot = x1 + x2 - one 
       z4 = x1/ztot 
       z5 = x2/ztot

       TCLim%KinInv = (/-s14,-s15,s45,z4,z5,zero/)


       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)


       TCLim%AmpMom(:,1) = p1- p4 - p5 
       TCLim%AmpMom(:,2) = p2


       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

       ! --------------  TC  + C54 ----------

       s14 = spart *x1*x3 
       s15 = spart *x2*x3
       s45 = spart *x1*x2*x3*x4**2/four/NF1(xx3,one,lambda)
       s1_45 = s14 * (x1+x2)/x1

!       z45 = (x1+x2)/(x1+x2-one)
!       z5 = x1/(x1+x2)
       zz1 = x2/(x1+x2)
       zz2 = (x1+x2)/(x1+x2-one)

       TCC54Lim%KinInv = (/-s1_45,zero,s45, zz1, zz2,zero/)

       TCC54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    spart = MV2/(one - x1) 
    if (spart.gt.sh) TCS5Lim%flag = .true.

    sqrts = sqrt(spart)
    Deltan = one - x1*(one - eta41) 
    Deltad = one - x1*(one - eta42)

    if ( ( TCS5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) TCS5Lim%Flag = .true.
    else 
       TCS5Lim%flag = .true.
    endif

    if (TCS5Lim%flag .eqv. .false.) then 

       TCS5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)

!       z4 = one/x1 
!       z5 = x2/x1
       z4 = x1/(x1-one)
       z5 = x2/(x1-one)

       s14 = spart *x1*x3*x4p 
       s15 = spart * x2*x3
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       TCS5Lim%KinInv = (/-s14,-s15,s45,z4,z5,zero/)

       pVV = p1 + p2 - p4 
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1- p4  
       TCS5Lim%AmpMom(:,2) = p2

       TCS5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2



       ! --------   TRIPLE COLLINEAR KINEMATICS  + C54 + SOFT5 ------------

       s14 = spart *x1*x3 
       s15 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)

       z45 = x1/(x1-one)
       z5 = x2/x1

       TCC54S5Lim%KinInv = (/-s14,-s15, s45, z45,z5, zero/)

       TCC54S5Lim%wgt = one/8.0_dp/pi * kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C54 FLM 
    !-------------------------------------------------------------------------------------------------------------


    phi4 = twopi*yr(4) 
    cos4 = one - two*x3 
    sin4 = sqrt(one-cos4**2)

    eta41 = x3
    eta51 = x3
    eta42 = (one - eta41) 
    eta52 = (one - eta51)
    eta54 = zero 

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C54Lim%flag = .true.
    if (C54Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C54Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( C54Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C54Lim%Flag = .true.
    else 
       C54Lim%flag = .true.
    endif

    if (C54Lim%flag .eqv. .false.) then 

       C54Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)


       s45 = spart *x1 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)
       z5 = x1/(x1 + x2) 

       C54Lim%KinInv = (/zero,zero,s45,zero,z5,zero/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C54Lim)

       C54Lim%AmpMom(:,1) = p1
       C54Lim%AmpMom(:,2) = p2
       C54Lim%AmpMom(:,3) = p4 + p5 

       C54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart & 
            * spart/MV2

    endif


  end subroutine kinematics_nnlo_4151d

  !-- we factor out an extra 1/(8*pisq)**2 from the phase space
  !-- which multiplied by gs**4 in the ME reconstructs asontwopi**2
  subroutine kinematics_nnlo_4152(MV2,Y,yr,&
       HardProc,S5Lim,C4Lim,C4S5Lim,&
       C5Lim,C5S5Lim,C4C5Lim,C4C5S5Lim) 
    real(dp), intent(in) ::  yr(9), MV2, Y
    type(KinConfig)      :: HardProc,S5Lim,C4Lim,C4S5Lim, C5Lim, C5S5Lim, C4C5Lim, C4C5S5Lim
    !---
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: xx3,x1, x2,x3, x4, lambda,  cos4, sin4, phi4, phi5, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4151, s14, s15, s45, z4, z5, ztot, weight4152
    real(dp) :: spart,sqrts, s25, s12, s24, yy, Deltas

    !    statement functions

    weight4152(eta41,eta42,eta51,eta52,eta54)  = eta42*eta51*eta54/(eta54+eta41+eta52)


    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)         !1
    call initialize_config(S5Lim)            !2
    call initialize_config(C5Lim)            !3 
    call initialize_config(C5S5Lim)          !4 
    call initialize_config(C4Lim)            !5
    call initialize_config(C4S5Lim)          !6
    call initialize_config(C4C5Lim)          !7
    call initialize_config(C4C5S5Lim)        !8

    !    MV2 and Y are generated

    MV =  sqrt(MV2)

    x1 =   yr(1) 
    x2 =   yr(2) 
    x3 =   yr(3) 
    x4 =   yr(5) 


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = one -  two*x3 
    sin4 = sqrt(one-cos4**2)

    phi5 = twopi*yr(6)
    cos5 = -one + two*x4
    sin5 = sqrt(one-cos5**2) 

    eta42 = (one - x3)
    eta52 = x4
    eta41 = x3 
    eta51 = (one - x4) 
    eta54 =   half*( one - sin4*sin5*cos(phi4-phi5) -cos4*cos5)

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (HardProc%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) HardProc%Flag = .true.
    else
       HardProc%flag = .true.
    endif


    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 



       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3: S5 FLM   ! angles the same as in hard // do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! --------------- SOFT 5 ----------


    spart = MV2/(one - x1)

    if (spart.gt.sh) S5Lim%flag = .true.
    sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41)
    Deltad = one - x1*(one - eta42)

    if ( (S5Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) S5Lim%Flag = .true.
    else 
       S5Lim%flag = .true.
    endif

    if (S5Lim%flag .eqv. .false.) then 

       S5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)


       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p4 

       ! --      note momentum swap below; this is because of the definition of the Eikonal_s5

       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = p4 
       S5Lim%LimMom(:,4) = p5 

       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

       ! ---------  C52 + S5 --------------

       phi4 = twopi*yr(4) 
       cos4 = one - two*x3 
       sin4 = sqrt(one-cos4**2)

       eta41 = x3
       eta51 = one
       eta42 = (one-x3) 
       eta52 = zero
       eta54 = (one-x3) 


       s25 = spart*x2*x4
       z5 = x2

       C5S5Lim%KinInv = (/zero,s25,zero,zero,z5,zero/)

       C5S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3:   C41  
    !-------------------------------------------------------------------------------------------------------------

    eta42 = one
    eta52 = x4 
    eta41 = zero  
    eta51 = (one - x4) 
    eta54 =   eta51


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C4Lim%flag  = .true.
    if (C4Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4Lim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( C4Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C4Lim%Flag = .true.
    else 
       C4Lim%flag = .true.
    endif

    if (C4Lim%flag .eqv. .false.) then 

       C4Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one,  zero, zero, one  /)
       p5 = half*sqrts*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)

       s14 = spart *x1*x3 
       z4 = x1/(x1-one)

       C4Lim%KinInv = (/s14,zero,zero,z4,zero,zero/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4Lim)

       C4Lim%AmpMom(:,1) = p1 - p4 
       C4Lim%AmpMom(:,2) = p2 
       C4Lim%AmpMom(:,3) = p5

       C4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2 

    endif


    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C41  + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------


    spart = MV2/(one - x1) 
    if (spart.gt.sh) C4S5Lim%flag = .true.

    sqrts = sqrt(spart)
    Deltan = one - x1*(one - eta41) 
    Deltad = one - x1*(one - eta42)

    if ( ( C4S5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C4S5Lim%Flag = .true.
    else 
       C4S5Lim%flag = .true.
    endif

    if (C4S5Lim%flag .eqv. .false.) then 

       C4S5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)

       s14 = spart *x1*x3 
       s15 = spart*x2*(one-x4) 
       s25 = spart*x2*x4
       s12 = two*scr(p1,p2) 

       z4 = x1/(x1 - one)  

       C4S5Lim%KinInv = (/s14,s15,zero,z4,zero,zero/)
       C4S5Lim%EikInv = (/s12,s15,s25/)

       pVV = p1 + p2 - p4 

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4S5Lim)

       C4S5Lim%AmpMom(:,1) = p1 - p4 
       C4S5Lim%AmpMom(:,2) = p2 


       C4S5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ---------   C41   + C52 + SOFT5 ------------


       eta41 = zero
       eta51 = one
       eta42 = one
       eta52 = zero
       eta54 = one           

       z4 = x1/(x1-one)
       z5 = x2

       s14 = spart *x1*x3 
       s25 = spart *x2*x4

       C4C5S5Lim%KinInv = (/-s14,-s25,zero,z4,z5,zero/)

       C4C5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: C52 FLM 
    !-------------------------------------------------------------------------------------------------------------


    phi4 = twopi*yr(4) 
    cos4 = one - two*x3 
    sin4 = sqrt(one-cos4**2)

    eta41 = x3
    eta51 = one
    eta42 = (one-x3) 
    eta52 = zero
    eta54 = (one-x3) 


    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C5Lim%flag = .true.
    if (C5Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C5Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( C5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C5Lim%Flag = .true.
    else 
       C5Lim%flag = .true.
    endif

    if (C5Lim%flag .eqv. .false.) then 

       C5Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1*(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)
       p5 = half*sqrts*x2*(/one, zero, zero, -one /)

       s25 = spart*x2*x4
       z5 = one/(one - x2) 

       C5Lim%LimMom(:,1) = (/zero,cos(phi5),sin(phi5),zero/)
       C5Lim%LimMom(:,2) = (/zero,zero,zero,-one/)

       C5Lim%KinInv = (/zero,-s25,zero,zero,z5,zero/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C5Lim)

       C5Lim%AmpMom(:,1) = p1 
       C5Lim%AmpMom(:,2) = p2 - p5 
       C5Lim%AmpMom(:,3) = p4 

       C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif




    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C41  + C52
    !-------------------------------------------------------------------------------------------------------------

    eta41 = zero
    eta51 = one
    eta42 = one
    eta52 = zero
    eta54 = one  

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C4C5Lim%flag = .true.
    if (C4C5Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4C5Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( ( C4C5Lim%flag .eqv. .false. ).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C4C5Lim%Flag = .true.
    else 
       C4C5Lim%flag = .true.
    endif

    if (C4C5Lim%flag .eqv. .false.) then 

       C4C5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)
       p5 = half*sqrts*x2*(/one,zero,zero,-one /)

       s14 = spart *x1*x3 
       s25 = spart *x2*x4
       z4 = x1/(x1-one)
       z5 = one/(one - x2) 

       C4C5Lim%KinInv = (/-s14,-s25,zero,z4,z5,zero/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4C5Lim)

       C4C5Lim%AmpMom(:,1) = p1 - p4 
       C4C5Lim%AmpMom(:,2) = p2 - p5

       C4C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


    endif


  end subroutine kinematics_nnlo_4152

end module mod_kinematicsNNLO_EunordAord
    

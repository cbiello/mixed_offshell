module mod_kinematicsNNLO_EordAord
  use mod_types; use mod_consts_dp; use common_def
  use mod_process
  use mod_auxfunctions; use mod_procfunctions
  use mod_kinematics_gen; use mod_kinematics_lept
  implicit none
  private

  public :: kinematics_nnlo_4151a,kinematics_nnlo_4151b,kinematics_nnlo_4151c,kinematics_nnlo_4151d,kinematics_nnlo_4251
  public :: kinematics_nnlo_4252a,kinematics_nnlo_4252b,kinematics_nnlo_4252c,kinematics_nnlo_4252d,kinematics_nnlo_4152

contains 
  
  !-- we factor out an extra 1/(8*pisq)**2 from the phase space
  !-- which multiplied by gs**4 in the ME reconstructs asontwopi**2
  subroutine kinematics_nnlo_4151a(MV2,Y,yr,&
       HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
       C5Lim,C5DSLim,C5S5Lim,C5DSS5Lim,TCC5Lim,TCC5DSLim,TCC5S5Lim,TCC5DSS5Lim)       
    real(dp), intent(in) ::  MV2,Y, yr(9)
    type(KinConfig)      ::  HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C5Lim,C5DSLim,C5DSS5Lim,C5S5Lim,TCC5Lim,TCC5DSLim,TCC5DSS5Lim,TCC5S5Lim
    !---
    real(dp) :: p1(4), p2(4), p4(4), p5(4)
    real(dp) :: MV, kallenF, sin5, cos5, pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4151, s14, s15, s45, z4, z5, ztot, weight4151
    real(dp) :: spart,sqrts, yy, Deltas
    real(dp) :: ncoll(4),nperp(4)
    real(dp) :: mysin4lim,mysin5lim,kt1kt2sq,z45,s4_15,zz1,zz2
    real(dp) :: xx3

    !    statement functions
    weight4151(eta41,eta42,eta51,eta52,eta54)  = eta42*eta52*( one + eta41/(eta54 + eta42 + eta51) & 
         + eta51/(eta54 + eta41 + eta52) )

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(DSLim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(DSS5Lim)        !4
    call initialize_config(TCLim)          !5
    call initialize_config(TCDSLim)        !6
    call initialize_config(TCS5Lim)        !7
    call initialize_config(TCDSS5Lim)      !8
    call initialize_config(C5Lim)          !9
    call initialize_config(C5DSLim)        !10
    call initialize_config(C5S5Lim)        !11
    call initialize_config(C5DSS5Lim)      !12
    call initialize_config(TCC5Lim)        !13
    call initialize_config(TCC5DSLim)      !14
    call initialize_config(TCC5S5Lim)      !15
    call initialize_config(TCC5DSS5Lim)    !16

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

    eta41 = x3
    eta51 = half*x3*x4
    eta42 = (one-x3) 
    eta52 = (one - half*x3*x4)
    eta54 = x3*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif
    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: DS FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! ------- DOUBLE SOFT and DOUBLE SOFT + S5 ------------

    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    DSLim%PartFrac = (/xi1,xi2/)

    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
    p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
         sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

    pVV = p1 + p2

    call get_lept_mom(MV,pVV,yr(7:9),kallenF,DSLim)

    DSLim%AmpMom(:,1) = p1
    DSLim%AmpMom(:,2) = p2

    DSLim%LimMom(:,1) = p1
    DSLim%LimMom(:,2) = p2
    DSLim%LimMom(:,3) = p4 
    DSLim%LimMom(:,4) = p5

    DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2  
    ! the weight for DS + S5 (identical to above, but avoids confusion)
    DSS5Lim%wgt= DSLim%wgt


    ! ------   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT  -------

    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero

    p1=  half*sqrts*(/one,zero,zero,one/)
    p2 = half*sqrts*(/one,zero,zero,-one/)
    p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)
    p5 = half*sqrts*x1*x2*(/one, zero,zero,one /)

    ztot = - one 
    z4 = x1/ztot 
    z5 = x1*x2/ztot

    s14 = spart *x1*x3 
    s15 = spart *x1*x2*(half*x3*x4)
    s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)


    TCDSLim%KinInv = (/s14,s15,s45,z4,z5/)

    pVV = p1 + p2

    call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCDSLim)

    TCDSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         *one/two/spart &
         *spart/MV2

    ! ------   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT  + S5 -------

    ztot = - one 
    z4 = x1/ztot
    z5 = x1*x2/ztot
    s14 = spart *x1*x3 
    s15 = spart *x1*x2*(half*x3*x4)
    s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

    TCDSS5Lim%KinInv = (/s14,s15,s45,z4,z5/)

    TCDSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         *one/two/spart &
         *spart/MV2

    ! --------  C5 DS FLM  ----------

    eta41 = x3
    eta51 = zero
    eta42 = (one-x3)
    eta52 = one
    eta54 = x3/NF1(x3,zero,lambda)


    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)
    p5 = half*sqrts*x1*x2*(/one, zero, zero, one /)

    ztot = - one
    z4 = x1/ztot
    z5 = x1*x2/ztot

    s14 = spart *x1*x3
    s15 = spart *x1*x2*(half*x3*x4)
    s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

    C5DSLim%KinInv = (/s14,s15,s45,z4,z5/)

    pVV = p1 + p2


    call get_lept_mom(MV,pVV,yr(7:9),kallenF,C5DSLim)

    C5DSLim%LimMom(:,1) = p1 
    C5DSLim%LimMom(:,2) = p2
    C5DSLim%LimMom(:,3) = p4 

    C5DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         *one/two/spart &
         * spart/MV2


    ! ----------- C5+DS+S5 --------------

    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

    s15 = spart*x1*x2*(half*x3*x4)
    z5 = x1*x2

    C5DSS5Lim%KinInv = (/zero,s15,zero,zero,z5/)

    pVV = p1 + p2

    call get_lept_mom(MV,pVV,yr(7:9),kallenF,C5DSS5Lim)


    C5DSS5Lim%LimMom(:,1) = p1
    C5DSS5Lim%LimMom(:,2) = p2
    C5DSS5Lim%LimMom(:,3) = p4 

    C5DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2


    ! -------   TRIPLE COLLINEAR KINEMATICS  + C5 + Double soft  --------------
    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero


    s14 = spart *x1*x3 
    s15 = spart *x1*x2*(half*x3*x4)
    s45 = spart *x1**2*x2*x3/NF1(xx3,zero,lambda)
    ztot =-one 
    z4 = x1/ztot 
    z5 = x1*x2/ztot

    TCC5DSLim%KinInv = (/s14,s15,s45,z4,z5/)

    TCC5DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(xx3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2

    !---------   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT5 +DS ------------


    z4 = x1 
    z5 = x1*x2
    s14 = spart *x1*x3 
    s15 = spart *x1*x2*(half*x3*x4)
    s45 = spart *x1**2*x2*x3/NF1(xx3,zero,lambda)

    TCC5DSS5Lim%KinInv = (/s14,s15,s45,z4,z5/)

    TCC5DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(xx3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 + C5S5   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------

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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
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
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2   


       ! kin inv and weight for C5 +  S5
       eta41 = x3
       eta51 = zero
       eta42 = (one-x3)
       eta52 = one
       eta54 = x3/NF1(x3,zero,lambda)

       s15 = spart*x1*x2*(half*x3*x4)
       z5 = x1*x2
       C5S5Lim%KinInv = (/zero,s15,zero,zero,z5/)

       C5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1**3*x2 & 
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


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one,zero,zero,one /)

       s14 = spart *x1*x3 
       s15 = spart *x1*x2*(half*x3*x4)
       s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

       ztot = x1*(one +x2) - one 
       z4 = x1/ztot 
       z5 = x1*x2/ztot

       TCLim%KinInv = (/s14,s15,s45,z4,z5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)

       TCLim%AmpMom(:,1) = p1 - p4 - p5
       TCLim%AmpMom(:,2) = p2

       if(2.0_dp*(one - half*x4)**2/NF1(xx3,half*x4,lambda)-two*(one+half*x4).gt.0) then
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(xx3,half*x4,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-half*x4)/NF1(xx3,half*x4,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif

       mysin4lim = two!*sqrt(x3)(1+O(x3))
       mysin5lim = two * sqrt(x4/two)!*sqrt(x3)(1+O(x3))

       TCLim%LimMom(:,1) = mysin4lim*(/ one, cos(phi4), sin(phi4), cos4 /)
       TCLim%LimMom(:,2) = mysin5lim*(/ one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)
       TCLim%LimMom(:,3) = (/ one, zero,zero,one/)

       
       ! Esum = E3+E4+E5

       ! z3 = E3/Esum
       ! z4 = E4/Esum
       ! z5 = E5/Esum
       
       ! !-- kti = ki - zi kp_sum, k3=0

       ! kt3(:) = (-z3*z4*kp4(:) - z3*z5*kp5(:))*Esum
       ! kt4(:) = ((one-z4)*z4*kp4(:) - z4*z5*kp5(:))*Esum
       ! kt5(:) = ((one-z5)*z5*kp5(:) - z5*z4*kp4(:))*Esum


       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2


       ! --------- TRIPLE COLLINEAR + C5  -------

       s14 = spart *x1*x3 
       s15 = spart *x1*x2*(half*x3*x4)
       s4_15 = s14 * (one-x1*x2)

       zz1 = x1*x2/(x1*x2-one)
       zz2 = (one-x1*x2)/(one-x1-x1*x2)

       kt1kt2sq = (-one+two*lambda)**2 * (one-xx3)

       TCC5Lim%KinInv = (/s15,s4_15,zz1,zz2,kt1kt2sq/)

       TCC5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1**3*x2 & 
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
       p5 = half*sqrts*x1*x2*(/one,zero,zero,one /)


       z4 = x1 
       z5 = x1*x2

       s14 = spart *x1*x3 
       s15 = spart *x1*x2*(half*x3*x4)
       s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

       TCS5Lim%KinInv = (/s14,s15,s45,z4,z5/)

       pVV = p1 + p2 - p4 

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1 - p4 
       TCS5Lim%AmpMom(:,2) = p2

       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            *spart/MV2


       ! -------   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT5 ---------
       z4 = x1 
       z5 = x1*x2

       s14 = spart *x1*x3 
       s15 = spart *x1*x2*(half*x3*x4)
       s45 = spart *x1**2*x2*x3/NF1(xx3,zero,lambda)

       TCC5S5Lim%KinInv = (/s14,s15,s45,z4,z5/)

       TCC5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
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

    eta41 = x3
    eta51 = zero
    eta42 = (one-x3) 
    eta52 = one
    eta54 = x3/NF1(x3,zero,lambda)


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C5Lim%flag = .true.
    if (C5Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C5Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, zero, zero, one /)

       s15 = spart*x1*x2*(half*x3*x4)
       z5 = one/(one - x1 * x2) 

       C5Lim%KinInv = (/zero,s15,zero,zero,z5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C5Lim)

       C5Lim%AmpMom(:,1) = p1 - p5
       C5Lim%AmpMom(:,2) = p2
       C5Lim%AmpMom(:,3) = p4 

       !-- spin-correlation vectors

       cosphi45 = -one + two * lambda
       sinphi45 = sqrt(one-cosphi45**2)

       nperp = (/zero,(cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            (sin(phi4)*cosphi45+cos(phi4)*sinphi45),zero/)

       ncoll = (/zero,zero,zero,one/)

       C5Lim%LimMom(:,1) = nperp
       C5Lim%LimMom(:,2) = ncoll

       C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif





  end subroutine kinematics_nnlo_4151a


  subroutine kinematics_nnlo_4151b(MV2,Y,yr,&
       HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
       C54Lim,C54DSLim,C54S5Lim,C54DSS5Lim,TCC54Lim,TCC54DSLim,TCC54S5Lim,TCC54DSS5Lim) 
    real(dp), intent(in) ::  yr(9), MV2, Y
    type(KinConfig)      ::  HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C54Lim,C54DSLim,C54DSS5Lim,C54S5Lim,TCC54Lim,TCC54DSLim,TCC54DSS5Lim,TCC54S5Lim
    !---
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, Deltas, w4151, s14, s15, s45, z4, z5, ztot, weight4151
    real(dp) :: spart,sqrts, x4p, nperp(4), s12, s24, z45, ncoll(1:4), yy 
    real(dp) :: kt1kt2sq,s1_45,z1,zz1,zz2
    real(dp) :: xx3

    !    statement functions
    weight4151(eta41,eta42,eta51,eta52,eta54)  = eta42*eta52*( one + eta41/(eta54 + eta42 + eta51) & 
         + eta51/(eta54 + eta41 + eta52) )

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(DSLim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(DSS5Lim)        !4
    call initialize_config(TCLim)          !5
    call initialize_config(TCDSLim)        !6
    call initialize_config(TCS5Lim)        !7
    call initialize_config(TCDSS5Lim)      !8
    call initialize_config(C54Lim)          !9
    call initialize_config(C54DSLim)        !10
    call initialize_config(C54S5Lim)        !11
    call initialize_config(C54DSS5Lim)      !12
    call initialize_config(TCC54Lim)        !13
    call initialize_config(TCC54DSLim)      !14
    call initialize_config(TCC54S5Lim)      !15
    call initialize_config(TCC54DSS5Lim)    !16

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


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,x4p,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart  & 
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: DS FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate ; weight is the same 
    !-------------------------------------------------------------------------------------------------------------

    !   DS + DS S5 FLM   ! angles etc. are the same ; weight is the same 
    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    DSLim%PartFrac = (/xi1,xi2/)

    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
    p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
         sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

    pVV = p1 + p2
    call get_lept_mom(MV,pVV,yr(7:9),kallenF,DSLim)

    DSLim%AmpMom(:,1) = p1
    DSLim%AmpMom(:,2) = p2

    DSLim%LimMom(:,1) = p1
    DSLim%LimMom(:,2) = p2
    DSLim%LimMom(:,3) = p4 
    DSLim%LimMom(:,4) = p5

    DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    ! -----   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT   --------

    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero 

    s14 = spart *x1*x3 
    s15 = spart *x1*x2*x3*x4p
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)
    ztot = - one 
    z4 = x1/ztot 
    z5 = x1*x2/ztot

    TCDSLim%KinInv = (/s14,s15,s45,z4,z5/)

    TCDSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    !-------------------------------------------------------------------------------------------------------------
    !   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT  + S5
    !-------------------------------------------------------------------------------------------------------------

    s14 = spart *x1*x3 
    s15 = spart *x1*x2*x3*x4p
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)
    ztot = - one 
    z4 = x1/ztot
    z5 = x1*x2/ztot

    TCDSS5Lim%KinInv = (/s14,s15,s45,z4,z5/)

    TCDSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    ! --------  C54 DS FLM  ----------
    eta41 = x3
    eta51 = x3
    eta42 = (one - eta41) 
    eta52 = (one - eta51)
    eta54 = zero

    z5 = x2/(one+x2) 

    s12 = two*scr(p1,p2) 
    s14 = spart *x1*x3 
    s24 = spart * x1* (  one - x3)
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)

    kt1kt2sq = lambda 
    C54DSLim%KinInv = (/s12,s14,s24,s45,z5,kt1kt2sq/)

    C54DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,one,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


    !------------  C54 +DS+S5 ---------------

    s12 = two*scr(p1,p2) 
    s14 = spart *x1*x3 
    s24 = spart * x1* (  one - x3)
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)
    z5 = x2 

    C54DSS5Lim%KinInv = (/s12,s45,s14,s24,z5/)
    C54DSS5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,one,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    ! ---------   TC  + C54 + Double soft  ----------

    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero 

    !--NEW
    s14 = spart *x1*x3 
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
    s1_45 = s14 * (one+x2)
    z5 = one/(one+x2)
    z1 = -one/(x1+x1*x2)
    kt1kt2sq = lambda * (one-xx3)

    TCC54DSLim%KinInv = (/s45,s1_45,z5,z1,kt1kt2sq/)

    TCC54DSLim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,one,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    !-------------------------------------------------------------------------------------------------------------
    !   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT5 +DS
    !-------------------------------------------------------------------------------------------------------------

    s14 = spart *x1*x3 
    s15 = spart *x1*x2*x3
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
    z45 = one - x1 
    z5 = x2
    TCC54DSS5Lim%KinInv = (/s14,s15, s45, z45, z5/)
    TCC54DSS5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,one,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2  


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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
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
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ------  C54 + S5 ---------

       eta41 = x3
       eta51 = x3
       eta42 = (one - eta41) 
       eta52 = (one - eta51)
       eta54 = zero

       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)
       z5 = x2

       C54S5Lim%KinInv = (/zero,s45,zero,zero,z5/)

       C54S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1**3*x2 & 
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
    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = x2*p4

       s14 = spart *x1*x3 
       s15 = spart *x1*x2*x3*x4p
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       ztot = x1*(one +x2) - one 
       z4 = x1/ztot 
       z5 = x1*x2/ztot

       TCLim%KinInv= (/s14,s15,s45,z4,z5/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)

       TCLim%AmpMom(:,1) = p1 - p4 - p5
       TCLim%AmpMom(:,2) = p2

       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! -------   TC  + C54 -------------

       s14 = spart *x1*x3 
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
       s1_45 = s14 * (one+x2)
       zz1 = x2/(one+x2)
       zz2 = x1*(one+x2)/(x1*(one+x2)-one)
       kt1kt2sq = lambda * (one-xx3)

       TCC54Lim%KinInv = (/s45,s1_45,zz1,zz2,kt1kt2sq/)

       TCC54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1**3*x2 & 
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

       z4 = x1 
       z5 = x1*x2
       s14 = spart *x1*x3 
       s15 = spart *x1*x2*x3*x4p
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       TCS5Lim%KinInv = (/s14,s15,s45,z4,z5/)
       pVV = p1 + p2 - p4 
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1 - p4 
       TCS5Lim%AmpMom(:,2) = p2


       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ------   TRIPLE COLLINEAR KINEMATICS  + C54 + SOFT5 -------

       z45 = one - x1

       s14 = spart *x1*x3 
       s15 = spart *x1*x2*x3
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)

       TCC54S5Lim%KinInv = (/s14,s15, s45, z45, zero/)

       TCC54S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1**3*x2 & 
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

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C54Lim%flag = .true.
    if (C54Lim%flag.eqv. .false.) spart = MV2/Deltas
    if (C54Lim%flag.eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = x2*p4 

       ncoll = (/ zero,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)

       nperp  = (/ zero,  -cos4*sqrt(lambda)*cos(phi4) - sin(phi4)*sqrt(one-lambda) , & 
            -cos4*sqrt(lambda)*sin(phi4) + cos(phi4)*sqrt(one-lambda), & 
            sin4*sqrt(lambda) /)

       pVV = p1 + p2 - p4 - p5

       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)

       z5 = x2/(one + x2) 

       C54Lim%KinInv = (/zero,zero,s45,zero,z5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C54Lim)


       C54Lim%AmpMom(:,1) = p1
       C54Lim%AmpMom(:,2) = p2
       C54Lim%AmpMom(:,3) = p4 + p5 

       C54Lim%LimMom(:,1) = nperp
       C54Lim%LimMom(:,2) = ncoll

       C54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif

  end subroutine kinematics_nnlo_4151b


  subroutine kinematics_nnlo_4151c(MV2,Y,yr, &
       HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
       C4Lim,C4DSLim,C4S5Lim,C4DSS5Lim,TCC4Lim,TCC4DSLim,TCC4S5Lim,TCC4DSS5Lim)
    real(dp), intent(in) ::  MV2, Y, yr(9)
    type(KinConfig) ::  HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C4Lim,C4DSLim,C4DSS5Lim,C4S5Lim,TCC4Lim,TCC4DSLim,TCC4DSS5Lim,TCC4S5Lim
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4) 
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5,  pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4151, s14, s15, s45, z4, z5, ztot, weight4151
    real(dp) :: spart,sqrts, s12, s25, z45, yy, Deltas 
    real(dp) :: nperp(4), ncoll(4), mycphi, z1, z2, s14_5
    real(dp) :: s5_14,zz1,zz2,kt1kt2sq
    real(dp) :: xx3

    !    statement functions
    weight4151(eta41,eta42,eta51,eta52,eta54)  = eta42*eta52*( one + eta41/(eta54 + eta42 + eta51) & 
         + eta51/(eta54 + eta41 + eta52) )

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(DSLim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(DSS5Lim)        !4
    call initialize_config(TCLim)          !5
    call initialize_config(TCDSLim)        !6
    call initialize_config(TCS5Lim)        !7
    call initialize_config(TCDSS5Lim)      !8
    call initialize_config(C4Lim)          !9
    call initialize_config(C4DSLim)        !10
    call initialize_config(C4S5Lim)        !11
    call initialize_config(C4DSS5Lim)      !12
    call initialize_config(TCC4Lim)        !13
    call initialize_config(TCC4DSLim)      !14
    call initialize_config(TCC4S5Lim)      !15
    call initialize_config(TCC4DSS5Lim)    !16

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

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2


    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: DS FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! ------- DOUBLE SOFT and DOUBLE SOFT + S5 ------------

    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    DSLIM%PartFrac = (/xi1,xi2/)

    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
    p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
         sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

    pVV = p1 + p2


    call get_lept_mom(MV,pVV,yr(7:9),kallenF,DSLim)

    DSLim%AmpMom(:,1) = p1
    DSLim%AmpMom(:,2) = p2

    DSLim%LimMom(:,1) = p1
    DSLim%LimMom(:,2) = p2
    DSLim%LimMom(:,3) = p4 
    DSLim%LimMom(:,4) = p5

    w4151 = eta42*eta52*( one + eta41/(eta54 + eta42 + eta51) & 
         + eta51/(eta54 + eta41 + eta52) )


    DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    ! identical to DS but avoids confusion 
    DSS5Lim%wgt = DSLim%wgt


    ! -----   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT   ----------

    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero

    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)
    p5 = half*sqrts*x1*x2*(/one,zero,zero,one /)

    s14 = spart *x1*(half*x3*x4) 
    s15 = spart *x1*x2*x3
    s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

    ztot = - one 
    z4 = x1/ztot 
    z5 = x1*x2/ztot

    TCDSLim%KinInv = (/s14,s15,s45,z4,z5/)

    pVV = p1 + p2

    call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCDSLim)

    TCDSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         *one/two/spart &
         * spart/MV2

    ! -----   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT  + S5 --------

    ztot = - one 
    z4 = x1/ztot
    z5 = x1*x2/ztot
    s14 = spart *x1*(half*x3*x4) 
    s15 = spart *x1*x2*x3
    s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

    TCDSS5Lim%KinInv = (/s14,s15,s45,z4,z5/)

    TCDSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2



    ! -------  C41 DS FLM  ----------------


    eta41 = zero
    eta51 = x3
    eta42 = one 
    eta52 = (one - eta51)
    eta54 = x3/NF1(x3,zero,lambda)


    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one, zero, zero, one /)

    s12 = two*scr(p1,p2) 
    s14 = spart *x1*(half*x3*x4) 
    s15 = spart *x1*x2*x3
    s25 = spart *x1*x2*(one-x3)

    z4 = x1

    C4DSLim%KinInv = (/s12,s15,s25,s14,z4/)

    C4DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         *one/two/spart & 
         * spart/MV2


    ! -------  C41+DS+S5 ----------

    s12 = two*scr(p1,p2) 
    s14 = spart *x1*(half*x3*x4) 
    s15 = spart *x1*x2*x3
    s25 = spart *x1*x2*(one-x3)

    z4 = x1
    C4DSS5Lim%KinInv = (/s12,s15,s25,s14,z4/)

    C4DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         *one/two/spart & 
         * spart/MV2


    ! -------    TRIPLE COLLINEAR KINEMATICS  + C41 + Double soft  --------------


    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero 

    s14 = spart *x1*(half * x3 * x4)  
    s15 = spart *x1*x2* x3
    s45 = spart *x1**2*x2*x3/NF1(xx3,zero,lambda)

    z4 = x1
    s25 = spart *x1*x2

    TCC4DSLim%KinInv = (/s12,s14,s15,s25,z4/)

    TCC4DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(xx3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2



    ! -----------   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT5 +DS --------------

    s14 = spart *x1*(half * x3 * x4)  
    s15 = spart *x1*x2* x3
    s45 = spart *x1**2*x2*x3/NF1(xx3,zero,lambda)

    z4 = one - x1 
    z5 = x1*x2

    TCC4DSS5Lim%KinInv= (/s14,s15,s45,z4,z5/)

    TCC4DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(xx3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
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
            * x1**3*x2 & 
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


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one,zero,zero,one /)

       s14 = spart *x1*(half*x3*x4) 
       s15 = spart *x1*x2*x3
       s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

       ztot = x1*(one +x2) - one 
       z4 = x1/ztot 
       z5 = x1*x2/ztot

       TCLim%KinInv = (/s14,s15, s45, z4, z5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)

       TCLim%AmpMom(:,1) = p1 - p4 - p5
       TCLim%AmpMom(:,2) = p2

       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            * spart/MV2


       ! ---------   TC  + C41 ----------------

       s14 = spart *x1*(half * x3 * x4)  
       s15 = spart *x1*x2* x3
       s5_14 = s15 * (one-x1)

       zz1 = x1/(x1-one)
       zz2 = (one-x1)/(one-x1-x1*x2)

       kt1kt2sq = (-one+two*lambda)**2 * (one-xx3)

       TCC4Lim%KinInv = (/s14,s5_14,zz1,zz2,kt1kt2sq/)


       TCC4Lim%wgt  = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1**3*x2 & 
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
       s15 = spart *x1*x2*x3
       s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)


       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)
       p5 = half*sqrts*x1*x2*(/one,zero,zero,one /)

       z4 = x1 
       z5 = x1*x2

       TCS5Lim%KinInv = (/s14,s15,s45,z4,z5/)

       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1 - p4 
       TCS5Lim%AmpMom(:,2) = p2

       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            * spart/MV2


       ! -------   TRIPLE COLLINEAR KINEMATICS  + C4 + SOFT5 -----------

       s14 = spart *x1*(half * x3 * x4)  
       s15 = spart *x1*x2* x3
       s45 = spart *x1**2*x2*x3/NF1(xx3,zero,lambda)

       z4 = one - x1 
       z5 = x1*x2

       TCC4S5Lim%KinInv = (/s14,s15,s45,z4,z5/)


       TCC4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ----------  C41 + S5 


       eta41 = zero
       eta51 = x3
       eta42 = one 
       eta52 = (one - eta51)
       eta54 = x3/NF1(x3,zero,lambda)


       z4 = one - x1
       s12 = two*scr(p1,p2) 
       s14 = spart *x1*(half*x3*x4) 
       s15 = spart *x1*x2*x3
       s25 = spart *x1*x2*(one-x3)

       C4S5Lim%KinInv = (/s12,s15,s25,s14,z4/)

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4S5Lim)

       C4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1**3*x2 & 
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


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C4Lim%flag  = .true.
    if (C4Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       s14 = spart*x1*(half*x3*x4)
       z4 = one/(one - x1 ) 

       nperp = (/zero,cos(phi4),sin(phi4),zero/)
       ncoll = (/zero,zero,zero,one/)

       C4Lim%LimMom(:,1) = nperp
       C4Lim%LimMom(:,2) = ncoll

       C4Lim%KinInv = (/zero,s14,zero,zero,z4/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4Lim)

       C4Lim%AmpMom(:,1) = p1 - p4
       C4Lim%AmpMom(:,2) = p2
       C4Lim%AmpMom(:,3) = p5 


       C4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif


  end subroutine kinematics_nnlo_4151c


  subroutine kinematics_nnlo_4151d(MV2,Y,yr,&
       HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
       C54Lim,C54DSLim,C54S5Lim,C54DSS5Lim,TCC54Lim,TCC54DSLim,TCC54S5Lim,TCC54DSS5Lim)
    real(dp), intent(in) ::  yr(9), MV2, Y
    type(KinConfig)      ::  HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C54Lim,C54DSLim,C54DSS5Lim,C54S5Lim,TCC54Lim,TCC54DSLim,TCC54DSS5Lim,TCC54S5Lim
    !---
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV,  p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4151, s14, s15, s45, z4, z5, ztot, weight4151
    real(dp) :: spart,sqrts, x4p, nperp(4), s12, s24, z45, ncoll(1:4), yy, Deltas 
    real(dp) :: kt1kt2sq,z1,s1_45,zz1,zz2
    real(dp) :: xx3

    !    statement functions
    weight4151(eta41,eta42,eta51,eta52,eta54)  = eta42*eta52*( one + eta41/(eta54 + eta42 + eta51) & 
         + eta51/(eta54 + eta41 + eta52) )

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(DSLim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(DSS5Lim)        !4
    call initialize_config(TCLim)          !5
    call initialize_config(TCDSLim)        !6
    call initialize_config(TCS5Lim)        !7
    call initialize_config(TCDSS5Lim)      !8
    call initialize_config(C54Lim)          !9
    call initialize_config(C54DSLim)        !10
    call initialize_config(C54S5Lim)        !11
    call initialize_config(C54DSS5Lim)      !12
    call initialize_config(TCC54Lim)        !13
    call initialize_config(TCC54DSLim)      !14
    call initialize_config(TCC54S5Lim)      !15
    call initialize_config(TCC54DSS5Lim)    !16

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

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,x4p,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: DS FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate ; weight is the same 
    !-------------------------------------------------------------------------------------------------------------


    ! --------- DOUBLE SOFT  and DOUBLE SOFT + S5--------------
    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    DSLim%PartFrac = (/xi1,xi2/)

    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
    p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
         sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

    pVV = p1 + p2
    call get_lept_mom(MV,pVV,yr(7:9),kallenF,DSLim)


    DSLim%AmpMom(:,1) = p1
    DSLim%AmpMom(:,2) = p2

    DSLim%LimMom(:,1) = p1
    DSLim%LimMom(:,2) = p2
    DSLim%LimMom(:,3) = p4 
    DSLim%LimMom(:,4) = p5

    DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2


    DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2 

    ! ------   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT  ------------

    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero 

    s14 = spart *x1*x3*x4p 
    s15 = spart *x1*x2*x3
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

    ztot = - one 
    z4 = x1/ztot 
    z5 = x1*x2/ztot

    TCDSLim%KinInv  = (/s14,s15,s45,z4,z5/)


    TCDSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    ! ------   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT  + S5 ---------
    ztot = - one 
    z4 = x1/ztot
    z5 = x1*x2/ztot
    s14 = spart *x1*x3*x4p 
    s15 = spart *x1*x2*x3
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

    TCDSS5Lim%KinInv = (/s14,s15,s45,z4,z5/)

    TCDSS5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    !--------  C54 DS FLM  ----------

    eta41 = x3
    eta51 = x3
    eta42 = (one - eta41) 
    eta52 = (one - eta51)
    eta54 = zero 

    z5 = x2/(one+x2) 
    s12 = two*scr(p1,p2) 
    s14 = spart *x1*x3 
    s24 = spart * x1* (  one - x3)
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)

    kt1kt2sq = lambda 
    C54DSLim%KinInv = (/s12,s14,s24,s45,z5,kt1kt2sq/)

    C54DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,one,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


    ! ------  C54 +DS+S5 ---------
    z5 = x2
    s12 = two*scr(p1,p2) 
    s14 = spart *x1*x3 
    s24 = spart * x1* (  one - x3)
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)

    C54DSS5Lim%KinInv = (/s12,s45,s14,s24,z5/)
    C54DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,one,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2


    ! ----   TRIPLE COLLINEAR KINEMATICS  + C54 + Double soft  ---------

    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero 

    !--NEW
    s14 = spart *x1*x3 
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
    s1_45 = s14 * (one+x2)
    z5 = one/(one+x2)
    z1 = -one/(x1+x1*x2)
    kt1kt2sq = lambda * (one-xx3)

    TCC54DSLim%KinInv = (/s45,s1_45,z5,z1,kt1kt2sq/)


    TCC54DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,one,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    ! -------   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT5 +DS ----------

    z45 = one - x1
    z5 = zero 
    s14 = spart *x1*x3 
    s15 = spart *x1*x2*x3
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)

    TCC54DSS5Lim%KinInv = (/s14,s15, s45, z45, z5/)
    TCC54DSS5Lim%wgt = one/8.0_dp/pi * kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,one,lambda) & 
         * x1**3*x2 & 
         * weight4151(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 FLM   ! angles the same as in hard, do not regenerate; weight is the same 
    !-------------------------------------------------------------------------------------------------------------

    eta41 = x3*x4p
    eta51 = x3
    eta42 = (one-eta41) 
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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
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
            * x1**3*x2 & 
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
    Deltad = one - x1*(one-eta42)

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


       z5 = x2
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)
       C54S5Lim%KinInv = (/zero,s45,zero,zero,z5/)

       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C54S5Lim)
       C54S5Lim%AmpMom(:,1) = p1
       C54S5Lim%AmpMom(:,2) = p2
       C54S5Lim%AmpMom(:,3) = p4



       C54S5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1**3*x2 & 
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

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one,zero,zero,one /)

       s14 = spart *x1*x3*x4p 
       s15 = spart *x1*x2*x3
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)
       ztot = x1*(one +x2) - one 
       z4 = x1/ztot 
       z5 = x1*x2/ztot

       TCLim%KinInv = (/s14,s15, s45, z4, z5/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)


       TCLim%AmpMom(:,1) = p1 - p4 - p5
       TCLim%AmpMom(:,2) = p2


       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

       ! --------------  TC  + C54 ----------

       s14 = spart *x1*x3 
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
       s1_45 = s14 * (one+x2)
       zz1 = x2/(one+x2)
       zz2 = x1*(one+x2)/(x1*(one+x2)-one)
       kt1kt2sq = lambda * (one-xx3)

       TCC54Lim%KinInv = (/s45,s1_45,zz1,zz2,kt1kt2sq/)


       TCC54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1**3*x2 & 
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

       z4 = x1 
       z5 = x1*x2

       s14 = spart *x1*x3*x4p 
       s15 = spart *x1*x2*x3
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       TCS5Lim%KinInv = (/s14,s15,s45,z4,z5/)

       pVV = p1 + p2 - p4 
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1 - p4 
       TCS5Lim%AmpMom(:,2) = p2

       TCS5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2



       ! --------   TRIPLE COLLINEAR KINEMATICS  + C54 + SOFT5 ------------

       s14 = spart *x1*x3 
       s15 = spart *x1*x2*x3
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
       z45 = one - x1
       z5 = zero 

       TCC54S5Lim%KinInv = (/s14,s15, s45, z45, zero/)

       TCC54S5Lim%wgt = one/8.0_dp/pi * kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1**3*x2 & 
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

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C54Lim%flag = .true.
    if (C54Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C54Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = x2*p4 

       ncoll = (/ zero,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)

       !-- note a sign change in some terms nperp relative to sector b

       nperp  = (/ zero,  -cos4*sqrt(lambda)*cos(phi4) + sin(phi4)*sqrt(one-lambda) , & 
            -cos4*sqrt(lambda)*sin(phi4) - cos(phi4)*sqrt(one-lambda), & 
            sin4*sqrt(lambda) /)    


       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)
       z5 = x2/(one + x2) 

       C54Lim%KinInv = (/zero,zero,s45,zero,z5/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C54Lim)

       C54Lim%AmpMom(:,1) = p1
       C54Lim%AmpMom(:,2) = p2
       C54Lim%AmpMom(:,3) = p4 + p5 
       C54Lim%LimMom(:,1) = nperp
       C54Lim%LimMom(:,2) = ncoll

       C54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1**3*x2 & 
            * weight4151(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart & 
            * spart/MV2

    endif


  end subroutine kinematics_nnlo_4151d


  subroutine kinematics_nnlo_4251(MV2,Y,yr,&
       HardProc,DSLim,S5Lim,DSS5Lim,C4Lim,C4DSLim,C4S5Lim,C4DSS5Lim,&
       C5Lim,C5DSLim,C5S5Lim,C5DSS5Lim,C4C5Lim,C4C5DSLim,C4C5S5Lim,C4C5DSS5Lim)
    real(dp), intent(in) ::  yr(9), MV2, Y
    type(KinConfig)      :: HardProc,DSLim,DSS5Lim,S5Lim,C4Lim,C4DSLim,C4DSS5Lim,C4S5Lim,&
         C5Lim,C5DSLim,C5DSS5Lim,C5S5Lim,C4C5Lim,C4C5DSLim,C4C5S5Lim,C4C5DSS5Lim
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos4, sin4, phi4, phi5, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4151, s14, s15, s45, z4, z5, ztot, weight4251
    real(dp) :: spart,sqrts, s25, s12, s24, yy, Deltas
    real(dp) :: ncoll(4),nperp(4)
    real(dp) :: xx3

    !    statement functions
    weight4251(eta41,eta42,eta51,eta52,eta54)  = eta41*eta52*eta54/(eta54 + eta51 + eta42)

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(DSLim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(DSS5Lim)        !4
    call initialize_config(C4Lim)          !5
    call initialize_config(C4DSLim)        !6
    call initialize_config(C4S5Lim)        !7
    call initialize_config(C4DSS5Lim)      !8
    call initialize_config(C5Lim)          !9
    call initialize_config(C5DSLim)        !10
    call initialize_config(C5S5Lim)        !11
    call initialize_config(C5DSS5Lim)      !12
    call initialize_config(C4C5Lim)        !13
    call initialize_config(C4C5DSLim)      !14
    call initialize_config(C4C5S5Lim)      !15
    call initialize_config(C4C5DSS5Lim)    !16

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


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 


       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
            * weight4251(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: DS FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! -------------- DOUBLE SOFT AND DOUBLE SOFT + S5 --------------------

    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    DSLim%PartFrac = (/xi1,xi2/)

    p1 = half*sqrts*(/one,zero,zero,one/) 
    p2 = half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
    p5 = half*sqrts*x1*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)

    pVV = p1 + p2
    call get_lept_mom(MV,pVV,yr(7:9),kallenF,DSLim)

    DSLim%AmpMom(:,1) = p1
    DSLim%AmpMom(:,2) = p2

    DSLim%LimMom(:,1) = p1
    DSLim%LimMom(:,2) = p2
    DSLim%LimMom(:,3) = p4 
    DSLim%LimMom(:,4) = p5

    DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4251(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2

    DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4251(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2


    ! ------------   C42  +  DOUBLE SOFT   -----------

    eta41 = one
    eta51 = x4 
    eta42 = zero  
    eta52 = (one - x4) 
    eta54 =   eta52

    s24 = spart *x1*x3 
    s25 = spart*x1*x2*(one-x4) 
    s15 = spart*x1*x2*x4
    s12 = two*scr(p1,p2) 
    z4 = x1 

    C4DSLim%KinInv = (/s12,s15,s25,s24,z4/)

    C4DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4251(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


    ! -------------   C42 +  DOUBLE SOFT  + S5 ---------------

    z4 = x1

    s24 = spart *x1*x3 
    s25 = spart*x1*x2*(one-x4) 
    s15 = spart*x1*x2*x4
    s12 = two*scr(p1,p2) 

    C4DSS5Lim%KinInv = (/s12,s15,s25,s24,z4/)

    C4DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4251(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2




    ! ------------  C51 DS FLM  ------------

    phi4 = twopi*yr(4) 
    cos4 = -one + two*x3 
    sin4 = sqrt(one-cos4**2)

    eta42 = x3
    eta52 = one
    eta41 = (one-x3) 
    eta51 = zero
    eta54 = (one-x3) 

    s12 = two*scr(p1,p2) 
    s24 = spart*x1*x3
    s14 = spart*x1*(one - x3) 

    s15 = spart*x1*x2*(x4)
    z5 = x1*x2

    C5DSLim%KinInv = (/s12,s14,s24,s15,z5/)

    C5DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4251(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


    ! -----------  C51+DS+S5 -----------------

    s12 = two*scr(p1,p2) 
    s24 = spart*x1*x3
    s14 = spart*x1*(one - x3) 

    s15 = spart*x1*x2*(x4)
    z5 = x1*x2

    C5DSS5Lim%KinInv = (/s12,s14,s24,s15,z5/)

    C5DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4251(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


    ! -----------   C42  + C51 + Double soft  ----------


    eta42 = zero
    eta52 = one
    eta41 = one
    eta51 = zero
    eta54 = one

    s24 = spart *x1*x3 
    s15 = spart *x1*x2*x4
    z4 = x1 
    z5 = x1*x2

    C4C5DSLim%KinInv = (/s24,s15,zero,z4,z5/)

    C4C5DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4251(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


    ! ----------   C41  + C52 + SOFT5 +DS ------------

    z4 = x1 
    z5 = x1*x2
    s24 = spart *x1*x3 
    s15 = spart *x1*x2*x4

    C4C5DSS5Lim%KinInv = (/s24,s15,zero,z4,z5/)

    C4C5DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4251(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2






    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3: S5 FLM   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------


    ! --------------- SOFT 5 ----------

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
       p5 = half*sqrts*x1*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)


       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p4 

       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = p4 
       S5Lim%LimMom(:,4) = p5 

       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
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


       s15 = spart*x1*x2*x4
       z5 = x1*x2

       C5S5Lim%KinInv = (/zero,s15,zero,zero,z5/)

       C5S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
            * weight4251(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif



    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3:   C42  sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------



    eta41 = one
    eta51 = x4 
    eta42 = zero  
    eta52 = (one - x4) 
    eta54 =   eta52


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C4Lim%flag  = .true.
    if (C4Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4Lim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)


       s24 = spart *x1*x3
       z4 = one/(one-x1)
       
       C4Lim%KinInv = (/s24,zero,zero, z4, zero/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4Lim)

       C4Lim%AmpMom(:,1) = p1 
       C4Lim%AmpMom(:,2) = p2 - p4
       C4Lim%AmpMom(:,3) = p5

       nperp = (/zero,cos(phi4),sin(phi4),zero/)
       ncoll = (/zero,zero,zero,-one/)

       C4Lim%LimMom(:,1) = nperp
       C4Lim%LimMom(:,2) = ncoll

       C4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
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
       s25 = spart*x1*x2*(one-x4) 
       s15 = spart*x1*x2*x4
       s12 = two*scr(p1,p2) 

       z4 = one - x1 

       C4S5Lim%KinInv = (/s12,s15,s25,s24,z4/)

       pVV = p1 + p2 - p4 
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4S5Lim)

       C4S5Lim%AmpMom(:,1) = p1 
       C4S5Lim%AmpMom(:,2) = p2 - p4


       C4S5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
            * weight4251(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ---------   C42   + C51 + SOFT5 ------------


       eta42 = zero
       eta52 = one
       eta41 = one
       eta51 = zero
       eta54 = one           

       z4 = one/(one-x1) 
       z5 = x1*x2

       s24 = spart *x1*x3 
       s15 = spart *x1*x2*x4

       C4C5S5Lim%KinInv = (/s24,s15,zero,z4,z5/)

       C4C5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
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


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C5Lim%flag = .true.
    if (C5Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C5Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, zero, zero, one /)

       s15 = spart*x1*x2*x4
       z5 = one/(one - x1 * x2) 

       C5Lim%KinInv = (/zero,s15,zero,zero,z5/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C5Lim)

       C5Lim%AmpMom(:,1) = p1 - p5 
       C5Lim%AmpMom(:,2) = p2 
       C5Lim%AmpMom(:,3) = p4 

       nperp = (/zero,cos(phi5),sin(phi5),zero/)
       ncoll = (/zero,zero,zero,one/)

       C5Lim%LimMom(:,1) = nperp
       C5Lim%LimMom(:,2) = ncoll

       C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
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

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C4C5Lim%flag = .true.
    if (C4C5Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4C5Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one,zero,zero,one /)

       s24 = spart *x1*x3 
       s15 = spart *x1*x2*x4
       z4 = one/(one - x1)
       z5 = one/(one - x1*x2) 

       C4C5Lim%KinInv = (/s24,s15,zero,z4,z5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4C5Lim)

       C4C5Lim%AmpMom(:,1) = p1 - p5 
       C4C5Lim%AmpMom(:,2) = p2 - p4

       C4C5Lim%LimMom(:,1) = (/zero,cos(phi4),sin(phi4),zero/) !-- nperp
       C4C5Lim%LimMom(:,2) = (/zero,zero,zero,-one/) !-- ndir
       C4C5Lim%LimMom(:,3) = (/zero,cos(phi5),sin(phi5),zero/) !-- nperp
       C4C5Lim%LimMom(:,4) = (/zero,zero,zero,one/) !-- ndir

       C4C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
            * weight4251(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


    endif


  end subroutine kinematics_nnlo_4251


  !-------------------------------------------------------------------------------------------------

  !-- symmetric sectors below

  subroutine kinematics_nnlo_4252a(MV2,Y,yr,&
       HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
       C5Lim,C5DSLim,C5S5Lim,C5DSS5Lim,TCC5Lim,TCC5DSLim,TCC5S5Lim,TCC5DSS5Lim)       
    real(dp), intent(in) ::  MV2,Y, yr(9)
    type(KinConfig)      ::  HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C5Lim,C5DSLim,C5DSS5Lim,C5S5Lim,TCC5Lim,TCC5DSLim,TCC5DSS5Lim,TCC5S5Lim
    !---
    real(dp) :: p1(4), p2(4), p4(4), p5(4)
    real(dp) :: MV, kallenF, sin5, cos5, pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4252, S24, s25, s45, z4, z5, ztot, weight4252
    real(dp) :: spart,sqrts, yy, Deltas
    real(dp) :: ncoll(4),nperp(4)
    real(dp) :: mysin4lim,mysin5lim,mycphi,z45
    real(dp) :: zz1,zz2,s4_25,kt1kt2sq
    real(dp) :: xx3

    !    statement functions
    weight4252(eta41,eta42,eta51,eta52,eta54)  = ETA41*ETA51*( one + ETA42/(eta54 + ETA41 + ETA52) & 
         + ETA52/(eta54 + ETA42 + ETA51) )

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(DSLim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(DSS5Lim)        !4
    call initialize_config(TCLim)          !5
    call initialize_config(TCDSLim)        !6
    call initialize_config(TCS5Lim)        !7
    call initialize_config(TCDSS5Lim)      !8
    call initialize_config(C5Lim)          !9
    call initialize_config(C5DSLim)        !10
    call initialize_config(C5S5Lim)        !11
    call initialize_config(C5DSS5Lim)      !12
    call initialize_config(TCC5Lim)        !13
    call initialize_config(TCC5DSLim)      !14
    call initialize_config(TCC5S5Lim)      !15
    call initialize_config(TCC5DSS5Lim)    !16

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

    ETA42 = x3
    ETA52 = half*x3*x4
    ETA41 = (one-x3) 
    ETA51 = (one - half*x3*x4)
    eta54 = x3*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif
    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: DS FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! ------- DOUBLE SOFT and DOUBLE SOFT + S5 ------------

    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    DSLim%PartFrac = (/xi1,xi2/)

    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
    p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
         sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

    pVV = p1 + p2

    call get_lept_mom(MV,pVV,yr(7:9),kallenF,DSLim)

    DSLim%AmpMom(:,1) = p1
    DSLim%AmpMom(:,2) = p2

    DSLim%LimMom(:,1) = p1
    DSLim%LimMom(:,2) = p2
    DSLim%LimMom(:,3) = p4 
    DSLim%LimMom(:,4) = p5

    DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2  
    ! the weight for DS + S5 (identical to above, but avoids confusion)
    DSS5Lim%wgt= DSLim%wgt


    ! ------   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT  -------

    ETA42 = zero
    ETA52 = zero
    ETA41 = one
    ETA51 = one
    eta54 = zero

    p1=  half*sqrts*(/one,zero,zero,one/)
    p2 = half*sqrts*(/one,zero,zero,-one/)
    p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)
    p5 = half*sqrts*x1*x2*(/one, zero,zero,one /)

    ztot = - one 
    z4 = x1/ztot 
    z5 = x1*x2/ztot

    S24 = spart *x1*x3 
    S25 = spart *x1*x2*(half*x3*x4)
    s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)


    TCDSLim%KinInv = (/S24,S25,s45,z4,z5/)

    pVV = p1 + p2

    call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCDSLim)

    TCDSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         *one/two/spart &
         *spart/MV2

    ! ------   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT  + S5 -------

    ztot = - one 
    z4 = x1/ztot
    z5 = x1*x2/ztot
    S24 = spart *x1*x3 
    s25 = spart *x1*x2*(half*x3*x4)
    s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

    TCDSS5Lim%KinInv = (/S24,s25,s45,z4,z5/)

    TCDSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         *one/two/spart &
         *spart/MV2

    ! --------  C5 DS FLM  ----------

    ETA42 = x3
    ETA52 = zero
    ETA41 = (one-x3)
    ETA51 = one
    eta54 = x3/NF1(x3,zero,lambda)


    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)
    p5 = half*sqrts*x1*x2*(/one, zero, zero, one /)

    ztot = - one
    z4 = x1/ztot
    z5 = x1*x2/ztot

    S24 = spart *x1*x3
    s25 = spart *x1*x2*(half*x3*x4)
    s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

    C5DSLim%KinInv = (/S24,s25,s45,z4,z5/)

    pVV = p1 + p2


    call get_lept_mom(MV,pVV,yr(7:9),kallenF,C5DSLim)

    C5DSLim%LimMom(:,1) = p1 
    C5DSLim%LimMom(:,2) = p2
    C5DSLim%LimMom(:,3) = p4 

    C5DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         *one/two/spart &
         * spart/MV2


    ! ----------- C5+DS+S5 --------------

    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

    s25 = spart*x1*x2*(half*x3*x4)
    z5 = x1*x2

    C5DSS5Lim%KinInv = (/zero,s25,zero,zero,z5/)

    pVV = p1 + p2

    call get_lept_mom(MV,pVV,yr(7:9),kallenF,C5DSS5Lim)


    C5DSS5Lim%LimMom(:,1) = p1
    C5DSS5Lim%LimMom(:,2) = p2
    C5DSS5Lim%LimMom(:,3) = p4 

    C5DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2


    ! -------   TRIPLE COLLINEAR KINEMATICS  + C5 + Double soft  --------------
    ETA42 = zero
    ETA52 = zero
    ETA41 = one
    ETA51 = one
    eta54 = zero


    S24 = spart *x1*x3 
    s25 = spart *x1*x2*(half*x3*x4)
    s45 = spart *x1**2*x2*x3/NF1(xx3,zero,lambda)
    ztot =-one 
    z4 = x1/ztot 
    z5 = x1*x2/ztot

    TCC5DSLim%KinInv = (/S24,s25,s45,z4,z5/)

    TCC5DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(xx3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2

    !---------   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT5 +DS ------------


    z4 = x1 
    z5 = x1*x2
    S24 = spart *x1*x3 
    s25 = spart *x1*x2*(half*x3*x4)
    s45 = spart *x1**2*x2*x3/NF1(xx3,zero,lambda)

    TCC5DSS5Lim%KinInv = (/S24,s25,s45,z4,z5/)

    TCC5DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(xx3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 + C5S5   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------

    ! --------- SOFT 5 ---------
    ETA42 = x3
    ETA52 = half*x3*x4
    ETA41 = (one-x3)
    ETA51 = (one - half*x3*x4)
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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
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
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2   


       ! kin inv and weight for C5 +  S5
       ETA42 = x3
       ETA52 = zero
       ETA41 = (one-x3)
       ETA51 = one
       eta54 = x3/NF1(x3,zero,lambda)

       s25 = spart*x1*x2*(half*x3*x4)
       z5 = x1*x2
       C5S5Lim%KinInv = (/zero,s25,zero,zero,z5/)

       C5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif
    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    ! --------- TRIPLE COLLINEAR -------

    ETA42 = zero
    ETA52 = zero
    ETA41 = one
    ETA51 = one
    eta54 = zero 


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one,zero,zero, -one /)

       S24 = spart *x1*x3 
       s25 = spart *x1*x2*(half*x3*x4)
       s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

       ztot = x1*(one +x2) - one 
       z4 = x1/ztot 
       z5 = x1*x2/ztot
       TCLim%KinInv = (/S24,s25, s45, z4, z5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)

       TCLim%AmpMom(:,1) = p1
       TCLim%AmpMom(:,2) = p2 - p4 - p5

       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2


       ! --------- TRIPLE COLLINEAR + C5  -------

       s24 = spart *x1*x3 
       s25 = spart *x1*x2*(half*x3*x4)
       s4_25 = s24 * (one-x1*x2)

       zz1 = x1*x2/(x1*x2-one)
       zz2 = (one-x1*x2)/(one-x1-x1*x2)

       kt1kt2sq = (-one+two*lambda)**2 * (one-xx3)

       TCC5Lim%KinInv = (/s25,s4_25,zz1,zz2,kt1kt2sq/)

       TCC5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1**3*x2 & 
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
       p5 = half*sqrts*x1*x2*(/one,zero,zero,-one /)


       z4 = x1 
       z5 = x1*x2

       S24 = spart *x1*x3 
       s25 = spart *x1*x2*(half*x3*x4)
       s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

       TCS5Lim%KinInv = (/S24,s25,s45,z4,z5/)

       pVV = p1 + p2 - p4 

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1
       TCS5Lim%AmpMom(:,2) = p2 - p4 

       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            *spart/MV2


       ! -------   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT5 ---------
       z4 = x1 
       z5 = x1*x2

       S24 = spart *x1*x3 
       s25 = spart *x1*x2*(half*x3*x4)
       s45 = spart *x1**2*x2*x3/NF1(xx3,zero,lambda)

       TCC5S5Lim%KinInv = (/S24,s25,s45,z4,z5/)

       TCC5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C5 FLM 
    !-------------------------------------------------------------------------------------------------------------


    ETA42 = x3
    ETA52 = zero
    ETA41 = (one-x3)
    ETA51 = one
    eta54 = x3/NF1(x3,zero,lambda)

    phi4 = twopi*yr(4) 
    cos4 = -one + two*x3 
    sin4 = sqrt(one-cos4**2)

    ETA42 = x3
    ETA52 = zero
    ETA41 = (one-x3) 
    ETA51 = one
    eta54 = x3/NF1(x3,zero,lambda)


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C5Lim%flag = .true.
    if (C5Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C5Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, zero, zero, -one /)

       s25 = spart*x1*x2*(half*x3*x4)
       z5 = one/(one - x1 * x2) 

       C5Lim%KinInv = (/zero,s25,zero,zero,z5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C5Lim)

       C5Lim%AmpMom(:,1) = p1
       C5Lim%AmpMom(:,2) = p2 - p5
       C5Lim%AmpMom(:,3) = p4

              !-- spin-correlation vectors

       cosphi45 = -one + two * lambda
       sinphi45 = sqrt(one-cosphi45**2)

       nperp = (/zero,(cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            (sin(phi4)*cosphi45+cos(phi4)*sinphi45),zero/)

       ncoll = (/zero,zero,zero,-one/)

       C5Lim%LimMom(:,1) = nperp
       C5Lim%LimMom(:,2) = ncoll

       C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif


  end subroutine kinematics_nnlo_4252a


  subroutine kinematics_nnlo_4252b(MV2,Y,yr,&
       HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
       C54Lim,C54DSLim,C54S5Lim,C54DSS5Lim,TCC54Lim,TCC54DSLim,TCC54S5Lim,TCC54DSS5Lim) 
    real(dp), intent(in) ::  yr(9), MV2, Y
    type(KinConfig)      ::  HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C54Lim,C54DSLim,C54DSS5Lim,C54S5Lim,TCC54Lim,TCC54DSLim,TCC54DSS5Lim,TCC54S5Lim
    !---
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, Deltas, w4252, s14, s25, s45, z4, z5, ztot, weight4252
    real(dp) :: spart,sqrts, x4p, nperp(4), s12, s24, z45, ncoll(1:4), yy 
    real(dp) :: kt1kt2sq,zz1,zz2,s2_45,z1,z2
    real(dp) :: xx3

    !    statement functions
    weight4252(eta41,eta42,eta51,eta52,eta54)  = ETA41*ETA51*( one + ETA42/(eta54 + ETA41 + ETA52) & 
         + ETA52/(eta54 + ETA42 + ETA51) )


    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(DSLim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(DSS5Lim)        !4
    call initialize_config(TCLim)          !5
    call initialize_config(TCDSLim)        !6
    call initialize_config(TCS5Lim)        !7
    call initialize_config(TCDSS5Lim)      !8
    call initialize_config(C54Lim)          !9
    call initialize_config(C54DSLim)        !10
    call initialize_config(C54S5Lim)        !11
    call initialize_config(C54DSS5Lim)      !12
    call initialize_config(TCC54Lim)        !13
    call initialize_config(TCC54DSLim)      !14
    call initialize_config(TCC54S5Lim)      !15
    call initialize_config(TCC54DSS5Lim)    !16

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

    ETA42 = x3
    ETA52 = x3*x4p
    ETA41 = (one-ETA42) 
    ETA51 = (one - ETA52)
    eta54 = x3*x4**2/four/NF1(x3,x4p,lambda)


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,x4p,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart  & 
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: DS FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate ; weight is the same 
    !-------------------------------------------------------------------------------------------------------------

    !   DS + DS S5 FLM   ! angles etc. are the same ; weight is the same 
    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    DSLim%PartFrac = (/xi1,xi2/)

    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
    p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
         sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

    pVV = p1 + p2
    call get_lept_mom(MV,pVV,yr(7:9),kallenF,DSLim)

    DSLim%AmpMom(:,1) = p1
    DSLim%AmpMom(:,2) = p2

    DSLim%LimMom(:,1) = p1
    DSLim%LimMom(:,2) = p2
    DSLim%LimMom(:,3) = p4 
    DSLim%LimMom(:,4) = p5

    DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    ! -----   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT   --------

    ETA42 = zero
    ETA52 = zero
    ETA41 = one
    ETA51 = one
    eta54 = zero 

    S24 = spart *x1*x3 
    S25 = spart *x1*x2*x3*x4p
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)
    ztot = - one 
    z4 = x1/ztot 
    z5 = x1*x2/ztot

    TCDSLim%KinInv = (/S24,S25,s45,z4,z5/)

    TCDSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    !-------------------------------------------------------------------------------------------------------------
    !   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT  + S5
    !-------------------------------------------------------------------------------------------------------------

    S24 = spart *x1*x3 
    S25 = spart *x1*x2*x3*x4p
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

    ztot = - one 
    z4 = x1/ztot
    z5 = x1*x2/ztot

    TCDSS5Lim%KinInv = (/S24,S25,s45,z4,z5/)

    TCDSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    ! --------  C54 DS FLM  ----------
    ETA42 = x3
    ETA52 = x3
    ETA41 = (one - ETA42) 
    ETA51 = (one - ETA52)
    eta54 = zero

    z5 = x2/(one+x2) 

    s12 = two*scr(p1,p2) 
    s24 = spart *x1*x3 
    s14 = spart * x1* (  one - x3)
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)

    kt1kt2sq = lambda !* (one-xx3)
    C54DSLim%KinInv = (/s12,s24,s14,s45,z5,kt1kt2sq/)

    C54DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,one,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


    !------------  C54 +DS+S5 ---------------

    s12 = two*scr(p1,p2) 
    S24 = spart *x1*x3 
    S14 = spart * x1* (  one - x3)
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)
    z5 = x2 

    C54DSS5Lim%KinInv = (/s12,s45,S24,S14,z5/)
    C54DSS5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,one,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    ! ---------   TC  + C54 + Double soft  ----------

    ETA42 = zero
    ETA52 = zero
    ETA41 = one
    ETA51 = one
    eta54 = zero 

    !--NEW
    s24 = spart *x1*x3 
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
    s2_45 = s24 * (one+x2)
    z5 = one/(one+x2)
    z2 = -one/(x1+x1*x2)
    kt1kt2sq = lambda * (one-xx3)

    TCC54DSLim%KinInv = (/s45,s2_45,z5,z2,kt1kt2sq/)

    TCC54DSLim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,one,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    !-------------------------------------------------------------------------------------------------------------
    !   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT5 +DS
    !-------------------------------------------------------------------------------------------------------------

    S24 = spart *x1*x3 
    S25 = spart *x1*x2*x3
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
    z45 = one - x1
    z5 = x2
    TCC54DSS5Lim%KinInv = (/S24,S25, s45, z45, z5/)
    TCC54DSS5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,one,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2  


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 FLM + C54 S5 FLM   ! angles the same as in hard, do not regenerate; weight is the same 
    !-------------------------------------------------------------------------------------------------------------

    ! ---------- SOFT 5------------
    ETA42 = x3
    ETA52 = x3*x4p
    ETA41 = (one-ETA42) 
    ETA51 = (one - ETA52)
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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
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
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ------  C54 + S5 ---------

       ETA42 = x3
       ETA52 = x3
       ETA41 = (one - ETA42) 
       ETA51 = (one - ETA52)
       eta54 = zero

       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)
       z5 = x2

       C54S5Lim%KinInv = (/zero,s45,zero,zero,z5/)

       C54S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif






    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    ! --------- TRIPLE COLLINEAR -------

    ETA42 = zero
    ETA52 = zero
    ETA41 = one
    ETA51 = one
    eta54 = zero 
    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = x2*p4

       S24 = spart *x1*x3 
       S25 = spart *x1*x2*x3*x4p
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       ztot = x1*(one +x2) - one 
       z4 = x1/ztot 
       z5 = x1*x2/ztot

       TCLim%KinInv= (/S24,S25,s45,z4,z5/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)

       TCLim%AmpMom(:,1) = p1
       TCLim%AmpMom(:,2) = p2 - p4 - p5

       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! -------   TC  + C54 -------------

       s24 = spart *x1*x3 
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
       s2_45 = s24 * (one+x2)
       zz1 = x2/(one+x2)
       zz2 = x1*(one+x2)/(x1*(one+x2)-one)
       kt1kt2sq = lambda * (one-xx3)

       TCC54Lim%KinInv = (/s45,s2_45,zz1,zz2,kt1kt2sq/)

       TCC54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1**3*x2 & 
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

       z4 = x1 
       z5 = x1*x2
       S24 = spart *x1*x3 
       S25 = spart *x1*x2*x3*x4p
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       TCS5Lim%KinInv = (/S24,S25,s45,z4,z5/)
       pVV = p1 + p2 - p4 
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1
       TCS5Lim%AmpMom(:,2) = p2 - p4 


       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ------   TRIPLE COLLINEAR KINEMATICS  + C54 + SOFT5 -------

       z45 = one - x1

       S24 = spart *x1*x3 
       S25 = spart *x1*x2*x3
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)

       TCC54S5Lim%KinInv = (/S24,S25, s45, z45, zero/)

       TCC54S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif



    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C54 FLM 
    !-------------------------------------------------------------------------------------------------------------

    cos4 = -one + two*x3 
    sin4 = sqrt(one-cos4**2)

    ETA42 = x3
    ETA52 = x3
    ETA41 = (one - ETA42) 
    ETA51 = (one - ETA52)
    eta54 = zero

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C54Lim%flag = .true.
    if (C54Lim%flag.eqv. .false.) spart = MV2/Deltas
    if (C54Lim%flag.eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = x2*p4 

       ncoll = (/ zero,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)

       !-- be careful, different from 14,15
       nperp  = (/ zero,  -cos4*sqrt(lambda)*cos(phi4) + sin(phi4)*sqrt(one-lambda) , & 
            -cos4*sqrt(lambda)*sin(phi4) - cos(phi4)*sqrt(one-lambda), & 
            sin4*sqrt(lambda) /)

       pVV = p1 + p2 - p4 - p5

       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)

       z5 = x2/(one + x2) 

       C54Lim%KinInv = (/zero,zero,s45,zero,z5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C54Lim)


       C54Lim%AmpMom(:,1) = p1
       C54Lim%AmpMom(:,2) = p2
       C54Lim%AmpMom(:,3) = p4 + p5 

       C54Lim%LimMom(:,1) = nperp
       C54Lim%LimMom(:,2) = ncoll

       C54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif


  end subroutine kinematics_nnlo_4252b


  subroutine kinematics_nnlo_4252c(MV2,Y,yr, &
       HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
       C4Lim,C4DSLim,C4S5Lim,C4DSS5Lim,TCC4Lim,TCC4DSLim,TCC4S5Lim,TCC4DSS5Lim)
    real(dp), intent(in) ::  MV2, Y, yr(9)
    type(KinConfig) ::  HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C4Lim,C4DSLim,C4DSS5Lim,C4S5Lim,TCC4Lim,TCC4DSLim,TCC4DSS5Lim,TCC4S5Lim
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4) 
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5,  pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4252, S24, S25, s45, z4, z5, ztot, weight4252
    real(dp) :: spart,sqrts, s12, S15, z45, yy, Deltas
    real(dp) :: nperp(4), ncoll(4), mycphi, z1, z2, s14_5
    real(dp) :: s5_24,zz1,zz2,kt1kt2sq
    real(dp) :: xx3
        
    !    statement functions
    weight4252(eta41,eta42,eta51,eta52,eta54)  = ETA41*ETA51*( one + ETA42/(eta54 + ETA41 + ETA52) & 
         + ETA52/(eta54 + ETA42 + ETA51) )

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(DSLim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(DSS5Lim)        !4
    call initialize_config(TCLim)          !5
    call initialize_config(TCDSLim)        !6
    call initialize_config(TCS5Lim)        !7
    call initialize_config(TCDSS5Lim)      !8
    call initialize_config(C4Lim)          !9
    call initialize_config(C4DSLim)        !10
    call initialize_config(C4S5Lim)        !11
    call initialize_config(C4DSS5Lim)      !12
    call initialize_config(TCC4Lim)        !13
    call initialize_config(TCC4DSLim)      !14
    call initialize_config(TCC4S5Lim)      !15
    call initialize_config(TCC4DSS5Lim)    !16

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

    ETA42 = half*x3*x4
    ETA52 = x3
    ETA41 = (one-ETA42) 
    ETA51 = (one - ETA52)
    eta54 = x3*(one-half*x4)**2/NF1(x3,half*x4,lambda)


    !  incoming momenta 

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2


    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: DS FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! ------- DOUBLE SOFT and DOUBLE SOFT + S5 ------------

    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    DSLIM%PartFrac = (/xi1,xi2/)

    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
    p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
         sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

    pVV = p1 + p2


    call get_lept_mom(MV,pVV,yr(7:9),kallenF,DSLim)

    DSLim%AmpMom(:,1) = p1
    DSLim%AmpMom(:,2) = p2

    DSLim%LimMom(:,1) = p1
    DSLim%LimMom(:,2) = p2
    DSLim%LimMom(:,3) = p4 
    DSLim%LimMom(:,4) = p5

    DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - half*x4)*two/NF1(x3,half*x4,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    ! identical to DS but avoids confusion 
    DSS5Lim%wgt = DSLim%wgt


    ! -----   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT   ----------

    ETA42 = zero
    ETA52 = zero
    ETA41 = one
    ETA51 = one
    eta54 = zero

    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one, zero, zero, one  /)
    p5 = half*sqrts*x1*x2*(/one,zero,zero,one /)

    S24 = spart *x1*(half*x3*x4) 
    S25 = spart *x1*x2*x3
    s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

    ztot = - one 
    z4 = x1/ztot 
    z5 = x1*x2/ztot

    TCDSLim%KinInv = (/S24,S25,s45,z4,z5/)

    pVV = p1 + p2

    call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCDSLim)

    TCDSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         *one/two/spart &
         * spart/MV2

    ! -----   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT  + S5 --------

    ztot = - one 
    z4 = x1/ztot
    z5 = x1*x2/ztot

    S24 = spart *x1*(half*x3*x4) 
    S25 = spart *x1*x2*x3
    s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

    TCDSS5Lim%KinInv = (/S24,S25,s45,z4,z5/)

    TCDSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2



    ! -------  C41 DS FLM  ----------------


    ETA42 = zero
    ETA52 = x3
    ETA41 = one 
    ETA51 = (one - ETA52)
    eta54 = x3/NF1(x3,zero,lambda)


    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one, zero, zero, one /)

    s12 = two*scr(p1,p2) 
    S24 = spart *x1*(half*x3*x4) 
    S25 = spart *x1*x2*x3
    S15 = spart *x1*x2*(one-x3)

    z4 = x1

    C4DSLim%KinInv = (/s12,S25,S15,S24,z4/)

    C4DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         *one/two/spart & 
         * spart/MV2


    ! -------  C41+DS+S5 ----------

    s12 = two*scr(p1,p2) 
    S24 = spart *x1*(half*x3*x4) 
    S25 = spart *x1*x2*x3
    S15 = spart *x1*x2*(one-x3)

    z4 = x1
    C4DSS5Lim%KinInv = (/s12,S25,S15,S24,z4/)

    C4DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         *one/two/spart & 
         * spart/MV2


    ! -------    TRIPLE COLLINEAR KINEMATICS  + C41 + Double soft  --------------


    ETA42 = zero
    ETA52 = zero
    ETA41 = one
    ETA51 = one
    eta54 = zero 

    S24 = spart *x1*(half * x3 * x4)  
    S25 = spart *x1*x2* x3
    s45 = spart *x1**2*x2*x3/NF1(xx3,zero,lambda)

    z4 = x1
    S15 = spart *x1*x2

    TCC4DSLim%KinInv = (/s12,S24,S25,S15,z4/)

    TCC4DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(xx3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2



    ! -----------   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT5 +DS --------------

    S24 = spart *x1*(half * x3 * x4)  
    S25 = spart *x1*x2* x3
    s45 = spart *x1**2*x2*x3/NF1(xx3,zero,lambda)

    z4 = one - x1 
    z5 = x1*x2

    TCC4DSS5Lim%KinInv= (/S24,S25,s45,z4,z5/)

    TCC4DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(xx3,zero,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 FLM   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------

    ETA42 = half*x3*x4
    ETA52 = x3
    ETA41 = (one-ETA42) 
    ETA51 = (one - ETA52)
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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
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
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    ! --------- TRIPLE COLLINEAR -------

    ETA42 = zero
    ETA52 = zero
    ETA41 = one
    ETA51 = one
    eta54 = zero 


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one,zero,zero,-one /)

       S24 = spart *x1*(half*x3*x4) 
       S25 = spart *x1*x2*x3
       s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)

       ztot = x1*(one +x2) - one 
       z4 = x1/ztot 
       z5 = x1*x2/ztot

       TCLim%KinInv = (/S24,S25, s45, z4, z5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)

       TCLim%AmpMom(:,1) = p1
       TCLim%AmpMom(:,2) = p2 - p4 - p5

       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            * spart/MV2


       ! ---------   TC  + C41 ----------------

       s24 = spart *x1*(half * x3 * x4)  
       s25 = spart *x1*x2* x3
       s5_24 = s25 * (one-x1)
       
       zz1 = x1/(x1-one)
       zz2 = (one-x1)/(one-x1-x1*x2)

       kt1kt2sq = (-one+two*lambda)**2 * (one-xx3)

       TCC4Lim%KinInv = (/s24,s5_24,zz1,zz2,kt1kt2sq/)

       TCC4Lim%wgt  = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1**3*x2 & 
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

       S24 = spart *x1*(half*x3*x4) 
       S25 = spart *x1*x2*x3
       s45 = spart *x1**2 * x2 * x3 * (one-half*x4)**2/NF1(xx3,half*x4,lambda)


       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, zero, zero, -one  /)
       p5 = half*sqrts*x1*x2*(/one,zero,zero,-one /)

       z4 = x1 
       z5 = x1*x2

       TCS5Lim%KinInv = (/S24,S25,s45,z4,z5/)

       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1
       TCS5Lim%AmpMom(:,2) = p2 - p4 

       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - half*x4)*two/NF1(xx3,half*x4,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            * spart/MV2


       ! -------   TRIPLE COLLINEAR KINEMATICS  + C4 + SOFT5 -----------

       S24 = spart *x1*(half * x3 * x4)  
       S25 = spart *x1*x2* x3
       s45 = spart *x1**2*x2*x3/NF1(xx3,zero,lambda)

       z4 = one - x1 
       z5 = x1*x2

       TCC4S5Lim%KinInv = (/S24,S25,s45,z4,z5/)


       TCC4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(xx3,zero,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ----------  C41 + S5 


       ETA42 = zero
       ETA52 = x3
       ETA41 = one 
       ETA51 = (one - ETA52)
       eta54 = x3/NF1(x3,zero,lambda)


       z4 = one - x1
       s12 = two*scr(p1,p2) 
       S24 = spart *x1*(half*x3*x4) 
       S25 = spart *x1*x2*x3
       S15 = spart *x1*x2*(one-x3)

       C4S5Lim%KinInv = (/s12,S25,S15,S24,z4/)

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4S5Lim)

       C4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1**3*x2 & 
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


    ETA42 = zero
    ETA52 = x3
    ETA41 = one 
    ETA51 = (one - ETA52)
    eta54 = x3/NF1(x3,zero,lambda)


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C4Lim%flag  = .true.
    if (C4Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       S24 = spart*x1*(half*x3*x4)
       z4 = one/(one - x1 )

       nperp = (/zero,cos(phi4),sin(phi4),zero/)
       ncoll = (/zero,zero,zero,-one/)

       C4Lim%LimMom(:,1) = nperp
       C4Lim%LimMom(:,2) = ncoll

       C4Lim%KinInv = (/zero,S24,zero,zero,z4/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4Lim)

       C4Lim%AmpMom(:,1) = p1
       C4Lim%AmpMom(:,2) = p2 - p4
       C4Lim%AmpMom(:,3) = p5 


       C4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif



  end subroutine kinematics_nnlo_4252c


  subroutine kinematics_nnlo_4252d(MV2,Y,yr,&
       HardProc,DSLim,S5Lim,DSS5Lim,TCLim,TCDSLim,TCS5Lim,TCDSS5Lim,&
       C54Lim,C54DSLim,C54S5Lim,C54DSS5Lim,TCC54Lim,TCC54DSLim,TCC54S5Lim,TCC54DSS5Lim)
    real(dp), intent(in) ::  yr(9), MV2, Y
    type(KinConfig)      ::  HardProc,DSLim,DSS5Lim,S5Lim,TCLim,TCDSLim,TCDSS5Lim,TCS5Lim,&
         C54Lim,C54DSLim,C54DSS5Lim,C54S5Lim,TCC54Lim,TCC54DSLim,TCC54DSS5Lim,TCC54S5Lim
    !---
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV,  p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos4, sin4, phi4, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4252, S24, S25, s45, z4, z5, ztot, weight4252
    real(dp) :: spart,sqrts, x4p, nperp(4), s12, S14, z45, ncoll(1:4), yy, Deltas 
    real(dp) :: kt1kt2sq,z2,s2_45,zz1,zz2
    real(dp) :: xx3

    !    statement functions
    weight4252(eta41,eta42,eta51,eta52,eta54)  = ETA41*ETA51*( one + ETA42/(eta54 + ETA41 + ETA52) & 
         + ETA52/(eta54 + ETA42 + ETA51) )

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(DSLim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(DSS5Lim)        !4
    call initialize_config(TCLim)          !5
    call initialize_config(TCDSLim)        !6
    call initialize_config(TCS5Lim)        !7
    call initialize_config(TCDSS5Lim)      !8
    call initialize_config(C54Lim)          !9
    call initialize_config(C54DSLim)        !10
    call initialize_config(C54S5Lim)        !11
    call initialize_config(C54DSS5Lim)      !12
    call initialize_config(TCC54Lim)        !13
    call initialize_config(TCC54DSLim)      !14
    call initialize_config(TCC54S5Lim)      !15
    call initialize_config(TCC54DSS5Lim)    !16

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

    ETA42 = x3*x4p
    ETA52 = x3
    ETA41 = (one-ETA42) 
    ETA51 = (one - ETA52)
    eta54 = x3*x4**2/four/NF1(x3,x4p,lambda)

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
            sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,x4p,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: DS FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate ; weight is the same 
    !-------------------------------------------------------------------------------------------------------------


    ! --------- DOUBLE SOFT  and DOUBLE SOFT + S5--------------
    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    DSLim%PartFrac = (/xi1,xi2/)

    p1= half*sqrts*(/one,zero,zero,one/) 
    p2 =half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
    p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
         sin5*(sin(phi4)*cosphi45+cos(phi4)*sinphi45 ), cos5 /)

    pVV = p1 + p2
    call get_lept_mom(MV,pVV,yr(7:9),kallenF,DSLim)


    DSLim%AmpMom(:,1) = p1
    DSLim%AmpMom(:,2) = p2

    DSLim%LimMom(:,1) = p1
    DSLim%LimMom(:,2) = p2
    DSLim%LimMom(:,3) = p4 
    DSLim%LimMom(:,4) = p5

    DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2


    DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2 

    ! ------   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT  ------------

    ETA42 = zero
    ETA52 = zero
    ETA41 = one
    ETA51 = one
    eta54 = zero 

    S24 = spart *x1*x3*x4p 
    S25 = spart *x1*x2*x3
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

    ztot = - one 
    z4 = x1/ztot 
    z5 = x1*x2/ztot

    TCDSLim%KinInv  = (/S24,S25,s45,z4,z5/)


    TCDSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    ! ------   TRIPLE COLLINEAR KINEMATICS +  DOUBLE SOFT  + S5 ---------
    ztot = - one 
    z4 = x1/ztot
    z5 = x1*x2/ztot

    S24 = spart *x1*x3*x4p 
    S25 = spart *x1*x2*x3
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

    TCDSS5Lim%KinInv = (/S24,S25,s45,z4,z5/)

    TCDSS5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,x4p,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    !--------  C54 DS FLM  ----------

    ETA42 = x3
    ETA52 = x3
    ETA41 = (one - ETA42) 
    ETA51 = (one - ETA52)
    eta54 = zero 

    z5 = x2/(one+x2) 
    s12 = two*scr(p1,p2) 
    s24 = spart *x1*x3 
    s14 = spart * x1* (  one - x3)
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)

    kt1kt2sq = lambda !* (one-xx3)
    C54DSLim%KinInv = (/s12,s24,s14,s45,z5,kt1kt2sq/)

    C54DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,one,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


    ! ------  C54 +DS+S5 ---------
    z5 = x2
    s12 = two*scr(p1,p2) 
    S24 = spart *x1*x3 
    S14 = spart * x1* (  one - x3)
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)

    C54DSS5Lim%KinInv = (/s12,s45,S24,S14,z5/)
    C54DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(x3,one,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2


    ! ----   TRIPLE COLLINEAR KINEMATICS  + C54 + Double soft  ---------

    ETA42 = zero
    ETA52 = zero
    ETA41 = one
    ETA51 = one
    eta54 = zero 

    !--NEW
    s24 = spart *x1*x3 
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
    s2_45 = s24 * (one+x2)
    z5 = one/(one+x2)
    z2 = -one/(x1+x1*x2)
    kt1kt2sq = lambda * (one-xx3)

    TCC54DSLim%KinInv = (/s45,s2_45,z5,z2,kt1kt2sq/)


    TCC54DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,one,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    ! -------   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT5 +DS ----------

    z45 = one - x1
    z5 = zero 
    S24 = spart *x1*x3 
    S25 = spart *x1*x2*x3
    s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)

    TCC54DSS5Lim%KinInv = (/S24,S25, s45, z45, z5/)
    TCC54DSS5Lim%wgt = one/8.0_dp/pi * kallenF * (spart/four)**2 & 
         * x3*x4/NF1(xx3,one,lambda) & 
         * x1**3*x2 & 
         * weight4252(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 FLM   ! angles the same as in hard, do not regenerate; weight is the same 
    !-------------------------------------------------------------------------------------------------------------

    ETA42 = x3*x4p
    ETA52 = x3
    ETA41 = (one-ETA42) 
    ETA51 = (one - ETA52)
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
       p5 = half*sqrts*x1*x2*(/one, sin5*( cos(phi4)*cosphi45 - sin(phi4)*sinphi45), & 
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
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2 
    endif


    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C54 + S5 
    !-------------------------------------------------------------------------------------------------------------


    ETA42 = x3
    ETA52 = x3
    ETA41 = (one - ETA42) 
    ETA51 = (one - ETA52)
    eta54 = zero 
    phi4 = twopi*yr(4) 
    cos4 = -one + two*x3 
    sin4 = sqrt(one-cos4**2)

    spart = MV2/(one - x1)
    if (spart.gt.sh) C54S5Lim%flag = .true.

    sqrts = sqrt(spart)
    Deltan = one - x1*(one - eta41)
    Deltad = one - x1*(one-eta42)

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


       z5 = x2
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)
       C54S5Lim%KinInv = (/zero,s45,zero,zero,z5/)

       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C54S5Lim)
       C54S5Lim%AmpMom(:,1) = p1
       C54S5Lim%AmpMom(:,2) = p2
       C54S5Lim%AmpMom(:,3) = p4



       C54S5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2
    endif



    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    ETA42 = zero
    ETA52 = zero
    ETA41 = one
    ETA51 = one
    eta54 = zero 

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) TCLim%flag = .true.
    if (TCLim%flag .eqv. .false.) spart = MV2/Deltas
    if (TCLim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one,zero,zero,-one /)

       S24 = spart *x1*x3*x4p 
       S25 = spart *x1*x2*x3
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)
       ztot = x1*(one +x2) - one 
       z4 = x1/ztot 
       z5 = x1*x2/ztot

       TCLim%KinInv = (/S24,S25, s45, z4, z5/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCLim)


       TCLim%AmpMom(:,1) = p1
       TCLim%AmpMom(:,2) = p2 - p4 - p5


       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

       ! --------------  TC  + C54 ----------

       s24 = spart *x1*x3 
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
       s2_45 = s24 * (one+x2)
       zz1 = x2/(one+x2)
       zz2 = x1*(one+x2)/(x1*(one+x2)-one)
       kt1kt2sq = lambda * (one-xx3)

       TCC54Lim%KinInv = (/s45,s2_45,zz1,zz2,kt1kt2sq/)

       TCC54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1**3*x2 & 
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

       z4 = x1 
       z5 = x1*x2

       S24 = spart *x1*x3*x4p 
       S25 = spart *x1*x2*x3
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,x4p,lambda)

       TCS5Lim%KinInv = (/S24,S25,s45,z4,z5/)

       pVV = p1 + p2 - p4 
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,TCS5Lim)

       TCS5Lim%AmpMom(:,1) = p1
       TCS5Lim%AmpMom(:,2) = p2 - p4 

       TCS5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,x4p,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2



       ! --------   TRIPLE COLLINEAR KINEMATICS  + C54 + SOFT5 ------------

       S24 = spart *x1*x3 
       S25 = spart *x1*x2*x3
       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(xx3,one,lambda)
       z45 = one - x1
       z5 = zero 

       TCC54S5Lim%KinInv = (/S24,S25, s45, z45, zero/)

       TCC54S5Lim%wgt = one/8.0_dp/pi * kallenF * (spart/four)**2 & 
            * x3*x4/NF1(xx3,one,lambda) & 
            * x1**3*x2 & 
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

    ETA42 = x3
    ETA52 = x3
    ETA41 = (one - ETA42) 
    ETA51 = (one - ETA52)
    eta54 = zero 

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C54Lim%flag = .true.
    if (C54Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C54Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = x2*p4 

       ncoll = (/ zero,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)

       !-- note a sign change in some terms nperp relative to sector b
       nperp  = (/ zero,  -cos4*sqrt(lambda)*cos(phi4) - sin(phi4)*sqrt(one-lambda) , & 
            -cos4*sqrt(lambda)*sin(phi4) + cos(phi4)*sqrt(one-lambda), & 
            sin4*sqrt(lambda) /)


       s45 = spart *x1**2 * x2 * x3 * x4**2/four/NF1(x3,one,lambda)
       z5 = x2/(one + x2) 

       C54Lim%KinInv = (/zero,zero,s45,zero,z5/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C54Lim)

       C54Lim%AmpMom(:,1) = p1
       C54Lim%AmpMom(:,2) = p2
       C54Lim%AmpMom(:,3) = p4 + p5 
       C54Lim%LimMom(:,1) = nperp
       C54Lim%LimMom(:,2) = ncoll

       C54Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*x4/NF1(x3,one,lambda) & 
            * x1**3*x2 & 
            * weight4252(eta41,eta42,eta51,eta52,eta54) &
            *one/two/spart & 
            * spart/MV2

    endif


  end subroutine kinematics_nnlo_4252d


  ! todo: finish this
  subroutine kinematics_nnlo_4152(MV2,Y,yr,&
       HardProc,DSLim,S5Lim,DSS5Lim,C4Lim,C4DSLim,C4S5Lim,C4DSS5Lim,&
       C5Lim,C5DSLim,C5S5Lim,C5DSS5Lim,C4C5Lim,C4C5DSLim,C4C5S5Lim,C4C5DSS5Lim)
    real(dp), intent(in) ::  yr(9), MV2, Y
    type(KinConfig)      :: HardProc,DSLim,DSS5Lim,S5Lim,C4Lim,C4DSLim,C4DSS5Lim,C4S5Lim,&
         C5Lim,C5DSLim,C5DSS5Lim,C5S5Lim,C4C5Lim,C4C5DSLim,C4C5S5Lim,C4C5DSS5Lim
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos4, sin4, phi4, phi5, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4252, S24, S25, s45, z4, z5, ztot, weight4152
    real(dp) :: spart,sqrts, nperp(4), ncoll(4), S15, s12, S14, yy, Deltas
    real(dp) :: xx3

    !    statement functions
    weight4152(eta41,eta42,eta51,eta52,eta54)  = ETA42*ETA51*eta54/(eta54 + ETA52 + ETA41)

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(DSLim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(DSS5Lim)        !4
    call initialize_config(C4Lim)          !5
    call initialize_config(C4DSLim)        !6
    call initialize_config(C4S5Lim)        !7
    call initialize_config(C4DSS5Lim)      !8
    call initialize_config(C5Lim)          !9
    call initialize_config(C5DSLim)        !10
    call initialize_config(C5S5Lim)        !11
    call initialize_config(C5DSS5Lim)      !12
    call initialize_config(C4C5Lim)        !13
    call initialize_config(C4C5DSLim)      !14
    call initialize_config(C4C5S5Lim)      !15
    call initialize_config(C4C5DSS5Lim)    !16

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

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = one -  two*x3 
    sin4 = sqrt(one-cos4**2)

    phi5 = twopi*yr(6)
    cos5 = -one + two*x4
    sin5 = sqrt(one-cos5**2) 

    ETA42 = (one - x3)
    ETA52 = x4
    ETA41 = x3 
    ETA51 = (one - x4) 
    eta54 =   half*( one - sin4*sin5*cos(phi4-phi5) -cos4*cos5)


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) HardProc%flag = .true.
    if (HardProc%flag .eqv. .false.) spart = MV2/Deltas
    if (HardProc%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4 
       HardProc%AmpMom(:,4) = p5 


       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: DS FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! -------------- DOUBLE SOFT AND DOUBLE SOFT + S5 --------------------

    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    DSLim%PartFrac = (/xi1,xi2/)

    p1 = half*sqrts*(/one,zero,zero,one/) 
    p2 = half*sqrts*(/one,zero,zero,-one/) 
    p4 = half*sqrts*x1 *(/ one,  sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
    p5 = half*sqrts*x1*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)

    pVV = p1 + p2
    call get_lept_mom(MV,pVV,yr(7:9),kallenF,DSLim)

    DSLim%AmpMom(:,1) = p1
    DSLim%AmpMom(:,2) = p2

    DSLim%LimMom(:,1) = p1
    DSLim%LimMom(:,2) = p2
    DSLim%LimMom(:,3) = p4 
    DSLim%LimMom(:,4) = p5

    DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4152(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2

    DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4152(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart &
         * spart/MV2


    ! ------------   C42  +  DOUBLE SOFT   -----------

    ETA42 = one
    ETA52 = x4 
    ETA41 = zero  
    ETA51 = (one - x4) 
    eta54 =   ETA51

    S14 = spart *x1*x3 
    S15 = spart*x1*x2*(one-x4) 
    S25 = spart*x1*x2*x4
    s12 = two*scr(p1,p2) 
    z4 = x1 

    C4DSLim%KinInv = (/s12,S25,S15,S14,z4/)

    C4DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4152(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


    ! -------------   C42 +  DOUBLE SOFT  + S5 ---------------

    z4 = x1

    S14 = spart *x1*x3 
    S15 = spart*x1*x2*(one-x4) 
    S25 = spart*x1*x2*x4
    s12 = two*scr(p1,p2) 

    C4DSS5Lim%KinInv = (/s12,S25,S15,S14,z4/)

    C4DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4152(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2




    ! ------------  C51 DS FLM  ------------

    phi4 = twopi*yr(4) 
    cos4 = one - two*x3 
    sin4 = sqrt(one-cos4**2)

    ETA41 = x3
    ETA51 = one
    ETA42 = (one-x3) 
    ETA52 = zero
    eta54 = (one-x3) 

    s12 = two*scr(p1,p2) 
    S14 = spart*x1*x3
    S24 = spart*x1*(one - x3) 

    S25 = spart*x1*x2*(x4)
    z5 = x1*x2

    C5DSLim%KinInv = (/s12,S24,S14,S25,z5/)

    C5DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4152(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


    ! -----------  C51+DS+S5 -----------------

    s12 = two*scr(p1,p2) 
    S14 = spart*x1*x3
    S24 = spart*x1*(one - x3) 

    S25 = spart*x1*x2*(x4)
    z5 = x1*x2

    C5DSS5Lim%KinInv = (/s12,S24,S14,S25,z5/)

    C5DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4152(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


    ! -----------   C42  + C51 + Double soft  ----------


    ETA41 = zero
    ETA51 = one
    ETA42 = one
    ETA52 = zero
    eta54 = one

    S14 = spart *x1*x3 
    S25 = spart *x1*x2*x4
    z4 = x1 
    z5 = x1*x2

    C4C5DSLim%KinInv = (/S14,S25,zero,z4,z5/)

    C4C5DSLim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4152(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2


    ! ----------   C41  + C52 + SOFT5 +DS ------------

    z4 = x1 
    z5 = x1*x2
    S14 = spart *x1*x3 
    S25 = spart *x1*x2*x4

    C4C5DSS5Lim%KinInv = (/S14,S25,zero,z4,z5/)

    C4C5DSS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1**3*x2 & 
         * weight4152(eta41,eta42,eta51,eta52,eta54) &
         * one/two/spart & 
         * spart/MV2






    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3: S5 FLM   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------


    ! --------------- SOFT 5 ----------

    phi4 = twopi*yr(4) 
    cos4 = one -  two*x3 
    sin4 = sqrt(one-cos4**2)

    phi5 = twopi*yr(6)
    cos5 = -one + two*x4
    sin5 = sqrt(one-cos5**2) 

    ETA42 = (one - x3)
    ETA52 = x4
    ETA41 = x3 
    ETA51 = (one - x4) 
    eta54 =   half*( one - sin4*sin5*cos(phi4-phi5) -cos4*cos5)    
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
       p5 = half*sqrts*x1*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)


       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p4 

       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = p4 
       S5Lim%LimMom(:,4) = p5 

       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart &
            * spart/MV2


       ! ---------  C51 + S5 --------------


       phi4 = twopi*yr(4) 
       cos4 = one - two*x3 
       sin4 = sqrt(one-cos4**2)

       ETA41 = x3
       ETA51 = one
       ETA42 = (one-x3) 
       ETA52 = zero
       eta54 = (one-x3) 


       S25 = spart*x1*x2*x4
       z5 = x1*x2

       C5S5Lim%KinInv = (/zero,S25,zero,zero,z5/)

       C5S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif



    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3:   C41  sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------



    ETA42 = one
    ETA52 = x4 
    ETA41 = zero  
    ETA51 = (one - x4) 
    eta54 =   ETA51


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C4Lim%flag  = .true.
    if (C4Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4Lim%flag .eqv. .false.) sqrts = sqrt(spart)


    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p4 = half*sqrts*x1 *(/ one,  zero, zero, +one  /)
       p5 = half*sqrts*x1*x2*(/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)


       S14 = spart *x1*x3 
       z4 = one/(one - x1)

       C4Lim%KinInv = (/S14,zero,zero, z4, zero/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4Lim)

       C4Lim%AmpMom(:,1) = p1 - p4 
       C4Lim%AmpMom(:,2) = p2
       C4Lim%AmpMom(:,3) = p5

       nperp = (/zero,cos(phi4),sin(phi4),zero/)
       ncoll = (/zero,zero,zero,one/)
       
       C4Lim%LimMom(:,1) = nperp
       C4Lim%LimMom(:,2) = ncoll

       C4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
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
       p4 = half*sqrts*x1 *(/ one, zero, zero, +one  /)

       S14 = spart *x1*x3 
       S15 = spart*x1*x2*(one-x4) 
       S25 = spart*x1*x2*x4
       s12 = two*scr(p1,p2) 

       z4 = one - x1 

       C4S5Lim%KinInv = (/s12,S25,S15,S14,z4/)

       pVV = p1 + p2 - p4 
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4S5Lim)

       C4S5Lim%AmpMom(:,1) = p1 - p4 
       C4S5Lim%AmpMom(:,2) = p2


       C4S5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


       ! ---------   C41   + C52 + SOFT5 ------------


       ETA41 = zero
       ETA51 = one
       ETA42 = one
       ETA52 = zero
       eta54 = one           

       z4 = one/(one-x1) 
       z5 = x1*x2

       S14 = spart *x1*x3 
       S25 = spart *x1*x2*x4

       C4C5S5Lim%KinInv = (/S14,S25,zero,z4,z5/)

       C4C5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
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

    ETA41 = x3
    ETA51 = one
    ETA42 = (one-x3) 
    ETA52 = zero
    eta54 = (one-x3) 


    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C5Lim%flag = .true.
    if (C5Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C5Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p5 = half*sqrts*x1*x2*(/one, zero, zero, -one /)

       S25 = spart*x1*x2*x4
       z5 = one/(one - x1 * x2) 

       C5Lim%KinInv = (/zero,S25,zero,zero,z5/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C5Lim)

       C5Lim%AmpMom(:,1) = p1
       C5Lim%AmpMom(:,2) = p2 - p5  
       C5Lim%AmpMom(:,3) = p4
       
       nperp = (/zero,cos(phi5),sin(phi5),zero/)
       ncoll = (/zero,zero,zero,-one/)

       C5Lim%LimMom(:,1) = nperp
       C5Lim%LimMom(:,2) = ncoll

       C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2

    endif




    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C41  + C52
    !-------------------------------------------------------------------------------------------------------------

    ETA41 = zero
    ETA51 = one
    ETA42 = one
    ETA52 = zero
    eta54 = one  

    Deltas = one - x1*(one+x2) + x1**2*x2*eta54
    if (Deltas.lt.MV2/sh) C4C5Lim%flag = .true.
    if (C4C5Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4C5Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x1*x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x1*x2*(one - eta52)  

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
       p4 = half*sqrts*x1 *(/ one, zero, zero, +one  /)
       p5 = half*sqrts*x1*x2*(/one,zero,zero,-one /)

       S14 = spart *x1*x3 
       S25 = spart *x1*x2*x4
       z4 = one/(one - x1)
       z5 = one/(one - x1*x2) 

       C4C5Lim%KinInv = (/S14,S25,zero,z4,z5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4C5Lim)

       C4C5Lim%AmpMom(:,1) = p1 - p4 
       C4C5Lim%AmpMom(:,2) = p2 - p5

       C4C5Lim%LimMom(:,1) = (/zero,cos(phi4),sin(phi4),zero/) !-- nperp
       C4C5Lim%LimMom(:,2) = (/zero,zero,zero,one/) !-- ndir
       C4C5Lim%LimMom(:,3) = (/zero,cos(phi5),sin(phi5),zero/) !-- nperp
       C4C5Lim%LimMom(:,4) = (/zero,zero,zero,-one/) !-- ndir

       C4C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1**3*x2 & 
            * weight4152(eta41,eta42,eta51,eta52,eta54) &
            * one/two/spart & 
            * spart/MV2


    endif

  end subroutine kinematics_nnlo_4152


end module mod_kinematicsNNLO_EordAord


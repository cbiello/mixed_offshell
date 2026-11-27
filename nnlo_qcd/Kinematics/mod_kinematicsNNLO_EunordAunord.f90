module mod_kinematicsNNLO_EunordAunord
  use mod_types; use mod_consts_dp; use common_def
  use mod_process
  use mod_auxfunctions
  use mod_kinematics_gen; use mod_kinematics_lept

  implicit none
  private

  public :: kinematics_nnlo

contains

  subroutine kinematics_nnlo(MV2,Y,yr,&
       HardProc,&
       C41Lim,C51Lim,C42Lim,C52Lim,&
       C4152Lim,C4251Lim)
    !
    real(dp), intent(in) ::  MV2,Y,yr(9)
    !
    type(KinConfig) :: HardProc
    type(KinConfig) :: C41Lim,C51Lim,C42Lim,C52Lim
    type(KinConfig) :: C4152Lim,C4251Lim
    !--
    real(dp) :: eta41,eta51,eta42,eta52,eta54
    real(dp) :: MV
    real(dp) :: x1,x2,x3,x4
    real(dp) :: cos4,sin4,cos5,sin5,phi4,phi5,cphi4,sphi4,cphi5,sphi5
    real(dp) :: deltas,deltan,deltad
    real(dp) :: spart,sqrts,xi1,xi2,yy
    real(dp) :: pVV(4),p1(4),p2(4),p4(4),p5(4)
    real(dp) :: kallenF
    real(dp) :: z1,z2,z4,z5,s14,s15,s24,s25

    !-- initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(C41Lim)         !2
    call initialize_config(C51Lim)         !3
    call initialize_config(C42Lim)         !4
    call initialize_config(C52Lim)         !5
    call initialize_config(C4152Lim)       !6
    call initialize_config(C4251Lim)       !7

    !--

    MV =  sqrt(MV2)

    x1 =   yr(1)
    x2 =   yr(2) 
    x3 =   yr(3) 
    x4 =   yr(5) 

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = one - two*x3 
    sin4 = sqrt(one-cos4**2)
    cphi4 = cos(phi4)
    sphi4 = sin(phi4)

    phi5 = twopi*yr(6)
    cos5 = one - two*x4 
    sin5 = sqrt(one-cos5**2) 
    cphi5 = cos(phi5)
    sphi5 = sin(phi5)

    eta41 = x3
    eta51 = x4
    eta42 = one - x3 
    eta52 = one - x4
    eta54 = half*(one - sin4*sin5*cos(phi4-phi5) - cos4*cos5)

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

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1*(/one,sin4*cphi4,sin4*sphi4,cos4/)
       p5 = half*sqrts*x2*(/one,sin5*cphi5,sin5*sphi5,cos5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,HardProc)
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4
       HardProc%AmpMom(:,4) = p5

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: C41 FLM 
    !-------------------------------------------------------------------------------------------------------------

    eta41 = zero
    eta51 = x4
    eta42 = one
    eta52 = one - x4
    eta54 = eta51

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C41Lim%flag = .true.
    if (C41Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C41Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (C41Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C41Lim%Flag = .true.
    else
       C41Lim%flag = .true.
    endif

    if (C41Lim%flag .eqv. .false.) then 

       C41Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)  
       p4 = half*sqrts*x1*(/one,zero,zero,one/)
       p5 = half*sqrts*x2*(/one,sin5*cphi5,sin5*sphi5,cos5/)
       
       pVV = p1 + p2 - p4 - p5

       s14 = spart * x1 * x3

       z1 = one/(one-x1)

       C41Lim%LimMom(:,1) = (/zero,cphi4,sphi4,zero/)
       C41Lim%LimMom(:,2) = (/zero,zero,zero,one/)

       C41Lim%KinInv = (/s14,zero,zero,z1,zero/)
       
       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C41Lim)

       C41Lim%AmpMom(:,1) = p1-p4
       C41Lim%AmpMom(:,2) = p2
       C41Lim%AmpMom(:,3) = p5

       C41Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 &
            * one/two/spart &
            * spart/MV2

    endif
       

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: C51 FLM 
    !-------------------------------------------------------------------------------------------------------------

    eta41 = x3
    eta51 = zero
    eta42 = one - x3
    eta52 = one
    eta54 = eta41

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C51Lim%flag = .true.
    if (C51Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C51Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (C51Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C51Lim%Flag = .true.
    else
       C51Lim%flag = .true.
    endif

    if (C51Lim%flag .eqv. .false.) then 

       C51Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)  
       p4 = half*sqrts*x1*(/one,sin4*cphi4,sin4*sphi4,cos4/)
       p5 = half*sqrts*x2*(/one,zero,zero,one/)
       
       pVV = p1 + p2 - p4 - p5

       s15 = spart * x2 * x4

       z1 = one/(one-x2)

       C51Lim%LimMom(:,1) = (/zero,cphi5,sphi5,zero/)
       C51Lim%LimMom(:,2) = (/zero,zero,zero,one/)

       C51Lim%KinInv = (/s15,zero,zero,z1,zero/)
       
       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C51Lim)

       C51Lim%AmpMom(:,1) = p1-p5
       C51Lim%AmpMom(:,2) = p2
       C51Lim%AmpMom(:,3) = p4

       C51Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: C42 FLM 
    !-------------------------------------------------------------------------------------------------------------

    eta41 = one
    eta51 = x4
    eta42 = zero
    eta52 = one - x4
    eta54 = eta52

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C42Lim%flag = .true.
    if (C42Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C42Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (C42Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C42Lim%Flag = .true.
    else
       C42Lim%flag = .true.
    endif

    if (C42Lim%flag .eqv. .false.) then 

       C42Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)  
       p4 = half*sqrts*x1*(/one,zero,zero,-one/)
       p5 = half*sqrts*x2*(/one,sin5*cphi5,sin5*sphi5,cos5/)
       
       pVV = p1 + p2 - p4 - p5

       s24 = spart * x1 * (one-x3)

       z2 = one/(one-x1)

       C42Lim%LimMom(:,1) = (/zero,cphi4,sphi4,zero/)
       C42Lim%LimMom(:,2) = (/zero,zero,zero,-one/)

       C42Lim%KinInv = (/s24,zero,zero,z2,zero/)
       
       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C42Lim)

       C42Lim%AmpMom(:,1) = p1
       C42Lim%AmpMom(:,2) = p2-p4
       C42Lim%AmpMom(:,3) = p5

       C42Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 &
            * one/two/spart &
            * spart/MV2

    endif
       

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: C52 FLM 
    !-------------------------------------------------------------------------------------------------------------

    eta41 = x3
    eta51 = one
    eta42 = one - x3
    eta52 = zero
    eta54 = eta42

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C52Lim%flag = .true.
    if (C52Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C52Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (C52Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C52Lim%Flag = .true.
    else
       C52Lim%flag = .true.
    endif

    if (C52Lim%flag .eqv. .false.) then 

       C52Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)  
       p4 = half*sqrts*x1*(/one,sin4*cphi4,sin4*sphi4,cos4/)
       p5 = half*sqrts*x2*(/one,zero,zero,-one/)
       
       pVV = p1 + p2 - p4 - p5

       s25 = spart * x2 * (one-x4)
       z2 = one/(one - x2)

       C52Lim%LimMom(:,1) = (/zero,cphi5,sphi5,zero/)
       C52Lim%LimMom(:,2) = (/zero,zero,zero,-one/)

       C52Lim%KinInv = (/s25,zero,zero,z2,zero/)
       
       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C52Lim)

       C52Lim%AmpMom(:,1) = p1
       C52Lim%AmpMom(:,2) = p2-p5
       C52Lim%AmpMom(:,3) = p4

       C52Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C4152 FLM 
    !-------------------------------------------------------------------------------------------------------------

    eta41 = zero
    eta51 = one
    eta42 = one
    eta52 = zero
    eta54 = one

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C4152Lim%flag = .true.
    if (C4152Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4152Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (C4152Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C4152Lim%Flag = .true.
    else
       C4152Lim%flag = .true.
    endif

    if (C4152Lim%flag .eqv. .false.) then 

       C4152Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)  
       p4 = half*sqrts*x1*(/one,zero,zero,+one/)
       p5 = half*sqrts*x2*(/one,zero,zero,-one/)
       
       pVV = p1 + p2 - p4 - p5

       s14 = spart * x1 * x3
       s25 = spart * x2 * (one-x4)

       z1 = one/(one - x1)
       z2 = one/(one - x2)

       C4152Lim%LimMom(:,1) = (/zero,cphi4,sphi4,zero/) !-- nperp
       C4152Lim%LimMom(:,2) = (/zero,zero,zero,one/) !-- ndir
       C4152Lim%LimMom(:,3) = (/zero,cphi5,sphi5,zero/) !-- nperp
       C4152Lim%LimMom(:,4) = (/zero,zero,zero,-one/) !-- ndir

       C4152Lim%KinInv = (/s14,s25,zero,z1,z2/)
       
       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4152Lim)

       C4152Lim%AmpMom(:,1) = p1-p4
       C4152Lim%AmpMom(:,2) = p2-p5

       C4152Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C4251 FLM 
    !-------------------------------------------------------------------------------------------------------------

    eta41 = one
    eta51 = zero
    eta42 = zero
    eta52 = one
    eta54 = one

    Deltas = one - (x1+x2) + x1*x2*eta54
    if (Deltas.lt.MV2/sh) C4251Lim%flag = .true.
    if (C4251Lim%flag .eqv. .false.) spart = MV2/Deltas
    if (C4251Lim%flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)  

    if ( (C4251Lim%flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) C4251Lim%Flag = .true.
    else
       C4251Lim%flag = .true.
    endif

    if (C4251Lim%flag .eqv. .false.) then 

       C4251Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)  
       p4 = half*sqrts*x1*(/one,zero,zero,-one/)
       p5 = half*sqrts*x2*(/one,zero,zero,+one/)
       
       pVV = p1 + p2 - p4 - p5

       s24 = spart * x1 * (one-x3)
       s15 = spart * x2 * x4

       z2 = one/(one - x1)
       z1 = one/(one - x2)

       C4251Lim%LimMom(:,1) = (/zero,cphi4,sphi4,zero/) !-- nperp
       C4251Lim%LimMom(:,2) = (/zero,zero,zero,-one/) !-- ndir
       C4251Lim%LimMom(:,3) = (/zero,cphi5,sphi5,zero/) !-- nperp
       C4251Lim%LimMom(:,4) = (/zero,zero,zero,+one/) !-- ndir

       C4251Lim%KinInv = (/s24,s15,zero,z2,z1/)
       
       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),kallenF,C4251Lim)

       C4251Lim%AmpMom(:,1) = p1-p5
       C4251Lim%AmpMom(:,2) = p2-p4

       C4251Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 &
            * one/two/spart &
            * spart/MV2

    endif

  end subroutine kinematics_nnlo

end module mod_kinematicsNNLO_EunordAunord

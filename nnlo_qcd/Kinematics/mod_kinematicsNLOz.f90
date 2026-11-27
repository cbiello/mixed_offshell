module mod_kinematicsNLOz
  use mod_types; use mod_consts_dp; use common_def
  use mod_auxfunctions
  use mod_process
  use mod_kinematics_gen
  use mod_kinematics_lept
  implicit none
  private

  public :: kinematics_fullONLO_nloz
  public :: kinematics_c42only_nloz,kinematics_c41only_nloz

contains 

  subroutine checkandgetxis(MV2,Y,sh,z1,z2,x1,eta41,eta42,flag,spart,sqrts,xi1,xi2)
    real(dp), intent(in) :: MV2,Y,sh,z1,z2,x1,eta41,eta42
    logical, intent(out) :: flag 
    real(dp), intent(out) :: xi1,xi2,spart,sqrts 
    real(dp) :: Deltas,Deltan,Deltad,yy
    
    flag = .false.
    xi1 = two 
    xi2 = two 
    spart = zero 
    sqrts = zero 

    Deltas = z1*z2 - z1*x1*eta41 - z2*x1*eta42
    
    
    if (Deltas.lt.MV2/sh) flag = .true.
    if (flag .eqv. .false.) spart = MV2/Deltas
    if (flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = z1 - x1*eta42
    Deltad = z2 - x1*eta41
    if ( (flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) flag = .true.
    else
       flag = .true.
    endif
    return 

  end subroutine checkandgetxis


  !------------------------------------------------------------------------------------------------------


  !-- we factor out an extra 1/(8*pisq) from the phase space
  !-- which multiplied by gs**2 in the ME reconstructs asontwopi
  subroutine kinematics_fullONLO_nloz(MV2,Y,z1,z2,yr,HardProc,SoftLim,Col1Lim,Col2Lim,SCol1Lim,SCol2Lim)
    real(dp), intent(in) :: MV2,Y,yr(6),z1,z2
    type(KinConfig)   :: HardProc,SoftLim,Col1Lim,Col2Lim,SCol1Lim,SCol2Lim
    real(dp) :: x1,x2,x3,cos4,sin4,phi4,eta41,eta42,spart,sqrts,xi1,xi2,kallenF,MV
    real(dp) :: p1(4),p2(4),p4(4),pVV(4),p124(4),p4coll(4)
    real(dp) :: s14,s24,z41,z42

    ! -- initialization of xi-values 

    call initialize_config(HardProc)
    call initialize_config(SoftLim)
    call initialize_config(Col1Lim)
    call initialize_config(Col2Lim)
    call initialize_config(SCol1Lim)
    call initialize_config(SCol2Lim)

    MV = sqrt(MV2) 

    x1 =   yr(1) 
    x2 =   yr(2) 
    x3 =   yr(3) 
    cos4 = one - two*x2 
    sin4 = sqrt(one-cos4**2)
    phi4 = twopi*x3 

    !-------------------------------------------------------------------------------------------------------------
    !   I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    cos4 = one - two*x2 
    sin4 = sqrt(one-cos4**2)
    phi4 = twopi*x3 

    eta41 = x2
    eta42 = one - x2


    call checkandgetxis(MV2,Y,sh,z1,z2,x1,eta41,eta42,HardProc%flag,spart,sqrts,xi1,xi2)


    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac= (/xi1,xi2/)

       p1=  half*sqrts*z1*(/one,zero,zero,one/) 
       p2 = half*sqrts*z2*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

       pVV = p1 + p2 - p4

       call get_lept_mom(MV,pVV,yr(4:6),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4

       HardProc%LimMom(:,1) = p1
       HardProc%LimMom(:,2) = p2
       HardProc%LimMom(:,3) = p4

       ! use for spin correlated limits
       HardProc%LimMom(:,4) = (/zero,zero,zero,one/)*sin4-(/zero,cos(phi4),sin(phi4),zero/)*cos4 
       HardProc%LimMom(:,5) = (/zero,sin4*cos(phi4),sin4*sin(phi4),cos4/)

       HardProc%wgt =  one/8.0_dp/pi*kallenF * ( x1 ) * ( spart/two) &
            *one/(two*spart*z1*z2) & ! -- flux
            *spart/MV2 

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   S4 FL(1,2,4,5)  !  momenta are identical to the hard piece
    !-------------------------------------------------------------------------------------------------------------

    call checkandgetxis(MV2,Y,sh,z1,z2,zero,eta41,eta42,SoftLim%flag,spart,sqrts,xi1,xi2)

    if (SoftLim%flag .eqv. .false.) then 

       SoftLim%PartFrac = (/xi1,xi2/)

       p1=  half*sqrts*z1*(/one,zero,zero,one/) 
       p2 = half*sqrts*z2*(/one,zero,zero,-one/)
       p4 = half*sqrts*x1*(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

       pVV = p1 + p2 
       call get_lept_mom(MV,pVV,yr(4:6),kallenF,SoftLim)      

       SoftLim%AmpMom(:,1) = p1
       SoftLim%AmpMom(:,2) = p2
       SoftLim%LimMom(:,1) = p1
       SoftLim%LimMom(:,2) = p2
       SoftLim%LimMom(:,3) = p4

       SoftLim%wgt =  one/8.0_dp/pi*kallenF * ( x1 ) * ( spart/two) &
            *one/(two*spart*z1*z2) & ! -- flux
            *spart/MV2 

    endif

    !----------------------------------------------------
    ! C41 * FLM 
    !-----------------------------------------------------

    eta41 = zero
    eta42 = one

    call checkandgetxis(MV2,Y,sh,z1,z2,x1,eta41,eta42,Col1Lim%flag,spart,sqrts,xi1,xi2)

    if (Col1Lim%flag .eqv. .false.) then 

       Col1Lim%PartFrac = (/xi1,xi2/)

       p1=  half*sqrts*z1*(/one,zero,zero,one/) 
       p2 = half*sqrts*z2*(/one,zero,zero,-one/)
       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

       p4coll = half*sqrts*x1 *(/ one,zero,zero,one /)
       pVV = p1 + p2 - p4coll
       call get_lept_mom(MV,pVV,yr(4:6),kallenF,Col1Lim)          

       s14 = spart*z1*x1*x2
       z41 = one - x1/z1

       Col1Lim%KinInv = (/z41,s14/) 

       Col1Lim%AmpMom(:,1) = p1 - p4coll 
       Col1Lim%AmpMom(:,2) = p2

       Col1Lim%LimMom(:,1) = p1 - p4coll 
       Col1Lim%LimMom(:,2) = p2
       Col1Lim%LimMom(:,3) = p4 


      Col1Lim%wgt =  one/8.0_dp/pi*kallenF * ( x1 ) * ( spart/two) &
            *one/(two*spart*z1*z2) &  ! -- flux
            *spart/MV2 


    endif

    !--------------------------------
    ! C42 FLM 
    !--------------------------------

    eta41 = one
    eta42 = zero

    call checkandgetxis(MV2,Y,sh,z1,z2,x1,eta41,eta42,Col2Lim%flag,spart,sqrts,xi1,xi2)
    if (Col2Lim%flag .eqv. .false.) then 

       Col2Lim%PartFrac = (/xi1,xi2/)

       p1=  half*sqrts*z1*(/one,zero,zero,one/) 
       p2 = half*sqrts*z2*(/one,zero,zero,-one/)
       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

       s24 = spart*z2*x1*(one-x2)
       z42 = one - x1/z2

       Col2Lim%KinInv = (/z42,s24/) 
       p4coll = half*sqrts*x1 *(/ one,zero,zero,-one /)

       pVV = p1 + p2 - p4coll
       call get_lept_mom(MV,pVV,yr(4:6),kallenF,Col2Lim)          

       Col2Lim%AmpMom(:,1) = p1 
       Col2Lim%AmpMom(:,2) = p2-p4coll

       Col2Lim%LimMom(:,1) = p1 
       Col2Lim%LimMom(:,2) = p2-p4coll
       Col2Lim%LimMom(:,3) = p4 


       Col2Lim%wgt =  one/8.0_dp/pi*kallenF * ( x1 ) * ( spart/two) &
            *one/(two*spart*z1*z2) & ! -- flux
            *spart/MV2 

         endif


    !----------------------------------------------------
    ! S4 * C41 * FLM 
    !----------------------------------------------------

    eta41 = zero
    eta42 = one

    call checkandgetxis(MV2,Y,sh,z1,z2,zero,eta41,eta42,SCol1Lim%flag,spart,sqrts,xi1,xi2)

    if (SCol1Lim%flag .eqv. .false.) then 

       SCol1Lim%PartFrac = (/xi1,xi2/)

       p1=  half*sqrts*z1*(/one,zero,zero,one/) 
       p2 = half*sqrts*z2*(/one,zero,zero,-one/)
       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

       pVV = p1 + p2 
       call get_lept_mom(MV,pVV,yr(4:6),kallenF,SCol1Lim)          
      
       SCol1Lim%AmpMom(:,1) = p1
       SCol1Lim%AmpMom(:,2) = p2

       SCol1Lim%LimMom(:,1) = p1
       SCol1Lim%LimMom(:,2) = p2
       SCol1Lim%LimMom(:,3) = p4

       s14 = spart*x1*z1*x2
       SCol1Lim%KinInv(1) = s14

       SCol1Lim%wgt =  one/8.0_dp/pi*kallenF * ( x1 ) * ( spart/two) &
            *one/(two*spart*z1*z2) & ! -- flux
           *spart/MV2 

    endif

    !--------------------------------
    ! S4 *C42 FLM 
    !--------------------------------

    eta41 = one
    eta42 = zero 

    call checkandgetxis(MV2,Y,sh,z1,z2,zero,eta41,eta42,SCol2Lim%flag,spart,sqrts,xi1,xi2)

    if (SCol2Lim%flag .eqv. .false.) then 

       SCol2Lim%PartFrac = (/xi1,xi2/)

       p1=  half*sqrts*z1*(/one,zero,zero,one/) 
       p2 = half*sqrts*z2*(/one,zero,zero,-one/)
       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

       pVV = p1 + p2 
       call get_lept_mom(MV,pVV,yr(4:6),kallenF,SCol2Lim)          

       SCol2Lim%AmpMom(:,1) = p1 
       SCol2Lim%AmpMom(:,2) = p2

       SCol2Lim%LimMom(:,1) = p1
       SCol2Lim%LimMom(:,2) = p2
       SCol2Lim%LimMom(:,3) = p4

       s24 = spart*x1*z2*(one-x2)
       SCol2Lim%KinInv(1) = s24

       SCol2Lim%wgt =  one/8.0_dp/pi*kallenF * ( x1 ) * ( spart/two) &
            *one/(two*spart*z1*z2) & ! -- flux
            *spart/MV2 

    endif

  end subroutine kinematics_fullONLO_nloz

  !-- identical to above but only has C42 limits
  subroutine kinematics_c42only_nloz(MV2,Y,z1,z2,yr,HardProc,Col2Lim)
    real(dp), intent(in) :: MV2,Y,yr(6),z1,z2
    type(KinConfig)   :: HardProc,Col2Lim
    real(dp) :: sqrts,spart
    real(dp) :: p1(4),p2(4),p4(4),pVV(4),p4coll(4),x1,x2,x3,kallenF
    real(dp) :: MV,cos4,sin4,phi4,eta41,eta42,xi1,xi2,s24,z42

    ! -- initialization of xi-values 

    call initialize_config(HardProc)
    call initialize_config(Col2Lim)

    !   MV2 and Y are generated 
    MV = sqrt(MV2) 

    x1 =   yr(1) 
    x2 =   yr(2) 
    x3 =   yr(3) 

    cos4 = one - two*x2 
    sin4 = sqrt(one-cos4**2)
    phi4 = twopi*x3 

    !-------------------------------------------------------------------------------------------------------------
    !   I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    cos4 = one - two*x2 
    sin4 = sqrt(one-cos4**2)
    phi4 = twopi*x3 


    eta41 = x2
    eta42 = one - x2


   call checkandgetxis(MV2,Y,sh,z1,z2,x1,eta41,eta42,HardProc%flag,spart,sqrts,xi1,xi2)


    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac= (/xi1,xi2/)

       p1=  half*sqrts*z1*(/one,zero,zero,one/) 
       p2 = half*sqrts*z2*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)
       pVV = p1 + p2 - p4

       call get_lept_mom(MV,pVV,yr(4:6),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4


       HardProc%LimMom(:,1) = p1
       HardProc%LimMom(:,2) = p2
       HardProc%LimMom(:,3) = p4

       ! use for spin correlated limits
       HardProc%LimMom(:,4) = (/zero,zero,zero,one/)*sin4-(/zero,cos(phi4),sin(phi4),zero/)*cos4 
       HardProc%LimMom(:,5) = (/zero,sin4*cos(phi4),sin4*sin(phi4),cos4/)

       HardProc%wgt =  one/8.0_dp/pi*kallenF * ( x1 ) * ( spart/two) &
            *one/(two*spart*z1*z2) & ! -- flux
            *spart/MV2 

    endif

    !--------------------------------
    ! C42 FLM 
    !--------------------------------

    eta41 = one
    eta42 = zero 

    call checkandgetxis(MV2,Y,sh,z1,z2,x1,eta41,eta42,Col2Lim%flag,spart,sqrts,xi1,xi2)

    if (Col2Lim%flag .eqv. .false.) then 

       Col2Lim%PartFrac = (/xi1,xi2/)

       p1=  half*sqrts*z1*(/one,zero,zero,one/) 
       p2 = half*sqrts*z2*(/one,zero,zero,-one/) 
       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

       s24 = two*scr(p2,p4) 
       z42 = p2(1)/(p2(1) - p4(1) )

       p4coll = half*sqrts*x1 *(/ one,zero,zero,-one /)

       pVV = p1 + p2 - p4coll
       call get_lept_mom(MV,pVV,yr(4:6),kallenF,Col2Lim)          

       Col2Lim%AmpMom(:,1) = p1 
       Col2Lim%AmpMom(:,2) = p2-p4coll

       Col2Lim%LimMom(:,1) = p1
       Col2Lim%LimMom(:,2) = p2-p4coll
       Col2Lim%LimMom(:,3) = p4 

       s24 = spart*x1*z2*(one-x2)
       z42 = one-x1/z2

       Col2Lim%KinInv = (/z42,s24/) 

       Col2Lim%wgt =  one/8.0_dp/pi*kallenF * ( x1 ) * ( spart/two) &
            *one/(two*spart*z1*z2) & ! -- flux
            *spart/MV2 

    endif


  end subroutine kinematics_c42only_nloz

  !-- identical to qqb but only has C41 limits
  subroutine kinematics_c41only_nloz(MV2,Y,z1,z2,yr,HardProc,Col1Lim)
    real(dp), intent(in) :: MV2,Y,yr(6),z1,z2
    type(KinConfig)   :: HardProc,Col1Lim
    real(dp) :: sqrts,spart
    real(dp) :: p1(4),p2(4),p4(4),pVV(4),p4coll(4),x1,x2,x3,kallenF
    real(dp) :: MV,cos4,sin4,phi4,eta41,eta42,xi1,xi2,s14,z41

    ! -- initialization of xi-values 

    call initialize_config(HardProc)
    call initialize_config(Col1Lim)

    !   MV2 and Y are generated 
    MV = sqrt(MV2) 

    x1 =   yr(1) 
    x2 =   yr(2) 
    x3 =   yr(3) 

    cos4 = one - two*x2 
    sin4 = sqrt(one-cos4**2)
    phi4 = twopi*x3 

    !-------------------------------------------------------------------------------------------------------------
    !   I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    cos4 = one - two*x2 
    sin4 = sqrt(one-cos4**2)
    phi4 = twopi*x3 


    eta41 = x2
    eta42 = one - x2


    call checkandgetxis(MV2,Y,sh,z1,z2,x1,eta41,eta42,HardProc%flag,spart,sqrts,xi1,xi2)


    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac= (/xi1,xi2/)

       p1=  half*sqrts*z1*(/one,zero,zero,one/) 
       p2 = half*sqrts*z2*(/one,zero,zero,-one/)
       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)
       pVV = p1 + p2 - p4

       call get_lept_mom(MV,pVV,yr(4:6),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p4


       HardProc%LimMom(:,1) = p1
       HardProc%LimMom(:,2) = p2
       HardProc%LimMom(:,3) = p4

       ! use for spin correlated limits
       HardProc%LimMom(:,4) = (/zero,zero,zero,one/)*sin4-(/zero,cos(phi4),sin(phi4),zero/)*cos4 
       HardProc%LimMom(:,5) = (/zero,sin4*cos(phi4),sin4*sin(phi4),cos4/)

       HardProc%wgt =  one/8.0_dp/pi*kallenF * ( x1 ) * ( spart/two) &
            *one/(two*spart*z1*z2) & ! -- flux
            *spart/MV2 

    endif


    !----------------------------------------------------
    ! C41 * FLM 
    !-----------------------------------------------------

    eta41 = zero	
    eta42 = one	

    call checkandgetxis(MV2,Y,sh,z1,z2,x1,eta41,eta42,Col1Lim%flag,spart,sqrts,xi1,xi2)

    if (Col1Lim%flag .eqv. .false.) then 

       Col1Lim%PartFrac = (/xi1,xi2/)

       p1=  half*sqrts*z1*(/one,zero,zero,one/) 
       p2 = half*sqrts*z2*(/one,zero,zero,-one/)
       p4 = half*sqrts*x1 *(/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

       p4coll = half*sqrts*x1 *(/ one,zero,zero,one /)
       pVV = p1 + p2 - p4coll
       call get_lept_mom(MV,pVV,yr(4:6),kallenF,Col1Lim)          

       Col1Lim%AmpMom(:,1) = p1 - p4coll
       Col1Lim%AmpMom(:,2) = p2

       Col1Lim%LimMom(:,1) = p1 - p4coll 
       Col1Lim%LimMom(:,2) = p2
       Col1Lim%LimMom(:,3) = p4

       s14 = spart*x1*z1*x2
       z41 = one-x1/z1

       Col1Lim%KinInv = (/z41,s14/) 

       Col1Lim%wgt =  one/8.0_dp/pi*kallenF * ( x1 ) * ( spart/two) &
            *one/(two*spart*z1*z2) &  ! -- flux
            *spart/MV2 

    endif

  end subroutine kinematics_c41only_nloz



end module mod_kinematicsNLOz


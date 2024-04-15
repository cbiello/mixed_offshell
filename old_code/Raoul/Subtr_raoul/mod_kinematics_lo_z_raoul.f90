module mod_kinematics_lo_z_raoul
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_process
  use mod_auxfunctions
  use mod_kinematics_lo
  use mod_kinematics_gen
  use mod_kinematics_lept
  implicit none

!!  integer, parameter :: LOdim_vegas = tauy_max + 2
!!  !
!!  integer, parameter :: kLO_min = tauy_max + 1
!!  integer, parameter :: kLO_max = tauy_max + 2
!!  integer, parameter :: kLO_max_full = kLO_max + 1
  
  public :: kinematics_lo_zzb

contains 
  

  subroutine kinematics_lo_zzb(yr,z1,z2,LOConfig)
    real(dp), intent(in) :: yr(kLO_max_full),z1,z2
    type(KinConfig)   :: LOConfig
    real(dp) :: tau,ylab,jac,kallenf
    real(dp) :: mv,mv2,spart,sqrts,xi1,xi2
    real(dp) :: s13,s23,s14,s24,s34,c,y,ometa34
    real(dp) :: p1(4),p2(4),p12(4),nlept(3,2)


    call initialize_config(LOConfig)
    
    call get_tauy(yr(1:tauy_max),tau,ylab,jac)

    mv2 = tau*sh
    mv = sqrt(mv2)

    spart = mv2/z1/z2
    sqrts = sqrt(spart)
    xi1 = sqrt(spart/sh)*sqrt(z2/z1)*exp(+ylab)
    xi2 = sqrt(spart/sh)*sqrt(z1/z2)*exp(-ylab)

    if (xi1 .gt. one .or. xi2 .gt. one) LOConfig%flag = .true.
    
    LOconfig%PartFrac(1:2) = (/xi1,xi2/)

    p1 = sqrts/two*z1*(/one,zero,zero, one/) 
    p2 = sqrts/two*z2*(/one,zero,zero,-one/) 
    p12(:) = p1(:)+p2(:)
    
    call get_lept_mom(MV,p12,yr(kLO_min:kLO_max_full),nlept,kallenF,LOConfig)
    LOConfig%AmpMom(:,1) = p1
    LOConfig%AmpMom(:,2) = p2

    !1/(two*mv2) is flux factor
    LOConfig%wgt  = one/8.0_dp/pi*kallenF/(two*mv2)*jac/z1/z2

    LOConfig%npart  = 4 
    
    LOConfig%mu2ref = mv2 !-- for hoppet

    ! From Federico
    c = sqrt(dot_product(nlept(:,1)+nlept(:,2),nlept(:,1)+nlept(:,2)))
    y = c**2/((one + (one + c))*((one - c) + one))


    LOConfig%Lim_Ei(1) = LOConfig%AmpMom(1,1)/z1
    LOConfig%Lim_Ei(2) = LOConfig%AmpMom(1,2)/z2
    LOConfig%Lim_Ei(3) = LOConfig%AmpMom(1,3)
    LOConfig%Lim_Ei(4) = LOConfig%AmpMom(1,4)

    s13 = two*scr( LOConfig%AmpMom(1:4,1),LOConfig%AmpMom(1:4,3) )/z1
    s23 = two*scr( LOConfig%AmpMom(1:4,2),LOConfig%AmpMom(1:4,3) )/z2
    s14 = two*scr( LOConfig%AmpMom(1:4,1),LOConfig%AmpMom(1:4,4) )/z1
    s24 = two*scr( LOConfig%AmpMom(1:4,2),LOConfig%AmpMom(1:4,4) )/z2
    s34 = two*scr( LOConfig%AmpMom(1:4,3),LOConfig%AmpMom(1:4,4) )

    LOConfig%Lim_etaij(3,1) = s13/four/LOConfig%Lim_Ei(1)/LOConfig%Lim_Ei(3)
    LOConfig%Lim_etaij(3,2) = s23/four/LOConfig%Lim_Ei(2)/LOConfig%Lim_Ei(3)
    LOConfig%Lim_etaij(4,1) = s14/four/LOConfig%Lim_Ei(1)/LOConfig%Lim_Ei(4)
    LOConfig%Lim_etaij(4,2) = s24/four/LOConfig%Lim_Ei(2)/LOConfig%Lim_Ei(4)
!    LOConfig%Lim_etaij(3,4) = s34/four/LOConfig%Lim_Ei(3)/LOConfig%Lim_Ei(4)
    LOConfig%Lim_etaij(3,4) = one/(one+y)
    ometa34 = y/(one+y)   ! 1-eta34 = y/(1+y) -- this is for numerical stability
    
!!    print *, "z1,z2",z1,z2
!!    print *, "s12,s34",two*scr(LOConfig%AmpMom(1:4,1),LOConfig%AmpMom(1:4,2)),s34
!!    print *, "E3,E4",LOConfig%Lim_Ei(3),LOConfig%Lim_Ei(4)
!!    print *, "p1",LOConfig%AmpMom(1:4,1)
!!    print *, "p2",LOConfig%AmpMom(1:4,2)
!!    print *, "p3",LOConfig%AmpMom(1:4,3)
!!    print *, "p4",LOConfig%AmpMom(1:4,4)
!!    print *, "mom  cons", LOConfig%AmpMom(1:4,1)+LOConfig%AmpMom(1:4,2)-(LOConfig%AmpMom(1:4,3)+LOConfig%AmpMom(1:4,4))
!!    print *, scr(LOConfig%AmpMom(1:4,3),LOConfig%AmpMom(1:4,3)),scr(LOConfig%AmpMom(1:4,4),LOConfig%AmpMom(1:4,4))
!!    pause


    LOConfig%Lim_KinInv = [s13,s23,s14,s24,s34,ometa34,zero]   !--this is a horrible place to put ometa34, but fine...
    
  end subroutine kinematics_lo_zzb
  
end module mod_kinematics_lo_z_raoul





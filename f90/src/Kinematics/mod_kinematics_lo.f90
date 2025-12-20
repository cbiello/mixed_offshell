module mod_kinematics_lo
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_process
  use mod_auxfunctions
  use mod_kinematics_gen
  use mod_kinematics_lept
  implicit none
  integer, public, parameter :: LOdim_vegas = tauy_max + 2
  !
  integer, public, parameter :: kLO_min = tauy_max + 1
  integer, public, parameter :: kLO_max = tauy_max + 2
  integer, public, parameter :: kLO_max_full = kLO_max + 1
  
  public :: kinematics_lo

contains 
  
  subroutine kinematics_lo(yr,LOConfig)
    real(dp), intent(in) :: yr(kLO_max_full)
    type(KinConfig)   :: LOConfig
    real(dp) :: tau,ylab,jac,kallenf
    real(dp) :: mv,mv2
    real(dp) :: p1(4),p2(4),p12(4),nlept(3,2)

    call initialize_config(LOConfig)
    
    call get_tauy(yr(1:tauy_max),tau,ylab,jac)
    LOconfig%PartFrac = sqrt(tau)*[exp(+ylab),exp(-ylab)]

    mv2 = tau*sh
    mv = sqrt(mv2)

    p1 = MV/two*(/one,zero,zero, one/) 
    p2 = MV/two*(/one,zero,zero,-one/) 
    p12(:) = p1(:)+p2(:)
    
    call get_lept_mom(MV,p12,yr(kLO_min:kLO_max_full),nlept,kallenF,LOConfig)
    LOConfig%AmpMom(:,1) = p1
    LOConfig%AmpMom(:,2) = p2

    !1/(two*mv2) is flux factor
    LOConfig%wgt  = one/8.0_dp/pi*kallenF/(two*mv2)*jac

    LOConfig%npart  = 4 
    
    allocate(LOConfig%part(LOConfig%npart))

    LOConfig%mu2ref = mv2 !-- for hoppet

    LOConfig%Lim_etaij(1,2) = one
    LOConfig%Lim_etaij(1,3) = half*(one-nlept(3,1))
    LOConfig%Lim_etaij(1,4) = half*(one-nlept(3,2))

    LOConfig%Lim_Ei(1) = MV/two

  end subroutine kinematics_lo
  
end module mod_kinematics_lo





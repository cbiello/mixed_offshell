module mod_kinematicsLO
  use mod_types; use mod_consts_dp; use common_def
  use mod_process
  use mod_auxfunctions
  use mod_kinematics_gen
  use mod_kinematics_lept
  implicit none
  private
  
  public :: kinematics_lo

contains 
  
  subroutine kinematics_lo(xi1,xi2,yr,LOConfig) 
    real(dp), intent(in) :: xi1,xi2,yr(3)
    type(KinConfig)   :: LOConfig
    real(dp) :: p1(4),p2(4)
    real(dp) :: m6,m7,mv,mv2,p12(4),kallenF

    LOConfig%flag = .false.

    MV =  sqrt(sh*xi1*xi2)

    MV2 = MV**2

    p1= MV/two*(/one,zero,zero,one/) 
    p2 =MV/two*(/one,zero,zero,-one/) 

    p12(:) = p1(:)+p2(:)
    
    call get_lept_mom(MV,p12,yr(1:3),kallenF,LOConfig)
    LOConfig%AmpMom(:,1) = p1
    LOConfig%AmpMom(:,2) = p2

    LOConfig%LimMom(:,1) = p1
    LOConfig%LimMom(:,2) = p2

    !1/(two*MV) is flux factor
    LOConfig%wgt  = one/8.0_dp/pi*kallenF/(two*MV2)

  end subroutine kinematics_lo

end module mod_kinematicsLO



module mod_kinematics_gen
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  implicit none
  private

  public :: get_tauy,get_tauyxi1xi2

contains

  subroutine get_tauy(xx,tau,y_lab,jac) !-- tau = Minv^2/s_had
    real(dp), intent(in) :: xx(1:tauy_max)
    real(dp), intent(out) :: tau,y_lab,jac
    
#if (_PhaseSpaceType == 1) 
    call get_tauy_flat(xx,tau,y_lab,jac)      !-- flat
#elif (_PhaseSpaceType == 2) 
    call get_tauy_BW(xx,mz,gaz,tau,y_lab,jac) !-- BW, z mass
#elif (_PhaseSpaceType == 3)
    call get_tauy_OS(xx,mh,tau,y_lab,jac)     !-- on-shell, h mass
#elif (_PhaseSpaceType == 4)
    call get_tauy_OS_dec(xx,mz,gaz,tau,y_lab,jac)     !-- on-shell, z mass
#elif (_PhaseSpaceType == 5)
    call get_tauy_OS_dec(xx,mz,gaz,tau,y_lab,jac)     !-- on-shell, z mass    
#elif (_PhaseSpaceType == 6)
    call get_tauy_OS_dec(xx,mw,gaw,tau,y_lab,jac)     !-- on-shell, w mass   
#endif
    return
  end subroutine get_tauy

  subroutine get_tauyxi1xi2(xx,z,zbar,tau,y_lab,xi1,xi2,jac) !-- tau = Minv^2/s_had
    real(dp), intent(in) :: xx(1:tauy_max),z,zbar
    real(dp), intent(out) :: tau,y_lab,jac,xi1,xi2
#if (_PhaseSpaceType == 1) 
    call get_tauyxi1xi2_flat(xx,z,zbar,tau,y_lab,xi1,xi2,jac)       !-- flat
#elif (_PhaseSpaceType == 2) 
    call get_tauyxi1xi2_BW(xx,z,zbar,mz,gaz,tau,y_lab,xi1,xi2,jac)  !-- BW, z mass
#elif (_PhaseSpaceType == 3)
    call get_tauyxi1xi2_OS(xx,z,zbar,mh,tau,y_lab,xi1,xi2,jac)      !-- on-shell, h mass
#elif (_PhaseSpaceType == 4)
    call get_tauyxi1xi2_OS_dec(xx,z,zbar,mz,gaz,tau,y_lab,xi1,xi2,jac)      !-- on-shell, z mass, with decay
#elif (_PhaseSpaceType == 5)
    call get_tauyxi1xi2_OS_dec(xx,z,zbar,mz,gaz,tau,y_lab,xi1,xi2,jac)      !-- on-shell, z mass, with decay 
#elif (_PhaseSpaceType == 6)
    call get_tauyxi1xi2_OS_dec(xx,z,zbar,mw,gaw,tau,y_lab,xi1,xi2,jac)      !-- on-shell, w mass, with decay 
#endif
    return
  end subroutine get_tauyxi1xi2

  !--------------------------------------------------------------------------------

  !-- here myq2min is not the minimum generated value
  subroutine get_tauy_flat(xx,tau,y_lab,jac) !-- tau = Minv^2/s_had
    real(dp), intent(in) :: xx(2)
    real(dp), intent(out) :: tau,y_lab,jac
    real(dp) :: lt,xx1tmp

    tau = myq2min*(one-xx(1))/(myq2min+sh*xx(1))
    lt = log(tau)

    !-- 1/2 Log[tau] < y_lab < 1/2 Log[1/tau]
    y_lab = lt * (half-xx(2))

    jac = -lt
    jac = jac*(tau*sh+myq2min)**2/(myq2min*(sh+myq2min))

    return

  end subroutine get_tauy_flat

  subroutine get_tauy_BW(xx,M,Ga,tau,y_lab,jac) !-- tau = Minv^2/s_had
    real(dp), intent(in) :: xx(2),M,Ga
    real(dp), intent(out) :: tau,y_lab,jac
    real(dp)              :: Msq,ximin,ximax,xi,sqrttau,ymax
    
    Msq=M**2
    ximin = +atan((myq2min-Msq)/(M*Ga))
    ximax = +atan((sh-Msq)/(M*Ga))
    xi = Msq + Ga*M*tan(ximin + (ximax-ximin)*xx(1))
    tau = xi/sh
    sqrttau = sqrt(tau)
    jac = (ximax-ximin)/(M*Ga*sh)*((xi-Msq)**2+Msq*Ga**2)

    ymax = -half * log(tau)
    y_lab = ymax * (two * xx(2) - one)
    jac = jac * two * ymax

    return

  end subroutine get_tauy_BW

  subroutine get_tauy_OS(xx,M,tau,y_lab,jac) !-- tau = Minv^2/s_had
    real(dp), intent(in) :: xx(1),M
    real(dp), intent(out) :: tau,y_lab,jac
    real(dp)              :: Msq,ximin,ximax,xi,sqrttau,ymax
    
    msq=M**2
    tau = msq/sh

    ymax = -half * log(tau)
    y_lab = ymax * (two * xx(1) - one)
    jac = two * ymax 

    return

  end subroutine get_tauy_OS

  subroutine get_tauy_OS_dec(xx,M,Ga,tau,y_lab,jac) !-- tau = Minv^2/s_had
    real(dp), intent(in) :: xx(1),M,Ga
    real(dp), intent(out) :: tau,y_lab,jac
    real(dp)              :: Msq,ximin,ximax,xi,sqrttau,ymax
    
    msq=M**2
    tau = msq/sh

    ymax = -half * log(tau)
    y_lab = ymax * (two * xx(1) - one)
    jac = two * ymax 

    jac = jac * M * Ga/sh * Pi

    return

  end subroutine get_tauy_OS_dec

  !--------------------------------------------------------------------------------

  subroutine get_tauyxi1xi2_flat(xx,z,zbar,tau,y_lab,xi1,xi2,jac) !-- tau = Minv^2/s_had
    real(dp), intent(in) :: xx(2),z,zbar
    real(dp), intent(out) :: tau,y_lab,jac,xi1,xi2
    real(dp) :: lt

    tau = myq2min*(one-xx(1))/(myq2min+sh*xx(1))
    lt = log(tau)

    !-- 1/2 Log[tau] < y_lab < 1/2 Log[1/tau]
    y_lab = lt * (half-xx(2))

    jac = -lt
    jac = jac*(tau*sh+myq2min)**2/(myq2min*(sh+myq2min))

    ! now get xi1 & xi2 at the same time
    xi1 = sqrt(tau)*exp(+y_lab)/z
    xi2 = sqrt(tau)*exp(-y_lab)/zbar

    jac = jac/z/zbar

    return

  end subroutine get_tauyxi1xi2_flat

  subroutine get_tauyxi1xi2_BW(xx,z,zbar,M,Ga,tau,y_lab,xi1,xi2,jac) !-- tau = Minv^2/s_had
    real(dp), intent(in) :: xx(2),M,Ga,z,zbar
    real(dp), intent(out) :: tau,y_lab,jac,xi1,xi2
    real(dp)              :: Msq,ximin,ximax,xi,sqrttau,ymax
    
    Msq=M**2
    ximin = +atan((-Msq)/(M*Ga))
    ximax = +atan((sh-Msq)/(M*Ga))
    xi = Msq + Ga*M*tan(ximin + (ximax-ximin)*xx(1))
    tau = xi/sh
    sqrttau = sqrt(tau)
    jac = (ximax-ximin)/(M*Ga*sh)*((xi-Msq)**2+Msq*Ga**2)

    ymax = -half * log(tau)
    y_lab = ymax * (two * xx(2) - one)
    jac = jac * two * ymax


    ! now get xi1 & xi2 at the same time -- as before
    xi1 = sqrt(tau)*exp(+y_lab)/z
    xi2 = sqrt(tau)*exp(-y_lab)/zbar

    jac = jac/z/zbar

    return

  end subroutine get_tauyxi1xi2_BW

  subroutine get_tauyxi1xi2_OS(xx,z,zbar,M,tau,y_lab,xi1,xi2,jac) !-- tau = Minv^2/s_had
    real(dp), intent(in) :: xx(1),M,z,zbar
    real(dp), intent(out) :: tau,y_lab,jac,xi1,xi2
    real(dp)              :: Msq,ximin,ximax,xi,sqrttau,ymax
    
    msq=M**2
    tau = msq/sh

    ymax = -half * log(tau)
    y_lab = ymax * (two * xx(1) - one)
    jac = two * ymax 

    ! now get xi1 & xi2 at the same time -- as before
    xi1 = sqrt(tau)*exp(+y_lab)/z
    xi2 = sqrt(tau)*exp(-y_lab)/zbar

    jac = jac/z/zbar

    return

  end subroutine get_tauyxi1xi2_OS

  subroutine get_tauyxi1xi2_OS_dec(xx,z,zbar,M,Ga,tau,y_lab,xi1,xi2,jac) !-- tau = Minv^2/s_had
    real(dp), intent(in) :: xx(1),M,Ga,z,zbar
    real(dp), intent(out) :: tau,y_lab,jac,xi1,xi2
    real(dp)              :: Msq,ximin,ximax,xi,sqrttau,ymax
    
    msq=M**2
    tau = msq/sh

    ymax = -half * log(tau)
    y_lab = ymax * (two * xx(1) - one)
    jac = two * ymax 

    ! now get xi1 & xi2 at the same time -- as before
    xi1 = sqrt(tau)*exp(+y_lab)/z
    xi2 = sqrt(tau)*exp(-y_lab)/zbar

    jac = jac/z/zbar

    jac = jac * M * Ga/sh * Pi

    return

  end subroutine get_tauyxi1xi2_OS_dec

end module mod_kinematics_gen

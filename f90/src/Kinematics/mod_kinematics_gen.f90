module mod_kinematics_gen
  use mod_types
  use mod_consts_dp
  use mod_proc_parms
  implicit none
  !-- dimension for vegas etc
  integer, public, parameter :: tauy_max = 2
  private

  public :: get_tauy

contains

  subroutine get_tauy(xx,tau,y_lab,jac) !-- tau = Minv^2/s_had
    real(dp), intent(in) :: xx(1:tauy_max)
    real(dp), intent(out) :: tau,y_lab,jac

    if(taumode == 1) then
      call get_tauy_newflat(xx,tau,y_lab,jac)   !-- flat, from min and max
    else if(taumode == 2) then
      call get_tauy_offshell(xx,tau,y_lab,jac)  !-- flat, from min and max
    else
      write(6,*) 'Invalid invariant mass generation mode. Options: 1/2. Aborting'
      stop
    endif

    return
  end subroutine get_tauy

  !--------------------------------------------------------------------------------

  !-- here q2min is indeed the minimum generated value
  subroutine get_tauy_newflat(xx,tau,y_lab,jac) !-- tau = Minv^2/s_had
    real(dp), intent(in) :: xx(2)
    real(dp), intent(out) :: tau,y_lab,jac
    real(dp) :: b,bmax,half_e,q2,ymax

    bmax = sqrt(one-q2min_tech/q2max_tech)
    b = bmax * xx(1)
    q2 = q2min_tech/(one-b**2)
    tau = q2/sh
    half_e = half*sqrt(q2)

    ymax = -half * log(tau)
    y_lab = ymax * (two * xx(2) - one)

    jac = four * tau * b * bmax * ymax/(one-b**2)

    return

  end subroutine get_tauy_newflat

  subroutine get_tauy_offshell(xx,tau,y_lab,jac)
    use mod_parms
    implicit none
    real(dp), intent(in) :: xx(2)
    real(dp), intent(out) :: tau,y_lab,jac
    real(dp) :: taumax, taumin
    real(dp) :: r, y, msh, a, b, wgt, ymax
    integer :: nu, p

    nu = 2
    p = nu - 1
    r = xx(1)
    y = xx(2)
    msh = mzsq/sh

    !-- invariant mass
    taumax = q2max_tech/sh
    taumin = q2min_tech/sh
    a = 1/(taumax - msh)**p
    b = 1/(taumin - msh)**p
    tau = one/(r*a + (one-r)*b) + msh
    wgt = (-one/p)*(tau-msh)**nu * (a - b)

    !-- rapidity
    ymax = -half * log(tau)
    y_lab = ymax * (two * y - one)

    jac = 2 * ymax * wgt

    return

  end subroutine get_tauy_offshell

  ! !-- careful about this routine, the max is not q2max but shad
  ! subroutine get_tauy_BW(xx,M,Ga,tau,y_lab,jac) !-- tau = Minv^2/s_had
  !   real(dp), intent(in) :: xx(2),M,Ga
  !   real(dp), intent(out) :: tau,y_lab,jac
  !   real(dp)              :: Msq,ximin,ximax,xi,sqrttau,ymax


  !   Msq=M**2
  !   ximin = +atan((q2min-Msq)/(M*Ga))
  !   ximax = +atan((sh-Msq)/(M*Ga))
  !   xi = Msq + Ga*M*tan(ximin + (ximax-ximin)*xx(1))
  !   tau = xi/sh
  !   sqrttau = sqrt(tau)
  !   jac = (ximax-ximin)/(M*Ga*sh)*((xi-Msq)**2+Msq*Ga**2)

  !   ymax = -half * log(tau)
  !   y_lab = ymax * (two * xx(2) - one)
  !   jac = jac * two * ymax

  !   return

  ! end subroutine get_tauy_BW

  ! subroutine get_tauy_OS(xx,M,tau,y_lab,jac) !-- tau = Minv^2/s_had
  !   real(dp), intent(in) :: xx(1),M
  !   real(dp), intent(out) :: tau,y_lab,jac
  !   real(dp)              :: msq,ymax
    
  !   msq=M**2
  !   tau = msq/sh

  !   ymax = -half * log(tau)
  !   y_lab = ymax * (two * xx(1) - one)
  !   jac = two * ymax 

  !   return

  ! end subroutine get_tauy_OS

end module mod_kinematics_gen

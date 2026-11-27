module mod_lumi
  use mod_types; use common_def; use mod_consts_dp
  implicit none
  private

  public :: get_pdf

contains

  subroutine get_pdf(z1,z2,mu,f1,f2)
    real(dp), intent(in) :: z1,z2,mu
    real(dp), intent(out) :: f1(-6:6),f2(-6:6)
    real(dp15) :: scale,x1,x2,f1_dp(-6:6),f2_dp(-6:6)

    f1 = zero
    f2 = zero

    x1 = real(z1,kind=dp15)
    x2 = real(z2,kind=dp15)
    
    scale = real(mu,kind=dp15)
    
    call evolvePDF(x1,scale,f1_dp)
    f1 = f1_dp/x1
    
    call evolvePDF(x2,scale,f2_dp)
    f2 = f2_dp/x2

    return

  end subroutine get_pdf

end module mod_lumi


module mod_limits_rv_dy_qg
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions; use mod_procfunctions
  use mod_splitting
  implicit none

  private

  public :: col_rv

contains

  !-- 1 -> PLOOP, 2 -> PTREE
  !-- Lmu = ln(musq/q2)
  function col_rv(Lmu,z,s14,eta41)
    real(dp), intent(in)  :: Lmu,z,s14,eta41
    real(dp)              :: col_rv(2)
    real(dp)              :: omz,lomz,lz,leta,Li2omz
    real(dp) :: ptree,pnew

    omz = one-z
    lz = log(z)
    lomz = log(omz)
    leta = log(eta41)
    Li2omz = dilog2(omz)

    pnew  = Pqgt_new_i(z)
    ptree  = Pqg_i(z)

    col_rv(1) = (Ca-Cf)*pnew + &
         ((leta/6.0_dp - Lmu/3.0_dp + lomz/6.0_dp - lz/3.0_dp)* &
         nf + Cf*(-2*Li2omz - 2*lomz**2 + 3*lz - lz**2 + lomz*(-3 + 2*lz) +  &
         leta*(-2*lomz + 2*lz) + (5*pisq)/3.0_dp) +  &
         Ca*(2*Li2omz + (11*Lmu)/6.0_dp - lomz/6.0_dp + lomz**2 +  &
         leta*(-5/3.0_dp + lomz) + (11*lz)/6.0_dp -  &
         (5*pisq)/6.0_dp))*ptree
    
    col_rv(2) = ptree

    col_rv = col_rv * two/s14

  end function col_rv
  
end module mod_limits_rv_dy_qg
    
  

module mod_limits_rv_dy_ns
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions; use mod_procfunctions
  use mod_splitting
  implicit none

  private

  public :: col_rv,soft_rv,softcol_rv

contains

  !-- 1 -> PLOOP, 2 -> PTREE
  !-- Lmu = ln(musq/q2)
  function col_rv(Lmu,z,s14,eta41)
    real(dp), intent(in)  :: Lmu,z,s14,eta41
    real(dp)              :: col_rv(2)
    real(dp)              :: omz,lomz,lz,leta,Li2omz
    real(dp) :: pnew,ptree

    omz = one-z
    lz = log(z)
    lomz = log(omz)
    leta = log(eta41)
    Li2omz = dilog2(omz)

    pnew  = Pqqt_new_i(z)
    ptree = Pqq_i(z)

    col_rv(1) = tag_ns_gg*(Ca-Cf) * pnew + &
         (Cf*(2*Li2omz + 3*lz + 2*leta*lz + 2*lomz*lz - lz**2)*tag_ns_gg +  &
         (leta/6.0_dp - Lmu/3.0_dp + lomz/3.0_dp - lz/3.0_dp)* &
         nf*tag_ns_nf + Ca*tag_ns_gg*(-2*Li2omz + (11*Lmu)/6.0_dp +  &
         leta*(-5/3.0_dp - lomz) - (10*lomz)/3.0_dp - lomz**2 + &
         (11*lz)/6.0_dp + (5*pisq)/6.0_dp)*tag_ns_gg)*ptree

    col_rv(2) = ptree*tag_ns_gg

    col_rv = col_rv * two/s14

  end function col_rv
  
  function soft_rv(musq,s12,s14,s24)
    real(dp), intent(in)  :: musq,s12,s14,s24
    real(dp)              :: soft_rv(2)
    real(dp)              :: eik,Lmu,l14,l24

    Lmu = log(musq/s12)
    l14 = log(s14/s12)
    l24 = log(s24/s12)

    eik= four*Cf*s12/s14/s24

    soft_rv(1) = (11.0_dp/6.0_dp*Ca*tag_ns_gg-two/three*TR*nf*tag_ns_nf)*(-l14/2.0_dp + Lmu - l24/2.0_dp) + &
         Ca*tag_ns_gg*(l14*(-3.0_dp/4.0_dp - l24) - (3*l24)/4.0_dp + (5*zeta2))

    soft_rv(2) = one*tag_ns_gg

    soft_rv = soft_rv * eik

  end function soft_rv

  !-- c41 + soft
  function softcol_rv(musq,s12,s14,E4oE2)
    real(dp), intent(in)  :: musq,s12,s14,E4oE2
    real(dp)              :: softcol_rv(2)
    real(dp)              :: eik,lmu,l14,l24

    lmu = log(musq/s12)
    l14 = log(s14/s12)
    l24 = log(E4oE2)

    eik= four*Cf/s14/E4oE2

    softcol_rv(1) = (11.0_dp/6.0_dp*Ca*tag_ns_gg-two/three*TR*Nf*tag_ns_nf)*(-l14/2.0_dp + Lmu - l24/2.0_dp) + &
         Ca*tag_ns_gg*(l14*(-3.0_dp/4.0_dp - l24) - (3*l24)/4.0_dp + (5*zeta2))
    softcol_rv(2) = one*tag_ns_gg

    softcol_rv = softcol_rv * eik

  end function softcol_rv

end module mod_limits_rv_dy_ns
    
  

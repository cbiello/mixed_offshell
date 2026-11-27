!-- AP and calligraphic P' splitting functions, for subtractions
module mod_splitting_aux
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions; use mod_procfunctions
  implicit none

  private

  !-- notation: [-1 --> delta, 0 --> reg including plus, +1 --> plus]
  public :: PggAP,PqqAP,PqgAP,PgqAP
  public :: calPp_gg,calPp_qq,calPp_qg,calPp_gq

contains

  function PggAP(z) result(res)
    real(dp) :: z
    real(dp) :: res(-1:1)

    res(-1) = b0

    res( 0) = two*Ca*((one-z)/z+z/(one-z)+z*(one-z))

    res(+1) = two*Ca/(one-z)
    
  end function PggAP

  function PqqAP(z) result(res)
    real(dp) :: z
    real(dp) :: res(-1:1)

    res(-1) = half*3*Cf

    res( 0) = Cf * (one+z**2)/(one-z)

    res(+1) = two*Cf/(one-z)

  end function PqqAP

  function PqgAP(z) result(res)
    real(dp) :: z
    real(dp) :: res(-1:1)

    res = 0

    res(0) = tr * (z**2 + (one-z)**2)

  end function PqgAP

  function PgqAP(z) result(res)
    real(dp) :: z
    real(dp) :: res(-1:1)

    res = 0
    
    res(0) = Cf * (one+(one-z)**2)/z

  end function PgqAP

  !--

  function calPp_gg(z) result(res)
    real(dp) :: z
    real(dp) :: res(-1:1)

    res(-1) = zero

    res( 0) = -two*Ca*((one-z)/z+z/(one-z)+z*(one-z))*log((one-z)**2/z)

    res(+1) = -four*Ca*log(one-z)/(one-z)

  end function calPp_gg

  function calPp_qq(z) result(res)
    real(dp) :: z
    real(dp) :: res(-1:1)

    res(-1) = zero

    res( 0) = -Cf*(one+z**2)/(one-z)*log((one-z)**2/z) - Cf * (one-z)

    res(+1) = -four*Cf*log(one-z)/(one-z)

  end function calPp_qq

  function calPp_qg(z) result(res)
    real(dp) :: z
    real(dp) :: res(-1:1)

    res(-1) = zero

    res( 0) = -tr*((z**2 + (one-z)**2)*(log((one-z)**2/z)-one) + one)

    res(+1) = zero

  end function CalPp_qg

  function calPp_gq(z) result(res)
    real(dp) :: z
    real(dp) :: res(-1:1)

    res = 0
    
    res(0) = -Cf * (one+(one-z)**2)/z*log((one-z)**2/z) - Cf * z

  end function CalPp_gq


end module mod_splitting_aux

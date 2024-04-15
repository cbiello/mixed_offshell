module mod_consts_dp
  use mod_types
  implicit none
  private

  ! -------------------------------------------------------------------------

  !-- numerical constants 
  real(dp), public, parameter :: pi    = 3.141592653589793238462643383279502884197_dp
  real(dp), public, parameter :: twopi = 6.283185307179586476925286766559005768394_dp
  real(dp), parameter, public :: pisq  = 9.869604401089358618834490999876151135314_dp
  !
  real(dp), public, parameter :: zeta2 = 1.644934066848226436472415166646025189219_dp
  real(dp), public, parameter :: zeta3 = 1.202056903159594285399738161511449990765_dp
  !
  real(dp), parameter, public :: ln2 = log(2.0_dp)
  real(dp), parameter, public :: li4half = 0.5174790616738993863307581618988629456224_dp ! Li_4[1/2]
  !
  real(dp), public, parameter :: zero  = 0.0_dp
  real(dp), public, parameter :: half  = 0.5_dp
  real(dp), public, parameter :: one   = 1.0_dp
  real(dp), public, parameter :: two   = 2.0_dp
  real(dp), public, parameter :: three = 3.0_dp
  real(dp), public, parameter :: four  = 4.0_dp
  real(dp), public, parameter :: five  = 5.0_dp
  real(dp), public, parameter :: six   = 6.0_dp
  !
  real(dp), public, parameter :: sqrt2 = 1.4142135623730950488016887242096980785696718753769_dp 
  !
  complex(dp), public, parameter :: ci    = cmplx(0.0_dp,1.0_dp,kind=dp)
  complex(dp), public, parameter :: czero = cmplx(0.0_dp,0.0_dp,kind=dp)
  complex(dp), public, parameter :: cone  = cmplx(1.0_dp,0.0_dp,kind=dp)
  complex(dp), public, parameter :: ctwo  = cmplx(2.0_dp,0.0_dp,kind=dp)

  !-- units
  real(dp), public, parameter :: GeVtoFb = 0.389379E12_dp
  real(dp), public, parameter :: GeVtoPb = 0.389379E9_dp
  real(dp), public, parameter :: GeVtoNb = 0.389379E6_dp

  !-- length of string for names
  integer, public, parameter :: strlen = 100
  
end module mod_consts_dp

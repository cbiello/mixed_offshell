module mod_limvals
  use mod_types
  real(dp), public :: limval_nlo(4)
  real(dp), public :: limval_nlo_is(6)
  real(dp), public :: limval_nnlo(16)
  real(dp), public :: limval_onlo(18)   !-- ONLO S sectors: per kinematic group g, slot 6(g-1)+k, k: H, S, C1 (C), C2, SC1 (SC), SC2

end module mod_limvals
  

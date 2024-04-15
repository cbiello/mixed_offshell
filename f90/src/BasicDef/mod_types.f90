!----------------------------------------------------------------------
! Defines kind parameters
module mod_types
  implicit none
  integer, parameter  :: dp = selected_real_kind(15)
  !integer, parameter  :: dp = selected_real_kind(30)
  integer, parameter  :: sp = selected_real_kind(6)
  integer, parameter  :: qp = selected_real_kind(30)
  
  integer, parameter  :: dp15 = selected_real_kind(15)
  integer, parameter  :: vp   = selected_real_kind(15)

end module mod_types

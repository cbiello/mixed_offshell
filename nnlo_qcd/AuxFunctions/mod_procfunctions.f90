module mod_procfunctions
  use mod_types; use mod_consts_dp; use common_def
  implicit none
  private

  public :: get_prefactor

contains

  subroutine get_prefactor(muR,as_ord,pref)
    real(dp), intent(in)    :: muR
    integer,intent(in)      :: as_ord
    real(dp), intent(out)   :: pref
    real(dp)                :: as,asontwopi,mu
    real(dp15)              :: alphasPDF

    as = alphasPDF(real(muR,kind=dp15))
    asontwopi = as/(two*pi)

    pref = asontwopi**as_ord * GeVtoFb

  end subroutine get_prefactor

end module mod_procfunctions

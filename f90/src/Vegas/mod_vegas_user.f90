!-- the routine vegas_user(..) is called after each iteration
!-- and can be used to do implement additional feature
!-- at the moment it is used to process histograms
module mod_vegas_user
  use mod_types
  use mod_mpi_common
  implicit none

  private

  public :: vegas_user

contains

  subroutine vegas_user(init,it,ti,res,tsi,err,chi2a)
    use mod_histo
    implicit none
    integer, intent(in) :: init,it
    real(vp) :: ti,res,tsi,err,chi2a

    call prepare_and_write(init,it,ti,res,tsi,err,chi2a)

    return
  end subroutine vegas_user

end module mod_vegas_user

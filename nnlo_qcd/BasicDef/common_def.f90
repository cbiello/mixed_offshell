module common_def
  use mod_types
  implicit none
  private

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! --- maximum number of histograms
  integer, public, parameter :: nobs_max = 50

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! --- generic
  real(dp), public :: sh
  real(dp), public :: musq,musq_in,xmuR,xmuF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dp), public :: buff,buffz,bufftc,buffrr,buffrv
  real(dp), public :: onet
  integer, public  :: buffrr_pow = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, public :: maxiternumber

  character*100, public :: outputfile !-- the file with histograms
  character*200, public :: histotmpfile !-- copy intermediate results for histograms not to loose them

  !-- command-line vegas parameters
  integer, public :: vegasNc0 = -1
  integer, public :: vegasNc1 = -1 
  integer, public :: vegasIt0 = -1
  integer, public :: vegasIt1 = -1
  integer, public :: seed
  character, public :: runtype = "w" !-- w(armup)/r(un)

  character, public :: pdffile*(50)

  logical, public, parameter :: withmax = .true.
  logical, public, parameter :: tcanalytic = .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !-- various dimensions
#if (_PhaseSpaceType == 3) || (_PhaseSpaceType == 4) || (_PhaseSpaceType == 5)  || (_PhaseSpaceType == 6)
  integer, public, parameter :: tauy_max = 1
  !
  integer, public, parameter :: kLO_min = 2
  integer, public, parameter :: kLO_max = 3
  integer, public, parameter :: kLO_max_full = kLO_max + 1 !-- random numbe
  !
  integer, public, parameter :: LOdim = 3  
  !
  integer, public, parameter :: xE4 = 2
  integer, public, parameter :: xE5 = 3
  integer, public, parameter :: xRHO1 = 4
  integer, public, parameter :: xRHO2 = 6
  !
  integer, public, parameter :: xE_nlo = 2
  integer, public, parameter :: xRHO_nlo = 3
  !
#else
  integer, public, parameter :: tauy_max = 2
  !
  integer, public, parameter :: kLO_min = 3
  integer, public, parameter :: kLO_max = 4
  integer, public, parameter :: kLO_max_full = kLO_max + 1 !-- random number
  !
  integer, public, parameter :: LOdim = 4
  !
  integer, public, parameter :: xE4 = 3
  integer, public, parameter :: xE5 = 4
  integer, public, parameter :: xRHO1 = 5
  integer, public, parameter :: xRHO2 = 7
  !
  integer, public, parameter :: xE_nlo = 3
  integer, public, parameter :: xRHO_nlo = 4
  !
#endif
  integer, public, parameter :: kNLO_min = kLO_min
  integer, public, parameter :: kNLO_max = kLO_max + 3
  integer, public, parameter :: kNLO_max_full = kNLO_max + 1 !-- random number
  !
  integer, public, parameter :: kNNLO_min = kLO_min
  integer, public, parameter :: kNNLO_max = kLO_max + 6
  integer, public, parameter :: kNNLO_max_full = kNNLO_max + 1 !-- random number

end  module common_def


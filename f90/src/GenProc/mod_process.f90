module mod_process
  use mod_types
  implicit none
  integer, parameter :: nmax = 6, zmax = 3, kininvmax = 7, kininvmax2 = 8

  type :: KinConfig
     real(dp) :: AmpMom(4,nmax)
     real(dp) :: LimMom(4,nmax)
     real(dp) :: Lim_etaij(nmax,nmax)
     real(dp) :: Lim_Ei(nmax)
     real(dp) :: Lim_sij(nmax,nmax)
     real(dp) :: Lim_z(zmax)
     real(dp) :: Lim_KinInv(kininvmax)
     real(dp) :: Lim_KinInv2(kininvmax2)
     integer  :: ids(nmax) !-- is is 0, only filled the one we use, use 1 for fs quark
     real(dp) :: PartFrac(2)
     real(dp) :: mu
     real(dp) :: mur(1),muf(1) !-- central,lower,upper
     real(dp) :: mu2ref !-- reference scale^2, for Lmu terms in the subtraction
     real(dp) :: wgt
     integer  :: npart  !-- how many particles in this configuration
     logical  :: flag    !-- technical flag from kinematic generation
     logical  :: makecut !-- actual fiducial cuts
     integer, allocatable :: part(:) ! part = (id_1, id_2, id_3, ... id_6)
                           !         id_3 = ( +/- id_el, +/- id_nu )
                           ! in Parms/mod_parms.f90 
                           ! integer, public, parameter :: id_el = 11 !-- + is always particle

  end type KinConfig

contains

  subroutine initialize_config(AKinConfig)
    !set all properties of KinConfig to initial values (mostly zero)
    type(KinConfig) :: AKinConfig
    
    AKinConfig%AmpMom = 0
    AKinConfig%LimMom = 0
    AKinConfig%Lim_etaij = 0
    AKinConfig%Lim_Ei = 0
    AKinConfig%Lim_sij = 0
    AKinConfig%Lim_z = 0
    AKinConfig%Lim_KinInv = 0
    AKinConfig%PartFrac = 2
    AKinConfig%mu  = 0
    AKinConfig%mu2ref  = -99
    AKinConfig%mur  = 0
    AKinConfig%muf  = 0
    AKinConfig%wgt = 0
    AKinConfig%ids = -99
    AKinConfig%flag    = .false.
    AKinConfig%makecut = .true.

  end subroutine initialize_config
    
end module mod_process



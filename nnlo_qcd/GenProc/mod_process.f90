module mod_process
  use mod_types
  implicit none

  type PartChannel
     integer :: ch(1:2)
  end type PartChannel

  type :: KinConfig
     real(dp),allocatable :: AmpMom(:,:)
     real(dp),allocatable :: LimMom(:,:)
     real(dp),allocatable :: KinInv(:)
     real(dp),allocatable :: EikInv(:)
     real(dp) :: PartFrac(2)
     real(dp) :: wgt
     logical :: flag
     integer :: NumLept
     integer :: numpartch
     type(PartChannel),allocatable  :: AllPartCh(:)

  end type KinConfig

  type KinSet
     type(KinConfig),allocatable :: AllKinConfig(:)
  end type KinSet

!!!                                                                                                                                                                                                         
  type(KinSet),public,target :: MyKinSet


contains


  subroutine initialize_config(AKinConfig)
    !set all properties of KinConfig to initial values (mostly zero)
    type(KinConfig) :: AKinConfig

    AKinConfig%AmpMom = 0.0_dp
    AKinConfig%LimMom = 0.0_dp
    AKinConfig%KinInv = 0.0_dp
!    AKinConfig%EikInv = 0.0_dp
    AKinConfig%PartFrac = 2.0_dp

    AKinConfig%wgt = 2.0_dp
    AKinConfig%flag = .false.
    
  end subroutine initialize_config
    
end module mod_process



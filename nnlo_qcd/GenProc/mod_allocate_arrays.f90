module mod_allocate_arrays
  use mod_types; use common_def; use mod_consts_dp
  use mod_process
  implicit none

contains

  subroutine allocate_arrays_lo()
    print *, "allocating arrays for LO"

    allocate(MyKinSet%AllKinConfig(1))

    allocate(MyKinSet%AllKinConfig(1)%LimMom(4,4))
    allocate(MyKinSet%AllKinConfig(1)%AmpMom(4,4))
    allocate(MyKinSet%AllKinConfig(1)%KinInv(1))

    MyKinSet%AllKinConfig(1)%NumLept=2

  end subroutine allocate_arrays_lo

  subroutine allocate_arrays_nlo(sec)
    character, intent(in) :: sec*(10)
    integer :: j

    print *, "allocating arrays for NLO"

    if (sec .eq. 'r') then

       print *, 'allocate NLO r for'

       allocate(MyKinSet%AllKinConfig(4)) 
       !
       allocate(MyKinSet%AllKinConfig(1)%LimMom(4,5)) !-- hard
       allocate(MyKinSet%AllKinConfig(1)%AmpMom(4,5))
       allocate(MyKinSet%AllKinConfig(1)%KinInv(1))   !-- q2
       !
       allocate(MyKinSet%AllKinConfig(2)%LimMom(4,5)) !-- soft
       allocate(MyKinSet%AllKinConfig(2)%AmpMom(4,4))
       allocate(MyKinSet%AllKinConfig(2)%KinInv(3))   !-- s12,s14,s24
       !
       allocate(MyKinSet%AllKinConfig(3)%LimMom(4,5)) !-- C1
       allocate(MyKinSet%AllKinConfig(3)%AmpMom(4,4))
       allocate(MyKinSet%AllKinConfig(3)%KinInv(3))   !-- z,s14,q2
       !
       allocate(MyKinSet%AllKinConfig(4)%LimMom(4,5)) !-- C2
       allocate(MyKinSet%AllKinConfig(4)%AmpMom(4,4))
       allocate(MyKinSet%AllKinConfig(4)%KinInv(3))   !-- z,s14,q2

       do j=1,4
          MyKinSet%AllKinConfig(j)%NumLept=2
       enddo

    elseif (sec .eq. 's' .or. sec .eq. 'v') then

       print *, 'allocate NLO virtual/subtractions'

       allocate(MyKinSet%AllKinConfig(1))
       !
       allocate(MyKinSet%AllKinConfig(1)%LimMom(4,4))    
       allocate(MyKinSet%AllKinConfig(1)%AmpMom(4,4))    
       allocate(MyKinSet%AllKinConfig(1)%KinInv(1))

       MyKinSet%AllKinConfig(1)%NumLept=2

    else

       print *, 'wrong sector'
       stop

    endif

  end subroutine allocate_arrays_nlo

  subroutine allocate_arrays_nnlo(sec,ch)
    character, intent(in) :: sec*(10),ch*(10)
    integer :: j,alloc_counter

    print *, "allocating arrays for NNLO"

    alloc_counter = 0 

    if (ch .eq. 'ns' .or. ch .eq. 's') then

       if (sec(1:2) .eq. 'rr') then
          include "./include/allocate_rr_dy_ns.f90"
       elseif (sec(1:2) .eq. 'rv') then
          include "./include/allocate_rv_dy.f90"
       elseif (sec .eq. 'sub124') then
          include "./include/allocate_sub124_dy_ns.f90"
       elseif (sec .eq. 'sub12') then
          include "./include/allocate_vv_dy_ns.f90"
       elseif (sec .eq. 'subv') then
          include "./include/allocate_vv_dy_ns.f90"
       elseif (sec .eq. 'vv') then
          include "./include/allocate_vv_dy_ns.f90"
       else

          print *, 'wrong allocation in ns/s'
          stop

       endif

    elseif (ch .eq. 'qg' .or. ch .eq. 'gq') then

       if (sec(1:2) .eq. 'rr') then
          include "./include/allocate_rr_dy_gq.f90"
       elseif (sec(1:2) .eq. 'rv') then
          include "./include/allocate_rv_dy.f90"
       elseif (sec .eq. 'sub124') then
          include "./include/allocate_sub124_dy_gq.f90"
       elseif (sec .eq. 'sub12') then
          include "./include/allocate_vv_dy_ns.f90"
       elseif (sec .eq. 'subv') then
          include "./include/allocate_vv_dy_ns.f90"
       else

          print *, 'wrong allocation in gq/qg'
          stop

       endif

    elseif (ch .eq. 'gg') then

       if (sec(1:2) .eq. 'rr') then
          include "./include/allocate_rr_dy_gg.f90"
       elseif (sec .eq. 'sub12') then
          include "./include/allocate_vv_dy_ns.f90"
       elseif (sec .eq. 'sub124') then
          include "./include/allocate_sub124_dy_gg.f90"
       else
          

          print *, 'wrong allocation in gg'
          stop

       endif
       
    else

       print *, 'allocating wrong channel'
       stop

    endif 

  end subroutine allocate_arrays_nnlo

end module mod_allocate_arrays



module mod_mpi_common
    
  integer, public, allocatable :: mpi_stat(:)
  integer, public, parameter :: mpi_root_process = 0
  logical, public :: root_process = .true.
  integer, public :: mpi_size = 1
  integer, public :: mpi_rank = 0
  integer, public :: ierr
  
  !-- time
  real :: mpi_overhead = 0.0
  real :: mpi_worker_time
  
end module mod_mpi_common

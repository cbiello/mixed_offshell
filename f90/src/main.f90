program main
  use mod_chooser
  use mod_mpi_common
#if ( _MPIbuild )
  use mpi
#endif
  !-- GNU
  real :: start,finish
    
#if ( _MPIbuild )
  call mpi_init( ierr )
  call mpi_comm_rank( mpi_comm_world, mpi_rank, ierr )
  call mpi_comm_size( mpi_comm_world, mpi_size, ierr )
  if ( mpi_rank .ne. mpi_root_process ) root_process = .false.
  allocate( mpi_stat(mpi_status_size) )
#endif

  if (root_process) then
     print *, '####################################'
     print *, '###  QCD-EW dilepton production  ###'
     print *, '####################################'
#if ( _MPIbuild )
     write(*,'(/,a,i0,a,/)') '# MPI build (# processes: ', mpi_size, ')'
#endif
  end if
     
  call cpu_time(start)

  call chooser()

#if ( _MPIbuild )
  deallocate( mpi_stat )
  call mpi_finalize( ierr )
#endif
  
  call cpu_time(finish)

  if (root_process) then
     print *, '#'
     print *, '# Time = ',finish-start
     print *, '#'
  endif

end program main

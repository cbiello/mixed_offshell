module mod_proc_histo
  use mod_parser
  use mod_consts_dp
  use mod_parms
  use mod_vegas_parms
  use mod_proc_parms
  use mod_mpi_common
  implicit none
  private
  
  public :: print_proc_histo,reset_iter_proc
#if (_MPIbuild)
  public :: mpi_proc_reduce
#endif

contains

  !-- at the end of each iteration, print stuff
  subroutine print_proc_histo(idev)
    integer, intent(in) :: idev
    character(2) :: str_units
    
    write(idev,*) '# ', adjustl(trim(command_line()))
    write(idev,*) '# '
#if defined(GIT_VERSION)
    write(idev,'(A,A)') " # Git version: ", GIT_VERSION
#else
    write(idev,'(A)') " # Git version: unknown"
#endif
#if defined(ANALYSIS)
    write(idev,'(A,A)') " # Analysis: ", ANALYSIS
#else
    write(idev,'(A)') " # Analysis: unknown"
#endif    
    write(idev,*) '# '
    
    call print_vegas_parms(idev)
    call print_ew_parms(ew_scheme,idev)
    call print_qcd_parms(idev)
    call print_proc_parms(idev)
    call print_counters(idev,iter=.true.,acc=.true.)
    
    if (units==GeVtoFb) str_units = 'fb'
    if (units==GeVtoPb) str_units = 'pb'
    if (units==GeVtoNb) str_units = 'nb'
    write(idev,*) '# Units: ', str_units
    write(idev,*) '# '
    
  end subroutine print_proc_histo

  !-- at the end of each iteration, reset stuff for the next iteration
  subroutine reset_iter_proc()

    acc_icount_nan    = acc_icount_nan    + icount_nan   
    acc_failed_points = acc_failed_points + failed_points
    acc_icount_good   = acc_icount_good   + icount_good  
    acc_icount_acc    = acc_icount_acc    + icount_acc   
    
    call reset_counters(all=.false.) !-- all user counters

  end subroutine reset_iter_proc

#if (_MPIbuild)
  !-- gather information from all mpi workers
  subroutine mpi_proc_reduce()
    use mpi
    integer :: icount_nan_tot,icount_good_tot,icount_acc_tot,failed_points_tot
    integer :: acc_icount_nan_tot,acc_icount_good_tot,acc_icount_acc_tot,acc_failed_points_tot

    call mpi_reduce( icount_nan,    icount_nan_tot, 1, mpi_integer, mpi_sum, mpi_root_process, mpi_comm_world, ierr )
    call mpi_reduce( icount_good,   icount_good_tot, 1, mpi_integer, mpi_sum, mpi_root_process, mpi_comm_world, ierr )
    call mpi_reduce( icount_acc,    icount_acc_tot, 1, mpi_integer, mpi_sum, mpi_root_process, mpi_comm_world, ierr )
    call mpi_reduce( failed_points, failed_points_tot, 1, mpi_integer, mpi_sum, mpi_root_process, mpi_comm_world, ierr )

    call mpi_reduce( acc_icount_nan,    acc_icount_nan_tot, 1, mpi_integer, mpi_sum, mpi_root_process, mpi_comm_world, ierr )
    call mpi_reduce( acc_icount_good,   acc_icount_good_tot, 1, mpi_integer, mpi_sum, mpi_root_process, mpi_comm_world, ierr )
    call mpi_reduce( acc_icount_acc,    acc_icount_acc_tot, 1, mpi_integer, mpi_sum, mpi_root_process, mpi_comm_world, ierr )
    call mpi_reduce( acc_failed_points, acc_failed_points_tot, 1, mpi_integer, mpi_sum, mpi_root_process, mpi_comm_world, ierr )

    if (root_process) then
       icount_nan    = icount_nan_tot
       icount_good   = icount_good_tot
       icount_acc    = icount_acc_tot
       failed_points = failed_points_tot
       !
       acc_icount_nan    = acc_icount_nan_tot
       acc_icount_good   = acc_icount_good_tot
       acc_icount_acc    = acc_icount_acc_tot
       acc_failed_points = acc_failed_points_tot
    else
       call reset_iter_proc() !-- reset everything user may want to reset
    end if
  end subroutine mpi_proc_reduce
#endif

end module mod_proc_histo

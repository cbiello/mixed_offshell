module mod_histo
  !-- for the time being, no cumulative histogram implemented
  !-- careful, changed step determination, to be tested
  use mod_types
  use mod_consts_dp
  use mod_auxfunctions
  use mod_proc_histo
  !-- one below for outputfile
  use mod_proc_parms
  use mod_mpi_common
  !--
  implicit none
  private
  integer, parameter :: nkin_max = 18 !-- how many different kinematics points per event
  integer, parameter :: nstring  = 15 !-- length of observable name
  real(dp), parameter :: tiny=1.0e-30_dp
  integer, save :: n_histo   = 0 !-- total number of histograms
  integer, save :: evalcount = 0 !-- old code -> set to 1, I think it was a bug
  logical, public, save :: nohistos = .false. !-- to avoid during histograms during warmup etc
  integer, public, parameter :: nweights = 3 !TODO: do this at the preprocessor level, with #define myweights = 3; neweights = myweights
  !--
  type histo
     character(nstring) :: name
     real(dp)      :: low,up,step
     integer       :: nbin
     !
     real(dp), allocatable :: tmp(:,:) !-- first: #weight. second: bin number
     real(dp), allocatable :: res(:,:) !-- results per iteration
     real(dp), allocatable :: err(:,:) 
     real(dp), allocatable :: acc_res(:,:) !-- accumulated results
     real(dp), allocatable :: acc_err(:,:)
     !
     integer       :: used_bin(nkin_max)
     integer       :: n_of_used_bins
     logical       :: discard_event
     !
     type(histo), pointer :: next_histo
  end type histo
  real(dp), public, allocatable, save :: obs(:)
  !--
  type(histo), target, save :: head
  type(histo), pointer,save :: current

  public :: histo
  public :: init_histo,setup_histo,new_histo,clear_histo
  public :: open_histo,fill_histo,close_histo
  public :: prepare_and_write
  
contains

  subroutine init_histo()

    current => head
    nullify(current%next_histo)

  end subroutine init_histo

  subroutine setup_histo(idev_in)
    integer, optional, intent(in) :: idev_in
    integer :: idev,i

    !-- allocate the observable vector
    allocate(obs(n_histo))

    !-- print information
    if (present(idev_in)) then
       idev = idev_in
    else
       idev = 6
    endif

    if (root_process) then
       current => head
       write(idev,*) '# '
       write(idev,*) '# number of histograms initialised = ',n_histo
       write(idev,*) '# name, lower bound, upper bound, number of bins'
       write(idev,*) '# last two bins: overflow, underflow in this order'
       write(idev,*) '# '
       do i = 1,n_histo
          current => current%next_histo
          write(idev,*) '# ',current%name, current%low, current%up, current%nbin
       enddo
       write(idev,*) '# '
    endif
       
    call clear_histo()
    
  end subroutine setup_histo

  !--
  
  subroutine new_histo(namein,low,up,step,cml)
    character(len=*), intent(in) :: namein
    real(dp), intent(in) :: low,up,step
    logical, optional, intent(in) :: cml
    type(histo), pointer :: temp
    logical :: cml_histo
    integer :: nbin

    if (present(cml)) then
       cml_histo = cml
    else
       cml_histo = .false.
    endif

    if (cml_histo) then
       write(6,*) 'Cumulative histograms not yet implemented'
       stop
    endif
    
    allocate(temp)

    nbin = floor((up-low)/step) !-- approximate by closest lower integer

    temp%name=trim(namein)
    !
    temp%low = low
    temp%up  = up
    temp%step = (up-low)/nbin
    !
    temp%nbin = nbin
    !
    allocate(temp%tmp(nweights,1:nbin+2)) !-- last two bins -> overflow/underflow
    allocate(temp%res(nweights,1:nbin+2))
    allocate(temp%err(nweights,1:nbin+2))
    allocate(temp%acc_res(nweights,1:nbin+2))
    allocate(temp%acc_err(nweights,1:nbin+2))
    !
    temp%tmp = zero
    temp%res = zero
    temp%err = zero
    temp%acc_res = zero
    temp%acc_err = zero
    !
    temp%n_of_used_bins = 0
    temp%discard_event = .false.
    !
    nullify(temp%next_histo)

    current%next_histo => temp
    current => current%next_histo
    
    n_histo = n_histo + 1

  end subroutine new_histo

  subroutine clear_histo()
    integer :: i 

    current => head

    do i = 1,n_histo

       current => current%next_histo

       current%tmp = zero
       current%res = zero
       current%err = zero
       current%acc_res = zero
       current%acc_err = zero
       !
       current%n_of_used_bins = 0
       current%discard_event = .false.

    enddo

  end subroutine clear_histo

  !-- 
  
  subroutine open_histo()
    integer :: i

    if (nohistos) return
    
    current => head
    
    do i = 1,n_histo
       !
       current => current%next_histo
       !
       current%n_of_used_bins = 0
       current%discard_event = .false.
       !
    enddo

  end subroutine open_histo

  subroutine fill_histo(res_diff,vegasweight)
    real(dp), intent(in)   :: res_diff(:)
    real(dp15), intent(in) :: vegasweight
    real(dp) :: myweight(nweights)
    integer :: i,mybin

    if (nohistos) return

    myweight = res_diff*vegasweight

    current => head
    
    do i = 1,n_histo
       
       current => current%next_histo

       if ( is_nan(obs(i)) .or. is_nan(myweight) ) then

          current%discard_event = .true.
          
       else

          !-- find the bin.
          !-- use floor: floor(1.789) = 1, floor(-1.789) = -2
          !-- careful with int: int(-1.789) = -1 and this gets the underflow wrong
          mybin = floor((obs(i)-current%low)/current%step) + 1
          if (mybin .gt. current%nbin+1) mybin = current%nbin + 1 !-- overflow
          if (mybin .le. zero)           mybin = current%nbin + 2 !-- underflow
          !
          current%n_of_used_bins = current%n_of_used_bins + 1
          current%used_bin(current%n_of_used_bins) = mybin

          !-- bin the event
          current%tmp(:,mybin) = current%tmp(:,mybin) + myweight(:)

       endif

    enddo

  end subroutine fill_histo

  subroutine close_histo()
    integer :: i,j,ibin

    if (nohistos) return
    
    current => head
    
    do i = 1,n_histo
       !
       current => current%next_histo
       !
       if (.not. current%discard_event) then
          !
          do j = 1, current%n_of_used_bins !-- only update modified bins
             ibin = current%used_bin(j)
             current%res(:,ibin) = current%res(:,ibin) + current%tmp(:,ibin)
             current%err(:,ibin) = current%err(:,ibin) + current%tmp(:,ibin)**2
             current%tmp(:,ibin) = zero
          enddo
          !
       endif
       !
    enddo

    evalcount = evalcount + 1

  end subroutine close_histo

  !--

  !-- do not care about optimization, this routine is only called once per iteration
  subroutine prepare_and_write(init,it,it_int,tot_int,it_err,tot_err,chi2a)
    !
    integer, intent(in)  :: init,it
    real(vp), intent(in) :: it_int,tot_int,it_err,tot_err,chi2a
    integer, parameter :: idev = 66
    logical, parameter :: writetofile = .true.
    !
    integer :: i,j,k
    real(dp) :: res(nweights),err(nweights),tmp(nweights)
    real(dp) :: acc_res(nweights),acc_err(nweights),tmp_err(nweights)

    if (nohistos) return
 
#if ( _MPIbuild )
    call mpi_gather_histos()
    if (.not.root_process) return
#endif   
    
    current => head

    do i = 1,n_histo
       !
       current => current%next_histo

       do j = 1,current%nbin+1 

          res(:) = current%res(:,j)
          err(:) = current%err(:,j)

          if (err(1) == 0) cycle !-- nothing to do if everything was untouched, take the first weight as a driver

          !-- do the same as vegas
          tmp = sqrt(err*evalcount)
          err = sqrt((tmp-res)*(tmp+res)/(evalcount-1))

          if (err(1) <= zero) err = tiny

          acc_res(:) = current%acc_res(:,j)
          acc_err(:) = current%acc_err(:,j)

          !-- the first time, there is no existing error
          if (acc_err(1) == 0) then

             acc_res = res
             acc_err = err

          else
          
             tmp_err = one/acc_err**2 + one/err**2
             acc_res = (acc_res/acc_err**2 + res/err**2)/tmp_err
             
             acc_err = one/sqrt(tmp_err)

          endif

          current%acc_res(:,j) = acc_res
          current%acc_err(:,j) = acc_err
          
       enddo
       !
       current%res = zero
       current%err = zero
       !
    enddo
    
    !-- write to file
    if (writetofile) open(unit=idev,file=outputfile,status='unknown')
    !
    call print_proc_histo(idev)
    !
    if (init==0) write(idev,*) '# warmup run, iteration', it
    if (init==1) write(idev,*) '# final run,  iteration', it
    write(idev,*) '# '
    write(idev,*) '# Iteration result:    ', it_int
    write(idev,*) '# Iteration error:     ', it_err
    write(idev,*) '#'
    write(idev,*) '# Accumulated result:  ', tot_int
    write(idev,*) '# Accumulated error:   ', tot_err
    write(idev,*) '# Chi2/Iteration: ', chi2a

    current => head

    do i = 1,n_histo

       current => current%next_histo

       write(idev,*) '#' !-- two spaces --> for Gavin's mergeidx.pl
       write(idev,*) '#' 

       write(idev,*) '# ----------Obs=',current%name,'----------------'
       write(idev,*) '# observable, bin number, left edge, result, error'

       do j = 1,current%nbin+1

          !-- TODO: arrange format in a better way, to avoid all this cut and paste
          if(nweights==1) then
             write(idev,"(I2,A,X,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)",advance='yes') i,' ', &
                  j, ' ',current%low+current%step*(j-1),' ',current%acc_res(1,j),' ', current%acc_err(1,j),' '
          else
             write(idev,"(I2,A,X,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)",advance='no') i,' ', &
                  j, ' ',current%low+current%step*(j-1),' ',current%acc_res(1,j),' ', current%acc_err(1,j),' '
             do k = 2,nweights-1
                write(idev,"(A,X,2X,1PE23.16,A,2X,1PE23.16,A)",advance='no') ' ',current%acc_res(k,j),' ', current%acc_err(k,j),' '
             enddo
             write(idev,"(A,X,2X,1PE23.16,A,2X,1PE23.16,A)",advance='yes') ' ',current%acc_res(nweights,j),' ', current%acc_err(nweights,j),' '
          endif

       enddo
      
    enddo

    if (writetofile) close(idev)

    call check_sumbin()

    !-- reset required variables
    evalcount = 0 !-- reset counter
    call reset_iter_proc() !-- reset everything user may want to reset
  
  end subroutine prepare_and_write

  subroutine check_sumbin()
    use mod_parser
    use mod_parms
    use mod_vegas_parms
    use mod_proc_parms
    !
    integer :: i,j
    real(dp) :: sumbin,sumerr

    if (nohistos) return

    current => head

    do i = 1,n_histo
       !
       current => current%next_histo

       sumbin = 0
       sumerr = 0

       do j = 1,current%nbin+1 

          sumbin = sumbin + current%acc_res(1,j)
          sumerr = sumerr + sqrt(current%acc_err(1,j)**2)
          
       enddo
       !
       write(6,*) '# Histo: ', current%name, ' Σ  -> ', sumbin, 'pm', sumerr
       !
    enddo
    
  end subroutine check_sumbin

#if ( _MPIbuild )
  subroutine mpi_gather_histos()
    use mpi
    implicit none
    real(dp), allocatable :: res(:,:), err(:,:)
    integer :: i, nbins, evalcount_tot
    real :: start, finish
    
    call cpu_time(start)

    !-- see comment below for its usage
    call mpi_reduce( evalcount, evalcount_tot, 1, mpi_integer, mpi_sum, mpi_root_process, mpi_comm_world, ierr )
    call mpi_proc_reduce()
    if (root_process) then
       evalcount = evalcount_tot
    else
       evalcount = 0 !-- evalcount set to 0 only for root_process in the main prepare_and_write
    endif
    
    current => head
    
    do i = 1,n_histo
       !
       current => current%next_histo

       nbins = current%nbin + 2

       allocate(res(nweights,nbins))
       allocate(err(nweights,nbins))
          
       !-- take current%res, send it to mpi_root_process doing the operation mpi_sum in mpi_double precission
       !-- only moves memory around, puts it into somthing of size nbins, which is res in 2nd position
       !-- res and err do not need to be defined here, they only need to be defined for root process
       call mpi_reduce( current%res, res, nbins*nweights, mpi_double_precision, mpi_sum, mpi_root_process, mpi_comm_world, ierr )
       call mpi_reduce( current%err, err, nbins*nweights, mpi_double_precision, mpi_sum, mpi_root_process, mpi_comm_world, ierr )

       if (root_process) then
          current%res = res
          current%err = err
       else
          current%res = zero !-- see comment above for evalcount
          current%err = zero
       endif

       deallocate( res )
       deallocate( err )

    enddo

    call cpu_time(finish)
    mpi_overhead = mpi_overhead + (finish - start)
    
    return
  end subroutine mpi_gather_histos
#endif
  
end module mod_histo

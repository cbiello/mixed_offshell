! This module contains the vegas lepage algortihm rewritten in Fortran 90
! It is heavily inspired by the Fortran 77 version from Numerical Recipies
! konstantin.asteriadis@t-online.de
! 
! - it is design to work together with the ramdom number generators in
!   mod_random.f90
! - main difference to Numercial Recipies version is that it can be compiled
!   using openMPI
! - MPI reference http://condor.cc.ku.edu/~grobe/docs/intro-MPI.shtml

module mod_vegas_mpi
    use mod_types, except => vp
    use mod_random
    use mod_vegas_parms
    use mod_vegas_user
    use mod_vegas_common
    use mod_histo
    
    use mod_mpi_common
#if ( _MPIbuild )
    use mpi
#endif
    
    implicit none
    
    private
    
    public :: vegas, vegas_integrate
    
!     integer, private, parameter :: dp = selected_real_kind(15)
    
!-- public vegas parameter
!     integer, parameter, public :: mxdim = 32
!     integer, parameter, public :: mxcomp = 16
    integer, parameter :: ncomp = 1

!-- private vegas parameter
!     integer, private, parameter :: ndmx = 50
    real(vp), private, parameter :: alph = 1.5_vp
    integer, parameter :: mds = 1 !-- extra care with our tweaks of vega if we use this with mds=0 and grid
    
    !-- public control parameters
    integer, public :: ncall, itmx, nprn, ndim
    logical, public :: readin, writeout, read_rand
    character(len=72), public :: ingridfile, outgridfile
    real(vp), public :: xl(mxdim), xu(mxdim), acc
    
!-- vegas internal variables that need to be saved (e.g. grid, ...)
!     real(vp), public :: grid(50,mxdim), si, si2, swgt, schi
!     integer, private :: ndo, iteration
    real(vp), private :: max_weight = 0.0_vp 
    
    !-- parameter needed by worker
    type, public :: worker_parameter
        integer :: ndim, nd, ng, npg, modus, userdata
        real(vp) :: xl(mxdim), xu(mxdim), xjacc
        real(vp) :: dxg, grid(50,mxdim), xnd, dx(mxdim)
    end type worker_parameter
    
    !-- output of worker
    type, private :: worker_output
        real(vp) :: ti, tsi
        real(vp) :: fb(3), f2b(3)
        real(vp) :: d(50,mxdim), di(50,mxdim)
        integer :: ia(mxdim)
        logical :: empty(3)
    end type worker_output
    
contains

  !-- wrapper for go_vegas_gen and grids
  subroutine vegas_integrate(ndim_in,integrand,vg_result,vg_error,vg_chi2)
    integer, intent(in) :: ndim_in
    real(dp15), intent(out) :: vg_result,vg_error,vg_chi2
    interface
       function integrand(xx,ff,weight)
         use mod_types
         integer    :: integrand
         real(dp15) :: xx(30),ff(1),weight
       end function integrand
    end interface

    readin    = readgrid
    writeout  = writegrid
    read_rand = .false.
    
    if (writegrid) then
       write (6,*) '# writing out vegas grid'
       call go_vegas_gen(ndim_in,integrand,vg_result,vg_error,vg_chi2,grid_out=gridfile)
    elseif (readgrid) then
       write (6,*) '# reading in vegas grid'
       call go_vegas_gen(ndim_in,integrand,vg_result,vg_error,vg_chi2,grid_in=gridfile)
    else
       call go_vegas_gen(ndim_in,integrand,vg_result,vg_error,vg_chi2)
    endif
    
  end subroutine vegas_integrate

  subroutine go_vegas_gen(ndim_in,integrand,res,err,chi2,grid_in,grid_out,io_grid)
    integer, intent(in)   :: ndim_in
    real(vp), intent(out) :: res,err,chi2
    integer :: iogrid
    character(*), intent(in), optional :: grid_in,grid_out
    integer, intent(in), optional      :: io_grid
    integer, parameter :: verbose = 0
    !
    interface
       function integrand(xx,ff,weight)
         use mod_types
         integer    :: integrand
         real(dp15) :: xx(30),ff(1),weight
       end function integrand
    end interface
    
    if (present(io_grid)) then
       iogrid = io_grid
    else
       iogrid = 66
    endif

    call init_rand(seed)

    ndim = ndim_in
    if (verbose.eq.0) then
       nprn = 1
    else
       nprn = 0
    endif

    xl = 0._vp
    xu = 1._vp
    acc = -1._vp
    
    if (present(grid_in)) then

       !-- read from the grid
       nohistos = .false.
       call clear_histo()
       !
       !call my_vegas_lepage(integrand,idum,1,ndim,VegasNc1,VegasIt1,verbose,res,err,chi2,grid_in=grid_in,io_grid=iogrid)
       itmx = VegasIt1
       ncall = VegasNc1
       ingridfile = grid_in
       call vegas(1,integrand,res,err,chi2)

    else

       !-- warmup
       nohistos = .true.
       call clear_histo()
       !
       itmx = VegasIt0
       ncall = VegasNc0
       
       if(present(grid_out)) outgridfile = grid_out
       !call my_vegas_lepage(integrand,idum,0,ndim,VegasNc0,VegasIt0,verbose,res,err,chi2,grid_out=grid_out,io_grid=iogrid)
       call vegas(0,integrand,res,err,chi2)

       !-- full run
       call init_rand(seed) !-- re-initialise random number to have identical results when reading-in grid
       nohistos = .false.
       call clear_histo()
       !
       itmx = VegasIt1
       ncall = VegasNc1
       call vegas(1,integrand,res,err,chi2)
       !call my_vegas_lepage(integrand,idum,1,ndim,VegasNc1,VegasIt1,verbose,res,err,chi2)
       
    endif

    return
    
  end subroutine go_vegas_gen
  
  
    subroutine vegas( init, integrand, avgi, sd, chi2a )
      implicit none
        
        interface
            function integrand(x,f,wgt)
                use mod_types
                integer :: integrand
                real(dp15) :: x(30), f(1), wgt
            end function integrand
        end interface
        
        integer, intent(in) :: init
        real(vp), intent(out) :: avgi, sd, chi2a
        
!         integer, intent(in) :: ndim
        
        !-- local variables
        integer :: nd, ng, k, ndm, npg
        real(vp) :: dx(mxdim), calls, xjacc, xnd, weight
        real(vp) :: dxg, dv2g, xin(50), dr, xn, rc, xo, r(50)
        real(vp) :: dt(mxdim), d(50,mxdim), di(50,mxdim)
        real(vp) :: ti, ti2, tsi
        integer :: it
        
        !-- container for vegas_worker
        type(worker_parameter) :: param
        type(worker_output) :: res

        !-- for random number generator
        integer :: idum, iy, iv(32)
        
        !-- for do loops
        integer :: i, j
        
        real :: efficiency
        real :: it_start, it_finish, it_time, single_core_esstimate, total_worker_time, best_case
        
        integer :: mpi_calls, mpi_start
#if ( _MPIbuild )
        type(worker_output), allocatable :: mpi_res(:)
        integer :: mpi_rest, mpi_worker_output
#endif
        
        !-- check input parameter
        if ( ndim .gt. mxdim ) goto 101     
        
    !--------------------------------------------------
    ! MPI initialisation
    !--------------------------------------------------
        
#if ( _MPIbuild )
        call mpi_create_worker_output_type( mpi_worker_output )
        if ( root_process ) then
            allocate( mpi_res(mpi_size) )
        end if
#endif
        
    !--------------------------------------------------
    !  VEGAS initialisation
    !--------------------------------------------------

        modus = mds
        
        if ( init .le. 0 ) then
            ndo = 1
            grid(1,1:ndim) = 1.0_vp
        end if
        
        !-- entry of vegas1_Lepage
        !-- initialises  cumulative  variables but not grid
        
        if ( init .le. 1 ) then
            iteration = 0
            si = 0.0_vp
            si2 = 0.0_vp
            swgt = 0.0_vp
            schi = 0.0_vp
        end if
        
        !-- entry of vegas2_Lepage
        !-- no initialisation at all
        
        if ( init .le. 2 ) then
            nd = ndmx
            ng = 1
            
            if ( modus .ne. 0 ) then
                ng = int( (dble(ncall) / 2d0)**(1d0 / dble(ndim)) )
                modus = 1
                if ( ( 2 * ng - ndmx ) .ge. 0 ) then
                    modus = -1
                    npg = ng / ndmx + 1
                    nd = ng / npg
                    ng = npg * nd
                end if
            end if
            
            k = ng**ndim
            npg = ncall / k
            if ( npg .lt. 2 ) npg = 2
            calls = dble(npg*k)
            dxg = 1.0_vp / ng

            dv2g = (calls * dxg**ndim)**2 
            dv2g = dv2g / dble(npg) / dble(npg) / dble( npg - 1.0_vp )

            xnd = dble(nd)
            ndm = nd - 1
            dxg = dxg * xnd
            
            !-- compute integration volume and jaccobian
            dx(1:ndim) = xu(1:ndim) - xl(1:ndim) 
            xjacc = product( dx(1:ndim) ) / calls
        end if

        if (init .gt. 2) then
           print *, 'wrong vegas initialisation'
           nd  = 0
           npg = 0
           ng  = 0
           xnd   = 0._vp
           xjacc = 0._vp
           dxg   = 0._vp
           stop
        endif
        
    !--------------------------------------------------
    ! read-in grid
    !--------------------------------------------------
        
        if ( readin) then
            if ( root_process ) then
                write(6,'(/,a,a)') 'load grid from ', ingridfile
                call flush(6)
            end if
            
            open( unit = 11, file = ingridfile, status = 'unknown')
            
            !-- read status of random number generator in buffer
            !-- and initialize if requested
            read(11,204) idum, iy, iv
            if ( read_rand ) then
                call init_rand( idum, iy, iv )
                write(6,'(/,a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                write(6,'(a)')   '!!!  WARNING: RANDOM NUMBER GENERATOR SET TO STATUS FROM GRID FILE  !!!'
                write(6,'(a)')   '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            end if
            
            !-- read grid
            do j = 1, ndim, 1
                read(11,203) ( grid(i,j), i = 1, nd )
            end do
            
            close(11)
            
            ndo = nd
            readin = .false.
        endif
        
    !--------------------------------------------------
    ! rebin preserving bin density
    !--------------------------------------------------
        
        if ( nd .ne. ndo ) then
            rc = ndo / xnd
            call rebin( nd, grid, rc )
            ndo = nd
        end if
        
    !--------------------------------------------------
    ! potentially print out run parameter and info
    !--------------------------------------------------
        
        if ( root_process ) then
            if ( nprn .ge. 0 ) write(6,200) ndim, calls, iteration, itmx, acc, modus, nd, (xl(j), xu(j), j=1,ndim)
            call flush(6)
        end if

    !--------------------------------------------------
    ! save parameter
    !--------------------------------------------------
        
        param%ndim = ndim
        param%nd = nd
        param%ng = ng
        param%npg = npg
        param%grid = grid
        param%modus = modus
        param%userdata = 0
            
        param%xl = xl
        param%xu = xu
        param%xjacc = xjacc
        param%dxg = dxg
        param%xnd = xnd
        param%dx = dx
        
    !--------------------------------------------------
    ! MPI: split integration into mpi_nproc parts
    !--------------------------------------------------
    
#if ( _MPIbuild )
        mpi_calls = int(calls) / mpi_size
        mpi_rest = mod(int(calls), mpi_size)
        mpi_start = mpi_rank * mpi_calls + min( mpi_rank, mpi_rest ) + 1
        if ( mpi_rank .lt. mpi_rest ) mpi_calls = mpi_calls + 1
#else
        mpi_calls = int(calls)
        mpi_start = 1
#endif
        
    !--------------------------------------------------
    ! main integration loop
    ! - entry of vegas3_Lepage (not implemented)
    !--------------------------------------------------
        
        do it = iteration + 1, itmx, 1
            call cpu_time( it_start )
        
            !-- update grid in parameter container
            param%grid(:,1:ndim) = grid(:,1:ndim)
        
#if ( _MPIbuild )
            !-- broadcast new grid to all processes
            call mpi_bcast_grid( param%grid, ndim )
            !-- "initialize random number generator"
            call random_offset( (mpi_start - 1) * ndim )
#endif
            
            !-- run wegas worker
            call vegas_worker(integrand, mpi_start, mpi_calls, param, res )
            
#if ( _MPIbuild )
            !-- "finalize random number generator"
            call random_offset( (int(calls) - (mpi_start - 1) - mpi_calls) * ndim )
            !-- collect results
            call mpi_reduce_results( res, param, mpi_worker_output )
#endif
            
            if ( root_process ) then
                !-- save results from worker
                ti = res%ti
                tsi = res%tsi
                
                d(1:nd,1:ndim) = res%d(1:nd,1:ndim)
                di(1:nd,1:ndim) = res%di(1:nd,1:ndim)
                
                ! compute final results for this iteration

                tsi = tsi * dv2g
                ti2 = ti**2
                
                weight = ti2 / tsi
                si = si + ti * weight
                si2 = si2 + ti2
                swgt = swgt + weight
                schi = schi + ti2 * weight

                avgi = si / swgt
                sd = swgt * dble(it) / si2
                chi2a = sd * (schi / swgt - avgi * avgi) / (dble(it) - .999d0)

                sd = dsqrt(1d0/ sd)
                
                ! eventually print final result of this iteration
                if ( nprn .ne. 0 ) then
                    tsi = dsqrt(tsi)
                
                    ! iteration: current iteration, ti: integral, tsi: std.dev.
                    ! avgi: accum. intgreal, sd: accum. std.dev., max_weight: max weight, chi2a: chi^2
                    write(6,201) it, ti, avgi, tsi, sd, max_weight, chi2a
                    write(15,'(2X,I3,4E20.8,I0)') it, ti, tsi, avgi, sd
                    call flush(6)
                    call flush(15)

                    if( nprn .lt. 0) then
                        do j = 1, ndim
                            write(6,202) j, (grid(i,j), di(i,j), d(i,j), i=1,nd)
                        end do
                    end if
                end if
        
                ! refine grid
                do j = 1, ndim
                    xo = d(1,j)
                    xn = d(2,j)
                    d(1,j) = (xo + xn) / 2.0_vp
                    dt(j) = d(1,j)
                    
                    do i = 2, ndm
                        d(i,j) = xo + xn
                        xo = xn
                        xn = d(i+1,j)
                        d(i,j) = (d(i,j) + xn) / 3.0_vp
                        dt(j) = dt(j) + d(i,j)
                    end do
                    
                    d(nd,j) = (xn + xo) / 2.0_vp
                    dt(j) = dt(j) + d(nd,j)
                end do

                do j = 1, ndim
                    rc = 0.0_vp
                    
                    do i = 1, nd
                        r(i) = 0.0_vp
                        if ( d(i,j) .gt. 0.0_vp ) then
                            xo = dt(j) / d(i,j)
                            r(i) = ((xo-1.0_vp) / xo / dlog(xo))**alph
                        end if
                        rc = rc + r(i)
                    end do
                    
                    rc = rc / xnd
                    
                    k = 0
                    xn = 0.0_vp
                    dr = xn
                    i = k
                    
25                  k=k+1
                    dr=dr+r(k)
                    xo=xn
                    xn=grid(k,j)
26                  if(rc.gt.dr) goto 25
        
                    i=i+1
                    dr=dr-rc
                    xin(i)=xn-(xn-xo)*dr/r(k)
                    if(i.lt.ndm) goto 26
                
                    grid(1:ndm,j) = xin(1:ndm)
                    grid(nd,j) = 1.0_vp
                 end do

            end if
            
            !-- call user function after iteration
            call vegas_user(init,it,ti,avgi,tsi,sd,chi2a)
#if ( _MPIbuild )
            call mpi_bcast_real( avgi )
#endif
            
            !-- some statistics about MPI efficiency
            call cpu_time( it_finish )
            it_time = it_finish - it_start
#if ( _MPIbuild )
            call mpi_reduce( mpi_worker_time, total_worker_time, 1, mpi_real, mpi_sum, mpi_root_process, mpi_comm_world, ierr ) ! put into container
            
            if ( root_process ) then
                single_core_esstimate = (it_time - mpi_overhead - mpi_worker_time) + total_worker_time ! esstimated time for single core
                best_case = single_core_esstimate / mpi_size
                efficiency = best_case / it_time * 100.0
                write(*,'(/,a,f12.3,a)') ' Wall time:      ', it_time / 60.0, ' minutes'
                write(*,'(a,f12.3,a)')   ' 1CPU esstimate: ', single_core_esstimate / 60.0, ' minutes'
                write(*,'(a,f12.3,a)')   ' MPI overhead:   ', mpi_overhead / 60.0, ' minutes'
                write(*,'(a,f12.3,a,/)') ' MPI efficiency: ', efficiency, ' %'
                mpi_overhead = 0.0
            end if
#endif
            
            !-- stop earlier if desired precision goal is already reached
            if ( acc * dabs(avgi) .ge. sd ) exit 
        end do

    !--------------------------------------------------
    !    write-out grid if necessary
    !--------------------------------------------------
    
        if ( root_process .and. writeout ) then
            open( unit = 11, file = outgridfile, status = 'unknown' )
            write(6,'(/,a,a)') 'write grid to ', outgridfile
            call flush(6)
            call get_rand_stat( idum, iy, iv )
            write(11,204) idum, iy, iv
            do j = 1, ndim
                write(11,203) ( grid(i,j), i  =1, nd )
            enddo
            close(11)
        end if
        
    !--------------------------------------------------
    ! finalize
    !--------------------------------------------------
        
#if ( _MPIbuild )
        if ( root_process ) then
            deallocate( mpi_res )
        end if
#endif
        
        return
        
    !--------------------------------------------------
    ! error messages
    !--------------------------------------------------
        
101     continue
        write(*,'(/,a,i0,a,i0/)') '[error][vegas] ndim = ', ndim, ' greater then mxdim = ', mxdim
        stop
        
    !--------------------------------------------------
    ! write formats
    !--------------------------------------------------
        
200     format(/ 1X,' Input parameters for vegas:  ndim=',i3, &
        '   ncall=',f10.0/28x,'  it=',i5,'    itmx=',i5/28x, &
        '  acc=',g9.3/28x,'  modus=',i3,'     nd=',i4/28x, &
        '  (xl,xu)=',(t40,'( ',g12.6,' , ',g12.6,' )'))
201     format(/'************* Integration by Vegas (iteration ',i3, &
        ') **************' / '*',63x,'*'/, &
        '*  integral  = ',g14.8,2x, &
        ' accum. integral = ',g14.8,'*'/, &
        '*  std. dev. = ',g14.8,2x, &
        ' accum. std. dev = ',g14.8,'*'/, &
        '*   max. wt. = ',g14.6,35x,'*'/,'*',63x,'*'/, &
        '**************   chi**2/iteration = ', &
        g10.4,'   ****************' /)
202     format(1X,' data for axis',i2,/,' ',6x,'x',7x,'  delt i ', &
        2x,'conv','ce   ',11x,'x',7x,'  delt i ',2x,'conv','ce  ' &
        ,11x,'x',7x,'   delt i ',2x,'conv','CE  ',/, &
        (1X,' ',3g12.4,5x,3g12.4,5x,3g12.4))
203     format(/(5z16))
204     format(/(5z16.16))
    end subroutine vegas

    subroutine rebin( nd, grid, rc )
        implicit none
        
        integer, intent(in) :: nd
        real(vp), intent(inout) :: grid(50,mxdim)
        real(vp), intent(in) :: rc
        
        integer :: j, i, k
        real(vp) :: dr, xo, xn, xin(50)

        do j = 1, ndim, 1
            k = 1
            dr = 1.0_vp
            xo = 0.0_vp
            xn = grid(k,j)
            
            do i = 1, nd - 1, 1
                do while ( rc .gt. dr )
                    k = k + 1
                    dr = dr + 1.0_vp
                    xo = xn
                    xn = grid(k,j)
                end do

                dr = dr - rc
                xin(i) = xn - (xn - xo) * dr
            end do

            grid(1:nd-1,j) = xin(1:nd-1)
            grid(nd,j) = 1.0_vp
        end do
        
        return
    end subroutine rebin
    
!----------------------------------------------------------------------------------
!   MPI send and recieve
!----------------------------------------------------------------------------------
    
#if ( _MPIbuild )
    subroutine mpi_bcast_grid( grid, ndim )
        implicit none
        real(vp), intent(inout) :: grid(50,mxdim)
        integer, intent(in) :: ndim
        real :: start, finish
        call cpu_time(start)
        call mpi_bcast( grid, 50 * ndim, mpi_double_precision, 0, mpi_comm_world, ierr )
        call cpu_time(finish)
        mpi_overhead = mpi_overhead + (finish - start)
        return
    end subroutine mpi_bcast_grid

    subroutine mpi_bcast_real( avgi )
        implicit none
        real(vp), intent(in) :: avgi
        real :: start, finish
        call cpu_time(start)
        call mpi_bcast( avgi, 1, mpi_double_precision, 0, mpi_comm_world, ierr )
        call cpu_time(finish)
        mpi_overhead = mpi_overhead + (finish - start)
        return
    end subroutine mpi_bcast_real
    
    ! subroutine mpi_send_result( res )
    !     implicit none
    !     
    !     type(worker_output), intent(in) :: res
    !     
    !     call mpi_send( res%ti, 1, mpi_double_precision, mpi_root_process, 2000, mpi_comm_world, ierr)
    !     call mpi_send( res%tsi, 1, mpi_double_precision, mpi_root_process, 2000, mpi_comm_world, ierr)
    !     call mpi_send( res%fb, 3, mpi_double_precision, mpi_root_process, 2000, mpi_comm_world, ierr)
    !     call mpi_send( res%f2b, 3, mpi_double_precision, mpi_root_process, 2000, mpi_comm_world, ierr)
    !     call mpi_send( res%ia, mxdim, mpi_integer, mpi_root_process, 2000, mpi_comm_world, ierr)
    !     call mpi_send( res%empty, 3, mpi_logical, mpi_root_process, 2000, mpi_comm_world, ierr)
    !     call mpi_send( res%d, 50*ndim, mpi_double_precision, mpi_root_process, 2000, mpi_comm_world, ierr)
    !     call mpi_send( res%di, 50*ndim, mpi_double_precision, mpi_root_process, 2000, mpi_comm_world, ierr)
    !     
    !     return
    ! end subroutine mpi_send_result
    
    ! subroutine mpi_receive_result( id, res )
    !     implicit none
        
    !     integer, intent(in) :: id
    !     type(worker_output), intent(out) :: res
    !     integer :: n
    !     real :: start, finish
        
    !     call cpu_time(start)
    !     call mpi_recv( res%ti, 1, mpi_double_precision, id, 2000, mpi_comm_world, mpi_stat, ierr)
    !     call mpi_recv( res%tsi, 1, mpi_double_precision, id, 2000, mpi_comm_world, mpi_stat, ierr)
    !     call mpi_recv( res%fb, 3, mpi_double_precision, id, 2000, mpi_comm_world, mpi_stat, ierr)
    !     call mpi_recv( res%f2b, 3, mpi_double_precision, id, 2000, mpi_comm_world, mpi_stat, ierr)
    !     call mpi_recv( res%ia, mxdim, mpi_integer, id, 2000, mpi_comm_world, mpi_stat, ierr)
    !     call mpi_recv( res%empty, 3, mpi_logical, id, 2000, mpi_comm_world, mpi_stat, ierr)
    !     call mpi_recv( res%d, 50*ndim, mpi_double_precision, id, 2000, mpi_comm_world, mpi_stat, ierr)
    !     call mpi_recv( res%di, 50*ndim, mpi_double_precision, id, 2000, mpi_comm_world, mpi_stat, ierr)
    !     call cpu_time(finish)
        
    !     mpi_overhead = mpi_overhead + (finish - start)
        
    !     return
    ! end subroutine mpi_receive_result
    
    subroutine mpi_create_worker_output_type( newtype )
        implicit none
        
!     !-- output of worker
!     type, private :: worker_output
!         real(vp) :: ti, tsi
!         real(vp) :: fb(3), f2b(3)
!         real(vp) :: d(50,mxdim), di(50,mxdim)
!         integer :: ia(mxdim)
!         logical :: empty(3)
!     end type worker_output
        
        integer, intent(out) :: newtype
        type(worker_output) :: dummy
        integer, parameter :: entries = 8
        integer :: n, types(entries), blocklengths(entries)
        integer(mpi_address_kind) :: base, displacements(entries)
        
        !-- variable types
        types(1:6) = mpi_double_precision
        types(7) = mpi_integer
        types(8) = mpi_logical
        
        !-- blocklengths
        blocklengths(1:2) = 1
        blocklengths(3:4) = 3
        blocklengths(5:6) = 50 * mxdim
        blocklengths(7) = mxdim
        blocklengths(8) = 3
        
        !-- displacements
        call mpi_get_address( dummy, base, ierr )
        call mpi_get_address( dummy%ti, displacements(1), ierr )
        call mpi_get_address( dummy%tsi, displacements(2), ierr )
        call mpi_get_address( dummy%fb, displacements(3), ierr )
        call mpi_get_address( dummy%f2b, displacements(4), ierr )
        call mpi_get_address( dummy%d, displacements(5), ierr )
        call mpi_get_address( dummy%di, displacements(6), ierr )
        call mpi_get_address( dummy%ia, displacements(7), ierr )
        call mpi_get_address( dummy%empty, displacements(8), ierr )
        
        do n = 1, entries, 1
            displacements(n) = mpi_aint_diff( displacements(n), base )
        end do
        
        !-- create datatype
        call mpi_type_create_struct(entries, blocklengths, displacements, types, newtype, ierr)
        call mpi_type_commit( newtype, ierr )
        
        return
    end subroutine mpi_create_worker_output_type
    
    subroutine mpi_reduce_results( res, param, mpi_worker_output )
        implicit none
        
        integer, intent(in) :: mpi_worker_output
        type(worker_parameter), intent(in) :: param
        type(worker_output), intent(inout) :: res
        type(worker_output), allocatable :: res_all(:)
        real :: start, finish
        
        call cpu_time(start)
        
        allocate( res_all(mpi_size) )
        call mpi_gather( res, 1, mpi_worker_output, res_all, 1, mpi_worker_output, mpi_root_process, mpi_comm_world, ierr )
        if ( root_process ) call combine( param, res_all, res )
        deallocate( res_all )
        
        call cpu_time(finish)
        mpi_overhead = mpi_overhead + (finish - start)
        
        return
    end subroutine mpi_reduce_results
    
!             if ( .not. root_process ) then
!                 call mpi_send_result( res )
!             else
!                 mpi_res(1) =  res
!                 do j = 1, mpi_size - 1, 1
!                     call mpi_recieve_result( j, mpi_res(j+1) )
!                 end do
!                 call combine( param, mpi_res, res )
!             end if
            

    
    
!     subroutine mpi_create_worker_output_op( newop )
!         implicit none
!         
!         integer, intent(out) :: newop
!         
!         !-- create operation to combine consecutive newtype objects
!         call mpi_op_create( combine_worker_output, .false., newop, ierr )
!         
!         return
!     end subroutine mpi_create_worker_output_op
!     
!     subroutine combine_worker_output( in, inout, length, datatype )
!         implicit none
!         
!         integer, intent(in) :: length, datatype
!         type(worker_output), intent(in) :: in(length)
!         type(worker_output), intent(inout) :: inout(length)
!         
!         return
!     end subroutine combine_worker_output
#endif

!  SUBROUTINE MPI_User_function(invec, inoutvec, len, datatype)
! USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
! TYPE(C_PTR), VALUE :: invec, inoutvec
! INTEGER :: len
! TYPE(MPI_Datatype) :: datatype
! SUBROUTINE USER_FUNCTION(INVEC, INOUTVEC, LEN, DATATYPE)
! <type> INVEC(LEN), INOUTVEC(LEN)
! INTEGER LEN, DATATYPE 


!----------------------------------------------------------------------------------
!   vegas_worker(...)
!   performs vegas main integration loop from kg_start(ndim) vector to 
!   kg_stop(ndim) vector all parameters are passed via worker_parameter type and 
!   results via worker_output type
!----------------------------------------------------------------------------------

    subroutine vegas_worker(integrand, offset, ncalls, param, res )
        implicit none
        
        interface
            function integrand(x,f,wgt)
                use mod_types
                integer :: integrand
                real(dp15) :: x(30), f(1), wgt
            end function integrand
        end interface
        
        integer, intent(in) :: offset
        integer, intent(in) :: ncalls
        type(worker_parameter), intent(in) :: param
        type(worker_output), intent(out) :: res
        
        integer :: k, j, ia(mxdim), ffdum
        real(vp) :: f, fb, f2b, weight, ff(1), f2
        real(vp) :: x, xn, xx(mxdim), rc, xo
        integer :: nd, ndim, kg(mxdim), k_min, ng, npg, todo
        logical :: lcomp, rcomp
        real :: start, finish
        
        call cpu_time(start)
        
        ng = param%ng
        npg = param%npg
        nd = param%nd
        ndim = param%ndim
        todo = ncalls
        
        res%fb = 0.0_vp
        res%f2b = 0.0_vp
        res%ti = 0.0_vp
        res%tsi = 0.0_vp
        res%d(1:nd,1:ndim) = 0.0_vp
        res%di(1:nd,1:ndim) = 0.0_vp
        res%empty = .true.
        
        !-- compute initial kg and k_min for first call
        call int_to_array( ng, (offset - 1) / npg, kg(1:ndim) )
        k_min = 1 + mod(offset - 1, npg)
        lcomp = .true.
        if ( k_min .ne. 1 ) lcomp = .false.

        do while ( .true.  )
            fb = 0.0_vp
            f2b = 0.0_vp
            rcomp = .false.
            
            do k = k_min, param%npg
                weight = param%xjacc
                
                do j = 1, ndim, 1
                    call random_real( x )
                    xn = (dble(kg(j)) - x) * param%dxg + 1.0_vp
                    ia(j) = int(xn)
                
                    if ( ia(j) .gt. 1 ) then
                        xO = param%grid(ia(j),j) - param%grid(ia(j) - 1,j)
                        rc = param%grid(ia(j) - 1,j) + (xn - dble(ia(j))) * xo
                    else 
                        xo = param%grid(ia(j),j)
                        rc = (xn - dble(ia(j))) * xo
                    end if
                    
                    xx(j) = param%xl(j) + rc * param%dx(j)
                    weight = weight * xo * param%xnd
                end do

                f = weight

                ffdum = integrand(xx,ff,weight)
                f = f * ff(1)

                f2 = f**2
                fb = fb + f
                f2b = f2b + f2
                
                do j = 1, ndim
                    res%di(ia(j),j) = res%di(ia(j),j) + f
                    if ( param%modus .ge. 0 ) res%d(ia(j),j) = res%d(ia(j),j) + f2
                end do
                
                todo = todo - 1
                if ( k .eq. param%npg ) rcomp = .true.
                if ( todo .eq. 0 ) exit
            end do

            if ( lcomp .and. rcomp ) then
                f2b = sqrt( f2b * param%npg )
                f2b = (f2b - fb) * (f2b + fb)
            
                res%ti = res%ti + fb
                res%tsi = res%tsi + f2b
                
                if ( param%modus .lt. 0 ) then
                    do j = 1, ndim
                        res%d(ia(j),j) = res%d(ia(j),j) + f2b
                    end do
                end if
            else if ( (.not. lcomp) .and. rcomp ) then
                res%fb(1) = fb
                res%f2b(1) = f2b
                res%ia = ia
                res%empty(1) = .false.
            else if ( lcomp .and. (.not. rcomp) ) then
                res%fb(2) = fb
                res%f2b(2) = f2b
                res%empty(2) = .false.
            else
                res%fb(3) = fb
                res%f2b(3) = f2b
                res%empty(3) = .false.
            end if
            
            k_min = 1
            lcomp = .true.
            if ( todo .eq. 0 ) exit
            call increment( param%ng, kg(1:ndim) )
        end do

        call cpu_time(finish)
        mpi_worker_time = finish - start
        
        return
    end subroutine vegas_worker
    
    subroutine combine( param, results, res )
        implicit none
        
        type(worker_parameter), intent(in) :: param
        type(worker_output), intent(in) :: results(:)
        type(worker_output), intent(out) :: res
        integer :: n
        
        res = results(1)
        do n = 2, size(results,1),  1
            call combine2( param, res, results(n) )
        end do
        
        return
    end subroutine combine
    
    subroutine combine2( param, res1, res2 )
        implicit none
        
        type(worker_parameter), intent(in) :: param
        type(worker_output), intent(inout) :: res1
        type(worker_output), intent(in) :: res2
        
        real(vp) :: fb, f2b
        integer :: j
        
        res1%ti = res1%ti + res2%ti
        res1%tsi = res1%tsi + res2%tsi
        res1%d = res1%d + res2%d
        res1%di = res1%di + res2%di

        if ( .not. res1%empty(2) .and. .not. res2%empty(1) ) then
            fb = res1%fb(2) + res2%fb(1)
            f2b = res1%f2b(2) + res2%f2b(1)
            
            f2b = sqrt( f2b * param%npg )
            f2b = (f2b - fb) * (f2b + fb)
        
            res1%ti = res1%ti + fb
            res1%tsi = res1%tsi + f2b

            if ( param%modus .lt. 0 ) then
                do j = 1, param%ndim
                    res1%d(res2%ia(j),j) = res1%d(res2%ia(j),j) + f2b
                end do
            end if
            
            if ( res2%empty(2) ) then
                res1%empty(2) = .true.
            else
                res1%empty(2) = .false.
                res1%fb(2) = res2%fb(2)
                res1%f2b(2) = res2%f2b(2)
            end if
        else if (all(res1%empty) .and. res2%empty(1) .and. (.not. res2%empty(2)) ) then
            res1%empty(2) = .false.
            res1%fb(2) = res2%fb(2)
            res1%f2b(2) = res2%f2b(2)
        else if ( .not. res1%empty(2) .and. .not. res2%empty(3) ) then
            res1%fb(2) = res1%fb(2) + res2%fb(3)
            res1%f2b(2) = res1%f2b(2) + res2%f2b(3)
         else if (all(res1%empty) .and. all(res2%empty)) then
            continue
         else
            print *, 'not covered in combine2',res1%empty,res2%empty
        end if
        
        return
    end subroutine combine2
    
!----------------------------------------------------------------------------------
!   functions to manipulate numbers in the format (/ dn, ..., d2, d1 /)
!   where dn are individual digits to a given base; note that digits, start
!   at 1 not 0; hence in e.g. binary (base = 2) 101101 is (/ 2, 1, 2, 2, 1, 2 /)
!----------------------------------------------------------------------------------
    
    subroutine increment( base, number )
        implicit none
        
        integer, intent(in) :: base
        integer, intent(inout) :: number(:)
        integer :: n
        
        n = size(number,1)
        do n = size(number,1), 1, -1
            number(n) = number(n) + 1
            if ( number(n) .le. base ) return
            number(n) = 1
        end do
        
        return
    end subroutine increment

    subroutine int_to_array( base, i, number )
        implicit none
        
        integer, intent(in) :: base, i
        integer, intent(out) :: number(:)
        integer :: n, ii
        
        ii = i
        number = 1
        do n = size(number,1), 1, -1
            number(n) = mod(ii, base) + 1
            ii = ii / base
            if ( ii .eq. 0 ) exit
        end do
        
        return
    end subroutine int_to_array
    
!     subroutine time_string( seconds, time_string )
!         implicit none
!         
!         real, intent(in) :: seconds
!         character(len=*), intent(out) :: time_string
!         
!         real :: sec
!         integer :: minute, hour
!         character(len=100) :: buffer
!         
!         
!         return
!     end subroutine time_string
    
  end module mod_vegas_mpi

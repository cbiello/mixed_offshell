! random number generator
! (C) Copr. 1986-92 Numerical Recipes Software ]2w.1,r1..
! rewritten in modern fortran

module mod_random
    use mod_mpi_common
    
    implicit none
    
    private
    
    integer, private, parameter :: dp = selected_real_kind(15)
    
    integer, private, parameter :: ia = 16807
    integer, private, parameter :: im = 2147483647
    integer, private, parameter :: iq = 127773
    integer, private, parameter :: ir = 2836
    
    integer,  public, parameter :: ntab = 32
    integer, private, parameter :: ndiv = int(1 + (im - 1._dp) / ntab)
    
    real(dp), private, parameter :: am = 1.0_dp / im
    real(dp), private, parameter :: eps = 3.0E-16_dp
    real(dp), private, parameter :: rnmx = 1.0_dp - eps
    
    integer, private :: idum
    integer, private :: iv(ntab) = 0
    integer, private :: iy = 0
    
    public :: random_real
    public :: init_rand
    public :: get_rand_stat, random_offset
    
    logical, private :: initialized = .false.
    
    interface init_rand
        module procedure init_rand_seed
        module procedure init_rand_stat
    end interface
    
contains
    
! >--------------------------------------------------------------------------
! | ...
! >--------------------------------------------------------------------------
    
    subroutine random_real( harvest )
        implicit none
        
        real(dp), intent(out) :: harvest
        integer :: j, k, l
        
        if (.not. initialized) call init_rand(1)
        
        k = idum / iq
        l = mod(idum, iq)
        idum = ia * l - ir * k
        if (idum .lt. 0) idum = idum + im
        
        j = 1 + iy / ndiv
        iy = iv(j)
        iv(j) = idum
        
        harvest = min(am * real(iy, kind = dp), rnmx)
        
        return
    end subroutine random_real
    
! >--------------------------------------------------------------------------
! | module procedure of subroutine init_rand
! | - init_rand_seed for ...
! | - init_rand_stat for ...
! >--------------------------------------------------------------------------
    
    subroutine init_rand_seed( seed )
        implicit none
        
        integer, intent(in) :: seed 
        integer :: j, k, l

        if (seed.eq.0) then
           idum = 19
        elseif (seed .lt. 0) then
           if ( root_process ) write(*,'(/,a,/)') '[random][warning] wrong seed, seed has to be an integer greater 0, default = 19 used'
           idum = 19 !1
        else
           idum = seed
        end if
        
        do j = ntab + 8, 1, -1
          k = idum / iq
          l = mod(idum, iq)
          idum = ia * l - ir * k
          if (idum .lt. 0) idum = idum + im
          if (j .le. ntab) iv(j) = idum
        end do
        
        iy = iv(1)
        
        initialized = .true.
        
        return
    end subroutine init_rand_seed
    
    subroutine init_rand_stat( idum_in, iy_in, iv_in )
        implicit none
        
        integer, intent(in) :: idum_in, iy_in, iv_in(ntab)

        idum = idum_in
        iy = iy_in
        iv = iv_in
        
        initialized = .true.
        
        return
    end subroutine init_rand_stat
    
    subroutine random_offset( offset )
        implicit none
        
        integer, intent(in) :: offset
        real(dp) :: xx
        integer :: n
        real :: start, finish
        
        call cpu_time(start)
        
        do n = 1, offset, 1
            call random_real( xx )
        end do
        
        call cpu_time(finish)
        mpi_overhead = mpi_overhead + (finish - start)
        
        return
    end subroutine random_offset
    
! >--------------------------------------------------------------------------
! | ...
! >--------------------------------------------------------------------------
    
    subroutine get_rand_stat( idum_out, iy_out, iv_out )
        implicit none
        
        integer, intent(out) :: idum_out, iy_out, iv_out(ntab)
        
        if (.not. initialized) call init_rand(1)
        
        idum_out = idum
        iy_out = iy
        iv_out = iv
        
        return
    end subroutine get_rand_stat

end module mod_random

module mod_vegas_common
    implicit none
    public
     
    integer, parameter :: vp = selected_real_kind(15) ! internal precision
    integer, parameter :: mxdim = 32
    integer, parameter :: mxcomp = 16
    integer, parameter :: ndmx = 50

    real(vp), save :: grid(50,mxdim), si, si2, swgt, schi
    integer,  save :: ndo, iteration, modus
    
    type :: vegas_data
        character :: runtype
        integer :: iteration, ndo, modus
        real(vp) :: si, si2
        real(vp) :: swgt, schi
    end type vegas_data
    
end module mod_vegas_common

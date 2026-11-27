module mod_vegas
  use mod_types; use common_def; use mod_consts_dp
  use mod_histograms
  implicit none
  private

  public :: go_vegas_gen

contains

  subroutine go_vegas_gen(ndim_in,inseed,integrand, VG_Result, VG_Error, VG_Chi2, gridfile)
    implicit none

    interface
       function integrand(ndim,xx,ncomp,ff,userdata,weight,iternumber)
         use mod_types
         integer :: ndim, ncomp, integrand, userdata, iternumber
         real(dp15) :: xx(30), ff(1), weight
       end function integrand
    end interface

    integer, parameter :: ncompmy = 1
    integer :: ndim_in,userdata,inseed
    real(dp15) :: VG_Result,VG_Error,VG_Chi2
    character :: gridfile*(70)
    include './vegas_lepage/vegas_common.f'
    
    ndim = ndim_in
    userdata = 0

    if (inseed.eq.0) then
       idum = -19
    else
       idum = -inseed
    endif

    evalcount = 1

    !-- vegas initialization
    xl(1:mxdim) = zero
    xu(1:mxdim) = one
    acc = -one
    nprn = one
    readin = .false.
    writeout = .false.
    outgridfile(1:70) = gridfile(1:70)

    if (gridfile.ne.'') writeout = .true.

    !-- warmup
    itmx = vegasIt0
    ncall = vegasNc0
    runtype = 'w'
    call vegas_lepage(integrand,vg_result,vg_error,vg_chi2)

    !-- real run
    itmx = vegasIt1
    ncall = vegasNc1
    runtype = 'r'

    !-- reset histograms
    HistoGlobalArray(:,:) = zero 
    HistoGlobalArray_aux(:,:) = zero  
    ErrorGlobalArray(:,:) = 1d9 
    ErrorGlobalArray_aux(:,:) = 1d9
    HistoGlobalArray_Markus(:,:) = zero 

    !-- never save the grid after warming-up
    writeout = .false.
    call vegas1_lepage(integrand,vg_result,vg_error,vg_chi2)

    return

  end subroutine go_vegas_gen

end module mod_vegas
    

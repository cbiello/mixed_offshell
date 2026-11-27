module mod_chooser
  use mod_types; use common_def; use mod_consts_dp
  use mod_consts_debug
  use mod_process
  use mod_allocate_arrays
  use mod_histograms
  !
  use mod_for_chooser
  use mod_for_chooser_proc
  !
  use mod_lim_proc
  use mod_run_proc
  use mod_tc_proc
  implicit none

contains

  subroutine chooser()
    integer :: numargs,countarg,i,j,iunit,ndim
    character :: arg*(300)
    real(dp15) :: vg_result,vg_error,vg_chi2
    character :: pdffile_lo*(50),pdffile_nlo*(50),pdffile_nnlo*(50)
    character :: gridfile*(100),corr*(8),pdford*(4),sec*(10),ch*(10),runname*(70),chxmuR*(4),chxmuF*(4)
    real(dp) :: mu
#if (_DebugFeatures == 0)
    character :: globalsec*(10) !-- to compile without global variables defined
#endif
    integer :: sizeseed
    integer, allocatable :: iseed(:)

    !-- some default values
    sh = 13000.0_dp**2
    !sh = 14000.0_dp**2

    myq2min = (1.0_dp)**2
    myq2max = sh

    !-- default scale choice
    mu = 100.0_dp
    musq_in = mu**2
    xMuR = 1
    xMuF = 1

    !-- PDF files. Which one is used is determined by pdford
    pdffile_lo = 'NNPDF31_lo_as_0130.LHgrid'
#if (_luxpdf==0)
    pdffile_nlo = 'NNPDF31_nlo_as_0118.LHgrid'
    pdffile_nnlo = 'NNPDF31_nnlo_as_0118.LHgrid'
#elif (_luxpdf==1)
    pdffile_nlo = 'NNPDF31_nlo_as_0118_luxqed.LHgrid'
    pdffile_nnlo = 'NNPDF31_nnlo_as_0118_luxqed.LHgrid'
#elif (_luxpdf==2)
    pdffile_nlo = 'NNPDF31_nlo_as_0118_luxqed.LHgrid'
    pdffile_nnlo = 'smpdf20tolsmpdf_NNPDF31_nnlo_as_0118_luxqed'
#endif
    pdford='nnlo'

    !-- default buffers
    buff = 1E-9_dp
    bufftc=1E-12_dp
    buffrr=-one
    buffrv=1E-10_dp
    buffz = 1E-9_dp
    onet = one-two*buff

    histotmpfile = ''

    sec = ''

    limcheck = 10**(-six)

    !-- vegas parameters
    vegasIt0 = 5
    vegasIt1 = 10
    vegasNc0 = 100000
    vegasNc1 = 200000
    seed = 0

    runname = "test"
    corr="lo"
    ch=""

    call set_default_proc()

    !-------------------------------------------------------------------

    !-- read user-parameters from command-line
#if defined(__GFORTRAN__)
    numargs = COMMAND_ARGUMENT_COUNT () !-- gfortran
#elif defined(__INTEL_COMPILER)
    numargs = nargs() - 1 !-- ifort
#else
    print *, 'wrong compiler'
    stop
#endif

    !-------------------------------------------------------------------

    countarg = 0

    call read_parms_local(numargs,countarg,corr,pdford,sec,ch,runname)
    call read_parms_gen(numargs,countarg)
    call read_parms_proc(numargs,countarg)

    if (countarg.ne.numargs) then
       print *, 'wrong command line'
       stop
    endif

    !-------------------------------------------------------------------

    write(chxMuR,'(F4.2)') xMuR
    write(chxMuF,'(F4.2)') xMuF

    if (corr.eq.'nnlo'.and.buffrr.eq.(-one)) call set_buff_default(sec,ch)

    !-------------------------------------------------------------------

    outputfile = './histo/'//trim(runname)//"_"//trim(corr)

    if (trim(sec) .ne. '') outputfile = trim(outputfile)//"_"//trim(sec)
    if (trim(ch)  .ne. '') outputfile = trim(outputfile)//"_"//trim(ch)

    outputfile = trim(outputfile)//"_xMuR_"//trim(chxMuR)//"_xMuF_"//trim(chxMuF)//'.histo'

    !-------------------------------------------------------------------

    if (pdford.eq.'lo') then
       pdffile = pdffile_lo
    elseif (pdford.eq.'nlo') then
       pdffile = pdffile_nlo
    elseif (pdford.eq.'nnlo') then
       pdffile = pdffile_nnlo
    endif

    !-------------------------------------------------------------------

    call InitPDFset(pdffile)
    call InitPDF(0)

    call init_histo_w()
    !-- random number for lepton decay rapidity
    call random_seed(size=sizeseed)
    allocate (iseed(sizeseed))
    iseed = 0
    iseed(1) = 123456789
    call random_seed(put=iseed)
    deallocate(iseed)

    !-------------------------------------------------------------------

    gridfile = ''

    !-- printout
    iunit = 6
#include "./printout.f90"
#include "./printout_proc.f90"
    write(iunit,*)
    write(iunit,*) '#runname         ', runname
    write(iunit,*) '#histogram file  ', outputfile
    write(iunit,*) ''

    !-- allocate arrays
    if (corr.eq.'lo') then
       call allocate_arrays_lo()
    elseif (corr.eq.'nlo' .or. corr.eq.'limnlo') then
       call allocate_arrays_nlo(sec)
    elseif (corr.eq.'nnlo' .or. corr.eq.'limnnlo') then
       call allocate_arrays_nnlo(sec,ch)
       !     elseif (corr.eq.'checkamp') then
       !        print *, 'checking amplitudes'
       !     elseif (corr.eq.'tcsplit') then
       !        print *, 'no allocation required'
    elseif (corr .eq. 'tcsplit') then
       print *, "Doing triple collinear splitting -- no allocation required"
    else
       print *, 'wrong correction'
       stop
    endif

    call run_tc_proc(corr,sec,ch,gridfile,vg_result,vg_error,vg_chi2)
    call run_limits_proc(corr,sec,ch,gridfile,vg_result,vg_error,vg_chi2)
    call run_sector_proc(corr,sec,ch,gridfile,vg_result,vg_error,vg_chi2)

  end subroutine chooser

end module mod_chooser

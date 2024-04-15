module mod_chooser
  use mod_mpi_common
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_parser
  use mod_sects_list  
  use mod_vegas_parms
  use mod_hoppet_tools
  use mod_ol_interface
  !--
  use mod_histo
  use mod_cut_histo
  !--
  use mod_checks
  !--
  use mod_vegas_parms
  use mod_parms
  use mod_proc_parms
  !--
  use mod_run_proc
  !--
  implicit none
  private

  public :: chooser

contains

  subroutine chooser()
    character(30) :: input_file
    logical       :: parse_error,do_check
    integer       :: infile_dev = 77
    integer       :: outdev     = 6
    !--
    real(dp15) :: vg_result,vg_error,vg_chi2
    !
    integer :: sizeseed
    integer, allocatable :: iseed(:)
    integer :: i,k

    if (log_val_opt('-h')) call print_usage(outdev)
    if (log_val_opt('-list')) call list_corr_ch_sec(outdev)
    
    if (root_process) then
       write(6,*) '# ', adjustl(trim(command_line()))
       write(6,*) '# '
    endif

    !----------------------------------------------------------------
    !--- begin parser
    !----------------------------------------------------------------
    
    !-- check if there is an input file, and parse it
    input_file = trim(string_val_opt('-f',''))
    if (input_file .ne. '') then
       call parse_file(input_file,infile_dev,outdev)
    endif

    !-- set up all the variables
    call get_vegas_parms()
    call get_parms()
    call get_proc_parms()
    
    !-- run type, helper for these in ./Chooser
    corr = trim(string_val_opt('-corr','lo'))
    ch   = trim(string_val_opt('-ch','ns'))
    sec  = trim(string_val_opt('-sec','na'))

    !-- run checks
    do_check = log_val_opt('-check',.false.)
    
    !-- check that parsing went through correctly
    parse_error =  ta_CheckAllArgsUsed(outdev)
    if (.not. parse_error) then
       print *, 'problem with command line'
       stop
    endif

    !----------------------------------------------------------------
    !--- end parser
    !----------------------------------------------------------------
    
    !----------------------------------------------------------------
    !--- set up LHAPDF 
    !----------------------------------------------------------------
    if ( .not. root_process ) call SetLHAPARM('SILENT')
    call InitPDFset(pdfset)
    call InitPDF(pdfmem)

    !----------------------------------------------------------------
    !--- set up hoppet
    !----------------------------------------------------------------
    call init_hoppet()
    
    !----------------------------------------------------------------
    !--- set up QCD loop
    !----------------------------------------------------------------
    !call qlinit()

    !----------------------------------------------------------------
    !--- set up and printout parameters
    !----------------------------------------------------------------
    call set_outfile(pref,outputfile,gridfile)

    call set_ew_parms_cms(ew_scheme,outdev)
    call set_qcd_parms()
    call set_proc_parms()

    if (root_process) then
       write(6,*) '#'
       call print_vegas_parms(outdev)
       !
       call print_ew_parms(ew_scheme,outdev)
       call print_qcd_parms(outdev)
       !
       call print_proc_parms(outdev)
    endif

    
    
    !----------------------------------------------------------------
    !--- set up Open Loops
    !----------------------------------------------------------------
    call initialise_ol(1)
    
    !----------------------------------------------------------------
    !--- set up histograms
    !----------------------------------------------------------------
    call init_histo()
    call init_user_histo()
    call setup_histo()

    !----------------------------------------------------------------
    !-- setup random number generator, for lepton rapidity
    !----------------------------------------------------------------
    call random_seed(size=sizeseed)
    allocate (iseed(sizeseed))
    do i = mpi_rank + 1, mpi_rank + sizeseed
       k = mod(i,sizeseed) + 1
       iseed(k) = int(123456789*sin(twopi*k*(mpi_rank+1)/sizeseed/mpi_size/sqrt(pi)))
    enddo
    call random_seed(put=iseed)
    deallocate(iseed)
  
    !----------------------------------------------------------------
    !--- run checks if needed
    !----------------------------------------------------------------
    if (do_check) then
       call do_checks()
       return
    endif
        
    !----------------------------------------------------------------
    !--- run vegas
    !----------------------------------------------------------------

    call reset_counters(all=.true.)
    call run_proc(vg_result,vg_error,vg_chi2)

    !----------------------------------------------------------------
    !--- print final results
    !----------------------------------------------------------------
    
    call print_final(outdev,vg_result,vg_error,vg_chi2)

  end subroutine chooser

  !--

  subroutine print_final(idev,vg_result,vg_error,vg_chi2)
    integer, intent(in)  :: idev
    real(dp15), intent(in) :: vg_result,vg_error,vg_chi2
    character(2) :: str_units

    if(.not.root_process) return    
    
    if (units==GeVtoFb) str_units = 'fb'
    if (units==GeVtoPb) str_units = 'pb'
    if (units==GeVtoNb) str_units = 'nb'

    write(idev,*) '################################################### '
    write(idev,*) '#'
    write(idev,*) '# ', adjustl(trim(command_line()))
    write(idev,*) '# '
#if defined(GIT_VERSION)
    write(idev,'(A,A)') " # Git version: ", GIT_VERSION
#else
    write(idev,'(A)') " # Git version: unknown"
#endif
#if defined(ANALYSIS)
    write(idev,'(A,A)') " # Analysis file: ", ANALYSIS
#else
    write(idev,'(A)') " # Analysis file unknown"
#endif    
    write(idev,*) '# '

    write(6,*) '#'
    call print_vegas_parms(idev)
    !
    call print_ew_parms(ew_scheme,idev)
    call print_qcd_parms(idev)
    !
    call print_proc_parms(idev)
    !
    write(idev,*) '################################################### '
    write(idev,*) '# '
    write(idev,*) '# Result:         ', vg_result, str_units
    write(idev,*) '# Error:          ', vg_error, str_units
    write(idev,*) '# Chi2/Iteration: ', vg_chi2
    write(idev,*) '# '
    write(idev,*) '################################################### '

    call print_counters(idev,iter=.false.,acc=.true.)
    
  end subroutine print_final
  
  subroutine print_usage(idev)
    integer :: idev

    if(.not.root_process) return    
    
    write(idev,*) 'Usage: '
    write(idev,*) ''
    write(idev,*) 'logical variables: '
    write(idev,*) 'command line: -logical --> true; in file: -logical true'
    write(idev,*) ''
    write(idev,*) "logical ``witharg'' variables: "
    write(idev,*) '-logical true in both command-line and file'
    write(idev,*)
    write(idev,*) ' -h -> print this guide'
    write(idev,*) ' -list sec/corr -> list all sectors, in computer-readable format' 
    write(idev,*) ''
    write(idev,*) ' -f input.DAT -> read input from input.DAT, overruled by command line'
    write(idev,*) ' -histodir histo -> selects histo as directory for output histograms'
    call help_vegas_parms(idev)
    call help_parms(idev)
    call help_proc_parms(idev)
    call help_run_proc(idev)
    write(idev,*) ' -check --> call check wrapper and stop'
    write(idev,*) ''
    stop
    
  end subroutine print_usage

  subroutine list_corr_ch_sec(idev)
    integer, intent(in) :: idev
    character(30) :: printlist

    if(.not.root_process) return

    printlist = trim(string_val_opt('-list','channels'))

    if(printlist.eq.'ch') then
       call print_out_channels_list(idev)
    else if(printlist.eq.'corr') then
       call print_out_corrections_list(idev)
    else if(printlist.eq.'sec') then
       call print_out_sectors_list() !-- no idev because it is interactive, it must be 6
    endif
    
    stop
    
  end subroutine list_corr_ch_sec
  
end module mod_chooser


!-- I/O module, mostly taken from Gavin, less fancy but with more consistency checks
module mod_parser
  use mod_types
  implicit none
  private

  integer, parameter :: max_arg_len  = 180
  integer, parameter :: max_line_len = 600
  integer, parameter :: max_args = 100

  integer, parameter :: verbose = 0
  
  !-- argument tracking block
  logical, allocatable, save :: argsused(:)
  logical, save :: started = .false.
  integer, save :: narg
 
  !-- input file block
  integer, parameter :: iodev  = 66
  integer, parameter :: outdev = 6
  !
  character(max_arg_len)  :: infile_opt_name(max_args),infile_opt_value(max_args)
  integer :: infile_n_values

  !--

  public :: command_line
  public :: log_val_opt,string_val_opt,int_val_opt,real_val_opt,log_val_opt_witharg
  public :: ta_CheckAllArgsUsed
  public :: parse_file
  
contains

  !-- gavin's functions, with extra option of reading input from file
  !-- simplified to use built-in read function to convert data typee
  !-- simplified to remove lcl_xxx

  !-- print the command line
  function command_line()
    implicit none
    character(len=max_line_len) ::  command_line
    !----------------------------------------------------------------------
    character(len=max_arg_len) :: string
    integer :: i, n
    
    n = iargc()
    call getarg(0,command_line)
    do i = 1, n
       call getarg(i,string)
       command_line = trim(command_line)//' '//trim(string)
    end do
  end function command_line

  !-- does not expect any argument after -opt
  function log_val_opt(opt,default_val)
    implicit none
    logical                       :: log_val_opt
    character(len=*), intent(in)  :: opt
    logical, intent(in), optional :: default_val
     !----------------------------------------------------------------------
    character(len=30) :: noopt
    integer :: i, j
    logical :: is_in_file

    !-- look into the file first
    is_in_file = .false.
    call iargc_file(opt,i,is_in_file)
    if (is_in_file) then
       read(infile_opt_value(i),*,err=111) log_val_opt
    endif
 
    
    i = len_trim(opt)
    noopt = opt(1:1)//'no'//opt(2:i)

    i = iargc_opt(trim(opt))
    j = iargc_opt(trim(noopt))
    
    if (i > 0 .neqv. j > 0) then
       log_val_opt = (i>0)
    elseif (is_in_file) then
       return
    else
       if (present(default_val)) then
          log_val_opt = default_val
       else
          log_val_opt = .false.
       end if
    end if

    return
    
111 write(outdev,*) 'problem reading variable ',trim(opt), ' from file'
    stop
    
  end function log_val_opt

  !-- functions below modified to read from file as well, if specified
  function string_val_opt(opt,default_val)
    implicit none
    character(len=max_arg_len)             :: string_val_opt
    character(len=*), intent(in)           :: opt
    character(len=*), intent(in), optional :: default_val
    integer :: i,n
    logical :: is_in_file

    !-- look into the file first
    call iargc_file(opt,i,is_in_file)
    if (is_in_file) string_val_opt = infile_opt_value(i)
    
    i = iargc_opt(opt)
    n = iargc()

    if (i >= n .or. i < 0) then
       if (is_in_file) then
          if(verbose >= 2) write(outdev,*) '# ',trim(opt),' set to ',trim(string_val_opt)
          return
       else if (present(default_val)) then
          string_val_opt = default_val
       else
          write(outdev,*) 'String value for option ',trim(opt),' has been&
               & requested'
          write(outdev,*) 'but that option is not present, and no default value&
               & was provided'
          write(outdev,*) 'use -h for help'
          stop
       end if
    else
       call getarg(i+1,string_val_opt)
       call ta_RegisterUsed(i+1)
       if(verbose >= 2) write(outdev,*) '# ',trim(opt),' set to ',trim(string_val_opt)
    end if
    
  end function string_val_opt

  !-- functions below modified to read from file as well, if specified
  function int_val_opt(opt,default_val)
    implicit none
    integer                                :: int_val_opt
    character(len=max_arg_len)             :: string_val_opt
    character(len=*), intent(in)           :: opt
    integer, intent(in), optional          :: default_val
    integer :: i,n
    logical :: is_in_file

    !-- look into the file first
    call iargc_file(opt,i,is_in_file)
    if (is_in_file) then
       read(infile_opt_value(i),*,err=111) int_val_opt
    endif
    
    i = iargc_opt(opt)
    n = iargc()

    if (i >= n .or. i < 0) then
       if (is_in_file) then
          if(verbose >= 2) write(outdev,*) '# ',trim(opt),' set to ',int_val_opt
          return
       else if (present(default_val)) then
          int_val_opt = default_val
       else
          write(outdev,*) 'String value for option ',trim(opt),' has been&
               & requested'
          write(outdev,*) 'but that option is not present, and no default value&
               & was provided'
          write(outdev,*) 'use -h for help'
          stop
       end if
    else
       call getarg(i+1,string_val_opt)
       read(string_val_opt,*,err=112) int_val_opt
       call ta_RegisterUsed(i+1)
       if(verbose >= 2) write(outdev,*) '# ',trim(opt),' set to ',int_val_opt
    end if

    return
    
111 write(outdev,*) 'problem reading variable ',trim(opt), ' from file'
    stop

112 write(outdev,*) 'problem reading variable ',opt
    stop
    
  end function int_val_opt

  !-- functions below modified to read from file as well, if specified
  function real_val_opt(opt,default_val)
    implicit none
    real(dp)                               :: real_val_opt
    character(len=max_arg_len)             :: string_val_opt
    character(len=*), intent(in)           :: opt
    real(dp), intent(in), optional         :: default_val
    integer :: i,n
    logical :: is_in_file

    !-- look into the file first
    call iargc_file(opt,i,is_in_file)
    if (is_in_file) then
       read(infile_opt_value(i),*,err=111) real_val_opt
    endif
    
    i = iargc_opt(opt)
    n = iargc()

    if (i >= n .or. i < 0) then
       if (is_in_file) then
          if(verbose >= 2) write(outdev,*) '# ',trim(opt),' set to ',real_val_opt
          return
       else if (present(default_val)) then
          real_val_opt = default_val
       else
          write(outdev,*) 'String value for option ',trim(opt),' has been&
               & requested'
          write(outdev,*) 'but that option is not present, and no default value&
               & was provided'
          write(outdev,*) 'use -h for help'
          stop
       end if
    else
       call getarg(i+1,string_val_opt)
       read(string_val_opt,*,err=112) real_val_opt
       call ta_RegisterUsed(i+1)
       if(verbose >= 2) write(outdev,*) '# ',trim(opt),' set to ',real_val_opt
    end if

    return

111 write(outdev,*) 'problem reading variable ',trim(opt), ' from file'
    stop

112 write(outdev,*) 'problem reading variable ',opt
    stop
    
  end function real_val_opt

  !-- functions below modified to read from file as well, if specified
  !-- same as log_val_opt, but requires argument
  function log_val_opt_witharg(opt,default_val)
    implicit none
    logical                                :: log_val_opt_witharg
    character(len=max_arg_len)             :: string_val_opt
    character(len=*), intent(in)           :: opt
    logical, intent(in), optional          :: default_val
    integer :: i,n
    logical :: is_in_file

    !-- look into the file first
    call iargc_file(opt,i,is_in_file)
    if (is_in_file) then
       read(infile_opt_value(i),*,err=111) log_val_opt_witharg
    endif
    
    i = iargc_opt(opt)
    n = iargc()

    if (i >= n .or. i < 0) then
       if (is_in_file) then
          if(verbose >= 2) write(outdev,*) '# ',trim(opt),' set to ',log_val_opt_witharg
          return
       else if (present(default_val)) then
          log_val_opt_witharg = default_val
       else
          write(outdev,*) 'String value for option ',trim(opt),' has been&
               & requested'
          write(outdev,*) 'but that option is not present, and no default value&
               & was provided'
          write(outdev,*) 'use -h for help'
          stop
       end if
    else
       call getarg(i+1,string_val_opt)
       read(string_val_opt,*,err=112) log_val_opt_witharg
       call ta_RegisterUsed(i+1)
       if(verbose >= 2) write(outdev,*) '# ',trim(opt),' set to ',log_val_opt_witharg
    end if

    return

111 write(outdev,*) 'problem reading variable ',trim(opt), ' from file'
    stop

112 write(outdev,*) 'problem reading variable ',opt
    stop
    
  end function log_val_opt_witharg
  
  !--------------------
  
  !-- generic functions
  !-- taken from gavin, but simplified by removing lcl_xxx interfaces

  !-- returns the position of opt in the command line
  !-- if more than one instance is present, error message
  function iargc_opt(opt)
    implicit none
    integer                      :: iargc_opt
    character(len=*), intent(in) :: opt
    !------------------------------------------------------------
    integer :: i, n
    character(len=max_arg_len) :: string
    
    n = iargc()
    do i = 1, n
       call getarg(i,string)
       if (trim(string) == trim(opt)) exit
    end do
    ! Check to see if argument found
    if (i > n) then
       !-- make sure that iargc_opt + reasonable number remains negative
       iargc_opt = -100
    else
       iargc_opt = i
       call ta_RegisterUsed(iargc_opt)
    end if
  end function iargc_opt

  !-- same as iargc, but returns position in the input file
  !-- match the last entry
  subroutine iargc_file(opt,imatch,matched)
    implicit none
    integer                      :: imatch
    character(len=*), intent(in) :: opt
    logical                      :: matched
    !------------------------------------------------------------
    integer :: i

    matched = .false.
    
    do i = 1, infile_n_values
       if (trim(infile_opt_name(i)) == trim(opt)) then
          matched = .true.
          imatch = i
       endif
    end do

    return

  end subroutine iargc_file
  
  !--------------------
  
  !-- tracking, taken from gavin
  subroutine ta_RegisterUsed(iarg)
    integer, intent(in)  :: iarg
    if (.not. started) then
       narg = iargc()
       allocate(argsused(narg))
       argsused = .false.
       started = .true.
    end if
    
    if (iarg > narg .or. iarg < 1) then
       write(outdev,*) 'In ta_RegisterUsed tried to set illegal arg', iarg
       stop
    end if
    argsused(iarg) = .true.
  end subroutine ta_RegisterUsed

  function ta_CheckAllArgsUsed(ErrDev) result(check)
    logical :: check
    integer, optional :: ErrDev
    integer :: i
    character(len=max_arg_len) :: argi

    if (narg > 0 ) then 
       check = all(argsused)
    elseif (narg == 0 .and. iargc() .ne. 0) then
       write(ErrDev,*) 'problem with command line'
       write(ErrDev,*) 'run with -h to see options'
       stop
    else
       check = .true. 
       return
    end if
    if ((.not. check) .and. present(ErrDev)) then
       write(ErrDev,*) 'ERROR: the following args were not recognized'
       do i = 1, narg
          if (.not. argsused(i)) then
             call getarg(i, argi)
             write(ErrDev,*) trim(argi)
          end if
       end do
       stop
    end if
  end function ta_CheckAllArgsUsed

  !---------------------
    
  subroutine parse_file(filename,idev,outdev)
    character(*) :: filename
    integer :: idev,outdev,io
    character(max_line_len) :: line
    character(max_arg_len)  :: val
    integer :: i,index
    logical :: itexists

    inquire(file=filename,exist=itexists)
    if (.not. itexists) then
       write (6,*) 'file ',trim(filename), ' does not exist'
       stop
    endif
    
    !--
    
    infile_opt_name  = ''
    infile_opt_value = ''
    infile_n_values  = 0
    
    !--

    if(verbose >= 1) write(outdev,*) '#'
    if(verbose >= 1) write(outdev,*) '# reading from ',trim(filename)
    if(verbose >= 1) write(outdev,*) '# careful, input from file not checked for inconsistencies'
    if(verbose >= 1) write(outdev,*) '# careful, input file overruled by command line'
    if(verbose >= 1) write(outdev,*) '#'
    
    !--
    
    i = 0
    
    open(unit=idev,file=filename,status='old')
    
    do
       read(idev,'(A)',iostat=io) line
       if (io/=0) exit
       
       if (line(1:1) == '#') cycle
       if (adjustl(trim(line)).eq.'') cycle

       i = i+1

       if (i .gt. max_args) then

          write(outdev,*) 'too many lines in ',filename
          write(outdev,*) 'change max_args in mod_parser.f90'
          stop

       endif
       
       line = adjustl(line)
       line = trim(line)
       index = scan(line,' ' )

       infile_opt_name(i)  = line(1:index-1)
       val          = line(index+1:)
       
       if (val == '') then
          
          write(outdev,*) 'problem in ',trim(filename)
          write(outdev,*) 'no option passed for ', infile_opt_name(i)
          stop
          
       else

          index = scan(val,'#')
          if (index.ne.0) val = val(1:index-1)

          infile_opt_value(i) = val

       endif
       
    enddo

    close(idev)

    infile_n_values = i
    
    return
    
  end subroutine parse_file
  
end module mod_parser

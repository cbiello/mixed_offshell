program mysummer
  implicit none
  
  character *(3) :: operation

  call GetArg(1,operation)

  if (operation.ne.'add' .and. operation.ne.'avg' .and. operation.ne.'nrm' .and. &
       operation.ne.'pls' .and. operation.ne.'rbn' .and. operation.ne.'sbn' &
       .and. operation.ne.'mul' .and. operation .ne. 'pck' .and. operation .ne. 'rmv' &
       .and. operation .ne. 'mrg' .and. operation .ne. 'kfc' .and. operation .ne. 'rb2' &
       .and. operation .ne. 'sym' .and. operation .ne. 'ad3') then
     write(6,*) 'usage: ./mysummer.x add filein_1 filein_2 ... filein_n file_out'
     write(6,*) 'usage: ./mysummer.x ad3 filein_1 filein_2 ... filein_n file_out --> 3-column files'
     write(6,*) 'usage: ./mysummer.x mul n filein file_out'
     write(6,*) '       ./mysummer.x pls filein_1 filein_2 ... filein_n file_out'
     write(6,*) '       ./mysummer.x avg filein_1 filein_2 ... filein_n file_out'
     write(6,*) '       ./mysummer.x nrm filein file_out'
     write(6,*) '       ./mysummer.x rbn nhisto nrebin filein file_out'
     write(6,*) '       ./mysummer.x sbn nhisto filein file_out'
     write(6,*) '       ./mysummer.x pck nhisto filein file_out'
     write(6,*) '       ./mysummer.x rmv nhisto minval maxval filein file_out'
     write(6,*) '       ./mysummer.x mrg nhisto minval maxval filein1 filein2 file_out'
     write(6,*) '       ./mysummer.x kfc nhisto filein1 filein2 file_out'
     write(6,*) '       ./mysummer.x rb2 nhisto nrebin xmin xmax filein file_out'
     write(6,*) '       ./mysummer.x sym nhisto filein file_out'
     stop
  endif

  if (operation.eq.'ad3') call mysummer_add3()
  if (operation.eq.'add') call mysummer_add()
  if (operation.eq.'add') call mysummer_add()
  if (operation.eq.'mul') call mysummer_mul()
  if (operation.eq.'pls') call mysummer_pls()
  if (operation.eq.'avg') call mysummer_avg()
  if (operation.eq.'rbn') call mysummer_rebin()
  if (operation.eq.'sbn') call mysummer_sumbin()
  if (operation.eq.'pck') call mysummer_pck()
  if (operation.eq.'rmv') call mysummer_rmvbins()
  if (operation.eq.'mrg') call mysummer_mrg()
  if (operation.eq.'kfc') call mysummer_kfactor()
  if (operation.eq.'rb2') call mysummer_rebin_limited()
  if (operation.eq.'sym') call mysummer_symmetrize()

end program mysummer

!-----------------------------------------------------------------------------------------------
! REWRITTEN
!-----------------------------------------------------------------------------------------------

!-- sum n histograms, 3 columns
subroutine mysummer_add3()
    implicit none
    
    character *(300) :: line
    character *(120) :: file_in(100), buffer
    character(len=:), allocatable :: file_out
    integer :: i, j, NumArgs, nfiles, nhisto, nbin
    real(8) :: bin(500), value(3,500), error(3,500)
    real(8) :: sumvalue(3), sumerror(3)
    character(len=*), parameter :: fmt1 = "(I2,A,I3,A,2X,1PE11.4,A,2X,1PE23.16,A,2X,1PE23.16,A)"
    character(len=1), parameter :: space = ' '
  
    character(len=15), dimension(:), allocatable :: id
    logical :: is_comment, flag
    
    nfiles = COMMAND_ARGUMENT_COUNT() - 2
    allocate( id(nfiles) )
    
    call print_logo()
    write(6,'(/,a,i0,a,/)') 'Adding ', nfiles, ' files together'
  
    do i = 1, nfiles
        call GET_COMMAND_ARGUMENT( i + 1, file_in(i) )
        file_in(i) = trim(file_in(i))
        write(*,'(a5,i3,a2,a120)') 'file ', i,': ', file_in(i)
        open( unit = 21 + i, file = file_in(i), status = 'old' )
    enddo
    
    call GET_COMMAND_ARGUMENT( nfiles + 2, buffer )
    file_out = buffer
    
    write(*,'(a,a,/)') 'file out: ', file_out
    open( unit = 21, file = file_out, status = 'unknown' )
    
    do while (.true.) ! loop over different histograms
        do i = 1, nfiles, 1
            call move_curser_to_next_histogram( 21 + i, id(i), flag )
            if ( flag ) goto 10
        end do
        
        if ( .not. all(id .eq. id(1)) ) goto 102
        write(21,'(a1)') '#'
        write(21,'(a1)') '#'
        write(21,'(a16,a15)') '# ----------Obs=', id(1)
        
        do while (.true.) ! loop over entries of one histogram
            read( 21 + 1,'(a)', end = 10 ) line
            backspace( 21 + 1 )
            if ( is_comment( line ) ) exit
             
            sumvalue = 0d0
            sumerror = 0d0

            do i = 1, nfiles, 1
                read( 21 + i, *, end = 10) nhisto, nbin, bin(i), &
                    value(1,i), error(1,i), value(2,i), error(2,i), value(3,i), error(3,i)

                do j=1,3
                   if ( abs(value(j,i)) .le. 1d-30 ) error(j,i) = 0d0
                   sumvalue(j) = sumvalue(j) + value(j,i)
                   sumerror(j) = sumerror(j) + error(j,i)**2
                enddo
            end do
            
            sumerror = dsqrt(sumerror)
            write(21,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A,1PE23.16,A,2X,1PE23.16,A,1PE23.16,A,2X,1PE23.16,A)") &
                 nhisto, space, nbin, space, bin(1),space,&
                 sumvalue(1),space,sumerror(1),space,sumvalue(2),space,sumerror(2),space,sumvalue(3),space,sumerror(3),space
        end do
    end do

10  continue
    return
    
102 continue
    write(*,'(a)') 'error: histograms in the files are not compatible (different histogram names)'
    return
end subroutine mysummer_add3


!-- sum n histograms
subroutine mysummer_add()
    implicit none
    
    character *(300) :: line
    character *(120) :: file_in(100), buffer
    character(len=:), allocatable :: file_out
    integer :: i, NumArgs, nfiles, nhisto, nbin
    real(8) :: bin(200), value(200), error(200)
    real(8) :: sumvalue, sumerror
    character(len=*), parameter :: fmt1 = "(I2,A,I3,A,2X,1PE11.4,A,2X,1PE23.16,A,2X,1PE23.16,A)"
    character(len=1), parameter :: space = ' '
  
    character(len=15), dimension(:), allocatable :: id
    logical :: is_comment, flag
    
    nfiles = COMMAND_ARGUMENT_COUNT() - 2
    allocate( id(nfiles) )
    
    call print_logo()
    write(6,'(/,a,i0,a,/)') 'Adding ', nfiles, ' files together'
  
    do i = 1, nfiles
        call GET_COMMAND_ARGUMENT( i + 1, file_in(i) )
        file_in(i) = trim(file_in(i))
        write(*,'(a5,i3,a2,a120)') 'file ', i,': ', file_in(i)
        open( unit = 21 + i, file = file_in(i), status = 'old' )
    enddo
    
    call GET_COMMAND_ARGUMENT( nfiles + 2, buffer )
    file_out = buffer
    
    write(*,'(a,a,/)') 'file out: ', file_out
    open( unit = 21, file = file_out, status = 'unknown' )
    
    do while (.true.) ! loop over different histograms
        do i = 1, nfiles, 1
            call move_curser_to_next_histogram( 21 + i, id(i), flag )
            if ( flag ) goto 10
        end do
        
        if ( .not. all(id .eq. id(1)) ) goto 102
        write(21,'(a1)') '#'
        write(21,'(a1)') '#'
        write(21,'(a16,a15)') '# ----------Obs=', id(1)
        
        do while (.true.) ! loop over entries of one histogram
            read( 21 + 1,'(a)', end = 10 ) line
            backspace( 21 + 1 )
            if ( is_comment( line ) ) exit
             
            sumvalue = 0d0
            sumerror = 0d0

            do i = 1, nfiles, 1
                read( 21 + i, *, end = 10) nhisto, nbin, bin(i), value(i), error(i)
                
                if ( abs(value(i)) .le. 1d-30 ) error(i) = 0d0
                sumvalue = sumvalue + value(i)
                sumerror = sumerror + error(i)**2
            end do
            
            sumerror = dsqrt(sumerror)
            write(21,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)") nhisto, space, nbin, &
                space, bin(1),space,sumvalue,space,sumerror,space
        end do
    end do

10  continue
    return
    
102 continue
    write(*,'(a)') 'error: histograms in the files are not compatible (different histogram names)'
    return
end subroutine mysummer_add


!-- average n histograms
subroutine mysummer_avg()
    implicit none
    
    character(len=300) :: line
    character *(120) :: file_in(1000), buffer
    character(len=:), allocatable :: file_out
    integer :: i, NumArgs, nfiles, nhisto, nbin, icount
    real(8) :: bin(1000), value(1000), error(1000)
    real(8) :: sumvalue, sumerror
    real(8) :: sumvalue2, sumerror2, sumerror3,sumerror4
    character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  
    character(len=1), parameter :: space = ' '
    
    character(len=15), dimension(:), allocatable :: id
    logical :: is_comment, flag
    
    nfiles = ( command_argument_count() + 1 ) - 3
    allocate( id(nfiles) )
    
    call print_logo()
    write(6,'(/,a,i0,a,/)') 'Averaging ', nfiles, ' files'
    
    do i = 1, nfiles
        call getarg( i + 1, file_in(i) )
        file_in(i) = trim(file_in(i))
        write(*,'(a5,i3,a2,a120)') 'file ', i,': ', file_in(i)
        open( unit = 11 + i, file = file_in(i), status = 'old' )
    enddo
    
    call getarg( nfiles + 2, buffer )
    file_out = buffer
    
    write(*,'(a,a,/)') 'file out: ', file_out
    open( unit = 11, file = file_out, status = 'unknown' )
    
    do while (.true.) ! loop over different histograms
        do i = 1, nfiles, 1
            call move_curser_to_next_histogram( 11 + i, id(i), flag )
            if ( flag ) goto 20
        end do
        
        if ( .not. all(id .eq. id(1)) ) goto 102
        write(11,'(a1)') '#'
        write(11,'(a1)') '#'
        write(11,'(a16,a15)') '# ----------Obs=', id(1)
        
        do while (.true.) ! loop over entries of one histogram
            read( 11 + 1,'(a)', end = 20 ) line
            backspace( 11 + 1 )
            if ( is_comment( line ) ) exit
                
            sumvalue = 0d0
            sumerror = 0d0
            sumvalue2 = 0d0
            sumerror2 = 0d0
            sumerror3 = 0d0
            sumerror4 = 0d0
            icount = 0
        
            do i = 1, nfiles, 1
                read( 11 + i, *, end = 20) nhisto, nbin, bin(i), value(i), error(i)
                
                icount = icount + 1
                !-- normal combination
                sumvalue = sumvalue + value(i)
                sumerror = sumerror + error(i)**2
                !-- gaussian combination
                !-- avoid bins with zero events and zero error
                if((value(i) .ne. 0d0) .and. (error(i) .ne. 0d0)) then
                  sumvalue2 = sumvalue2 + value(i)/error(i)**2
                  sumerror2 = sumerror2 + 1d0/error(i)**2
                endif
                !-- error is just the spread, gaussian
                sumerror3 = sumerror3 + value(i)**2/error(i)**2
            end do
            
            !-- normal combination
            sumerror = sumerror/icount
            sumvalue = sumvalue/icount
            sumerror = dsqrt(sumerror)
            
            !-- gaussian combination
            sumvalue2 = sumvalue2/sumerror2
            !-- error is just the spread, gaussian central value
            sumerror3 = sumerror3/sumerror2
            sumerror3 = dsqrt(abs(sumerror3-sumvalue2**2)) !-- abs should not be there, safety precaution
            !-- original gaussian error
            sumerror2 = 1d0/sumerror2
            sumerror2 = dsqrt(sumerror2)

            if (sumvalue.eq.0d0 .and. sumerror2 .gt. 1000d0) sumerror2 = 0d0
            if (sumvalue.eq.0d0 .and. sumerror3 .gt. 1000d0) sumerror3 = 0d0
            if (sumvalue.lt.1d-50) sumvalue = 0d0
            if (.not.sumvalue2.le.0d0 .and. .not.sumvalue2.gt.0d0 .and. sumerror2.eq.0d0) sumvalue2 = 0d0 !-- remove NaN from 0 error

            !-- gaussian
            if ( nfiles .le. 5 ) then
                write(11,fmt="(I2,A,X,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)") nhisto, space, &
                    nbin, space, bin(1),space,sumvalue2,space,sumerror2,space
            else
                write(11,fmt="(I2,A,X,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)") nhisto, space, &
                    nbin, space, bin(1),space,sumvalue2,space,sumerror2,space
            endif
        end do
    end do
    
20  continue
    return
    
102 continue
    write(*,'(a)') 'error: histograms in the files are not compatible (different histogram names)'
    return
end subroutine mysummer_avg

!-- rebin n histograms
subroutine mysummer_rebin()
    implicit none
    character *(1) :: dummy
    character *(150) :: file_in, file_out, charnrebin, charnhisto
    integer :: i, NumArgs, ihisto,nhisto, nrebin, nbin, ibin
    real(8) :: bin, value, error
    real(8) :: sumvalue, sumerror
    integer, parameter :: tothisto = 4
    character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)" 
    
    character(len=120) :: buffer
    character(len=:), allocatable :: input
    integer :: buffer_length, ios, sumnbin
    logical :: is_comment
    real(8) :: sumbin
        
    call print_logo()
    write(6,'(/,a,/)') 'Rebinning'

    call GetArg(2,charnhisto) !-- which histogram to rebin
    read(charnhisto,*) nhisto
    call GetArg(3,charnrebin) !-- how many bin to combine
    read(charnrebin,*) nrebin
    call GetArg(4,file_in) 
    call GetArg(5,file_out)

    file_in = trim(file_in)
    file_out = trim(file_out)
    
    write(6,'(a,i0)') 'histogram number: ', nhisto
    write(6,'(a,i0)') 'number of bins to combine: ', nrebin
    write(6,'(a,a)') 'file in: ', file_in
    write(6,'(a,a,/)') 'file out: ', file_out

    open (unit=12,file=file_in,status='old')
    open (unit=11,file=file_out,status='unknown')
    
    dummy = ' '
    sumnbin = 0
    
    do while (.true.)
        !-- get full line
        input = ''
        do while ( .true. )
            read(12,'(a)', end = 50, advance = 'no', iostat = ios, size = buffer_length ) buffer
            input = input // buffer(1:buffer_length)
            if ( is_iostat_eor(ios) ) exit
        end do
        
        !-- pass to output file if comment
        if ( is_comment(input) ) then
            write(11,'(a)') input
            cycle
        end if
        
        !-- read input as histogram bin
        backspace(unit=12)
        read(12,*, end = 50) ihisto, nbin, bin, value, error
        
        !-- do nothing if wrong histogram
        if ( ihisto .ne. nhisto ) then
            write(11,'(a)') input
            cycle
        end if
        
        sumnbin = sumnbin + 1
        sumbin = bin
        sumvalue = value
        sumerror = error**2
        
        do ibin = 1, nrebin - 1, 1
            read(12,*, iostat = ios) ihisto, nbin, bin, value, error
            if ( is_iostat_end(ios) ) then
                backspace(unit=12)
                exit
            end if
            if ( ios .ne. 0 ) exit
            if ( ihisto .ne. nhisto ) then
                backspace(unit=12)
                exit
            end if
            sumvalue = sumvalue + value
            sumerror = sumerror + error**2
        end do
        
        sumerror = dsqrt(sumerror)
        write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16)") ihisto, &
            dummy, sumnbin, dummy, sumbin, dummy, sumvalue, dummy, sumerror
    end do
    
20  continue
    return
    
50  continue
    close(11)
    close(12)
    return
end subroutine mysummer_rebin


!-- rebin n histograms, but apply only to bins with xmin < x < xmax
subroutine mysummer_rebin_limited()
    implicit none
    character *(1) :: dummy
    character *(150) :: file_in, file_out, charnrebin, charnhisto,charxmin,charxmax
    integer :: i, NumArgs, ihisto, iihisto,nhisto, nrebin, nbin, ibin
    real(8) :: bin, value, error,xmin,xmax
    real(8) :: sumvalue, sumerror
    integer, parameter :: tothisto = 4
    character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  

    character(len=120) :: buffer
    character(len=:), allocatable :: input
    integer :: buffer_length, ios, sumnbin
    logical :: is_comment
    real(8) :: sumbin
    
    call GetArg(2,charnhisto) !-- which histogram to rebin
    read(charnhisto,*) nhisto
    call GetArg(3,charnrebin) !-- how many bin to combine
    read(charnrebin,*) nrebin
    
    call GetArg(4,charxmin) !-- how many bin to combine
    read(charxmin,*) xmin
    call GetArg(5,charxmax) !-- how many bin to combine
    read(charxmax,*) xmax
    call GetArg(6,file_in) 
    call GetArg(7,file_out)

    file_in = trim(file_in)
    file_out = trim(file_out)

    call print_logo()
    write(6,'(/,a)') 'Rebinning range'
    write(6,'(a)',advance = 'no') 'xmin: '
    write(6,*) xmin
    write(6,'(a)',advance = 'no') 'xmax: '
    write(6,*) xmax
    write(6,*) ''
    
    write(6,'(a,i0)') 'histogram number: ', nhisto
    write(6,'(a,i0)') 'number of bins to combine: ', nrebin
    write(6,'(a,a)') 'file in: ', file_in
    write(6,'(a,a,/)') 'file out: ', file_out
    
    open (unit=12,file=file_in,status='old')
    open (unit=11,file=file_out,status='unknown')
    
    dummy = ' '
    sumnbin = 0
    
    do while (.true.)
        !-- get full line
        input = ''
        do while ( .true. )
            read(12,'(a)', end = 50, advance = 'no', iostat = ios, size = buffer_length ) buffer
            input = input // buffer(1:buffer_length)
            if ( is_iostat_eor(ios) ) exit
        end do
        
        !-- pass to output file if comment
        if ( is_comment(input) ) then
            write(11,'(a)') input
            cycle
        end if
        
        !-- read input as histogram bin
        backspace(unit=12)
        read(12,*, end = 50) ihisto, nbin, bin, value, error
        
        !-- do nothing if wrong histogram
        if ( ihisto .ne. nhisto ) then
            write(11,'(a)') input
            cycle
        end if
        
        !-- correct histogram
        sumnbin = sumnbin + 1
        
        !-- out of range
        if ( bin .lt. xmin .or. bin .gt. xmax ) then
            write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16)") ihisto, &
                dummy, sumnbin, dummy, bin, dummy, value, dummy, error
            cycle
        end if
        
        sumbin = bin
        sumvalue = value
        sumerror = error**2
        
        do ibin = 1, nrebin - 1, 1
            read(12,*, iostat = ios) ihisto, nbin, bin, value, error
            if ( is_iostat_end(ios) ) then
                backspace(unit=12)
                exit
            end if
            if ( ios .ne. 0 ) exit
            if ( ihisto .ne. nhisto .or. bin .lt. xmin .or. bin .gt. xmax ) then
                backspace(unit=12)
                exit
            end if
            sumvalue = sumvalue + value
            sumerror = sumerror + error**2
        end do
        
        sumerror = dsqrt(sumerror)
        write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16)") ihisto, &
            dummy, sumnbin, dummy, sumbin, dummy, sumvalue, dummy, sumerror
    end do
    
20  continue
    return
    
50  continue
    close(11)
    close(12)
    return
end subroutine mysummer_rebin_limited


!-- rebin n histograms
subroutine mysummer_symmetrize()
    implicit none
    character *(1) :: dummy
    character *(150) :: file_in, file_out, charnhisto
    integer :: i, NumArgs, ihisto,nhisto
    real(8) :: bin, value, error
    real(8) :: sumvalue, sumerror
    integer, parameter :: tothisto = 4
    character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)" 
    
    character(len=120) :: buffer
    character(len=:), allocatable :: input
    integer :: buffer_length, ios, sumnbin, nbins, nbin
    logical :: is_comment, flag
    real(8) :: sumbin
        
    real(8), allocatable :: bins(:), values(:), errors(:)
    real(8) :: symvalue, symerror
        
    call print_logo()
    write(6,'(/,a,/)') 'Symmetrizing'

    call GetArg(2,charnhisto) !-- which histogram to symmetrize
    read(charnhisto,*) nhisto
    call GetArg(3,file_in) 
    call GetArg(4,file_out)

    file_in = trim(file_in)
    file_out = trim(file_out)
    
    write(6,'(a,i0)') 'histogram number: ', nhisto
    write(6,'(a,a)') 'file in: ', file_in
    write(6,'(a,a,/)') 'file out: ', file_out

    open (unit=12,file=file_in,status='old')
    open (unit=11,file=file_out,status='unknown')
    
    dummy = ' '
    sumnbin = 0
    
    do while (.true.)
        !-- get full line
        input = ''
        do while ( .true. )
            read(12,'(a)', end = 50, advance = 'no', iostat = ios, size = buffer_length ) buffer
            input = input // buffer(1:buffer_length)
            if ( is_iostat_eor(ios) ) exit
        end do
        
        !-- pass to output file if comment
        if ( is_comment(input) ) then
            write(11,'(a)') input
            cycle
        end if
        
        !-- read input as histogram bin
        backspace(unit=12)
        read(12,*, end = 50) ihisto, nbin, bin, value, error
        
        !-- do nothing if wrong histogram
        if ( ihisto .ne. nhisto ) then
            write(11,'(a)') input
            cycle
        end if
        
        !----------------------------------------------------------------
        ! Symmetrize histogram nhisto
        !----------------------------------------------------------------
        
        !-- get number of bins
        nbins = 1
        do while (.true.)
            read(12,*, iostat = ios) ihisto, nbin, bin, value, error
            if ( ios .ne. 0 ) exit
            if ( ihisto .ne. nhisto ) exit
            nbins = nbins + 1
        end do
        
        !-- move curser to beginning of nhisto
        call MOVE_TO_BEGINNING( 12, nhisto, flag )
        
        !-- allocate arrays to store full histogram
        allocate( bins(nbins) )
        allocate( values(nbins) )
        allocate( errors(nbins) )
        
        !-- read full histogram
        do i = 1, nbins, 1
            read(12,*) ihisto, nbin, bins(i), values(i), errors(i)
        end do
        
        !-- write symmetrized
        do i = 1, nbins - 1, 1
            symvalue = ( values(i) + values(nbins - i) ) / 2.0
            symerror = sqrt( (errors(i)**2 + errors(nbins - i)**2) / 2.0 )
            write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16)") nhisto, &
                dummy, i, dummy, bins(i), dummy, symvalue, dummy, symerror
        end do
        
        !-- write overflow / underflow bin
        write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16)") nhisto, &
            dummy, nbins, dummy, bins(nbins), dummy, values(nbins), dummy, errors(nbins)
        
        !-- deallocate and continue
        deallocate( bins )
        deallocate( values )
        deallocate( errors )
    end do
    
50  continue
    close(11)
    close(12)
    return
end subroutine mysummer_symmetrize

!-----------------------------------------------------------------------------------------------
! OLD
!-----------------------------------------------------------------------------------------------

!-- multiply by a constant histograms
subroutine mysummer_mul()
  implicit none
  character *(1) :: dummy
  character *(120) :: file_in, file_out, chnum
  integer :: i, NumArgs, n_in_files, nhisto, nbin
  real(8) :: bin, value, error, num
  character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  

  n_in_files = (COMMAND_ARGUMENT_COUNT()+1)-4

  if (n_in_files.ne.1) then
     write(6,*) 'only one file at the time in input'
     stop
  endif

  write(6,*) 'Multiply all histograms by a factor'

  call GetArg(2,chnum)
  read(chnum,*) num
  print *, 'multiply by ', num

  call GetArg(3,file_in)
  file_in = adjustl(trim(file_in))
  print *, 'file in: ', file_in
  
  open (unit=12,file=file_in,status='old')
  
  call GetArg(4,file_out)
  file_out = adjustl(trim(file_out))
  print *, 'file out: ', file_out

  open (unit=11,file=file_out,status='unknown')

  do while (.true.)
     read(12,*,end=35) dummy
     if(dummy(1:1).eq."#") cycle
     backspace(unit=12)
     read(12,*,end=35) nhisto, nbin, bin, value, error
     value = value * num
     error = error * num
     if(dummy(1:1).eq."#") cycle
     dummy = " "
     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)") nhisto, &
          dummy, nbin, dummy, bin,dummy,value,dummy,error,dummy
 
  enddo

35 continue

end subroutine mysummer_mul

!-- put a single file in a suitable form
subroutine mysummer_nrm()
  implicit none
  character *(1) :: dummy
  character *(70) :: file_in, file_out
  integer :: i, NumArgs, n_in_files, nhisto, nbin, icount
  real(8) :: bin, value, error
  character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  

  n_in_files = (COMMAND_ARGUMENT_COUNT()+1)-3

  if (n_in_files.ne.1) then
     write(6,*) 'only one file at the time in input'
     stop
  endif

  write(6,*) 'Put in a form suitable for combination'

  call GetArg(2,file_in)
  print *, 'file in: ', file_in
  file_in = trim(file_in)

  open (unit=12,file=file_in,status='old')
  
  call GetArg(3,file_out)
  file_out = trim(file_out)
  print *, 'file out: ', file_out

  open (unit=11,file=file_out,status='unknown')

  do while (.true.)
     read(12,*,end=30) dummy
     if(dummy(1:1).eq."#") cycle
     backspace(unit=12)
     read(12,*,end=30) nhisto, nbin, bin, value, error
     if(dummy(1:1).eq."#") cycle
     dummy = " "
     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)") nhisto, &
          dummy, nbin, dummy, bin,dummy,value,dummy,error,dummy
 
  enddo

30 continue

end subroutine mysummer_nrm

  
!-- sum n histograms to show poles cancellation -> divide by absolute value, ie res = (RR + RV + VV )/(1/3(|RR|+|RV|+|VV|))
subroutine mysummer_pls()
  implicit none
  character *(1) :: dummy
  character *(70) :: file_in(50), file_out
  integer :: i, NumArgs, n_in_files, nhisto, nbin
  real(8) :: bin(2), value(2), error(2)
  real(8) :: sumvalue, sumerror, meanabsvalue
  character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  

  n_in_files = (COMMAND_ARGUMENT_COUNT()+1)-3

  write(6,*) 'For poles: sum contributions, with meaningful normalization:'
  write(6,*) 'res = (RR+RV+VV+...)/( (|RR|+|RV|+|VV|+...)/N)'

  do i = 1, n_in_files
     call GetArg(i+1,file_in(i))
     file_in(i) = trim(file_in(i))
     write(*,fmt="(A5,I2,A2,A40)") 'file ', i,': ', file_in(i)
     open (unit=11+i,file=file_in(i),status='old')
  enddo
  
  call GetArg(n_in_files+2,file_out)
  file_out = trim(file_out)
  print *, 'file out: ', file_out

  open (unit=11,file=file_out,status='unknown')

  do while (.true.)
     meanabsvalue = 0d0
     sumvalue = 0d0
     sumerror = 0d0
     do i = 1, n_in_files
        read(11+i,*,end=40) dummy
        if(dummy(1:1).eq."#") cycle
        backspace(unit=11+i)
        read(11+i,*,end=40) nhisto, nbin, bin(i),value(i),error(i)
        sumvalue = sumvalue + value(i)
        meanabsvalue = meanabsvalue + abs(value(i))/n_in_files
        sumerror = sumerror + error(i)**2
     enddo
     if(dummy(1:1).eq."#") cycle

     sumerror = dsqrt(sumerror)

     !-- divide by the mean to see relative cancellation
     sumvalue = sumvalue/meanabsvalue
     sumerror = sumerror/meanabsvalue
     dummy = " "
     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)") nhisto, dummy, &
          nbin, dummy, bin(1),dummy,sumvalue,dummy,sumerror,dummy
  enddo

40 continue

end subroutine mysummer_pls




!-- sum the n histogram
subroutine mysummer_sumbin()
  implicit none
  character *(1) :: dummy
  character *(100) :: file_in, file_out, charnhisto
  integer :: i, NumArgs, ihisto, nhisto, nbin(1500), ibin
  real(8) :: bin(1500), value, error, binmin, binmax
  real(8) :: sumvalue, sumerror, summom1, errmom1, summom2, errmom2
  real(8) :: lastvalue, lasterror
  integer, parameter :: tothisto = 4
  character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  

  write(6,*) 'Sum over bins'

  call GetArg(2,charnhisto) !-- which histogram to sum (-1 -> all of them)
  read(charnhisto,*) nhisto
  call GetArg(3,file_in) 

  file_in = trim(file_in)

  print *, 'file in: ', file_in
  print *, 'nhisto: ', nhisto

  open (unit=12,file=file_in,status='old')

  sumvalue = 0d0
  sumerror = 0d0
  summom1 = 0d0
  errmom1 = 0d0
  summom2 = 0d0
  errmom2 = 0d0
  ibin = 1
  do while (.true.)
     read(12,*,end=60) dummy
     if(dummy(1:1).eq."#") cycle
     backspace(unit=12)

     read(12,*,end=60) ihisto, nbin(ibin), bin(ibin),value,error
     if (ihisto.eq.nhisto) then
        binmin = bin(1)
        sumvalue = sumvalue + value
        if (value .ne. 0d0) then
           lastvalue = value
           lasterror = error
           sumerror = sumerror + error**2
        endif
        
        !summom1 = summom1 + value*bin(ibin)
        !errmom1 = errmom1 + (error*bin(ibin))**2
        !summom2 = summom2 + value*(3d0 * bin(ibin)**2 - 1d0)
        !errmom2 = errmom2 + (error*(3d0 * bin(ibin)**2 - 1d0))**2

        binmax = bin(ibin)
        ibin = ibin + 1
     endif
     
     if(dummy(1:1).eq."#") cycle

  enddo

60 continue

  dummy = " "
  sumerror = dsqrt(sumerror)
  errmom1 = dsqrt(errmom1)
  errmom2 = dsqrt(errmom2)
  write(*,*) 'sum: ', sumvalue
  write(*,*) 'error: ', sumerror
  write(*,*) 'sum without last bin: ', sumvalue-lastvalue
  write(*,*) 'error without last bin: ', dsqrt(sumerror**2-lasterror**2)
  !write(*,*) 'first moment: ', summom1
  !write(*,*) 'error: ', errmom1
  !write(*,*) 'second moment: ', summom2
  !write(*,*) 'error: ', errmom2
  return




end subroutine mysummer_sumbin

!-- pick out a given histogram and write only it to a file 
subroutine mysummer_pck()
  implicit none
  character *(1) :: dummy
  character *(100) :: file_in, file_out, charnhisto
  integer :: i, NumArgs, n_in_files, nhisto, nbin, thehisto
  real(8) :: bin, value, error, num
  character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  

  n_in_files = (COMMAND_ARGUMENT_COUNT()+1)-4

  if (n_in_files.ne.1) then
     write(6,*) 'only one file at the time in input'
     stop
  endif

  write(6,*) 'Pick out a given histogram'

  call GetArg(2,charnhisto) !-- which histogram to pick out
  read(charnhisto,*) thehisto

  call GetArg(3,file_in)
  file_in = adjustl(trim(file_in))
  print *, 'file in: ', file_in
  
  open (unit=12,file=file_in,status='old')
  
  call GetArg(4,file_out)
  file_out = adjustl(trim(file_out))
  print *, 'file out: ', file_out

  open (unit=11,file=file_out,status='unknown')

  do while (.true.)
     read(12,*,end=35) dummy
     if(dummy(1:1).eq."#") cycle
     backspace(unit=12)
     read(12,*,end=35) nhisto, nbin, bin, value, error
     if(dummy(1:1).eq."#") cycle
     dummy = " "
     if (nhisto .eq. thehisto) then
        write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)") nhisto, &
             dummy, nbin, dummy, bin,dummy,value,dummy,error,dummy
     endif
 
  enddo

35 continue

end subroutine mysummer_pck

!-- pick out a given histogram and write only it to a file 
subroutine mysummer_rmvbins()
  implicit none
  character *(1) :: dummy
  character *(100) :: file_in, file_out, charnhisto, charminval, charmaxval
  integer :: i, NumArgs, n_in_files, nhisto, nbin, thehisto
  real(8) :: bin, value, error, num, minval, maxval
  character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  

  n_in_files = (COMMAND_ARGUMENT_COUNT()+1)-6

  if (n_in_files.ne.1) then
     write(6,*) 'only one file at the time in input'
     stop
  endif



  call GetArg(2,charnhisto) !-- which histogram to pick out
  read(charnhisto,*) thehisto

  call GetArg(3,charminval) !-- which histogram to pick out
  read(charminval,*) minval

  call GetArg(4,charmaxval) !-- which histogram to pick out
  read(charmaxval,*) maxval

  write(6,*) 'Removing bins with x-value <', minval, "and x-value > ",maxval

  call GetArg(5,file_in)
  file_in = adjustl(trim(file_in))
  print *, 'file in: ', file_in
  
  open (unit=12,file=file_in,status='old')
  
  call GetArg(6,file_out)
  file_out = adjustl(trim(file_out))
  print *, 'file out: ', file_out

  open (unit=11,file=file_out,status='unknown')

  do while (.true.)
     read(12,*,end=35) dummy
     if(dummy(1:1).eq."#") cycle
     backspace(unit=12)
     read(12,*,end=35) nhisto, nbin, bin, value, error
     if(dummy(1:1).eq."#") cycle
     dummy = " "
!     if (nhisto .ne. thehisto) cycle
     !     if (bin .lt. minval .or. bin .gt. maxval) cycle
     if (nhisto .eq. thehisto) then
        if (bin .lt. minval .or. bin .gt. maxval) cycle
     endif
     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)") nhisto, &
          dummy, nbin, dummy, bin,dummy,value,dummy,error,dummy
!          write(11,fmt="(1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"),  &
!               bin,dummy,value,dummy,error,dummy
  
 
  enddo

35 continue

end subroutine mysummer_rmvbins

!
!-- pick out a given histogram and write only it to a file 
subroutine mysummer_mrg()
  implicit none
  character *(1) :: dummy1,dummy2
  character *(100) :: file_in1,file_in2, file_out, charnhisto, charminval, charmaxval
  integer :: i, NumArgs, n_in_files,nhisto1,nhisto2,nbin1,nbin2, thehisto
  real(8) :: bin1,bin2,binout,value1,value2,valueout,error1,error2,errorout, num, minval, maxval
  character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  

  n_in_files = (COMMAND_ARGUMENT_COUNT()+1)-6

  if (n_in_files.ne.2) then
     write(6,*) 'two files at the time in input'
     stop
  endif



  call GetArg(2,charnhisto) !-- which histogram to pick out
  read(charnhisto,*) thehisto

  call GetArg(3,charminval) !-- min value at which to start merg
  read(charminval,*) minval

  call GetArg(4,charmaxval) !-- min value at which to start merg
  read(charmaxval,*) maxval

  call GetArg(5,file_in1)
  file_in1 = adjustl(trim(file_in1))
  print *, 'file in 1: ', file_in1

  call GetArg(6,file_in2)
  file_in2 = adjustl(trim(file_in2))
  print *, 'file in 2: ', file_in2
  close(12)
  close(13)
  open (unit=12,file=file_in1,status='old')
  open (unit=13,file=file_in2,status='old')

  write(6,*) 'Merging: using file1 for bins with x-value <', minval, "and x-value > ",maxval
  write(6,*) 'Using and file2 for bins with', minval, "<  x-value < ",maxval
   
  call GetArg(7,file_out)
  file_out = adjustl(trim(file_out))
  print *, 'file out: ', file_out
  open (unit=11,file=file_out,status='unknown')

  do while (.true.)
     read(12,*,end=35) dummy1
     read(13,*,end=35) dummy2
     if(dummy1(1:1).eq."#") then
        if (dummy2(1:1) .ne. "#") then
           print *, "error"
           stop
        else
           cycle
        endif
     endif
     backspace(unit=12)
     backspace(unit=13)
     read(12,*,end=35) nhisto1, nbin1, bin1, value1, error1
     read(13,*,end=35) nhisto2, nbin2, bin2, value2, error2
     if(dummy1(1:1).eq."#") then
        if (dummy2(1:1) .ne. "#") then
           print *, "error"
           stop
        else
           cycle
        endif
     endif
     if (nhisto1 .ne. nhisto2 .or. bin1 .ne. bin2 .or. nbin1 .ne. nbin2) then
        print *, "error"
        stop
     endif
     dummy1 = " "
     dummy2 = " "
     if (nhisto1 .eq. thehisto) then
        if (bin1 .lt. minval .or. bin1 .gt. maxval) then
           binout   = bin1
           valueout = value1
           errorout = error1
        else
           binout   = bin2
           valueout = value2
           errorout = error2
        endif
     else
        binout   = bin1
        valueout = value1
        errorout = error1
     endif
     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)") nhisto1, &
          dummy1, nbin1, dummy1, binout,dummy1,valueout,dummy1,errorout,dummy1
  
 
  enddo
  close(11)
  close(12)
  close(13)

35 continue
end subroutine mysummer_mrg

! get a kfactor. Order of inputs: LO,NLO
  subroutine mysummer_kfactor()
    implicit none
    character *(1) :: dummy1,dummy2
  character *(100) :: file_lo,file_nlo, file_out, charnhisto, charminval, charmaxval
  integer :: i, NumArgs, n_in_files,nhisto1,nhisto2,nbin1,nbin2, thehisto
  real(8) :: bin1,bin2,binout,value1,value2,valueout,error1,error2,errorout, num, minval, maxval,kfactor,errk
  character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)" 

  n_in_files = (COMMAND_ARGUMENT_COUNT()+1)-4

  write(6,*) 'Getting kfactor'
  print *, n_in_files
  

  call GetArg(2,charnhisto) !-- which histogram to pick out
  read(charnhisto,*) thehisto

  call GetArg(3,file_lo)
  file_lo = adjustl(trim(file_lo))
  print *, 'file in 1: ', file_lo

  call GetArg(4,file_nlo)
  file_nlo = adjustl(trim(file_nlo))
  print *, 'file in 2: ', file_nlo
  close(12)
  close(13)
  open (unit=12,file=file_lo,status='old')
  open (unit=13,file=file_nlo,status='old')
   
  call GetArg(5,file_out)
  file_out = adjustl(trim(file_out))
  print *, 'file out: ', file_out
  open (unit=11,file=file_out,status='unknown')

  do while (.true.)
     read(12,*,end=35) dummy1
     read(13,*,end=35) dummy2
     if(dummy1(1:1).eq."#") then
        if (dummy2(1:1) .ne. "#") then
           print *, "error"
           stop
        else
           cycle
        endif
     endif
     backspace(unit=12)
     backspace(unit=13)
     read(12,*,end=35) nhisto1, nbin1, bin1, value1, error1
     read(13,*,end=35) nhisto2, nbin2, bin2, value2, error2
     if(dummy1(1:1).eq."#") then
        if (dummy2(1:1) .ne. "#") then
           print *, "error"
           stop
        else
           cycle
        endif
     endif
     if (nhisto1 .ne. nhisto2 .or. bin1 .ne. bin2 .or. nbin1 .ne. nbin2) then
        print *, "error"
        stop
     endif
     dummy1 = " "
     dummy2 = " "
     if (abs(value1) .lt. 1d-30 .and. abs(value2) .le. 1d-30) then
        kfactor=1d0
        errk=0d0
     else
        kfactor=value2/value1
        errk=value2/value1*sqrt((error1/value1)**2 + (error2/value2)**2)
     endif
     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)") nhisto1, &
          dummy1, nbin1, dummy1, bin1,dummy1,kfactor,dummy1,errk,dummy1
  
 
  enddo
  close(11)
  close(12)
  close(13)

35 continue
end subroutine mysummer_kfactor








! Konstantin (added 21.4.2021 because Fabrizio was not happy)

! subroutine extract_histo()
!     implicit none
!     
!     
!     
!     return
! end subroutine extract_histo

! pure subroutine get_histo_bin( unit, ihisto, nbin, bin, value, error, flag )
!     implicit none
!     
!     integer, intent(in) :: unit
!     integer, intent(out) :: ihisto, nbin
!     real(8), intent(out) :: bin, value, error
!     logical, intent(out) :: flag
!     
!     !-- get full line
!     input = ''
!     do while ( .true. )
!         read(12,'(a)', end = 50, advance = 'no', iostat = ios, size = buffer_length ) buffer
!         input = input // buffer(1:buffer_length)
!         if ( is_iostat_eor(ios) ) exit
!     end do
!     
!     !-- pass to output file if comment
!     if ( is_comment(input) ) then
!         write(11,'(a)') input
!         cycle
!     end if
!     
!     !-- read input as histogram bin
!     backspace(unit=12)
!     read(12,*, end = 50) ihisto, nbin(1), bin(1), value, error
!     
!     
!     return
! end subroutine get_histo_bin


pure function is_comment(string)
    implicit none
    
    logical :: is_comment
    character(len=*), intent(in) :: string
    integer :: n, m
    n = len(string)
    
    is_comment = .false.
    do m = 1, n, 1
        if ( string(m:m) .eq. ' ' ) cycle
        if ( string(m:m) .eq. '#' ) exit
        return
    end do
    
    is_comment = .true.
    return
end function is_comment

! pure function is_id(string)
!     implicit none
!     
!     logical :: is_id
!     character(len=*), intent(in) :: string
!     integer :: n, m
!     n = len(string)
!     
!     is_id = .false.
!     do m = 1, n, 1
!         if ( string(m:m) .eq. ' ' ) cycle
!         if ( string(m:m) .eq. '#' ) exit
!         return
!     end do
!     
!     if ( n - m .lt. 16 ) return
!     if ( string(m:m + 15) .eq. '# ----------Obs=' ) is_id = .true.
!     
!     return
! end function is_id

pure function get_id(string)
    implicit none
    
    character(len=15) :: get_id
    character(len=*), intent(in) :: string
    integer :: n, m
    n = len(string)
    
    get_id = ''
    do m = 1, n, 1
        if ( string(m:m) .eq. ' ' ) cycle
        if ( string(m:m) .eq. '#' ) exit
        return
    end do
    
    if ( n - m .lt. 16 ) return
    if ( string(m:m + 15) .eq. '# ----------Obs=' ) get_id = string(m + 16:m + 30)
    
    return
end function get_id


subroutine move_curser_to_next_histogram( file_unit, name, flag )
    implicit none
    
    integer, intent(in) :: file_unit
    character(len=*), intent(out) :: name
    logical, intent(out) :: flag
    
    integer, parameter :: mxlength = 300
    character(len=mxlength) :: line
    integer :: n
    
    do while ( .true. )
        read( file_unit, '(a)', end = 101 ) line
        do n = 1, mxlength, 1
            if ( line(n:n) .eq. ' ' ) cycle
            if ( line(n:n) .eq. '#' ) exit
            goto 102
        end do
        if ( mxlength - n .lt. 16 ) cycle
        if ( line(n:n + 15) .eq. '# ----------Obs=' ) name = line(n + 16:n + 30)
    end do
    
101 continue
    flag = .true.
    return
    
102 continue
    flag = .false.
    backspace( file_unit )
    return
end subroutine move_curser_to_next_histogram




subroutine MOVE_TO_BEGINNING( file_unit, nhisto, flag )
    implicit none
    
    integer, intent(in) :: file_unit, nhisto
    logical, intent(out) :: flag
    
    character(len=120) :: buffer
    character(len=:), allocatable :: input
    integer :: buffer_length, ios
    logical :: is_comment

    integer :: ihisto, nbin
    real(8) :: bin, value, error
    
    !-- go to beginning of file
    REWIND(unit=file_unit)
    
    do while (.true.)
        !-- get full line
        input = ''
        do while ( .true. )
            read(12,'(a)', end = 50, advance = 'no', iostat = ios, size = buffer_length ) buffer
            input = input // buffer(1:buffer_length)
            if ( is_iostat_eor(ios) ) exit
        end do
        
        !-- cycle if comment
        if ( is_comment(input) ) cycle
        
        !-- read input as histogram bin
        backspace(unit=file_unit)
        read(file_unit,*) ihisto, nbin, bin, value, error
        
        !-- cycle if wrong histo
        if ( ihisto .ne. nhisto ) cycle

        !-- return if correct histo
        BACKSPACE(unit=file_unit)
        flag = .false.
        return
    end do
    
50  continue
    flag = .true.
    return
end subroutine MOVE_TO_BEGINNING





















subroutine print_logo()
    implicit none
    return
    print *, ' ______'
    print *, '|  __  |'
    print *, ' \ \ |_|  SUMmer - manipulation of histograms'
    print *, '  \ \'
    print *, '  / / _   by Konstantin Asteriadis (konstantin.asteriadis@t-online.de)'
    print *, ' / /_| |  Version: beta.2021.09.01'
    print *, '|______|'
    return
end subroutine



















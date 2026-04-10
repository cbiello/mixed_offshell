program mysummer
  implicit none
  character *(3) :: operation

  call GetArg(1,operation)

  if (operation.ne.'add') then
     write(6,*) 'usage: ./mysummer.x add filein_1 filein_2 ... filein_n file_out'
     stop
  endif

  if (operation.eq.'add') call mysummer_add_scalvar()

end program mysummer

!-- sum n histograms including scale variations
subroutine mysummer_add_scalvar()
  implicit none
  character*(1) :: dummy
  character*(140) :: file_in(80), file_out
  integer :: i, j, n_in_files, nhisto, nbin
  real(8) :: bin(200)
  real(8) :: myvalue(3,200), error(3,200)
  real(8) :: sumvalue(3), sumerror(3)

  n_in_files = (COMMAND_ARGUMENT_COUNT()+1)-3

  write(6,*) 'Adding together files'
  do i = 1, n_in_files
     call GetArg(i+1,file_in(i))
     file_in(i) = trim(file_in(i))
     write(*,"(A5,I3,A2,A120)") 'file ', i,': ', file_in(i)
     open(unit=21+i,file=file_in(i),status='old')
  enddo

  call GetArg(n_in_files+2,file_out)
  file_out = trim(file_out)
  print *, 'file out: ', file_out

  open(unit=21,file=file_out,status='unknown')

  do while (.true.)

     sumvalue = 0d0
     sumerror = 0d0

     do i = 1, n_in_files
        read(21+i,*,end=10) dummy
        if(dummy(1:1).eq."#") cycle
        backspace(unit=21+i)

        read(21+i,*,end=10) nhisto, nbin, bin(i), &
             myvalue(1,i), error(1,i), &
             myvalue(2,i), error(2,i), &
             myvalue(3,i), error(3,i)

        do j = 1,3
           sumvalue(j) = sumvalue(j) + myvalue(j,i)
           sumerror(j) = sumerror(j) + error(j,i)**2
        enddo

     enddo

     if(dummy(1:1).eq."#") cycle

     do j = 1,3
        sumerror(j) = dsqrt(sumerror(j))
     enddo

     dummy = " "

     write(21,"(I2,A,I3,A,2X,1PE10.3,A,6(2X,1PE23.16,A))") &
          nhisto, dummy, nbin, dummy, bin(1), dummy, &
          sumvalue(1), dummy, sumerror(1), dummy, &
          sumvalue(2), dummy, sumerror(2), dummy, &
          sumvalue(3), dummy, sumerror(3), dummy

  enddo

10 return

end subroutine mysummer_add_scalvar

!-- sum n histograms
subroutine mysummer_add()
  implicit none
  character *(1) :: dummy
  character *(140) :: file_in(80), file_out
  integer :: i, NumArgs, n_in_files, nhisto, nbin
  real(8) :: bin(200), myvalue(200), error(200)
  real(8) :: sumvalue, sumerror
  character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  

  n_in_files = (COMMAND_ARGUMENT_COUNT()+1)-3

  write(6,*) 'Adding together files'
  do i = 1, n_in_files
     call GetArg(i+1,file_in(i))
     file_in(i) = trim(file_in(i))
     write(*,fmt="(A5,I3,A2,A120)"), 'file ', i,': ', file_in(i)
     open (unit=21+i,file=file_in(i),status='old')
  enddo


  
  call GetArg(n_in_files+2,file_out)
  file_out = trim(file_out)
  print *, 'file out: ', file_out

  open (unit=21,file=file_out,status='unknown')

  do while (.true.)
     sumvalue = 0d0
     sumerror = 0d0
     do i = 1, n_in_files
        read(21+i,*,end=10) dummy
        if(dummy(1:1).eq."#") cycle
        backspace(unit=21+i)
        read(21+i,*,end=10), nhisto, nbin, bin(i),myvalue(i),error(i)
        sumvalue = sumvalue + myvalue(i)
        sumerror = sumerror + error(i)**2
     enddo
     if(dummy(1:1).eq."#") cycle

     sumerror = dsqrt(sumerror)
     dummy = " "
     write(21,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), nhisto, dummy, nbin, &
          dummy, bin(1),dummy,sumvalue,dummy,sumerror,dummy

  enddo

10 return

end subroutine mysummer_add



!-- sum n-1 histograms from first one
subroutine mysummer_sbt()
  implicit none
  character *(1) :: dummy
  character *(120) :: file_in(80), file_out
  integer :: i, NumArgs, n_in_files, nhisto, nbin
  real(8) :: bin(200), myvalue(200), error(200)
  real(8) :: sumvalue, sumerror,isign
  character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  

  n_in_files = (COMMAND_ARGUMENT_COUNT()+1)-3

  write(6,*) 'Subtracting files from first file'
  do i = 1, n_in_files
     call GetArg(i+1,file_in(i))
     file_in(i) = trim(file_in(i))
     write(*,fmt="(A5,I3,A2,A70)"), 'file ', i,': ', file_in(i)
     open (unit=21+i,file=file_in(i),status='old')
  enddo


  
  call GetArg(n_in_files+2,file_out)
  file_out = trim(file_out)
  print *, 'file out: ', file_out

  open (unit=21,file=file_out,status='unknown')

  do while (.true.)
     sumvalue = 0d0
     sumerror = 0d0
     do i = 1, n_in_files
        if (i .eq. 1) isign=1d0
        if (i .gt. 1) isign=-1d0
        read(21+i,*,end=10) dummy
        if(dummy(1:1).eq."#") cycle
        backspace(unit=21+i)
        
        read(21+i,*,end=10), nhisto, nbin, bin(i),myvalue(i),error(i)
        sumvalue = sumvalue + isign*myvalue(i)
        sumerror = sumerror + error(i)**2
     enddo
     if(dummy(1:1).eq."#") cycle

     sumerror = dsqrt(sumerror)
     dummy = " "
     write(21,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), nhisto, dummy, nbin, &
          dummy, bin(1),dummy,sumvalue,dummy,sumerror,dummy
     
  enddo

10 return

end subroutine mysummer_sbt



!-- multiply by a constant histograms
subroutine mysummer_mul()
  implicit none
  character *(1) :: dummy
  character *(100) :: file_in, file_out, chnum
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
     read(12,*,end=35), nhisto, nbin, bin, value, error
     value = value * num
     error = error * num
     if(dummy(1:1).eq."#") cycle
     dummy = " "
     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), nhisto, &
          dummy, nbin, dummy, bin,dummy,value,dummy,error,dummy
 
  enddo

35 continue

end subroutine mysummer_mul

!-- average n histograms
subroutine mysummer_avg()
  implicit none
  character *(1) :: dummy
  character *(150) :: file_in(1000), file_out
  integer :: i, NumArgs, n_in_files, nhisto, nbin, icount
  real(8) :: bin(1000), value(1000), error(1000)
  real(8) :: sumvalue, sumerror
  real(8) :: sumvalue2, sumerror2, sumerror3,sumerror4
  character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  

  n_in_files = (COMMAND_ARGUMENT_COUNT()+1)-3

  write(6,*) 'Averaging files'

  do i = 1, n_in_files
     call GetArg(i+1,file_in(i))
     file_in(i) = trim(file_in(i))
     write(*,fmt="(A5,I2,A2,A120)"), 'file ', i,': ', file_in(i)
     open (unit=11+i,file=file_in(i),status='old')
  enddo
  
  call GetArg(n_in_files+2,file_out)
  file_out = trim(file_out)
  print *, 'file out: ', file_out

  open (unit=11,file=file_out,status='unknown')

  do while (.true.)
     sumvalue = 0d0
     sumerror = 0d0
     sumvalue2 = 0d0
     sumerror2 = 0d0
     sumerror3 = 0d0
     sumerror4 = 0d0
     icount = 0
     do i = 1, n_in_files
        read(11+i,*,end=20) dummy
        if(dummy(1:1).eq."#") cycle
        backspace(unit=11+i)
        read(11+i,*,end=20), nhisto, nbin, bin(i),value(i),error(i)
        icount = icount + 1
        !-- normal combination
        sumvalue = sumvalue + value(i)
        sumerror = sumerror + error(i)**2
        !-- gaussian combination
        sumvalue2 = sumvalue2 + value(i)/error(i)**2
        sumerror2 = sumerror2 + 1d0/error(i)**2
        !-- error is just the spread, gaussian
        sumerror3 = sumerror3 + value(i)**2/error(i)**2
     enddo
     if(dummy(1:1).eq."#") cycle

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

     dummy= " "
     if (sumvalue.eq.0d0 .and. sumerror2 .gt. 1000d0) sumerror2 = 0d0
     if (sumvalue.eq.0d0 .and. sumerror3 .gt. 1000d0) sumerror3 = 0d0
     if (sumvalue.lt.1d-50) sumvalue = 0d0
     if (.not.sumvalue2.le.0d0 .and. .not.sumvalue2.gt.0d0 .and. sumerror2.eq.0d0) sumvalue2 = 0d0 !-- remove NaN from 0 error
     !-- gaussian
     if (n_in_files.le.5) then
!        write(11,fmt="(I2,A,X,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), nhisto, dummy, &
        !             nbin, dummy, bin(1),dummy,sumvalue2,dummy,sumerror2,dummy
        write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), nhisto, dummy, &
             nbin, dummy, bin(1),dummy,sumvalue2,dummy,sumerror2,dummy
     else
!       write(11,fmt="(I2,A,X,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), nhisto, dummy, &
        !            nbin, dummy, bin(1),dummy,sumvalue2,dummy,sumerror2,dummy
        write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), nhisto, dummy, &
             nbin, dummy, bin(1),dummy,sumvalue2,dummy,sumerror2,dummy
     endif

  enddo

20 return

end subroutine mysummer_avg


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
     read(12,*,end=30), nhisto, nbin, bin, value, error
     if(dummy(1:1).eq."#") cycle
     dummy = " "
     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), nhisto, &
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
     write(*,fmt="(A5,I2,A2,A40)"), 'file ', i,': ', file_in(i)
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
        read(11+i,*,end=40), nhisto, nbin, bin(i),value(i),error(i)
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
     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), nhisto, dummy, &
          nbin, dummy, bin(1),dummy,sumvalue,dummy,sumerror,dummy
  enddo

40 continue

end subroutine mysummer_pls


!-- rebin n histograms
subroutine mysummer_rebin()
  implicit none
  character *(1) :: dummy
  character *(100) :: file_in, file_out, charnrebin, charnhisto
  integer :: i, NumArgs, ihisto, nhisto, nrebin, nbin(100), ibin, ihisto_dum, ihisto_tmp
  real(8) :: bin(100), value, error
  real(8) :: sumvalue, sumerror
  integer, parameter :: tothisto = 4
  character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  

  write(6,*) 'Rebinning'

  call GetArg(2,charnhisto) !-- which histogram to rebin (-1 -> all of them)
  read(charnhisto,*) nhisto
  call GetArg(3,charnrebin) !-- how many bin to combine
  read(charnrebin,*) nrebin
  call GetArg(4,file_in) 
  call GetArg(5,file_out)

  file_in = trim(file_in)
  file_out = trim(file_out)

  print *, 'number of bins to combine: ', nrebin
  print *, 'file in: ', file_in
  print *, 'file out: ', file_out
  print *, 'nhisto: ', nhisto
  print *, 'nrebin: ', nrebin

  open (unit=12,file=file_in,status='old')
  open (unit=11,file=file_out,status='unknown')

  nrebin = nrebin-1

  do while (.true.)
     read(12,*,end=50) dummy
     if(dummy(1:1).eq."#") cycle
     backspace(unit=12)

     read(12,*,end=50), ihisto, nbin(1), bin(1),value,error
     if (ihisto.eq.nhisto) then

        sumvalue = value
        sumerror = error**2
        do ibin = 1, nrebin
           read(12,*,end=50), ihisto_tmp, nbin(ibin+1), bin(ibin+1),value,error
           ! this seems to help if the total number of bins is not a multiple of the number of bins that we're rebinning
           if (ihisto_tmp .ne. nhisto) then
              backspace(unit=12)
              cycle
           endif
           ihisto = ihisto_tmp
           if(dummy(1:1).eq."#") cycle

           sumvalue = sumvalue + value
           sumerror = sumerror + error**2
        enddo
! RR modified to remove dividing by number of bins
        sumvalue = sumvalue         !/(nrebin+1d0)
        sumerror = sumerror         !/(nrebin+1d0)
        sumerror = dsqrt(sumerror)
     else
        sumvalue = value
        sumerror = error
     endif
     
     if(dummy(1:1).eq."#") cycle

     dummy = " "

     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), ihisto, &
          dummy, nbin(1), dummy, bin(1),dummy,sumvalue,dummy,sumerror,dummy

  enddo

50 return

end subroutine mysummer_rebin


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

     read(12,*,end=60), ihisto, nbin(ibin), bin(ibin),value,error
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
  write(*,*) 'ratio w last bin/wo last bin', sumvalue/(sumvalue-lastvalue)
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
     read(12,*,end=35), nhisto, nbin, bin, value, error
     if(dummy(1:1).eq."#") cycle
     dummy = " "
     if (nhisto .eq. thehisto) then
        write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), nhisto, &
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
     read(12,*,end=35), nhisto, nbin, bin, value, error
     if(dummy(1:1).eq."#") cycle
     dummy = " "
     if (bin .lt. minval .or. bin .gt. maxval) cycle
     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), nhisto, &
          dummy, nbin, dummy, bin,dummy,value,dummy,error,dummy
  
 
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
     read(12,*,end=35), nhisto1, nbin1, bin1, value1, error1
     read(13,*,end=35), nhisto2, nbin2, bin2, value2, error2
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
     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), nhisto1, &
          dummy1, nbin1, dummy1, binout,dummy1,valueout,dummy1,errorout,dummy1
  
 
  enddo
  close(11)
  close(12)
  close(13)

35 continue

end subroutine mysummer_mrg

! -- produces the same histogram mirrored about x=0
subroutine mysummer_mrr()
  implicit none
  character *(1) :: dummy
  character *(100) :: file_in, file_out, charnhisto,charnbins
  integer :: i, NumArgs, n_in_files, nhisto, nbin, thehisto,binnum,mynhisto(1:300),mynbin(1:300),cntr
  real(8) :: mybin(1:300),myvalue(1:300), myerror(1:300)
  real(8) :: bin, value, error, num
  character(len=*),parameter :: fmt1 = "(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"  

  n_in_files = (COMMAND_ARGUMENT_COUNT()+1)-5

  if (n_in_files.ne.1) then
     write(6,*) 'only one file at the time in input'
     stop
  endif

  write(6,*) 'Pick out a given histogram'

  call GetArg(2,charnhisto) !-- which histogram to pick out
  read(charnhisto,*) thehisto

  call GetArg(3,charnbins) !-- which histogram to pick out
  read(charnbins,*) binnum

  call GetArg(4,file_in)
  file_in = adjustl(trim(file_in))
  print *, 'file in: ', file_in
  
  open (unit=12,file=file_in,status='old')
  
  call GetArg(5,file_out)
  file_out = adjustl(trim(file_out))
  print *, 'file out: ', file_out

  open (unit=11,file=file_out,status='unknown')

  do cntr=1,binnum
     read(12,*,end=35), mynhisto(cntr), mynbin(cntr), mybin(cntr), myvalue(cntr), myerror(cntr)
  enddo

!  do while (.true.)
!     read(12,*,end=35) dummy
!     if(dummy(1:1).eq."#") cycle
!     backspace(unit=12)
!     read(12,*,end=35), nhisto, nbin, bin, value, error
!     if(dummy(1:1).eq."#") cycle
!     if (bin .eq. 0) cycle
!     dummy = " "
!     if (nhisto .eq. thehisto) then
!        write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), nhisto, &
!             dummy, nbin, dummy, -bin,dummy,value,dummy,error,dummy
!     endif
! 
!  enddo

  dummy = " "
  do cntr=binnum,2,-1
     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"), mynhisto(cntr), &
          dummy, mynbin(cntr), dummy, -mybin(cntr),dummy,myvalue(cntr),dummy,myerror(cntr),dummy
  enddo
  
  do  cntr=1,binnum
     write(11,fmt="(I2,A,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)"),mynhisto(cntr), &
          dummy, mynbin(cntr), dummy, mybin(cntr),dummy,myvalue(cntr),dummy,myerror(cntr),dummy
  enddo
  

35 continue

end subroutine mysummer_mrr


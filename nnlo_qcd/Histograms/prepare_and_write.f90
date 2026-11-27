  !-- prepare histograms and write them down to file

  do i=1,Nobs*Numofblocks
     do j=1,Histo(i)%Nbins+1

        HistoGlobalArray_tmp(i,j) = HistoGlobalArray(i,j)
        ErrorGlobalArray_tmp(i,j) = ErrorGlobalArray(i,j)

        Error_aux(i,j) = sqrt(Intaversq(i,j) - Intaver(i,j)**2/evalcount)  

        if (Error_aux(i,j).eq.zero) Error_aux(i,j) = ErrorGlobalArray(i,j) 

        HistoGlobalArray_tmp(i,j) = (Intaver(i,j)/Error_aux(i,j)**2  & 
             + HistoGlobalArray_tmp(i,j)/ErrorGlobalArray_tmp(i,j)**2 )/(one/Error_aux(i,j)**2+one/ErrorGlobalArray_tmp(i,j)**2)

        ErrorGlobalArray_tmp(i,j) = one/sqrt(one/ErrorGlobalArray_tmp(i,j)**2 + one/Error_aux(i,j)**2)

     enddo
  enddo
  
  call get_nan(HistoGlobalArray_tmp(1,1),flag_nan)
  if (flag_nan) return

  !--------------------------------------------------------------------------------

  !-- write down our histograms
  call writeHisto(outputfile,HistoGlobalArray_tmp,ErrorGlobalArray_tmp,myiternumber)
  print *, "Output file: ", outputfile


  !-- intermediate runs
  write(ch_iternumber,'(I3)') myiternumber
!  tmpfile = trim(adjustl(outputfile))//trim(adjustl(ch_iternumber))//runtype
!  print *, 'Temporary file: ', tmpfile
!  call writeHisto(tmpfile,HistoGlobalArray_tmp,ErrorGlobalArray_tmp,myiternumber,.true.)


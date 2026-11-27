  !-- check that the total integral agrees, as well as sum of bins of all integrals
  print *, 'evalcount ', evalcount - 1
  print *, 'HistoGlobalArray_tmp(1,1)', HistoGlobalArray_tmp(1,1)
  print *, 'ErrorGlobalArray_tmp(1,1)', ErrorGlobalArray_tmp(1,1)
  print *, ''
  print *, 'Sum of bins, non-cumulative'
  do j = 1, Nnoncum
     mytest = zero
     mytest_err = zero
     do i = 1,Histo(j)%nbins+1
        call get_nan(HistoGlobalArray_tmp(j,i),flag_nan)
        if (.not. flag_nan) then
           mytest = mytest + HistoGlobalArray_tmp(j,i)
           if (HistoGlobalArray_tmp(j,i).ne.zero) mytest_err = mytest_err + ErrorGlobalArray_tmp(j,i)**2
        endif
     enddo
     write(*,'(A,I2,0PE14.4,A,0PE14.4)') ' histo', j, mytest,' pm',sqrt(mytest_err)
  enddo
  print *, ''
  print *, 'Cumulative, first bin'
  do j = Nnoncum+1, Nobs
     write(*,'(A,I2,0PE14.4,A,0PE14.4)') ' Chisto', j, HistoGlobalArray_tmp(j,1), &
           ' pm',ErrorGlobalArray_tmp(j,1)
  enddo

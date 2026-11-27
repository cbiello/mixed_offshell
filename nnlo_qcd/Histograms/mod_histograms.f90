module mod_histograms
  use mod_types; use mod_consts_dp; use common_def
  use mod_auxfunctions
  use mod_kinfunctions
  implicit none
#include "include_nobs.f90"
  !
  integer, public, parameter :: NumOfBlocks = 1  ! replicas of histograms, to check pole cancellation. Here we don't do it
  integer, public, parameter :: Nbinmax = 300  ! maximal number of bins 
  !
  integer, public :: evalcount ! count the number of evaluations
  integer, public :: sumcount = 0
  integer, public :: evalcountmax  ! maximum number of evaluations 
  
  real(dp), public :: auxweight = one ! a useful parameter to account for 
                                      ! some global factors that are used 
                                      ! in finalresults instead of FLO
  !
  real(dp), public :: HistoGlobalArray(Nobs*NumOfBlocks,Nbinmax) = zero ! second argument is just large enough 
  real(dp), public :: HistoGlobalArray_aux(Nobs*NumOfBlocks,Nbinmax) = zero  ! second argument is just large enough 
  real(dp), public :: ErrorGlobalArray(Nobs*NumOfBlocks,Nbinmax) = 1d9 ! second argument is just large enough 
  real(dp), public :: ErrorGlobalArray_aux(Nobs*NumOfBlocks,Nbinmax) ! second argument is just large enough 
  !
  real(dp), public :: Intaver(Nobs*NumOfBlocks,Nbinmax)  ! second argument is just large enough 
  real(dp), public :: Intaver_aux(Nobs*NumOfBlocks,Nbinmax)  ! second argument is just large enough 
  !
  real(dp), public :: Intaversq(Nobs*NumOfBlocks,Nbinmax)  ! second argument is just large enough 
  real(dp), public :: Intaversq_new(Nobs*NumOfBlocks,Nbinmax)  ! second argument is just large enough 
  real(dp), public :: Error_aux(Nobs*NumOfBlocks,Nbinmax)  ! second argument is just large enough 
  !
  real(dp), public :: HistoGlobalArray_Markus(Nobs*NumOfBlocks,Nbinmax) = zero ! second argument is just large enough 

  type :: Histogram
     character *15 :: obs_code
     real(dp) :: LowLimit,UpLimit,Step
     integer :: Nbins
  end type Histogram

  type(Histogram), public :: histo(Nobs)

  private

  public :: init_histo_w
  public :: cut_histo_w
  public :: find_bin,HistoFill
  public :: writeHisto
  
contains

#include "./cut_histo_w.f90"

  subroutine find_bin(obs,bin_array,weight_array)
    real(dp), intent(in) :: obs(nobs)
    integer,  intent(out) :: bin_array(Nobs)
    real(dp), intent(out) :: weight_array(Nobs)
    real(dp) :: binlow,binhigh,hstep
    integer :: i,j,nbins_local

    bin_array = 0
    weight_array = one

    !----- now finding the bins where we should put them

    !--- first do normal histograms
    do i=1,Nnoncum
       
       hstep = Histo(i)%Step
       nbins_local  = Histo(i)%Nbins
       binlow = Histo(i)%LowLimit

       j = (Obs(i)-Histo(i)%LowLimit)/Histo(i)%Step + 1

       if (j.gt.nbins_local + 1) j = nbins_local + 1 !-- overflow
       if (j.le.0) j = nbins_local + 1 !-- underflow in the overflow for now

       bin_array(i) = j

    enddo

    !--- then do the cumulative ones
    do i=Nnoncum+1,Nobs

       hstep = Histo(i)%Step
       nbins_local  = Histo(i)%Nbins
       binlow = Histo(i)%LowLimit

       do j=1,nbins_local
              
          binhigh = binlow + hstep
          
          if ( (Obs(i).lt.binhigh))  then !-- cumulative histogram
             bin_array(i) = j
             exit 
          endif

          binlow = binhigh
          
       enddo

       !-- overflow bin
       if (bin_array(i).eq.0) bin_array(i) = nbins_local + 1

    enddo
    
  end subroutine find_bin

  subroutine HistoFill(STATUS,bin_array,weight_array,res,vegasweight,iternumber,histoblock) 
    integer, intent(in) ::  bin_array(Nobs)
    real(dp), intent(in) :: weight_array(Nobs)
    integer, intent(in), optional :: histoblock
    integer :: iternumber
    real(dp15) :: vegasweight
    real(dp), intent(in) :: res 
    integer, intent(in) :: STATUS
    integer :: i,j,ishift
    real(dp) :: HistoGlobalArray_tmp(Nobs*NumOfBlocks,NbinMax),ErrorGlobalArray_tmp(Nobs*NumOfBlocks,NbinMax)
    character :: tmpfile*(100),ch_iternumber*(3)
    real(dp) :: mytest,mytest_err
    integer :: myiternumber = 1
    integer :: iskip = 20
    logical :: flag_nan

    call get_nan(res*real(vegasweight,kind=dp)*auxweight,flag_nan)
    if (flag_nan) return

    !-- shift to the right histogram block
    if (.not. present(histoblock)) then
       ishift = 0
    else
       ishift = Nobs * (histoblock - 1)
    endif

    if (evalcount.eq.1.and.STATUS.eq.0) then 
       Intaver = zero 
       Intaversq = zero
       Intaversq_new = zero
       Intaver_aux = zero 
    endif
    
    !-- use ishift to go to the correct block
    !-- first do the non-cumulative ones
    do i=1,Nnoncum
       if (bin_array(i).ne.0) then 
          Intaver_aux(i+ishift,bin_array(i)) = Intaver_aux(i+ishift,bin_array(i))  & !--standard
               + res*real(vegasweight,kind=dp)*auxweight*weight_array(i)   ! res/evalcountmax*auxweight
       endif
    enddo

    !-- then do the cumulative
    do i=Nnoncum+1,Nobs
       if (bin_array(i).ne.0) then 
          do j = 1, bin_array(i) !-- cumulative
             Intaver_aux(i+ishift,j) = Intaver_aux(i+ishift,j)  & !--cumulative
                  + res*real(vegasweight,kind=dp)*auxweight*weight_array(i)   ! res/evalcountmax*auxweight
          enddo
        endif
    enddo

    if (STATUS.eq.10) then 
            Intaversq = Intaversq + Intaver_aux**2/auxweight 
            Intaver = Intaver + Intaver_aux
    endif 

    if (STATUS.eq.10) evalcount = evalcount + 1 
    if (STATUS.eq.10) Intaver_aux = zero 

  end subroutine HistoFill

  subroutine WriteHisto(TheFile,HistoArray,ErrorArray,iternumber,DontCopyButMove)
    character*100, intent(in) :: TheFile
    real(dp), intent(in) :: HistoArray(Nobs*NumOfBlocks,Nbinmax), ErrorArray(Nobs*NumOfBlocks,Nbinmax)
    integer, intent(in) :: iternumber
    logical, intent(in), optional :: DontCopyButMove
    integer :: i, j, k, nbins_local, iunit
    real(dp) :: hstep, binlow
    character*1500 :: mystring
    
    
    open(unit=12,file=TheFile)
    iunit = 12
#include "./../printout.f90"
#include "./../printout_proc.f90"
    
    do k = 1, NumOfBlocks

       write(12,*) '# ----------Iteration number ',iternumber,runtype,'----------------'   
       write(12,*)  '# ----------Block = ',k,'----------------'   


       do i=1,Nobs   

          write(12,*) '# ----------Obs=',Histo(i)%Obs_code,'----------------'   
          write(12,*) '# bin number--', 'value of the variable--', 'value of the function--'
          write(12,*) '# last bin is overflow '
          
          hstep = Histo(i)%Step
          nbins_local  = Histo(i)%Nbins+1
          binlow = Histo(i)%LowLimit
          
          do j=1,nbins_local
             
             write(12,"(I2,A,X,I3,A,2X,1PE10.3,A,2X,1PE23.16,A,2X,1PE23.16,A)") i+(k-1)*nobs,' ', & 
                  j, ' ',binlow+hstep*(j-1),' ',HistoArray(i+(k-1)*nobs,j),' ', ErrorArray(i+(k-1)*nobs,j),' '
          
          enddo
       
       enddo
    
    enddo

    include "dump_const.f90"
    include "dump_cut.f90"

    close(12)

    !-- for CERN cluster, move histogram files to a fixed localtion during run-time
    if (histotmpfile.ne.'') then

       if ( present(DontCopyButMove) ) then
          if (DontCopyButMove) then
             mystring = "mv "//trim(adjustl(thefile))//" "//trim(adjustl(histotmpfile))//"/."
             call system(mystring)
          else
             mystring = "cp "//trim(adjustl(thefile))//" "//trim(adjustl(histotmpfile))//"/."
             call system(mystring)
          endif
       else
          mystring = "cp "//trim(adjustl(thefile))//" "//trim(adjustl(histotmpfile))//"/."
          call system(mystring)
       endif

    endif

  end subroutine WriteHisto

end module mod_histograms



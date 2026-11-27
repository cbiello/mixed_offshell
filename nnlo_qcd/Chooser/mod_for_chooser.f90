module mod_for_chooser
  use mod_types; use common_def; use mod_consts_dp
  use mod_consts_debug
  implicit none
  private

  public :: read_parms_local, read_parms_gen

contains

  subroutine read_parms_local(numargs,countarg,corr,pdford,sec,ch,runname)
    integer, intent(in) :: numargs
    integer :: countarg, i
    character :: arg*(300)
    character, intent(out) :: corr*(8),pdford*(4),sec*(10),ch*(10),runname*(70)

    do i = 1, numargs

       call getarg(i,arg)

       !-- correction, sectors...

       if (arg(1:5).eq.'corr=') then
          read(arg(6:13),*) corr
          corr=trim(adjustl(corr))
          countarg = countarg + 1
       elseif (arg(1:7).eq.'pdford=') then
          read(arg(8:11),*) pdford
          pdford=trim(adjustl(pdford))
          countarg = countarg + 1
       elseif (arg(1:4).eq.'sec=') then
          read(arg(5:14),*) sec
          sec=trim(adjustl(sec))
          countarg = countarg + 1
       elseif (arg(1:3).eq.'ch=') then
          read(arg(4:13),*) ch
          ch=trim(adjustl(ch))
          countarg = countarg + 1
       elseif(arg(1:8).eq.'runname=') then
          read(arg(9:68),*) runname
          runname=trim(adjustl(runname))
          countarg = countarg + 1

       endif

    enddo

  end subroutine read_parms_local

  subroutine read_parms_gen(numargs,countarg)
    integer, intent(in) :: numargs
    integer :: countarg, i
    character :: arg*(300)
    real(dp) :: mu,myqmin,myqmax

    do i = 1, numargs

       call getarg(i,arg)

       !-- vegas parameters
       if(arg(1:9) .eq.'vegasNc0=') then
          read(arg(10:20),*) vegasNc0
          countarg = countarg + 1
       elseif(arg(1:9) .eq.'vegasNc1=') then
          read(arg(10:20),*) vegasNc1
          countarg = countarg + 1
       elseif(arg(1:9) .eq.'vegasIt0=') then
          read(arg(10:20),*) vegasIt0
          countarg = countarg + 1
       elseif(arg(1:9) .eq.'vegasIt1=') then
          read(arg(10:20),*) vegasIt1
          countarg = countarg + 1
       elseif(arg(1:5) .eq.'seed=') then
          read(arg(6:20),*) seed
          countarg = countarg + 1

       elseif(arg(1:3).eq.'mu=') then
          read(arg(4:12),*) mu
          musq_in = mu**2
          countarg = countarg + 1
       elseif(arg(1:5).eq.'xMuR=') then
          read(arg(6:18),*) xMuR
          countarg = countarg + 1
       elseif(arg(1:5).eq.'xMuF=') then
          read(arg(6:18),*) xMuF
          countarg = countarg + 1
       elseif(arg(1:5).eq.'qmin=') then
          read(arg(6:12),*) myqmin
          countarg = countarg + 1
          myq2min =myqmin**2
       elseif(arg(1:5).eq.'qmax=') then
          read(arg(6:12),*) myqmax
          countarg = countarg + 1
          myq2max =myqmax**2

       !-- buffers
       elseif(arg(1:7).eq.'bufftc=') then
          read(arg(8:12),*) bufftc
          bufftc = 10**(bufftc)
          countarg = countarg + 1
       elseif(arg(1:7).eq.'buffrr=') then
          read(arg(8:12),*) buffrr
          buffrr = 10**(buffrr)
          countarg = countarg + 1
       elseif(arg(1:7).eq.'buffrv=') then
          read(arg(8:12),*) buffrv
          buffrv = 10**(buffrv)
          countarg = countarg + 1
       elseif(arg(1:5).eq.'buff=') then
          read(arg(6:10),*) buff
          buff = 10**(buff)
          onet = one-two*buff
          countarg = countarg + 1

       !-- limit check
       elseif(arg(1:9).eq.'limcheck=') then
          read(arg(10:14),*) limcheck
          limcheck = 10**(limcheck)
          countarg = countarg + 1

       endif

    enddo

    return

  end subroutine read_parms_gen


end module mod_for_chooser

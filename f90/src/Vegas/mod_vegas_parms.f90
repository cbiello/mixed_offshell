module mod_vegas_parms
  use mod_types
  use mod_consts_dp
  use mod_parser
  implicit none
  integer, public, save :: vegasNc0 = 10000
  integer, public, save :: vegasNc1 = 50000
  integer, public, save :: vegasIt0 = 5
  integer, public, save :: vegasIt1 = 10
  integer, public, save :: seed = 0
  !
  character(strlen), public :: gridfile
  logical, save, public :: writegrid,readgrid
  private

  public :: get_vegas_parms, help_vegas_parms, print_vegas_parms

contains

  subroutine get_vegas_parms()
    
    vegasIt0 = int_val_opt('-vegasIt0',5)
    vegasNc0 = int_val_opt('-vegasNc0',10000)
    !
    vegasIt1 = int_val_opt('-vegasIt1',10)
    vegasNc1 = int_val_opt('-vegasNc1',50000)
    !
    seed = int_val_opt('-seed',0)

    gridfile   = trim(string_val_opt('-gridfile','auto'))
    writegrid  = log_val_opt('-writegrid',.false.)
    readgrid   = log_val_opt('-readgrid',.false.)
    if (writegrid .and. readgrid) then
       print *, 'cannot both write and read grid'
       stop
    endif


  end subroutine get_vegas_parms

  subroutine help_vegas_parms(idev)
    integer, intent(in) :: idev

    write(idev,*) ''
    write(idev,*) ' -vegasIt0 5'
    write(idev,*) ' -vegasNc0 10000'
    write(idev,*) ' -vegasIt1 10'
    write(idev,*) ' -vegasNc1 50000'
    write(idev,*) ' -seed 0'
    write(idev,*) ''
    write(idev,*) ' -gridfile auto --> automatic name for the grid file'
    write(idev,*) ' -writegrid false'
    write(idev,*) ' -readgrid false'
    
  end subroutine help_vegas_parms
  
  subroutine print_vegas_parms(outdev)
    integer, intent(in) :: outdev

    write(outdev,*) '# vegasIt0 = ', vegasIt0
    write(outdev,*) '# vegasNc0 = ', vegasNc0
    write(outdev,*) '#'
    write(outdev,*) '# vegasIt1 = ', vegasIt1
    write(outdev,*) '# vegasNc1 = ', vegasNc1
    write(outdev,*) '#'
    write(outdev,*) '# seed = ', seed
    write(outdev,*) '#'
    write(outdev,*) '# gridfile = ', gridfile
    write(outdev,*) '# writegrid = ', writegrid
    write(outdev,*) '# readgrid  = ', readgrid
    write(outdev,*) '#'

  end subroutine print_vegas_parms

end module mod_vegas_parms

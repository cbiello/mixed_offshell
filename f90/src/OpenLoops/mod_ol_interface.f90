module mod_ol_interface
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_mpi_common
  implicit none
  integer,  public :: OL_id(9), OL_rr_id(12)
  real(dp15), public, save :: as_ol,aem_ol,gsq_ol,eesq_ol
  real(dp), public, parameter :: mu_ol = 100._dp
  private

  public :: initialise_ol
  
contains

  subroutine initialise_ol(ol_verb)
    use openloops
    integer, intent(in) :: ol_verb !-- this should be 1 to printout 
    integer :: stab_mode = 11
    integer :: istart_rr

    if (.not. root_process) then
       call set_parameter("no_splash",1)
       call set_parameter("verbose",0)
    else
       call set_parameter("verbose",ol_verb)
       call set_parameter("parameters_write",ol_verb)
    endif

    call set_parameter("expert_mode",1)
    call set_parameter("stability_mode",stab_mode)
    
    !-- Set OL physical parameters
    !-- we use the complex mass scheme for W,Z only
    !-- --> set GaH = GaTop = 0
    call set_parameter("mass(5)",real(mb,dp15))      ! b, by default 0
    call set_parameter("mass(6)",real(mt,dp15))      ! t
    call set_parameter("mass(23)",real(mz_in,dp15))  ! Z
    call set_parameter("mass(24)",real(mw_in,dp15))  ! W
    call set_parameter("mass(25)",real(mh,dp15))     ! H
    
    call set_parameter("width(5)",real(0._dp,dp15))   ! top
    call set_parameter("width(23)",real(GaZ,dp15))    ! Z
    call set_parameter("width(24)",real(GaW,dp15))    ! W
    call set_parameter("width(25)",real(0._dp,dp15))  ! H

    call set_parameter("alphas",real(0.118_dp,dp15))  !-- by default, OL will keep this fixed
    if (cm_scheme) then
       call set_parameter("complex_mass_scheme",1)
    else
       call set_parameter("complex_mass_scheme",0)
    endif
           
    !-- EW input schemes
    if (ew_scheme.eq.'a0') then           !-- alpha(0)
       call set_parameter("ew_scheme", 0) 
       call set_parameter("alpha_qed_0",real(alpha_0,dp15))
       call set_parameter("alpha_qed_mz",real(alpha_mz,dp15))
       call set_parameter("onshell_photons_lsz",0)
       call set_parameter("offshell_photons_lsz",0)
    elseif (ew_scheme == 'Gmu') then
       call set_parameter("ew_scheme", 1) !-- Gmu 
       call set_parameter("Gmu",real(Gf,dp15))
    elseif (ew_scheme == 'amz') then
       call set_parameter("ew_scheme", 2) !-- alpha(mz)
       call set_parameter("alpha_qed_0",real(alpha_0,dp15))
       call set_parameter("alpha_qed_mz",real(alpha_mz,dp15))
    endif

    
    !-- initialise only the processes that you are using, for collier chaching 
    !-- register the various amplitudes
    if (corr(1:3).eq.'nlo' .and. ch.eq.'ns') then

       if (corr.eq.'nloqcd' .and. sec.eq.'v') then
          call register_red_qcd(1)
       elseif (corr.eq.'nloewk' .and. sec.eq.'v') then
          call register_red_ew(1)
       endif

    elseif (corr.eq.'nloewk' .and. sec.eq.'v' .and. ch.eq.'aa') then
       call register_red_ew_aa(1)

    elseif (corr.eq.'nnlo' .and. sec.eq.'vvfc') then
       call register_red_ew(1)

    elseif (corr.eq.'nnlo' .and. sec(1:2).eq.'rv') then

       if (ch.eq.'ns') then
          if (sec .eq.'rvewk_is') then
             call register_red_ew(1)
             call set_orders_g_virtew()
#if (_Vcharge == 0)
             OL_id(4) = register_process("1 -1 -> 11 -11 21",11)   ! d db > e-e+ g  + d <--> db from crossing, so that I keep the same cache
             OL_id(5) = register_process("2 -2 -> 11 -11 21",11)   ! u ub > e-e+ g  + u <--> ub
             OL_id(6) = register_process("5 -5 -> 11 -11 21",11)   ! b bb > e-e+ g -- only needed for EW
#elseif (_Vcharge == 1)
             OL_id(3) = register_process("-1 2 -> 12 -11 21",11)
             OL_id(4) = register_process("2 -1 -> 12 -11 21",11)
             ! b t > absent
#elseif (_Vcharge == -1)
             OL_id(3) = register_process("1 -2 -> 11 -12 21",11)
             OL_id(4) = register_process("-2 1 -> 11 -12 21",11)
             ! b t > absent             
#endif
          elseif (sec .eq.'rvewknf_is') then
             call set_parameter("approximation","ewnf") !-- select 1L-ewk nf amplitudes from OpenLoops
             call register_red_ew(1)
             call set_orders_g_virtew()
#if (_Vcharge == 0)
             OL_id(4) = register_process("1 -1 -> 11 -11 21",11)   ! d db > e-e+ g  + d <--> db from crossing, so that I keep the same cache
             OL_id(5) = register_process("2 -2 -> 11 -11 21",11)   ! u ub > e-e+ g  + u <--> ub
#elseif (_Vcharge == 1)
             OL_id(3) = register_process("-1 2 -> 12 -11 21",11)
             OL_id(4) = register_process("2 -1 -> 12 -11 21",11)
#elseif (_Vcharge == -1)
             OL_id(3) = register_process("1 -2 -> 11 -12 21",11)
             OL_id(4) = register_process("-2 1 -> 11 -12 21",11)
#endif       
          elseif (sec(1:5).eq.'rvqcd') then
             call register_red_qcd(1)
             call set_orders_a_virtqcd()
#if (_Vcharge == 0)
             OL_id(3) = register_process("1 -1 -> 11 -11 22",11)   ! d db > e-e+ a  + d <--> db from crossing, so that I keep the same cache
             OL_id(4) = register_process("2 -2 -> 11 -11 22",11)   ! u ub > e-e+ a  + u <--> ub
#elseif (_Vcharge == 1)
             OL_id(3) = register_process("-1 2 -> 12 -11 22",11)
             OL_id(4) = register_process("1 -2 -> 12 -11 22",11)
#elseif (_Vcharge == -1)
             OL_id(3) = register_process("1 -2 -> 11 -12 22",11)
             OL_id(4) = register_process("-2 1 -> 11 -12 22",11)
#endif
          else
             print *, 'unrecognized option for ol, nnlo, ns'
          endif

       elseif (ch.eq.'gq' .and. sec .eq.'rvewk_is') then
          call register_red_ew(1)
          call set_orders_g_virtew()
#if (_Vcharge == 0)
          OL_id(4) = register_process("21 -1 -> 11 -11 -1",11) !-- ask federico about crossing
          OL_id(5) = register_process("21 -2 -> 11 -11 -2",11)
          OL_id(6) = register_process("21 -5 -> 11 -11 -5",11)
          !
          OL_id(7) = register_process("21 1 -> 11 -11 1",11) !-- ask federico about crossing
          OL_id(8) = register_process("21 2 -> 11 -11 2",11)
          OL_id(9) = register_process("21 5 -> 11 -11 5",11)
#elseif (_Vcharge == 1)
          OL_id(3) = register_process("21 2 -> 12 -11 1",11)
          OL_id(4) = register_process("21 -1 -> 12 -11 -2",11)
          ! g b > W t has a different signature
#elseif (_Vcharge == -1)
          OL_id(3) = register_process("21 -2 -> 11 -12 -1",11)
          OL_id(4) = register_process("21 1 -> 11 -12 2",11)
          ! g b > W t has a different signature
#endif
       elseif (ch.eq.'qg' .and. sec .eq.'rvewk_is') then
          call register_red_ew(1)
          call set_orders_g_virtew()
#if (_Vcharge == 0)
          OL_id(4) = register_process("1 21  -> 11 -11 1",11) !-- ask federico about crossing
          OL_id(5) = register_process("2 21  -> 11 -11 2",11)
          OL_id(6) = register_process("5 21  -> 11 -11 5",11)
          !
          OL_id(7) = register_process("-1 21 -> 11 -11 -1",11) !-- ask federico about crossing
          OL_id(8) = register_process("-2 21 -> 11 -11 -2",11)
          OL_id(9) = register_process("-5 21 -> 11 -11 -5",11)
#elseif (_Vcharge == 1)
          OL_id(3) = register_process("2 21 -> 12 -11 1",11)
          OL_id(4) = register_process("-1 21 -> 12 -11 -2",11)
#elseif (_Vcharge == -1)
          OL_id(3) = register_process("-2 21 -> 11 -12 -1",11)
          OL_id(4) = register_process("1 21 -> 11 -12 2",11)
#endif
       elseif (ch.eq.'gq' .and. sec .eq.'rvewknf_is') then
          call set_parameter("approximation","ewnf") !-- select 1L-ewk nf amplitudes from OpenLoops
          call register_red_ew(1)
          call set_orders_g_virtew()
#if (_Vcharge == 0)
          OL_id(4) = register_process("21 -1 -> 11 -11 -1",11) !-- ask federico about crossing
          OL_id(5) = register_process("21 -2 -> 11 -11 -2",11)
          !
          OL_id(7) = register_process("21 1 -> 11 -11 1",11) !-- ask federico about crossing
          OL_id(8) = register_process("21 2 -> 11 -11 2",11)
#elseif (_Vcharge == 1)
          OL_id(3) = register_process("21 2 -> 12 -11 1",11)
          OL_id(4) = register_process("21 -1 -> 12 -11 -2",11)
#elseif (_Vcharge == -1)
          OL_id(3) = register_process("21 -2 -> 11 -12 -1",11)
          OL_id(4) = register_process("21 1 -> 11 -12 2",11)
#endif    
       elseif (ch.eq.'qg' .and. sec .eq.'rvewknf_is') then
          call set_parameter("approximation","ewnf") !-- select 1L-ewk nf amplitudes from OpenLoops
          call register_red_ew(1)
          call set_orders_g_virtew()
#if (_Vcharge == 0)
          OL_id(4) = register_process("1 21  -> 11 -11 1",11) !-- ask federico about crossing
          OL_id(5) = register_process("2 21  -> 11 -11 2",11)
          !
          OL_id(7) = register_process("-1 21 -> 11 -11 -1",11) !-- ask federico about crossing
          OL_id(8) = register_process("-2 21 -> 11 -11 -2",11)
#elseif (_Vcharge == 1)
          OL_id(3) = register_process("2 21 -> 12 -11 1",11)
          OL_id(4) = register_process("-1 21 -> 12 -11 -2",11)
#elseif (_Vcharge == -1)
          OL_id(3) = register_process("-2 21 -> 11 -12 -1",11)
          OL_id(4) = register_process("1 21 -> 11 -12 2",11)
#endif
       elseif (ch.eq.'aq' .and. sec(1:5).eq.'rvqcd') then
          call register_red_qcd(1)
          call set_orders_a_virtqcd()
#if (_Vcharge == 0)
          OL_id(3) = register_process("22 -1 -> 11 -11 -1",11) !-- ask federico about crossing
          OL_id(4) = register_process("22 -2 -> 11 -11 -2",11)
          !
          OL_id(5) = register_process("22 1 -> 11 -11 1",11) !-- ask federico about crossing
          OL_id(6) = register_process("22 2 -> 11 -11 2",11)
#elseif (_Vcharge == 1)
          OL_id(3) = register_process("22 2 -> 12 -11 1",11)
          OL_id(4) = register_process("22 -1 -> 12 -11 -2",11)
#elseif (_Vcharge == -1)
          OL_id(3) = register_process("22 -2 -> 11 -12 -1",11)
          OL_id(4) = register_process("22 1 -> 11 -12 2",11)
#endif
       elseif (ch.eq.'qa' .and. sec(1:5).eq.'rvqcd') then
          call register_red_qcd(1)
          call set_orders_a_virtqcd()
#if (_Vcharge == 0)
          OL_id(3) = register_process("1 22  -> 11 -11 1",11) !-- ask federico about crossing
          OL_id(4) = register_process("2 22  -> 11 -11 2",11)
          !          
          OL_id(5) = register_process("-1 22 -> 11 -11 -1",11) !-- ask federico about crossing
          OL_id(6) = register_process("-2 22 -> 11 -11 -2",11)
#elseif	(_Vcharge == 1)
          OL_id(3) = register_process("2 22 -> 12 -11 1",11)
          OL_id(4) = register_process("-1 22 -> 12 -11 -2",11)
#elseif (_Vcharge == -1)
          OL_id(3) = register_process("-2 22 -> 11 -12 -1",11)
          OL_id(4) = register_process("1 22 -> 11 -12 2",11)
#endif
       endif

    elseif (corr.eq.'nnlo' .and. sec(1:2).eq.'rr') then
      call set_orders_rr()
      istart_rr = 1

      if(ch.eq.'ns_qqb_w') then

#if (_Vcharge == 0)
        OL_rr_id(istart_rr)   = register_process(" 2 -2 -> 11 -11 1 -1",1) !-- d  db -> e- e+ u ub
        OL_rr_id(istart_rr+1) = register_process("-2  2 -> 11 -11 1 -1",1) !-- db d  -> e- e+ u ub
        OL_rr_id(istart_rr+2) = register_process(" 1 -1 -> 11 -11 2 -2",1) !-- u  ub -> e- e+ d db
        OL_rr_id(istart_rr+3) = register_process("-1  1 -> 11 -11 2 -2",1) !-- ub u  -> e- e+ d db
#else
        print*, 'TO IMPLEMENT OL rr for W in the 4q case'
        stop
#endif

      elseif(ch.eq.'ns_qqp_w') then

#if (_Vcharge == 0)
        OL_rr_id(istart_rr+4) = register_process(" 2  1 -> 11 -11  2  1",1) !-- d  u  -> e- e+ d  u
        OL_rr_id(istart_rr+5) = register_process("-2 -1 -> 11 -11 -2 -1",1) !-- db ub -> e- e+ db ub
        OL_rr_id(istart_rr+6) = register_process(" 1  2 -> 11 -11  2  1",1) !-- u  d  -> e- e+ d  u
        OL_rr_id(istart_rr+7) = register_process("-1 -2 -> 11 -11 -2 -1",1) !-- ub db -> e- e+ db ub
#else
        print*, 'TO IMPLEMENT OL rr for W in the 4q case'
        stop
#endif
        
      elseif(ch.eq.'ns_qqpb_w') then

#if (_Vcharge == 0)
        OL_rr_id(istart_rr+8)  = register_process(" 2 -1 -> 11 -11  2 -1",1)  !-- d  ub -> e- e+ d  ub
        OL_rr_id(istart_rr+9)  = register_process("-2  1 -> 11 -11 -2  1",1)  !-- db u  -> e- e+ db u
        OL_rr_id(istart_rr+10) = register_process("-1  2 -> 11 -11  2 -1",1)  !-- ub d  -> e- e+ d  ub
        OL_rr_id(istart_rr+11) = register_process(" 1 -2 -> 11 -11 -2  1",1)  !-- u  db -> e- e+ db u
#else
        print*, 'TO IMPLEMENT OL rr for W in the 4q case'
        stop
#endif
        
      endif

   elseif (corr .eq. 'nnlo' .and. sec(1:4) .eq. 'subv') then
      if (sec .eq. 'subv_qcd') then
         call register_red_qcd(1)
      elseif (sec .eq. 'subv_ewk') then
         call register_red_ew(1)
      endif

   elseif (corr .eq. 'nnlo' .and. sec .eq. 's_vewk') then
      if (ch.eq.'ns' .or. ch.eq.'gq' .or. ch.eq.'qg') call register_red_ew(1)

    endif

    !-- set OL scale once and for all
    call set_parameter("mu",real(mu_ol,kind=dp15))
    
    call start()

    !-- get back the couplings
    call get_parameter("alpha_s",as_ol)
    call get_parameter("alpha",aem_ol)

    gsq_ol  = as_ol  * four * pi
    eesq_ol = aem_ol * four * pi
    
  contains

    !-- q qb -> e- e+ [NLO QCD]
    subroutine register_red_qcd(istart)
      integer, intent(in) :: istart
      call set_parameter("order_ew", 2)  !-- tree-level
      call set_parameter("order_qcd", 0) !-- tree-level
      call set_parameter("loop_order_ew", 2) !-- one-loop
      call set_parameter("loop_order_qcd",1) !-- one-loop
      !
#if (_Vcharge == 0)
      OL_id(istart  ) = register_process("1 -1 -> 11 -11",11) !-- d db -> e- e+
      OL_id(istart+1) = register_process("2 -2 -> 11 -11",11) !-- u ub -> e- e+ 
#endif

#if (_Vcharge == 1)
      OL_id(istart  ) = register_process("-1 2 -> 12 -11",11) !--
      OL_id(istart+1) = register_process("2 -1 -> 12 -11",11) !-- crossed channel
#endif

#if (_Vcharge == -1)
      OL_id(istart  ) = register_process("1 -2 -> 11 -12",11) !--
      OL_id(istart+1) = register_process("-2 1 -> 11 -12",11) !-- crossed channel
#endif
      
    end subroutine register_red_qcd

    !-- quark-quark > lep lep photon [NLO EW]
    subroutine register_red_ew(istart)
      integer, intent(in) :: istart
      call set_parameter("order_ew", 2)  !-- tree-level
      call set_parameter("order_qcd", 0) !-- tree-level
      call set_parameter("loop_order_ew", 3) !-- one-loop
      call set_parameter("loop_order_qcd",0) !-- one-loop

#if (_Vcharge == 0)
      !-- q qb -> e- e+ 
      OL_id(istart  ) = register_process("1 -1  -> 11 -11",11) !-- d db -> e- e+
      OL_id(istart+1) = register_process("2 -2  -> 11 -11",11) !-- u ub -> e- e+
      OL_id(istart+2) = register_process("5 -5  -> 11 -11",11) !-- b bb -> e- e+
#elseif (_Vcharge == 1)
      !-- q q' > W+ > v l+
      OL_id(istart) = register_process("-1 2  -> 12 -11",11) !-- db u -> ve e+
      OL_id(istart+1) = register_process("2 -1  -> 12 -11",11) !-- u db -> ve e+
#elseif (_Vcharge == -1)
      !-- q q' > W- > l- v~
      OL_id(istart) = register_process("1 -2  -> 11 -12",11) !-- d ub -> e- ve~ 
      OL_id(istart+1) = register_process("-2 1  -> 11 -12",11) !-- ub d -> e- ve~
#endif
      
    end subroutine register_red_ew

    !-- q qpb -> l nu [NLO EW]
    subroutine register_red_ewreal(iw)
      integer, intent(in) :: iw
      call set_parameter("order_ew", 3)  !-- tree-level
      !call set_parameter("order_qcd", 0) !-- tree-level --CB: turned off in order to avoid problem in the process registration

#if (_Vcharge == 0)
         OL_id(1  ) = register_process("1 -1  -> 11 -11 22",11) !-- d db -> e- e+ a
         OL_id(2  ) = register_process("2 -2  -> 11 -11 22",11) !-- u ub -> e- e+ a
#endif

#if (_Vcharge == 1)
         OL_id(1  ) = register_process("2 -1  -> -11 12 22",11) !-- u db -> e+ nu a
         OL_id(2  ) = register_process("-1 2  -> -11 12 22",11) !-- db u -> e+ nu a
#endif

#if (_Vcharge == -1)
         OL_id(1  ) = register_process("1 -2  -> 11 -12 22",11) !-- d ub -> e+ nu a
         OL_id(2  ) = register_process("-2 1  -> 11 -12 22",11) !-- ub d -> e+ nu a
#endif
      
    end subroutine register_red_ewreal

    
    !-- a a -> e- e+ [NLO EW]
    subroutine register_red_ew_aa(istart)
      integer, intent(in) :: istart
      call set_parameter("order_ew", 2)  !-- tree-level
      call set_parameter("order_qcd", 0) !-- tree-level
      call set_parameter("loop_order_ew", 3) !-- one-loop
      call set_parameter("loop_order_qcd",0) !-- one-loop

#if (_Vcharge == 0)
      if (ew_scheme.eq.'a0') then
         OL_id(istart) = register_process("2002 2002 -> 11 -11",11)   !-- a a   -> e- e+
      else
         OL_id(istart) = register_process("-2002 -2002 -> 11 -11",11) !-- a* a* -> e- e+
      endif
#endif

    end subroutine register_red_ew_aa

    !-- extra gluon emission, ew loop
    subroutine set_orders_g_virtew()
      call set_parameter("order_ew", 2)  !-- tree-level
      call set_parameter("order_qcd", 1) !-- tree-level
      call set_parameter("loop_order_ew", 3)
      call set_parameter("loop_order_qcd",1)
    end subroutine set_orders_g_virtew

    !-- extra photon emisssion, qcd loop
    subroutine set_orders_a_virtqcd()
      call set_parameter("order_ew", 3)  !-- tree-level
      call set_parameter("order_qcd", 0) !-- tree-level
      call set_parameter("loop_order_ew", 3)
      call set_parameter("loop_order_qcd",1)
    end subroutine set_orders_a_virtqcd

    !-- double real emission with four external quarks.
    subroutine set_orders_rr()
      call set_parameter("order_ew", 3) !-- tree-level
      call set_parameter("order_qcd",1) !-- tree-level
    end subroutine set_orders_rr

  end subroutine initialise_ol
  
end module mod_ol_interface

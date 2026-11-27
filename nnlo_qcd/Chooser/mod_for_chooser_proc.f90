module mod_for_chooser_proc
  use mod_types; use common_def; use mod_consts_dp
  use mod_consts_debug
  private

  public :: set_buff_default, set_default_proc, read_parms_proc

contains

  subroutine set_default_proc()
    
    !-- vN tags
    vn_taga  = one
    vn_tagb  = one
    vn_tagc  = one
    vn_tagab = one
    vn_tagac = one
    vn_tagbc = one
    vn_tagcd = one
    vn_tagce = one
    vn_tagcf = one

    tag_ns_gg = one
    tag_ns_nf = one
    tag_ns_id = one
    tag_ns_s = one

    tag_s_s = one
    tag_s_v = one
    
    ! -- tag for NLO channels
    do_nlo_qqb = 1
    do_nlo_qg  = 1

    ! -- tag for NNLO channels
    do_nnlo_qqb_ns = 1
    do_nnlo_qqb_ot = 1

    ! -- abelian/non abelian
    tagab  = one
    tagnab = one

    ! -- triple collinear
    tag_tcpls = one
    tag_tcdel = one
    tag_tcreg = one
    tag_nontc = one

    tag_12    = one
    tag_z12   = one
    tag_z1zb2 = one
  
  end subroutine set_default_proc

  !-- 

  subroutine read_parms_proc(numargs,countarg)
    integer, intent(in) :: numargs
    integer :: countarg, i
    character :: arg*(300)

    do i = 1, numargs

       call getarg(i,arg)

       !-- tags
       if(arg(1:10).eq.'tag_ns_gg=') then
          read(arg(11:14),*) tag_ns_gg
          countarg = countarg + 1
       elseif(arg(1:10).eq.'tag_ns_nf=') then
          read(arg(11:14),*) tag_ns_nf
          countarg = countarg + 1
       elseif(arg(1:10).eq.'tag_ns_id=') then
          read(arg(11:14),*) tag_ns_id
          countarg = countarg + 1
       elseif(arg(1:9).eq.'tag_ns_s=') then
          read(arg(10:14),*) tag_ns_s
          countarg = countarg + 1
       elseif(arg(1:8).eq.'tag_s_s=') then
          read(arg(9:12),*) tag_s_s
          countarg = countarg + 1
       elseif(arg(1:8).eq.'tag_s_v=') then
          read(arg(9:12),*) tag_s_v
          countarg = countarg + 1
       elseif(arg(1:8).eq.'vn_taga=') then
          read(arg(9:12),*) vn_taga
          countarg = countarg + 1
       elseif(arg(1:8).eq.'vn_tagb=') then
          read(arg(9:12),*) vn_tagb
          countarg = countarg + 1
       elseif(arg(1:8).eq.'vn_tagc=') then
          read(arg(9:12),*) vn_tagc
          countarg = countarg + 1
       elseif(arg(1:9).eq.'vn_tagab=') then
          read(arg(10:13),*) vn_tagab
          countarg = countarg + 1
       elseif(arg(1:9).eq.'vn_tagac=') then
          read(arg(10:13),*) vn_tagac
          countarg = countarg + 1
       elseif(arg(1:9).eq.'vn_tagbc=') then
          read(arg(10:13),*) vn_tagbc
          countarg = countarg + 1
       elseif(arg(1:9).eq.'vn_tagcd=') then
          read(arg(10:13),*) vn_tagcd
          countarg = countarg + 1
       elseif(arg(1:9).eq.'vn_tagce=') then
          read(arg(10:13),*) vn_tagce
          countarg = countarg + 1
       elseif(arg(1:9).eq.'vn_tagcf=') then
          read(arg(10:13),*) vn_tagcf
          countarg = countarg + 1

       !-- NLO channels tag
       elseif(arg(1:11).eq.'do_nlo_qqb=') then
          read(arg(12:15),*) do_nlo_qqb
          countarg = countarg + 1
       elseif(arg(1:10).eq.'do_nlo_qg=') then
          read(arg(11:14),*) do_nlo_qg
          countarg = countarg + 1

       !-- NNLO channels tag
       elseif(arg(1:15).eq.'do_nnlo_qqb_ns=') then
          read(arg(16:17),*) do_nnlo_qqb_ns
          countarg = countarg + 1
       elseif(arg(1:15).eq.'do_nnlo_qqb_ot=') then
          read(arg(16:17),*) do_nnlo_qqb_ot
          countarg = countarg + 1

       ! -- abelian/non abelian
       elseif(arg(1:6).eq.'tagab=') then
          read(arg(7:8),*) tagab
          countarg = countarg + 1
       elseif(arg(1:7).eq.'tagnab=') then
          read(arg(8:9),*) tagnab
          countarg = countarg + 1
          
       !-- TC tags
       elseif(arg(1:10).eq.'tag_tcpls=') then
          read(arg(11:14),*) tag_tcpls
          countarg = countarg + 1
       elseif(arg(1:10).eq.'tag_tcdel=') then
          read(arg(11:14),*) tag_tcdel
          countarg = countarg + 1
       elseif(arg(1:10).eq.'tag_tcreg=') then
          read(arg(11:14),*) tag_tcreg
          countarg = countarg + 1
       elseif(arg(1:10).eq.'tag_nontc=') then
          read(arg(11:14),*) tag_nontc
          countarg = countarg + 1
       elseif(arg(1:7).eq.'tag_12=') then
          read(arg(8:11),*) tag_12
          countarg = countarg + 1
       elseif(arg(1:8).eq.'tag_z12=') then
          read(arg(9:12),*) tag_z12
          countarg = countarg + 1
       elseif(arg(1:10).eq.'tag_z1zb2=') then
          read(arg(11:14),*) tag_z1zb2
          countarg = countarg + 1
       endif

    enddo

       ! now recompute the EW parameters for the input mw
       call compute_ew_params()


    return

  end subroutine read_parms_proc

 subroutine compute_ew_params()
    implicit none

    mwsq = mw**2
    cw = mw/mz
    cosW2 = mwsq/mzsq
    sinW2 = one - cosW2
    alpha_ew = (sqrt2*Gf*mwsq*(one - mwsq/mzsq))/pi
    sw = sqrt(sinW2)
    eesq  = 4*pi*alpha_ew
    ee    = sqrt(eesq)
    gwsq  = eesq/sinW2


    ! couplings --  quarks
    Vup  = 1.0_dp/2.0_dp - 4.0_dp/3.0_dp*sinW2
    Aup  = 1.0_dp/2.0_dp
    Vdn  = -1.0_dp/2.0_dp + 2.0_dp/3.0_dp*sinW2
    Adn  = -1.0_dp/2.0_dp
    Lup  = Vup + Aup
    Rup  = Vup - Aup
    Ldn  = Vdn + Adn
    Rdn  = Vdn - Adn
    cLup = Lup/(two*cw*sw) !-- coupling in units of e
                                                                                                                                                                                                             
    cRup = Rup/(two*cw*sw)
    cLdn = Ldn/(two*cw*sw)
    cRdn = Rdn/(two*cw*sw)

!     couplings -- leptons

    Vel  = -half + 2.0_dp*sinW2
    Ael  = -half
    Lel  = Vel + Ael
    Rel  = Vel - Ael
    Lnu  = Vnu + Anu
    Rnu  = Vnu - Anu
    cLel = Lel/(two*cw*sw)
    cRel = Rel/(two*cw*sw)
    cLnu = Lnu/(two*cw*sw)
    cRnu = Rnu/(two*cw*sw)

    cAWW  = one
    cZWW  = -cw/sw
    cLWud  = 1/sqrt2/sw*Vud
    cLWnue = 1/sqrt2/sw*Vnue

    cL_lep = cLel
    cR_lep = cRel


  end subroutine compute_ew_params


  subroutine set_buff_default(sec,ch)
    !-- these are default settings with buffrr_pow = 1
    character, intent(in) :: sec*(10),ch*(10)

    if (sec .ne. 'sub12' .and. sec .ne. 'sub12' .and. sec .ne. 'subv') then  ! these sectors don't need buffrr
       if (buffrr_pow .ne. 1) then
          print *, 'no default buffer for buffrr_pow != 1 implemented'
          stop
       endif
    endif

    if (ch .eq. 'ns' .or. ch .eq. 's' .or. ch .eq. 'qg' .or. ch .eq. 'gq' ) then

       if ( sec .eq. 'rr4151a' .or. sec .eq. 'rr4252a' ) buffrr = 1E-10_dp
       if ( sec .eq. 'rr4151c' .or. sec .eq. 'rr4252c' ) buffrr = 1E-10_dp

       if ( sec .eq. 'rr4151b' .or. sec .eq. 'rr4252b' ) buffrr = 1E-6_dp
       if ( sec .eq. 'rr4151d' .or. sec .eq. 'rr4252d' ) buffrr = 1E-7_dp

       if ( sec .eq. 'rr4152' .or. sec .eq. 'rr4251'   ) buffrr = 1E-11_dp

    endif

    !-- gg
    if (sec.eq.'ggrr') buffrr = 1E-9_dp

    !-- 2q
    if (sec.eq.'2q1415a' .or. &
        sec.eq.'2q1415b' .or. &
        sec.eq.'2q2425a' .or. &
        sec.eq.'2q2425b' ) then
       
       buffrr = 1E-7_dp

    endif

    if (buffrr.eq.-one) then
       buffrr = 1E-8_dp
       print *, 'buffer not set for sector ', sec
       print *, 'buffer set to default ', buffrr
    endif

  end subroutine set_buff_default

end module mod_for_chooser_proc

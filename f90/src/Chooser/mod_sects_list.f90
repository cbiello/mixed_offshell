module mod_sects_list
  use mod_types
  implicit none
  private

  public :: print_out_channels_list, print_out_corrections_list, print_out_sectors_list

  !-- number of partonic channels
  integer, parameter :: n_channels = 14

  !-- number of higher-order corrections
  integer, parameter :: n_corrections = 4

  character(len=10), allocatable :: ch_list(:)
  character(len=10), allocatable :: corr_list(:)

  character(len=20), allocatable :: ns_lo(:),ns_nloqcd(:),ns_nloewk(:),ns_nnlo(:)
  character(len=20), allocatable :: aa_lo(:),aa_nloewk(:)
  character(len=20), allocatable :: qg_nloqcd(:),qg_nnlo(:),gq_nloqcd(:),gq_nnlo(:)
  character(len=20), allocatable :: qa_nloewk(:),qa_nnlo(:),aq_nloewk(:),aq_nnlo(:)
  character(len=20), allocatable :: ag_nnlo(:),ga_nnlo(:)

  character(len=20), allocatable :: ns_ga_nnlo(:),ns_qqb_nnlo(:),ns_qq_nnlo(:)
  character(len=20), allocatable :: ns_qqb_w_nnlo(:),ns_qqp_w_nnlo(:),ns_qqpb_w_nnlo(:)

contains

  subroutine print_out_channels_list(idev)
    integer, intent(in) :: idev
    integer :: i

    call initialise_ch_corr_sec_list(1)
    do i = 1, size(ch_list)
      write(idev,*) '-ch ', ch_list(i)
    enddo
    call finalise_ch_corr_sec_list()

  end subroutine print_out_channels_list

  subroutine print_out_corrections_list(idev)
    integer, intent(in) :: idev
    integer :: i

    call initialise_ch_corr_sec_list(2)
    do i = 1, size(corr_list)
      write(idev,*) '-corr ', corr_list(i)
    enddo
    call finalise_ch_corr_sec_list()

  end subroutine print_out_corrections_list

  subroutine print_out_sectors_list()
    integer :: i, outid, sel_lst, sel_ch, sel_corr
    character(len=32) :: iString

    write(6,*) 'List sectors by [1/2]: '
    write(6,*) '1: corrections'
    write(6,*) '2: channels   '
    read *, sel_lst

    write(6,*) ''

    outid = 7
    
    if(sel_lst==1) then
      write(6,*) "Select perturbative order. Possible choices: "
      call initialise_ch_corr_sec_list(2)
      write(6,*) '0  : all'
      do i = 1, size(corr_list)
        if(i < 10) then
          write(iString,*) i
          write(6,*) trim(adjustl(iString)), '  : ' , adjustl(corr_list(i))
        else
          write(iString,*) i
          write(6,*) trim(adjustl(iString)), ' : ' , adjustl(corr_list(i))
        endif
      enddo

      read *, sel_corr

      call initialise_ch_corr_sec_list(3)

      if(sel_corr .eq. 1) then
        open(outid, file = "lo_sectors_list.dat", status = 'replace')
        call write_sects_lo(outid)
        write(6,*) 'List of sectors written to file: lo_sectors_list.dat'
      elseif(sel_corr.eq. 2) then
        open(outid, file = "nloqcd_sectors_list.dat", status = 'replace')
        call write_sects_nloqcd(outid)
        write(6,*) 'List of sectors written to file: nloqcd_sectors_list.dat'
      elseif(sel_corr.eq. 3) then
        open(outid, file = "nloewk_sectors_list.dat", status = 'replace')
        call write_sects_nloewk(outid)
        write(6,*) 'List of sectors written to file: nloewk_sectors_list.dat'
      elseif(sel_corr.eq. 4) then
        open(outid, file = "nnlo_sectors_list.dat", status = 'replace')
        call write_sects_nnlo(outid)
        write(6,*) 'List of sectors written to file: nnlo_sectors_list.dat'
      elseif(sel_corr.eq. 0) then
        open(outid, file = "sectors_list.dat", status = 'replace')
        call write_sects_lo(outid)
        call write_sects_nloqcd(outid)
        call write_sects_nloewk(outid)
        call write_sects_nnlo(outid)
        write(6,*) 'List of sectors written to file: sectors_list.dat'
      endif

    else if(sel_lst==2) then
      write(6,*) "Select partonic channel. Possible choices: "
      call initialise_ch_corr_sec_list(1)
      write(6,*) '0  : all'
      do i = 1, size(ch_list)
        if(i < 10) then
          write(iString,*) i
          write(6,*) trim(adjustl(iString)), '  : ' , adjustl(ch_list(i))
        else
          write(iString,*) i
          write(6,*) trim(adjustl(iString)), ' : ' , adjustl(ch_list(i))
        endif
      enddo

      read *, sel_ch

      outid = 7
      open(outid, file = "sectors_list.dat", status = 'replace')
      call initialise_ch_corr_sec_list(3)

      if(sel_ch .eq. 1) then
        call write_sects_ns(outid)
      elseif(sel_ch .eq. 2) then
        call write_sects_aa(outid)
      elseif(sel_ch .eq. 3) then
        call write_sects_qg(outid)
      elseif(sel_ch .eq. 4) then
        call write_sects_gq(outid)
      elseif(sel_ch .eq. 5) then
        call write_sects_qa(outid)
      elseif(sel_ch .eq. 6) then
        call write_sects_aq(outid)
      elseif(sel_ch .eq. 7) then
        call write_sects_ag(outid)
      elseif(sel_ch .eq. 8) then
        call write_sects_ga(outid)
      elseif(sel_ch .eq. 9) then
        call write_sects_ns_ga(outid)
      elseif(sel_ch .eq. 10) then
        call write_sects_ns_qqb(outid)
      elseif(sel_ch .eq. 11) then
        call write_sects_ns_qq(outid)
      elseif(sel_ch .eq. 12) then
        call write_sects_ns_qqb_w(outid)
      elseif(sel_ch .eq. 13) then
        call write_sects_ns_qqp_w(outid)
      elseif(sel_ch .eq. 14) then
        call write_sects_ns_qqpb_w(outid)
      elseif(sel_ch .eq. 0) then
        call write_sects_ns(outid)        ! 1
        call write_sects_aa(outid)        ! 2
        call write_sects_qg(outid)        ! 3
        call write_sects_gq(outid)        ! 4
        call write_sects_qa(outid)        ! 5
        call write_sects_aq(outid)        ! 6
        call write_sects_ag(outid)        ! 7
        call write_sects_ga(outid)        ! 8
        call write_sects_ns_ga(outid)     ! 9
        call write_sects_ns_qqb(outid)    ! 10
        call write_sects_ns_qq(outid)     ! 11
        call write_sects_ns_qqb_w(outid)  ! 12
        call write_sects_ns_qqp_w(outid)  ! 13
        call write_sects_ns_qqpb_w(outid) ! 14
      endif
      write(6,*) 'List of sectors written to file: sectors_list.dat'
    else
      write(6,*) 'Invalid selection. Aborting'
    endif

    call finalise_ch_corr_sec_list()

    close(outid)

  end subroutine print_out_sectors_list

  !---------------------------------------------------------------------------!

  subroutine initialise_ch_corr_sec_list(list_id)
    integer, intent(in) :: list_id

    if(list_id == 1) then
      allocate (ch_list(n_channels))
      ch_list = [character(len=10) :: &
        "ns",       "aa",       "qg",       "gq",       "qa",       "aq",       & !-- 6
        "ag",       "ga",       "ns_ga",    "ns_qqb",   "ns_qq",    "ns_qqb_w", & !-- 12
        "ns_qqp_w", "ns_qqpb_w"]

    elseif(list_id == 2) then
      allocate (corr_list(n_corrections))
      corr_list = [character(len=10) :: "lo", "nloqcd", "nloewk", "nnlo"]

    else if(list_id == 3) then
      allocate (ns_lo(1))
      allocate (ns_nloqcd(3))
      allocate (ns_nloewk(5))
      allocate (ns_nnlo(16))
      ns_lo     = [character(len=20) :: 'na']
      ns_nloqcd = [character(len=20) :: 'r_is','s','v']
      ns_nloewk = [character(len=20) :: 'r_is','s','v','r_fs_53','r_fs_54']
      ns_nnlo   = [character(len=20) :: 'vvnf','vvfc','vvxf','rvewk_is','rvqcd_is','rvqcd_fs_53','rvqcd_fs_54',&
        's_vewk','s_vqcd','s_oqcd','s_oewk_is','s_oewk_fs_53','s_oewk_fs_54','s_12','s_qqb','s_qq']

      allocate (aa_lo(1))
      allocate (aa_nloewk(4))
      aa_lo     = [character(len=20) :: 'na']
      aa_nloewk = [character(len=20) :: 's','v','r_fs_53','r_fs_54']

      allocate (qg_nloqcd(2))
      allocate (qg_nnlo(14))
      qg_nloqcd  = [character(len=20) :: 'r_is','s']
      qg_nnlo    = [character(len=20) :: 'rvewk_is', &
        'rr_5262a','rr_5262b','rr_5262c','rr_5262d','rr_5261','rr_5263','rr_5264', &
        's_vewk','s_oqcd','s_oewk_is','s_oewk_fs_53','s_oewk_fs_54','s_12']

      allocate (gq_nloqcd(2))
      allocate (gq_nnlo(14))
      gq_nloqcd  = [character(len=20) :: 'r_is','s']
      gq_nnlo    = [character(len=20) :: 'rvewk_is', &
        'rr_5161a','rr_5161b','rr_5161c','rr_5161d','rr_5162','rr_5163','rr_5164', &
        's_vewk','s_oqcd','s_oewk_is','s_oewk_fs_53','s_oewk_fs_54','s_12']

      allocate (qa_nloewk(2))
      allocate (qa_nnlo(15))
      qa_nloewk  = [character(len=20) :: 'r_is','s']
      qa_nnlo    = [character(len=20) :: 'rvqcd_is', &
        'rr_5161a','rr_5161b','rr_5161c','rr_5161d','rr_5262a','rr_5262b','rr_5262c','rr_5262d','rr_5162','rr_5261', &
        's_vqcd','s_oqcd','s_oewk','s_12']

      allocate (aq_nloewk(2))
      allocate (aq_nnlo(15))
      aq_nloewk  = [character(len=20) :: 'r_is','s']
      aq_nnlo    = [character(len=20) :: 'rvqcd_is', &
        'rr_5161a','rr_5161b','rr_5161c','rr_5161d','rr_5262a','rr_5262b','rr_5262c','rr_5262d','rr_5162','rr_5261', &
        's_vqcd','s_oqcd','s_oewk','s_12']

      allocate (ag_nnlo(11))
      ag_nnlo = [character(len=20) :: & 
        'rr_5161','rr_5262a','rr_5262b','rr_5262c','rr_5262d','rr_5162','rr_5261', &
        's_oqcd','s_oewk','s_12']

      allocate (ga_nnlo(11))
      ga_nnlo = [character(len=20) :: &
        'rr_5161a','rr_5161b','rr_5161c','rr_5161d','rr_5262','rr_5162','rr_5261', &
        's_oqcd','s_oewk','s_12']

      allocate (ns_ga_nnlo(10))
      ns_ga_nnlo = [character(len=20) :: &
        'rr_5161a','rr_5161c','rr_5262a','rr_5262c','rr_5162', 'rr_5163','rr_5164','rr_5261','rr_5263','rr_5264']

      allocate (ns_qqb_nnlo(4))
      ns_qqb_nnlo = [character(len=20) :: 'rr_5161a','rr_5161c','rr_5262a','rr_5262c']

      allocate (ns_qq_nnlo(4))
      ns_qq_nnlo = [character(len=20) :: 'rr_5161a','rr_5161c','rr_5262a','rr_5262c']

      allocate (ns_qqb_w_nnlo(1))
      ns_qqb_w_nnlo  = [character(len=20) :: 'rr']
      allocate (ns_qqp_w_nnlo(1))
      ns_qqp_w_nnlo  = [character(len=20) :: 'rr']
      allocate (ns_qqpb_w_nnlo(1))
      ns_qqpb_w_nnlo = [character(len=20) :: 'rr']

    endif

  end subroutine initialise_ch_corr_sec_list

  subroutine finalise_ch_corr_sec_list()
    if(allocated(ch_list))    deallocate(ch_list)
    if(allocated(corr_list))  deallocate(corr_list)

    if(allocated(ns_lo))      deallocate(ns_lo)
    if(allocated(ns_nloqcd))  deallocate(ns_nloqcd)
    if(allocated(ns_nloewk))  deallocate(ns_nloewk)
    if(allocated(ns_nnlo))    deallocate(ns_nnlo)

    if(allocated(aa_lo))      deallocate(aa_lo)
    if(allocated(aa_nloewk))  deallocate(aa_nloewk)

    if(allocated(qg_nloqcd))  deallocate(qg_nloqcd)
    if(allocated(qg_nnlo))    deallocate(qg_nnlo)
    if(allocated(gq_nloqcd))  deallocate(gq_nloqcd)
    if(allocated(gq_nnlo))    deallocate(gq_nnlo)

    if(allocated(qa_nloewk))  deallocate(qa_nloewk)
    if(allocated(qa_nnlo))    deallocate(qa_nnlo)
    if(allocated(aq_nloewk))  deallocate(aq_nloewk)
    if(allocated(aq_nnlo))    deallocate(aq_nnlo)

    if(allocated(ag_nnlo))    deallocate(ag_nnlo)
    if(allocated(ga_nnlo))    deallocate(ga_nnlo)

    if(allocated(ns_ga_nnlo))  deallocate(ns_ga_nnlo)
    if(allocated(ns_qqb_nnlo)) deallocate(ns_qqb_nnlo)
    if(allocated(ns_qq_nnlo))  deallocate(ns_qq_nnlo)

    if(allocated(ns_qqb_w_nnlo))  deallocate(ns_qqb_w_nnlo)
    if(allocated(ns_qqp_w_nnlo))  deallocate(ns_qqp_w_nnlo)
    if(allocated(ns_qqpb_w_nnlo)) deallocate(ns_qqpb_w_nnlo)

  end subroutine finalise_ch_corr_sec_list

  !---------------------------------------------------------------------------!

  subroutine write_sects_lo(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- LO ns
    do i = 1, size(ns_lo)
      write(outid,*) '-corr lo -ch ns -sec ', ns_lo(i)
    enddo
    !-- LO aa
    do i = 1, size(aa_lo)
      write(outid,*) '-corr lo -ch aa -sec ', aa_lo(i)
    enddo
  end subroutine write_sects_lo

  subroutine write_sects_nloqcd(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NLO QCD ns
    do i = 1, size(ns_nloqcd)
      write(outid,*) '-corr nloqcd -ch ns -sec ', ns_nloqcd(i)
    enddo
    !-- NLO QCD qg
    do i = 1, size(qg_nloqcd)
      write(outid,*) '-corr nloqcd -ch qg -sec ', qg_nloqcd(i)
    enddo
    !-- NLO QCD gq
    do i = 1, size(gq_nloqcd)
      write(outid,*) '-corr nloqcd -ch gq -sec ', gq_nloqcd(i)
    enddo
  end subroutine write_sects_nloqcd

  subroutine write_sects_nloewk(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NLO EWK ns
    do i = 1, size(ns_nloewk)
      write(outid,*) '-corr nloewk -ch ns -sec ', ns_nloewk(i)
    enddo
    !-- NLO EWK aa
    do i = 1, size(aa_nloewk)
      write(outid,*) '-corr nloewk -ch aa -sec ', aa_nloewk(i)
    enddo
    !-- NLO EWK qa
    do i = 1, size(qa_nloewk)
      write(outid,*) '-corr nloewk -ch qa -sec ', qa_nloewk(i)
    enddo
    !-- NLO EWK aq
    do i = 1, size(aq_nloewk)
      write(outid,*) '-corr nloewk -ch aq -sec ', aq_nloewk(i)
    enddo
  end subroutine write_sects_nloewk

  subroutine write_sects_nnlo(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NNLO QCDxEW ns
    do i = 1, size(ns_nnlo)
      write(outid,*) '-corr nnlo -ch ns -sec ', ns_nnlo(i)
    enddo
    !-- NNLO QCDxEW qg
    do i = 1, size(qg_nnlo)
      write(outid,*) '-corr nnlo -ch qg -sec ', qg_nnlo(i)
    enddo
    !-- NNLO QCDxEW gq
    do i = 1, size(gq_nnlo)
      write(outid,*) '-corr nnlo -ch gq -sec ', gq_nnlo(i)
    enddo
    !-- NNLO QCDxEW qa
    do i = 1, size(qa_nnlo)
      write(outid,*) '-corr nnlo -ch qa -sec ', qa_nnlo(i)
    enddo
    !-- NNLO QCDxEW aq
    do i = 1, size(aq_nnlo)
      write(outid,*) '-corr nnlo -ch aq -sec ', aq_nnlo(i)
    enddo
    !-- NNLO QCDxEW ag
    do i = 1, size(ag_nnlo)
      write(outid,*) '-corr nnlo -ch ag -sec ', ag_nnlo(i)
    enddo
    !-- NNLO QCDxEW ga
    do i = 1, size(ga_nnlo)
      write(outid,*) '-corr nnlo -ch ga -sec ', ga_nnlo(i)
    enddo
    !-- NNLO QCDxEW
    do i = 1, size(ns_ga_nnlo)
      write(outid,*) '-corr nnlo -ch ns_ga -sec ', ns_ga_nnlo(i)
    enddo
    !-- NNLO QCDxEW
    do i = 1, size(ns_qqb_nnlo)
      write(outid,*) '-corr nnlo -ch ns_qqb -sec ', ns_qqb_nnlo(i)
    enddo
    do i = 1, size(ns_qq_nnlo)
      write(outid,*) '-corr nnlo -ch ns_qq -sec ', ns_qq_nnlo(i)
    enddo
    !-- NNLO QCDxEW
    do i = 1, size(ns_qqb_w_nnlo)
      write(outid,*) '-corr nnlo -ch ns_qqb_w -sec ', ns_qqb_w_nnlo(i)
    enddo
    !-- NNLO QCDxEW
    do i = 1, size(ns_qqp_w_nnlo)
      write(outid,*) '-corr nnlo -ch ns_qqp_w -sec ', ns_qqp_w_nnlo(i)
    enddo
    !-- NNLO QCDxEW
    do i = 1, size(ns_qqpb_w_nnlo)
      write(outid,*) '-corr nnlo -ch ns_qqpb_w -sec ', ns_qqpb_w_nnlo(i)
    enddo
  end subroutine write_sects_nnlo

  !---------------------------------------------------------------------------!

  !-- ns
  subroutine write_sects_ns(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- LO
    do i = 1, size(ns_lo)
      write(outid,*) '-corr lo -ch ns -sec ', ns_lo(i)
    enddo
    !-- NLO QCD
    do i = 1, size(ns_nloqcd)
      write(outid,*) '-corr nloqcd -ch ns -sec ', ns_nloqcd(i)
    enddo
    !-- NLO EWK
    do i = 1, size(ns_nloewk)
      write(outid,*) '-corr nloewk -ch ns -sec ', ns_nloewk(i)
    enddo
    !-- NNLO QCDxEW
    do i = 1, size(ns_nnlo)
      write(outid,*) '-corr nnlo -ch ns -sec ', ns_nnlo(i)
    enddo
  end subroutine write_sects_ns

  !-- aa
  subroutine write_sects_aa(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- LO
    do i = 1, size(aa_lo)
      write(outid,*) '-corr lo -ch aa -sec ', aa_lo(i)
    enddo
    !-- NLO EWK
    do i = 1, size(aa_nloewk)
      write(outid,*) '-corr nloewk -ch aa -sec ', aa_nloewk(i)
    enddo
  end subroutine write_sects_aa

  !-- qg
  subroutine write_sects_qg(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NLO QCD
    do i = 1, size(qg_nloqcd)
      write(outid,*) '-corr nloqcd -ch qg -sec ', qg_nloqcd(i)
    enddo
    !-- NNLO QCDxEW
    do i = 1, size(qg_nnlo)
      write(outid,*) '-corr nnlo -ch qg -sec ', qg_nnlo(i)
    enddo
  end subroutine write_sects_qg

  !-- gq
  subroutine write_sects_gq(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NLO QCD
    do i = 1, size(gq_nloqcd)
      write(outid,*) '-corr nloqcd -ch gq -sec ', gq_nloqcd(i)
    enddo
    !-- NNLO QCDxEW
    do i = 1, size(gq_nnlo)
      write(outid,*) '-corr nnlo -ch gq -sec ', gq_nnlo(i)
    enddo
  end subroutine write_sects_gq

  !-- qa
  subroutine write_sects_qa(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NLO EWK
    do i = 1, size(qa_nloewk)
      write(outid,*) '-corr nloewk -ch qa -sec ', qa_nloewk(i)
    enddo
    !-- NNLO QCDxEW
    do i = 1, size(qa_nnlo)
      write(outid,*) '-corr nnlo -ch qa -sec ', qa_nnlo(i)
    enddo
  end subroutine write_sects_qa

  !-- aq
  subroutine write_sects_aq(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NLO EWK
    do i = 1, size(aq_nloewk)
      write(outid,*) '-corr nloewk -ch aq -sec ', aq_nloewk(i)
    enddo
    !-- NNLO QCDxEW
    do i = 1, size(aq_nnlo)
      write(outid,*) '-corr nnlo -ch aq -sec ', aq_nnlo(i)
    enddo
  end subroutine write_sects_aq

  !-- ag
  subroutine write_sects_ag(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NNLO QCDxEW
    do i = 1, size(ag_nnlo)
      write(outid,*) '-corr nnlo -ch ag -sec ', ag_nnlo(i)
    enddo
  end subroutine write_sects_ag

  !-- ga
  subroutine write_sects_ga(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NNLO QCDxEW
    do i = 1, size(ga_nnlo)
      write(outid,*) '-corr nnlo -ch ga -sec ', ga_nnlo(i)
    enddo
  end subroutine write_sects_ga

  !-- qqx(qxq) -> l-l+ g a
  subroutine write_sects_ns_ga(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NNLO QCDxEW
    do i = 1, size(ns_ga_nnlo)
      write(outid,*) '-corr nnlo -ch ns_ga -sec ', ns_ga_nnlo(i)
    enddo
  end subroutine write_sects_ns_ga

  !-- qqx(qxq) -> l-l+ q qx
  subroutine write_sects_ns_qqb(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NNLO QCDxEW
    do i = 1, size(ns_qqb_nnlo)
      write(outid,*) '-corr nnlo -ch ns_qqb -sec ', ns_qqb_nnlo(i)
    enddo
  end subroutine write_sects_ns_qqb

  !-- qqx(qxq) -> l-l+ q qx
  subroutine write_sects_ns_qq(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NNLO QCDxEW
    do i = 1, size(ns_qq_nnlo)
      write(outid,*) '-corr nnlo -ch ns_qq -sec ', ns_qq_nnlo(i)
    enddo
  end subroutine write_sects_ns_qq

  !-- qqx(qxq) -> l-l+ q qx (W exchange)
  subroutine write_sects_ns_qqb_w(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NNLO QCDxEW
    do i = 1, size(ns_qqb_w_nnlo)
      write(outid,*) '-corr nnlo -ch ns_qqb_w -sec ', ns_qqb_w_nnlo(i)
    enddo
  end subroutine write_sects_ns_qqb_w

  !-- qqx(qxq) -> l-l+ q qx (W exchange)
  subroutine write_sects_ns_qqp_w(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NNLO QCDxEW
    do i = 1, size(ns_qqp_w_nnlo)
      write(outid,*) '-corr nnlo -ch ns_qqp_w -sec ', ns_qqp_w_nnlo(i)
    enddo
  end subroutine write_sects_ns_qqp_w

  !-- qqx(qxq) -> l-l+ q qx (W exchange)
  subroutine write_sects_ns_qqpb_w(outid)
    integer, intent(in) :: outid
    integer :: i
    !-- NNLO QCDxEW
    do i = 1, size(ns_qqpb_w_nnlo)
      write(outid,*) '-corr nnlo -ch ns_qqpb_w -sec ', ns_qqpb_w_nnlo(i)
    enddo
  end subroutine write_sects_ns_qqpb_w

end module mod_sects_list

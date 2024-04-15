module mod_hoppet_nnlo
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF, except => Cf, except => Ca, except => tr
  use mod_types, except => dp
  use consts_dp
  use mod_parms
  use mod_proc_parms
  use convolution
  use convolution_communicator
  use mod_hoppet_tools
  use mod_auxfunctions
  implicit none
  private
  !-- res = xPij + xPij_Lmu * log(musq/mvsq)
  !-- could put Lmu term inside xPij, but would not work for dynamic mvsq
  type(pdf_table), public, save :: xPij,xPij_Lmu,xPij_Lmu2
  type(pdf_table), public, save :: xPij_a,xPij_a_Lmu,xPij_a_Lmu2

  type(pdf_table), public, save :: xPij_A_1,xPij_A_1_Lmu,xPij_A_2,xPij_A_2_Lmu
  type(pdf_table), public, save :: xPij_B_1,xPij_B_1_Lmu,xPij_B_2
  type(pdf_table), public, save :: xPij_C_1,xPij_C_2,xPij_C_2_Lmu
  type(pdf_table), public, save :: xPij_D_1,xPij_D_2

  public :: init_xPij_nnlo

contains

  subroutine init_xPij_nnlo()
    integer :: iQ,iflav

    type(grid_conv) :: mySub_ns,mySub_ns_Lmu
    type(grid_conv) :: mySub_qq,mySub_qq_Lmu,mySub_qq_b,mySub_qq_b_Lmu
    type(grid_conv) :: mySub_TCqq,mySub_TCqq_Lmu,mySub_TCqq_Lmu2
    type(grid_conv) :: mySub_ns_qqb,mySub_ns_qqb_Lmu
    type(grid_conv) :: mySub_ns_qq,mySub_ns_qq_Lmu

    type(grid_conv) :: mySub_gq,mySub_gq_Lmu,mySub_gq_b,mySub_gq_b_Lmu
    type(grid_conv) :: mySub_TCgq,mySub_TCgq_Lmu,mySub_TCgq_Lmu2

    type(grid_conv) :: mySub_TCaa,mySub_TCaa_Lmu,mySub_TCaa_Lmu2

    type(grid_conv) :: mySub_TCag,mySub_TCag_Lmu,mySub_TCag_Lmu2

    !-- use:
    !-- call AllocGridConv(grid,myPij)
    !-- call InitGridConv(grid,myPij,function_to_load)
    !-- call AddWithCoeff(myPij_1,myPij_2)               --> myPij_1 = myPij_1 + myPij_2
    !-- call AddWithCoeff(myPij_1,myPij_2,const)         --> myPij_1 = myPij_1 + myPij_2*const
    !-- call SetToConvolution(myPij_1x2,myPij_1,myPij_2) --> myPij_1x2 = myPij_1 \otimes myPij_2

    !!=======================================================================!!
    !!                  Grids allocation and initialisation                  !!
    !!=======================================================================!!

    !-- PDFs Tables
    call init_pdfs(xPij)
    call init_pdfs(xPij_Lmu)
    call init_pdfs(xPij_Lmu2)
    call init_pdfs(xPij_a)
    call init_pdfs(xPij_a_Lmu)
    call init_pdfs(xPij_a_Lmu2)

    call init_pdfs(xPij_A_1)
    call init_pdfs(xPij_A_1_Lmu)
    call init_pdfs(xPij_A_2)
    call init_pdfs(xPij_A_2_Lmu)

    call init_pdfs(xPij_B_1)
    call init_pdfs(xPij_B_1_Lmu)
    call init_pdfs(xPij_B_2)

    call init_pdfs(xPij_C_1)
    call init_pdfs(xPij_C_2)
    call init_pdfs(xPij_C_2_Lmu)

    call init_pdfs(xPij_D_1)
    call init_pdfs(xPij_D_2)

    !-- used for single boosts
    xPij%tab = 0
    xPij_Lmu%tab = 0
    xPij_Lmu2%tab = 0
    xPij_a%tab = 0
    xPij_a_Lmu%tab = 0
    xPij_a_Lmu2%tab = 0
    !-- used for double boosts
    xPij_A_1%tab = 0
    xPij_A_1_Lmu%tab = 0
    xPij_A_2%tab = 0
    xPij_A_2_Lmu%tab = 0
    xPij_B_1%tab = 0
    xPij_B_1_Lmu%tab = 0
    xPij_B_2%tab = 0
    xPij_C_1%tab = 0
    xPij_C_2%tab = 0
    xPij_C_2_Lmu%tab = 0
    xPij_D_1%tab = 0
    xPij_D_2%tab = 0


    !!-----------------------------------------------------------------------!!

    iQ = 0

    !!=======================================================================!!
    !!                              qqb channel                              !!
    !!=======================================================================!!

    if(ch.eq.'ns') then

      if(sec.eq.'s_12') then

        ! ---------------- Triple Collinear, TC_gq \otimes q ---------------- !
        call AllocGridConv(grid,mySub_TCqq)
        call AllocGridConv(grid,mySub_TCqq_Lmu)
        call AllocGridConv(grid,mySub_TCqq_Lmu2)
        call InitGridConv(grid,mySub_TCqq,mySub_TCqq_func)
        call InitGridConv(grid,mySub_TCqq_Lmu,mySub_TCqq_Lmu_func)
        call InitGridConv(grid,mySub_TCqq_Lmu2,mySub_TCqq_Lmu2_func)

        ! -------------------------- Pqq \otimes q -------------------------- !
        call AllocGridConv(grid,mySub_qq)
        call AllocGridConv(grid,mySub_qq_Lmu)
        call InitGridConv(grid,mySub_qq,PqqNLO_func)
        call InitGridConv(grid,mySub_qq_Lmu,PqqNLO_Lmu_func)

        call AllocGridConv(grid,mySub_qq_b)
        call AllocGridConv(grid,mySub_qq_b_Lmu)
        call InitGridConv(grid,mySub_qq_b,PqqNLO_b_func)
        call InitGridConv(grid,mySub_qq_b_Lmu,PqqNLO_b_Lmu_func)

        do iQ = 0, PDFs%nQ
          do iflav = -5,5
            if (iflav == 0) cycle
            xPij_A_1%tab(:,iflav,iQ)     = mySub_qq .conv. PDFs%tab(:,iflav,iQ)
            xPij_A_1_Lmu%tab(:,iflav,iQ) = mySub_qq_Lmu .conv. PDFs%tab(:,iflav,iQ)
            xPij_A_2%tab(:,iflav,iQ)     = xPij_A_1%tab(:,iflav,iQ)
            xPij_A_2_Lmu%tab(:,iflav,iQ) = xPij_A_1_Lmu%tab(:,iflav,iQ)

            xPij_B_1%tab(:,iflav,iQ)     = mySub_qq_b .conv. PDFs%tab(:,iflav,iQ)
            xPij_B_1_Lmu%tab(:,iflav,iQ) = mySub_qq_b_Lmu .conv. PDFs%tab(:,iflav,iQ)
            xPij_B_2%tab(:,iflav,iQ)     = xPij_A_2_Lmu%tab(:,iflav,iQ)

            xPij_C_1%tab(:,iflav,iQ)     = xPij_A_1_Lmu%tab(:,iflav,iQ)
            xPij_C_2%tab(:,iflav,iQ)     = xPij_B_1%tab(:,iflav,iQ)
            xPij_C_2_Lmu%tab(:,iflav,iQ) = xPij_B_1_Lmu%tab(:,iflav,iQ)

            xPij_D_1%tab(:,iflav,iQ) = xPij_B_1_Lmu%tab(:,iflav,iQ)
            xPij_D_2%tab(:,iflav,iQ) = xPij_C_2_Lmu%tab(:,iflav,iQ)

            xPij%tab(:,iflav,iQ)      = mySub_TCqq      .conv. PDFs%tab(:,iflav,iQ)
            xPij_Lmu%tab(:,iflav,iQ)  = mySub_TCqq_Lmu  .conv. PDFs%tab(:,iflav,iQ)
            xPij_Lmu2%tab(:,iflav,iQ) = mySub_TCqq_Lmu2 .conv. PDFs%tab(:,iflav,iQ)
          enddo
        enddo


      elseif(sec(1:3).eq.'s_v') then

        !-- q(qb) qb(q): both s_vqcd and s_vewk
        call AllocGridConv(grid,mySub_ns)
        call AllocGridConv(grid,mySub_ns_Lmu)
        call InitGridConv(grid,mySub_ns,mySub_ns_func)
        call InitGridConv(grid,mySub_ns_Lmu,mySub_ns_Lmu_func)

        do iQ = 0, PDFs%nQ
           do iflav = -5,5
              if (iflav == 0) cycle
              xPij%tab(:,iflav,iQ)     = mySub_ns .conv. PDFs%tab(:,iflav,iQ)
              xPij_Lmu%tab(:,iflav,iQ) = mySub_ns_Lmu .conv. PDFs%tab(:,iflav,iQ)
           enddo
        enddo

      elseif(sec.eq.'s_qqb') then

        call AllocGridConv(grid,mySub_ns_qqb)
        call AllocGridConv(grid,mySub_ns_qqb_Lmu)
        call InitGridConv(grid,mySub_ns_qqb,mySub_ns_qqb_func)
        call InitGridConv(grid,mySub_ns_qqb_Lmu,mySub_ns_qqb_Lmu_func)

        do iQ = 0, PDFs%nQ
           do iflav = -5,5
              if (iflav == 0) cycle
              xPij%tab(:,iflav,iQ) = mySub_ns_qqb .conv. PDFs%tab(:,iflav,iQ)
              xPij_Lmu%tab(:,iflav,iQ) = mySub_ns_qqb_Lmu .conv. PDFs%tab(:,iflav,iQ)
           enddo
        enddo

      elseif(sec.eq.'s_qq') then

        call AllocGridConv(grid,mySub_ns_qq)
        call AllocGridConv(grid,mySub_ns_qq_Lmu)
        call InitGridConv(grid,mySub_ns_qq,mySub_ns_qq_func)
        call InitGridConv(grid,mySub_ns_qq_Lmu,mySub_ns_qq_Lmu_func)

        do iQ = 0, PDFs%nQ
           do iflav = -5,5
              if (iflav == 0) cycle
              xPij%tab(:,iflav,iQ) = mySub_ns_qq .conv. PDFs%tab(:,-iflav,iQ)
              xPij_Lmu%tab(:,iflav,iQ) = mySub_ns_qq_Lmu .conv. PDFs%tab(:,-iflav,iQ)
           enddo
        enddo

      endif

    !!=======================================================================!!
    !!                            gq/qg channels                             !!
    !!=======================================================================!!
    elseif(ch.eq.'gq' .or. ch.eq.'qg') then

      if(sec(1:6).eq.'s_vewk') then

        !-- Pgq \otimes q
        call AllocGridConv(grid,mySub_gq)
        call AllocGridConv(grid,mySub_gq_Lmu)
        call InitGridConv(grid,mySub_gq,mySub_gq_func)
        call InitGridConv(grid,mySub_gq_Lmu,mySub_gq_Lmu_func)
  
        do iQ = 0, PDFs%nQ
          do iflav = -5,5
            if (iflav == 0) cycle
               xPij%tab(:,iflav,iQ) = mySub_gq .conv. PDFs%tab(:,0,iQ)
               xPij_Lmu%tab(:,iflav,iQ) = mySub_gq_Lmu .conv. PDFs%tab(:,0,iQ)
          enddo
        enddo

      elseif(sec.eq.'s_12') then

        ! -------------------------- Pgq \otimes q -------------------------- !
        call AllocGridConv(grid,mySub_gq)
        call AllocGridConv(grid,mySub_gq_Lmu)
        call InitGridConv(grid,mySub_gq,mySub_gq_func)
        call InitGridConv(grid,mySub_gq_Lmu,mySub_gq_Lmu_func)

        call AllocGridConv(grid,mySub_gq_b)
        call AllocGridConv(grid,mySub_gq_b_Lmu)
        call InitGridConv(grid,mySub_gq_b,mySub_gq_b_func)
        call InitGridConv(grid,mySub_gq_b_Lmu,mySub_gq_b_Lmu_func)

        ! -------------------------- Pqq \otimes q -------------------------- !
        call AllocGridConv(grid,mySub_qq)
        call AllocGridConv(grid,mySub_qq_Lmu)
        call InitGridConv(grid,mySub_qq,PqqNLO_func)
        call InitGridConv(grid,mySub_qq_Lmu,PqqNLO_Lmu_func)

        call AllocGridConv(grid,mySub_qq_b)
        call AllocGridConv(grid,mySub_qq_b_Lmu)
        call InitGridConv(grid,mySub_qq_b,PqqNLO_b_func)
        call InitGridConv(grid,mySub_qq_b_Lmu,PqqNLO_b_Lmu_func)

        ! ---------------- Triple Collinear, TC_gq \otimes q ---------------- !
        call AllocGridConv(grid,mySub_TCgq)
        call AllocGridConv(grid,mySub_TCgq_Lmu)
        call AllocGridConv(grid,mySub_TCgq_Lmu2)
        call InitGridConv(grid,mySub_TCgq,mySub_TCgq_func)
        call InitGridConv(grid,mySub_TCgq_Lmu,mySub_TCgq_Lmu_func)
        call InitGridConv(grid,mySub_TCgq_Lmu2,mySub_TCgq_Lmu2_func)

        do iQ = 0, PDFs%nQ
          do iflav = -5,5
            if (iflav == 0) cycle
            xPij_A_1%tab(:,iflav,iQ)     = mySub_gq .conv. PDFs%tab(:,0,iQ)
            xPij_A_1_Lmu%tab(:,iflav,iQ) = mySub_gq_Lmu .conv. PDFs%tab(:,0,iQ)
            xPij_A_2%tab(:,iflav,iQ)     = mySub_qq .conv. PDFs%tab(:,iflav,iQ)
            xPij_A_2_Lmu%tab(:,iflav,iQ) = mySub_qq_Lmu .conv. PDFs%tab(:,iflav,iQ)

            xPij_B_1%tab(:,iflav,iQ)     = mySub_gq_b .conv. PDFs%tab(:,0,iQ)
            xPij_B_1_Lmu%tab(:,iflav,iQ) = mySub_gq_b_Lmu .conv. PDFs%tab(:,0,iQ)
            xPij_B_2%tab(:,iflav,iQ)     = xPij_A_2_Lmu%tab(:,iflav,iQ)

            xPij_C_1%tab(:,iflav,iQ)     = xPij_A_1_Lmu%tab(:,iflav,iQ)
            xPij_C_2%tab(:,iflav,iQ)     = mySub_qq_b .conv. PDFs%tab(:,iflav,iQ)
            xPij_C_2_Lmu%tab(:,iflav,iQ) = mySub_qq_b_Lmu .conv. PDFs%tab(:,iflav,iQ)

            xPij_D_1%tab(:,iflav,iQ) = xPij_B_1_Lmu%tab(:,iflav,iQ)
            xPij_D_2%tab(:,iflav,iQ) = xPij_C_2_Lmu%tab(:,iflav,iQ)

            xPij%tab(:,iflav,iQ) = mySub_TCgq .conv. PDFs%tab(:,0,iQ)
            xPij_Lmu%tab(:,iflav,iQ) = mySub_TCgq_Lmu .conv. PDFs%tab(:,0,iQ)
            xPij_Lmu2%tab(:,iflav,iQ) = mySub_TCgq_Lmu2 .conv. PDFs%tab(:,0,iQ)
          enddo
        enddo

      endif

    !!=======================================================================!!
    !!                            aq/qa channels                             !!
    !!=======================================================================!!
    elseif(ch.eq.'aq' .or. ch.eq.'qa') then

      if(sec.eq.'s_vqcd') then

        ! -------------------------- Paq \otimes q -------------------------- !
        call AllocGridConv(grid,mySub_gq)
        call AllocGridConv(grid,mySub_gq_Lmu)
        call InitGridConv(grid,mySub_gq,mySub_gq_func)
        call InitGridConv(grid,mySub_gq_Lmu,mySub_gq_Lmu_func)
  
        do iQ = 0, PDFs%nQ
          do iflav = -5,5
            if (iflav == 0) cycle
            xPij%tab(:,iflav,iQ) = mySub_gq .conv. PDFs%tab(:,6,iQ)
            xPij_Lmu%tab(:,iflav,iQ) = mySub_gq_Lmu .conv. PDFs%tab(:,6,iQ)
          enddo
        enddo

      elseif(sec.eq.'s_12') then

        ! -------------------------- Pqa \otimes a -------------------------- !
        call AllocGridConv(grid,mySub_TCaa)
        call AllocGridConv(grid,mySub_TCaa_Lmu)
        call AllocGridConv(grid,mySub_TCaa_Lmu2)
        call InitGridConv(grid,mySub_TCaa,mySub_TCaa_func)
        call InitGridConv(grid,mySub_TCaa_Lmu,mySub_TCaa_Lmu_func)
        call InitGridConv(grid,mySub_TCaa_Lmu2,mySub_TCaa_Lmu2_func)

        do iQ = 0, PDFs%nQ
          do iflav = -5,5,2
             if (iflav==0) cycle
             xPij_a%tab(:,6,iQ) = xPij_a%tab(:,6,iQ) + Qdn2 * (mySub_TCaa .conv. PDFs%tab(:,iflav,iQ))
             xPij_a_Lmu%tab(:,6,iQ) = xPij_a_Lmu%tab(:,6,iQ) + Qdn2 * (mySub_TCaa_Lmu .conv. PDFs%tab(:,iflav,iQ))
             xPij_a_Lmu2%tab(:,6,iQ) = xPij_a_Lmu2%tab(:,6,iQ) + Qdn2 * (mySub_TCaa_Lmu2 .conv. PDFs%tab(:,iflav,iQ))
          enddo
          do iflav = -4,4,2
             if (iflav==0) cycle
             xPij_a%tab(:,6,iQ) = xPij_a%tab(:,6,iQ) + Qup2 * (mySub_TCaa .conv. PDFs%tab(:,iflav,iQ))
             xPij_a_Lmu%tab(:,6,iQ) = xPij_a_Lmu%tab(:,6,iQ) + Qup2 * (mySub_TCaa_Lmu .conv. PDFs%tab(:,iflav,iQ))
             xPij_a_Lmu2%tab(:,6,iQ) = xPij_a_Lmu2%tab(:,6,iQ) + Qup2 * (mySub_TCaa_Lmu2 .conv. PDFs%tab(:,iflav,iQ))
          enddo
        enddo

        ! -------------------------- Paq \otimes q -------------------------- !
        call AllocGridConv(grid,mySub_gq)
        call AllocGridConv(grid,mySub_gq_Lmu)
        call InitGridConv(grid,mySub_gq,mySub_gq_func)
        call InitGridConv(grid,mySub_gq_Lmu,mySub_gq_Lmu_func)

        call AllocGridConv(grid,mySub_gq_b)
        call AllocGridConv(grid,mySub_gq_b_Lmu)
        call InitGridConv(grid,mySub_gq_b,mySub_gq_b_func)
        call InitGridConv(grid,mySub_gq_b_Lmu,mySub_gq_b_Lmu_func)

        ! -------------------------- Pqq \otimes q -------------------------- !
        call AllocGridConv(grid,mySub_qq)
        call AllocGridConv(grid,mySub_qq_Lmu)
        call InitGridConv(grid,mySub_qq,PqqNLO_func)
        call InitGridConv(grid,mySub_qq_Lmu,PqqNLO_Lmu_func)

        call AllocGridConv(grid,mySub_qq_b)
        call AllocGridConv(grid,mySub_qq_b_Lmu)
        call InitGridConv(grid,mySub_qq_b,PqqNLO_b_func)
        call InitGridConv(grid,mySub_qq_b_Lmu,PqqNLO_b_Lmu_func)

        ! ---------------- Triple Collinear, TC_aq \otimes q ---------------- !
        call AllocGridConv(grid,mySub_TCgq)
        call AllocGridConv(grid,mySub_TCgq_Lmu)
        call AllocGridConv(grid,mySub_TCgq_Lmu2)
        call InitGridConv(grid,mySub_TCgq,mySub_TCgq_func)
        call InitGridConv(grid,mySub_TCgq_Lmu,mySub_TCgq_Lmu_func)
        call InitGridConv(grid,mySub_TCgq_Lmu2,mySub_TCgq_Lmu2_func)

        do iQ = 0, PDFs%nQ
          do iflav = -5,5
            if (iflav == 0) cycle
            xPij_A_1%tab(:,iflav,iQ)     = mySub_gq .conv. PDFs%tab(:,6,iQ)
            xPij_A_1_Lmu%tab(:,iflav,iQ) = mySub_gq_Lmu .conv. PDFs%tab(:,6,iQ)
            xPij_A_2%tab(:,iflav,iQ)     = mySub_qq .conv. PDFs%tab(:,iflav,iQ)
            xPij_A_2_Lmu%tab(:,iflav,iQ) = mySub_qq_Lmu .conv. PDFs%tab(:,iflav,iQ)

            xPij_B_1%tab(:,iflav,iQ)     = mySub_gq_b .conv. PDFs%tab(:,6,iQ)
            xPij_B_1_Lmu%tab(:,iflav,iQ) = mySub_gq_b_Lmu .conv. PDFs%tab(:,6,iQ)
            xPij_B_2%tab(:,iflav,iQ)     = xPij_A_2_Lmu%tab(:,iflav,iQ)

            xPij_C_1%tab(:,iflav,iQ)     = xPij_A_1_Lmu%tab(:,iflav,iQ)
            xPij_C_2%tab(:,iflav,iQ)     = mySub_qq_b .conv. PDFs%tab(:,iflav,iQ)
            xPij_C_2_Lmu%tab(:,iflav,iQ) = mySub_qq_b_Lmu .conv. PDFs%tab(:,iflav,iQ)

            xPij_D_1%tab(:,iflav,iQ) = xPij_B_1_Lmu%tab(:,iflav,iQ)
            xPij_D_2%tab(:,iflav,iQ) = xPij_C_2_Lmu%tab(:,iflav,iQ)

            xPij%tab(:,iflav,iQ) = mySub_TCgq .conv. PDFs%tab(:,6,iQ)
            xPij_Lmu%tab(:,iflav,iQ) = mySub_TCgq_Lmu .conv. PDFs%tab(:,6,iQ)
            xPij_Lmu2%tab(:,iflav,iQ) = mySub_TCgq_Lmu2 .conv. PDFs%tab(:,6,iQ)
          enddo
        enddo

      endif

    !!=======================================================================!!
    !!                            ag/ga channels                             !!
    !!=======================================================================!!
    elseif(ch.eq.'ag' .or. ch.eq.'ga') then

      if(sec.eq.'s_12') then

        ! ---------------- Triple Collinear, TC_gqa \otimes a --------------- !
        call AllocGridConv(grid,mySub_TCag)
        call AllocGridConv(grid,mySub_TCag_Lmu)
        call AllocGridConv(grid,mySub_TCag_Lmu2)
        call InitGridConv(grid,mySub_TCag,mySub_TCag_func)
        call InitGridConv(grid,mySub_TCag_Lmu,mySub_TCag_Lmu_func)
        call InitGridConv(grid,mySub_TCag_Lmu2,mySub_TCag_Lmu2_func)

        do iQ = 0, PDFs%nQ
          xPij%tab(:,6,iQ) = mySub_TCag .conv. PDFs%tab(:,0,iQ)
          xPij_Lmu%tab(:,6,iQ) = mySub_TCag_Lmu .conv. PDFs%tab(:,0,iQ)
          xPij_Lmu2%tab(:,6,iQ) = mySub_TCag_Lmu2 .conv. PDFs%tab(:,0,iQ)
        enddo

        ! -------------------------- P[a,g]q \otimes q -------------------------- !
        call AllocGridConv(grid,mySub_gq)
        call AllocGridConv(grid,mySub_gq_Lmu)
        call InitGridConv(grid,mySub_gq,mySub_gq_func)
        call InitGridConv(grid,mySub_gq_Lmu,mySub_gq_Lmu_func)

        call AllocGridConv(grid,mySub_gq_b)
        call AllocGridConv(grid,mySub_gq_b_Lmu)
        call InitGridConv(grid,mySub_gq_b,mySub_gq_b_func)
        call InitGridConv(grid,mySub_gq_b_Lmu,mySub_gq_b_Lmu_func)

        do iQ = 0, PDFs%nQ
          do iflav = -5,5
            if (iflav == 0) cycle
            xPij_A_1%tab(:,iflav,iQ)     = mySub_gq .conv. PDFs%tab(:,6,iQ)
            xPij_A_1_Lmu%tab(:,iflav,iQ) = mySub_gq_Lmu .conv. PDFs%tab(:,6,iQ)
            xPij_A_2%tab(:,iflav,iQ)     = mySub_gq .conv. PDFs%tab(:,0,iQ)
            xPij_A_2_Lmu%tab(:,iflav,iQ) = mySub_gq_Lmu .conv. PDFs%tab(:,0,iQ)

            xPij_B_1%tab(:,iflav,iQ)     = mySub_gq_b .conv. PDFs%tab(:,6,iQ)
            xPij_B_1_Lmu%tab(:,iflav,iQ) = mySub_gq_b_Lmu .conv. PDFs%tab(:,6,iQ)
            xPij_B_2%tab(:,iflav,iQ)     = xPij_A_2_Lmu%tab(:,iflav,iQ)

            xPij_C_1%tab(:,iflav,iQ)     = xPij_A_1_Lmu%tab(:,iflav,iQ)
            xPij_C_2%tab(:,iflav,iQ)     = mySub_gq_b .conv. PDFs%tab(:,0,iQ)
            xPij_C_2_Lmu%tab(:,iflav,iQ) = mySub_gq_b_Lmu .conv. PDFs%tab(:,0,iQ)

            xPij_D_1%tab(:,iflav,iQ) = xPij_B_1_Lmu%tab(:,iflav,iQ)
            xPij_D_2%tab(:,iflav,iQ) = xPij_C_2_Lmu%tab(:,iflav,iQ)

          enddo
        enddo
        
      elseif(sec.eq.'s_oqcd') then

        ! -------------------------- Paq \otimes q -------------------------- !
        call AllocGridConv(grid,mySub_gq)
        call AllocGridConv(grid,mySub_gq_Lmu)
        call InitGridConv(grid,mySub_gq,mySub_gq_func)
        call InitGridConv(grid,mySub_gq_Lmu,mySub_gq_Lmu_func)
  
        do iQ = 0, PDFs%nQ
          do iflav = -5,5
            if (iflav == 0) cycle
            xPij%tab(:,iflav,iQ) = mySub_gq .conv. PDFs%tab(:,6,iQ)
            xPij_Lmu%tab(:,iflav,iQ) = mySub_gq_Lmu .conv. PDFs%tab(:,6,iQ)
          enddo
        enddo

      endif

    endif

  end subroutine init_xPij_nnlo

  !!*************************************************************************!!

  !---------------------------------------------------------------------------!
  function mySub_ns_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,omx,lx,lomx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       omx = one-x
       lx = log(x)
       lomx = log(omx)
       !-- reg
       res = (-two*(one+x)*lomx-(one+x**2)/omx*lx+(one-x))
       !-- plus
       res = res + (four*lomx)/omx
    end select

    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)       
       omx = one-x
       lomx = log(omx)
       res = res -((four*lomx)/omx)
    case(cc_DELTA)
       res = 2*zeta2   !-- half -> one per leg
    end select

    if (cc_piece /= cc_DELTA) res = res * x

  end function mySub_ns_func

  function mySub_ns_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
      xb = one-x
      res = -(one+x**2)/xb
    end select

    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
      xb = one-x
      res = res + 2/xb
    case(cc_DELTA)
      res = -1.5_dp     !-- half -> one per leg
    end select

    if (cc_piece /= cc_DELTA) res = res * x

  end function mySub_ns_Lmu_func

  function PqqNLO_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,omx,lx,lomx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       omx = one-x
       lx = log(x)
       lomx = log(omx)
       !-- reg
       res = omx-two*(one+x)*lomx+(one+x-2/omx)*lx
       !-- plus
       res = res + (four*lomx)/omx
    end select

    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)       
       omx = one-x
       lomx = log(omx)
       res = res - (four*lomx)/omx
    case(cc_DELTA)
    end select

    if (cc_piece /= cc_DELTA) res = res * x

  end function PqqNLO_func

  function PqqNLO_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,omx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       omx = one-x
       !-- reg
       res = one + x
       !-- plus
       res = res - 2/omx
    end select

    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)       
       omx = one-x
       res = res + 2/omx
    case(cc_DELTA)
    end select

    if (cc_piece /= cc_DELTA) res = res * x

  end function PqqNLO_Lmu_func

  function PqqNLO_b_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,omx,lx,lomx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       omx = one-x
       lx = log(x)
       lomx = log(omx)

       res = omx-two*(one+x)*lomx+(one+x-2/omx)*lx+(four*lomx)/omx
       res = res * lx

    end select

    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)       
    case(cc_DELTA)
    end select

    if (cc_piece /= cc_DELTA) res = res * x

  end function PqqNLO_b_func

  function PqqNLO_b_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,omx,lx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       omx = one-x
       lx = log(x)

       res = one+x-2/omx
       res = res * lx

    end select

    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
    case(cc_DELTA)
    end select

    if (cc_piece /= cc_DELTA) res = res * x

  end function PqqNLO_b_Lmu_func

  function mySub_TCqq_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,omx,opx
    real(dp)             :: lx,lomx,li2x,li3x,li3xb
    real(dp)             :: D0,D1,D2,D3

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       opx = one+x
       omx = one-x
       lomx = log(omx)
       D0 = one/omx
       D1 = lomx*D0
       D2 = lomx*D1
       D3 = lomx*D2
       lx = log(x)
       li2x  = real(dilog2(x),kind=dp) 
       li3x  = real(trilog(x),kind=dp) 
       li3xb = real(trilog(omx),kind=dp) 
       !-- reg
       res = (-16*lomx**2*lx)/omx-(4*lx**3)/omx+12*li3x*opx &
           - 6*li3xb*opx-8*lomx**3*opx-(13*lx**3*opx)/6 &
           + li2x*(-6*omx-2*lomx*opx-4*lx*opx)+lomx**2*(-4*omx+lx*opx)+(4*lx*(-3*lx+4*lomx*lx&
           + 2*(li2x-zeta2)))/omx+6*(3*omx+(5*lomx*opx)/3-(8*lx*opx)/3)*zeta2-28*opx*zeta3&
           + (12*lx**2*((-4*lx)/9+1._dp))/omx+lx**2*(-x+2._dp)+lx**2*(-4*lomx*opx+3*lx*opx-4*(x&
           + 2._dp))+lomx*(4*lx**2*opx-8*lx*(x+2._dp)-3*(3*x-4._dp))+lx*(-9*x+5._dp)&
           + lx*((-24*lomx)/omx-(24*lomx**2)/omx+12*lomx**2*opx-lx**2*opx+(8*zeta2)/omx-4*opx*zeta2&
           + 4*lx*(x+2._dp)+lomx*(-8*lx*opx+8*(x+2._dp))+2*(5*x-6._dp))+2*(8*x-9._dp) &
           + (10*lx**3-4*li3xb+12*lomx*(2*lx-li2x+zeta2)+8*lx*(li2x+(3*lomx**2)/4&
           + 4*zeta2)-32*(li3x-zeta3))/omx
       !-- plus
       res = res + 16*(D3 - zeta2*D1 + 2*zeta3*D0)
    end select

    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       omx = one-x
       lomx = log(omx)
       omx = one-x
       lomx = log(omx)
       D0 = one/omx
       D1 = lomx*D0
       D2 = lomx*D1
       D3 = lomx*D2
       !-- subtract plus
       res = res - 16*(D3 - zeta2*D1 + 2*zeta3*D0)
    case(cc_DELTA)
       res = -2*pisq**2/45 !-- half per leg
    end select

    if (cc_piece /= cc_DELTA) res = res * x

  end function mySub_TCqq_func

  function mySub_TCqq_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,omx,opx,lx,lomx,li2x
    real(dp)             :: D0,D1,D2

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       opx = one+x
       omx = one-x
       lomx = log(omx)
       D0 = one/omx
       D1 = lomx*D0
       D2 = lomx*D1
       lx = log(x)
       li2x = real(dilog2(x),kind=dp) 
       !-- reg
       res = 4*(-3*lx+4*lomx*lx+2*(li2x-zeta2))/omx -8*lx**2/omx &
           + 2*(-4*(two+x) + 12/omx + 8*lomx/omx + (3*lx-4*lomx)*opx)*lx &
           - 4*zeta2*opx + 2*(-6._dp + 5*x) + 12*opx*lomx**2 + 4*(two + x)*lx &
           - opx*lx**2 + 8*lomx*(two+x-opx*lx)
       !-- plus
       res = res + 8*(zeta2*D0-3*D1-3*D2)
    end select

    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       omx = one-x
       lomx = log(omx)
       D0 = one/omx
       D1 = lomx*D0
       D2 = lomx*D1
       !-- subtract plus
       res = res - 8*(zeta2*D0-3*D1-3*D2)
    case(cc_DELTA)
       res = -2*(pisq + 8*zeta3) !-- half per leg
    end select

    if (cc_piece /= cc_DELTA) res = res * x

  end function mySub_TCqq_Lmu_func

  function mySub_TCqq_Lmu2_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,omx,opx,lx,lomx
    real(dp)             :: D0,D1
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       omx = one-x
       lomx = log(omx)
       D0 = one/omx
       D1 = lomx*D0
       lx = log(x)
       opx = one + x
       !-- reg
       res = -4*(lx/omx + two + x + opx*lomx) + 3*opx*lx
       !-- plus
       res = res + 12*D0 + 8*D1
    end select

    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       omx = one-x
       lomx = log(omx)
       D0 = one/omx
       D1 = lomx*D0
       !-- subtract plus
       res = res - 12*D0 - 8*D1
    case(cc_DELTA)
       res = 4.5_dp - 4*zeta2 !-- half per leg
    end select

    if (cc_piece /= cc_DELTA) res = res * x

  end function mySub_TCqq_Lmu2_func

  !---------------------------------------------------------------------------!
  function mySub_ns_qqb_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xp,xb
    real(dp) :: lx,lxp,lxb,li2x,li2mx,li3x,li3xb,li3mx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      xp = one+x
      lx = log(x)
      lxp = log(xp)
      lxb = log(xb)
      li2x = real(dilog2(x),kind=dp) 
      li2mx = real(dilog2(-x),kind=dp)
      li3x = trilog(x)
      li3xb = trilog(xb)
      li3mx = trilog(-x)

      !-- reg
      res = 7*x+6*lx*lxb*xb+(8*li2mx+8*lx*lxp)*xp+1._dp+((16*li3mx+18*li3x+12*li3xb-8*li2mx*lx&
      -10*li2x*lx-lx**3/3+8*li2x*lxb-5*lx**2*lxb+8*lx*lxb**2-4*lx*zeta2-8*lxb*zeta2-&
      6*zeta3)*(x**2+1._dp))/xb-(lx**2*(2*x-1._dp)*(2*x-5._dp))/(2*xb)+(lx*(-11*x+27*x**2-&
      6._dp))/xb+4*lxb*(7*x-8._dp)-(2*zeta2*(-6*x+x**2+11._dp))/xb-(2*li2x*(6*x+x**2-&
      13._dp))/xb

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)

    case(cc_DELTA)
    end select

    if (cc_piece /= cc_DELTA) res = res * x

  end function mySub_ns_qqb_func

  function mySub_ns_qqb_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb
    real(dp) :: lx,lxb,li2x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      lx = log(x)
      lxb = log(xb)
      li2x = real(dilog2(x),kind=dp)

      !-- reg
      res = ((-4*li2x+2*lx**2-4*lx*lxb+4*zeta2)*(x**2+1._dp))/xb-&
        (2*lx*(2*x**2-5._dp))/xb-2*(7*x-8._dp)

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x 

  end function mySub_ns_qqb_Lmu_func

  !---------------------------------------------------------------------------!
  function mySub_ns_qq_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xp,xb
    real(dp) :: lx,lxp,lxb,li2x,li2mx,li3x,li3xb,li3mx,li3xoxp,li3xbxp
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      xp = one+x
      lx = log(x)
      lxp = log(xp)
      lxb = log(xb)
      li2x = real(dilog2(x),kind=dp) 
      li2mx = real(dilog2(-x),kind=dp)
      li3x = trilog(x)
      li3xb = trilog(xb)
      li3mx = trilog(-x)
      li3xoxp = trilog(x/xp)
      li3xbxp = trilog(xb*xp)

      !-- reg
      res = 4*lx**2+(-12*li2mx-12*lx*lxp)*xp+((-36*li3mx-16*li3x+16*li3xb-8*li3xbxp-24*li3xoxp+&
      12*li2mx*lx+4*li2x*lx+lx**3/3-16*li2mx*lxb+6*lx**2*lxp-16*lx*lxb*lxp-12*lx*lxp**2+&
      4*lxp**3+12*lx*zeta2-8*lxb*zeta2-12*lxp*zeta2+10*zeta3)*(x**2+1._dp))/xp+2*zeta2*(-x+&
      3._dp)-4*li2x*(x+3._dp)+lx*(19*x+11._dp)+xb*(16*lxb-4*lx*lxb+15._dp)

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_ns_qq_func

  function mySub_ns_qq_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xp,xb
    real(dp) :: lx,lxp,li2mx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      xp = one+x
      lx = log(x)
      lxp = log(xp)
      li2mx = real(dilog2(-x),kind=dp)

      !-- reg
      res = -8*xb-4*lx*xp+((8*li2mx-2*lx**2+8*lx*lxp+4*zeta2)*(x**2+1._dp))/xp

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)

    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_ns_qq_Lmu_func

  !---------------------------------------------------------------------------!

  function mySub_gq_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb
    real(dp) :: lx,lxb
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      lx = log(x)
      lxb = log(xb)

      !-- reg
      res = 2*xb*x + (xb**2 + x**2)*(2*lxb - lx)

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_gq_func

  function mySub_gq_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      !-- reg
      res = -(xb**2 + x**2)

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_gq_Lmu_func

  !-- includes extra log(z)
  function mySub_gq_b_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb
    real(dp) :: lx,lxb
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      lx = log(x)
      lxb = log(xb)

      !-- reg
      res = 2*xb*x + (xb**2 + x**2)*(2*lxb - lx)
      res = res * lx

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_gq_b_func

  function mySub_gq_b_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb,lx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      lx = log(x)
      !-- reg
      res = -(xb**2 + x**2)
      res = res * lx

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_gq_b_Lmu_func

  function mySub_TCgq_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb
    real(dp) :: lx,lxb,li2x,li2xb,li3x,li3xb
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      lx = log(x)
      lxb = log(xb)
      li2x  = real(dilog2(x),kind=dp) 
      li2xb = real(dilog2(xb),kind=dp) 
      li3x  = trilog(x)
      li3xb = trilog(xb)

      !-- reg
      res = 3*li2x-8*li2x*lxb*(x-1._dp)**2+(2*li2xb*lx-3*li2xb*lxb+(5*lx**2*lxb)/2+&
      (11*lxb**3)/6)*(-2*x+2*x**2+1._dp)+2*lx*zeta2*(-2*x+4*x**2+1._dp)+(lx**3*(-6*x+4*x**2&
      +3._dp))/12+(lx**2*(60*x-48*x**2-5._dp))/8+(li2x*lx+2*lxb*zeta2)*(-10*x+6*x**2+&
      5._dp)+lxb**2*(21*x-17*x**2-7._dp)+lx*lxb*(-22*x+18*x**2+7._dp)+li3x*(18*x-14*x**2&
      -9._dp)+li3xb*(18*x-10*x**2-9._dp)+2*zeta3*(-18*x+16*x**2+9._dp)+(zeta2*(-26*x+&
      18*x**2+11._dp))/2+(lx*lxb**2*(34*x-26*x**2-17._dp))/2+(lx*(117*x-92*x**2-28._dp))/4&
      +(lxb*(-109*x+88*x**2+32._dp))/2+(255*x-196*x**2-69._dp)/4

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)

    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x

  end function mySub_TCgq_func

  function mySub_TCgq_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb
    real(dp) :: lx,lxb,li2x,li2xb
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      lx = log(x)
      lxb = log(xb)
      li2x  = real(dilog2(x),kind=dp) 
      li2xb = real(dilog2(xb),kind=dp) 

      !-- reg
      res = 4*li2x*(x-1._dp)**2+lx*(6*x-4*x**2+1._dp)+(2*li2xb+6*lx*lxb-3*lxb**2)*(-2*x+2*x**2&
      +1._dp)+lxb*(-8*x+4*x**2+1._dp)+(lx**2*(6*x-8*x**2-3._dp))/2+2*zeta2*(6*x-4*x**2-&
      3._dp)+(41*x-30*x**2-16._dp)/2

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_TCgq_Lmu_func

  function mySub_TCgq_Lmu2_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb
    real(dp) :: lx,lxb
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      lx = log(x)
      lxb = log(xb)

      !-- reg
      res = (lx*(2*x-4*x**2-1._dp))/2+lxb*(-2*x+2*x**2+1._dp)+(-8*x+12*x**2+5._dp)/4

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_TCgq_Lmu2_func

  !---------------------------------------------------------------------------!
  function mySub_TCaa_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb
    real(dp) :: lx,lxb,li2x,li3x,li3xb
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      lx = log(x)
      lxb = log(xb)
      li2x = real(dilog2(x),kind=dp) 
      li3x = trilog(x)
      li3xb = trilog(xb)

      !-- reg
      res = (lx**2*(-four - 5*x))/8 + li2x*(-8._dp - 2*lx*(-two + x) + 4*lxb*(-two + x) - 5*x) + &
      (lx**3*(two - x))/12 + 2*li3x*(-two + x) + 4*li3xb*(-two + x) + &
      (13*lxb**3*(two - 2*x + x**2))/(6*x) - (108._dp - 146*x + 29*x**2)/(4*x) + &
      lxb**2*(-((lx*(6._dp - 2*x + x**2))/x) - (42._dp - 58*x + 27*x**2)/(4*x)) + &
      ((30._dp - 22*x + 27*x**2)*zeta2)/(2*x) + lx*((-20._dp + 27*x + 13*x**2)/(4*x) + &
      (2*(four - 6*x + 3*x**2)*zeta2)/x) + lxb*((2*lx*(-one + x)*(-3._dp + 2*x))/x + &
      (lx**2*(two - 2*x + x**2))/x - (-36._dp + 46*x + 3*x**2)/(2*x) - &
      (2*(10._dp - 14*x + 7*x**2)*zeta2)/x) + (2*(8._dp - 6*x + 3*x**2)*zeta3)/x

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_TCaa_func

  function mySub_TCaa_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb
    real(dp) :: lx,lxb,li2x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      lx = log(x)
      lxb = log(xb)
      li2x = real(dilog2(x),kind=dp) 

      !-- reg
      res = lx*(-two - 3*x) + (lx**2*(two - x))/2 - 2*li2x*(-two + x) - &
      (3*lxb**2*(two - 2*x + x**2))/x + (-10._dp + 15*x + 4*x**2)/(2*x) + &
      lxb*((2*(-one + x)*(-3._dp + 2*x))/x + (2*lx*(two - 2*x + x**2))/x) + &
      (2*(four - 6*x + 3*x**2)*zeta2)/x

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_TCaa_Lmu_func

  function mySub_TCaa_Lmu2_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb
    real(dp) :: lx,lxb
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      lx = log(x)
      lxb = log(xb)

      !-- reg
      res = (four-x)/4 + (two-2*x+x**2)*lxb/x + (two-x)*lx/2

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_TCaa_Lmu2_func

  !---------------------------------------------------------------------------!
  function mySub_TCag_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb
    real(dp) :: lx,lxb,li2x,li3x,li3xb
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      lx = log(x)
      lxb = log(xb)
      li2x = real(dilog2(x),kind=dp) 
      li3x = trilog(x)
      li3xb = trilog(xb)

      !-- reg
      res = 4*li2x*(one + 3*x) - (4*lx*(-10._dp - 33*x + 6*x**2 + 19*x**3))/(9*x) + &
      (4*lxb*(-lx + lxb)*(four + 7*x + 4*x**2)*xb)/(3*x) + &
      (lx**2*(8._dp + 11*x + 8*x**2)*xb)/(6*x) - (4*lxb*(20._dp + 77*x + 38*x**2)*xb)/(9.*x) + &
      (2*(103._dp + 631*x + 211*x**2)*xb)/(27*x) + &
      (8*(-two - 3*x - 3*x**2 + 2*x**3)*zeta2)/(3*x) + &
      (one + x)*(-8*li3x - 16*li3xb - 8*li2x*(-lx + 2*lxb) + 8*(-lx + 2*lxb)*zeta2 + &
      (lx**3 - 24*lx*lxb**2 + 24*zeta3)/3)

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_TCag_func

  function mySub_TCag_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb,lx,lxb,li2x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      lx   = log(x)
      lxb  = log(xb)
      li2x = real(dilog2(x),kind=dp) 

      !-- reg
      res = (-4*lx*(-two - 3*x - 3*x**2 + 2*x**3))/(3*x) - &
      (4*lxb*(four + 7*x + 4*x**2)*xb)/(3*x) + (2*(20._dp + 77*x + 38*x**2)*xb)/(9*x) +&
      (one + x)*(8*li2x + 2*lx**2 - 8*zeta2)

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)

    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x

  end function mySub_TCag_Lmu_func

  function mySub_TCag_Lmu2_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x,xb,lx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

      xb = one-x
      lx = log(x)

      !-- reg
      res = 2*lx*(one + x) + ((four + 7*x + 4*x**2)*xb)/(3*x)

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)

    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_TCag_Lmu2_func

end module mod_hoppet_nnlo

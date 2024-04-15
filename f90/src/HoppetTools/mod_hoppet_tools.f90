module mod_hoppet_tools
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use mod_types, except => dp
  use consts_dp
  use convolution
  implicit none
  private
  type(grid_def),  public, save :: grid
  type(pdf_table), public, save :: PDFs
  !
  public :: init_hoppet,init_pdfs
  public :: get_pdf_hoppet_qed
  
  
contains

  subroutine init_hoppet
    interface
       subroutine evolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine evolvePDF
    end interface
    !--------------------------------------------------
    !-- test against LHAPDF
    !-- real(dp) :: hf(-6:6),lf(-6:7)
    !-- end test
    
    call init_grid()
    call init_pdfs(PDFs)

    call FillPdfTable_LHAPDF(PDFs,evolvePDF_photon_is_six)

    !-- test against LHAPDF
    ! call EvalPdfTable_xQ(PDFs,0.123_dp,91._dp,hf)
    ! call evolvePDFphoton(0.123_dp,91._dp,lf(-6:6),lf(7))
    ! do iq=-6,6
    !    if (iq.eq.6) then
    !       print *, 'photon',lf(7),hf(6)
    !    else
    !       print *, 'iq',iq,lf(iq),hf(iq)
    !    endif
    ! enddo
    ! stop
    !-- end test
    
  end subroutine init_hoppet

  !-- put the photon in position 6
  subroutine evolvePDF_photon_is_six(x,Q,res)
    use types; implicit none
    real(dp), intent(in)  :: x,Q
    real(dp), intent(out) :: res(*)
    real(dp) :: resA

    !-- the above does not work on linux with new LHAPDF
    !call evolvePDFphoton(x,Q,res(1:13),res(14))
    !res(13) = res(14)

    call evolvePDFphoton(x,Q,res(1:13),resA)
    res(13) = resA

  end subroutine evolvePDF_photon_is_six
  
  subroutine init_grid()
    type(grid_def) ::  gdarray(4) ! grid
    integer :: pdf_fit_order
    real(dp) :: dy, ymax, dlnlnQ
    ! build the PDF grid 
    dy = 0.10_dp
    dlnlnQ = dy/4.0_dp
    ymax = 15.0_dp
    pdf_fit_order = -5

    call SetDefaultConvolutionEps(1E-5_dp)

    call InitGridDef(gdarray(4),dy/27.0_dp, 0.2_dp, order=pdf_fit_order)
    call InitGridDef(gdarray(3),dy/9.0_dp,  0.5_dp, order=pdf_fit_order)
    call InitGridDef(gdarray(2),dy/3.0_dp,  2.0_dp, order=pdf_fit_order)
    call InitGridDef(gdarray(1),dy,         ymax  , order=pdf_fit_order)
    call InitGridDef(grid,gdarray(1:4),locked=.true.)
    
  end subroutine init_grid

  subroutine init_pdfs(myPDFs,Qmax)
    real(dp), intent(in), optional :: Qmax
    type(pdf_table) :: myPDFs
    real(dp) :: dy,dlnlnQ

    !-- as in the grid
    dy = 0.10_dp
    dlnlnQ = dy/4.0_dp

    if (present(Qmax)) then
       call AllocPdfTable(grid, myPDFs, 1d0, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    else
       call AllocPdfTable(grid, myPDFs, 1d0, 2d4, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    end if

  end subroutine init_pdfs

  !--

  subroutine get_pdf_hoppet_qed(myPDFs1,myPDFs2,x1,x2,muF,f1,f2,myPDFs1_Lmu,myPDFs2_Lmu,mu2ref)
    type(pdf_table),           intent(in) :: myPDFs1,myPDFs2
    real(dp), intent(in)  :: x1,x2,muF
    real(dp), intent(out) :: f1(-6:7),f2(-6:7)
    type(pdf_table), optional, intent(in) :: myPDFs1_Lmu(:),myPDFs2_Lmu(:)
    real(dp), optional, intent(in) :: mu2ref
    real(dp) :: f1_Lmu(-6:7),f2_Lmu(-6:7),logfact
    integer :: k

    call EvalPdfTable_xQ(myPDFs1,x1,muF,f1)
    call EvalPdfTable_xQ(myPDFs2,x2,muF,f2)

    call fix_pdf(x1,f1)
    call fix_pdf(x2,f2)

    logfact = one
    if (present(myPDFs1_Lmu)) then
       do k=1, size(myPDFs1_Lmu)
          call EvalPdfTable_xQ(myPDFs1_Lmu(k),x1,muF,f1_Lmu)
          call fix_pdf(x1,f1_Lmu)
          logfact = logfact*log(muF**2/mu2ref)
          f1 = f1 + f1_Lmu*logfact
       enddo
    endif

    logfact = one
    if (present(myPDFs2_Lmu)) then
       do k=1, size(myPDFs2_Lmu)
         call EvalPdfTable_xQ(myPDFs2_Lmu(k),x2,muF,f2_Lmu)
         call fix_pdf(x2,f2_Lmu)
         logfact = logfact*log(muF**2/mu2ref)
         f2 = f2 + f2_Lmu*logfact
       enddo
    endif

     contains

       subroutine fix_pdf(xi,fi)
         real(dp), intent(in)    :: xi
         real(dp), intent(inout) :: fi(-6:7)

         fi = fi/xi

         !-- we are using the top as a placeholder for the photon
         fi(7) = fi(6)

       end subroutine fix_pdf
       
  end subroutine get_pdf_hoppet_qed

end module mod_hoppet_tools

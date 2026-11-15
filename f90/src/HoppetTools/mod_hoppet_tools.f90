module mod_hoppet_tools
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use mod_types, except => dp
  use consts_dp
  use convolution
  use convolution_communicator
  use mod_auxfunctions, only: dilog2
  use mod_parms, only: pdfscale,rge,rgescale,aem,Qsq_IS
  implicit none
  private
  type(grid_def),  public, save :: grid
  type(pdf_table), public, save :: PDFs
  !
  public :: init_hoppet,init_pdfs,rge_shift_pdfs
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
    use types; use mod_parms, only: cptest; implicit none
    real(dp), intent(in)  :: x,Q
    real(dp), intent(out) :: res(*)
    real(dp) :: resA

    !-- the above does not work on linux with new LHAPDF
    !call evolvePDFphoton(x,Q,res(1:13),res(14))
    !res(13) = res(14)

    call evolvePDFphoton(x,merge(pdfscale,Q,pdfscale > 0),res(1:13),resA)
    if (cptest) res(1:13) = res(13:1:-1)
    res(13) = resA

  end subroutine evolvePDF_photon_is_six
  
  !-- RGE test: quark PDFs -> PDFs + sgn(Q) rgescale (splitting function x PDFs), sgn = -1, 0, 1 at Q = pdfscale/2, 1, 2
  subroutine rge_shift_pdfs()
    type(grid_conv) :: K
    real(dp) :: as2pi,a2pi,alphasPDF,sgn
    real(dp), allocatable :: add(:),orig(:,:)
    integer  :: iQ,q

    if (pdfscale <= 0) stop 'rge_shift_pdfs: -rge needs -pdfscale'
    as2pi = alphasPDF(pdfscale)/twopi
    a2pi  = aem/twopi
    allocate(add(lbound(PDFs%tab,1):ubound(PDFs%tab,1)),orig(lbound(PDFs%tab,1):ubound(PDFs%tab,1),-5:5))

    select case (rge)
    case ('qg','qa');     call InitGridConv(grid,K,rge_Pqg)
    case ('qg11','qa11'); call InitGridConv(grid,K,rge_Pqg11)
    case ('qq','qqs');    call InitGridConv(grid,K,rge_Pqq)
    case ('qq11');        call InitGridConv(grid,K,rge_Pqq11)
    case ('qb11');        call InitGridConv(grid,K,rge_Pqqb11)
    case ('qq11d')
    case default;  stop 'rge_shift_pdfs: unknown -rge'
    end select

    do iQ = 0,PDFs%nQ
       sgn = rgescale*log(PDFs%Q_vals(iQ)/pdfscale)/log(two)
       orig = PDFs%tab(:,-5:5,iQ)
       if (rge == 'qg' .or. rge == 'qg11') add = K .conv. PDFs%tab(:,0,iQ)
       if (rge == 'qa' .or. rge == 'qa11') add = K .conv. PDFs%tab(:,6,iQ)   !-- photon
       do q = -5,5
          if (q == 0) cycle
          if (rge == 'qg')   PDFs%tab(:,q,iQ) = PDFs%tab(:,q,iQ) + sgn*as2pi*half*add
          if (rge == 'qg11') PDFs%tab(:,q,iQ) = PDFs%tab(:,q,iQ) + sgn*as2pi*a2pi*Qsq_IS(q)*add
          if (rge == 'qq')   PDFs%tab(:,q,iQ) = PDFs%tab(:,q,iQ) + sgn*a2pi*Qsq_IS(q)*(K .conv. PDFs%tab(:,q,iQ))
          !-- photon -> q: TR -> Nc Q^2; QCD Pqq: Q^2 -> CF
          if (rge == 'qa')   PDFs%tab(:,q,iQ) = PDFs%tab(:,q,iQ) + sgn*a2pi*3*Qsq_IS(q)*add
          if (rge == 'qa11') PDFs%tab(:,q,iQ) = PDFs%tab(:,q,iQ) + sgn*as2pi*a2pi*8*Qsq_IS(q)*add
          if (rge == 'qqs')  PDFs%tab(:,q,iQ) = PDFs%tab(:,q,iQ) + sgn*as2pi*4/three*(K .conv. PDFs%tab(:,q,iQ))
          if (rge == 'qq11') PDFs%tab(:,q,iQ) = PDFs%tab(:,q,iQ) + sgn*as2pi*a2pi*8/three*Qsq_IS(q)*(K .conv. PDFs%tab(:,q,iQ))
          if (rge == 'qb11') PDFs%tab(:,q,iQ) = PDFs%tab(:,q,iQ) + sgn*as2pi*a2pi*8/three*Qsq_IS(q)*(K .conv. orig(:,-q))
          !-- delta(1-x) part of qq11 alone: minus the scale slope of the two-loop sector
          if (rge == 'qq11d') PDFs%tab(:,q,iQ) = PDFs%tab(:,q,iQ) + sgn*as2pi*a2pi*8/three*Qsq_IS(q)*(3/8._dp - pi**2/2 + 6*zeta3)*orig(:,q)
       enddo
    enddo

  end subroutine rge_shift_pdfs

  function rge_Pqg(y) result(res)
    real(dp), intent(in) :: y
    real(dp) :: res,x
    x = exp(-y)
    res = zero
    if (cc_piece == cc_REAL .or. cc_piece == cc_REALVIRT) res = (x**2 + (one-x)**2)*x
  end function rge_Pqg

  !-- (1+x^2)/(1-x)_+ + 3/2 delta(1-x)
  function rge_Pqq(y) result(res)
    real(dp), intent(in) :: y
    real(dp) :: res,x
    x = exp(-y)
    res = zero
    select case (cc_piece)
    case (cc_REAL);     res = (one+x**2)/(one-x)*x
    case (cc_VIRT);     res = -two/(one-x)*x
    case (cc_REALVIRT); res = ((one+x**2)/(one-x) - two/(one-x))*x
    case (cc_DELTA);    res = 1.5_dp
    end select
  end function rge_Pqq

  !-- P_qq^(1,1)/(2 CF Q^2): CF^2 part of the two-loop non-singlet splitting function
  function rge_Pqq11(y) result(res)
    real(dp), intent(in) :: y
    real(dp) :: res,x,lx,l1
    x = exp(-y)
    res = zero
    select case (cc_piece)
    case (cc_REAL,cc_REALVIRT)
       lx = log(x); l1 = log(one-x)
       res = (-(2*lx*l1 + 1.5_dp*lx)*(one+x**2)/(one-x) - (1.5_dp + 3.5_dp*x)*lx - half*(one+x)*lx**2 - 5*(one-x))*x
    case (cc_DELTA)
       res = 3/8._dp - pi**2/2 + 6*zeta3
    end select
  end function rge_Pqq11

  !-- P_qqb^(1,1)/(2 CF Q^2) = 2 p_qq(-x) S2(x) + 2 (1+x) log(x) + 4 (1-x)
  function rge_Pqqb11(y) result(res)
    real(dp), intent(in) :: y
    real(dp) :: res,x,lx,s2
    x = exp(-y)
    res = zero
    if (cc_piece == cc_REAL .or. cc_piece == cc_REALVIRT) then
       lx = log(x)
       s2 = -2*real(dilog2(-x),kind=dp) - 2*lx*log(one+x) + half*lx**2 - pi**2/6
       res = (2*(one+x**2)/(one+x)*s2 + 2*(one+x)*lx + 4*(one-x))*x
    endif
  end function rge_Pqqb11

  !-- P_qg^(1,1)/Q^2 (de Florian, Sborlini, Rodrigo), normalisation as/2pi a/2pi
  function rge_Pqg11(y) result(res)
    real(dp), intent(in) :: y
    real(dp) :: res,x,lx,l1,lr
    x = exp(-y)
    res = zero
    if (cc_piece == cc_REAL .or. cc_piece == cc_REALVIRT) then
       lx = log(x); l1 = log(one-x); lr = l1 - lx
       res = half*half*(4 - 9*x - (one-4*x)*lx - (one-2*x)*lx**2 + 4*l1 &
            + (2*lr**2 - 4*lr - two*pi**2/3 + 10)*(x**2 + (one-x)**2))*x
    endif
  end function rge_Pqg11

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
          logfact = logfact*log(merge(pdfscale,muF,rge /= 'none')**2/mu2ref)
          f1 = f1 + f1_Lmu*logfact
       enddo
    endif

    logfact = one
    if (present(myPDFs2_Lmu)) then
       do k=1, size(myPDFs2_Lmu)
         call EvalPdfTable_xQ(myPDFs2_Lmu(k),x2,muF,f2_Lmu)
         call fix_pdf(x2,f2_Lmu)
         logfact = logfact*log(merge(pdfscale,muF,rge /= 'none')**2/mu2ref)
         f2 = f2 + f2_Lmu*logfact
       enddo
    endif

     contains

       subroutine fix_pdf(xi,fi)
         real(dp), intent(in)    :: xi
         real(dp), intent(inout) :: fi(-6:7)

         fi(-6:6) = fi(-6:6)/xi

         !-- we are using the top as a placeholder for the photon
         fi(7) = fi(6)

       end subroutine fix_pdf
       
  end subroutine get_pdf_hoppet_qed

end module mod_hoppet_tools

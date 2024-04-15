module mod_hoppet_nnlo_fc
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
  type(pdf_table), public, save :: xPij,xPij_Lmu
  type(pdf_table), public, save :: xPij_2,xPij_2_Lmu

  public :: init_xPij_nnlo_fc
  
contains

  subroutine init_xPij_nnlo_fc()
    integer :: iQ,iflav
    type(grid_conv) :: mySub_ns_qqb,mySub_ns_qqb_Lmu
    type(grid_conv) :: mySub_ns_qq,mySub_ns_qq_Lmu
    type(grid_conv) :: mySub_qg,mySub_qg_Lmu

    call init_pdfs(xPij)
    call init_pdfs(xPij_Lmu)
    !
    call init_pdfs(xPij_2)
    call init_pdfs(xPij_2_Lmu)

    !-- use:
    !-- call AllocGridConv(grid,myPij)
    !-- call InitGridConv(grid,myPij,function_to_load)
    !-- call AddWithCoeff(myPij_1,myPij_2)                --> myPij_1 = myPij_1 + myPij_2
    !-- call AddWithCoeff(myPij_1,myPij_2,const)          --> myPij_1 = myPij_1 + myPij_2*const
    !-- call SetToConvolutionn(myPij_1x2,myPij_1,myPij_2) --> myPij_1x2 = myPij_1 \otimes myPij_2
    
    call AllocGridConv(grid,mySub_ns_qqb)
    call AllocGridConv(grid,mySub_ns_qqb_Lmu)
    !
    call AllocGridConv(grid,mySub_ns_qq)
    call AllocGridConv(grid,mySub_ns_qq_Lmu)
    !
    call AllocGridConv(grid,mySub_qg)
    call AllocGridConv(grid,mySub_qg_Lmu)
    
    call InitGridConv(grid,mySub_ns_qqb,mySubNNLO_ns_qqb_func)
    call InitGridConv(grid,mySub_ns_qqb_Lmu,mySubNNLO_ns_qqb_Lmu_func)
    !
    call InitGridConv(grid,mySub_ns_qq,mySubNNLO_ns_qq_func)
    call InitGridConv(grid,mySub_ns_qq_Lmu,mySubNNLO_ns_qq_Lmu_func)
    !
    call InitGridConv(grid,mySub_qg,mySub_qg_func)
    call InitGridConv(grid,mySub_qg_Lmu,mySub_qg_Lmu_func)
    
    xPij%tab = 0
    xPij_Lmu%tab = 0

    xPij_2%tab = 0
    xPij_2_Lmu%tab = 0

    do iQ = 0, PDFs%nQ

       !-- :,iflav,iQ
       !xPij%tab(:,:,iQ) = PDFs%tab(:,:,iQ)

       do iflav = -5,5
          if (iflav == 0) cycle

          if (ch.eq.'ns' .and. sec.eq.'fcsub12_qqb') then
             xPij%tab(:,iflav,iQ) = mySub_ns_qqb .conv. PDFs%tab(:,iflav,iQ)
             xPij_Lmu%tab(:,iflav,iQ) = mySub_ns_qqb_Lmu .conv. PDFs%tab(:,iflav,iQ)
          endif

          if (ch.eq.'ns' .and. sec.eq.'fcsub12_qq') then
             xPij%tab(:,iflav,iQ) = mySub_ns_qq .conv. PDFs%tab(:,-iflav,iQ)
             xPij_Lmu%tab(:,iflav,iQ) = mySub_ns_qq_Lmu .conv. PDFs%tab(:,-iflav,iQ)
          endif

          if (ch.eq.'ag' .and. sec.eq.'fcsub_onloqcd') then
             !-- Pqa \otimes a
             xPij%tab(:,iflav,iQ) = mySub_qg .conv. PDFs%tab(:,6,iQ)
             xPij_Lmu%tab(:,iflav,iQ) = mySub_qg_Lmu .conv. PDFs%tab(:,6,iQ)
          endif
          
       enddo

    enddo

  end subroutine init_xPij_nnlo_fc

  !--
  
  function mySubNNLO_ns_qqb_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: omx,lx,pqq,lomx,lopx
    real(dp) :: Li2mx,Li2omx,Li3x,Li3mx,Li3omx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

       omx = one-x
       lx = log(x)
       pqq = (one+x**2)/omx
       lomx = log(omx)
       lopx = log(one+x)

       Li2mx  = real(dilog2(-x),kind=dp)
       Li2omx = real(dilog2(omx),kind=dp)

       Li3x   = trilog(x)
       Li3mx  = trilog(-x)
       Li3omx = trilog(omx)

       !-- reg
       res = (1 + 7*x + 8*Li2mx*(1 + x) + (2*pisq*(1 + x))/3._dp + &
            (2*Li2omx*(-13 + 6*x + x**2))/omx - (lx**2*(5 - 12*x + 4*x**2))/(2._dp*omx) + &
            lx*(8*lopx*(1 + x) - (6 + 11*x - 27*x**2)/omx) + &
            lomx*(4*(-8 + 7*x) - (lx*(20 - 8*x**2))/omx) + &
            pqq*(16*Li3mx + 12*Li3omx + 18*Li3x - lx**3/3._dp + lomx*(-8*Li2omx + 5*lx**2) + &
            lx*(-8*Li2mx + 10*Li2omx - (7*pisq)/3._dp) - 6*zeta3))
       !print *, 'x',x,res
       
       !-- plus

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySubNNLO_ns_qqb_func

  function mySubNNLO_ns_qqb_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: omx,pqq,lx,Li2omx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

       omx = one-x
       pqq = (one+x**2)/omx
       
       lx = log(x)
       Li2omx = real(dilog2(omx),kind=dp)

       res = (16 + (4*Li2omx + 2*lx**2)*pqq - 14*x - (2*lx*(-5 + 2*x**2))/omx)
       !print *, 'res,Lmu2',x,res

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySubNNLO_ns_qqb_Lmu_func

  !--

  function mySubNNLO_ns_qq_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lx,pqqmx,lomx,lopx
    real(dp) :: Li2mx,Li2omx,Li3x,Li3mx,Li3omx,Li3xoopx,Li3omx2
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

       lx = log(x)
       pqqmx = (one+x**2)/(one+x)
       lomx = log(one-x)
       lopx = log(one+x)
       !
       Li2mx  = real(dilog2(-x),kind=dp)
       Li2omx = real(dilog2(one-x),kind=dp)
       !
       Li3x   = trilog(x)
       Li3mx  = trilog(-x)
       Li3omx = trilog(one-x)
       Li3xoopx = trilog(x/(one+x))
       Li3omx2  = trilog(one-x**2)

       !-- reg
       res =  4*lx**2 + pisq*(-1 - x) - 15*(-1 + x) - 12*Li2mx*(1 + x) + 4*Li2omx*(3 + x) + &
            lx*(11 + 19*x - 12*lopx*(1 + x)) + &
            lomx*(-16*(-1 + x) + 8*lx*(1 + x)) + &
            pqqmx*(-36*Li3mx + 16*Li3omx - 8*Li3omx2 - 16*Li3x - 24*Li3xoopx + 4*lopx**3 + &
            6*lopx*lx**2 + lx**3/3._dp + &
            lomx*(-16*Li2mx - 16*lopx*lx - 4*lx**2 - (4*pisq)/3._dp) - 2*lopx*pisq + &
            lx*(12*Li2mx - 4*Li2omx - 12*lopx**2 + (8*pisq)/3._dp) + 10*zeta3)
       !print *, 'x',x,res
       
       !-- plus

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySubNNLO_ns_qq_func

  function mySubNNLO_ns_qq_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: pqqmx,lx,lopx,Li2mx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

       pqqmx = (one+x**2)/(one+x)
       
       lx = log(x)
       lopx = log(one+x)
       Li2mx = real(dilog2(-x),kind=dp)

       res = ((8*Li2mx + 8*lopx*lx - 2*lx**2 + (2*pisq)/3._dp)*pqqmx + 8*(-1 + x) - &
            4*lx*(1 + x))
       !print *, 'res,Lmu2',x,res

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySubNNLO_ns_qq_Lmu_func

  !-- copied from mod_hoppet_nlo

    function mySub_qg_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: omx,lx,lomx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

       omx = one-x
       lx = log(x)
       lomx = log(omx)

       !-- reg
       res = ((x**2+omx**2)*(two*lomx-lx)+two*x*omx)

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_qg_func

  !-- -Pqg
  function mySub_qg_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: omx,lx,lomx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

       omx = one-x
       lx = log(x)
       lomx = log(omx)

       !-- reg
       res = -(x**2+omx**2)

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_qg_Lmu_func

end module mod_hoppet_nnlo_fc

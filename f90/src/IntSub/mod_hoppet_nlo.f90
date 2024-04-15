module mod_hoppet_nlo
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF, except => Cf, except => Ca, except => tr
  use mod_types, except => dp
  use consts_dp
  use mod_parms
  use mod_proc_parms
  use convolution
  use convolution_communicator
  use mod_hoppet_tools
  implicit none
  private
  !-- res = xPij + xPij_Lmu * log(musq/mvsq)
  !-- could put Lmu term inside xPij, but would not work for dynamic mvsq
  type(pdf_table), public, save :: xPij,xPij_Lmu
  type(pdf_table), public, save :: xPij_2,xPij_2_Lmu

  public :: init_xPij_nlo
  
contains

  subroutine init_xPij_nlo()
    integer :: iQ,iflav
    type(grid_conv) :: mySub_ns,mySub_ns_Lmu
    type(grid_conv) :: mySub_qg,mySub_qg_Lmu
    type(grid_conv) :: mySub_gq,mySub_gq_Lmu

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

    call AllocGridConv(grid,mySub_ns)
    call AllocGridConv(grid,mySub_ns_Lmu)
    !
    call AllocGridConv(grid,mySub_qg)
    call AllocGridConv(grid,mySub_qg_Lmu)
    !
    call AllocGridConv(grid,mySub_gq)
    call AllocGridConv(grid,mySub_gq_Lmu)

    call InitGridConv(grid,mySub_ns,mySub_ns_func)
    call InitGridConv(grid,mySub_ns_Lmu,mySub_ns_Lmu_func)
    !
    call InitGridConv(grid,mySub_qg,mySub_qg_func)
    call InitGridConv(grid,mySub_qg_Lmu,mySub_qg_Lmu_func)
    !
    call InitGridConv(grid,mySub_gq,mySub_gq_func)
    call InitGridConv(grid,mySub_gq_Lmu,mySub_gq_Lmu_func)

    xPij%tab = 0
    xPij_Lmu%tab = 0

    xPij_2%tab = 0
    xPij_2_Lmu%tab = 0

    do iQ = 0, PDFs%nQ

       !-- :,iflav,iQ
       !xPij%tab(:,:,iQ) = PDFs%tab(:,:,iQ)

       do iflav = -5,5
          if (iflav == 0) cycle
          if (ch.eq.'ns') then
             xPij%tab(:,iflav,iQ) = mySub_ns .conv. PDFs%tab(:,iflav,iQ)
             xPij_Lmu%tab(:,iflav,iQ) = mySub_ns_Lmu .conv. PDFs%tab(:,iflav,iQ)
          elseif (ch.eq.'gq' .or. ch.eq.'qg') then
             xPij%tab(:,iflav,iQ) = mySub_qg .conv. PDFs%tab(:,0,iQ)
             xPij_Lmu%tab(:,iflav,iQ) = mySub_qg_Lmu .conv. PDFs%tab(:,0,iQ)
          elseif (ch.eq.'aq' .or. ch.eq.'qa') then
             !-- Pqa \otimes a
             xPij%tab(:,iflav,iQ) = mySub_qg .conv. PDFs%tab(:,6,iQ)
             xPij_Lmu%tab(:,iflav,iQ) = mySub_qg_Lmu .conv. PDFs%tab(:,6,iQ)
          endif
       enddo

       !-- Paq \otimes q
       if (ch.eq.'aq' .or. ch.eq.'qa') then
          do iflav = -5,5,2
             xPij_2%tab(:,6,iQ) = xPij_2%tab(:,6,iQ) + Qdn2 * (mySub_gq .conv. PDFs%tab(:,iflav,iQ))
             xPij_2_Lmu%tab(:,6,iQ) = xPij_2_Lmu%tab(:,6,iQ) + Qdn2 * (mySub_gq_Lmu .conv. PDFs%tab(:,iflav,iQ))
          enddo
          do iflav = -4,4,2
             if (iflav==0) cycle
             xPij_2%tab(:,6,iQ) = xPij_2%tab(:,6,iQ) + Qup2 * (mySub_gq .conv. PDFs%tab(:,iflav,iQ))
             xPij_2_Lmu%tab(:,6,iQ) = xPij_2_Lmu%tab(:,6,iQ) + Qup2 * (mySub_gq_Lmu .conv. PDFs%tab(:,iflav,iQ))
          enddo
       endif

       !-- to add scale variation term
       !xPij_Lmu%tab(:,:,iQ) = PDFs%tab(:,:,iQ)*log(mz/PDFs%Q_vals(iQ))

    enddo

  end subroutine init_xPij_nlo

  function mySub_ns_func(y) result(res)
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

       res = (2*zeta2) !-- half -> one per leg
       
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_ns_func

  !--

  !-- -Pqq
  function mySub_ns_Lmu_func(y) result(res)
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

       res = -(one+x**2)/omx

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
       omx = one-x

       res = res + two/omx

    case(cc_DELTA) 

       res = -three/two
       
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_ns_Lmu_func

  !-- 
  
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

  !--
  
  function mySub_gq_func(y) result(res)
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
       res = ((one+omx**2)/x)*(two*lomx-lx)+x

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_gq_func

  !-- -Pgq
  function mySub_gq_Lmu_func(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: omx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)

       omx = one-x

       !-- reg
       res = -(one+omx**2)/x

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       
    case(cc_DELTA) 
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function mySub_gq_Lmu_func
  
end module mod_hoppet_nlo

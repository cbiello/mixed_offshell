module mod_lumi
  use mod_types
  use mod_consts_dp
  implicit none
  private

  public :: get_pdf_qed

contains

  subroutine get_pdf_qed(z1,z2,muf,f1,f2)
    real(dp), intent(in)  :: z1,z2,muf
    real(dp), intent(out) :: f1(:),f2(:)
    real(dp15) :: x1,x2,muf_dp15,f1_dp15(-6:7),f2_dp15(-6:7)

    f1 = zero
    f2 = zero

    x1 = real(z1,kind=dp15)
    x2 = real(z2,kind=dp15)

    muf_dp15 = real(muf,kind=dp15)

    call evolvePDFphoton(x1,muf_dp15,f1_dp15,f1_dp15(7))
    f1 = f1_dp15/x1

    call evolvePDFphoton(x2,muf_dp15,f2_dp15,f2_dp15(7))
    f2 = f2_dp15/x2

    return

  end subroutine get_pdf_qed

  ! !-- use this one if you don't need QED PDFs
  ! subroutine get_pdf(z1,z2,muf,f1,f2)
  !   real(dp), intent(in)  :: z1,z2,muf
  !   real(dp), intent(out) :: f1(:),f2(:)
  !   real(dp15) :: x1,x2,muf_dp15,f1_dp15(-6:6),f2_dp15(-6:6)

  !   f1 = zero
  !   f2 = zero

  !   x1 = real(z1,kind=dp15)
  !   x2 = real(z2,kind=dp15)

  !   muf_dp15 = real(muf,kind=dp15)

  !   call evolvePDF(x1,muf_dp15,f1_dp15)
  !   f1 = f1_dp15/x1

  !   call evolvePDF(x2,muf_dp15,f2_dp15)
  !   f2 = f2_dp15/x2

  !   return

  ! end subroutine get_pdf
  
end module mod_lumi


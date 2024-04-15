module mod_eikonals
  use mod_types
  use mod_consts_dp
  implicit none

  public :: get_qcd_eik, get_qed_eik

  private

contains

  subroutine get_qcd_eik(casimir,etas,array_pos,g_pos,eik)
    real(dp), intent(in)  :: casimir, etas(:,:)
    real(dp), intent(out) :: eik
    integer, intent(in) :: array_pos(:), g_pos
    integer  :: i, j, m, n

    eik = zero
    do i = 1, size(array_pos)
      do j = i+1, size(array_pos)
        m = array_pos(i)
        n = array_pos(j)
        eik = eik + &
              casimir*etas(m,n)/etas(m,g_pos)/etas(n,g_pos)
      enddo
    enddo

  end subroutine get_qcd_eik

  subroutine get_qed_eik(charges,etas,array_pos,a_pos,eik)
    real(dp), intent(in)  :: charges(:,:)
    real(dp), intent(in)  :: etas(:,:)
    real(dp), intent(out) :: eik(4)
    integer,  intent(in)  :: array_pos(:), a_pos
    integer  :: i, j, m, n

    eik = zero
    do i = 1, size(array_pos)
      do j = i+1, size(array_pos)
        m = array_pos(i)
        n = array_pos(j)
        eik = eik - &
              charges(:,i)*charges(:,j)*etas(m,n)/etas(m,a_pos)/etas(n,a_pos)
      enddo
    enddo

  end subroutine get_qed_eik

end module mod_eikonals

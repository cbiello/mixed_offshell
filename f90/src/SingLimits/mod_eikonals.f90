module mod_eikonals
  use mod_types
  use mod_consts_dp
  use mod_parms
  implicit none

  public :: get_qcd_eik, get_qed_eik, get_qed_eik_gen

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



  subroutine get_qed_eik_gen(ampl,etas,array_pos,a_pos,res_out)
    real(dp), intent(in)  :: ampl(-5:7,-5:7)
    real(dp), intent(in)  :: etas(:,:)
    integer,  intent(in)  :: array_pos(:), a_pos
    real(dp), intent(out) :: res_out(-5:7,-5:7)
    real(dp)              :: eik
    real(dp)              :: Qvec(5), Qq(-5:5)
    integer  :: i, j, m, n, al, be

    if (size(array_pos) .ne. 4) then
       print *, "incorrect size of array in QED eikonal"
       stop
    endif



    eik = zero
    res_out = zero
    do i = 1, size(array_pos)
       do j = i+1, size(array_pos)
          m = array_pos(i)
          n = array_pos(j)
          eik =  -etas(m,n)/etas(m,a_pos)/etas(n,a_pos)
          do al = -5,5
             do be = -5,5
                Qvec = [ -Q_IS(al), -Q_IS(be), Q3, Q4, - Qq(al)  - Qq(be) - Q3 - Q4]                
                res_out(al,be) = res_out(al,be) + ampl(al,be)*eik*Qvec(m)*Qvec(n)

             enddo
          enddo
          
      enddo
    enddo

  end subroutine get_qed_eik_gen

end module mod_eikonals

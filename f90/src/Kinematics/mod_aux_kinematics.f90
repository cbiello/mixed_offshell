module mod_aux_kinematics
  use mod_types
  use mod_consts_dp
  use mod_proc_parms
  use mod_process
  implicit none

  public :: etaij, fill_etas
  public :: NF1
  public :: get_part_s_xi_nlo,get_part_s_xi_nnlo
  public :: get_part_s_xi_nlo_z

  private

contains

  function NF1(x3,x4,lambda)
    implicit none
    real(dp), intent(in) :: x3,x4,lambda
    real(dp) :: NF1

    if (x4.eq.zero) then 
       NF1 = one 
       return 
    elseif(x4.eq.one) then 
       NF1 = four*lambda*(one-x3) 
       return 
    else
       NF1 = one + x4*(one-two*x3) & 
            - two*(one-two*lambda)*sqrt(abs(x4*(one-x3)*(one-x3*x4)))
    endif

  end function NF1

  function etaij(ni,nj)
    real(dp), intent(in) :: ni(4), nj(4)
    real(dp) :: etaij
    etaij = half*(one-dot_product(ni(2:4),nj(2:4)))
  end function etaij

  !-- Consider if a more efficient version if possible
  !-- for higher multiplicities, e.g. using masks
  subroutine fill_etas(proc,ns,eta_pos)
    type(KinConfig), intent(inout) :: proc
    real(dp), intent(in) :: ns(:,:)
    integer, intent(in) :: eta_pos(:)
    integer :: i,m,n

    do i = 1, size(eta_pos)
      n = eta_pos(i)
      do m = 1, size(ns,2)
        if(m == n) cycle
        proc%Lim_etaij(m,n) = etaij(ns(:,m),ns(:,n))
        proc%Lim_etaij(n,m) = proc%Lim_etaij(m,n)
      enddo
    enddo

  end subroutine

  subroutine get_part_s_xi_nlo(mv2,y,x1,eta51,eta52,spart,sqrts,xi1,xi2,flag)
    real(dp), intent(in)  :: mv2,y,x1,eta51,eta52
    logical,  intent(out) :: flag
    real(dp), intent(out) :: sqrts,spart
    real(dp), intent(out) :: xi1,xi2
    real(dp) :: yy,Deltan,Deltad

    flag = .true.

    spart = mv2/(one - x1)
    if (spart.gt.sh) return

    sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta51)
    Deltad = one - x1*(one - eta52)

    if ((Deltan.gt.zero).and.(Deltad.gt.zero)) then

       yy = Y - half*log(Deltan/Deltad)
       xi1 = sqrt(spart/sh)*exp( yy)
       xi2 = sqrt(spart/sh)*exp(-yy)

       if (xi1.gt.one .or. xi2.gt.one) return

    else

       return

    endif

    flag = .false.

  end subroutine get_part_s_xi_nlo

  subroutine get_part_s_xi_nnlo(Q2,Y,x1,x2,eta65,eta51,eta61,eta52,eta62,spart,sqrts,xi1,xi2,flag)
    real(dp), intent(in)  :: Q2,Y,x1,x2,eta65,eta51,eta61,eta52,eta62
    real(dp), intent(out) :: spart,sqrts
    real(dp), intent(out) :: xi1,xi2
    logical,  intent(out) :: flag
    real(dp) :: yy,Deltas,Deltan,Deltad

    flag = .true.

    Deltas = one - (x1+x2) + x1*x2*eta65
    if (Deltas.lt.Q2/sh) return

    spart = Q2/Deltas
    sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta51) - x2*(one - eta61)
    Deltad = one - x1*(one - eta52) - x2*(one - eta62)

    if ((Deltan.gt.zero) .and. (Deltad.gt.zero)) then 

       yy = Y - half*log(Deltan/Deltad)
       xi1 = sqrt(spart/sh)*exp( yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)

       if (xi1.gt.one .or. xi2.gt.one) return

    else

       return

    endif

    flag = .false.

  end subroutine get_part_s_xi_nnlo

  !-- To be used when p1/p2 are boosted
  subroutine get_part_s_xi_nlo_z(mv2,y,x1,z1,z2,eta51,eta52,spart,sqrts,xi1,xi2,flag)
    real(dp), intent(in)  :: mv2,y,x1,z1,z2,eta51,eta52
    logical,  intent(out) :: flag
    real(dp), intent(out) :: sqrts,spart
    real(dp), intent(out) :: xi1,xi2
    real(dp) :: yy,Deltan,Deltad,Deltas

    flag = .true.

    Deltas = z1*z2 - z1*x1*eta51 - z2*x1*eta52

    if (Deltas < mv2/sh) return

    spart = mv2/Deltas
    sqrts = sqrt(spart)

    Deltan = z1 - x1*eta52
    Deltad = z2 - x1*eta51

    if ((Deltan.gt.zero).and.(Deltad.gt.zero)) then

       yy = Y - half*log(Deltan/Deltad)
       xi1 = sqrt(spart/sh)*exp( yy)
       xi2 = sqrt(spart/sh)*exp(-yy)

       if (xi1.gt.one .or. xi2.gt.one) return

    else

       return

    endif

    flag = .false.

  end subroutine get_part_s_xi_nlo_z

end module mod_aux_kinematics

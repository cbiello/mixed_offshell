module mod_aux_kinematicsNNLO_EunordAord_ga
  use mod_types
  use mod_consts_dp
  use mod_proc_parms
  implicit none
  private

  !-- interf not needed, identical to nlo
  !-- sector_nnlo_minimal -> partition [fact]
  public :: get_flag_nnlo, get_flag_nlo, get_flag_nnlo_interf
  public :: sector_nnlo_minimal

  !-- get flag: get_part_xi --> 1 for NLO [not really needed for nlo but used if x1 or x2 = 0], 1 for NNLO
  !-- add NF1 here

  
  
contains

  function NF1(x3,x4,x5) 
    implicit none  
    real(dp), intent(in) :: x3,x4,x5
    real(dp) :: NF1
    
    if (x4.eq.zero) then 
       NF1 = one 
       return 
    elseif(x4.eq.one) then 
       NF1 = four*x5*(one-x3) 
       return 
    else
       NF1 = one + x4*(one-two*x3) & 
            - two*(one-two*x5)*sqrt(abs(x4*(one-x3)*(one-x3*x4)))
    endif
    
  end function NF1

  subroutine get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,flag)
    real(dp), intent(in)  :: MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41
    logical,  intent(out) :: flag
    real(dp), intent(out) :: spart,sqrts
    real(dp), intent(out) :: xi1,xi2
    real(dp) :: yy,Deltas,Deltan,Deltad

    Deltas = one - (x1+x2) + x1*x2*eta54

    if (Deltas.lt.MV2/sh) flag = .true.
    if (flag .eqv. .false.) spart = MV2/Deltas
    if (flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) - x2*(one - eta51)
    Deltad = one - x1*(one - eta42) - x2*(one - eta52)

    if ( (flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad)

       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)

       if (xi1.gt.one .or. xi2.gt.one) flag = .true.
    else
       flag = .true.
    endif

  end subroutine get_flag_nnlo

  subroutine get_flag_nlo(xi,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,flag)
    real(dp), intent(in)  :: xi,MV2,Y,x2,eta51,eta52
    logical,  intent(out) :: flag
    real(dp), intent(out) :: sqrts,spart
    real(dp), intent(out) :: xi1,xi2
    real(dp) :: yy,Deltan,Deltad

    spart = MV2/(one - xi)      
    if (spart.gt.sh) flag = .true.  
    sqrts = sqrt(spart)

    Deltan = one - xi*(one - eta51)
    Deltad = one - xi*(one - eta52)

    if ( (flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)

       if (xi1.gt.one .or. xi2.gt.one) flag = .true.
    else 
       flag = .true.
    endif

  end subroutine get_flag_nlo

  subroutine get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,flag)
    real(dp), intent(in)  :: sh,MV2,Y,x1,eta41,eta42
    real(dp), intent(out) :: xi1,xi2,sqrts,spart
    logical, intent(out)  :: flag
    real(dp) :: yy,Deltas,Deltan,Deltad

    Deltas = one - x1 
    if (Deltas.lt.MV2/sh) flag = .true.
    if (flag .eqv. .false.) spart = MV2/Deltas
    if (flag.eqv. .false.) sqrts = sqrt(spart)

    Deltan = one - x1*(one - eta41) 
    Deltad = one - x1*(one - eta42) 

    if ( (flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then 
       yy = Y - half*log(Deltan/Deltad) 
       xi1 = sqrt(spart/sh)*exp(yy) 
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) flag = .true.
    else
       flag = .true.
    endif

    return
  end subroutine get_flag_nnlo_interf

  subroutine sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,i,j,weight)
    real(dp), intent(in) :: eta41,eta42,eta51,eta52,eta56,eta57
    integer, intent(in) :: i,j
    real(dp), intent(out) :: weight
    real(dp) :: denom,wpart

    denom = (eta51+eta52)*eta56*eta57+(eta56+eta57)*eta51*eta52

    if (j.eq.1) wpart = eta52*eta56*eta57
    if (j.eq.2) wpart = eta51*eta56*eta57
    if (j.eq.6) wpart = eta51*eta52*eta57
    if (j.eq.7) wpart = eta51*eta52*eta56


    if (i.eq.1) weight = eta42*wpart/denom
    if (i.eq.2) weight = eta41*wpart/denom

    if (i.eq.4) weight = wpart/denom

    !print *, 'switching off damping'
    !weight = 1
    
    return

  end subroutine sector_nnlo_minimal

end module mod_aux_kinematicsNNLO_EunordAord_ga


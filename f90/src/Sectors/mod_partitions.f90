module mod_partitions
  use mod_types
  use mod_consts_dp
  use mod_proc_parms
  use mod_process
  implicit none

  public :: partition_nlo_qed,partition_nloAA_qed
  public :: partition_nnlo_qcd,partition_nnlo_fact

  private

contains

  !-- split into initial state (i=12), and two final states (i=3,4)
  !-- partitions: 1/(eta51*eta52),1/eta53,1/eta54
  subroutine partition_nlo_qed(proc,damp,i)
    type(KinConfig), intent(in) :: proc
    integer, intent(in)         :: i
    real(dp), intent(out)       :: damp
    real(dp) :: eta5i(4)
    real(dp) :: den

    eta5i(:) = proc%Lim_etaij(1:4,5)
    
    den = eta5i(1)*eta5i(2)*eta5i(3) + eta5i(1)*eta5i(2)*eta5i(4) + eta5i(3)*eta5i(4)

    if (i == 12) then
       damp = eta5i(3)*eta5i(4)
    elseif (i==3) then
       damp = eta5i(1)*eta5i(2)*eta5i(4)
    elseif (i==4) then
       damp = eta5i(1)*eta5i(2)*eta5i(3)
    endif

    damp = damp/den

  end subroutine partition_nlo_qed

  !-- split into two final states (i=3,4)
  !-- partitions: 1/(eta51*eta52),1/eta53,1/eta54
  subroutine partition_nloAA_qed(proc,damp,i)
    type(KinConfig), intent(in) :: proc
    integer, intent(in)         :: i
    real(dp), intent(out)       :: damp
    real(dp) :: eta5i(4)
    real(dp) :: den

    eta5i(:) = proc%Lim_etaij(1:4,5)
    
    den = eta5i(3) + eta5i(4)

    if (i==3) then
       damp = eta5i(4)
    elseif (i==4) then
       damp = eta5i(3)
    endif

    damp = damp/den

  end subroutine partition_nloAA_qed

  
  !! Written for channel q(1) qb(2) -> l-(3) l+(4) g(5) a(6)
  !! eta_{15} -> eta_gq , eta_{25} -> eta_gqb
  !! eta_{16} -> eta_aq , eta_{26} -> eta_aqb
  !! eta_{36} -> eta_al-, eta_{46} -> eta_al+
  !! iconf_qcd/qed are the crossings from the process above
  !! i_qcd is the position in iconf_qcd of the particle emitting the gluon
  !! i_qed is the position in iconf_qed of the particle emitting the photon
  subroutine partition_nnlo_fact(proc,part,iconf_qcd,i_qcd,iconf_qed,i_qed)
    type(KinConfig), intent(in) :: proc
    real(dp), intent(out) :: part
    integer, optional, intent(in) :: iconf_qcd(:), i_qcd
    integer, optional, intent(in) :: iconf_qed(:), i_qed
    real(dp) :: eta_51, eta_52
    real(dp) :: eta_61, eta_62, eta_63, eta_64
    real(dp) :: denom, part_ew, part_qcd

    part_ew  = one
    part_qcd = one

    if(present(i_qed)) then
      eta_61 = proc%Lim_etaij(iconf_qed(1),iconf_qed(5))
      eta_62 = proc%Lim_etaij(iconf_qed(2),iconf_qed(5))
      eta_63 = proc%Lim_etaij(iconf_qed(3),iconf_qed(5))
      eta_64 = proc%Lim_etaij(iconf_qed(4),iconf_qed(5))

      denom = (eta_61+eta_62)*eta_63*eta_64+(eta_63+eta_64)*eta_61*eta_62
      if (i_qed .eq. 1) part_ew = eta_62 * eta_63 * eta_64 / denom
      if (i_qed .eq. 2) part_ew = eta_61 * eta_63 * eta_64 / denom
      if (i_qed .eq. 3) part_ew = eta_61 * eta_62 * eta_64 / denom
      if (i_qed .eq. 4) part_ew = eta_61 * eta_62 * eta_63 / denom
    endif

    if(present(i_qcd)) then
      eta_51 = proc%Lim_etaij(iconf_qcd(1),iconf_qcd(3))
      eta_52 = proc%Lim_etaij(iconf_qcd(2),iconf_qcd(3))
      denom = eta_51 + eta_52
      if (i_qcd .eq. 1) part_qcd = eta_52/denom
      if (i_qcd .eq. 2) part_qcd = eta_51/denom
    endif

    part = part_ew * part_qcd

  end subroutine partition_nnlo_fact

  !! Partition function defined in (3.12) of 1902.02081
  !! iconf(:,1) = positions of initial-state particles
  !! iconf(:,2) = positions final-state (radiated) particles
  subroutine partition_nnlo_qcd(proc, part, iconf, i, j)
    type(KinConfig), intent(in) :: proc
    real(dp), intent(out) :: part
    integer, intent(in)   :: i, j, iconf(:,:)
    real(dp) :: eta51, eta52, eta61, eta62, eta65, denom1, denom2

    eta51 = proc%Lim_etaij(iconf(1,1),iconf(1,2))
    eta52 = proc%Lim_etaij(iconf(2,1),iconf(1,2))
    eta61 = proc%Lim_etaij(iconf(1,1),iconf(2,2))
    eta62 = proc%Lim_etaij(iconf(2,1),iconf(2,2))
    eta65 = proc%Lim_etaij(iconf(1,2),iconf(2,2))

    denom1 = eta65 + eta52 + eta61
    denom2 = eta65 + eta51 + eta62

    !! w_5161
    if((i.eq.1) .and. (j.eq.1)) then
      part = eta52*eta62*(one + eta51/denom1 + eta61/denom2)

    !! w_5262
    elseif((i.eq.2) .and. (j.eq.2)) then
      part = eta51*eta61*(one + eta62/denom1 + eta52/denom2)

    !! w_5162
    elseif((i.eq.1) .and. (j.eq.2)) then
      part = eta52*eta61*eta65/denom2

    !! w_5261
    elseif((i.eq.2) .and. (j.eq.1)) then
      part = eta51*eta62*eta65/denom1

    else
      print *,'In mod_partitions/partition_nnlo_qcd'
      print *,'Calling partition_nnlo_qcd with illegal values of (i,j). Aborting'
      stop
    endif

  end subroutine partition_nnlo_qcd

end module mod_partitions

!-- Double-collinear phase space, ISR.
!-- Damping factor not included
module mod_kinematics_nnlo_dc_ii
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_process
  use mod_auxfunctions
  use mod_kinematics_gen
  use mod_kinematics_lept
  use mod_aux_kinematics
  use mod_proc_parms
  implicit none

  integer, public, parameter :: NNLOdim_vegas = tauy_max + 2 + 3 + 3
  !
  integer, public, parameter :: kNNLO_min = tauy_max + 1
  integer, public, parameter :: kNNLO_max = tauy_max + 2 + 3 + 3
  integer, public, parameter :: kNNLO_max_full = kNNLO_max + 2
  !
  integer, public, parameter :: xE5 = kNNLO_min
  integer, public, parameter :: xE6 = kNNLO_min + 1
  integer, public, parameter :: xRHO5 = kNNLO_min + 2
  integer, public, parameter :: xRHO6 = kNNLO_min + 4
  !!
  integer, parameter :: xPHI5   = kNNLO_min + 3
  integer, parameter :: xPHI6   = kNNLO_min + 5
  !
  integer, parameter :: xLMIN  = kNNLO_min + 6
  integer, parameter :: xLMAX  = kNNLO_min + 8
  !
  integer, parameter :: npart_lo = 4, npart_nlo = 5, npart_nnlo = 6

  public :: kinematics_nnlo_ii_5i6j, kinematics_nnlo_ii_5i6i

  private

contains

  !!=========================================================================!!
  !! Kinematics for IS double-collinear configuration                        !!
  !! {i,j} = {1,2} or {i,j} = {2,1}                                          !!
  !!                                                                         !!
  !! opt_etas: indices k for eta_kl, l=1,...,n.                              !!
  !! eta_ij, where i,j refer to initial-state particles and                  !!
  !! radiated particles, are always computed by default                      !!
  !!=========================================================================!!
  subroutine kinematics_nnlo_ii_5i6j(yr, i, j, &
    HardProc, S5Lim, C5Lim, C5S5Lim, S6Lim, C6Lim, C6S6Lim,  &
    S5S6Lim, C5C6Lim, C5S6Lim, C6S5Lim,     &
    C5S5S6Lim, C6S5S6Lim, C5C6S5Lim, C5C6S6Lim, C5C6S5S6Lim, &
    opt_etas)
    real(dp),        intent(in)    :: yr(:)
    integer,         intent(in)    :: i, j
    type(KinConfig), intent(inout) :: HardProc
    type(KinConfig), optional, intent(inout) :: S5Lim, C5Lim, S6Lim, C6Lim
    type(KinConfig), optional, intent(inout) :: S5S6Lim, C5S5Lim, C6S6Lim
    type(KinConfig), optional, intent(inout) :: C5C6Lim, C5S6Lim, C6S5Lim
    type(KinConfig), optional, intent(inout) :: C5S5S6Lim, C6S5S6Lim
    type(KinConfig), optional, intent(inout) :: C5C6S5Lim, C5C6S6Lim
    type(KinConfig), optional, intent(inout) :: C5C6S5S6Lim
    integer, optional, intent(in) :: opt_etas(:)
    !---
    real(dp) :: mv, mv2, tau, ylab, jac
    real(dp) :: p1(4), p2(4), p5(4), p6(4), pv(4)
    real(dp) :: sgn, kallenF
    real(dp) :: x1, x2, x3, x4, x5, x6
    real(dp) :: cos5i, sin5i, cos6j, sin6j
    real(dp) :: phi5, phi6, cos_phi5, cos_phi6, sin_phi5, sin_phi6
    real(dp) :: xi1, xi2, spart, sqrts
    real(dp) :: nlept(3,2), ns(4,6)
    real(dp) :: eta5i, eta6i, eta5j, eta6j
    real(dp) :: eta51, eta61, eta52, eta62, eta65

    !-- process-dependent part
    call get_tauy(yr(1:tauy_max),tau,ylab,jac)
    mv2 = tau*sh
    mv  = sqrt(mv2)

    !-- Initialization of kinematic configurations
    call initialize_config(HardProc)                             !  1
    if(present(S5Lim))       call initialize_config(S5Lim)       !  2
    if(present(C5S5Lim))     call initialize_config(C5S5Lim)     !  3
    if(present(S6Lim))       call initialize_config(S6Lim)       !  4
    if(present(C6S6Lim))     call initialize_config(C6S6Lim)     !  5
    if(present(C5Lim))       call initialize_config(C5Lim)       !  6
    if(present(C6Lim))       call initialize_config(C6Lim)       !  7
    if(present(S5S6Lim))     call initialize_config(S5S6Lim)     !  8
    if(present(C5S5S6Lim))   call initialize_config(C5S5S6Lim)   !  9
    if(present(C6S5S6Lim))   call initialize_config(C6S5S6Lim)   ! 10
    if(present(C5S6Lim))     call initialize_config(C5S6Lim)     ! 11
    if(present(C5C6S6Lim))   call initialize_config(C5C6S6Lim)   ! 12
    if(present(C6S5Lim))     call initialize_config(C6S5Lim)     ! 13
    if(present(C5C6S5Lim))   call initialize_config(C5C6S5Lim)   ! 14
    if(present(C5C6Lim))     call initialize_config(C5C6Lim)     ! 15
    if(present(C5C6S5S6Lim)) call initialize_config(C5C6S5S6Lim) ! 16

    x1 = yr(xE5)   ! E_5
    x2 = yr(xE6)   ! E_6
    x3 = yr(xRHO5) ! \theta_{5i}
    x4 = yr(xRHO6) ! \theta_{6j}
    x5 = yr(xPHI5) ! \phi_{5}
    x6 = yr(xPHI6) ! related to \lambda

    !! Choose sign of \cos\theta depending on emitter 
    !! {i,j} = {1,2} or {i,j} = {2,1}
    if(i .eq. 1) then
      sgn = +one
    elseif(i .eq. 2) then
      sgn = -one
    else
      print *, 'wrong emitter'
      sgn = zero
      stop
    endif

    !!-----------------------------------------------------------------------!!
    !!                      STRUCTURE 1: Hard Process.                       !!
    !!-----------------------------------------------------------------------!!

    !! Azimuthal variables
    phi5 = twopi*x5 
    phi6 = twopi*x6
    cos_phi5 = cos(phi5)
    sin_phi5 = sin(phi5)
    cos_phi6 = cos(phi6)
    sin_phi6 = sin(phi6)

    !! Polar variables
    eta5i = x3
    eta5j = one - x3
    eta6i = one - x4
    eta6j = x4

    cos5i = sgn*( one - 2*x3)
    sin5i = two*sqrt(x3*(one-x3))

    cos6j = sgn*(-one + 2*x4)
    sin6j = two*sqrt(x4*(one-x4))

    eta65 = half*(one - sin5i*sin6j*cos(phi5-phi6) - cos5i*cos6j)

    !! set eta variables in the limit
    HardProc%Lim_etaij(1,2) = one
    HardProc%Lim_etaij(i,5) = eta5i
    HardProc%Lim_etaij(j,5) = eta5j
    HardProc%Lim_etaij(i,6) = eta6i
    HardProc%Lim_etaij(j,6) = eta6j
    HardProc%Lim_etaij(5,6) = eta65

    eta51 = HardProc%Lim_etaij(1,5)
    eta52 = HardProc%Lim_etaij(2,5)
    eta61 = HardProc%Lim_etaij(1,6)
    eta62 = HardProc%Lim_etaij(2,6)

    ns(:,1) = [one,zero,zero, one]
    ns(:,2) = [one,zero,zero,-one]

    ns(:,5) = [one,sin5i*cos_phi5,sin5i*sin_phi5,cos5i]
    ns(:,6) = [one,sin6j*cos_phi6,sin6j*sin_phi6,cos6j]

    !! Compute partonic fractions {xi1, xi2}, partonic s and its square root
    !! and check if point is physical
    call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
      spart,sqrts,xi1,xi2,HardProc%flag)

    if (.not. HardProc%flag) then

      HardProc%npart = npart_nnlo
       
      HardProc%PartFrac = (/xi1,xi2/)

      p1 = half*sqrts*ns(:,1)
      p2 = half*sqrts*ns(:,2)
      p5 = half*sqrts*x1*ns(:,5)
      p6 = half*sqrts*x2*ns(:,6)

      pv = p1 + p2 - p5 - p6

      call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, HardProc)
      ns(:,3) = (/one,nlept(:,1)/)
      ns(:,4) = (/one,nlept(:,2)/)

      HardProc%AmpMom(:,1) = p1
      HardProc%AmpMom(:,2) = p2
      HardProc%AmpMom(:,5) = p5
      HardProc%AmpMom(:,6) = p6

      if(present(opt_etas)) call fill_etas(HardProc,ns,opt_etas)

      HardProc%wgt = one/8.0_dp/pi * kallenF &
           * (spart/two)**2 & 
           * x1*x2          & 
           * one/two/spart  &
           * spart/mv2      &
           * jac

    endif

    !!-----------------------------------------------------------------------!!
    !!                        STRUCTURE 2: S5 FLM.                           !!
    !!-----------------------------------------------------------------------!!

    if(present(S5Lim)) then

      !! Angles same as in hard configuration, thus not recomputed
      call get_part_s_xi_nlo(mv2,ylab,x2,eta61,eta62,spart,sqrts,xi1,xi2,S5Lim%flag)

      !! ------------------------------- S5 -------------------------------- !!

      if (.not. S5Lim%flag) then 

        S5Lim%npart = npart_nlo

        S5Lim%PartFrac = (/xi1,xi2/)

        ns(:,6) = [one,sin6j*cos_phi6,sin6j*sin_phi6,cos6j]

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p6 = half*sqrts*x2*ns(:,6)

        pv = p1 + p2 - p6
        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, S5Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        S5Lim%AmpMom(:,1) = p1
        S5Lim%AmpMom(:,2) = p2
        S5Lim%AmpMom(:,5) = p6

        !! E5
        S5Lim%Lim_KinInv(1) = x1*half*sqrts

        !! in the following: set eta_ij, ij = {21,51,52,61,62,65}
        S5Lim%Lim_etaij(1,2)   = one
        S5Lim%Lim_etaij(1:2,5) = HardProc%Lim_etaij(1:2,5)
        S5Lim%Lim_etaij(1:2,6) = HardProc%Lim_etaij(1:2,6)
        S5Lim%Lim_etaij(5,6)   = HardProc%Lim_etaij(5,6)
        S5Lim%Lim_etaij(6,5)   = HardProc%Lim_etaij(5,6)

        if(present(opt_etas)) call fill_etas(S5Lim,ns,opt_etas)

        S5Lim%wgt = one/8.0_dp/pi * kallenF &
           * (spart/two)**2 & 
           * x1*x2          & 
           * one/two/spart  &
           * spart/mv2      &
           * jac

    !! -------------------------------- C5j + S5 --------------------------- !!

        if(present(C5S5Lim)) then

          C5S5Lim%PartFrac = S5Lim%PartFrac

          ns(:,5) = [one,zero,zero,sgn]

          C5S5Lim%Lim_etaij = S5Lim%Lim_etaij
          C5S5Lim%Lim_etaij(i,5) = zero
          C5S5Lim%Lim_etaij(j,5) = one
          C5S5Lim%Lim_etaij(5,6) = one - x4

          C5S5Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)

          C5S5Lim%Lim_z(1) = x1

          if(present(opt_etas)) call fill_etas(C5S5Lim,ns,opt_etas)

          C5S5Lim%wgt =  one/8.0_dp/pi * kallenF &
             * (spart/two)**2 & 
             * x1*x2          & 
             * one/two/spart  & 
             * spart/mv2      &
             * jac

        endif
      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                        STRUCTURE 3: S6 FLM.                           !!
    !!-----------------------------------------------------------------------!!

    if(present(S6Lim)) then

      !! Angles same as in hard configuration, thus not recomputed
      call get_part_s_xi_nlo(mv2,ylab,x1,eta51,eta52,spart,sqrts,xi1,xi2,S6Lim%flag)

      !! ------------------------------- S6 -------------------------------- !!

      if (.not. S6Lim%flag) then 

        S6Lim%npart = npart_nlo
         
        S6Lim%PartFrac = (/xi1,xi2/)

        ns(:,5) = [one,sin5i*cos_phi5,sin5i*sin_phi5,cos5i]

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = half*sqrts*x1*ns(:,5)

        pv = p1 + p2 - p5
        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, S6Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        S6Lim%AmpMom(:,1) = p1
        S6Lim%AmpMom(:,2) = p2
        S6Lim%AmpMom(:,5) = p5

        !! E6
        S6Lim%Lim_KinInv(1) = x2*half*sqrts

        !! in the following: set eta_ij, ij = {21,51,52,61,62,65}
        S6Lim%Lim_etaij(1,2)   = one
        S6Lim%Lim_etaij(1:2,5) = HardProc%Lim_etaij(1:2,5)
        S6Lim%Lim_etaij(1:2,6) = HardProc%Lim_etaij(1:2,6)
        S6Lim%Lim_etaij(5,6)   = HardProc%Lim_etaij(5,6)
        S6Lim%Lim_etaij(6,5)   = HardProc%Lim_etaij(5,6)

        if(present(opt_etas)) call fill_etas(S6Lim,ns,opt_etas)

        S6Lim%wgt = one/8.0_dp/pi * kallenF &
           * (spart/two)**2 & 
           * x1*x2          & 
           * one/two/spart  &
           * spart/mv2      &
           * jac

      !! ------------------------------- C6j + S6 -------------------------- !!

        if(present(C6S6Lim)) then

          C6S6Lim%PartFrac = S6Lim%PartFrac

          ns(:,6) = [one,zero,zero,-sgn]

          C6S6Lim%Lim_etaij = S6Lim%Lim_etaij
          C6S6Lim%Lim_etaij(i,6) = one
          C6S6Lim%Lim_etaij(j,6) = zero
          C6S6Lim%Lim_etaij(5,6) = one - x3

          C6S6Lim%Lim_sij(j,6) = spart*x2*HardProc%Lim_etaij(j,6)

          C6S6Lim%Lim_z(2) = x2

          if(present(opt_etas)) call fill_etas(C6S6Lim,ns,opt_etas)

          C6S6Lim%wgt =  one/8.0_dp/pi * kallenF &
             * (spart/two)**2 & 
             * x1*x2          & 
             * one/two/spart  & 
             * spart/mv2      &
             * jac

        endif
      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                        STRUCTURE 4: C5i FLM.                          !!
    !!-----------------------------------------------------------------------!!

    if(present(C5Lim)) then

      !! x3 -> 0
      C5Lim%Lim_etaij(1,2) = one
      C5Lim%Lim_etaij(i,5) = zero
      C5Lim%Lim_etaij(j,5) = one
      C5Lim%Lim_etaij(i,6) = one - x4
      C5Lim%Lim_etaij(j,6) = x4
      C5Lim%Lim_etaij(5,6) = one - x4

      eta51 = C5Lim%Lim_etaij(1,5)
      eta52 = C5Lim%Lim_etaij(2,5)
      eta61 = C5Lim%Lim_etaij(1,6)
      eta62 = C5Lim%Lim_etaij(2,6)
      eta65 = one - x4

      call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
        spart,sqrts,xi1,xi2,C5Lim%flag)

      if (.not. C5Lim%flag) then 

        C5Lim%npart = npart_nlo         
         
        C5Lim%PartFrac = (/xi1,xi2/)

        ns(:,5) = [one,zero,zero,sgn]
        ns(:,6) = [one,sin6j*cos_phi6,sin6j*sin_phi6,cos6j]

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = half*sqrts*x1*ns(:,5)
        p6 = half*sqrts*x2*ns(:,6)

        if(i .eq. 1) then
          p1 = p1 - p5
        elseif(i .eq. 2) then
          p2 = p2 - p5
        endif

        pv = p1 + p2 - p6
        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C5Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        C5Lim%AmpMom(:,1) = p1
        C5Lim%AmpMom(:,2) = p2
        C5Lim%AmpMom(:,5) = p6

        C5Lim%Lim_z(1) = one/(one - x1)
        C5Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)

        if(present(opt_etas)) call fill_etas(C5Lim,ns,opt_etas)

        C5Lim%wgt = one/8.0_dp/pi * kallenF &
             * (spart/two)**2 & 
             * x1*x2          & 
             * one/two/spart  & 
             * spart/MV2      &
             * jac

      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                        STRUCTURE 5: C6j FLM.                          !!
    !!-----------------------------------------------------------------------!!

    if(present(C6Lim)) then

      !! x4 -> 0
      C6Lim%Lim_etaij(1,2) = one
      C6Lim%Lim_etaij(i,5) = x3
      C6Lim%Lim_etaij(j,5) = one - x3
      C6Lim%Lim_etaij(i,6) = one
      C6Lim%Lim_etaij(j,6) = zero
      C6Lim%Lim_etaij(5,6) = one - x3

      eta51 = C6Lim%Lim_etaij(1,5)
      eta52 = C6Lim%Lim_etaij(2,5)
      eta61 = C6Lim%Lim_etaij(1,6)
      eta62 = C6Lim%Lim_etaij(2,6)
      eta65 = one - x3

      call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
        spart,sqrts,xi1,xi2,C6Lim%flag)

      if (.not. C6Lim%flag) then

        C6Lim%npart = npart_nlo

        C6Lim%PartFrac = (/xi1,xi2/)

        ns(:,5) = [one,sin5i*cos_phi5,sin5i*sin_phi5,cos5i]
        ns(:,6) = [one,zero,zero,-sgn]

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = half*sqrts*x1*ns(:,5)
        p6 = half*sqrts*x2*ns(:,6)

        if(j .eq. 1) then
          p1 = p1 - p6
        elseif(j .eq. 2) then
          p2 = p2 - p6
        endif

        pv = p1 + p2 - p5
        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C6Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        C6Lim%AmpMom(:,1) = p1
        C6Lim%AmpMom(:,2) = p2
        C6Lim%AmpMom(:,5) = p5

        C6Lim%Lim_z(2) = one/(one - x2)
        C6Lim%Lim_sij(j,6) = spart*x2*HardProc%Lim_etaij(j,6)

        if(present(opt_etas)) call fill_etas(C6Lim,ns,opt_etas)

        C6Lim%wgt = one/8.0_dp/pi * kallenF &
             * (spart/two)**2 & 
             * x1*x2          & 
             * one/two/spart  & 
             * spart/MV2      &
             * jac

      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                       STRUCTURE 6: S5 + S6 FLM.                       !!
    !!-----------------------------------------------------------------------!!

    if(present(S5S6Lim)) then

      !! ----------------------------- S5 + S6 ----------------------------- !!

      S5S6Lim%npart = npart_lo       
       
      sqrts = mv
      spart = mv2
      S5S6Lim%PartFrac = sqrt(tau)*[exp(+ylab),exp(-ylab)]

      ns(:,5) = [one,sin5i*cos_phi5,sin5i*sin_phi5,cos5i]
      ns(:,6) = [one,sin6j*cos_phi6,sin6j*sin_phi6,cos6j]

      p1 = half*sqrts*ns(:,1)
      p2 = half*sqrts*ns(:,2)

      pv = p1 + p2
      call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, S5S6Lim)
      ns(:,3) = (/one,nlept(:,1)/)
      ns(:,4) = (/one,nlept(:,2)/)

      S5S6Lim%AmpMom(:,1) = p1 
      S5S6Lim%AmpMom(:,2) = p2

      S5S6Lim%Lim_etaij(1,2)   = one
      S5S6Lim%Lim_etaij(1:2,5) = HardProc%Lim_etaij(1:2,5)
      S5S6Lim%Lim_etaij(1:2,6) = HardProc%Lim_etaij(1:2,6)
      S5S6Lim%Lim_etaij(5,6)   = HardProc%Lim_etaij(5,6)
      S5S6Lim%Lim_etaij(6,5)   = HardProc%Lim_etaij(5,6)

      S5S6Lim%Lim_KinInv(1) = x1*half*sqrts
      S5S6Lim%Lim_KinInv(2) = x2*half*sqrts

      if(present(opt_etas)) call fill_etas(S5S6Lim,ns,opt_etas)

      S5S6Lim%wgt = one/8.0_dp/pi * kallenF &
           * (spart/two)**2 &
           * x1*x2          &
           * one/two/spart  &
           * spart/MV2      &
           * jac

      !! -------------------------- C5i + S5 + S6 -------------------------- !!

      if(present(C5S5S6Lim)) then

        C5S5S6Lim%PartFrac = C5S5S6Lim%PartFrac

        ns(:,5) = [one,zero,zero,sgn]
        ns(:,6) = [one,sin6j*cos_phi6,sin6j*sin_phi6,cos6j]

        C5S5S6Lim%Lim_etaij = S5S6Lim%Lim_etaij
        C5S5S6Lim%Lim_etaij(i,5) = zero
        C5S5S6Lim%Lim_etaij(j,5) = one
        C5S5S6Lim%Lim_etaij(5,6) = one - x4

        C5S5S6Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
        C5S5S6Lim%Lim_z(1) = x1

        if(present(opt_etas)) call fill_etas(C5S5S6Lim,ns,opt_etas)

        C5S5S6Lim%wgt =  one/8.0_dp/pi * kallenF &
           * (spart/two)**2 & 
           * x1*x2          & 
           * one/two/spart  & 
           * spart/mv2      &
           * jac
      endif

      !! -------------------------- C6j + S5 + S6 -------------------------- !!

      if(present(C6S5S6Lim)) then

        C6S5S6Lim%PartFrac = C6S5S6Lim%PartFrac

        ns(:,5) = [one,sin5i*cos_phi5,sin5i*sin_phi5,cos5i]
        ns(:,6) = [one,zero,zero,-sgn]

        C6S5S6Lim%Lim_etaij = S5S6Lim%Lim_etaij
        C6S5S6Lim%Lim_etaij(i,6) = one
        C6S5S6Lim%Lim_etaij(j,6) = zero
        C6S5S6Lim%Lim_etaij(5,6) = one - x3

        C6S5S6Lim%Lim_sij(j,6) = spart*x2*HardProc%Lim_etaij(j,6)
        C6S5S6Lim%Lim_z(2) = x2

        if(present(opt_etas)) call fill_etas(C6S5S6Lim,ns,opt_etas)

        C6S5S6Lim%wgt =  one/8.0_dp/pi * kallenF &
           * (spart/two)**2 & 
           * x1*x2          & 
           * one/two/spart  & 
           * spart/mv2      &
           * jac
      endif

    endif

    !!-----------------------------------------------------------------------!!
    !!                      STRUCTURE 7: C5i + S6 FLM.                       !!
    !!-----------------------------------------------------------------------!!

    if(present(C5S6Lim)) then

      !! x3 -> 0
      C5S6Lim%Lim_etaij(1,2) = one
      C5S6Lim%Lim_etaij(i,5) = zero
      C5S6Lim%Lim_etaij(j,5) = one
      C5S6Lim%Lim_etaij(i,6) = one - x4
      C5S6Lim%Lim_etaij(j,6) = x4
      C5S6Lim%Lim_etaij(5,6) = one - x4
      C5S6Lim%Lim_etaij(6,5) = one - x4

      eta51 = C5S6Lim%Lim_etaij(1,5)
      eta52 = C5S6Lim%Lim_etaij(2,5)
      eta61 = C5S6Lim%Lim_etaij(1,6)
      eta62 = C5S6Lim%Lim_etaij(2,6)
      eta65 = one - x4

      call get_part_s_xi_nlo(mv2,ylab,x1,eta51,eta52,spart,sqrts,xi1,xi2,C5S6Lim%flag)

      if (.not. C5S6Lim%flag) then 

        C5S6Lim%npart = npart_lo
         
      !! ---------------------------- C5i + S6 ----------------------------- !!

        C5S6Lim%PartFrac = (/xi1,xi2/)

        ns(:,5) = [one,zero,zero,sgn]
        ns(:,6) = [one,sin6j*cos_phi6,sin6j*sin_phi6,cos6j]

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = half*sqrts*x1*ns(:,5)

        if(i .eq. 1) then
          p1 = p1 - p5
        elseif(i .eq. 2) then
          p2 = p2 - p5
        endif

        pv = p1 + p2
        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C5S6Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        C5S6Lim%AmpMom(:,1) = p1 
        C5S6Lim%AmpMom(:,2) = p2

        !! E6
        C5S6Lim%Lim_KinInv(1) = x2*half*sqrts

        C5S6Lim%Lim_z(1)  = one/(one - x1)
        C5S6Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)

        if(present(opt_etas)) call fill_etas(C5S6Lim,ns,opt_etas)

        C5S6Lim%wgt= one/8.0_dp/pi * kallenF &
             * (spart/two)**2 & 
             * x1*x2          &
             * one/two/spart  & 
             * spart/MV2      &
             * jac

      !! -------------------------- C5i + C6j + S6 ------------------------- !!

        if(present(C5C6S6Lim)) then

          C5C6S6Lim%PartFrac = C5S6Lim%PartFrac

          ns(:,6) = [one,zero,zero,-sgn]

          C5C6S6Lim%Lim_etaij = C5S6Lim%Lim_etaij
          C5C6S6Lim%Lim_etaij(i,6) = one
          C5C6S6Lim%Lim_etaij(j,6) = zero
          C5C6S6Lim%Lim_etaij(5,6) = one - x3

          C5C6S6Lim%Lim_z(1) = one/(one-x1)
          C5C6S6Lim%Lim_z(2) = x2

          C5C6S6Lim%Lim_sij(i,5) = C5S6Lim%Lim_sij(i,5)
          C5C6S6Lim%Lim_sij(j,6) = spart*x2*HardProc%Lim_etaij(j,6)

          C5C6S6Lim%wgt = one/8.0_dp/pi * kallenF &
               * (spart/two)**2 & 
               * x1*x2          & 
               * one/two/spart  & 
               * spart/MV2      &
               * jac
        endif
      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                      STRUCTURE 8: C6j + S5 FLM.                       !!
    !!-----------------------------------------------------------------------!!

    if(present(C6S5Lim)) then

      !! x4 -> 0
      C6S5Lim%Lim_etaij(1,2) = one
      C6S5Lim%Lim_etaij(i,5) = x3
      C6S5Lim%Lim_etaij(j,5) = one - x3
      C6S5Lim%Lim_etaij(i,6) = one
      C6S5Lim%Lim_etaij(j,6) = zero
      C6S5Lim%Lim_etaij(5,6) = one - x3
      C6S5Lim%Lim_etaij(6,5) = one - x3

      eta51 = C6S5Lim%Lim_etaij(1,5)
      eta52 = C6S5Lim%Lim_etaij(2,5)
      eta61 = C6S5Lim%Lim_etaij(1,6)
      eta62 = C6S5Lim%Lim_etaij(2,6)
      eta65 = one - x3

      call get_part_s_xi_nlo(mv2,ylab,x2,eta61,eta62,spart,sqrts,xi1,xi2,C6S5Lim%flag)

      if (.not. C6S5Lim%flag) then 

        C6S5Lim%npart = npart_lo         
         
      !! ------------------------------- C6j + S5 -------------------------- !!

        C6S5Lim%PartFrac = (/xi1,xi2/)

        ns(:,5) = [one,sin5i*cos_phi5,sin5i*sin_phi5,cos5i]
        ns(:,6) = [one,zero,zero,-sgn]

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p6 = half*sqrts*x2*ns(:,6)

        if(j .eq. 1) then
          p1 = p1 - p6
        elseif(j .eq. 2) then
          p2 = p2 - p6
        endif

        pv = p1 + p2
        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C6S5Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        C6S5Lim%AmpMom(:,1) = p1 
        C6S5Lim%AmpMom(:,2) = p2

        !! E5
        C6S5Lim%Lim_KinInv(1) = x1*half*sqrts

        C6S5Lim%Lim_z(2)  = one/(one - x2)
        C6S5Lim%Lim_sij(j,6) = spart*x2*HardProc%Lim_etaij(j,6)

        if(present(opt_etas)) call fill_etas(C6S5Lim,ns,opt_etas)

        C6S5Lim%wgt= one/8.0_dp/pi * kallenF &
             * (spart/two)**2 &
             * x1*x2          &
             * one/two/spart  &
             * spart/MV2      &
             * jac

      !! -------------------------- C6j + C5i + S5 ------------------------- !!

        if(present(C5C6S5Lim)) then

          C5C6S5Lim%PartFrac = C6S5Lim%PartFrac

          ns(:,5) = [one,zero,zero,sgn]

          C5C6S5Lim%Lim_etaij = C6S5Lim%Lim_etaij
          C5C6S5Lim%Lim_etaij(i,5) = zero
          C5C6S5Lim%Lim_etaij(j,5) = one
          C5C6S5Lim%Lim_etaij(5,6) = one - x4

          C5C6S5Lim%Lim_z(1) = x1
          C5C6S5Lim%Lim_z(2) = one/(one-x2)

          C5C6S5Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
          C5C6S5Lim%Lim_sij(j,6) = C6S5Lim%Lim_sij(j,6)

          C5C6S5Lim%wgt = one/8.0_dp/pi * kallenF &
               * (spart/two)**2 & 
               * x1*x2          & 
               * one/two/spart  & 
               * spart/MV2      &
               * jac
        endif
      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                      STRUCTURE 9: C5i + C6j FLM.                      !!
    !!-----------------------------------------------------------------------!!
    if(present(C5C6Lim)) then

      !! x3,x4 -> 0
      C5C6Lim%Lim_etaij(1,2) = one
      C5C6Lim%Lim_etaij(i,5) = zero
      C5C6Lim%Lim_etaij(j,5) = one
      C5C6Lim%Lim_etaij(i,6) = one
      C5C6Lim%Lim_etaij(j,6) = zero
      C5C6Lim%Lim_etaij(5,6) = one
      C5C6Lim%Lim_etaij(6,5) = one

      eta51 = C5C6Lim%Lim_etaij(1,5)
      eta52 = C5C6Lim%Lim_etaij(2,5)
      eta61 = C5C6Lim%Lim_etaij(1,6)
      eta62 = C5C6Lim%Lim_etaij(2,6)
      eta65 = one

      call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
        spart,sqrts,xi1,xi2,C5C6Lim%flag)

      if (.not. C5C6Lim%flag) then 

        C5C6Lim%npart = npart_lo         
         
        C5C6Lim%PartFrac = (/xi1,xi2/)

        ns(:,5) = [one,zero,zero, sgn]
        ns(:,6) = [one,zero,zero,-sgn]

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = half*sqrts*x1*ns(:,5)
        p6 = half*sqrts*x2*ns(:,6)

        if(i .eq. 1) then
          p1 = p1 - p5
          p2 = p2 - p6
        elseif(i .eq. 2) then
          p1 = p1 - p6
          p2 = p2 - p5
        endif

        pv = p1 + p2
        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C5C6Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        C5C6Lim%AmpMom(:,1) = p1
        C5C6Lim%AmpMom(:,2) = p2

        C5C6Lim%Lim_z(1) = one/(one - x1)
        C5C6Lim%Lim_z(2) = one/(one - x2)
        C5C6Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
        C5C6Lim%Lim_sij(j,6) = spart*x2*HardProc%Lim_etaij(j,6)

        if(present(opt_etas)) call fill_etas(C5C6Lim,ns,opt_etas)

        C5C6Lim%wgt = one/8.0_dp/pi * kallenF &
             * (spart/two)**2 & 
             * x1*x2          & 
             * one/two/spart  & 
             * spart/MV2      &
             * jac
      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                STRUCTURE 10: S5 + C5i + S6 + C6j FLM.                 !!
    !!-----------------------------------------------------------------------!!
    if(present(C5C6S5S6Lim)) then

      !! x3,x4 -> 0
      C5C6S5S6Lim%Lim_etaij(1,2) = one
      C5C6S5S6Lim%Lim_etaij(i,5) = zero
      C5C6S5S6Lim%Lim_etaij(j,5) = one
      C5C6S5S6Lim%Lim_etaij(i,6) = one
      C5C6S5S6Lim%Lim_etaij(j,6) = zero
      C5C6S5S6Lim%Lim_etaij(5,6) = one
      C5C6S5S6Lim%Lim_etaij(6,5) = one

      sqrts = mv
      spart = mv2
      C5C6S5S6Lim%PartFrac = sqrt(tau)*[exp(+ylab),exp(-ylab)]

      ns(:,5) = [one,zero,zero, sgn]
      ns(:,6) = [one,zero,zero,-sgn]

      p1 = half*sqrts*ns(:,1)
      p2 = half*sqrts*ns(:,2)

      pv = p1 + p2
      call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C5C6S5S6Lim)
      ns(:,3) = (/one,nlept(:,1)/)
      ns(:,4) = (/one,nlept(:,2)/)

      C5C6S5S6Lim%AmpMom(:,1) = p1 
      C5C6S5S6Lim%AmpMom(:,2) = p2

      C5C6S5S6Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
      C5C6S5S6Lim%Lim_sij(j,6) = spart*x2*HardProc%Lim_etaij(j,6)

      if(present(opt_etas)) call fill_etas(C5C6S5S6Lim,ns,opt_etas)

      C5C6S5S6Lim%wgt = one/8.0_dp/pi * kallenF &
           * (spart/two)**2 & 
           * x1*x2          & 
           * one/two/spart  & 
           * spart/MV2      &
           * jac

    endif

  end subroutine kinematics_nnlo_ii_5i6j

  !!=========================================================================!!
  !! Kinematics for IS double-collinear configuration                        !!
  !! Collinearity to the same leg i = {1,2}                                  !!
  !! j is the other leg {2,1}                                                !!
  !!                                                                         !!
  !! opt_etas: indices k for eta_kl, l=1,...,n.                              !!
  !! eta_ij, where i,j refer to initial-state particles and                  !!
  !! radiated particles, are always computed by default                      !!
  !!=========================================================================!!
  subroutine kinematics_nnlo_ii_5i6i(yr, i, j, &
       HardProc, C5Lim, C6Lim)
    use mod_kinematics_gen
    use mod_kinematics_lept
    use mod_aux_kinematics
    !
    real(dp),        intent(in)    :: yr(:)
    integer,         intent(in)    :: i, j
    type(KinConfig), intent(inout) :: HardProc
    type(KinConfig), optional, intent(inout) :: C5Lim, C6Lim
    !---
    real(dp) :: mv, mv2, tau, ylab, jac
    real(dp) :: p1(4), p2(4), p5(4), p6(4), pv(4)
    real(dp) :: sgn, kallenF
    real(dp) :: x1, x2, x3, x4, x5, x6
    real(dp) :: cos5i, sin5i, cos6i, sin6i
    real(dp) :: phi5, phi6, cos_phi5, cos_phi6, sin_phi5, sin_phi6
    real(dp) :: xi1, xi2, spart, sqrts
    real(dp) :: nlept(3,2), ns(4,6)
    real(dp) :: eta5i, eta6i, eta5j, eta6j
    real(dp) :: eta51, eta61, eta52, eta62, eta65
    
    !-- process-dependent part
    call get_tauy(yr(1:tauy_max),tau,ylab,jac)
    mv2 = tau*sh
    mv  = sqrt(mv2)
    
    !-- Initialization of kinematic configurations
    call initialize_config(HardProc)                             !  1
    if(present(C5Lim))       call initialize_config(C5Lim)       !  2
    if(present(C6Lim))       call initialize_config(C6Lim)       !  3

    x1 = yr(xE5)   ! E_5
    x2 = yr(xE6)   ! E_6
    x3 = yr(xRHO5) ! \theta_{5i}
    x4 = yr(xRHO6) ! \theta_{6i}
    x5 = yr(xPHI5) ! \phi_{5}
    x6 = yr(xPHI6) ! related to \lambda
    
    !! Choose sign of \cos\theta depending on emitter 
    !! {i,j} = {1,2} or {i,j} = {2,1}
    if(i .eq. 1) then
       sgn = +one
    elseif(i .eq. 2) then
       sgn = -one
    else
       print *, 'wrong emitter'
       sgn = zero
       stop
    endif
    
    !!-----------------------------------------------------------------------!!
    !!                      STRUCTURE 1: Hard Process.                       !!
    !!-----------------------------------------------------------------------!!
    
    !! Azimuthal variables
    phi5 = twopi*x5 
    phi6 = twopi*x6
    cos_phi5 = cos(phi5)
    sin_phi5 = sin(phi5)
    cos_phi6 = cos(phi6)
    sin_phi6 = sin(phi6)
    
    !! Polar variables
    eta5i = x3
    eta5j = one - x3
    eta6i = x4
    eta6j = one - x4
    
    cos5i = sgn*( one - 2*x3)
    sin5i = two*sqrt(x3*(one-x3))
    
    cos6i = sgn*( one - 2*x4)
    sin6i = two*sqrt(x4*(one-x4))
    
    eta65 = half*(one - sin5i*sin6i*cos(phi5-phi6) - cos5i*cos6i)
    
    !! set eta variables in the limit
    HardProc%Lim_etaij(1,2) = one
    HardProc%Lim_etaij(i,5) = eta5i
    HardProc%Lim_etaij(j,5) = eta5j
    HardProc%Lim_etaij(i,6) = eta6i
    HardProc%Lim_etaij(j,6) = eta6j
    HardProc%Lim_etaij(5,6) = eta65
    
    eta51 = HardProc%Lim_etaij(1,5)
    eta52 = HardProc%Lim_etaij(2,5)
    eta61 = HardProc%Lim_etaij(1,6)
    eta62 = HardProc%Lim_etaij(2,6)
    
    ns(:,1) = [one,zero,zero, one]
    ns(:,2) = [one,zero,zero,-one]
    
    ns(:,5) = [one,sin5i*cos_phi5,sin5i*sin_phi5,cos5i]
    ns(:,6) = [one,sin6i*cos_phi6,sin6i*sin_phi6,cos6i]
    
    !! Compute partonic fractions {xi1, xi2}, partonic s and its square root
    !! and check if point is physical
    call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
         spart,sqrts,xi1,xi2,HardProc%flag)
    
    if (.not. HardProc%flag) then

       HardProc%npart = npart_nnlo
       
       HardProc%PartFrac = (/xi1,xi2/)
       
       p1 = half*sqrts*ns(:,1)
       p2 = half*sqrts*ns(:,2)
       p5 = half*sqrts*x1*ns(:,5)
       p6 = half*sqrts*x2*ns(:,6)
       
       pv = p1 + p2 - p5 - p6
       
       call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, HardProc)
       ns(:,3) = (/one,nlept(:,1)/)
       ns(:,4) = (/one,nlept(:,2)/)
       
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,5) = p5
       HardProc%AmpMom(:,6) = p6
       
       HardProc%wgt = one/8.0_dp/pi * kallenF &
            * (spart/two)**2 & 
            * x1*x2          & 
            * one/two/spart  &
            * spart/mv2      &
            * jac
       
    endif

    !!-----------------------------------------------------------------------!!
    !!                        STRUCTURE 2: C5i FLM.                          !!
    !!-----------------------------------------------------------------------!!
    
    if(present(C5Lim)) then
       
       !! x3 -> 0
       C5Lim%Lim_etaij(1,2) = one
       C5Lim%Lim_etaij(i,5) = zero
       C5Lim%Lim_etaij(j,5) = one
       C5Lim%Lim_etaij(i,6) = x4
       C5Lim%Lim_etaij(j,6) = one - x4
       C5Lim%Lim_etaij(5,6) = x4
       
       eta51 = C5Lim%Lim_etaij(1,5)
       eta52 = C5Lim%Lim_etaij(2,5)
       eta61 = C5Lim%Lim_etaij(1,6)
       eta62 = C5Lim%Lim_etaij(2,6)
       eta65 = x4
       
       call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
            spart,sqrts,xi1,xi2,C5Lim%flag)
       
       if (.not. C5Lim%flag) then 
          
          C5Lim%npart = npart_nlo         
          
          C5Lim%PartFrac = (/xi1,xi2/)
          
          ns(:,5) = [one,zero,zero,sgn]
          ns(:,6) = [one,sin6i*cos_phi6,sin6i*sin_phi6,cos6i]
          
          p1 = half*sqrts*ns(:,1)
          p2 = half*sqrts*ns(:,2)
          p5 = half*sqrts*x1*ns(:,5)
          p6 = half*sqrts*x2*ns(:,6)
          
          if(i .eq. 1) then
             p1 = p1 - p5
          elseif(i .eq. 2) then
             p2 = p2 - p5
          endif
          
          pv = p1 + p2 - p6
          call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C5Lim)
          ns(:,3) = (/one,nlept(:,1)/)
          ns(:,4) = (/one,nlept(:,2)/)
          
          C5Lim%AmpMom(:,1) = p1
          C5Lim%AmpMom(:,2) = p2
          C5Lim%AmpMom(:,5) = p6
          
          C5Lim%Lim_z(1) = one/(one - x1)
          C5Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
          
          C5Lim%wgt = one/8.0_dp/pi * kallenF &
               * (spart/two)**2 & 
               * x1*x2          & 
               * one/two/spart  & 
               * spart/MV2      &
               * jac
          
       endif
    endif
    
    !!-----------------------------------------------------------------------!!
    !!                        STRUCTURE 5: C6i FLM.                          !!
    !!-----------------------------------------------------------------------!!
    
    if(present(C6Lim)) then
       
       !! x4 -> 0
       C6Lim%Lim_etaij(1,2) = one
       C6Lim%Lim_etaij(i,5) = x3
       C6Lim%Lim_etaij(j,5) = one - x3
       C6Lim%Lim_etaij(i,6) = zero
       C6Lim%Lim_etaij(j,6) = one
       C6Lim%Lim_etaij(5,6) = x3

       eta51 = C6Lim%Lim_etaij(1,5)
       eta52 = C6Lim%Lim_etaij(2,5)
       eta61 = C6Lim%Lim_etaij(1,6)
       eta62 = C6Lim%Lim_etaij(2,6)
       eta65 = x3
       
       call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
            spart,sqrts,xi1,xi2,C6Lim%flag)
       
       if (.not. C6Lim%flag) then
          
          C6Lim%npart = npart_nlo
          
          C6Lim%PartFrac = (/xi1,xi2/)
          
          ns(:,5) = [one,sin5i*cos_phi5,sin5i*sin_phi5,cos5i]
          ns(:,6) = [one,zero,zero,sgn]
          
          p1 = half*sqrts*ns(:,1)
          p2 = half*sqrts*ns(:,2)
          p5 = half*sqrts*x1*ns(:,5)
          p6 = half*sqrts*x2*ns(:,6)
          
          if(i .eq. 1) then
             p1 = p1 - p6
          elseif(i .eq. 2) then
             p2 = p2 - p6
          endif

          pv = p1 + p2 - p5
          call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C6Lim)
          ns(:,3) = (/one,nlept(:,1)/)
          ns(:,4) = (/one,nlept(:,2)/)
          
          C6Lim%AmpMom(:,1) = p1
          C6Lim%AmpMom(:,2) = p2
          C6Lim%AmpMom(:,5) = p5
          
          C6Lim%Lim_z(2) = one/(one - x2)
          C6Lim%Lim_sij(i,6) = spart*x2*HardProc%Lim_etaij(i,6)
          
          C6Lim%wgt = one/8.0_dp/pi * kallenF &
               * (spart/two)**2 & 
               * x1*x2          & 
               * one/two/spart  & 
               * spart/MV2      &
               * jac
          
       endif
    end if

  end subroutine kinematics_nnlo_ii_5i6i

end module mod_kinematics_nnlo_dc_ii

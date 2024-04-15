!-- Double-collinear phase space, IS + FS radiation.
!-- Damping factor not included
module mod_kinematics_nnlo_dc_if
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
  integer, parameter :: xPHI5 = kNNLO_min + 3
  integer, parameter :: xPHI6 = kNNLO_min + 5
  !!
  integer, parameter :: xPHIe = kNNLO_min + 6 ! e -> FS emitter
  integer, parameter :: xRHOe = kNNLO_min + 7 ! e -> FS emitter
  !
  integer, parameter :: npart_lo = 4, npart_nlo = 5, npart_nnlo = 6

  public :: kinematics_nnlo_if_5i6k

  private

contains

  !!=========================================================================!!
  !! Kinematics for IS-FS double-collinear configuration                     !!
  !! i,j -> initial-state emitter                                            !!
  !! k,l -> final-state emitter                                              !!
  !! e.g. i \in {1,2} and k \in {3,4}                                        !!
  !!                                                                         !!
  !! opt_etas: indices k for eta_kl, l=1,...,n.                              !!
  !! eta_ij, where i,j refer to initial-state particles and                  !!
  !! radiated particles, are always computed by default                      !!
  !!=========================================================================!!
  subroutine kinematics_nnlo_if_5i6k(yr, i, j, k, l, &
    HardProc, S6Lim, C6Lim, C6S6Lim,    &
    S5Lim, S5S6Lim, C6S5Lim, C6S5S6Lim, &
    C5Lim, C5S6Lim, C5C6Lim, C5C6S6Lim, &
    C5S5Lim, C5S5S6Lim, C5C6S5Lim, C5C6S5S6Lim, &
    opt_etas)
    real(dp),        intent(in)    :: yr(:)
    integer,         intent(in)    :: i, j, k, l
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
    real(dp) :: sgn
    real(dp) :: x1, x2, x3, x4, x5, x6, x7, x8
    real(dp) :: cos5i, sin5i, cos6k, sin6k, cosk, sink
    real(dp) :: phi5, phi6, phik
    real(dp) :: cos_phi5, cos_phi6, cos_phik, sin_phi5, sin_phi6, sin_phik
    real(dp) :: xi1, xi2, spart, sqrts
    real(dp), dimension(4) :: a1, a2, pk, pl
    real(dp) :: Emax, E6, Ek, Jf, Q2, n6aux(3), ns(4,6)
    real(dp) :: eta5i, eta5j
    real(dp) :: eta51, eta52

    !-- process-dependent part
    call get_tauy(yr(1:tauy_max),tau,ylab,jac)
    mv2 = tau*sh
    mv  = sqrt(mv2)

    !-- Initialization of kinematic configurations
    call initialize_config(HardProc)                             !  1
    if(present(S6Lim))       call initialize_config(S6Lim)       !  2
    if(present(C6Lim))       call initialize_config(C6Lim)       !  3
    if(present(C6S6Lim))     call initialize_config(C6S6Lim)     !  4
    if(present(S5Lim))       call initialize_config(S5Lim)       !  5
    if(present(S5S6Lim))     call initialize_config(S5S6Lim)     !  6
    if(present(C6S5Lim))     call initialize_config(C6S5Lim)     !  7
    if(present(C6S5S6Lim))   call initialize_config(C6S5S6Lim)   !  8
    if(present(C5Lim))       call initialize_config(C5Lim)       !  9
    if(present(C5S6Lim))     call initialize_config(C5S6Lim)     ! 10
    if(present(C5C6Lim))     call initialize_config(C5C6Lim)     ! 11
    if(present(C5C6S6Lim))   call initialize_config(C5C6S6Lim)   ! 12
    if(present(C5S5Lim))     call initialize_config(C5S5Lim)     ! 13
    if(present(C5S5S6Lim))   call initialize_config(C5S5S6Lim)   ! 14
    if(present(C5C6S5Lim))   call initialize_config(C5C6S5Lim)   ! 15
    if(present(C5C6S5S6Lim)) call initialize_config(C5C6S5S6Lim) ! 16

    x1 = yr(xE5)    ! E_5
    x2 = yr(xE6)    ! E_6
    x3 = yr(xRHO5)  ! \theta_{5i}
    x4 = yr(xRHO6)  ! \theta_{6k}, measured wrt the emitter
    x5 = yr(xPHI5)  ! \phi_{5}
    x6 = yr(xPHI6)  ! \phi_{6}, measured around the emitter
    x7 = yr(xPHIe)  ! azimuthal angle of the FS emitter
    x8 = yr(xRHOe)  ! polar angle of the FS emitter

    !! Choose sign of \cos\theta depending on IS emitter, i=1,2
    if(i .eq. 1) then
      sgn = +one
    elseif(i .eq. 2) then
      sgn = -one
    else
      print *, 'wrong emitter'
      sgn = zero
      stop
    endif

    !-- IS
    phi5 = twopi*x5
    cos_phi5 = cos(phi5)
    sin_phi5 = sin(phi5)

    eta5i = x3
    eta5j = one - x3
    cos5i = sgn*( one - 2*x3)
    sin5i = two*sqrt(x3*(one-x3))

    !-- FS, radiator
    phik = twopi*x7
    cos_phik = cos(phik)
    sin_phik = sin(phik)

    cosk = one - 2*x8
    sink = two*sqrt(x8*(one-x8))

    !-- FS, radiated
    phi6 = twopi*x6
    cos_phi6 = cos(phi6)
    sin_phi6 = sin(phi6)

    cos6k = one - 2*x4
    sin6k = two*sqrt(x4*(one-x4))

    !-- unit vectors
    a1 = [one, -cos_phik*cosk,-sin_phik*cosk, sink]
    a2 = [one, -sin_phik, cos_phik, zero]
    ns(:,k) = [one,cos_phik*sink,sin_phik*sink,cosk] !-- emitter

    n6aux = cos6k*ns(2:4,k) + sin6k*(cos_phi6*a1(2:4) + sin_phi6*a2(2:4))

    ns(:,1) = [one,zero,zero, one]
    ns(:,2) = [one,zero,zero,-one]
    ns(:,5) = [one,sin5i*cos_phi5,sin5i*sin_phi5,cos5i]
    ns(:,6) = [one,n6aux(1),n6aux(2),n6aux(3)]

    !!-----------------------------------------------------------------------!!
    !!                  STRUCTURE 1: Hard Process + S6 + C6                  !!
    !!-----------------------------------------------------------------------!!

    !! ------------------------------ Hard 5 ------------------------------- !!

    !! set eta variables in the limit
    HardProc%Lim_etaij(1,2) = one
    HardProc%Lim_etaij(i,5) = eta5i
    HardProc%Lim_etaij(j,5) = eta5j

    eta51 = HardProc%Lim_etaij(1,5)
    eta52 = HardProc%Lim_etaij(2,5)

    call get_part_s_xi_nlo(mv2,ylab,x1,eta51,eta52,spart,sqrts,xi1,xi2,HardProc%flag)

    if (HardProc%flag) then

       if(present(C6Lim))   C6Lim%flag = .true.
       if(present(S6Lim))   S6Lim%flag = .true.
       if(present(C6S6Lim)) C6S6Lim%flag = .true.

    else

      HardProc%npart = npart_nnlo

      !! ------------------------------ Hard 6 ------------------------------- !!

      HardProc%PartFrac = (/xi1,xi2/)
      Emax = sqrts*half
      E6 = x2*Emax

      p1 = Emax*ns(:,1)
      p2 = Emax*ns(:,2)
      p5 = x1*Emax*ns(:,5)
      p6 = x2*Emax*ns(:,6)

      pv = p1 + p2 - p5
      Q2 = scr(pv,pv)

      Jf = 2*scr(pv,ns(:,k)) - 4*Emax*x2*x4
      Ek = (Q2 - 2*E6*scr(pv,ns(:,6)))/Jf

      pk = Ek * ns(:,k)
      pl = pv - pk - p6
      ns(:,l) = pl/pl(1)

      !-- Check valid point according to energy of final-state particles
      if(Ek < zero .or. pl(1) < zero) then

        HardProc%flag = .true.

      else

        HardProc%AmpMom(:,1) = p1
        HardProc%AmpMom(:,2) = p2
        HardProc%AmpMom(:,5) = p5
        HardProc%AmpMom(:,6) = p6
        HardProc%AmpMom(:,k) = pk
        HardProc%AmpMom(:,l) = pl

        HardProc%Lim_etaij(1,6) = (one - ns(4,6))*half
        HardProc%Lim_etaij(2,6) = (one + ns(4,6))*half
        HardProc%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,6)))*half

        if(present(opt_etas)) call fill_etas(HardProc,ns,opt_etas)

        HardProc%Lim_etaij(k,6) = x4

        HardProc%wgt = one/(2*spart)      &
          * x1 * spart**2/mv2             &
          * (x2/pi) * (Ek/Jf) * (Emax**2) &
          * jac

      endif
    !! -------------------------------- S6 --------------------------------- !!
      if(present(S6Lim)) then

        S6Lim%npart = npart_nlo

        S6Lim%PartFrac = HardProc%PartFrac

        Jf = 2*scr(pv,ns(:,k))
        Ek = Q2/Jf

        pk = Ek * ns(:,k)
        pl = pv - pk
        ns(:,l) = pl/pl(1)

        !-- Check valid point according to energy of final-state particles
        if(Ek < zero .or. pl(1) < zero) then

          S6Lim%flag = .true.

        else

          Emax = half*sqrts

          S6Lim%AmpMom(:,1) = p1
          S6Lim%AmpMom(:,2) = p2
          S6Lim%AmpMom(:,5) = p5
          S6Lim%AmpMom(:,k) = pk
          S6Lim%AmpMom(:,l) = pl

          S6Lim%Lim_KinInv(1) = x2*Emax

          S6Lim%Lim_etaij(1,2) = one
          S6Lim%Lim_etaij(i,5) = HardProc%Lim_etaij(i,5)
          S6Lim%Lim_etaij(j,5) = HardProc%Lim_etaij(j,5)
          S6Lim%Lim_etaij(1,6) = (one - ns(4,6))*half
          S6Lim%Lim_etaij(2,6) = (one + ns(4,6))*half
          S6Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,6)))*half
          S6Lim%Lim_etaij(6,5) = S6Lim%Lim_etaij(5,6)

          if(present(opt_etas)) call fill_etas(S6Lim,ns,opt_etas)

          S6Lim%Lim_etaij(k,6) = x4

          S6Lim%wgt = one/(2*spart)         &
            * x1 * spart**2/mv2             &
            * (x2/pi) * (Ek/Jf) * (Emax**2) &
            * jac

    !! ------------------------------ S6 + C6k ----------------------------- !!
          if(present(C6S6Lim)) then

            C6S6Lim%PartFrac = HardProc%PartFrac
            ns(:,6) = ns(:,k)

            C6S6Lim%Lim_etaij(1,2) = one
            C6S6Lim%Lim_etaij(i,5) = HardProc%Lim_etaij(i,5)
            C6S6Lim%Lim_etaij(j,5) = HardProc%Lim_etaij(j,5)
            C6S6Lim%Lim_etaij(1,6) = (one - ns(4,k))*half
            C6S6Lim%Lim_etaij(2,6) = (one + ns(4,k))*half
            C6S6Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,k)))*half

            if(present(opt_etas)) call fill_etas(C6S6Lim,ns,opt_etas)

            C6S6Lim%Lim_etaij(k,6) = zero

            C6S6Lim%Lim_z(2) = E6/Ek
            C6S6Lim%Lim_sij(k,6) = 4*E6*Ek*x4

            C6S6Lim%wgt = S6Lim%wgt

          endif
        endif
      endif

    !! -------------------------------- C6k --------------------------------- !!
      if(present(C6Lim)) then

        C6Lim%PartFrac = HardProc%PartFrac
        ns(:,6) = ns(:,k)

        Jf = 2*scr(pv,ns(:,k))
        Ek = (Q2 - 2*E6*scr(pv,ns(:,k)))/Jf
        pk = (Ek + E6) * ns(:,k)

        pl = pv - pk
        ns(:,l) = pl/pl(1)

        !-- Check valid point according to energy of final-state particles
        if(Ek < zero .or. pl(1) < zero) then

           C6Lim%flag = .true.

        else

          C6Lim%npart = npart_nlo

          C6Lim%AmpMom(:,1) = p1
          C6Lim%AmpMom(:,2) = p2
          C6Lim%AmpMom(:,5) = p5
          C6Lim%AmpMom(:,k) = pk
          C6Lim%AmpMom(:,l) = pl

          C6Lim%Lim_etaij(1,2) = one
          C6Lim%Lim_etaij(i,5) = HardProc%Lim_etaij(i,5)
          C6Lim%Lim_etaij(j,5) = HardProc%Lim_etaij(j,5)
          C6Lim%Lim_etaij(1,6) = (one - ns(4,k))*half
          C6Lim%Lim_etaij(2,6) = (one + ns(4,k))*half
          C6Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,k)))*half

          if(present(opt_etas)) call fill_etas(C6Lim,ns,opt_etas)

          C6Lim%Lim_etaij(k,6) = zero

          C6Lim%Lim_z(2) = E6/(E6+Ek)
          C6Lim%Lim_sij(k,6) = 4*E6*Ek*x4

          C6Lim%wgt = one/(2*spart)         &
            * x1 * spart**2/mv2             &
            * (x2/pi) * (Ek/Jf) * (Emax**2) &
            * jac

        endif
      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                     STRUCTURE 2: Soft 5 + S6 + C6                     !!
    !!-----------------------------------------------------------------------!!

    if(present(S5Lim)) then
    !! ------------------------------ Soft 5 ------------------------------- !!

      ns(:,6) = [one,n6aux(1),n6aux(2),n6aux(3)]

      S5Lim%Lim_etaij(1,2) = HardProc%Lim_etaij(1,2)
      S5Lim%Lim_etaij(i,5) = HardProc%Lim_etaij(i,5)
      S5Lim%Lim_etaij(j,5) = HardProc%Lim_etaij(j,5)

    !! ------------------------------ Hard 6 ------------------------------- !!

      S5Lim%PartFrac = sqrt(tau)*[exp(+ylab),exp(-ylab)]
      sqrts = mv
      spart = mv2
      Emax = sqrts*half

      p1 = Emax*ns(:,1)
      p2 = Emax*ns(:,2)
      p5 = x1*Emax*ns(:,5)

      S5Lim%Lim_KinInv(1) = x1*Emax

      pv = p1 + p2
      Q2 = spart !-- scr(pv,pv)

      E6 = x2*Emax
      p6 = E6*ns(:,6)

      Jf = 2*scr(pv,ns(:,k)) - 4*Emax*x2*x4
      Ek = (Q2 - 2*E6*scr(pv,ns(:,6)))/Jf

      pk = Ek * ns(:,k)
      pl = pv - pk - p6
      ns(:,l) = pl/pl(1)

      !-- Check valid point according to energy of final-state particles
      if(Ek < zero .or. pl(1) < zero) then

        S5Lim%flag = .true.

      else

        S5Lim%npart = npart_nlo

        S5Lim%AmpMom(:,1) = p1
        S5Lim%AmpMom(:,2) = p2
        S5Lim%AmpMom(:,5) = p6
        S5Lim%AmpMom(:,k) = pk
        S5Lim%AmpMom(:,l) = pl

        S5Lim%Lim_etaij(1,6) = (one - ns(4,6))*half
        S5Lim%Lim_etaij(2,6) = (one + ns(4,6))*half
        S5Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,6)))*half

        if(present(opt_etas)) call fill_etas(S5Lim,ns,opt_etas)

        S5Lim%Lim_etaij(k,6) = x4

        S5Lim%wgt = one/(2*spart)         &
          * x1 * spart**2/mv2             &
          * (x2/pi) * (Ek/Jf) * (Emax**2) &
          * jac

      endif
    !! -------------------------------- S6 --------------------------------- !!
      if(present(S5S6Lim)) then

        S5S6Lim%PartFrac = S5Lim%PartFrac

        Jf = 2*scr(pv,ns(:,k))
        Ek = Q2/Jf

        pk = Ek * ns(:,k)
        pl = pv - pk
        ns(:,l) = pl/pl(1)

        !-- Check valid point according to energy of final-state particles
        if(Ek < zero .or. pl(1) < zero) then

          S5S6Lim%flag = .true.

        else

          S5S6Lim%npart = npart_lo

          S5S6Lim%AmpMom(:,1) = p1
          S5S6Lim%AmpMom(:,2) = p2
          S5S6Lim%AmpMom(:,k) = pk
          S5S6Lim%AmpMom(:,l) = pl

          S5S6Lim%Lim_KinInv(1) = S5Lim%Lim_KinInv(1)
          S5S6Lim%Lim_KinInv(2) = E6

          S5S6Lim%Lim_etaij(1,2) = one
          S5S6Lim%Lim_etaij(i,5) = S5Lim%Lim_etaij(i,5)
          S5S6Lim%Lim_etaij(j,5) = S5Lim%Lim_etaij(j,5)
          S5S6Lim%Lim_etaij(1,6) = (one - ns(4,6))*half
          S5S6Lim%Lim_etaij(2,6) = (one + ns(4,6))*half
          S5S6Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,6)))*half
          S5S6Lim%Lim_etaij(6,5) = S5S6Lim%Lim_etaij(5,6)

          if(present(opt_etas)) call fill_etas(S5S6Lim,ns,opt_etas)

          S5S6Lim%Lim_etaij(k,6) = x4

          S5S6Lim%wgt = one/(2*spart)       &
            * x1 * spart**2/mv2             &
            * (x2/pi) * (Ek/Jf) * (Emax**2) &
            * jac

    !! ------------------------------ S6 + C6k ----------------------------- !!
          if(present(C6S5S6Lim)) then

            C6S5S6Lim%PartFrac = S5Lim%PartFrac
            ns(:,6) = ns(:,k)

            C6S5S6Lim%Lim_etaij(1,2) = one
            C6S5S6Lim%Lim_etaij(i,5) = S5Lim%Lim_etaij(i,5)
            C6S5S6Lim%Lim_etaij(j,5) = S5Lim%Lim_etaij(j,5)
            C6S5S6Lim%Lim_etaij(1,6) = (one - ns(4,k))*half
            C6S5S6Lim%Lim_etaij(2,6) = (one + ns(4,k))*half
            C6S5S6Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,k)))*half

            if(present(opt_etas)) call fill_etas(C6S5S6Lim,ns,opt_etas)

            C6S5S6Lim%Lim_etaij(k,6) = zero

            C6S5S6Lim%Lim_KinInv(1) = S5Lim%Lim_KinInv(1)
            C6S5S6Lim%Lim_z(2) = E6/Ek
            C6S5S6Lim%Lim_sij(k,6) = 4*E6*Ek*x4

            C6S5S6Lim%wgt = S5S6Lim%wgt

          endif
      endif
    endif

    !! -------------------------------- C6k --------------------------------- !!
      if(present(C6S5Lim)) then

        C6S5Lim%PartFrac = S5Lim%PartFrac
        ns(:,6) = ns(:,k)

        Jf = 2*scr(pv,ns(:,k))
        Ek = (Q2 - 2*E6*scr(pv,ns(:,k)))/Jf
        pk = (Ek + E6) * ns(:,k)

        pl = pv - pk
        ns(:,l) = pl/pl(1)

        !-- Check valid point according to energy of final-state particles
        if(Ek < zero .or. pl(1) < zero) then

          C6S5Lim%flag = .true.

        else

          C6S5Lim%npart = npart_lo

          C6S5Lim%AmpMom(:,1) = p1
          C6S5Lim%AmpMom(:,2) = p2
          C6S5Lim%AmpMom(:,k) = pk
          C6S5Lim%AmpMom(:,l) = pl

          C6S5Lim%Lim_etaij(1,2) = one
          C6S5Lim%Lim_etaij(i,5) = S5Lim%Lim_etaij(i,5)
          C6S5Lim%Lim_etaij(j,5) = S5Lim%Lim_etaij(j,5)
          C6S5Lim%Lim_etaij(1,6) = (one - ns(4,k))*half
          C6S5Lim%Lim_etaij(2,6) = (one + ns(4,k))*half
          C6S5Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,k)))*half

          if(present(opt_etas)) call fill_etas(C6S5Lim,ns,opt_etas)

          C6S5Lim%Lim_etaij(k,6) = zero

          C6S5Lim%Lim_KinInv(1) = S5Lim%Lim_KinInv(1)
          C6S5Lim%Lim_z(2) = E6/(E6+Ek)
          C6S5Lim%Lim_sij(k,6) = 4*E6*Ek*x4

          C6S5Lim%wgt = one/(2*spart)       &
            * x1 * spart**2/mv2             &
            * (x2/pi) * (Ek/Jf) * (Emax**2) &
            * jac

        endif
      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                  STRUCTURE 3: Collinear 5 + S6 + C6                   !!
    !!-----------------------------------------------------------------------!!

    if(present(C5Lim)) then
    !! ------------------------------ Coll 5 ------------------------------- !!

      ns(:,6) = [one,n6aux(1),n6aux(2),n6aux(3)]

      !-- x3 -> 0
      C5Lim%Lim_etaij(1,2) = one
      C5Lim%Lim_etaij(i,5) = zero
      C5Lim%Lim_etaij(j,5) = one

      eta51 = C5Lim%Lim_etaij(1,5)
      eta52 = C5Lim%Lim_etaij(2,5)

      ns(:,5) = [one,zero,zero,sgn]

      call get_part_s_xi_nlo(mv2,ylab,x1,eta51,eta52,spart,sqrts,xi1,xi2,C5Lim%flag)

      if(C5Lim%flag) then

         if(present(C5S6Lim))   C5S6Lim%flag=.true.
         if(present(C5C6Lim))   C5C6Lim%flag=.true.
         if(present(C5C6S6Lim)) C5C6S6Lim%flag=.true.

      else
    !! ------------------------------ Hard 6 ------------------------------- !!

        C5Lim%npart = npart_nlo

        C5Lim%PartFrac = (/xi1,xi2/)
        Emax = sqrts*half

        C5Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
        C5Lim%Lim_z(1) = x1/(x1-one)

        p1 = Emax*ns(:,1)
        p2 = Emax*ns(:,2)
        p5 = x1*Emax*ns(:,5)

        if(i .eq. 1) then
          p1 = p1 - p5
        elseif(i .eq. 2) then
          p2 = p2 - p5
        endif

        pv = p1 + p2
        Q2 = scr(pv,pv)

        E6 = x2*Emax
        p6 = E6*ns(:,6)

        Jf = 2*scr(pv,ns(:,k)) - 4*Emax*x2*x4
        Ek = (Q2 - 2*E6*scr(pv,ns(:,6)))/Jf

        pk = Ek * ns(:,k)
        pl = pv - pk - p6
        ns(:,l) = pl/pl(1)

        !-- Check valid point according to energy of final-state particles
        if(Ek < zero .or. pl(1) < zero) then

          C5Lim%flag = .true.

        else

          C5Lim%AmpMom(:,1) = p1
          C5Lim%AmpMom(:,2) = p2
          C5Lim%AmpMom(:,5) = p6
          C5Lim%AmpMom(:,k) = pk
          C5Lim%AmpMom(:,l) = pl

          C5Lim%Lim_etaij(1,6) = (one - ns(4,6))*half
          C5Lim%Lim_etaij(2,6) = (one + ns(4,6))*half
          C5Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,6)))*half

          if(present(opt_etas)) call fill_etas(C5Lim,ns,opt_etas)

          C5Lim%Lim_etaij(k,6) = x4

          C5Lim%wgt = one/(2*spart)         &
            * x1 * spart**2/mv2             &
            * (x2/pi) * (Ek/Jf) * (Emax**2) &
            * jac

        endif
    !! -------------------------------- S6 --------------------------------- !!
        if(present(C5S6Lim)) then

          C5S6Lim%PartFrac = C5Lim%PartFrac

          Jf = 2*scr(pv,ns(:,k))
          Ek = Q2/Jf

          pk = Ek * ns(:,k)
          pl = pv - pk
          ns(:,l) = pl/pl(1)

          !-- Check valid point according to energy of final-state particles
          if(Ek < zero .or. pl(1) < zero) then

            C5S6Lim%flag = .true.

          else

            C5S6Lim%npart = npart_lo

            C5S6Lim%AmpMom(:,1) = p1
            C5S6Lim%AmpMom(:,2) = p2
            C5S6Lim%AmpMom(:,k) = pk
            C5S6Lim%AmpMom(:,l) = pl

            C5S6Lim%Lim_KinInv(1) = E6

            C5S6Lim%Lim_etaij(1,2) = one
            C5S6Lim%Lim_etaij(i,5) = C5Lim%Lim_etaij(i,5)
            C5S6Lim%Lim_etaij(j,5) = C5Lim%Lim_etaij(j,5)
            C5S6Lim%Lim_etaij(1,6) = (one - ns(4,6))*half
            C5S6Lim%Lim_etaij(2,6) = (one + ns(4,6))*half
            C5S6Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,6)))*half
            C5S6Lim%Lim_etaij(6,5) = C5S6Lim%Lim_etaij(5,6)

            if(present(opt_etas)) call fill_etas(C5S6Lim,ns,opt_etas)

            C5S6Lim%Lim_sij(i,5) = C5Lim%Lim_sij(i,5)
            C5S6Lim%Lim_z(1) = C5Lim%Lim_z(1)

            C5S6Lim%Lim_etaij(k,6) = x4

            C5S6Lim%wgt = one/(2*spart)       &
              * 4 * x1 * spart/mv2            &
              * (x2/pi) * (Ek/Jf) * (Emax**2) &
              * spart/4 * jac

    !! ------------------------------ S6 + C6k ----------------------------- !!
            if(present(C5C6S6Lim)) then

              C5C6S6Lim%PartFrac = C5Lim%PartFrac
              ns(:,6) = ns(:,k)

              C5C6S6Lim%Lim_etaij(1,2) = one
              C5C6S6Lim%Lim_etaij(i,5) = C5Lim%Lim_etaij(i,5)
              C5C6S6Lim%Lim_etaij(j,5) = C5Lim%Lim_etaij(j,5)
              C5C6S6Lim%Lim_etaij(1,6) = (one - ns(4,k))*half
              C5C6S6Lim%Lim_etaij(2,6) = (one + ns(4,k))*half
              C5C6S6Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,k)))*half

              if(present(opt_etas)) call fill_etas(C5C6S6Lim,ns,opt_etas)

              C5C6S6Lim%Lim_etaij(k,6) = zero

              C5C6S6Lim%Lim_z(1) = C5Lim%Lim_z(1)
              C5C6S6Lim%Lim_z(2) = E6/Ek
              C5C6S6Lim%Lim_sij(i,5) = C5Lim%Lim_sij(i,5)
              C5C6S6Lim%Lim_sij(k,6) = 4*E6*Ek*x4

              C5C6S6Lim%wgt = C5S6Lim%wgt

            endif
          endif
        endif

    !! -------------------------------- C6k --------------------------------- !!
        if(present(C5C6Lim)) then

          C5C6Lim%PartFrac = C5Lim%PartFrac
          ns(:,6) = ns(:,k)

          Jf = 2*scr(pv,ns(:,k))
          Ek = (Q2 - 2*E6*scr(pv,ns(:,k)))/Jf
          pk = (Ek + E6) * ns(:,k)

          pl = pv - pk
          ns(:,l) = pl/pl(1)

          !-- Check valid point according to energy of final-state particles
          if(Ek < zero .or. pl(1) < zero) then

            C5C6Lim%flag = .true.

          else

            C5C6Lim%npart = npart_lo

            C5C6Lim%AmpMom(:,1) = p1
            C5C6Lim%AmpMom(:,2) = p2
            C5C6Lim%AmpMom(:,k) = pk
            C5C6Lim%AmpMom(:,l) = pl

            C5C6Lim%Lim_etaij(1,2) = one
            C5C6Lim%Lim_etaij(i,5) = C5Lim%Lim_etaij(i,5)
            C5C6Lim%Lim_etaij(j,5) = C5Lim%Lim_etaij(j,5)
            C5C6Lim%Lim_etaij(1,6) = (one - ns(4,k))*half
            C5C6Lim%Lim_etaij(2,6) = (one + ns(4,k))*half
            C5C6Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,k)))*half

            if(present(opt_etas)) call fill_etas(C5C6Lim,ns,opt_etas)

            C5C6Lim%Lim_etaij(k,6) = zero

            C5C6Lim%Lim_z(1) = C5Lim%Lim_z(1)
            C5C6Lim%Lim_z(2) = E6/(E6+Ek)

            C5C6Lim%Lim_sij(i,5) = C5Lim%Lim_sij(i,5)
            C5C6Lim%Lim_sij(k,6) = 4*E6*Ek*x4

            C5C6Lim%wgt = one/(2*spart)       &
              * x1 * spart**2/mv2             &
              * (x2/pi) * (Ek/Jf) * (Emax**2) &
              * jac

          endif
        endif
      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                      STRUCTURE 4: C5S5 + S6 + C6                      !!
    !!-----------------------------------------------------------------------!!

    if(present(C5S5Lim)) then
    !! ------------------------------- C5S5 -------------------------------- !!

      ns(:,6) = [one,n6aux(1),n6aux(2),n6aux(3)]

      C5S5Lim%Lim_etaij(1,2) = one
      C5S5Lim%Lim_etaij(i,5) = zero
      C5S5Lim%Lim_etaij(j,5) = one

      ns(:,5) = [one,zero,zero,sgn]

    !! ------------------------------ Hard 6 ------------------------------- !!

      C5S5Lim%npart = npart_nlo

      C5S5Lim%PartFrac = sqrt(tau)*[exp(+ylab),exp(-ylab)]
      sqrts = mv
      spart = mv2
      Emax = sqrts*half

      C5S5Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)

      p1 = Emax*ns(:,1)
      p2 = Emax*ns(:,2)

      pv = p1 + p2
      Q2 = spart !-- scr(pv,pv)

      E6 = x2*Emax
      p6 = E6*ns(:,6)

      Jf = 2*scr(pv,ns(:,k)) - 4*Emax*x2*x4
      Ek = (Q2 - 2*E6*scr(pv,ns(:,6)))/Jf

      pk = Ek * ns(:,k)
      pl = pv - pk - p6
      ns(:,l) = pl/pl(1)

      !-- Check valid point according to energy of final-state particles
      if(Ek < zero .or. pl(1) < zero) then

        C5S5Lim%flag = .true.

      else

        C5S5Lim%AmpMom(:,1) = p1
        C5S5Lim%AmpMom(:,2) = p2
        C5S5Lim%AmpMom(:,5) = p6
        C5S5Lim%AmpMom(:,k) = pk
        C5S5Lim%AmpMom(:,l) = pl

        C5S5Lim%Lim_etaij(1,6) = (one - ns(4,6))*half
        C5S5Lim%Lim_etaij(2,6) = (one + ns(4,6))*half
        C5S5Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,6)))*half

        if(present(opt_etas)) call fill_etas(C5S5Lim,ns,opt_etas)

        C5S5Lim%Lim_etaij(k,6) = x4

        C5S5Lim%Lim_sij(k,5) = spart*x1*x3
        C5S5Lim%Lim_z(1) = x1

        C5S5Lim%wgt = one/(2*spart)       &
          * x1 * spart**2/mv2             &
          * (x2/pi) * (Ek/Jf) * (Emax**2) &
          * jac

      endif
    !! -------------------------------- S6 --------------------------------- !!
      if(present(C5S5S6Lim)) then

        C5S5S6Lim%PartFrac = C5S5Lim%PartFrac

        Jf = 2*scr(pv,ns(:,k))
        Ek = Q2/Jf

        pk = Ek * ns(:,k)
        pl = pv - pk
        ns(:,l) = pl/pl(1)

        !-- Check valid point according to energy of final-state particles
        if(Ek < zero .or. pl(1) < zero) then

          C5S5S6Lim%flag = .true.

        else

          C5S5S6Lim%npart = npart_lo

          C5S5S6Lim%AmpMom(:,1) = p1
          C5S5S6Lim%AmpMom(:,2) = p2
          C5S5S6Lim%AmpMom(:,k) = pk
          C5S5S6Lim%AmpMom(:,l) = pl

          C5S5S6Lim%Lim_KinInv(2) = E6

          C5S5S6Lim%Lim_etaij(1,2) = one
          C5S5S6Lim%Lim_etaij(i,5) = C5S5Lim%Lim_etaij(i,5)
          C5S5S6Lim%Lim_etaij(j,5) = C5S5Lim%Lim_etaij(j,5)
          C5S5S6Lim%Lim_etaij(1,6) = (one - ns(4,6))*half
          C5S5S6Lim%Lim_etaij(2,6) = (one + ns(4,6))*half
          C5S5S6Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,6)))*half
          C5S5S6Lim%Lim_etaij(6,5) = C5S5S6Lim%Lim_etaij(5,6)

          if(present(opt_etas)) call fill_etas(C5S5S6Lim,ns,opt_etas)

          C5S5S6Lim%Lim_etaij(k,6) = x4

          C5S5S6Lim%Lim_sij(i,5) = C5S5Lim%Lim_sij(i,5)

          C5S5S6Lim%Lim_z(1) = x1

          C5S5S6Lim%wgt = one/(2*spart)     &
            * x1 * spart**2/mv2             &
            * (x2/pi) * (Ek/Jf) * (Emax**2) &
            * jac

    !! ------------------------------ S6 + C6k ----------------------------- !!
          if(present(C5C6S5S6Lim)) then

            C5C6S5S6Lim%PartFrac = C5S5Lim%PartFrac
            ns(:,6) = ns(:,k)

            C5C6S5S6Lim%Lim_etaij(1,2) = one
            C5C6S5S6Lim%Lim_etaij(i,5) = C5S5Lim%Lim_etaij(i,5)
            C5C6S5S6Lim%Lim_etaij(j,5) = C5S5Lim%Lim_etaij(j,5)
            C5C6S5S6Lim%Lim_etaij(1,6) = (one - ns(4,k))*half
            C5C6S5S6Lim%Lim_etaij(2,6) = (one + ns(4,k))*half
            C5C6S5S6Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,k)))*half

            if(present(opt_etas)) call fill_etas(C5C6S5S6Lim,ns,opt_etas)

            C5C6S5S6Lim%Lim_etaij(k,6) = zero

            C5C6S5S6Lim%Lim_sij(i,5) = C5S5Lim%Lim_sij(i,5)
            C5C6S5S6Lim%Lim_z(2) = E6/Ek
            C5C6S5S6Lim%Lim_sij(k,6) = 4*E6*Ek*x4

            C5C6S5S6Lim%wgt = C5S5S6Lim%wgt

          endif
        endif
      endif

    !! -------------------------------- C6k --------------------------------- !!
      if(present(C5C6S5Lim)) then

        C5C6S5Lim%PartFrac = C5S5Lim%PartFrac
        ns(:,6) = ns(:,k)

        Jf = 2*scr(pv,ns(:,k))
        Ek = (Q2 - 2*E6*scr(pv,ns(:,k)))/Jf
        pk = (Ek + E6) * ns(:,k)

        pl = pv - pk
        ns(:,l) = pl/pl(1)

        !-- Check valid point according to energy of final-state particles
        if(Ek < zero .or. pl(1) < zero) then

          C5C6S5Lim%flag = .true.

        else

          C5C6S5Lim%npart = npart_lo

          C5C6S5Lim%AmpMom(:,1) = p1
          C5C6S5Lim%AmpMom(:,2) = p2
          C5C6S5Lim%AmpMom(:,k) = pk
          C5C6S5Lim%AmpMom(:,l) = pl

          C5C6S5Lim%Lim_etaij(1,2) = one
          C5C6S5Lim%Lim_etaij(i,5) = C5S5Lim%Lim_etaij(i,5)
          C5C6S5Lim%Lim_etaij(j,5) = C5S5Lim%Lim_etaij(j,5)
          C5C6S5Lim%Lim_etaij(1,6) = (one - ns(4,k))*half
          C5C6S5Lim%Lim_etaij(2,6) = (one + ns(4,k))*half
          C5C6S5Lim%Lim_etaij(5,6) = (one - dot_product(ns(2:4,5),ns(2:4,k)))*half

          if(present(opt_etas)) call fill_etas(C5C6S5Lim,ns,opt_etas)

          C5C6S5Lim%Lim_etaij(k,6) = zero

          C5C6S5Lim%Lim_sij(i,5) = C5S5Lim%Lim_sij(i,5)
          C5C6S5Lim%Lim_z(2) = E6/(E6+Ek)
          C5C6S5Lim%Lim_sij(k,6) = 4*E6*Ek*x4

          C5C6S5Lim%Lim_z(1) = x1

          C5C6S5Lim%wgt = one/(2*spart)     &
            * x1 * spart**2/mv2             &
            * (x2/pi) * (Ek/Jf) * (Emax**2) &
            * jac

        endif
      endif
    endif

  end subroutine kinematics_nnlo_if_5i6k

end module mod_kinematics_nnlo_dc_if

!-- Triple-collinear phase space, ISR, energy unordered, angular ordering (a,b,c,d).
!-- Damping factor not included
module mod_kinematics_nnlo_tc_is_eu_abcd
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
  integer, parameter :: xLAMBDA = kNNLO_min + 5
  integer, parameter :: xSGN56  = kNNLO_min + 9
  !
  integer, parameter :: xLMIN  = kNNLO_min + 6
  integer, parameter :: xLMAX  = kNNLO_min + 8

  public :: kinematics_nnlo_abcd_5i6iac
  public :: kinematics_nnlo_abcd_5i6ibd

  private

contains

  !! i = {1,2} collinear direction
  !! j = {2,1} opposite direction
  !! opt_etas: indices k for eta_kl, l=1,...,n. 
  !! eta_ij, where i,j refer to initial-state particles and 
  !! radiated particles, are always computed by default
  subroutine kinematics_nnlo_abcd_5i6iac(yr, i, j, sec, &
    HardProc, S5Lim, S6Lim, TCLim, TCC5Lim, TCC6Lim,    &
    TCS5Lim, TCC6S5Lim, TCS6Lim, TCC5S6Lim, TCC5S5Lim,  &
    C5Lim, C6Lim, C5S6Lim, C6S5Lim, C5S5Lim, &
    opt_etas)
    real(dp),        intent(in) :: yr(:)
    integer,         intent(in) :: i, j
    integer,         intent(in) :: sec
    type(KinConfig), intent(inout) :: HardProc
    type(KinConfig), optional, intent(inout) :: S5Lim, S6Lim, TCLim, TCC5Lim, TCC6Lim
    type(KinConfig), optional, intent(inout) :: TCS5Lim, TCC6S5Lim, TCS6Lim, TCC5S6Lim
    type(KinConfig), optional, intent(inout) :: TCC5S5Lim
    type(KinConfig), optional, intent(inout) :: C5Lim, C6Lim, C5S6Lim, C6S5Lim, C5S5Lim
    integer, optional, intent(in) :: opt_etas(:)
    !--
    real(dp) :: mv, mv2, tau, ylab, jac
    real(dp) :: p1(4), p2(4), p5(4), p6(4), pv(4)
    real(dp) :: sgn, sgn56, kallenF
    real(dp) :: x1, x2, x3, x4, x5, x6
    real(dp) :: cos5i, sin5i, cos6i, sin6i, lambda, phi5, NN
    real(dp) :: cos_phi5, sin_phi5, cos_phi56, sin_phi56
    real(dp) :: xi1, xi2, spart, sqrts
    real(dp) :: nlept(3,2)
    real(dp), dimension(4) :: a, b, t, e3
    real(dp) :: ns(4,6)
    real(dp) :: eta5i, eta6i, eta5j, eta6j
    real(dp) :: eta51, eta61, eta52, eta62, eta65

    call get_tauy(yr(1:tauy_max),tau,ylab,jac)
    mv2 = tau*sh
    mv  = sqrt(mv2)

    x1 = yr(xE5)     ! E_5
    x2 = yr(xE6)     ! E_6
    x3 = yr(xRHO5)   ! \theta_{5i}
    x4 = yr(xRHO6)   ! \theta_{6j}
    x5 = yr(xPHI5)   ! \phi_{5}
    x6 = yr(xLAMBDA) ! related to \lambda

    !! Choose sign of \cos\theta depending on is emitter
    !! i = 1 or i = 2
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

    call initialize_config(HardProc)

    !! Azimuthal variables
    phi5 = twopi*x5
    cos_phi5 = cos(phi5)
    sin_phi5 = sin(phi5)

    !! Polar variables
    if(sec .eq. 1) then       !-- sec = 1 -> a
      eta5i = x3
      eta6i = x3*x4/2
    else if (sec .eq. 3) then !-- sec = 3 -> c
      eta5i = x3*x4/2
      eta6i = x3
    else
      print *, 'wrong sector'
      eta5i = -one
      eta6i = -one
      stop
    end if
    eta5j = one - eta5i
    eta6j = one - eta6i

    cos5i = sgn*(one - 2*eta5i)
    sin5i = two*sqrt(eta5i*(one-eta5i))

    cos6i = sgn*(one - 2*eta6i)
    sin6i = two*sqrt(eta6i*(one-eta6i))

    if(yr(xSGN56) .ge. half) then
      sgn56 = +one
    else
      sgn56 = -one
    end if

    t  = [one,      zero,    zero,zero]
    a  = [zero,-sin_phi5,cos_phi5,zero]
    b  = [zero, cos_phi5,sin_phi5,zero]
    e3 = [zero,     zero,    zero, one]

    lambda = sin(x6*pi/2)**2

    NN = NF1(x3,x4*half,lambda)

    sin_phi56 = sgn56*two*sqrt(lambda*(one-lambda))*(one-half*x4)/NN

    if((one-half*x4)**2/NN - (one+half*x4*(one-two*x3)).gt.0) then
      cos_phi56 = -sqrt(abs(one-sin_phi56**2))
    else
      cos_phi56 =  sqrt(abs(one-sin_phi56**2))
    endif

    ns(:,1) = t + e3
    ns(:,2) = t - e3
    ns(:,5) = t + cos5i*e3 + sin5i*b
    ns(:,6) = t + cos6i*e3 + sin6i*(cos_phi56*b + sin_phi56*a)

    eta65 = x3*(one-half*x4)**2/NN

    HardProc%Lim_etaij(1,2) = one
    HardProc%Lim_etaij(i,5) = eta5i
    HardProc%Lim_etaij(i,6) = eta6i
    HardProc%Lim_etaij(j,5) = eta5j
    HardProc%Lim_etaij(j,6) = eta6j
    HardProc%Lim_etaij(5,6) = eta65

    eta51 = HardProc%Lim_etaij(1,5)
    eta52 = HardProc%Lim_etaij(2,5)
    eta61 = HardProc%Lim_etaij(1,6)
    eta62 = HardProc%Lim_etaij(2,6)

    !! Compute partonic fractions {xi1, xi2}, partonic s and its square root
    !! and check if point is physical
    call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
      spart,sqrts,xi1,xi2,HardProc%flag)

    if (.not. HardProc%flag) then

      HardProc%PartFrac = (/xi1,xi2/)

      p1 = half*sqrts*ns(:,1)
      p2 = half*sqrts*ns(:,2)

      p5 = x1*half*sqrts*ns(:,5)
      p6 = x2*half*sqrts*ns(:,6)

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
           * one/two/spart  &
           * (spart/4)**2   &
           * 4*x1*x2*x3*(one-x4/2)/NN/2 &
           * spart/mv2      &
           * jac

    endif

    !!-----------------------------------------------------------------------!!
    !!                         STRUCTURE 2: S5 FLM.                          !!
    !!-----------------------------------------------------------------------!!

    if(present(S5Lim)) then

      call initialize_config(S5Lim)

      !! Angles same as in hard configuration, thus not recomputed
      call get_part_s_xi_nlo(mv2,ylab,x2,eta61,eta62,spart,sqrts,xi1,xi2,S5Lim%flag)

      if (.not. S5Lim%flag) then

        S5Lim%PartFrac = (/xi1,xi2/)

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p6 = x2*half*sqrts*ns(:,6)

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
             * one/two/spart  &
             * (spart/4)**2   &
             * 4*x1*x2*x3*(one-x4/2)/NN/2 &
             * spart/mv2      &
             * jac

      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                         STRUCTURE 3: S6 FLM.                          !!
    !!-----------------------------------------------------------------------!!

    if(present(S6Lim)) then

      call initialize_config(S6Lim)

      !! Angles same as in hard configuration, thus not recomputed
      call get_part_s_xi_nlo(mv2,ylab,x1,eta51,eta52,spart,sqrts,xi1,xi2,S6Lim%flag)

      if (.not. S6Lim%flag) then

        S6Lim%PartFrac = (/xi1,xi2/)

        p1 = half*sqrts*ns(:,1) 
        p2 = half*sqrts*ns(:,2)

        p5 = x1*half*sqrts*ns(:,5)

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

        if(present(opt_etas)) call fill_etas(S6Lim,ns,opt_etas)

        S6Lim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * 4*x1*x2*x3*(one-x4/2)/NN/2 &
             * spart/mv2      &
             * jac

      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                       STRUCTURE 4: TC FLM.                            !!
    !!-----------------------------------------------------------------------!!

    if(present(TCLim)) then

      call initialize_config(TCLim)

      !! x3 -> 0
      eta65 = zero
      TCLim%Lim_etaij(1,2) = one
      TCLim%Lim_etaij(i,5) = zero
      TCLim%Lim_etaij(i,6) = zero
      TCLim%Lim_etaij(j,5) = one
      TCLim%Lim_etaij(j,6) = one
      TCLim%Lim_etaij(5,6) = eta65

      eta51 = TCLim%Lim_etaij(1,5)
      eta52 = TCLim%Lim_etaij(2,5)
      eta61 = TCLim%Lim_etaij(1,6)
      eta62 = TCLim%Lim_etaij(2,6)

      call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
        spart,sqrts,xi1,xi2,TCLim%flag)

      if (.not. TCLim%flag) then

      !! ---------------------------------- TC ---------------------------- !!

        TCLim%PartFrac = (/xi1,xi2/)

        ns(:,5) = t + sgn*e3
        ns(:,6) = t + sgn*e3

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = x1*half*sqrts*ns(:,5)
        p6 = x2*half*sqrts*ns(:,6)

        if(i .eq. 1) then
          p1 = p1 - p5 - p6
        elseif(i .eq. 2) then
          p2 = p2 - p5 - p6
        endif

        pv = p1 + p2

        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, TCLim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        TCLim%AmpMom(:,1) = p1
        TCLim%AmpMom(:,2) = p2

        TCLim%Lim_z(1) = x1/(x1 + x2 - one)
        TCLim%Lim_z(2) = x2/(x1 + x2 - one)
        TCLim%Lim_z(3) = one/(one - x1 - x2)

        !! set sij in the limit. sij = 2*pi.pj.
        !! Only those that are needed are computed
        TCLim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
        TCLim%Lim_sij(i,6) = spart*x2*HardProc%Lim_etaij(i,6)
        TCLim%Lim_sij(5,6) = spart*x1*x2*HardProc%Lim_etaij(5,6)

        if(present(opt_etas)) call fill_etas(TCLim,ns,opt_etas)

        TCLim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * 4*x1*x2*x3*(one-x4/2)/NN/2 &
             * spart/mv2      &
             * jac

      !! ----------------------------- TC + C5 ----------------------------- !!
        if(present(TCC5Lim)) then

#if(_withchecks == 1)
          if(sec == 1) then
            print *,'In mod_kinematics_nnlo_tc_is_eu_abcd/kinematics_nnlo_abcd_5i6iac'
            print *,'Calling TCC5 limit in sector a. Aborting'
            stop
          endif
#endif

          call initialize_config(TCC5Lim)

          TCC5Lim%Lim_z(1) = x1/(x1 - one)
          TCC5Lim%Lim_z(2) = x2/(x1 + x2 - one)

          !! si5 + si6 - s56 in the x4 -> 0 limit
          TCC5Lim%Lim_KinInv(1) = (one-x1)*spart*x2*x3

          TCC5Lim%wgt = one/8.0_dp/pi * kallenF &
               * one/two/spart  &
               * (spart/4)**2   &
               * 4*x1*x2*x3/NF1(x3,zero,lambda)/2 &
               * spart/mv2      &
               * jac

        endif

      !! ----------------------------- TC + C6 ----------------------------- !!
        if(present(TCC6Lim)) then

#if(_withchecks == 1)
          if(sec == 3) then
            print *,'In mod_kinematics_nnlo_tc_is_eu_abcd/kinematics_nnlo_abcd_5i6iac'
            print *,'Calling TCC6 limit in sector c. Aborting'
            stop
          endif
#endif

          call initialize_config(TCC6Lim)

          TCC6Lim%Lim_z(1) = x1/(x1 + x2 - one)
          TCC6Lim%Lim_z(2) = x2/(x2 - one)

          !! si5 + si6 - s56 in the x4 -> 0 limit
          TCC6Lim%Lim_KinInv(1) = (one-x2)*spart*x1*x3

          TCC6Lim%wgt = one/8.0_dp/pi * kallenF &
               * one/two/spart  &
               * (spart/4)**2   &
               * 4*x1*x2*x3/NF1(x3,zero,lambda)/2 &
               * spart/mv2      &
               * jac

        endif

      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                     STRUCTURE 5: TC + S5 FLM.                         !!
    !!-----------------------------------------------------------------------!!

    if(present(TCS5Lim)) then

      call initialize_config(TCS5Lim)

      !! x3 -> 0
      eta65 = zero
      TCS5Lim%Lim_etaij(1,2) = one
      TCS5Lim%Lim_etaij(i,5) = zero
      TCS5Lim%Lim_etaij(i,6) = zero
      TCS5Lim%Lim_etaij(j,5) = one
      TCS5Lim%Lim_etaij(j,6) = one
      TCS5Lim%Lim_etaij(5,6) = eta65

      eta51 = TCS5Lim%Lim_etaij(1,5)
      eta52 = TCS5Lim%Lim_etaij(2,5)
      eta61 = TCS5Lim%Lim_etaij(1,6)
      eta62 = TCS5Lim%Lim_etaij(2,6)

      call get_part_s_xi_nlo(mv2,ylab,x2,eta61,eta62,spart,sqrts,xi1,xi2,TCS5Lim%flag)

      !! ------------------------------- TC + S5 --------------------------- !!

      if (.not. TCS5Lim%flag) then

        TCS5Lim%PartFrac = (/xi1,xi2/)

        ns(:,5) = t + sgn*e3
        ns(:,6) = t + sgn*e3

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = x1*half*sqrts*ns(:,5)
        p6 = x2*half*sqrts*ns(:,6)

        if(i .eq. 1) then
          p1 = p1 - p6
        elseif(i .eq. 2) then
          p2 = p2 - p6
        endif

        pv = p1 + p2

        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, TCS5Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        TCS5Lim%AmpMom(:,1) = p1
        TCS5Lim%AmpMom(:,2) = p2

        !!TODO: adopt a similar form in TCS6
        TCS5Lim%Lim_z(1) = x1/(x2 - one)
        TCS5Lim%Lim_z(2) = x2/(x2 - one)
        TCS5Lim%Lim_z(3) = one/(one - x2)

        !! set sij in the limit. sij = 2*pi.pj.
        !! Only those that are needed are computed
        TCS5Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
        TCS5Lim%Lim_sij(i,6) = spart*x2*HardProc%Lim_etaij(i,6)
        TCS5Lim%Lim_sij(5,6) = spart*x1*x2*HardProc%Lim_etaij(5,6)

        if(present(opt_etas)) call fill_etas(TCS5Lim,ns,opt_etas)

        TCS5Lim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * 4*x1*x2*x3*(one-x4/2)/NN/2 &
             * spart/mv2      &
             * jac

    !! --------------------------- TC + S5 + C5 -------------------------- !!
        if(present(TCC5S5Lim)) then

#if(_withchecks == 1)
          if(sec == 1) then
            print *,'In mod_kinematics_nnlo_tc_is_eu_abcd/kinematics_nnlo_abcd_5i6iac'
            print *,'Calling TCC5S5 limit in sector a. Aborting'
            stop
          endif
#endif

          call initialize_config(TCC5S5Lim)

          TCC5S5Lim%Lim_z(1) = x1
          TCC5S5Lim%Lim_z(2) = one/(one-x2)

          TCC5S5Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
          TCC5S5Lim%Lim_sij(i,6) = spart*x2*HardProc%Lim_etaij(i,6)

          TCC5S5Lim%wgt = one/8.0_dp/pi * kallenF &
               * one/two/spart  &
               * (spart/4)**2   &
               * 4*x1*x2*x3/NF1(x3,zero,lambda)/2 &
               * spart/mv2      &
               * jac

        endif

    !! --------------------------- TC + S5 + C6 -------------------------- !!
        if(present(TCC6S5Lim)) then

#if(_withchecks == 1)
          if(sec == 3) then
            print *,'In mod_kinematics_nnlo_tc_is_eu_abcd/kinematics_nnlo_abcd_5i6iac'
            print *,'Calling TCC6S5 limit in sector c. Aborting'
            stop
          endif
#endif

          call initialize_config(TCC6S5Lim)

          !! s56 in the x4 -> 0 limit
          TCC6S5Lim%Lim_sij(5,6) = spart*x1*x2*x3

          TCC6S5Lim%wgt = one/8.0_dp/pi * kallenF &
               * one/two/spart  &
               * (spart/4)**2   &
               * 4*x1*x2*x3/NF1(x3,zero,lambda)/2 &
               * spart/mv2      &
               * jac

        endif
      endif
    endif

    !!----------------------------------------------------------------------!!
    !!                    STRUCTURE 6: TC + S6 FLM.                         !!
    !!----------------------------------------------------------------------!!

    if(present(TCS6Lim)) then

      call initialize_config(TCS6Lim)

      !! x3 -> 0
      eta65 = zero
      TCS6Lim%Lim_etaij(1,2) = one
      TCS6Lim%Lim_etaij(i,5) = zero
      TCS6Lim%Lim_etaij(i,6) = zero
      TCS6Lim%Lim_etaij(j,5) = one
      TCS6Lim%Lim_etaij(j,6) = one
      TCS6Lim%Lim_etaij(5,6) = eta65

      eta51 = TCS6Lim%Lim_etaij(1,5)
      eta52 = TCS6Lim%Lim_etaij(2,5)
      eta61 = TCS6Lim%Lim_etaij(1,6)
      eta62 = TCS6Lim%Lim_etaij(2,6)

      call get_part_s_xi_nlo(mv2,ylab,x1,eta51,eta52,spart,sqrts,xi1,xi2,TCS6Lim%flag)

      !! ------------------------------- TC + S6 --------------------------- !!

      if (.not. TCS6Lim%flag) then

        TCS6Lim%PartFrac = (/xi1,xi2/)
        
        ns(:,5) = t + sgn*e3
        ns(:,6) = t + sgn*e3
        
        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = x1*half*sqrts*ns(:,5)
        p6 = x2*half*sqrts*ns(:,6)
        
        if(i .eq. 1) then
          p1 = p1 - p5
        elseif(i .eq. 2) then
          p2 = p2 - p5
        endif
        
        pv = p1 + p2
        
        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, TCS6Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)
        
        TCS6Lim%AmpMom(:,1) = p1
        TCS6Lim%AmpMom(:,2) = p2
        
        TCS6Lim%Lim_z(1) = x2/x1
        TCS6Lim%Lim_z(3) = one/(one-x1)
        
        !! set sij in the limit. sij = 2*pi.pj.
        !! Only those that are needed are computed
        TCS6Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
        TCS6Lim%Lim_sij(i,6) = spart*x2*HardProc%Lim_etaij(i,6)
        TCS6Lim%Lim_sij(5,6) = spart*x1*x2*HardProc%Lim_etaij(5,6)
        
        if(present(opt_etas)) call fill_etas(TCS6Lim,ns,opt_etas)
        
        TCS6Lim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * 4*x1*x2*x3*(one-x4/2)/NN/2 &
             * spart/mv2      &
             * jac

    !! --------------------------- TC + S6 + C5 -------------------------- !!
        if(present(TCC5S6Lim)) then

#if(_withchecks == 1)
          if(sec == 1) then
            print *,'In mod_kinematics_nnlo_tc_is_eu_abcd/kinematics_nnlo_abcd_5i6iac'
            print *,'Calling TCC5S6 limit in sector a. Aborting'
            stop
          endif
#endif

          call initialize_config(TCC5S6Lim)

          !! s56 in the x4 -> 0 limit
          TCC5S6Lim%Lim_sij(5,6) = spart*x1*x2*x3

          TCC5S6Lim%wgt = one/8.0_dp/pi * kallenF &
               * one/two/spart  &
               * (spart/4)**2   &
               * 4*x1*x2*x3/NF1(x3,zero,lambda)/2 &
               * spart/mv2      &
               * jac

        endif
      endif
    endif

    !!----------------------------------------------------------------------!!
    !!                       STRUCTURE 7: C5 FLM.                            !!
    !!----------------------------------------------------------------------!!

    !! x4 -> 0. Only sector 'c' relevant in this case
    if(present(C5Lim)) then

#if(_withchecks == 1)
      if(sec == 1) then
        print *,'In mod_kinematics_nnlo_tc_is_eu_abcd/kinematics_nnlo_abcd_5i6iac'
        print *,'Calling C5 limit in sector a. Aborting'
        stop
      endif
#endif

      call initialize_config(C5Lim)

      NN = NF1(x3,zero,lambda)
      eta65 = x3/NN

      C5Lim%Lim_etaij(1,2) = one
      C5Lim%Lim_etaij(i,5) = zero
      C5Lim%Lim_etaij(i,6) = x3
      C5Lim%Lim_etaij(j,5) = one
      C5Lim%Lim_etaij(j,6) = one - x3
      C5Lim%Lim_etaij(5,6) = eta65

      eta51 = C5Lim%Lim_etaij(1,5)
      eta52 = C5Lim%Lim_etaij(2,5)
      eta61 = C5Lim%Lim_etaij(1,6)
      eta62 = C5Lim%Lim_etaij(2,6)

      call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
        spart,sqrts,xi1,xi2,C5Lim%flag)

      if (.not. C5Lim%flag) then 

        C5Lim%PartFrac = (/xi1,xi2/)

        sin_phi56 = sgn56*two*sqrt(lambda*(one-lambda))/NN

        if(one-two*lambda .ge. zero) then
          cos_phi56 = -sqrt(abs(one-sin_phi56**2))
        else
          cos_phi56 =  sqrt(abs(one-sin_phi56**2))
        endif

        ns(:,5) = t + sgn*e3
        ns(:,6) = t + cos6i*e3 + sin6i*(cos_phi56*b + sin_phi56*a)

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = x1*half*sqrts*ns(:,5)
        p6 = x2*half*sqrts*ns(:,6)

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

        C5Lim%Lim_z(1) = one/(one-x1)

        !! set sij in the limit. sij = 2*pi.pj.
        !! Only those that are needed are computed
        C5Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
        C5Lim%Lim_sij(j,5) = spart*x1*HardProc%Lim_etaij(j,5)

        if(present(opt_etas)) call fill_etas(C5Lim,ns,opt_etas)

        !! x4 -> 0
        C5Lim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * 4*x1*x2*x3/NF1(x3,zero,lambda)/2 &
             * spart/mv2      &
             * jac

      endif
    endif

    !!----------------------------------------------------------------------!!
    !!                       STRUCTURE 8: C6 FLM.                            !!
    !!----------------------------------------------------------------------!!

    !! x4 -> 0. Only sector 'a' relevant in this case
    if(present(C6Lim)) then

#if(_withchecks == 1)
      if(sec == 3) then
        print *,'In mod_kinematics_nnlo_tc_is_eu_abcd/kinematics_nnlo_abcd_5i6iac'
        print *,'Calling C6 limit in sector c. Aborting'
        stop
      endif
#endif

      call initialize_config(C6Lim)

      NN = NF1(x3,zero,lambda)
      eta65 = x3/NN

      C6Lim%Lim_etaij(1,2) = one
      C6Lim%Lim_etaij(i,5) = x3
      C6Lim%Lim_etaij(i,6) = zero
      C6Lim%Lim_etaij(j,5) = one - x3
      C6Lim%Lim_etaij(j,6) = one
      C6Lim%Lim_etaij(5,6) = eta65

      eta51 = C6Lim%Lim_etaij(1,5)
      eta52 = C6Lim%Lim_etaij(2,5)
      eta61 = C6Lim%Lim_etaij(1,6)
      eta62 = C6Lim%Lim_etaij(2,6)

      call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
        spart,sqrts,xi1,xi2,C6Lim%flag)

      if (.not. C6Lim%flag) then

        C6Lim%PartFrac = (/xi1,xi2/)

        sin_phi56 = sgn56*two*sqrt(lambda*(one-lambda))/NN

        if(one-two*lambda .ge. zero) then
          cos_phi56 = -sqrt(abs(one-sin_phi56**2))
        else
          cos_phi56 =  sqrt(abs(one-sin_phi56**2))
        endif

        ns(:,5) = t + cos5i*e3 + sin5i*b
        ns(:,6) = t + sgn*e3

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = x1*half*sqrts*ns(:,5)
        p6 = x2*half*sqrts*ns(:,6)

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

        C6Lim%Lim_z(2) = one/(one-x2)

        !! set sij in the limit. sij = 2*pi.pj.
        !! Only those that are needed are computed
        C6Lim%Lim_sij(i,6) = spart*x2*HardProc%Lim_etaij(i,6)
        C6Lim%Lim_sij(j,6) = spart*x2*HardProc%Lim_etaij(j,6)

        if(present(opt_etas)) call fill_etas(C6Lim,ns,opt_etas)

        !! x4 -> 0
        C6Lim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * 4*x1*x2*x3/NF1(x3,zero,lambda)/2 &
             * spart/mv2      &
             * jac

      endif
    endif

    !!----------------------------------------------------------------------!!
    !!                       STRUCTURE 9: S6 + C5 FLM.                       !!
    !!----------------------------------------------------------------------!!

    !! x4 -> 0. Only sector 'c' relevant in this case
    if(present(C5S6Lim)) then

#if(_withchecks == 1)
      if(sec == 1) then
        print *,'In mod_kinematics_nnlo_tc_is_eu_abcd/kinematics_nnlo_abcd_5i6iac'
        print *,'Calling C5S6 limit in sector a. Aborting'
        stop
      endif
#endif

      call initialize_config(C5S6Lim)

      NN = NF1(x3,zero,lambda)
      eta65 = x3/NN

      C5S6Lim%Lim_etaij(1,2) = one
      C5S6Lim%Lim_etaij(i,5) = zero
      C5S6Lim%Lim_etaij(i,6) = x3
      C5S6Lim%Lim_etaij(j,5) = one
      C5S6Lim%Lim_etaij(j,6) = one - x3
      C5S6Lim%Lim_etaij(5,6) = eta65

      eta51 = C5S6Lim%Lim_etaij(1,5)
      eta52 = C5S6Lim%Lim_etaij(2,5)
      eta61 = C5S6Lim%Lim_etaij(1,6)
      eta62 = C5S6Lim%Lim_etaij(2,6)

      call get_part_s_xi_nlo(mv2,ylab,x1,eta51,eta52,spart,sqrts,xi1,xi2,C5S6Lim%flag)

      !! ------------------------------- S6 + C5 --------------------------- !!

      if (.not. C5S6Lim%flag) then 

        C5S6Lim%PartFrac = (/xi1,xi2/)

        sin_phi56 = sgn56*two*sqrt(lambda*(one-lambda))/NN

        if(one-two*lambda .ge. zero) then
          cos_phi56 = -sqrt(abs(one-sin_phi56**2))
        else
          cos_phi56 =  sqrt(abs(one-sin_phi56**2))
        endif

        ns(:,5) = t + sgn*e3
        ns(:,6) = t + cos6i*e3 + sin6i*(cos_phi56*b + sin_phi56*a)

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = x1*half*sqrts*ns(:,5)
        p6 = x2*half*sqrts*ns(:,6)

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

        C5S6Lim%Lim_z(1) = one/(one-x1)

        !! set sij in the limit. sij = 2*pi.pj.
        !! Only those that are needed are computed
        C5S6Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)

        !! E6
        C5S6Lim%Lim_KinInv(1) = x2*half*sqrts
   
        if(present(opt_etas)) call fill_etas(C5S6Lim,ns,opt_etas)

        !! x4 -> 0
        C5S6Lim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * 4*x1*x2*x3/NF1(x3,zero,lambda)/2 &
             * spart/mv2      &
             * jac

      endif
    endif

    !!----------------------------------------------------------------------!!
    !!                      STRUCTURE 10: S5 + C6 FLM.                      !!
    !!----------------------------------------------------------------------!!

    !! x4 -> 0. Only sector 'a' relevant in this case
    if(present(C6S5Lim)) then

#if(_withchecks == 1)
      if(sec == 3) then
        print *,'In mod_kinematics_nnlo_tc_is_eu_abcd/kinematics_nnlo_abcd_5i6iac'
        print *,'Calling C6S5 limit in sector c. Aborting'
        stop
      endif
#endif

      call initialize_config(C6S5Lim)

      NN = NF1(x3,zero,lambda)
      eta65 = x3/NN

      C6S5Lim%Lim_etaij(1,2) = one
      C6S5Lim%Lim_etaij(i,5) = x3
      C6S5Lim%Lim_etaij(i,6) = zero
      C6S5Lim%Lim_etaij(j,5) = one - x3
      C6S5Lim%Lim_etaij(j,6) = one
      C6S5Lim%Lim_etaij(5,6) = eta65
      C6S5Lim%Lim_etaij(6,5) = eta65

      eta51 = C6S5Lim%Lim_etaij(1,5)
      eta52 = C6S5Lim%Lim_etaij(2,5)
      eta61 = C6S5Lim%Lim_etaij(1,6)
      eta62 = C6S5Lim%Lim_etaij(2,6)

      call get_part_s_xi_nlo(mv2,ylab,x2,eta61,eta62,spart,sqrts,xi1,xi2,C6S5Lim%flag)

      !! ------------------------------- S6 + C5 --------------------------- !!

      if (.not. C6S5Lim%flag) then 

        C6S5Lim%PartFrac = (/xi1,xi2/)

        sin_phi56 = sgn56*two*sqrt(lambda*(one-lambda))/NN

        if(one-two*lambda .ge. zero) then
          cos_phi56 = -sqrt(abs(one-sin_phi56**2))
        else
          cos_phi56 =  sqrt(abs(one-sin_phi56**2))
        endif

        ns(:,5) = t + cos5i*e3 + sin5i*b
        ns(:,6) = t + sgn*e3

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = x1*half*sqrts*ns(:,5)
        p6 = x2*half*sqrts*ns(:,6)

        if(i .eq. 1) then
          p1 = p1 - p6
        elseif(i .eq. 2) then
          p2 = p2 - p6
        endif

        pv = p1 + p2

        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C6S5Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        C6S5Lim%AmpMom(:,1) = p1
        C6S5Lim%AmpMom(:,2) = p2

        C6S5Lim%Lim_z(2) = one/(one-x2)

        !! set sij in the limit. sij = 2*pi.pj.
        !! Only those that are needed are computed
        C6S5Lim%Lim_sij(i,6) = spart*x2*HardProc%Lim_etaij(i,6)

        !! E5
        C6S5Lim%Lim_KinInv(1) = x1*half*sqrts

        if(present(opt_etas)) call fill_etas(C6S5Lim,ns,opt_etas)

        !! x4 -> 0
        C6S5Lim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * 4*x1*x2*x3/NF1(x3,zero,lambda)/2 &
             * spart/mv2      &
             * jac

      endif
    endif

    !!----------------------------------------------------------------------!!
    !!                      STRUCTURE 11: S5 + C5 FLM.                      !!
    !!----------------------------------------------------------------------!!

    !! x4 -> 0. Only sector 'c' relevant in this case
    if(present(C5S5Lim)) then

#if(_withchecks == 1)
      if(sec == 1) then
        print *,'In mod_kinematics_nnlo_tc_is_eu_abcd/kinematics_nnlo_abcd_5i6iac'
        print *,'Calling C5S5 limit in sector a. Aborting'
        stop
      endif
#endif

      call initialize_config(C5S5Lim)

      NN = NF1(x3,zero,lambda)
      eta65 = x3/NN

      C5S5Lim%Lim_etaij(1,2) = one
      C5S5Lim%Lim_etaij(i,5) = zero
      C5S5Lim%Lim_etaij(i,6) = x3
      C5S5Lim%Lim_etaij(j,5) = one
      C5S5Lim%Lim_etaij(j,6) = one - x3
      C5S5Lim%Lim_etaij(5,6) = eta65

      eta51 = C5S5Lim%Lim_etaij(1,5)
      eta52 = C5S5Lim%Lim_etaij(2,5)
      eta61 = C5S5Lim%Lim_etaij(1,6)
      eta62 = C5S5Lim%Lim_etaij(2,6)

      call get_part_s_xi_nlo(mv2,ylab,x2,eta61,eta62,spart,sqrts,xi1,xi2,C5S5Lim%flag)

      !! ------------------------------- S5 + C5 --------------------------- !!

      if (.not. C5S5Lim%flag) then 

        C5S5Lim%PartFrac = (/xi1,xi2/)

        sin_phi56 = sgn56*two*sqrt(lambda*(one-lambda))/NN

        if(one-two*lambda .ge. zero) then
          cos_phi56 = -sqrt(abs(one-sin_phi56**2))
        else
          cos_phi56 =  sqrt(abs(one-sin_phi56**2))
        endif

        ns(:,5) = t + sgn*e3
        ns(:,6) = t + cos6i*e3 + sin6i*(cos_phi56*b + sin_phi56*a)

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p6 = x2*half*sqrts*ns(:,6)

        pv = p1 + p2 - p6

        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C5S5Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        C5S5Lim%AmpMom(:,1) = p1
        C5S5Lim%AmpMom(:,2) = p2
        C5S5Lim%AmpMom(:,5) = p6

        C5S5Lim%Lim_z(1) = x1
        C5S5Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)

        if(present(opt_etas)) call fill_etas(C5S5Lim,ns,opt_etas)

        !! x4 -> 0
        C5S5Lim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * 4*x1*x2*x3/NF1(x3,zero,lambda)/2 &
             * spart/mv2      &
             * jac

      endif
    endif

  end subroutine kinematics_nnlo_abcd_5i6iac

  !!*************************************************************************!!

  !! i = {1,2} collinear direction
  !! j = {2,1} opposite direction
  !! opt_etas: indices k for eta_kl, l=1,...,n. 
  !! eta_ij, where i,j refer to initial-state particles and 
  !! radiated particles, are always computed by default
  subroutine kinematics_nnlo_abcd_5i6ibd(yr, i, j, sec, &
    HardProc, S5Lim, S6Lim, TCLim, TCC6Lim,  &
    TCS5Lim, TCC6S5Lim, TCS6Lim, TCC6S6Lim,  &
    C6Lim, C6S5Lim, C6S6Lim, opt_etas)
    real(dp),        intent(in) :: yr(:)
    integer,         intent(in) :: i, j
    integer,         intent(in) :: sec
    type(KinConfig), intent(inout) :: HardProc
    type(KinConfig), optional, intent(inout) :: S5Lim, S6Lim, TCLim, TCC6Lim
    type(KinConfig), optional, intent(inout) :: TCS5Lim, TCC6S5Lim, TCS6Lim, TCC6S6Lim
    type(KinConfig), optional, intent(inout) :: C6Lim, C6S5Lim, C6S6Lim
    integer, optional, intent(in) :: opt_etas(:)
    !--
    real(dp) :: mv, mv2, tau, ylab, jac
    real(dp) :: p1(4), p2(4), p5(4), p6(4), pv(4)
    real(dp) :: sgn, sgn56, kallenF
    real(dp) :: x1, x2, x3, x4, x5, x6
    real(dp) :: cos5i, sin5i, cos6i, sin6i, lambda, phi5, NN
    real(dp) :: cos_phi5, sin_phi5, cos_phi56, sin_phi56
    real(dp) :: xi1, xi2, spart, sqrts
    real(dp) :: nlept(3,2)
    real(dp), dimension(4) :: a, b, t, e3
    real(dp) :: ns(4,6), x4_half
    real(dp) :: eta5i, eta6i, eta5j, eta6j
    real(dp) :: eta51, eta61, eta52, eta62, eta65

    call get_tauy(yr(1:tauy_max),tau,ylab,jac)
    mv2 = tau*sh
    mv  = sqrt(mv2)

    x1 =   yr(xE5)     ! E_5
    x2 =   yr(xE6)     ! E_6
    x3 =   yr(xRHO5)   ! \theta_{5i}
    x4 =   yr(xRHO6)   ! \theta_{6j}
    x5 =   yr(xPHI5)   ! \phi_{5}
    x6 =   yr(xLAMBDA) ! related to \lambda

    x4_half = x4*half

    !! Choose sign of \cos\theta depending on is emitter
    !! i = 1 or i = 2
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

    call initialize_config(HardProc)

    !! Azimuthal variables
    phi5 = twopi*x5
    cos_phi5 = cos(phi5)
    sin_phi5 = sin(phi5)

    !! Polar variables
    if(sec .eq. 2) then       !-- sec = 2 -> b
      eta5i = x3
      eta6i = x3*(one-x4_half)
    else if (sec .eq. 4) then !-- sec = 4 -> d
      eta5i = x3*(one-x4_half)
      eta6i = x3
    else
      print *, 'wrong sector'
      eta5i = -one
      eta6i = -one
      stop
    end if
   
    eta5j = one - eta5i
    eta6j = one - eta6i

    cos5i = sgn*(one - 2*eta5i)
    sin5i = two*sqrt(eta5i*(one-eta5i))

    cos6i = sgn*(one - 2*eta6i)
    sin6i = two*sqrt(eta6i*(one-eta6i))

    if(yr(xSGN56) .ge. half) then
      sgn56 = +one
    else
      sgn56 = -one
    end if

    t  = [one,      zero,    zero,zero]
    a  = [zero,-sin_phi5,cos_phi5,zero]
    b  = [zero, cos_phi5,sin_phi5,zero]
    e3 = [zero,     zero,    zero, one]

    lambda = sin(x6*pi/2)**2

    NN = NF1(x3,one-x4_half,lambda)

    sin_phi56 = sgn56*two*sqrt(lambda*(one-lambda))*x4_half/NN

    if((x4_half**2/NN - (one + (one - 2*x3)*(one-x4_half))).gt.0) then
      cos_phi56 = -sqrt(abs(one-sin_phi56**2))
    else
      cos_phi56 =  sqrt(abs(one-sin_phi56**2))
    endif

    ns(:,1) = t + e3
    ns(:,2) = t - e3
    ns(:,5) = t + cos5i*e3 + sin5i*b
    ns(:,6) = t + cos6i*e3 + sin6i*(cos_phi56*b + sin_phi56*a)

    eta65 = x3*x4_half**2/NN

    HardProc%Lim_etaij(1,2) = one
    HardProc%Lim_etaij(i,5) = eta5i
    HardProc%Lim_etaij(i,6) = eta6i
    HardProc%Lim_etaij(j,5) = eta5j
    HardProc%Lim_etaij(j,6) = eta6j
    HardProc%Lim_etaij(5,6) = eta65

    eta51 = HardProc%Lim_etaij(1,5)
    eta52 = HardProc%Lim_etaij(2,5)
    eta61 = HardProc%Lim_etaij(1,6)
    eta62 = HardProc%Lim_etaij(2,6)

    !! Compute partonic fractions {xi1, xi2}, partonic s and its square root
    !! and check if point is physical
    call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
      spart,sqrts,xi1,xi2,HardProc%flag)

    if (.not. HardProc%flag) then

      HardProc%PartFrac = (/xi1,xi2/)

      p1 = half*sqrts*ns(:,1)
      p2 = half*sqrts*ns(:,2)

      p5 = x1*half*sqrts*ns(:,5)
      p6 = x2*half*sqrts*ns(:,6)

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
           * one/two/spart  &
           * (spart/4)**2   &
           * x1*x2*x3*x4/NF1(x3,one-x4_half,lambda) &
           * spart/mv2      &
           * jac

    endif

    !!-----------------------------------------------------------------------!!
    !!                         STRUCTURE 2: S5 FLM.                          !!
    !!-----------------------------------------------------------------------!!

    if(present(S5Lim)) then

      call initialize_config(S5Lim)

      !! Angles same as in hard configuration, thus not recomputed
      call get_part_s_xi_nlo(mv2,ylab,x2,eta61,eta62,spart,sqrts,xi1,xi2,S5Lim%flag)

      if (.not. S5Lim%flag) then

        S5Lim%PartFrac = (/xi1,xi2/)

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p6 = x2*half*sqrts*ns(:,6)

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
           * one/two/spart  &
           * (spart/4)**2   &
           * x1*x2*x3*x4/NF1(x3,one-x4_half,lambda) &
           * spart/mv2      &
           * jac

      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                         STRUCTURE 3: S6 FLM.                          !!
    !!-----------------------------------------------------------------------!!

    if(present(S6Lim)) then

      call initialize_config(S6Lim)

      !! Angles same as in hard configuration, thus not recomputed
      call get_part_s_xi_nlo(mv2,ylab,x1,eta51,eta52,spart,sqrts,xi1,xi2,S6Lim%flag)

      if (.not. S6Lim%flag) then

        S6Lim%PartFrac = (/xi1,xi2/)

        p1 = half*sqrts*ns(:,1) 
        p2 = half*sqrts*ns(:,2)

        p5 = x1*half*sqrts*ns(:,5)

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

        if(present(opt_etas)) call fill_etas(S6Lim,ns,opt_etas)

        S6Lim%wgt = one/8.0_dp/pi * kallenF &
           * one/two/spart  &
           * (spart/4)**2   &
           * x1*x2*x3*x4/NF1(x3,one-x4_half,lambda) &
           * spart/mv2      &
           * jac

      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                       STRUCTURE 4: TC FLM.                            !!
    !!-----------------------------------------------------------------------!!

    if(present(TCLim)) then

      call initialize_config(TCLim)

      !! x3 -> 0
      eta65 = zero
      TCLim%Lim_etaij(1,2) = one
      TCLim%Lim_etaij(i,5) = zero
      TCLim%Lim_etaij(i,6) = zero
      TCLim%Lim_etaij(j,5) = one
      TCLim%Lim_etaij(j,6) = one
      TCLim%Lim_etaij(5,6) = eta65

      eta51 = TCLim%Lim_etaij(1,5)
      eta52 = TCLim%Lim_etaij(2,5)
      eta61 = TCLim%Lim_etaij(1,6)
      eta62 = TCLim%Lim_etaij(2,6)

      call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
        spart,sqrts,xi1,xi2,TCLim%flag)

      if (.not. TCLim%flag) then

    !! ---------------------------------- TC ---------------------------- !!

        TCLim%PartFrac = (/xi1,xi2/)

        ns(:,5) = t + sgn*e3
        ns(:,6) = t + sgn*e3

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = x1*half*sqrts*ns(:,5)
        p6 = x2*half*sqrts*ns(:,6)

        if(i .eq. 1) then
          p1 = p1 - p5 - p6
        elseif(i .eq. 2) then
          p2 = p2 - p5 - p6
        endif

        pv = p1 + p2

        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, TCLim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        TCLim%AmpMom(:,1) = p1
        TCLim%AmpMom(:,2) = p2

        TCLim%Lim_z(1) = x1/(x1 + x2 - one)
        TCLim%Lim_z(2) = x2/(x1 + x2 - one)
        TCLim%Lim_z(3) = one/(one - x1 - x2)

        !! set sij in the limit. sij = 2*pi.pj.
        !! Only those that are needed are computed
        TCLim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
        TCLim%Lim_sij(i,6) = spart*x2*HardProc%Lim_etaij(i,6)
        TCLim%Lim_sij(5,6) = spart*x1*x2*HardProc%Lim_etaij(5,6)

        if(present(opt_etas)) call fill_etas(TCLim,ns,opt_etas)

        TCLim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * x1*x2*x3*x4/NF1(x3,one-x4_half,lambda) &
             * spart/mv2      &
             * jac

      !! ----------------------------- TC + C6 ----------------------------- !!
        if(present(TCC6Lim)) then

          call initialize_config(TCC6Lim)

          TCC6Lim%Lim_z(1) = x1/(x1 + x2)
          TCC6Lim%Lim_z(2) = x2/(x1 + x2)
          TCC6Lim%Lim_z(3) = one/(one - x1 - x2)

          !! si5 + si6 - s56 in the x4 -> 0 limit
          TCC6Lim%Lim_KinInv(1) = -spart*(x1 + x2)*x3
          TCC6Lim%Lim_sij(5,6) = spart*x1*x2*x3*x4_half**2/NF1(x3,one,lambda)

          TCC6Lim%wgt = one/8.0_dp/pi * kallenF &
               * one/two/spart  &
               * (spart/4)**2   &
               * x1*x2*x3*x4/NF1(x3,one,lambda) &
               * spart/mv2      &
               * jac

        endif

      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                    STRUCTURE 5: TC + S5 FLM.                          !!
    !!-----------------------------------------------------------------------!!

    if(present(TCS5Lim)) then

      call initialize_config(TCS5Lim)

      !! x3 -> 0
      eta65 = zero
      TCS5Lim%Lim_etaij(1,2) = one
      TCS5Lim%Lim_etaij(i,5) = zero
      TCS5Lim%Lim_etaij(i,6) = zero
      TCS5Lim%Lim_etaij(j,5) = one
      TCS5Lim%Lim_etaij(j,6) = one
      TCS5Lim%Lim_etaij(5,6) = eta65

      eta51 = TCS5Lim%Lim_etaij(1,5)
      eta52 = TCS5Lim%Lim_etaij(2,5)
      eta61 = TCS5Lim%Lim_etaij(1,6)
      eta62 = TCS5Lim%Lim_etaij(2,6)

      call get_part_s_xi_nlo(mv2,ylab,x2,eta61,eta62,spart,sqrts,xi1,xi2,TCS5Lim%flag)

      !! ------------------------------- TC + S5 --------------------------- !!

      if (.not. TCS5Lim%flag) then

        TCS5Lim%PartFrac = (/xi1,xi2/)

        ns(:,5) = t + sgn*e3
        ns(:,6) = ns(:,5)

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p6 = x2*half*sqrts*ns(:,6)

        if(i .eq. 1) then
          p1 = p1 - p6
        elseif(i .eq. 2) then
          p2 = p2 - p6
        endif

        pv = p1 + p2

        call get_lept_mom(mv, pv,  yr(xLMIN:xLMAX), nlept, kallenF, TCS5Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        TCS5Lim%AmpMom(:,1) = p1
        TCS5Lim%AmpMom(:,2) = p2

        TCS5Lim%Lim_z(1) = x1/(x2 - one)
        TCS5Lim%Lim_z(2) = x2/(x2 - one)
        TCS5Lim%Lim_z(3) = one/(one - x2)

        !! set sij in the limit. sij = 2*pi.pj.
        !! Only those that are needed are computed
        TCS5Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
        TCS5Lim%Lim_sij(i,6) = spart*x2*HardProc%Lim_etaij(i,6)
        TCS5Lim%Lim_sij(5,6) = spart*x1*x2*x3*x4_half**2/NF1(x3,one-x4_half,lambda)

        if(present(opt_etas)) call fill_etas(TCS5Lim,ns,opt_etas)

        TCS5Lim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * x1*x2*x3*x4/NF1(x3,one-x4_half,lambda) &
             * spart/mv2      &
             * jac

      !! --------------------------- TC + S5 + C6 -------------------------- !!

        if(present(TCC6S5Lim)) then

          TCC6S5Lim%Lim_sij(i,5) = spart*x1*x3
          TCC6S5Lim%Lim_sij(i,6) = spart*x2*x3
          TCC6S5Lim%Lim_sij(5,6) = spart*x1*x2*x3*x4_half**2/NF1(x3,one,lambda)

          TCC6S5Lim%wgt = one/8.0_dp/pi * kallenF &
               * one/two/spart  &
               * (spart/4)**2   &
               * x1*x2*x3*x4/NF1(x3,one,lambda) &
               * spart/mv2      &
               * jac
        endif

      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                    STRUCTURE 6: TC + S6 FLM.                          !!
    !!-----------------------------------------------------------------------!!

    if(present(TCS6Lim)) then

      call initialize_config(TCS6Lim)

      !! x3 -> 0
      eta65 = zero
      TCS6Lim%Lim_etaij(1,2) = one
      TCS6Lim%Lim_etaij(i,5) = zero
      TCS6Lim%Lim_etaij(i,6) = zero
      TCS6Lim%Lim_etaij(j,5) = one
      TCS6Lim%Lim_etaij(j,6) = one
      TCS6Lim%Lim_etaij(5,6) = eta65

      eta51 = TCS6Lim%Lim_etaij(1,5)
      eta52 = TCS6Lim%Lim_etaij(2,5)
      eta61 = TCS6Lim%Lim_etaij(1,6)
      eta62 = TCS6Lim%Lim_etaij(2,6)

      call get_part_s_xi_nlo(mv2,ylab,x1,eta51,eta52,spart,sqrts,xi1,xi2,TCS6Lim%flag)

      !! ------------------------------- TC + S6 --------------------------- !!

      if (.not. TCS6Lim%flag) then

        TCS6Lim%PartFrac = (/xi1,xi2/)

        ns(:,5) = t + sgn*e3
        ns(:,6) = ns(:,5)

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = x1*half*sqrts*ns(:,5)

        if(i .eq. 1) then
          p1 = p1 - p5
        elseif(i .eq. 2) then
          p2 = p2 - p5
        endif

        pv = p1 + p2

        call get_lept_mom(mv, pv,  yr(xLMIN:xLMAX), nlept, kallenF, TCS6Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        TCS6Lim%AmpMom(:,1) = p1
        TCS6Lim%AmpMom(:,2) = p2

        TCS6Lim%Lim_z(1) = x1/(x2 - one)
        TCS6Lim%Lim_z(2) = x2/(x2 - one)
        TCS6Lim%Lim_z(3) = one/(one-x1)

        !! set sij in the limit. sij = 2*pi.pj.
        !! Only those that are needed are computed
        TCS6Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
        TCS6Lim%Lim_sij(i,6) = spart*x2*HardProc%Lim_etaij(i,6)
        TCS6Lim%Lim_sij(5,6) = spart*x1*x2*HardProc%Lim_etaij(5,6)

        if(present(opt_etas)) call fill_etas(TCS6Lim,ns,opt_etas)

        TCS6Lim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * x1*x2*x3*x4/NF1(x3,one-x4_half,lambda) &
             * spart/mv2      &
             * jac

      !! --------------------------- TC + S6 + C6 -------------------------- !!

        if(present(TCC6S6Lim)) then

          TCC6S6Lim%Lim_z(3) = x1/(x1-one)

          TCC6S6Lim%Lim_sij(i,6) = spart*x2*x3
          TCC6S6Lim%Lim_sij(5,6) = spart*x1*x2*x3*x4_half**2/NF1(x3,one,lambda)

          TCC6S6Lim%wgt = one/8.0_dp/pi * kallenF &
               * one/two/spart  &
               * (spart/4)**2   &
               * x1*x2*x3*x4/NF1(x3,one,lambda) &
               * spart/mv2      &
               * jac
        endif

      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                       STRUCTURE 7: C6 FLM.                            !!
    !!-----------------------------------------------------------------------!!

    if(present(C6Lim)) then

      call initialize_config(C6Lim)

      !! x4 -> 0
      eta65 = zero
      C6Lim%Lim_etaij(1,2) = one
      C6Lim%Lim_etaij(i,5) = x3
      C6Lim%Lim_etaij(i,6) = x3
      C6Lim%Lim_etaij(j,5) = one-x3
      C6Lim%Lim_etaij(j,6) = one-x3
      C6Lim%Lim_etaij(5,6) = zero

      eta51 = C6Lim%Lim_etaij(1,5)
      eta52 = C6Lim%Lim_etaij(2,5)
      eta61 = C6Lim%Lim_etaij(1,6)
      eta62 = C6Lim%Lim_etaij(2,6)

      call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
        spart,sqrts,xi1,xi2,C6Lim%flag)

      if (.not. C6Lim%flag) then 

        C6Lim%PartFrac = (/xi1,xi2/)

        cos5i = sgn*(one - 2*x3)
        sin5i = two*sqrt(x3*(one-x3))

        ns(:,5) = t + cos5i*e3 + sin5i*b
        ns(:,6) = ns(:,5)

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = x1*half*sqrts*ns(:,5)
        p6 = x2*half*sqrts*ns(:,6)

        p5 = p5 + p6

        pv = p1 + p2 - p5

        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C6Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        C6Lim%AmpMom(:,1) = p1
        C6Lim%AmpMom(:,2) = p2
        C6Lim%AmpMom(:,5) = p5

        C6Lim%Lim_z(1) = x1/(x1+x2)
        C6Lim%Lim_z(2) = x2/(x1+x2)

        !! set sij in the limit. sij = 2*pi.pj.
        !! Only those that are needed are computed
        C6Lim%Lim_sij(5,6) = spart*x1*x2*x3*x4_half**2/NF1(x3,one,lambda)

        if(present(opt_etas)) call fill_etas(C6Lim,ns,opt_etas)

        C6Lim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * x1*x2*x3*x4/NF1(x3,one,lambda) &
             * spart/mv2      &
             * jac

      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                       STRUCTURE 9: S5 + C6 FLM.                       !!
    !!-----------------------------------------------------------------------!!

    if(present(C6S5Lim)) then

      call initialize_config(C6S5Lim)

      !! x4 -> 0
      eta65 = zero
      C6S5Lim%Lim_etaij(i,5) = x3
      C6S5Lim%Lim_etaij(i,6) = x3
      C6S5Lim%Lim_etaij(j,5) = one-x3
      C6S5Lim%Lim_etaij(j,6) = one-x3
      C6S5Lim%Lim_etaij(5,6) = zero

      eta51 = C6S5Lim%Lim_etaij(1,5)
      eta52 = C6S5Lim%Lim_etaij(2,5)
      eta61 = C6S5Lim%Lim_etaij(1,6)
      eta62 = C6S5Lim%Lim_etaij(2,6)

      call get_part_s_xi_nlo(mv2,ylab,x2,eta61,eta62,spart,sqrts,xi1,xi2,C6S5Lim%flag)

      if (.not. C6S5Lim%flag) then

        C6S5Lim%PartFrac = (/xi1,xi2/)

        cos5i = sgn*(one - 2*x3)
        sin5i = two*sqrt(x3*(one-x3))

        ns(:,5) = t + cos5i*e3 + sin5i*b
        ns(:,6) = ns(:,5)

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p6 = x2*half*sqrts*ns(:,6)

        pv = p1 + p2 - p6

        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C6S5Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        C6S5Lim%AmpMom(:,1) = p1
        C6S5Lim%AmpMom(:,2) = p2
        C6S5Lim%AmpMom(:,5) = p6

        C6S5Lim%Lim_z(1) = x2/x1
        C6S5Lim%Lim_z(2) = x1/x2
        C6S5Lim%Lim_sij(5,6) = spart*x1*x2*x3*x4_half**2/NF1(x3,one,lambda)

        if(present(opt_etas)) call fill_etas(C6S5Lim,ns,opt_etas)

        C6S5Lim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * x1*x2*x3*x4/NF1(x3,one,lambda) &
             * spart/mv2      &
             * jac

      endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                      STRUCTURE 10: S6 + C6 FLM.                       !!
    !!-----------------------------------------------------------------------!!

    if(present(C6S6Lim)) then

      call initialize_config(C6S6Lim)

      !! x4 -> 0
      eta65 = zero
      C6S6Lim%Lim_etaij(i,5) = x3
      C6S6Lim%Lim_etaij(i,6) = x3
      C6S6Lim%Lim_etaij(j,5) = one-x3
      C6S6Lim%Lim_etaij(j,6) = one-x3
      C6S6Lim%Lim_etaij(5,6) = zero

      eta51 = C6S6Lim%Lim_etaij(1,5)
      eta52 = C6S6Lim%Lim_etaij(2,5)
      eta61 = C6S6Lim%Lim_etaij(1,6)
      eta62 = C6S6Lim%Lim_etaij(2,6)

      call get_part_s_xi_nlo(mv2,ylab,x1,eta51,eta52,spart,sqrts,xi1,xi2,C6S6Lim%flag)

      if (.not. C6S6Lim%flag) then 

        C6S6Lim%PartFrac = (/xi1,xi2/)

        cos5i = sgn*(one - 2*x3)
        sin5i = two*sqrt(x3*(one-x3))

        ns(:,5) = t + cos5i*e3 + sin5i*b
        ns(:,6) = ns(:,5)

        p1 = half*sqrts*ns(:,1)
        p2 = half*sqrts*ns(:,2)
        p5 = x1*half*sqrts*ns(:,5)

        pv = p1 + p2 - p5

        call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C6S6Lim)
        ns(:,3) = (/one,nlept(:,1)/)
        ns(:,4) = (/one,nlept(:,2)/)

        C6S6Lim%AmpMom(:,1) = p1
        C6S6Lim%AmpMom(:,2) = p2
        C6S6Lim%AmpMom(:,5) = p5

        C6S6Lim%Lim_z(1) = x2/x1
        C6S6Lim%Lim_sij(5,6) = spart*x1*x2*x3*x4_half**2/NF1(x3,one,lambda)

        if(present(opt_etas)) call fill_etas(C6S6Lim,ns,opt_etas)

        C6S6Lim%wgt = one/8.0_dp/pi * kallenF &
             * one/two/spart  &
             * (spart/4)**2   &
             * x1*x2*x3*x4/NF1(x3,one,lambda) &
             * spart/mv2      &
             * jac

      endif
    endif

  end subroutine kinematics_nnlo_abcd_5i6ibd

end module mod_kinematics_nnlo_tc_is_eu_abcd

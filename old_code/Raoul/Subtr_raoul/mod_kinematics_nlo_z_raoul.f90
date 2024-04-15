!-- does not include damping factors: for nlo qcd, we do not partition
!-- also, do not set particle types because we don't know whether we are emitting
!   photons or gluons
module mod_kinematics_nlo_z_raoul
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_process
  use mod_auxfunctions
  use mod_kinematics_gen
  use mod_kinematics_lept
  implicit none
  integer, public, parameter :: NLOdim_vegas = tauy_max + 2 + 3
  !
  integer, public, parameter :: kNLO_min = tauy_max + 1
  integer, public, parameter :: kNLO_max = tauy_max + 2 + 3
  integer, public, parameter :: kNLO_max_full = kNLO_max + 1
  !
  integer, public, parameter :: xE   = kNLO_min
  integer, public, parameter :: xRHO = kNLO_min + 1
  private

  public :: kinematics_nlo_z_is
  public :: kinematics_nlo_z_fs

contains

  subroutine checkandgetxis(MV2,Y,sh,z1,z2,x1,eta41,eta42,flag,spart,sqrts,xi1,xi2)
    real(dp), intent(in) :: MV2,Y,sh,z1,z2,x1,eta41,eta42
    logical, intent(out) :: flag
    real(dp), intent(out) :: xi1,xi2,spart,sqrts
    real(dp) :: Deltas,Deltan,Deltad,yy

    flag = .false.
    xi1 = two
    xi2 = two
    spart = zero
    sqrts = zero

    Deltas = z1*z2 - z1*x1*eta41 - z2*x1*eta42


    if (Deltas.lt.MV2/sh) flag = .true.
    if (flag .eqv. .false.) spart = MV2/Deltas
    if (flag .eqv. .false.) sqrts = sqrt(spart)

    Deltan = z1 - x1*eta42
    Deltad = z2 - x1*eta41
    if ( (flag .eqv. .false.).and.(Deltan.gt.zero).and.(Deltad.gt.zero) ) then       
       yy = Y - half*log(Deltan/Deltad)
       xi1 = sqrt(spart/sh)*exp(yy)
       xi2 = sqrt(spart/sh)*exp(-yy)
       if (xi1.gt.one .or. xi2.gt.one) flag = .true.
    else
       flag = .true.
    endif
    return

  end subroutine checkandgetxis

  
  !-- both 51 and 52
  subroutine kinematics_nlo_z_is(yr,z1,z2,HardProc,C1Lim,C2Lim,SLim,SC1Lim,SC2Lim,compute_etas)
    real(dp), intent(in) :: yr(kNLO_max_full),z1,z2
    type(KinConfig) :: HardProc
    type(KinConfig), optional :: C1Lim,SC1Lim,C2Lim,SC2Lim,SLim
    logical,  intent(in) :: compute_etas
    real(dp) :: tau,ylab,jac
    real(dp) :: x1,x2,sin5,cos5,phi5
    real(dp) :: xi1,xi2
    real(dp) :: slocal,sqrts,mv,mv2
    real(dp) :: p1(4),p2(4),p5(4),pv(4)
    real(dp) :: kallenF,ni(3,5)
    real(dp) :: z,s51,s52,eta5i(4)
    integer :: i,j

    call initialize_config(HardProc)

    call get_tauy(yr(1:tauy_max),tau,ylab,jac)

    !-- variable assignment
    x1   = yr(kNLO_min)         !-- E5 3
    x2   = yr(kNLO_min+1)       !-- eta51 4 
    phi5 = twopi*yr(kNLO_min+2) !-- 5

    cos5 = one-two*x2
    sin5 = two*sqrt(abs(x2*(one-x2)))

    MV2 = tau*sh
    MV = sqrt(MV2)

    eta5i(1) = x2
    eta5i(2) = one-x2

    !-- HardProc
    call checkandgetxis(MV2,ylab,sh,z1,z2,x1,eta5i(1),eta5i(2),HardProc%flag,slocal,sqrts,xi1,xi2)


    if (HardProc%flag .eqv. .false.) then
       
       HardProc%PartFrac = [xi1,xi2]

       p1 = half*sqrts*z1*[one,zero,zero, one]
       p2 = half*sqrts*z2*[one,zero,zero,-one]
       p5 = half*sqrts*x1*[one,sin5*cos(phi5),sin5*sin(phi5),cos5]
       
       pv  = p1+p2-p5
!       mv2 = slocal*(one-x1)
!       mv  = sqrt(mv2)
       call get_lept_mom(mv,pv,yr(kNLO_min+3:kNLO_max_full),ni(:,3:4),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,5) = p5

       !-- flux: 1/2/(z1*z2*slocal)
       !-- xi1xi2jac = spart/MV2
       ! -- -> wgt = kallenF/32/pi * x1/z1/z2 *slocal/MV2
       HardProc%wgt = (kallenF/32.0_dp/pi) * (x1/z1/z2) * slocal/MV2  * jac
       HardProc%npart = 5

       HardProc%Lim_Ei(1:5) = [sqrts/two,sqrts/two,HardProc%AmpMom(1,3),HardProc%AmpMom(1,4),HardProc%AmpMom(1,5)]

       if (compute_etas) then

          ni(:,1) = [0,0, 1]
          ni(:,2) = [0,0,-1]
          ni(:,5) = [sin5*cos(phi5),sin5*sin(phi5),cos5]
          do i = 3,4
             HardProc%Lim_etaij(1,i) = half*(one-dot_product(ni(:,i),ni(:,1)))
             HardProc%Lim_etaij(2,i) = half*(one-dot_product(ni(:,i),ni(:,2)))
             HardProc%Lim_etaij(5,i) = half*(one-dot_product(ni(:,i),ni(:,5)))
             HardProc%Lim_etaij(i,5) = HardProc%Lim_etaij(5,i) 
          enddo
          HardProc%Lim_etaij(3,4) = half*(one-dot_product(ni(:,3),ni(:,4)))

          HardProc%Lim_etaij(1,5) = x2
          HardProc%Lim_etaij(2,5) = one-x2

       endif
       
    endif

    !-- C1
    if (present(C1Lim)) then

       call initialize_config(C1Lim)

       eta5i(1) = zero
       eta5i(2) = one

       call checkandgetxis(MV2,ylab,sh,z1,z2,x1,eta5i(1),eta5i(2),C1Lim%flag,slocal,sqrts,xi1,xi2)


       if (C1Lim%flag .eqv. .false.) then

       
          C1Lim%PartFrac = [xi1,xi2]
          
          p1 = half*sqrts*(z1-x1)*[one,zero,zero,one]
          p2 = half*sqrts*z2*[one,zero,zero,-one]
          
          pv  = p1+p2
          call get_lept_mom(mv,pv,yr(kNLO_min+3:kNLO_max_full),ni(:,3:4),kallenF,C1Lim)
          
          C1Lim%AmpMom(:,1) = p1
          C1Lim%AmpMom(:,2) = p2

          C1Lim%LimMom(:,1) = [zero,cos(phi5),sin(phi5),zero] !-- nperp
          C1Lim%LimMom(:,2) = [zero,zero,zero,one]             !-- ndir

          s51 = slocal*x1*x2*z1
          z   = x1/z1 !-- use Pqg instead of Pqq to get the soft limit in x and not in 1-(1-x)
          C1Lim%Lim_KinInv(1:3) = [z,s51,x2]

          C1Lim%Lim_Ei(1:5) = [sqrts/two,sqrts/two,C1Lim%AmpMom(1,3),C1Lim%AmpMom(1,4),half*sqrts*x1]

          C1Lim%Lim_etaij(1,5) = x2
          C1Lim%Lim_etaij(2,5) = one-x2

          C1Lim%wgt = (kallenF/32.0_dp/pi) * (x1/z1/z2) * slocal/MV2  * jac          
          C1Lim%npart = 4

       if (compute_etas) then

          ni(:,1) = [0,0, 1]
          ni(:,2) = [0,0,-1]
          do i = 3,4
             C1Lim%Lim_etaij(1,i) = half*(one-dot_product(ni(:,i),ni(:,1)))
             C1Lim%Lim_etaij(2,i) = half*(one-dot_product(ni(:,i),ni(:,2)))
          enddo
          C1Lim%Lim_etaij(3,4) = half*(one-dot_product(ni(:,3),ni(:,4)))

       endif

          
       endif
    endif

    !-- C2
    if (present(C2Lim)) then

       call initialize_config(C2Lim)

       eta5i(1) = one
       eta5i(2) = zero

       call checkandgetxis(MV2,ylab,sh,z1,z2,x1,eta5i(1),eta5i(2),C2Lim%flag,slocal,sqrts,xi1,xi2)


       if (C2Lim%flag .eqv. .false.) then

          C2Lim%PartFrac = [xi1,xi2]
          
          slocal = sh*xi1*xi2
          sqrts  = sqrt(slocal)
          
          p1 = half*sqrts*z1*[one,zero,zero,one]
          p2 = half*sqrts*(z2-x1)*[one,zero,zero,-one]
          
          pv  = p1+p2
!          mv2 = slocal*(one-x1)
!          mv  = sqrt(mv2)
          call get_lept_mom(mv,pv,yr(kNLO_min+3:kNLO_max_full),ni(:,3:4),kallenF,C2Lim)
          
          C2Lim%AmpMom(:,1) = p1
          C2Lim%AmpMom(:,2) = p2
          
          s52 = slocal*x1*(one-x2)*z2
          z   = x1/z2 !-- Use Pqq, see comment above
          C2Lim%Lim_KinInv(1:3) = [z,s52,one-x2]

          C2Lim%LimMom(:,1) = [zero,cos(phi5),sin(phi5),zero] !-- nperp
          C2Lim%LimMom(:,2) = [zero,zero,zero,-one]             !-- ndir

          C2Lim%Lim_etaij(1,5) = x2
          C2Lim%Lim_etaij(2,5) = one-x2

          C2Lim%Lim_Ei(1:5) = [sqrts/two,sqrts/two,C2Lim%AmpMom(1,3),C2Lim%AmpMom(1,4),half*sqrts*x1]
          C2Lim%wgt = (kallenF/32.0_dp/pi) * (x1/z1/z2) * slocal/MV2  * jac
          C2Lim%npart = 4

       if (compute_etas) then

          ni(:,1) = [0,0, 1]
          ni(:,2) = [0,0,-1]
          do i = 3,4
             C2Lim%Lim_etaij(1,i) = half*(one-dot_product(ni(:,i),ni(:,1)))
             C2Lim%Lim_etaij(2,i) = half*(one-dot_product(ni(:,i),ni(:,2)))
          enddo

          C2Lim%Lim_etaij(3,4) = half*(one-dot_product(ni(:,3),ni(:,4)))

       endif
          
       endif
    endif
       
    !-- S, SC1, SC2
    if (present(SLim)) then

      call initialize_config(SLim)

      eta5i(1) = x2
      eta5i(2) = one-x2
      

       call checkandgetxis(MV2,ylab,sh,z1,z2,zero,eta5i(1),eta5i(2),SLim%flag,slocal,sqrts,xi1,xi2)

       if (SLim%flag .eqv. .false.) then

          SLim%PartFrac = [xi1,xi2]

          p1 = half*sqrts*z1*[one,zero,zero, one]
          p2 = half*sqrts*z2*[one,zero,zero,-one]
          
          pv  = p1+p2
!          mv  = sqrts
          call get_lept_mom(mv,pv,yr(kNLO_min+3:kNLO_max_full),ni(:,3:4),kallenF,SLim)
          
          SLim%AmpMom(:,1) = p1
          SLim%AmpMom(:,2) = p2
          
          !-- prepare the etas
          ni(:,1) = [0,0, 1]
          ni(:,2) = [0,0,-1]
          ni(:,5) = [sin5*cos(phi5),sin5*sin(phi5),cos5]
          do i = 1,5
             do j = i+1,5
                if (i.lt.3 .and. j.eq. 5) cycle
                SLim%Lim_etaij(i,j) = half*(one-dot_product(ni(:,i),ni(:,j)))
             enddo
          enddo
          
          !-- now the ones with 5 and the collinear direction
          SLim%Lim_etaij(1,5) = x2
          SLim%Lim_etaij(2,5) = one-x2

          SLim%Lim_KinInv(1) = half*sqrts*x1 !-- e5

          SLim%Lim_Ei(1:5) = [sqrts/two,sqrts/two,SLim%AmpMom(1,3),SLim%AmpMom(1,4),half*sqrts*x1]
          
          SLim%wgt = (kallenF/32.0_dp/pi) * (x1/z1/z2) * slocal/MV2  * jac
          SLim%npart = 4

       if (compute_etas) then

          ni(:,1) = [0,0, 1]
          ni(:,2) = [0,0,-1]
          do i = 3,4
             SLim%Lim_etaij(1,i) = half*(one-dot_product(ni(:,i),ni(:,1)))
             SLim%Lim_etaij(2,i) = half*(one-dot_product(ni(:,i),ni(:,2)))
             SLim%Lim_etaij(5,i) = half*(one-dot_product(ni(:,i),ni(:,5)))
             SLim%Lim_etaij(i,5) = SLim%Lim_etaij(5,i)
          enddo

          SLim%Lim_etaij(3,4) = half*(one-dot_product(ni(:,3),ni(:,4)))
          SLim%Lim_etaij(1,5) = x2
          SLim%Lim_etaij(2,5) = one-x2

       endif

       endif
    endif

    !-- SC1
     if (present(SC1Lim)) then

        call initialize_config(SC1Lim)

        eta5i(1) = zero
        eta5i(2) = one

        call checkandgetxis(MV2,ylab,sh,z1,z2,zero,eta5i(1),eta5i(2),SC1Lim%flag,slocal,sqrts,xi1,xi2)


        if (SC1Lim%flag .eqv. .false.) then

           SC1Lim%PartFrac = [xi1,xi2]

           p1 = half*sqrts*z1*[one,zero,zero,one]
           p2 = half*sqrts*z2*[one,zero,zero,-one]

           pv  = p1+p2
           call get_lept_mom(mv,pv,yr(kNLO_min+3:kNLO_max_full),ni(:,3:4),kallenF,SC1Lim)

           SC1Lim%AmpMom(:,1) = p1
           SC1Lim%AmpMom(:,2) = p2

           SC1Lim%Lim_KinInv(1:2) = [half*sqrts*x1,x2] !-- E5,eta
           SC1Lim%Lim_Ei(1:5) = [sqrts/two,sqrts/two,SC1Lim%AmpMom(1,3),SC1Lim%AmpMom(1,4),half*sqrts*x1]

           SC1Lim%Lim_etaij(1,5) = x2
           SC1Lim%Lim_etaij(2,5) = one-x2

           SC1Lim%wgt = (kallenF/32.0_dp/pi) * (x1/z1/z2) * slocal/MV2  * jac
           SC1Lim%npart = 4

       if (compute_etas) then

          ni(:,1) = [0,0, 1]
          ni(:,2) = [0,0,-1]
          do i = 3,4
             SC1Lim%Lim_etaij(1,i) = half*(one-dot_product(ni(:,i),ni(:,1)))
             SC1Lim%Lim_etaij(2,i) = half*(one-dot_product(ni(:,i),ni(:,2)))
          enddo
          SC1Lim%Lim_etaij(3,4) = half*(one-dot_product(ni(:,3),ni(:,4)))

       endif

        endif
     endif


     !-- SC2
     if (present(SC2Lim)) then

        call initialize_config(SC2Lim)

        eta5i(1) = one
        eta5i(2) = zero

        call checkandgetxis(MV2,ylab,sh,z1,z2,zero,eta5i(1),eta5i(2),SC2Lim%flag,slocal,sqrts,xi1,xi2)


        if (SC2Lim%flag .eqv. .false.) then

           SC2Lim%PartFrac = [xi1,xi2]

           slocal = sh*xi1*xi2
           sqrts  = sqrt(slocal)

           p1 = half*sqrts*z1*[one,zero,zero,one]
           p2 = half*sqrts*z2*[one,zero,zero,-one]

           pv  = p1+p2
           call get_lept_mom(mv,pv,yr(kNLO_min+3:kNLO_max_full),ni(:,3:4),kallenF,SC2Lim)

           SC2Lim%AmpMom(:,1) = p1
           SC2Lim%AmpMom(:,2) = p2


           SC2Lim%Lim_KinInv(1:2) = [half*sqrts*x1,one-x2] !-- E5,eta
           SC2Lim%Lim_Ei(1:5) = [sqrts/two,sqrts/two,SC2Lim%AmpMom(1,3),SC2Lim%AmpMom(1,4),half*sqrts*x1]

           SC2Lim%Lim_etaij(1,5) = x2
           SC2Lim%Lim_etaij(2,5) = one-x2

           SC2Lim%wgt = (kallenF/32.0_dp/pi) * (x1/z1/z2) * slocal/MV2  * jac
           SC2Lim%npart = 4

       if (compute_etas) then

          ni(:,1) = [0,0, 1]
          ni(:,2) = [0,0,-1]
          do i = 3,4
             SC2Lim%Lim_etaij(1,i) = half*(one-dot_product(ni(:,i),ni(:,1)))
             SC2Lim%Lim_etaij(2,i) = half*(one-dot_product(ni(:,i),ni(:,2)))
          enddo
          SC2Lim%Lim_etaij(3,4) = half*(one-dot_product(ni(:,3),ni(:,4)))

       endif

        endif
     endif

   end subroutine kinematics_nlo_z_is

  !-- Both sectors. icoll: lepton collinear to emission; jother: other lepton
  subroutine kinematics_nlo_z_fs(yr,z1,z2,icoll,jother,HardProc,CLim,SCLim,SLim)
    real(dp), intent(in) :: yr(kNLO_max_full),z1,z2
    integer,  intent(in) :: icoll, jother
    type(KinConfig) :: HardProc,CLim,SCLim,SLim
    real(dp) :: tau,ylab,jac
    real(dp) :: xi1,xi2
    real(dp) :: x1,x2,sin5,cos5,phi5,eta51,eta52
    real(dp) :: slocal,sqrts,mv,mv2,flux
    real(dp) :: p1(4),p2(4),pv(4),p5(4),q(4),pi_aux(4),pi(4),pj(4)
    real(dp) :: cosi,sini,cphi,sphi
    real(dp) :: nicoll_ta(3),nicoll_tb(3),ni(3,5),e5,ei,ej,s_qi
    real(dp) :: si5,z
    integer  :: i,j

    call initialize_config(HardProc)
    call initialize_config(CLim)
    call initialize_config(SCLim)
    call initialize_config(SLim)

    call get_tauy(yr(1:tauy_max),tau,ylab,jac)

    !-- variable assignment
    x1   = yr(kNLO_min)         !-- E5 
    x2   = yr(kNLO_min+1)       !-- eta5
    phi5 = twopi*yr(kNLO_min+2) 
    
    cos5 = one-two*x2
    sin5 = two*sqrt(abs(x2*(one-x2)))

    eta51 = x2
    eta52 = one-x2


    MV2 = tau*sh
    MV = sqrt(MV2)

!    call checkandgetxis(MV2,ylab,sh,z1,z2,x1,eta51,eta52,all_conf_flag,slocal,sqrts,xi1,xi2)   ! only need to do this once for all configs

!    if (all_conf_flag) then

    slocal = mv2/z1/z2
    sqrts = sqrt(slocal)
    xi1 = sqrt(slocal/sh)*sqrt(z2/z1)*exp(+ylab)
    xi2 = sqrt(slocal/sh)*sqrt(z1/z2)*exp(-ylab)

    if (xi1 .gt. one .or. xi2 .gt. one) then

       HardProc%flag = .true.
       CLim%flag     = .true.
       SCLim%flag    = .true.
       SLim%flag     = .true.

    else

       flux = one/(two*mv)
       
       p1 = half*sqrts*z1*[one,zero,zero, one]
       p2 = half*sqrts*z2*[one,zero,zero,-one]
       pv = p1+p2

       ni(:,1) = [zero,zero, one]
       ni(:,2) = [zero,zero,-one]
       
       !-- randomize lepton axis
       cosi = one-two*yr(kNLO_min+3)
       sini = sqrt(abs(one-cosi**2))
       cphi = cos(twopi*yr(kNLO_min+4))
       sphi = sin(twopi*yr(kNLO_min+4))

       ni(:,icoll) = [sini*cphi,sini*sphi, cosi]
       nicoll_ta   = [cosi*cphi,cosi*sphi,-sini]
       nicoll_tb   = [-sphi,cphi,0._dp]
       pi_aux(:) = [one,ni(:,icoll)]

       ni(:,5) = cos5*ni(:,icoll) + sin5*(cos(phi5)*nicoll_ta(:) + sin(phi5)*nicoll_tb(:))

       !-- hard kinematics
       e5 = half*sqrts*x1
       p5 = e5*[one,ni(:,5)]
       q  = pv - p5
       s_qi = scr(q,pi_aux)
       ei = scr(q,q)/two/s_qi
       pi = ei*pi_aux
       pj = q - pi
       ej = pj(1)
       
       !-- hard directions             
       if (ei.lt.zero .or. ej.lt.zero) then

          HardProc%flag = .true.

       else

          HardProc%PartFrac = [xi1,xi2]
       
          HardProc%AmpMom(:,1) = p1
          HardProc%AmpMom(:,2) = p2
          HardProc%AmpMom(:,icoll)  = pi
          HardProc%AmpMom(:,jother) = pj
          HardProc%AmpMom(:,5) = p5

          HardProc%Lim_Ei(1:5) = [sqrts/two,sqrts/two,HardProc%AmpMom(1,3),HardProc%AmpMom(1,4),HardProc%AmpMom(1,5)]

          !-- prepare etas for the damping
          ni(:,jother) = 1/ej*pj(2:4)
          do i = 1,4
             if (i==icoll) then
                HardProc%Lim_etaij(i,5) = x2
             else
                HardProc%Lim_etaij(i,5) = half*(one-dot_product(ni(:,i),ni(:,5)))
             endif
             HardProc%Lim_etaij(1,i) = half*(one-dot_product(ni(:,i),ni(:,1)))
             HardProc%Lim_etaij(2,i) = half*(one-dot_product(ni(:,i),ni(:,2)))
          enddo
          HardProc%Lim_etaij(3,4) = half*(one-dot_product(ni(:,3),ni(:,4)))

          HardProc%wgt = one/(two*mv) &
               * (mv/two)**2 * x1/z1/z2 &
               * ei/s_qi/twopi &
               * two * flux * jac


          HardProc%npart = 5
          
       endif

       !-- C
       p5 = e5*pi_aux
       q  = pv - p5
       s_qi = scr(q,pi_aux)
       ei = scr(q,q)/two/s_qi
       pi = ei*pi_aux
       pj = q - pi
       ej = pj(1)

       if (ei.lt.zero .or. ej.lt.zero) then

          CLim%flag = .true.

       else

          CLim%PartFrac = [xi1,xi2]
       
          CLim%AmpMom(:,1) = p1
          CLim%AmpMom(:,2) = p2
          CLim%AmpMom(:,icoll)  = pi+p5
          CLim%AmpMom(:,jother) = pj

          CLim%wgt = one/(two*mv) &
               * (mv/two)**2 * x1/z1/z2 &
!               * (1-x1)/two/twopi &     
!               * ei/s_qi/twopi &
               *two*(z1*z2 - x1/two*(z1+z2-cosi*(z1-z2)))/(z1+z2-cosi*(z1-z2))**2/twopi &
               * two * flux * jac
          CLim%npart = 4

          !-- prepare etas for the ONLO subtraction functions
          ni(:,jother) = 1/ej*pj(2:4)
          do i = 1,4
             if (i==icoll) then
                CLim%Lim_etaij(i,5) = x2
             else
                CLim%Lim_etaij(i,5) = half*(one-dot_product(ni(:,i),ni(:,5)))
             endif
             CLim%Lim_etaij(1,i) = half*(one-dot_product(ni(:,i),ni(:,1)))
             CLim%Lim_etaij(2,i) = half*(one-dot_product(ni(:,i),ni(:,2)))
          enddo
          CLim%Lim_etaij(3,4) = half*(one-dot_product(ni(:,3),ni(:,4)))


          si5 = four*ei*e5*x2
          z = e5/(e5+ei) !-- to be used with Pqg, so that sing is 1/x1
                         !-- and not 1/[1-(1-x1)] which leads to precision loss

          CLim%Lim_KinInv(1:2) = [z,si5]
          CLim%Lim_Ei(1:5) = [sqrts/two,sqrts/two,CLim%AmpMom(1,3),CLim%AmpMom(1,4),e5]
          
       endif

       !-- S and SC
       q  = pv
       s_qi = scr(q,pi_aux)
       ei = scr(q,q)/two/s_qi
       pi = ei*pi_aux
       pj = q - pi
       ej = pj(1)


       if (ei.lt.zero .or. ej.lt.zero) then

          SLim%flag = .true.
          SCLim%flag = .true.

       else

          !-- S
          SLim%PartFrac = [xi1,xi2]
       
          SLim%AmpMom(:,1) = p1
          SLim%AmpMom(:,2) = p2
          SLim%AmpMom(:,icoll)  = pi
          SLim%AmpMom(:,jother) = pj

          !-- prepare ij and damping
          ni(:,jother) = 1/ej*pj(2:4)
          do i = 1,4
             do j = i+1,5
                if (i==icoll .and. j==5) then
                   SLim%Lim_etaij(i,j) = x2
                else
                   SLim%Lim_etaij(i,j) = half*(one-dot_product(ni(:,i),ni(:,j)))
                endif
             enddo
             SLim%Lim_etaij(1,i) = half*(one-dot_product(ni(:,i),ni(:,1)))
             SLim%Lim_etaij(2,i) = half*(one-dot_product(ni(:,i),ni(:,2)))
          enddo
          SLim%Lim_etaij(3,4) = half*(one-dot_product(ni(:,3),ni(:,4)))

          SLim%Lim_Ei(1:5) = [sqrts/two,sqrts/two,SLim%AmpMom(1,3),SLim%AmpMom(1,4),e5]
                    
          SLim%wgt = one/(two*mv) &
               * (mv/two)**2 * x1/z1/z2 &
               * ei/s_qi/twopi &
               * two * flux * jac

          SLim%npart = 4

          SLim%Lim_KinInv(1) = e5

          !-- SC
          SCLim%wgt = one/(two*mv) &
               * (mv/two)**2 * x1/z1/z2 &
!               * one/two/twopi &
               * two*z1*z2/(z1+z2-cosi*(z1-z2))**2/twopi &
               * two * flux * jac

          !-- pass e5 and eta5i
          SCLim%Lim_KinInv(1:2) = [e5,x2]

          SCLim%PartFrac = SLim%PartFrac
          SCLim%AmpMom = SLim%AmpMom
          SCLim%Lim_Ei = SLim%Lim_Ei
          SCLim%Lim_etaij = SLim%Lim_etaij
          SCLim%npart = SLim%npart

       endif

       
    endif

  end subroutine kinematics_nlo_z_fs
  
end module mod_kinematics_nlo_z_raoul

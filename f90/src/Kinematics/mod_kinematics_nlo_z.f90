!-- does not include damping factors: for nlo qcd, we do not partition
!-- also, do not set particle types because we don't know whether we are emitting
!   photons or gluons
module mod_kinematics_nlo_z
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

  public :: kinematics_nlo_z_is, kinematics_nlo_z_fs

contains
  
  !===========================================================================!
  !                            Boosted Kinematics                             !
  !===========================================================================!

  !-- both 51 and 52
  subroutine kinematics_nlo_z_is(yr,z1,z2,HardProc,C1Lim,C2Lim,SLim,SC1Lim,SC2Lim,compute_etas)
    use mod_aux_kinematics
    implicit none
    real(dp), intent(in) :: yr(kNLO_max_full),z1,z2
    type(KinConfig) :: HardProc
    logical, intent(in) :: compute_etas
    type(KinConfig), optional :: C1Lim,SC1Lim,C2Lim,SC2Lim,SLim
    real(dp) :: tau,ylab,jac
    real(dp) :: x1,x2,sin5,cos5,phi5
    real(dp) :: xi1,xi2,eta51,eta52,y,c
    real(dp) :: spart,sqrts,mv,mv2
    real(dp) :: p1(4),p2(4),p5(4),pv(4)
    real(dp) :: kallenF,ni(3,5)
    integer :: i,j

    call get_tauy(yr(1:tauy_max),tau,ylab,jac)
    mv2 = tau*sh
    mv  = sqrt(mv2)

    x1   = yr(kNLO_min)         !-- E5     3
    x2   = yr(kNLO_min+1)       !-- eta51  4
    phi5 = twopi*yr(kNLO_min+2) !-- phi5   5

    cos5 = one-two*x2
    sin5 = two*sqrt(abs(x2*(one-x2)))

    !!-----------------------------------------------------------------------!!
    !!                      STRUCTURE 1: Hard Process.                       !!
    !!-----------------------------------------------------------------------!!

    call initialize_config(HardProc)

    HardProc%Lim_etaij(1,2) = one
    HardProc%Lim_etaij(1,5) = x2
    HardProc%Lim_etaij(2,5) = one-x2

    eta51 = HardProc%Lim_etaij(1,5)
    eta52 = HardProc%Lim_etaij(2,5)

    call get_part_s_xi_nlo_z(mv2,ylab,x1,z1,z2,eta51,eta52,spart,sqrts,&
      xi1,xi2,HardProc%flag)

    if (.not. HardProc%flag) then

       HardProc%npart = 5
       HardProc%PartFrac = [xi1,xi2]

       p1 = half*sqrts*z1*[one,zero,zero, one]
       p2 = half*sqrts*z2*[one,zero,zero,-one]
       p5 = half*sqrts*x1*[one,sin5*cos(phi5),sin5*sin(phi5),cos5]

       pv  = p1+p2-p5

       call get_lept_mom(mv,pv,yr(kNLO_min+3:kNLO_max_full),ni(:,3:4),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,5) = p5

       HardProc%Lim_Ei(1) = half*sqrts
       HardProc%Lim_Ei(2) = half*sqrts*x1
       HardProc%Lim_Ei(3) = HardProc%AmpMom(1,3)
       HardProc%Lim_Ei(4) = HardProc%AmpMom(1,4)

       !-- PS: kallenF/8/pi * [x1*s/2 * s/mv2]
       !-- flux: 1/(2*s)
       !-- from FLM: 1/z1/z2
       HardProc%wgt = jac * ((kallenF/32._dp/pi) * (x1*spart/mv2))/z1/z2

       if (compute_etas) then
          ni(:,5) = [sin5*cos(phi5),sin5*sin(phi5),cos5]
          do i = 3,4
             HardProc%Lim_etaij(i,5) = half*(one-dot_product(ni(:,i),ni(:,5)))
          enddo
          c = sqrt(dot_product(ni(:,3)+ni(:,4),ni(:,3)+ni(:,4)))
          y = c**2/((one + (one + c))*((one - c) + one))
          HardProc%Lim_etaij(1,3) = half*(one-ni(3,3))
          HardProc%Lim_etaij(2,3) = half*(one+ni(3,3))
          HardProc%Lim_etaij(1,4) = half*(one-ni(3,4))
          HardProc%Lim_etaij(2,4) = half*(one+ni(3,4))
          HardProc%Lim_etaij(3,4) = one/(one+y)
          HardProc%Lim_KinInv(5)  = y/(one+y) ! 1-eta34
       endif

    endif

    !!-----------------------------------------------------------------------!!
    !!                      STRUCTURE 2: Coll 51 Limit                       !!
    !!-----------------------------------------------------------------------!!
    if (present(C1Lim)) then

       call initialize_config(C1Lim)

       C1Lim%Lim_etaij(1,2) = one
       C1Lim%Lim_etaij(1,5) = zero
       C1Lim%Lim_etaij(2,5) = one

       eta51 = C1Lim%Lim_etaij(1,5)
       eta52 = C1Lim%Lim_etaij(2,5)

       call get_part_s_xi_nlo_z(mv2,ylab,x1,z1,z2,eta51,eta52,spart,sqrts,&
         xi1,xi2,C1Lim%flag)

       if (.not. C1Lim%flag) then

          C1Lim%npart = 4
          C1Lim%PartFrac = [xi1,xi2]

          spart = sh*xi1*xi2
          sqrts  = sqrt(spart)

          p1 = half*sqrts*z1*(one-x1/z1)*[one,zero,zero,one]
          p2 = half*sqrts*z2*[one,zero,zero,-one]

          pv  = p1+p2
          call get_lept_mom(mv,pv,yr(kNLO_min+3:kNLO_max_full),ni(:,3:4),kallenF,C1Lim)

          C1Lim%Lim_Ei(1) = half*sqrts
          C1Lim%Lim_Ei(2) = half*sqrts*x1
          C1Lim%Lim_Ei(3) = C1Lim%AmpMom(1,3)
          C1Lim%Lim_Ei(4) = C1Lim%AmpMom(1,4)

          C1Lim%AmpMom(:,1) = p1
          C1Lim%AmpMom(:,2) = p2

          !-- use Pqg instead of Pqq to get the soft limit in x and not in 1-(1-x)
          C1Lim%Lim_z(1) = x1/z1
          C1Lim%Lim_sij(1,5) = spart*z1*x1*x2

          C1Lim%wgt = jac * ((kallenF/32._dp/pi) * (x1*spart/mv2))/z1/z2

          C1Lim%Lim_KinInv(4) = x2

          if (compute_etas) then
             ni(:,5) = [zero,zero,one]
             do i = 3,4
                C1Lim%Lim_etaij(i,5) = half*(one-dot_product(ni(:,i),ni(:,5)))
             enddo
             c = sqrt(dot_product(ni(:,3)+ni(:,4),ni(:,3)+ni(:,4)))
             y = c**2/((one + (one + c))*((one - c) + one))
             C1Lim%Lim_etaij(1,3) = half*(one-ni(3,3))
             C1Lim%Lim_etaij(2,3) = half*(one+ni(3,3))
             C1Lim%Lim_etaij(1,4) = half*(one-ni(3,4))
             C1Lim%Lim_etaij(2,4) = half*(one+ni(3,4))
             C1Lim%Lim_etaij(3,4) = one/(one+y)
             C1Lim%Lim_KinInv(5)  = y/(one+y) ! 1-eta34
          endif

       endif
    endif

    !!-----------------------------------------------------------------------!!
    !!                      STRUCTURE 3: Coll 52 Limit                       !!
    !!-----------------------------------------------------------------------!!
    if (present(C2Lim)) then

       call initialize_config(C2Lim)

       C2Lim%Lim_etaij(1,2) = one
       C2Lim%Lim_etaij(1,5) = one
       C2Lim%Lim_etaij(2,5) = zero

       eta51 = C2Lim%Lim_etaij(1,5)
       eta52 = C2Lim%Lim_etaij(2,5)

       call get_part_s_xi_nlo_z(mv2,ylab,x1,z1,z2,eta51,eta52,spart,sqrts,&
         xi1,xi2,C2Lim%flag)

       if (.not. C2Lim%flag) then

          C2Lim%npart = 4
          C2Lim%PartFrac = [xi1,xi2]

          spart = sh*xi1*xi2
          sqrts  = sqrt(spart)

          p1 = half*sqrts*z1*[one,zero,zero,one]
          p2 = half*sqrts*z2*(one-x1/z2)*[one,zero,zero,-one]

          pv  = p1+p2
          call get_lept_mom(mv,pv,yr(kNLO_min+3:kNLO_max_full),ni(:,3:4),kallenF,C2Lim)

          C2Lim%Lim_Ei(1) = half*sqrts
          C2Lim%Lim_Ei(2) = half*sqrts*x1
          C2Lim%Lim_Ei(3) = C2Lim%AmpMom(1,3)
          C2Lim%Lim_Ei(4) = C2Lim%AmpMom(1,4)

          C2Lim%AmpMom(:,1) = p1
          C2Lim%AmpMom(:,2) = p2

          !-- use Pqg instead of Pqq to get the soft limit in x and not in 1-(1-x)
          C2Lim%Lim_z(2) = x1/z2
          C2Lim%Lim_sij(2,5) = spart*z2*x1*(one-x2)

          C2Lim%wgt = jac * ((kallenF/32._dp/pi) * (x1*spart/mv2))/z1/z2

          C2Lim%Lim_KinInv(4) = one - x2

          if (compute_etas) then
             ni(:,5) = [zero,zero,-one]
             do i = 3,4
                C2Lim%Lim_etaij(i,5) = half*(one-dot_product(ni(:,i),ni(:,5)))
             enddo
             c = sqrt(dot_product(ni(:,3)+ni(:,4),ni(:,3)+ni(:,4)))
             y = c**2/((one + (one + c))*((one - c) + one))
             C2Lim%Lim_etaij(1,3) = half*(one-ni(3,3))
             C2Lim%Lim_etaij(2,3) = half*(one+ni(3,3))
             C2Lim%Lim_etaij(1,4) = half*(one-ni(3,4))
             C2Lim%Lim_etaij(2,4) = half*(one+ni(3,4))
             C2Lim%Lim_etaij(3,4) = one/(one+y)
             C2Lim%Lim_KinInv(5)  = y/(one+y) ! 1-eta34
          endif

       endif
    endif
       
    !!-----------------------------------------------------------------------!!
    !!                      STRUCTURE 4: Soft 5 Limit                        !!
    !!-----------------------------------------------------------------------!!
    if (present(SLim)) then

       call initialize_config(SLim)

       SLim%Lim_etaij(1,2) = one
       SLim%Lim_etaij(1,5) = x2
       SLim%Lim_etaij(2,5) = one-x2

       eta51 = SLim%Lim_etaij(1,5)
       eta52 = SLim%Lim_etaij(2,5)

       call get_part_s_xi_nlo_z(mv2,ylab,zero,z1,z2,eta51,eta52,spart,sqrts,&
         xi1,xi2,SLim%flag)

       if (.not. SLim%flag) then

          SLim%npart = 4
          SLim%PartFrac = [xi1,xi2]
          
          spart = sh*xi1*xi2
          sqrts  = sqrt(spart)

          p1 = half*sqrts*z1*[one,zero,zero, one]
          p2 = half*sqrts*z2*[one,zero,zero,-one]

          pv  = p1+p2
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

          SLim%Lim_Ei(1) = half*sqrts
          SLim%Lim_Ei(2) = half*sqrts*x1

          SLim%wgt = jac * ((kallenF/32._dp/pi) * (x1*spart/mv2))/z1/z2
          
          !-- now prepare the soft-collinear
          if (present(SC1Lim)) then

             call initialize_config(SC1Lim)
             SC1Lim%npart = SLim%npart

             SC1Lim%Lim_etaij(1,5) = zero
             SC1Lim%Lim_etaij(2,5) = one

             SC1Lim%wgt = jac * ((kallenF/32._dp/pi) * (x1*spart/mv2))/z1/z2

             SC1Lim%Lim_KinInv(4) = x2

          endif

          if (present(SC2Lim)) then
             
             call initialize_config(SC2Lim)
             SC2Lim%npart = SLim%npart

             SC2Lim%Lim_etaij(1,5) = one
             SC2Lim%Lim_etaij(2,5) = zero

             SC2Lim%wgt = jac * ((kallenF/32._dp/pi) * (x1*spart/mv2))/z1/z2

             SC2Lim%Lim_KinInv(4) = one - x2

          endif
       
       endif

    endif
       
  end subroutine kinematics_nlo_z_is

  !-- Both sectors. icoll: lepton collinear to emission; jother: other lepton
  subroutine kinematics_nlo_z_fs(yr,icoll,jother,z1,z2,HardProc,CLim,SLim,SCLim)
    implicit none
    real(dp), intent(in) :: yr(kNLO_max_full),z1,z2
    integer,  intent(in) :: icoll, jother
    type(KinConfig) :: HardProc,CLim,SLim,SCLim
    real(dp) :: tau,ylab,jac
    real(dp) :: xi1,xi2
    real(dp) :: x1,x2,sin5,cos5,phi5
    real(dp) :: slocal,sqrts,mv,mv2,flux
    real(dp) :: p1(4),p2(4),pv(4),p5(4),q(4),pi_aux(4),pi(4),pj(4)
    real(dp) :: cosi,sini,cphi,sphi
    real(dp) :: nicoll_ta(3),nicoll_tb(3),ni(3,5),e5,ei,ej,s_qi
    real(dp) :: si5,z
    integer  :: i,j

    call initialize_config(HardProc)
    call initialize_config(CLim)
    call initialize_config(SLim)
    call initialize_config(SCLim)

    call get_tauy(yr(1:tauy_max),tau,ylab,jac)

    xi1 = sqrt(tau)*exp(+ylab)/z1
    xi2 = sqrt(tau)*exp(-ylab)/z2

    if (xi1.gt.one .or. xi2.gt.one) then

       HardProc%flag = .true.
       CLim%flag     = .true.
       SLim%flag     = .true.
       SCLim%flag    = .true.

    else

       !-- variable assignment
       x1   = yr(kNLO_min)         !-- E5 
       x2   = yr(kNLO_min+1)       !-- eta5
       phi5 = twopi*yr(kNLO_min+2) 
       
       cos5 = one-two*x2
       sin5 = two*sqrt(abs(x2*(one-x2)))

       mv2 = tau * sh
       mv  = sqrt(mv2)

       slocal = mv2/z1/z2
       sqrts  = sqrt(slocal)

       flux = one/2/slocal

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
       pi_aux(:)   = [one,ni(:,icoll)]

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

    !!-----------------------------------------------------------------------!!
    !!                      STRUCTURE 1: Hard Process.                       !!
    !!-----------------------------------------------------------------------!!

       if (ei.lt.zero .or. ej.lt.zero) then

          HardProc%flag = .true.

       else

          HardProc%npart = 5
          HardProc%PartFrac = [xi1,xi2]

          HardProc%AmpMom(:,1) = p1
          HardProc%AmpMom(:,2) = p2
          HardProc%AmpMom(:,icoll)  = pi
          HardProc%AmpMom(:,jother) = pj
          HardProc%AmpMom(:,5) = p5

          !-- prepare etas for the damping
          ni(:,jother) = 1/ej*pj(2:4)
          do i = 1,4
             if (i==icoll) then
                HardProc%Lim_etaij(i,5) = x2
             else
                HardProc%Lim_etaij(i,5) = half*(one-dot_product(ni(:,i),ni(:,5)))
             endif
          enddo

          HardProc%Lim_etaij(1,2) = one
          HardProc%Lim_etaij(1,5) = half*(one - ni(3,5))
          HardProc%Lim_etaij(2,5) = half*(one + ni(3,5))
          HardProc%Lim_Ei(1) = half*sqrts

          HardProc%wgt = slocal/4 * ei/s_qi/twopi * x1 &
                       * flux/z1/z2 * jac

       endif

    !!-----------------------------------------------------------------------!!
    !!                      STRUCTURE 2: Collinear Limit                     !!
    !!-----------------------------------------------------------------------!!
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
          CLim%npart = 4

          CLim%AmpMom(:,1) = p1
          CLim%AmpMom(:,2) = p2
          CLim%AmpMom(:,icoll)  = pi+p5
          CLim%AmpMom(:,jother) = pj

          si5 = four*ei*e5*x2
          z = e5/(e5+ei) !-- to be used with Pqg, so that sing is 1/x1
                         !-- and not 1/[1-(1-x1)] which leads to precision loss

          CLim%Lim_KinInv(1:2) = [z,si5]
          CLim%Lim_Ei(1) = half*sqrts

          CLim%Lim_etaij(1,2) = one
          CLim%Lim_etaij(1,5) = half*(one - ni(3,icoll))
          CLim%Lim_etaij(2,5) = half*(one + ni(3,icoll))

          CLim%Lim_z(1) = z
          CLim%Lim_sij(icoll,5) = si5

          CLim%wgt = slocal/4 * ei/s_qi/twopi * x1 &
                   * flux/z1/z2 * jac
       endif

    !!-----------------------------------------------------------------------!!
    !!                        STRUCTURE 3: Soft Limit                        !!
    !!-----------------------------------------------------------------------!!
       q  = pv
       s_qi = scr(q,pi_aux)
       ei = scr(q,q)/two/s_qi
       pi = ei*pi_aux
       pj = q - pi
       ej = pj(1)

       if (ei.lt.zero .or. ej.lt.zero) then

          SLim%flag = .true.

       else

          !-- S
          SLim%npart = 4
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
          enddo

          SLim%wgt = slocal/4 * ei/s_qi/twopi * x1 &
                   * flux/z1/z2 * jac

          SLim%Lim_KinInv(1) = e5
          SLim%Lim_Ei(1) = half*sqrts
          SLim%Lim_Ei(2) = E5

          !-- SC
          SCLim%wgt = slocal/4 * ei/s_qi/twopi * x1 &
                    * flux/z1/z2 * jac

          SCLim%Lim_etaij(1,2) = one
          SCLim%Lim_etaij(1,5) = half*(one - ni(3,5))
          SCLim%Lim_etaij(2,5) = half*(one + ni(3,5))
          SCLim%Lim_etaij(icoll,5) = x2

          !-- pass e5 and eta5i
          SCLim%Lim_KinInv(1:2) = [e5,x2]

          SCLim%Lim_Ei(1) = half*sqrts
          SCLim%Lim_Ei(2) = E5

       endif

    endif
       
  end subroutine kinematics_nlo_z_fs

end module mod_kinematics_nlo_z

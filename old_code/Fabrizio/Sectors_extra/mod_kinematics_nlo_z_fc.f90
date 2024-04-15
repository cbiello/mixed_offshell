!-- does not include damping factors: for nlo qcd, we do not partition
!-- also, do not set particle types because we don't know whether we are emitting
!   photons or gluons
module mod_kinematics_nlo_z_fc
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

  public :: kinematics_nlo_is_z, kinematics_nlo_fs_z

contains

  !-- tau: mv2/sh, ylab: y_lab
  subroutine get_xi(tau,ylab,z1,z2,x1,x2,xi1,xi2,xi1xi2jac,flag)
    real(dp), intent(in)  :: tau,ylab,z1,z2,x1,x2
    real(dp), intent(out) :: xi1,xi2,xi1xi2jac
    logical, intent(out)  :: flag
    real(dp) :: DS,DP,DM

    flag = .true.
    xi1 = two
    xi2 = two
    xi1xi2jac = -one
    
    DS = z1*z2-z1*x1*x2-z2*x1*(one-x2)
    if (DS.lt.tau) return
    
    DP = z2-x1*x2
    if (DP.lt.zero) return
    
    DM = z1-x1*(one-x2)
    if (DM.lt.zero) return

    xi1 = sqrt(tau/DS) * exp(+ylab) * sqrt(DP/DM)
    xi2 = sqrt(tau/DS) * exp(-ylab) * sqrt(DM/DP)

    if (xi1.gt.one .or. xi2.gt.one) return

    flag = .false.
    xi1xi2jac = one/DS !-- from (xi1,xi2) to (tau,y)

  end subroutine get_xi

  subroutine fill_kininv2(slocal,ni,proc)
    real(dp), intent(in) :: slocal,ni(3,5)
    type(KinConfig) :: proc
    real(dp) :: spart,e3,e4,eta13,eta14,eta23,eta24,eta34
    
    spart = slocal
    e3    = proc%AmpMom(1,3)
    e4    = proc%AmpMom(1,4)
    eta13 = half*(one-dot_product(ni(:,1),ni(:,3)))
    eta14 = half*(one-dot_product(ni(:,1),ni(:,4)))
    eta23 = half*(one-dot_product(ni(:,2),ni(:,3)))
    eta24 = half*(one-dot_product(ni(:,2),ni(:,4)))
    eta34 = half*(one-dot_product(ni(:,3),ni(:,4))) !-- not so sure how Federico is getting this

    proc%Lim_KinInv2 = [spart,e3,e4,eta13,eta14,eta23,eta24,eta34]

  end subroutine fill_kininv2
  
  !-- both 51 and 52
  subroutine kinematics_nlo_is_z(yr,z1,z2,HardProc,C1Lim,C2Lim,SLim,SC1Lim,SC2Lim,compute_etas)
    real(dp), intent(in) :: yr(kNLO_max_full),z1,z2
    type(KinConfig) :: HardProc
    type(KinConfig), optional :: C1Lim,SC1Lim,C2Lim,SC2Lim,SLim
    logical,  intent(in) :: compute_etas
    real(dp) :: tau,ylab,jac
    real(dp) :: x1,x2,sin5,cos5,phi5
    real(dp) :: xi1,xi2,xi1xi2jac
    real(dp) :: slocal,sqrts,mv,mv2
    real(dp) :: p1(4),p2(4),p5(4),pv(4)
    real(dp) :: kallenF,ni(3,5)
    real(dp) :: z,s51,s52
    integer :: i,j
    logical :: xiflag

    call initialize_config(HardProc)

    call get_tauy(yr(1:tauy_max),tau,ylab,jac)

    !-- variable assignment
    x1   = yr(kNLO_min)         !-- E5 3
    x2   = yr(kNLO_min+1)       !-- eta51 4 
    phi5 = twopi*yr(kNLO_min+2) !-- 5

    cos5 = one-two*x2
    sin5 = two*sqrt(abs(x2*(one-x2)))
    
    ni(:,1) = [0,0, 1]
    ni(:,2) = [0,0,-1]

    mv2 = tau*sh
    mv  = sqrt(mv2)
    
    !-- HardProc
    call get_xi(tau,ylab,z1,z2,x1,x2,xi1,xi2,xi1xi2jac,xiflag)
    if (xiflag) then

       HardProc%flag = .true.

    else

       HardProc%PartFrac = [xi1,xi2]

       slocal = sh*xi1*xi2
       sqrts  = sqrt(slocal)

       p1 = half*sqrts*z1*[one,zero,zero, one]
       p2 = half*sqrts*z2*[one,zero,zero,-one]
       p5 = half*sqrts*x1*[one,sin5*cos(phi5),sin5*sin(phi5),cos5]

       pv  = p1+p2-p5
       call get_lept_mom(mv,pv,yr(kNLO_min+3:kNLO_max_full),ni(:,3:4),kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,5) = p5

       !-- PS: kallenF/8/pi * [x1*s/2]
       !-- flux: 1/2/s -->
       !-- wgt = kallenF*x1/32/pi
       HardProc%wgt = kallenF*x1/32._dp/pi * xi1xi2jac * jac/z1/z2
       HardProc%npart = 5

       if (compute_etas) then

          ni(:,5) = [sin5*cos(phi5),sin5*sin(phi5),cos5]
          do i = 3,4
             HardProc%Lim_etaij(i,5) = half*(one-dot_product(ni(:,i),ni(:,5)))
          enddo
         
       endif

       !-- always need these for the damping factors
       HardProc%Lim_etaij(1,5) = x2
       HardProc%Lim_etaij(2,5) = one-x2

       !-- fill variables for subtraction in kininv2
       call fill_kininv2(slocal,ni,HardProc)

       !-- use kininv[5-7] to fill e5,eta51L,eta52L
       HardProc%Lim_KinInv(5:7) = [half*sqrts*x1,x2,one-x2]
      
    endif

    !-- C1
    if (present(C1Lim)) then

       call initialize_config(C1Lim)

       call get_xi(tau,ylab,z1,z2,x1,zero,xi1,xi2,xi1xi2jac,xiflag)
       if (xiflag) then
          
          C1Lim%flag = .true.
          
       else
          
          C1Lim%PartFrac = [xi1,xi2]
          
          slocal = sh*xi1*xi2
          sqrts  = sqrt(slocal)
          
          p1 = half*sqrts*(z1-x1)*[one,zero,zero,one]
          p2 = half*sqrts*z2*[one,zero,zero,-one]
          
          pv  = p1+p2
          call get_lept_mom(mv,pv,yr(kNLO_min+3:kNLO_max_full),ni(:,3:4),kallenF,C1Lim)
          
          C1Lim%AmpMom(:,1) = p1
          C1Lim%AmpMom(:,2) = p2

          C1Lim%LimMom(:,1) = [zero,cos(phi5),sin(phi5),zero] !-- nperp
          C1Lim%LimMom(:,2) = [zero,zero,zero,one]            !-- ndir

          s51 = slocal*x1*x2*z1
          z   = x1/z1 !-- use Pqg instead of Pqq to get the soft limit in x and not in 1-(1-x)
          C1Lim%Lim_KinInv(1:3) = [z,s51,x2]
          
          C1Lim%wgt = kallenF*x1/32._dp/pi * xi1xi2jac * jac/z1/z2
          C1Lim%npart = 4

          !-- always need these for the damping factors
          C1Lim%Lim_etaij(1,5) = zero
          C1Lim%Lim_etaij(2,5) = one

          !-- for ONLO qg subtractions
          C1Lim%Lim_etaij(3,5) = half*(one-ni(3,3))
          C1Lim%Lim_etaij(4,5) = half*(one-ni(3,4))
          
          !-- fill variables for subtraction
          call fill_kininv2(slocal,ni,C1Lim)

          !-- use kininv[6-7] to fill e5,eta51L,eta52L
          C1Lim%Lim_KinInv(5:7) = [half*sqrts*x1,x2,one]
         
       endif
    endif

    !-- C2
    if (present(C2Lim)) then

       call initialize_config(C2Lim)

       call get_xi(tau,ylab,z1,z2,x1,one,xi1,xi2,xi1xi2jac,xiflag)
       if (xiflag) then
          
          C2Lim%flag = .true.
          
       else
          
          C2Lim%PartFrac = [xi1,xi2]
          
          slocal = sh*xi1*xi2
          sqrts  = sqrt(slocal)
          
          p1 = half*sqrts*z1*[one,zero,zero,one]
          p2 = half*sqrts*(z2-x1)*[one,zero,zero,-one]
          
          pv  = p1+p2
          call get_lept_mom(mv,pv,yr(kNLO_min+3:kNLO_max_full),ni(:,3:4),kallenF,C2Lim)
          
          C2Lim%AmpMom(:,1) = p1
          C2Lim%AmpMom(:,2) = p2
          
          s52 = slocal*x1*(one-x2)*z2
          z   = x1/z2 !-- Use Pqq, see comment above
          C2Lim%Lim_KinInv(1:3) = [z,s52,one-x2]

          C2Lim%LimMom(:,1) = [zero,cos(phi5),sin(phi5),zero] !-- nperp
          C2Lim%LimMom(:,2) = [zero,zero,zero,-one]           !-- ndir
          
          C2Lim%wgt = kallenF*x1/32._dp/pi * xi1xi2jac * jac/z1/z2
          C2Lim%npart = 4

          !-- always need these for the damping factors
          C2Lim%Lim_etaij(1,5) = one
          C2Lim%Lim_etaij(2,5) = zero

          !-- for ONLO qg subtractions
          C2Lim%Lim_etaij(3,5) = half*(one+ni(3,3))
          C2Lim%Lim_etaij(4,5) = half*(one+ni(3,4))

          !-- fill variables for subtraction
          call fill_kininv2(slocal,ni,C2Lim)

          !-- use kininv[5-7] to fill e5,eta51L,eta52L
          C2Lim%Lim_KinInv(5:7) = [half*sqrts*x1,one,one-x2]
          
       endif
    endif
       
    !-- S, SC1, SC2
    if (present(SLim)) then

       call initialize_config(SLim)

       call get_xi(tau,ylab,z1,z2,zero,x2,xi1,xi2,xi1xi2jac,xiflag)
       if (xiflag) then
          
          SLim%flag = .true.
          
       else
          
          SLim%PartFrac = [xi1,xi2]
          
          slocal = sh*xi1*xi2
          sqrts  = sqrt(slocal)
          
          p1 = half*sqrts*z1*[one,zero,zero, one]
          p2 = half*sqrts*z2*[one,zero,zero,-one]
          
          pv  = p1+p2
          call get_lept_mom(mv,pv,yr(kNLO_min+3:kNLO_max_full),ni(:,3:4),kallenF,SLim)
          
          SLim%AmpMom(:,1) = p1
          SLim%AmpMom(:,2) = p2
          
          !-- prepare the etas
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
          
          SLim%wgt = kallenF*x1/32._dp/pi * xi1xi2jac * jac/z1/z2
          SLim%npart = 4

          !-- always need these for the damping factors
          SLim%Lim_etaij(1,5) = x2
          SLim%Lim_etaij(2,5) = one-x2

          !-- fill variables for subtraction
          call fill_kininv2(slocal,ni,SLim)

          !-- use kininv[6-7] to fill eta51L,eta52L
          SLim%Lim_KinInv(6:7) = [x2,one-x2]

          
          !-- now prepare the soft-collinar
          if (present(SC1Lim)) then

             call initialize_config(SC1Lim)

             SC1Lim%Lim_KinInv(1:2) = [half*sqrts*x1,x2] !-- E5,eta
             SC1Lim%wgt = kallenF*x1/32._dp/pi * xi1xi2jac * jac/z1/z2
             SC1Lim%npart = SLim%npart

             !-- always need these for the damping factors
             SC1Lim%Lim_etaij(1,5) = zero
             SC1Lim%Lim_etaij(2,5) = one

             SC1Lim%Lim_KinInv2 = SLim%Lim_KinInv2
             SC1Lim%Lim_KinInv(4:7) = [one,zero,x2,one-x2]
             
          endif

          if (present(SC2Lim)) then
             
             call initialize_config(SC2Lim)

             SC2Lim%Lim_KinInv(1:2) = [half*sqrts*x1,one-x2] !-- E5,eta
             SC2Lim%wgt = kallenF*x1/32._dp/pi * xi1xi2jac * jac/z1/z2
             SC2Lim%npart = SLim%npart

             !-- always need these for the damping factors
             SC2Lim%Lim_etaij(1,5) = one
             SC2Lim%Lim_etaij(2,5) = zero

             SC2Lim%Lim_KinInv2 = SLim%Lim_KinInv2
             SC2Lim%Lim_KinInv(4:7) = [zero,one,x2,one-x2]

          endif
       
       endif

    endif
       
  end subroutine kinematics_nlo_is_z

  !-- Both sectors. icoll: lepton collinear to emission; jother: other lepton
  subroutine kinematics_nlo_fs_z(yr,z1,z2,icoll,jother,HardProc,CLim,SCLim,SLim)
    real(dp), intent(in) :: yr(kNLO_max_full),z1,z2
    integer,  intent(in) :: icoll, jother
    type(KinConfig) :: HardProc,CLim,SCLim,SLim
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
    call initialize_config(SCLim)
    call initialize_config(SLim)

    call get_tauy(yr(1:tauy_max),tau,ylab,jac)

    xi1 = sqrt(tau) * exp(+ylab)/z1
    xi2 = sqrt(tau) * exp(-ylab)/z2

    if (xi1.gt.one .or. xi2.gt.one) then

       HardProc%flag = .true.
       CLim%flag     = .true.
       SCLim%flag    = .true.
       SLim%flag     = .true.

    else

       !-- variable assignment
       x1   = yr(kNLO_min)         !-- E5 
       x2   = yr(kNLO_min+1)       !-- eta5
       phi5 = twopi*yr(kNLO_min+2) 
       
       cos5 = one-two*x2
       sin5 = two*sqrt(abs(x2*(one-x2)))

       slocal = sh*xi1*xi2
       sqrts  = sqrt(slocal)

       mv  = sqrts
       mv2 = slocal
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

          !-- prepare etas for the damping
          ni(:,jother) = 1/ej*pj(2:4)
          do i = 1,4
             if (i==icoll) then
                HardProc%Lim_etaij(i,5) = x2
             else
                HardProc%Lim_etaij(i,5) = half*(one-dot_product(ni(:,i),ni(:,5)))
             endif
          enddo

          HardProc%wgt = one/(two*mv) &
               * (mv/two)**2 * x1 &
               * ei/s_qi/twopi &
               * two * flux * jac/z1/z2
          
          HardProc%npart = 5

          !-- fill variables for subtraction in kininv2
          HardProc%Lim_KinInv2(1) = slocal

          !-- use kininv[6-7] to fill eta51L,eta52L
          HardProc%Lim_KinInv(6:7) = [HardProc%Lim_etaij(1,5),HardProc%Lim_etaij(2,5)]
          
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
               * (mv/two)**2 * x1 &
               * ei/s_qi/twopi &
               * two * flux * jac/z1/z2
          
          CLim%npart = 4

          si5 = four*ei*e5*x2
          z = e5/(e5+ei) !-- to be used with Pqg, so that sing is 1/x1
                         !-- and not 1/[1-(1-x1)] which leads to precision loss

          CLim%Lim_KinInv(1:2) = [z,si5]

          !-- fill variables for subtraction in kininv2
          CLim%Lim_KinInv2(1) = slocal

          !-- use kininv[6-7] to fill eta51L,eta52L
          CLim%Lim_KinInv(6) = half*(one-dot_product(ni(:,1),ni(:,icoll)))
          CLim%Lim_KinInv(7) = half*(one-dot_product(ni(:,2),ni(:,icoll)))
          
       endif

       !-- S and CS
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
                    
          SLim%wgt = one/(two*mv) &
               * (mv/two)**2 * x1 &
               * ei/s_qi/twopi &
               * two * flux * jac/z1/z2
          
          SLim%npart = 4

          SLim%Lim_KinInv(1) = e5

          !-- fill variables for subtraction in kininv2
          SLim%Lim_KinInv2(1) = slocal          
          
          !-- use kininv[6-7] to fill eta51L,eta52L
          SLim%Lim_KinInv(6:7) = [SLim%Lim_etaij(1,5),SLim%Lim_etaij(2,5)]
          
          !-- SC
          SCLim%wgt = one/(two*mv) &
               * (mv/two)**2 * x1 &
               * ei/s_qi/twopi &
               * two * flux * jac/z1/z2

          !-- pass e5 and eta5i
          SCLim%Lim_KinInv(1:2) = [e5,x2]

          !-- use kininv[6-7] to fill eta51L,eta52L
          SCLim%Lim_KinInv(6) = half*(one-dot_product(ni(:,1),ni(:,icoll)))
          SCLim%Lim_KinInv(7) = half*(one-dot_product(ni(:,2),ni(:,icoll)))

       endif

       
    endif
       
  end subroutine kinematics_nlo_fs_z
  
end module mod_kinematics_nlo_z_fc

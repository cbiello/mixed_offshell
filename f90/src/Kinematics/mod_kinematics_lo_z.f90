module mod_kinematics_lo_z
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_process
  use mod_auxfunctions
  use mod_kinematics_gen
  use mod_kinematics_lo
  use mod_kinematics_lept

  public :: kinematics_lo_z

contains

  subroutine kinematics_lo_z(yr,z1,z2,LOConfig)
    real(dp), intent(in) :: yr(kLO_max_full)
    type(KinConfig)   :: LOConfig
    real(dp) :: z1,z2
    real(dp) :: tau,ylab,jac,kallenf
    real(dp) :: mv,mv2,EClocal
    real(dp) :: p1(4),p2(4),p12(4),nlept(3,2)
    real(dp) :: spart,e3,e4,eta13,eta14,eta23,eta24,eta34,n1(3),n2(3)
    real(dp) :: y,c

    call initialize_config(LOConfig)

    call get_tauy(yr(1:tauy_max),tau,ylab,jac)
    LOconfig%PartFrac = sqrt(tau)*[exp(+ylab)/z1,exp(-ylab)/z2]

    if (any(LOConfig%partFrac >= one)) then
       LOConfig%flag = .true.
       return
    endif

    mv2 = tau*sh
    mv = sqrt(mv2)
    spart = mv2/z1/z2
    EClocal = half*sqrt(spart)

    p1 = EClocal*z1*[one,zero,zero, one] 
    p2 = EClocal*z2*[one,zero,zero,-one] 
    p12 = p1 + p2

    call get_lept_mom(MV,p12,yr(kLO_min:kLO_max_full),nlept,kallenF,LOConfig)
    LOConfig%AmpMom(:,1) = p1
    LOConfig%AmpMom(:,2) = p2

    !1/(two*mv2) is flux factor, extra z1,z2 from the jacobian
    LOConfig%wgt  = one/8.0_dp/pi*kallenF/(two*mv2)*jac/z1/z2

    LOConfig%npart  = 4

    allocate(LOConfig%part(LOConfig%npart))
    
    LOConfig%mu2ref = mv2 !-- for hoppet

    !-- extra parameters I need
    n1 = [zero,zero,+one]
    n2 = [zero,zero,-one]

    E3    = LOConfig%AmpMom(1,3)
    E4    = LOConfig%AmpMom(1,4)
    eta13 = half*(one-dot_product(n1,nlept(:,1)))
    eta14 = half*(one-dot_product(n1,nlept(:,2)))
    eta23 = half*(one-dot_product(n2,nlept(:,1)))
    eta24 = half*(one-dot_product(n2,nlept(:,2)))
    eta34 = half*(one-dot_product(nlept(:,1),nlept(:,2)))

    LOConfig%Lim_Ei(1) = EClocal
    LOConfig%Lim_Ei(2) = EClocal
    LOConfig%Lim_Ei(3) = E3
    LOConfig%Lim_Ei(4) = E4

    LOConfig%Lim_etaij(1,3) = half*(one-nlept(3,1))
    LOConfig%Lim_etaij(2,3) = half*(one+nlept(3,1))
    LOConfig%Lim_etaij(1,4) = half*(one-nlept(3,2))
    LOConfig%Lim_etaij(2,4) = half*(one+nlept(3,2))

    c = sqrt(dot_product(nlept(:,1)+nlept(:,2),nlept(:,1)+nlept(:,2)))
    y = c**2/((one + (one + c))*((one - c) + one))
    LOConfig%Lim_etaij(3,4) = one/(one+y)
    LOConfig%Lim_KinInv(5)  = y/(one+y) ! 1-eta34, stable version

    LOConfig%Lim_KinInv2 = [spart,e3,e4,eta13,eta14,eta23,eta24,eta34]

  end subroutine kinematics_lo_z

end module mod_kinematics_lo_z

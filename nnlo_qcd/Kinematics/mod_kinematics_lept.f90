module mod_kinematics_lept
  use mod_types; use mod_consts_dp; use common_def
  use mod_auxfunctions
  use mod_process
  implicit none
  private

  public :: get_lept_mom

contains

  subroutine get_lept_mom(MV,pV,yr,jac,ProcConfig)
    real(dp), intent(in) :: MV,pV(1:4),yr(3)
    type(KinConfig) :: ProcConfig
    real(dp), intent(out) :: jac
    integer         :: numlept,numleptrand,numpart,i
    logical :: flag

    jac=zero
    numlept = ProcConfig%NumLept
    numpart = size(ProcConfig%AmpMom,dim=2)-numlept

#if (_PhaseSpaceType == 3)
    jac = 16*pisq/sh !-- 16 Pi: we should have 2Pi delta(mhsq-sh x1 x2) 
                     !--        we write this as 1/8 Pi --> 16 pisq
                     !--        1/sh is jacobian for the delta. Works for on-shell Higgs
    ProcConfig%AmpMom(:,numpart+1) = pV
#else
    call get_DY_lept(MV,zero,zero,pV,yr(:),jac,ProcConfig%AmpMom(1:4,numpart+1:numpart+numlept))
#endif



  end subroutine get_lept_mom

  !----------------------------------------------------------------------


  ! ------------------------------------------------------------------------------------
  !-- DY

  !-- the lepton momenta are generated in the CMS frame
  subroutine  get_DY_lept(MV,m6,m7,p124,yR,kallenF,plept)
    real(dp), intent(in) :: MV, m6,m7,p124(1:4),yR(1:3)
    real(dp), intent(out) :: plept(1:4,1:2),kallenF
    real(dp) :: MV2,E6,E7,p67,cos6,sin6,phi6,p6aux(1:4),p7aux(1:4),x6,x7
    real(dp) :: pV(4),phi, pVt,pVt2,pt(4),Delta,pta,YV,ptmax,y

    x6 = yR(1)
    x7 = yR(2)

    pt = zero
    MV2 = MV**2

    !   kinematics of the V

    pV = p124

    YV = half*log( abs ( ( pV(1) + pV(4))/(pV(1) - pV(4) ) ) )

    pVt2 = pV(2)**2+pV(3)**2
    pVt = sqrt(pVt2)

    !   lepton
    phi = two*pi*x7

    ptmax = half*MV2/( sqrt(MV2 + pVt2) - ( pV(2)*cos(phi) + pV(3)*sin(phi))  )

    pta = ptmax*x6

    pt =  pta*(/zero, cos(phi), sin(phi), zero/)

    Delta = (  MV2 + two*( pV(2)*pt(2) + pV(3)*pt(3) ) )/two/pta/sqrt(MV2 + pVt2)

    if (yR(3).ge.half) then
       y = YV - acosh(Delta)
    else
       y = YV + acosh(Delta)
    endif

    plept(1:4,1) = (/pta*cosh(y),pta*cos(phi),pta*sin(phi),pta*sinh(y)/)
    plept(1:4,2) = pV - plept(1:4,1)

    kallenF = two * ptmax/sqrt(mV2 + pVt**2)/(abs(sinh(YV-y)))

  end subroutine get_DY_lept


end module mod_kinematics_lept

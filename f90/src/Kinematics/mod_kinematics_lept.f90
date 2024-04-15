module mod_kinematics_lept
  use mod_types
  use mod_consts_dp
  use mod_process
  implicit none
  private

  public :: get_lept_mom

contains
  
  subroutine get_lept_mom(mv,pV,yr,nlept,jac,ProcConfig)
    real(dp), intent(in) :: mv,pV(1:4),yr(3)
    type(KinConfig) :: ProcConfig
    real(dp), intent(out) :: jac,nlept(1:3,1:2)

    jac=zero

    call get_DY_lept(mv,pV,yr(:),nlept,jac,ProcConfig%AmpMom(1:4,3:4))

  end subroutine get_lept_mom

  !----------------------------------------------------------------------


  ! ------------------------------------------------------------------------------------
  !-- DY

  !-- the lepton momenta are generated in the CMS frame
  subroutine  get_DY_lept(mv,pv,yR,nlept,kallenF,plept)
    real(dp), intent(in) :: mv,pv(1:4),yR(1:3)
    real(dp), intent(out) :: plept(1:4,1:2),kallenF,nlept(1:3,1:2)
    real(dp) :: x1,x2
    real(dp) :: mv2,phi,pVt,pVt2,pt(4),Delta,pta,yv,ptmax,y
    real(dp) :: jacpt

    x1 = yR(1)
    x2 = yR(2)

    pt = zero
    mv2 = mv**2

    !-- kinematics of the V
    yv = half*log(abs((pV(1)+pV(4))/(pV(1)-pV(4))))

    pVt2 = pV(2)**2+pV(3)**2
    pVt = sqrt(pVt2)

    !-- leptons
    phi = twopi*x2

    ptmax = half*mv2/(sqrt(mv2+pVt2)-(pV(2)*cos(phi)+pV(3)*sin(phi)))

    pta   = ptmax*x1
    jacpt = ptmax

    pt =  pta*(/zero,cos(phi),sin(phi),zero/)

    Delta = (MV2+two*(pV(2)*pt(2)+pV(3)*pt(3)))/two/pta/sqrt(MV2+pVt2)

    if (yR(3).ge.half) then
       y = yv - acosh(Delta)
    else
       y = yv + acosh(Delta)
    endif

    plept(1:4,1) = (/pta*cosh(y),pta*cos(phi),pta*sin(phi),pta*sinh(y)/)
    plept(1:4,2) = pV - plept(1:4,1)

    nlept(1:3,1) = (/plept(2,1), plept(3,1), plept(4,1)/)/plept(1,1)
    nlept(1:3,2) = (/plept(2,2), plept(3,2), plept(4,2)/)/plept(1,2)

    kallenF = two * jacpt/sqrt(mV2 + pVt**2)/(abs(sinh(YV-y)))

  end subroutine get_DY_lept
  
end module mod_kinematics_lept


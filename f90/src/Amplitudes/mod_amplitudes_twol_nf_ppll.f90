!-- amplitudes with extra jet: factor out gs**2 or ee**2 per extra gluon/photon emission
module mod_amplitudes_twol_nf_ppll
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_auxfunctions
  use mod_coupl
  use handyG
  implicit none
  private
  complex(dp), save :: dZe,dZaa,SiTAA_light_Z,dSiTAA_heavy_0
  complex(dp), save :: dZaz,dZza,dZzz,dZmz2,dZmw2,deltaR11
  real(dp), parameter :: munf = 200._dp !-- results will not depend on munf
  complex(kind=prec), parameter :: zp=(zero,zero), op=(one,zero), tp=(two,zero)

  !-- res_tree_emitted-particles_initial-state
  !-- q stands for both q and qb
  !-- res_tree   --> 0 -> q e- e+ qb
  
  public :: ew_renorm_nf_2L,res_twol_nf_qqb

contains

  subroutine ew_renorm_nf_2L()
    implicit none
    complex(dp15) :: SiTAZz,dSiTZZz,SiTZZz,SiTWWw,dSiTWWw,SiTWW0
    real(dp) :: mw2,mz2,mu2

    !-- All counterterms are computed using real values for the masses
    !-- The couplings, e.g. cw and sw are complex, following OpenLoop convention.
    !-- The only exceptions are in the Z,W mass counterterms dZmz2 and dZmw2
    !-- where there's an extra term coming from the complex-mass scheme

    mz2 = real(cms_mzsq,kind=dp)
    mw2 = real(cms_mwsq,kind=dp)
    mu2 = munf**2

    !-- SiTAZ(mz2)
    SiTAZz  = xn * CF * (-Qup * vcu * mtsq * F(mz2,mtsq,mu2) &
                        -(nup*Qup*vcu + ndn*Qdn*vcd) * mz2 * A(mz2,mu2))

    !-- SiTZZ(mz2)
    SiTZZz  = xn * CF * (acu**2 * mtsq * X(mz2,mtsq,mu2) + (acu**2 + vcu**2) * mtsq * F(mz2,mtsq,mu2) &
                      + (nup * (acu**2 + vcu**2) + ndn * (acd**2 + vcd**2)) * mz2 * A(mz2,mu2))

    !-- SiTWW(mw2)
    SiTWWw  = xn * CF * (mtsq*E(mw2,mtsq,mu2) + nlq * mw2 * (55._dp/12 - 4*zeta3+log(mu2/mw2)+ci*pi))/cms_sw2

    !-- SiTWW(0)
    SiTWW0  = xn * CF * mtsq/cms_sw2 * (5*pisq/12 + 11*log(mu2/mtsq)/2 + 3*log(mu2/mtsq)**2 + 35._dp/8)

    !-- d/ds SiTZZ(s) |s=mz2
    dSiTZZz = xn * CF * (acu**2 * mtsq * dX(mz2,mtsq) + (acu**2 + vcu**2) * mtsq * dF(mz2,mtsq,mu2) &
                      + (nup * (acu**2 + vcu**2) + ndn * (acd**2 + vcd**2)) * (A(mz2,mu2)-four))

    !-- d/ds SiTWW(s) |s=mw2
    dSiTWWw = xn * CF * (one/cms_sw2*dE(mw2,mtsq,mu2)+nlq/cms_sw2*(43._dp/12 + ci*pi + log(mu2/mw2) - 4*zeta3))

    dSiTAA_heavy_0 = xn * CF * (Qup2 * (15._dp + 4*log(mu2/mtsq)))
    SiTAA_light_Z  = xn * CF * (nup*Qup2 + ndn*Qdn2) * A(mz2,mu2) * mz2

    dZaa = - dSiTAA_heavy_0 - SiTAA_light_Z/mz2

    !-- Eqs.(2.37,2.38) of 2009.02229
    dZza  = zero
    dZaz  = -2/mz2 * SiTAZz
    dZzz  = -dSiTZZz
    dZmz2 = SiTZZz + (cms_mzsq - mz2)*dSiTZZz
    dZmw2 = SiTWWw + (cms_mwsq - mw2)*dSiTWWw

    DeltaR11 = - dZaa - cms_cw2/cms_sw2*(dZmz2/cms_mzsq - dZmw2/cms_mwsq) + (SiTWW0-dZmw2)/cms_mwsq

    if(ew_scheme.eq.'Gmu') then
      dZe = -dZaa/2 - DeltaR11/2
    elseif(ew_scheme.eq.'amz') then
      dZe = -dZaa/2
    endif

    dZe = real(dZe,kind=dp)

  end subroutine ew_renorm_nf_2L

  !----------------------------------------------------------------------
  !-- 4-point amplitudes
  !----------------------------------------------------------------------

  !-- res(1,1) = d dx -> l lb
  !-- res(2,1) = dx d -> l lb
  !-- res(1,2) = u ux -> l lb
  !-- res(2,2) = ux u -> l lb
  subroutine res_twol_nf_qqb(p,res0,res2)
    implicit none
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res0(2,2),res2(2,2)
    integer  :: i
    real(dp15) :: ss,tt,uu
    complex(dp15) :: K11(3)
    real(dp) :: Qq(2),Iq(2)
    complex(dp) :: mz2,vq(2),aq(2)
    complex(dp) :: cvv(2),cva(2),cav(2),caa(2)
    complex(dp) :: HDY00(4,2,2)
    complex(dp) :: HDY11(4,2,2)

    Qq = [Qdn,Qup]
    Iq = [I3dn,I3up]
    vq = [vcd,vcu]
    aq = [acd,acu]

    ss = two*scr(p(:,1),p(:,2))
    mz2 = cms_mzsq

    call evaluate_K11s(ss,K11)

    do i = 1,2

      cvv(i) = -Qel*Qq(i) * K11(1)/ss - (Qq(i)*vce + Qel*vq(i))/(ss-mz2) * K11(2) - vq(i)*vce/(ss-mz2) * K11(3)

      cav(i) = aq(i)*Qel/(ss-mz2) * K11(2) + aq(i)*vce/(ss-mz2) * K11(3)

      cva(i) = ace*Qq(i)/(ss-mz2) * K11(2) + ace*vq(i)/(ss-mz2) * K11(3)

      caa(i) = -aq(i)*ace/(ss-mz2) * K11(3)
    end do

    !! qqx channel
    tt = -two*scr(p(:,1),p(:,3))
    uu = -ss-tt

    do i = 1, 2
      !-- -+-+
      HDY00(1,i,1) = -2*Qq(i)*uu/ss-(one/6/cms_cw2+Iq(i)/cms_sw2-2*Qq(i))*uu/(ss-mz2)
      HDY11(1,i,1) = 2*uu*(cvv(i)-cav(i)-cva(i)+caa(i))
      !-- +--+
      HDY00(2,i,1) = 2*Qq(i)*tt/ss - Qq(i)*(two - 1/cms_cw2)*tt/(ss - mz2)
      HDY11(2,i,1) = -2*tt*(cvv(i)+cav(i)-cva(i)-caa(i))
      !-- -++-
      HDY00(3,i,1) = 2*Qq(i)*tt/ss + (one/3/cms_cw2 - 2*Qq(i))*tt/(ss - mz2)
      HDY11(3,i,1) = -2*tt*(cvv(i)-cav(i)+cva(i)-caa(i))
      !-- +-+-
      HDY00(4,i,1) = -2*Qq(i)*uu/ss - 2*Qq(i)*(cms_sw2/cms_cw2)*uu/(ss - mz2)
      HDY11(4,i,1) = 2*uu*(cvv(i)+cav(i)+cva(i)+caa(i))
    end do

    !! qxq channel
    tt = -two*scr(p(:,1),p(:,4))
    uu = -ss-tt
    do i = 1, 2
      !-- -+-+
      HDY00(1,i,2) = -2*Qq(i)*uu/ss-(one/6/cms_cw2+Iq(i)/cms_sw2-2*Qq(i))*uu/(ss-mz2)
      HDY11(1,i,2) = 2*uu*(cvv(i)-cav(i)-cva(i)+caa(i))
      !-- +--+
      HDY00(2,i,2) = 2*Qq(i)*tt/ss - Qq(i)*(two - 1/cms_cw2)*tt/(ss - mz2)
      HDY11(2,i,2) = -2*tt*(cvv(i)+cav(i)-cva(i)-caa(i))
      !-- -++-
      HDY00(3,i,2) = 2*Qq(i)*tt/ss + (one/3/cms_cw2 - 2*Qq(i))*tt/(ss - mz2)
      HDY11(3,i,2) = -2*tt*(cvv(i)-cav(i)+cva(i)-caa(i))
      !-- +-+-
      HDY00(4,i,2) = -2*Qq(i)*uu/ss - 2*Qq(i)*(cms_sw2/cms_cw2)*uu/(ss - mz2)
      HDY11(4,i,2) = 2*uu*(cvv(i)+cav(i)+cva(i)+caa(i))
    end do

    res0(1,1) = real(dot_product(HDY00(:,1,1),HDY00(:,1,1)),kind=dp)
    res0(2,1) = real(dot_product(HDY00(:,1,2),HDY00(:,1,2)),kind=dp)
    res0(1,2) = real(dot_product(HDY00(:,2,1),HDY00(:,2,1)),kind=dp)
    res0(2,2) = real(dot_product(HDY00(:,2,2),HDY00(:,2,2)),kind=dp)

    !-- HDY11 are expanded in as/4/pi,aw/4/pi
    !-- thus: 2*Re[(as/4/pi)*(aw/4/pi)*A0*A2] = (as/2/pi)*(aw/2/pi)*Re[A0*A2]/2
    res2(1,1) = real(dot_product(HDY00(:,1,1),HDY11(:,1,1)),kind=dp)/2
    res2(2,1) = real(dot_product(HDY00(:,1,2),HDY11(:,1,2)),kind=dp)/2
    res2(1,2) = real(dot_product(HDY00(:,2,1),HDY11(:,2,1)),kind=dp)/2
    res2(2,2) = real(dot_product(HDY00(:,2,2),HDY11(:,2,2)),kind=dp)/2

    res0 = res0 * xn * eesq2 * aveqq
    res2 = res2 * xn * eesq2 * aveqq

  end subroutine res_twol_nf_qqb

  subroutine evaluate_K11s(s,K11)
     implicit none
     real(dp15), intent(in) :: s
     complex(dp15), intent(out) :: K11(3)
     real(dp) :: mu2
     complex(dp15) :: As,Fs,Xs
     complex(dp15) :: SiTAA,SiTAZ,SiTZZ

     mu2 = munf**2

     As = A(s,mu2)
     Fs = F(s,mtsq,mu2)
     Xs = X(s,mtsq,mu2)

     !-- Evaluate VV self-energies
     SiTAA = xn * CF * (Qup2 * mtsq * Fs + (nup*Qup2 + ndn*Qdn2) * s * As)
     SiTAZ = xn * CF * (-Qup * vcu * mtsq * Fs - (nup*Qup*vcu + ndn*Qdn*vcd) * s * As)
     SiTZZ = xn * CF * (acu**2 * mtsq * Xs + (acu**2 + vcu**2) * mtsq * Fs &
                     + (nup * (acu**2 + vcu**2) + ndn * (acd**2 + vcd**2)) * s * As)

     K11(1) = SiTAA/s - 2*dZe
     K11(2) = cms_cw/cms_sw * (dZmz2/cms_mzsq - dZmw2/cms_mwsq) - SiTAZ/s
     K11(3) = (SiTZZ-dZmz2)/(s-cms_mzsq) - 2*dZe + (cms_cw2-cms_sw2)/cms_sw2 * (dZmz2/cms_mzsq - dZmw2/cms_mwsq)

  end subroutine evaluate_K11s

  function A(s,mu2)
    implicit none
    real(dp), intent(in) :: s,mu2
    complex(dp) :: A

    A = 55._dp/3 + 4*ci*pi - 16*zeta3 + 4*log(mu2/s)

  end function
  
  function F(s,mt2,mu2)
    implicit none
    real(dp), intent(in) :: s,mt2,mu2
    complex(dp) :: F
    complex(dp) :: xt,t(8)

    if (s.lt.four*mt2) then
       xt = (s + ci*sqrt(-s*(s-4*mt2)))/(2*mt2)
    else
       xt = (s + sqrt(s*(s-4*mt2)))/(2*mt2)
    endif

    t(1) = G([op], xt)
    t(2) = G([op,op], xt)
    t(3) = G([zp,op], xt)
    t(4) = G([tp,op], xt)
    t(5) = G([zp,op,op], xt)
    t(6) = G([op,zp,op], xt)
    t(7) = G([op,tp,op], xt)
    t(8) = G([tp,op,op], xt)

    F = -(55*xt**4+104*xt**2*(xt-1._dp)+32*(t(3)+2*t(4))*xt*(xt-2._dp)*(2*xt+xt**2-2._dp)+&
      32*(t(5)-t(6)-2*t(7)+2*t(8))*(8*xt-4*xt**2+xt**4-4._dp)-12*t(1)*xt*(xt-2._dp)*(6*xt+xt**2-&
      6._dp)-16*t(2)*(xt-1._dp)*(-19*xt+4*xt**2+6*xt**3+7._dp))/(3*xt**2*(-xt+1._dp))+&
      (4*xt**2*log(mu2/mt2))/(xt-1._dp)

  end function

  function E(s,mt2,mu2)
    implicit none
    real(dp), intent(in) :: s,mt2,mu2
    complex(dp) :: E
    real(dp) :: y
    real(kind=prec):: ze=zero,on=one

    y = s/mt2

    E = -(((-one + y)*(pisq*(-two - 2*y + 4*y**2) &
        + 3*(-5._dp - 9*y + 6*y**2))*G([on],y))/(18*y**2)) &
        - (2*(-one + y + 2*y**2)*G([ze,on], y))/(3*y) &
        + ((-one + y)**2*(five + 4*y)*G([on,on], y))/(3*y**2) &
        + (two - 2/(3*y**2) - (4*y)/3)*G([ze,on,on], y) &
        + two/3*(-3._dp + 1/y**2 + 2*y)*G([on,ze,on], y) &
        + (-60._dp + 213*y + 330*y**2 + 2*pisq*(-four + 13*y + 8*y**2) &
        + 36*y*(11._dp + 2*y)*log(mu2/mt2) &
        + 216*y*log(mu2/mt2)**2)/(72*y)

  end function

  function X(s,mt2,mu2)
    implicit none
    real(dp), intent(in) :: s,mt2,mu2
    complex(dp) :: X
    complex(dp) :: xt,t(8)

    if (s.lt.four*mt2) then
       xt = (s + ci*sqrt(-s*(s-4*mt2)))/(2*mt2)
    else
       xt = (s + sqrt(s*(s-4*mt2)))/(2*mt2)
    endif

    t(1) = G([op], xt)
    t(2) = G([op,op], xt)
    t(3) = G([zp,op], xt)
    t(4) = G([tp,op], xt)
    t(5) = G([zp,op,op], xt)
    t(6) = G([op,zp,op], xt)
    t(7) = G([op,tp,op], xt)
    t(8) = G([tp,op,op], xt)

    X = 88*log(mu2/mt2)+48*log(mu2/mt2)**2+2*(2*Pi**2+(8*(xt-1._dp))/xt**2-11._dp+(4*(xt-&
      2._dp)*(-2*xt+9*xt**2+2._dp)*t(1))/xt**3+(8*(4*xt+2*xt**2-16*xt**3+9*xt**4-&
      2._dp)*t(2))/xt**4-(32*(xt-2._dp)*(t(3)+2*t(4)))/xt-(32*(-2*xt+xt**2+&
      2._dp)*(t(5)-t(6)-2*t(7)+2*t(8)))/xt**2)

  end function

  function dF(s,mt2,mu2)
    implicit none
    real(dp), intent(in) :: s,mt2,mu2
    complex(dp) :: dF
    complex(dp) :: xt,t(8)

    if (s.lt.four*mt2) then
       xt = (s + ci*sqrt(-s*(s-4*mt2)))/(2*mt2)
    else
       xt = (s + sqrt(s*(s-4*mt2)))/(2*mt2)
    endif

    t(1) = G([op], xt)
    t(2) = G([op,op], xt)
    t(3) = G([zp,op], xt)
    t(4) = G([tp,op], xt)
    t(5) = G([zp,op,op], xt)
    t(6) = G([op,zp,op], xt)
    t(7) = G([op,tp,op], xt)
    t(8) = G([tp,op,op], xt)

    dF = (4*log(mu2/mt2))/mt2+((-72*(xt-1._dp))/xt**2+43._dp-(4*(xt-2._dp)*(-10*xt+3*xt**2+&
      10._dp)*t(1))/xt**3+(16*(xt-1._dp)*(-136*xt+135*xt**2-55*xt**3+14*xt**4+&
      44._dp)*t(2))/(xt**4*(xt-2._dp)**2)-(64*(xt-1._dp)*(xt-2._dp)*(t(3)+&
      2*t(4)))/xt**3+(32*(-8*xt+4*xt**2+xt**4+4._dp)*(t(5)-t(6)-2*t(7)+&
      2*t(8)))/xt**4)/(3*mt2)

  end function

  function dX(s,mt2)
    implicit none
    real(dp), intent(in) :: s,mt2
    complex(dp) :: dX
    complex(dp) :: xt,t(8)

    if (s.lt.four*mt2) then
       xt = (s + ci*sqrt(-s*(s-4*mt2)))/(2*mt2)
    else
       xt = (s + sqrt(s*(s-4*mt2)))/(2*mt2)
    endif

    t(1) = G([op], xt)
    t(2) = G([op,op], xt)
    t(3) = G([zp,op], xt)
    t(4) = G([tp,op], xt)
    t(5) = G([zp,op,op], xt)
    t(6) = G([op,zp,op], xt)
    t(7) = G([op,tp,op], xt)
    t(8) = G([tp,op,op], xt)

    dX = (8*(xt-1._dp)*(xt*(-4*xt+9*xt**2+4._dp)-2*(xt-2._dp)*(-2*xt+3*xt**2+2._dp)*t(1)-&
      (8*(xt-1._dp)*(-8*xt-3*xt**2+19*xt**3-12*xt**4+3*xt**5+4._dp)*t(2))/(xt*(xt-&
      2._dp)**2)+8*xt**2*(xt-2._dp)*(t(3)+2*t(4))-16*xt*(xt-1._dp)*(t(5)-t(6)&
      -2*t(7)+2*t(8))))/(mt2*xt**5)

  end function

  function dE(s,mt2,mu2)
    implicit none
    real(dp), intent(in) :: s,mt2,mu2
    complex(dp) :: dE
    real(dp) :: y
    real(kind=prec):: ze=zero,on=one

    y = s/mt2

    dE = (60._dp + 8*pisq + 54*y + 4*pisq*y + 129*y**2)/(36*y**2) &
       - ((-one + y)*(15._dp + 2*pisq + 12*y + 2*pisq*y + 9*y**2 &
       + 2*pisq*y**2)*G([on], y))/(9*y**3) - (2*(two + y)*G([ze,on], y))/(3*y**2) &
       + (2*(-one + y)*(two + y)*G([on,on], y))/y**3 &
       - (4*(-one + y)*(one + y + y**2)*G([ze,on,on], y))/(3*y**3) &
       + (4*(-one + y)*(one + y + y**2)*G([on,ze,on], y))/(3*y**3) + log(mu2/mt2)

  end function

end module mod_amplitudes_twol_nf_ppll

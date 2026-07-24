module mod_int_sub_nnlo
  use types
  use consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_auxfunctions
  implicit none
  private

  public :: fin_ns_vqcd
  public :: Pqqbqqb, Pqqbqqb_Lmu
  public :: calG_QqQl,calG_QqQl_L,calG_CF,calG_Q2
  public :: calG_ONLOQCD_ns

  interface calG_CF
    module procedure calG_CF_hard, calG_CF_coll
  end interface

contains

  !!=========================================================================!!

  function fin_ns_vqcd(etasb,Qq,Ql) result(res)
     implicit none
     real(dp), intent(in) :: etasb(:)
     real(dp), intent(in) :: Qq(:),Ql
     real(dp) :: a(2),res(4)
     real(dp) :: li2eta13,li2eta14
     integer  :: i

     li2eta13 = real(dilog2(etasb(1)),kind=dp)
     li2eta14 = real(dilog2(etasb(2)),kind=dp)

     a(1) = 4*(li2eta13 - li2eta14) - 6*log(etasb(1)/etasb(2))
     a(2) = 13._dp - 4*zeta2

     res = zero
     do i = 1, size(Qq)
       res(i) = res(i) + a(1)*Qq(i)*Ql
       res(i) = res(i) + a(2)*Ql**2
     end do

  end function fin_ns_vqcd

  function Pqqbqqb(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res, xp,xb
    real(dp) :: lx,lxp,lxb,li2x,li2mx,li3x,li3xb,li3mx

    xb = one-x
    xp = one+x
    lx = log(x)
    lxp = log(xp)
    lxb = log(xb)
    li2x = real(dilog2(x),kind=dp) 
    li2mx = real(dilog2(-x),kind=dp)
    li3x = trilog(x)
    li3xb = trilog(xb)
    li3mx = trilog(-x)

    !-- reg
    res = 7*x+6*lx*lxb*xb+(8*li2mx+8*lx*lxp)*xp+1._dp+((16*li3mx+18*li3x+12*li3xb-8*li2mx*lx&
      -10*li2x*lx-lx**3/3+8*li2x*lxb-5*lx**2*lxb+8*lx*lxb**2-4*lx*zeta2-8*lxb*zeta2-&
      6*zeta3)*(x**2+1._dp))/xb-(lx**2*(2*x-1._dp)*(2*x-5._dp))/(2*xb)+(lx*(-11*x+27*x**2-&
      6._dp))/xb+4*lxb*(7*x-8._dp)-(2*zeta2*(-6*x+x**2+11._dp))/xb-(2*li2x*(6*x+x**2-13._dp))/xb

  end function Pqqbqqb

  function Pqqbqqb_Lmu(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: res,xb
    real(dp) :: lx,lxb,li2x

    xb = one-x
    lx = log(x)
    lxb = log(xb)
    li2x = real(dilog2(x),kind=dp)

    !-- reg
    res = ((-4*li2x+2*lx**2-4*lx*lxb+4*zeta2)*(x**2+1._dp))/xb-(2*lx*(2*x**2-5._dp))/xb-2*(7*x-8._dp)

  end function Pqqbqqb_Lmu

  function calG_QqQl(proc,z,Qcharges,Qlept,icoll,chrg)
    use mod_process, only: KinConfig
    implicit none
    type(KinConfig),   intent(in) :: proc
    real(dp),          intent(in) :: z,Qcharges(:),Qlept
    integer, optional, intent(in) :: icoll, chrg
    real(dp) :: calG_QqQl(size(Qcharges))
    real(dp) :: EC,E3,E4
    real(dp) :: etai3,etai4,eta13,eta14,eta23,eta24,eta34,eta34b
    real(dp) :: logEC3,logEC4,logE34,li213,li214,li223,li224,li234
    real(dp) :: a,b,c
    integer  :: mychrg

    EC = proc%Lim_Ei(1)
    E3 = proc%Lim_Ei(3)
    E4 = proc%Lim_Ei(4)

    eta13  = proc%Lim_etaij(1,3)
    eta23  = proc%Lim_etaij(2,3)
    eta14  = proc%Lim_etaij(1,4)
    eta24  = proc%Lim_etaij(2,4)
    eta34  = proc%Lim_etaij(3,4)
    eta34b = proc%Lim_KinInv(5)

    logEC3 = log(EC/E3)
    logEC4 = log(EC/E4)
    logE34 = log(E3/E4)

    li234 = real(dilog2(eta34b),kind=dp)
    li213 = real(dilog2(one-eta13),kind=dp)
    li214 = real(dilog2(one-eta14),kind=dp)
    li223 = real(dilog2(one-eta23),kind=dp)
    li224 = real(dilog2(one-eta24),kind=dp)

    a = 4*zeta2
    b = 2*(li213-li214-li223+li224) &
      + (3._dp + 2*logEC3)*log(eta13/eta23) + (3._dp + 2*logEC4)*log(eta24/eta14)            ! eq 2.67, Geq^ij * two for some reason
    c = 13._dp - 4*zeta2 + logE34**2 + (3._dp + 2*logEC3 + 2*logEC4)*log(eta34) + 2*li234    ! eq 2.67 of paper, G_e^2

    mychrg = 1
    if(present(chrg)) mychrg = chrg

    if(present(icoll)) then
      etai3 = proc%Lim_etaij(icoll,3)
      etai4 = proc%Lim_etaij(icoll,4)
      b = b - 2*log(z) * log(E3/E4 * etai3/etai4) * mychrg
      a = a - 4*zeta2
    endif

    calG_QqQl = a * Qcharges*Qcharges + b * Qcharges*Qlept + c*Qlept*Qlept 

  end function calG_QqQl

  function calG_QqQl_L(proc,Qcharges,Qlept)
    use mod_process, only: KinConfig
    implicit none
    type(KinConfig),   intent(in) :: proc
    real(dp),          intent(in) :: Qcharges(:),Qlept
    real(dp) :: calG_QqQl_L(size(Qcharges))
    real(dp) :: EC,E3,E4
    real(dp) :: eta13,eta14,eta23,eta24,eta34,eta34b
    real(dp) :: logEC3,logEC4,logE34,li213,li214,li223,li224,li234
    real(dp) :: b,c

    EC = proc%Lim_Ei(1)
    E3 = proc%Lim_Ei(3)
    E4 = proc%Lim_Ei(4)

    eta13  = proc%Lim_etaij(1,3)
    eta23  = proc%Lim_etaij(2,3)
    eta14  = proc%Lim_etaij(1,4)
    eta24  = proc%Lim_etaij(2,4)
    eta34  = proc%Lim_etaij(3,4)
    eta34b = proc%Lim_KinInv(5)

    logEC3 = log(EC/E3)
    logEC4 = log(EC/E4)
    logE34 = log(E3/E4)

    li234 = real(dilog2(eta34b),kind=dp)
    li213 = real(dilog2(one-eta13),kind=dp)
    li214 = real(dilog2(one-eta14),kind=dp)
    li223 = real(dilog2(one-eta23),kind=dp)
    li224 = real(dilog2(one-eta24),kind=dp)

    b = 2*(li213-li214-li223+li224) &
      + (3._dp + 2*logEC3)*log(eta13/eta23) + (3._dp + 2*logEC4)*log(eta24/eta14)
    c = 13._dp - 4*zeta2 + logE34**2 + (3._dp + 2*logEC3 + 2*logEC4)*log(eta34) + 2*li234

    calG_QqQl_L = -3 * (b * Qcharges*Qlept + c*Qlept*Qlept)

  end function calG_QqQl_L

  !-- Eq. (484) of the notes, hard configuration
  function calG_CF_hard(EC,Ek,etaik,etajk,damp1,damp2)
  real(dp), intent(in) :: EC,Ek,etaik,etajk,damp1,damp2
  real(dp) :: calG_CF_hard,logEcEk,li2jk

  logEcEk = log(EC/Ek)
  li2jk = real(dilog2(one-etajk),kind=dp)

  calG_CF_hard = 13/two - pisq + logEcEk**2 + (3._dp + 2*logEcEk)*log(4*etajk) -&
    (1.5_dp + 2*logEcEk)*(damp1*log(etaik/(one-etaik)) + &
    damp2*log(etajk/(one-etajk))) + 2*li2jk

  end function calG_CF_hard

  !-- Eq. (484) of the notes, collinear limit
  function calG_CF_coll(EC,Ek,etaik,etajk,n)
  real(dp), intent(in) :: EC,Ek,etaik,etajk
  integer, intent(in) :: n !-- collinear direction, n=1,2
  real(dp) :: calG_CF_coll,logEcEk

  logEcEk = log(EC/Ek)
  calG_CF_coll = 0

  if(n == 1) then

    calG_CF_coll = 13/two - pisq + logEcEk**2 + (3._dp + 2*logEcEk)*log(4._dp) -&
      (1.5_dp + 2*logEcEk)*log(etaik)

  else if(n == 2) then

    calG_CF_coll = 13/two - 4*zeta2 + logEcEk**2 + (3._dp + 2*logEcEk)*log(4._dp) +&
      1.5_dp*log(etajk)

  endif

  end function calG_CF_coll

  function calG_Q2(proc,Qcharges,Qlept,icoll,jother,coll_5i)
    use mod_process, only: KinConfig
    implicit none
    type(KinConfig), intent(in) :: proc
    real(dp),        intent(in) :: Qcharges(:),Qlept
    integer,         intent(in) :: icoll,jother
    logical,         intent(in) :: coll_5i
    real(dp) :: calG_Q2(size(Qcharges))
    real(dp) :: EC,E3,E4,E5
    real(dp) :: etaj3,etaj4,eta5i,eta5j,eta53,eta54,eta34,eta34b
    real(dp) :: logEC3,logEC4,logE34,logEC5,li234,li2j3,li2j4,li25j,li253,li254
    real(dp) :: a,b,c

    EC = proc%Lim_Ei(1)
    E5 = proc%Lim_Ei(2)
    E3 = proc%Lim_Ei(3)
    E4 = proc%Lim_Ei(4)

    eta5i  = proc%Lim_etaij(icoll,5)
    eta5j  = proc%Lim_etaij(jother,5)
    etaj3  = proc%Lim_etaij(jother,3)
    etaj4  = proc%Lim_etaij(jother,4)
    eta34  = proc%Lim_etaij(3,4)
    eta53  = proc%Lim_etaij(3,5)
    eta54  = proc%Lim_etaij(4,5)
    eta34b = proc%Lim_KinInv(5)

    logEC3 = log(EC/E3)
    logEC4 = log(EC/E4)
    logEC5 = log(EC/E5)
    logE34 = log(E3/E4)

    li2j3 = real(dilog2(one-etaj3),kind=dp)
    li2j4 = real(dilog2(one-etaj4),kind=dp)
    li25j = real(dilog2(one-eta5j),kind=dp)
    li253 = real(dilog2(one-eta53),kind=dp)
    li254 = real(dilog2(one-eta54),kind=dp)
    li234 = real(dilog2(eta34b),kind=dp)

    !-- collinear limit 5||i
    if(coll_5i) then
      a = 13/two - pisq + logEC5**2 - (1.5_dp + 2*logEC5)*log(eta5i) &
        + (3._dp + 2*logEC5)*log(4*eta5j) + 2*li25j
    !-- hard configuration
    else
      a = 13/two - pisq + logEC5**2 - (1.5_dp + 2*logEC5)*log(eta5i/(one-eta5i)) &
        + (3._dp + 2*logEC5)*log(4*eta5j) + 2*li25j
    endif

    b = 2*(li253-li254-li2j3+li2j4) + 3._dp*(log(eta53/eta54)-log(etaj3/etaj4)) &
      - 2*logEC3*log(etaj3/eta53) + 2*logEC4*log(etaj4/eta54) &
      + 2*logEC5*(logE34 + log(eta53/eta54))
    c = 13._dp - 4*zeta2 + logE34**2 + (3._dp + 2*logEC3 + 2*logEC4)*log(eta34) + 2*li234

    calG_Q2 = a * Qcharges*Qcharges + b * Qcharges*Qlept + c*Qlept*Qlept

  end function calG_Q2


  function calG_ONLOQCD_ns(proc,ampl,lmu,Qlept,ilept)
    use mod_process, only: KinConfig
    ! lmu  = log(mu^2/4/EC^2)
    implicit none
    type(KinConfig),   intent(in) :: proc
    real(dp),          intent(in) :: Qlept,ampl(-5:7,-5:7),lmu(:)
    integer, intent(in)           :: ilept    ! this is the label of the charged lepton
!    integer, optional, intent(in) :: icoll
    real(dp) :: calG_ONLOQCD_ns(-5:7,-5:7,size(lmu))
    real(dp) :: EC,El
    real(dp) :: eta1l,eta2l,lECoEl,li2_ometa1l, li2_ometa2l
    real(dp) :: leta1l,leta2l,li2_1l,li2_2l
    real(dp) :: coeff_Qlsq(size(lmu)), coeff_QqQl(size(lmu)), coeff_QqbQl(size(lmu)),coeff_QqQqb(size(lmu))
    integer  :: al,be

    calG_ONLOQCD_ns = zero
    
    EC = proc%Lim_Ei(1)
    El = proc%Lim_Ei(ilept)   ! for CC, specify the final state lepton

    eta1l  = proc%Lim_etaij(1,ilept)
    eta2l  = proc%Lim_etaij(2,ilept)


    lECoEl = log(EC/El)

    leta1l = log(eta1l)
    leta2l = log(eta2l)
    li2_ometa1l = real(dilog2(one-eta1l),kind=dp)
    li2_ometa2l = real(dilog2(one-eta2l),kind=dp)

    ! mu-independent part
    coeff_QqQqb = -6.0_dp*LECoEl - two*pisq/three
    coeff_Qlsq = 13.0_dp/two - two*LECoEl**2 - two*pisq/three
    coeff_QqQl  = -three*LECoEl*(one + LECoEl) - leta1l*(three + two*LECoEl) - two*li2_ometa1l + pisq/three
    coeff_QqbQl = -three*LECoEl*(one + LECoEl) - leta2l*(three + two*LECoEl) - two*li2_ometa2l + pisq/three

    ! mu-dependent part
    coeff_QqQqb = coeff_QqQqb + three*lmu
    coeff_QqQl  = coeff_QqQl  + lmu*(three + two*lECoEl)
    coeff_QqbQl = coeff_QqbQl + lmu*(three + two*lECoEl)
    coeff_Qlsq  = coeff_Qlsq  + lmu*(three/two + two*lECoEl)


    do al = -5,5
       do be = -5,5
          calG_ONLOQCD_ns(al,be, :) = ampl(al,be)* ( Qlept**2 * coeff_Qlsq(:) + &
               Qlept*Q_IS(al)*coeff_QqQl(:) + Qlept*Q_IS(be)*coeff_QqbQl(:) + Q_IS(al)*Q_IS(be)*coeff_QqQqb(:) )
       enddo
    enddo


  end function calG_ONLOQCD_ns


  function calG_ONLOQCD_gq(proc,ampl,mu,Qlept,ilept,ig,iq)
    use mod_process, only: KinConfig
    ! set ig,iq=1,2 for gq channel; ig,iq=2,1 for qg channel
    implicit none
    type(KinConfig),   intent(in) :: proc
    real(dp),          intent(in) :: Qlept,ampl(-5:7,-5:7),mu(:)
    integer, intent(in)           :: ilept    ! this is the label of the charged lepton
    integer, intent(in)           :: ig,iq    ! this is the label of the IS gluon and IS quark
    !    integer, optional, intent(in) :: icoll
    real(dp) :: calG_ONLOQCD_gq(-5:7,-5:7,size(mu))
    real(dp)    :: Emax,Eq,Eg,El,E5,LE5oEmax,LE5oEl,LEloEq,LEmaxsqoEqEl,LEmaxsqoElE5,LEmaxsqoEqE5,Lmusqo4E2sq(size(mu))
    real(dp)    :: etal5,etaql,etaq5,etag5,li2ometal5,li2ometaql,li2ometaq5,Qqp
    real(dp)    :: coeff_Qqsq(size(mu)),coeff_Qqpsq(size(mu)),coeff_QqQqp(size(mu))
    integer  :: al,be

    calG_ONLOQCD_gq = zero
    
    Eg = proc%Lim_Ei(ig)       ! energy of IS gluon
    Eq = proc%Lim_Ei(iq)       ! energy of IS quark
    El = proc%Lim_Ei(ilept)   ! for CC, specify the final state lepton
    E5 = proc%Lim_Ei(5)       ! energy of FS quark
    Emax = Eq

    LE5oEmax = log(E5/Emax)
    LE5oEl = log(E5/El)
    LEloEq = log(El/Eq)
    LEmaxsqoEqEl = log(Emax**2/Eq/El)
    LEmaxsqoElE5 = log(Emax**2/El/E5)
    LEmaxsqoEqE5 = log(Emax**2/Eq/E5)
    Lmusqo4E2sq = log(mu**2/four/Eq**2)

    etal5  = proc%Lim_etaij(ilept,5)
    etaql  = proc%Lim_etaij(iq,ilept)
    etaq5  = proc%Lim_etaij(iq,5)
    etag5  = proc%Lim_etaij(ig,5)

    li2ometal5 = real(dilog2(one-etal5),kind=dp)
    li2ometaql = real(dilog2(one-etaql),kind=dp)
    li2ometaq5 = real(dilog2(one-etaq5),kind=dp)


    ! mu-independent part
    coeff_Qqsq =  13.0_dp/two + LEloEq**2 + two*li2ometaql - pisq + (three + two*LEmaxsqoEqEl)*Log(etaql)
    coeff_Qqpsq = 13.0_dp + LE5oEl**2 + two*li2ometal5 - two*pisq/three + 3*ln2 - 4*LE5oEmax*ln2 + (three + two*LEmaxsqoElE5)*Log(etal5) + (three/two - two*LE5oEmax)*(Log(one - etag5) - Log(etag5))
    coeff_QqQqp = -13.0_dp + two*LE5oEl*LEloEq - two*li2ometal5 + two*li2ometaq5 - two*li2ometaql + two*pisq/three - &
         (three + two*LEmaxsqoElE5)*Log(etal5) - (three + two*LEmaxsqoEqEl)*Log(etaql) + (three + two*LEmaxsqoEqE5)*Log(etaq5)

    ! mu-dependent part
    coeff_Qqsq = coeff_Qqsq - three*Lmusqo4E2sq/two


    do al = -5,5
       Qqp = Q_IS(al) - Qlept
       calG_ONLOQCD_gq(0,al,:) = ampl(0,al)* ( Q_IS(al)**2 * coeff_Qqsq(:)  + Qqp*Q_IS(al) * coeff_QqQqp(:) + Qqp**2 * coeff_Qqpsq(:))
    enddo


  end function calG_ONLOQCD_gq

end module mod_int_sub_nnlo

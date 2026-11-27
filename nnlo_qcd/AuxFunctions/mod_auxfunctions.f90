module mod_auxfunctions
  use mod_types; use mod_consts_dp; use common_def
  implicit none
  private

  public :: get_NaN, NF1, boost, boostz, scr, spinoru
  public :: aa, bb
  public :: setscale
  public :: lnrat, L0, L1, Lsm1, dilog2, trilog
  public :: F2me4F,Lsm1_2me,klog

contains

  !--- spinor products
  subroutine spinoru(j,p,zza,zzb,ss1) 
    implicit none 
    integer, intent(in) :: j
    complex(dp), intent(out)  :: zza(j,j),zzb(j,j)
    real(dp), intent(out) :: ss1(j,j)
    real(dp), intent(in) :: p(4,j)
    integer :: i1,i2

    zza = czero
    zzb = czero
    ss1 = zero

    do i1=1,j
       do i2=i1+1,j

          zza(i1,i2) = aa(p(:,i1),p(:,i2))
          zzb(i1,i2) = bb(p(:,i1),p(:,i2)) 
          zza(i2,i1) = -zza(i1,i2)
          zzb(i2,i1) = -zzb(i1,i2)

          ss1(i1,i2) = zza(i1,i2)*zzb(i2,i1)
          ss1(i2,i1) = ss1(i1,i2)

       enddo
    enddo

    return 
  end subroutine spinoru


  ! -- spinor helicity routines 

  ! | p > spinor

  function spa(p) 
    real(dp), intent(in) :: p(4)
    complex(dp) :: spa(4)

    spa = u0(p,1)

  end function spa


  ! | p ] spinor

  function spb(p) 
    real(dp), intent(in) :: p(4)
    complex(dp) :: spb(4)

    spb = u0(p,-1)

  end function spb


  ! < p | spinor

  function bspa(p) 
    real(dp), intent(in) :: p(4)
    complex(dp) :: bspa(4)

    bspa = ubar0(p,-1)

  end function bspa


  ! [ p | spinor

  function bspb(p) 
    real(dp), intent(in) :: p(4)
    complex(dp) :: bspb(4)

    bspb = ubar0(p,1)

  end function bspb



  !--- spinor products 


  ! <1 2 >

  function aa(p2,p1)
    real(dp), intent(in) :: p1(4), p2(4) 
    complex(dp) :: aa

    aa = sum(bspa(p2)*spa(p1))

  end function aa


  ! [1 2] 

  function bb(p2,p1)
    real(dp), intent(in) :: p1(4), p2(4) 
    complex(dp) :: bb

    bb = sum(bspb(p2)*spb(p1))

  end function bb



  !--- Spinors

  ! -- ubar spinor, massless
  function ubar0(p,i)
    real(dp), intent(in) :: p(4)
    integer, intent(in) :: i
    complex(dp) :: ubar0(4)
    real(dp)    :: p0,px,py,pz, theta, phi, rrr
    complex(dp) :: cp0
    logical :: flag_nan

    flag_nan = .false.

    p0=p(1)
    px=p(2)
    py=p(3)
    pz=p(4)


    cp0 = cmplx(p0,kind=dp)

    if (p0.eq.zero) then
       write(6,*) 'error in ubar -> p0=0'
       !       pause
    elseif(px.eq.zero.and.py.eq.zero) then
       if ((pz/p0).gt.0.0_dp) theta = zero
       if ((pz/p0).lt.0.0_dp) theta = pi
       phi = zero
    elseif(px.eq.zero.and.py.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = py/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = asin(rrr)
    elseif(py.eq.zero.and.px.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    else
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    endif


    call get_NaN(theta,flag_nan)
    if (flag_nan) then
       write(6,*) 'ubar-th', p0,px,py,pz
       stop
    endif


    call get_NaN(phi,flag_nan)
    if (flag_nan) then
       write(6,*) 'ubar-phi', p0,px,py,pz
       stop
    endif


    if (i.eq.1) then
       ubar0(1)=czero
       ubar0(2)=czero
       ubar0(3)=sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
       ubar0(4)=sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),-sin(phi),kind=dp)
    elseif(i.eq.-1) then
       ubar0(1)= sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),sin(phi),kind=dp)
       ubar0(2)=-sqrt2*sqrt(cp0)*cmplx(abs(cos(theta/two)),0.0_dp,kind=dp)
       ubar0(3)=czero
       ubar0(4)=czero
    else
       stop 'ubar0: i out of range'
    endif

  end function ubar0



  ! -- u0  spinor, massless
  function u0(p,i)
    real(dp), intent(in) :: p(4)
    integer, intent(in) :: i
    complex(dp) :: u0(4)
    real(dp)    :: p0,px,py,pz, theta, phi, rrr
    complex(dp) :: cp0
    logical :: flag_nan

    flag_nan = .false.

    p0=p(1)
    px=p(2)
    py=p(3)
    pz=p(4)

    cp0 = cmplx(p0,kind=dp)

    if (p0.eq.zero) then
       write(6,*) 'error in v0 -> p0=0'
       !       pause
    elseif(px.eq.zero.and.py.eq.zero) then
       if ((pz/p0).gt.0.0_dp) theta = zero
       if ((pz/p0).lt.0.0_dp) theta = pi
       phi = zero
    elseif(px.eq.zero.and.py.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       if (py/p0.gt.zero)  phi = pi/two
       if (py/p0.lt.zero)  phi = three*pi/two
    elseif(py.eq.zero.and.px.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       if (px/p0.gt.zero)  phi = zero
       if (px/p0.lt.zero)  phi = pi
    else
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    endif

    call get_NaN(theta,flag_nan)
    if (flag_nan) then
       write(6,*) 'th-v0', p0,px,py,pz
       write(6,*) 'ratio', pz/p0, acos(pz/p0)
       stop
    endif


    call get_NaN(phi,flag_nan)
    if (flag_nan) then
       write(6,*) 'ph-v0', p0,px,py,pz, px/p0/sin(theta)
       stop
    endif

    if (i.eq.-1) then
       u0(1)=czero
       u0(2)=czero
       u0(3)=sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),-sin(phi),kind=dp)
       u0(4)=-sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
    elseif (i.eq.1) then
       u0(1)= sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
       u0(2)= sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),sin(phi),kind=dp)
       u0(3)=czero
       u0(4)=czero
    else
       stop 'u0: i out of range'
    endif

  end function u0


  ! -- v0  spinor, massless
  function v0(p,i)
    real(dp), intent(in) :: p(4)
    integer, intent(in) :: i
    complex(dp) :: v0(4)
    real(dp)    :: p0,px,py,pz, theta, phi, rrr
    complex(dp) :: cp0
    logical :: flag_nan

    flag_nan = .false.

    p0=real(p(1),kind=dp)
    px=real(p(2),kind=dp)
    py=real(p(3),kind=dp)
    pz=real(p(4),kind=dp)

    cp0 = cmplx(p0,kind=dp)

    if (p0.eq.zero) then
       write(6,*) 'error in v0 -> p0=0'
       !       pause
    elseif(px.eq.zero.and.py.eq.zero) then
       if (pz/p0.gt.0.0_dp) theta = zero
       if (pz/p0.lt.0.0_dp) theta = pi
       phi = zero
    elseif(px.eq.zero.and.py.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       if (py/p0.gt.zero)  phi = pi/two
       if (py/p0.lt.zero)  phi = three*pi/two
    elseif(py.eq.zero.and.px.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       if (px/p0.gt.zero)  phi = zero
       if (px/p0.lt.zero)  phi = pi
    else
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    endif

    call get_NaN(theta,flag_nan)
    if (flag_nan) then
       write(6,*) 'th-v0', p0,px,py,pz
       write(6,*) 'ratio', pz/p0, acos(pz/p0)
       stop
    endif

    call get_NaN(phi,flag_nan)
    if (flag_nan) then
       write(6,*) 'ph-v0', p0,px,py,pz, px/p0/sin(theta)
       stop
    endif

    if (i.eq.1) then
       v0(1)=czero
       v0(2)=czero
       v0(3)=sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),-sin(phi),kind=dp)
       v0(4)=-sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
    elseif (i.eq.-1) then
       v0(1)= sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
       v0(2)= sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),sin(phi),kind=dp)
       v0(3)=czero
       v0(4)=czero
    else
       stop 'v0: i out of range'
    endif

  end function v0



  function scr(p1,p2)   !scalar product of real vectors
    real(dp), intent(in) :: p1(4), p2(4)
    real(dp) :: scr
    scr = p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)-p1(4)*p2(4)
  end function scr



  ! take momenta p_in in frame in which particle one is at rest with mass
  ! "mass" and convert to frame in which particle one has fourvector p1
  subroutine boost(mass,p1,p_in,p_out)
    use mod_types; use mod_consts_dp; use common_def
    implicit none
    real(dp), intent(in)  ::  mass,p1(4),p_in(4)
    real(dp), intent(out) ::  p_out(4)
    ! ------------------------------------
    real(dp) ::  gam,beta(2:4),bdotp
    integer :: j,k

    gam=p1(1)/mass
    bdotp=zero
    do j=2,4
       beta(j)=-p1(j)/p1(1)
       bdotp=bdotp+p_in(j)*beta(j)
    enddo
    p_out(1)=gam*(p_in(1)-bdotp)
    do k=2,4
       p_out(k)=p_in(k)+gam*beta(k)*(gam/(gam+one)*bdotp-p_in(1))
    enddo
  end subroutine boost


  !-- boost along the z-axis from the com frame to the lab frame according to
  !-- p1 = sqrt(s*xa*xb)/2 * (1,0,0,+1) --> sqrt(s)/2 * xa * (1,0,0,+1)
  !-- p2 = sqrt(s*xa*xb)/2 * (1,0,0,-1) --> sqrt(s)/2 * xb * (1,0,0,-1)
  subroutine boostz(xa,xb,npart,pin,pout)
    real(dp), intent(in) :: xa,xb
    integer, intent(in) :: npart
    real(dp), intent(in) :: pin(4,npart)
    real(dp), intent(out) :: pout(4,npart)
    real(dp) :: beta,gamma
    integer :: i

    beta = (xb-xa)/(xb+xa)
    gamma = (xa+xb)/two/sqrt(xa*xb)

    pout = pin

    do i = 1,npart
       pout(1,i) = gamma * (pin(1,i)-beta*pin(4,i))
       pout(4,i) = gamma * (pin(4,i)-beta*pin(1,i))
    enddo

    return

  end subroutine boostz


  function NF1(x3,x4,x5) 
    implicit none  
    real(dp), intent(in) :: x3,x4,x5
    real(dp) :: NF1

    if (x4.eq.zero) then 
     NF1 = one 
     return 
    elseif(x4.eq.one) then 
     NF1 = four*x5*(one-x3) 
     return 
    else
    NF1 = 1.0_dp + x4*(1.0_dp-2.0_dp*x3) & 
         - 2.0_dp*(1.0_dp-2.0_dp*x5)*sqrt(abs(x4*(1.0_dp-x3)*(1.0_dp-x3*x4)))
    endif 

  end function NF1


  subroutine get_NaN(value,flag_nan)
    implicit none
    real(dp), intent(in)  :: value
    logical, intent(out) :: flag_nan

    flag_nan = .false.

    if (.not.value.le.0.0_dp .and. .not.value.gt.0.0_dp ) then
       flag_nan = .true.
    endif

  end subroutine get_NaN

  !-- useful integrals

  function klog(x)
    real(dp), intent(in)  :: x
    complex(dp) :: klog
    
    if (x.gt.zero) then 
       klog = log(x) 
    else
       klog = log(abs(x)) - ci*pi
    endif
    
  end function klog

  function Lnrat(x,y) 
    real(dp), intent(in) :: x,y
    complex(dp) :: Lnrat
    
    Lnrat = klog(x) - klog(y)
    
  end function Lnrat

  function L0(r1,r2) 
    implicit none
    real(dp), intent(in) :: r1,r2
    complex(dp) :: L0
    
    L0 = Lnrat(r1,r2)/(one-r1/r2)
    
  end function L0
  
  function L1(r1,r2) 
    implicit none
    real(dp), intent(in) :: r1,r2
    complex(dp) :: L1
    
    L1 = (L0(r1,r2)+one)/(one-r1/r2)
    
  end function L1

  function Lsm1(r1,r2,r3,r4) 
    implicit none
    real(dp), intent(in) :: r1,r2,r3,r4
    complex(dp) :: Lsm1
    
    Lsm1 = kdilog(r1,r2)+ kdilog(r3,r4) & 
         +Lnrat(r1,r2)*Lnrat(r3,r4)-pisq/6.0_dp
    
  end function Lsm1

  function kdilog(r1,r2)  ! this is really just dilog(1-r1/r2) but with 
    implicit none           ! particular imaginary part 
    real(dp), intent(in) :: r1, r2
    real(dp) :: x 
    complex(dp) :: kdilog 
    
    x = r1/r2
    
    if (x.gt.zero) then 
       kdilog = dilog2(one-x)
    else
       kdilog = pisq/6.0_dp - dilog2(x) & 
            - log(one-x)*Lnrat(r1,r2) 
    endif
    
  end function kdilog

  FUNCTION dilog2(xIn)
    IMPLICIT NONE
    real(dp), intent(in):: xIn
    real(dp):: x,z,z2
    complex(dp):: dilog2, Li2tmp,Const
    real(dp):: Fact
    real(dp), parameter:: Pi26=1.64493406684822643647241516665_dp
    real(dp), parameter:: Pi23=3.28986813369645287294483033329_dp
    real(dp), parameter:: B1= -0.25_dp
    real(dp), parameter:: B2=  2.7777777777777778d-02
    real(dp), parameter:: B4= -2.7777777777777778d-04
    real(dp), parameter:: B6=  4.7241118669690098d-06
    real(dp), parameter:: B8= -9.1857730746619635d-08
    real(dp), parameter:: B10= 1.8978869988970999d-09
    real(dp), parameter:: B12=-4.0647616451442255d-11
    real(dp), parameter:: B14= 8.9216910204564525d-13
    real(dp), parameter:: B16=-1.9939295860721075d-14
    real(dp), parameter:: B18= 4.5189800296199181d-16
    
    Const=cmplx(0.0_dp,0.0_dp)
    Fact =1.0_dp
    
    x = xIn
    if ( x.gt.1.0_dp.and.x.lt.2.0_dp  ) then
       Fact = -1.0_dp
       Const = Pi23-0.5_dp*log(x)**2-pi*log(x)*ci
       x = 1.0_dp/x
       write(6,*) 'error in dilog2'
       stop
    elseif(x.gt.2.0_dp) then
       Fact = -1.0_dp
       Const = Pi23-0.5_dp*log(x)**2-pi*log(x)*ci
       x = 1.0_dp/x
    elseif ( x.lt.(-1.0_dp) ) then
       Fact = -1.0_dp
       Const =-Pi26-0.5_dp*log(-x)**2
       x = 1.0_dp/x
    elseif ( x.eq.1.0_dp) then
       dilog2 = Pi26
       return
    endif
    
    if ( x.gt.0.5_dp ) then
       Fact = -1.0_dp*Fact
       Const = Const + Pi26 - log(x)*log(1.0_dp-x)
       x = 1.0_dp-x
    endif
    
    z = -log(1.0_dp-x)
    Li2tmp = z
    z2 = z*z
    Li2tmp = Li2tmp + B1  * z2
    z = z*z2
    Li2tmp = Li2tmp + B2  * z
    z = z*z2
    Li2tmp = Li2tmp + B4  * z
    z = z*z2
    Li2tmp = Li2tmp + B6  * z
    z = z*z2
    Li2tmp = Li2tmp + B8  * z
    z = z*z2
    Li2tmp = Li2tmp + B10 * z
    z = z*z2
    Li2tmp = Li2tmp + B12 * z
    z = z*z2
    Li2tmp = Li2tmp + B14 * z
    z = z*z2
    Li2tmp = Li2tmp + B16 * z
    z = z*z2
    Li2tmp = Li2tmp + B18 * z
    
    dilog2 = Fact*Li2tmp + Const
    
    return
  END FUNCTION dilog2

  function trilog(z)
    real(dp) :: z, trilog
    complex(8) :: xli3
    
    trilog = real(xli3(dcmplx(z,0.0d0)),kind=dp)
    
  end function trilog

  function F2me4F(s,t,Psq,Qsq)
    complex(dp)    :: F2me4F
    real(dp), intent(in)   :: s,t,Psq,Qsq

    F2me4F = -two*( kdilog(Psq,s)+kdilog(Psq,t)+kdilog(Qsq,s)+kdilog(Qsq,t) &
         -kdilog(Psq*Qsq,s*t) + half*(lnrat(s,t))**2 )

  end function F2me4F


  function Lsm1_2me(s,t,m1sq,m3sq)
    implicit none
    !---- Taken from MCFM v8       
    !---- formula taken from G.~Duplancic and B~Nizic [arXiv:hep-ph/0006249 v2]
    !---- Eq 71
    !---- Lsm1_2me notation follows from 
    !----  Z.~Bern, L.~J.~Dixon and D.~A.~Kosower,
    !----  %``Dimensionally regulated pentagon integrals,''
    !----  Nucl.\ Phys.\ B {\bf 412}, 751 (1994)
    !----  [arXiv:hep-ph/9306240].
    !----  %%CITATION = HEP-PH 9306240;%%
    !----  Eqs. (I.13)
    !---- analytic continuation has been checked by calculating numerically.
    
    complex(dp):: Lsm1_2me
    real(dp), intent(in)  :: s,t,m1sq,m3sq
    integer:: j
    real(dp):: arg(4),omarg(4),f2me,htheta
    complex(dp):: Li2(4),wlog(4)
    !--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
    htheta(s)=half+half*sign(one,s)
    
    f2me=(s+t-m1sq-m3sq)/(s*t-m1sq*m3sq)
    
    arg(1)=f2me*s
    arg(2)=f2me*t
    arg(3)=f2me*m1sq
    arg(4)=f2me*m3sq
    
    do j=1,4
       omarg(j)=one-arg(j)
       wlog(j)=log(abs(arg(j)))    &
            +ci*pi*cone*(htheta(-arg(j))*sign(one,f2me))
       if (omarg(j) > one) then 
          Li2(j)=cone*(pisq/6.0_dp-dilog2(arg(j))) &
               -wlog(j)*cone*(log(omarg(j)))
       else
          Li2(j)=cone*(dilog2(omarg(j)))
       endif
    enddo
    Lsm1_2me=Li2(1)+Li2(2)-Li2(3)-Li2(4)
    
  end function Lsm1_2me

  subroutine setscale(mu,muF,muR)
    ! takes calculated scale mu and multiplies by input scale variation factor to get mupdf and musq
    real(dp), intent(in)   :: mu
    real(dp), intent(out)  :: muF,muR

    muF = mu * xMuF
    muR = mu * xMuR

  end subroutine setscale

end module mod_auxfunctions

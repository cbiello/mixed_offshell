module mod_auxfunctions
  use mod_types
  use mod_consts_dp
  implicit none
  private

  public :: is_nan
  public :: mypause
  public :: boost,boostz
  public :: spinoru,aa,bb
  public :: scr
  public :: dilog2
  public :: trilog

  interface is_nan
     module procedure is_nan_single
     module procedure is_nan_vec
  end interface is_nan
  
contains

  subroutine boost(mass,p1,p_in,p_out)
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
    integer,  intent(in) :: npart
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
  
  !-- intrinsic isnan does not catch \pm infinity
  !-- here use that NaN is unordered, i.e. NaN \= NaN, NaN is not greater or smaller than NaN
  !-- not tested on ifort
  !-- 
  function is_nan_single(myvalue) result(is_nan)
    real(dp), intent(in)  :: myvalue
    logical :: is_nan
    
    is_nan = .false.

    if (.not. abs(myvalue) <= huge(myvalue) ) is_nan = .true.

  end function is_nan_single

  !-- perhaps merge these two in an interface
  function is_nan_vec(myvec) result(is_nan)
    real(dp), intent(in)  :: myvec(:)
    logical :: is_nan
    integer :: i
    
    is_nan = .false.

    do i = 1,size(myvec)
       if (.not. abs(myvec(i)) <= huge(myvec(i)) ) then
          is_nan = .true.
          return
       endif
    enddo

  end function is_nan_vec

  subroutine mypause()
    write(*,*) 'Press a key to continue'
    read(*,*)
  end subroutine mypause

  function scr(p1,p2)
    real(dp) :: scr,p1(4),p2(4)

    scr = p1(1)*p2(1)-dot_product(p1(2:4),p2(2:4))

  end function scr

  !-- Li2
  function dilog2(xIn)
    implicit none
    real(dp), intent(in):: xIn
    real(dp):: x,z,z2
    complex(dp):: dilog2, Li2tmp,Const
    real(dp):: Fact
    real(dp), parameter:: Pi26=1.64493406684822643647241516665_dp
    real(dp), parameter:: Pi23=3.28986813369645287294483033329_dp
    real(dp), parameter:: B1= -0.25_dp
    real(dp), parameter:: B2=  2.7777777777777778E-02_dp
    real(dp), parameter:: B4= -2.7777777777777778E-04_dp
    real(dp), parameter:: B6=  4.7241118669690098E-06_dp
    real(dp), parameter:: B8= -9.1857730746619635E-08_dp
    real(dp), parameter:: B10= 1.8978869988970999E-09_dp
    real(dp), parameter:: B12=-4.0647616451442255E-11_dp
    real(dp), parameter:: B14= 8.9216910204564525E-13_dp
    real(dp), parameter:: B16=-1.9939295860721075E-14_dp
    real(dp), parameter:: B18= 4.5189800296199181E-16_dp
    
    Const=cmplx(0.0_dp,0.0_dp,kind=dp)
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
  end function dilog2

  ! Li3
  function trilog(z)
    real(dp) :: z, trilog
    complex(8) :: xli3

    trilog = real(xli3(dcmplx(z,0.0d0)),kind=dp)

  end function trilog


  
  !-- spinor routines

  subroutine spinoru(j,p,zza,zzb,ss1) 
    implicit none 
    integer, intent(in) :: j
    complex(dp), intent(out) :: zza(:,:),zzb(:,:)
    real(dp), intent(out) :: ss1(:,:)
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

          ss1(i1,i2) = real(zza(i1,i2)*zzb(i2,i1),kind=dp)
          ss1(i2,i1) = ss1(i1,i2)

       enddo
    enddo

    return 

  end subroutine spinoru

  !-- <1 2>
  function aa(p2,p1)
    real(dp), intent(in) :: p1(:), p2(:) 
    complex(dp) :: aa
    
    aa = sum(bspa(p2)*spa(p1))
    
  end function aa

  !-- [1 2] 
  function bb(p2,p1)
    real(dp), intent(in) :: p1(:), p2(:) 
    complex(dp) :: bb

    bb = sum(bspb(p2)*spb(p1))

  end function bb

  !--

  !-- | p >
  function spa(p) 
    real(dp), intent(in) :: p(:)
    complex(dp) :: spa(4)

    spa = u0(p,1)

  end function spa

  !-- | p ]
  function spb(p) 
    real(dp), intent(in) :: p(:)
    complex(dp) :: spb(4)

    spb = u0(p,-1)

  end function spb

  !-- < p | 
  function bspa(p) 
    real(dp), intent(in) :: p(:)
    complex(dp) :: bspa(4)

    bspa = ubar0(p,-1)

  end function bspa

  !--[ p | spinor
  function bspb(p) 
    real(dp), intent(in) :: p(:)
    complex(dp) :: bspb(4)

    bspb = ubar0(p,1)

  end function bspb

  !--

  ! -- u0  spinor, massless
  function u0(p,i)
    real(dp), intent(in) :: p(:)
    integer, intent(in) :: i
    complex(dp) :: u0(4)
    real(dp)    :: p0,px,py,pz, theta, phi, rrr
    complex(dp) :: cp0

    p0=p(1)
    px=p(2)
    py=p(3)
    pz=p(4)

    cp0 = cmplx(p0,kind=dp)

    if (p0.eq.zero) then
       write(6,*) 'error in v0 -> p0=0'
       theta = zero
       phi   = zero
       stop
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

    if (is_NaN(theta)) then
       write(6,*) 'th-v0', p0,px,py,pz
       write(6,*) 'ratio', pz/p0, acos(pz/p0)
       stop
    endif

    if (is_NaN(phi)) then
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

  function ubar0(p,i)
    real(dp), intent(in) :: p(:)
    integer, intent(in) :: i
    complex(dp) :: ubar0(4)
    real(dp)    :: p0,px,py,pz, theta, phi, rrr
    complex(dp) :: cp0

    p0=p(1)
    px=p(2)
    py=p(3)
    pz=p(4)

    cp0 = cmplx(p0,kind=dp)

    if (p0.eq.zero) then
       write(6,*) 'error in ubar -> p0=0'
       theta = zero
       phi   = zero
       stop
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

    if (is_NaN(theta)) then
       write(6,*) 'ubar-th', p0,px,py,pz
       stop
    endif

    if (is_NaN(phi)) then
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

end module mod_auxfunctions

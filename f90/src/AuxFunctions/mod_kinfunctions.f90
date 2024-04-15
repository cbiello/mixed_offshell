module mod_kinfunctions
  use mod_types
  use mod_consts_dp
  use mod_auxfunctions
  implicit none
  real(dp), parameter :: ybuff = 1E-8_dp    !-- buffer for rapidity
  real(dp), parameter :: phibuff = 1E-14_dp !-- buffer for azimuth
  real(dp), parameter :: yinf  = 99._dp
  private

  public :: get_pt,get_y,get_eta,get_dy,get_r,get_dphi,get_costh_star
  public :: ktjet,pt_order

contains

  function get_pt(p)
    real(dp) :: get_pt,p(4)

    get_pt = sqrt(p(2)**2+p(3)**2)

    return

  end function get_pt

  function get_y(p)
    real(dp) :: get_y,p(4)
    real(dp) :: pp,pm

    pp = p(1)+p(4)
    pm = p(1)-p(4)

    if (abs(pm).lt.ybuff) then
       get_y =  yinf
    elseif (abs(pp).lt.ybuff) then
       get_y = -yinf
    else
       get_y = half * log(pp/pm)
    endif

    return

  end function get_y
  
  function get_eta(p)
    real(dp) :: get_eta,p(4)
    real(dp) :: pmod,pp,pm

    pmod = sqrt(p(2)**2+p(3)**2+p(4)**2)
    pp = pmod+p(4)
    pm = pmod-p(4)

    if (abs(pm).lt.ybuff) then
       get_eta =  yinf
    elseif (abs(pp).lt.ybuff) then
       get_eta = -yinf
    else
       get_eta = half * log(pp/pm)
    endif

    return

  end function get_eta

  function get_dy(p1,p2)
    real(dp) :: get_dy,p1(4),p2(4)
    real(dp) :: pp1,pp2,pm1,pm2,pnum,pden

    pp1 = p1(1)+p1(4)
    pp2 = p2(1)+p2(4)

    pm1 = p1(1)-p1(4)
    pm2 = p2(1)-p2(4)

    pnum = pp1*pm2
    pden = pm1*pp2

    if (abs(pden).lt.ybuff) then
       get_dy =  yinf
    elseif (abs(pnum).lt.ybuff) then
       get_dy = -yinf
    else
       get_dy = half * log(pnum/pden)
    endif

    return

  end function get_dy
     
  function get_r(p1,p2)
    real(dp) :: get_r,p1(4),p2(4)
    real(dp) :: dy,dphi

    dy = get_dy(p1,p2)
    dphi = get_dphi(p1,p2)

    get_r = sqrt(dy**2+dphi**2)

    return

  end function get_r

  function get_dphi(p1,p2)
    real(dp) :: get_dphi,p1(4),p2(4)
    real(dp) :: pt1,pt2,cosdphi

    pt1 = get_pt(p1)
    pt2 = get_pt(p2)

    cosdphi = (p1(2)*p2(2)+p1(3)*p2(3))/pt1/pt2

    if (cosdphi.gt.one) then
       get_dphi = zero
    elseif (cosdphi.lt.-one) then
       get_dphi = pi
    elseif (abs(cosdphi+one) < phibuff) then
       get_dphi = pi
    else
       get_dphi = acos(cosdphi)
    endif

    return

  end function get_dphi

  function get_costh_star(p1,p2)
    real(dp) :: get_costh_star,p1(4),p2(4)
    real(dp) :: p1p,p1m,p2p,p2m
    real(dp) :: pll(4)
    real(dp) :: mll2,ptll,pzll
    
    p1p = (p1(1) + p1(4))/sqrt2
    p1m = (p1(1) - p1(4))/sqrt2
    p2p = (p2(1) + p2(4))/sqrt2
    p2m = (p2(1) - p2(4))/sqrt2
    
    pll  = p1 + p2
    mll2 = scr(pll,pll)
    ptll = get_pt(pll)
    pzll = pll(4)
    
    get_costh_star = two*(p1p*p2m - p1m*p2p)/sqrt(mll2*(mll2 + ptll**2))
    get_costh_star = get_costh_star*pzll/abs(pzll)

    return
    
  end function get_costh_star

  !--

  subroutine ktjet(jetR, ptjetcut, pin, flin, nin, pjet, fljet, njet)
    integer, intent(in) :: nin
    real(dp), intent(in) :: jetR,ptjetcut
    real(dp), intent(in) :: pin(4,nin)
    integer,intent(in):: flin(1:2,nin)
    integer, intent(out) :: njet
    real(dp), intent(out) :: pjet(4,nin)
    integer,intent(out):: fljet(1:2,nin)
    real(dp) :: pti,dbi,ptj,dbj
    real(dp) :: yi,yj,dy,d_cosphi,dphi,r12sq,dij
    real(dp) :: dmin,rsq
    real(dp) :: prehard(4,nin),prejet(4,nin)
    integer :: prefljet(1:2,nin),flhard(1:2,nin)
    integer :: nh,nj,ijmin(2)
    integer :: i,j,le
    real(dp), parameter :: ktclass = -one   ! antikt

    rsq = jetR**2

    le = 0

    njet = 0
    pjet = zero

    do i=1,nin
       if (pin(1,i).ne.zero.and.(abs(pin(4,i)).ne.abs(pin(1,i)))) then
          le = le+1
          prehard(:,le) = pin(:,i)
          flhard(:,le) = flin(:,i)
       endif
    enddo

    nj = 0
    nh = le
    prejet = zero
    prefljet= 0

    !   recombine

    do while (nh.gt.0) ! recombination loop

       pti = sqrt(prehard(2,nh)**2+prehard(3,nh)**2)
       dbi = pti**(2*ktclass)
       dmin = dbi
       ijmin = (/nh,0/)

       do i = 1, nh-1 ! compute db, dij

          pti = sqrt(prehard(2,i)**2+prehard(3,i)**2)
          dbi = pti**(2*ktclass)

          if (dbi.lt.dmin) then
             dmin = dbi
             ijmin = (/i,0/)
          endif

          do j = i+1,nh
             ptj = sqrt(prehard(2,j)**2+prehard(3,j)**2)
             dbj = ptj**(2*ktclass)

             yi = half * log((prehard(1,i)+prehard(4,i))/(prehard(1,i)-prehard(4,i)))
             yj = half * log((prehard(1,j)+prehard(4,j))/(prehard(1,j)-prehard(4,j)))

             dy = yi-yj

             ! cos(a-b) = cos(a)cos(b)+sin(a)sin(b)
             d_cosphi = (prehard(2,i)*prehard(2,j)+prehard(3,i)*prehard(3,j))/pti/ptj
             if (d_cosphi.gt.one) then
                dphi = zero
             elseif (d_cosphi.lt.-one) then
                dphi = pi
             else
                dphi = acos(d_cosphi) ! 0<dphi<pi
             endif

             r12sq = dy**2 + dphi**2

             dij = min(dbi,dbj) * r12sq/rsq

             if (dij.lt.dmin) then
                dmin = dij
                ijmin = (/min(i,j),max(i,j)/)
             endif

          enddo

       enddo ! compute db, dij

       if (ijmin(2).eq.0) then ! no recombination, 1 new jet candidate
          nj = nj+1
          prejet(:,nj) = prehard(:,ijmin(1))
          prefljet(:,nj) = flhard(:,ijmin(1))

          do i = ijmin(1),nh-1
             prehard(:,i) = prehard(:,i+1)
             flhard(:,i) = flhard(:,i+1)
          enddo

       else ! recombine, no new jet candidates

          prehard(:,ijmin(1)) = prehard(:,ijmin(1))+prehard(:,ijmin(2))
          flhard(:,ijmin(1)) = flhard(:,ijmin(1)) + flhard(:,ijmin(2))

          do i = ijmin(2),nh-1
             prehard(:,i) = prehard(:,i+1)
             flhard(:,i) = flhard(:,i+1)
          enddo

       endif

       nh = nh - 1

    enddo ! recombination loop


    do i = 1, nj
       pti = sqrt(prejet(2,i)**2+prejet(3,i)**2)
       if (pti.gt.ptjetcut) then
          njet = njet + 1
          pjet(:,njet) = prejet(:,i)
          fljet(:,njet) = prefljet(:,i)
       endif

    enddo

    return

  end subroutine ktjet

  !#######################################################################
  !##### this subroutine is for ordering jets according to pT
  !#######################################################################
  subroutine pT_order(N,Mom)
    implicit none
    integer :: N
    !!!TODO: should the following be real(dp) ?
    real(8) :: Mom(1:4,1:N),Mom_Tmp(1:4,1:N),pTList(1:N)
    integer :: i,MomOrder(1:N)


    if(N.lt.1) return
    do i=1,N
       pTList(i) = get_PT(real(Mom(1:4,i),kind=dp))
       MomOrder(i) = i
    enddo

    call BubleSort(N,pTList(1:N),MomOrder(1:N))

    Mom_Tmp(1:4,1:N) = Mom(1:4,1:N)
    do i=1,N
       Mom(1:4,i) = Mom_Tmp(1:4,MomOrder(i))
    enddo

  end subroutine pT_order

  subroutine BubleSort(N,X, IY)
    IMPLICIT NONE
    integer n
    real(8) x(1:n)
    integer iy(1:n)
    real(8) temp
    integer i, j, jmax, itemp

    jmax=n-1
    do i=1,n-1
       temp=1d38
       do j=1,jmax
          if(x(j).gt.x(j+1)) cycle
          temp=x(j)
          x(j)=x(j+1)
          x(j+1)=temp
          itemp=iy(j)
          iy(j)=iy(j+1)
          iy(j+1)=itemp
       enddo
       if(temp.eq.1d38) return
       jmax=jmax-1
    enddo

    return
  end subroutine BubleSort
 
end module mod_kinfunctions

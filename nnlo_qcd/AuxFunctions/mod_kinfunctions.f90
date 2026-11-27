module mod_kinfunctions
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  implicit none
  real(dp), parameter :: ybuff = 1E-8_dp !-- buffer for rapidity
  private

  public :: get_pt,get_y,get_eta,get_dy,get_r,get_dphi

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
       get_y = 99.0_dp
    elseif (abs(pp).lt.ybuff) then
       get_y = -99.0_dp
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
       get_eta = 99.0_dp
    elseif (abs(pp).lt.ybuff) then
       get_eta = -99.0_dp
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
       get_dy = 99.0_dp
    elseif (abs(pnum).lt.ybuff) then
       get_dy = -99.0_dp
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
    else
       get_dphi = acos(cosdphi)
    endif

    return

  end function get_dphi

end module mod_kinfunctions

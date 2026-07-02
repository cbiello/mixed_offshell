!-- In this file all splittings are taken without any colour/charge factors
module mod_splittings_bare
  use mod_types
  use mod_consts_dp
  implicit none
  private

  public :: Pqq,Pqg,Pgq_spav
  public :: Pgaq,Pgqqb_spav_ab,Pqbqq_id_symm
  public :: Pgaq_ds
  public :: Pgaq_s2
  public :: PgqAP,PqgAP,PqqAP
  public :: calPp_qq,calPp_gq,calPp_qg
  public :: Pqq_NLO,PqqAP_0_R

  
contains

  !-- Double collinear
  function Pqq(z)
    real(dp) :: z,Pqq

    Pqq = (one+z**2)/(one-z) ! * Cf

  end function Pqq

  function Pqg(z)
    real(dp) :: z,Pqg

    Pqg = (one+(one-z)**2)/z ! * Cf

  end function Pqg

  function Pgq_spav(z)
    real(dp) :: z,Pgq_spav

    Pgq_spav = (one-z)**2 + z**2 ! * Tr

  end function Pgq_spav
  
  !-- Triple collinear
  function Pgaq(s12,s13,s23,z1,z2,z3)
    real(dp) :: Pgaq
    real(dp), intent(in) :: s12,s13,s23,z1,z2,z3
    real(dp) :: s123

    s123 = s12+s13+s23

    !-- s12 is not needed
    pgaq = Pggq_ab_unsymm(s13,s23,s123,z1,z2,z3) &
         + Pggq_ab_unsymm(s23,s13,s123,z2,z1,z3)

    pgaq = pgaq * four/s123**2

    return

  end function Pgaq

  !! hep-ph/9908523 equation (61)
  !-- s12 is not needed
  function Pggq_ab_unsymm(s13,s23,s123,z1,z2,z3)
    real(dp) :: Pggq_ab_unsymm
    real(dp), intent(in) :: s13,s23,s123,z1,z2,z3

    Pggq_ab_unsymm = s123**2/two/s13/s23*z3*(one+z3**2)/z1/z2 &
         + s123/s13* (z3*(one-z1)+(one-z2)**3)/z1/z2 &
         - s23/s13

  end function Pggq_ab_unsymm

  !! Eq. (64) of hep-ph/9908523 averaged over vector boson polarisation,
  !! i.e. P^{\mu\nu(ab)}_{g_1 q_2 \bar{q}_3} d_{\mu\nu}/(d-2)
  function Pgqqb_spav_ab(s12,s13,s23,z1,z2,z3)
    real(dp) :: Pgqqb_spav_ab
    real(dp), intent(in) :: s12,s13,s23,z1,z2,z3
    real(dp) :: s123, A, B, C

    s123 = s12 + s13 + s23

    A = -2*(one + s123*(one/s12+one/s13)*(one+z1) &
            - (one+z1**2)*s123**2/s12/s13)

    !-- piece proportional to 1/(1-ep)
    B = 2*s123*((z1+z2)/s12 + (z1+z3)/s13 - s123*(z1+2*z2*z3)/s12/s13)

    !-- piece proportional to (1-ep)
    C = two + s12/s13 + s13/s12

    Pgqqb_spav_ab = (A + B + C) * 4/s123**2 ! * Tr * Cf

    return

  end function Pgqqb_spav_ab

  function Pqbqq_id_unsymm(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pqbqq_id_unsymm
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3
    
    Pqbqq_id_unsymm = two*s23/s12 &
         +s123/s12*( (one+z1**2)/(one-z2)-two*z2/(one-z3) ) &
         -s123**2/s12/s13*z1/two*( (one+z1**2)/(one-z2)/(one-z3) ) 
    
    !Pqbqq_id_unsymm = Pqbqq_id_unsymm * Cf * (Cf-half*Ca)

    return

  end function Pqbqq_id_unsymm

  function Pqbqq_id_symm(s12,s13,s23,z1,z2,z3)
    real(dp) :: Pqbqq_id_symm
    real(dp), intent(in) :: s12,s13,s23,z1,z2,z3
    real(dp) :: s123,p1,p2

    s123 = s12 + s13 + s23

    p1 = Pqbqq_id_unsymm(s12,s13,s23,s123,z1,z2,z3)
    p2 = Pqbqq_id_unsymm(s13,s12,s23,s123,z1,z3,z2)

    Pqbqq_id_symm = p1+p2

    Pqbqq_id_symm = Pqbqq_id_symm * four/s123**2
    
    return

  end function Pqbqq_id_symm
  
  !-- Triple-collinear + double-soft
  function Pgaq_ds(s13,s23,z1,z2)
    real(dp), intent(in) :: s13,s23,z1,z2
    real(dp) :: Pgaq_ds
    real(dp) :: z12,s123,t123
    
    z12 = z1+z2
    
    s123 = s13+s23
    t123 = two * (z1*s23-z2*s13)/z12
    
    !-- ``abelian-like'' piece
    Pgaq_ds = 16.0_dp/s13/s23/z1/z2 ! * Cf * Q**2
    
    return
    
  end function Pgaq_ds

  !-- Triple-collinear + single-soft
  !-- z2~s12~s23->0 limit of Pgaq
  !-- z1 = e1/(e1+e3), z2 = e2/(e1+e3)
  function Pgaq_s2(s13,s23,z1,z2)
    real(dp), intent(in) :: s13,s23,z1,z2
    real(dp) :: Pgaq_s2
    
    !-- ``abelian-like'' piece
    Pgaq_s2 = 8._dp*(two-4._dp*z1+3._dp*z1**2-z1**3)/z1/s13
    Pgaq_s2 = Pgaq_s2/z2/s23
    
    return
    
  end function Pgaq_s2

!! notation for Altarelli-Parisi functions: [-1 --> delta, 0 --> reg including plus, +1 --> plus]

  function PqqAP(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res(-1:1)

    res(-1) = half*three

    res( 0) =  (one+z**2)/(one-z)

    res(+1) = two/(one-z)

  end function PqqAP

  function PqgAP(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res(-1:1)

    res = 0

    res(0) =z**2 + (one-z)**2

  end function PqgAP

function PgqAP(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res(-1:1)

    res = 0

    res(0) = one+(one-z)**2/z

  end function PgqAP


  function calPp_qq(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res(-1:1)

    res(-1) = zero

    res( 0) = -(one+z**2)/(one-z)*log((one-z)**2/z) - (one-z)

    res(+1) = -four*log(one-z)/(one-z)

  end function calPp_qq



 function calPp_qg(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res(-1:1)

    res(-1) = zero

    res( 0) = -(z**2 + (one-z)**2)*(log((one-z)**2/z)-one) - one

    res(+1) = zero

  end function CalPp_qg

  function calPp_gq(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: res(-1:1)

    res = 0

    res(0) = -(one+(one-z)**2)/z*log((one-z)**2/z) - z

  end function CalPp_gq

  !-- only scale indepedent term
  function Pqq_NLO(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: omz,lomz,res(-1:1)

    omz = one - z
    lomz = log(omz)

    res(-1) = zero

    res( 0) = omz - 2*(one+z)*lomz !-- reg
    res( 0) = res(0) + 4*lomz/omz  !-- plus

    res(+1) = 4*lomz/omz           !-- plus

  end function Pqq_NLO

  function PqqAP_0_R(z) result(res)
    real(dp), intent(in) :: z
    real(dp) :: omz,res(-1:1)

    omz = one-z

    res(-1) = zero

    res( 0) = one + z        !-- reg
    res( 0) = res(0) - 2/omz !-- plus

    res(+1) = -2/omz         !-- "minus" component of the plus distribution

  end function PqqAP_0_R

end module mod_splittings_bare

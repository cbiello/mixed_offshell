module mod_limits_dy_qg
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions; use mod_procfunctions
  use mod_splitting
  implicit none

  private

  !-- soft
  public :: eik_s5

  !-- soft-collinear
  public :: Pggq_s5

contains

  function Eik_s5(p1,p2,p4,p5)
    real(dp), intent(in) :: p1(4),p2(4),p4(4),p5(4)
    real(dp) :: Eik_s5
    real(dp) :: ES12_5, ES14_5, ES24_5, ES12_4 
    real(dp) :: s14, s15, s24, s25, s12, s45

    !   sij's >are defined without factors of two here 

    s14 = -scr(p1,p4) 
    s15 = -scr(p1,p5) 
    s24 = -scr(p2,p4) 
    s25 = -scr(p2,p5) 
    s45 = scr(p4,p5) 
    s12 = scr(p1,p2) 

    ES12_5 = s12/s15/s25
    ES14_5 = s14/s15/s45
    ES24_5 = s24/s25/s45 

    Eik_s5 = (two*Cf-Ca) * ES14_5 + Ca * ( ES12_5 + ES24_5 ) 

    return 

  end function Eik_s5

  function Pggq_s5(kinlist) 
    real(dp) :: Pggq_s5 
    real(dp), intent(in) :: kinlist(6)
    real(dp) :: scr24,scr25,scr45,z2,z4,z5,z

    scr24 = kinlist(1)*half
    scr25 = kinlist(2)*half
    scr45 = kinlist(3)*half


    z4 = kinlist(4)
    z5 = kinlist(5)

    z2 = one-z4
    z  = z2

    Pggq_s5 = Ca/(scr25*scr45) + Ca*(z2/z5/scr25/scr24) + (two*Cf-Ca)*(z4/z5/scr45/scr24)
    Pggq_s5 = Pggq_s5 * Pqg_flip(z)

  end function Pggq_s5

end module mod_limits_dy_qg

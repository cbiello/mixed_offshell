module mod_limits_dy_ns
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions; use mod_procfunctions
  use mod_splitting
  implicit none

  private

  !-- soft
  public :: eik_s4,eik_s5,eik_ds,eik_dss5

  !-- soft-collinear
  public :: Pqq_c5_ds, Pqqg_s5
  public :: Pgg_c54_ds, Pgq_c54_ds

contains
  
  function Eik_s4(p)
    real(dp), intent(in) :: p(4,5)
    real(dp) :: Eik_s4
    real(dp) :: p1(4),p2(4),p4(4)

    p1 = p(:,1)
    p2 = p(:,2)
    p4 = p(:,3)

    Eik_s4 = two*Cf*scr(p1,p2)/scr(p1,p4)/scr(p2,p4)

    return 

  end function Eik_s4

  function Eik_s5(p1,p2,p4,p5)
    real(dp), intent(in) :: p1(4),p2(4),p4(4),p5(4)
    real(dp) :: Eik_s5
    real(dp) :: ES12_5, ES14_5, ES24_5, ES12_4 
    real(dp) :: s14, s15, s24, s25, s12, s45

    !   sij's are defined without factors of two here 

    s14 = -scr(p1,p4) 
    s15 = -scr(p1,p5) 
    s24 = -scr(p2,p4) 
    s25 = -scr(p2,p5) 
    s45 = scr(p4,p5) 
    s12 = scr(p1,p2) 

    ES12_5 = s12/s15/s25
    ES14_5 = s14/s15/s45
    ES24_5 = s24/s25/s45 

    Eik_s5 = (two*Cf-Ca)*ES12_5 + Ca*(ES14_5 + ES24_5 ) 

    return 

  end function Eik_s5

  function Eik_ds(p1,p2,p4,p5)
    real(dp), intent(in) :: p1(4), p2(4), p4(4), p5(4)
    real(dp) :: Eik_ds
    real(dp) :: ES12_4, ES12_5, ES12_45 
    real(dp) :: s14, s15, s24, s25, s12, s45, ES11_45, ES22_45
    real(dp) :: I11,I22,I12

    !-- note division by two, because of how things are defined in Catani's double-soft limits

    s14 = -scr(p1,p4) 
    s15 = -scr(p1,p5) 
    s24 = -scr(p2,p4) 
    s25 = -scr(p2,p5) 
    s45 = scr(p4,p5) 
    s12 = scr(p1,p2) 


    ES12_4 = s12/s14/s24
    ES12_5 = s12/s15/s25

    ES12_45 = one/s45**2 * ( s14*s25 + s15*s24 )/( ( s14+s15 ) *( s24+s25 ) )  &
         - s12**2/( two*s14*s25*s15*s24 ) * ( two  - ( s14*s25 + s15*s24 )/( (s14 + s15)* ( s24+s25 ) ) ) &
         + s12/two/s45 * ( two/s14/s25 + two/s24/s15  & 
         - one/( s14+s15 )/(s24 + s25 ) *( four + ( s14*s25 + s15*s24 )**2/s14/s25/s15/s24  ) )

    ES11_45 = two*s14*s15/s45**2/(s14+s15)**2

    ES22_45 = two*s24*s25/s45**2/(s24+s25)**2 


    Eik_ds  = 4.0_dp*Cf**2 * ES12_4 * ES12_5 * tag_ns_gg

    Eik_ds = Eik_ds + Ca*Cf*( two*ES12_45 - ES11_45 - ES22_45 ) * tag_ns_gg

    !-- nf contribution

    I12 = -s12/(s45*(s14+s15)*(s24+s25))-(s14-s15)*(s24-s25)/(two*s45**2*(s14+s15)*(s24+s25))
    I11 = -(s14-s15)**2/(two*s45**2*(s14+s15)**2)
    I22 = -(s24-s25)**2/(two*s45**2*(s24+s25)**2)

    Eik_ds = Eik_ds + two * nf * Cf*tr*(I11+I22-two*I12) * tag_ns_nf 


    return 

  end function Eik_ds

  function Eik_dss5(p1,p2,p4,p5)
    real(dp), intent(in) :: p1(4), p2(4), p4(4), p5(4)
    real(dp) :: Eik_dss5
    real(dp) :: ES12_5, ES14_5, ES24_5, ES12_4 
    real(dp) :: s14, s15, s24, s25, s12, s45

    !   sij's are defined without factors of two here 

    s14 = -scr(p1,p4) 
    s15 = -scr(p1,p5) 
    s24 = -scr(p2,p4) 
    s25 = -scr(p2,p5) 
    s45 = scr(p4,p5) 
    s12 = scr(p1,p2) 

    ES12_5 = s12/s15/s25
    ES14_5 = s14/s15/s45
    ES24_5 = s24/s25/s45 

    ES12_4 = s12/s14/s24

    Eik_dss5 = ( (two*Cf-Ca)*ES12_5 + Ca*(ES14_5 + ES24_5 ) )*two*Cf*ES12_4 

    return 

  end function Eik_dss5

!   !-- soft-collinear limits

  function Pqq_c5_ds(listkin,p1,p2,p4) 
    real(dp), intent(in) :: listkin(5),p1(4),p2(4),p4(4) 
    real(dp) :: Pqq_c5_ds
    real(dp) :: s15,z5 

    s15 = listkin(2)
    z5 = listkin(5)

    Pqq_c5_ds = -two/z5/s15 * four*Cf**2 *scr(p1,p2)/scr(p4,p1)/scr(p4,p2)

  end function Pqq_c5_ds

  function Pgg_c54_ds(listkin)
    real(dp), intent(in) :: listkin(6)
    real(dp) :: Pgg_c54_ds
    real(dp) :: s12,s14,s24,s45,z5,kt1kt2sq
    real(dp) :: e4_over_e45

    s12 = listkin(1)
    s14 = listkin(2)
    s24 = listkin(3)
    
    s45 = listkin(4)
    z5  = listkin(5)

    kt1kt2sq = listkin(6)

    e4_over_e45 = one-z5

    Pgg_c54_ds = two/s45 * ( pgg_spav(z5) + (kt1kt2sq-half)*(four*Ca*z5*(one-z5)) ) &
         * four * Cf * e4_over_e45**2 * s12/s14/s24

  end function Pgg_c54_ds

  function Pgq_c54_ds(listkin)
    real(dp), intent(in) :: listkin(6)
    real(dp) :: Pgq_c54_ds
    real(dp) :: s12,s14,s24,s45,z5,kt1kt2sq
    real(dp) :: e4_over_e45

    s12 = listkin(1)
    s14 = listkin(2)
    s24 = listkin(3)
    
    s45 = listkin(4)
    z5  = listkin(5)

    kt1kt2sq = listkin(6)

    e4_over_e45 = one-z5

    Pgq_c54_ds = two/s45 * ( pgq_spav(z5) + (kt1kt2sq-half)*(-four*tr*z5*(one-z5)) ) &
         * four * Cf * e4_over_e45**2 * s12/s14/s24

  end function Pgq_c54_ds

  function Pqqg_s5(kinlist) 
    real(dp) :: Pqqg_s5 
    real(dp), intent(in) :: kinlist(5)
    real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, z, Pqqonemz

    s14 = kinlist(1) 
    s15 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 

    z = one-z4

    Pqqonemz = one + z**2

    Pqqg_s5 = Cf*four*((two*Cf-Ca)*one/s14/s15/z4/z5 + Ca*(one/s15/s45/z4 + one/s14/s45/z5))*Pqqonemz

    Pqqg_s5 = Pqqg_s5/z

  end function Pqqg_s5

end module mod_limits_dy_ns

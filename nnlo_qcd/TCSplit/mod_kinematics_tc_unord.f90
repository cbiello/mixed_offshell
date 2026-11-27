module mod_kinematics_tc_unord
  use mod_types; use mod_consts_dp
  use mod_auxfunctions
  implicit none
  private
  
  public :: kinematics_tcl_unord_A, kinematics_tcl_unord_B, kinematics_tcl_unord_C, kinematics_tcl_unord_D

contains 

  ! triple-collinear kinematics for qg initial state; the collinear direction is ``2''
  
  subroutine kinematics_tcl_unord_A(z,yr,smallangle,full,fullc51,soft5,soft5c5)
    real(dp), intent(in) :: z,yr(4)
    real(dp), intent(out) :: full(6),fullc51(6),soft5(6),soft5c5(6)
    real(dp), intent(out) :: smallangle
    real(dp) :: x3,x4,lambda,E2,E4,E5,r,s24,s25,s45,rho24,rho25,rho45,z1,z4,z5

    r =  yr(1) 
    x4 = yr(2) 
    lambda = yr(3) 
    x3 = yr(4)

    E2 = one

    !  full triple collinear kinematics

    E4 = E2 * (one-z)*(one-r)
    E5 = E2 * (one-z)*r

    z4 = -(one-z)*(one-r)/z
    z5 = -(one-z)*r/z
    z1 = one/z

    rho24 = two ! *x3 will give the full  result 
    rho25 = x4  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    smallangle = rho25

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    full = (/s24,s25,s45,z4,z5,z1/) 

    !  full triple collinear, strongly ordered angles (x4->0)  

    rho24 = two ! *x3 will give the full  result 
    rho25 = one  ! *x3*x4 will give the full result 
    rho45 = two  ! x4 =0

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    fullc51 = (/s24,s25,s45,z4,z5,z1/) 

    !  ---- soft 5

    E4 = E2*(one-z)
    E5 = E2*(one-z)

    z4 = -(one-z)/z
    z5 = -(one-z)/z  ! r removed 
    z1 = one/z

    rho24 = two ! *x3 will give the full  result 
    rho25 = x4  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    soft5 = (/s24,s25,s45,z4,z5,z1/) 

    !  --- soft 5 collinear 5 

    rho24 = two ! *x3 will give the full  result 
    rho25 = one  ! *x3*x4 will give the full result 
    rho45 = two

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    soft5c5 = (/s24,s25,s45,z4,z5,z1/) 

  end subroutine kinematics_tcl_unord_A


  subroutine kinematics_tcl_unord_B(z,yr,smallangle,full,fullc54,soft5,soft5c5)
    real(dp), intent(in) :: z,yr(4)
    real(dp), intent(out) :: full(6),fullc54(6),soft5(6),soft5c5(6)
    real(dp), intent(out) :: smallangle
    real(dp) :: x3,x4,lambda,E2,E4,E5,r,s24,s25,s45,rho24,rho25,rho45,z1,z4,z5

    r =  yr(1) 
    x4 = yr(2) 
    lambda = yr(3) 
    x3 = yr(4)

    E2 = one

    !  full triple collinear kinematics ----------------------

    E4 = E2*(one-z)*(one - r)
    E5 = E2*(one-z)*r

    z4 = -(one-z)*(one-r)/z
    z5 = -(one-z)*r/z
    z1 = one/z

    rho24 = two ! *x3 will give the full  result 
    rho25 = two*(one-half*x4)  ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)

    smallangle = rho45

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    full = (/s24,s25,s45,z4,z5,z1/) 

    !   4 || 5

    rho24 = two ! *x3 will give the full  result 
    rho25 = two  ! *x3*x4 will give the full result 
    rho45 = one  ! x4 =0

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    fullc54 = (/s24,s25,s45,z4,z5,z1/) 

    !  soft 5  -------------------------------------

    E4 = E2*(one-z)
    E5 = E2*(one-z)

    z4 = -(one-z)/z
    z5 = -(one-z)/z  ! r removed 
    z1 = one/z

    rho24 = two ! *x3 will give the full  result 
    rho25 = two*(one-half*x4)  ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    soft5 = (/s24,s25,s45,z4,z5,z1/) 

    !  soft 5 and 5 || 4

    rho24 = two ! *x3 will give the full  result 
    rho25 = two  ! *x3 will give the full result 
    rho45 = one

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    soft5c5 = (/ s24,s25,s45,z4,z5,z1/) 

  end subroutine kinematics_tcl_unord_B


  subroutine kinematics_tcl_unord_C(z,yr,smallangle,full,fullc41,soft5,soft5c41)
    real(dp), intent(in) :: z,yr(4)
    real(dp), intent(out) :: full(6),fullc41(6),soft5(6),soft5c41(6)
    real(dp), intent(out) :: smallangle
    real(dp) :: x3,x4,lambda,E2,E4,E5,r,s24,s25,s45,rho24,rho25,rho45,z1,z4,z5

    r =  yr(1) 
    x4 = yr(2) 
    lambda = yr(3) 
    x3 = yr(4)

    E2 = one

    !---------------------------------------------------------  full triple collinear kinematics ;

    E4 = E2*(one-z)*(one-r)
    E5 = E2*(one-z)*r

    z4 = -(one-z)*(one-r)/z
    z5 = -(one-z)*r/z
    z1 = one/z

    rho24 = x4 ! *x3 will give the full  result 
    rho25 = two  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    smallangle = rho24

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    full = (/s24,s25,s45,z4,z5,z1/)    

    !-----------------------------------------------------------  full triple collinear, collinear 4

    rho24 = one ! *x3 will give the full  result 
    rho25 = two  ! *x3*x4 will give the full result 
    rho45 = two  ! x4 =0

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    fullc41 = (/s24,s25,s45,z4,z5,z1/) 

    !   soft 5

    E4 = E2*(one-z)   
    E5 = E2*(one-z)

    z4 = -(one-z)/z
    z5 = -(one-z)/z  ! r removed 
    z1 = one/z

    rho24 = x4 ! *x3 will give the full  result 
    rho25 = two  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    soft5 = (/s24,s25,s45,z4,z5,z1/) 

    !   soft 5 and 4 || 1 

    rho24 = one ! *x3 will give the full  result 
    rho25 = two  ! *x3 will give the full result 
    rho45 = two

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    soft5c41 = (/s24,s25,s45,z4,z5,z1/) 

  end subroutine kinematics_tcl_unord_C


  subroutine kinematics_tcl_unord_D(z,yr,smallangle,full,fullc54,soft5,soft5c5) 
    real(dp), intent(in) ::  z,yr(4)
    real(dp), intent(out) :: full(6),fullc54(6),soft5(6),soft5c5(6)
    real(dp), intent(out) :: smallangle
    real(dp) :: x3,x4,lambda,E2,E4,E5,r,s24,s25,s45,rho24,rho25,rho45,z1,z4,z5

    r =  yr(1) 
    x4 = yr(2) 
    lambda = yr(3) 
    x3 = yr(4)

    E2 = one

    !  full triple collinear kinematics ----------------------

    E4 = E2*(one-z)*(one - r)
    E5 = E2*(one-z)*r

    z4 = -(one-z)*(one-r)/z
    z5 = -(one-z)*r/z
    z1 = one/z

    rho24 = two*(one-half*x4)  ! *x3 will give the full  result 
    rho25 = two  ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)

    smallangle = rho45

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    full = (/s24,s25,s45,z4,z5,z1/) 

    !  full triple collinear, 4||5 ---------------------------------

    rho24 = two ! *x3 x4->0 will give the full  result 
    rho25 = two  ! *x3  will give the full result 
    rho45 = one  ! x4 =0

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    fullc54 = (/s24,s25,s45,z4,z5,z1/)

    !  full kinematics,  with E5 (r->0) ->0 -------------------------------------

    E4 = E2*(one-z)
    E5 = E2*(one-z)

    z4 = -(one-z)/z
    z5 = -(one-z)/z  ! r removed 
    z1 = one/z

    rho24 = two*(one-half*x4) ! *x3 will give the full  result 
    rho25 = two   ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    soft5 = (/s24,s25,s45,z4,z5,z1/) 

    !  full kinematics,  with E5->0  and collinear 5||4

    rho24 = two ! *x3 will give the full  result 
    rho25 = two  ! *x3 will give the full result 
    rho45 = one

    s24 = two*E2*E4*rho24
    s25 = two*E2*E5*rho25
    s45 = two*E4*E5*rho45

    soft5c5 = (/s24,s25,s45,z4,z5,z1/) 

  end subroutine kinematics_tcl_unord_D


end module mod_kinematics_tc_unord



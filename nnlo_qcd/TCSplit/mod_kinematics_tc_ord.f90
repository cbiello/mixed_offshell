module mod_kinematics_tc_ord
  use mod_types; use mod_consts_dp
  use mod_auxfunctions
  implicit none
  private

  public :: kinematics_tcl_a, kinematics_tcl_b, kinematics_tcl_c, kinematics_tcl_d

contains

  subroutine kinematics_tcl_a(z,yr,smallangle,full,fullc51,soft5,soft5c5,dsoft,dsord,dsoftc5,dsordc5,dsoftz,dsoftzc5)
    real(dp), intent(in) :: z, yr(4)
    real(dp), intent(out) :: full(6),fullc51(6),soft5(6),soft5c5(6),dsoft(6),dsord(6),dsoftc5(6),dsordc5(6),dsoftz(6),dsoftzc5(6)
    real(dp), intent(out) :: smallangle
    real(dp) :: x3,x4,lambda,E1,E4,E5,r,s14,s15,s45,rho14,rho15,rho45,z1,z4,z5

    r =  yr(1) 
    x4 = yr(2) 
    lambda = yr(3) 
    x3 = yr(4)

    E1 = one

    !  full triple collinear kinematics ;

    E4 = E1 * (one-z)*(one-half*r)
    E5 = E1 * (one-z)*half*r

    z4 = -(one-z)/z*(one-r/two)
    z5 = -(one-z)/z*r/two
    z1 = one/z

    rho14 = two ! *x3 will give the full  result 
    rho15 = x4  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    smallangle = rho15

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    full = (/ s14,s15,s45,z4,z5,z1/) 


    !  full triple collinear, strongly ordered angles (x4->0)  ---------------------------------

    ! energies & z's are the same as before 

    rho14 = two ! *x3 will give the full  result 
    rho15 = one  ! *x3*x4 will give the full result 
    rho45 = two  ! x4 =0

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    fullc51 = (/ s14,s15,s45,z4,z5,z1/) 


    !  ---- soft 5

    E4 = E1*(one-z)
    E5 = E1*(one-z)

    z4 = -(one-z)/z
    z5 = -(one-z)/z
    z1 = one/z


    rho14 = two ! *x3 will give the full  result 
    rho15 = x4  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    soft5 = (/ s14,s15,s45,z4,z5,z1/) 


    !  --- soft 5 collinear 5 


    rho14 = two ! *x3 will give the full  result 
    rho15 = one  ! *x3 will give the full result 
    rho45 = two

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    soft5c5 = (/ s14,s15,s45,z4,z5,z1/) 


    !  double soft   (xi ( not z)  limit )

    E4 = E1 * (one- half*r) 
    E5 = E1 *half * r

    z4 = -(one-half*r)
    z5 = -half*r
    z1 = one

    rho14 = two ! *x3 will give the full  result 
    rho15 = x4  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoft = (/ s14,s15,s45,z4,z5,z1/) 


    !  double soft  collinear 5 

    rho14 = two ! *x3 will give the full  result 
    rho15 = one  ! *x3 will give the full result 
    rho45 = two  

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoftc5 = (/ s14,s15,s45,z4,z5,z1/) 



    !  double soft   z 

    E4 = E1  
    E5 = E1 * r

    z4 = -one
    z5 = -r
    z1 = one

    rho14 = two ! *x3 will give the full  result 
    rho15 = x4  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoftz = (/ s14,s15,s45,z4,z5,z1/) 


    !  double soft  collinear 5 

    rho14 = two ! *x3 will give the full  result 
    rho15 = one  ! *x3 will give the full result 
    rho45 = two  

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoftzc5 = (/ s14,s15,s45,z4,z5,z1/) 




    !  double soft (xi), strongly ordered (r/2->0)


    E4 = E1
    E5 = E1

    z4 = -one
    z5 = -one
    z1 = one

    rho14 = two ! *x3 will give the full  result 
    rho15 = x4  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsord = (/ s14,s15,s45,z4,z5,z1/)    

    !  double soft (xi), strongly ordered (r/2 -> 0) & 5 collinear 

    rho14 = two ! *x3 will give the full  result 
    rho15 = one  ! *x3 will give the full result 
    rho45 = two

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsordc5 = (/ s14,s15,s45,z4,z5,z1/)    

  end subroutine kinematics_tcl_a



  subroutine kinematics_tcl_b(z,yr,smallangle,full,fullc54,soft5,soft5c5,dsoft,dsord,dsoftc5,dsordc5,dsoftz,dsoftzc5)
    real(dp), intent(in) :: z, yr(4)
    real(dp), intent(out) :: full(6), fullc54(6), soft5(6), soft5c5(6), dsoft(6) ,dsord(6), dsoftc5(6), dsordc5(6)
    real(dp), intent(out) :: dsoftz(6), dsoftzc5(6)
    real(dp), intent(out) :: smallangle
    real(dp) :: x3, x4, lambda, E1, E4, E5, r, s14,s15,s45,rho14,rho15,rho45,z1,z4,z5

    r =  yr(1) 
    x4 = yr(2) 
    lambda = yr(3) 
    x3 = yr(4)

    E1 = one

    !  full triple collinear kinematics ----------------------

    E4 = E1*(one-z)*(one - half*r)
    E5 = E1*(one-z)*(half*r)

    z4 = -(one-z)/z*(one-half*r)
    z5 = -(one-z)/z*(half*r)
    z1 = one/z

    rho14 = two ! *x3 will give the full  result 
    rho15 = two*(one-half*x4)  ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)

    smallangle = rho45

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    full = (/ s14,s15,s45,z4,z5,z1/) 


    !    4 || 5

    rho14 = two ! *x3 will give the full  result 
    rho15 = two  ! *x3*x4 will give the full result 
    rho45 = one  ! x4 =0

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    fullc54 = (/ s14,s15,s45,z4,z5,z1/) 



    !  soft 5  -------------------------------------

    E4 = E1*(one-z)
    E5 = E1*(one-z)

    z4 = -(one-z)/z  
    z5 = -(one-z)/z
    z1 = one/z


    rho14 = two ! *x3 will give the full  result 
    rho15 = two*(one-half*x4)  ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45


    soft5 = (/ s14,s15,s45,z4,z5,z1/) 


    !  soft 5 and 5 || 4

    rho14 = two ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = one

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    soft5c5 = (/ s14,s15,s45,z4,z5,z1/) 


    !  double soft 


    E4 = E1*(one - half*r)
    E5 = E1*(half*r)

    z4 = -(one-half*r)
    z5 = -half*r
    z1 = one


    rho14 = two ! *x3 will give the full  result 
    rho15 = two*(one-half*x4)  ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)


    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoft = (/ s14,s15,s45,z4,z5,z1/) 


    !  double soft triple collinear kinematics  5||4

    rho14 = two ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = one

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoftc5 = (/ s14,s15,s45,z4,z5,z1/) 




    !  double soft z


    E4 = E1
    E5 = E1*r

    z4 = -one
    z5 = -r
    z1 = one


    rho14 = two ! *x3 will give the full  result 
    rho15 = two*(one-half*x4)  ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)


    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoftz = (/ s14,s15,s45,z4,z5,z1/) 


    !  double soft  z triple collinear kinematics  5||4

    rho14 = two ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = one

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoftzc5 = (/ s14,s15,s45,z4,z5,z1/) 



    !  double soft strongly ordered triple collinear kinematics  -----

    E4 = E1
    E5 = E1  

    z4 = -one
    z5 = -one 
    z1 = one

    rho14 = two ! *x3 will give the full  result 
    rho15 = two*(one-half*x4)  ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45


    dsord = (/ s14,s15,s45,z4,z5,z1/)    


    !  double soft strongly ordered triple collinear angle ordered c5||4


    rho14 = two ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = one

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsordc5 = (/ s14,s15,s45,z4,z5,z1/)    


  end subroutine kinematics_tcl_b



  !  sector qqb_C

  subroutine kinematics_tcl_c(z,yr,smallangle,full,fullc41,soft5,soft5c41,dsoft,dsord,dsoftc41,dsordc41,dsoftz,dsoftzc41)
    real(dp), intent(in) :: z,  yr(4)
    real(dp), intent(out) :: full(6), fullc41(6), soft5(6), soft5c41(6), dsoft(6) ,dsord(6), dsoftc41(6), dsordc41(6) 
    real(dp), intent(out) :: smallangle, dsoftz(6), dsoftzc41(6)
    real(dp) :: x3, x4, lambda, E1, E4, E5, r, s14,s15,s45,rho14,rho15,rho45,z1,z4,z5

    r =  yr(1) 
    x4 = yr(2) 
    lambda = yr(3) 
    x3 = yr(4)

    E1 = one

    !  full triple collinear kinematics ;

    E4 = E1*(one-z)*(one-half*r)
    E5 = E1*(one-z)*half*r

    z4 = -(one-z)/z*(one-half*R)
    z5 = -(one-z)/z*(half*r)
    z1 = one/z

    rho14 = x4 ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    smallangle = rho14

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    full = (/ s14,s15,s45,z4,z5,z1/) 

    !  full triple collinear, collinear 4

    rho14 = one ! *x3 will give the full  result 
    rho15 = two  ! *x3*x4 will give the full result 
    rho45 = two  ! x4 =0

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    fullc41 = (/ s14,s15,s45,z4,z5,z1/) 


    !   soft 5

    E4 = E1*(one-z)   
    E5 = E1*(one-z)

    z4 = -(one-z)/z
    z5 = -(one-z)/z
    z1 = one/z

    rho14 = x4 ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    soft5 = (/ s14,s15,s45,z4,z5,z1/) 

    !   soft 5 and 4 || 1 

    rho14 = one ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = two

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    soft5c41 = (/ s14,s15,s45,z4,z5,z1/) 


    !  double soft triple collinear kinematics  -------------------------------------------------

    E4 = E1*(one - half*r)
    E5 = E1*half*r

    z4 = -(one-half*r)
    z5 = -half*r
    z1 = one

    rho14 = x4 ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoft = (/ s14,s15,s45,z4,z5,z1/) 


    !  double soft triple c4 || 1


    rho14 = one ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = two  

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoftc41 = (/ s14,s15,s45,z4,z5,z1/) 



    !  double soft z  triple collinear kinematics  -------------------------------------------------

    E4 = E1
    E5 = E1*r

    z4 = -one
    z5 = -r
    z1 = one

    rho14 = x4 ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoftz = (/ s14,s15,s45,z4,z5,z1/) 

    !  double soft triple c4 || 1

    rho14 = one ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = two  

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoftzc41 = (/ s14,s15,s45,z4,z5,z1/) 



    !  double soft strongly ordered triple collinear kinematics  -----

    E4 = E1
    E5 = E1

    z4 = -one
    z5 = -one
    z1 = one


    rho14 = x4 ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = two*(one-half*x4)**2/NF1(x3,half*x4,lambda)

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsord = (/ s14,s15,s45,z4,z5,z1/)    


    !  double soft strongly ordered triple collinear angle ordered c4 -------


    rho14 = one ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = two

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsordc41 = (/ s14,s15,s45,z4,z5,z1/)    

  end subroutine kinematics_tcl_c

  !---------------------------------------------------------------------------------------------------------
  ! sector qqb_D
  !--------------------------------------------------------------------------------------------------------

  subroutine kinematics_tcl_d(z,yr,smallangle,full,fullc54,soft5,soft5c5,dsoft,dsord,dsoftc5,dsordc5,dsoftz,dsoftzc5)
    real(dp), intent(in) ::  z, yr(4)
    real(dp), intent(out) :: full(6), fullc54(6), soft5(6), soft5c5(6), dsoft(6) ,dsord(6), dsoftc5(6), dsordc5(6) 
    real(dp), intent(out) :: smallangle, dsoftz(6), dsoftzc5(6)
    real(dp) :: x3, x4, lambda, E1, E4, E5, r, s14,s15,s45,rho14,rho15,rho45,z1,z4,z5

    r =  yr(1) 
    x4 = yr(2) 
    lambda = yr(3) 
    x3 = yr(4)

    E1 = one

    !  full triple collinear kinematics ----------------------

    E4 = E1*(one-z)*(one - half*r)
    E5 = E1*(one-z)*half*r

    z4 = -(one-z)/z*(one-half*r)
    z5 = -(one-z)/z*(half*r)
    z1 = one/z

    rho14 = two*(one-half*x4)  ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)

    smallangle = rho45

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45


    full = (/ s14,s15,s45,z4,z5,z1/) 

    !  full triple collinear, 4||5 ---------------------------------

    rho14 = two ! *x3 x4->0 will give the full  result 
    rho15 = two  ! *x3  will give the full result 
    rho45 = one  ! x4 =0

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    fullc54 = (/ s14,s15,s45,z4,z5,z1/) 


    !  full kinematics,  with E5 (r->0) ->0 -------------------------------------

    E4 = E1*(one-z)
    E5 = E1*(one-z)

    z4 = -(one-z)/z
    z5 = -(one-z)/z
    z1 = one/z


    rho14 = two*(one-half*x4) ! *x3 will give the full  result 
    rho15 = two   ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)


    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45


    soft5 = (/ s14,s15,s45,z4,z5,z1/) 

    !  full kinematics,  with E5->0  and collinear 5||4

    rho14 = two ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = one

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    soft5c5 = (/ s14,s15,s45,z4,z5,z1/) 

    !  double soft triple collinear kinematics  -------------------------------------------------

    E4 = E1*( one -half*r)
    E5 = E1*( half*r )

    z4 = -(one-half*r)
    z5 = -half*r
    z1 = one

    rho14 = two*(one-half*x4) ! *x3 will give the full  result 
    rho15 = two   ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoft = (/ s14,s15,s45,z4,z5,z1/) 


    !  double soft triple collinear kinematics  strong angle c5-------


    rho14 = two ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = one

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoftc5 = (/ s14,s15,s45,z4,z5,z1/) 



    !  double soft z triple collinear kinematics  -------------------------------------------------

    E4 = E1
    E5 = E1 *r 

    z4 = -one
    z5 = -r
    z1 = one

    rho14 = two*(one-half*x4) ! *x3 will give the full  result 
    rho15 = two   ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoftz = (/ s14,s15,s45,z4,z5,z1/) 


    !  double soft triple collinear kinematics  strong angle c5-------


    rho14 = two ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = one

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsoftzc5 = (/ s14,s15,s45,z4,z5,z1/) 




    !  double soft strongly ordered triple collinear kinematics  -----

    E4 = E1
    E5 = E1

    z4 = -one
    z5 = -one
    z1 = one

    rho14 = two*(one-half*x4) ! *x3 will give the full  result 
    rho15 = two   ! *x3 will give the full result 
    rho45 = two*(half*x4)**2/NF1(x3,one-half*x4,lambda)

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsord = (/ s14,s15,s45,z4,z5,z1/)    


    !  double soft strongly ordered triple collinear angle ordered c5||4

    rho14 = two ! *x3 will give the full  result 
    rho15 = two  ! *x3 will give the full result 
    rho45 = one

    s14 = two*E1*E4*rho14
    s15 = two*E1*E5*rho15
    s45 = two*E4*E5*rho45

    dsordc5 = (/ s14,s15,s45,z4,z5,z1/)    

  end subroutine kinematics_tcl_d

end module mod_kinematics_tc_ord

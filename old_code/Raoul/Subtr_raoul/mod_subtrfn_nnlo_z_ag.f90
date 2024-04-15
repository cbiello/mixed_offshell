module mod_subtrfn_nnlo_z_ag
! integrated subtraction terms for a g->qqb  channel
! cf. Notes from Chiara and Federica
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_auxfunctions
  use mod_splittings_bare
  implicit none

  private 

  public :: ag_qqb_sub_aqq_zi,ag_qqb_sub_qgq_zi
  public :: ag_qqb_sub_qqbz1zbar2,ag_qqb_sub_aa1z2
  
contains

 
  function ag_qqb_sub_qgq_zi(Ec,mu,z,Qq)
    ! Integrated subtraction terms to multiply FLMqqbg[z.1_q,2_g,3,4|5_q ], cf. Eq. 522
    real(dp), intent(in)   :: Ec,mu,z,Qq
    real(dp)               :: ag_qqb_sub_qgq_zi
    real(dp)               :: Lmu,omz,lomz,Pqg(-1:1)

    Lmu = log(mu/two/Ec)
    omz = one-z
    lomz = log(omz)

    Pqg = PqgAP(z)
   
! scale dependent pieces
    ag_qqb_sub_qgq_zi = -two*Pqg(0)*Lmu

! scale independent pieces
    ag_qqb_sub_qgq_zi = ag_qqb_sub_qgq_zi + &
         two*(omz*z + lomz*Pqg(0)) 

    ag_qqb_sub_qgq_zi = ag_qqb_sub_qgq_zi * xn * Qq**2
    


  end function ag_qqb_sub_qgq_zi



  function ag_qqb_sub_aqq_zi(Ec,mu,z,eta52,eta52log)
    ! Integrated subtraction terms to multiply FLMqqbg[z.1_q,2_g,3,4|5_q ], cf. Eq. 522
    real(dp), intent(in)   :: Ec,mu,z,eta52,eta52log
    real(dp)               :: ag_qqb_sub_aqq_zi
    real(dp)               :: Lmu,omz,lomz,Pqg(-1:1),wt62

    Lmu = log(mu/two/Ec)
    omz = one-z
    lomz = log(omz)
    wt62 = one-eta52**2

    Pqg = PqgAP(z)
   
! scale dependent pieces
    ag_qqb_sub_aqq_zi = -two*Pqg(0)*Lmu

! scale independent pieces
    ag_qqb_sub_aqq_zi = ag_qqb_sub_aqq_zi + &
         two*(omz*z + lomz*Pqg(0)) +&
         Pqg(0)*log(eta52log/two)*wt62

    ag_qqb_sub_aqq_zi = ag_qqb_sub_aqq_zi * tr
    


  end function ag_qqb_sub_aqq_zi

   function ag_qqb_sub_qqbz1zbar2(mu,z1,z2,Ec,Qq)
     ! Integrated subtraction terms to multiply FLM[z1,zbar2]
     real(dp), intent(in) :: mu,z1,z2,Ec,Qq
     real(dp)             :: ag_qqb_sub_qqbz1zbar2
     real(dp)             :: Lmu,lomz1,lomz2


     Lmu = log(mu/two/Ec) + log(z1)/two + log(z2)/two
     lomz1 = log(one-z1)
     lomz2 = log(one-z2)


     ! Lmu^2
     ag_qqb_sub_qqbz1zbar2 = four*(one - two*z1 + two*z1**2)*(one - two*z2 + two*z2**2)*Lmu**2

     ! Lmu
     ag_qqb_sub_qqbz1zbar2 = ag_qqb_sub_qqbz1zbar2 + &
          ( -four*(z1 - z1**2 + z2 - four*z1*z2 + four*z1**2*z2 - z2**2 + four*z1*z2**2 - four*z1**2*z2**2 + &
          (one - two*z1 + two*z1**2)*(one - two*z2 + two*z2**2)*lomz1 + (one - two*z1 + two*z1**2)*(one - two*z2 + two*z2**2)*lomz2) )*Lmu

     ! mu independent
     ag_qqb_sub_qqbz1zbar2 = ag_qqb_sub_qqbz1zbar2 + &
          four*(-((-one + z1)*z1) + (one - two*z1 + two*z1**2)*lomz1)*(-((-one + z2)*z2) + (one - two*z2 + two*z2**2)*lomz2)

     ! overall factor Nc*TR*Qq^2
     ag_qqb_sub_qqbz1zbar2 = xn*TR*Qq**2 * ag_qqb_sub_qqbz1zbar2
     

   end function ag_qqb_sub_qqbz1zbar2



   function ag_qqb_sub_aa1z2(mu,z,Ec,Qq)
     ! Integrated subtraction terms to multiply FLM[z1,zbar2]
     real(dp), intent(in) :: mu,z,Ec,Qq
     real(dp)             :: ag_qqb_sub_aa1z2
     real(dp)             :: Lmu,lomz,lz,li2z,li3z,li3omz


     Lmu = log(mu/two/Ec) + log(z )/two

     lz = log(z)
     lomz = log(one-z)
     li2z = real(dilog2(z),kind=dp)
     li3omz = trilog(one-z)
     li3z   =  trilog(+z)



     ! Lmu^2
     ag_qqb_sub_aa1z2 = (four*(three + four/z - three*z - four*z**2 + 6.0_dp*lz*(one + z)))/three * Lmu**2

     ! Lmu
     ag_qqb_sub_aa1z2 = ag_qqb_sub_aa1z2 + & 
          Lmu * four*(20.0_dp + 57.0_dp*z + 9.0_dp*lz*z - 9.0_dp*lz**2*z - 6.0_dp*pisq*z - 39.0_dp*z**2 + 27.0_dp*lz*z**2 - 9.0_dp*lz**2*z**2 - 6.0_dp*pisq*z**2 - 38.0_dp*z**3 + &
          36.0_dp*li2z*z*(one + z) + 6.0_dp*lomz*(-four - three*z + three*z**2 + four*z**3))/(9.0_dp*z)

     ! mu independent
     ag_qqb_sub_aa1z2 = ag_qqb_sub_aa1z2 + & 
          (-864.0_dp*li3omz*z*(one + z) - 432.0_dp*li3z*z*(one + z) + 18.0_dp*lz**3*z*(one + z) + 216.0_dp*li2z*z*(one + three*z) - 27.0_dp*lz**2*z*(5.0_dp + 11.0_dp*z) - &
          72.0_dp*lomz**2*(-one + z)*(four + 7.0_dp*z + four*z**2) + 108.0_dp*lz*z*(one + three*z - four*lomz**2*(one + z)) + &
          24.0_dp*lomz*(-20.0_dp + (-57.0_dp + 6.0_dp*pisq)*z + (39.0_dp + 6.0_dp*pisq)*z**2 + 38.0_dp*z**3 - 36.0_dp*li2z*z*(one + z)) + &
          two*(206.0_dp + 1056.0_dp*z - 840.0_dp*z**2 - 422.0_dp*z**3 + 12.0_dp*pisq*(-two - three*z - three*z**2 + two*z**3)) + 432.0_dp*z*(one + z)*zeta3)/(54.0_dp*z)
     

     ! overall factor TR*Qq^2
     ag_qqb_sub_aa1z2 = TR*Qq**2 * ag_qqb_sub_aa1z2
     

   end function ag_qqb_sub_aa1z2
   
 

end module mod_subtrfn_nnlo_z_ag

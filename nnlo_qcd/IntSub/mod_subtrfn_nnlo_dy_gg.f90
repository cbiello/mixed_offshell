module mod_subtrfn_nnlo_dy_gg
! integrated subtraction terms for gg->qqb channel
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  implicit none
  private 
  
  public :: gg_sub_qgz124,gg_sub_qqbz1zb2

contains

  
  function gg_sub_qgz124(shat,mu2,z)
    real(dp), intent(in)      :: shat,mu2,z
    real(dp)                  :: gg_sub_qgz124
    real(dp)                  :: Lmu,Pqg

    Lmu = log(mu2/shat)
    
    ! scale dependent term
    gg_sub_qgz124 = -tr*(one-two*z+two*z**2)*Lmu
    ! scale independent term
    gg_sub_qgz124 = gg_sub_qgz124 + two*tr*z*(one-z) &
         +two*tr*log(one-z)*(one-two*z+two*z**2)

!    print *, "z->",z, ",s->",shat,",musq->",mu2
!    print *, "gg_sub_qgz124",gg_sub_qgz124
!    stop


  end function gg_sub_qgz124

  
  function gg_sub_qqbz1zb2(shat,mu2,z,zbar)
    real(dp), intent(in) :: shat,mu2,z,zbar
    real(dp)             :: gg_sub_qqbz1zb2
    real(dp)             :: Lmu,coeffLsq,coeffL,coeffnoL

    Lmu = log(mu2/shat) 
    ! L^2
    coeffLsq = tr**2*(one-two*z+two*z**2)*(one-two*zbar+two*zbar**2)
    ! L
    coeffL = -two*tr**2*(one-two*z+two*z**2)*(one-two*zbar+two*zbar**2)*(log(one-z)+log(one-zbar)) +&
         two*tr**2*( (z**2-z)*(one-two*zbar)**2 -zbar*(one-zbar) )
    ! scale independent
    coeffnoL = tr**2*( (two*z*(one-z) + two*log(one-z)*(one-two*z+two*z**2))* &
         (two*zbar*(one-zbar) + two*log(one-zbar)*(one-two*zbar+two*zbar**2)) )
    
         
    gg_sub_qqbz1zb2 = coeffLsq * Lmu**2 + coeffL * Lmu + coeffnoL
!    print *, "z->",z, ",zbar->",zbar,",s->",shat,",musq->",mu2
!    print *, coeffnoL,coeffL,coeffLsq
!    print *, "gg_sub_qqbz1zb2",gg_sub_qqbz1zb2
!    stop

  end function gg_sub_qqbz1zb2

end module mod_subtrfn_nnlo_dy_gg

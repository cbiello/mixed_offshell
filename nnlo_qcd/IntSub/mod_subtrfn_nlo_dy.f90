module mod_subtrfn_nlo_dy
! integrated subtraction terms for gg->ff channel, where ff=gg or qqb
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions
  use mod_splitting_aux

  implicit none

  private

  public :: qqb_g_sub_qqbz12,gq_q_sub_qqbz12

contains

  function qqb_g_sub_qqbz12(q2,mu2,z)
! -- output indices: -1=delta, 0=reg+pls, 1=pls    
    real(dp), intent(in)  :: q2,mu2,z
    real(dp)              :: qqb_g_sub_qqbz12(-1:1) 
    real(dp)              :: Lmu

    Lmu = log(mu2/q2)

    qqb_g_sub_qqbz12 = zero

    qqb_g_sub_qqbz12 = -calPp_qq(z) - PqqAP(z) * Lmu

    qqb_g_sub_qqbz12(-1) = qqb_g_sub_qqbz12(-1) + Cf/three*pisq   !-- remainder from virtual definition. 
                                                              !-- Half w.r.t. notes because multiplied by #legs

  end function qqb_g_sub_qqbz12

  function gq_q_sub_qqbz12(q2,mu2,z)
! -- output indices: -1=delta, 0=reg+pls, 1=pls    
    real(dp), intent(in)  :: q2,mu2,z
    real(dp)              :: gq_q_sub_qqbz12(-1:1)
    real(dp)              :: Lmu

    Lmu = log(mu2/q2)

    gq_q_sub_qqbz12 = zero

    gq_q_sub_qqbz12 = -calPp_qg(z) - PqgAP(z) * Lmu

  end  function gq_q_sub_qqbz12

end module mod_subtrfn_nlo_dy

    


    
    

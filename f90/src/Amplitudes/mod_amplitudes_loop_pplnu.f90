!-- loop amplituds, all from OpenLoops
!-- factorise out as/twopi [aem/twopi] and gs^2[eesq] per qcd[ew] loop and gluon[photon] emission
module mod_amplitudes_loop_pplnu
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_ol_interface
  use mod_amplitudes_tree_ppnul
  use mod_coupl
  use mod_auxfunctions
  use openloops
  implicit none
  private

  !-- for qcdloop amplitudes, use mine
  public :: res_qcdloop_qqb_wp

contains


  !-----------------------------------------------------------------
  !--- amplitudes without extra radiation
  !-----------------------------------------------------------------
  
  !-- res(1,:) = q qb -> e- e [dn,up]
  !-- res(2,:) = qb q -> e- e [dn,up]
  subroutine res_qcdloop_qqb_wp(p,res0,res1fin)
    real(dp), intent(in)  :: p(4,4)
    real(dp), intent(out) :: res0(2,2),res1fin(1,2)

    call res_tree_qqb_w(p,res0)
    res1fin(1,:) = -8._dp*Cf*res0(1,:)
    

  end subroutine res_qcdloop_qqb_wp

  
end module mod_amplitudes_loop_pplnu

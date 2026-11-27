module mod_consts_debug
  use mod_types; use common_def
  implicit none
  private

  integer, public :: callcount

  real(dp), public :: limcheck
  integer, public :: icount
  
  real(dp), public :: FintNLO_CheckLim(4)

  real(dp), public :: Fint_fullONLO_rv_CheckLim(6)
  real(dp), public :: Fint_c41only_rv_CheckLim(2)
  real(dp), public :: Fint_c42only_rv_CheckLim(2)

  character, public :: globalsec*(10)


  real(dp), public :: FintNNLO_CheckLim_wp_ns(16)
  real(dp), public :: FintNNLO_CheckLim_wp_qg(8)
  real(dp), public :: FintNNLO_CheckLim_wp_gq(8)
  real(dp), public :: FintNNLO_CheckLim_wp_gg(7)


  real(dp), public :: Fint_gg_ff_sub_GGG_z124_CheckLim(6),&
       Fint_gg_ff_sub_GGG_1z24_CheckLim(6),&
       Fint_gg_ff_sub_GGG_124_CheckLim(6)

  real(dp), public :: Fint_gg_ff_sub_QGQ_z124_CheckLim(6),&
       Fint_gg_ff_sub_GQQ_1z24_CheckLim(6)

  real(dp), public :: Fint_qg_qg_sub_QGQ_z124_CheckLim(2),&
       Fint_qg_qg_sub_QGQ_1z24_CheckLim(2), &
       Fint_qg_qg_sub_QGQ_124_CheckLim(2),&
       Fint_qg_qg_sub_QQBG_1z24_CheckLim(6)


  real(dp), public :: Fint_gq_qg_sub_GQQ_z124_CheckLim(2),&
       Fint_gq_qg_sub_GQQ_1z24_CheckLim(2), &
       Fint_gq_qg_sub_GQQ_124_CheckLim(2),&
       Fint_gq_qg_sub_QQBG_z124_CheckLim(6)

  real(dp), public :: Fint_gg_sub_z124_CheckLim(2),Fint_gg_sub_1z24_CheckLim(2)


  real(dp), public :: Fint_qqb_ns_sub_QQBG_z124_CheckLim(6),&
       Fint_qqb_ns_sub_QQBG_1z24_CheckLim(6),&
       Fint_qqb_ns_sub_QQBG_124_CheckLim(6),&
       Fint_qqb_ns_sub_GQQ_z124_CheckLim(2),&
       Fint_qqb_ns_sub_QGQ_1z24_CheckLim(2)

end module mod_consts_debug

  ! -- vn tags
  real(dp), public :: vn_taga   !-- tag_ns_nf
  real(dp), public :: vn_tagb   !-- finite [in the ns_nf tag if run without tag]
  real(dp), public :: vn_tagc   !-- whatever_singlet
  real(dp), public :: vn_tagab  !-- zero after integration
  real(dp), public :: vn_tagac  !-- ns_id
  real(dp), public :: vn_tagbc  !-- finite [in ns_id]
  real(dp), public :: vn_tagcd  !-- finite [in whatever_singlet]
  real(dp), public :: vn_tagce  !-- s_v
  real(dp), public :: vn_tagcf  !-- finite [in s_v]

  !-- non-singlet tags
  real(dp), public :: tag_ns_gg
  real(dp), public :: tag_ns_nf
  real(dp), public :: tag_ns_id
  real(dp), public :: tag_ns_s

  !-- singlet tags
  real(dp), public :: tag_s_s !-- pure singlet
  real(dp), public :: tag_s_v !-- vector

  ! -- tag for NLO channels
  integer, public :: do_nlo_qqb
  integer, public :: do_nlo_qg

  ! -- tag for NNLO channels
  integer, public :: do_nnlo_qqb_ns
  integer, public :: do_nnlo_qqb_ot

  ! -- abelian, non ablelian
  real(dp), public :: tagab 
  real(dp), public :: tagnab 

  ! -- 
  real(dp), public :: myq2min, myq2max

  ! -- tag for triple collinear splitting function
  real(dp), public :: tag_tcpls
  real(dp), public :: tag_tcdel
  real(dp), public :: tag_tcreg
  real(dp), public :: tag_nontc
  ! 
  real(dp), public :: tag_12
  real(dp), public :: tag_z12
  real(dp), public :: tag_z1zb2

  ! -- coupling for Z decay
  real(dp), public, save :: cL_lep
  real(dp), public, save :: cR_lep
  real(dp), public, parameter :: Q_lep = Qel


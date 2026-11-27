module mod_consts_dp
  use mod_types; use common_def
  implicit none
  private

  ! -------------------------------------------------------------------------

  ! --- numerical constants 
  real(dp), public, parameter :: pi =&
       & 3.141592653589793238462643383279502884197_dp
  real(dp), public, parameter :: twopi =&
       & 6.283185307179586476925286766559005768394_dp
  real(dp), public, parameter :: zeta2 =&
       & 1.644934066848226436472415166646025189219_dp
  real(dp), public, parameter :: zeta3 =&
       & 1.202056903159594285399738161511449990765_dp
  real(dp), parameter, public :: pisq =&
       & 9.869604401089358618834490999876151135314_dp
  real(dp), parameter, public :: eulergamma =&
       & 0.577215664901532860606512090082402431042_dp
 real(dp), parameter, public :: li4half = &              ! Li_4[1/2]
       0.5174790616738993863307581618988629456224_dp 
  real(dp), public, parameter :: half  = 0.5_dp, two = 2.0_dp
  real(dp), public, parameter :: zero  = 0.0_dp, one = 1.0_dp
  real(dp), public, parameter :: mone = -1._dp 
  real(dp), public, parameter :: three = 3._dp 
  real(dp), public, parameter :: four  = 4._dp 
  real(dp), public, parameter :: five  = 5._dp 
  real(dp), public, parameter :: six  =  6._dp 
  real(dp), public, parameter :: sqrt2 = &
       &1.4142135623730950488016887242096980785696718753769_dp 
  real(dp), public, parameter :: msqrt2 = &
       &-1.4142135623730950488016887242096980785696718753769_dp 
  real(dp), public, parameter :: sqrt3 = &
       &1.7320508075688772935274463415058723669428_dp
  real(dp), public, parameter :: log2=log(two)
  complex(dp), public, parameter :: ci = cmplx(zero,1.0_dp,kind=dp)
  complex(dp), public, parameter :: czero = cmplx(0.0_dp,0.0_dp,kind=dp)
  complex(dp), public, parameter :: cone = cmplx(1.0_dp,0.0_dp,kind=dp)

  ! --- color factors
  real(dp), public, parameter :: CF = four/three
  real(dp), public, parameter :: CA = three
  real(dp), public, parameter :: Tr = half
  real(dp), public, parameter :: xn = 3.0_dp
  real(dp), public, parameter :: xnsq = 9.0_dp
  real(dp), public, parameter :: V = xnsq-one !-- for MCFM routines

  ! --- averaging factors
  real(dp), public, parameter :: aveqq = one/four/xnsq
  real(dp), public, parameter :: avegg = one/four/(xnsq-one)**2
  real(dp), public, parameter :: aveqg = one/four/(xnsq-one)/xn

  ! --- couplings and parameters  
  real(dp), public, parameter :: mt = 173.2_dp
  real(dp), public, parameter :: mtsq = mt**2
  real(dp), public, parameter :: mb = zero
  real(dp), public, parameter :: mbsq = zero
  real(dp), public, parameter :: mh = 125.0_dp
  real(dp), public, parameter :: mhsq = mh**2
  real(dp), public, parameter :: mz = 91.1876_dp 
  real(dp), public, parameter :: mzsq = mz**2
  real(dp), public, save :: mw = 80.398_dp
!  real(dp), public, parameter :: mwsq = mw**2
  real(dp), public, parameter :: GaW = 2.1054_dp !-- MCFM W width
  real(dp), public, parameter :: GaZ = 2.4952_dp !-- MCFM Z width
  real(dp), public, parameter :: GaH = 0.0041650_dp !-- MCFM H width
  real(dp), public, parameter :: Gf = 1.16639d-5
 real(dp), public, save :: sinW2
  real(dp), public, save :: gwsq
  real(dp), public, save :: eesq
  real(dp), public, save :: alpha_ew
  real(dp), public, save :: mwsq
  real(dp), public, parameter :: vev = one/sqrt(Gf*sqrt2)
  real(dp), public, save :: cosW2
  real(dp), public, save :: sw
  real(dp), public, save :: cw
  real(dp), public, save :: ee

  !--
  real(dp), public, parameter :: Qup =  2.0_dp/3.0_dp
  real(dp), public, parameter :: Qdn = -1.0_dp/3.0_dp
  real(dp), public, parameter :: Qel = -one
  real(dp), public, parameter :: Vnu = half
  real(dp), public, parameter :: Anu = half
  real(dp), public, parameter :: Qw = one
  real(dp), public, save :: Vup
  real(dp), public, save :: Aup
  real(dp), public, save :: Vdn
  real(dp), public, save :: Adn
  real(dp), public, save :: Lup
  real(dp), public, save :: Rup
  real(dp), public, save :: Ldn
  real(dp), public, save :: Rdn
  real(dp), public, save :: cLup
  real(dp), public, save :: cRup
  real(dp), public, save :: cLdn
  real(dp), public, save :: cRdn

  real(dp), public, save :: Vel
  real(dp), public, save :: Ael
  real(dp), public, save :: Lel
  real(dp), public, save :: Rel
  real(dp), public, save :: Lnu
  real(dp), public, save :: Rnu
  real(dp), public, save :: cLel
  real(dp), public, save :: cRel
  real(dp), public, save :: cLnu
  real(dp), public, save :: cRnu
  !--  
  !Mixing matrix elements
  real(dp), public, parameter :: Vud  = one
  real(dp), public, parameter :: Vnue = one
  !--
  real(dp), public, save :: cLWud
  real(dp), public, save :: cLWnue
  !-- 
  real(dp), public, save :: cAWW
  real(dp), public, save :: cZWW


  !--
  real(dp), public, parameter :: GeVtoFb=0.389379d12

  !-- quark name definition, according to LHAPDF
  integer, public, parameter :: gl = 0
  integer, public, parameter :: dn = 1
  integer, public, parameter :: up = 2
  integer, public, parameter :: st = 3
  integer, public, parameter :: ch = 4
  integer, public, parameter :: bt = 5
  integer, public, parameter :: tp = 6 

  ! --- active flavors and QCD beta function
  integer,  public, parameter :: nup = 2
  integer,  public, parameter :: ndn = 3     ! b quarks dealt with separately in pdf so this should be set to 2 for everything to work
  integer,  public, parameter :: nf = nup + ndn 
  integer,  public, parameter :: nup_pdf = nup
  integer,  public, parameter :: ndn_pdf = ndn
  integer,  public, parameter :: nf_pdf  = nup_pdf+ndn_pdf
  real(dp), public, parameter :: b0=11.0_dp/6.0_dp * Ca - 2.0_dp/3.0_dp * tr * nf
  real(dp), public, parameter :: b1=(17.0_dp*Ca**2 - 5.0_dp*Ca*nf - 3.0_dp*Cf*nf)/6.0_dp

  ! -- color structures
  real(dp), public, parameter :: CaCf = CA*CF
  real(dp), public, parameter :: CfCf = CF**2
  real(dp), public, parameter :: CaCa = CA**2
  real(dp), public, parameter :: CaTr = Ca/two
  real(dp), public, parameter :: CfTr = Cf/two
  real(dp), public, parameter :: CfNf = CF*nf
  real(dp), public, parameter :: CaNf = CA*nf
  real(dp), public, parameter :: NfNf = nf*nf

  !-- process-specific constants
  include 'inc_procconsts_dp.f90'

end module mod_consts_dp

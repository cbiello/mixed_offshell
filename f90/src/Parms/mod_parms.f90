module mod_parms
  use mod_types
  use mod_consts_dp
  use mod_parser
  implicit none

  private

  !--------------------------------------------------------------------------
  !-- EW
  !--------------------------------------------------------------------------

  !-- variables
  real(dp), public, save :: mt,mb,mw,mz,mh
  real(dp), public, save :: mz_in,mw_in !-- for complex mass-scheme
  real(dp), public, save :: gaw,gaz,gah
  real(dp), public, save :: Gf,vev,alpha_0,alpha_mz,aem,eesq,eesq2,gwsq,sw2,cw2,ee,cw,sw,gw
  real(dp), public, save :: cLup,cRup,cLdn,cRdn !-- coupling in units of ee
  real(dp), public, save :: cLnu,cRnu,cLel,cRel
  complex(dp), public, save :: acu,vcu,acd,vcd,ace,vce,acn,vcn !-- for complex mass-scheme
  real(dp),    public, save :: acur,vcur,acdr,vcdr             !-- real-values, needed for UV renorm

  real(dp), public, save :: mtsq,mbsq,mwsq,mzsq,mhsq
  complex(dp), public, save :: mzsq_prop,mwsq_prop !-- msq - i*m*ga

  !chiara
  real(dp), public, save :: cLWud,cLWnue
  real(dp), public, save :: cLW

  ! raoul added

  ! warning: all testing done with VCKM = 1, so diagonal elements only!
  real(dp), public, parameter :: VCKM_ud = one
  real(dp), public, parameter :: VCKM_us = zero
  real(dp), public, parameter :: VCKM_ub = zero
  real(dp), public, parameter :: VCKM_cd = zero
  real(dp), public, parameter :: VCKM_cs = one
  real(dp), public, parameter :: VCKM_cb = zero
  real(dp), public, parameter :: VCKM_td = zero
  real(dp), public, parameter :: VCKM_ts = zero
  real(dp), public, parameter :: VCKM_tb = one
  real(dp), public, parameter, dimension(3,3) :: VCKM(3,3)=reshape( [VCKM_ud,VCKM_cd,VCKM_td,VCKM_us,VCKM_cs,VCKM_ts,VCKM_ub,VCKM_cb,VCKM_tb],[3,3])
  
  

  
  !-- same for complex mass scheme
  character(30), public, save :: ew_scheme
  logical, public, save :: cm_scheme
  logical, public, save :: only_up
  integer, public, save :: flag_down

  complex(dp), public, save :: cms_gwsq,cms_sw2,cms_cw2,cms_cw,cms_sw,cms_gw
  complex(dp), public, save :: cms_cLup,cms_cRup,cms_cLdn,cms_cRdn
  complex(dp), public, save :: cms_cLnu,cms_cRnu,cms_cLel,cms_cRel
  complex(dp), public, save :: cms_mwsq,cms_mzsq

  complex(dp), public, save :: lfav

  !chiara
  complex(dp), public, save :: cms_cLWud,cms_cLWnue
  complex(dp), public, save :: cms_cLW


  !-- fixed parameters
  real(dp), public, parameter :: I3up =  half
  real(dp), public, parameter :: I3dn = -half

  real(dp), public, parameter :: I3nu =  half
  real(dp), public, parameter :: I3el = -half

  real(dp), public, parameter :: Qup =  2.0_dp/3.0_dp
  real(dp), public, parameter :: Qdn = -1.0_dp/3.0_dp
  real(dp), public, parameter :: Qup2 = Qup**2, Qdn2=Qdn**2
  
  real(dp), public, parameter :: Qnu = 0._dp
  real(dp), public, parameter :: Qel = -1.0_dp
  real(dp), public, parameter :: Qnu2 = Qnu**2, Qel2=Qel**2

  !Mixing matrix elements
  real(dp), public, parameter :: Vnue = one

  !--------------------------------------------------------------------------
  !-- QCD
  !--------------------------------------------------------------------------

  !-- variables
  integer, public, save :: nup,ndn,nlq,nle,nne,nf
  real(dp), public, save :: b0,b1
  real(dp), public, save :: gamma_a !-- photon anomalous dimension

  !-- parameters
  real(dp), public, parameter :: Cf = four/three
  real(dp), public, parameter :: Ca = three
  real(dp), public, parameter :: tr = half
  real(dp), public, parameter :: xn = 3.0_dp
  real(dp), public, parameter :: xnsq = xn**2
  real(dp), public, parameter :: v = xnsq-one

  !--------------------------------------------------------------------------
  !-- averages
  !--------------------------------------------------------------------------

  real(dp), public, parameter :: aveqq = one/four/xnsq
  real(dp), public, parameter :: aveqg = one/four/xn/v
  real(dp), public, parameter :: avegg = one/four/v**2

  real(dp), public, parameter :: aveqa = one/four/xn
  real(dp), public, parameter :: aveaa = one/four
  real(dp), public, parameter :: avega = one/four/v
  
  !--------------------------------------------------------------------------
  !-- human numbers
  integer, public, parameter :: id_q = 1 !-- generic quark
  integer, public, parameter :: id_b = 5 !-- b-quark
  integer, public, parameter :: id_g = 21
  integer, public, parameter :: id_a = 22
  integer, public, parameter :: id_z = 23
  integer, public, parameter :: id_w = 24 !-- + is always particle
  integer, public, parameter :: id_h = 25
  !
  integer, public, parameter :: id_el = 11 !-- + is always particle
  integer, public, parameter :: id_mu = 13
  !chiara
  integer, public, parameter :: id_nue = 12 !-- neutrino e
  integer, public, parameter :: id_numu = 14 !-- neutrino mu
  integer, public, parameter :: id_nutau = 15 !-- neutrino tau 
  !--------------------------------------------------------------------------

  ! raoul added
  !-- charges for final state leptons
  real(dp), public, save  :: Q3
  real(dp), public, save  :: Q4
  !-- generation labels for quarks
  real(dp), public, parameter  :: Qgeneration(-5:5)=[3,2,2,1,1,0,1,1,2,2,3]
  real(dp), public, parameter  :: Q_IS(-5:7) = [-Qdn,-Qup,-Qdn,-Qup,-Qdn,zero,Qdn,Qup,Qdn,Qup,Qdn,Qup,zero]
  real(dp), public, parameter  :: Qsq_IS(-5:7) = [Qdn2,Qup2,Qdn2,Qup2,Qdn2,zero,Qdn2,Qup2,Qdn2,Qup2,Qdn2,Qup2,zero]
  ! ----- 
  public :: get_parms,help_parms
  public :: set_qcd_parms,set_ew_parms_cms
  public :: print_qcd_parms,print_ew_parms
  
contains

  subroutine get_parms()
    
    !-- QCD
    nup = int_val_opt('-nup',2)
    ndn = int_val_opt('-ndn',3)
    nlq = int_val_opt('-nlq',2)

    !-- Number of light leptons
    nle = int_val_opt('-nle',3)

    !-- Number of neutrinos
    nne = int_val_opt('-nne',3)

    !-- masses
    mt = real_val_opt('-mt',173.2_dp)
    mb = real_val_opt('-mb',zero)
    mz_in = real_val_opt('-mz',91.1876_dp) !-- for complex mass scheme 
    mw_in = real_val_opt('-mw',80.398_dp)  !-- for complex mass scheme
    mh = real_val_opt('-mh',125._dp)

    !-- widths
    gaw = real_val_opt('-gaw',2.1054_dp)
    gaz = real_val_opt('-gaz',2.4952_dp)
    gah = real_val_opt('-gah',0.004165_dp)

    !-- EW scheme
    ew_scheme = trim(string_val_opt('-ew_scheme','Gmu'))

    !-- complex mass scheme
    cm_scheme = log_val_opt_witharg('-cm_scheme',.true.)

    !-- select only up-quarks
    only_up = log_val_opt_witharg('-only_up',.false.)
    if(only_up) then
      flag_down = 0
    else
      flag_down = 1
    endif

    !-- couplings
    Gf  = real_val_opt('-Gf',1.16639E-5_dp)        !-- for Gmu scheme
    alpha_0  = real_val_opt('-alpha_0',1/137._dp)  !-- for alpha(0) scheme
    alpha_mz = real_val_opt('-alpha_mz',1/128._dp) !-- for alpha(0) scheme

    ! couplings of final state leptons
#if (_Vcharge == 0)
    Q3 = -one
    Q4 = +one
#elif  (_Vcharge == -1)
    Q3 = -one
    Q4 = zero
#elif  (_Vcharge == +1)
    Q3 = zero
    Q4 = +one
#endif
    
    
  end subroutine get_parms

  subroutine help_parms(idev)
    integer, intent(in) :: idev

    write(idev,*) ''
    write(idev,*) ' -nup 2'
    write(idev,*) ' -ndn 3'
    write(idev,*) ''
    write(idev,*) ' -ew_scheme Gmu'
    write(idev,*) ' -cm_scheme true --> complex mass scheme'
    write(idev,*) ''
    write(idev,*) ' -only_up true --> only contribution from up-type quarks, i.e. u,c'
    write(idev,*) ''
    write(idev,*) ' -m[t/b/z/w/h] 173.2 --> for mz,mw, this is mz_in, mw_in'
    write(idev,*) ' -ga[w/z/h] 2.1054'
    write(idev,*) ' -Gf 1.16639E-5'
    
  end subroutine help_parms
  
  subroutine set_qcd_parms()

    nf = nup + ndn

    b0 = 11.0_dp/6.0_dp * Ca-2.0_dp/3.0_dp*tr*nf
    b1 = (17.0_dp*Ca**2-5.0_dp*Ca*nf-3.0_dp*Cf*nf)/6.0_dp
        
  end subroutine set_qcd_parms

  subroutine print_qcd_parms(outdev)
    integer, intent(in) :: outdev

    write(outdev,*) '# Input nup  = ', nup
    write(outdev,*) '# Input ndn  = ', ndn
    write(outdev,*) '# Input nle  = ', nle
    write(outdev,*) '# Input nne  = ', nne
    write(outdev,*) '# Derived nf = ', nf
    write(outdev,*) '#'

  end subroutine print_qcd_parms

  subroutine set_ew_parms_cms(my_ew_scheme,outdev)
    character(*), intent(in) :: my_ew_scheme
    integer, intent(in) :: outdev

    real(dp) :: swMG

    mz = mz_in
    mw = mw_in

    mzsq_prop = mz**2 - ci * mz * gaz
    mwsq_prop = mw**2 - ci * mw * gaw
    
    !-- compute other parameters
    if (cm_scheme) then
       cms_mzsq = mzsq_prop
       cms_mwsq = mwsq_prop

       mz = mz/sqrt(one+gaz**2/mz**2)
       mw = mw/sqrt(one+gaw**2/mw**2)

       mzsq = mz**2
       mwsq = mw**2
    else
       cms_mzsq = mz**2
       cms_mwsq = mw**2

       mzsq = mz**2
       mwsq = mw**2
    endif
       
    mhsq = mh**2
    mtsq = mt**2
    mbsq = mb**2

    if (my_ew_scheme.eq.'Gmu') then     !-- G_mu

       vev  = one/sqrt(Gf*sqrt2)
       cms_gwsq = four*sqrt2*cms_mwsq*Gf 
       cms_cw2  = cms_mwsq/cms_mzsq
       cms_sw2  = one - cms_cw2
       aem  = sqrt2*Gf*abs(cms_mwsq*cms_sw2)/pi
       eesq = aem * four * pi
              
    elseif (my_ew_scheme.eq.'amz') then !-- aem(mz)

       aem = alpha_mz
       eesq = aem * four * pi
       cms_cw2  = cms_mwsq/cms_mzsq
       cms_sw2  = one-cms_mwsq/cms_mzsq
       Gf   = aem*pi/sqrt2/abs(cms_mwsq*cms_sw2)
       vev  = one/sqrt(Gf*sqrt2)
       cms_gwsq = four*sqrt2*cms_mwsq*Gf

    elseif (my_ew_scheme.eq.'a0') then  !-- alpha(0)

       aem = alpha_0
       eesq = aem * four * pi
       cms_cw2 = cms_mwsq/cms_mzsq
       cms_sw2 = one - cms_cw2
       Gf   = aem*pi/sqrt2/abs(cms_mwsq*cms_sw2)
       vev  = one/sqrt(Gf*sqrt2)
       cms_gwsq = four*sqrt2*cms_mwsq*Gf

    else
       write(outdev,*) 'ew scheme not implemented'
       stop       
    endif

    cms_sw = sqrt(cms_sw2)
    cms_cw = sqrt(cms_cw2)
    cms_gw = sqrt(cms_gwsq)

    ee = sqrt(eesq)
    eesq2 = eesq**2

    !--
    
    cms_cLup = (I3up-Qup*cms_sw2)/cms_cw/cms_sw
    cms_cRup = -Qup*cms_sw/cms_cw
    
    cms_cLdn = (I3dn-Qdn*cms_sw2)/cms_cw/cms_sw
    cms_cRdn = -Qdn*cms_sw/cms_cw
    
    cms_cLel = (I3el-Qel*cms_sw2)/cms_cw/cms_sw
    cms_cRel = -Qel*cms_sw/cms_cw
    
    cms_cLnu = I3nu/cms_cw/cms_sw
    cms_cRnu = zero

    !chiara
    cms_cLWud = 1/sqrt2/cms_sw!*Vud
    cms_cLW   = 1/sqrt2/cms_sw
    cms_cLWnue = 1/sqrt2/cms_sw*Vnue

    !-- now set the real values
    gwsq = real(cms_gwsq,dp)
    sw2  = real(cms_sw2,dp)
    cw2  = real(cms_cw2,dp)

    sw = real(cms_sw,dp)
    cw = real(cms_cw,dp)
    gw = real(cms_gw,dp)

    cLup = real(cms_cLup,dp)
    cRup = real(cms_cRup,dp)
    cLdn = real(cms_cLdn,dp)
    cRdn = real(cms_cRdn,dp)
    cLel = real(cms_cLel,dp)
    cRel = real(cms_cRel,dp)
    cLnu = real(cms_cLnu,dp)
    cRnu = real(cms_cRnu,dp)

    !chiara
    !--
    swMG=sqrt( 0.22224648578577766_dp) 
    !cLWud  = 1/sqrt2/swMG!*Vud
    !cLW    = 1/sqrt2/swMG
    !cLWnue = 1/sqrt2/swMG*Vnue
    cLWud  = real(cms_cLWud,dp)
    cLW  = real(cms_cLW,dp)
    cLWnue = real(cms_cLWnue,dp)


    !-- axial and vector couplings, complex-mass scheme
    acd = I3dn/(2*cms_cw*cms_sw)
    acu = I3up/(2*cms_cw*cms_sw)
    ace = I3el/(2*cms_cw*cms_sw)
    acn = I3nu/(2*cms_cw*cms_sw)
    vcd = (I3dn - 2*cms_sw2*Qdn)/(2*cms_cw*cms_sw)
    vcu = (I3up - 2*cms_sw2*Qup)/(2*cms_cw*cms_sw)
    vce = (I3el - 2*cms_sw2*Qel)/(2*cms_cw*cms_sw)
    vcn = (I3nu - 2*cms_sw2*Qnu)/(2*cms_cw*cms_sw)

    !-- Ligh-fermion a/v couplings. 
    !-- Multiplicity for closed fermion loops with exernal Z bosons
    lfav = xn * (nup*(acu**2 + vcu**2) + ndn*(acd**2 + vcd**2)) &
         + nle*(ace**2 + vce**2) + nne*(acn**2 + vcn**2)

    !-- axial and vector couplings, real values
    acdr = I3dn/(2*cw*sw)
    acur = I3up/(2*cw*sw)
    vcdr = (I3dn - 2*sw2*Qdn)/(2*cw*sw)
    vcur = (I3up - 2*sw2*Qup)/(2*cw*sw)

    !-- photon anomalous dimension
    gamma_a = -two/3 * (xn*(nup*Qup2+ndn*Qdn2) + nle*Qel2)

  end subroutine set_ew_parms_cms
  
  subroutine print_ew_parms(my_ew_scheme,outdev)
    character(*), intent(in) :: my_ew_scheme
    integer, intent(in) :: outdev

    write(outdev,*) '# PROCESS  = ', _Vcharge
    write(outdev,*) '#'

    write(outdev,*) '# complex mass scheme = ', cm_scheme
    write(outdev,*) '#'

    write(outdev,*) '# only_up = ', only_up
    write(outdev,*) '#'
    
    if (my_ew_scheme.eq.'Gmu') then
       write(outdev,*) '# Using the Gmu scheme'
       write(outdev,*) '# Input Gf      = ', Gf
       write(outdev,*) '# Derived 1/aem = ', one/aem
       write(outdev,*) '# '
    elseif (my_ew_scheme.eq.'amz') then
       write(outdev,*) '# Using the alpha(mz) scheme'
       write(outdev,*) '# Input 1/aem = ', one/aem
       write(outdev,*) '# Derived Gf  = ', Gf
       write(outdev,*) '# '
    elseif (my_ew_scheme.eq.'a0') then
       write(outdev,*) '# Using the alpha(0) scheme'
       write(outdev,*) '# Input 1/aem = ', one/aem
       write(outdev,*) '# Derived Gf  = ', Gf
       write(outdev,*) '# '
    else
       write(outdev,*) 'ew scheme not implemented'
       stop
    endif

    write(outdev,*) '# Input mz_in = ', mz_in
    write(outdev,*) '# Input mw_in = ', mw_in
    write(outdev,*) '# Input mz = ', mz
    write(outdev,*) '# Input mw = ', mw
    write(outdev,*) '# Input mh = ', mh
    write(outdev,*) '# Input mb = ', mb
    write(outdev,*) '# Input mt = ', mt
    write(outdev,*) '# Input Gf = ', Gf
    write(outdev,*) '# Input gaz = ', gaz
    write(outdev,*) '# Input gaw = ', gaw
    write(outdev,*) '# Input gah = ', gah
    write(outdev,*) '#'
    write(outdev,*) '# Derived vev ', vev
    write(outdev,*) '# Derived gw  ', gw
    write(outdev,*) '# Derived sw2 ', sw2
    write(outdev,*) '# Derived cw2 ', cw2
    write(outdev,*) '#'

  end subroutine print_ew_parms
  
end module mod_parms

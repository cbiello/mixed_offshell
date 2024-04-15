!-- phase space for q[1] qb[2] -> l[3] lb[4] g[5] a[6]
module mod_kinematicsNNLO_EunordAord_ga
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_process
  use mod_auxfunctions
  use mod_aux_kinematicsNNLO_EunordAord_ga
  use mod_kinematics_gen
  use mod_kinematics_lept
  implicit none
  integer, public, parameter :: NNLOdim_vegas = tauy_max + 2 + 3 + 3
  !
  integer, public, parameter :: kNNLO_min = tauy_max + 1
  integer, public, parameter :: kNNLO_max = tauy_max + 2 + 3 + 3
  integer, public, parameter :: kNNLO_max_full = kNNLO_max + 2
  !
  integer, public, parameter :: xE4 = kNNLO_min
  integer, public, parameter :: xE5 = kNNLO_min + 1
  integer, public, parameter :: xRHO1 = kNNLO_min + 2
  integer, public, parameter :: xRHO2 = kNNLO_min + 4
  logical, parameter :: withmax = .true.
  private

  public ::  kinematics_nnlo_5161a,kinematics_nnlo_5161c
  public ::  kinematics_nnlo_5262a, kinematics_nnlo_5262c
  public ::  kinematics_nnlo_5261, kinematics_nnlo_5162

  public :: kinematics_nnlo_5163, kinematics_nnlo_5164
  public :: kinematics_nnlo_5263, kinematics_nnlo_5264




contains 

  !-- we factor out an extra 1/(8*pisq)**2 from the phase space
  !-- which multiplied by gs**4 in the ME reconstructs asontwopi**2

  subroutine kinematics_nnlo_5161a(mv2,y,yr,&
       HardProc,S5Lim,S6Lim,S5S6Lim,TCLim,TCS5Lim,TCS6Lim,TCS5S6Lim,&
       C6Lim,C6S5Lim,C6S6Lim,C6S5S6Lim,TCC6Lim,TCC6S5Lim,TCC6S6Lim,TCC6S5S6Lim)
    real(dp), intent(in) ::  mv2,y,yr(10)
    type(KinConfig)      ::  HardProc,S6Lim,TCLim,TCS6Lim,C6Lim,C6S6Lim,TCC6Lim,TCC6S6Lim
    type(KinConfig)      ::  S5Lim,S5S6Lim,TCS5Lim,C6S5Lim,TCC6S5Lim,TCS5S6Lim,TCC6S5S6Lim,C6S5S6Lim
    !---
    real(dp) :: mv
    real(dp) :: x1,x2,x3,x4,x5,lambda,sinphi45sign
    real(dp) :: cos41,sin41,phi41,cos51,sin51,sinphi45,cosphi45
    real(dp) :: mysin4lim,mysin5lim
    real(dp) :: eta41,eta51,eta42,eta52,eta54,eta56,eta57
    real(dp) :: s12,s41,s51,s42,s52,s54,s4_15,s5_14
    real(dp) :: z,zz,z4,z5,zz1,zz2,ztot,kt1kt2sq,nperp(4),ncoll(4)
    real(dp) :: spart,sqrts,xi1,xi2,yy
    real(dp) :: p1(4),p2(4),p4(4),p5(4),pV(4),E4,E5,n4(3),n5(3)
    real(dp) :: kallenF,nlept(3,2),nlept1(3),nlept2(3)
    real(dp) :: weight

    call initialize_config(HardProc)       !1
    call initialize_config(S5Lim)          !2
    call initialize_config(S6Lim)          !3
    call initialize_config(S5S6Lim)        !4
    call initialize_config(TCLim)          !5
    call initialize_config(TCS5Lim)        !6
    call initialize_config(TCS6Lim)        !7
    call initialize_config(TCS5S6Lim)      !8
    call initialize_config(C6Lim)          !9
    call initialize_config(C6S5Lim)        !10
    call initialize_config(C6S6Lim)        !11
    call initialize_config(C6S5S6Lim)      !12
    call initialize_config(TCC6Lim)        !13
    call initialize_config(TCC6S5Lim)      !14
    call initialize_config(TCC6S6Lim)      !15
    call initialize_config(TCC6S5S6Lim)    !16

    mv =  sqrt(mv2)

    x1 = yr(1) 
    x2 = yr(2) 
    x3 = yr(3) 
    x4 = yr(5)
    x5 = yr(6)

    lambda = sin(half*pi*x5)**2

    !-- the sign of sinphi45 is not well defined -> randomize it

    if (yr(10) .ge. half) then
       sinphi45sign = +one
    else
       sinphi45sign = -one
    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    cos41 = one - two*x3
    sin41 = sqrt(one - cos41**2)
    phi41 = twopi*yr(4)

    cos51 = one - two*x3*x4
    sin51 = sqrt(one - cos51**2)

    if(2.0_dp*(one - x4)**2/NF1(x3,x4,lambda)-two*(one+x4*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    sinphi45 = sinphi45sign*sinphi45 

    !-- eta variables
    eta41 = x3
    eta42 = one - eta41
    eta51 = x3*x4
    eta52 = one - eta51
    eta54 = x3*(one-x4)**2/NF1(x3,x4,lambda)

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,HardProc%flag)

    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin41*cos(phi41),sin41*sin(phi41),cos41/)
       n5 = (/sin51*(cos(phi41)*cosphi45-sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45),cos51/)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,HardProc)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,5) = p4 ! gluon 
       HardProc%AmpMom(:,6) = p5 ! photon

       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)

       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight & !-- partition
            * one/two/spart &
            * spart/MV2

       HardProc%ids = [0,0,id_el,-id_el,id_g,id_a]
       HardProc%npart = 6

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: S5 S6 * FLM  : both gluon and photon soft
    !-------------------------------------------------------------------------------------------------------------

    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    S5S6Lim%PartFrac = (/xi1,xi2/)

    p1=  half*sqrts*(/one,zero,zero, one/) 
    p2 = half*sqrts*(/one,zero,zero,-one/) 
    E4 = half*sqrts*x1
    E5 = half*sqrts*x2
    n4 = (/sin41*cos(phi41),sin41*sin(phi41),cos41/)
    n5 = (/sin51*(cos(phi41)*cosphi45-sin(phi41)*sinphi45), & 
         sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45),cos51/)
    p4 = E4*(/one,n4(1),n4(2),n4(3)/)
    p5 = E5*(/one,n5(1),n5(2),n5(3)/)

    pV = p1 + p2

    call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,S5S6Lim)

    S5S6Lim%AmpMom(:,1) = p1
    S5S6Lim%AmpMom(:,2) = p2

    S5S6Lim%LimMom(:,1) = p1
    S5S6Lim%LimMom(:,2) = p2
    S5S6Lim%LimMom(:,3:4) = S5S6Lim%AmpMom(:,3:4)
    S5S6Lim%LimMom(:,5) = p4  !gluon
    S5S6Lim%LimMom(:,6) = p5  !photon

    nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
    nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)

    eta56 = half*(one-dot_product(n5,nlept1))
    eta57 = half*(one-dot_product(n5,nlept2))

    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

    S5S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
         * x1*x2 & 
         * weight &
         * one/two/spart &
         * spart/MV2

    S5S6Lim%ids(1:4) = [0,0,id_el,-id_el]
    S5S6Lim%npart = 4

    ! ------   TRIPLE COLLINEAR KINEMATICS + S5S6  -------
    ! TC = C156 S5 S6 : both gluon and photon collinear to 1,both gluon and photon soft
    ! x3 = 0, x4 = 0, x2 = 0, x1 = 0

    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero

    p1=  half*sqrts*(/one,zero,zero,one/)
    p2 = half*sqrts*(/one,zero,zero,-one/)
    E4 = half*sqrts*x1
    E5 = half*sqrts*x2
    n4 = (/zero,zero,one/)
    n5 = (/zero,zero,one/)
    p4 = E4*(/one,n4(1),n4(2),n4(3)/)
    p5 = E5*(/one,n5(1),n5(2),n5(3)/)

    !see Kosta. thesis z5 = E5/(E4+E5-E1), z4 = E4/(E4+E5-E1), z = 1-z4-z5 -> see eq.B23
    ztot = - one 
    z4 = x1/ztot 
    z5 = x2/ztot

    !sij = 2pi.pj = 2Ei.Ej.rhoij = 4Ei.Ej.etaij
    s41 = spart*x1*x3
    s51 = spart*x2*x3*x4
    s54 = spart*x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

    TCS5S6Lim%Lim_KinInv(1:5) = (/s41,s51,s54,z4,z5/)

    pV = p1 + p2

    call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,TCS5S6Lim)

    nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
    nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

    eta56 = half*(one-dot_product(n5,nlept1))
    eta57 = half*(one-dot_product(n5,nlept2))

    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

    ! triple collinear operators do not modify the unresolved phase space
    ! sector functions are affected by the limit
    TCS5S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
         * x1*x2 & 
         * weight &
         * one/two/spart &
         * spart/MV2

    TCS5S6Lim%ids = [0,0,id_el,-id_el,0,0]
    TCS5S6Lim%npart = 4
    
    ! --------  C5 S5S6 FLM  ----------
    ! x4 ->0

    eta41 = x3
    eta51 = zero
    eta42 = (one-x3)
    eta52 = one
    eta54 = x3/NF1(x3,zero,lambda)

    p1 = half*sqrts*(/one,zero,zero,one/)
    p2 = half*sqrts*(/one,zero,zero,-one/)

    E4 = half*sqrts*x1
    E5 = half*sqrts*x2
    n4 = (/sin41*cos(phi41),  sin41*sin(phi41), cos41 /)
    n5 = (/zero,zero,one/)

    p4 = E4*(/one,n4(1),n4(2),n4(3)/)
    p5 = E5*(/one,n5(1),n5(2),n5(3)/)

    !print*, 'momentum of the photon for C5S5 S4', p5

    ztot = - one
    z4 = x1/ztot
    z5 = x2/ztot

    !sij = 2pi.pj = 2Ei.Ej.rhoij = 4Ei.Ej.etaij
    s41 = spart*x1*x3
    s51 = spart*x2*x3*x4
    s54 = spart*x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

    C6S5S6Lim%Lim_KinInv(1:5) = (/s41,s51,s54,z4,z5/)
    pV = p1 + p2 

    call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,C6S5S6Lim)

    C6S5S6Lim%LimMom(:,1) = p1
    C6S5S6Lim%LimMom(:,2) = p2
    C6S5S6Lim%LimMom(:,3) = p4


    nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
    nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

    eta56 = half*(one-dot_product(n5,nlept1))
    eta57 = half*(one-dot_product(n5,nlept2))

    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

    C6S5S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1*x2 & 
         * weight &       
         * one/two/spart &
         * spart/MV2

    C6S5S6Lim%ids = [0,0,id_el,-id_el,0,0]
    C6S5S6Lim%npart = 4
    
    ! -------   TRIPLE COLLINEAR KINEMATICS  + C5 + S5S6 --------------
    ! x3 = 0, x4 = 0, x2 = 0, x1 = 0

    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero

    !sij = 2pi.pj = 2Ei.Ej.rhoij = 4Ei.Ej.etaij
    s41 = spart*x1*x3
    s51 = spart*x2*x3*x4
    s54 = spart*x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

    !see Kosta. thesis z5 = E5/(E4+E5-E1), z4 = E4/(E4+E5-E1), z = 1-z4-z5 -> see eq.B23
    ztot =-one 
    z4 = x1/ztot 
    z5 = x2/ztot

    TCC6S5S6Lim%Lim_KinInv(1:5) = (/s41,s51,s54,z4,z5/)
    pV = p1 + p2

    call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,TCS5S6Lim)  !! get leptons momenta

    nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
    nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

    eta56 = half*(one-dot_product(n5,nlept1))
    eta57 = half*(one-dot_product(n5,nlept2))

    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

    TCC6S5S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1*x2 & 
         * weight &
         * one/two/spart &
         * spart/MV2

    TCC6S5S6Lim%ids = [0,0,id_el,-id_el,0,0]
    TCC6S5S6Lim%npart = 4

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5, gluon is soft, angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------

    eta41 = x3
    eta51 = x3*x4
    eta42 = (one - x3)
    eta52 = (one - x3*x4)
    eta54 = x3*(one-x4)**2/NF1(x3,x4,lambda)

    call get_flag_nlo(x2,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,S5Lim%flag)

    if (S5Lim%flag .eqv. .false.) then 
       S5Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin41*cos(phi41),sin41*sin(phi41),cos41/)
       n5 = (/sin51*(cos(phi41)*cosphi45-sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45),cos51/)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pV = p1 + p2 - p5
       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,5) = p5 !-- photon

       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = p4

       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)

       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

       S5Lim%ids = [0,0,id_el,-id_el,id_a,0]
       S5Lim%npart = 5

    endif

    ! --------- C6 S5 ---------
    ! C61 S5 : photon collinear, gluon soft
    ! x4 = 0, x1 = 0

    eta41 = x3
    eta42 = one - eta41
    eta51 = zero
    eta52 = one
    eta54 = x3/NF1(x3,zero,lambda)

    call get_flag_nlo(x2,MV2,y,x2,eta51,eta52,xi1,xi2,spart,sqrts,C6S5Lim%flag)

    if (C6S5Lim%flag .eqv. .false.) then 

       C6S5Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/)
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin41*cos(phi41),sin41*sin(phi41),cos41/)
       n5 = (/zero,zero,one/)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pV = p1 + p2 - p5

       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,C6S5Lim)

       C6S5Lim%AmpMom(:,1) = p1-p5
       C6S5Lim%AmpMom(:,2) = p2

       C6S5Lim%LimMom(:,1) =  p1
       C6S5Lim%LimMom(:,2) =  p2
       C6S5Lim%LimMom(:,3) =  p4

       ! see Kosta. thesis z = (E1-E5)/E1, z4 = E4/(E4+E5-E1) eq.B12
       s51 = spart*x2*(x3*x4)
       z5 = one/(one - x2) 

       C6S5Lim%Lim_KinInv(1:5) = (/zero,s51,zero,zero,z5/)

       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)

       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

       C6S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

       C6S5Lim%ids = [0,0,id_el,-id_el,0,0]
       C6S5Lim%npart = 4
       
    endif

    ! --------- TC C6 S5 ---------
    ! x3 = 0, x4 = 0, x1 = 0, eta61 << eta51

    eta41 = zero
    eta42 = one - eta41
    eta51 = zero
    eta52 = one - eta51
    eta54 = x3/NF1(zero,zero,lambda)

    call get_flag_nlo(x2,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,TCC6S5Lim%flag)

    if (TCC6S5Lim%flag .eqv. .false.) then 

       TCC6S5Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero,zero,one/)
       n5 = (/zero,zero,one/)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pV = p1 + p2 - p5

       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,TCC6S5Lim)

       TCC6S5Lim%AmpMom(:,1) = p1-p5
       TCC6S5Lim%AmpMom(:,2) = p2

       TCC6S5Lim%LimMom(:,1) =  p1
       TCC6S5Lim%LimMom(:,2) =  p2
       TCC6S5Lim%LimMom(:,3) =  p4

       z5 = one/(one - x2) 
       s41 = spart*x1*x3
       s51 = spart*x2*x3*x4

       TCC6S5Lim%Lim_KinInv(1:5) = (/s41,s51,zero,z4,z5/)

       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)

       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

       TCC6S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

       TCC6S5Lim%ids = [0,0,id_el,-id_el,0,0]
       TCC6S5Lim%npart = 4

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S6 + C6 S6   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------

    ! --------- SOFT 6 ---------
    !  S6: photon is soft

    eta41 = x3
    eta51 = x3*x4
    eta42 = (one-x3)
    eta52 = (one - x3*x4)
    eta54 = x3*(one-x4)**2/NF1(x3,x4,lambda)

    call get_flag_nlo(x1,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts,S6Lim%flag)

    if (S6Lim%flag .eqv. .false.) then 

       S6Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin41*cos(phi41),sin41*sin(phi41),cos41/)
       n5 = (/sin51*(cos(phi41)*cosphi45-sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45),cos51/)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pV = p1 + p2 - p4
       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,S6Lim)

       S6Lim%AmpMom(:,1) = p1
       S6Lim%AmpMom(:,2) = p2
       S6Lim%AmpMom(:,5) = p4

       S6Lim%LimMom(:,1) = p1
       S6Lim%LimMom(:,2) = p2
       S6Lim%LimMom(:,3) = S6Lim%AmpMom(:,3) 
       S6Lim%LimMom(:,4) = S6Lim%AmpMom(:,4) 
       S6Lim%LimMom(:,5) = p5 !photon

       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)

       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

       S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2   

       S6Lim%ids = [0,0,id_el,-id_el,id_g,0]
       S6Lim%npart = 5

       ! kin inv and weight for C6 +  S6
       ! C61 S6 : photon collinear, photon soft
       ! x4 = 0, x2 = 0

       eta41 = x3
       eta51 = zero
       eta42 = (one-x3)
       eta52 = one
       eta54 = x3/NF1(x3,zero,lambda)

       s51 = spart*x2*(x3*x4)
       z5 = x2
       C6S6Lim%Lim_KinInv(1:5) = (/zero,s51,zero,zero,z5/)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

       C6S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

       C6S6Lim%ids = [0,0,id_el,-id_el,id_g,0]
       C6S6Lim%npart = 5

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    ! --------- TRIPLE COLLINEAR -------
    eta41 = zero
    eta42 = one - eta41
    eta51 = zero
    eta52 = one - eta51
    eta54 = zero

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,TCLim%flag)

    if (TCLim%flag .eqv. .false.) then

       TCLim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero,zero,one/)
       n5 = (/zero,zero,one/)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       !see Kosta. thesis z5 = E5/(E4+E5-E1), z4 = E4/(E4+E5-E1), z = 1-z4-z5 -> see eq.B23
       z5 = x2/(x1 + x2 - one)
       z4 = x1/(x1 + x2 - one)
       z = one/(one - x1- x2)

       s41 = spart*x1*x3
       s51 = spart*x2*x3*x4
       s54 = spart*x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

       TCLim%Lim_KinInv(1:7) = (/s41,s51,s54,z4,z5,x1,x2/)

       pV = p1 + p2 - p4 - p5 

       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,TCLim)  !! get leptons momenta

       TCLim%AmpMom(:,1) = p1 - p4 - p5
       TCLim%AmpMom(:,2) = p2

       if(2.0_dp*(one - x4)**2/NF1(x3,x4,lambda)-two*(one+x4).gt.0) then
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif
       sinphi45 = sinphi45sign*sinphi45     ! randomize sign of sinphi45

       mysin4lim = two             !*sqrt(x3)(1+O(x3))
       mysin5lim = two * sqrt(x4)  !*sqrt(x3)(1+O(x3))

       TCLim%LimMom(:,1) = mysin4lim*(/one,cos(phi41),sin(phi41),cos41/)
       TCLim%LimMom(:,2) = mysin5lim*(/one,sin51*(cos(phi41)*cosphi45-sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45),cos51/)
       TCLim%LimMom(:,3) = (/one,zero,zero,one/)

       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)

       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

       TCLim%ids = [0,0,id_el,-id_el,0,0]
       TCLim%npart = 4

       ! --------- TRIPLE COLLINEAR + C6  -------
       ! TC = C156 C16 : both gluon and photon collinear to 1, photon super collinear
       ! x3 = 0, x4 = 0, eta61 << eta51

       s41 = spart * x1*x3 
       s51 = spart * x2*(x3*x4)
       s4_15 = s41 * (one-x2)

       zz1 = x2/(x2-one)
       zz2 = (one-x2)/(one-x1-x2)

       kt1kt2sq = (-one+two*lambda)**2 * (one-x3)

       TCC6Lim%Lim_KinInv(1:5) = (/s51,s4_15,zz1,zz2,kt1kt2sq/)

       TCC6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

       TCC6Lim%ids = [0,0,id_el,-id_el,0,0]
       TCC6Lim%npart = 4
       
    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 5
    !-------------------------------------------------------------------------------------------------------------

    ! -------- TRIPLE COLLINEAR + SOFT 5

    ! TC = C156 S5 : both gluon and photon collinear to 1, gluon soft
    ! x3 = 0, x4 = 0, x1 = 0

    eta41 = zero
    eta42 = one - eta41
    eta51 = zero
    eta52 = one - eta51
    eta54 = zero

    call get_flag_nlo(x2,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,TCS5Lim%flag)

    if (TCS5Lim%flag .eqv. .false.) then 

       TCS5Lim%PartFrac = (/xi1,xi2/)

       p1=  half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero,zero,one/)
       n5 = (/zero,zero,one/)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       z4 = x1 
       z5 = one/(one-x2)

       s41 = spart*x1*x3
       s51 = spart*x2*x3*x4
       s54 = spart*x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

       TCS5Lim%Lim_KinInv(1:5) = (/s41,s51,s54,z4,z5/)

       pV = p1 + p2 - p5 

       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,TCS5Lim) 

       TCS5Lim%AmpMom(:,1) = p1 - p5
       TCS5Lim%AmpMom(:,2) = p2

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

       TCS5Lim%ids = [0,0,id_el,-id_el,0,0]
       TCS5Lim%npart = 4

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 6
    !-------------------------------------------------------------------------------------------------------------

    ! -------- TRIPLE COLLINEAR + SOFT 6
    ! TC = C156 S6 : both gluon and photon collinear to 1, photon soft
    ! x3 = 0, x4 = 0, x2 = 0

    call get_flag_nlo(x1,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts,TCS6Lim%flag)

    if (TCS6Lim%flag .eqv. .false.) then

       TCS6Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero,zero,one/)
       n5 = (/zero,zero,one/)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       z4 = x1 
       z5 = x2

       !sij = 2pi.pj = 2Ei.Ej.rhoij = 4Ei.Ej.etaij
       s41 = spart*x1*x3
       s51 = spart*x2*x3*x4
       s54 = spart*x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

       TCS6Lim%Lim_KinInv(1:5) = (/s41,s51,s54,z4,z5/)

       pV = p1 + p2 - p4

       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,TCS6Lim)  !! get leptons momenta

       TCS6Lim%AmpMom(:,1) = p1 - p4
       TCS6Lim%AmpMom(:,2) = p2

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

       TCS6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

       TCS6Lim%ids = [0,0,id_el,-id_el,0,0]
       TCS6Lim%npart = 4

       ! -------   TRIPLE COLLINEAR KINEMATICS  + C6 + SOFT6 ---------
       z4 = one/(one-x1) 
       z5 = x2

       s41 = spart *x1*x3 
       s51 = spart *x2*(x3*x4)
       s54 = spart *x1*x2*x3/NF1(x3,zero,lambda)

       TCC6S6Lim%Lim_KinInv(1:5) = (/s41,s51,s54,z4,z5/)

       TCC6S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

       TCC6S6Lim%ids = [0,0,id_el,-id_el,0,0]
       TCC6S6Lim%npart = 4
       
    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C6 FLM 
    !-------------------------------------------------------------------------------------------------------------

    phi41 = twopi*yr(4) 
    cos41 = one - two*x3 
    sin41 = sqrt(one-cos41**2)

    eta41 = x3
    eta51 = zero
    eta42 = (one-x3) 
    eta52 = one
    eta54 = x3/NF1(x3,zero,lambda)

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,C6Lim%flag)

    if (C6Lim%flag .eqv. .false.) then 

       C6Lim%PartFrac = (/xi1,xi2/)
       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin41*cos(phi41), sin41*sin(phi41), cos41 /)
       n5 = (/zero,zero,one/)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       s51 = spart*x2*(x3*x4)
       z5 = one/(one - x2) 

       C6Lim%Lim_KinInv(1:5) = (/zero,s51,zero,zero,z5/)

       pV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,C6Lim)

       C6Lim%AmpMom(:,1) = p1 - p5
       C6Lim%AmpMom(:,2) = p2
       C6Lim%AmpMom(:,5) = p4

       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)

       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

       C6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

       C6Lim%ids = [0,0,id_el,-id_el,id_g,0]
       C6Lim%npart = 5

    endif

  end subroutine kinematics_nnlo_5161a

  subroutine kinematics_nnlo_5161c(MV2,Y,yr, &
       HardProc,S5Lim,S6Lim,S5S6Lim,TCLim,TCS5Lim,TCS6Lim,TCS5S6Lim,&
       C5Lim,C5S5Lim,C5S6Lim,C5S5S6Lim,TCC5Lim,TCC5S5Lim,TCC5S6Lim,TCC5S5S6Lim)
    real(dp), intent(in) ::  mv2,y,yr(10)
    type(KinConfig)      ::  HardProc,S6Lim,TCLim,TCS6Lim,C5Lim,C5S6Lim,TCC5Lim,TCC5S6Lim
    type(KinConfig)      ::  S5Lim,S5S6Lim,TCS5Lim,C5S5Lim,TCC5S5Lim,TCS5S6Lim,TCC5S5S6Lim,C5S5S6Lim
    !
    real(dp) :: mv
    real(dp) :: x1,x2,x3,x4,x5,lambda,sinphi45sign
    real(dp) :: cos41,sin41,phi41,cos51,sin51,sinphi45,cosphi45
    real(dp) :: mysin4lim,mysin5lim
    real(dp) :: eta41,eta51,eta42,eta52,eta54,eta56,eta57
    real(dp) :: s12,s41,s51,s42,s52,s54,s4_15,s5_14
    real(dp) :: z,zz,z4,z5,zz1,zz2,ztot,kt1kt2sq,nperp(4),ncoll(4)
    real(dp) :: spart,sqrts,xi1,xi2,yy
    real(dp) :: p1(4),p2(4),p4(4),p5(4),pV(4),E4,E5,n4(3),n5(3)
    real(dp) :: kallenF,nlept(3,2),nlept1(3),nlept2(3)
    real(dp) :: weight

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S5Lim)          !2
    call initialize_config(S6Lim)          !3
    call initialize_config(S5S6Lim)        !4
    call initialize_config(TCLim)          !5
    call initialize_config(TCS5Lim)        !6
    call initialize_config(TCS6Lim)        !7
    call initialize_config(TCS5S6Lim)      !8
    call initialize_config(C5Lim)          !9
    call initialize_config(C5S5Lim)        !10
    call initialize_config(C5S6Lim)        !11
    call initialize_config(C5S5S6Lim)      !12
    call initialize_config(TCC5Lim)        !13
    call initialize_config(TCC5S5Lim)      !14
    call initialize_config(TCC5S6Lim)      !15
    call initialize_config(TCC5S5S6Lim)    !16

    HardProc%npart = 6
    S5Lim%npart = 5
    S6Lim%npart = 5
    S5S6Lim%npart = 4
    TCLim%npart = 4
    TCS5Lim%npart = 4
    TCS6Lim%npart = 4
    TCS5S6Lim%npart = 4
    C5Lim%npart = 5
    C5S5Lim%npart = 5
    C5S6Lim%npart = 4
    C5S5S6Lim%npart = 4
    TCC5Lim%npart = 4 
    TCC5S5Lim%npart = 4
    TCC5S6Lim%npart = 4
    TCC5S5S6Lim%npart = 4

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a] 
    S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_a] 
    S6Lim%ids(1:5) = [0,0,id_el,-id_el,id_g] 
    S5S6Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCLim%ids(1:4) = [0,0,id_el,-id_el] 
    TCS5Lim%ids(1:4) = [0,0,id_el,-id_el]
    TCS6Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCS5S6Lim%ids(1:4) = [0,0,id_el,-id_el]
    C5Lim%ids(1:5) = [0,0,id_el,-id_el,id_a] 
    C5S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]
    C5S6Lim%ids(1:4) = [0,0,id_el,-id_el] 
    C5S5S6Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCC5Lim%ids(1:4) = [0,0,id_el,-id_el]  
    TCC5S5Lim%ids(1:4) = [0,0,id_el,-id_el]
    TCC5S6Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCC5S5S6Lim%ids(1:4) = [0,0,id_el,-id_el]
    
    
    !-- MV2 and Y are generated
    MV =  sqrt(MV2)
    
    x1 = yr(1) 
    x2 = yr(2) 
    x3 = yr(3) 
    x4 = yr(5)
    
    if (yr(10) .ge. half) then
       sinphi45sign = +one
    else
       sinphi45sign = -one
    endif
    
    lambda = sin(half*pi*yr(6))**2

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi41 = twopi*yr(4)
    
    cos41 = one - two*x3*x4
    sin41 = sqrt(one-cos41**2)
    
    cos51 = one - two*x3 
    sin51 = sqrt(one-cos51**2)
    
    if(2.0_dp*(one - x4)**2/NF1(x3,x4,lambda)-two*(one+x4*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif
    
    sinphi45 = sinphi45sign*sinphi45     ! randomize sign of sinphi45
    
    eta41 = x3*x4
    eta51 = x3
    eta42 = (one-eta41) 
    eta52 = (one - eta51)
    eta54 = x3*(one-x4)**2/NF1(x3,x4,lambda)
    
    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,HardProc%flag)
    
    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       
       n4 = (/sin41*cos(phi41),sin41*sin(phi41),cos41/)
       n5 = (/sin51*(cos(phi41)*cosphi45-sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45),cos51/)
       
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)
       
       pV = p1 + p2 - p4 - p5       
       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,HardProc)

       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)

       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,5) = p4 ! gluon 
       HardProc%AmpMom(:,6) = p5 ! photon

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: S5S6 FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! -------  S5 S6 ------------

    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    S5S6Lim%PartFrac = (/xi1,xi2/)

     p1 = half*sqrts*(/one,zero,zero,one/) 
     p2 = half*sqrts*(/one,zero,zero,-one/)

     E4 = half*sqrts*x1
     E5 = half*sqrts*x2

     n4 = (/sin41*cos(phi41),sin41*sin(phi41),cos41/)
     n5 = (/sin51*(cos(phi41)*cosphi45-sin(phi41)*sinphi45), & 
          sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45),cos51/)

     p4 = E4*(/one,n4(1),n4(2),n4(3)/)
     p5 = E5*(/one,n5(1),n5(2),n5(3)/)

     pV = p1 + p2
     call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,S5S6Lim)
     
     nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
     nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)
     
     eta56 = half*(one-dot_product(n5,nlept1))
     eta57 = half*(one-dot_product(n5,nlept2))

     call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)
     
     S5S6Lim%AmpMom(:,1) = p1
     S5S6Lim%AmpMom(:,2) = p2
     
     S5S6Lim%LimMom(:,1) = p1
     S5S6Lim%LimMom(:,2) = p2
     S5S6Lim%LimMom(:,3) = S5S6Lim%AmpMom(:,3)
     S5S6Lim%LimMom(:,4) = S5S6Lim%AmpMom(:,4)
     S5S6Lim%LimMom(:,5) = p4 
     S5S6Lim%LimMom(:,6) = p5
     
     S5S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
          * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
          * x1*x2 & 
          * weight &
          * one/two/spart & 
          * spart/MV2
     
     ! -----   TRIPLE COLLINEAR KINEMATICS +  S5S6   ----------

     eta41 = zero
     eta51 = zero
     eta42 = one
     eta52 = one
     eta54 = zero
     
     p1 = half*sqrts*(/one,zero,zero,one/) 
     p2 = half*sqrts*(/one,zero,zero,-one/)
     E4 = half*sqrts*x1
     E5 = half*sqrts*x2
     n4 = (/zero, zero, one /)
     n5 = (/zero, zero, one /)
     p4 = E4*(/one,n4(1),n4(2),n4(3)/)
     p5 = E5*(/one,n5(1),n5(2),n5(3)/)
     
     pV = p1 + p2
     call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,TCS5S6Lim)
     
     nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
     nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)
     
     eta56 = half*(one-dot_product(n5,nlept1))
     eta57 = half*(one-dot_product(n5,nlept2))

     call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

     s41 = spart * x1*(x3*x4) 
     s51 = spart * x2*x3
     s54 = spart * x1 * x2 * x3 * (one-x4)**2/NF1(x3,x4,lambda)
     
     ztot = - one 
     z4 = x1/ztot 
     z5 = x2/ztot
     
     TCS5S6Lim%Lim_KinInv(1:5) = (/s41,s51,s54,z4,z5/)
     
     TCS5S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
          * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
          * x1*x2 & 
          * weight &
          * one/two/spart &
          * spart/MV2
     
     ! -------  C51 S5S6 FLM  ----------------
     
     phi41 = twopi*yr(4)
     
     cos41 = one 
     sin41 = zero
     
     cos51 = one - two*x3 
     sin51 = sqrt(one-cos51**2)
     
     eta41 = zero
     eta51 = x3
     eta42 = one 
     eta52 = (one - eta51)
     eta54 = x3/NF1(x3,zero,lambda)

     if(one-two*lambda .ge. zero) then
        sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
        cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
     else
        sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
        cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
     endif
     
     sinphi45 = sinphi45sign*sinphi45     ! randomize sign of sinphi45
     
     p1 = half*sqrts*(/one,zero,zero,one/) 
     p2 = half*sqrts*(/one,zero,zero,-one/)
     
     E4 = half*sqrts*x1
     E5 = half*sqrts*x2
     
     n4 = (/zero, zero, one/)
     n5 = (/sin51*(cos(phi41)*cosphi45-sin(phi41)*sinphi45), & 
          sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45),cos51/)
     
     p4 = E4*(/one,n4(1),n4(2),n4(3)/)
     p5 = E5*(/one,n5(1),n5(2),n5(3)/)
     
     pV = p1 + p2
     call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,C5S5S6Lim)
    
     nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
     nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)

     eta56 = half*(one-dot_product(n5,nlept1))
     eta57 = half*(one-dot_product(n5,nlept2))
     
     call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

     C5S5S6Lim%LimMom(:,1) =   p1 
     C5S5S6Lim%LimMom(:,2) =   p2
     C5S5S6Lim%LimMom(:,3) =  C5S5S6Lim%AmpMom(:,3)
     C5S5S6Lim%LimMom(:,4) =  C5S5S6Lim%AmpMom(:,4)
     C5S5S6Lim%LimMom(:,5) =   p5
     
     s12 = two*scr(p1,p2) 
     s41 = spart *x1*(x3*x4) 
     s51 = spart *x2*x3
     s52 = spart *x2*(one-x3)
     
     z4 = x1
     
     C5S5S6Lim%Lim_KinInv(1:5) = (/s12,s51,s52,s41,z4/)
     
     C5S5S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
          * x3*two/NF1(x3,zero,lambda) & 
          * x1*x2 & 
          * weight &
          * one/two/spart & 
          * spart/MV2
     
     ! -------    TRIPLE COLLINEAR KINEMATICS  + C5 + S5 + S6  --------------
     
     eta41 = zero
     eta51 = zero
     eta42 = one
     eta52 = one
     eta54 = zero
     
     p1 = half*sqrts*(/one,zero,zero,one/)
     p2 = half*sqrts*(/one,zero,zero,-one/)
     
     E4 = half*sqrts*x1
     E5 = half*sqrts*x2

     n4 = (/zero,zero,one/)
     n5 = (/zero,zero,one/)
     
     p4 = E4*(/one,n4(1),n4(2),n4(3)/)
     p5 = E5*(/one,n5(1),n5(2),n5(3)/)
     
     eta56 = half*(one-dot_product(n5,nlept1))
     eta57 = half*(one-dot_product(n5,nlept2))

     call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)

    s41 = spart *x1*(x3 * x4)  
    s51 = spart *x2*x3
    s54 = spart *x1*x2*x3/NF1(x3,zero,lambda)

    ztot =-one 
    z4 = x1/ztot 
    z5 = x2/ztot

    TCC5S5S6Lim%Lim_KinInv(1:5) = (/s41,s51,s54,z4,z5/)
    
    TCC5S5S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1*x2 & 
         * weight &
         * one/two/spart & 
         * spart/MV2
    
    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S6 FLM   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------
    !  S6 : photon is soft
    
    phi41 = twopi*yr(4)
    
    cos41 = one - two*x3*x4
    sin41 = sqrt(one-cos41**2)
    
    cos51 = one - two*x3 
    sin51 = sqrt(one-cos51**2)
    
    if(2.0_dp*(one - x4)**2/NF1(x3,x4,lambda)-two*(one+x4*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif
    
    sinphi45 = sinphi45sign*sinphi45     ! randomize sign of sinphi45
    
    eta41 = x3*x4
    eta51 = x3
    eta42 = (one-eta41) 
    eta52 = (one - eta51)
    eta54 = x3*(one-x4)**2/NF1(x3,x4,lambda)
    
    call get_flag_nlo(x1,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts,S6Lim%flag)
    
    if (S6Lim%flag .eqv. .false.) then 

       S6Lim%PartFrac = (/xi1,xi2/)
       
       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       
       n4 = (/sin41*cos(phi41),sin41*sin(phi41),cos41/)
       n5 = (/sin51*(cos(phi41)*cosphi45-sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45),cos51/)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)
       
       pV = p1 + p2 - p4
       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,S6Lim)
       
       S6Lim%AmpMom(:,1) = p1
       S6Lim%AmpMom(:,2) = p2
       S6Lim%AmpMom(:,5) = p4

       S6Lim%LimMom(:,1) = p1
       S6Lim%LimMom(:,2) = p2
       S6Lim%LimMom(:,3) = S6Lim%AmpMom(:,3) 
       S6Lim%LimMom(:,4) = S6Lim%AmpMom(:,4) 
       S6Lim%LimMom(:,5) = p5 !photon
       
       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)

       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))
       
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)
       
       S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2
       
    endif
    
    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5  and S5 + C5  FLM   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------
    !  S5: gluon is soft

    phi41 = twopi*yr(4)
    cos41 = one - two*x3*x4
    sin41 = sqrt(one-cos41**2)
    
    cos51 = one - two*x3 
    sin51 = sqrt(one-cos51**2)
    
    if(2.0_dp*(one - x4)**2/NF1(x3,x4,lambda)-two*(one+x4*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif
    
    sinphi45 = sinphi45sign*sinphi45     ! randomize sign of sinphi45
    
    eta41 = x3*x4
    eta51 = x3
    eta42 = (one-eta41) 
    eta52 = (one - eta51)
    eta54 = x3*(one-x4)**2/NF1(x3,x4,lambda)
    
    call get_flag_nlo(x2,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,S5Lim%flag)
    
    if (S5Lim%flag .eqv. .false.) then 
       S5Lim%PartFrac = (/xi1,xi2/)
       
       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       
       n4 = (/sin41*cos(phi41),sin41*sin(phi41),cos41/)
       n5 = (/sin51*(cos(phi41)*cosphi45-sin(phi41)*sinphi45),& 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45),cos51/)
       
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pV = p1 + p2 - p5
       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,S5Lim)
       
       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,5) = p5 ! photon
       
       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = p4
       
       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)
       
       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))
       
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)
       
       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2
       
    endif
   
    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    ! --------- TRIPLE COLLINEAR -------
    
    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero
    
    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,TCLim%flag)
    
    if (TCLim%flag .eqv. .false.) then 

       TCLim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       
       n4 = (/zero,zero,one/)
       n5 = (/zero,zero,one/)
       
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)
       
       pV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,TCLim) 
       
       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)
       
       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))
       
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)
       
       s41 = spart *x1*(x3*x4) 
       s51 = spart *x2*x3
       s54 = spart *x1* x2 * x3 * (one-x4)**2/NF1(x3,x4,lambda)
       
       ztot = x1 +x2 - one 
       z4 = x1/ztot 
       z5 = x2/ztot
       
       TCLim%Lim_KinInv(1:7) = (/s41,s51,s54,z4,z5,x1,x2/)
       
       TCLim%AmpMom(:,1) = p1 - p4 - p5
       TCLim%AmpMom(:,2) = p2
       
       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2
       
       
       ! ---------   TC  + C5 ----------------
       
       s41 =  spart *x1*(x3*x4)  
       s51 =  spart *x2*x3
       s5_14 = s51 *(one-x1)
       
       zz1 = x1/(x1-one)
       zz2 = (one-x1)/(one-x1-x2)
       
       kt1kt2sq = (-one+two*lambda)**2 * (one-x3)
       
       TCC5Lim%Lim_KinInv(1:5) = (/s41,s5_14,zz1,zz2,kt1kt2sq/)
       
       TCC5Lim%wgt  = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 6
    !-------------------------------------------------------------------------------------------------------------

    ! -------------- TRIPLE COLLINEAR KINEMATICS + SOFT 6  -------------
    
    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero 
    
    call get_flag_nlo(x1,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts,TCS6Lim%flag)
    
    if (TCS6Lim%flag .eqv. .false.) then 

       TCS6Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       
       n4 = (/zero,zero,one/)
       n5 = (/zero,zero,one/)
       
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)
       
       pV = p1 + p2 - p4
       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,TCS6Lim)
       
       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)
       
       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))
       
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)
       
       s41 = spart *x1*(x3*x4) 
       s51 = spart *x2*x3
       s54 = spart *x1 * x2 * x3 * (one-x4)**2/NF1(x3,x4,lambda)
       
       z4 = x1 
       z5 = x2
       
       TCS6Lim%Lim_KinInv(1:5) = (/s41,s51,s54,z4,z5/)
       TCS6Lim%AmpMom(:,1) = p1 - p4 
       TCS6Lim%AmpMom(:,2) = p2
       
       TCS6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2
       
       ! -------   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT6 -----------
       
       s41 = spart *x1*(x3 * x4)  
       s51 = spart *x2* x3
       s54 = spart *x1*x2*x3/NF1(x3,zero,lambda)
       
       
       z4 = one - x1 
       z5 = x2
       
       TCC5S6Lim%Lim_KinInv(1:5) = (/s41,s51,s54,z4,z5/)
       
       TCC5S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2


       ! ----------  C51 + S6
       
       cos51 = one - two*x3 
       sin51 = sqrt(one-cos51**2) 
       
       if(one-two*lambda .ge. zero) then
          sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif
       
       sinphi45 = sinphi45sign*sinphi45     ! randomize sign of sinphi45
       
       eta41 = zero
       eta51 = x3
       eta42 = one 
       eta52 = (one - eta51)
       eta54 = x3/NF1(x3,zero,lambda)
       
       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       
       n4 = (/zero,zero,one/)
       n5 = (/sin51*(cos(phi41)*cosphi45-sin(phi41)*sinphi45),& 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45),cos51/)
       
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)
       
       pV = p1 + p2 - p4
       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,C5S6Lim)
       
       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))
       
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)
       
       C5S6Lim%AmpMom(:,1) = p1 - p4
       C5S6Lim%AmpMom(:,2) = p2
       
       C5S6Lim%LimMom(:,1) = p1 
       C5S6Lim%LimMom(:,2) = p2
       C5S6Lim%LimMom(:,3) = C5S6Lim%AmpMom(:,3)
       C5S6Lim%LimMom(:,4) = C5S6Lim%AmpMom(:,4) 
       C5S6Lim%LimMom(:,5) = p5 !photon
       
       z4 = one/(one - x1) 
       s12 = two*scr(p1,p2) 
       s41 = spart *x1*(x3*x4) 
       s51 = spart *x2*x3
       s52 = spart *x2*(one-x3)
       
       C5S6Lim%Lim_KinInv(1:5) = (/s12,s51,s52,s41,z4/)
       
       C5S6Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2
       
    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 4
    !-------------------------------------------------------------------------------------------------------------
    
    ! -------------- TRIPLE COLLINEAR KINEMATICS + SOFT 4  -------------
    
    eta41 = zero
    eta51 = zero
    eta42 = one
    eta52 = one
    eta54 = zero 
    
    call get_flag_nlo(x2,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts,TCS5Lim%flag)

    if (TCS5Lim%flag  .eqv. .false.) then 

       TCS5Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero,zero,one/)
       n5 = (/zero,zero,one/)
       
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)
       
       pV = p1 + p2 - p5
       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,TCS5Lim)  !! get leptons momenta
       
       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)
       
       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))
       
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)
       
       s41 = spart *x1*(x3*x4) 
       s51 = spart *x2*x3
       s54 = spart *x1 * x2 * x3 * (one-x4)**2/NF1(x3,x4,lambda)
       
       z4 = x1 
       z5 = one/(one-x2)
       
       TCS5Lim%Lim_KinInv(1:5) = (/s41,s51,s54,z4,z5/)
       TCS5Lim%AmpMom(:,1) = p1 - p5 
       TCS5Lim%AmpMom(:,2) = p2
       
       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2
       
       !---------> TCC5S5

       s41 = spart *x1*(x3*x4) 
       s51 = spart *x2*x3
       s54 = spart *x1 * x2 * x3 /NF1(x3,zero,lambda)
       
       z4 = x1 
       z5 = one/(one-x2)
       
       TCC5S5Lim%Lim_KinInv(1:5) = (/s41,s51,s54,z4,z5/)

       TCC5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2
       
    endif
    

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C5 FLM 
    !-------------------------------------------------------------------------------------------------------------

    !-------- These angular relations hold for C5 and for C5S5 below
       
    cos51 = one - two*x3 
    sin51 = sqrt(one-cos51**2) 
    
    if(one-two*lambda .ge. zero) then
       sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif
    
    sinphi45 = sinphi45sign*sinphi45     ! randomize sign of sinphi45
    
    eta41 = zero
    eta51 = x3
    eta42 = one 
    eta52 = (one - eta51)
    eta54 = x3/NF1(x3,zero,lambda)
    
    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,C5Lim%flag)
    
    if ( C5Lim%flag .eqv. .false.) then

       C5Lim%PartFrac = (/xi1,xi2/)
       
       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       
       n4 = (/zero,zero,one/)
       n5 = (/sin51*(cos(phi41)*cosphi45-sin(phi41)*sinphi45),& 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45),cos51/)
       
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)
       
       pV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,C5Lim)
       
       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)
       
       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))
       
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)
       
       s41 = spart*x1*(x3*x4)
       z4 = one/(one - x1) 
       
       nperp = (/zero,cos(phi41),sin(phi41),zero/)
       ncoll = (/zero,zero,zero,one/)
       
       C5Lim%AmpMom(:,1) = p1 - p4
       C5Lim%AmpMom(:,2) = p2
       C5Lim%AmpMom(:,5) = p5
       
       C5Lim%LimMom(:,1) = nperp
       C5Lim%LimMom(:,2) = ncoll
       
       C5Lim%Lim_KinInv(1:5) = (/zero,s41,zero,zero,z4/)
       
       C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2
       
    endif
    
    !---- C5 S5

    call get_flag_nlo(x2,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,C5S5Lim%flag)
        
    if (C5S5Lim%flag .eqv. .false.) then 

       C5S5Lim%PartFrac = (/xi1,xi2/)
       
       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       
       n4 = (/zero,zero,one /)
       n5 = (/sin51*(cos(phi41)*cosphi45-sin(phi41)*sinphi45),& 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45),cos51/)
       
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pV = p1 + p2 - p5
       call get_lept_mom(MV,pV,yr(7:9),nlept,kallenF,C5S5Lim)
       
       nlept1 = (/nlept(1,1),nlept(2,1),nlept(3,1)/)
       nlept2 = (/nlept(1,2),nlept(2,2),nlept(3,2)/)
       
       eta56 = half*(one-dot_product(n5,nlept1))
       eta57 = half*(one-dot_product(n5,nlept2))
       
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,1,weight)
       
       C5S5Lim%AmpMom(:,1) = p1
       C5S5Lim%AmpMom(:,2) = p2
       C5S5Lim%AmpMom(:,3) = C5S5Lim%AmpMom(:,3)
       C5S5Lim%AmpMom(:,4) = C5S5Lim%AmpMom(:,4)
       C5S5Lim%AmpMom(:,5) = p5
       
       s41 = spart*x1*(x3*x4)
       z4 = x1
       C5S5Lim%Lim_KinInv(1:5) = (/zero,s41,zero,zero,z4/)
       
       C5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2
       
    endif
     
  end subroutine kinematics_nnlo_5161c

  subroutine kinematics_nnlo_5262a(MV2,Y,yr,&
       HardProc,S4Lim,S5Lim,S4S5Lim,TCLim,TCS4Lim,TCS5Lim,TCS4S5Lim,&
       C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,TCC5Lim,TCC5S4Lim,TCC5S5Lim,TCC5S4S5Lim)
    real(dp), intent(in) ::  MV2,Y,yr(10)
    type(KinConfig)      ::  HardProc,S5Lim,TCLim,TCS5Lim,C5Lim,C5S5Lim,TCC5Lim,TCC5S5Lim
    type(KinConfig)      ::  S4Lim, S4S5Lim,TCS4Lim,C5S4Lim,TCC5S4Lim, TCS4S5Lim, TCC5S4S5Lim, C5S4S5Lim
    !---
    real(dp) :: p1(4), p2(4), p4(4), p5(4)
    real(dp) :: MV, kallenF, sin15, cos51, pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos41, sin41, phi41, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w1415, s45, z4, z5, ztot, weight5262, z45, ncoll(1:4), nperp(1:4)
    real(dp) :: spart,sqrts, yy, Deltas, acosphi45, asinphi45
    real(dp) :: s4_25,zz1,zz2,kt1kt2sq
    real(dp) :: sinphi45sign
    real(dp) :: n4(3),n5(3), nlept1(3),nlept2(3), nlept(3,2), zz, z, E4, E5, eta56, eta57, x5
    real(dp) :: s51, s41, s54,sin51, weight, plepHard1(4), plepHard2(4)

    real(dp) :: mysin4lim, mysin5lim, s4_15
    real(dp) :: s24, s5_24, s42, s52

    ! initialize all the properties of all kinematic configs
    ! p4 gluon
    ! p5 photon

   call initialize_config(HardProc)       !1
    call initialize_config(S4Lim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(S4S5Lim)        !4
    call initialize_config(TCLim)          !5
    call initialize_config(TCS4Lim)        !6
    call initialize_config(TCS5Lim)        !7
    call initialize_config(TCS4S5Lim)      !8
    call initialize_config(C5Lim)          !9
    call initialize_config(C5S4Lim)        !10
    call initialize_config(C5S5Lim)        !11
    call initialize_config(C5S4S5Lim)      !12
    call initialize_config(TCC5Lim)        !13
    call initialize_config(TCC5S4Lim)      !14
    call initialize_config(TCC5S5Lim)      !15
    call initialize_config(TCC5S4S5Lim)    !16

    HardProc%npart = 6
    S4Lim%npart = 5
    S5Lim%npart = 5
    S4S5Lim%npart = 4
    TCLim%npart = 4
    TCS4Lim%npart = 4
    TCS5Lim%npart = 4
    TCS4S5Lim%npart = 4
    C5Lim%npart = 5
    C5S4Lim%npart = 4
    C5S5Lim%npart = 5
    C5S4S5Lim%npart = 4
    TCC5Lim%npart = 4
    TCC5S4Lim%npart = 4
    TCC5S5Lim%npart = 4
    TCC5S4S5Lim%npart = 4

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a]
    S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]  
    S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]
    S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]
    TCLim%ids(1:4) = [0,0,id_el,-id_el]
    TCS4Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCS5Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCS4S5Lim%ids(1:4) = [0,0,id_el,-id_el] 
    C5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g] 
    C5S4Lim%ids(1:4) = [0,0,id_el,-id_el] 
    C5S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g] 
    C5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCC5Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCC5S4Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCC5S5Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCC5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el] 

    MV =  sqrt(MV2)

    x1 = yr(1) 
    x2 = yr(2) 
    x3 = yr(3) 
    x4 = yr(5)
    x5 = yr(6)


    !print*, 'sector 5262 a'
    !print*, 'x1', x1
    !print*, 'x2', x2
    !print*, 'x3', x3
    !print*, 'x4', x4

    

    !x1 = 1d-6
    !x2 = 1d-6
    !x3 = 1d-12
    !x4 = 1d-12

    lambda = sin(half*pi*x5)**2

    if (yr(10) .ge. half) then
       sinphi45sign = +one
    else
       sinphi45sign = -one
    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    cos41 = -one + two*x3
    sin41 = sqrt(one - cos41**2)
    phi41 = two*pi*yr(4)

    cos51 = -one + two*x3*x4
    sin51 = sqrt(one - cos51**2)


    if (x4.eq.zero) then 
       if(one-two*lambda .ge. zero) then
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif
    elseif(2.0_dp*(one - x4)**2/NF1(x3,x4,lambda)-two*(one+x4*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    sinphi45 = sinphi45sign*sinphi45 

    eta42 = x3
    eta41 = one - x3
    eta52 = x3*x4
    eta51 = one - x3*x4
    eta54 = x3*(one-x4)**2/NF1(x3,x4,lambda)

    
    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,HardProc%flag)

    !print*, 'HardProc%flag', HardProc%flag

    if (HardProc%flag .eqv. .false.) then 

       HardProc%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin41*cos(phi41),  sin41*sin(phi41), cos41 /)
       n5 = (/sin51*( cos(phi41)*cosphi45 - sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45 ), cos51 /)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)


       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,HardProc)

       plepHard1 = HardProc%AmpMom(:,5)
       plepHard2 = HardProc%AmpMom(:,6)

       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       !HardProc%AmpMom(:,3) = plepHard1
       !HardProc%AmpMom(:,4) = plepHard2
       HardProc%AmpMom(:,5) = p4 ! gluon 
       HardProc%AmpMom(:,6) = p5 ! photon

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

        call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &     
            * one/two/spart &
            * spart/MV2


    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: S4 S5 * FLM  : both gluon and photon soft
    !-------------------------------------------------------------------------------------------------------------

    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    S4S5Lim%PartFrac = (/xi1,xi2/)

    p1=  half*sqrts*(/one,zero,zero,one/) 
    p2 = half*sqrts*(/one,zero,zero,-one/) 
    E4 = half*sqrts*x1
    E5 = half*sqrts*x2
    n4 = (/sin41*cos(phi41),  sin41*sin(phi41), cos41 /)
    n5 = (/sin51*( cos(phi41)*cosphi45 - sin(phi41)*sinphi45), & 
         sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45 ), cos51 /)
    p4 = E4*(/one,n4(1),n4(2),n4(3)/)
    p5 = E5*(/one,n5(1),n5(2),n5(3)/)

    pVV = p1 + p2

    call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,S4S5Lim)

    S4S5Lim%AmpMom(:,1) = p1
    S4S5Lim%AmpMom(:,2) = p2
    !S4S5Lim%AmpMom(:,3)  lepton
    !S4S5Lim%AmpMom(:,4)  antilepton

    S4S5Lim%LimMom(:,1) = p1
    S4S5Lim%LimMom(:,2) = p2
    S4S5Lim%LimMom(:,3) = S4S5Lim%AmpMom(:,3)
    S4S5Lim%LimMom(:,4) = S4S5Lim%AmpMom(:,4)
    S4S5Lim%LimMom(:,5) = p4  !gluon
    S4S5Lim%LimMom(:,6) = p5  !photon

    nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
    nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

    eta56 = half*(one-sc3(n5,nlept1))
    eta57 = half*(one-sc3(n5,nlept2))

    !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

    
    S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
         * x1*x2 & 
         * weight &
         * one/two/spart &
         * spart/MV2


    ! ------   TRIPLE COLLINEAR KINEMATICS +  S4S5  -------
    ! TC = C245 S5 S4 : both gluon and photon collinear to 2,both gluon and photon soft
    ! x3 = 0, x4 = 0, x2 = 0, x1 = 0

    eta42 = zero
    eta52 = zero
    eta41 = one
    eta51 = one
    eta54 = zero

    p1=  half*sqrts*(/one,zero,zero,one/)
    p2 = half*sqrts*(/one,zero,zero,-one/)
    E4 = half*sqrts*x1
    E5 = half*sqrts*x2
    n4 = (/zero, zero, - one/)
    n5 = (/zero, zero, - one/)
    p4 = E4*(/one,n4(1),n4(2),n4(3)/)
    p5 = E5*(/one,n5(1),n5(2),n5(3)/)

    ztot = - one 
    z4 = x1/ztot 
    z5 = x2/ztot

    s42 = spart*x1*x3
    s52 = spart*x2*x3*x4
    s54 = spart*x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

    TCS4S5Lim%Lim_KinInv(1:5) = (/s42,s52,s54,z4,z5/)

    pVV = p1 + p2

    call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,TCS4S5Lim)

    nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
    nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

    eta56 = half*(one-sc3(n5,nlept1))
    eta57 = half*(one-sc3(n5,nlept2))

    !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

    TCS4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
         * x1*x2 & 
         * weight &
         * one/two/spart &
         * spart/MV2

    ! --------  C5 S4S5 FLM  ----------
    ! x4 ->0

    eta42 = x3
    eta52 = zero
    eta41 = (one-x3)
    eta51 = one
    eta54 = x3/NF1(x3,zero,lambda)

    p1 = half*sqrts*(/one,zero,zero,one/)
    p2 = half*sqrts*(/one,zero,zero,-one/)

    E4 = half*sqrts*x1
    E5 = half*sqrts*x2
    n4 = (/sin41*cos(phi41),  sin41*sin(phi41), cos41 /)
    n5 = (/zero,zero, - one/)

    p4 = E4*(/one,n4(1),n4(2),n4(3)/)
    p5 = E5*(/one,n5(1),n5(2),n5(3)/)

    ztot = - one
    z4 = x1/ztot
    z5 = x2/ztot

    s42 = spart*x1*x3
    s52 = spart*x2*x3*x4
    s54 = spart*x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

    C5S4S5Lim%Lim_KinInv(1:5) = (/s42,s52,s54,z4,z5/)
    pVV = p1 + p2 

    call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C5S4S5Lim)

    C5S4S5Lim%LimMom(:,1) = p1
    C5S4S5Lim%LimMom(:,2) = p2
    C5S4S5Lim%LimMom(:,3) = p4


    nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
    nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

    eta56 = half*(one-sc3(n5,nlept1))
    eta57 = half*(one-sc3(n5,nlept2))

    !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

    C5S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1*x2 & 
         * weight &       
         * one/two/spart &
         * spart/MV2

    ! -------   TRIPLE COLLINEAR KINEMATICS  + C5 + S4S5 --------------
    ! x3 = 0, x4 = 0, x2 = 0, x1 = 0

    eta42 = zero
    eta52 = zero
    eta41 = one
    eta51 = one
    eta54 = zero

    s42 = spart*x1*x3
    s52 = spart*x2*x3*x4
    s54 = spart*x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

    ztot =-one 
    z4 = x1/ztot 
    z5 = x2/ztot

    TCC5S4S5Lim%Lim_KinInv(1:5) = (/s42,s52,s54,z4,z5/)
    pVV = p1 + p2

    call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,TCS4S5Lim)  !! get leptons momenta

    nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
    nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

    eta56 = half*(one-sc3(n5,nlept1))
    eta57 = half*(one-sc3(n5,nlept2))

    !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

    TCC5S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1*x2 & 
         * weight &
         * one/two/spart &
         * spart/MV2


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S4    ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------

    ! --------- SOFT 4 ---------
    !  S4 : gluon is soft

    eta42 = x3
    eta52 = x3*x4
    eta41 = one - eta42
    eta51 = one - eta52
    eta54 = x3*(one-x4)**2/NF1(x3,x4,lambda)

    call get_flag_nlo(x2,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,S4Lim%flag)

    if (S4Lim%flag .eqv. .false.) then 
       S4Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin41*cos(phi41),  sin41*sin(phi41), cos41 /)
       n5 = (/sin51*( cos(phi41)*cosphi45 - sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45 ), cos51 /)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p5
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,S4Lim)

       S4Lim%AmpMom(:,1) = p1
       S4Lim%AmpMom(:,2) = p2
       !S4Lim%AmpMom(:,3) = S4Lim%AmpMom(:,4)
       !S4Lim%AmpMom(:,4) = S4Lim%AmpMom(:,5)
       S4Lim%AmpMom(:,5) = p5 ! photon

       S4Lim%LimMom(:,1) = p1
       S4Lim%LimMom(:,2) = p2
       S4Lim%LimMom(:,3) = p4

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       S4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

    endif

    ! --------- C5 S4 ---------
    ! C52 S4 : photon collinear, gluon soft
    ! x4 = 0, x1 -> 0


    eta42 = x3
    eta52 = zero
    eta41 = one - x3
    eta51 = one 
    eta54 = x3/NF1(x3,zero,lambda)

    call get_flag_nlo(x2,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,C5S4Lim%flag)

    if (C5S4Lim%flag .eqv. .false.) then 
       C5S4Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/)
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin41*cos(phi41),  sin41*sin(phi41), cos41 /)
       n5 = (/zero,zero,-one/)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p5

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C5S4Lim)

       C5S4Lim%AmpMom(:,1) = p1
       C5S4Lim%AmpMom(:,2) = p2-p5

       C5S4Lim%LimMom(:,1) =  p1
       C5S4Lim%LimMom(:,2) =  p2
       C5S4Lim%LimMom(:,3) =  p4

       s52 = spart*x2*(x3*x4)
       z5 = one/(one - x2) 

       C5S4Lim%Lim_KinInv(1:5) = (/zero,s52,zero,zero,z5/)

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       C5S4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

    endif



    ! --------- TC C5 S4 ---------
    ! x3 = 0, x4 = 0, x1 = 0, x4 << x3

    eta42 = zero
    eta52 = zero
    eta41 = one 
    eta51 = one 
    eta54 = x3/NF1(zero,zero,lambda)

    call get_flag_nlo(x2,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,TCC5S4Lim%flag)

    !call get_flag_nlo(x2,MV2,Y,x2,eta52,eta51,yy,xi1,xi2,spart,sqrts,TCC5S4Lim%flag)
   

    if (TCC5S4Lim%flag .eqv. .false.) then 
       TCC5S4Lim%PartFrac = (/xi1,xi2/)
       
       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero,zero, -one/)
       n5 = (/zero,zero, -one/)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p5

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,TCC5S4Lim)

       TCC5S4Lim%AmpMom(:,1) = p1 
       TCC5S4Lim%AmpMom(:,2) = p2 - p5

       TCC5S4Lim%LimMom(:,1) =  p1
       TCC5S4Lim%LimMom(:,2) =  p2
       TCC5S4Lim%LimMom(:,3) =  p4

       z5 = one/(one - x2)
       s42 = spart*x1*x3
       s52 = spart*x2*x3*x4

       TCC5S4Lim%Lim_KinInv(1:5) = (/s42,s52,zero,z4,z5/)

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       TCC5S4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 + C5 S5   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------

    ! --------- SOFT 5 ---------
    !  S5 : photon is soft


    eta42 = x3
    eta52 = x3*x4
    eta41 = one- x3
    eta51 = one -  x3*x4
    eta54 = x3*(one-x4)**2/NF1(x3,x4,lambda)

    call get_flag_nlo(x1,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts,S5Lim%flag)

    if (S5Lim%flag .eqv. .false.) then 
       S5Lim%PartFrac = (/xi1,xi2/)


       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin41*cos(phi41), sin41*sin(phi41), cos41/)
       n5 = (/sin51*( cos(phi41)*cosphi45 - sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45 ), cos51 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p4

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       !S5Lim%AmpMom(:,3) =  S5Lim%AmpMom(:,4) ! lepton
       !S5Lim%AmpMom(:,4) =  S5Lim%AmpMom(:,5) ! antilepton
       S5Lim%AmpMom(:,5) = p4

       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = S5Lim%AmpMom(:,3) 
       S5Lim%LimMom(:,4) = S5Lim%AmpMom(:,4) 
       S5Lim%LimMom(:,5) = p5 !photon

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

      
       ! kin inv and weight for C5 +  S5
       ! C52 S5 : photon collinear, photon soft
       ! x4 = 0, x2 -> 0

       eta42 = x3
       eta52 = zero
       eta41 = one- eta42
       eta51 = one - eta52
       eta54 = x3/NF1(x3,zero,lambda)

       n5 = (/zero,zero,-one/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       s52 = spart*x2*(x3*x4)
       z5 = x2
       C5S5Lim%Lim_KinInv(1:5) =  (/zero,s52,zero,zero,z5/)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       C5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    ! --------- TRIPLE COLLINEAR -------

    eta42 = zero
    eta52 = zero
    eta41 = one 
    eta51 = one 
    eta54 = zero

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,TCLim%flag)

    if (TCLim%flag .eqv. .false.) then 
       TCLim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero,zero,-one/)
       n5 = (/zero,zero,-one/)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       z5 = x2/(x1 + x2 - one)
       z4 = x1/(x1 + x2 - one)
       z = one/(one - x1- x2)

       s42 = spart*x1*x3
       s52 = spart*x2*x3*x4
       s54 = spart*x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

       TCLim%Lim_KinInv(1:7) = (/s42,s52,s54,z4,z5,x1,x2/)

       pVV = p1 + p2 - p4 - p5 

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,TCLim) 

       TCLim%AmpMom(:,1) = p1 - p4 - p5
       TCLim%AmpMom(:,2) = p2

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))
       
       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)


       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2



       ! --------- TRIPLE COLLINEAR + C5  -------
       !  TC = C245 C25 : both gluon and photon collinear to 1, photon super collinear
       ! x3 = 0, x4 = 0, x4 << x3

       s42 = spart * x1*x3 
       s52 = spart * x2*(x3*x4)
       s4_25 = s42 * (one-x2)

       zz1 = x2/(x2-one)
       zz2 = (one-x2)/(one-x1-x2)

       kt1kt2sq = (-one+two*lambda)**2 * (one-x3)

       TCC5Lim%Lim_KinInv(1:5) = (/s52,s4_25,zz1,zz2,kt1kt2sq/)

       TCC5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 4
    !-------------------------------------------------------------------------------------------------------------

    ! -------- TRIPLE COLLINEAR + SOFT 4

    ! TC = C245 S4 : both gluon and photon collinear to 2, gluon soft
    ! x3 = 0, x4 = 0, x1 = 0

    eta42 = zero
    eta52 = zero
    eta41 = one 
    eta51 = one 
    eta54 = zero

    call get_flag_nlo(x2,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,TCS4Lim%flag)

    if (TCS4Lim%flag .eqv. .false.) then 
       TCS4Lim%PartFrac = (/xi1,xi2/)

       p1=  half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero, zero, -one/)
       n5 = (/zero, zero, -one/)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       z4 = x1 
       z5 = one/(one-x2)

       s42 = spart*x1*x3
       s52 = spart*x2*x3*x4
       s54 = spart*x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

       TCS4Lim%Lim_KinInv(1:5) = (/s42,s52,s54,z4,z5/)

       pVV = p1 + p2 - p5 

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,TCS4Lim) 

       TCS4Lim%AmpMom(:,1) = p1
       TCS4Lim%AmpMom(:,2) = p2 - p5

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       TCS4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    ! -------- TRIPLE COLLINEAR + SOFT 5
    ! TC = C245 S5 : both gluon and photon collinear to 2, photon soft
    ! x3 = 0, x4 = 0, x2 = 0

    call get_flag_nlo(x1,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts,TCS5Lim%flag)

    if (TCS5Lim%flag .eqv. .false.) then 
       TCS5Lim%PartFrac = (/xi1,xi2/)


       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero,zero,-one/)
       n5 = (/zero,zero,-one/)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)
       
       z4 = x1 
       z5 = x2

       s42 = spart*x1*x3
       s52 = spart*x2*x3*x4
       s54 = spart*x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

       TCS5Lim%Lim_KinInv(1:5) = (/s42,s52,s54,z4,z5/)

       pVV = p1 + p2 - p4

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,TCS5Lim) 

       TCS5Lim%AmpMom(:,1) = p1 
       TCS5Lim%AmpMom(:,2) = p2 - p4

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

!    print*, 'nlept1', nlept1
!    print*, 'nlept2', nlept2

!print*, 'weight', weight 
!print*, 'eta56', eta56
!print*, 'eta57', eta57
!print*, 'n5', n5
!print*, 'eta41', eta41
!print*, 'eta42', eta42
!print*, 'eta51', eta51
!print*, 'eta52', eta52
!print*, 'nlept1', nlept1
!print*, 'nlept2', nlept2

       ! -------   TRIPLE COLLINEAR KINEMATICS  + C5 + SOFT5 ---------

       z4 = one/(one-x1) 
       z5 = x2

       s42 = spart *x1*x3 
       s52 = spart *x2*(x3*x4)
       s45 = spart *x1*x2*x3/NF1(x3,zero,lambda)

       TCC5S5Lim%Lim_KinInv(1:5) = (/s42,s52,s45,z4,z5/)

       TCC5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C5 FLM 
    !-------------------------------------------------------------------------------------------------------------


    phi41 = twopi*yr(4) 
    cos41 = - one + two*x3 
    sin41 = sqrt(one-cos41**2)

    eta42 = x3
    eta52 = zero
    eta41 = (one-x3) 
    eta51 = one
    eta54 = x3/NF1(x3,zero,lambda)

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,C5Lim%flag)
    
    if (C5Lim%flag .eqv. .false.) then 

       C5Lim%PartFrac = (/xi1,xi2/)
       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin41*cos(phi41), sin41*sin(phi41), cos41 /)
       n5 = (/zero, zero, -one /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/) 
    
       s52 = spart*x2*(x3*x4)
       z5 = one/(one - x2) 

       C5Lim%Lim_KinInv(1:5) = (/zero,s52,zero,zero,z5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C5Lim)

       C5Lim%AmpMom(:,1) = p1 
       C5Lim%AmpMom(:,2) = p2- p5
       !C5Lim%AmpMom(:,3) = C5Lim%AmpMom(:,4)
       !C5Lim%AmpMom(:,4) = C5Lim%AmpMom(:,5)
       C5Lim%AmpMom(:,5) = p4

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))


       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2


    endif


  end subroutine kinematics_nnlo_5262a

    
  subroutine kinematics_nnlo_5262c(MV2,Y,yr, &
       HardProc,S4Lim,S5Lim,S4S5Lim,TCLim,TCS4Lim,TCS5Lim,TCS4S5Lim,&
       C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,TCC4Lim,TCC4S4Lim,TCC4S5Lim,TCC4S4S5Lim)
    real(dp), intent(in) ::  MV2, Y, yr(10)
    type(KinConfig) ::  HardProc,S4Lim,S4S5Lim,S5Lim,TCLim,TCS4Lim,TCS4S5Lim,TCS5Lim, &
         C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,TCC4Lim,TCC4S4Lim,TCC4S5Lim,TCC4S4S5Lim
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4) 
    real(dp) :: E6,E7, MV, p67, kallenF, sin51, cos51,  pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos41, sin41, phi41, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4151, s14, s15, s45, z4, z5, ztot
    real(dp) :: spart,sqrts, s12, s25, z45, yy, Deltas 
    real(dp) :: nperp(4), ncoll(4), mycphi, z1, z2, s14_5
    real(dp) :: s5_14,zz1,zz2,kt1kt2sq
    real(dp) :: mysign
    real(dp) :: weight,nlept(3,2),nlept1(3), nlept2(3), eta56, eta57, n5(3), E4, E5, n4(3)
    real(dp) :: s24, s5_24, s42, s52

    call initialize_config(HardProc)       !1
    call initialize_config(S4Lim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(S4S5Lim)        !4
    call initialize_config(TCLim)          !5
    call initialize_config(TCS4Lim)        !6
    call initialize_config(TCS5Lim)        !7
    call initialize_config(TCS4S5Lim)      !8
    call initialize_config(C4Lim)          !9
    call initialize_config(C4S4Lim)        !10
    call initialize_config(C4S5Lim)        !11
    call initialize_config(C4S4S5Lim)      !12
    call initialize_config(TCC4Lim)        !13
    call initialize_config(TCC4S4Lim)      !14
    call initialize_config(TCC4S5Lim)      !15
    call initialize_config(TCC4S4S5Lim)    !16

    HardProc%npart = 6
    S4Lim%npart = 5
    S5Lim%npart = 5
    S4S5Lim%npart = 4
    TCLim%npart = 4
    TCS4Lim%npart = 4
    TCS5Lim%npart = 4
    TCS4S5Lim%npart = 4
    C4Lim%npart = 5
    C4S4Lim%npart = 5
    C4S5Lim%npart = 4
    C4S4S5Lim%npart = 4
    TCC4Lim%npart = 4 
    TCC4S4Lim%npart = 4
    TCC4S5Lim%npart = 4
    TCC4S4S5Lim%npart = 4

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a] 
    S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a] 
    S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g] 
    S4S5Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCLim%ids(1:4) = [0,0,id_el,-id_el] 
    TCS4Lim%ids(1:4) = [0,0,id_el,-id_el]
    TCS5Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCS4S5Lim%ids(1:4) = [0,0,id_el,-id_el]
    C4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a] 
    C4S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]
    C4S5Lim%ids(1:4) = [0,0,id_el,-id_el] 
    C4S4S5Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCC4Lim%ids(1:4) = [0,0,id_el,-id_el]  
    TCC4S4Lim%ids(1:4) = [0,0,id_el,-id_el]
    TCC4S5Lim%ids(1:4) = [0,0,id_el,-id_el] 
    TCC4S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]

    
    MV =  sqrt(MV2)

    x1 = yr(1) 
    x2 = yr(2) 
    x3 = yr(3) 
    x4 = yr(5) 

    !x1 = 1d-6
    !x2 = 1d-6
    !x3 = 1d-12
    !x4 = 1d-12

    if (yr(10) .ge. half) then
       mysign = +one
    else
       mysign = -one
    endif

    lambda = sin( half * pi * yr(6) )**2

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi41 = twopi*yr(4)

    cos41 = - one + two*x3*x4
    sin41 = sqrt(one-cos41**2)

    cos51 = - one + two*x3 
    sin51 = sqrt(one-cos51**2)

    if (x4.eq.zero) then 
       if(one-two*lambda .ge. zero) then
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif
    elseif(2.0_dp*(one - x4)**2/NF1(x3,x4,lambda)-two*(one+x4*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    sinphi45 = mysign*sinphi45     ! randomize sign of sinphi45

    eta42 = x3*x4
    eta52 = x3
    eta41 = one - x3*x4
    eta51 = one -  x3
    eta54 = x3*(one-x4)**2/NF1(x3,x4,lambda)

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,HardProc%flag)

    if (HardProc%flag .eqv. .false.) then 
       HardProc%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin41*cos(phi41),  sin41*sin(phi41), cos41 /)
       n5 = (/sin51*( cos(phi41)*cosphi45 - sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45 ), cos51 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)


       pVV = p1 + p2 - p4 - p5       
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,HardProc)

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       !reshuffle momenta for the amplitude
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       !HardProc%AmpMom(:,3) = HardProc%AmpMom(:,5) !lepton 1
       !HardProc%AmpMom(:,4) = HardProc%AmpMom(:,6) !lepton 2
       HardProc%AmpMom(:,5) = p4 ! gluon 
       HardProc%AmpMom(:,6) = p5 ! photon

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2
       
    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: S4S5 FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! -------  S4 S5 ------------

    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    S4S5Lim%PartFrac = (/xi1,xi2/)

    p1 = half*sqrts*(/one,zero,zero,one/) 
    p2 = half*sqrts*(/one,zero,zero,-one/)
    E4 = half*sqrts*x1
    E5 = half*sqrts*x2
    n4 = (/sin41*cos(phi41),  sin41*sin(phi41), cos41 /)
    n5 = (/sin51*( cos(phi41)*cosphi45 - sin(phi41)*sinphi45), & 
         sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45 ), cos51 /)
    p4 = E4*(/one,n4(1),n4(2),n4(3)/)
    p5 = E5*(/one,n5(1),n5(2),n5(3)/)

    pVV = p1 + p2

    call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,S4S5Lim)

    nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
    nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

    eta56 = half*(one-sc3(n5,nlept1))
    eta57 = half*(one-sc3(n5,nlept2))

    !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

    S4S5Lim%AmpMom(:,1) = p1
    S4S5Lim%AmpMom(:,2) = p2
    !S4S5Lim%AmpMom(:,3) = plept1
    !S4S5Lim%AmpMom(:,4) = plept2

    S4S5Lim%LimMom(:,1) = p1
    S4S5Lim%LimMom(:,2) = p2
    S4S5Lim%LimMom(:,3) = S4S5Lim%AmpMom(:,3)
    S4S5Lim%LimMom(:,4) = S4S5Lim%AmpMom(:,4)
    S4S5Lim%LimMom(:,5) = p4 
    S4S5Lim%LimMom(:,6) = p5

    S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
         * x1*x2 & 
         * weight &
         * one/two/spart & 
         * spart/MV2

    ! -----   TRIPLE COLLINEAR KINEMATICS +  S4S5   ----------

    eta42 = zero
    eta52 = zero
    eta41 = one
    eta51 = one
    eta54 = zero

    p1 = half*sqrts*(/one,zero,zero,one/) 
    p2 = half*sqrts*(/one,zero,zero,-one/)
    E4 = half*sqrts*x1
    E5 = half*sqrts*x2
    n4 = (/zero, zero, -one /)
    n5 = (/zero, zero, -one /)
    p4 = E4*(/one,n4(1),n4(2),n4(3)/)
    p5 = E5*(/one,n5(1),n5(2),n5(3)/)

    pVV = p1 + p2

    call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,TCS4S5Lim)

    nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
    nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

    eta56 = half*(one-sc3(n5,nlept1))
    eta57 = half*(one-sc3(n5,nlept2))

    !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

    s24 = spart * x1*(x3*x4) 
    s25 = spart * x2*x3
    s45 = spart * x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

    ztot = - one 
    z4 = x1/ztot 
    z5 = x2/ztot

    TCS4S5Lim%Lim_KinInv(1:5) = (/s24,s25,s45,z4,z5/)

    TCS4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
         * x1*x2 & 
         * weight &
         * one/two/spart &
         * spart/MV2

    ! -------  C42 S4S5 FLM  ----------------

    !momenta are not stored:
    !for the amplitude we will use
    !S4S5Lim%AmpMom(:,1) = p1
    !S4S5Lim%AmpMom(:,2) = p2
    !S4S5Lim%AmpMom(:,3) = plept1
    !S4S5Lim%AmpMom(:,4) = plept2
    !for the soft photon kernel we will use
    !S4S5Lim%LimMom(:,1) = p1
    !S4S5Lim%LimMom(:,2) = p2
    !S4S5Lim%LimMom(:,3) = S4S5Lim%AmpMom(:,3)
    !S4S5Lim%LimMom(:,4) = S4S5Lim%AmpMom(:,4)
    !S4S5Lim%LimMom(:,5) = p4  -> we do not need the gluon momentum for the photon kernel
    !S4S5Lim%LimMom(:,6) = p5

    !from here

    phi41 = twopi*yr(4)
    
    cos41 = - one 
    sin41 = zero

    cos51 = -one + two*x3 
    sin51 = sqrt(one-cos51**2)
    
    eta52 = x3
    eta42 = zero
    eta51 = (one-x3)
    eta41 = one
    eta54 = x3/NF1(x3,zero,lambda)

    if(one-two*lambda .ge. zero) then
       sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    sinphi45 = mysign*sinphi45     ! randomize sign of sinphi45
    
       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero, zero, - one/)
       n5 = (/sin51*( cos(phi41)*cosphi45 - sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45 ), cos51 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

    !to here
    

    ztot = - one
    z4 = x1/ztot
    z5 = x2/ztot
    
    pVV = p1 + p2 

    call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C4S4S5Lim)
   
   
    C4S4S5Lim%LimMom(:,1) =   p1 
    C4S4S5Lim%LimMom(:,2) =   p2
    C4S4S5Lim%LimMom(:,3) =  C4S4S5Lim%AmpMom(:,3)
    C4S4S5Lim%LimMom(:,4) =  C4S4S5Lim%AmpMom(:,4)
    C4S4S5Lim%LimMom(:,5) =   p5

    nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
    nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

    eta56 = half*(one-sc3(n5,nlept1))
    eta57 = half*(one-sc3(n5,nlept2))

    !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

    s12 = two*scr(p1,p2) 
    s24 = spart *x1*(x3*x4) 
    s25 = spart *x2*x3
    s15 = spart *x2*(one-x3)

    z4 = x1

    C4S4S5Lim%Lim_KinInv(1:5) = (/s12,s25,s15,s24,z4/)

    C4S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1*x2 & 
         * weight &
         * one/two/spart & 
         * spart/MV2


!!$    eta42 = zero
!!$    eta52 = x3
!!$    eta41 = one 
!!$    eta51 = one - x3
!!$    eta54 = x3/NF1(x3,zero,lambda)
!!$
!!$    p1 = half*sqrts*(/one,zero,zero,one/) 
!!$    p2 = half*sqrts*(/one,zero,zero,-one/)
!!$    E4 = half*sqrts*x1
!!$    E5 = half*sqrts*x2
!!$    n4 = (/zero, zero, -one/)
!!$    n5 = (/sin51*( cos(phi41)*cosphi45 - sin(phi41)*sinphi45), & 
!!$         sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45 ), cos51 /)
!!$    p4 = E4*(/one,n4(1),n4(2),n4(3)/)
!!$    p5 = E5*(/one,n5(1),n5(2),n5(3)/)
!!$
!!$    eta56 = half*(one-sc3(n5,nlept1))
!!$    eta57 = half*(one-sc3(n5,nlept2))
!!$
!!$    !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
!!$    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)
!!$
!!$    s12 = two*scr(p1,p2) 
!!$    s24 = spart *x1*(x3*x4) 
!!$    s25 = spart *x2*x3
!!$    s15 = spart *x2*(one-x3)
!!$
!!$    z4 = x1
!!$
!!$    C4S4S5Lim%Lim_KinInv = (/s12,s25,s15,s24,z4/)
!!$
!!$    C4S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
!!$         * x3*two/NF1(x3,zero,lambda) & 
!!$         * x1*x2 & 
!!$         * weight &
!!$         * one/two/spart & 
!!$         * spart/MV2

    ! -------    TRIPLE COLLINEAR KINEMATICS  + C4 + S4 + S5  --------------

    !momenta are not stored:
    !for the amplitude we will use
    !S4S5Lim%AmpMom(:,1) = p1
    !S4S5Lim%AmpMom(:,2) = p2
    !S4S5Lim%AmpMom(:,3) = plept1
    !S4S5Lim%AmpMom(:,4) = plept2

    eta42 = zero
    eta52 = zero
    eta41 = one 
    eta51 = one 
    eta54 = zero

    p1 = half*sqrts*(/one,zero,zero,one/)
    p2 = half*sqrts*(/one,zero,zero,-one/)

    E4 = half*sqrts*x1
    E5 = half*sqrts*x2
    n4 = (/zero,zero,-one/)
    n5 = (/zero,zero,-one/)

    p4 = E4*(/one,n4(1),n4(2),n4(3)/)
    p5 = E5*(/one,n5(1),n5(2),n5(3)/)

    eta56 = half*(one-sc3(n5,nlept1))
    eta57 = half*(one-sc3(n5,nlept2))

    !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)


    s12 = spart 
    
    s24 = spart * x1 * ( x3 * x4)  
    s25 = spart * x2 * x3
    s45 = spart * x1 * x2 * x3/NF1(x3,zero,lambda)

    z4 = -x1
    s15 = spart *x2

    TCC4S4S5Lim%Lim_KinInv(1:5) = (/s12,s24,s25,s15,z4/)


    TCC4S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
         * x3*two/NF1(x3,zero,lambda) & 
         * x1*x2 & 
         * weight &
         * one/two/spart & 
         * spart/MV2

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S5 FLM   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------
    ! --------- SOFT 5 ---------
    !  S5 : photon is soft


    phi41 = twopi*yr(4)

    cos41 = - one + two*x3*x4
    sin41 = sqrt(one-cos41**2)

    cos51 = - one + two*x3 
    sin51 = sqrt(one-cos51**2)

    if (x4.eq.zero) then 
       if(one-two*lambda .ge. zero) then
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
          cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
       else
          sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
          cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
       endif
    elseif(2.0_dp*(one - x4)**2/NF1(x3,x4,lambda)-two*(one+x4*(one-two*x3)).gt.0) then
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))*(one-x4)/NF1(x3,x4,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    sinphi45 = mysign*sinphi45     ! randomize sign of sinphi45


    eta42 = x3*x4
    eta52 = x3
    eta41 = one-x3*x4
    eta51 = one - x3
    eta54 = x3*(one-x4)**2/NF1(x3,x4,lambda)

    call get_flag_nlo(x1,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts,S5Lim%flag)

    if (S5Lim%flag .eqv. .false.) then 
       S5Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin41*cos(phi41), sin41*sin(phi41), cos41/)
       n5 = (/sin51*( cos(phi41)*cosphi45 - sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45 ), cos51 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p4

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       !S5Lim%AmpMom(:,3) = S5Lim%AmpMom(:,4) ! lepton
       !S5Lim%AmpMom(:,4) = S5Lim%AmpMom(:,5) ! antilepton
       S5Lim%AmpMom(:,5) = p4

       S5Lim%LimMom(:,1) = p1
       S5Lim%LimMom(:,2) = p2
       S5Lim%LimMom(:,3) = S5Lim%AmpMom(:,3) 
       S5Lim%LimMom(:,4) = S5Lim%AmpMom(:,4) 
       S5Lim%LimMom(:,5) = p5 !photon

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 3: S4  and S4 + C4  FLM   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------


    ! --------- SOFT 4 ---------
    !  S4 : gluon is soft

    eta42 = x3*x4
    eta52 = x3
    eta41 = one - x3*x4
    eta51 = one -  x3
    eta54 = x3*(one-x4)**2/NF1(x3,x4,lambda)

    call get_flag_nlo(x2,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,S4Lim%flag)

    if (S4Lim%flag .eqv. .false.) then 
       S4Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin41*cos(phi41),  sin41*sin(phi41), cos41 /)
       n5 = (/sin51*( cos(phi41)*cosphi45 - sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45 ), cos51 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p5

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,S4Lim)

       S4Lim%AmpMom(:,1) = p1
       S4Lim%AmpMom(:,2) = p2
       !S4Lim%AmpMom(:,3) = S4Lim%AmpMom(:,4)
       !S4Lim%AmpMom(:,4) = S4Lim%AmpMom(:,5)
       S4Lim%AmpMom(:,5) = p5 ! photon

       S4Lim%LimMom(:,1) = p1
       S4Lim%LimMom(:,2) = p2
       S4Lim%LimMom(:,3) = p4

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       S4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 4: TRIPLE COLLINEAR KINEMATICS sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    ! --------- TRIPLE COLLINEAR -------

    eta42 = zero
    eta52 = zero
    eta41 = one
    eta51 = one
    eta54 = zero

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,TCLim%flag)

    if (TCLim%flag .eqv. .false.) then 
       TCLim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero,zero,-one/)
       n5 = (/zero,zero,-one/)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,TCLim) 

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       s24 = spart *x1*(x3*x4) 
       s25 = spart *x2*x3
       s45 = spart *x1* x2 * x3 * (one-x4)**2/NF1(x3,x4,lambda)

       ztot = x1 +x2 - one 
       z4 = x1/ztot 
       z5 = x2/ztot

       TCLim%Lim_KinInv(1:7) = (/s24,s25, s45, z4, z5,x1,x2/)

       TCLim%AmpMom(:,1) = p1 - p4 - p5
       TCLim%AmpMom(:,2) = p2

       TCLim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2


       ! ---------   TC  + C4 ----------------

       s24 =  spart *x1*(x3*x4)  
       s25 =  spart *x2*x3
       s5_24 = s25 *(one-x1)

       zz1 = x1/(x1-one)
       zz2 = (one-x1)/(one-x1-x2)

       kt1kt2sq = (-one+two*lambda)**2 * (one-x3)

       TCC4Lim%Lim_KinInv(1:5) = (/s24,s5_24,zz1,zz2,kt1kt2sq/)

       TCC4Lim%wgt  = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    ! -------------- TRIPLE COLLINEAR KINEMATICS + SOFT 5  -------------

    eta42 = zero
    eta52 = zero
    eta41 = one 
    eta51 = one 
    eta54 = zero 

    call get_flag_nlo(x1,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts,TCS5Lim%flag)

    if (TCS5Lim%flag .eqv. .false.) then 
       TCS5Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero,zero,-one/)
       n5 = (/zero,zero,-one/)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,TCS5Lim)

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       s24 = spart *x1*(x3*x4) 
       s25 = spart *x2*x3
       s45 = spart *x1*x2*x3*(one-x4)**2/NF1(x3,x4,lambda)

       z4 = x1 
       z5 = x2

       TCS5Lim%Lim_KinInv(1:5) = (/s24,s25,s45,z4,z5/)
       TCS5Lim%AmpMom(:,1) = p1 - p4 
       TCS5Lim%AmpMom(:,2) = p2

       TCS5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2


       ! -------   TRIPLE COLLINEAR KINEMATICS  + C4 + SOFT5 -----------

       s24 = spart *x1*(x3 * x4)  
       s25 = spart *x2* x3
       s45 = spart *x1*x2*x3/NF1(x3,zero,lambda)

       z4 = one - x1 
       z5 = x2

       TCC4S5Lim%Lim_KinInv(1:5) = (/s24,s25,s45,z4,z5/)

       TCC4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2


       ! ----------  C42 + S5

       phi41 = twopi*yr(4)

       cos41 = - one 
       sin41 = zero

       cos51 = - one + two*x3 
       sin51 = sqrt(one-cos51**2)

       
    if(one-two*lambda .ge. zero) then
       sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    sinphi45 = mysign*sinphi45     ! randomize sign of sinphi45

       eta42 = zero
       eta52 = x3
       eta41 = one 
       eta51 = one - x3
       eta54 = x3/NF1(x3,zero,lambda)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero, zero, -one/)
       n5 = (/sin51*( cos(phi41)*cosphi45 - sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45 ), cos51 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)
      

       pVV = p1 + p2 - p4

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C4S5Lim)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)
       
       
       C4S5Lim%AmpMom(:,1) = p1 
       C4S5Lim%AmpMom(:,2) = p2 - p4
       !C4S5Lim%AmpMom(:,3) = C4S5Lim%AmpMom(:,4) ! lepton
       !C4S5Lim%AmpMom(:,4) = C4S5Lim%AmpMom(:,5) ! antilepton

       C4S5Lim%LimMom(:,1) = p1
       C4S5Lim%LimMom(:,2) = p2 - p4
       C4S5Lim%LimMom(:,3) = C4S5Lim%AmpMom(:,3)
       C4S5Lim%LimMom(:,4) = C4S5Lim%AmpMom(:,4) 
       C4S5Lim%LimMom(:,5) = p5 !photon

       z4 = one/(one - x1)  
       s12 = two*scr(p1,p2) 
       s24 = spart *x1*(x3*x4) 
       s25 = spart *x2*x3
       s15 = spart *x2*(one-x3)

       C4S5Lim%Lim_KinInv(1:5) = (/s12,s25,s15,s24,z4/)
       
       C4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2
       
    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: TRIPLE COLLINEAR KINEMATICS + SOFT 4
    !-------------------------------------------------------------------------------------------------------------

    ! -------------- TRIPLE COLLINEAR KINEMATICS + SOFT 4  -------------

    eta42 = zero
    eta52 = zero
    eta41 = one
    eta51 = one
    eta54 = zero 

    call get_flag_nlo(x2,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts,TCS4Lim%flag)

    if (TCS4Lim%flag  .eqv. .false.) then 
       TCS4Lim%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero,zero,-one/)
       n5 = (/zero,zero,-one/)

       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p5
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,TCS4Lim)  !! get leptons momenta

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       s24 = spart *x1*(x3*x4) 
       s25 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 * (one-x4)**2/NF1(x3,x4,lambda)

       z4 = x1 
       z5 = one/(one-x2)

       TCS4Lim%Lim_KinInv(1:5) = (/s24,s25,s45,z4,z5/)
       TCS4Lim%AmpMom(:,1) = p1 - p5 
       TCS4Lim%AmpMom(:,2) = p2

       TCS4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*(one - x4)*two/NF1(x3,x4,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

       !---------> TCC4S4

       s24 = spart *x1*(x3*x4) 
       s25 = spart *x2*x3
       s45 = spart *x1 * x2 * x3 /NF1(x3,zero,lambda)

       z4 = x1 
       z5 = one/(one-x2)

       TCC4S4Lim%Lim_KinInv(1:5) = (/s24,s25,s45,z4,z5/)

       TCC4S4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C4 FLM 
    !-------------------------------------------------------------------------------------------------------------

    !-------- These angular relations hold for C4 and for C4S4 below


    phi41 = twopi*yr(4)

       cos41 = - one 
       sin41 = zero

       cos51 = - one + two*x3 
       sin51 = sqrt(one-cos51**2) 

    if(one-two*lambda .ge. zero) then
       sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
       cosphi45 = (-1.0_dp)*sqrt(abs(1.0_dp-sinphi45**2)) 
    else
       sinphi45= two*sqrt(lambda*(one-lambda))/NF1(x3,zero,lambda)
       cosphi45 = sqrt(abs(1.0_dp-sinphi45**2)) 
    endif

    sinphi45 = mysign*sinphi45     ! randomize sign of sinphi45

    eta42 = zero
    eta52 = x3
    eta41 = one 
    eta51 = one - x3
    eta54 = x3/NF1(x3,zero,lambda)

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,C4Lim%flag)

    if ( C4Lim%flag .eqv. .false.) then
       C4Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)

       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero, zero, -one /)
       n5 = (/sin51*( cos(phi41)*cosphi45 - sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45 ), cos51 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C4Lim)

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       s24 = spart*x1*(x3*x4)
       z4 = one/(one - x1 ) 

       nperp = (/zero,cos(phi41),sin(phi41),zero/)
       ncoll = (/zero,zero,zero,-one/)

       C4Lim%AmpMom(:,1) = p1 
       C4Lim%AmpMom(:,2) = p2 - p4
       !C4Lim%AmpMom(:,3) = C4Lim%AmpMom(:,4)
       !C4Lim%AmpMom(:,4) = C4Lim%AmpMom(:,5)
       C4Lim%AmpMom(:,5) = p5

       C4Lim%LimMom(:,1) = nperp
       C4Lim%LimMom(:,2) = ncoll

       C4Lim%Lim_KinInv(1:5) = (/zero,s24,zero,zero,z4/)

       C4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2

     

    endif
    
    !---- C4 S4

    call get_flag_nlo(x2,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,C4S4Lim%flag)

    if (C4S4Lim%flag .eqv. .false.) then 
       C4S4Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/zero, zero, -one /)
       n5 = (/sin51*( cos(phi41)*cosphi45 - sin(phi41)*sinphi45), & 
            sin51*(sin(phi41)*cosphi45+cos(phi41)*sinphi45 ), cos51 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p5

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C4S4Lim)

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,2, weight)

       C4S4Lim%AmpMom(:,1) = p1
       C4S4Lim%AmpMom(:,2) = p2
       !C4S4Lim%AmpMom(:,3) = C4S4Lim%AmpMom(:,4)
       !C4S4Lim%AmpMom(:,4) = C4S4Lim%AmpMom(:,5)
       C4S4Lim%AmpMom(:,5) = p5

       s24 = spart*x1*(x3*x4)
       z4 = x1
       C4S4Lim%Lim_KinInv(1:5) = (/zero,s24,zero,zero,z4/)

       C4S4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/four)**2 & 
            * x3*two/NF1(x3,zero,lambda) & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

!!$       print*, 'C4S4'
!!$       print*, 'weight', weight
!!$         print*, 'p1', p1
!!$             print*, 'p2', p2
!!$             print*, 'p4', p4
!!$             print*, 'p5', p5
!!$
!!$             print*, 'n5', n5
!!$       print*, 'pVV', pvv

    endif


  end subroutine kinematics_nnlo_5262c


  
  subroutine kinematics_nnlo_5162(MV2,Y,yr,&
       HardProc,S4Lim,S5Lim,S4S5Lim,C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,&
       C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim)
    real(dp), intent(in) ::  yr(10), MV2, Y
    type(KinConfig)      :: HardProc,S4Lim,S4S5Lim,S5Lim,C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,&
         C5Lim,C5S4Lim,C5S4S5Lim,C5S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos4, sin4, phi4, phi5, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4151, s14, s15, s45, z4, z5, ztot, weight5261
    real(dp) :: spart,sqrts, s25, s12, s24, yy, Deltas
    real(dp) :: ncoll(4),nperp(4)
    real(dp) :: xx3
    real(dp) :: mysign
    real(dp) :: n4(3), n5(3), nlept(3,2), E4, E5, eta56, eta57, nlept1(3), nlept2(3), weight
    real(dp) :: n5test(4), nlept1test(4), nlept2test(4)

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S4Lim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(S4S5Lim)        !4
    call initialize_config(C4Lim)          !5
    call initialize_config(C4S4Lim)        !6
    call initialize_config(C4S5Lim)        !7
    call initialize_config(C4S4S5Lim)      !8
    call initialize_config(C5Lim)          !9
    call initialize_config(C5S4Lim)        !10
    call initialize_config(C5S5Lim)        !11
    call initialize_config(C5S4S5Lim)      !12
    call initialize_config(C4C5Lim)        !13
    call initialize_config(C4C5S4Lim)      !14
    call initialize_config(C4C5S5Lim)      !15
    call initialize_config(C4C5S4S5Lim)    !16

    !-- this is not ideal, but for a check
    
    HardProc%npart = 6    
    S4Lim%npart = 5           
    S4S5Lim%npart = 4         
    S5Lim%npart = 5           
    C4Lim%npart = 5           
    C4S4Lim%npart = 5         
    C4S4S5Lim%npart = 4       
    C4S5Lim%npart =  4        
    C5Lim%npart =  5          
    C5S4Lim%npart = 4         
    C5S4S5Lim%npart = 4       
    C5S5Lim%npart = 5         
    C4C5Lim%npart = 4         
    C4C5S4Lim%npart = 4       
    C4C5S5Lim%npart = 4       
    C4C5S4S5Lim%npart = 4    

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a]     
    S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]            
    S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]            
    C4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]            
    C4S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]          
    C4S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4S5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]            
    C5S4Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C5S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]          
    C4C5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C4C5S4Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4C5S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4C5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]     

    !    MV2 and Y are generated

    MV =  sqrt(MV2)


    x1 =   yr(1) 
    x2 =   yr(2) 
    x3 =   yr(3) 
    x4 =   yr(5) 

    !x1 = 1d-6
    !x2 = 1d-6
    !x3 = 1d-12
    !x4 = 1d-12

    if (yr(10) .ge. half) then
       mysign = +one
    else
       mysign = -one
    endif

    if (withmax) then
       xx3 = x3
    else
       xx3 = zero
    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = one -  two*x3 
    sin4 = sqrt(one-cos4**2)

    phi5 = twopi*yr(6)
    cos5 = -one + two*x4
    sin5 = sqrt(one-cos5**2) 

    eta42 = (one - x3)
    eta52 = x4
    eta41 = x3 
    eta51 = (one - x4) 
    eta54 =   half*( one - sin4*sin5*cos(phi4-phi5) -cos4*cos5)

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,HardProc%flag)

    if (HardProc%flag .eqv. .false.) then 
       HardProc%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       n5 = (/sin5*cos(phi5),  sin5*sin(phi5), cos5 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p4 - p5       
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,HardProc)

       !nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       !nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       !eta56 = half*(one-sc3(n5,nlept1))
       !eta57 = half*(one-sc3(n5,nlept2))

       nlept1test = (/one, nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2test = (/one, nlept(1,2), nlept(2,2), nlept(3,2)/)
       n5test = (/one, sin5*cos(phi5),  sin5*sin(phi5), cos5/)

        eta56 = half*scr(nlept1test,n5test)
        eta57 = half*scr(nlept2test,n5test)
        
       

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2,weight)

       !reshuffle momenta for the amplitude
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       !HardProc%AmpMom(:,3) = HardProc%AmpMom(:,5) !lepton 1
       !HardProc%AmpMom(:,4) = HardProc%AmpMom(:,6) !lepton 2
       HardProc%AmpMom(:,5) = p4 ! gluon 
       HardProc%AmpMom(:,6) = p5 ! photon

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2


    endif



    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: S4S5 FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! -------------- DOUBLE SOFT AND DOUBLE SOFT + Collinear --------------------

    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    S4S5Lim%PartFrac = (/xi1,xi2/)

    p1 = half*sqrts*(/one,zero,zero,one/) 
    p2 = half*sqrts*(/one,zero,zero,-one/)
    E4 = half*sqrts*x1
    E5 = half*sqrts*x2
    n4 = (/sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
    n5 = (/sin5*cos(phi5),  sin5*sin(phi5), cos5 /)
    p4 = E4*(/one,n4(1),n4(2),n4(3)/)
    p5 = E5*(/one,n5(1),n5(2),n5(3)/)


    pVV = p1 + p2
    call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,S4S5Lim)

    !nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
    !nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

    !eta56 = half*(one-sc3(n5,nlept1))
    !eta57 = half*(one-sc3(n5,nlept2))

     nlept1test = (/one, nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2test = (/one, nlept(1,2), nlept(2,2), nlept(3,2)/)
       n5test = (/one, sin5*cos(phi5),  sin5*sin(phi5), cos5/)

        eta56 = half*scr(nlept1test,n5test)
        eta57 = half*scr(nlept2test,n5test)
       

    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2,weight)

    S4S5Lim%AmpMom(:,1) = p1
    S4S5Lim%AmpMom(:,2) = p2
    !S4S5Lim%AmpMom(:,3) = S4S5Lim%AmpMom(:,3)  lepton
    !S4S5Lim%AmpMom(:,4) = S4S5Lim%AmpMom(:,4)  antilepton

    S4S5Lim%LimMom(:,1) = p1
    S4S5Lim%LimMom(:,2) = p2
    S4S5Lim%LimMom(:,3) = S4S5Lim%AmpMom(:,3)
    S4S5Lim%LimMom(:,4) = S4S5Lim%AmpMom(:,4)
    S4S5Lim%LimMom(:,5) = p4  !gluon
    S4S5Lim%LimMom(:,6) = p5  !photon

    S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1*x2 & 
         * weight &
         * one/two/spart &
         * spart/MV2

    

    ! ------------   C4  +  S4S5   -----------

    eta42 = one
    eta52 = x4 
    eta41 = zero  
    eta51 = (one - x4) 
    eta54 = eta52

    s14 = spart *x1*x3 
    s15 = spart*x2*(one-x4) 
    s25 = spart*x2*x4
    s12 = two*scr(p1,p2) 
    z4 = x1 

    C4S4S5Lim%Lim_KinInv(1:5) = (/s12,s25,s15,s14,z4/)

    C4S4S5Lim%LimMom(:,1) = p1
    C4S4S5Lim%LimMom(:,2) = p2
    C4S4S5Lim%LimMom(:,3) = S4S5Lim%LimMom(:,3)
    C4S4S5Lim%LimMom(:,4) = S4S5Lim%LimMom(:,4)
    C4S4S5Lim%LimMom(:,5) = p5

       nlept1test = (/one, nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2test = (/one, nlept(1,2), nlept(2,2), nlept(3,2)/)
       n5test = (/one, sin5*cos(phi5),  sin5*sin(phi5), cos5/)

        eta56 = half*scr(nlept1test,n5test)
        eta57 = half*scr(nlept2test,n5test)

    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2,weight)
    
    C4S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1*x2 & 
         * weight &
         * one/two/spart & 
         * spart/MV2


    ! ------------  C51 S4S5 FLM  ------------

    phi4 = twopi*yr(4) 
    cos4 = one - two*x3 
    sin4 = sqrt(one-cos4**2)

    eta41 = x3
    eta51 = one
    eta42 = (one-x3) 
    eta52 = zero
    eta54 = (one-x3)

    n5 = (/zero, zero, -one/)

    s12 = two*scr(p1,p2) 
    s14 = spart*x1*x3
    s24 = spart*x1*(one - x3) 

    s25 = spart*x2*(x4)
    z5 = x2

    C5S4S5Lim%Lim_KinInv(1:5) = (/s12,s14,s24,s25,z5/)

    !eta56 = half*(one-sc3(n5,nlept1))
    !eta57 = half*(one-sc3(n5,nlept2))

    nlept1test = (/one, nlept(1,1), nlept(2,1), nlept(3,1)/)
    nlept2test = (/one, nlept(1,2), nlept(2,2), nlept(3,2)/)
    n5test = (/one, zero, zero, -one/)

        eta56 = half*scr(nlept1test,n5test)
        eta57 = half*scr(nlept2test,n5test)

    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2,weight)

    C5S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1*x2 & 
         * weight &
         * one/two/spart & 
         * spart/MV2



    ! -----------   C4  + C5 + S4S5  ----------

    eta41 = zero
    eta51 = one
    eta42 = one
    eta52 = zero
    eta54 = one

    s14 = spart *x1*x3 
    s25 = spart *x2*x4
    z4 = x1 
    z5 = x2

    C4C5S4S5Lim%Lim_KinInv(1:5) = (/s14,s25,zero,z4,z5/)

    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2,weight)
    
    C4C5S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1*x2 & 
         * weight &
         * one/two/spart & 
         * spart/MV2



    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3: S5 FLM   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------


    ! --------------- SOFT 5 ----------

    phi4 = twopi*yr(4) 
    cos4 = one -  two*x3 
    sin4 = sqrt(one-cos4**2)

    phi5 = twopi*yr(6)
    cos5 = -one + two*x4
    sin5 = sqrt(one-cos5**2) 

    eta42 = (one - x3)
    eta52 = x4
    eta41 = x3 
    eta51 = (one - x4) 
    eta54 =   half*( one - sin4*sin5*cos(phi4-phi5) -cos4*cos5)

    call get_flag_nlo(x1,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts,S5Lim%flag)

    if (S5Lim%flag .eqv. .false.) then 

       S5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/ sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       n5 = (/sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)


       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       !S5Lim%AmpMom(:,3) = S5Lim%AmpMom(:,4) ! lepton
       !S5Lim%AmpMom(:,4) = S5Lim%AmpMom(:,5) ! antilepton
       S5Lim%AmpMom(:,5) = p4

       S5Lim%LimMom(:,1:4) = S5Lim%AmpMom(:,1:4)
       S5Lim%LimMom(:,5) = p5 !photon

       !nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       !nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       !eta56 = half*(one-sc3(n5,nlept1))
       !eta57 = half*(one-sc3(n5,nlept2))

       nlept1test = (/one, nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2test = (/one, nlept(1,2), nlept(2,2), nlept(3,2)/)
       n5test = (/one, sin5*cos(phi5),  sin5*sin(phi5), cos5/)

        eta56 = half*scr(nlept1test,n5test)
        eta57 = half*scr(nlept2test,n5test)

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2, weight)


       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2


       ! ---------  C52 + S5 --------------

       !n5 = (/zero,zero,-one/)

       eta41 = x3
       eta51 = one
       eta42 = (one-x3) 
       eta52 = zero
       eta54 = (one-x3)

        !eta56 = half*(one-sc3(n5,nlept1))
        !eta57 = half*(one-sc3(n5,nlept2))
        
       nlept1test = (/one, nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2test = (/one, nlept(1,2), nlept(2,2), nlept(3,2)/)
       n5test = (/one, zero,zero,-one/)

        eta56 = half*scr(nlept1test,n5test)
        eta57 = half*scr(nlept2test,n5test)

      

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2, weight)

       s25 = spart*x2*x4
       z5 = x2

       C5S5Lim%Lim_KinInv(1:5) = (/zero,s25,zero,zero,z5/)

       C5S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2


    endif

    ! --------------- SOFT 4 ----------

    phi4 = twopi*yr(4) 
    cos4 = one -  two*x3 
    sin4 = sqrt(one-cos4**2)

    phi5 = twopi*yr(6)
    cos5 = -one + two*x4
    sin5 = sqrt(one-cos5**2) 

    eta42 = (one - x3)
    eta52 = x4
    eta41 = x3 
    eta51 = (one - x4) 
    eta54 =   half*( one - sin4*sin5*cos(phi4-phi5) -cos4*cos5)


    call get_flag_nlo(x2,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,S4Lim%flag)

    if (S4Lim%flag .eqv. .false.) then 

       S4Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/ sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       n5 = (/sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p5
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,S4Lim)

       !nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       !nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       !eta56 = half*(one-sc3(n5,nlept1))
       !eta57 = half*(one-sc3(n5,nlept2))

       nlept1test = (/one, nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2test = (/one, nlept(1,2), nlept(2,2), nlept(3,2)/)
       n5test = (/one, sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)

        eta56 = half*scr(nlept1test,n5test)
        eta57 = half*scr(nlept2test,n5test)
        
       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2, weight)


       S4Lim%AmpMom(:,1) = p1
       S4Lim%AmpMom(:,2) = p2
       !S4Lim%AmpMom(:,3) = S4Lim%AmpMom(:,4) ! lepton
       !S4Lim%AmpMom(:,4) = S4Lim%AmpMom(:,5) ! antilepton
       S4Lim%AmpMom(:,5) = p5


       S4Lim%LimMom(:,1) = p1
       S4Lim%LimMom(:,2) = p2
       S4Lim%LimMom(:,3) = p4 
       S4Lim%LimMom(:,4) = p5 

       S4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

       ! ---------  C4 + S4 --------------

       eta41 = zero
       eta51 = one-x4
       eta42 = one 
       eta52 = x4
       eta54 = (one-x4)      

       s14 = spart*x1*x3
       z4 = x1

       C4S4Lim%Lim_KinInv(1:5) = (/zero,s14,zero,zero,z4/)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2, weight)


       C4S4Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2

    endif



    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3:   C41  sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    eta42 = one
    eta52 = x4 
    eta41 = zero  
    eta51 = (one - x4) 
    eta54 =   eta51

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,C4Lim%flag)

    if (C4Lim%flag .eqv. .false.) then 

       C4Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/ zero, zero, one /)
       n5 = (/sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)


       s14 = spart *x1*x3
       z4 = one/(one-x1)

       C4Lim%Lim_KinInv(1:5) = (/s14,zero,zero, z4, zero/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C4Lim)

       !nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       !nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       !eta56 = half*(one-sc3(n5,nlept1))
       !eta57 = half*(one-sc3(n5,nlept2))

       nlept1test = (/one, nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2test = (/one, nlept(1,2), nlept(2,2), nlept(3,2)/)
       n5test = (/one, sin5*cos(phi5),  sin5*sin(phi5), cos5/)

        eta56 = half*scr(nlept1test,n5test)
        eta57 = half*scr(nlept2test,n5test)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2, weight)

       C4Lim%AmpMom(:,1) = p1 -p4
       C4Lim%AmpMom(:,2) = p2 
       !C4Lim%AmpMom(:,3) = C4Lim%AmpMom(:,4)
       !C4Lim%AmpMom(:,4) = C4Lim%AmpMom(:,5)
       C4Lim%AmpMom(:,5) = p5

       C4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2

    endif



    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C41  + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------



    call get_flag_nlo(x1,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts, C4S5Lim%flag)

    if (C4S5Lim%flag .eqv. .false.) then 

       C4S5Lim%PartFrac = (/xi1,xi2/)


       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/ zero, zero, one /)
       n5 = (/sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       s15 = spart*x2*(one-x4) 
       s25 = spart*x2*x4
       s12 = two*scr(p1,p2) 

       z4 = one/(one - x1)
       s14 = spart *x1*x3

       C4S5Lim%Lim_KinInv(1:5) = (/s12,s25,s15,s14,z4/)

       pVV = p1 + p2 - p4

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C4S5Lim)

      ! nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
      ! nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

      ! eta56 = half*(one-sc3(n5,nlept1))
       ! eta57 = half*(one-sc3(n5,nlept2))


       nlept1test = (/one, nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2test = (/one, nlept(1,2), nlept(2,2), nlept(3,2)/)
       n5test = (/one, sin5*cos(phi5),  sin5*sin(phi5), cos5/)

        eta56 = half*scr(nlept1test,n5test)
        eta57 = half*scr(nlept2test,n5test)

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2, weight)


       C4S5Lim%AmpMom(:,1) = p1 - p4
       C4S5Lim%AmpMom(:,2) = p2 
       !C4S5Lim%AmpMom(:,3) = C4S5Lim%AmpMom(:,4) ! lepton
       !C4S5Lim%AmpMom(:,4) = C4S5Lim%AmpMom(:,5) ! antilepton

       C4S5Lim%LimMom(:,1) = p1- p4
       C4S5Lim%LimMom(:,2) = p2 
       C4S5Lim%LimMom(:,3) = C4S5Lim%AmpMom(:,3)
       C4S5Lim%LimMom(:,4) = C4S5Lim%AmpMom(:,4) 
       C4S5Lim%LimMom(:,5) = p5 !photon

       C4S5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2

   
       ! ---------   C41   + C52 + S5 ------------

       eta41 = zero
       eta51 = one
       eta42 = one
       eta52 = zero
       eta54 = one           

       !n5 = (/zero,zero,-one/)
       
       z4 = one/(one-x1) 
       z5 = x2

       s14 = spart *x1*x3 
       s25 = spart *x2*x4

       C4C5S5Lim%Lim_KinInv(1:5) = (/s14,s25,zero,z4,z5/)

       !eta56 = half*(one-sc3(n5,nlept1))
       !eta57 = half*(one-sc3(n5,nlept2))

        nlept1test = (/one, nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2test = (/one, nlept(1,2), nlept(2,2), nlept(3,2)/)
       n5test = (/one, zero,zero,-one/)

        eta56 = half*scr(nlept1test,n5test)
        eta57 = half*scr(nlept2test,n5test)

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2, weight)

       C4C5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2

    endif




    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: C52 FLM 
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = one - two*x3 
    sin4 = sqrt(one-cos4**2)

    eta41 = x3
    eta51 = one
    eta42 = (one-x3) 
    eta52 = zero
    eta54 = (one-x3)

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,C5Lim%flag)
  
    if (C5Lim%flag .eqv. .false.) then 

       C5Lim%PartFrac = (/xi1,xi2/)


       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/ sin4*cos(phi4), sin4*sin(phi4), cos4 /)
       n5 = (/zero, zero, -one /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)


       s25 = spart*x2*x4
       z5 = one/(one - x2) 

       C5Lim%Lim_KinInv(1:5) = (/zero,s25,zero,zero,z5/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C5Lim)

       C5Lim%AmpMom(:,1) = p1 
       C5Lim%AmpMom(:,2) = p2  - p5      
       !C5Lim%AmpMom(:,3) = C5Lim%AmpMom(:,4) !lepton
       !C5Lim%AmpMom(:,4) = C5Lim%AmpMom(:,5) !anti lepton
       C5Lim%AmpMom(:,5) = p4

       
      ! print*,'++++++++++++++++++'
      ! print*, 'C5Lim%AmpMom(:,1)', C5Lim%AmpMom(:,1)
      ! print*, 'C5Lim%AmpMom(:,2)', C5Lim%AmpMom(:,2)
      ! print*, 'C5Lim%AmpMom(:,3)', C5Lim%AmpMom(:,3)
      ! print*, 'C5Lim%AmpMom(:,4)', C5Lim%AmpMom(:,4)
      ! print*, 'C5Lim%AmpMom(:,5)', C5Lim%AmpMom(:,5)
      ! print*,'++++++++++++++++++'
       
       !eta56 = half*(one-sc3(n5,nlept1))
       !eta57 = half*(one-sc3(n5,nlept2))

       nlept1test = (/one, nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2test = (/one, nlept(1,2), nlept(2,2), nlept(3,2)/)
       n5test = (/one, zero, zero, -one /)

        eta56 = half*scr(nlept1test,n5test)
        eta57 = half*scr(nlept2test,n5test)

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2, weight)

       C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2



    !print*, '****************'
    !print*, 'x1', x1
    !print*, 'x2', x2
    !print*, 'x3', x3
    !print*, 'x4', x4
    !print*, 'C5Lim%wgt', C5Lim%wgt
    !print*, 'weight', weight 
    !print*, 'eta56', eta56
   !print*, 'eta57', eta57
  !print*, 'eta41', eta41
 !print*, 'eta42', eta42 
  ! print*, 'eta41', eta51
 !print*, 'eta42', eta52

!print*, 

!print*, 'eta54', eta54
    
    
 !   print*, '****************'
    endif




    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C41  + C52
    !-------------------------------------------------------------------------------------------------------------

    eta41 = zero
    eta51 = one
    eta42 = one
    eta52 = zero
    eta54 = one

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,C4C5Lim%flag)

    if (C4C5Lim%flag .eqv. .false.) then 

       C4C5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/ zero, zero, one/)
       n5 = (/zero, zero, -one /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       s14 = spart *x1*x3 
       s25 = spart *x2*x4
       z4 = one/(one - x1)
       z5 = one/(one - x2) 

       C4C5Lim%Lim_KinInv(1:5) = (/s14,s25,zero,z4,z5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C4C5Lim)

       C4C5Lim%AmpMom(:,1) = p1 - p4 
       C4C5Lim%AmpMom(:,2) = p2 - p5

       !eta56 = half*(one-sc3(n5,nlept1))
       !eta57 = half*(one-sc3(n5,nlept2))

       nlept1test = (/one, nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2test = (/one, nlept(1,2), nlept(2,2), nlept(3,2)/)
       n5test = (/one, zero, zero, -one /)

        eta56 = half*scr(nlept1test,n5test)
        eta57 = half*scr(nlept2test,n5test)
       

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2, weight)

       C4C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2


    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: C5 + S4 
    !-------------------------------------------------------------------------------------------------------------

    phi4 =  twopi*yr(4) 
    cos4 = one - two*x3 
    sin4 =  sqrt(one-cos4**2)

    eta41 = x3
    eta51 = one
    eta42 = (one-x3) 
    eta52 = zero
    eta54 = (one-x3)


    call get_flag_nnlo(MV2,Y,zero,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,C5S4Lim%flag)

    if (C5S4Lim%flag .eqv. .false.) then 

       C5S4Lim%PartFrac = (/xi1,xi2/)


       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin4*cos(phi4), sin4*sin(phi4), cos4/)
       n5 = (/zero, zero, -one /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       s12 = two*scr(p1,p2)
       s24 = spart*x1*(one-x3) 
       s14 = spart*x1*x3
       s25 = spart*x2*x4
       z5 = one/(one - x2) 

       C5S4Lim%Lim_KinInv(1:5) = (/s12,s14,s24,s25,z5/)

       pVV = p1 + p2 - p5
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C5S4Lim)

       C5S4Lim%AmpMom(:,1) = p1
       C5S4Lim%AmpMom(:,2) = p2  - p5 
       !C5S4Lim%AmpMom(:,3) = C5S4Lim%AmpMom(:, 4)
       !C5S4Lim%AmpMom(:,4) = C5S4Lim%AmpMom(:, 5)

       !eta56 = half*(one-sc3(n5,nlept1))
       !eta57 = half*(one-sc3(n5,nlept2))

       nlept1test = (/one, nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2test = (/one, nlept(1,2), nlept(2,2), nlept(3,2)/)
       n5test = (/one,zero, zero, -one  /)

        eta56 = half*scr(nlept1test,n5test)
        eta57 = half*scr(nlept2test,n5test)

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2, weight)

       C5S4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2


       !------ C5 C4 S4

       eta42 = one 

       z4 = x1
       s24 = spart*x1
       s14 = spart*x1*x3
       s25 = spart*x2*x4
       z5 = one/(one - x2) 

       C4C5S4Lim%Lim_KinInv(1:5) = (/z4,s24,s14,s25,z5/)

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,2, weight)


       C4C5S4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2



    endif



  end subroutine kinematics_nnlo_5162






  

  subroutine kinematics_nnlo_5261(MV2,Y,yr,&
       HardProc,S4Lim,S5Lim,S4S5Lim,C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,&
       C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim)
    real(dp), intent(in) ::  yr(10), MV2, Y
    type(KinConfig)      :: HardProc,S4Lim,S4S5Lim,S5Lim,C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,&
         C5Lim,C5S4Lim,C5S4S5Lim,C5S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: x1, x2,x3, x4, lambda,  cos4, sin4, phi4, phi5, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52, eta54, sinphi45, cosphi45, phi45, MV2aux
    real(dp) :: Deltan, Deltad, w4151, s14, s15, s45, z4, z5, ztot, weight5261
    real(dp) :: spart,sqrts, s25, s12, s24, yy, Deltas
    real(dp) :: ncoll(4),nperp(4)
    real(dp) :: xx3
    real(dp) :: mysign
    real(dp) :: n4(3), n5(3), nlept(3,2), E4, E5, eta56, eta57, nlept1(3), nlept2(3), weight

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S4Lim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(S4S5Lim)        !4
    call initialize_config(C4Lim)          !5
    call initialize_config(C4S4Lim)        !6
    call initialize_config(C4S5Lim)        !7
    call initialize_config(C4S4S5Lim)      !8
    call initialize_config(C5Lim)          !9
    call initialize_config(C5S4Lim)        !10
    call initialize_config(C5S5Lim)        !11
    call initialize_config(C5S4S5Lim)      !12
    call initialize_config(C4C5Lim)        !13
    call initialize_config(C4C5S4Lim)      !14
    call initialize_config(C4C5S5Lim)      !15
    call initialize_config(C4C5S4S5Lim)    !16

    HardProc%npart = 6    
    S4Lim%npart = 5           
    S4S5Lim%npart = 4         
    S5Lim%npart = 5           
    C4Lim%npart = 5           
    C4S4Lim%npart = 5         
    C4S4S5Lim%npart = 4       
    C4S5Lim%npart =  4        
    C5Lim%npart =  5          
    C5S4Lim%npart = 4         
    C5S4S5Lim%npart = 4       
    C5S5Lim%npart = 5         
    C4C5Lim%npart = 4         
    C4C5S4Lim%npart = 4       
    C4C5S5Lim%npart = 4       
    C4C5S4S5Lim%npart = 4    

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a]     
    S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]            
    S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]            
    C4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]            
    C4S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]          
    C4S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4S5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]            
    C5S4Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C5S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]          
    C4C5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C4C5S4Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4C5S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4C5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]     

    
    !    MV2 and Y are generated

    MV =  sqrt(MV2)


    x1 =   yr(1) 
    x2 =   yr(2) 
    x3 =   yr(3) 
    x4 =   yr(5) 

    !x1 = 1d-6
    !x2 = 1d-6
    !x3 = 1d-12
    !x4 = 1d-12

    if (yr(10) .ge. half) then
       mysign = +one
    else
       mysign = -one
    endif

    if (withmax) then
       xx3 = x3
    else
       xx3 = zero
    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 1: I * FLM  
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = -one +  two*x3 
    sin4 = sqrt(one-cos4**2)

    phi5 = twopi*yr(6)
    cos5 = one - two*x4
    sin5 = sqrt(one-cos5**2) 

    eta41 = (one - x3)
    eta51 = x4
    eta42 = x3 
    eta52 = (one - x4) 
    eta54 =   half*( one - sin4*sin5*cos(phi4-phi5) -cos4*cos5)

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,HardProc%flag)

    if (HardProc%flag .eqv. .false.) then 
       HardProc%PartFrac = (/xi1,xi2/)

       p1 = half*sqrts*(/one,zero,zero,one/) 
       p2 = half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       n5 = (/sin5*cos(phi5),  sin5*sin(phi5), cos5 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p4 - p5       
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,HardProc)

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1,weight)

       !reshuffle momenta for the amplitude
       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       !HardProc%AmpMom(:,3) = HardProc%AmpMom(:,5) !lepton 1
       !HardProc%AmpMom(:,4) = HardProc%AmpMom(:,6) !lepton 2
       HardProc%AmpMom(:,5) = p4 ! gluon 
       HardProc%AmpMom(:,6) = p5 ! photon

       HardProc%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2


    endif



    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: S4S5 FL(1,2,4,5)  ! angles etc. identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! -------------- DOUBLE SOFT AND DOUBLE SOFT + Collinear --------------------

    spart = MV2
    sqrts = sqrt(spart)
    yy = Y 

    xi1 = sqrt(spart/sh)*exp(yy) 
    xi2 = sqrt(spart/sh)*exp(-yy)

    S4S5Lim%PartFrac = (/xi1,xi2/)

    p1 = half*sqrts*(/one,zero,zero,one/) 
    p2 = half*sqrts*(/one,zero,zero,-one/)
    E4 = half*sqrts*x1
    E5 = half*sqrts*x2
    n4 = (/sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
    n5 = (/sin5*cos(phi5),  sin5*sin(phi5), cos5 /)
    p4 = E4*(/one,n4(1),n4(2),n4(3)/)
    p5 = E5*(/one,n5(1),n5(2),n5(3)/)


    pVV = p1 + p2
    call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,S4S5Lim)

    nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
    nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

    eta56 = half*(one-sc3(n5,nlept1))
    eta57 = half*(one-sc3(n5,nlept2))

    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1,weight)

    S4S5Lim%AmpMom(:,1) = p1
    S4S5Lim%AmpMom(:,2) = p2
    !S4S5Lim%AmpMom(:,3) = S4S5Lim%AmpMom(:,3)  lepton
    !S4S5Lim%AmpMom(:,4) = S4S5Lim%AmpMom(:,4)  antilepton

    S4S5Lim%LimMom(:,1) = p1
    S4S5Lim%LimMom(:,2) = p2
    S4S5Lim%LimMom(:,3) = S4S5Lim%AmpMom(:,3)
    S4S5Lim%LimMom(:,4) = S4S5Lim%AmpMom(:,4)
    S4S5Lim%LimMom(:,5) = p4  !gluon
    S4S5Lim%LimMom(:,6) = p5  !photon

    S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1*x2 & 
         * weight &
         * one/two/spart &
         * spart/MV2


    ! ------------   C4  +  S4S5   -----------

    eta41 = one
    eta51 = x4 
    eta42 = zero  
    eta52 = (one - x4) 
    eta54 = eta52

    s24 = spart *x1*x3 
    s25 = spart*x2*(one-x4) 
    s15 = spart*x2*x4
    s12 = two*scr(p1,p2) 
    z4 = x1 

    C4S4S5Lim%Lim_KinInv(1:5) = (/s12,s15,s25,s24,z4/)

    C4S4S5Lim%LimMom(:,1) = p1
    C4S4S5Lim%LimMom(:,2) = p2
    C4S4S5Lim%LimMom(:,3) = S4S5Lim%LimMom(:,3)
    C4S4S5Lim%LimMom(:,4) = S4S5Lim%LimMom(:,4)
    C4S4S5Lim%LimMom(:,5) = p5

    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1,weight)
    
    C4S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1*x2 & 
         * weight &
         * one/two/spart & 
         * spart/MV2


    ! ------------  C51 S4S5 FLM  ------------

    phi4 = twopi*yr(4) 
    cos4 = -one + two*x3 
    sin4 = sqrt(one-cos4**2)

    eta42 = x3
    eta52 = one
    eta41 = (one-x3) 
    eta51 = zero
    eta54 = (one-x3)

    n5 = (/zero, zero, one/)

    s12 = two*scr(p1,p2) 
    s24 = spart*x1*x3
    s14 = spart*x1*(one - x3) 

    s15 = spart*x2*(x4)
    z5 = x2

    C5S4S5Lim%Lim_KinInv(1:5) = (/s12,s14,s24,s15,z5/)

    eta56 = half*(one-sc3(n5,nlept1))
    eta57 = half*(one-sc3(n5,nlept2))

    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1,weight)

    C5S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1*x2 & 
         * weight &
         * one/two/spart & 
         * spart/MV2



    ! -----------   C4  + C5 + S4S5  ----------

    eta42 = zero
    eta52 = one
    eta41 = one
    eta51 = zero
    eta54 = one

    s24 = spart *x1*x3 
    s15 = spart *x2*x4
    z4 = x1 
    z5 = x2

    C4C5S4S5Lim%Lim_KinInv(1:5) = (/s24,s15,zero,z4,z5/)

    call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1,weight)
    
    C4C5S4S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
         * x1*x2 & 
         * weight &
         * one/two/spart & 
         * spart/MV2



    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3: S5 FLM   ! angles the same as in hard 
    !-------------------------------------------------------------------------------------------------------------


    ! --------------- SOFT 5 ----------

    phi4 = twopi*yr(4) 
    cos4 = -one +  two*x3 
    sin4 = sqrt(one-cos4**2)

    phi5 = twopi*yr(6)
    cos5 = one - two*x4
    sin5 = sqrt(one-cos5**2) 

    eta41 = (one - x3)
    eta51 = x4
    eta42 = x3 
    eta52 = (one - x4) 
    eta54 =   half*( one - sin4*sin5*cos(phi4-phi5) -cos4*cos5)

    call get_flag_nlo(x1,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts,S5Lim%flag)

    if (S5Lim%flag .eqv. .false.) then 

       S5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/ sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       n5 = (/sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)


       pVV = p1 + p2 - p4
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,S5Lim)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       !S5Lim%AmpMom(:,3) = S5Lim%AmpMom(:,4) ! lepton
       !S5Lim%AmpMom(:,4) = S5Lim%AmpMom(:,5) ! antilepton
       S5Lim%AmpMom(:,5) = p4

       S5Lim%LimMom(:,1:4) = S5Lim%AmpMom(:,1:4)
       S5Lim%LimMom(:,5) = p5 !photon

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1, weight)


       S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2


       ! ---------  C5 + S5 --------------

       n5 = (/zero,zero,one/)

       eta42 = x3
       eta52 = one
       eta41 = (one-x3) 
       eta51 = zero
       eta54 = (one-x3)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1, weight)

       s15 = spart*x2*x4
       z5 = x2

       C5S5Lim%Lim_KinInv(1:5) = (/zero,s15,zero,zero,z5/)

       C5S5Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2


    endif

    ! --------------- SOFT 4 ----------

    phi4 = twopi*yr(4) 
    cos4 = -one +  two*x3 
    sin4 = sqrt(one-cos4**2)

    phi5 = twopi*yr(6)
    cos5 = one - two*x4
    sin5 = sqrt(one-cos5**2) 

    eta41 = (one - x3)
    eta51 = x4
    eta42 = x3 
    eta52 = (one - x4) 
    eta54 =   half*( one - sin4*sin5*cos(phi4-phi5) -cos4*cos5)


    call get_flag_nlo(x2,MV2,Y,x2,eta51,eta52,xi1,xi2,spart,sqrts,S4Lim%flag)

    if (S4Lim%flag .eqv. .false.) then 

       S4Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/ sin4*cos(phi4),  sin4*sin(phi4), cos4 /)
       n5 = (/sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       pVV = p1 + p2 - p5
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,S4Lim)

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_5262(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1, weight)


       S4Lim%AmpMom(:,1) = p1
       S4Lim%AmpMom(:,2) = p2
       !S4Lim%AmpMom(:,3) = S4Lim%AmpMom(:,4) ! lepton
       !S4Lim%AmpMom(:,4) = S4Lim%AmpMom(:,5) ! antilepton
       S4Lim%AmpMom(:,5) = p5


       S4Lim%LimMom(:,1) = p1
       S4Lim%LimMom(:,2) = p2
       S4Lim%LimMom(:,3) = p4 
       S4Lim%LimMom(:,4) = p5 

       S4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart &
            * spart/MV2

       ! ---------  C4 + S4 --------------

       eta42 = zero
       eta52 = one-x4
       eta41 = one 
       eta51 = x4
       eta54 = (one-x4)      

       s24 = spart*x1*x3
       z4 = x1

       C4S4Lim%Lim_KinInv(1:5) = (/zero,s24,zero,zero,z4/)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1, weight)


       C4S4Lim%wgt =  one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2

    endif



    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3:   C42  sij = 2*pi*pj
    !-------------------------------------------------------------------------------------------------------------

    eta41 = one
    eta51 = x4 
    eta42 = zero  
    eta52 = (one - x4) 
    eta54 =   eta52

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,C4Lim%flag)

    if (C4Lim%flag .eqv. .false.) then 

       C4Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/ zero, zero, -one /)
       n5 = (/sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)


       s24 = spart *x1*x3
       z4 = one/(one-x1)

       C4Lim%Lim_KinInv(1:5) = (/s24,zero,zero, z4, zero/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C4Lim)

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1, weight)

       C4Lim%AmpMom(:,1) = p1 
       C4Lim%AmpMom(:,2) = p2 - p4
       !C4Lim%AmpMom(:,3) = C4Lim%AmpMom(:,4)
       !C4Lim%AmpMom(:,4) = C4Lim%AmpMom(:,5)
       C4Lim%AmpMom(:,5) = p5

       C4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2 

    endif



    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C42  + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------



    call get_flag_nlo(x1,MV2,Y,x2,eta41,eta42,xi1,xi2,spart,sqrts, C4S5Lim%flag)

    if (C4S5Lim%flag .eqv. .false.) then 

       C4S5Lim%PartFrac = (/xi1,xi2/)


       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/ zero, zero, -one /)
       n5 = (/sin5*cos(phi5),  sin5*sin(phi5),  cos5 /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       s25 = spart*x2*(one-x4) 
       s15 = spart*x2*x4
       s12 = two*scr(p1,p2) 

       z4 = one/(one - x1)
       s24 = spart *x1*x3

       C4S5Lim%Lim_KinInv(1:5) = (/s12,s15,s25,s24,z4/)

       pVV = p1 + p2 - p4

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C4S5Lim)

       nlept1 = (/nlept(1,1), nlept(2,1), nlept(3,1)/)
       nlept2 = (/nlept(1,2), nlept(2,2), nlept(3,2)/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1, weight)


       C4S5Lim%AmpMom(:,1) = p1 
       C4S5Lim%AmpMom(:,2) = p2 - p4
       !C4S5Lim%AmpMom(:,3) = C4S5Lim%AmpMom(:,4) ! lepton
       !C4S5Lim%AmpMom(:,4) = C4S5Lim%AmpMom(:,5) ! antilepton

       C4S5Lim%LimMom(:,1) = p1
       C4S5Lim%LimMom(:,2) = p2 - p4
       C4S5Lim%LimMom(:,3) = C4S5Lim%AmpMom(:,3)
       C4S5Lim%LimMom(:,4) = C4S5Lim%AmpMom(:,4) 
       C4S5Lim%LimMom(:,5) = p5 !photon

       C4S5Lim%wgt= one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2

   
       ! ---------   C4   + C51 + S5 ------------

       eta42 = zero
       eta52 = one
       eta41 = one
       eta51 = zero
       eta54 = one           

       n5 = (/zero,zero,one/)
       
       z4 = one/(one-x1) 
       z5 = x2

       s24 = spart *x1*x3 
       s15 = spart *x2*x4

       C4C5S5Lim%Lim_KinInv(1:5) = (/s24,s15,zero,z4,z5/)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1, weight)

       C4C5S5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2

    endif




    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: C51 FLM 
    !-------------------------------------------------------------------------------------------------------------

    phi4 = twopi*yr(4) 
    cos4 = -one + two*x3 
    sin4 = sqrt(one-cos4**2)

    eta42 = x3
    eta52 = one
    eta41 = (one-x3) 
    eta51 = zero
    eta54 = (one-x3)

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,C5Lim%flag)
  
    if (C5Lim%flag .eqv. .false.) then 

       C5Lim%PartFrac = (/xi1,xi2/)


       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/ sin4*cos(phi4), sin4*sin(phi4), cos4 /)
       n5 = (/zero, zero, one /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)


       s15 = spart*x2*x4
       z5 = one/(one - x2) 

       C5Lim%Lim_KinInv(1:5) = (/zero,s15,zero,zero,z5/)

       pVV = p1 + p2 - p4 - p5
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C5Lim)

       C5Lim%AmpMom(:,1) = p1 - p5 
       C5Lim%AmpMom(:,2) = p2       
       !C5Lim%AmpMom(:,3) = C5Lim%AmpMom(:,4) !lepton
       !C5Lim%AmpMom(:,4) = C5Lim%AmpMom(:,5) !anti lepton
       C5Lim%AmpMom(:,5) = p4

       
       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1, weight)

       C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2

    endif




    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C42  + C51
    !-------------------------------------------------------------------------------------------------------------

    eta42 = zero
    eta52 = one
    eta41 = one
    eta51 = zero
    eta54 = one

    call get_flag_nnlo(MV2,Y,x1,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,C4C5Lim%flag)

    if (C4C5Lim%flag .eqv. .false.) then 

       C4C5Lim%PartFrac = (/xi1,xi2/)

       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/ zero, zero, -one/)
       n5 = (/zero, zero, one /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       s24 = spart *x1*x3 
       s15 = spart *x2*x4
       z4 = one/(one - x1)
       z5 = one/(one - x2) 

       C4C5Lim%Lim_KinInv(1:5) = (/s24,s15,zero,z4,z5/)

       pVV = p1 + p2 - p4 - p5

       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C4C5Lim)

       C4C5Lim%AmpMom(:,1) = p1 - p5 
       C4C5Lim%AmpMom(:,2) = p2 - p4

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1, weight)

       C4C5Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2


    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: C5 + S4 
    !-------------------------------------------------------------------------------------------------------------

    phi4 =  twopi*yr(4) 
    cos4 = -one + two*x3 
    sin4 =  sqrt(one-cos4**2)

    eta42 = x3
    eta52 = one
    eta41 = (one-x3) 
    eta51 = zero
    eta54 = (one-x3)


    call get_flag_nnlo(MV2,Y,zero,x2,eta54,eta51,eta52,eta42,eta41,xi1,xi2,spart,sqrts,C5S4Lim%flag)

    if (C5S4Lim%flag .eqv. .false.) then 

       C5S4Lim%PartFrac = (/xi1,xi2/)


       p1= half*sqrts*(/one,zero,zero,one/) 
       p2 =half*sqrts*(/one,zero,zero,-one/)
       E4 = half*sqrts*x1
       E5 = half*sqrts*x2
       n4 = (/sin4*cos(phi4), sin4*sin(phi4), cos4/)
       n5 = (/zero, zero, one /)
       p4 = E4*(/one,n4(1),n4(2),n4(3)/)
       p5 = E5*(/one,n5(1),n5(2),n5(3)/)

       s12 = two*scr(p1,p2)
       s14 = spart*x1*(one-x3) 
       s24 = spart*x1*x3
       s15 = spart*x2*x4
       z5 = one/(one - x2) 

       C5S4Lim%Lim_KinInv(1:5) = (/s12,s14,s24,s15,z5/)

       pVV = p1 + p2 - p5
       call get_lept_mom(MV,pVV,yr(7:9),nlept,kallenF,C5S4Lim)

       C5S4Lim%AmpMom(:,1) = p1 - p5 
       C5S4Lim%AmpMom(:,2) = p2
       !C5S4Lim%AmpMom(:,3) = C5S4Lim%AmpMom(:, 4)
       !C5S4Lim%AmpMom(:,4) = C5S4Lim%AmpMom(:, 5)

       eta56 = half*(one-sc3(n5,nlept1))
       eta57 = half*(one-sc3(n5,nlept2))

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1, weight)

       C5S4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2


       !------ C5 C4 S4

       eta41 = one 

       z4 = x1
       s14 = spart*x1
       s24 = spart*x1*x3
       s15 = spart*x2*x4
       z5 = one/(one - x2) 

       C4C5S4Lim%Lim_KinInv(1:5) = (/z4,s14,s24,s15,z5/)

       !call sector_nnlo_4151(eta41,eta42,eta51,eta52,eta56,eta57,weight)
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,1, weight)


       C4C5S4Lim%wgt = one/8.0_dp/pi*kallenF * (spart/two)**2 & 
            * x1*x2 & 
            * weight &
            * one/two/spart & 
            * spart/MV2



    endif






  end subroutine kinematics_nnlo_5261
  

  
  subroutine kinematics_nnlo_5163(MV2,Y,yr,&
       HardProc,S4Lim,S5Lim,S4S5Lim,C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,&
       C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim)
    real(dp), intent(in) ::   MV2, Y,yr(10)
    type(KinConfig)      :: HardProc,S4Lim,S4S5Lim,S5Lim,C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,&
         C5Lim,C5S4Lim,C5S4S5Lim,C5S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim
    real(dp) :: MV, kallenF, xi1, xi2, spart, sqrts
    real(dp) :: x1, x2,x3, x4, lambda
    real(dp) :: cos4, sin4, phi4, cos5, sin5, phi5, phi6, cos6, sin6, cosphi6, sinphi6
    real(dp) :: n4(4), n6aux(3), n6t1(3), n6t2(3), n6(4), n5aux(3), n5(4), n1(4), n2(4), n4c(4), n5c(4), n7(4)
    real(dp) :: eta41, eta42, eta56, eta51, eta52, C4eta41, C4eta42, C5eta56, C5eta51, C5eta52, eta57
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4), p6(4), p7(4), Q(4), E5, E6, Qsq
    real(dp) :: s14, s15, s12, s24, s65, z4, z5, sector 
    integer :: i

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S4Lim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(S4S5Lim)        !4
    call initialize_config(C4Lim)          !5
    call initialize_config(C4S4Lim)        !6
    call initialize_config(C4S5Lim)        !7
    call initialize_config(C4S4S5Lim)      !8
    call initialize_config(C5Lim)          !9
    call initialize_config(C5S4Lim)        !10
    call initialize_config(C5S5Lim)        !11
    call initialize_config(C5S4S5Lim)      !12
    call initialize_config(C4C5Lim)        !13
    call initialize_config(C4C5S4Lim)      !14
    call initialize_config(C4C5S5Lim)      !15
    call initialize_config(C4C5S4S5Lim)    !16

    HardProc%npart = 6    
    S4Lim%npart = 5           
    S4S5Lim%npart = 4         
    S5Lim%npart = 5           
    C4Lim%npart = 5           
    C4S4Lim%npart = 5         
    C4S4S5Lim%npart = 4       
    C4S5Lim%npart =  4        
    C5Lim%npart =  5          
    C5S4Lim%npart = 4         
    C5S4S5Lim%npart = 4       
    C5S5Lim%npart = 5         
    C4C5Lim%npart = 4         
    C4C5S4Lim%npart = 4       
    C4C5S5Lim%npart = 4       
    C4C5S4S5Lim%npart = 4    

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a]     
    S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]            
    S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]            
    C4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]            
    C4S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]          
    C4S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4S5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]            
    C5S4Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C5S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]          
    C4C5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C4C5S4Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4C5S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4C5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]     
    
    !    MV2 and Y are generated

    MV =  sqrt(MV2)

    x1 =   yr(1)  !energy of the gluon
    x2 =   yr(2)  !energy of the photon
    x3 =   yr(3)  !angle of the gluon
    x4 =   yr(5)  !angle of the photon

   ! x1 = 1d-6
   ! x2 = 1d-6
   ! x3 = 1d-12
   ! x4 = 1d-12

    !gluon
    cos4 = one - two*x3         
    sin4 = sqrt(one-cos4**2)
    phi4 = twopi*yr(4)           

    n4 = (/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)

    !photon: coordinates
    cos5 = one-two*x4
    sin5 = sqrt(one-cos5**2)
    phi5 = twopi*yr(6)

    !lepton: direction
    phi6 = twopi*yr(7)
    cos6 = one -  two*yr(8) 
    sin6 = sqrt(one-cos6**2)
    cosphi6 = cos(phi6)
    sinphi6 = sin(phi6)

    !randomize 1 axis
    n6aux = (/sin6*cosphi6, sin6*sinphi6, cos6/)
    n6t1 = (/-cos6*cosphi6,-cos6*sinphi6,sin6/) 
    n6t2 = (/-sinphi6,cosphi6, zero/)

    n6 = (/one,n6aux(1),n6aux(2),n6aux(3)/)

    !photon: direction
    n5aux = cos5*n6aux+sin5*(cos(phi5)*n6t1+sin(phi5)*n6t2)

    n5 = (/one, n5aux(1), n5aux(2), n5aux(3)/)

    !quarks
    n1 = (/one,zero,zero,one/)
    n2 = (/one,zero,zero,-one/)

    !eta variables 
    eta41 = x3
    eta42 = one - x3
    eta56 = x4
    eta51 = half*scr(n5,n1)
    eta52 = half*scr(n5,n2)

    !limits
    n4c = n1
    n5c = n6

    C4eta41 = zero
    C4eta42 = one
    
    C5eta56 = zero
    C5eta51 = half*scr(n5c,n1)
    C5eta52 = half*scr(n5c,n2)

!***************************
!      I.FLM(1,2,6,7)
!***************************

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,HardProc%flag)

    if (HardProc%flag .eqv. .false.) then 

       HardProc%ids = [0,0,id_el,-id_el,id_g,id_a]
       HardProc%npart = 6
       
       !-- q[1] qb[2] -> XX g[4] a[5] l[6] lb[7]
       !-- n5c = n6
       !-- v -> hard1 hard2 g1 [g can be coll to 1]
       !-- qt = pv - pg1
       !-- so pV = p6+p7+pa = p1 + p2 - p4
       !--    Qt  = pv - p5
       !-- chiara calls pv ``q''
       !-- E6 = scr(Qt,Qt)/two/(qt.n6)
       !--    = scr(Q-p5,Q-p5)/2/([Q-p5].n6)
       !--    = half*(q2 - 2 q.p5)/(q.n6-p5.n6)
       
       p1 = half*sqrts*n1
       p2 = half*sqrts*n2
       p4 = half*sqrts*x1*n4
       p5 = half*sqrts*x2*n5

       Q = p1 + p2 - p4  !momentum of the intermediate boson 
       Qsq = scr(Q,Q)

       ! generate the decay of the boson in the full kinematics     

       E6 = half*(two*scr(Q,p5)-Qsq)/(scr(n6,p5)-scr(n6,Q))
       p6 = E6*n6

       p7 = Q - p6 - p5
       n7 = p7(1:4)/p7(1)

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  HardProc%flag = .true. 
    if (HardProc%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,n6)-scr(p5,n6))

       eta57 = half*scr(n5,n7)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,6, sector)



       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p6
       HardProc%AmpMom(:,4) = p7
       HardProc%AmpMom(:,5) = p4
       HardProc%AmpMom(:,6) = p5


       HardProc%PartFrac=(/xi1,xi2/)

       HardProc%wgt = two/pi*kallenF &
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) & ! -- flux
            *sector &
            *spart/MV2 

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: S4 S5 + collinear all angle identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! -------------- S4 S5

    call get_flag_nnlo_interf(sh,MV2,Y,zero,eta41,eta42,spart,sqrts,xi1,xi2,S4S5Lim%flag)

    if (S4S5Lim%flag .eqv. .false.) then

       S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]
       S4S5Lim%npart = 4

       
    p1 = half*sqrts*n1
    p2 = half*sqrts*n2
    p4 = half*sqrts*x1*n4
    p5 = half*sqrts*x2*n5

    Q = p1 + p2  
    Qsq = scr(Q,Q)

    E6 = half*Qsq/scr(n6,Q)
    p6 = E6*n6

    p7 = Q - p6
    n7 = p7(1:4)/p7(1)

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  S4S5Lim%flag = .true. 
    if (S4S5Lim%flag .eqv. .false.) then

       kallenF = E6 / scr(n6,Q)

       eta57 = half*scr(n5,n7)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,6, sector)


       S4S5Lim%AmpMom(:,1) = p1
       S4S5Lim%AmpMom(:,2) = p2
       S4S5Lim%AmpMom(:,3) = p6
       S4S5Lim%AmpMom(:,4) = p7

       S4S5Lim%LimMom(:,1:4) = S4S5Lim%AmpMom(:,1:4)
       S4S5Lim%LimMom(:,5) = p4                        
       S4S5Lim%LimMom(:,6) = p5

       S4S5Lim%PartFrac=(/xi1,xi2/)

       S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            * spart/MV2

       ! ------------   C4  +  S4 + S5  -----------

       p4 = half*sqrts*x1*n4c

       s14 = spart*x1*x3 
       z4 = x1 

       C4S4S5Lim%Lim_KinInv(1:5) = (/zero,zero,zero,s14,z4/)

       C4S4S5Lim%LimMom(:,1) = p1 - p4
       C4S4S5Lim%LimMom(:,2) = p2
       C4S4S5Lim%LimMom(:,3) = p6
       C4S4S5Lim%LimMom(:,4) = p7                           
       C4S4S5Lim%LimMom(:,5) = p5

       call sector_nnlo_minimal(C4eta41,C4eta42,eta51,eta52,eta56,eta57,1,6, sector)


       C4S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector&
            *spart/MV2

       ! ------------  C5 S4S5 FLM  ------------

       eta57 = half*scr(n5c,n7)

       s12 = two*scr(p1,p2) 
       s14 = spart*x1*x3
       s24 = spart*x1*(one - x3) 

       E5 = half*sqrts*x2
       s65 = four*E6*E5*x4
       z5 = E5/E6

       C5S4S5Lim%Lim_KinInv(1:5) = (/s12,s14,s24,s65,z5/)
    
       call sector_nnlo_minimal(eta41,eta42,C5eta51,C5eta52,C5eta56,eta57,1,6, sector)

       C5S4S5Lim%wgt = two/pi*kallenF &
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector&
            *spart/MV2
       

       ! -----------   C4  + C5 + S4S5  ----------

       s14 = spart *x1*x3 
       z4 = x1 
       s65 = four*E6*E5*x4
       z5 = E5/E6

       C4C5S4S5Lim%Lim_KinInv(1:5) = (/s14,s65,zero,z4,z5/)

       call sector_nnlo_minimal(C4eta41,C4eta42,C5eta51,C5eta52,C5eta56,eta57,1,6, sector)

       C4C5S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            *spart/MV2

    endif

    ! --------------- SOFT 4 ----------
    ! gluon soft

    call get_flag_nnlo_interf(sh,MV2,Y,zero,eta41,eta42,spart,sqrts,xi1,xi2,S4Lim%flag)

    if (S4Lim%flag .eqv. .false.) then

       S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]
       S4Lim%npart = 5

       
    p1 = half*sqrts*n1
    p2 = half*sqrts*n2
    p4 = half*sqrts*x1*n4
    p5 = half*sqrts*x2*n5

    Q = p1 + p2  
    Qsq = scr(Q,Q)

    E6 = half*(two*scr(Q,p5)-Qsq)/(scr(n6,p5)-scr(n6,Q))
    p6 = E6*n6

    p7 = Q - p6 - p5
    n7 = p7(1:4)/p7(1)

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  S4Lim%flag = .true. 
    if (S4Lim%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,n6)-scr(p5,n6))
       
       eta57 = half*scr(n5,n7)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,6, sector)

       S4Lim%PartFrac = (/xi1,xi2/)

       S4Lim%AmpMom(:,1) = p1
       S4Lim%AmpMom(:,2) = p2
       S4Lim%AmpMom(:,3) = p6
       S4Lim%AmpMom(:,4) = p7
       S4Lim%AmpMom(:,5) = p5 

       S4Lim%LimMom(:,1) = p1
       S4Lim%LimMom(:,2) = p2
       S4Lim%LimMom(:,3) = p4

       S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)&
            /(two*spart) * sector &
            *spart/MV2

       ! ---------  C4 + S4 --------------

       s14 = spart*x1*x3
       z4 = x1

       C4S4Lim%Lim_KinInv(1:5) = (/zero,s14,zero,zero,z4/)

       call sector_nnlo_minimal(C4eta41,C4eta42,eta51,eta52,eta56,eta57,1,6, sector)

       C4S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)&
            /(two*spart) * sector &
            *spart/MV2

    endif

    ! --------------- SOFT 5 ----------

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,S5Lim%flag)

    if (S5Lim%flag .eqv. .false.) then

       S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]
       S5Lim%npart = 5

       
    p1 = half*sqrts*n1
    p2 = half*sqrts*n2
    p4 = half*sqrts*x1*n4
    p5 = half*sqrts*x2*n5

    Q = p1 + p2 - p4  !momentum of the intermediate boson 
    Qsq = scr(Q,Q)

    E6 = half*Qsq/scr(n6,Q)
    p6 = E6*n6

    p7 = Q - p6
    n7 = p7(1:4)/p7(1)

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero) S5Lim%flag = .true.
    if (S5Lim%flag .eqv. .false.) then

       kallenF = E6 / scr(Q,n6)

       eta57 = half*scr(n5,n7)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,6, sector)

       S5Lim%PartFrac = (/xi1, xi2/)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p6
       S5Lim%AmpMom(:,4) = p7
       S5Lim%AmpMom(:,5) = p4 

       S5Lim%LimMom(:,1:4) = S5Lim%AmpMom(:,1:4)
       S5Lim%LimMom(:,5) = p5

       S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            * spart/MV2

       ! ---------  C5  + S5 --------------

       eta57 = half* (one-sc3(n5,n7))

       E5 = half* sqrts*x2      
       s65 = four*E6*E5*x4
       z5 = E5/E6

       C5S5Lim%Lim_KinInv(1:5) = (/zero,s65,zero,zero,z5/)

       call sector_nnlo_minimal(eta41,eta42,C5eta51,C5eta52,C5eta56,eta57,1,6, sector)

       C5S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            *spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3:   C4
    !-------------------------------------------------------------------------------------------------------------

    call get_flag_nnlo_interf(sh,MV2,Y,x1,C4eta41,C4eta42,spart,sqrts,xi1,xi2,C4Lim%flag)

    if (C4Lim%flag .eqv. .false.) then 

       C4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]
       C4Lim%npart = 5

       
    p1 = half*sqrts*n1
    p2 = half*sqrts*n2
    p4 = half*sqrts*x1*n4c
    p5 = half*sqrts*x2*n5

    Q = p1 + p2 - p4
    Qsq = scr(Q,Q)

    E6 = half*(two*scr(Q,p5)-Qsq)/(scr(n6,p5)-scr(n6,Q))
    p6 = E6*n6

    p7 = Q - p6 - p5
    n7 = p7(1:4)/p7(1)

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C4Lim%flag = .true. 
    if (C4Lim%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,n6)-scr(p5,n6))

       eta57 = half*scr(n5,n7)

       call sector_nnlo_minimal(C4eta41,C4eta42,eta51,eta52,eta56,eta57,1,6, sector)

       C4Lim%PartFrac=(/xi1,xi2/)

       C4Lim%AmpMom(:,1) = p1 - p4 
       C4Lim%AmpMom(:,2) = p2
       C4Lim%AmpMom(:,3) = p6
       C4Lim%AmpMom(:,4) = p7
       C4Lim%AmpMom(:,5) = p5

       s14 = spart*x1*x3 
       z4 = one/(one - x1)

       C4Lim%Lim_KinInv(1:5) = (/s14,zero,zero, z4, zero/)

       C4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C41  + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    call get_flag_nnlo_interf(sh,MV2,Y,x1,C4eta41,C4eta42,spart,sqrts,xi1,xi2,C4S5Lim%flag)


    if (C4S5Lim%flag .eqv. .false.) then 

       C4S5Lim%ids(1:4) = [0,0,id_el,-id_el]
       C4S5Lim%npart = 4

       
    p1 = half*sqrts*n1
    p2 = half*sqrts*n2
    p4 = half*sqrts*x1*n4c
    p5 = half*sqrts*x2*n5

    Q = p1 + p2 - p4
    Qsq = scr(Q,Q)

    E6 = half*Qsq/scr(n6,Q)
    p6 = E6*n6

    p7 = Q - p6
    n7 = p7(1:4)/p7(1)

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C4S5Lim%flag = .true. 
    if (C4S5Lim%flag .eqv. .false.) then

       kallenF = E6 / scr(Q,n6)

       eta57 = half*scr(n5,n7)

       call sector_nnlo_minimal(C4eta41,C4eta42,eta51,eta52,eta56,eta57,1,6, sector)

       C4S5Lim%PartFrac=(/xi1,xi2/)

       C4S5Lim%AmpMom(:,1) = p1 - p4 
       C4S5Lim%AmpMom(:,2) = p2
       C4S5Lim%AmpMom(:,3) = p6
       C4S5Lim%AmpMom(:,4) = p7

       C4S5Lim%LimMom(:,1:4) = C4S5Lim%AmpMom(:,1:4)
       C4S5Lim%LimMom(:,5) = p5

       s14 = spart*x1*x3 
       z4 = one/(one - x1)

       C4S5Lim%Lim_KinInv(1:5) = (/zero,zero,zero,s14,z4/)

       C4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2


       ! ---------   C4   + C5 + SOFT5 ------------

       E5 = half*sqrts*x2
       z5 =  E5/E6

       s15 = spart *x2*(one-x4) 
       s14 = spart *x1*x3 
       s65 = four*E6*E5*x4

       C4C5S5Lim%Lim_KinInv(1:5) = (/s14,s65,zero,z4,z5/)
       
       eta57 = half* (one-sc3(n5c,n7))

       call sector_nnlo_minimal(C4eta41,C4eta42,C5eta51,C5eta52,C5eta56,eta57,1,6, sector)

       C4C5S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: C5  FLM 
    !-------------------------------------------------------------------------------------------------------------

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,C5Lim%flag)

    if (C5Lim%flag .eqv. .false.) then 

       C5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]
       C5Lim%npart = 5

       
    p1 = half*sqrts*n1
    p2 = half*sqrts*n2
    p4 = half*sqrts*x1*n4
    p5 = half*sqrts*x2*n5c

    Q = p1 + p2 - p4   
    Qsq = scr(Q,Q)

    E6 = half*(two*scr(Q,p5)-Qsq)/(scr(n6,p5)-scr(n6,Q))
    p6 = E6*n6

    p7 = Q - p6 - p5
    n7 = p7(1:4)/p7(1)

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C5Lim%flag = .true. 
    if (C5Lim%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,n6)-scr(p5,n6))

       eta57 = half*scr(n5,n7)
       
       call sector_nnlo_minimal(eta41,eta42,C5eta51,C5eta52,C5eta56,eta57,1,6, sector)

       C5Lim%PartFrac = (/xi1,xi2/)

       C5Lim%AmpMom(:,1) = p1
       C5Lim%AmpMom(:,2) = p2 
       C5Lim%AmpMom(:,3) = p6 + p5
       C5Lim%AmpMom(:,4) = p7
       C5Lim%AmpMom(:,5) = p4

       !-- photon 5, lepton 6, antilepton 7
       !-- coll lim: p5.p6 = e5*rho56*e6[x2=0]

       E5 = half*sqrts*x2
       s65 = four*E6*E5*x4
       z5 = E6/(E6 + E5) 

       
       C5Lim%Lim_KinInv(1:5) = (/zero,s65,zero,zero,z5/)
       C5Lim%PartFrac=(/xi1,xi2/)

       C5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C41  + C52
    !-------------------------------------------------------------------------------------------------------------

    call get_flag_nnlo_interf(sh,MV2,Y,x1,C4eta41,C4eta42,spart,sqrts,xi1,xi2,C4C5Lim%flag)

    if (C4C5Lim%flag .eqv. .false.) then 

       C4C5Lim%ids(1:4) = [0,0,id_el,-id_el]
       C4C5Lim%npart = 4

       
    p1 = half*sqrts*n1
    p2 = half*sqrts*n2
    p4 = half*sqrts*x1*n4c
    p5 = half*sqrts*x2*n5c

    Q = p1 + p2 - p4  
    Qsq = scr(Q,Q)

    E6 = half*(two*scr(Q,p5)-Qsq)/(scr(n6,p5)-scr(n6,Q))
    p6 = E6*n6

    p7 = Q - p6 - p5
    n7 = p7(1:4)/p7(1)

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C4C5Lim%flag = .true. 
    if (C4C5Lim%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,n6)-scr(p5,n6))

       eta57 = half*scr(n5c,n7)

       call sector_nnlo_minimal(C4eta41,C4eta42,C5eta51,C5eta52,C5eta56,eta57,1,6, sector)

       C4C5Lim%PartFrac = (/xi1,xi2/)

       C4C5Lim%AmpMom(:,1) = p1 - p4
       C4C5Lim%AmpMom(:,2) = p2 
       C4C5Lim%AmpMom(:,3) = p6 + p5
       C4C5Lim%AmpMom(:,4) = p7


       E5 = half*sqrts*x2
       s14 = spart *x1*x3 
       s65 = four*E6*E5*x4
       z4 = one/(one - x1)
       z5 = E6/(E6 + E5) 

       C4C5Lim%Lim_KinInv(1:5) = (/s14,s65,zero,z4,z5/)

       C4C5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C5S4  
    !-------------------------------------------------------------------------------------------------------------

    call get_flag_nnlo_interf(sh,MV2,Y,zero,eta41,eta42,spart,sqrts,xi1,xi2,C5S4Lim%flag)

    if (C5S4Lim%flag .eqv. .false.) then 

       C5S4Lim%ids(1:4) = [0,0,id_el,-id_el]
       C5S4Lim%npart = 4

       
    p1 = half*sqrts*n1
    p2 = half*sqrts*n2
    p4 = half*sqrts*x1*n4
    p5 = half*sqrts*x2*n5c

    Q = p1 + p2
    Qsq = scr(Q,Q)

    E6 = half*(two*scr(Q,p5)-Qsq)/(scr(n6,p5)-scr(n6,Q))
    p6 = E6*n6

    p7 = Q - p6 - p5
    n7 = p7(1:4)/p7(1)

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C5S4Lim%flag = .true. 
    if (C5S4Lim%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,n6)-scr(p5,n6))

       eta57 = half*scr(n5c,n7)

       call sector_nnlo_minimal(eta41,eta42,C5eta51,C5eta52,C5eta56,eta57,1,6, sector)

       C5S4Lim%PartFrac = (/xi1,xi2/)

       C5S4Lim%AmpMom(:,1) = p1
       C5S4Lim%AmpMom(:,2) = p2 
       C5S4Lim%AmpMom(:,3) = p6 + p5
       C5S4Lim%AmpMom(:,4) = p7


       E5 = half*sqrts*x2
       s14 = spart *x1*x3
       s24 = spart *x1*(one-x3) 
       s65 = four*E6*E5*x4
       s12 = two*scr(p1,p2) 

       z5 = E6/(E6 + E5) 

       C5S4Lim%Lim_KinInv(1:5) = (/s12,s14,s24,s65,z5/)

       C5S4Lim%PartFrac=(/xi1,xi2/)

       C5S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

       ! ---------   C4   + C5 + S4  ------------

       call sector_nnlo_minimal(C4eta41,C4eta42,C5eta51,C5eta52,C5eta56,eta57,1,6, sector)

       z4 = x1 
       z5 =  E6/(E6 + E5)

       s14 = spart *x1*x3 
       s65 = four*E5*E6*x4

       C4C5S4Lim%Lim_KinInv(1:5) = (/s14,s65,zero,z4,z5/)

     
       C4C5S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

    endif


  end subroutine kinematics_nnlo_5163

  


   subroutine kinematics_nnlo_5164(MV2,Y,yr,&
       HardProc,S4Lim,S5Lim,S4S5Lim,C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,&
       C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim)
    real(dp), intent(in) ::   MV2, Y,yr(10)
    type(KinConfig)      :: HardProc,S4Lim,S4S5Lim,S5Lim,C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,&
         C5Lim,C5S4Lim,C5S4S5Lim,C5S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: x1, x2,x3, x4, cos4, sin4, phi4, phi5, xi1, xi2, phi7, cos7, sin7
    real(dp) :: eta41, eta42, eta51, eta52
    real(dp) :: s14, s15, z4, z5
    real(dp) :: spart,sqrts, s12, s75,s24
    real(dp) :: n6(3), n7(3), n5(3), p6(4),p7(4),Q(4), Qsq, E5, eta56, eta57
    real(dp) :: n4(3), n3(3), n1(3), n2(3)
    real(dp) :: cosphi7, sinphi7, n7t1(3), n7t2(3)
    real(dp) :: n1aux,n2aux, sector, s65
    real(dp) :: p1aux(4), p2aux(4), p4aux(4), p4auxc(4), p4c(4), p5aux(4), p5auxc(4), p6aux(4), p7aux(4)
    integer :: i

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S4Lim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(S4S5Lim)        !4
    call initialize_config(C4Lim)          !5
    call initialize_config(C4S4Lim)        !6
    call initialize_config(C4S5Lim)        !7
    call initialize_config(C4S4S5Lim)      !8
    call initialize_config(C5Lim)          !9
    call initialize_config(C5S4Lim)        !10
    call initialize_config(C5S5Lim)        !11
    call initialize_config(C5S4S5Lim)      !12
    call initialize_config(C4C5Lim)        !13
    call initialize_config(C4C5S4Lim)      !14
    call initialize_config(C4C5S5Lim)      !15
    call initialize_config(C4C5S4S5Lim)    !16

    HardProc%npart = 6    
    S4Lim%npart = 5           
    S4S5Lim%npart = 4         
    S5Lim%npart = 5           
    C4Lim%npart = 5           
    C4S4Lim%npart = 5         
    C4S4S5Lim%npart = 4       
    C4S5Lim%npart =  4        
    C5Lim%npart =  5          
    C5S4Lim%npart = 4         
    C5S4S5Lim%npart = 4       
    C5S5Lim%npart = 5         
    C4C5Lim%npart = 4         
    C4C5S4Lim%npart = 4       
    C4C5S5Lim%npart = 4       
    C4C5S4S5Lim%npart = 4    

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a]     
    S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]            
    S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]            
    C4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]            
    C4S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]          
    C4S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4S5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]            
    C5S4Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C5S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]          
    C4C5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C4C5S4Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4C5S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4C5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]     
    
    
    !    MV2 and Y are generated

    MV =  sqrt(MV2)

    x1 =   yr(1)  !energy of the gluon
    x2 =   yr(2)  !energy of the photon
    x3 =   yr(3)  !angle of the gluon
    x4 =   yr(5)  !angle of the photon


  !  print*, (yr(i),i=1,10)


    !x1 = 1d-6
    !x2 = 1d-6
    !x3 = 1d-12
    !x4 = 1d-12


    !gluon
    cos4 = one - two*x3         
    sin4 = sqrt(one-cos4**2)
    phi4 = twopi*yr(4)           

    eta41 = x3
    eta42 = one - x3

    p4aux = (/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)
    p4auxC = (/ one, zero, zero, one /)
    
    !photon: coordinates
    cos5 = one-two*x4
    sin5 = sqrt(one-cos5**2)
    phi5 = twopi*yr(6)

    !lepton 7: direction
    phi7 = twopi*yr(7)
    cos7 = one -  two*yr(8) 
    sin7 = sqrt(one-cos7**2)
    cosphi7 = cos(phi7)
    sinphi7 = sin(phi7)

    !randomize 1 axis
    n7 = (/sin7*cosphi7, sin7*sinphi7, cos7/)
    n7t1 = (/-cos7*cosphi7,-cos7*sinphi7,sin7/) 
    n7t2 = (/-sinphi7,cosphi7, zero/)

    p7aux = (/one,n7(1),n7(2),n7(3)/)

    !photon: direction
    n5 = cos5*n7+sin5*(cos(phi5)*n7t1+sin(phi5)*n7t2)
    p5aux = (/one,n5(1),n5(2),n5(3)/)

    eta57 = x4

    p5auxC = p7aux

    !quarks
    n1 = (/zero,zero,one/)
    n2 = (/zero,zero,-one/)

    p1aux = (/one,zero,zero,one/)
    p2aux = (/one,zero,zero,-one/)

  

    !  full kinematics


    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,HardProc%flag)

    if (HardProc%flag .eqv. .false.) then 

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(two*scr(Q,p5)-Qsq)/(scr(p7aux,p5)-scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 - p5
       E6 = p6(1)
       p6aux = p6(1:4)/E6


    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  HardProc%flag = .true. 

    if (HardProc%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux)- scr(p5,p7aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta56 = half*scr(p5aux, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)


       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p6
       HardProc%AmpMom(:,4) = p7
       HardProc%AmpMom(:,5) = p4
       HardProc%AmpMom(:,6) = p5

       HardProc%PartFrac=(/xi1,xi2/)

       HardProc%wgt = kallenF/two/pi &
                     *(spart/two)*x2 & 
                     * (spart/two)*x1 &  
            * sector  &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: S4 S5 + collinear all angle identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! -------------- S4 S5

    call get_flag_nnlo_interf(sh,MV2,Y,zero,eta41,eta42,spart,sqrts,xi1,xi2,S4S5Lim%flag)

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2   !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(Qsq)/(scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7
       E6 = p6(1)
       p6aux = p6(1:4)/E6

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  S4S5Lim%flag = .true. 

    if (S4S5Lim%flag .eqv. .false.) then

        kallenF = E7 / (scr(Q,p7aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta56 = half*scr(p5aux, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)

       S4S5Lim%AmpMom(:,1) = p1
       S4S5Lim%AmpMom(:,2) = p2
       S4S5Lim%AmpMom(:,3) = p6
       S4S5Lim%AmpMom(:,4) = p7

       S4S5Lim%LimMom(:,1:4) = S4S5Lim%AmpMom(:,1:4)
       S4S5Lim%LimMom(:,5) = p4                        
       S4S5Lim%LimMom(:,6) = p5

       S4S5Lim%PartFrac=(/xi1,xi2/)

       S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            * spart/MV2

       ! ------------   C4  +  S4 + S5  -----------

       eta41 = zero
       eta42 = one
       eta57 = x4

       p4c = half*sqrts*x1*p4auxC

       s14 = spart*x1*x3 
       z4 = x1 

       C4S4S5Lim%Lim_KinInv(1:5) = (/zero,zero,zero,s14,z4/)

       C4S4S5Lim%LimMom(:,1) = p1 - p4c
       C4S4S5Lim%LimMom(:,2) = p2
       C4S4S5Lim%LimMom(:,3) = p6
       C4S4S5Lim%LimMom(:,4) = p7                           
       C4S4S5Lim%LimMom(:,5) = p5

      
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)

       C4S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector&
            *spart/MV2

       ! ------------  C5 S4S5 FLM  ------------
       
       eta41 = x3
       eta42 = one - x3
       eta57 = zero

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta56 = half*scr(p5auxC, p6aux)
    

    s12 = two*scr(p1,p2) 
    s14 = spart*x1*x3
    s24 = spart*x1*(one - x3) 

    s75 = four*E7*E5*x4
    z5 = E5/E7

    C5S4S5Lim%Lim_KinInv(1:5) = (/s12,s14,s24,s75,z5/)
    
     call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)

       C5S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector&
            *spart/MV2
       

       ! -----------   C4  + C5 + S4S5  ----------

       eta41 = zero
       eta42 = one
       eta57 = zero

       s14 = spart *x1*x3 
       z4 = x1 
       s75 = four*E7*E5*x4
       z5 = E5/E7

       C4C5S4S5Lim%Lim_KinInv(1:5) = (/s14,s75,zero,z4,z5/)
    
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)

       C4C5S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            *spart/MV2

    endif

    ! --------------- SOFT 4 ----------
    ! gluon soft

    eta41 = x3
    eta42 = one - x3
    eta57 = x4
   
    call get_flag_nnlo_interf(sh,MV2,Y,zero,eta41,eta42,spart,sqrts,xi1,xi2,S4Lim%flag)

    if (S4Lim%flag .eqv. .false.) then

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(two*scr(Q,p5)-Qsq)/(scr(p7aux,p5)-scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 - p5
       E6 = p6(1)
       p6aux = p6(1:4)/E6


    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  S4Lim%flag = .true. 

    if (S4Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux)- scr(p5,p7aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta56 = half*scr(p5aux, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)


       S4Lim%PartFrac = (/xi1,xi2/)


       S4Lim%AmpMom(:,1) = p1
       S4Lim%AmpMom(:,2) = p2
       S4Lim%AmpMom(:,3) = p6
       S4Lim%AmpMom(:,4) = p7
       S4Lim%AmpMom(:,5) = p5 


       S4Lim%LimMom(:,1) = p1
       S4Lim%LimMom(:,2) = p2
       S4Lim%LimMom(:,3) = p4

       S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)&
            /(two*spart) * sector &
            *spart/MV2

       ! ---------  C4 + S4 --------------

       eta41 = zero
       eta42 = one

       s14 = spart*x1*x3
       z4 = x1

       C4S4Lim%Lim_KinInv(1:5) = (/zero,s14,zero,zero,z4/)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)

       C4S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)&
            /(two*spart) * sector &
            *spart/MV2


    endif

    ! --------------- SOFT 5 ----------

    eta41 = x3
    eta42 = one - x3
    eta56 = x4

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,S5Lim%flag)


    if (S5Lim%flag .eqv. .false.) then

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(Qsq)/(scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 
       E6 = p6(1)
       p6aux = p6(1:4)/E6

       if ( p7(1).lt.zero.or.p6(1).lt.zero) S5Lim%flag = .true.

    endif

    if (S5Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta56 = half*scr(p5aux, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)


       S5Lim%PartFrac = (/xi1, xi2/)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p6
       S5Lim%AmpMom(:,4) = p7
       S5Lim%AmpMom(:,5) = p4 

       S5Lim%LimMom(:,1:4) = S5Lim%AmpMom(:,1:4)
       S5Lim%LimMom(:,5) = p5

       S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            * spart/MV2


       ! ---------  C5  + S5 --------------

       eta57 = zero

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta56 = half*scr(p5auxC, p6aux)

       s75 = four*E7*E5*x4
       z5 = E5/E7

       C5S5Lim%Lim_KinInv(1:5) = (/zero,s75,zero,zero,z5/)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)

       C5S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            *spart/MV2


    endif

    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3:   C4
    !-------------------------------------------------------------------------------------------------------------

    eta41 = zero
    eta42 = one
    eta57 = x4
    
    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,C4Lim%flag)


    if (C4Lim%flag .eqv. .false.) then

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4auxC
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(two*scr(Q,p5)-Qsq)/(scr(p7aux,p5)-scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 - p5
       E6 = p6(1)
       p6aux = p6(1:4)/E6

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C4Lim%flag = .true. 

    if (C4Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux)- scr(p5,p7aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta56 = half*scr(p5aux, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)

       C4Lim%PartFrac=(/xi1,xi2/)

       C4Lim%AmpMom(:,1) = p1 - p4 
       C4Lim%AmpMom(:,2) = p2
       C4Lim%AmpMom(:,3) = p6
       C4Lim%AmpMom(:,4) = p7
       C4Lim%AmpMom(:,5) = p5

       s14 = spart*x1*x3 
       z4 = one/(one - x1)

       C4Lim%Lim_KinInv(1:5) = (/s14,zero,zero, z4, zero/)


       C4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2


       
    endif


    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C41  + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    eta41 = zero
    eta42 = one

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,C4S5Lim%flag)


    if (C4S5Lim%flag .eqv. .false.) then 

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4auxC
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(Qsq)/(scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 
       E6 = p6(1)
       p6aux = p6(1:4)/E6

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C4S5Lim%flag = .true. 

    if (C4S5Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta56 = half*scr(p5aux, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)

       C4S5Lim%PartFrac=(/xi1,xi2/)

       C4S5Lim%AmpMom(:,1) = p1 - p4 
       C4S5Lim%AmpMom(:,2) = p2
       C4S5Lim%AmpMom(:,3) = p6
       C4S5Lim%AmpMom(:,4) = p7

       C4S5Lim%LimMom(:,1:4) = C4S5Lim%AmpMom(:,1:4)
       C4S5Lim%LimMom(:,5) = p5

       s14 = spart*x1*x3 
       z4 = one/(one - x1)

       C4S5Lim%Lim_KinInv(1:5) = (/zero,zero,zero,s14,z4/)


       C4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

       ! ---------   C4   + C5 + SOFT5 ------------

       eta57 = zero

       z5 =  E5/E7

       s15 = spart *x2*(one-x4) 
       s14 = spart *x1*x3 
       s75 = four*E7*E5*x4

       C4C5S5Lim%Lim_KinInv(1:5) = (/s14,s75,zero,z4,z5/)

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta56 = half*scr(p5auxC, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)

       C4C5S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2


    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: C5  FLM 
    !-------------------------------------------------------------------------------------------------------------

    eta41 = x3
    eta42 = one - x3
    eta57 = zero

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,C5Lim%flag)


    if (C5Lim%flag .eqv. .false.) then 

       p1= half*sqrts*p1aux
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 = E5*p5auxC

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(two*scr(Q,p5)-Qsq)/(scr(p7aux,p5)-scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 - p5
       E6 = p6(1)
       p6aux = p6(1:4)/E6
       

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C5Lim%flag = .true. 

    if (C5Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux)- scr(p5,p7aux))

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta56 = half*scr(p5auxC, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)

       C5Lim%PartFrac = (/xi1,xi2/)

       C5Lim%AmpMom(:,1) = p1
       C5Lim%AmpMom(:,2) = p2 
       C5Lim%AmpMom(:,3) = p6 
       C5Lim%AmpMom(:,4) = p7 + p5
       C5Lim%AmpMom(:,5) = p4


       s75 = four*E7*E5*x4
       z5 = E7/(E7 + E5) 

       
       C5Lim%Lim_KinInv(1:5) = (/zero,s75,zero,zero,z5/)
       HardProc%PartFrac=(/xi1,xi2/)


       C5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C41  + C52
    !-------------------------------------------------------------------------------------------------------------

    eta41 = zero
    eta42 = one 
    eta57 = zero

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,C4C5Lim%flag)


    if (C4C5Lim%flag .eqv. .false.) then
       
       p1= half*sqrts*p1aux
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4auxC
       E5 = half*sqrts*x2
       p5 =  E5*p5auxC

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(two*scr(Q,p5)-Qsq)/(scr(p7aux,p5)-scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 - p5
       E6 = p6(1)
       p6aux = p6(1:4)/E6

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C4C5Lim%flag = .true. 

    if (C4C5Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux)- scr(p5,p7aux))

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta56 = half*scr(p5auxC, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)

       C4C5Lim%PartFrac = (/xi1,xi2/)

       C4C5Lim%AmpMom(:,1) = p1 - p4
       C4C5Lim%AmpMom(:,2) = p2 
       C4C5Lim%AmpMom(:,3) = p6 
       C4C5Lim%AmpMom(:,4) = p7 + p5


       s14 = spart *x1*x3 
       s65 = four*E7*E5*x4
       z4 = one/(one - x1)
       z5 = E7/(E7 + E5) 

       C4C5Lim%Lim_KinInv(1:5) = (/s14,s65,zero,z4,z5/)


       C4C5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            * spart/MV2

    
    endif


    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C5S4  
    !-------------------------------------------------------------------------------------------------------------

    eta41 = x3
    eta42 = one - x3
    eta57 = zero

    call get_flag_nnlo_interf(sh,MV2,Y,zero,eta41,eta42,spart,sqrts,xi1,xi2,C5S4Lim%flag)


    if (C5S4Lim%flag .eqv. .false.) then 


       p1= half*sqrts*p1aux
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5auxC

       Q = p1 + p2   !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(two*scr(Q,p5)-Qsq)/(scr(p7aux,p5)-scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 - p5
       E6 = p6(1)
       p6aux = p6(1:4)/E6

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C5S4Lim%flag = .true. 

    if (C5S4Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux)- scr(p5,p7aux))

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta56 = half*scr(p5auxC, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)

       C5S4Lim%PartFrac = (/xi1,xi2/)

       C5S4Lim%AmpMom(:,1) = p1
       C5S4Lim%AmpMom(:,2) = p2 
       C5S4Lim%AmpMom(:,3) = p6 
       C5S4Lim%AmpMom(:,4) = p7 + p5


       s14 = spart *x1*x3
       s24 = spart *x1*(one-x3) 
       s75 = four*E7*E5*x4
       s12 = two*scr(p1,p2) 

       z5 = E7/(E7 + E5) 

       C5S4Lim%Lim_KinInv(1:5) = (/s12,s14,s24,s75,z5/)

       C5S4Lim%PartFrac=(/xi1,xi2/)


       C5S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

       ! ---------   C4   + C5 + S4  ------------

       eta41 = zero
       eta42 = one 
       eta57 = zero
       
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,1,7, sector)

       z4 = x1 
       z5 =  E7/(E7 + E5)

       s14 = spart *x1*x3 
       s75 = four*E5*E7*x4

       C4C5S4Lim%Lim_KinInv(1:5) = (/s14,s75,zero,z4,z5/)

       C4C5S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2


    endif


  end subroutine kinematics_nnlo_5164







  subroutine kinematics_nnlo_5263(MV2,Y,yr,&
       HardProc,S4Lim,S5Lim,S4S5Lim,C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,&
       C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim)
    real(dp), intent(in) ::   MV2, Y,yr(10)
    type(KinConfig)      :: HardProc,S4Lim,S4S5Lim,S5Lim,C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,&
         C5Lim,C5S4Lim,C5S4S5Lim,C5S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: x1, x2,x3, x4, cos4, sin4, phi4, phi5, xi1, xi2, phi6, cos6, sin6
    real(dp) :: eta41, eta42, eta51, eta52
    real(dp) :: s14, s25, z4, z5
    real(dp) :: spart,sqrts, s12, s75,s24
    real(dp) :: n6(3), n7(3), n5(3), p6(4),p7(4),Q(4), Qsq, E5, eta56, eta57
    real(dp) :: n4(3), n3(3), n1(3), n2(3)
    real(dp) :: cosphi6, sinphi6, n6t1(3), n6t2(3)
    real(dp) :: n1aux,n2aux, sector, s65
    real(dp) :: p1aux(4), p2aux(4), p4aux(4), p4auxc(4), p4c(4), p5aux(4), p5auxc(4), p6aux(4), p7aux(4)
    integer :: i

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S4Lim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(S4S5Lim)        !4
    call initialize_config(C4Lim)          !5
    call initialize_config(C4S4Lim)        !6
    call initialize_config(C4S5Lim)        !7
    call initialize_config(C4S4S5Lim)      !8
    call initialize_config(C5Lim)          !9
    call initialize_config(C5S4Lim)        !10
    call initialize_config(C5S5Lim)        !11
    call initialize_config(C5S4S5Lim)      !12
    call initialize_config(C4C5Lim)        !13
    call initialize_config(C4C5S4Lim)      !14
    call initialize_config(C4C5S5Lim)      !15
    call initialize_config(C4C5S4S5Lim)    !16

    HardProc%npart = 6    
    S4Lim%npart = 5           
    S4S5Lim%npart = 4         
    S5Lim%npart = 5           
    C4Lim%npart = 5           
    C4S4Lim%npart = 5         
    C4S4S5Lim%npart = 4       
    C4S5Lim%npart =  4        
    C5Lim%npart =  5          
    C5S4Lim%npart = 4         
    C5S4S5Lim%npart = 4       
    C5S5Lim%npart = 5         
    C4C5Lim%npart = 4         
    C4C5S4Lim%npart = 4       
    C4C5S5Lim%npart = 4       
    C4C5S4S5Lim%npart = 4    

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a]     
    S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]            
    S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]            
    C4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]            
    C4S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]          
    C4S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4S5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]            
    C5S4Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C5S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]          
    C4C5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C4C5S4Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4C5S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4C5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]     
    
    
    !    MV2 and Y are generated

    MV =  sqrt(MV2)

    x1 =   yr(1)  !energy of the gluon
    x2 =   yr(2)  !energy of the photon
    x3 =   yr(3)  !angle of the gluon: wrt quark 2
    x4 =   yr(5)  !angle of the photon: wrt lepton 6


  !  print*, (yr(i),i=1,10)


    !x1 = 1d-6
    !x2 = 1d-6
    !x3 = 1d-12
    !x4 = 1d-12


    !gluon
    cos4 = - one + two*x3         
    sin4 = sqrt(one-cos4**2)
    phi4 = twopi*yr(4)           

    eta42 = x3
    eta41 = one - x3

    p4aux = (/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)
    p4auxC = (/ one, zero, zero, - one /)
    
    !photon: coordinates
    cos5 = one-two*x4
    sin5 = sqrt(one-cos5**2)
    phi5 = twopi*yr(6)

    !lepton 6: direction
    phi6 = twopi*yr(7)
    cos6 = one -  two*yr(8) 
    sin6 = sqrt(one-cos6**2)
    cosphi6 = cos(phi6)
    sinphi6 = sin(phi6)

    !randomize 1 axis
    n6 = (/sin6*cosphi6, sin6*sinphi6, cos6/)
    n6t1 = (/-cos6*cosphi6,-cos6*sinphi6,sin6/) 
    n6t2 = (/-sinphi6,cosphi6, zero/)

    p6aux = (/one,n6(1),n6(2),n6(3)/)

    !photon: direction
    n5 = cos5*n6+sin5*(cos(phi5)*n6t1+sin(phi5)*n6t2)
    p5aux = (/one,n5(1),n5(2),n5(3)/)

    eta56 = x4

    p5auxC = p6aux

    !quarks
    n1 = (/zero,zero,one/)
    n2 = (/zero,zero,-one/)

    p1aux = (/one,zero,zero,one/)
    p2aux = (/one,zero,zero,-one/)

  

    !  full kinematics


    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,HardProc%flag)

    if (HardProc%flag .eqv. .false.) then 

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E6 = half*(two*scr(Q,p5)-Qsq)/(scr(p6aux,p5)-scr(p6aux,Q))

       p6 = E6*p6aux

       p7 = Q - p6 - p5
       E7 = p7(1)
       p7aux = p7(1:4)/E7


    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  HardProc%flag = .true. 

    if (HardProc%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,p6aux)- scr(p5,p6aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta57 = half*scr(p5aux, p7aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6,sector)


       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p6
       HardProc%AmpMom(:,4) = p7
       HardProc%AmpMom(:,5) = p4
       HardProc%AmpMom(:,6) = p5

       HardProc%PartFrac=(/xi1,xi2/)

       HardProc%wgt = kallenF/two/pi &
                     *(spart/two)*x2 & 
                     * (spart/two)*x1 &  
            * sector  &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: S4 S5 + collinear all angle identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! -------------- S4 S5

    call get_flag_nnlo_interf(sh,MV2,Y,zero,eta41,eta42,spart,sqrts,xi1,xi2,S4S5Lim%flag)

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2   !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E6 = half*(Qsq)/(scr(p6aux,Q))

       p6 = E6*p6aux

       p7 = Q - p6
       E7 = p7(1)
       p7aux = p7(1:4)/E7

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  S4S5Lim%flag = .true. 

    if (S4S5Lim%flag .eqv. .false.) then

        kallenF = E6 / (scr(Q,p6aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta57 = half*scr(p5aux, p7aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)

       S4S5Lim%AmpMom(:,1) = p1
       S4S5Lim%AmpMom(:,2) = p2
       S4S5Lim%AmpMom(:,3) = p6
       S4S5Lim%AmpMom(:,4) = p7

       S4S5Lim%LimMom(:,1:4) = S4S5Lim%AmpMom(:,1:4)
       S4S5Lim%LimMom(:,5) = p4                        
       S4S5Lim%LimMom(:,6) = p5

       S4S5Lim%PartFrac=(/xi1,xi2/)

       S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            * spart/MV2

       ! ------------   C4  +  S4 + S5  -----------

       eta42 = zero
       eta41 = one
       eta56 = x4

       p4c = half*sqrts*x1*p4auxC

       s24 = spart*x1*x3 
       z4 = x1 

       C4S4S5Lim%Lim_KinInv(1:5) = (/zero,zero,zero,s24,z4/)

       C4S4S5Lim%LimMom(:,1) = p1 
       C4S4S5Lim%LimMom(:,2) = p2 - p4c
       C4S4S5Lim%LimMom(:,3) = p6
       C4S4S5Lim%LimMom(:,4) = p7                           
       C4S4S5Lim%LimMom(:,5) = p5

      
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)

       C4S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector&
            *spart/MV2

       ! ------------  C5 S4S5 FLM  ------------
       
       eta42 = x3
       eta41 = one - x3
       eta56 = zero

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta57 = half*scr(p5auxC, p7aux)
    

    s12 = two*scr(p1,p2) 
    s24 = spart*x1*x3
    s14 = spart*x1*(one - x3) 

    s65 = four*E6*E5*x4
    z5 = E5/E6

    C5S4S5Lim%Lim_KinInv(1:5) = (/s12,s24,s14,s65,z5/)
    
     call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)

       C5S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector&
            *spart/MV2
       

       ! -----------   C4  + C5 + S4S5  ----------

       eta42 = zero
       eta41 = one
       eta56 = zero

       s24 = spart *x1*x3 
       z4 = x1 
       s65 = four*E6*E5*x4
       z5 = E5/E6

       C4C5S4S5Lim%Lim_KinInv(1:5) = (/s24,s65,zero,z4,z5/)
    
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)

       C4C5S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            *spart/MV2

    endif

    ! --------------- SOFT 4 ----------
    ! gluon soft

    eta42 = x3
    eta41 = one - x3
    eta56 = x4
   
    call get_flag_nnlo_interf(sh,MV2,Y,zero,eta41,eta42,spart,sqrts,xi1,xi2,S4Lim%flag)

    if (S4Lim%flag .eqv. .false.) then

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E6 = half*(two*scr(Q,p5)-Qsq)/(scr(p6aux,p5)-scr(p6aux,Q))

       p6 = E6*p6aux

       p7 = Q - p6 - p5
       E7 = p7(1)
       p7aux = p7(1:4)/E7


    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  S4Lim%flag = .true. 

    if (S4Lim%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,p6aux)- scr(p5,p6aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta57 = half*scr(p5aux, p7aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)


       S4Lim%PartFrac = (/xi1,xi2/)


       S4Lim%AmpMom(:,1) = p1
       S4Lim%AmpMom(:,2) = p2
       S4Lim%AmpMom(:,3) = p6
       S4Lim%AmpMom(:,4) = p7
       S4Lim%AmpMom(:,5) = p5 


       S4Lim%LimMom(:,1) = p1
       S4Lim%LimMom(:,2) = p2
       S4Lim%LimMom(:,3) = p4

       S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)&
            /(two*spart) * sector &
            *spart/MV2

       ! ---------  C4 + S4 --------------

       eta42 = zero
       eta41 = one

       s24 = spart*x1*x3
       z4 = x1

       C4S4Lim%Lim_KinInv(1:5) = (/zero,s24,zero,zero,z4/)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)

       C4S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)&
            /(two*spart) * sector &
            *spart/MV2


    endif

    ! --------------- SOFT 5 ----------

    eta42 = x3
    eta41 = one - x3
    eta56 = x4

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,S5Lim%flag)


    if (S5Lim%flag .eqv. .false.) then

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E6 = half*(Qsq)/(scr(p6aux,Q))

       p6 = E6*p6aux

       p7 = Q - p6 
       E7 = p7(1)
       p7aux = p7(1:4)/E7

       if ( p7(1).lt.zero.or.p6(1).lt.zero) S5Lim%flag = .true.

    endif

    if (S5Lim%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,p6aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta57 = half*scr(p5aux, p7aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)


       S5Lim%PartFrac = (/xi1, xi2/)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p6
       S5Lim%AmpMom(:,4) = p7
       S5Lim%AmpMom(:,5) = p4 

       S5Lim%LimMom(:,1:4) = S5Lim%AmpMom(:,1:4)
       S5Lim%LimMom(:,5) = p5

       S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            * spart/MV2


       ! ---------  C5  + S5 --------------

       eta56 = zero

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta57 = half*scr(p5auxC, p7aux)

       s65 = four*E6*E5*x4
       z5 = E5/E6

       C5S5Lim%Lim_KinInv(1:5) = (/zero,s65,zero,zero,z5/)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)

       C5S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            *spart/MV2


    endif

    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3:   C4
    !-------------------------------------------------------------------------------------------------------------

        
    eta42 = zero
    eta41 = one
    eta56 = x4
        
    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,C4Lim%flag)


    if (C4Lim%flag .eqv. .false.) then

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4auxC
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E6 = half*(two*scr(Q,p5)-Qsq)/(scr(p6aux,p5)-scr(p6aux,Q))

       p6 = E6*p6aux

       p7 = Q - p6 - p5
       E7 = p7(1)
       p7aux = p7(1:4)/E7

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C4Lim%flag = .true. 

    if (C4Lim%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,p6aux)- scr(p5,p6aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta57 = half*scr(p5aux, p7aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)

       C4Lim%PartFrac=(/xi1,xi2/)

       C4Lim%AmpMom(:,1) = p1 
       C4Lim%AmpMom(:,2) = p2 - p4 
       C4Lim%AmpMom(:,3) = p6
       C4Lim%AmpMom(:,4) = p7
       C4Lim%AmpMom(:,5) = p5

       s24 = spart*x1*x3 
       z4 = one/(one - x1)

       C4Lim%Lim_KinInv(1:5) = (/s24,zero,zero, z4, zero/)


       C4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2


       
    endif


    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C41  + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    eta42 = zero
    eta41 = one

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,C4S5Lim%flag)


    if (C4S5Lim%flag .eqv. .false.) then 

       p1 = half*sqrts*p1aux 
       p2 = half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4auxC
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E6 = half*(Qsq)/(scr(p6aux,Q))

       p6 = E6*p6aux

       p7 = Q - p6 
       E7 = p7(1)
       p7aux = p7(1:4)/E7

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C4S5Lim%flag = .true. 

    if (C4S5Lim%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,p6aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta57 = half*scr(p5aux, p7aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)

       C4S5Lim%PartFrac=(/xi1,xi2/)

       C4S5Lim%AmpMom(:,1) = p1  
       C4S5Lim%AmpMom(:,2) = p2 - p4
       C4S5Lim%AmpMom(:,3) = p6
       C4S5Lim%AmpMom(:,4) = p7

       C4S5Lim%LimMom(:,1:4) = C4S5Lim%AmpMom(:,1:4)
       C4S5Lim%LimMom(:,5) = p5

       s24 = spart*x1*x3 
       z4 = one/(one - x1)

       C4S5Lim%Lim_KinInv(1:5) = (/zero,zero,zero,s24,z4/)


       C4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

       ! ---------   C4   + C5 + SOFT5 ------------

       eta56 = zero

       z5 =  E5/E6

       s25 = spart *x2*(one-x4) 
       s24 = spart *x1*x3 
       s65 = four*E6*E5*x4

       C4C5S5Lim%Lim_KinInv(1:5) = (/s24,s65,zero,z4,z5/)

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta57 = half*scr(p5auxC, p7aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)

       C4C5S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2


    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: C5  FLM 
    !-------------------------------------------------------------------------------------------------------------

    eta42 = x3
    eta41 = one - x3
    eta56 = zero

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,C5Lim%flag)


    if (C5Lim%flag .eqv. .false.) then 

       p1= half*sqrts*p1aux
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 = E5*p5auxC

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E6 = half*(two*scr(Q,p5)-Qsq)/(scr(p6aux,p5)-scr(p6aux,Q))

       p6 = E6*p6aux

       p7 = Q - p6 - p5
       E7 = p7(1)
       p7aux = p7(1:4)/E7
       

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C5Lim%flag = .true. 

    if (C5Lim%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,p6aux)- scr(p5,p6aux))

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta57 = half*scr(p5auxC, p7aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)

       C5Lim%PartFrac = (/xi1,xi2/)

       C5Lim%AmpMom(:,1) = p1
       C5Lim%AmpMom(:,2) = p2 
       C5Lim%AmpMom(:,3) = p6 + p5
       C5Lim%AmpMom(:,4) = p7 
       C5Lim%AmpMom(:,5) = p4


       s65 = four*E6*E5*x4
       z5 = E6/(E6 + E5) 

       
       C5Lim%Lim_KinInv(1:5) = (/zero,s65,zero,zero,z5/)
       HardProc%PartFrac=(/xi1,xi2/)


       C5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C41  + C52
    !-------------------------------------------------------------------------------------------------------------

    eta42 = zero
    eta41 = one 
    eta56 = zero

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,C4C5Lim%flag)


    if (C4C5Lim%flag .eqv. .false.) then
       
       p1= half*sqrts*p1aux
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4auxC
       E5 = half*sqrts*x2
       p5 =  E5*p5auxC

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E6 = half*(two*scr(Q,p5)-Qsq)/(scr(p6aux,p5)-scr(p6aux,Q))

       p6 = E6*p6aux

       p7 = Q - p6 - p5
       E7 = p7(1)
       p7aux = p7(1:4)/E7

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C4C5Lim%flag = .true. 

    if (C4C5Lim%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,p6aux)- scr(p5,p6aux))

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta57 = half*scr(p5auxC, p7aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)

       C4C5Lim%PartFrac = (/xi1,xi2/)

       C4C5Lim%AmpMom(:,1) = p1 
       C4C5Lim%AmpMom(:,2) = p2 - p4
       C4C5Lim%AmpMom(:,3) = p6 + p5
       C4C5Lim%AmpMom(:,4) = p7 


       s24 = spart *x1*x3 
       s65 = four*E6*E5*x4
       z4 = one/(one - x1)
       z5 = E6/(E6 + E5) 

       C4C5Lim%Lim_KinInv(1:5) = (/s24,s65,zero,z4,z5/)


       C4C5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            * spart/MV2

    
    endif


    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C5S4  
    !-------------------------------------------------------------------------------------------------------------

    eta42 = x3
    eta41 = one - x3
    eta56 = zero

    call get_flag_nnlo_interf(sh,MV2,Y,zero,eta41,eta42,spart,sqrts,xi1,xi2,C5S4Lim%flag)


    if (C5S4Lim%flag .eqv. .false.) then 


       p1= half*sqrts*p1aux
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5auxC

       Q = p1 + p2   !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E6 = half*(two*scr(Q,p5)-Qsq)/(scr(p6aux,p5)-scr(p6aux,Q))

       p6 = E6*p6aux

       p7 = Q - p6 - p5
       E7 = p7(1)
       p7aux = p7(1:4)/E7

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C5S4Lim%flag = .true. 

    if (C5S4Lim%flag .eqv. .false.) then

       kallenF = E6 / (scr(Q,p6aux)- scr(p5,p6aux))

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta57 = half*scr(p5auxC, p7aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)

       C5S4Lim%PartFrac = (/xi1,xi2/)

       C5S4Lim%AmpMom(:,1) = p1
       C5S4Lim%AmpMom(:,2) = p2 
       C5S4Lim%AmpMom(:,3) = p6 + p5
       C5S4Lim%AmpMom(:,4) = p7 


       s24 = spart *x1*x3
       s14 = spart *x1*(one-x3) 
       s65 = four*E6*E5*x4
       s12 = two*scr(p1,p2) 

       z5 = E6/(E6 + E5) 

       C5S4Lim%Lim_KinInv(1:5) = (/s12,s24,s14,s65,z5/)

       C5S4Lim%PartFrac=(/xi1,xi2/)


       C5S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

       ! ---------   C4   + C5 + S4  ------------

       eta42 = zero
       eta41 = one 
       eta56 = zero
       
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,6, sector)

       z4 = x1 
       z5 =  E6/(E6 + E5)

       s24 = spart *x1*x3 
       s65 = four*E5*E6*x4

       C4C5S4Lim%Lim_KinInv(1:5) = (/s24,s65,zero,z4,z5/)

       C4C5S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2


    endif


  end subroutine kinematics_nnlo_5263





   subroutine kinematics_nnlo_5264(MV2,Y,yr,&
       HardProc,S4Lim,S5Lim,S4S5Lim,C4Lim,C4S4Lim,C4S5Lim,C4S4S5Lim,&
       C5Lim,C5S4Lim,C5S5Lim,C5S4S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim)
    real(dp), intent(in) ::   MV2, Y,yr(10)
    type(KinConfig)      :: HardProc,S4Lim,S4S5Lim,S5Lim,C4Lim,C4S4Lim,C4S4S5Lim,C4S5Lim,&
         C5Lim,C5S4Lim,C5S4S5Lim,C5S5Lim,C4C5Lim,C4C5S4Lim,C4C5S5Lim,C4C5S4S5Lim
    real(dp) :: p1(4), p2(4), p3(4), p4(4), p5(4)
    real(dp) :: E6,E7, MV, p67, kallenF, sin5, cos5, pVV(4)
    real(dp) :: x1, x2,x3, x4, cos4, sin4, phi4, phi5, xi1, xi2, phi7, cos7, sin7
    real(dp) :: eta41, eta42, eta51, eta52
    real(dp) :: s14, s25, z4, z5
    real(dp) :: spart,sqrts, s12, s75,s24
    real(dp) :: n6(3), n7(3), n5(3), p6(4),p7(4),Q(4), Qsq, E5, eta56, eta57
    real(dp) :: n4(3), n3(3), n1(3), n2(3)
    real(dp) :: cosphi7, sinphi7, n7t1(3), n7t2(3)
    real(dp) :: n1aux,n2aux, sector, s65
    real(dp) :: p1aux(4), p2aux(4), p4aux(4), p4auxc(4), p4c(4), p5aux(4), p5auxc(4), p6aux(4), p7aux(4)
    integer :: i

    ! initialize all the properties of all kinematic configs
    call initialize_config(HardProc)       !1
    call initialize_config(S4Lim)          !2
    call initialize_config(S5Lim)          !3
    call initialize_config(S4S5Lim)        !4
    call initialize_config(C4Lim)          !5
    call initialize_config(C4S4Lim)        !6
    call initialize_config(C4S5Lim)        !7
    call initialize_config(C4S4S5Lim)      !8
    call initialize_config(C5Lim)          !9
    call initialize_config(C5S4Lim)        !10
    call initialize_config(C5S5Lim)        !11
    call initialize_config(C5S4S5Lim)      !12
    call initialize_config(C4C5Lim)        !13
    call initialize_config(C4C5S4Lim)      !14
    call initialize_config(C4C5S5Lim)      !15
    call initialize_config(C4C5S4S5Lim)    !16

        HardProc%npart = 6    
    S4Lim%npart = 5           
    S4S5Lim%npart = 4         
    S5Lim%npart = 5           
    C4Lim%npart = 5           
    C4S4Lim%npart = 5         
    C4S4S5Lim%npart = 4       
    C4S5Lim%npart =  4        
    C5Lim%npart =  5          
    C5S4Lim%npart = 4         
    C5S4S5Lim%npart = 4       
    C5S5Lim%npart = 5         
    C4C5Lim%npart = 4         
    C4C5S4Lim%npart = 4       
    C4C5S5Lim%npart = 4       
    C4C5S4S5Lim%npart = 4    

    HardProc%ids(1:6) = [0,0,id_el,-id_el,id_g,id_a]     
    S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]            
    S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]            
    C4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]            
    C4S4Lim%ids(1:5) = [0,0,id_el,-id_el,id_a]          
    C4S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4S5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]            
    C5S4Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C5S5Lim%ids(1:5) = [0,0,id_el,-id_el,id_g]          
    C4C5Lim%ids(1:4) = [0,0,id_el,-id_el]          
    C4C5S4Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4C5S5Lim%ids(1:4) = [0,0,id_el,-id_el]        
    C4C5S4S5Lim%ids(1:4) = [0,0,id_el,-id_el]     

    

    !    MV2 and Y are generated
    
    MV =  sqrt(MV2)

    x1 =   yr(1)  !energy of the gluon
    x2 =   yr(2)  !energy of the photon
    x3 =   yr(3)  !angle of the gluon: wrt quark 2
    x4 =   yr(5)  !angle of the photon: wrt lepton 6


  !  print*, (yr(i),i=1,10)


    !x1 = 1d-6
    !x2 = 1d-6
    !x3 = 1d-12
    !x4 = 1d-12


    !gluon
    cos4 = - one + two*x3         
    sin4 = sqrt(one-cos4**2)
    phi4 = twopi*yr(4)           

    eta42 = x3
    eta41 = one - x3

    p4aux = (/ one, sin4*cos(phi4), sin4*sin(phi4), cos4 /)
    p4auxC = (/ one, zero, zero, - one /)
    
    !photon: coordinates
    cos5 = one-two*x4
    sin5 = sqrt(one-cos5**2)
    phi5 = twopi*yr(6)

    !lepton 7: direction
    phi7 = twopi*yr(7)
    cos7 = one -  two*yr(8) 
    sin7 = sqrt(one-cos7**2)
    cosphi7 = cos(phi7)
    sinphi7 = sin(phi7)

    !randomize 1 axis
    n7 = (/sin7*cosphi7, sin7*sinphi7, cos7/)
    n7t1 = (/-cos7*cosphi7,-cos7*sinphi7,sin7/) 
    n7t2 = (/-sinphi7,cosphi7, zero/)

    p7aux = (/one,n7(1),n7(2),n7(3)/)

    !photon: direction
    n5 = cos5*n7+sin5*(cos(phi5)*n7t1+sin(phi5)*n7t2)
    p5aux = (/one,n5(1),n5(2),n5(3)/)

    eta57 = x4

    p5auxC = p7aux

    !quarks
    n1 = (/zero,zero,one/)
    n2 = (/zero,zero,-one/)

    p1aux = (/one,zero,zero,one/)
    p2aux = (/one,zero,zero,-one/)

  

    !  full kinematics


    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,HardProc%flag)

    if (HardProc%flag .eqv. .false.) then 

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(two*scr(Q,p5)-Qsq)/(scr(p7aux,p5)-scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 - p5
       E6 = p6(1)
       p6aux = p6(1:4)/E6


    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  HardProc%flag = .true. 

    if (HardProc%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux)- scr(p5,p7aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta56 = half*scr(p5aux, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7,sector)


       HardProc%AmpMom(:,1) = p1
       HardProc%AmpMom(:,2) = p2
       HardProc%AmpMom(:,3) = p6
       HardProc%AmpMom(:,4) = p7
       HardProc%AmpMom(:,5) = p4
       HardProc%AmpMom(:,6) = p5

       HardProc%PartFrac=(/xi1,xi2/)

       HardProc%wgt = kallenF/two/pi &
                     *(spart/two)*x2 & 
                     * (spart/two)*x1 &  
            * sector  &
            * one/two/spart &
            * spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 2: S4 S5 + collinear all angle identical to the hard piece, do not regenerate 
    !-------------------------------------------------------------------------------------------------------------

    ! -------------- S4 S5

    call get_flag_nnlo_interf(sh,MV2,Y,zero,eta41,eta42,spart,sqrts,xi1,xi2,S4S5Lim%flag)

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2   !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(Qsq)/(scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7
       E6 = p6(1)
       p6aux = p6(1:4)/E6

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  S4S5Lim%flag = .true. 

    if (S4S5Lim%flag .eqv. .false.) then

        kallenF = E7 / (scr(Q,p7aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta56 = half*scr(p5aux, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)

       S4S5Lim%AmpMom(:,1) = p1
       S4S5Lim%AmpMom(:,2) = p2
       S4S5Lim%AmpMom(:,3) = p6
       S4S5Lim%AmpMom(:,4) = p7

       S4S5Lim%LimMom(:,1:4) = S4S5Lim%AmpMom(:,1:4)
       S4S5Lim%LimMom(:,5) = p4                        
       S4S5Lim%LimMom(:,6) = p5

       S4S5Lim%PartFrac=(/xi1,xi2/)

       S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            * spart/MV2

       ! ------------   C4  +  S4 + S5  -----------

       eta42 = zero
       eta41 = one
       eta57 = x4

       p4c = half*sqrts*x1*p4auxC

       s24 = spart*x1*x3 
       z4 = x1 

       C4S4S5Lim%Lim_KinInv(1:5) = (/zero,zero,zero,s24,z4/)

       C4S4S5Lim%LimMom(:,1) = p1 
       C4S4S5Lim%LimMom(:,2) = p2 - p4c
       C4S4S5Lim%LimMom(:,3) = p6
       C4S4S5Lim%LimMom(:,4) = p7                           
       C4S4S5Lim%LimMom(:,5) = p5

      
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)

       C4S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector&
            *spart/MV2

       ! ------------  C5 S4S5 FLM  ------------
       
       eta42 = x3
       eta41 = one - x3
       eta57 = zero

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta56 = half*scr(p5auxC, p6aux)
    

    s12 = two*scr(p1,p2) 
    s24 = spart*x1*x3
    s14 = spart*x1*(one - x3) 

    s75 = four*E7*E5*x4
    z5 = E5/E7

    C5S4S5Lim%Lim_KinInv(1:5) = (/s12,s24,s14,s75,z5/)
    
     call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)

       C5S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector&
            *spart/MV2
       

       ! -----------   C4  + C5 + S4S5  ----------

       eta42 = zero
       eta41 = one
       eta57 = zero

       s24 = spart *x1*x3 
       z4 = x1 
       s75 = four*E7*E5*x4
       z5 = E5/E7

       C4C5S4S5Lim%Lim_KinInv(1:5) = (/s24,s75,zero,z4,z5/)
    
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)

       C4C5S4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            *spart/MV2

    endif

    ! --------------- SOFT 4 ----------
    ! gluon soft

    eta42 = x3
    eta41 = one - x3
    eta57 = x4
   
    call get_flag_nnlo_interf(sh,MV2,Y,zero,eta41,eta42,spart,sqrts,xi1,xi2,S4Lim%flag)

    if (S4Lim%flag .eqv. .false.) then

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(two*scr(Q,p5)-Qsq)/(scr(p7aux,p5)-scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 - p5
       E6 = p6(1)
       p6aux = p6(1:4)/E6


    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  S4Lim%flag = .true. 

    if (S4Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux)- scr(p5,p7aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta56 = half*scr(p5aux, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)


       S4Lim%PartFrac = (/xi1,xi2/)


       S4Lim%AmpMom(:,1) = p1
       S4Lim%AmpMom(:,2) = p2
       S4Lim%AmpMom(:,3) = p6
       S4Lim%AmpMom(:,4) = p7
       S4Lim%AmpMom(:,5) = p5 


       S4Lim%LimMom(:,1) = p1
       S4Lim%LimMom(:,2) = p2
       S4Lim%LimMom(:,3) = p4

       S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)&
            /(two*spart) * sector &
            *spart/MV2

       ! ---------  C4 + S4 --------------

       eta42 = zero
       eta41 = one

       s24 = spart*x1*x3
       z4 = x1

       C4S4Lim%Lim_KinInv(1:5) = (/zero,s24,zero,zero,z4/)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)

       C4S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)&
            /(two*spart) * sector &
            *spart/MV2


    endif

    ! --------------- SOFT 5 ----------

    eta42 = x3
    eta41 = one - x3
    eta57 = x4

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,S5Lim%flag)


    if (S5Lim%flag .eqv. .false.) then

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(Qsq)/(scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 
       E6 = p6(1)
       p6aux = p6(1:4)/E6

       if ( p7(1).lt.zero.or.p6(1).lt.zero) S5Lim%flag = .true.

    endif

    if (S5Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta56 = half*scr(p5aux, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)


       S5Lim%PartFrac = (/xi1, xi2/)

       S5Lim%AmpMom(:,1) = p1
       S5Lim%AmpMom(:,2) = p2
       S5Lim%AmpMom(:,3) = p6
       S5Lim%AmpMom(:,4) = p7
       S5Lim%AmpMom(:,5) = p4 

       S5Lim%LimMom(:,1:4) = S5Lim%AmpMom(:,1:4)
       S5Lim%LimMom(:,5) = p5

       S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            * spart/MV2


       ! ---------  C5  + S5 --------------

       eta57 = zero

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta56 = half*scr(p5auxC, p6aux)

       s75 = four*E7*E5*x4
       z5 = E5/E7

       C5S5Lim%Lim_KinInv(1:5) = (/zero,s75,zero,zero,z5/)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)

       C5S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four) /(two*spart) *sector &
            *spart/MV2


    endif

    !-------------------------------------------------------------------------------------------------------------
    !  STRUCTURE 3:   C4
    !-------------------------------------------------------------------------------------------------------------

        
    eta42 = zero
    eta41 = one
    eta57 = x4
        
    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,C4Lim%flag)


    if (C4Lim%flag .eqv. .false.) then

       p1= half*sqrts*p1aux 
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4auxC
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(two*scr(Q,p5)-Qsq)/(scr(p7aux,p5)-scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 - p5
       E6 = p6(1)
       p6aux = p6(1:4)/E6

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C4Lim%flag = .true. 

    if (C4Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux)- scr(p5,p7aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta56 = half*scr(p5aux, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)

       C4Lim%PartFrac=(/xi1,xi2/)

       C4Lim%AmpMom(:,1) = p1 
       C4Lim%AmpMom(:,2) = p2 - p4 
       C4Lim%AmpMom(:,3) = p6
       C4Lim%AmpMom(:,4) = p7
       C4Lim%AmpMom(:,5) = p5

       s24 = spart*x1*x3 
       z4 = one/(one - x1)

       C4Lim%Lim_KinInv(1:5) = (/s24,zero,zero, z4, zero/)


       C4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2


       
    endif


    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C41  + SOFT 5 
    !-------------------------------------------------------------------------------------------------------------

    eta42 = zero
    eta41 = one

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,C4S5Lim%flag)


    if (C4S5Lim%flag .eqv. .false.) then 

       p1 = half*sqrts*p1aux 
       p2 = half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4auxC
       E5 = half*sqrts*x2
       p5 =  E5*p5aux

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(Qsq)/(scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 
       E6 = p6(1)
       p6aux = p6(1:4)/E6

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C4S5Lim%flag = .true. 

    if (C4S5Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux))

       eta51 = half*scr(p5aux, p1aux)
       eta52 = half*scr(p5aux, p2aux)
       eta56 = half*scr(p5aux, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)

       C4S5Lim%PartFrac=(/xi1,xi2/)

       C4S5Lim%AmpMom(:,1) = p1  
       C4S5Lim%AmpMom(:,2) = p2 - p4
       C4S5Lim%AmpMom(:,3) = p6
       C4S5Lim%AmpMom(:,4) = p7

       C4S5Lim%LimMom(:,1:4) = C4S5Lim%AmpMom(:,1:4)
       C4S5Lim%LimMom(:,5) = p5

       s24 = spart*x1*x3 
       z4 = one/(one - x1)

       C4S5Lim%Lim_KinInv(1:5) = (/zero,zero,zero,s24,z4/)


       C4S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

       ! ---------   C4   + C5 + SOFT5 ------------

       eta57 = zero

       z5 =  E5/E7

       s25 = spart *x2*(one-x4) 
       s24 = spart *x1*x3 
       s75 = four*E7*E5*x4

       C4C5S5Lim%Lim_KinInv(1:5) = (/s24,s75,zero,z4,z5/)

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta56 = half*scr(p5auxC, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)

       C4C5S5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2


    endif

    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 5: C5  FLM 
    !-------------------------------------------------------------------------------------------------------------

    eta42 = x3
    eta41 = one - x3
    eta57 = zero

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,C5Lim%flag)


    if (C5Lim%flag .eqv. .false.) then 

       p1= half*sqrts*p1aux
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 = E5*p5auxC

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(two*scr(Q,p5)-Qsq)/(scr(p7aux,p5)-scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 - p5
       E6 = p6(1)
       p6aux = p6(1:4)/E6
       

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C5Lim%flag = .true. 

    if (C5Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux)- scr(p5,p7aux))

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta56 = half*scr(p5auxC, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)

       C5Lim%PartFrac = (/xi1,xi2/)

       C5Lim%AmpMom(:,1) = p1
       C5Lim%AmpMom(:,2) = p2 
       C5Lim%AmpMom(:,3) = p6 
       C5Lim%AmpMom(:,4) = p7 + p5
       C5Lim%AmpMom(:,5) = p4


       s75 = four*E7*E5*x4
       z5 = E7/(E7 + E5) 

       
       C5Lim%Lim_KinInv(1:5) = (/zero,s75,zero,zero,z5/)
       HardProc%PartFrac=(/xi1,xi2/)


       C5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

    endif


    !-------------------------------------------------------------------------------------------------------------
    !   STRUCTURE 6: C41  + C52
    !-------------------------------------------------------------------------------------------------------------

    eta42 = zero
    eta41 = one 
    eta57 = zero

    call get_flag_nnlo_interf(sh,MV2,Y,x1,eta41,eta42,spart,sqrts,xi1,xi2,C4C5Lim%flag)


    if (C4C5Lim%flag .eqv. .false.) then
       
       p1= half*sqrts*p1aux
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4auxC
       E5 = half*sqrts*x2
       p5 =  E5*p5auxC

       Q = p1 + p2 - p4  !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(two*scr(Q,p5)-Qsq)/(scr(p7aux,p5)-scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 - p5
       E6 = p6(1)
       p6aux = p6(1:4)/E6

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C4C5Lim%flag = .true. 

    if (C4C5Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux)- scr(p5,p7aux))

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta56 = half*scr(p5auxC, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)

       C4C5Lim%PartFrac = (/xi1,xi2/)

       C4C5Lim%AmpMom(:,1) = p1 
       C4C5Lim%AmpMom(:,2) = p2 - p4
       C4C5Lim%AmpMom(:,3) = p6 
       C4C5Lim%AmpMom(:,4) = p7 + p5


       s24 = spart *x1*x3 
       s75 = four*E7*E5*x4
       z4 = one/(one - x1)
       z5 = E7/(E7 + E5) 

       C4C5Lim%Lim_KinInv(1:5) = (/s24,s75,zero,z4,z5/)


       C4C5Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            * spart/MV2

    
    endif


    !-------------------------------------------------------------------------------------------------------------
    ! STRUCTURE 4:  C5S4  
    !-------------------------------------------------------------------------------------------------------------

    eta42 = x3
    eta41 = one - x3
    eta57 = zero

    call get_flag_nnlo_interf(sh,MV2,Y,zero,eta41,eta42,spart,sqrts,xi1,xi2,C5S4Lim%flag)


    if (C5S4Lim%flag .eqv. .false.) then 


       p1= half*sqrts*p1aux
       p2 =half*sqrts*p2aux 

       p4 = half*sqrts*x1*p4aux
       E5 = half*sqrts*x2
       p5 =  E5*p5auxC

       Q = p1 + p2   !momentum of the intermediate boson 

       Qsq = scr(Q,Q)

       E7 = half*(two*scr(Q,p5)-Qsq)/(scr(p7aux,p5)-scr(p7aux,Q))

       p7 = E7*p7aux

       p6 = Q - p7 - p5
       E6 = p6(1)
       p6aux = p6(1:4)/E6

    endif

    if ( p7(1).lt.zero.or.p6(1).lt.zero)  C5S4Lim%flag = .true. 

    if (C5S4Lim%flag .eqv. .false.) then

       kallenF = E7 / (scr(Q,p7aux)- scr(p5,p7aux))

       eta51 = half*scr(p5auxC, p1aux)
       eta52 = half*scr(p5auxC, p2aux)
       eta56 = half*scr(p5auxC, p6aux)

       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)

       C5S4Lim%PartFrac = (/xi1,xi2/)

       C5S4Lim%AmpMom(:,1) = p1
       C5S4Lim%AmpMom(:,2) = p2 
       C5S4Lim%AmpMom(:,3) = p6 
       C5S4Lim%AmpMom(:,4) = p7 + p5


       s24 = spart *x1*x3
       s14 = spart *x1*(one-x3) 
       s75 = four*E7*E5*x4
       s12 = two*scr(p1,p2) 

       z5 = E7/(E7 + E5) 

       C5S4Lim%Lim_KinInv(1:5) = (/s12,s24,s14,s75,z5/)

       C5S4Lim%PartFrac=(/xi1,xi2/)


       C5S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2

       ! ---------   C4   + C5 + S4  ------------

       eta42 = zero
       eta41 = one 
       eta57 = zero
       
       call sector_nnlo_minimal(eta41,eta42,eta51,eta52,eta56,eta57,2,7, sector)

       z4 = x1 
       z5 =  E7/(E7 + E5)

       s24 = spart *x1*x3 
       s75 = four*E5*E7*x4

       C4C5S4Lim%Lim_KinInv(1:5) = (/s24,s75,zero,z4,z5/)

       C4C5S4Lim%wgt = two/pi*kallenF&
            *( x1 ) * ( spart/four) &
            *( x2 ) * ( spart/four)/(two*spart) &
            *sector &
            *spart/MV2


    endif


  end subroutine kinematics_nnlo_5264
  
  
  

  
  function NF1(x3,x4,x5) 
    implicit none  
    real(dp),intent(in) :: x3,x4,x5
    real(dp) :: NF1
    
    if (x4.eq.zero) then 
       NF1 = one 
       return 
    elseif(x4.eq.one) then 
       NF1 = four*x5*(one-x3) 
       return 
    else
       NF1 = 1.0_dp + x4*(1.0_dp-2.0_dp*x3) & 
            - 2.0_dp*(1.0_dp-2.0_dp*x5)*sqrt(abs(x4*(1.0_dp-x3)*(1.0_dp-x3*x4)))
    endif
    
  end function NF1

  function sc3(a,b)
    real(dp) :: sc3,a(3),b(3)

    sc3 = dot_product(a,b)

  end function sc3
  
end module mod_kinematicsNNLO_EunordAord_ga
   
    
!!!!! TO DO:

! 1) check which invariants do we actually need for the kernels -> see appendix B in Kostantinos thesis
! 2) check which partitioning are actually affected by the limits
! 3) check which phase space measures are actually affected by the limits,and why -> see footnote6 in mixed
!    QCD.QED W paper
! 4) check if combinations of limits are more conveniently generated than limits by their own
! 5) find a way to consistently generate final state lepton momenta -> see Kostantin thesis
! 6) check if sector functions are correctly generated by their own
! 7) check where %LimMom are actually used and why -> verify that I'm doing the right thing,consistently in the code
! 8) verify limits on energy fractions z!!!!!

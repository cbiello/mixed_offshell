module mod_splitting
  use mod_types; use common_def; use mod_consts_dp
  use mod_auxfunctions; use mod_procfunctions
  implicit none

  private

  !-- flip -> for non-diagonal channels, change color factor and flip overall sign

  !-- tree-level
  public :: Pgg_spav, Pgq_spav, Pqq, Pqg
  public :: Pgg_spav_i, Pgq_spav_i, Pqq_i, Pqg_i, Pgq_spav_flip, Pqg_flip
  public :: comp_hell_mat_pgg, comp_hell_mat_pgq, comp_hell_mat_pgq_flip
  !--
  !-- triple collinear
  public :: Pggg_spav, Pgqqb_spav, Pgqqb_spav_flip, Pggq, Pggq_flip, PqbqR, Pqbqq_id_symm
  public :: Pggg_ds, Pggq_ds, Pqqbg_ds, PqqbQ_ds
  public :: Pggg_dsso, Pggq_dsso
  public :: Ptcso_spav_gg_g, Ptcso_spav_qg_q, Ptcso_spav_qqb_g, Ptcso_spav_qqb_g_flip
  public :: Ptcso_gg_q, Ptcso_gg_q_flip, Ptcso_qbq_r
  public :: Ptcso_gg_g_ds, Ptcso_qqb_g_ds
  public :: Ptcso_gg_q_ds, Ptcso_qqb_q_ds
  !--
  !-- for RV
  public :: Pgg_new_spav_i, Pqqt_new_i, Pqgt_new_i 

  public :: zipp, zipp2, comp_hell_mat

contains

  !----------------------------------------------------------
  !-- final state LO splitting functions
  !----------------------------------------------------------

  function Pgg_spav(z)
    real(dp) :: z,Pgg_spav
    
    Pgg_spav = two*Ca*((one-z)/z + z/(one-z) + z*(one-z))

  end function Pgg_spav

  function Pgq_spav(z)
    real(dp) :: z,Pgq_spav

    Pgq_spav = tr * (one-two*z+two*z**2)

  end function Pgq_spav

  function Pqq(z)
    real(dp) :: z,Pqq

    Pqq = Cf * ((one+z**2)/(one-z))

  end function Pqq

  function Pqg(z)
    real(dp) :: z,Pqg

    Pqg = Cf * (one+(one-z)**2)/z

  end function Pqg

  !----------------------------------------------------------
  !-- final state off-diagonal LO splitting functions
  !-- to be used with initial state crossing, 
  !-- in unphysical z region
  !-- 
  !-- flip color factor and add extra minus sign
  !-- to compensate for fermion crossing --> 
  !-- use as diagonal ones
  !----------------------------------------------------------

  function Pgq_spav_flip(z)
    real(dp) :: z,Pgq_spav_flip

    Pgq_spav_flip = -Cf * (one-two*z+two*z**2)

  end function Pgq_spav_flip

  function Pqg_flip(z)
    real(dp) :: z,Pqg_flip

    Pqg_flip = -tr * (one+(one-z)**2)/z

  end function Pqg_flip

  !----------------------------------------------------------
  !-- initial state LO splitting functions
  !-- repsect to the notes, here there is an extra 1/z because I'm stupid
  !----------------------------------------------------------

  function Pgg_spav_i(z)
    real(dp) :: z,Pgg_spav_i
    
    Pgg_spav_i = two*Ca*((one-z)/z + z/(one-z) + z*(one-z))
    Pgg_spav_i = Pgg_spav_i/z

  end function Pgg_spav_i

  function Pgq_spav_i(z)
    real(dp) :: z,Pgq_spav_i

    Pgq_spav_i = Cf * (one+(one-z)**2)/z
    Pgq_spav_i = Pgq_spav_i/z

  end function Pgq_spav_i

  function Pqq_i(z)
    real(dp) :: z,Pqq_i

    Pqq_i = Cf * ((one+z**2)/(one-z))
    Pqq_i = Pqq_i/z

  end function Pqq_i

  function Pqg_i(z)
    real(dp) :: z,Pqg_i

    Pqg_i = tr * (one-two*z*(one-z))
    Pqg_i = Pqg_i/z

  end function Pqg_i

  !----------------------------------------------------------
  !-- spin-correlated splitting functions
  !----------------------------------------------------------
  
  function comp_hell_mat_pgg(z,nperp,ndir) ! nperp: correlation vector, ndir: collinear direction
    real(dp) :: z,nperp(3),ndir(3)
    complex(dp) :: comp_hell_mat_pgg(-1:1,-1:1)
    real(dp) :: t1(4),t2(4),pp(4),pn(4)
    complex(dp) :: foo,ept_p1,ept_m1
    real(dp) :: P_dmunu,P_kappa

    comp_hell_mat_pgg = zero

    P_dmunu = two*Ca*(z/(one-z)+(one-z)/z)
    P_kappa = four*Ca*z*(one-z)

    t1 = half*(/one,+nperp(1),+nperp(2),+nperp(3)/)
    t2 = half*(/one,-nperp(1),-nperp(2),-nperp(3)/)
    
    pp(:) = (/one,+ndir(1),+ndir(2),+ndir(3)/)
    pn(:) = (/one,-ndir(1),-ndir(2),-ndir(3)/)
    
    foo = sqrt2*aa(pn,pp)
    
    ept_p1 = (bb(pp,t1)*aa(t1,pn)-bb(pp,t2)*aa(t2,pn))/foo
    
    foo = sqrt2*bb(pn,pp)
    
    ept_m1 = -(aa(pp,t1)*bb(t1,pn)-aa(pp,t2)*bb(t2,pn))/foo

    comp_hell_mat_pgg(+1,+1) = P_dmunu+half*P_kappa
    comp_hell_mat_pgg(-1,-1) = conjg(comp_hell_mat_pgg(1,1))
    comp_hell_mat_pgg(+1,-1) = P_kappa * ept_m1 * conjg(ept_p1)
    comp_hell_mat_pgg(-1,+1) = conjg(comp_hell_mat_pgg(1,-1))
    
    return
    
  end function comp_hell_mat_pgg

  !-- g -> q(z) qb
  function comp_hell_mat_pgq(z,nperp,ndir) ! nperp: correlation vector, ndir: collinear direction
    real(dp) :: z,nperp(3),ndir(3)
    complex(dp) :: comp_hell_mat_pgq(-1:1,-1:1)
    real(dp) :: t1(4),t2(4),pp(4),pn(4)
    complex(dp) :: foo,ept_p1,ept_m1
    real(dp) :: P_dmunu,P_kappa

    comp_hell_mat_pgq = czero

    P_dmunu = tr
    P_kappa = -four*tr*z*(one-z)

    t1 = half*(/one,+nperp(1),+nperp(2),+nperp(3)/)
    t2 = half*(/one,-nperp(1),-nperp(2),-nperp(3)/)
    
    pp(:) = (/one,+ndir(1),+ndir(2),+ndir(3)/)
    pn(:) = (/one,-ndir(1),-ndir(2),-ndir(3)/)
    
    foo = sqrt2*aa(pn,pp)
    
    ept_p1 = (bb(pp,t1)*aa(t1,pn)-bb(pp,t2)*aa(t2,pn))/foo
    
    foo = sqrt2*bb(pn,pp)
    
    ept_m1 = -(aa(pp,t1)*bb(t1,pn)-aa(pp,t2)*bb(t2,pn))/foo

    comp_hell_mat_pgq(+1,+1) = P_dmunu+half*P_kappa
    comp_hell_mat_pgq(-1,-1) = conjg(comp_hell_mat_pgq(1,1))
    comp_hell_mat_pgq(+1,-1) = P_kappa * ept_m1 * conjg(ept_p1)
    comp_hell_mat_pgq(-1,+1) = conjg(comp_hell_mat_pgq(1,-1))

    return
    
  end function comp_hell_mat_pgq

  !-- g -> q(z) qb
  function comp_hell_mat_pgq_flip(z,nperp,ndir) ! nperp: correlation vector, ndir: collinear direction
    real(dp) :: z,nperp(3),ndir(3)
    complex(dp) :: comp_hell_mat_pgq_flip(-1:1,-1:1)
    real(dp) :: t1(4),t2(4),pp(4),pn(4)
    complex(dp) :: foo,ept_p1,ept_m1
    real(dp) :: P_dmunu,P_kappa

    comp_hell_mat_pgq_flip = czero

    !-- flip sign and change color factor
    P_dmunu = -Cf
    P_kappa = +four*Cf*z*(one-z)

    t1 = half*(/one,+nperp(1),+nperp(2),+nperp(3)/)
    t2 = half*(/one,-nperp(1),-nperp(2),-nperp(3)/)
    
    pp(:) = (/one,+ndir(1),+ndir(2),+ndir(3)/)
    pn(:) = (/one,-ndir(1),-ndir(2),-ndir(3)/)
    
    foo = sqrt2*aa(pn,pp)
    
    ept_p1 = (bb(pp,t1)*aa(t1,pn)-bb(pp,t2)*aa(t2,pn))/foo
    
    foo = sqrt2*bb(pn,pp)
    
    ept_m1 = -(aa(pp,t1)*bb(t1,pn)-aa(pp,t2)*bb(t2,pn))/foo

    comp_hell_mat_pgq_flip(+1,+1) = P_dmunu+half*P_kappa
    comp_hell_mat_pgq_flip(-1,-1) = conjg(comp_hell_mat_pgq_flip(1,1))
    comp_hell_mat_pgq_flip(+1,-1) = P_kappa * ept_m1 * conjg(ept_p1)
    comp_hell_mat_pgq_flip(-1,+1) = conjg(comp_hell_mat_pgq_flip(1,-1))

    return
    
  end function comp_hell_mat_pgq_flip

  subroutine zipp(me2spin,hell_matrix,res)
    complex(dp), intent(in) :: me2spin(-1:1,-1:1), hell_matrix(-1:1,-1:1)
    real(dp), intent(out) :: res
    integer :: i, j
    
    res = zero 
    
    do i=-1,1,2
    do j=-1,1,2
          
       res = res + real(me2spin(i,j)*hell_matrix(i,j),kind=dp)
          
    enddo
    enddo
    
    return  
    
  end subroutine zipp

  subroutine zipp2(hell_matrix_a,me2spin,hell_matrix_b,res)
    complex(dp), intent(in) :: hell_matrix_a(-1:1,-1:1),me2spin(-1:1,-1:1,-1:1,-1:1),hell_matrix_b(-1:1,-1:1)
    real(dp), intent(out) :: res
    integer :: ia,ja,ib,jb
    
    res = zero 
    
    do ia=-1,1,2
    do ja=-1,1,2
    do ib=-1,1,2
    do jb=-1,1,2
          
       res = res + real(hell_matrix_a(ia,ja)*me2spin(ia,ja,ib,jb)*hell_matrix_b(ib,jb),kind=dp)
          
    enddo
    enddo
    enddo
    enddo
    
    return  
    
  end subroutine zipp2

  !-- this is for nmu nnu * Mmunu, without splitting function
  function comp_hell_mat(nperp,ndir) ! nperp: correlation vector, ndir: collinear direction
    real(dp) :: nperp(3),ndir(3)
    complex(dp) :: comp_hell_mat(-1:1,-1:1)
    real(dp) :: t1(4), t2(4), pp(4), pn(4)
    complex(dp) :: foo,ept_p1,ept_m1
    
    integer :: iep
    
    comp_hell_mat = zero
    
    t1 = half*(/one,+nperp(1),+nperp(2),+nperp(3)/)
    t2 = half*(/one,-nperp(1),-nperp(2),-nperp(3)/)
    
    pp(:) = (/one,+ndir(1),+ndir(2),+ndir(3)/)
    pn(:) = (/one,-ndir(1),-ndir(2),-ndir(3)/)
    
    foo = sqrt2*aa(pn,pp)
    
    ept_p1 = (bb(pp,t1)*aa(t1,pn)-bb(pp,t2)*aa(t2,pn))/foo
    
    foo = sqrt2*bb(pn,pp)
    
    ept_m1 = -(aa(pp,t1)*bb(t1,pn)-aa(pp,t2)*bb(t2,pn))/foo

    comp_hell_mat(1,1) = ept_p1 * conjg(ept_p1) 
    comp_hell_mat(-1,-1) = conjg(comp_hell_mat(1,1))
    comp_hell_mat(1,-1) = ept_m1 * conjg(ept_p1) 
    comp_hell_mat(-1,1) = conjg(comp_hell_mat(1,-1))
           
    return
    
  end function comp_hell_mat

  !---

  !----------------------------------------------------------
  !-- triple-collinear splitting functions
  !----------------------------------------------------------
  
  function Pggg_spav(s12,s13,s23,z1,z2,z3)
    real(dp) :: Pggg_spav
    real(dp), intent(in) :: s12,s13,s23,z1,z2,z3
    real(dp) :: s123

    s123 = s12+s13+s23

    pggg_spav = Pggg_spav_unsymm(s12, s13, s23, s123, z1, z2, z3) + &
         Pggg_spav_unsymm(s13, s12, s23, s123, z1, z3, z2) + &
         Pggg_spav_unsymm(s12, s23, s13, s123, z2, z1, z3) + &
         Pggg_spav_unsymm(s23, s12, s13, s123, z2, z3, z1) + &
         Pggg_spav_unsymm(s13, s23, s12, s123, z3, z1, z2) + &
         Pggg_spav_unsymm(s23, s13, s12, s123, z3, z2, z1)

    pggg_spav = pggg_spav * four/s123**2

    return

  end function Pggg_spav

  function Pggg_spav_unsymm(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pggg_spav_unsymm
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3
    real(dp) :: t123

    t123 = tijk(s12,s13,s23,z1,z2,z3)

    Pggg_spav_unsymm = t123**2/four/s12**2 + three/four + &
         s123/s12 * (four*(z1*z2-one)/(one-z3) + (z1*z2-two)/z3 + three/two + five/two*z3 &
         +(one-z3*(one-z3))**2/z3/z1/(one-z1) ) &
         +s123**2/s12/s13 * ( z1*z2*(one-z2)*(one-two*z3)/z3/(one-z3) + z2*z3 - two &
         + z1*(one+two*z1)/two + (one+two*z1*(one+z1))/two/(one-z2)/(one-z3) + &
         (one-two*z1*(one-z1))/two/z2/z3 )

    Pggg_spav_unsymm = Pggg_spav_unsymm * Ca**2

    return

  end function Pggg_spav_unsymm

  !-

  function Pgqqb_spav(s12,s13,s23,z1,z2,z3)
    real(dp) :: Pgqqb_spav
    real(dp), intent(in) :: s12,s13,s23,z1,z2,z3
    real(dp) :: s123

    s123 = s12+s13+s23

    Pgqqb_spav = tr * Cf * Pgqqb_spav_ab(s12,s13,s23,s123,z1,z2,z3) + &
         tr * Ca * (Pgqqb_spav_nab_unsymm(s12,s13,s23,s123,z1,z2,z3) + &
         Pgqqb_spav_nab_unsymm(s13,s12,s23,s123,z1,z3,z2))

    Pgqqb_spav = Pgqqb_spav * four/s123**2

  end function Pgqqb_spav

  function Pgqqb_spav_flip(s12,s13,s23,z1,z2,z3)
    real(dp) :: Pgqqb_spav_flip
    real(dp), intent(in) :: s12,s13,s23,z1,z2,z3
    real(dp) :: s123

    s123 = s12+s13+s23

    Pgqqb_spav_flip = Cf * Cf * Pgqqb_spav_ab(s12,s13,s23,s123,z1,z2,z3) + &
         Cf * Ca * (Pgqqb_spav_nab_unsymm(s12,s13,s23,s123,z1,z2,z3) + &
         Pgqqb_spav_nab_unsymm(s13,s12,s23,s123,z1,z3,z2))

    Pgqqb_spav_flip = -Pgqqb_spav_flip * four/s123**2

  end function Pgqqb_spav_flip

  function Pgqqb_spav_ab(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pgqqb_spav_ab
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3
    real(dp) :: overomep

    !-- piece proportional to 1/(1-ep)
    overomep = -two * s123**2/s12/s13 * (z1+two*z2*z3) &
         + s123/s12*two*(z1+z2) + s123/s13*two*(z1+z3)

    Pgqqb_spav_ab = -two - s23*(one/s12+one/s13) + two*s123**2/s12/s13*(one+z1**2) &
         -s123/s12*(one+two*z1) -s123/s13*(one+two*z1) + overomep

    return

  end function Pgqqb_spav_ab

  function Pgqqb_spav_nab_unsymm(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pgqqb_spav_nab_unsymm
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3
    real(dp) :: t231,overomep

    t231 = tijk(s23,s12,s13,z2,z3,z1)

    !-- piece proportional to 1/(1-ep)
    overomep = -s123**2/s13/s23*z3**2*(one-z3-two*z1*z2)/z1/(one-z1) &
         -s123/s13*(one-z2)**2*z2/z1/(one-z1) &
         +s123/two/s23*(z1*(z3-z2)**2-two*z2*z3*(one+z1))/z1/(one-z1) &
         +s123**2/two/s12/s13*(z1+two*z2*z3)

    Pgqqb_spav_nab_unsymm = -t231**2/four/s23**2 + s123**2/two/s13/s23*z3*( &
         ((one-z1)**3-z1**3)/z1/(one-z1) ) &
         +s123/two/s13*(one-z2)*(one+one/z1/(one-z1)) &
         +s123/two/s23*(one+z1**3)/z1/(one-z1) &
          -one/four - s123**2/two/s12/s13*(one+z1**2) + overomep

    return

  end function Pgqqb_spav_nab_unsymm
  
  function Pggq(s12,s13,s23,z1,z2,z3)
    real(dp) :: Pggq
    real(dp), intent(in) :: s12,s13,s23,z1,z2,z3
    real(dp) :: s123

    s123 = s12+s13+s23

    pggq = Cf*Cf*tagab*(Pggq_ab_unsymm(s12,s13,s23,s123,z1,z2,z3) + &
         Pggq_ab_unsymm(s12,s23,s13,s123,z2,z1,z3)) + &
         Cf*Ca*tagnab*(Pggq_nab_unsymm(s12,s13,s23,s123,z1,z2,z3) + &
         Pggq_nab_unsymm(s12,s23,s13,s123,z2,z1,z3))

    pggq = pggq * four/s123**2

    return

  end function Pggq

  function Pggq_flip(s12,s13,s23,z1,z2,z3)
    real(dp) :: Pggq_flip
    real(dp), intent(in) :: s12,s13,s23,z1,z2,z3
    real(dp) :: s123

    s123 = s12+s13+s23

    pggq_flip = tr*Cf*tagab*(Pggq_ab_unsymm(s12,s13,s23,s123,z1,z2,z3) + &
         Pggq_ab_unsymm(s12,s23,s13,s123,z2,z1,z3)) + &
         tr*Ca*tagnab*(Pggq_nab_unsymm(s12,s13,s23,s123,z1,z2,z3) + &
         Pggq_nab_unsymm(s12,s23,s13,s123,z2,z1,z3))

    pggq_flip = -pggq_flip * four/s123**2

    return

  end function Pggq_flip

  function Pggq_ab_unsymm(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pggq_ab_unsymm
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3

    Pggq_ab_unsymm = s123**2/two/s13/s23*z3*(one+z3**2)/z1/z2 &
         +s123/s13* (z3*(one-z1)+(one-z2)**3)/z1/z2 &
         - s23/s13

    return

  end function Pggq_ab_unsymm

  function Pggq_nab_unsymm(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pggq_nab_unsymm
    real(dp) :: s12,s13,s23,s123,z1,z2,z3,t123
    real(dp) :: bit

    t123 = tijk(s12,s13,s23,z1,z2,z3)

    bit = (z1*(two-two*z1+z1**2)-z2*(six-six*z2+z2**2))/z2/(one-z3)

    Pggq_nab_unsymm = t123**2/four/s12**2 + one/four + &
         s123**2/two/s12/s13 * ( (one+z3**2)/z2 &
         +(z2**2+two*(one-z2))/(one-z3) ) &
         -s123**2/four/s13/s23*z3 * ( (one+z3**2)/z1/z2 ) &
         +s123/two/s12 * ( bit ) &
         +s123/two/s13 * ( ((one-z2)**3+z3**2-z2)/z2/(one-z3) &
         -(z3*(one-z1)+(one-z2)**3)/z1/z2 )

    return

  end function Pggq_nab_unsymm

  !-- different quarks
  function PqbqR(s12,s13,s23,z1,z2,z3)
    real(dp) :: PqbqR
    real(dp), intent(in) :: s12,s13,s23,z1,z2,z3
    real(dp) :: s123,t123,bit

    t123 = tijk(s12,s13,s23,z1,z2,z3)

    s123 = s12 + s13 + s23

    bit = z1+z2-s12/s123

    PqbqR = half * Cf * tr * s123/s12 * (&
         -t123**2/s12/s123 + (four*z3+(z1-z2)**2)/(z1+z2) + bit)
    
    PqbqR = PqbqR * four/s123**2

    return

  end function PqbqR

  function Pqbqq_id_unsymm(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pqbqq_id_unsymm
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3

    Pqbqq_id_unsymm = two*s23/s12 &
         +s123/s12*( (one+z1**2)/(one-z2)-two*z2/(one-z3) ) &
         -s123**2/s12/s13*z1/two*( (one+z1**2)/(one-z2)/(one-z3) )

    Pqbqq_id_unsymm = Pqbqq_id_unsymm * Cf * (Cf-half*Ca)

    return

  end function Pqbqq_id_unsymm

  function Pqbqq_id_symm(s12,s13,s23,z1,z2,z3)
    real(dp) :: Pqbqq_id_symm
    real(dp), intent(in) :: s12,s13,s23,z1,z2,z3
    real(dp) :: s123,p1,p2

    s123 = s12 + s13 + s23

    p1 = Pqbqq_id_unsymm(s12,s13,s23,s123,z1,z2,z3)
    p2 = Pqbqq_id_unsymm(s13,s12,s23,s123,z1,z3,z2)

    Pqbqq_id_symm = p1+p2

    Pqbqq_id_symm = Pqbqq_id_symm * four/s123**2
    
    return

  end function Pqbqq_id_symm

  !---

  function tijk(sij,sik,sjk,zi,zj,zk)
    real(dp) :: tijk
    real(dp), intent(in) :: sij,sik,sjk,zi,zj,zk

    tijk = two * (zi*sjk-zj*sik)/(zi+zj) + (zi-zj)/(zi+zj)*sij

    return

  end function tijk

  !--

  !-- includes (2/s123)^2 prefactor
  !-- E1,E2 << E3
  function Pggg_ds(s12,s13,s23,z1,z2)
    real(dp), intent(in) :: s12,s13,s23,z1,z2
    real(dp) :: Pggg_ds
    real(dp) :: z12,s123,t123

    z12 = z1+z2

    s123 = s13+s23
    t123 = two * (z1*s23-z2*s13)/z12

    !-- ``abelian-like'' piece
    Pggg_ds = Ca**2*16.0_dp/s13/s23/z1/z2

    !-- ``non-abelian-like'' piece
    Pggg_ds = Pggg_ds + Ca**2 * four*( t123**2/(two*s12**2*s123**2) + &
         (one/(z2*s13)+one/(z1*s23))/s12 + &
         (one/s13+one/s23)/(s12*z12) - one/(s13*s23*z1*z2) &
         +(z1/z2+z2/z1-6.0_dp)/(s123*s12*z12) &
         +(one/(z2*s13)+one/(z1*s23))/(s123*z12) &
         -(one/s13+one/s23)/(s123*z1*z2) )

    return

  end function Pggg_ds

  function Pggq_ds(s12,s13,s23,z1,z2)
    real(dp), intent(in) :: s12,s13,s23,z1,z2
    real(dp) :: Pggq_ds
    real(dp) :: z12,s123,t123

    z12 = z1+z2

    s123 = s13+s23
    t123 = two * (z1*s23-z2*s13)/z12

    !-- ``abelian-like'' piece
    Pggq_ds = Cf*Cf*16.0_dp/s13/s23/z1/z2

    !-- ``non-abelian-like'' piece
    Pggq_ds = Pggq_ds + Ca*Cf * four*( t123**2/(two*s12**2*s123**2) + &
         (one/(z2*s13)+one/(z1*s23))/s12 + &
         (one/s13+one/s23)/(s12*z12) - one/(s13*s23*z1*z2) &
         +(z1/z2+z2/z1-6.0_dp)/(s123*s12*z12) &
         +(one/(z2*s13)+one/(z1*s23))/(s123*z12) &
         -(one/s13+one/s23)/(s123*z1*z2) )

    return

  end function Pggq_ds
    
  function Pqqbg_ds(s12,s13,s23,z1,z2)
    real(dp), intent(in) :: s12,s13,s23,z1,z2
    real(dp) :: Pqqbg_ds
    real(dp) :: z12,s123,t123

    z12 = z1+z2

    s123 = s13+s23
    t123 = two * (z1*s23-z2*s13)/z12

    Pqqbg_ds = two * Ca * tr / (s123*s12) * &
         (-t123**2/s12/s123 + four/z12)

  end function Pqqbg_ds

  function PqqbQ_ds(s12,s13,s23,z1,z2)
    real(dp), intent(in) :: s12,s13,s23,z1,z2
    real(dp) :: PqqbQ_ds
    real(dp) :: z12,s123,t123

    z12 = z1+z2

    s123 = s13+s23
    t123 = two * (z1*s23-z2*s13)/z12

    PqqbQ_ds = two * Cf * tr / (s123*s12) * &
         (-t123**2/s12/s123 + four/z12)

  end function PqqbQ_ds

  !-- E1 << E2 << E3
  function Pggg_dsso(s12,s13,s23,z1,z2)
    real(dp), intent(in) :: s12,s13,s23,z1,z2
    real(dp) :: Pggg_dsso

    !-- abelian-like piece
    Pggg_dsso = Ca**2*16.0_dp/s23/s13/z2/z1

    !-- non-abelian-like piece
    Pggg_dsso = Pggg_dsso + Ca*Ca*(8.0_dp)*(- one/s23/s13/z2/z1 + one/s13/s12/z2 + one/s23/s12/z1)

    return

  end function Pggg_dsso

  function Pggq_dsso(s12,s13,s23,z1,z2)
    real(dp), intent(in) :: s12,s13,s23,z1,z2
    real(dp) :: Pggq_dsso

    !-- abelian piece
    Pggq_dsso = Cf**2*16.0_dp/s23/s13/z2/z1

    !-- non-abelian piece
    Pggq_dsso = Pggq_dsso + Cf*Ca*(8.0_dp)*(- one/s23/s13/z2/z1 + one/s13/s12/z2 + one/s23/s12/z1)

    return

  end function Pggq_dsso

  !-- strongly collinear ordered
  !-- s12 << s13, s23
  function Ptcso_spav_gg_g(s12,s12_3,zz1,zz2,kt1kt2sq)
    real(dp) :: Ptcso_spav_gg_g
    real(dp), intent(in) :: s12,s12_3,zz1,zz2,kt1kt2sq

    Ptcso_spav_gg_g = pgg_spav(zz2)*pgg_spav(zz1) &
         + 8.0_dp * Ca**2 * zz1*(one-zz1)*(one-zz2)/zz2*(kt1kt2sq - half)

    Ptcso_spav_gg_g = Ptcso_spav_gg_g * (two/s12)*(two/s12_3)

    return

  end function Ptcso_spav_gg_g

  function Ptcso_gg_q(s12,s12_3,zz1,zz2,kt1kt2sq)
    real(dp) :: Ptcso_gg_q
    real(dp), intent(in) :: s12,s12_3,zz1,zz2,kt1kt2sq

    Ptcso_gg_q = pqg(zz2)*pgg_spav(zz1) &
         + 8.0_dp * Ca*Cf * zz1*(one-zz1)*(one-zz2)/zz2*(kt1kt2sq - half)

    Ptcso_gg_q = Ptcso_gg_q * (two/s12)*(two/s12_3)

    return

  end function Ptcso_gg_q

  function Ptcso_gg_q_flip(s12,s12_3,zz1,zz2,kt1kt2sq)
    real(dp) :: Ptcso_gg_q_flip
    real(dp), intent(in) :: s12,s12_3,zz1,zz2,kt1kt2sq

    Ptcso_gg_q_flip = pqg_flip(zz2)*pgg_spav(zz1) &
         - 8.0_dp * Ca*tr * zz1*(one-zz1)*(one-zz2)/zz2*(kt1kt2sq - half)

    Ptcso_gg_q_flip = Ptcso_gg_q_flip * (two/s12)*(two/s12_3)

    return

  end function Ptcso_gg_q_flip

  function Ptcso_spav_qg_q(s12,s12_3,zz1,zz2)
    real(dp) :: Ptcso_spav_qg_q
    real(dp), intent(in) :: s12,s12_3,zz1,zz2

    Ptcso_spav_qg_q = Pgq_spav(zz2)*Pqq(zz1) * (two/s12)*(two/s12_3)

    return

  end function Ptcso_spav_qg_q

  function Ptcso_spav_qqb_g(s12,s12_3,zz1,zz2,kt1kt2sq)
    real(dp) :: Ptcso_spav_qqb_g
    real(dp), intent(in) :: s12,s12_3,zz1,zz2,kt1kt2sq

    Ptcso_spav_qqb_g = pgg_spav(zz2)*pgq_spav(zz1) &
         - 8.0_dp * Ca*tr * zz1*(one-zz1)*(one-zz2)/zz2*(kt1kt2sq - half)

    Ptcso_spav_qqb_g = Ptcso_spav_qqb_g * (two/s12)*(two/s12_3)

    return

  end function Ptcso_spav_qqb_g

  function Ptcso_spav_qqb_g_flip(s12,s12_3,zz1,zz2,kt1kt2sq)
    real(dp) :: Ptcso_spav_qqb_g_flip
    real(dp), intent(in) :: s12,s12_3,zz1,zz2,kt1kt2sq

    Ptcso_spav_qqb_g_flip = pgg_spav(zz2)*pgq_spav_flip(zz1) &
         + 8.0_dp * Ca*Cf * zz1*(one-zz1)*(one-zz2)/zz2*(kt1kt2sq - half)

    Ptcso_spav_qqb_g_flip = Ptcso_spav_qqb_g_flip * (two/s12)*(two/s12_3)

    return

  end function Ptcso_spav_qqb_g_flip

  function Ptcso_qbq_r(s12,s12_3,zz1,zz2,kt1kt2sq)
    real(dp) :: Ptcso_qbq_r
    real(dp), intent(in) :: s12,s12_3,zz1,zz2,kt1kt2sq

    Ptcso_qbq_r = Pqg(zz2)*Pgq_spav(zz1) &
         - 8 * Cf * tr * zz1 * (one-zz1) * (one-zz2)/zz2 * (kt1kt2sq - half)

    Ptcso_qbq_r = Ptcso_qbq_r * (two/s12)*(two/s12_3)

    return

  end function Ptcso_qbq_r

  !-- s12<<s13,s23, e12 << e3
  function Ptcso_gg_g_ds(s12,s12_3,z,z3,kt1kt2sq)
    real(dp), intent(in) :: s12,s12_3,z,z3,kt1kt2sq
    real(dp) :: Ptcso_gg_g_ds

    Ptcso_gg_g_ds = two*Ca*z3 * two/s12 * two/s12_3 * (&
         Pgg_spav(z)+four*Ca*z*(one-z)*(kt1kt2sq-half))

    return

  end function Ptcso_gg_g_ds

  function Ptcso_gg_q_ds(s12,s12_3,z,z3,kt1kt2sq)
    real(dp), intent(in) :: s12,s12_3,z,z3,kt1kt2sq
    real(dp) :: Ptcso_gg_q_ds

    Ptcso_gg_q_ds = two*Cf*z3 * two/s12 * two/s12_3 * (&
         Pgg_spav(z)+four*Ca*z*(one-z)*(kt1kt2sq-half))

    return

  end function Ptcso_gg_q_ds

  !-- s12<<s13,s23, e12 << e3
  function Ptcso_qqb_g_ds(s12,s12_3,z,z3,kt1kt2sq)
    real(dp), intent(in) :: s12,s12_3,z,z3,kt1kt2sq
    real(dp) :: Ptcso_qqb_g_ds

    Ptcso_qqb_g_ds = two*Ca*z3 * two/s12 * two/s12_3 * (&
         Pgq_spav(z)-four*tr*z*(one-z)*(kt1kt2sq-half))

    return

  end function Ptcso_qqb_g_ds

  function Ptcso_qqb_q_ds(s12,s12_3,z,z3,kt1kt2sq)
    real(dp), intent(in) :: s12,s12_3,z,z3,kt1kt2sq
    real(dp) :: Ptcso_qqb_q_ds

    Ptcso_qqb_q_ds = two*Cf*z3 * two/s12 * two/s12_3 * (&
         Pgq_spav(z)-four*tr*z*(one-z)*(kt1kt2sq-half))

    return

  end function Ptcso_qqb_q_ds

  !----------------------------------------------------------
  !-- extra structures for RV, initial state
  !-- in the notes, are defined without the 1/z
  !-- here I add it, see comment above
  !----------------------------------------------------------

  function Pgg_new_spav_i(z)
    real(dp) :: z,Pgg_new_spav_i

    Pgg_new_spav_i = Ca * z/z

  end function Pgg_new_spav_i

  !-- tilde version
  function Pqqt_new_i(z)
    real(dp) :: z,Pqqt_new_i

    Pqqt_new_i = -Cf * z/z
    
  end function Pqqt_new_i

  !-- tilde version
  function Pqgt_new_i(z)
    real(dp) :: z,Pqgt_new_i

    Pqgt_new_i = tr * z/z
    
  end function Pqgt_new_i

end module mod_splitting

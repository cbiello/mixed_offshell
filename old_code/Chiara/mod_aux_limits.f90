module mod_aux_limits
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_auxfunctions
  implicit none
  private

  public :: Pqq,Pqg,Pqg_i
  public :: Pgaq,Pgaq_ds
  public :: comp_hell_mat_pgq,zipp
  
contains

  function Pqq(z,myesq)
    real(dp) :: z,Pqq
    real(dp), intent(in) :: myesq

    Pqq = myesq * ((one+z**2)/(one-z))

  end function Pqq

  function Pqg(z,myesq)
    real(dp) :: z,Pqg
    real(dp), intent(in) :: myesq
    
    Pqg = myesq * (one+(one-z)**2)/z
    
  end function Pqg

  function Pqg_i(z,myesq)
    real(dp) :: z,Pqg_i
    real(dp), intent(in) :: myesq
    
    Pqg_i = myesq * (z**2+(one-z)**2)
    
  end function Pqg_i

  !-- g -> q(z) qb
  function comp_hell_mat_pgq(z,nperp,ndir,myesq) ! nperp: correlation vector, ndir: collinear direction
    real(dp), intent(in) :: z,nperp(3),ndir(3),myesq
    complex(dp) :: comp_hell_mat_pgq(-1:1,-1:1)
    real(dp) :: t1(4),t2(4),pp(4),pn(4)
    complex(dp) :: foo,ept_p1,ept_m1
    real(dp) :: P_dmunu,P_kappa

    comp_hell_mat_pgq = czero

    P_dmunu = myesq
    P_kappa = -four*myesq*z*(one-z)

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

  !--
  
  function Pgaq(s12,s13,s23,z1,z2,z3)
    real(dp) :: Pgaq
    real(dp), intent(in) :: s12,s13,s23,z1,z2,z3
    real(dp) :: s123
    real(dp), parameter :: esq = one

    s123 = s12+s13+s23

    !print *, 'removed because of esq, careful'
    
    pgaq = Cf*esq*(Pggq_ab_unsymm(s12,s13,s23,s123,z1,z2,z3) + &
         Pggq_ab_unsymm(s12,s23,s13,s123,z2,z1,z3)) 
    
    pgaq = pgaq * four/s123**2

    return

  end function Pgaq

  function Pggq_ab_unsymm(s12,s13,s23,s123,z1,z2,z3)
    real(dp) :: Pggq_ab_unsymm
    real(dp), intent(in) :: s12,s13,s23,s123,z1,z2,z3

    Pggq_ab_unsymm = s123**2/two/s13/s23*z3*(one+z3**2)/z1/z2 &
         +s123/s13* (z3*(one-z1)+(one-z2)**3)/z1/z2 &
         - s23/s13

    return

  end function Pggq_ab_unsymm
  
  function Pgaq_ds(s12,s13,s23,z1,z2)
    real(dp), intent(in) :: s12,s13,s23,z1,z2
    real(dp) :: Pgaq_ds
    real(dp) :: z12,s123,t123
    real(dp), parameter :: esq = one
    
    z12 = z1+z2
    
    s123 = s13+s23
    t123 = two * (z1*s23-z2*s13)/z12
    
    !-- ``abelian-like'' piece
    !print *, 'removed because of esq, careful'
    Pgaq_ds = Cf*esq*16.0_dp/s13/s23/z1/z2
    
    return
    
  end function Pgaq_ds
  
end module mod_aux_limits

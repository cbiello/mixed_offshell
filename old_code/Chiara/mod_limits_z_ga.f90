module mod_limits_z_ga
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_auxfunctions
  use mod_aux_limits
  implicit none
  real(dp), parameter :: esq = one
  !private

  public :: myEik_g,get_Eik_ph_z_ga,get_Eik_ds_z_ga
  public :: get_Eik_ph_z_ga_etas
  public :: AP_nlo_ph,AP_nnlo_ga
  public :: get_Pgaq_ds_z_ga,get_Pgaq_C_z_ga,get_Pgaq_S_z_ga, &
       get_Pgaq_SC_z_ga,get_Pggq_SC_z_ga
  public :: product_for_pdf,product_for_pdf_matrix

contains

  function myEik_g(p1,p2,pg,charge)
    real(dp), intent(in) :: p1(4),p2(4),pg(4)
    real(dp), intent(in) :: charge
    real(dp) :: myEik_g
    
    myEik_g = two*charge*scr(p1,p2)/scr(p1,pg)/scr(p2,pg)
    
    return 
  end function myEik_g

  !-- this is a duplicate of the one below, but
  !-- a) it has explicit lepton charge
  !-- b) should not loose precision
  !-- trick: all-outgoing, part +q antipart -q
  !--        symmetric under crossing
  !--
  !-- particle id is always the one of HardProc
  !-- subroutine get_Eik_ph_z_ga_etas(charge,etas,eik,photon_i)
  !--            charge(:,i)*charge(:,j)*etas(i,j)/etas(i,photon_i)/etas(j,photon_i)
  !-- eta(2,1) do i=1,n-1,do j = i+1,n
  subroutine get_Eik_ph_z_ga_etas(etas,eik)
    real(dp), intent(in)  :: etas(5,5)
    real(dp), intent(out) :: eik(4)
    real(dp) :: charge(4,4)
    integer  :: i,j

    charge(:,1) = [-Qdn,-Qup, Qdn, Qup]
    charge(:,2) = [ Qdn, Qup,-Qdn,-Qup]
    charge(:,3) =  Q_lep
    charge(:,4) = -Q_lep

    eik = 0
    do i = 1,4
       do j = i+1,4
          eik = eik - &
               charge(:,i)*charge(:,j)*etas(i,j)/etas(i,5)/etas(j,5)
       enddo
    enddo

  end subroutine get_Eik_ph_z_ga_etas

  !-- same as above, but more general
  !-- pass [q1,q2,q3,q4,always]
  subroutine get_qed_eik(charges,etas,array_pos,a_pos,eik)
    real(dp), intent(in)  :: charges(:,:)
    real(dp), intent(in)  :: etas(:,:)
    real(dp), intent(out) :: eik(4)
    integer, intent(in) :: array_pos(:), a_pos
    integer  :: i, j, m, n
    
    eik = zero
    do i = 1, size(array_pos)
       do j = i+1, size(array_pos)
          m = array_pos(i)
          n = array_pos(j)
          eik = eik - &
               charges(:,i)*charges(:,j)*etas(m,n)/etas(m,a_pos)/etas(n,a_pos)
       enddo
    enddo
    
  end subroutine get_qed_eik

  
  subroutine get_Eik_ph_z_ga(p1,p2,p3,p4,p6,eik)
    real(dp), intent(in)  :: p1(4), p2(4), p3(4), p4(4), p6(4)
    real(dp), intent(out) :: eik(4)
    real(dp)              :: eikg, eikph(4)
    real(dp)              :: eikvec4(4),eikvec5(4)
    real(dp)              :: momw4(4),momw5(4)
    real(dp)              :: charge(4)
    real(dp), parameter :: esq = one

    charge = [Qdn,Qup,-Qdn,-Qup]
    
    !p6 photon

    !print *, 'removed because of esq, careful'
    
    eikph(:) = (charge(:)**2*myEik_g(p1,p2,p6,esq) &
         - charge(:)*myEik_g(p1,p3,p6,esq)&
         + charge(:)*myEik_g(p1,p4,p6,esq)&
         + charge(:)*myEik_g(p2,p3,p6,esq)&
         - charge(:)*myEik_g(p2,p4,p6,esq)&
         + myEik_g(p3,p4,p6,esq))
    
    ! d(1) db(2) 
    eik(1) = eikph(1)
    
    ! u(1) ub(2) 
    eik(2) = eikph(2)
    
    ! db(1) d(2) 
    eik(3) = eikph(3)
    
    ! ub(1) u(2) 
    eik(4) = eikph(4)

    return 

  end subroutine get_Eik_ph_z_ga

  subroutine get_Eik_ds_z_ga(p1,p2,p3,p4,p5,p6,eik)
    real(dp), intent(in)  :: p1(4), p2(4), p3(4), p4(4), p5(4), p6(4)
    real(dp), intent(out) :: eik(4)
    real(dp)              :: eikg, eikph(4)
    real(dp)              :: eikvec4(4),eikvec5(4)
    real(dp)              :: momw4(4),momw5(4)
    real(dp)              :: charge(4)
    real(dp), parameter :: esq = one

    charge = (/Qdn,Qup,-Qdn,-Qup/)
    
    !p5 gluon, p6 photon
    eikg = myEik_g(p1,p2,p5,Cf)

    !print *, 'removed because of esq, careful'
    
    eikph(:) = (-charge(:)**2*myEik_g(p1,p2,p6,esq) &
         + charge(:)*myEik_g(p1,p3,p6,esq)&
         - charge(:)*myEik_g(p1,p4,p6,esq)&
         - charge(:)*myEik_g(p2,p3,p6,esq)&
         + charge(:)*myEik_g(p2,p4,p6,esq)&
         - myEik_g(p3,p4,p6,esq))
    
    ! d(1) db(2) 
    eik(1) = eikg*eikph(1)
    
    ! u(1) ub(2) 
    eik(2) = eikg*eikph(2)
    
    ! db(1) d(2) 
    eik(3) = eikg*eikph(3)
    
    ! ub(1) u(2) 
    eik(4) = eikg*eikph(4)
    
    return 
    
  end subroutine get_Eik_ds_z_ga
  
  !--

  subroutine AP_nlo_ph(z,esq,split)
    real(dp), intent(in) :: z,esq
    real(dp), intent(out) :: split(2,2)
    real(dp) :: charge(2)

    charge = (/Qdn,Qup/)

    split(1,:) = charge(:)**2*Pqq(z,esq)
    split(2,:) = charge(:)**2*Pqq(z,esq)

  end subroutine AP_nlo_ph

  subroutine AP_nnlo_ga(s45,s14,s15,z4,z5,z1,split)

    real(dp), intent(in) :: s45,s14,s15,z4,z5,z1
    real(dp), intent(out) :: split(2,2)
    real(dp) :: charge(2)

    charge = (/Qdn,Qup/)
    
    !d dx
    split(1,1) = charge(1)**2*Pgaq(s45,-s14,-s15,z4,z5,z1)
    !u ux
    split(1,2) = charge(2)**2*Pgaq(s45,-s14,-s15,z4,z5,z1)

    !dx d
    split(2,1) = charge(1)**2*Pgaq(s45,-s14,-s15,z4,z5,z1)
    !ux u
    split(2,2) = charge(2)**2*Pgaq(s45,-s14,-s15,z4,z5,z1) 

  end subroutine AP_nnlo_ga

  subroutine get_Pgaq_ds_z_ga(s45,s14,s15,z4,z5,split)
    real(dp), intent(in) :: s45,s14,s15,z4,z5
    real(dp), intent(out) :: split(2,2)
    real(dp) :: charge(2)

    charge = (/Qdn,Qup/)
    
    !d dx
    split(1,1) = charge(1)**2*Pgaq_ds(s45,-s14,-s15,z4,z5)
    !u ux
    split(1,2) = charge(2)**2*Pgaq_ds(s45,-s14,-s15,z4,z5)
    
    !dx d
    split(2,1) = charge(1)**2*Pgaq_ds(s45,-s14,-s15,z4,z5)
    !ux u
    split(2,2) = charge(2)**2*Pgaq_ds(s45,-s14,-s15,z4,z5) 
    
  end subroutine get_Pgaq_ds_z_ga

  subroutine get_Pgaq_C_z_ga(s15,s4_15,z1,z2,esq,split)
    real(dp), intent(in) :: s15,s4_15,z1,z2,esq
    real(dp), intent(out) :: split(2,2)
    real(dp) :: charge(2)

    charge = (/Qdn,Qup/)
    
    !d dx
    split(1,1) = charge(1)**2*four/s15/s4_15*Pqg(z1,Cf)*Pqq(z2,esq)
    !u ux
    split(1,2) = charge(2)**2*four/s15/s4_15*Pqg(z1,Cf)*Pqq(z2,esq)
    
    !dx d
    split(2,1) = charge(1)**2*four/s15/s4_15*Pqg(z1,Cf)*Pqq(z2,esq)
    !ux u
    split(2,2) = charge(2)**2*four/s15/s4_15*Pqg(z1,Cf)*Pqq(z2,esq)
    
  end subroutine get_Pgaq_C_z_ga
  
  subroutine get_Pgaq_S_z_ga(kinlist,split)
    real(dp), intent(in) :: kinlist(5)
    real(dp), intent(out) :: split(2,2)
    real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, z, Pqqonemz
    real(dp) :: Pqqg_s5,charge(2)
    real(dp), parameter :: esq = one
    
    s14 = kinlist(1) 
    s15 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 
    
    z = one-z4
    
    Pqqonemz = one + z**2

    !print *, 'removed because of esq, careful'
    
    Pqqg_s5 = Cf*four*(two*esq*one/s14/s15/z4/z5)*Pqqonemz

    Pqqg_s5 = Pqqg_s5/z

    charge = (/Qdn,Qup/)
    
    !d dx
    split(1,1) = charge(1)**2*Pqqg_s5
    !u ux
    split(1,2) = charge(2)**2*Pqqg_s5

    !dx d
    split(2,1) = charge(1)**2*Pqqg_s5
    !ux u
    split(2,2) = charge(2)**2*Pqqg_s5

  end subroutine get_Pgaq_S_z_ga

  subroutine get_Pgaq_SC_z_ga(s14,s15,z4,z5,esq,split)
     real(dp), intent(in) :: s14,s15,z4,z5,esq
     real(dp), intent(out) :: split(2,2)
     real(dp) :: charge(2)

    charge = (/Qdn,Qup/)
    
    !d dx
    split(1,1) = charge(1)**2*four*esq/s15/z5*two/s14*Pqq(z4,Cf )
    !u ux
    split(1,2) = charge(2)**2*four*esq/s15/z5*two/s14*Pqq(z4,Cf )
    
    !dx d
    split(2,1) = charge(1)**2*four*esq/s15/z5*two/s14*Pqq(z4,Cf )
    !ux u
    split(2,2) = charge(2)**2*four*esq/s15/z5*two/s14*Pqq(z4,Cf )

  end subroutine get_Pgaq_SC_z_ga
  
  !--
  
  subroutine product_for_pdf(vec1,vec2,vec)
    real(dp), intent(in)  :: vec1(2,2),vec2(4)
    real(dp), intent(out) :: vec(2,2)
    
    vec(1,1) = vec1(1,1)*vec2(1)
    vec(1,2) = vec1(1,2)*vec2(2)
    vec(2,1) = vec1(2,1)*vec2(3)
    vec(2,2) = vec1(2,2)*vec2(4)
    
    return
  end subroutine product_for_pdf

  subroutine product_for_pdf_matrix(vec1,vec2,vec)
    real(dp), intent(in)  :: vec1(2,2),vec2(2,2)
    real(dp), intent(out) :: vec(2,2)
    
    vec(1,1) = vec1(1,1)*vec2(1,1)
    vec(1,2) = vec1(1,2)*vec2(1,2)
    vec(2,1) = vec1(2,1)*vec2(2,1)
    vec(2,2) = vec1(2,2)*vec2(2,2)
    
    return
  end subroutine product_for_pdf_matrix


  


  
!   function Eik_w_nlo(p)
!     real(dp), intent(in) :: p(4,5)
!     real(dp) :: Eik_w_nlo
!     real(dp) :: p1(4),p2(4),p4(4)

!     p1 = p(:,1)
!     p2 = p(:,2)
!     p4 = p(:,3)

!     Eik_w_nlo = two*scr(p1,p2)/(scr(p1,p4) + scr(p2,p4))**2

!     return 

!   end function Eik_w_nlo


!   function myEik_g_li(s12,s1g,s2g,charge)
!     real(dp), intent(in) :: s12,s1g,s2g
!     real(dp), intent(in) :: charge
!     real(dp) :: myEik_g_li

!     myEik_g_li = four*charge*s12/s1g/s2g

!     return 

!   end function myEik_g_li

!   subroutine dotfunc(vec1,vec2,vec)
!     real(dp), intent(in)  :: vec1(4),vec2(4)
!     real(dp), intent(out) :: vec(4)

!     vec(1) = vec1(1)*vec2(1)
!     vec(2) = vec1(2)*vec2(2)
!     vec(3) = vec1(3)*vec2(3)
!     vec(4) = vec1(4)*vec2(4)

!     return

!   end subroutine dotfunc

!   !chaged from here

  
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine get_Eik_ph_z_ga_gluon(p,res_S5,eik)
!     real(dp), intent(in)  :: p(4,6), res_S5(6,2)
!     real(dp), intent(out) :: eik(6,2)
!     real(dp)              :: eikph(3,4)
!     real(dp)              :: charge(4), p1(4), p2(4), p3(4), p4(4), p5(4), p6(4)

!     charge = (/- one/three, two/three, one/three, - two/three /)

!     !p(:,1) = (1,0,0,1)
!     !p(:,2) = (1,0,0,-1)
!     !p(:,3) = lepton
!     !p(:,4) = antilepton
!     !p(:,5) = outgoing fermion
!     !p(:,6) = photon

!      p1 = p(:,1) 
!      p2 = p(:,2)
!      p3 = p(:,3)
!      p4 = p(:,4)
!      p5 = p(:,5) 
!      p6 = p(:,6)

!     eikph(1,:) =  (charge(:)**2*myEik_g(p5,p2,p6,esq) &
!          - charge(:)*myEik_g(p5,p3,p6,esq)&
!          + charge(:)*myEik_g(p5,p4,p6,esq)&
!          + charge(:)*myEik_g(p2,p3,p6,esq)&
!          - charge(:)*myEik_g(p2,p4,p6,esq)&
!          + myEik_g(p3,p4,p6,esq))

!  if (eikph(1,1) .ne. eikph(1,1)) then
!           print*, 'myEik_g(p5,p2,p6,esq)', myEik_g(p5,p2,p6,esq)
!           print*, 'myEik_g(p5,p3,p6,esq)', myEik_g(p5,p3,p6,esq)
!           print*, 'myEik_g(p2,p3,p6,esq)', myEik_g(p2,p3,p6,esq)
!           print*, 'myEik_g(p2,p4,p6,esq)', myEik_g(p2,p4,p6,esq)
!           print*, 'myEik_g(p3,p4,p6,esq)', myEik_g(p3,p4,p6,esq)

!           print *, ''

! print*, 'p5', p5
! print*, 'p6', p6
!           print *, ''
!           print*, 'scr(p5,p6)',scr(p5,p6)

!           print *, ''
!    endif

!    eikph(2,:) =  (charge(:)**2*myEik_g(p5,p1,p6,esq) &
!          - charge(:)*myEik_g(p5,p3,p6,esq)&
!          + charge(:)*myEik_g(p5,p4,p6,esq)&
!          + charge(:)*myEik_g(p1,p3,p6,esq)&
!          - charge(:)*myEik_g(p1,p4,p6,esq)&
!          + myEik_g(p3,p4,p6,esq))


!    eikph(3,:) =  (charge(:)**2*myEik_g(p1,p2,p6,esq) &
!          - charge(:)*myEik_g(p1,p3,p6,esq)&
!          + charge(:)*myEik_g(p1,p4,p6,esq)&
!          + charge(:)*myEik_g(p2,p3,p6,esq)&
!          - charge(:)*myEik_g(p2,p4,p6,esq)&
!          + myEik_g(p3,p4,p6,esq))

! ! print*, 'eikph(3,1)', eikph(1,1) 
! ! print*, 

!     ! d(1) db(2)
!     eik(1,1) = eikph(3,1)*res_S5(1,1)

!     ! u(1) ub(2)
!     eik(1,2) = eikph(3,2)*res_S5(3,1)

!     ! db(1) d(2)
!     eik(2,1) = eikph(3,3)*res_S5(3,1)

!     ! ub(1) u(2)
!     eik(2,2) = eikph(3,4)*res_S5(3,1)



!     ! g(1) db(2)
!     eik(3,1) = eikph(1,1)*res_S5(3,1)

!     ! g(1) ub(2)
!     eik(3,2) = eikph(1,2)*res_S5(3,2)

!     ! g(1) d(2)
!     eik(4,1) = eikph(1,3)*res_S5(4,1)

!     ! g(1) u(2)
!     eik(4,2) = eikph(1,4)*res_S5(4,2)

!     ! d(1) g(2) 
!     eik(5,1) = eikph(2,3)*res_S5(5,1)

!     ! u(1) g(2)
!     eik(5,2) = eikph(2,4)*res_S5(5,2)

!     ! dx(1) g(2)
!     eik(6,1) = eikph(2,1)*res_S5(6,1)

!     ! ux(1) g(2)
!     eik(6,2) = eikph(2,2)*res_S5(6,2)

   
!     return

!   end subroutine get_Eik_ph_z_ga_gluon

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  subroutine get_Eik_ph_z_ga_gq(p,res_S5,eik)
!     real(dp), intent(in)  :: p(4,5), res_S5(2,2)
!     real(dp), intent(out) :: eik(2,2)
!     real(dp)              :: eikph(4)
!     real(dp)              :: charge(4), p1(4), p2(4), p3(4), p4(4), p5(4)

!     charge = (/- one/three, two/three, one/three, - two/three /)

!     !p(:,1) = outgoing fermion 
!     !p(:,2) = ingoing quark
!     !p(:,3) = lepton
!     !p(:,4) = antilepton
!     !p(:,5) = photon


!      p1 = p(:,1) 
!      p2 = p(:,2)
!      p3 = p(:,3)
!      p4 = p(:,4)
!      p5 = p(:,5) 
     
!      eikph(:) =  (charge(:)**2*myEik_g(p1,p2,p5,esq) &
!          - charge(:)*myEik_g(p1,p3,p5,esq)&
!          + charge(:)*myEik_g(p1,p4,p5,esq)&
!          + charge(:)*myEik_g(p2,p3,p5,esq)&
!          - charge(:)*myEik_g(p2,p4,p5,esq)&
!          + myEik_g(p3,p4,p5,esq))


!     ! g(1) db(2)
!     eik(1,1) = eikph(1)*res_S5(1,1)

!     ! g(1) ub(2)
!     eik(1,2) = eikph(2)*res_S5(1,2)

!     ! g(1) d(2)
!     eik(2,1) = eikph(3)*res_S5(2,1)

!     ! g(1) u(2)
!     eik(2,2) = eikph(4)*res_S5(2,2)

!     return

!   end subroutine get_Eik_ph_z_ga_gq
 
  
  




!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   subroutine product_for_pdf_gluon(vec1,vec2,vec)
!     real(dp), intent(in)  :: vec1(6,2),vec2(4)
!     real(dp), intent(out) :: vec(4,2)

!     vec(1,1) = vec1(3,1)*vec2(1)
!     vec(1,2) = vec1(3,2)*vec2(2)
!     vec(2,1) = vec1(4,1)*vec2(3)
!     vec(2,2) = vec1(4,2)*vec2(4)

!     vec(3,1) = vec1(5,1)*vec2(3)
!     vec(3,2) = vec1(5,2)*vec2(4)
!     vec(4,1) = vec1(6,1)*vec2(1)
!     vec(4,2) = vec1(6,2)*vec2(2)

!     return

!   end subroutine product_for_pdf_gluon
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    function Pggq_ds(s12,s13,s23,z1,z2)
    real(dp), intent(in) :: s12,s13,s23,z1,z2
    real(dp) :: Pggq_ds
    real(dp) :: z12,s123,t123

    z12 = z1+z2

    s123 = s13+s23
    t123 = two * (z1*s23-z2*s13)/z12

    !-- ``abelian-like'' piece
    Pggq_ds = Cf*esq*16.0_dp/s13/s23/z1/z2

    !-- ``non-abelian-like'' piece
    !Pggq_ds = Pggq_ds + Ca*Cf * four*( t123**2/(two*s12**2*s123**2) + &
    !     (one/(z2*s13)+one/(z1*s23))/s12 + &
    !     (one/s13+one/s23)/(s12*z12) - one/(s13*s23*z1*z2) &
    !     +(z1/z2+z2/z1-6.0_dp)/(s123*s12*z12) &
    !     +(one/(z2*s13)+one/(z1*s23))/(s123*z12) &
    !     -(one/s13+one/s23)/(s123*z1*z2) )

    return

  end function Pggq_ds



  subroutine get_Pggq_ds_z_ga(s45,s14,s15,z4,z5,split)
    real(dp), intent(in) :: s45,s14,s15,z4,z5
    real(dp), intent(out) :: split(2,2)
    real(dp) :: charge(2)

    charge = (/ -one/three, two/three /)
    
    !d dx
    split(1,1) = charge(1)**2*Pggq_ds(s45,-s14,-s15,z4,z5)
    !u ux
    split(1,2) = charge(2)**2*Pggq_ds(s45,-s14,-s15,z4,z5)

    !dx d
    split(2,1) = charge(1)**2*Pggq_ds(s45,-s14,-s15,z4,z5)
    !ux u
    split(2,2) = charge(2)**2*Pggq_ds(s45,-s14,-s15,z4,z5) 
    

  end subroutine get_Pggq_ds_z_ga


  subroutine get_Pggq_SC_z_ga(s14,s15,z4,z5,esq,split)
     real(dp), intent(in) :: s14,s15,z4,z5,esq
     real(dp), intent(out) :: split(2,2)
     real(dp) :: charge(2)

    charge = (/ -one/three, two/three /)
    
    !d dx
    split(1,1) = charge(1)**2*four*esq/s15/z5*two/s14*Pqq(z4,Cf )
    !u ux
    split(1,2) = charge(2)**2*four*esq/s15/z5*two/s14*Pqq(z4,Cf )

    !dx d
    split(2,1) = charge(1)**2*four*esq/s15/z5*two/s14*Pqq(z4,Cf )
    !ux u
    split(2,2) = charge(2)**2*four*esq/s15/z5*two/s14*Pqq(z4,Cf )
    

  end subroutine get_Pggq_SC_z_ga


  subroutine get_Pggq_S_z_ga(kinlist,split)
    real(dp), intent(in) :: kinlist(5)
    real(dp), intent(out) :: split(2,2)
    real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, z, Pqqonemz
    real(dp) :: Pqqg_s5,charge(2)
   
    s14 = kinlist(1) 
    s15 = kinlist(2) 
    s45 = kinlist(3) 
    z4 = kinlist(4) 
    z5 = kinlist(5) 

    z = one-z4

    Pqqonemz = one + z**2

    Pqqg_s5 = Cf*four*(two*esq*one/s14/s15/z4/z5)*Pqqonemz

    Pqqg_s5 = Pqqg_s5/z

    charge = (/ -one/three, two/three /)
    
    !d dx
    split(1,1) = charge(1)**2*Pqqg_s5
    !u ux
    split(1,2) = charge(2)**2*Pqqg_s5

    !dx d
    split(2,1) = charge(1)**2*Pqqg_s5
    !ux u
    split(2,2) = charge(2)**2*Pqqg_s5
      

  end subroutine get_Pggq_S_z_ga

  

!   subroutine get_Pggq_s5(kinlist,myesq,mytr,splitting)
!     real(dp) :: Pggq_s5
!     real(dp), intent(in) :: kinlist(6)
!     real(dp), intent(in) :: myesq,mytr
!     real(dp), intent(out) :: splitting
!     real(dp) :: scr14,scr15,scr45,z2,z4,z5,z

!     scr14 = kinlist(1)*half
!     scr15 = kinlist(2)*half
!     scr45 = kinlist(3)*half

!     z4 = kinlist(4)
!     z5 = kinlist(5)

!     z2 = one-z4
!     z  = z2

!     Pggq_s5 = Ca/(scr15*scr45) + Ca*(z2/z5/scr15/scr14) + (two*myesq-Ca)*(z4/z5/scr45/scr14)

!     !propto tr
!     splitting = Pggq_s5 * Pqg_flip(z,zero,mytr)

!     end subroutine get_Pggq_s5
  

!   !to here

!   subroutine dotfunc_ns_nlo(vec1,vec2,vec)
!     real(dp), intent(in)  :: vec1(12),vec2(4)
!     real(dp), intent(out) :: vec(4)

!     vec(1) = vec1(1)*vec2(1)
!     vec(2) = vec1(2)*vec2(2)
!     vec(3) = vec1(7)*vec2(3)
!     vec(4) = vec1(8)*vec2(4)

!     return

!   end subroutine dotfunc_ns_nlo

!   subroutine dotfunc_qg_nlo(vec1,vec2,vec)
!     real(dp), intent(in)  :: vec1(12),vec2(4)
!     real(dp), intent(out) :: vec(4)

!     vec(1) = vec1(9 )*vec2(1)
!     vec(2) = vec1(10)*vec2(2)
!     vec(3) = vec1(11)*vec2(3)
!     vec(4) = vec1(12)*vec2(4)

!     return

!   end subroutine dotfunc_qg_nlo

!   subroutine get_Eik_w_ns_li(s12,s1p,s2p,eik)
!     real(dp), intent(in)  :: s12,s1p,s2p
!     real(dp), intent(out) :: eik(4)
!     real(dp) :: term1,term2,term3,term4
!     real(dp) :: p12,p1p,p2p,p1w,p2w,pwp,pww

!     p12 = s12/two
!     p1p = s1p/two
!     p2p = s2p/two
!     p1w = s12/two
!     p2w = s12/two
!     pwp = (s1p + s2p)/two
!     pww = s12

!     term1 = + p12/p1p/p2p*two
!     term2 = + p1w/p1p/pwp*two
!     term3 = + p2w/p2p/pwp*two
!     term4 = + pww/pwp/pwp

!     ! u(1) db(2) -> W+
!     eik(1) = + Qup*Qdn*term1 &
!          + Qup*(Qup - Qdn)*term2 &
!          - Qdn*(Qup - Qdn)*term3 &
!          - (Qup - Qdn)*(Qup - Qdn)*term4

!     ! d(1) ub(2) -> W-
!     eik(2) = + Qup*Qdn*term1 &
!          + Qdn*(- Qup + Qdn)*term2 &
!          - Qup*(- Qup + Qdn)*term3 &
!          - (- Qup + Qdn)*(- Qup + Qdn)*term4

!     ! ub(1) d(2) -> W-
!     eik(3) = + Qup*Qdn*term1 &
!          - Qup*(- Qup + Qdn)*term2 &
!          + Qdn*(- Qup + Qdn)*term3 &
!          - (- Qup + Qdn)*(- Qup + Qdn)*term4

!     ! db(1) u(2) -> W+
!     eik(4) = + Qup*Qdn*term1 &
!          - Qdn*(Qup - Qdn)*term2 &
!          + Qup*(Qup - Qdn)*term3 &
!          - (Qup - Qdn)*(Qup - Qdn)*term4

!     return 

!   end subroutine get_Eik_w_ns_li

!   subroutine get_Eik_w_ns(p1,p2,pw,pp,eik)
!     real(dp), intent(in)  :: p1(4),p2(4),pw(4),pp(4)
!     real(dp), intent(out) :: eik(4)
!     real(dp) :: term1,term2,term3,term4
!     real(dp) :: p12,p1p,p2p,p1w,p2w,pwp,pww

!     p12 = scr(p1,p2)
!     p1p = scr(p1,pp)
!     p2p = scr(p2,pp)
!     p1w = scr(p1,pw)
!     p2w = scr(p2,pw)
!     pwp = scr(pw,pp)
!     pww = scr(pw,pw)

!     term1 = + p12/p1p/p2p*two
!     term2 = + p1w/p1p/pwp*two
!     term3 = + p2w/p2p/pwp*two
!     term4 = + pww/pwp/pwp

!     ! u(1) db(2) -> W+
!     eik(1) = + Qup*Qdn*term1 &
!          + Qup*(Qup - Qdn)*term2 &
!          - Qdn*(Qup - Qdn)*term3 &
!          - (Qup - Qdn)*(Qup - Qdn)*term4

!     ! d(1) ub(2) -> W-
!     eik(2) = + Qup*Qdn*term1 &
!          + Qdn*(- Qup + Qdn)*term2 &
!          - Qup*(- Qup + Qdn)*term3 &
!          - (- Qup + Qdn)*(- Qup + Qdn)*term4

!     ! ub(1) d(2) -> W-
!     eik(3) = + Qup*Qdn*term1 &
!          - Qup*(- Qup + Qdn)*term2 &
!          + Qdn*(- Qup + Qdn)*term3 &
!          - (- Qup + Qdn)*(- Qup + Qdn)*term4

!     ! db(1) u(2) -> W+
!     eik(4) = + Qup*Qdn*term1 &
!          - Qdn*(Qup - Qdn)*term2 &
!          + Qup*(Qup - Qdn)*term3 &
!          - (Qup - Qdn)*(Qup - Qdn)*term4

!     return 

!   end subroutine get_Eik_w_ns

!   subroutine get_Eik_C15S5(sinv,efrac,eik)
!     real(dp), intent(in)  :: sinv,efrac
!     real(dp), intent(out) :: eik(4)
!     real(dp) :: Qub,Qdb

!     Qub = - Qup
!     Qdb = - Qdn

!     ! u(1) db(2) -> W+
!     eik(1) = + four*Qup**2/sinv/efrac

!     ! d(1) ub(2) -> W-
!     eik(2) = + four*Qdn**2/sinv/efrac

!     ! ub(1) d(2) -> W-
!     eik(3) = + four*Qub**2/sinv/efrac

!     ! db(1) u(2) -> W+
!     eik(4) = + four*Qdb**2/sinv/efrac

!     return 

!   end subroutine get_Eik_C15S5

!   subroutine get_Eik_C25S5(sinv,efrac,eik)
!     real(dp), intent(in)  :: sinv,efrac
!     real(dp), intent(out) :: eik(4)
!     real(dp) :: Qub,Qdb

!     Qub = - Qup
!     Qdb = - Qdn

!     ! u(1) db(2) -> W+
!     eik(1) = + four*Qdb**2/sinv/efrac

!     ! d(1) ub(2) -> W-
!     eik(2) = + four*Qub**2/sinv/efrac

!     ! ub(1) d(2) -> W-
!     eik(3) = + four*Qdn**2/sinv/efrac

!     ! db(1) u(2) -> W+
!     eik(4) = + four*Qup**2/sinv/efrac

!     return 

!   end subroutine get_Eik_C25S5

!   subroutine Pqq1_w(z,split)
!     real(dp), intent(in)  :: z
!     real(dp), intent(out) :: split(4)
!     real(dp)              :: pqq
!     real(dp)              :: Qub,Qdb

!     Qub = - Qup
!     Qdb = - Qdn

!     pqq = ((one+z**2)/(one-z))

!     ! u(1) db(2) -> W+
!     split(1) = Qup**2*pqq

!     ! d(1) ub(2) -> W-
!     split(2) = Qdn**2*pqq

!     ! ub(1) d(2) -> W-
!     split(3) = Qub**2*pqq

!     ! db(1) u(2) -> W+
!     split(4) = Qdb**2*pqq    

!   end subroutine Pqq1_w

!   subroutine Pqq2_w(z,split)
!     real(dp), intent(in)  :: z
!     real(dp), intent(out) :: split(4)
!     real(dp)              :: pqq
!     real(dp)              :: Qub,Qdb

!     Qub = - Qup
!     Qdb = - Qdn

!     pqq = ((one+z**2)/(one-z))

!     ! u(1) db(2) -> W+
!     split(1) = Qdb**2*pqq

!     ! d(1) ub(2) -> W-
!     split(2) = Qub**2*pqq

!     ! ub(1) d(2) -> W-
!     split(3) = Qdn**2*pqq

!     ! db(1) u(2) -> W+
!     split(4) = Qup**2*pqq

!   end subroutine Pqq2_w

!   subroutine Pqq1_w_i(z,split)
!     real(dp), intent(in)  :: z
!     real(dp), intent(out) :: split(4)
!     real(dp)              :: pqq
!     real(dp)              :: Qub,Qdb

!     Qub = - Qup
!     Qdb = - Qdn

!     pqq = ((one+z**2)/(one-z))
!     pqq = pqq/z

!     ! u(1) db(2) -> W+
!     split(1) = Qup**2*pqq

!     ! d(1) ub(2) -> W-
!     split(2) = Qdn**2*pqq

!     ! ub(1) d(2) -> W-
!     split(3) = Qub**2*pqq

!     ! db(1) u(2) -> W+
!     split(4) = Qdb**2*pqq

!   end subroutine Pqq1_w_i

!   subroutine Pqq2_w_i(z,split)
!     real(dp), intent(in)  :: z
!     real(dp), intent(out) :: split(4)
!     real(dp)              :: pqq
!     real(dp)              :: Qub,Qdb

!     Qub = - Qup
!     Qdb = - Qdn

!     pqq = ((one+z**2)/(one-z))
!     pqq = pqq/z

!     ! u(1) db(2) -> W+
!     split(1) = Qdb**2*pqq

!     ! d(1) ub(2) -> W-
!     split(2) = Qub**2*pqq

!     ! ub(1) d(2) -> W-
!     split(3) = Qdn**2*pqq

!     ! db(1) u(2) -> W+
!     split(4) = Qup**2*pqq

!   end subroutine Pqq2_w_i


!   subroutine get_Eik_ds_w(p1,p2,p4,p5,eik)
!     real(dp), intent(in)  :: p1(4), p2(4), p4(4), p5(4)
!     real(dp), intent(out) :: eik(4)
!     real(dp)              :: eikg4, eikg5
!     real(dp)              :: eikvec4(4),eikvec5(4)
!     real(dp)              :: momw4(4),momw5(4)

!     !p4 gluon, p5 photon
!     eikg4 = myEik_g(p1,p2,p4,Cf)
!     momw5 = p1 + p2
!     call get_Eik_w_ns(p1,p2,momw5,p5,eikvec5)

!     !p4 gluon, p3 photon
!     eikg5 = myEik_g(p1,p2,p5,Cf)
!     momw4 = p1 + p2
!     call get_Eik_w_ns(p1,p2,momw4,p4,eikvec4)

!     ! u(1) db(2) -> W+
!     eik(1) = eikg4*eikvec5(1) + eikg5*eikvec4(1)

!     ! d(1) ub(2) -> W-
!     eik(2) = eikg4*eikvec5(2) + eikg5*eikvec4(2)

!     ! ub(1) d(2) -> W-
!     eik(3) = eikg4*eikvec5(3) + eikg5*eikvec4(3)

!     ! db(1) u(2) -> W+
!     eik(4) = eikg4*eikvec5(4) + eikg5*eikvec4(4)

!     return 

!   end subroutine get_Eik_ds_w

!   subroutine get_Eik_dss5_w(p1,p2,p4,p5,eik)
!     real(dp), intent(in)  :: p1(4), p2(4), p4(4), p5(4)
!     real(dp), intent(out) :: eik(4)
!     real(dp)              :: thiseik(4)

!     call get_Eik_ds_w(p1,p2,p4,p5,thiseik)
!     eik = thiseik

!     return 

!   end subroutine get_Eik_dss5_w

!   subroutine get_Eik_C15ds_w(sinv,efrac,p1,p2,p4,eik)
!     real(dp), intent(in)  :: p1(4), p2(4), p4(4)
!     real(dp), intent(in)  :: sinv,efrac
!     real(dp), intent(out) :: eik(4)
!     real(dp)              :: eikg4, eikg5
!     real(dp)              :: eikvec4(4),eikvec5(4)
!     real(dp)              :: momw4(4),momw5(4)

!     !p4 gluon, p5 photon
!     eikg4 = myEik_g(p1,p2,p4,Cf)
!     call get_Eik_C15S5(sinv,efrac,eikvec5)

!     !p5 gluon, p4 photon
!     eikg5 = four*Cf/sinv/efrac
!     momw4 = p1 + p2
!     call get_Eik_w_ns(p1,p2,momw4,p4,eikvec4)

!     ! u(1) db(2) -> W+
!     eik(1) = eikg4*eikvec5(1) + eikg5*eikvec4(1)

!     ! d(1) ub(2) -> W-
!     eik(2) = eikg4*eikvec5(2) + eikg5*eikvec4(2)

!     ! ub(1) d(2) -> W-
!     eik(3) = eikg4*eikvec5(3) + eikg5*eikvec4(3)

!     ! db(1) u(2) -> W+
!     eik(4) = eikg4*eikvec5(4) + eikg5*eikvec4(4)

!     return 

!   end subroutine get_Eik_C15ds_w

!   subroutine get_Eik_C15dss5_w(sinv,efrac,p1,p2,p4,eik)
!     real(dp), intent(in)  :: p1(4), p2(4), p4(4)
!     real(dp), intent(in)  :: sinv,efrac
!     real(dp), intent(out) :: eik(4)
!     real(dp)              :: thiseik(4)

!     call get_Eik_C15ds_w(sinv,efrac,p1,p2,p4,thiseik)
!     eik = thiseik

!     return 

!   end subroutine get_Eik_C15dss5_w

!   subroutine get_Eik_TC1ds_w(sinv1,sinv2,efrac1,efrac2,eik)
!     real(dp), intent(in)  :: sinv1,sinv2,efrac1,efrac2
!     real(dp), intent(out) :: eik(4)
!     real(dp)              :: eikg4, eikg5
!     real(dp)              :: eikvec4(4),eikvec5(4)

!     !p4 gluon, p5 photon
!     eikg4 = four*Cf/sinv1/efrac1
!     call get_Eik_C15S5(sinv2,efrac2,eikvec5)

!     !p5 gluon, p4 photon
!     eikg5 = four*Cf/sinv2/efrac2
!     call get_Eik_C15S5(sinv1,efrac1,eikvec4)

!     ! u(1) db(2) -> W+
!     eik(1) = eikg4*eikvec5(1) + eikg5*eikvec4(1)

!     ! d(1) ub(2) -> W-
!     eik(2) = eikg4*eikvec5(2) + eikg5*eikvec4(2)

!     ! ub(1) d(2) -> W-
!     eik(3) = eikg4*eikvec5(3) + eikg5*eikvec4(3)

!     ! db(1) u(2) -> W+
!     eik(4) = eikg4*eikvec5(4) + eikg5*eikvec4(4)

!     return 

!   end subroutine get_Eik_TC1ds_w

!   subroutine get_Eik_TC1dss5_w(sinv1,sinv2,efrac1,efrac2,eik)
!     real(dp), intent(in)  :: sinv1,sinv2,efrac1,efrac2
!     real(dp), intent(out) :: eik(4)
!     real(dp)              :: thiseik(4)

!     call get_Eik_TC1ds_w(sinv1,sinv2,efrac1,efrac2,thiseik)
!     eik = thiseik

!     return 

!   end subroutine get_Eik_TC1dss5_w

!   subroutine get_Eik_TC1C5ds_w(sinv1,sinv2,efrac1,efrac2,eik)
!     real(dp), intent(in)  :: sinv1,sinv2,efrac1,efrac2
!     real(dp), intent(out) :: eik(4)
!     real(dp)              :: thiseik(4)

!     call get_Eik_TC1ds_w(sinv1,sinv2,efrac1,efrac2,thiseik)
!     eik = thiseik

!     return 

!   end subroutine get_Eik_TC1C5ds_w

!   subroutine get_Eik_TC1C5dss5_w(sinv1,sinv2,efrac1,efrac2,eik)
!     real(dp), intent(in)  :: sinv1,sinv2,efrac1,efrac2
!     real(dp), intent(out) :: eik(4)
!     real(dp)              :: thiseik(4)

!     call get_Eik_TC1ds_w(sinv1,sinv2,efrac1,efrac2,thiseik)
!     eik = thiseik

!     return 

!   end subroutine get_Eik_TC1C5dss5_w

!   subroutine get_Eik_TC1s5_w(sinv1,sinv2,efrac1,efrac2,eik)
!     real(dp), intent(in)  :: sinv1,sinv2,efrac1,efrac2
!     real(dp), intent(out) :: eik(4)
!     real(dp)              :: splitg4, eikg5
!     real(dp)              :: splitvec4(4),eikvec5(4)
!     real(dp)              :: myz

!     myz = one-efrac1

!     !p4 gluon, p5 photon
!     splitg4 = Pqq_i(myz,Cf)/sinv1
!     call get_Eik_C15S5(sinv2,efrac2,eikvec5)

!     !p5 gluon, p4 photon
!     call Pqq1_w_i(myz,splitvec4)
!     splitvec4 = splitvec4/sinv1
!     eikg5 = four*Cf/sinv2/efrac2*(-one)

!     ! u(1) db(2) -> W+
!     eik(1) = splitg4*eikvec5(1) + eikg5*splitvec4(1)

!     ! d(1) ub(2) -> W-
!     eik(2) = splitg4*eikvec5(2) + eikg5*splitvec4(2)

!     ! ub(1) d(2) -> W-
!     eik(3) = splitg4*eikvec5(3) + eikg5*splitvec4(3)

!     ! db(1) u(2) -> W+
!     eik(4) = splitg4*eikvec5(4) + eikg5*splitvec4(4)

!     return 

!   end subroutine get_Eik_TC1s5_w

!   subroutine get_Eik_TC1C5s5_w(sinv1,sinv2,efrac1,efrac2,eik)
!     real(dp), intent(in)  :: sinv1,sinv2,efrac1,efrac2
!     real(dp), intent(out) :: eik(4)
!     real(dp)              :: thiseik(4)

!     call get_Eik_TC1s5_w(sinv1,sinv2,efrac1,efrac2,thiseik)
!     eik = thiseik

!     return 

!   end subroutine get_Eik_TC1C5s5_w

!   !   !-- soft-collinear limits

!   function Pgg_c54_ds(listkin)
!     real(dp), intent(in) :: listkin(6)
!     real(dp) :: Pgg_c54_ds
!     real(dp) :: s12,s14,s24,s45,z5,kt1kt2sq
!     real(dp) :: e4_over_e45

!     s12 = listkin(1)
!     s14 = listkin(2)
!     s24 = listkin(3)

!     s45 = listkin(4)
!     z5  = listkin(5)

!     kt1kt2sq = listkin(6)

!     e4_over_e45 = one-z5

!     !propto Ca
!     Pgg_c54_ds = two/s45 * ( pgg_spav(z5,Cf) + (kt1kt2sq-half)*(four*Ca*z5*(one-z5)) ) &
!          * four * Cf * e4_over_e45**2 * s12/s14/s24

!   end function Pgg_c54_ds

!   function Pgq_c54_ds(listkin)
!     real(dp), intent(in) :: listkin(6)
!     real(dp) :: Pgq_c54_ds
!     real(dp) :: s12,s14,s24,s45,z5,kt1kt2sq
!     real(dp) :: e4_over_e45

!     s12 = listkin(1)
!     s14 = listkin(2)
!     s24 = listkin(3)

!     s45 = listkin(4)
!     z5  = listkin(5)

!     kt1kt2sq = listkin(6)

!     e4_over_e45 = one-z5

!     !propto tr
!     Pgq_c54_ds = two/s45 * ( pgq_spav(z5,Cf) + (kt1kt2sq-half)*(-four*tr*z5*(one-z5)) ) &
!          * four * Cf * e4_over_e45**2 * s12/s14/s24

!   end function Pgq_c54_ds

!   function Pqqg_s5(kinlist) 
!     real(dp) :: Pqqg_s5 
!     real(dp), intent(in) :: kinlist(5)
!     real(dp) :: s14,s15,s45,z4,z5, z45, s13, s23, s12, s123, t123, z1, z2, z3, z, Pqqonemz

!     s14 = kinlist(1) 
!     s15 = kinlist(2) 
!     s45 = kinlist(3) 
!     z4 = kinlist(4) 
!     z5 = kinlist(5) 

!     z = one-z4

!     Pqqonemz = one + z**2

!     Pqqg_s5 = Cf*four*((two*esq-Ca)*one/s14/s15/z4/z5 + Ca*(one/s15/s45/z4 + one/s14/s45/z5))*Pqqonemz

!     Pqqg_s5 = Pqqg_s5/z

!   end function Pqqg_s5

    subroutine get_Pggq_C_z_ga(s15,s4_15,z1,z2,esq,split)
    real(dp), intent(in) :: s15,s4_15,z1,z2,esq
    real(dp), intent(out) :: split(2,2)
    real(dp) :: charge(2)

    charge = (/ -one/three, two/three /)
    
    !d dx
    split(1,1) = charge(1)**2*four/s15/s4_15*Pqg(z1,Cf)*Pqq(z2,esq)
    !u ux
    split(1,2) = charge(2)**2*four/s15/s4_15*Pqg(z1,Cf)*Pqq(z2,esq)

    !dx d
    split(2,1) = charge(1)**2*four/s15/s4_15*Pqg(z1,Cf)*Pqq(z2,esq)
    !ux u
    split(2,2) = charge(2)**2*four/s15/s4_15*Pqg(z1,Cf)*Pqq(z2,esq)

  end subroutine get_Pggq_C_z_ga

end module mod_limits_z_ga


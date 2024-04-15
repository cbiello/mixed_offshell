module mod_subtrfn_nlo_z
! integrated subtraction terms for NLO QCD and NLO EW
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_auxfunctions
  use mod_splittings_bare

  implicit none
  integer, public, parameter :: nintsub= 1 !-- number of scales used in integrated subtraction
                                           !-- todo-raoul: is this the best way to do this?
  integer, parameter        :: iQqsq=1
  integer, parameter        :: iQeQq=2
  integer, parameter        :: iQesq=3
  private


  public :: sub_g_qqb_qqb_z, sub_g_gq_qqb_z
  public :: sub_a_qqb_qqb_z, sub_a_aq_qqb_z, sub_a_aq_aa_z

  ! -- output indices: -1=delta, 0=reg+pls, 1=pls    
  ! -- naming convention: [sub]_[emitted-parton]_[initial-state-unreduced]_[initial-state-reduced]_[boosts]
  ! -- where [boosts] indicates if we have a unboosted, boosted, or doubly boosted initial state
  ! -- e.g. if we have qg->q+V and the reduced ME is qqb, we call this g_qg_qqb_z

contains


  !! functions for NLO QCD
  function sub_g_qqb_qqb_z(q2,mu2,z)

    real(dp), intent(in)  :: q2,mu2(:),z
    real(dp)              :: sub_g_qqb_qqb_z(-1:1,nintsub) 
    real(dp)              :: Lmu(size(mu2))
    integer               :: i

    if (size(mu2) .ne. nintsub) then
       print *, "Error: requested size of integrated subtraction array not the same as scale array"
       print *, size(mu2), nintsub
       stop
    endif

    Lmu = log(mu2/q2)
    
    sub_g_qqb_qqb_z = zero

    do i =1,size(mu2)
       sub_g_qqb_qqb_z(-1:1,i) = -Cf*calPp_qq(z) - Cf*PqqAP(z) * Lmu(i)    !todo-raoul: easier way to do this?
    enddo

    sub_g_qqb_qqb_z(-1,:) = sub_g_qqb_qqb_z(-1,:) + Cf/three*pisq   !-- remainder from virtual definition. 
                                                                    !-- Half w.r.t. notes because multiplied by #legs

  end function sub_g_qqb_qqb_z

  function sub_g_gq_qqb_z(q2,mu2,z)
    real(dp), intent(in)  :: q2,mu2(:),z
    real(dp)              :: sub_g_gq_qqb_z(-1:1,nintsub)
    real(dp)              :: Lmu(nintsub)
    integer               :: i


    if (size(mu2) .ne. nintsub) then
       print *, "Error: requested size of integrated subtraction array not the same as scale array"
       print *, size(mu2), nintsub
       stop
    endif

    Lmu = log(mu2/q2)

    sub_g_gq_qqb_z = zero

    do i=1,nintsub
       sub_g_gq_qqb_z(-1:1,i) = -tr*calPp_qg(z) - tr*PqgAP(z) * Lmu(i)
    enddo

  end  function sub_g_gq_qqb_z



  !! functions for NLO EW
  function sub_a_qqb_qqb_z(q2,mu2,z,eta31,eta32,eta41,eta42,Qq,Qe)
    real(dp), intent(in)  :: q2,mu2(:),z,eta31,eta32,eta41,eta42,Qq,Qe
    real(dp)              :: sub_a_qqb_qqb_z(-1:1,nintsub) 
    real(dp)              :: Lmu(nintsub), logeta31,logeta32,logeta41,logeta42
    integer               :: i


    if (size(mu2) .ne. nintsub) then
       print *, "Error: requested size of integrated subtraction array not the same as scale array"
       print *, size(mu2), nintsub
       stop
    endif

    Lmu = log(mu2/q2)

!    LsoEmsq = log(q2/Emax**2)   ! Log[s/Emax**2]
    logeta31 = log(eta31)
    logeta32 = log(eta32)
    logeta41 = log(eta41)
    logeta42 = log(eta42)
    
    sub_a_qqb_qqb_z = zero

    do i =1,nintsub
       
       sub_a_qqb_qqb_z(-1:1,i) = (-Qq**2*calPp_qq(z) - Qq**2*PqqAP(z) * Lmu(i))       ! plus+reg and plus


       ! delta, factor 1/2 because multiplied by number of legs in sector routine
       sub_a_qqb_qqb_z(-1,i) = sub_a_qqb_qqb_z(-1,i) + Qq**2*two/three*pisq/two
       sub_a_qqb_qqb_z(-1,i) = sub_a_qqb_qqb_z(-1,i) + Qe**2*(13.0_dp - two/three*pisq)/two
       sub_a_qqb_qqb_z(-1,i) = sub_a_qqb_qqb_z(-1,i) + Qe*Qq* &
            ( three*logeta31 - three*logeta41 - three*logeta32 + three*logeta42 +&
            two*real(dilog2(1 - eta31),kind=dp) - two*real(dilog2(1 - eta41),kind=dp) - &
            two*real(dilog2(1 - eta32),kind=dp) + two*real(dilog2(1 - eta42),kind=dp))/two

       
    enddo

  end function sub_a_qqb_qqb_z

  function sub_a_aq_qqb_z(q2,mu2,z)
    real(dp), intent(in)  :: q2,mu2(:),z
    real(dp)              :: sub_a_aq_qqb_z(-1:1,nintsub)
    real(dp)              :: Lmu(nintsub),omz
    integer               :: i

    if (size(mu2) .ne. nintsub) then
       print *, "Error: requested size of integrated subtraction array not the same as scale array"
       print *, size(mu2), nintsub
       stop
    endif


    Lmu = log(mu2/q2)
    omz = one-z

    sub_a_aq_qqb_z = zero
    
    do i=1,nintsub
       sub_a_aq_qqb_z(-1:1,i) = xn*(two*z*omz + (log(omz**2)-log(z)-Lmu(i))*PqgAP(z))
    enddo


  end  function sub_a_aq_qqb_z


  function sub_a_aq_aa_z(q2,mu2,z)
    real(dp), intent(in)  :: q2,mu2(:),z
    real(dp)              :: sub_a_aq_aa_z(-1:1,nintsub)
    real(dp)              :: Lmu(nintsub),omz
    integer               :: i

    if (size(mu2) .ne. nintsub) then
       print *, "Error: requested size of integrated subtraction array not the same as scale array"
       print *, size(mu2), nintsub
       stop
    endif


    Lmu = log(mu2/q2)
    omz = one-z
    
    do i=1,nintsub
       sub_a_aq_aa_z(-1:1,i) = z + (log(omz**2)-log(z)-Lmu(i))*PgqAP(z)
    enddo


  end  function sub_a_aq_aa_z
  
end module mod_subtrfn_nlo_z

    


    
    

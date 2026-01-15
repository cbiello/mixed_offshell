module mod_check_lim
  use mod_types
  use mod_consts_dp
!  use mod_wbfj_N2LOxLO_r_q_raoul
  use mod_limvals
  use mod_proc_parms
  use mod_kinematics_nlo
  use mod_xsects_nloqcd_r
  use mod_xsects_nloewk_r
  use mod_kinematics_nnlo_tc_is_eu_ac
  use mod_xsects_nnlo_rr_ns
  
  private

  integer, parameter :: limdepth_min = -3
  integer, parameter :: limdepth_max = -15
  integer, parameter :: number_limits_nlo = 3
  integer, parameter :: number_terms_nlo = 4
  integer, parameter :: number_limits_nlo_is = 5
  integer, parameter :: number_terms_nlo_is = 6
  integer, parameter :: number_limits_nnlo = 15
  integer, parameter :: number_terms_nnlo = 16
  logical, parameter :: verbose = .true.

  !  integer, public    :: cancelling_terms_nlo(4,3) = reshape([2,1,4,3, 3,4,1,2, 4,3,2,1],[4,3])

  

  public :: check_lim_nlo, check_lim_nlo_is, check_lim_nnlo


contains

  integer function check_lim_nlo(yRnd,ff,vegasweight) result(r)
    ! this checks NLO general, i.e. (I-S_m)(I-C_im), so 1 hard term + 3 subtr terms
!    integer, intent(in)    :: ndim
    real(dp15)    :: yRnd(30), vegasweight, ff(1)
!    interface
!       function f_to_be_checked(xx,ff,weight)
!         use mod_types
!         real(dp15)       :: xx(30), ff(1), weight
!       end  function f_to_be_checked
!    end interface
    real(dp15)             :: x(30),xlim(30)
    integer                :: rlim(number_limits_nlo)
    
    real(dp)               :: limmat(limdepth_min-limdepth_max+1,number_terms_nlo,number_limits_nlo),ssign
    real(dp)               :: limcancellation(limdepth_min-limdepth_max+1,number_limits_nlo,number_terms_nlo/2),loglimcancellation(limdepth_min-limdepth_max+1,number_limits_nlo,number_terms_nlo/2)
    real(dp)               :: target_scaling_nlo(number_limits_nlo,number_terms_nlo/2),scaling(limdepth_min-limdepth_max+1,number_limits_nlo,number_terms_nlo/2)
    integer                :: limdepth, i, ilim,iterm,i1,i2,ifile,icoll,jother
    integer                :: cancelling_terms_nlo(number_limits_nlo,number_terms_nlo/2,2)
    integer, save          :: num_checks = 1
    logical                :: sing_terms_nlo(number_limits_nlo,number_terms_nlo/2)

    call set_sing_terms_nlo(sing_terms_nlo)
    
    cancelling_terms_nlo(1,1,1:2)=[2,1]       ! first  pair for soft limit
    cancelling_terms_nlo(1,2,1:2)=[4,3]       ! second pair for soft limit

    cancelling_terms_nlo(2,1,1:2)=[3,1]       ! first  pair for coll limit
    cancelling_terms_nlo(2,2,1:2)=[4,2]       ! second pair for coll limit

    cancelling_terms_nlo(3,1,1:2)=[4,1]       ! first  pair for softcoll limit
    cancelling_terms_nlo(3,2,1:2)=[3,2]       ! second pair for softcoll limit

    target_scaling_nlo(1,:) = -one            ! soft limit
    target_scaling_nlo(2,:) = -half            ! coll limit
    target_scaling_nlo(3,:) = -half            ! softcoll limit
    
        
!    x(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    x = yRnd

    i  = 1
    r = 0
    ff(1) = zero

    icoll  = 3
    jother = 4
    
    limmat = zero
    limval_nlo = zero
    ! run without limits
    !    r = xsect_wbfj_n2loxlo_r(ndim,x,vegasweight)
    ! r = xsect_nloqcd_r_is_gq(xlim,ff,vegasweight)
    !call xsect_selector(xlim,ff,vegasweight,r)


    print *, "----- USING POINT # ", num_checks
    do limdepth = limdepth_min,limdepth_max,-1
       print *, "---------limdepth----------",limdepth

       ! soft limit
       if (verbose) print *, "soft limit"
       xlim = x
!       xlim(kNLO_min) = 10.0_dp**limdepth
       !       xlim(kNLO_min) = ( xlim(kNLO_min) - buff)/onet
       xlim = set_small_x(x,limdepth,kNLO_min)
       call xsect_selector(xlim,ff,vegasweight,rlim(1))

       if (verbose) print *, "limval", limval_nlo(1:number_terms_nlo)


       limmat(i,1:number_terms_nlo,1) = limval_nlo(1:number_terms_nlo)
!
       !       ! coll limit
       if (verbose) print *, "coll limit"
       xlim = x
!       xlim(kNLO_min+1) = 10.0_dp**limdepth
!       xlim(kNLO_min+1) = ( xlim(kNLO_min+1) - buff)/onet
       xlim = set_small_x(x,limdepth,kNLO_min+1)
       call xsect_selector(xlim,ff,vegasweight,rlim(2))

       if (verbose) print *, "limval", limval_nlo(1:number_terms_nlo)


       limmat(i,1:number_terms_nlo,2) = limval_nlo(1:number_terms_nlo)
!
       !       ! soft-coll limit
       if (verbose) print *, "softcoll limit"
       xlim = x
!       xlim(kNLO_min) = 10.0_dp**limdepth
!       xlim(kNLO_min+1) = 10.0_dp**limdepth
!       xlim(kNLO_min) = ( xlim(kNLO_min) - buff)/onet
       !       xlim(kNLO_min+1) = ( xlim(kNLO_min+1) - buff)/onet
       xlim = set_small_x(x,limdepth,kNLO_min,kNLO_min+1)

       call xsect_selector(xlim,ff,vegasweight,rlim(3))

       if (verbose) print *, "limval", limval_nlo(1:number_terms_nlo)
       
       limmat(i,1:number_terms_nlo,3) = limval_nlo(1:number_terms_nlo)

!       print *, "soft limit",limmat(i,:,1)
!       print *, "coll limit",limmat(i,:,2)
!       print *, "softcoll limit",limmat(i,:,3)
!       print *, "r",r

       !if (product(limmat(i,:,1))*product(limmat(i,:,2))*product(limmat(i,:,3))*product(limmat(i,:,4)) .eq. zero) then
!       if (product(rlim) .eq. zero) then
!          return
!       else
!
!
       ! check cancellation
          do ilim = 1, number_limits_nlo
             ssign = +one
             if (ilim .eq. 3) ssign = -one  ! SC
             do iterm = 1, number_terms_nlo/2
                i1 =  cancelling_terms_nlo(ilim,iterm,1)
                i2 =  cancelling_terms_nlo(ilim,iterm,2)
                if (.not. sing_terms_nlo(ilim,iterm)) cycle    ! don't check if there is no sing (sub)limit

!                cancelling_term_index = cancelling_terms_nlo(iterm,ilim)
                if (limmat(i,i1,ilim) .ne. zero) then
!                   print *, limmat(i,i1,ilim),limmat(i,i2,ilim)
                   limcancellation(i,ilim,iterm) = (limmat(i,i1,ilim) + ssign*limmat(i,i2,ilim))/limmat(i,i1,ilim)
                   loglimcancellation(i,ilim,iterm) = log(abs(limcancellation(i,ilim,iterm)))/log(10.0_dp)
                                   
!                   print *, loglimcancellation(i,ilim,iterm)
                else
                   limcancellation(i,ilim,iterm) = (limmat(i,i1,ilim) + ssign*limmat(i,i2,ilim))
                   loglimcancellation(i,ilim,iterm) = -1000.0_dp
                                   
!                   print *, ilim, iterm, (limmat(i,iterm,ilim) + ssign*limmat(i,cancelling_term_index,ilim))
                endif
                
                print *, i,ilim,iterm,limcancellation(i,ilim,iterm),loglimcancellation(i,ilim,iterm)
             enddo
          enddo
!       endif
                
       i = i+1
    enddo

    ifile = 1
    do ilim = 1,number_limits_nlo
       if (ilim .eq. 1) then
          print *, "SCALING SOFT LIMIT"
          
       elseif (ilim .eq. 2) then
          print *, "SCALING COLL LIMIT"
       elseif (ilim .eq. 3) then
          print *, "SCALING SOFTCOLL LIMIT"
       endif
          
       do iterm = 1, number_terms_nlo/2
          print *, iterm
          if (.not. sing_terms_nlo(ilim,iterm)) cycle
!          print *, ilim,iterm,loglimcancellation(:,ilim,iterm)
          i = 1
          do limdepth = limdepth_min,limdepth_max+1,-1
             !             write(101+ifile,*) limdepth,loglimcancellation(i,ilim,iterm)
             scaling(i,ilim,iterm) = loglimcancellation(i+1,ilim,iterm)-loglimcancellation(i,ilim,iterm)
             print *, limcancellation(i+1,ilim,iterm),limcancellation(i,ilim,iterm),scaling(i,ilim,iterm)
             i = i+1
          enddo
                          
          ifile = ifile+1
       enddo
       call evaluate_scaling(scaling(:,ilim,:),target_scaling_nlo(ilim,:),sing_terms_nlo(ilim,:))
    enddo


    
    print *, "num_checks", num_checks
    if (num_checks .eq. 5) stop

    num_checks = num_checks+1


        
    return
    
  end function check_lim_nlo


  
  integer function check_lim_nlo_is(yRnd,ff,vegasweight) result(r)
    ! this checks NLO sometimes used for initial state
    ! , i.e. (I-S_m) (I - C_1m - C_2m) + , so 1 hard term + 5 subtr terms
    real(dp15)  :: yRnd(30),vegasweight,ff(1),ff_dum(1)

    real(dp15)             :: x(30),xlim(30)
    integer                :: rlim(number_limits_nlo_is)
    
    real(dp)               :: limmat(limdepth_min-limdepth_max+1,number_terms_nlo_is,number_limits_nlo_is),ssign
    real(dp)               :: limcancellation(limdepth_min-limdepth_max+1,number_limits_nlo_is,number_terms_nlo_is/2),loglimcancellation(limdepth_min-limdepth_max+1,number_limits_nlo_is,number_terms_nlo_is/2)
    real(dp)               :: target_scaling_nlo(number_limits_nlo_is,number_terms_nlo_is/2),scaling(limdepth_min-limdepth_max+1,number_limits_nlo_is,number_terms_nlo_is/2)
    integer                :: limdepth, i, ilim,iterm,i1,i2,ifile
    integer                :: cancelling_terms_nlo(number_limits_nlo_is,number_terms_nlo_is/2,2)
    logical                :: sing_terms_nlo(number_limits_nlo_is,number_terms_nlo_is/2)
    integer, save          :: num_checks = 1


    call set_sing_terms_nlo_is(sing_terms_nlo)

    cancelling_terms_nlo(1,1,1:2)=[2,1]       ! first  pair for soft limit
    cancelling_terms_nlo(1,2,1:2)=[5,3]       ! second pair for soft limit
    cancelling_terms_nlo(1,3,1:2)=[6,4]       ! ... 

    cancelling_terms_nlo(2,1,1:2)=[3,1]       ! first  pair for coll1 limit
    cancelling_terms_nlo(2,2,1:2)=[5,2]       ! second pair for coll1 limit
    cancelling_terms_nlo(2,3,1:2)=[6,4]       ! dummy!
    sing_terms_nlo(2,3) = .false.

    cancelling_terms_nlo(3,1,1:2)=[4,1]       ! first  pair for coll2 limit
    cancelling_terms_nlo(3,2,1:2)=[6,2]       ! second pair for coll2 limit
    cancelling_terms_nlo(3,3,1:2)=[5,3]       ! dummy
    sing_terms_nlo(3,3) = .false.

    cancelling_terms_nlo(4,1,1:2)=[5,1]       ! first  pair for softcoll1 limit
    cancelling_terms_nlo(4,2,1:2)=[3,2]       ! second pair for softcoll1 limit
    cancelling_terms_nlo(4,3,1:2)=[6,4]       ! dummy!
!    sing_terms_nlo(4,3) = .false.
    
    cancelling_terms_nlo(5,1,1:2)=[6,1]       ! first  pair for softcoll2 limit
    cancelling_terms_nlo(5,2,1:2)=[4,2]       ! second pair for softcoll2 limit
    cancelling_terms_nlo(5,3,1:2)=[5,3]       ! dummy
    !    sing_terms_nlo(5,3) = .false.


    target_scaling_nlo(1,:) = -one            ! soft limit
    target_scaling_nlo(2,:) = -half            ! coll limit
    target_scaling_nlo(3,:) = -half            ! softcoll limit
    target_scaling_nlo(4,1:2) = -half            ! coll limit
    target_scaling_nlo(5,1:2) = -half            ! softcoll limit
    target_scaling_nlo(4,3) = -one            ! softcoll limit, on "other" IS parton --> soft  
    target_scaling_nlo(5,3) = -one            ! softcoll limit, on "other" IS parton --> soft

    x = yRnd
    i  = 1
    r = 0
    ff(1) = zero
    
    limmat = zero
    limval_nlo_is = zero
    ! run without limits
    ! todo call correct function here
    !call xsect_selector(x,ff,vegasweight,r)
        

    print *, "----- USING POINT # ", num_checks
    do limdepth = limdepth_min,limdepth_max,-1
       print *, "---------limdepth----------",limdepth
       if (verbose) print *, "soft limit"
       xlim = x
       xlim(kNLO_min) = 10.0_dp**limdepth
       xlim(kNLO_min) = ( xlim(kNLO_min) - buff)/onet
       call xsect_selector(xlim,ff,vegasweight,rlim(1))
              
       limmat(i,1:number_terms_nlo_is,1) = limval_nlo_is(1:number_terms_nlo_is)
       if (verbose) print *, "limval", limval_nlo_is(1:number_terms_nlo_is)
!
       !       ! coll1 limit
       if (verbose) print *, "coll1 limit"
       xlim = x
       xlim(kNLO_min+1) = 10.0_dp**limdepth
       xlim(kNLO_min+1) = ( xlim(kNLO_min+1) - buff)/onet

       call xsect_selector(xlim,ff,vegasweight,rlim(2))     


       limmat(i,1:number_terms_nlo_is,2) = limval_nlo_is(1:number_terms_nlo_is)
       if (verbose) print *, "limval", limval_nlo_is(1:number_terms_nlo_is)

!       ! coll2 limit
       xlim = x
       if (verbose) print *, "coll2 limit"
       xlim(kNLO_min+1) = one - 10.0_dp**limdepth
       xlim(kNLO_min+1) = ( xlim(kNLO_min+1) - buff)/onet

       call xsect_selector(xlim,ff,vegasweight,rlim(3))

       limmat(i,1:number_terms_nlo_is,3) = limval_nlo_is(1:number_terms_nlo_is)
       if (verbose) print *, "limval", limval_nlo_is(1:number_terms_nlo_is)
!
       !       ! soft-coll1 limit
       if (verbose) print *, "softcoll1 limit"
       xlim = x
       xlim(kNLO_min) = 10.0_dp**limdepth
       xlim(kNLO_min+1) = 10.0_dp**limdepth
       xlim(kNLO_min) = ( xlim(kNLO_min) - buff)/onet
       xlim(kNLO_min+1) = ( xlim(kNLO_min+1) - buff)/onet

       call xsect_selector(xlim,ff,vegasweight,rlim(4))
       
       limmat(i,1:number_terms_nlo_is,4) = limval_nlo_is(1:number_terms_nlo_is)
       if (verbose) print *, "limval", limval_nlo_is(1:number_terms_nlo_is)

       !       ! soft-coll2 limit
       if (verbose) print *, "softcoll2 limit"
       xlim = x
       xlim(kNLO_min) = 10.0_dp**limdepth
       xlim(kNLO_min+1) = one - 10.0_dp**limdepth
       xlim(kNLO_min) = ( xlim(kNLO_min) - buff)/onet
       xlim(kNLO_min+1) = ( xlim(kNLO_min+1) - buff)/onet

       call xsect_selector(xlim,ff,vegasweight,rlim(5))

       limmat(i,1:number_terms_nlo_is,5) = limval_nlo_is(1:number_terms_nlo_is)
       if (verbose) print *, "limval", limval_nlo_is(1:number_terms_nlo_is)
       

 !      print *, "soft limit",limmat(i,:,1)
 !      print *, "coll limit",limmat(i,:,2)
 !      print *, "softcoll limit",limmat(i,:,3)
 !      print *, "r",r

       !if (product(limmat(i,:,1))*product(limmat(i,:,2))*product(limmat(i,:,3))*product(limmat(i,:,4)) .eq. zero) then
!       if (product(limval_nlo_is) .eq. zero) then
!          print *, 'returning'
       !   return
!       else

!
       ! check cancellation
       do ilim = 1, number_limits_nlo_is
          if (verbose) print *, "limit",ilim
             ssign = +one
!             if (ilim .eq. 5 .or. ilim .eq. 6) ssign = -one  ! SC
             do iterm = 1, number_terms_nlo_is/2
                ssign = +one
                if ((ilim .eq. 4 .or. ilim .eq. 5) .and. (iterm .ne. number_terms_nlo_is/2)) ssign = -one  ! SC, but not for one coll limit!
                i1 =  cancelling_terms_nlo(ilim,iterm,1)
                i2 =  cancelling_terms_nlo(ilim,iterm,2)
                if (.not. sing_terms_nlo(ilim,iterm)) cycle    ! don't check if there is no sing (sub)limit

!                cancelling_term_index = cancelling_terms_nlo_is(iterm,ilim)
                if (limmat(i,i1,ilim) .ne. zero) then
                   limcancellation(i,ilim,iterm) = (limmat(i,i1,ilim) + ssign*limmat(i,i2,ilim))/limmat(i,i1,ilim)
                   loglimcancellation(i,ilim,iterm) = log(abs(limcancellation(i,ilim,iterm)))/log(10.0_dp)
                                   
!                   print *, ilim, iterm, (limmat(i,iterm,ilim) + ssign*limmat(i,cancelling_term_index,ilim))/limmat(i,iterm,ilim)
                else
!                   print *, limmat(i,i1,ilim),limmat(i,i2,ilim)
                   limcancellation(i,ilim,iterm) = (limmat(i,i1,ilim) + ssign*limmat(i,i2,ilim))
                   loglimcancellation(i,ilim,iterm) = -1000.0_dp
                                   
 !                  print *, ilim, iterm, (limmat(i,iterm,ilim) + ssign*limmat(i,cancelling_term_index,ilim))
                endif
                
!               print *, i,ilim,iterm,limcancellation(i,ilim,iterm),loglimcancellation(i,ilim,iterm)
             enddo
          enddo
!       endif
                
       i = i+1
    enddo
!    pause

    ifile = 1
    do ilim = 1,number_limits_nlo_is
       if (ilim .eq. 1) then
          print *, "SCALING SOFT LIMIT"
          
       elseif (ilim .eq. 2) then
          print *, "SCALING COLL1 LIMIT"
       elseif (ilim .eq. 3) then
          print *, "SCALING COLL2 LIMIT"
       elseif (ilim .eq. 4) then
          print *, "SCALING SOFTCOLL1 LIMIT"
       elseif (ilim .eq. 5) then
          print *, "SCALING SOFTCOLL2 LIMIT"
       endif
          
       do iterm = 1, number_terms_nlo_is/2
          if (.not. sing_terms_nlo(ilim,iterm)) cycle
                    
!          print *, ilim,iterm,loglimcancellation(:,ilim,iterm)
          i = 1
          do limdepth = limdepth_min,limdepth_max+2,-1
             scaling(i+1,ilim,iterm) = loglimcancellation(i+1,ilim,iterm)-loglimcancellation(i,ilim,iterm) 
             !print *, limcancellation(i+1,ilim,iterm),limcancellation(i,ilim,iterm),scaling(i+1,ilim,iterm)
             !print *, loglimcancellation(i+1,ilim,iterm)-loglimcancellation(i,ilim,iterm)
             !print *, limcancellation(i,ilim,iterm)
             !             write(101+ifile,*) limdepth,loglimcancellation(i,ilim,iterm)
             i = i+1
          enddo
                          
          ifile = ifile+1
       enddo
       call evaluate_scaling(scaling(:,ilim,:),target_scaling_nlo(ilim,:),sing_terms_nlo(ilim,:))
    enddo
    
    print *, "num_checks sdc", num_checks
    if (num_checks .eq. 20) stop

    num_checks = num_checks+1
        
    return
    
  end function check_lim_nlo_is
  


  integer function check_lim_nnlo(yRnd,ff,vegasweight) result(r)
!    integer, intent(in)    :: ndim
    real(dp15)              :: yRnd(30), vegasweight,ff(1)
    real(dp15)             :: xlim(30),x(30)
    integer                :: rlim(number_limits_nnlo)
    
    real(dp)               :: limmat(limdepth_min-limdepth_max+1,number_terms_nnlo,number_limits_nnlo),ssign(number_limits_nnlo)
    real(dp)               :: limcancellation(limdepth_min-limdepth_max+1,number_limits_nnlo,number_terms_nnlo/2),loglimcancellation(limdepth_min-limdepth_max+1,number_limits_nnlo,number_terms_nnlo/2)
    
    integer                :: limdepth, i, ilim,iterm,i1,i2,ifile
    integer                :: cancelling_terms_nnlo(number_limits_nnlo,number_terms_nnlo/2,2)
    logical                :: sing_terms_nnlo(number_limits_nnlo,number_terms_nnlo/2)
    integer, save          :: num_checks = 1

    cancelling_terms_nnlo(1,1,1:2)=[2,1]       ! first  pair for x1->0 limit
    cancelling_terms_nnlo(1,2,1:2)=[6,3]       ! second pair for x1->0 limit
    cancelling_terms_nnlo(1,3,1:2)=[7,4]       ! ...
    cancelling_terms_nnlo(1,4,1:2)=[8,5]       ! ...
    cancelling_terms_nnlo(1,5,1:2)=[12,9]      ! ...
    cancelling_terms_nnlo(1,6,1:2)=[13,10]     ! ...
    cancelling_terms_nnlo(1,7,1:2)=[14,11]     ! ...
    cancelling_terms_nnlo(1,8,1:2)=[16,15]     ! ...


    cancelling_terms_nnlo(2,1,1:2)=[3,1]       ! first  pair for x2->0 limit
    cancelling_terms_nnlo(2,2,1:2)=[6,2]       ! second pair for x2->0 limit
    cancelling_terms_nnlo(2,3,1:2)=[9,4]       ! ...
    cancelling_terms_nnlo(2,4,1:2)=[10,5]      ! ...
    cancelling_terms_nnlo(2,5,1:2)=[12,7]      ! ...
    cancelling_terms_nnlo(2,6,1:2)=[13,8]      ! ...
    cancelling_terms_nnlo(2,7,1:2)=[15,11]     ! ...
    cancelling_terms_nnlo(2,8,1:2)=[16,14]     ! ...
    

    cancelling_terms_nnlo(3,1,1:2)=[4,1]       ! first  pair for x3->0 limit
    cancelling_terms_nnlo(3,2,1:2)=[7,2]       ! second pair for x3->0 limit
    cancelling_terms_nnlo(3,3,1:2)=[9,3]       ! ...
    cancelling_terms_nnlo(3,4,1:2)=[11,5]      ! ...
    cancelling_terms_nnlo(3,5,1:2)=[12,6]      ! ...
    cancelling_terms_nnlo(3,6,1:2)=[14,8]      ! ...
    cancelling_terms_nnlo(3,7,1:2)=[15,10]     ! ...
    cancelling_terms_nnlo(3,8,1:2)=[16,13]     ! ...


    cancelling_terms_nnlo(4,1,1:2)=[5,1]       ! first  pair for x4->0 limit
    cancelling_terms_nnlo(4,2,1:2)=[8,2]       ! second pair for x4->0 limit
    cancelling_terms_nnlo(4,3,1:2)=[10,3]      ! ...
    cancelling_terms_nnlo(4,4,1:2)=[11,4]      ! ...
    cancelling_terms_nnlo(4,5,1:2)=[13,6]      ! ...
    cancelling_terms_nnlo(4,6,1:2)=[14,7]      ! ...
    cancelling_terms_nnlo(4,7,1:2)=[15,9]      ! ...
    cancelling_terms_nnlo(4,8,1:2)=[16,12]     ! ...


    cancelling_terms_nnlo(5,1,1:2)=[6,1]       ! first  pair for x1,x2->0 limit
    cancelling_terms_nnlo(5,2,1:2)=[3,2]       ! second pair for x1,x2->0 limit
    cancelling_terms_nnlo(5,3,1:2)=[12,4]      ! ...
    cancelling_terms_nnlo(5,4,1:2)=[13,5]      ! ...
    cancelling_terms_nnlo(5,5,1:2)=[10,8]      ! ...
    cancelling_terms_nnlo(5,6,1:2)=[9,7]       ! ...
    cancelling_terms_nnlo(5,7,1:2)=[16,11]     ! ...
    cancelling_terms_nnlo(5,8,1:2)=[15,14]     ! ...


    cancelling_terms_nnlo(6,1,1:2)=[7,1]       ! first  pair for x1,x3->0 limit
    cancelling_terms_nnlo(6,2,1:2)=[4,2]       ! second pair for x1,x3->0 limit
    cancelling_terms_nnlo(6,3,1:2)=[12,3]      ! ...
    cancelling_terms_nnlo(6,4,1:2)=[14,5]      ! ...
    cancelling_terms_nnlo(6,5,1:2)=[9,6]       ! ...
    cancelling_terms_nnlo(6,6,1:2)=[11,8]      ! ...
    cancelling_terms_nnlo(6,7,1:2)=[16,10]     ! ...
    cancelling_terms_nnlo(6,8,1:2)=[15,13]     ! ...

    cancelling_terms_nnlo(7,1,1:2)=[8,1]       ! first  pair for x1,x4->0 limit
    cancelling_terms_nnlo(7,2,1:2)=[5,2]       ! second pair for x1,x4->0 limit
    cancelling_terms_nnlo(7,3,1:2)=[13,3]      ! ...
    cancelling_terms_nnlo(7,4,1:2)=[14,4]      ! ...
    cancelling_terms_nnlo(7,5,1:2)=[10,6]      ! ...
    cancelling_terms_nnlo(7,6,1:2)=[11,7]      ! ...
    cancelling_terms_nnlo(7,7,1:2)=[16,9]      ! ...
    cancelling_terms_nnlo(7,8,1:2)=[15,12]     ! ...
    
    cancelling_terms_nnlo(8,1,1:2)=[9,1]       ! first  pair for x2,x3->0 limit
    cancelling_terms_nnlo(8,2,1:2)=[12,2]      ! second pair for x2,x3->0 limit
    cancelling_terms_nnlo(8,3,1:2)=[4,3]       ! ...
    cancelling_terms_nnlo(8,4,1:2)=[15,5]      ! ...
    cancelling_terms_nnlo(8,5,1:2)=[7,6]       ! ...
    cancelling_terms_nnlo(8,6,1:2)=[16,8]      ! ...
    cancelling_terms_nnlo(8,7,1:2)=[11,10]     ! ...
    cancelling_terms_nnlo(8,8,1:2)=[14,13]     ! ...

    cancelling_terms_nnlo(9,1,1:2)=[10,1]       ! first  pair for x2,x4->0 limit
    cancelling_terms_nnlo(9,2,1:2)=[13,2]       ! second pair for x2,x4->0 limit
    cancelling_terms_nnlo(9,3,1:2)=[5,3]        ! ...
    cancelling_terms_nnlo(9,4,1:2)=[15,4]       ! ...
    cancelling_terms_nnlo(9,5,1:2)=[8,6]        ! ...
    cancelling_terms_nnlo(9,6,1:2)=[16,7]       ! ...
    cancelling_terms_nnlo(9,7,1:2)=[11,9]       ! ...
    cancelling_terms_nnlo(9,8,1:2)=[14,12]      ! ...


    cancelling_terms_nnlo(10,1,1:2)=[11,1]       ! first  pair for x3,x4->0 limit
    cancelling_terms_nnlo(10,2,1:2)=[14,2]       ! second pair for x3,x4->0 limit
    cancelling_terms_nnlo(10,3,1:2)=[15,3]       ! ...
    cancelling_terms_nnlo(10,4,1:2)=[5,4]        ! ...
    cancelling_terms_nnlo(10,5,1:2)=[16,6]       ! ...
    cancelling_terms_nnlo(10,6,1:2)=[8,7]        ! ...
    cancelling_terms_nnlo(10,7,1:2)=[10,9]       ! ...
    cancelling_terms_nnlo(10,8,1:2)=[13,12]      ! ...

    cancelling_terms_nnlo(11,1,1:2)=[12,1]       ! first  pair for x1,x2,x3->0 limit
    cancelling_terms_nnlo(11,2,1:2)=[9,2]        ! second pair for x1,x2,x3->0 limit
    cancelling_terms_nnlo(11,3,1:2)=[7,3]        ! ...
    cancelling_terms_nnlo(11,4,1:2)=[6,4]        ! ...
    cancelling_terms_nnlo(11,5,1:2)=[16,5]       ! ...
    cancelling_terms_nnlo(11,6,1:2)=[15,8]       ! ...
    cancelling_terms_nnlo(11,7,1:2)=[14,10]      ! ...
    cancelling_terms_nnlo(11,8,1:2)=[13,11]      ! ...

    cancelling_terms_nnlo(12,1,1:2)=[13,1]       ! first  pair for x1,x2,x4->0 limit
    cancelling_terms_nnlo(12,2,1:2)=[10,2]       ! second pair for x1,x2,x4->0 limit
    cancelling_terms_nnlo(12,3,1:2)=[8,3]        ! ...
    cancelling_terms_nnlo(12,4,1:2)=[16,4]       ! ...
    cancelling_terms_nnlo(12,5,1:2)=[6,5]        ! ...
    cancelling_terms_nnlo(12,6,1:2)=[15,7]       ! ...
    cancelling_terms_nnlo(12,7,1:2)=[14,9]       ! ...
    cancelling_terms_nnlo(12,8,1:2)=[12,11]      ! ...

    cancelling_terms_nnlo(13,1,1:2)=[14,1]       ! first  pair for x1,x3,x4->0 limit
    cancelling_terms_nnlo(13,2,1:2)=[11,2]       ! second pair for x1,x3,x4->0 limit
    cancelling_terms_nnlo(13,3,1:2)=[16,3]       ! ...
    cancelling_terms_nnlo(13,4,1:2)=[8,4]        ! ...
    cancelling_terms_nnlo(13,5,1:2)=[7,5]        ! ...
    cancelling_terms_nnlo(13,6,1:2)=[15,6]       ! ...
    cancelling_terms_nnlo(13,7,1:2)=[13,9]       ! ...
    cancelling_terms_nnlo(13,8,1:2)=[12,10]      ! ...

    cancelling_terms_nnlo(14,1,1:2)=[15,1]       ! first  pair for x2,x3,x4->0 limit
    cancelling_terms_nnlo(14,2,1:2)=[16,2]       ! second pair for x2,x3,x4->0 limit
    cancelling_terms_nnlo(14,3,1:2)=[11,3]       ! ...
    cancelling_terms_nnlo(14,4,1:2)=[10,4]       ! ...
    cancelling_terms_nnlo(14,5,1:2)=[9,5]        ! ...
    cancelling_terms_nnlo(14,6,1:2)=[14,6]       ! ...
    cancelling_terms_nnlo(14,7,1:2)=[13,7]       ! ...
    cancelling_terms_nnlo(14,8,1:2)=[12,8]       ! ...

    cancelling_terms_nnlo(15,1,1:2)=[16,1]       ! first  pair for x1,x2,x3,x4->0 limit
    cancelling_terms_nnlo(15,2,1:2)=[15,2]       ! second pair for x1,x2,x3,x4->0 limit
    cancelling_terms_nnlo(15,3,1:2)=[14,3]       ! ...
    cancelling_terms_nnlo(15,4,1:2)=[13,4]       !  ...
    cancelling_terms_nnlo(15,5,1:2)=[12,5]       ! ...
    cancelling_terms_nnlo(15,6,1:2)=[11,6]       ! ...
    cancelling_terms_nnlo(15,7,1:2)=[10,7]       ! ...
    cancelling_terms_nnlo(15,8,1:2)=[9,8]        ! ...

    ssign = [one, one,one,one,-one,-one,-one,-one,-one,-one,one,one,one,one,-one]

    x = yRnd

    sing_terms_nnlo = .true.
    i = 1
    r = 0
    limmat = zero
    limval_nnlo = zero
    ! run without limits
    ! r = my_nnlo_function(ndim,x,vegasweight)


    print *, "----- USING POINT # ", num_checks
    do limdepth = limdepth_min,limdepth_max,-1
       do ilim = 1, number_limits_nnlo
          print *, "limit",ilim
          xlim = x

          select case (ilim)
          case(1)
             xlim = set_small_x(x,limdepth,xE5)
          case(2)
             xlim = set_small_x(x,limdepth,xE6)
          case(3)
             xlim = set_small_x(x,limdepth,xRHO5)
          case(4)
             xlim = set_small_x(x,limdepth,xRHO6)
          case(5)
             xlim = set_small_x(x,limdepth,xE5,xE6)
          case(6)
             xlim = set_small_x(x,limdepth,xE5,xRHO5)
          case(7)
             xlim = set_small_x(x,limdepth,xE5,xRHO6)
          case(8)
             xlim = set_small_x(x,limdepth,xE6,xRHO5)
          case(9)
             xlim = set_small_x(x,limdepth,xE6,xRHO6)
          case(10)
             xlim = set_small_x(x,limdepth,xRHO5,xRHO6)
          case(11)
             xlim = set_small_x(x,limdepth,xE5,xE6,xRHO5)
          case(12)
             xlim = set_small_x(x,limdepth,xE5,xE6,xRHO6)
          case(13)
             xlim = set_small_x(x,limdepth,xE5,xRHO5,xRHO6)
          case(14)
             xlim = set_small_x(x,limdepth,xE6,xRHO5,xRHO6)
          case(15)
             xlim = set_small_x(x,limdepth,xE5,xE6,xRHO5,xRHO6)
          end  select
          call xsect_selector(xlim,ff(1),vegasweight,rlim(ilim))
          limmat(i,1:number_terms_nnlo,ilim) = limval_nnlo(1:number_terms_nnlo)
          if (verbose) print *, limval_nnlo
! only proceed with check if all contributions are non-zero          
!          if (product(rlim) .eq. zero) then
!             return
!          else

! check cancellation            
!             do ilim = 1, number_limits_nnlo

          do iterm = 1, number_terms_nnlo/2
             print *, "term", iterm

             i1 =  cancelling_terms_nnlo(ilim,iterm,1)
             i2 =  cancelling_terms_nnlo(ilim,iterm,2)
             print *,i1,i2,limmat(i,i1,ilim),limmat(i,i2,ilim)
             if (.not. sing_terms_nnlo(ilim,iterm)) cycle    ! don't check if there is no sing (sub)limit
             
             if (abs(limmat(i,i1,ilim)) .gt. 1E-300_dp) then
                limcancellation(i,ilim,iterm) = (limmat(i,i1,ilim) + ssign(ilim)*limmat(i,i2,ilim))/limmat(i,i1,ilim)
                loglimcancellation(i,ilim,iterm) = log(abs(limcancellation(i,ilim,iterm)))/log(10.0_dp)
             else
                limcancellation(i,ilim,iterm) = (limmat(i,i1,ilim) + ssign(ilim)*limmat(i,i2,ilim))
                loglimcancellation(i,ilim,iterm) = -1000.0_dp
             endif
             print *, ilim,iterm, limcancellation(i,ilim,iterm)
!             pause
             
             !               print *, i,ilim,iterm,limcancellation(i,ilim,iterm),loglimcancellation(i,ilim,iterm)
          enddo


       enddo
!    endif
    
    i = i+1
 enddo

 
 do ilim = 1,number_limits_nnlo
    if (ilim .eq. 1) then
       print *, "SCALING SOFT m LIMIT"
    elseif (ilim .eq. 2) then
       print *, "SCALING SOFT n LIMIT"
    elseif (ilim .eq. 3) then
       print *, "SCALING COLL m LIMIT"
    elseif (ilim .eq. 4) then
          print *, "SCALING COLL n LIMIT"
       elseif (ilim .eq. 5) then
          print *, "SCALING SOFT m SOFT n LIMIT"
       elseif (ilim .eq. 6) then
          print *, "SCALING SOFT m COLL m LIMIT"
       elseif (ilim .eq. 7) then
          print *, "SCALING SOFT m COLL n LIMIT"
       elseif (ilim .eq. 8) then
          print *, "SCALING SOFT n COLL m LIMIT"
       elseif (ilim .eq. 9) then
          print *, "SCALING SOFT n COLL n LIMIT"
       elseif (ilim .eq. 10) then
          print *, "SCALING COLL m COLL n LIMIT"
       elseif (ilim .eq. 11) then
          print *, "SCALING SOFT m SOFT n COLL m LIMIT"
       elseif (ilim .eq. 12) then
          print *, "SCALING SOFT m SOFT n COLL n LIMIT"
       elseif (ilim .eq. 13) then
          print *, "SCALING SOFT m COLL m COLL n LIMIT"
       elseif (ilim .eq. 14) then
          print *, "SCALING  SOFT n COLL m COLL n LIMIT"
       elseif (ilim .eq. 15) then
          print *, "SCALING  SOFT n SOFT m COLL m COLL n LIMIT"
       endif
          
       do iterm = 1,  1 !number_terms_nnlo/2
          print *, iterm
!          print *, sing_terms_nlo(ilim,iterm)
          if (.not. sing_terms_nnlo(ilim,iterm)) cycle
                    
!          print *, ilim,iterm,loglimcancellation(:,ilim,iterm)
          i = 1
          do limdepth = limdepth_min,limdepth_max+2,-1
             print *, i
             print *, limcancellation(i+2,ilim,iterm),limcancellation(i,ilim,iterm),(loglimcancellation(i+2,ilim,iterm)-loglimcancellation(i,ilim,iterm))/two
             !print *, loglimcancellation(i+1,ilim,iterm)-loglimcancellation(i,ilim,iterm)
             !print *, limcancellation(i,ilim,iterm)
             !             write(101+ifile,*) limdepth,loglimcancellation(i,ilim,iterm)
             i = i+1
          enddo
          
          ifile = ifile+1
       enddo
       pause
    enddo

 
    print *, "num_checks", num_checks
    
    
    if (num_checks .eq. 20) stop

    num_checks = num_checks+1
       
    return


  end function check_lim_nnlo


 function set_small_x(x,limdepth,i1,i2,i3,i4) 
   real(dp15), intent(in)   :: x(30)
   integer, intent(in)    ::  limdepth
   integer, intent(in)    :: i1
   integer, intent(in),optional    :: i2,i3,i4
   real(dp15)  :: set_small_x(30)
   
   set_small_x = x
   set_small_x(i1) = 10.0_dp**limdepth
   set_small_x(i1) = ( set_small_x(i1) - buff)/onet
   
   if (present(i2)) then
      set_small_x(i2) = 10.0_dp**limdepth
      set_small_x(i2) = ( set_small_x(i2) - buff)/onet
      if (present(i3)) then
         set_small_x(i3) = 10.0_dp**limdepth
         set_small_x(i3) = ( set_small_x(i3) - buff)/onet
         if (present(i4)) then
            set_small_x(i4) = 10.0_dp**limdepth
            set_small_x(i4) = ( set_small_x(i4) - buff)/onet
         endif
      endif
   endif
   
 end function set_small_x


 subroutine xsect_selector(xlim,ff,vegasweight,r)
   real(dp15), intent(in)     :: xlim(30)
   integer, intent(out)       :: r
   real(dp15), intent(out)    :: vegasweight, ff(1)   

   if (corr .eq. 'nloqcd') then

      if (ch .eq. 'ns') then

         if (sec .eq. 'r_is') then
            r = xsect_nloqcd_r_is_ns(xlim,ff,vegasweight)
         else
            call err_unknown_sec()
         endif         
            
!      elseif (ch .eq. 'gq') then

      else
         call err_unknown_ch()
      endif
      
   elseif  (corr .eq. 'nloewk') then
      if (ch .eq. 'ns') then
      
         if (sec .eq. 'r_is') then
            r = xsect_nloewk_r_is_ns(xlim,ff,vegasweight)
         elseif (sec .eq. 'r_fs_53') then
            r = xsect_nloewk_r_fs_53_ns(xlim,ff,vegasweight)        
         else
            call err_unknown_sec()
         endif
      elseif (ch.eq.'aq') then
         if (sec.eq.'r_is') then
            r = xsect_nloewk_r_is_aq(xlim,ff,vegasweight)        
         else
            call err_unknown_sec()
         endif
      else
         call err_unknown_ch()
      endif
      
      
   elseif   (corr .eq. 'nnlo') then
      if (ch .eq.  'ns_ga') then
         if (sec .eq. 'rr_5161a') then
            r = xsect_nnlo_rr_5161a_ns_ga(xlim,ff,vegasweight)
         elseif (sec .eq. 'rr_5161c') then
            r = xsect_nnlo_rr_5161c_ns_ga(xlim,ff,vegasweight)
         elseif (sec .eq. 'rr_5262a') then
            r = xsect_nnlo_rr_5262a_ns_ga(xlim,ff,vegasweight)
         elseif (sec .eq. 'rr_5262c') then
            r = xsect_nnlo_rr_5262c_ns_ga(xlim,ff,vegasweight)
         else
            call err_unknown_sec()
         endif
      else
         call err_unknown_ch()
      endif
      
   else
      
      call err_unknown_corr()
   endif
   
 end subroutine xsect_selector


 subroutine set_sing_terms_nlo(sing_terms_nlo)    
   logical, intent(out)           :: sing_terms_nlo(number_limits_nlo,number_terms_nlo/2)

   sing_terms_nlo = .true.
   
   if (corr .eq. 'nloqcd') then
      if (ch .eq. 'gq') then
         if (sec .eq. 'r_is') then
            sing_terms_nlo(1,:) = .false.     ! no soft lim
            sing_terms_nlo(2,2) = .false.     ! no soft lim inside coll limit
            sing_terms_nlo(3,:) = .false.     ! no soft-coll lim
         endif
      endif
   endif
            
  
 end subroutine set_sing_terms_nlo

 subroutine set_sing_terms_nlo_is(sing_terms_nlo)    
   logical, intent(out)          :: sing_terms_nlo(number_limits_nlo_is,number_terms_nlo_is/2)


   sing_terms_nlo = .true.
   if (corr .eq. 'nloewk') then
      if (ch .eq. 'aq') then
         if (sec .eq. 'r_is') then
            sing_terms_nlo(1,:) = .false.     ! no soft lim
            sing_terms_nlo(2:3,2) = .false.     ! no soft lim inside coll limit
!            sing_terms_nlo(2:3,4:5) = .false.     ! no softcoll lim inside coll limit
            sing_terms_nlo(4:5,:) = .false.     ! no soft-coll lim
         endif
      endif
   elseif (corr .eq. 'nloqcd') then
      if (ch .eq. 'ns') then
         if (sec .eq. 'r_is') then
            sing_terms_nlo = .false.           ! only have C1 and C2
            sing_terms_nlo(2,1) = .true.
            sing_terms_nlo(3,1) = .true.
         endif

      endif
   endif

   
 end subroutine  set_sing_terms_nlo_is
   

 subroutine err_unknown_corr()

   print *, "error in limit check!"
   print *, "unknown correction:", corr
   stop

 end subroutine err_unknown_corr

  subroutine err_unknown_ch()

   print *, "error in limit check, corr = ", corr
   print *, "unknown channel:", ch   
   stop

 end subroutine err_unknown_ch

 subroutine err_unknown_sec()

   print *, "error in limit check, corr, sec  = ", corr, sec
   print *, "unknown sector:", sec
   stop

 end subroutine err_unknown_sec
   


  subroutine evaluate_scaling(scaling,target_scaling,sing_terms)
    real(dp), intent(in)   :: scaling(:,:),target_scaling(:)
    logical, intent(in)    :: sing_terms(:)
    integer                :: nlim,nterms,ilim,iterm,scaling_rating
    

    nterms = size(scaling,dim=2)
 
    if (nterms .ne. size(target_scaling,dim=1)) then
       print *, 'error', nterms, size(target_scaling,dim=1)
       stop
    endif
    do iterm = 1,nterms
       if (.not. sing_terms(iterm)) return
       call test_scaling(scaling(:,iterm), target_scaling(iterm),scaling_rating)
       print *, "scaling rating", scaling_rating
       select case(scaling_rating)
       case(5)
          print *, "scaling of term #", iterm, " is perfect"
       case(4)
          print *, "scaling of term #", iterm, " is excellent"
       case(3)
          print *, "scaling of term #", iterm, " is very good"
       case(2)
          print *, "scaling of term #", iterm, " is good"
       case(1)
          print *, "scaling of term #", iterm, " is ok"
       case(0)
          print *, "scaling of term #", iterm, " failed!"
          pause
       case default
          print *, "error getting scaling_rating"
          stop
       end select
    enddo


  end subroutine evaluate_scaling

  subroutine test_scaling(scaling,target_scaling,scaling_rating)
    real(dp), intent(in)   :: scaling(:),target_scaling
    integer, intent(out)  :: scaling_rating
    integer                :: ndepth, idepth,i
    real(dp)               :: tol1, tol2, tol3

    tol1 = 0.05_dp
    tol2 = 0.2_dp
    tol3 = 0.5_dp
    
    scaling_rating = -1
    ndepth = size(scaling,dim=1)
    scaling_rating  = 0
    
    do idepth = 1,ndepth-3
       if ( close_enough(scaling(idepth),target_scaling,tol1) .and.  close_enough(scaling(idepth+1),target_scaling,tol1) &
            .and. close_enough(scaling(idepth+2),target_scaling,tol1) .and. close_enough(scaling(idepth+3),target_scaling,tol1) ) then
          scaling_rating = 5
          return
       endif
    enddo

    do idepth = 1,ndepth-3
       if ( close_enough(scaling(idepth),target_scaling,tol2) .and.  close_enough(scaling(idepth+1),target_scaling,tol2) &
            .and. close_enough(scaling(idepth+2),target_scaling,tol2) .and. close_enough(scaling(idepth+3),target_scaling,tol2) ) then
          scaling_rating = 4
          return
       endif
    enddo
   
     do idepth = 1,ndepth-2
       if (close_enough(scaling(idepth),target_scaling,tol2) .and.  close_enough(scaling(idepth+1),target_scaling,tol2) &
            .and. close_enough(scaling(idepth+2),target_scaling,tol2)  ) then
          scaling_rating = 3
          return
       endif
    enddo
    
    do idepth = 1,ndepth-2
       if (close_enough(scaling(idepth),target_scaling,tol3) .and.  close_enough(scaling(idepth+1),target_scaling,tol3) &
            .and. close_enough(scaling(idepth+2),target_scaling,tol3)  ) then
          scaling_rating = 2
          return
       endif
    enddo

    do idepth = 1,ndepth-1
       if (close_enough(scaling(idepth),target_scaling,tol3) .and.  close_enough(scaling(idepth+1),target_scaling,tol3) ) then
          scaling_rating = 1
          return
       endif
    enddo

  end  subroutine test_scaling

  logical function close_enough(a,b,tol)
    real(dp), intent(in)   :: a,b,tol

    if ( abs((a-b)/a) .lt. abs(tol) ) then
       close_enough = .true.
    else
       close_enough = .false.
    end if

  end function close_enough




end module mod_check_lim
       
       

       


    

    

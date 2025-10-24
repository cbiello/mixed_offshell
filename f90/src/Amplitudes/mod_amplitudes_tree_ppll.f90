!-- amplitudes with extra jet: factor out gs**2 or ee**2 per extra gluon/photon emission
module mod_amplitudes_tree_ppll
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_auxfunctions
  use mod_coupl
  implicit none
  private

  !-- res_tree_emitted-particles_initial-state
  !-- q stands for both q and qb
  !-- res_tree   --> 0 -> q e- e+ qb
  !-- res_treeAA --> 0 -> e- a a e+
  
  public :: res_tree_qqb,res_treeAA_aa
  public :: res_tree_qqb_gen
  public :: res_tree_g_qqb_gen, res_tree_g_gq_gen, res_tree_g_qg_gen
  
  public :: res_tree_g_qqb,res_tree_a_qqb 
  public :: res_tree_a_aq,res_tree_a_qa
  public :: res_tree_g_gq,res_tree_g_qg
  public :: res_treeAA_a_aa

  public :: res_tree_ga_qqb 
  public :: res_tree_qqb_qqb,res_tree_qq_qq 
  public :: res_tree_ga_aq,res_tree_ga_qa
  public :: res_tree_ga_gq,res_tree_ga_qg
  public :: res_tree_ga_ag,res_tree_ga_ga 

  public :: master_amp_qgqb_llb
  
contains

  ! New from Raoul !

  ! conventions: 0 -> q(1) qb(2) l(3) lb(4)
  ! amp returned as -5:7 x -5:7 matrix
  ! used for NC, CC+, CC-

  

  !----------------------------------------------------------------------
  !-- 4-point amplitudes
  !----------------------------------------------------------------------

    subroutine res_tree_qqb_gen(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(-5:7,-5:7)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4),zb(4,4),coupl(1:3,-1:1,-1:1)
    complex(dp) :: amp
    real(dp) :: aa,bb
    integer  :: i1,i2,i3,i4,i,j,ii

    !-- 0 --> q(1)^- qb(2) [V --> l(3)^- lb(4)]
    !    amp(i1,i2,i3,i4) = two * za(i1,i3)*zb(i4,i2)/sprod(i1,i2)
        amp(i1,i2,i3,i4) = two * za(i1,i3)*zb(i4,i2)/sprod(i1,i2)

    call spinoru(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)

    
!    call get_coupl(sprod(1,2),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
    !         [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],coupl)

    call get_coupl_gen(sprod(1,2),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
         [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],coupl)

    res = zero
    
    aa = abs(amp(2,1,3,4))**2  ! qqb LL,RR + qbq LR, RL
    bb = abs(amp(1,2,3,4))**2  ! qqb LR, RL + qbq LL, RR

!    print *,  " d db LL",coupl(1,-1,-1,-1)
!    print *,  " d db LR",coupl(1,-1,-1,+1)
!    print *,  " d db RL",coupl(1,-1,+1,-1)
!    print *,  " d db RR",coupl(1,-1,+1,+1)
!
!
!    print *,  " u ub LL",coupl(2,-2,-1,-1)
!    print *,  " u ub LR",coupl(2,-2,-1,+1)
!    print *,  " u ub RL",coupl(2,-2,+1,-1)
!    print *,  " u ub RR",coupl(2,-2,+1,+1)

    do i = -5,5
       do j = -5,5

#if (_Vcharge == 0)          
          if ( i .eq. -j) then
             if (mod(abs(i),2) .eq. 0) ii=2
             if (mod(abs(i),2) .eq. 1) ii=1
          else
             cycle
          endif
#elif (_Vcharge == 1)
          if (i+j .eq. 1) then                                   ! this supposes a unit CKM matrix
             ii = 3
          else
             cycle
          endif
#elif (_Vcharge == -1)
          if (i+j .eq. -1) then                                   ! this supposes a unit CKM matrix
             ii = 3
          else
             cycle
          endif
#endif       
          if ( i .gt. 0 .and. j .lt. 0) then ! qqb
             res(i,j) =  aa * ( abs(coupl(ii,-1,-1))**2 + abs(coupl(ii,+1,+1))**2 ) + &
                  bb * ( abs(coupl(ii,-1,+1))**2 + abs(coupl(ii,+1,-1))**2 )
          elseif ( i .lt. 0 .and. j .gt. 0) then ! qbq
             res(i,j) =  bb * ( abs(coupl(ii,-1,-1))**2 + abs(coupl(ii,+1,+1))**2 ) + &
                  aa * ( abs(coupl(ii,-1,+1))**2 + abs(coupl(ii,+1,-1))**2 )
          else
             cycle
          endif
       enddo
    enddo

    print *, "eesq2",eesq2

    res = res * xn * eesq2 * aveqq


    return

  end subroutine res_tree_qqb_gen


  !----------------------------------------------------------------------
  !-- 5-point amplitudes
  !----------------------------------------------------------------------



  !-- generic amplitude for nlo QCD mission
  !-- iconf(5,:) are the required crossings from 0 -> q g qb l lb [e.g. qqb -> e-e+ g is 2,5,1,3,4]
  !-- ave is the averaging factor
  !-- always assume that leptons are 3 and 4, otherwise coupling is wrong
  !-- if need to move leptons around: see routine for QED emission
  subroutine res_tree_j_qcd_gen(p,iconf,ave,res)
    real(dp), intent(in)  :: p(:,:),ave
    integer, intent(in)   :: iconf(:,:)
    real(dp), intent(out) :: res(size(iconf,2),3)
    integer     :: i,hq,hg,hl
    real(dp)    :: sprod(5,5)
    complex(dp) :: za(5,5),zb(5,5),coupl(1:3,-1:1,-1:1),amp(-1:1,-1:1,-1:1)

    res = zero

    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)
    call get_coupl_gen(sprod(3,4),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
         [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],coupl)

    do i = 1,size(iconf,2)
       amp = master_amp_qgqb_llb(iconf(1,i),iconf(2,i),iconf(3,i),iconf(4,i),iconf(5,i),za,zb)

       do hl=-1,1,2
          do hg=-1,1,2
             do hq=-1,1,2
                res(i,:) = res(i,:) + abs(amp(hq,hg,hl)*coupl(:,hq,hl))**2
             enddo
          enddo
       enddo
       
    enddo
                    

    !-- eesq already in coupl
    res = ave * 8 * xn * Cf * eesq2 * res

  end subroutine res_tree_j_qcd_gen

  
  !-- qqb channel
  !-- q(p1) qb(p2) -> e-(p3) e+(p4) g(p5)
  !-- res(1,:) = q qb for dn and up
  !-- res(2,:) = qb q for dn and up
  subroutine res_tree_g_qqb_gen(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(-5:7,-5:7)
    integer               :: i,j,ii,jj
    real(dp)              :: res1(2,3)
    integer, parameter :: iconf(5,2) = reshape([2,5,1,3,4, 1,5,2,3,4],[5,2])

    call res_tree_j_qcd_gen(p,iconf,aveqq,res1)

    res = zero
    do i = -5,5
       do j = -5,5
#if (_Vcharge == 0)          
          if ( i .eq. -j) then
             if (mod(abs(i),2) .eq. 0) ii=2
             if (mod(abs(i),2) .eq. 1) ii=1
          else
             cycle
          endif
#elif (_Vcharge == 1)
          if (i+j .eq. 1) then                                   ! this supposes a unit CKM matrix
             ii = 3
          else
             cycle
          endif
#elif (_Vcharge == -1)
          if (i+j .eq. -1) then                                   ! this supposes a unit CKM matrix
             ii = 3
          else
             cycle
          endif
#endif       
!          if ( i .eq. -j) then
!             if (mod(abs(i),2) .eq. 0) jj=2
!             if (mod(abs(i),2) .eq. 1) jj=1
!             !          elseif(mod(abs(i*j),2) .eq. 1) then  ! one of i and j is even, the other odd -- this would be for any FC current
!             elseif (abs(i+j) .eq. 1) then                   ! this supposes a unit CKM matrix
!             ii=3
!          else
!             cycle
!          endif
          if (i .gt. 0 .and. j .lt. 0) then   ! qqb
             res(i,j) = res1(1,ii)
          elseif (i .lt. 0 .and. j .gt. 0) then ! qbq
             res(i,j) = res1(2,ii)
          endif
       enddo
    enddo


    
  end subroutine res_tree_g_qqb_gen

!  !-- qqb channel
!  !-- q(p1) qb(p2) -> e-(p3) e+(p4) γ(p5)
!  !-- res(1,:) = q qb for dn and up
!  !-- res(2,:) = qb q for dn and up
!  subroutine res_tree_a_qqb(p,res)
!    real(dp), intent(in)  :: p(:,:)
!    real(dp), intent(out) :: res(2,2)
!    integer, parameter :: iconf(5,2) = reshape([2,5,1,3,4, 1,5,2,3,4],[5,2])
!
!    call res_tree_j_qed(p,iconf,aveqq,res)
!
!  end subroutine res_tree_a_qqb

  !-- gq channel
  !-- g(p1) q(p2) -> e-(p3) e+(p4) q(p5)
  !-- res(1,:) = g qb for dn and up
  !-- res(2,:) = g q  for dn and up
  subroutine res_tree_g_gq_gen(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(-5:7,-5:7)
    real(dp)              :: res1(2,2)
    integer               :: i,j,ii,jj
    integer, parameter :: iconf(5,2) = reshape([2,1,5,3,4, 5,1,2,3,4],[5,2])

    call res_tree_j_qcd(p,iconf,aveqg,res1)


    res = zero
    do i = -5,5
       do j = -5,5
#if (_Vcharge == 0)                   
          jj = quark_type(j)
#else
          jj = 3
#endif          
          if (i .eq. 0 .and. j .gt. 0) then   ! gq
             res(i,j) = res1(2,jj)
          elseif (i .eq. 0 .and. j .lt. 0) then ! gqb
             res(i,j) = res1(1,jj)
          endif
       enddo
    enddo

  end subroutine res_tree_g_gq_gen

  !-- qg channel
  !-- q(p1) g(p2) -> e-(p3) e+(p4) q(p5)
  !-- res(1,:) = q  g for dn and up
  !-- res(2,:) = qb g  for dn and up
  subroutine res_tree_g_qg_gen(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(-5:7,-5:7)
    real(dp)              :: res1(2,2)
    integer               :: i,j,ii,jj
    integer, parameter :: iconf(5,2) = reshape([5,2,1,3,4, 1,2,5,3,4],[5,2])

    call res_tree_j_qcd(p,iconf,aveqg,res1)

    res = zero
    do i = -5,5
       do j = -5,5
#if (_Vcharge == 0)          
          ii = quark_type(i)
#else
          ii = 3
#endif
          
          if (i .gt. 0 .and. j .eq. 0) then   ! qg
             res(i,j) = res1(1,ii)
          elseif (i .lt. 0 .and. j .eq. 0) then ! qbg
             res(i,j) = res1(2,ii)
          endif
       enddo
    enddo
    
  end subroutine res_tree_g_qg_gen
  
!  !-- aq channel
!  !-- as before, but with g <-> a
!  subroutine res_tree_a_aq(p,res)
!    real(dp), intent(in)  :: p(:,:)
!    real(dp), intent(out) :: res(2,2)
!    integer, parameter :: iconf(5,2) = reshape([2,1,5,3,4, 5,1,2,3,4],[5,2])
!
!    call res_tree_j_qed(p,iconf,aveqa,res)
!
!  end subroutine res_tree_a_aq
!
!  !-- qa channel
!  !-- as above, but with g <-> a
!  subroutine res_tree_a_qa(p,res)
!    real(dp), intent(in)  :: p(:,:)
!    real(dp), intent(out) :: res(2,2)
!    integer, parameter :: iconf(5,2) = reshape([5,2,1,3,4, 1,2,5,3,4],[5,2])
!
!    call res_tree_j_qed(p,iconf,aveqa,res)
!
!  end subroutine res_tree_a_qa


  !----------------------------------------------------------------------
  !-- 4-point amplitudes
  !----------------------------------------------------------------------
  
  !-- res(1,1) = d dx -> l lb
  !-- res(1,2) = u ux -> l lb
  !-- res(2,1) = dx d -> l lb
  !-- res(2,2) = ux u -> l lb
  subroutine res_tree_qqb(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4),zb(4,4),coupl(1:2,-1:1,-1:1)
    complex(dp) :: amp
    real(dp) :: aa,bb,cc(1:2),dd(1:2)
    integer  :: i1,i2,i3,i4

    !-- 0 --> q(1)^- qb(2) [V --> l(3)^- lb(4)]
    amp(i1,i2,i3,i4) = two * za(i1,i3)*zb(i4,i2)/sprod(i1,i2)

    call spinoru(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)
    
    call get_coupl(sprod(1,2),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
         [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],coupl)

    aa = abs(amp(2,1,3,4))**2
    bb = abs(amp(1,2,3,4))**2
    cc(:) = abs(coupl(:,-1,-1))**2 + abs(coupl(:,+1,+1))**2
    dd(:) = abs(coupl(:,-1,+1))**2 + abs(coupl(:,+1,-1))**2

    res(1,:) = aa * cc(:) + bb * dd(:)
    res(2,:) = aa * dd(:) + bb * cc(:)

    res = res * xn * eesq2 * aveqq

    return

  end subroutine res_tree_qqb

  !-- a a -> l lb
  subroutine res_treeAA_aa(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res
    complex(dp) :: za(4,4),zb(4,4)
    real(dp)    :: sprod(4,4),me2

    call spinoru(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)

    call me2_qaaqb(3,1,2,4,za,zb,me2)

    res = me2 * aveaa * Q_lep**4 * eesq2

  end subroutine res_treeAA_aa

  !-- amplitude for 0 -> em[p1] a[p2] a[p3] ep[p4]
  !-- factor out Q_lep**4 * eesq**2
  subroutine me2_qaaqb(i1,i2,i3,i4,za,zb,res)
    integer, intent(in) :: i1,i2,i3,i4
    complex(dp), intent(in) :: za(4,4),zb(4,4)
    real(dp), intent(out) :: res
    complex(dp) :: amp(-1:1,-1:1)
    integer :: he,ha

    amp(-1,+1) = two*zb(i4,i2)**2/zb(i3,i1)/zb(i4,i3) !-- a[- + - +]
    amp(-1,-1) = two*zb(i4,i3)**2/zb(i2,i1)/zb(i4,i2) !-- a[- - + +]
    amp(+1,+1) = two*za(i3,i4)**2/za(i1,i2)/za(i2,i4) !-- a[+ + - -]
    amp(+1,-1) = two*za(i2,i4)**2/za(i1,i3)/za(i3,i4) !-- a[+ - + -]

    res = zero

    do ha=-1,1,2
       do he=-1,1,2
          res = res + real(amp(he,ha)*conjg(amp(he,ha)),dp)
       enddo
    enddo

  end subroutine me2_qaaqb


  !----------------------------------------------------------------------
  !-- 5-point amplitudes
  !----------------------------------------------------------------------

  !-- qqb channel
  !-- q(p1) qb(p2) -> e-(p3) e+(p4) g(p5)
  !-- res(1,:) = q qb for dn and up
  !-- res(2,:) = qb q for dn and up
  subroutine res_tree_g_qqb(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    integer, parameter :: iconf(5,2) = reshape([2,5,1,3,4, 1,5,2,3,4],[5,2])

    call res_tree_j_qcd(p,iconf,aveqq,res)

  end subroutine res_tree_g_qqb

  !-- qqb channel
  !-- q(p1) qb(p2) -> e-(p3) e+(p4) γ(p5)
  !-- res(1,:) = q qb for dn and up
  !-- res(2,:) = qb q for dn and up
  subroutine res_tree_a_qqb(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    integer, parameter :: iconf(5,2) = reshape([2,5,1,3,4, 1,5,2,3,4],[5,2])

    call res_tree_j_qed(p,iconf,aveqq,res)

  end subroutine res_tree_a_qqb

  !-- gq channel
  !-- g(p1) q(p2) -> e-(p3) e+(p4) q(p5)
  !-- res(1,:) = g qb for dn and up
  !-- res(2,:) = g q  for dn and up
  subroutine res_tree_g_gq(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    integer, parameter :: iconf(5,2) = reshape([2,1,5,3,4, 5,1,2,3,4],[5,2])

    call res_tree_j_qcd(p,iconf,aveqg,res)

  end subroutine res_tree_g_gq

  !-- qg channel
  !-- q(p1) g(p2) -> e-(p3) e+(p4) q(p5)
  !-- res(1,:) = q  g for dn and up
  !-- res(2,:) = qb g  for dn and up
  subroutine res_tree_g_qg(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    integer, parameter :: iconf(5,2) = reshape([5,2,1,3,4, 1,2,5,3,4],[5,2])

    call res_tree_j_qcd(p,iconf,aveqg,res)

  end subroutine res_tree_g_qg

  !-- aq channel
  !-- as before, but with g <-> a
  subroutine res_tree_a_aq(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    integer, parameter :: iconf(5,2) = reshape([2,1,5,3,4, 5,1,2,3,4],[5,2])

    call res_tree_j_qed(p,iconf,aveqa,res)

  end subroutine res_tree_a_aq

  !-- qa channel
  !-- as above, but with g <-> a
  subroutine res_tree_a_qa(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    integer, parameter :: iconf(5,2) = reshape([5,2,1,3,4, 1,2,5,3,4],[5,2])

    call res_tree_j_qed(p,iconf,aveqa,res)

  end subroutine res_tree_a_qa
  
  !-- master amplitudes below
  
  !-- generic amplitude for nlo QCD mission
  !-- iconf(5,:) are the required crossings from 0 -> q g qb l lb [e.g. qqb -> e-e+ g is 2,5,1,3,4]
  !-- ave is the averaging factor
  !-- always assume that leptons are 3 and 4, otherwise coupling is wrong
  !-- if need to move leptons around: see routine for QED emission
  subroutine res_tree_j_qcd(p,iconf,ave,res)
    real(dp), intent(in)  :: p(:,:),ave
    integer, intent(in)   :: iconf(:,:)
    real(dp), intent(out) :: res(size(iconf,2),2)
    integer     :: i,hq,hg,hl
    real(dp)    :: sprod(5,5)
    complex(dp) :: za(5,5),zb(5,5),coupl(1:2,-1:1,-1:1),amp(-1:1,-1:1,-1:1)

    res = zero

    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)
    call get_coupl(sprod(3,4),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
         [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],coupl)

    do i = 1,size(iconf,2)

       amp = master_amp_qgqb_llb(iconf(1,i),iconf(2,i),iconf(3,i),iconf(4,i),iconf(5,i),za,zb)

       do hl=-1,1,2
          do hg=-1,1,2
             do hq=-1,1,2
                res(i,:) = res(i,:) + abs(amp(hq,hg,hl)*coupl(:,hq,hl))**2
             enddo
          enddo
       enddo
       
    enddo

    !-- eesq already in coupl
    res = ave * 8 * xn * Cf * eesq2 * res

  end subroutine res_tree_j_qcd

  !-- generic amplitude for nlo QED mission
  !-- iconf(5,:) are the required crossings from 0 -> q a qb l lb [e.g. qqb -> e-e+ a is 2,5,1,3,4]
  !-- ave is the averaging factor
  subroutine res_tree_j_qed(p,iconf,ave,res)
    real(dp), intent(in)  :: p(:,:),ave
    integer, intent(in)   :: iconf(:,:)
    real(dp), intent(out) :: res(size(iconf,2),2)
    integer     :: i,hq,ha,hl
    real(dp)    :: sprod(5,5)
    complex(dp) :: za(5,5),zb(5,5),coupl_is(1:2,-1:1,-1:1),coupl_fs(1:2,-1:1,-1:1)
    complex(dp) :: amp_is(-1:1,-1:1,-1:1),amp_fs(-1:1,-1:1,-1:1),amp(2)
    integer :: ismin,ismax,fsmin,fsmax
    logical :: need

    res = zero
    ismin = -1; ismax = -1
    fsmin = -1; fsmax = -1

    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)

    do i = 1,size(iconf,2)

       !-- initial-state emission
       amp_is = master_amp_qgqb_llb(iconf(1,i),iconf(2,i),iconf(3,i),iconf(4,i),iconf(5,i),za,zb)
       call need_coupl(iconf(4,i),iconf(5,i),ismin,ismax,need)
       if (need) call get_coupl(sprod(ismin,ismax),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
         [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],coupl_is)

       !-- final_state emission
       amp_fs = master_amp_qgqb_llb(iconf(4,i),iconf(2,i),iconf(5,i),iconf(1,i),iconf(3,i),za,zb)
       call need_coupl(iconf(1,i),iconf(3,i),fsmin,fsmax,need)
       if (need) call get_coupl(sprod(fsmin,fsmax),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
         [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],coupl_fs)

       do hl=-1,1,2
          do ha=-1,1,2
             do hq=-1,1,2

                amp = [Qdn,Qup]*amp_is(hq,ha,hl)*coupl_is(:,hq,hl) + Q_lep*amp_fs(hl,ha,hq)*coupl_fs(:,hq,hl) 
                res(i,:) = res(i,:) + abs(amp)**2

             enddo
          enddo
       enddo
       
    enddo

    res = res * 8._dp * xn * ave * eesq2 ! * ee**2

  contains

    subroutine need_coupl(i1,i2,imin_ref,imax_ref,need)
      integer, intent(in)  :: i1,i2
      integer, intent(inout) :: imin_ref,imax_ref
      logical, intent(out) :: need
      integer :: imin,imax

      imin = min(i1,i2)
      imax = max(i1,i2)
      if (imin_ref.ne.imin .or. imax_ref.ne.imax) then 
         need = .true.
         imin_ref = imin
         imax_ref = imax
      else
         need = .false.
      endif
    end subroutine need_coupl

  end subroutine res_tree_j_qed

  !---

  !-- master amplitude for 0 -> q[p1] g[p2] qb[p3] l[p4] lb[p5]
  !-- helicity label: h1,h2,h4
  !-- amp[-h1,-h2,-h3] = -(amp[h1,h2,h3])^*
  !-- vertices are all gamma/sqrt2
  !-- propagator factor for the photon: 1/s45
  function master_amp_qgqb_llb(i1,i2,i3,i4,i5,za,zb) result(res)
    complex(dp) :: res(-1:1,-1:1,-1:1)
    integer, intent(in)     :: i1,i2,i3,i4,i5
    complex(dp), intent(in) :: za(5,5),zb(5,5)
    complex(dp) :: dena,denb

    dena = one/(za(i1,i2)*za(i2,i3)*za(i4,i5))
    denb = one/(zb(i1,i2)*zb(i2,i3)*zb(i5,i4))

    !-- basic amplitudes
    res(-1,+1,-1) = za(i1,i4)**2*dena
    res(-1,-1,-1) = zb(i5,i3)**2*denb

    !-- 4 <-> 5
    res(-1,+1,+1) = -za(i1,i5)**2*dena
    res(-1,-1,+1) = -zb(i4,i3)**2*denb

    !-- other helicities
    !-- careful, complex conjugate does not
    !-- work with crossing
    res(+1,-1,+1) = zb(i1,i4)**2*denb
    res(+1,+1,+1) = za(i5,i3)**2*dena
    
    res(+1,-1,-1) = -zb(i1,i5)**2*denb
    res(+1,+1,-1) = -za(i4,i3)**2*dena
    
  end function master_amp_qgqb_llb

  !-- routine for a[1] a[2] > e-[3] e+[4] a[5]
  subroutine res_treeAA_a_aa(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res
    complex(dp) :: za(5,5),zb(5,5)
    real(dp)    :: sprod(5,5)
    integer     :: i1,i2,i3,i4,i5
    real(dp)    :: me2

    !-- 0 -> em ep a a a, summed over pol, no averages
    me2(i1,i2,i3,i4,i5) = 16._dp * sprod(i1,i2) * ( &
         (sprod(i1,i3)**2 + sprod(i2,i3)**2)/(sprod(i1,i4)*sprod(i2,i4)*sprod(i1,i5)*sprod(i2,i5)) + &
         (sprod(i1,i4)**2 + sprod(i2,i4)**2)/(sprod(i1,i3)*sprod(i2,i3)*sprod(i1,i5)*sprod(i2,i5)) + &
         (sprod(i1,i5)**2 + sprod(i2,i5)**2)/(sprod(i1,i3)*sprod(i2,i3)*sprod(i1,i4)*sprod(i2,i4)) ) ! * eesq**3
    
    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)
    
    res = me2(3,4,1,2,5)*aveaa*eesq2*Q_lep2**3 ! * eesq

    return

  end subroutine res_treeAA_a_aa
    

  !----------------------------------------------------------------------
  !-- 6-point amplitudes
  !----------------------------------------------------------------------
  
  !-- qqb channel
  !-- q(p1) qb(p2) -> e-(p3) e+(p4) g(p5) a(p6)
  !-- res(1,:) -> q qb for dn and up
  !-- res(2,:) -> qb q for dn and up
  subroutine res_tree_ga_qqb(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    integer, parameter :: iconf1(6,2) = reshape([1,5,2,3,6,4, 2,5,1,3,6,4],[6,2])
    integer, parameter :: iconf2(6,2) = reshape([1,5,6,2,3,4, 2,5,6,1,3,4],[6,2])

    call res_tree_jj_ag(p,iconf1,iconf2,aveqq,res)

  end subroutine res_tree_ga_qqb

  !-- gq channel
  !-- g(p1) qb/q(p2) -> e-(p3) e+(p4) qb/q(p5) a(p6)
  !-- res(1,:) -> g qb for dn and up
  !-- res(2,:) -> g q  for dn and up
  subroutine res_tree_ga_gq(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    integer, parameter :: iconf1(6,2) = reshape([1,2,5,3,6,4, 5,2,1,3,6,4],[6,2])
    integer, parameter :: iconf2(6,2) = reshape([1,2,6,5,3,4, 5,2,6,1,3,4],[6,2])

    call res_tree_jj_ag(p,iconf1,iconf2,aveqg,res)

  end subroutine res_tree_ga_gq

  !-- qg channel
  !-- q(p1) g(p2) -> e-(p3) e+(p4) q(p5) a(p6)
  !-- res(1,:) -> q  g for dn and up
  !-- res(2,:) -> qb g  for dn and up
  subroutine res_tree_ga_qg(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    integer, parameter :: iconf1(6,2) = reshape([5,1,2,3,6,4, 2,1,5,3,6,4],[6,2])
    integer, parameter :: iconf2(6,2) = reshape([5,1,6,2,3,4, 2,1,6,5,3,4],[6,2])

    call res_tree_jj_ag(p,iconf1,iconf2,aveqg,res)

  end subroutine res_tree_ga_qg

  !-- aq channel
  !-- a(p1) qb/q(p2) -> e-(p3) e+(p4) g(p5) qb/q(p6)
  !-- res(1,:) -> a qb for dn and up
  !-- res(2,:) -> a q  for dn and up
  subroutine res_tree_ga_aq(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    integer, parameter :: iconf1(6,2) = reshape([1,5,6,3,2,4, 6,5,1,3,2,4],[6,2])
    integer, parameter :: iconf2(6,2) = reshape([1,5,2,6,3,4, 6,5,2,1,3,4],[6,2])

    call res_tree_jj_ag(p,iconf1,iconf2,aveqa,res)

  end subroutine res_tree_ga_aq

  !-- qa channel
  !-- q(p1) a(p2) -> e-(p3) e+(p4) g(p5) q(p6)
  !-- res(1,:) -> q  a for dn and up
  !-- res(2,:) -> qb a for dn and up
  subroutine res_tree_ga_qa(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    integer, parameter :: iconf1(6,2) = reshape([6,5,2,3,1,4, 2,5,6,3,1,4],[6,2])
    integer, parameter :: iconf2(6,2) = reshape([6,5,1,2,3,4, 2,5,1,6,3,4],[6,2])

    call res_tree_jj_ag(p,iconf1,iconf2,aveqa,res)

  end subroutine res_tree_ga_qa

  !-- ag channel
  !-- a(p1) g(p2) -> e-(p3) e+(p4) q(p5) qb(p6)
  !-- res(1,:) -> q qb for dn and up
  !-- res(2,:) -> qb q for dn and up
  subroutine res_tree_ga_ag(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res
    real(dp) :: res_tmp(1,2)
    integer, parameter :: iconf1(6,1) = reshape([5,1,6,3,2,4],[6,1])
    integer, parameter :: iconf2(6,1) = reshape([5,1,2,6,3,4],[6,1])

    call res_tree_jj_ag(p,iconf1,iconf2,avega,res_tmp)

    res = res_tmp(1,1) * ndn + res_tmp(1,2) * nup

  end subroutine res_tree_ga_ag

  !-- ga channel
  !-- g(p1) a(p2) -> e-(p3) e+(p4) q(p5) qb(p6)
  !-- res(1,:) -> q qb for dn and up
  !-- res(2,:) -> qb q for dn and up
  subroutine res_tree_ga_ga(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res
    real(dp) :: res_tmp(1,2)
    integer, parameter :: iconf1(6,1) = reshape([5,2,6,3,1,4],[6,1])
    integer, parameter :: iconf2(6,1) = reshape([5,2,1,6,3,4],[6,1])

    call res_tree_jj_ag(p,iconf1,iconf2,avega,res_tmp)

    res = res_tmp(1,1) * ndn + res_tmp(1,2) * nup

  end subroutine res_tree_ga_ga

  !-- master amplitude for gluon/photon emission, temporary
  !-- 3 and 4 must always be the leptons, otherwise couplings won't work
  subroutine res_tree_jj_ag(p,iconf1,iconf2,ave,res)
    real(dp), intent(in)  :: p(:,:), ave
    integer, intent(in)   :: iconf1(:,:), iconf2(:,:)
    real(dp), intent(out) :: res(size(iconf1,2),2)
    integer     :: i1,i2,i3,i4
    integer     :: j
    complex(dp) :: za(6,6), zb(6,6)
    real(dp)    :: sprod(6,6), sijk
    complex(dp) :: ampl_res1(-1:1,-1:1,-1:1,-1:1),ampl_res2(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: ampl_res_tot(-1:1,-1:1,-1:1,-1:1,1:2,1:2)
    complex(dp) :: ampl_res(-1:1,-1:1,-1:1,-1:1)
    complex(dp) :: coupl(1:2,-1:1,-1:1),coupl_fact(1:2,-1:1,-1:1)
    real(dp) :: charge(2)
    
    !contribution of the double ISR        
    call spinoru(6,(/-p(:,2),-p(:,1),p(:,3),p(:,4),p(:,5),p(:,6)/),za,zb,sprod)
    call get_coupl(sprod(3,4),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
         [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],coupl)

    !contribution of the factorised diagrams
    sijk = sprod(iconf1(4,1),iconf1(5,1)) &
         + sprod(iconf1(4,1),iconf1(6,1)) &
         + sprod(iconf1(5,1),iconf1(6,1))        
    call get_coupl(sijk,[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
         [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],coupl_fact)

    charge = [Qdn,Qup]

    res = zero
    
    do j = 1,size(iconf1,2)

       !contribution of the factorised diagrams
       !-- notation for the factorised: 0 -> j1(q,h1) j2(g,h2) j3(qb,-h1) + j4(q',h3) j5(g,h4) j6(qb',-h3)
       
       ampl_res = helamp_tree_jj_1(iconf1(1,j),iconf1(2,j),iconf1(3,j),&
            iconf1(4,j),iconf1(5,j),iconf1(6,j),&
            za,zb)
       ampl_res = ampl_res/sijk
       
       !contribution of the double ISR
       !-- notation for the ISR: 0 -> j1(q,h1) j2(g,h2) j3(g,h3) j4(qb,-h1) + j5(q',h4) j6(qb',-h4)

       ampl_res1 = helamp_tree_jj_2(iconf2(1,j),iconf2(2,j),iconf2(3,j),&
            iconf2(4,j),iconf2(5,j),iconf2(6,j),&
            za,zb,sprod)
       ampl_res1 = ampl_res1/sprod(3,4)
       
       ampl_res2 = helamp_tree_jj_2(iconf2(1,j),iconf2(3,j),iconf2(2,j),&
            iconf2(4,j),iconf2(5,j),iconf2(6,j),&
            za,zb,sprod)
       ampl_res2 = ampl_res2/sprod(3,4)

       do i4=-1,1,2
          do i3=-1,1,2
             do i2=-1,1,2
                do i1=-1,1,2
                   ampl_res_tot(i1,i2,i3,i4,j,:) = &
                        (ampl_res1(i1,i2,i3,i4) + ampl_res2(i1,i3,i2,i4))*coupl(:,i1,i4)*charge(:) &
                        + ampl_res(i1,i2,i4,i3)*coupl_fact(:,i1,i4)*Q_lep
                   res(j,:) = res(j,:) + real(ampl_res_tot(i1,i2,i3,i4,j,:)*conjg(ampl_res_tot(i1,i2,i3,i4,j,:)),kind=dp)
                enddo
             enddo
          enddo
       enddo
    enddo

    res = eesq2 * ave * four * xn * Cf * res

  end subroutine res_tree_jj_ag

    !-- notation: 0 -> j1(q,h1) j2(g,h2) j3(qb,-h1) + j4(q',h3) j5(g,h4) j6(qb',-h3)
  !-- factorised diagrams
  function helamp_tree_jj_1(j1, j2, j3, j4, j5, j6, za, zb)
    complex(dp) :: helamp_tree_jj_1(-1:1,-1:1,-1:1,-1:1)
    integer, intent(in)     :: j1, j2, j3, j4, j5, j6
    complex(dp), intent(in) :: za(6,6), zb(6,6)
    integer     :: i1, i2
    complex(dp) :: rs1_13, rs2_13, rs3_13, rs4_13
    complex(dp) :: rs1_31, rs2_31, rs3_31, rs4_31
    complex(dp) :: num1_13, num2_13, num3_13, num4_13
    complex(dp) :: num1_31, num2_31, num3_31, num4_31
    complex(dp) :: rs1_13_c, rs2_13_c, rs3_13_c, rs4_13_c
    complex(dp) :: rs1_31_c, rs2_31_c, rs3_31_c, rs4_31_c
    complex(dp) :: num1_13_c, num2_13_c, num3_13_c, num4_13_c
    complex(dp) :: num1_31_c, num2_31_c, num3_31_c, num4_31_c
    
    ! -- recurring structures
    rs1_13(i1,i2) = TWO * za(j4,j1) * zb(i1,i2) * za(i2,j4) * za(i1,j1)
    rs2_13(i1,i2) = TWO * zb(j3,j6) * za(i2,i1) * zb(i1,j3) * zb(i2,j6)
    rs3_13(i1,i2) = TWO * za(i1,j1) * zb(j6,i2) * za(i2,j1) * zb(i1,j6)
    rs4_13(i1,i2) = TWO * za(i2,j4) * zb(j3,i1) * za(i1,j4) * zb(i2,j3)
    
    rs1_31(i1,i2) = TWO * za(j4,j3) * zb(i1,i2) * za(i2,j4) * za(i1,j3)
    rs2_31(i1,i2) = TWO * zb(j1,j6) * za(i2,i1) * zb(i1,j1) * zb(i2,j6)
    rs3_31(i1,i2) = TWO * za(i1,j3) * zb(j6,i2) * za(i2,j3) * zb(i1,j6)
    rs4_31(i1,i2) = TWO * za(i2,j4) * zb(j1,i1) * za(i1,j4) * zb(i2,j1)
    
    rs1_13_c(i1,i2) = TWO * zb(j4,j1) * za(i1,i2) * zb(i2,j4) * zb(i1,j1)
    rs2_13_c(i1,i2) = TWO * za(j3,j6) * zb(i2,i1) * za(i1,j3) * za(i2,j6)
    rs3_13_c(i1,i2) = TWO * zb(i1,j1) * za(j6,i2) * zb(i2,j1) * za(i1,j6)
    rs4_13_c(i1,i2) = TWO * zb(i2,j4) * za(j3,i1) * zb(i1,j4) * za(i2,j3)
    
    rs1_31_c(i1,i2) = TWO * zb(j4,j3) * za(i1,i2) * zb(i2,j4) * zb(i1,j3)
    rs2_31_c(i1,i2) = TWO * za(j1,j6) * zb(i2,i1) * za(i1,j1) * za(i2,j6)
    rs3_31_c(i1,i2) = TWO * zb(i1,j3) * za(j6,i2) * zb(i2,j3) * za(i1,j6)
    rs4_31_c(i1,i2) = TWO * zb(i2,j4) * za(j1,i1) * zb(i1,j4) * za(i2,j1)
    
    num1_13 = rs1_13(j2,j5) + rs1_13(j2,j6) + rs1_13(j3,j5) + rs1_13(j3,j6)
    num2_13 = rs2_13(j1,j4) + rs2_13(j1,j5) + rs2_13(j2,j4) + rs2_13(j2,j5)
    num3_13 = rs3_13(j2,j4) + rs3_13(j2,j5) + rs3_13(j3,j4) + rs3_13(j3,j5)
    num4_13 = rs4_13(j1,j5) + rs4_13(j1,j6) + rs4_13(j2,j5) + rs4_13(j2,j6)
    
    num1_31 = rs1_31(j2,j5) + rs1_31(j2,j6) + rs1_31(j1,j5) + rs1_31(j1,j6)
    num2_31 = rs2_31(j3,j4) + rs2_31(j3,j5) + rs2_31(j2,j4) + rs2_31(j2,j5)
    num3_31 = rs3_31(j2,j4) + rs3_31(j2,j5) + rs3_31(j1,j4) + rs3_31(j1,j5)
    num4_31 = rs4_31(j3,j5) + rs4_31(j3,j6) + rs4_31(j2,j5) + rs4_31(j2,j6)
    
    !recurring structure for complex conjugate      
    
    num1_13_c = rs1_13_c(j2,j5) + rs1_13_c(j2,j6) + rs1_13_c(j3,j5) + rs1_13_c(j3,j6)
    num2_13_c = rs2_13_c(j1,j4) + rs2_13_c(j1,j5) + rs2_13_c(j2,j4) + rs2_13_c(j2,j5)
    num3_13_c = rs3_13_c(j2,j4) + rs3_13_c(j2,j5) + rs3_13_c(j3,j4) + rs3_13_c(j3,j5)
    num4_13_c = rs4_13_c(j1,j5) + rs4_13_c(j1,j6) + rs4_13_c(j2,j5) + rs4_13_c(j2,j6)
    
    num1_31_c = rs1_31_c(j2,j5) + rs1_31_c(j2,j6) + rs1_31_c(j1,j5) + rs1_31_c(j1,j6)
    num2_31_c = rs2_31_c(j3,j4) + rs2_31_c(j3,j5) + rs2_31_c(j2,j4) + rs2_31_c(j2,j5)
    num3_31_c = rs3_31_c(j2,j4) + rs3_31_c(j2,j5) + rs3_31_c(j1,j4) + rs3_31_c(j1,j5)
    num4_31_c = rs4_31_c(j3,j5) + rs4_31_c(j3,j6) + rs4_31_c(j2,j5) + rs4_31_c(j2,j6)
    
    ! Eq. (4.9) and Eq. (4.10)
    helamp_tree_jj_1(-1,+1,-1,+1) =   ONE / za(j1,j2) / za(j2,j3) / za(j4,j5) / za(j5,j6) * num1_13
    helamp_tree_jj_1(-1,-1,-1,-1) =   ONE / zb(j1,j2) / zb(j2,j3) / zb(j4,j5) / zb(j5,j6) * num2_13
    helamp_tree_jj_1(-1,+1,-1,-1) = - ONE / za(j1,j2) / za(j2,j3) / zb(j4,j5) / zb(j5,j6) * num3_13
    helamp_tree_jj_1(-1,-1,-1,+1) = - ONE / zb(j1,j2) / zb(j2,j3) / za(j4,j5) / za(j5,j6) * num4_13
    
    !conjugated amplitudes
    
    helamp_tree_jj_1(+1,-1,+1,-1) =   ONE / zb(j1,j2) / zb(j2,j3) / zb(j4,j5) / zb(j5,j6) * num1_13_c
    helamp_tree_jj_1(+1,+1,+1,+1) =   ONE / za(j1,j2) / za(j2,j3) / za(j4,j5) / za(j5,j6) * num2_13_c
    helamp_tree_jj_1(+1,-1,+1,+1) = - ONE / zb(j1,j2) / zb(j2,j3) / za(j4,j5) / za(j5,j6) * num3_13_c
    helamp_tree_jj_1(+1,+1,+1,-1) = - ONE / za(j1,j2) / za(j2,j3) / zb(j4,j5) / zb(j5,j6) * num4_13_c 
    !helamp_tree_jj_1(+1,-1,+1,-1) = conjg(helamp_tree_jj_1(-1,+1,-1,+1))
    !helamp_tree_jj_1(+1,+1,+1,+1) = conjg(helamp_tree_jj_1(-1,-1,-1,-1))
    !helamp_tree_jj_1(+1,-1,+1,+1) = conjg(helamp_tree_jj_1(-1,+1,-1,-1))
    !helamp_tree_jj_1(+1,+1,+1,-1) = conjg(helamp_tree_jj_1(-1,-1,-1,+1))
    
    helamp_tree_jj_1(+1,+1,-1,+1) = - ONE / za(j3,j2) / za(j2,j1) / za(j4,j5) / za(j5,j6) * num1_31
    helamp_tree_jj_1(+1,-1,-1,-1) = - ONE / zb(j3,j2) / zb(j2,j1) / zb(j4,j5) / zb(j5,j6) * num2_31
    helamp_tree_jj_1(+1,+1,-1,-1) =   ONE / za(j3,j2) / za(j2,j1) / zb(j4,j5) / zb(j5,j6) * num3_31
    helamp_tree_jj_1(+1,-1,-1,+1) =   ONE / zb(j3,j2) / zb(j2,j1) / za(j4,j5) / za(j5,j6) * num4_31
    
    !conjugated amplitudes
    
    helamp_tree_jj_1(-1,-1,+1,-1) = - ONE / zb(j3,j2) / zb(j2,j1) / zb(j4,j5) / zb(j5,j6) * num1_31_c
    helamp_tree_jj_1(-1,+1,+1,+1) = - ONE / za(j3,j2) / za(j2,j1) / za(j4,j5) / za(j5,j6) * num2_31_c
    helamp_tree_jj_1(-1,-1,+1,+1) =   ONE / zb(j3,j2) / zb(j2,j1) / za(j4,j5) / za(j5,j6) * num3_31_c
    helamp_tree_jj_1(-1,+1,+1,-1) =   ONE / za(j3,j2) / za(j2,j1) / zb(j4,j5) / zb(j5,j6) * num4_31_c
    !helamp_tree_jj_1(-1,-1,+1,-1) = conjg(helamp_tree_jj_1(+1,+1,-1,+1))
    !helamp_tree_jj_1(-1,+1,+1,+1) = conjg(helamp_tree_jj_1(+1,-1,-1,-1))
    !helamp_tree_jj_1(-1,-1,+1,+1) = conjg(helamp_tree_jj_1(+1,+1,-1,-1))
    !helamp_tree_jj_1(-1,+1,+1,-1) = conjg(helamp_tree_jj_1(+1,-1,-1,+1))
    
    return

  end function helamp_tree_jj_1

  !-- notation: 0 -> j1(q,h1) j2(g,h2) j3(g,h3) j4(qb,-h1) + j5(q',h4) j6(qb',-h4)
  !-- double ISR
  function helamp_tree_jj_2(j1, j2, j3, j4, j5, j6, za, zb, sprod)
    complex(dp) :: helamp_tree_jj_2(-1:1,-1:1,-1:1,-1:1)
    integer, intent(in)     :: j1, j2, j3, j4, j5, j6
    complex(dp), intent(in) :: za(6,6), zb(6,6)
    real(dp), intent(in)    :: sprod(6,6)
    integer     :: i1,i2,i3,i4
    real(dp)    :: s
    complex(dp) :: rf1, rf2, rf3, rf4, rf1rev
    complex(dp) :: rs1
    complex(dp) :: rs1_c, rf1_c, rf2_c, rf3_c, rf4_c
    complex(dp) :: rf1rev_c
    
    !-- recurring structures
    s(i1,i2,i3) = sprod(i1,i2) + sprod(i1,i3) + sprod(i2,i3) ! s123 = s12 + s13 + s23
    rs1(i1,i2,i3,i4) = za(i1,i2) * zb(i2,i4) + za(i1,i3) * zb(i3,i4) ! <i1|i2+i3|i4]
    
    rf1(i1,i2) = - TWO * za(j1,i1) * (zb(i2,j2) * za(j2,j1) + zb(i2,j3) * za(j3,j1) + zb(i2,j4) * za(j4,j1))
    rf2(i1,i2) = - TWO * zb(i2,j4) * (zb(j4,j1) * za(j1,i1) + zb(j4,j2) * za(j2,i1) + zb(j4,j3) * za(j3,i1))
    rf3(i1,i2) = - TWO * (za(i1,j1) * zb(j1,j3) + za(i1,j2) * zb(j2,j3)) &
         * (rs1(j4,j1,j2,j3) * zb(j4,i2) + zb(j3,i2) * s(j1,j2,j3))
    rf4(i1,i2) = - TWO * (za(j2,j3) * zb(j3,i2) + za(j2,j4) * zb(j4,i2)) * (za(i1,j1) &
         * rs1(j2,j3,j4,j1) + za(i1,j2) * s(j2,j3,j4))
    
    rf1rev(i1,i2) = - TWO * za(j4,i1) * (zb(i2,j3) * za(j3,j4) + zb(i2,j2) * za(j2,j4) + zb(i2,j1) * za(j1,j4))
    
    ! recurring structure for complex conjugate
    rs1_c(i1,i2,i3,i4) = zb(i1,i2) * za(i2,i4) + zb(i1,i3) * za(i3,i4) 
    rf1_c(i1,i2) = - TWO * zb(j1,i1) * (za(i2,j2) * zb(j2,j1) + za(i2,j3) * zb(j3,j1) + za(i2,j4) * zb(j4,j1))
    rf2_c(i1,i2) = - TWO * za(i2,j4) * (za(j4,j1) * zb(j1,i1) + za(j4,j2) * zb(j2,i1) + za(j4,j3) * zb(j3,i1))
    rf3_c(i1,i2) = - TWO * (zb(i1,j1) * za(j1,j3) + zb(i1,j2) * za(j2,j3)) &
         * (rs1_c(j4,j1,j2,j3) * za(j4,i2) + za(j3,i2) * s(j1,j2,j3))
    rf4_c(i1,i2) = - TWO * (zb(j2,j3) * za(j3,i2) + zb(j2,j4) * za(j4,i2)) * (zb(i1,j1) &
         * rs1_c(j2,j3,j4,j1) + zb(i1,j2) * s(j2,j3,j4))
    rf1rev_c(i1,i2) = - TWO * zb(j4,i1) * (za(i2,j3) * zb(j3,j4) + za(i2,j2) * zb(j2,j4) + za(i2,j1) * zb(j1,j4))
    
    !-- Eq. (4.11) and Eq. (4.5)
    helamp_tree_jj_2(-1,+1,+1,-1) = ONE / za(j1,j2) / za(j2,j3) /za(j3,j4) * rf1(j5,j6)
    helamp_tree_jj_2(-1,+1,+1,+1) = ONE / za(j1,j2) / za(j2,j3) /za(j3,j4) * rf1(j6,j5)
    
    
    helamp_tree_jj_2(+1,-1,-1,+1) = ONE / zb(j1,j2) / zb(j2,j3) /zb(j3,j4) * rf1_c(j5,j6)
    helamp_tree_jj_2(+1,-1,-1,-1) = ONE / zb(j1,j2) / zb(j2,j3) /zb(j3,j4) * rf1_c(j6,j5)
    !helamp_tree_jj_2(+1,-1,-1,+1) = conjg(helamp_tree_jj_2(-1,+1,+1,-1))
    !helamp_tree_jj_2(+1,-1,-1,-1) = conjg(helamp_tree_jj_2(-1,+1,+1,+1))
    
    !-- Eq. (4.11) (using "line reversal" symmetry) and Eq. (4.5)
    !-- Note: in constrast to Eq. (4.5) there is no additional minus sign (but I think this has no effect anyway)
    helamp_tree_jj_2(+1,+1,+1,-1) = ONE / za(j4,j3) / za(j3,j2) / za(j2,j1) * rf1rev(j5,j6) 
    helamp_tree_jj_2(+1,+1,+1,+1) = ONE / za(j4,j3) / za(j3,j2) / za(j2,j1) * rf1rev(j6,j5) 
    
    helamp_tree_jj_2(-1,-1,-1,+1) = ONE / zb(j4,j3) / zb(j3,j2) / zb(j2,j1) * rf1rev_c(j5,j6)
    helamp_tree_jj_2(-1,-1,-1,-1) = ONE / zb(j4,j3) / zb(j3,j2) / zb(j2,j1) * rf1rev_c(j6,j5)
    !helamp_tree_jj_2(-1,-1,-1,+1) = conjg(helamp_tree_jj_2(+1,+1,+1,-1))
    !helamp_tree_jj_2(-1,-1,-1,-1) = conjg(helamp_tree_jj_2(+1,+1,+1,+1))
    
    !-- Eq. (4.12) and Eq. (4.5)
    helamp_tree_jj_2(-1,-1,+1,-1) = zb(j1,j3) / s(j1,j2,j3) / zb(j1,j2) / zb(j2,j3) / rs1(j4,j2,j3,j1) * rf3(j5,j6) &
         + za(j2,j4) / s(j2,j3,j4) / za(j2,j3) / za(j3,j4) / rs1(j4,j2,j3,j1) * rf4(j5,j6)
    helamp_tree_jj_2(-1,-1,+1,+1) = zb(j1,j3) / s(j1,j2,j3) /zb(j1,j2) / zb(j2,j3) / rs1(j4,j2,j3,j1) * rf3(j6,j5) &
         + za(j2,j4) / s(j2,j3,j4) / za(j2,j3) / za(j3,j4) / rs1(j4,j2,j3,j1) * rf4(j6,j5)
    
    helamp_tree_jj_2(+1,+1,-1,+1) = za(j1,j3) / s(j1,j2,j3) / za(j1,j2) / za(j2,j3) / rs1_c(j4,j2,j3,j1) * rf3_c(j5,j6) &
         + zb(j2,j4) / s(j2,j3,j4) / zb(j2,j3) / zb(j3,j4) / rs1_c(j4,j2,j3,j1) * rf4_c(j5,j6)
    helamp_tree_jj_2(+1,+1,-1,-1) = za(j1,j3) / s(j1,j2,j3) /za(j1,j2) / za(j2,j3) / rs1_c(j4,j2,j3,j1) * rf3_c(j6,j5) &
         + zb(j2,j4) / s(j2,j3,j4) / zb(j2,j3) / zb(j3,j4) / rs1_c(j4,j2,j3,j1) * rf4_c(j6,j5)
    !helamp_tree_jj_2(+1,+1,-1,+1) = conjg(helamp_tree_jj_2(-1,-1,+1,-1))
    !helamp_tree_jj_2(+1,+1,-1,-1) = conjg(helamp_tree_jj_2(-1,-1,+1,+1))
    
    !-- Eq. (4.13) and Eq. (4.5)
    helamp_tree_jj_2(-1,+1,-1,-1) = zb(j2,j4)**3 / s(j2,j3,j4) / zb(j2,j3) / zb(j3,j4) / rs1(j1,j2,j3,j4) * rf1(j5,j6) &
         + za(j1,j3)**3 / s(j1,j2,j3) / za(j1,j2) / za(j2,j3) / rs1(j1,j2,j3,j4) * rf2(j5,j6)
    helamp_tree_jj_2(-1,+1,-1,+1) = zb(j2,j4)**3 / s(j2,j3,j4) / zb(j2,j3) / zb(j3,j4) / rs1(j1,j2,j3,j4) * rf1(j6,j5) &
         + za(j1,j3)**3 / s(j1,j2,j3) / za(j1,j2) / za(j2,j3) / rs1(j1,j2,j3,j4) * rf2(j6,j5)
    
    helamp_tree_jj_2(+1,-1,+1,+1) = za(j2,j4)**3 / s(j2,j3,j4) / za(j2,j3) / za(j3,j4) / rs1_c(j1,j2,j3,j4) * rf1_c(j5,j6) &
         + zb(j1,j3)**3 / s(j1,j2,j3) / zb(j1,j2) / zb(j2,j3) / rs1_c(j1,j2,j3,j4) * rf2_c(j5,j6)
    helamp_tree_jj_2(+1,-1,+1,-1) = za(j2,j4)**3 / s(j2,j3,j4) / za(j2,j3) / za(j3,j4) / rs1_c(j1,j2,j3,j4) * rf1_c(j6,j5) &
         + zb(j1,j3)**3 / s(j1,j2,j3) / zb(j1,j2) / zb(j2,j3) / rs1_c(j1,j2,j3,j4) * rf2_c(j6,j5)
    !helamp_tree_jj_2(+1,-1,+1,+1) = conjg(helamp_tree_jj_2(-1,+1,-1,-1))
    !helamp_tree_jj_2(+1,-1,+1,-1) = conjg(helamp_tree_jj_2(-1,+1,-1,+1))
    
    return

  end function helamp_tree_jj_2

  
  !----------------------------------------------------------------------------------------------------

  
  !-- qqb channel, qqb final state
  !-- q(p1) qb(p2) -> l(p3) lb(p4) q(p5) qb(p6)
  !-- res(1,:) = q qb for dn and up
  !-- res(2,:) = qb q for dn and up
  subroutine res_tree_qqb_qqb(p,res_jj)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res_jj(2,2)
    integer, parameter :: iconf(6,2) = reshape([1,2,5,6,3,4, 2,1,5,6,3,4],[6,2])

    call restree_jj_qqll(p,iconf,res_jj)

  end subroutine res_tree_qqb_qqb

  !-- qq channel, qq final state
  !-- q(p1) q(p2) -> l(p3) lb(p4) q(p5) q(p6)
  !-- res(1,:) = q  q  for dn and up
  !-- res(2,:) = qb qb for dn and up
  !-- it includes symmetry factor
  subroutine res_tree_qq_qq(p,res_jj)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res_jj(2,2)
    integer, parameter :: iconf(6,2) = reshape([1,5,6,2,3,4, 5,1,2,6,3,4],[6,2])

    call restree_jj_qqll(p,iconf,res_jj)
    res_jj = res_jj/two !-- 1/2! symmetry factor

  end subroutine res_tree_qq_qq

  !-- master amplitude for qqb -> qqb, for generic crossing
  !-- iconf: crossing from 0 -> qb q Q QB l lb
  !-- already includes aveqq, does not include any 1/2! symmetry factors
  subroutine restree_jj_qqll(p,myis,res_jj)
    real(dp), intent(in)  :: p(:,:)
    integer, intent(in)   :: myis(:,:)
    real(dp), intent(out) :: res_jj(size(myis,2),2)
    complex(dp) :: za(6,6),zb(6,6)
    real(dp) :: sprod(6,6)
    integer :: h1,h5
    complex(dp) :: amp_a(-1:1,-1:1,-1:1)
    complex(dp) :: amp_b(-1:1,-1:1,-1:1)
    complex(dp) :: amp_c(-1:1,-1:1,-1:1)
    complex(dp) :: amp_d(-1:1,-1:1,-1:1)
    complex(dp) :: amp_g(-1:1,-1:1,-1:1)
    complex(dp) :: amp_h(-1:1,-1:1,-1:1)
    real(dp) :: pref
    real(dp), parameter :: Q4(2) = (/Qdn**4,Qup**4/)*Q_lep**2
    real(dp), parameter :: Q3(2) = (/Qdn**3,Qup**3/)*Q_lep**3
    integer :: i,i1,i2,i3,i4,i5,i6
    complex(dp) :: c12_ql(1:2,-1:1,-1:1),c34_ql(1:2,-1:1,-1:1)
    complex(dp) :: c13_ql(1:2,-1:1,-1:1),c24_ql(1:2,-1:1,-1:1)
    complex(dp) :: c56_ql(1:2,-1:1,-1:1)
    complex(dp) :: c12_qq(1:2,-1:1,-1:1),c34_qq(1:2,-1:1,-1:1)
    complex(dp) :: c13_qq(1:2,-1:1,-1:1),c24_qq(1:2,-1:1,-1:1)

    res_jj = zero

    call spinoru(6,(/-p(:,1),-p(:,2),p(:,3:6)/),za,zb,sprod)

    !-- routine in the do-loop is for 0 -> qb q Q QB l bl

    do i = 1,size(myis,2)
       i1 = myis(1,i); i2 = myis(2,i); i3 = myis(3,i)
       i4 = myis(4,i); i5 = myis(5,i); i6 = myis(6,i)
    
       call master_amp_qqbQQBV_az(i2,i1,i3,i4,i5,i6,za,zb,sprod,amp_a)
       call master_amp_qqbQQBV_az(i3,i4,i2,i1,i5,i6,za,zb,sprod,amp_b)
       call master_amp_qqbQQBV_az(i3,i1,i2,i4,i5,i6,za,zb,sprod,amp_c)
       call master_amp_qqbQQBV_az(i2,i4,i3,i1,i5,i6,za,zb,sprod,amp_d)
       !                          
       call master_amp_qqbQQBV_az(i5,i6,i2,i1,i3,i4,za,zb,sprod,amp_g) 
       call master_amp_qqbQQBV_az(i5,i6,i3,i1,i2,i4,za,zb,sprod,amp_h)

       call get_coupl(sprod(i5,i6),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
            [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],c56_ql)
       
       call get_coupl(sprod(i1,i2),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
            [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],c12_ql)
       call get_coupl(sprod(i1,i2),[Qdn,Qup],[Qdn,Qup],[cms_cLdn,cms_cLup],[cms_cLdn,cms_cLup],&
            [cms_cRdn,cms_cRup],[cms_cRdn,cms_cRup],c12_qq)
       
       call get_coupl(sprod(i1,i3),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
            [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],c13_ql)
       call get_coupl(sprod(i1,i3),[Qdn,Qup],[Qdn,Qup],[cms_cLdn,cms_cLup],[cms_cLdn,cms_cLup],&
            [cms_cRdn,cms_cRup],[cms_cRdn,cms_cRup],c13_qq)
       
       call get_coupl(sprod(i2,i4),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
            [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],c24_ql)
       call get_coupl(sprod(i2,i4),[Qdn,Qup],[Qdn,Qup],[cms_cLdn,cms_cLup],[cms_cLdn,cms_cLup],&
            [cms_cRdn,cms_cRup],[cms_cRdn,cms_cRup],c24_qq)
       
       call get_coupl(sprod(i3,i4),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLdn,cms_cLup],[cms_cL_lep,cms_cL_lep],&
            [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],c34_ql)
       call get_coupl(sprod(i3,i4),[Qdn,Qup],[Qdn,Qup],[cms_cLdn,cms_cLup],[cms_cLdn,cms_cLup],&
            [cms_cRdn,cms_cRup],[cms_cRdn,cms_cRup],c34_qq)

       do h5 = -1,1,2
          do h1 = -1,1,2
             !-- minus sign: odd fermion flip
             res_jj(i,:) = res_jj(i,:) &
                  !-- (a+b) vs (c+d); c34 + c24 etc: ga + ag 
                  - abs(c56_ql(:,h1,h5))**2 * ( & 
                  two * real(amp_a(h1,h1,h5)*conjg(amp_c(h1,h1,h5))*(c34_qq(:,h1,h1)+conjg(c24_qq(:,h1,h1))),kind=dp) + &
                  two * real(amp_a(h1,h1,h5)*conjg(amp_d(h1,h1,h5))*(c34_qq(:,h1,h1)+conjg(c13_qq(:,h1,h1))),kind=dp) + &
                  two * real(amp_b(h1,h1,h5)*conjg(amp_c(h1,h1,h5))*(c12_qq(:,h1,h1)+conjg(c24_qq(:,h1,h1))),kind=dp) + &
                  two * real(amp_b(h1,h1,h5)*conjg(amp_d(h1,h1,h5))*(c12_qq(:,h1,h1)+conjg(c13_qq(:,h1,h1))),kind=dp)   &
                  ) &
                  !-- (c+d) vs g, g is purely electroweak
                  - ( &
                  two * real(amp_c(h1,h1,h5)*conjg(amp_g(h5,h1,h1))*(c56_ql(:,h1,h5)*conjg(c34_ql(:,h1,h5)*c12_ql(:,h1,h5))),kind=dp) + &
                  two * real(amp_d(h1,h1,h5)*conjg(amp_g(h5,h1,h1))*(c56_ql(:,h1,h5)*conjg(c34_ql(:,h1,h5)*c12_ql(:,h1,h5))),kind=dp)   &
                  ) &
                  !-- (a+b) vs h, h is purely electroweak
                  - ( &
                  two * real(amp_a(h1,h1,h5)*conjg(amp_h(h5,h1,h1))*(c56_ql(:,h1,h5)*conjg(c13_ql(:,h1,h5)*c24_ql(:,h1,h5))),kind=dp) + &
                  two * real(amp_b(h1,h1,h5)*conjg(amp_h(h5,h1,h1))*(c56_ql(:,h1,h5)*conjg(c13_ql(:,h1,h5)*c24_ql(:,h1,h5))),kind=dp)   &
                  )
          enddo
       enddo

    enddo
       
    !-- four**2: to compensate for the 1/sqrt2 in the vertex
    pref = 16._dp * (xn * Cf) * eesq2 * aveqq ! * eesq * gs**2
    res_jj = res_jj * pref
    
  end subroutine restree_jj_qqll
   
  !-- amplitude for 0 -> q(1) qb(2) [g* -> Q(3) QB(4)] qb(4) [V --> l(5) lb(6)]
  !-- V attached to the 12 line
  !-- index: h1,h3,h5,down/up
  ! -- elementary vertex: 1/sqrt2 * gamma(mu)
  subroutine master_amp_qqbQQBV_az(i1,i2,i3,i4,i5,i6,za,zb,sprod,amp)
    integer, intent(in) :: i1,i2,i3,i4,i5,i6
    complex(dp), intent(in) :: za(6,6),zb(6,6)
    real(dp), intent(in) :: sprod(6,6)
    complex(dp), intent(out) :: amp(-1:1,-1:1,-1:1)
    complex(dp) :: zab2
    real(dp) :: dena,denb
    
    zab2(i1,i2,i3,i4) = za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
    
    dena = (sprod(i1,i3)+sprod(i1,i4)+sprod(i3,i4))*sprod(i3,i4)*sprod(i5,i6) !-- s134*s34*s56
    denb = (sprod(i2,i3)+sprod(i2,i4)+sprod(i3,i4))*sprod(i3,i4)*sprod(i5,i6) !-- s234*s34*s56
    
    amp(-1,-1,-1) =  (za(i1,i3)*zb(i6,i2)*zab2(i5,i1,i3,i4)/dena-za(i1,i5)*zb(i4,i2)*zab2(i3,i2,i4,i6)/denb)
    amp(-1,-1,+1) =  (za(i1,i3)*zb(i5,i2)*zab2(i6,i1,i3,i4)/dena-za(i1,i6)*zb(i4,i2)*zab2(i3,i2,i4,i5)/denb)
    amp(-1,+1,-1) =  (za(i1,i4)*zb(i6,i2)*zab2(i5,i1,i4,i3)/dena-za(i1,i5)*zb(i3,i2)*zab2(i4,i2,i3,i6)/denb)
    amp(-1,+1,+1) =  (za(i1,i4)*zb(i5,i2)*zab2(i6,i1,i4,i3)/dena-za(i1,i6)*zb(i3,i2)*zab2(i4,i2,i3,i5)/denb)
    !
    ! amp(+1,-1,-1) = (-za(i2,i3)*zb(i6,i1)*zab2(i5,i2,i3,i4)/denb+za(i2,i5)*zb(i4,i1)*zab2(i3,i1,i4,i6)/dena)
    ! amp(+1,-1,+1) = (-za(i2,i3)*zb(i5,i1)*zab2(i6,i2,i3,i4)/denb+za(i2,i6)*zb(i4,i1)*zab2(i3,i1,i4,i5)/dena)
    ! amp(+1,+1,-1) = (-za(i2,i4)*zb(i6,i1)*zab2(i5,i2,i4,i3)/denb+za(i2,i5)*zb(i3,i1)*zab2(i4,i1,i3,i6)/dena)
    ! amp(+1,+1,+1) = (-za(i2,i4)*zb(i5,i1)*zab2(i6,i2,i4,i3)/denb+za(i2,i6)*zb(i3,i1)*zab2(i4,i1,i3,i5)/dena)
    !
    amp(+1,-1,-1) = conjg(amp(-1,+1,+1))
    amp(+1,-1,+1) = conjg(amp(-1,+1,-1))
    amp(+1,+1,-1) = conjg(amp(-1,-1,+1))
    amp(+1,+1,+1) = conjg(amp(-1,-1,-1))
    
    return

  end subroutine master_amp_qqbQQBV_az

end module mod_amplitudes_tree_ppll
     

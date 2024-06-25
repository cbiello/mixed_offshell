!-- amplitudes with extra jet: factor out gs**2 or ee**2 per extra gluon/photon emission
module mod_amplitudes_tree_ppnul
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_auxfunctions
  use mod_coupl
  use openloops
  use mod_ol_interface
  implicit none
  private

  !-- res_tree_emitted-particles_initial-state
  !-- res_tree   --> 0 -> q l nu qb
  public :: res_tree_qqb_w 
  public :: res_tree_g_qqb_w
  public :: res_tree_a_qqb_w
  
  !master amplitudes
  public :: res_tree_j_qcd
contains

  !----------------------------------------------------------------------
  !-- 4-point amplitudes
  !----------------------------------------------------------------------

  !-- f(1) + f(2) --> [W --> eb(3) + nu(4)]
  !-- summed / averaged

  !-- res(1,1)+ = u db -> W+ -> ep nue
  !-- res(1,2)+ = db u -> W+ -> ep nue
  !-- res(1,1)- = d ub -> W- -> em nueb
  !-- res(1,2)- = ub d -> W- -> em nueb

  subroutine res_tree_qqb_w(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    real(dp) :: sprod(4,4)
    complex(dp) :: za(4,4),zb(4,4),coupl(1:2,-1:1,-1:1)
    complex(dp) :: amp
    real(dp) :: aa,bb,cc(1:2),dd(1:2)
    integer  :: i1,i2,i3,i4

    ! u(1) db(2) -> W+ nue(4) ep(3)
    !res(1,1) = abs(amp(2,1,3,4)*coupl(1,-1,-1))**2

    ! d(1) ub(2) -> W- nueb(4) em(3)
    !res(1,2) = abs(amp(1,2,3,4)*coupl(2,-1,-1))**2

    ! ub(1) d(2) -> W- nueb(4) em(3)
    !res(2,1) = abs(amp(2,1,3,4)*coupl(2,-1,-1))**2

    ! db(1) u(2) -> W+ nue(4) ep(3)
    !res(2,2) = abs(amp(1,2,3,4)*coupl(1,-1,-1))**2

    !res = res * xn * aveqq

    !-- 0 --> q(1)^- qb(2) [V --> l(3)^- lb(4)]
    amp(i1,i2,i3,i4) = two * za(i1,i4)*zb(i3,i2)/sprod(i1,i2)
    !CB: cheanged in order to have the nu in the fourth position

    !Here we need to pass the new order of PS
    call spinoru(4,(/-p(:,1),-p(:,2),p(:,3),p(:,4)/),za,zb,sprod)

    call get_coupl(sprod(1,2),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLWud,cms_cLWud],[cms_cLWnue,cms_cLWnue],&
         [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],coupl,1)
    
    aa = abs(amp(2,1,3,4))**2
    bb = abs(amp(1,2,3,4))**2
    cc(:) = abs(coupl(:,-1,-1))**2 + abs(coupl(:,+1,+1))**2
    dd(:) = abs(coupl(:,-1,+1))**2 + abs(coupl(:,+1,-1))**2

    res(1,:) = aa * cc(:) + bb * dd(:)
    res(2,:) = aa * dd(:) + bb * cc(:)

    res = res * xn * eesq2 * aveqq

    return

  end subroutine res_tree_qqb_w


  !----------------------------------------------------------------------
  !-- 5-point amplitudes
  !----------------------------------------------------------------------

  !-- f(1) + f(2) --> [W --> nu(3) + eb(4)] + g(5)
  !-- summed / averaged
  !-- qqb channel
  !-- q(p1) qb(p2) -> nu(p3) eb(p4) g(p5)
  !-- res(1,:) = q qb for dn and up
  !-- res(2,:) = qb q for dn and up
  subroutine res_tree_g_qqb_w(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    integer, parameter :: iconf(5,2) = reshape([2,5,1,4,3, 1,5,2,4,3],[5,2])
    !CB changed iconf in order to switch nue<->e

    call res_tree_j_qcd(p,iconf,aveqq,res)

  end subroutine res_tree_g_qqb_w

  subroutine res_tree_a_qqb_w(p,res)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res(2,2)
    integer, parameter :: iconf(5,2) = reshape([2,5,1,4,3, 1,5,2,4,3],[5,2])
    !CB changed iconf in order to switch nue<->e

    call res_tree_j_qed(p,iconf,aveqq,res)

  end subroutine res_tree_a_qqb_w


  !-- master amplitudes below

  !-- generic amplitude for nlo QCD mission
  !-- iconf(5,:) are the required crossings from 0 -> q g qb l lb [e.g. q qb -> e-e+ g is 2,5,1,4,3]
  !-- ave is the averaging factor
  !-- always assume that leptons are 3 and 4, otherwise coupling is wrong
  subroutine res_tree_j_qcd(p,iconf,ave,res)
    real(dp), intent(in)  :: p(:,:),ave
    integer, intent(in)   :: iconf(:,:)
    real(dp), intent(out) :: res(size(iconf,2),2)
    integer     :: i,hq,hg,hl
    real(dp)    :: sprod(5,5)
    complex(dp) :: za(5,5),zb(5,5),coupl(1:2,-1:1,-1:1),amp(-1:1,-1:1,-1:1)

    res = zero

    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)
    !the spinors now are determined with the new order nue<->e
    call get_coupl(sprod(3,4),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLWud,cms_cLWud],[cms_cLWnue,cms_cLWnue],&
         [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],coupl,1)
 
 
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
  !-- iconf(5,:) are the required crossings from 0 -> q a qb l lb [e.g. qqb -> e-e+ a is 2,5,1,4,3]
  !-- ave is the averaging factor
  subroutine res_tree_j_qed(p,iconf,ave,res)
    real(dp), intent(in)  :: p(:,:),ave
    real(dp15) :: p_ol(4,5)
    integer, intent(in)   :: iconf(:,:)
    real(dp), intent(out) :: res(size(iconf,2),2)
    integer     :: i,hq,ha,hl
    real(dp)    :: sprod(5,5)
    complex(dp) :: za(5,5),zb(5,5),coupl_is(1:2,-1:1,-1:1),coupl_fs(1:2,-1:1,-1:1)
    complex(dp) :: amp_is1(-1:1,-1:1,-1:1),amp_is2(-1:1,-1:1,-1:1),amp_fs(-1:1,-1:1,-1:1),amp(2)
    integer :: ismin,ismax,fsmin,fsmax
    logical :: need
    real(dp) :: Wcorrup, Wcorrdown

    real(dp) :: mwMG, gawMG, mypropsq

    !CB MG settings
    mwMG = 80.419002445756163_dp
    gawMG = 2.0476000000000001_dp
    mypropsq = mwMG**2 - mwMG * gawMG

    res = zero
    ismin = -1; ismax = -1
    fsmin = -1; fsmax = -1

    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)

    !-- W emission as a correction factor of the initial-state emissions
    Wcorrup = Q_lep*sprod(1,5)/(sprod(1,2)-mypropsq)
    Wcorrdown = Q_lep*sprod(2,5)/(sprod(1,2)-mypropsq)

    do i = 1,size(iconf,2)

       !-- initial-state emission
       amp_is1 = master_amp_5pt_firstdiagram(iconf(1,i),iconf(2,i),iconf(3,i),iconf(4,i),iconf(5,i),za,zb)
       amp_is2 = master_amp_5pt_seconddiagram(iconf(1,i),iconf(2,i),iconf(3,i),iconf(4,i),iconf(5,i),za,zb)
       call need_coupl(iconf(4,i),iconf(5,i),ismin,ismax,need)
       !CB (11May): here we have to put 1 in order to call the coupling for charged bosons
       if(need) call get_coupl(sprod(ismin,ismax),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLWud,cms_cLWud],[cms_cLWnue,cms_cLWnue],&
         [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],coupl_is,1)
       
       !-- final_state emission
       amp_fs = master_amp_5pt_seconddiagram(iconf(4,i),iconf(2,i),iconf(5,i),iconf(1,i),iconf(3,i),za,zb)
       call need_coupl(iconf(1,i),iconf(3,i),fsmin,fsmax,need)
       if (need) call get_coupl(sprod(fsmin,fsmax),[Qdn,Qup],[Q_lep,Q_lep],[cms_cLWud,cms_cLWud],[cms_cLWnue,cms_cLWnue],&
         [cms_cRdn,cms_cRup],[cms_cR_lep,cms_cR_lep],coupl_fs,1)
       
       do hl=-1,1,2
          do ha=-1,1,2
             do hq=-1,1,2

                !amp = [Qdn-Wcorrdown,Qup-Wcorrup]*amp_is1(hq,ha,hl)*coupl_is(:,hq,hl) + [Qup-Wcorrup,Qdn-Wcorrdown]*amp_is2(hq,ha,hl)*coupl_is(:,hq,hl) + Q_lep*amp_fs(hl,ha,hq)*coupl_fs(:,hq,hl)
                amp = (Qdn+Wcorrdown)*amp_is1(hq,ha,hl)*coupl_is(:,hq,hl) + (Qup+Wcorrup)*amp_is2(hq,ha,hl)*coupl_is(:,hq,hl) - Q_lep*amp_fs(hl,ha,hq)*coupl_fs(:,hq,hl)

                if(i.eq.2) amp = (Qup+Wcorrup)*amp_is1(hq,ha,hl)*coupl_is(:,hq,hl) + (Qdn+Wcorrdown)*amp_is2(hq,ha,hl)*coupl_is(:,hq,hl) - Q_lep*amp_fs(hl,ha,hq)*coupl_fs(:,hq,hl)

                res(i,:) = res(i,:) + abs(amp)**2
 
             enddo
          enddo
       enddo
       
    enddo

    res = res * 8._dp * xn * ave * eesq2 ! * ee**2

    !res(:) gives both w+ and w- amplitudes, 
    !this is set in mod_ol_interface.f90
    ! --> for both w+ and w- we have only res(1,1) and res(1,2) 
    
    !CB: data la disperazione imposto OL
    p_ol(:,1:5)=p(:,1:5)
    call evaluate_tree(OL_id(1), p_ol, res(1,1))
    call evaluate_tree(OL_id(2), p_ol, res(1,2))


    print*, 'OL_id(1)= ', OL_id(1)


    res(1,1)=res(1,1)/eesq
    res(1,2)=res(1,2)/eesq
    !For now we use the OL amplitudes, in this way Chiara can use them for the checks of the subtraction.
    !Later (soon) I will fix the analytic amplitudes since it is a shame...
    
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
    res(-1,+1,-1) = za(i1,i4)**2*dena  !-- emission from i1
    res(-1,-1,-1) = zb(i5,i3)**2*denb  !-- emission from i3

    !-- 4 <-> 5
    res(-1,+1,+1) = -zero*za(i1,i5)**2*dena
    res(-1,-1,+1) = -zero*zb(i4,i3)**2*denb

    !-- other helicities
    !-- careful, complex conjugate does not
    !-- work with crossing
    res(+1,-1,+1) = zero*zb(i1,i4)**2*denb
    res(+1,+1,+1) = zero*za(i5,i3)**2*dena

    res(+1,-1,-1) = -zero*zb(i1,i5)**2*denb
    res(+1,+1,-1) = -zero*za(i4,i3)**2*dena

  end function master_amp_qgqb_llb

  !-- CB (25apr24)
  !-- master amplitude for 0 -> q[p1] g[p2] qb[p3] l[p4] lb[p5]
  !-- but the gluon g[p2] is radiated only from q[p3] and not form qb[p1]
  !-- these are useful for NLO EW W for final state radiations
  !-- helicity label: h1,h2,h4
  !-- amp[-h1,-h2,-h3] = -(amp[h1,h2,h3])^*
  !-- vertices are all gamma/sqrt2
  !-- propagator factor for the photon: 1/s45
  function master_amp_qgqb_llb_afin(i1,i2,i3,i4,i5,za,zb) result(res)
    complex(dp) :: res(-1:1,-1:1,-1:1)
    integer, intent(in)     :: i1,i2,i3,i4,i5
    complex(dp), intent(in) :: za(5,5),zb(5,5)
    complex(dp) :: dena,denb

    dena = one/(za(i1,i2)*za(i2,i3)*za(i4,i5))
    denb = one/(zb(i1,i2)*zb(i2,i3)*zb(i5,i4))

    !-- basic amplitudes
    !-- I inserted a zero just to switch off the photon radiation from the neutrino
    res(-1,+1,-1) = za(i1,i4)**2*dena
    res(-1,-1,-1) = zero*zb(i5,i3)**2*denb

    !-- 4 <-> 5
    res(-1,+1,+1) = -zero*za(i1,i5)**2*dena
    res(-1,-1,+1) = -zero*zb(i4,i3)**2*denb

    !-- other helicities
    !-- careful, complex conjugate does not
    !-- work with crossing
    res(+1,-1,+1) = zero*zb(i1,i4)**2*denb
    res(+1,+1,+1) = zero*za(i5,i3)**2*dena

    res(+1,-1,-1) = -zero*zb(i1,i5)**2*denb
    res(+1,+1,-1) = -zero*za(i4,i3)**2*dena

  end function master_amp_qgqb_llb_afin

  !-- CB (19may24)
  !-- I realised that, in the derivation of the master amp,
  !-- there was the hypotesis of same coupling for both initial radiations.
  !-- Since we don't have ux-u or dx-d but we have fermion lines with different
  !-- charges (u-dx and dx-u) I have to split the two initial state emisisons.
  !-- Let me introduce firstdiagram and seconddigram for this purpose.
  function master_amp_5pt_firstdiagram(i1,i2,i3,i4,i5,za,zb) result(res)
    complex(dp) :: res(-1:1,-1:1,-1:1)
    integer, intent(in)     :: i1,i2,i3,i4,i5
    complex(dp), intent(in) :: za(5,5),zb(5,5)
    complex(dp) :: dena,denb

    dena = one/(za(i1,i2)*za(i2,i3)*za(i4,i5))
    denb = one/(zb(i1,i2)*zb(i2,i3)*zb(i5,i4))

    !-- basic amplitudes
    !-- I inserted a zero just to switch off the photon radiation from the neutrino
    res(-1,+1,-1) = za(i1,i4)**2*dena
    res(-1,-1,-1) = zero*zb(i5,i3)**2*denb

    !-- 4 <-> 5
    res(-1,+1,+1) = -zero
    res(-1,-1,+1) = -zero
    res(+1,-1,+1) = zero
    res(+1,+1,+1) = zero
    res(+1,-1,-1) = -zero
    res(+1,+1,-1) = -zero

  end function master_amp_5pt_firstdiagram

  function master_amp_5pt_seconddiagram(i1,i2,i3,i4,i5,za,zb) result(res)
    complex(dp) :: res(-1:1,-1:1,-1:1)
    integer, intent(in)     :: i1,i2,i3,i4,i5
    complex(dp), intent(in) :: za(5,5),zb(5,5)
    complex(dp) :: dena,denb

    dena = one/(za(i1,i2)*za(i2,i3)*za(i4,i5))
    denb = one/(zb(i1,i2)*zb(i2,i3)*zb(i5,i4))

    !-- basic amplitudes
    !-- I inserted a zero just to switch off the photon radiation from the neutrino
    res(-1,+1,-1) = zero*za(i1,i4)**2*dena
    res(-1,-1,-1) = zb(i5,i3)**2*denb

    !-- 4 <-> 5
    res(-1,+1,+1) = -zero
    res(-1,-1,+1) = -zero
    res(+1,-1,+1) = zero
    res(+1,+1,+1) = zero
    res(+1,-1,-1) = -zero
    res(+1,+1,-1) = -zero

  end function master_amp_5pt_seconddiagram

  


end module mod_amplitudes_tree_ppnul
     

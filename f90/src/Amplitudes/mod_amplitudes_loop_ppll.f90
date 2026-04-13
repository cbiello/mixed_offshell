!-- loop amplituds, all from OpenLoops
!-- factorise out as/twopi [aem/twopi] and gs^2[eesq] per qcd[ew] loop and gluon[photon] emission
module mod_amplitudes_loop_ppll
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_ol_interface
  use mod_amplitudes_tree_ppll
  use mod_coupl
  use mod_auxfunctions
  use openloops
  implicit none
  private

  !-- general routines [work in progress]
  public :: res_qcdloop_qqb_gen
  public :: res_ewkloop_qqb_gen
  public :: ol_res_qcdloop_qqb_gen

  public :: res_ewkloop_g_qqb_gen
  public :: res_ewkloop_g_gq_gen
  public :: res_ewkloop_g_qg_gen
  public :: res_ewkloop_g_qqb_nf_gen
  public :: res_ewkloop_g_gq_nf_gen
  public :: res_ewkloop_g_qg_nf_gen
  public :: ol_res_qcdloop_a_qqb_gen, ol_res_qcdloop_a_aq_gen, ol_res_qcdloop_a_qa_gen

  !-- for qcdloop amplitudes, use mine
  public :: res_qcdloop_qqb
  public :: res_ewkloop_qqb,res_ewkloopAA_aa

  public :: res_ewkloop_g_qqb,res_ewkloop_g_gq,res_ewkloop_g_qg
  public :: res_qcdloop_a_qqb,res_qcdloop_a_aq,res_qcdloop_a_qa

  public :: res_ewkloop_qqb_nf,res_ewkloop_g_qqb_nf
  public :: res_ewkloop_g_gq_nf,res_ewkloop_g_qg_nf

  !-- OL QCD amplitudes, for checks
  public :: ol_res_qcdloop_qqb
  public :: ol_res_qcdloop_a_qqb,ol_res_qcdloop_a_aq,ol_res_qcdloop_a_qa

  public :: ew_renorm_nf_1L

  !-- EW counterterms from fermionic loop contributions
  complex(dp), save :: dZe,dZmz2,dZmw2
  real(dp), parameter :: munf = 200._dp !-- results will not depend on munf

contains

  subroutine ew_renorm_nf_1L()
    implicit none
    real(dp) :: mw2,mz2,mu2
    complex(dp) :: xtz,cw2,sw2
    complex(dp) :: SiTZZmz,dSiTZZ,SiTWWmw,dSiTWW,SiTW0

    mz2 = real(cms_mzsq,kind=dp)
    mw2 = real(cms_mwsq,kind=dp)
    mu2 = munf**2
    sw2 = cms_sw2
    cw2 = cms_cw2

    xtz = (mz2 + ci*sqrt(mz2*(4*mtsq-mz2)))/2/mtsq

    SiTZZmz = four/3 * mz2 * lfav * (five/3 + ci*pi + log(mu2/mz2)) &
            + four/3 * xn * (4*mtsq-mz2) * (one-xtz) * log(one-xtz) &
            * (acu**2 * (four-mz2/mtsq) - (two+mz2/mtsq)*vcu**2)/(xtz*(xtz-two)) &
            - four/3 * xn * (4*mtsq*(2*acu**2-vcu**2) - five/3*mz2*(acu**2+vcu**2) &
            + (6*acu**2*mtsq-mz2*(acu**2+vcu**2))*log(mu2/mtsq))

    dSiTZZ  = four/3 * lfav * (two/3 + ci*pi + log(mu2/mz2)) &
            + xn * (8*vcu**2 * (-mtsq/mz2 - mtsq**2/mz2**2 * 2 * xtz * log(one-xtz)/(xtz-two)) &
            + four/3 * (acu**2+vcu**2)*(two/3 * (-6._dp + 6*xtz+xtz**2)/xtz**2 + log(mu2/mtsq) &
            - (xtz-two)*(-two+2*xtz+xtz**2)*log(one-xtz)/xtz**3))

    SiTWWmw = mw2/3/sw2 * nle * (five/3 + ci*pi + log(mu2/mw2)) &
            + mw2/3/sw2 * xn * nlq * (five/3 + ci*pi + log(mu2/mw2)) &
            + mw2/3/sw2 * xn * ((-3*mtsq**2 - 6*mtsq*mw2 + 10*mw2**2)/(6*mw2**2) &
            + ((-3*mtsq + 2*mw2)*log(mu2/mtsq))/(2*mw2) &
            - ((mtsq + 2*mw2)*(mtsq**2 - 2*mtsq*mw2 + mw2**2)*log(one - mw2/mtsq))/(2*mw2**3))

    dSiTWW = (nle*(two + 3*ci*pi))/(9*sw2) &
           + ((6*mtsq**2 + 3*mtsq*mw2 + 2*mw2**2*(two + nlq*(two + 3*ci*pi)))*xn)/(18*mw2**2*sw2) &
           + (xn*log(mu2/mtsq))/(3*sw2) + (nle/(3*sw2) + (nlq*xn)/(3*sw2))*log(mu2/mw2) &
           + ((-one + mtsq**3/mw2**3)*xn*log(one - mw2/mtsq))/(3*sw2)

    SiTW0  = -xn*(mtsq + 2*mtsq*log(mu2/mtsq))/(4*sw2)

    dZmz2 = SiTZZmz + (cms_mzsq - mz2)*dSiTZZ
    dZmw2 = SiTWWmw + (cms_mwsq - mw2)*dSiTWW

    if(ew_scheme.eq.'Gmu') then
      dZe = (cw2/sw2*(dZmz2/cms_mzsq - dZmw2/cms_mwsq) - (SiTW0-dZmw2)/cms_mwsq)/2
      dZe = real(dZe,kind=dp)
    elseif(ew_scheme.eq.'amz') then
      print *,'a(mz) scheme for 1l-ew nf contributions not implemented yet'
    endif

  end subroutine

  !-----------------------------------------------------------------
  !--- amplitudes without extra radiation
  !-----------------------------------------------------------------

  ! conventions: 0 -> q(1) qb(2) l(3) lb(4)
  ! amp returned as -5:7 x -5:7 matrix
  ! used for NC, CC+, CC-

  subroutine res_qcdloop_qqb_gen(p,res0,res1fin)
    real(dp), intent(in)  :: p(4,4)
    real(dp), intent(out) :: res0(-5:7,-5:7),res1fin(-5:7,-5:7)

    call res_tree_qqb_gen(p,res0)
    res1fin = -8._dp*Cf*res0

  end subroutine res_qcdloop_qqb_gen


  ! old amplitudes for NC  
  !-- res(1,:) = q qb -> e- e [dn,up]
  !-- res(2,:) = qb q -> e- e [dn,up]
  subroutine res_qcdloop_qqb(p,res0,res1fin)
    real(dp), intent(in)  :: p(4,4)
    real(dp), intent(out) :: res0(2,2),res1fin(2,2)

    call res_tree_qqb(p,res0)
    res1fin = -8._dp*Cf*res0
    
  end subroutine res_qcdloop_qqb

  !-- res(1,:) = q qb -> e- e [dn,up,bot]
  !-- res(2,:) = qb q -> e- e [dn,up,bot]
  subroutine res_ewkloop_qqb(p,res0,res1fin)
    real(dp), intent(in)  :: p(4,4)
    real(dp), intent(out) :: res0(2,3),res1fin(2,3)
    real(dp15) :: p_ol(4,4)
    real(dp15) :: res1(0:2,2,3),acc(2,3),res0_ol(2,3)
    integer  :: ord(2,2),i,j
    real(dp) :: Lij(4,4),Lij2(4,4)

    p_ol(:,3:4) = p(:,3:4)
    
    ord(:,1) = [1,2] !-- q qb
    ord(:,2) = [2,1] !-- qb q

    do i = 1,2
       p_ol(:,1) = p(:,ord(1,i)); p_ol(:,2) = p(:,ord(2,i))
       do j = 1,3
          call evaluate_loop(OL_id(j),p_ol,res0_ol(i,j),res1(:,i,j),acc(i,j))
       enddo
    enddo
    
    res0 = res0_ol
    
    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(aem_ol/twopi)

    !-- now go to finite part
    call get_Lij(p,mu_ol,[1,2,3,4],Lij,Lij2)
    res1fin(1,:) = res1(0,1,:) - get_fin_ewk_qqbllb(2,1,3,4,Lij,Lij2)*res0(1,:)
    res1fin(2,:) = res1(0,2,:) - get_fin_ewk_qqbllb(1,2,3,4,Lij,Lij2)*res0(2,:)


#if(_withchecks == 1)
    if (checkpoles) then
       call check_poles('q qb->e- e+ ew',i1ewk_qqbllb(2,1,3,4,Lij,Lij2),res0(1,:),res1(:,1,:))
       call check_poles('qb q->e- e+ ew',i1ewk_qqbllb(1,2,3,4,Lij,Lij2),res0(2,:),res1(:,2,:))
    endif
#endif
    
  end subroutine res_ewkloop_qqb

  ! conventions: 0 -> q(1) qb(2) l(3) lb(4)
  ! amp returned as -5:7 x -5:7 matrix
  ! used for NC, CC+, CC-

  subroutine res_ewkloop_qqb_gen(p,res0loop,res1loopfin)
    real(dp), intent(in)  :: p(4,4)
    real(dp) :: res0(2,5),res1fin(2,5)
    real(dp15) :: p_ol(4,4)
    real(dp15) :: res1(0:2,2,5),acc(2,5),res0_ol(2,5)
    integer  :: ord(2,2),i,j
    real(dp) :: Lij(4,4),Lij2(4,4)
    real(dp) :: res1loopfin(-5:7,-5:7)
    real(dp) :: res0loop(-5:7,-5:7)
    !-- debug
    real(dp) :: finiteZdn,finiteZup
    real(dp) :: finiteWdu,finiteWud
    logical :: debug

    debug = .true.

    p_ol(:,3:4) = p(:,3:4)

    ord(:,1) = [1,2] !-- q qb
    ord(:,2) = [2,1] !-- qb q

#if (_Vcharge==0)
    do i = 1,2
       p_ol(:,1) = p(:,ord(1,i)); p_ol(:,2) = p(:,ord(2,i))
       do j = 1,3
          if(debug) then
             call evaluate_tree(OL_id(j),p_ol,res0_ol(i,j))
             res1(:,i,j)=0d0
          else
             call evaluate_loop(OL_id(j),p_ol,res0_ol(i,j),res1(:,i,j),acc(i,j))
          endif
       enddo
    enddo
#else
    do i = 1,1 ! no additional ord since we have defined all the process
               ! therefore we don't need to obtain the results via crossing
       p_ol(:,1) = p(:,ord(1,i)); p_ol(:,2) = p(:,ord(2,i))
       do j = 1,2
          if(debug) then
             call evaluate_tree(OL_id(j),p_ol,res0_ol(i,j))
             res1(:,1,j) = 0d0
          else
             call evaluate_loop(OL_id(j),p_ol,res0_ol(i,j),res1(:,i,j),acc(i,j))
          endif
       enddo
    enddo
#endif

    res0 = res0_ol

    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(aem_ol/twopi)
    
    !-- now go to finite part
    call get_Lij(p,mu_ol,[1,2,3,4],Lij,Lij2)
    
#if (_Vcharge == 0)
    
    finiteZdn = get_fin_ewk_qqbllb_gen(1,2,3,4,Lij,Lij2,Qdn,-Qdn,Q_lep,-Q_lep)
    finiteZup = get_fin_ewk_qqbllb_gen(1,2,3,4,Lij,Lij2,Qup,-Qup,Q_lep,-Q_lep)
    res1fin(1,1) = res1(0,1,1) - finiteZdn * res0(1,1)
    res1fin(1,2) = res1(0,1,2) - finiteZup * res0(1,2)
    res1fin(1,3) = res1(0,1,3) - finiteZdn * res0(1,3)

    finiteZdn = get_fin_ewk_qqbllb_gen(2,1,3,4,Lij,Lij2,Qdn,-Qdn,Q_lep,-Q_lep)
    finiteZup = get_fin_ewk_qqbllb_gen(2,1,3,4,Lij,Lij2,Qup,-Qup,Q_lep,-Q_lep)
    !--CB: equivalently (1,2,3,4,Lij,Lij2,-Qdn,Qdn,Q_lep,-Q_lep) and same for Qup
    res1fin(2,1) = res1(0,2,1) - finiteZdn * res0(2,1)
    res1fin(2,2) = res1(0,2,2) - finiteZup * res0(2,2)
    res1fin(2,3) = res1(0,2,3) - finiteZdn * res0(2,3)

#elif (_Vcharge == 1)

    !-1 2  -> 12 -11
    finiteWdu = get_fin_ewk_qqbllb_gen(1,2,3,4,Lij,Lij2,-Qdn,Qup,0d0,-Q_lep)
    res1fin(1,1) = res1(0,1,1) - finiteWdu * res0(1,1)
    !2 -1  -> 12 -11
    finiteWud = get_fin_ewk_qqbllb_gen(1,2,3,4,Lij,Lij2,Qup,-Qdn,0d0,-Q_lep)
    res1fin(1,2) = res1(0,1,2) - finiteWud * res0(1,2)

#elif (_Vcharge == -1)
    
    !1 -2  -> 11 -12
    finiteWdu = get_fin_ewk_qqbllb_gen(1,2,3,4,Lij,Lij2,Qdn,-Qup,Q_lep,0d0)
    res1fin(1,1) = res1(0,1,1) - finiteWdu * res0(1,1)
    !-2 1  -> 11 -12
    finiteWud = get_fin_ewk_qqbllb_gen(1,2,3,4,Lij,Lij2,-Qup,Qdn,Q_lep,0d0)
    res1fin(1,2) = res1(0,1,2) - finiteWud * res0(1,2)

#endif
    
    !--reshape the resfin matrix element to have the new shape
    res1loopfin=0
#if (_Vcharge == 0) 
    res1loopfin(-1,1) = res1fin(2,1)
    res1loopfin(-2,2) = res1fin(2,2)
    res1loopfin(-3,3) = res1fin(2,1)
    res1loopfin(-4,4) = res1fin(2,2)
    res1loopfin(-5,5) = res1fin(2,3)

    res1loopfin(1,-1) = res1fin(1,1)
    res1loopfin(2,-2) = res1fin(1,2)
    res1loopfin(3,-3) = res1fin(1,1)
    res1loopfin(4,-4) = res1fin(1,2)
    res1loopfin(5,-5) = res1fin(1,3)


    res0loop(-1,1) = res0(2,1)
    res0loop(-2,2) = res0(2,2)
    res0loop(-3,3) = res0(2,1)
    res0loop(-4,4) = res0(2,2)
    res0loop(-5,5) = res0(2,3)

    res0loop(1,-1) = res0(1,1)
    res0loop(2,-2) = res0(1,2)
    res0loop(3,-3) = res0(1,1)
    res0loop(4,-4) = res0(1,2)
    res0loop(5,-5) = res0(1,3)


#elif (_Vcharge == 1)
    res1loopfin(-1,2) = res1fin(1,1)
    res1loopfin(2,-1) = res1fin(1,2)
    res1loopfin(4,-3) = res1loopfin(2,-1)
    res1loopfin(-3,4) = res1loopfin(-1,2)

    res0loop(-1,2) = res0(1,1)
    res0loop(2,-1) = res0(1,2)
    res0loop(4,-3) = res0(1,2)
    res0loop(-3,4) = res0(1,1)

#elif (_Vcharge == -1)
    res1loopfin(1,-2) = res1fin(1,1)
    res1loopfin(-2,1) = res1fin(1,2)
    res1loopfin(-4,3) = res1loopfin(-2,1)
    res1loopfin(3,-4) = res1loopfin(1,-2)


    res0loop(1,-2) = res0(1,1)
    res0loop(-2,1) = res0(1,2)
    res0loop(-4,3) = res0(1,2)
    res0loop(3,-4) = res0(1,1)
#endif
  
  end subroutine res_ewkloop_qqb_gen


  !-- res = a a  -> e- e 
  !-- res(2,:) = qb q -> e- e [dn,up,bot]
  subroutine res_ewkloopAA_aa(p,res0,res1fin)
    real(dp), intent(in)  :: p(4,4)
    real(dp), intent(out) :: res0,res1fin
    real(dp15) :: p_ol(4,4)
    real(dp15) :: res1(0:2),acc,res0_ol
    real(dp) :: Lij(4,4),Lij2(4,4)
    logical :: debug

    debug=.true.

#if(_Vcharge==0)
    
    p_ol = p

    if(debug) then
       call evaluate_tree(OL_id(1),p_ol,res0_ol)
       res1=0d0
    else
       call evaluate_loop(OL_id(1),p_ol,res0_ol,res1,acc)
    endif

    
    res0 = real(res0_ol,kind=dp)
    
    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(aem_ol/twopi)

    !-- now go to finite part
    call get_Lij(p,mu_ol,[1,2],Lij,Lij2)
    res1fin = res1(0) - get_fin_ewk_aallb(1,2,Lij,Lij2)*res0

#else

    write(*,*), 'no aa process for CC DY'
    stop

#endif

#if(_withchecks == 1)
    if (checkpoles) then
       call check_poles_sing('q qb->e- e+ ew',i1ewk_aallb(1,2,Lij,Lij2),res0,res1)
   endif
#endif
    
  end subroutine res_ewkloopAA_aa

  !-----------------------------------------------------------------
  !--- amplitudes with extra radiation, ns
  !-----------------------------------------------------------------
  
  !-- res1loopfin(-5:7,-5:7) is a matrix in the flavour space
  !-- we use res1fin internally for historical (Z) reasons
  subroutine res_ewkloop_g_qqb_gen(p,res0,res1loopfin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res0(2,3), res1loopfin(-5:7,-5:7)
    real(dp) :: res1fin(2,3)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,3),acc(2,3),res0_ol(2,3)
    integer  :: ord(2,2),i,j
    real(dp) :: Lij(5,5),Lij2(5,5)
    real(dp) :: finiteZdn, finiteZup, finiteWdu, finiteWud

    p_ol(:,3:5) = p(:,3:5)
    
    ord(:,1) = [1,2] !-- q qb
    ord(:,2) = [2,1] !-- qb q

#if (_Vcharge==0)
    do i = 1,2 !for crossing symmetry
       p_ol(:,1) = p(:,ord(1,i)); p_ol(:,2) = p(:,ord(2,i))
       do j = 1,3
          call evaluate_loop(OL_id(j+3),p_ol,res0_ol(i,j),res1(:,i,j),acc(i,j))
       enddo
    enddo
#else
    do i = 1,1 !no needed crossing-symmetry
       p_ol(:,1) = p(:,ord(1,i)); p_ol(:,2) = p(:,ord(2,i))
       do j = 1,2
          call evaluate_loop(OL_id(j+2),p_ol,res0_ol(i,j),res1(:,i,j),acc(i,j))
       enddo
    enddo
#endif
    res0 = res0_ol

    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(two*aem_ol*as_ol)
    res0 = res0/(as_ol*four*pi)

    !-- now go to finite part
    call get_Lij(p,mu_ol,[1,2,3,4],Lij,Lij2)

#if (_Vcharge==0)

    finiteZdn = get_fin_ewk_qqbllb_gen(1,2,3,4,Lij,Lij2,Qdn,-Qdn,Q_lep,-Q_lep)
    finiteZup = get_fin_ewk_qqbllb_gen(1,2,3,4,Lij,Lij2,Qup,-Qup,Q_lep,-Q_lep)
    res1fin(1,1) = res1(0,1,1) - finiteZdn * res0(1,1)
    res1fin(1,2) = res1(0,1,2) - finiteZup * res0(1,2)
    res1fin(1,3) = res1(0,1,3) - finiteZdn * res0(1,3)

    finiteZdn = get_fin_ewk_qqbllb_gen(2,1,3,4,Lij,Lij2,Qdn,-Qdn,Q_lep,-Q_lep)
    !equivalently get_fin_ewk_qqbllb_gen(1,2,3,4,Lij,Lij2,-Qdn,Qdn,Q_lep,-Q_lep)
    !CB: I suggest to always use 1234 as order and touch the charges only
    finiteZup = get_fin_ewk_qqbllb_gen(2,1,3,4,Lij,Lij2,Qup,-Qup,Q_lep,-Q_lep)
    res1fin(2,1) = res1(0,2,1) - finiteZdn * res0(2,1)
    res1fin(2,2) = res1(0,2,2) - finiteZup * res0(2,2)
    res1fin(2,3) = res1(0,2,3) - finiteZdn * res0(2,3)
    
#elif (_Vcharge==1)

    !-1 2 -> 12 -11 21
    finiteWdu = get_fin_ewk_qqbllb_gen(1,2,3,4,Lij,Lij2,-Qdn,Qup,0d0,-Q_lep)
    res1fin(1,1) = res1(0,1,1) - finiteWdu * res0(1,1)
    !2 -1 -> 12 -11 21
    finiteWud = get_fin_ewk_qqbllb_gen(1,2,3,4,Lij,Lij2,Qup,-Qdn,0d0,-Q_lep)
    res1fin(1,2) = res1(0,1,2) - finiteWud * res0(1,2)
    
#elif (_Vcharge==-1)

    !1 -2 -> 11 -12 21
    finiteWdu = get_fin_ewk_qqbllb_gen(1,2,3,4,Lij,Lij2,Qdn,-Qup,Q_lep,0d0)
    res1fin(1,1) = res1(0,1,1) - finiteWdu * res0(1,1)
    !-2 1 -> 11 -12 21
    finiteWud = get_fin_ewk_qqbllb_gen(1,2,3,4,Lij,Lij2,-Qup,Qdn,Q_lep,0d0)
    res1fin(1,2) = res1(0,1,2) - finiteWud * res0(1,2)

#endif

    !--reshape the resfin matrix element to have the new shape
    res1loopfin=0
#if (_Vcharge == 0) 
    res1loopfin(-1,1) = res1fin(2,1)
    res1loopfin(-2,2) = res1fin(2,2)
    res1loopfin(-3,3) = res1fin(2,1)
    res1loopfin(-4,4) = res1fin(2,2)
    res1loopfin(-5,5) = res1fin(2,3)

    res1loopfin(1,-1) = res1fin(1,1)
    res1loopfin(2,-2) = res1fin(1,2)
    res1loopfin(3,-3) = res1fin(1,1)
    res1loopfin(4,-4) = res1fin(1,2)
    res1loopfin(5,-5) = res1fin(1,3)
#elif (_Vcharge == 1)
    res1loopfin(-1,2) = res1fin(1,1)
    res1loopfin(2,-1) = res1fin(1,2)
    res1loopfin(4,-3) = res1loopfin(2,-1)
    res1loopfin(-3,4) = res1loopfin(-1,2)
#elif (_Vcharge == -1)
    res1loopfin(1,-2) = res1fin(1,1)
    res1loopfin(-2,1) = res1fin(1,2)
    res1loopfin(-4,3) = res1loopfin(-2,1)
    res1loopfin(3,-4) = res1loopfin(1,-2)
#endif
    
  end subroutine res_ewkloop_g_qqb_gen
  
  !----- OLD CODE
  !-- res(1,:) = q qb -> e- e+ g [dn,up,bot]
  !-- res(2,:) = qb q -> e- e+ g [dn,up,bot]
  subroutine res_ewkloop_g_qqb(p,res0,res1fin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res0(2,3),res1fin(2,3)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,3),acc(2,3),res0_ol(2,3)
    integer  :: ord(2,2),i,j
    real(dp) :: Lij(5,5),Lij2(5,5)

    p_ol(:,3:5) = p(:,3:5)
    
    ord(:,1) = [1,2] !-- q qb
    ord(:,2) = [2,1] !-- qb q

    do i = 1,2
       p_ol(:,1) = p(:,ord(1,i)); p_ol(:,2) = p(:,ord(2,i))
       do j = 1,3
          call evaluate_loop(OL_id(j+3),p_ol,res0_ol(i,j),res1(:,i,j),acc(i,j))
       enddo
    enddo
    res0 = res0_ol

    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(two*aem_ol*as_ol)
    res0 = res0/(as_ol*four*pi)

    !-- now go to finite part
    call get_Lij(p,mu_ol,[1,2,3,4],Lij,Lij2)
    res1fin(1,:) = res1(0,1,:) - get_fin_ewk_qqbllb(2,1,3,4,Lij,Lij2)*res0(1,:)
    res1fin(2,:) = res1(0,2,:) - get_fin_ewk_qqbllb(1,2,3,4,Lij,Lij2)*res0(2,:)

#if(_withchecks == 1)
    if (checkpoles) then
       call check_poles('q qb->e- e+ g ew',i1ewk_qqbllb(2,1,3,4,Lij,Lij2),res0(1,:),res1(:,1,:))
       call check_poles('qb q->e- e+ g ew',i1ewk_qqbllb(1,2,3,4,Lij,Lij2),res0(2,:),res1(:,2,:))
    endif
#endif
    
  end subroutine res_ewkloop_g_qqb

  !!!!!!!!!!!!!!!!!!
  !TODO: QCD TO DO!!
  !!!!!!!!!!!!!!!!!!
  !-- res(1,:) = [d db, db d] -> e- e+ a
  !-- res(2,:) = [u ub, ub u] -> e- e+ a
  subroutine res_qcdloop_a_qqb(p,res0,res1fin)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res0(2,2),res1fin(2,2)
    integer, parameter :: iconf(5,2) = reshape([2,5,1,3,4, 1,5,2,3,4],[5,2])

    call res_qcdloop_j(p,iconf,aveqq,res0,res1fin)

  end subroutine res_qcdloop_a_qqb

  !-----------------------------------------------------------------
  !--- amplitudes with extra radiation, gq and qg
  !-----------------------------------------------------------------

  !-- gq amplitudes oneloop ewk
  !-- res1loopfin(-5:7,-5:7) is a matrix in the flavour space
  subroutine res_ewkloop_g_gq_gen(p,res0,res1loopfin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res0(2,3),res1loopfin(-5:7,-5:7)
    real(dp) :: res1fin(2,3)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,3),acc(2,3),res0_ol(2,3)
    integer  :: j
    real(dp) :: Lij(5,5),Lij2(5,5)
    real(dp) :: finiteZdn, finiteZup, finiteWdu, finiteWud

    p_ol(:,1:5) = p(:,1:5)

#if (_Vcharge == 0)
    !antiquarks in is
    do j = 1,3
       call evaluate_loop(OL_id(j+3),p_ol,res0_ol(1,j),res1(:,1,j),acc(1,j))
    enddo
    !quarks in is
    do j = 1,3
       call evaluate_loop(OL_id(j+6),p_ol,res0_ol(2,j),res1(:,2,j),acc(2,j))
    enddo
#else
    !up (or anti-up) in is
    call evaluate_loop(OL_id(3),p_ol,res0_ol(1,1),res1(:,1,1),acc(1,1))
    !anti-down (or down) in is
    call evaluate_loop(OL_id(4),p_ol,res0_ol(1,2),res1(:,1,2),acc(1,2))
#endif
    
    res0 = res0_ol
    
    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(two*aem_ol*as_ol)
    res0 = res0/(as_ol*four*pi)

    !-- now go to finite part
    call get_Lij(p,mu_ol,[2,3,4,5],Lij,Lij2)

#if (_Vcharge == 0)

    !21 -1 -> 11 -11 -1  |---> 1 -1 -> 11 -11
    finiteZdn = get_fin_ewk_qqbllb_gen(5,2,3,4,Lij,Lij2,Qdn,-Qdn,Q_lep,-Q_lep)    
    finiteZup = get_fin_ewk_qqbllb_gen(5,2,3,4,Lij,Lij2,Qup,-Qup,Q_lep,-Q_lep)
    res1fin(1,1) = res1(0,1,1) - finiteZdn * res0(1,1)
    res1fin(1,2) = res1(0,1,2) - finiteZup * res0(1,2)
    res1fin(1,3) = res1(0,1,3) - finiteZdn * res0(1,3)
    ! print*, 'finiteZdn= ', finiteZdn
    ! print*, 'finiteZup= ', finiteZup
    ! print*, 'old= ', get_fin_ewk_qqbllb(2,5,3,4,Lij,Lij2)

    !21 1 -> 11 -11 1  |---> -1 1 -> 11 -11
    finiteZdn = get_fin_ewk_qqbllb_gen(5,2,3,4,Lij,Lij2,-Qdn,Qdn,Q_lep,-Q_lep)
    finiteZup =	get_fin_ewk_qqbllb_gen(5,2,3,4,Lij,Lij2,-Qup,Qup,Q_lep,-Q_lep)
    res1fin(2,1) = res1(0,2,1) - finiteZdn * res0(2,1)
    res1fin(2,2) = res1(0,2,2) - finiteZup * res0(2,2)
    res1fin(2,3) = res1(0,2,3) - finiteZdn * res0(2,3)
    ! print*, 'finiteZdn= ', finiteZdn
    ! print*, 'finiteZup= ', finiteZup
    ! print*, 'old= ', get_fin_ewk_qqbllb(5,2,3,4,Lij,Lij2)
    
#elif (_Vcharge == 1)

    !21 2 -> 12 -11 1
    finiteWdu = get_fin_ewk_qqbllb_gen(5,2,3,4,Lij,Lij2,-Qdn,Qup,0d0,-Q_lep)
    res1fin(1,1) = res1(0,1,1) - finiteWdu * res0(1,1)
    !21 -1 -> 12 -11 -2
    finiteWud = get_fin_ewk_qqbllb_gen(5,2,3,4,Lij,Lij2,Qup,-Qdn,0d0,-Q_lep)
    res1fin(1,2) = res1(0,1,2) - finiteWud * res0(1,2)

#elif (_Vcharge == -1)

    !21 -2 -> 11 -12 -1
    finiteWdu = get_fin_ewk_qqbllb_gen(5,2,3,4,Lij,Lij2,Qdn,-Qup,Q_lep,0d0)
    res1fin(1,1) = res1(0,1,1) - finiteWdu * res0(1,1)
    !21 1 -> 11 -12 2
    finiteWud = get_fin_ewk_qqbllb_gen(5,2,3,4,Lij,Lij2,-Qup,Qdn,Q_lep,0d0)
    res1fin(1,2) = res1(0,1,2) - finiteWud * res0(1,2)
    
#endif

    !--reshape the resfin matrix element to have the new shape
    res1loopfin=0
#if (_Vcharge == 0) 
    res1loopfin(0,-1) = res1fin(1,1)
    res1loopfin(0,-2) = res1fin(1,2)
    res1loopfin(0,-3) = res1fin(1,1)
    res1loopfin(0,-4) = res1fin(1,2)
    res1loopfin(0,-5) = res1fin(1,3)
    !
    res1loopfin(0,1) = res1fin(2,1)
    res1loopfin(0,2) = res1fin(2,2)
    res1loopfin(0,3) = res1fin(2,1)
    res1loopfin(0,4) = res1fin(2,2)
    res1loopfin(0,5) = res1fin(2,3)
#elif (_Vcharge == 1)
    res1loopfin(0,2) = res1fin(1,1)
    res1loopfin(0,-1) = res1fin(1,2)
    res1loopfin(0,-3) = res1loopfin(0,-1)
    res1loopfin(0,4) = res1loopfin(0,2)
#elif (_Vcharge == -1)
    res1loopfin(0,-2) = res1fin(1,1)
    res1loopfin(0,1) = res1fin(1,2)
    res1loopfin(0,3) = res1loopfin(0,1)
    res1loopfin(0,-4) = res1loopfin(0,-2)
#endif
    
  end subroutine res_ewkloop_g_gq_gen

  !-- res(1,:) = g qb -> e- e+ qb [dn,up,bot]
  !-- res(2,:) = g q  -> e- e+ q  [dn,up,bot]
  subroutine res_ewkloop_g_gq(p,res0,res1fin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res0(2,3),res1fin(2,3)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,3),acc(2,3),res0_ol(2,3)
    integer  :: j
    real(dp) :: Lij(5,5),Lij2(5,5)

    p_ol(:,1:5) = p(:,1:5)
    
    !quark in is
    do j = 1,3
       call evaluate_loop(OL_id(j+3),p_ol,res0_ol(1,j),res1(:,1,j),acc(1,j))
    enddo
    !antiquarks in is
    do j = 1,3
       call evaluate_loop(OL_id(j+6),p_ol,res0_ol(2,j),res1(:,2,j),acc(2,j))
    enddo
    
    res0 = res0_ol
    
    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(two*aem_ol*as_ol)
    res0 = res0/(as_ol*four*pi)

    !-- now go to finite part
    call get_Lij(p,mu_ol,[2,3,4,5],Lij,Lij2)
    res1fin(1,:) = res1(0,1,:) - get_fin_ewk_qqbllb(2,5,3,4,Lij,Lij2)*res0(1,:)
    res1fin(2,:) = res1(0,2,:) - get_fin_ewk_qqbllb(5,2,3,4,Lij,Lij2)*res0(2,:)

#if(_withchecks == 1)
    if (checkpoles) then
       call check_poles('g qb->e- e+ qb ew',i1ewk_qqbllb(2,5,3,4,Lij,Lij2),res0(1,:),res1(:,1,:))
       call check_poles('g  q->e- e+ q  ew',i1ewk_qqbllb(5,2,3,4,Lij,Lij2),res0(2,:),res1(:,2,:))
    endif
#endif
    
  end subroutine res_ewkloop_g_gq

  !!!!!!!!!!!!!!!!!!!!!
  ! g_qg_gen
  ! res1loopfin(-5:7,-5:7) is a matrix in the flavour space 
  subroutine res_ewkloop_g_qg_gen(p,res0,res1loopfin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res0(2,3)
    real(dp)   :: res1fin(2,3)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,3),acc(2,3),res0_ol(2,3)
    integer  :: j
    real(dp) :: Lij(5,5),Lij2(5,5)
    real(dp), intent(out) :: res1loopfin(-5:7,-5:7)
    real(dp) :: finiteZdn, finiteZup, finiteWdu, finiteWud

    p_ol(:,1:5) = p(:,1:5)

#if (_Vcharge == 0)
    !quarks
    do j = 1,3
       call evaluate_loop(OL_id(j+3),p_ol,res0_ol(1,j),res1(:,1,j),acc(1,j))
    enddo
    !antiquarks
    do j = 1,3
       call evaluate_loop(OL_id(j+6),p_ol,res0_ol(2,j),res1(:,2,j),acc(2,j))
    enddo
#else
    !ol_id(3)= 2 21 -> 12 -11 1  OR  -1 21 -> 12 -11 -2
    !up (or anti-up) in is
    call evaluate_loop(OL_id(3),p_ol,res0_ol(1,1),res1(:,1,1),acc(1,1))
    !anti-down (or down) in is
    !ol_id(4)= -1 21 -> 12 -11 -2  OR  1 21 -> 11 -12 2
    call evaluate_loop(OL_id(4),p_ol,res0_ol(1,2),res1(:,1,2),acc(1,2))
#endif
    
    res0 = res0_ol
       
    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(two*aem_ol*as_ol)
    res0 = res0/(as_ol*four*pi)

    !-- now go to finite part
    call get_Lij(p,mu_ol,[1,3,4,5],Lij,Lij2)

#if (_Vcharge == 0)    
    !1 21 -> 11 -11 -1  |---> 1 -1 -> 11 -11
    finiteZdn = get_fin_ewk_qqbllb_gen(1,5,3,4,Lij,Lij2,Qdn,-Qdn,Q_lep,-Q_lep)    
    finiteZup = get_fin_ewk_qqbllb_gen(1,5,3,4,Lij,Lij2,Qup,-Qup,Q_lep,-Q_lep)
    res1fin(1,1) = res1(0,1,1) - finiteZdn * res0(1,1)
    res1fin(1,2) = res1(0,1,2) - finiteZup * res0(1,2)
    res1fin(1,3) = res1(0,1,3) - finiteZdn * res0(1,3)
    !print*, 'new dn= ', get_fin_ewk_qqbllb_gen(1,5,3,4,Lij,Lij2,Qdn,-Qdn,Q_lep,-Q_lep)
    !print*, 'new up= ', get_fin_ewk_qqbllb_gen(1,5,3,4,Lij,Lij2,Qup,-Qup,Q_lep,-Q_lep)
    !print*, 'old= ', get_fin_ewk_qqbllb(5,1,3,4,Lij,Lij2)

    !-1 21 -> 11 -11 1  |---> -1 1 -> 11 -11
    finiteZdn = get_fin_ewk_qqbllb_gen(1,5,3,4,Lij,Lij2,-Qdn,Qdn,Q_lep,-Q_lep)
    finiteZup =	get_fin_ewk_qqbllb_gen(1,5,3,4,Lij,Lij2,-Qup,Qup,Q_lep,-Q_lep)
    res1fin(2,1) = res1(0,2,1) - finiteZdn * res0(2,1)
    res1fin(2,2) = res1(0,2,2) - finiteZup * res0(2,2)
    res1fin(2,3) = res1(0,2,3) - finiteZdn * res0(2,3)
    !print*, 'finiteZdn= ', finiteZdn
    !print*, 'finiteZup= ', finiteZup
    !print*, 'old= ', get_fin_ewk_qqbllb(1,5,3,4,Lij,Lij2)
    
#elif (_Vcharge == 1)

    !2 21 -> 12 -11 1
    finiteWdu = get_fin_ewk_qqbllb_gen(1,5,3,4,Lij,Lij2,Qup,-Qdn,0d0,-Q_lep)
    res1fin(1,1) = res1(0,1,1) - finiteWdu * res0(1,1)
    !-1 21 -> 12 -11 -2
    finiteWud = get_fin_ewk_qqbllb_gen(1,5,3,4,Lij,Lij2,-Qdn,Qup,0d0,-Q_lep)
    res1fin(1,2) = res1(0,1,2) - finiteWud * res0(1,2)

#elif (_Vcharge == -1)

    !-2 21 -> 11 -12 -1
    finiteWdu = get_fin_ewk_qqbllb_gen(1,5,3,4,Lij,Lij2,-Qup,Qdn,Q_lep,0d0)
    res1fin(1,1) = res1(0,1,1) - finiteWdu * res0(1,1)
    !1 21 -> 11 -12 2
    finiteWud = get_fin_ewk_qqbllb_gen(1,5,3,4,Lij,Lij2,Qdn,-Qup,Q_lep,0d0)
    res1fin(1,2) = res1(0,1,2) - finiteWud * res0(1,2)
    
#endif

    ! FIlling the res1loopfin in flavour space
    res1loopfin=0
#if (_Vcharge == 0) 
    res1loopfin(1,0) = res1fin(1,1)
    res1loopfin(2,0) = res1fin(1,2)
    res1loopfin(3,0) = res1fin(1,1)
    res1loopfin(4,0) = res1fin(1,2)
    res1loopfin(5,0) = res1fin(1,3)
    !
    res1loopfin(-1,0) = res1fin(2,1)
    res1loopfin(-2,0) = res1fin(2,2)
    res1loopfin(-3,0) = res1fin(2,1)
    res1loopfin(-4,0) = res1fin(2,2)
    res1loopfin(-5,0) = res1fin(2,3)
#elif (_Vcharge == 1)
    res1loopfin(2,0) = res1fin(1,1)
    res1loopfin(-1,0) = res1fin(1,2)
    res1loopfin(-3,0) = res1loopfin(-1,0)
    res1loopfin(4,0) = res1loopfin(2,0)
#elif (_Vcharge == -1)
    res1loopfin(-2,0) = res1fin(1,1)
    res1loopfin(1,0) = res1fin(1,2)
    res1loopfin(3,0) = res1loopfin(1,0)
    res1loopfin(-4,0) = res1loopfin(-2,0)
#endif
    
  end subroutine res_ewkloop_g_qg_gen

  !-- res(1,:) = q g  -> e- e+ q  [dn,up,bot]
  !-- res(2,:) = qb g -> e- e+ qb [dn,up,bot]
  subroutine res_ewkloop_g_qg(p,res0,res1fin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res0(2,3),res1fin(2,3)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,3),acc(2,3),res0_ol(2,3)
    integer  :: j
    real(dp) :: Lij(5,5),Lij2(5,5)

    p_ol(:,1:5) = p(:,1:5)
    
    do j = 1,3
       call evaluate_loop(OL_id(j+3),p_ol,res0_ol(1,j),res1(:,1,j),acc(1,j))
    enddo

    do j = 1,3
       call evaluate_loop(OL_id(j+6),p_ol,res0_ol(2,j),res1(:,2,j),acc(2,j))
    enddo
    res0 = res0_ol
       
    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(two*aem_ol*as_ol)
    res0 = res0/(as_ol*four*pi)

    !-- now go to finite part
    call get_Lij(p,mu_ol,[1,3,4,5],Lij,Lij2)
    res1fin(1,:) = res1(0,1,:) - get_fin_ewk_qqbllb(5,1,3,4,Lij,Lij2)*res0(1,:)
    res1fin(2,:) = res1(0,2,:) - get_fin_ewk_qqbllb(1,5,3,4,Lij,Lij2)*res0(2,:)

#if(_withchecks == 1)
    if (checkpoles) then
       call check_poles('q  g->e- e+ q  ew',i1ewk_qqbllb(5,1,3,4,Lij,Lij2),res0(1,:),res1(:,1,:))
       call check_poles('qb g->e- e+ qb ew',i1ewk_qqbllb(1,5,3,4,Lij,Lij2),res0(2,:),res1(:,2,:))
    endif
#endif
    
  end subroutine res_ewkloop_g_qg

  !----------------------------------------------------------------------------
  !--- 1L-EWK amplitudes with corrections from closed fermionic loops only
  !----------------------------------------------------------------------------

  subroutine res_ewkloop_qqb_nf(p,res1fin)
    real(dp), intent(in)  :: p(4,4)
    real(dp), intent(out) :: res1fin(2,2)
    integer  :: i
    real(dp) :: s,t,u,mu2
    complex(dp) :: xt,lflog,SiTAA,SiTAZ,SiTZZ,KK(3)
    complex(dp) :: cvv(2),cva(2),cav(2),caa(2),vq(2),aq(2)
    real(dp) :: Qq(2),Iq(2)
    complex(dp) :: HDY00(4,2,2),HDY01(4,2,2)

    s =  2*scr(p(:,1),p(:,2))
    t = -2*scr(p(:,1),p(:,3))
    u = -s -t
    mu2 = munf**2

    Qq = [Qdn,Qup]
    Iq = [I3dn,I3up]
    vq = [vcd,vcu]
    aq = [acd,acu]

    if (s.lt.four*mtsq) then
      xt = (s + ci*sqrt(-(s*(-4*mtsq + s))))/(2*mtsq)
    else
      xt = (s + sqrt(s*(-4*mtsq + s)))/(2*mtsq)
    endif

    lflog = five/3 + ci*pi + log(mu2/s)

    SiTAA = four/3 * xn * Qup2*((mtsq*(-four + 4*xt + (5*xt**2)/3))/(xt-one) &
          + s*log(mu2/mtsq) - (mtsq*(-two + xt)*(-two + 2*xt + xt**2)*log(one - xt))/((-one + xt)*xt)) &
          + four/3 * s * lflog * (xn*ndn*Qdn2 + xn*nup*Qup2 + nle*Qel2)

    SiTAZ = -four/3 * vcu * xn * Qup * (4*mtsq + s*five/3 + s*log(mu2/mtsq) &
          + ((4*mtsq - s)*(2*mtsq + s)*(-one + xt)*log(one - xt))/(mtsq*(-two + xt)*xt)) &
          - four/3 * s * lflog * (xn*ndn*vcd*Qdn + xn*nup*vcu*Qup + nle*vce*Qel)

    SiTZZ = four/3 * s * lfav * lflog &
          - four/3 * xn * (acu**2*(8*mtsq - s*five/3)-(4*mtsq + s*five/3)*vcu**2&
          +(acu**2*(6*mtsq-s)-s*vcu**2)*log(mu2/mtsq)) &
          - four/3 * xn * (((four - s/mtsq)*(acu**2*(4*mtsq - s) &
          - (2*mtsq + s)*vcu**2)*(-one + xt)*log(one - xt))/((-two + xt)*xt))
    

    KK(1) = SiTAA/s - 2*dZe
    KK(2) = cms_cw/cms_sw * (dZmz2/cms_mzsq - dZmw2/cms_mwsq) - SiTAZ/s
    KK(3) = (SiTZZ-dZmz2)/(s-cms_mzsq) - 2*dZe + (cms_cw2-cms_sw2)/cms_sw2 * (dZmz2/cms_mzsq - dZmw2/cms_mwsq)

    do i = 1, 2
      cvv(i) = -Qel*Qq(i) * KK(1)/s - (Qq(i)*vce + Qel*vq(i))/(s-cms_mzsq) * KK(2) - vq(i)*vce/(s-cms_mzsq) * KK(3)
      cva(i) = ace*Qq(i)/(s-cms_mzsq) * KK(2) + ace*vq(i)/(s-cms_mzsq) * KK(3)
      cav(i) = aq(i)*Qel/(s-cms_mzsq) * KK(2) + aq(i)*vce/(s-cms_mzsq) * KK(3)
      caa(i) = -aq(i)*ace/(s-cms_mzsq) * KK(3)
    enddo

    !! qqx channel
    do i = 1, 2
      !-- -+-+
      HDY00(1,i,1) = -2*Qq(i)*u/s-(one/6/cms_cw2 + Iq(i)/cms_sw2-2*Qq(i))*u/(s-cms_mzsq)
      HDY01(1,i,1) =  2*u*(cvv(i)-cav(i)-cva(i)+caa(i))
      !-- +--+
      HDY00(2,i,1) =  2*Qq(i)*t/s - Qq(i)*(two - one/cms_cw2)*t/(s - cms_mzsq)
      HDY01(2,i,1) = -2*t*(cvv(i)+cav(i)-cva(i)-caa(i))
      !-- -++-
      HDY00(3,i,1) = 2*Qq(i)*t/s + (one/3/cms_cw2 - 2*Qq(i))*t/(s - cms_mzsq)
      HDY01(3,i,1) = -2*t*(cvv(i)-cav(i)+cva(i)-caa(i))
      !-- +-+-
      HDY00(4,i,1) = -2*Qq(i)*u/s - 2*Qq(i)*(cms_sw2/cms_cw2)*u/(s - cms_mzsq)
      HDY01(4,i,1) = 2*u*(cvv(i)+cav(i)+cva(i)+caa(i))
    enddo

    !! qxq channel
    t = -two*scr(p(:,1),p(:,4))
    u = -s-t
    do i = 1, 2
      !-- -+-+
      HDY00(1,i,2) = -2*Qq(i)*u/s-(one/6/cms_cw2 + Iq(i)/cms_sw2-2*Qq(i))*u/(s-cms_mzsq)
      HDY01(1,i,2) =  2*u*(cvv(i)-cav(i)-cva(i)+caa(i))
      !-- +--+
      HDY00(2,i,2) =  2*Qq(i)*t/s - Qq(i)*(two - one/cms_cw2)*t/(s - cms_mzsq)
      HDY01(2,i,2) = -2*t*(cvv(i)+cav(i)-cva(i)-caa(i))
      !-- -++-
      HDY00(3,i,2) = 2*Qq(i)*t/s + (one/3/cms_cw2 - 2*Qq(i))*t/(s - cms_mzsq)
      HDY01(3,i,2) = -2*t*(cvv(i)-cav(i)+cva(i)-caa(i))
      !-- +-+-
      HDY00(4,i,2) = -2*Qq(i)*u/s - 2*Qq(i)*(cms_sw2/cms_cw2)*u/(s - cms_mzsq)
      HDY01(4,i,2) = 2*u*(cvv(i)+cav(i)+cva(i)+caa(i))
    enddo

    res1fin(1,1) = real(dot_product(HDY00(:,1,1),HDY01(:,1,1)),kind=dp)
    res1fin(2,1) = real(dot_product(HDY00(:,1,2),HDY01(:,1,2)),kind=dp)
    res1fin(1,2) = real(dot_product(HDY00(:,2,1),HDY01(:,2,1)),kind=dp)
    res1fin(2,2) = real(dot_product(HDY00(:,2,2),HDY01(:,2,2)),kind=dp)

    res1fin = res1fin * xn * eesq2 * aveqq

  end subroutine res_ewkloop_qqb_nf

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! nf amps from OL !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine res_ewkloop_g_qqb_nf_gen(p,res1loopfin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res1loopfin(-5:7,-5:7)
    real(dp)   :: res1fin(2,2)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: ord(2,2),i,j

#if (_Vcharge == 0)
    p_ol(:,3:5) = p(:,3:5)
    ord(:,1) = [1,2] !-- q qb
    ord(:,2) = [2,1] !-- qb q
    do i = 1,2
       p_ol(:,1) = p(:,ord(1,i)); p_ol(:,2) = p(:,ord(2,i))
       do j = 1,2
          call evaluate_loop(OL_id(j+3),p_ol,res0_ol(i,j),res1(:,i,j),acc(i,j))
       enddo
    enddo
#else
    p_ol(:,1:5) = p(:,1:5)
    !-1 2 -> 12 -11 21 or 1 -2 -> 11 -12 21
    call evaluate_loop(OL_id(3),p_ol,res0_ol(1,1),res1(:,1,1),acc(1,1))
    !2 -1 -> 12 -11 21 or -2 1 -> 11 -12 21
    call evaluate_loop(OL_id(4),p_ol,res0_ol(1,2),res1(:,1,2),acc(1,2))    
#endif

    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(aem_ol/twopi)/(4*pi*as_ol)
    res1fin(1,:) = res1(0,1,:)
    res1fin(2,:) = res1(0,2,:)

    !-- filling the res1loopfin in flavour space
    res1loopfin=0
#if (_Vcharge == 0) 
    res1loopfin(1,-1) = res1fin(1,1)
    res1loopfin(2,-2) = res1fin(1,2)
    res1loopfin(3,-3) = res1fin(1,1)
    res1loopfin(4,-4) = res1fin(1,2)
    res1loopfin(5,-5) = res1fin(1,1)
    !
    res1loopfin(-1,1) = res1fin(2,1)
    res1loopfin(-2,2) = res1fin(2,2)
    res1loopfin(-3,3) = res1fin(2,1)
    res1loopfin(-4,4) = res1fin(2,2)
    res1loopfin(-5,5) = res1fin(2,1) !no diff in the nf part
#elif (_Vcharge == 1)
    res1loopfin(-1,2) = res1fin(1,1)
    res1loopfin(2,-1) = res1fin(1,2)
    res1loopfin(-3,4) = res1loopfin(-1,2)
    res1loopfin(4,-3) = res1loopfin(2,-1)
#elif (_Vcharge == -1)
    res1loopfin(1,-2) = res1fin(1,1)
    res1loopfin(-2,1) = res1fin(1,2)
    res1loopfin(3,-4) = res1loopfin(1,-2)
    res1loopfin(-4,3) = res1loopfin(-2,1)
#endif

  end subroutine res_ewkloop_g_qqb_nf_gen
  
  
  subroutine res_ewkloop_g_qqb_nf(p,res1fin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res1fin(2,2)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: ord(2,2),i,j

    p_ol(:,3:5) = p(:,3:5)
    
    ord(:,1) = [1,2] !-- q qb
    ord(:,2) = [2,1] !-- qb q

    do i = 1,2
       p_ol(:,1) = p(:,ord(1,i)); p_ol(:,2) = p(:,ord(2,i))
       do j = 1,2
          call evaluate_loop(OL_id(j+3),p_ol,res0_ol(i,j),res1(:,i,j),acc(i,j))
       enddo
    enddo

    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(aem_ol/twopi)/(4*pi*as_ol)

    res1fin(1,:) = res1(0,1,:)
    res1fin(2,:) = res1(0,2,:)

  end subroutine res_ewkloop_g_qqb_nf


  subroutine res_ewkloop_g_gq_nf_gen(p,res1loopfin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res1loopfin(-5:7,-5:7)
    real(dp)   :: res1fin(2,2)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: j

    p_ol(:,1:5) = p(:,1:5)

#if (_Vcharge == 0)
    !anti-quark
    do j = 1,2
       call evaluate_loop(OL_id(j+3),p_ol,res0_ol(1,j),res1(:,1,j),acc(1,j))
    enddo
    !quarks
    do j = 1,2
       call evaluate_loop(OL_id(j+6),p_ol,res0_ol(2,j),res1(:,2,j),acc(2,j))
    enddo
#else
    !21 2 -> 12 -11 1 or 21 -2 -> 11 -12 -1
    call evaluate_loop(OL_id(3),p_ol,res0_ol(1,1),res1(:,1,1),acc(1,1))
    !21 -1 -> 12 -11 -2 or 21 1 -> 11 -12 2
    call evaluate_loop(OL_id(4),p_ol,res0_ol(1,2),res1(:,1,2),acc(1,2))
#endif

    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(aem_ol/twopi)/(4*pi*as_ol)
    res1fin(1,:) = res1(0,1,:)
    res1fin(2,:) = res1(0,2,:)

    !-- filling the res1loopfin in flavour space
    res1loopfin=0
#if (_Vcharge == 0) 
    res1loopfin(0,-1) = res1fin(1,1)
    res1loopfin(0,-2) = res1fin(1,2)
    res1loopfin(0,-3) = res1fin(1,1)
    res1loopfin(0,-4) = res1fin(1,2)
    res1loopfin(0,-5) = res1fin(1,1)
    !
    res1loopfin(0,1) = res1fin(2,1)
    res1loopfin(0,2) = res1fin(2,2)
    res1loopfin(0,3) = res1fin(2,1)
    res1loopfin(0,4) = res1fin(2,2)
    res1loopfin(0,5) = res1fin(2,1) !no diff in the nf part
#elif (_Vcharge == 1)
    res1loopfin(0,2) = res1fin(1,1)
    res1loopfin(0,-1) = res1fin(1,2)
    res1loopfin(0,4) = res1loopfin(0,2)
    res1loopfin(0,-3) = res1loopfin(0,-1)
#elif (_Vcharge == -1)
    res1loopfin(0,-2) = res1fin(1,1)
    res1loopfin(0,1) = res1fin(1,2)
    res1loopfin(0,-4) = res1loopfin(0,-2)
    res1loopfin(0,3) = res1loopfin(0,1)
#endif    
    
  end subroutine res_ewkloop_g_gq_nf_gen

  
  subroutine res_ewkloop_g_gq_nf(p,res1fin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res1fin(2,2)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: j

    p_ol(:,1:5) = p(:,1:5)
    
    do j = 1,2
       call evaluate_loop(OL_id(j+3),p_ol,res0_ol(1,j),res1(:,1,j),acc(1,j))
    enddo

    do j = 1,2
       call evaluate_loop(OL_id(j+6),p_ol,res0_ol(2,j),res1(:,2,j),acc(2,j))
    enddo

    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(aem_ol/twopi)/(4*pi*as_ol)

    res1fin(1,:) = res1(0,1,:)
    res1fin(2,:) = res1(0,2,:)
    
  end subroutine res_ewkloop_g_gq_nf

  subroutine res_ewkloop_g_qg_nf_gen(p,res1loopfin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res1loopfin(-5:7,-5:7)
    real(dp)   :: res1fin(2,2)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: j

    p_ol(:,1:5) = p(:,1:5)

#if (_Vcharge == 0)
    !quark
    do j = 1,2
       call evaluate_loop(OL_id(j+3),p_ol,res0_ol(1,j),res1(:,1,j),acc(1,j))
    enddo
    !anti-quark
    do j = 1,2
       call evaluate_loop(OL_id(j+6),p_ol,res0_ol(2,j),res1(:,2,j),acc(2,j))
    enddo
#else
    !2 21 -> 12 -11 1 or -2 21 -> 11 -12 -1
    call evaluate_loop(OL_id(3),p_ol,res0_ol(1,1),res1(:,1,1),acc(1,1))
    !-1 21 -> 12 -11 -2 or 1 21 -> 11 -12 2
    call evaluate_loop(OL_id(4),p_ol,res0_ol(1,2),res1(:,1,2),acc(1,2))
#endif
       
    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(aem_ol/twopi)/(4*pi*as_ol)
    res1fin(1,:) = res1(0,1,:)
    res1fin(2,:) = res1(0,2,:)
    
    !-- filling the res1loopfin in flavour space
    res1loopfin=0
#if (_Vcharge == 0) 
    res1loopfin(1,0) = res1fin(1,1)
    res1loopfin(2,0) = res1fin(1,2)
    res1loopfin(3,0) = res1fin(1,1)
    res1loopfin(4,0) = res1fin(1,2)
    res1loopfin(5,0) = res1fin(1,1)
    !
    res1loopfin(-1,0) = res1fin(2,1)
    res1loopfin(-2,0) = res1fin(2,2)
    res1loopfin(-3,0) = res1fin(2,1)
    res1loopfin(-4,0) = res1fin(2,2)
    res1loopfin(-5,0) = res1fin(2,1) !no diff in the nf part
#elif (_Vcharge == 1)
    res1loopfin(2,0) = res1fin(1,1)
    res1loopfin(-1,0) = res1fin(1,2)
    res1loopfin(4,0) = res1loopfin(2,0)
    res1loopfin(-3,0) = res1loopfin(-1,0)
#elif (_Vcharge == -1)
    res1loopfin(-2,0) = res1fin(1,1)
    res1loopfin(1,0) = res1fin(1,2)
    res1loopfin(-4,0) = res1loopfin(-2,0)
    res1loopfin(3,0) = res1loopfin(1,0)
#endif    
    
  end subroutine res_ewkloop_g_qg_nf_gen

  
  subroutine res_ewkloop_g_qg_nf(p,res1fin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res1fin(2,2)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: j

    p_ol(:,1:5) = p(:,1:5)
    
    do j = 1,2
       call evaluate_loop(OL_id(j+3),p_ol,res0_ol(1,j),res1(:,1,j),acc(1,j))
    enddo

    do j = 1,2
       call evaluate_loop(OL_id(j+6),p_ol,res0_ol(2,j),res1(:,2,j),acc(2,j))
    enddo
       
    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(aem_ol/twopi)/(4*pi*as_ol)

    res1fin(1,:) = res1(0,1,:)
    res1fin(2,:) = res1(0,2,:)

  end subroutine res_ewkloop_g_qg_nf

  !-----------------------------------------------------------------
  !--- amplitudes with extra radiation, aq and qa
  !-----------------------------------------------------------------

  subroutine res_qcdloop_a_aq(p,res0,res1fin)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res0(2,2),res1fin(2,2)
    integer, parameter :: iconf(5,2) = reshape([2,1,5,3,4, 5,1,2,3,4],[5,2])

    call res_qcdloop_j(p,iconf,aveqa,res0,res1fin)

  end subroutine res_qcdloop_a_aq

  subroutine res_qcdloop_a_qa(p,res0,res1fin)
    real(dp), intent(in)  :: p(:,:)
    real(dp), intent(out) :: res0(2,2),res1fin(2,2)
    integer, parameter :: iconf(5,2) = reshape([5,2,1,3,4, 1,2,5,3,4],[5,2])

    call res_qcdloop_j(p,iconf,aveqa,res0,res1fin)

  end subroutine res_qcdloop_a_qa

  !-----------------------------------------------------------------
  !--- poles check
  !-----------------------------------------------------------------

#if(_withchecks == 1)
  
  subroutine check_poles(str,i1,res0,res1)
    character(*), intent(in) :: str
    real(dp), intent(in)   :: i1(-2:,:),res0(:)
    real(dp15), intent(in) :: res1(0:,:)
    real(dp) :: tmp
    integer  :: i,j

    print *, 'checking poles, ',str
    
    do i = 1,size(res0)
       do j = -2,-1
          tmp = (res1(-j,i)-i1(j,i)*res0(i))/res0(i)
          print *, i,j,res1(-j,i),i1(j,i)*res0(i),tmp
       enddo
       print *, ''
    enddo

    print *, 'finite part divided by tree-level'
    do i = 1,size(res0)
       tmp = (res1(-0,i)-i1(0,i)*res0(i))/res0(i)
       tmp = tmp - pisq/12._dp * res1(2,i)/res0(i) !-- fix normalization
       print *, i,0,res1(-0,i),i1(0,i)*res0(i),tmp
    enddo
    print *, ''
    
  end subroutine check_poles

  subroutine check_poles_sing(str,i1,res0,res1)
    character(*), intent(in) :: str
    real(dp), intent(in)   :: i1(-2:),res0
    real(dp15), intent(in) :: res1(0:)
    real(dp) :: tmp
    integer  :: j

    print *, 'checing poles, ',str
    
    do j = -2,-1
       tmp = (res1(-j)-i1(j)*res0)/res0
       print *, j,res1(-j),i1(j)*res0,tmp
    enddo

    print *, 'finite part divided by tree-level'
    tmp = (res1(-0)-i1(0)*res0)/res0
    tmp = tmp - pisq/12._dp * res1(2)/res0 !-- fix normalization
    print *, 0,res1(-0),i1(0)*res0,tmp
    print *, ''
    
  end subroutine check_poles_sing

  !-- 0 -> q[i1] qb[i2] l[i3] lb[i4]
  !-- returns poles for dn,up,
  function i1qcd(i1,i2,Lij,Lij2) result(res)
    integer, intent(in)   :: i1,i2
    real(dp), intent(in)  :: Lij(:,:),Lij2(:,:)
    real(dp)              :: res(-2:0,2)

    res(-2,:) = -two*Cf
    res(-1,:) = Cf*(-three + two*Lij(i1,i2))
    res( 0,:) = Cf*(zeta2 + three*Lij(i1,i2) - Lij2(i1,i2))
    
  end function i1qcd

  !-- 0 -> a[i1] a[i2] l[i3] lb[i4]
  !-- returns poles identical to QCD with Cf -> Q_lep2
  function i1ewk_aallb(i1,i2,Lij,Lij2) result(res)
    integer, intent(in)   :: i1,i2
    real(dp), intent(in)  :: Lij(:,:),Lij2(:,:)
    real(dp)              :: res(-2:0)

    res(-2) = -two*Q_lep2
    res(-1) = Q_lep2*(-three + two*Lij(i1,i2))
    res( 0) = Q_lep2*(zeta2 + three*Lij(i1,i2) - Lij2(i1,i2))

  end function i1ewk_aallb

  !-- 0 -> q[i1] qb[i2] l[i3] lb[i4]
  !-- returns poles for dn,up,bt
  function i1ewk_qqbllb(i1,i2,i3,i4,Lij,Lij2) result(res)
    integer, intent(in)   :: i1,i2,i3,i4
    real(dp), intent(in)  :: Lij(:,:),Lij2(:,:)
    real(dp)              :: res(-2:0,3)
    real(dp), parameter :: Qq(3) = [Qdn,Qup,Qdn]
    real(dp) :: bit1(3),bit2(3),bit3(3)

    bit1 = Q_lep**2 + Qq**2
    
    bit2 = Qq**2*Lij(i1,i2) + Q_lep**2*Lij(i3,i4) &
         -Qq*Q_lep*(Lij(i1,i3)+Lij(i2,i4)-Lij(i2,i3)-Lij(i1,i4))
    
    bit3 = Qq**2*Lij2(i1,i2) + Q_lep**2*Lij2(i3,i4) &
         -Qq*Q_lep*(Lij2(i1,i3)+Lij2(i2,i4)-Lij2(i2,i3)-Lij2(i1,i4))
   
    res(-2,:) = -2*bit1
    res(-1,:) = -3*bit1 + 2*bit2
    res( 0,:) = zeta2*bit1 + 3*bit2 - bit3
   
  end function i1ewk_qqbllb

#endif  

  !-----------------------------------------------------------------
  !--- generic routines
  !-----------------------------------------------------------------
  
  !-- add pisq to Lij2 where needed
  subroutine get_Lij(p,mu,is,Lij,Lij2)
    real(dp), intent(in)  :: p(:,:),mu
    integer,  intent(in)  :: is(:)
    real(dp), intent(out) :: Lij(size(p,2),size(p,2))
    real(dp), intent(out) :: Lij2(size(p,2),size(p,2))
    integer  :: i,j,n,i1,i2
    real(dp) :: musq

    musq = mu**2
    
    n = size(is)
    do i = 1,n-1
       do j = i+1,n
          i1 = is(i); i2 = is(j)
          Lij(i1,i2) = log(two*scr(p(:,i1),p(:,i2))/musq)
          if ( (i1.le.2.and.i2.le.2) .or. (i1.gt.2.and.i2.gt.2) ) then
             Lij2(i1,i2) = Lij(i1,i2)**2 - pisq
          else
             Lij2(i1,i2) = Lij(i1,i2)**2
          endif
          Lij(i2,i1)  = Lij(i1,i2)
          Lij2(i2,i1) = Lij2(i1,i2)
       enddo
    enddo
  end subroutine get_Lij

  !-- returns finite part for dn,up,bt
  !-- to be subtracted to ol output, as it is
  function get_fin_ewk_qqbllb(i1,i2,i3,i4,Lij,Lij2) result(res)
    integer, intent(in)   :: i1,i2,i3,i4
    real(dp), intent(in)  :: Lij(:,:),Lij2(:,:)
    real(dp)              :: res(3)
    real(dp), parameter :: Qq(3) = [Qdn,Qup,Qdn]
    real(dp) :: bit2(3),bit3(3)

    bit2 = Qq**2*Lij(i1,i2) + Q_lep**2*Lij(i3,i4) &
         -Qq*Q_lep*(Lij(i1,i3)+Lij(i2,i4)-Lij(i2,i3)-Lij(i1,i4))
    
    bit3 = Qq**2*Lij2(i1,i2) + Q_lep**2*Lij2(i3,i4) &
         -Qq*Q_lep*(Lij2(i1,i3)+Lij2(i2,i4)-Lij2(i2,i3)-Lij2(i1,i4))
    
    res(:) = 3*bit2 - bit3
   
  end function get_fin_ewk_qqbllb
  


  function get_fin_ewk_qqbllb_gen(i1,i2,i3,i4,Lij,Lij2,Qq,Qqbp,Ql,Qlb) result(res) 
  implicit none

  integer, intent(in)   :: i1,i2,i3,i4
  real(dp), intent(in)  :: Lij(:,:), Lij2(:,:)
  real(dp) :: res

  real(dp) :: Qq, Qqbp, Ql, Qlb
  real(dp) :: logbit, log2bit

  !---------------------------------------------------
  ! Compute logarithmic pieces
  !---------------------------------------------------
  logbit  = ewk_log_bit (Qq, Qqbp, Ql, Qlb, i1, i2, i3, i4, Lij)
  log2bit = ewk_log2_bit(Qq, Qqbp, Ql, Qlb, i1, i2, i3, i4, Lij2)

  res = 3.0_dp*logbit - log2bit

contains

  !---------------------------------------------------
  ! Single-log contribution
  !---------------------------------------------------
  function ewk_log_bit(Qq, Qqbp, Ql, Qlb, i1, i2, i3, i4, Lij) result(val)
    implicit none
    real(dp), intent(in) :: Qq, Qqbp, Ql, Qlb
    integer, intent(in)  :: i1, i2, i3, i4
    real(dp), intent(in) :: Lij(:,:)
    real(dp)             :: val

    val = - Qq *Qqbp*Lij(i1,i2) + Ql *Qq  *Lij(i1,i3) &
          + Qlb*Qq  *Lij(i1,i4) + Ql *Qqbp*Lij(i2,i3) &
          + Qlb*Qqbp*Lij(i2,i4) - Ql *Qlb *Lij(i3,i4)
  end function ewk_log_bit


  !---------------------------------------------------
  ! Double-log contribution
  !---------------------------------------------------
  function ewk_log2_bit(Qq, Qqbp, Ql, Qlb, i1, i2, i3, i4, Lij2) result(val)
    implicit none
    real(dp), intent(in) :: Qq, Qqbp, Ql, Qlb
    integer, intent(in)  :: i1, i2, i3, i4
    real(dp), intent(in) :: Lij2(:,:)
    real(dp)             :: val

    val = -  Qq *Qqbp*Lij2(i1,i2) + Ql *Qq  *Lij2(i1,i3) &
          + Qlb*Qq  *Lij2(i1,i4) + Ql *Qqbp*Lij2(i2,i3) &
          + Qlb*Qqbp*Lij2(i2,i4) - Ql *Qlb *Lij2(i3,i4)
  end function ewk_log2_bit

end function get_fin_ewk_qqbllb_gen


  !-- returns finite part for dn,up
  !-- to be subtracted to ol output, as it is
  function get_fin_qcd(i1,i2,Lij,Lij2) result(res)
    integer, intent(in)   :: i1,i2
    real(dp), intent(in)  :: Lij(:,:),Lij2(:,:)
    real(dp)              :: res(2)

    res(:) = Cf*(three*Lij(i1,i2) - Lij2(i1,i2))
    
  end function get_fin_qcd

  !-- returns finite part for dn,up
  function get_fin_ewk_aallb(i1,i2,Lij,Lij2) result(res)
    integer, intent(in)   :: i1,i2
    real(dp), intent(in)  :: Lij(:,:),Lij2(:,:)
    real(dp)              :: res

    res = 2 * gamma_a * Lij(i1,i2) - Lij2(i1,i2) + 3*Lij(i1,i2)

  end function get_fin_ewk_aallb

  !-----------------------------------------------------------------
  !--- Master amplitude for QCD corrections
  !-----------------------------------------------------------------

  !-- return finite part -> no mu dependence
  subroutine res_qcdloop_j(p,iconf,ave,res0,res1fin)
    real(dp), intent(in)  :: p(:,:),ave
    integer, intent(in)   :: iconf(:,:)
    real(dp), intent(out) :: res0(size(iconf,2),2),res1fin(size(iconf,2),2)
    integer     :: i,hq,ha,hl
    real(dp)    :: sprod(5,5)
    complex(dp) :: za(5,5),zb(5,5),coupl_is(1:2,-1:1,-1:1),coupl_fs(1:2,-1:1,-1:1)
    complex(dp) :: amp_is(1:2,-1:1,-1:1,-1:1),amp_fs(-1:1,-1:1,-1:1),amp0(2),amp1(2)
    integer :: ismin,ismax,fsmin,fsmax
    logical :: need

    res0 = zero
    res1fin = zero
    
    ismin = -1; ismax = -1
    fsmin = -1; fsmax = -1

    call spinoru(5,(/-p(:,1),-p(:,2),p(:,3),p(:,4),p(:,5)/),za,zb,sprod)

    do i = 1,size(iconf,2)

       !-- initial-state emission
       amp_is = master_qcdloop_amp_qgqb_llb(iconf(1,i),iconf(2,i),iconf(3,i),iconf(4,i),iconf(5,i),za,zb,sprod)
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

                amp0 = [Qdn,Qup]*amp_is(1,hq,ha,hl)*coupl_is(:,hq,hl) + Q_lep*amp_fs(hl,ha,hq)*coupl_fs(:,hq,hl)

                !-- finite part of FF: -8*Cf
                amp1 = [Qdn,Qup]*amp_is(2,hq,ha,hl)*coupl_is(:,hq,hl) + Q_lep*amp_fs(hl,ha,hq)*coupl_fs(:,hq,hl)*(-8._dp*Cf)

                res0(i,:) = res0(i,:) + abs(amp0)**2
                res1fin(i,:) = res1fin(i,:) + real(amp0*conjg(amp1),kind=dp) !-- no factor 2 because of 1L normalization

             enddo
          enddo
       enddo
       
    enddo

    res0 = res0 * 8._dp * xn * ave * eesq2 ! * ee**2
    res1fin = res1fin * 8._dp * xn * ave * eesq2 ! * ee**2 * asontwopi

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
    
  end subroutine res_qcdloop_j

  function master_qcdloop_amp_qgqb_llb(i1,i2,i3,i4,i5,za,zb,sprod) result(res)
    complex(dp) :: res(2,-1:1,-1:1,-1:1)
    integer, intent(in) :: i1,i2,i3,i4,i5
    complex(dp), intent(in) :: za(5,5),zb(5,5)
    real(dp), intent(in) :: sprod(5,5)

    !-- this could be done better, but good enough for now
    res(1:2,+1,+1,-1) = -A5NLOfin_THV(i1,i2,i3,i4,i5,za,zb,sprod)
    res(1:2,+1,-1,-1) = -A5NLOfin_THV(i3,i2,i1,i5,i4,zb,za,sprod)
    res(1:2,+1,+1,+1) = -A5NLOfin_THV(i1,i2,i3,i5,i4,za,zb,sprod)
    res(1:2,+1,-1,+1) = -A5NLOfin_THV(i3,i2,i1,i4,i5,zb,za,sprod)
    res(1:2,-1,+1,-1) =  A5NLOfin_THV(i3,i2,i1,i4,i5,za,zb,sprod)
    res(1:2,-1,-1,-1) =  A5NLOfin_THV(i1,i2,i3,i5,i4,zb,za,sprod)
    res(1:2,-1,+1,+1) =  A5NLOfin_THV(i3,i2,i1,i5,i4,za,zb,sprod)
    res(1:2,-1,-1,+1) =  A5NLOfin_THV(i1,i2,i3,i4,i5,zb,za,sprod)

  end function master_qcdloop_amp_qgqb_llb
  
  !-- 0 -> 1q+ 2g+ 3qb- 4l- 5lb+ with Ca -> 0
  function A5NLOfin_THV(j1,j2,j3,j4,j5,za,zb,sprod) result(res)
    complex(dp) :: res(2)
    integer, intent(in) :: j1,j2,j3,j4,j5
    complex(dp), intent(in) :: za(5,5),zb(5,5)
    real(dp), intent(in) :: sprod(5,5)

    res = A52fin_THV(j1,j3,j2,j4,j5,za,zb,sprod)
    res(2) = -two*Cf*res(2)
    
  end function A5NLOfin_THV

  !-- 0 -> 1q+ 2qb- 3g+ 4l- 5lb+, SLC
  function A52fin_THV(j1,j2,j3,j4,j5,za,zb,sprod) result(res)
    complex(dp) :: res(2)
    integer, intent(in) :: j1,j2,j3,j4,j5
    complex(dp), intent(in) :: za(5,5),zb(5,5)
    real(dp), intent(in) :: sprod(5,5)
    complex(dp) :: A5lom
    complex(dp) :: Fcc,Fsc,VccPlusVscFin_THV

    !-- -i * A5tree  
    A5lom=za(j2,j4)**2/(za(j2,j3)*za(j3,j1)*za(j4,j5))
    
    Fcc=-za(j2,j4)**2/(za(j2,j3)*za(j3,j1)*za(j4,j5)) &
         *Lsm1(-sprod(j1,j2),-sprod(j4,j5),-sprod(j1,j3), -sprod(j4,j5)) &
         +za(j2,j4)*(za(j1,j2)*za(j3,j4)-za(j1,j4)*za(j2,j3)) &
         /(za(j2,j3)*za(j1,j3)**2*za(j4,j5)) &
         *Lsm1(-sprod(j1,j2),-sprod(j4,j5),-sprod(j2,j3),-sprod(j4,j5)) &
         +two*zb(j1,j3)*za(j1,j4)*za(j2,j4)/(za(j1,j3)*za(j4,j5)) &
         *L0(-sprod(j2,j3),-sprod(j4,j5))/sprod(j4,j5)
    
    Fsc=za(j1,j4)**2*za(j2,j3)/(za(j1,j3)**3*za(j4,j5)) &
         *Lsm1(-sprod(j1,j2),-sprod(j4,j5),-sprod(j2,j3),-sprod(j4,j5)) &
         -half*(za(j4,j1)*zb(j1,j3))**2*za(j2,j3)/(za(j1,j3)*za(j4,j5)) &
         *L1(-sprod(j4,j5),-sprod(j2,j3))/sprod(j2,j3)**2 &
         +za(j1,j4)**2*za(j2,j3)*zb(j3,j1)/(za(j1,j3)**2*za(j4,j5)) &
         *L0(-sprod(j4,j5),-sprod(j2,j3))/sprod(j2,j3) &
         -za(j2,j1)*zb(j1,j3)*za(j4,j3)*zb(j3,j5)/za(j1,j3) &
         *L1(-sprod(j4,j5),-sprod(j1,j2))/sprod(j1,j2)**2 &
         -za(j2,j1)*zb(j1,j3)*za(j3,j4)*za(j1,j4)/(za(j1,j3)**2*za(j4,j5)) &
         *L0(-sprod(j4,j5),-sprod(j1,j2))/sprod(j1,j2) &
         -half*zb(j3,j5)*(zb(j1,j3)*zb(j2,j5)+zb(j2,j3)*zb(j1,j5)) &
         /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*zb(j4,j5))

    !-- last -1/2: FDH -> THV
    VccPlusVscFin_THV = -3.5_dp + 1.5_dp*lnrat(-sprod(j4,j5),-sprod(j1,j2)) - half

    res(1) = A5lom
    res(2) = -(VccPlusVscFin_THV*A5lom+Fcc+Fsc) !-- extra minus sign to have same phase as us

    return

  contains
    
    function klog(x)
      real(dp), intent(in)  :: x
      complex(dp) :: klog
      
      if (x.gt.zero) then 
         klog = log(x) 
      else
         klog = log(abs(x)) - ci*pi
      endif
      
    end function klog
    
    function Lnrat(x,y) 
      real(dp), intent(in) :: x,y
      complex(dp) :: Lnrat
      
      Lnrat = klog(x) - klog(y)
      
    end function Lnrat
    
    function L0(r1,r2) 
      implicit none
      real(dp), intent(in) :: r1,r2
      complex(dp) :: L0
      
      L0 = Lnrat(r1,r2)/(one-r1/r2)
      
    end function L0
    
    function L1(r1,r2) 
      implicit none
      real(dp), intent(in) :: r1,r2
      complex(dp) :: L1
      
      L1 = (L0(r1,r2)+one)/(one-r1/r2)
      
    end function L1
    
    function Lsm1(r1,r2,r3,r4) 
      implicit none
      real(dp), intent(in) :: r1,r2,r3,r4
      complex(dp) :: Lsm1
      
      Lsm1 = kdilog(r1,r2)+ kdilog(r3,r4) & 
           +Lnrat(r1,r2)*Lnrat(r3,r4)-pisq/6.0_dp
      
    end function Lsm1
    
    function kdilog(r1,r2)  ! this is really just dilog(1-r1/r2) but with 
      implicit none           ! particular imaginary part 
      real(dp), intent(in) :: r1, r2
      real(dp) :: x 
      complex(dp) :: kdilog 
      
      x = r1/r2
      
      if (x.gt.zero) then 
         kdilog = dilog2(one-x)
      else
         kdilog = pisq/6.0_dp - dilog2(x) & 
              - log(one-x)*Lnrat(r1,r2) 
      endif
      
    end function kdilog
    
  end function A52fin_THV

  
  !-----------------------------------------------------------------
  !--- OpenLoops checks for QCD amplitudes
  !-----------------------------------------------------------------

  !res1loopfin is a matrix in flavour space (-5:7)
  subroutine ol_res_qcdloop_qqb_gen(p,res0,res1loopfin)
    real(dp), intent(in)  :: p(4,4)
    real(dp), intent(out) :: res0(2,2), res1loopfin(-5:7,-5:7)
    real(dp)   :: res1fin(2,2)
    real(dp15) :: p_ol(4,4)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: ord(2,2),i,j
    real(dp) :: Lij(4,4),Lij2(4,4),myfin(1:2)

#if (_Vcharge == 0)
    
    p_ol(:,3:4) = p(:,3:4)
    ord(:,1) = [1,2] !-- q qb
    ord(:,2) = [2,1] !-- qb q

    do i = 1,2
       p_ol(:,1) = p(:,ord(1,i)); p_ol(:,2) = p(:,ord(2,i))
       do j = 1,2
          call evaluate_loop(OL_id(j),p_ol,res0_ol(i,j),res1(:,i,j),acc(i,j))
       enddo
    enddo

#else

    p_ol(:,1:4) = p(:,1:4)
    !q qb: 2 -1 -> 12 -11 OR 1 -2 -> 11 -12
    call evaluate_loop(OL_id(1),p_ol,res0_ol(1,1),res1(:,1,1),acc(1,1))
    !qb q: -1 2 -> 12 -11 OR -2 1 -> 11 -12
    call evaluate_loop(OL_id(2),p_ol,res0_ol(2,1),res1(:,2,1),acc(2,1))

#endif
    
    res0 = res0_ol
    !-- remove our couplings: aem/twopi
    res1 = res1/(as_ol/twopi)

    !-- go to ``fin''
    call get_Lij(p,mu_ol,[1,2],Lij,Lij2)
    myfin = get_fin_qcd(1,2,Lij,Lij2)

    res1fin(1,:) = res1(0,1,:) - myfin(:)*res0(1,:)
    res1fin(2,:) = res1(0,2,:) - myfin(:)*res0(2,:)

    !-- filling the res1loopfin in flavour space
    res1loopfin=0
#if (_Vcharge == 0) 
    res1loopfin(1,-1) = res1fin(1,1)
    res1loopfin(2,-2) = res1fin(1,2)
    res1loopfin(3,-3) = res1fin(1,1)
    res1loopfin(4,-4) = res1fin(1,2)
    res1loopfin(5,-5) = res1fin(1,1)
    !
    res1loopfin(-1,1) = res1fin(2,1)
    res1loopfin(-2,2) = res1fin(2,2)
    res1loopfin(-3,3) = res1fin(2,1)
    res1loopfin(-4,4) = res1fin(2,2)
    res1loopfin(-5,5) = res1fin(2,1) !no diff in the nf part
#elif (_Vcharge == 1)
    res1loopfin(2,-1) = res1fin(1,1)
    res1loopfin(-1,2) = res1fin(2,1)
    res1loopfin(4,-3) = res1loopfin(2,-1)
    res1loopfin(-3,4) = res1loopfin(-1,2)
#elif (_Vcharge == -1)
    res1loopfin(1,-2) = res1fin(1,1)
    res1loopfin(-2,1) = res1fin(2,1)
    res1loopfin(3,-4) = res1loopfin(1,-2)
    res1loopfin(-4,3) = res1loopfin(-2,1)
#endif    


  end subroutine ol_res_qcdloop_qqb_gen

  
  !-- res(1,:) = q qb -> e- e [dn,up]
  !-- res(2,:) = qb q -> e- e [dn,up]
  subroutine ol_res_qcdloop_qqb(p,res0,res1fin)
    real(dp), intent(in)  :: p(4,4)
    real(dp), intent(out) :: res0(2,2),res1fin(2,2)
    real(dp15) :: p_ol(4,4)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: ord(2,2),i,j
    real(dp) :: Lij(4,4),Lij2(4,4),myfin(1:2)

    p_ol(:,3:4) = p(:,3:4)
    
    ord(:,1) = [1,2] !-- q qb
    ord(:,2) = [2,1] !-- qb q

    do i = 1,2
       p_ol(:,1) = p(:,ord(1,i)); p_ol(:,2) = p(:,ord(2,i))
       do j = 1,2
          call evaluate_loop(OL_id(j),p_ol,res0_ol(i,j),res1(:,i,j),acc(i,j))
       enddo
    enddo
    res0 = res0_ol

    !-- remove our couplings: aem/twopi
    res1 = res1/(as_ol/twopi)

    !-- go to ``fin''
    call get_Lij(p,mu_ol,[1,2],Lij,Lij2)
    myfin = get_fin_qcd(1,2,Lij,Lij2)
    res1fin(1,:) = res1(0,1,:) - myfin(:)*res0(1,:)
    res1fin(2,:) = res1(0,2,:) - myfin(:)*res0(2,:)

#if(_withchecks == 1)
    if (checkpoles) then
       call check_poles('q qb->e- e+ qcd',i1qcd(1,2,Lij,Lij2),res0(1,:),res1(:,1,:))
       call check_poles('qb q->e- e+ qcd',i1qcd(1,2,Lij,Lij2),res0(2,:),res1(:,2,:))
    endif
#endif

  end subroutine ol_res_qcdloop_qqb

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! gen version of QCDloop correction to QED emission in qqb
  subroutine ol_res_qcdloop_a_qqb_gen(p,res0,res1loopfin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res0(2,2), res1loopfin(-5:7,-5:7)
    real(dp)   :: res1fin(2,2)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: ord(2,2),i,j
    real(dp) :: Lij(5,5),Lij2(5,5),myfin(1:2)

#if (_Vcharge ==0 )
    
    p_ol(:,3:5) = p(:,3:5)
    
    ord(:,1) = [1,2] !-- q qb
    ord(:,2) = [2,1] !-- qb q

    do i = 1,2
       p_ol(:,1) = p(:,ord(1,i)); p_ol(:,2) = p(:,ord(2,i))
       do j = 1,2
          call evaluate_loop(OL_id(j+2),p_ol,res0_ol(i,j),res1(:,i,j),acc(i,j))
       enddo
    enddo

#else

    p_ol(:,1:5) = p(:,1:5)
    !q qb: 2 -1 -> 12 -11 22 OR 1 -2 -> 11 -12 22
    call evaluate_loop(OL_id(3),p_ol,res0_ol(1,1),res1(:,1,1),acc(1,1))
    !qb q: -1 2 -> 12 -11 22 OR -2 1 -> 11 -12 22
    call evaluate_loop(OL_id(4),p_ol,res0_ol(2,1),res1(:,2,1),acc(2,1))

#endif


    res0 = res0_ol

    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(two*aem_ol*as_ol)
    res0 = res0/(aem_ol*four*pi)

    !-- go to ``fin''
    call get_Lij(p,mu_ol,[1,2],Lij,Lij2)
    myfin = get_fin_qcd(1,2,Lij,Lij2)
    res1fin(1,:) = res1(0,1,:) - myfin(:)*res0(1,:)
    res1fin(2,:) = res1(0,2,:) - myfin(:)*res0(2,:)

    !-- filling the res1loopfin in flavour space
    res1loopfin=0
#if (_Vcharge == 0) 
    res1loopfin(1,-1) = res1fin(1,1)
    res1loopfin(2,-2) = res1fin(1,2)
    res1loopfin(3,-3) = res1fin(1,1)
    res1loopfin(4,-4) = res1fin(1,2)
    res1loopfin(5,-5) = res1fin(1,1)
    !
    res1loopfin(-1,1) = res1fin(2,1)
    res1loopfin(-2,2) = res1fin(2,2)
    res1loopfin(-3,3) = res1fin(2,1)
    res1loopfin(-4,4) = res1fin(2,2)
    res1loopfin(-5,5) = res1fin(2,1)
#elif (_Vcharge == 1)
    res1loopfin(2,-1) = res1fin(1,1)
    res1loopfin(-1,2) = res1fin(2,1)
    res1loopfin(4,-3) = res1loopfin(2,-1)
    res1loopfin(-3,4) = res1loopfin(-1,2)
#elif (_Vcharge == -1)
    res1loopfin(1,-2) = res1fin(1,1)
    res1loopfin(-2,1) = res1fin(2,1)
    res1loopfin(3,-4) = res1loopfin(1,-2)
    res1loopfin(-4,3) = res1loopfin(-2,1)
#endif    


  end subroutine ol_res_qcdloop_a_qqb_gen
  
  !-- res(1,:) = [d db, db d] -> e- e+ a
  !-- res(2,:) = [u ub, ub u] -> e- e+ a
  subroutine ol_res_qcdloop_a_qqb(p,res0,res1fin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res0(2,2),res1fin(2,2)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: ord(2,2),i,j
    real(dp) :: Lij(5,5),Lij2(5,5),myfin(1:2)

    p_ol(:,3:5) = p(:,3:5)
    
    ord(:,1) = [1,2] !-- q qb
    ord(:,2) = [2,1] !-- qb q

    do i = 1,2
       p_ol(:,1) = p(:,ord(1,i)); p_ol(:,2) = p(:,ord(2,i))
       do j = 1,2
          call evaluate_loop(OL_id(j+2),p_ol,res0_ol(i,j),res1(:,i,j),acc(i,j))
       enddo
    enddo
    res0 = res0_ol

    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(two*aem_ol*as_ol)
    res0 = res0/(aem_ol*four*pi)

    !-- go to ``fin''
    call get_Lij(p,mu_ol,[1,2],Lij,Lij2)
    myfin = get_fin_qcd(1,2,Lij,Lij2)
    res1fin(1,:) = res1(0,1,:) - myfin(:)*res0(1,:)
    res1fin(2,:) = res1(0,2,:) - myfin(:)*res0(2,:)

#if(_withchecks == 1)
    if (checkpoles) then
       call check_poles('q qb->e- e+ γ qcd',i1qcd(1,2,Lij,Lij2),res0(1,:),res1(:,1,:))
       call check_poles('qb q->e- e+ γ qcd',i1qcd(1,2,Lij,Lij2),res0(2,:),res1(:,2,:))
    endif
#endif

  end subroutine ol_res_qcdloop_a_qqb

  !-- gen version of aq channel with QCDloop
  subroutine ol_res_qcdloop_a_aq_gen(p,res0,res1loopfin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res0(2,2), res1loopfin(-5:7,-5:7)
    real(dp)   :: res1fin(2,2)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: j
    real(dp) :: Lij(5,5),Lij2(5,5),myfin(1:2)

#if (_Vcharge == 0 )
    
    p_ol(:,1:5) = p(:,1:5)
    ! a qbar -> underline process q qb
    do j = 1,2
       call evaluate_loop(OL_id(j+2),p_ol,res0_ol(1,j),res1(:,1,j),acc(1,j))
    enddo
    ! a q
    do j = 1,2
       call evaluate_loop(OL_id(j+4),p_ol,res0_ol(2,j),res1(:,2,j),acc(2,j))
    enddo

#else

    p_ol(:,1:5) = p(:,1:5)
    !a qb: 22 -1 -> 12 -11 -2 OR 22 -2 -> 11 -12 -1
    !I fill the res1(:,1,:) since the underline process is q qb
    call evaluate_loop(OL_id(3),p_ol,res0_ol(1,1),res1(:,1,1),acc(1,1))
    !a q: 22 2 -> 12 -11 1 OR 22 1 -> 11 -12 2
    !I fill the res1(:,2,:) since the underline process is qb q
    call evaluate_loop(OL_id(4),p_ol,res0_ol(2,1),res1(:,2,1),acc(2,1))    

#endif
    
    res0 = res0_ol

    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(two*aem_ol*as_ol)
    res0 = res0/(aem_ol*four*pi)
    
    !-- go to ``fin''
    call get_Lij(p,mu_ol,[2,5],Lij,Lij2)
    myfin = get_fin_qcd(2,5,Lij,Lij2)
    res1fin(1,:) = res1(0,1,:) - myfin(:)*res0(1,:)
    res1fin(2,:) = res1(0,2,:) - myfin(:)*res0(2,:)

    !-- filling the res1loopfin in flavour space
    res1loopfin=0
#if (_Vcharge == 0) 
    res1loopfin(7,-1) = res1fin(1,1)
    res1loopfin(7,-2) = res1fin(1,2)
    res1loopfin(7,-3) = res1fin(1,1)
    res1loopfin(7,-4) = res1fin(1,2)
    res1loopfin(7,-5) = res1fin(1,1)
    !
    res1loopfin(7,1) = res1fin(2,1)
    res1loopfin(7,2) = res1fin(2,2)
    res1loopfin(7,3) = res1fin(2,1)
    res1loopfin(7,4) = res1fin(2,2)
    res1loopfin(7,5) = res1fin(2,1)
#elif (_Vcharge == 1)
    res1loopfin(7,-1) = res1fin(1,1)
    res1loopfin(7,2) = res1fin(2,1)
    res1loopfin(7,-3) = res1loopfin(7,-1)
    res1loopfin(7,4) = res1loopfin(7,2)
#elif (_Vcharge == -1)
    res1loopfin(7,-2) = res1fin(1,1)
    res1loopfin(7,1) = res1fin(2,1)
    res1loopfin(7,-4) = res1loopfin(7,-2)
    res1loopfin(7,3) = res1loopfin(7,1)
#endif    
    
  end subroutine ol_res_qcdloop_a_aq_gen

  
  !-- res(1,:) = a qb -> e- e+ qb [dn,up]
  !-- res(2,:) = a q  -> e- e+ q  [dn,up]
  subroutine ol_res_qcdloop_a_aq(p,res0,res1fin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res0(2,2),res1fin(2,2)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: j
    real(dp) :: Lij(5,5),Lij2(5,5),myfin(1:2)

    p_ol(:,1:5) = p(:,1:5)
    
    do j = 1,2
       call evaluate_loop(OL_id(j+2),p_ol,res0_ol(1,j),res1(:,1,j),acc(1,j))
    enddo

    do j = 1,2
       call evaluate_loop(OL_id(j+4),p_ol,res0_ol(2,j),res1(:,2,j),acc(2,j))
    enddo
    res0 = res0_ol

    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(two*aem_ol*as_ol)
    res0 = res0/(aem_ol*four*pi)

    !-- go to ``fin''
    call get_Lij(p,mu_ol,[2,5],Lij,Lij2)
    myfin = get_fin_qcd(2,5,Lij,Lij2)
    res1fin(1,:) = res1(0,1,:) - myfin(:)*res0(1,:)
    res1fin(2,:) = res1(0,2,:) - myfin(:)*res0(2,:)

#if(_withchecks == 1)
    if (checkpoles) then
       call check_poles('γ qb->e- e+ qb qcd',i1qcd(2,5,Lij,Lij2),res0(1,:),res1(:,1,:))
       call check_poles('γ  q->e- e+ q  qcd',i1qcd(2,5,Lij,Lij2),res0(2,:),res1(:,2,:))
    endif
#endif

  end subroutine ol_res_qcdloop_a_aq

  !-- gen version of a_qa for QCDloop with photon radiation
  subroutine ol_res_qcdloop_a_qa_gen(p,res0,res1loopfin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res0(2,2), res1loopfin(-5:7,-5:7)
    real(dp)   :: res1fin(2,2)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: j
    real(dp) :: Lij(5,5),Lij2(5,5),myfin(1:2)

#if (_Vcharge == 0)
    
    p_ol(:,1:5) = p(:,1:5)
    ! q a -> underline event is q qb
    do j = 1,2
       call evaluate_loop(OL_id(j+2),p_ol,res0_ol(1,j),res1(:,1,j),acc(1,j))
    enddo
    ! qb a -> underline event is qb q
    do j = 1,2
       call evaluate_loop(OL_id(j+4),p_ol,res0_ol(2,j),res1(:,2,j),acc(2,j))
    enddo

#else

    p_ol(:,1:5) = p(:,1:5)
    ! 2 22 -> 12 -11 1 or 1 22 -> 11 -12 2
    call evaluate_loop(OL_id(3),p_ol,res0_ol(1,1),res1(:,1,1),acc(1,1))
    ! -1 22 -> 12 -11 -2 or -2 22 -> 11 -12 -1
    call evaluate_loop(OL_id(4),p_ol,res0_ol(2,1),res1(:,2,1),acc(2,1))    

#endif

    
    res0 = res0_ol
    
    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(two*aem_ol*as_ol)
    res0 = res0/(aem_ol*four*pi)

    !-- go to ``fin''
    call get_Lij(p,mu_ol,[1,5],Lij,Lij2)
    myfin = get_fin_qcd(1,5,Lij,Lij2)
    res1fin(1,:) = res1(0,1,:) - myfin(:)*res0(1,:)
    res1fin(2,:) = res1(0,2,:) - myfin(:)*res0(2,:)

    !-- filling the res1loopfin in flavour space
    res1loopfin=0
#if (_Vcharge == 0) 
    res1loopfin(1,7) = res1fin(1,1)
    res1loopfin(2,7) = res1fin(1,2)
    res1loopfin(3,7) = res1fin(1,1)
    res1loopfin(4,7) = res1fin(1,2)
    res1loopfin(5,7) = res1fin(1,1)
    !
    res1loopfin(-1,7) = res1fin(2,1)
    res1loopfin(-2,7) = res1fin(2,2)
    res1loopfin(-3,7) = res1fin(2,1)
    res1loopfin(-4,7) = res1fin(2,2)
    res1loopfin(-5,7) = res1fin(2,1)
#elif (_Vcharge == 1)
    res1loopfin(2,7) = res1fin(1,1)
    res1loopfin(-1,7) = res1fin(2,1)
    res1loopfin(4,7) = res1loopfin(2,7)
    res1loopfin(-3,7) = res1loopfin(-1,7)
#elif (_Vcharge == -1)
    res1loopfin(1,7) = res1fin(1,1)
    res1loopfin(-2,7) = res1fin(2,1)
    res1loopfin(-4,7) = res1loopfin(-2,7)
    res1loopfin(3,7) = res1loopfin(1,7)
#endif    
    

  end subroutine ol_res_qcdloop_a_qa_gen
  
  !-- res(1,:) = q a  -> e- e+ q  [dn,up]
  !-- res(2,:) = qb a -> e- e+ qb [dn,up]
  subroutine ol_res_qcdloop_a_qa(p,res0,res1fin)
    real(dp), intent(in)  :: p(4,5)
    real(dp), intent(out) :: res0(2,2),res1fin(2,2)
    real(dp15) :: p_ol(4,5)
    real(dp15) :: res1(0:2,2,2),acc(2,2),res0_ol(2,2)
    integer  :: j
    real(dp) :: Lij(5,5),Lij2(5,5),myfin(1:2)

    p_ol(:,1:5) = p(:,1:5)
    
    do j = 1,2
       call evaluate_loop(OL_id(j+2),p_ol,res0_ol(1,j),res1(:,1,j),acc(1,j))
    enddo

    do j = 1,2
       call evaluate_loop(OL_id(j+4),p_ol,res0_ol(2,j),res1(:,2,j),acc(2,j))
    enddo
    res0 = res0_ol
    
    !-- remove our couplings: aem/twopi * (as*four*pi) = aem*as*two
    res1 = res1/(two*aem_ol*as_ol)
    res0 = res0/(aem_ol*four*pi)

    !-- go to ``fin''
    call get_Lij(p,mu_ol,[1,5],Lij,Lij2)
    myfin = get_fin_qcd(1,5,Lij,Lij2)
    res1fin(1,:) = res1(0,1,:) - myfin(:)*res0(1,:)
    res1fin(2,:) = res1(0,2,:) - myfin(:)*res0(2,:)

#if(_withchecks == 1)
    if (checkpoles) then
       call check_poles('q  γ->e- e+ q  qcd',i1qcd(1,5,Lij,Lij2),res0(1,:),res1(:,1,:))
       call check_poles('qb γ->e- e+ qb qcd',i1qcd(1,5,Lij,Lij2),res0(2,:),res1(:,2,:))
    endif
#endif

  end subroutine ol_res_qcdloop_a_qa
  
end module mod_amplitudes_loop_ppll

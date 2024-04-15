module mod_intsub_fc
  use mod_types
  use mod_consts_dp
  use mod_process
  use mod_parms
  use mod_proc_parms
  use mod_auxfunctions
  implicit none
  private

  public :: intsub_test_z,intsub_test_12 !-- nloqcd, ns
  public :: intsub_qqb_ga_onloewk_z2,intsub_qqb_ga_onloewk_1zb,intsub_qqb_ga_onloewk_12
  public :: intsub_qqb_ga_onloqcd_z2,intsub_qqb_ga_onloqcd_1zb,intsub_qqb_ga_onloqcd_12
  public :: intsub_qqb_ga_sub12_zzb,intsub_qqb_ga_sub12_z2,intsub_qqb_ga_sub12_1zb,intsub_qqb_ga_sub12_12

  !--

  public :: intsub_aq_onloqcd_z2
  public :: intsub_aq_onloewk_1zb,intsub_aq_onloewk_12
  public :: intsub_aq_sub12_zzb,intsub_aq_sub12_z2,intsub_aq_sub12_1zb

  !--

  public :: intsub_ag_onloqcd_z2,intsub_ag_onloqcd_z2_split
  public :: intsub_ag_sub12_zzb,intsub_ag_sub12_1zb

  !--

  public :: intsub_gq_onloqcd_1zb,intsub_gq_onloqcd_12
  
contains

  !-- Lmu2 = log(spart/mu2)

  !-- ONLOEWK below
  function intsub_qqb_ga_onloewk_z2(z,Lmu2,w11,kininv,iflag) result(res)
    real(dp) :: z,Lmu2(:),w11,kininv(7)
    real(dp) :: res(size(Lmu2))
    integer  :: iflag
    real(dp) :: eta51L,l51

    eta51L = kininv(6)
    l51    = log(eta51L)
    
    res = Cf*(PqqTNLO(z,Lmu2,iflag) + w11*l51*PqqAPR(z,iflag))
    
  end function intsub_qqb_ga_onloewk_z2

  function intsub_qqb_ga_onloewk_1zb(zb,Lmu2,w22,kininv,iflag) result(res)
    real(dp) :: zb,Lmu2(:),w22,kininv(7)
    real(dp) :: res(size(Lmu2))
    integer  :: iflag
    real(dp) :: eta52L,l52

    eta52L = kininv(7)
    l52    = log(eta52L)

    res = Cf*(PqqTNLO(zb,Lmu2,iflag) + w22*l52*PqqAPR(zb,iflag))
    
  end function intsub_qqb_ga_onloewk_1zb

  function intsub_qqb_ga_onloewk_12(Lmu2) result(res)
    real(dp) :: Lmu2(:)
    real(dp) :: res(size(Lmu2))

    res = Cf*(two/three*pisq + 3*Lmu2)

  end function intsub_qqb_ga_onloewk_12
  
  !-- ONLOQCD below
  function intsub_qqb_ga_onloqcd_z2(z,Lmu2,w11,kininv,iflag) result(res)
    real(dp) :: z,Lmu2(:),w11,kininv(7)
    real(dp) :: res(2,2,size(Lmu2))
    integer  :: iflag
    real(dp) :: eta51L,l51,restmp(size(Lmu2))

    eta51L = kininv(6)
    l51    = log(eta51L)
    
    restmp = PqqTNLO(z,Lmu2,iflag) + w11*l51*PqqAPR(z,iflag)

    res(1,1,:) = Qdn2*restmp(:)
    res(1,2,:) = Qup2*restmp(:)
    res(2,1,:) = Qdn2*restmp(:)
    res(2,2,:) = Qup2*restmp(:)
    
  end function intsub_qqb_ga_onloqcd_z2
  
  function intsub_qqb_ga_onloqcd_1zb(zb,Lmu2,w22,kininv,iflag) result(res)
    real(dp) :: zb,Lmu2(:),w22,kininv(7)
    real(dp) :: res(2,2,size(Lmu2))
    integer  :: iflag
    real(dp) :: eta52L,l52,restmp(size(Lmu2))

    eta52L = kininv(7)
    l52    = log(eta52L)

    restmp = PqqTNLO(zb,Lmu2,iflag) + w22*l52*PqqAPR(zb,iflag)

    res(1,1,:) = Qdn2*restmp(:)
    res(1,2,:) = Qup2*restmp(:)
    res(2,1,:) = Qdn2*restmp(:)
    res(2,2,:) = Qup2*restmp(:)
    
  end function intsub_qqb_ga_onloqcd_1zb

  function intsub_qqb_ga_onloqcd_12(Lmu2,kininv2) result(res)
    real(dp) :: Lmu2(:),kininv2(8)
    real(dp) :: res(2,2,size(Lmu2))
    real(dp) :: restmp_Qq2(size(Lmu2)),restmp_QqQl,restmp_Ql2,tmp(1)
    real(dp) :: spart,ec,e3,e4,eta13,eta14,eta23,eta24,eta34

    spart = kininv2(1)
    ec = half*sqrt(spart)
    e3 = kininv2(2)
    e4 = kininv2(3)
    eta13 = kininv2(4)
    eta14 = kininv2(5)
    eta23 = kininv2(6)
    eta24 = kininv2(7)
    eta34 = kininv2(8)

    restmp_Qq2  = two/three*pisq + 3*Lmu2
    
    tmp = Geq(e3,e4,[ec**2],eta13,eta14,eta23,eta24)
    restmp_QqQl = tmp(1)
    
    restmp_Ql2  = Ge2(e3,e4,ec,eta34)

    !--

    res(1,1,:) = (Qdn2*restmp_Qq2(:)+Qdn*Q_lep*restmp_QqQl+Q_lep2*restmp_Ql2)
    res(1,2,:) = (Qup2*restmp_Qq2(:)+Qup*Q_lep*restmp_QqQl+Q_lep2*restmp_Ql2)
    res(2,1,:) = (Qdn2*restmp_Qq2(:)-Qdn*Q_lep*restmp_QqQl+Q_lep2*restmp_Ql2)
    res(2,2,:) = (Qup2*restmp_Qq2(:)-Qup*Q_lep*restmp_QqQl+Q_lep2*restmp_Ql2)

  end function intsub_qqb_ga_onloqcd_12
  
  !--

  !-- Test below [NLO QCD]
  function intsub_test_z(z,Lmu2,iflag) result(res)
    real(dp) :: z,Lmu2(:)
    real(dp) :: res(2,2,size(Lmu2))
    real(dp) :: restmp(size(Lmu2))
    integer  :: iflag,i,j

    restmp = Cf*PqqTNLO(z,Lmu2,iflag)

    do i = 1,2
       do j = 1,2
          res(j,i,:) = restmp
       enddo
    enddo
    
  end function intsub_test_z

  function intsub_test_12(Lmu2) result(res)
    real(dp) :: Lmu2(:)
    real(dp) :: res(2,2,size(Lmu2))
    real(dp) :: restmp(size(Lmu2))
    integer  :: i,j

    restmp = two/three*Cf*pisq + three * Cf * Lmu2

    do i = 1,2
       do j = 1,2
          res(j,i,:) = restmp
       enddo
    enddo
    
  end function intsub_test_12

  !-- SUB12 below
  function intsub_qqb_ga_sub12_zzb(z,zb,Lmu2,iflagz,iflagzb) result(res)
    real(dp) :: z,zb,Lmu2(:)
    real(dp) :: res(2,2,size(Lmu2))
    real(dp) :: restmp(size(Lmu2))
    integer  :: iflagz,iflagzb

    restmp = PqqTNLO(z,Lmu2,iflagz)*PqqTNLO(zb,Lmu2,iflagzb)

    res(1,1,:) = two*Cf*Qdn2 * restmp !-- d db
    res(1,2,:) = two*Cf*Qup2 * restmp !-- u ub
    res(2,1,:) = two*Cf*Qdn2 * restmp !-- db d
    res(2,2,:) = two*Cf*Qup2 * restmp !-- ub u

  end function intsub_qqb_ga_sub12_zzb

  function intsub_qqb_ga_sub12_z2(z,mu2,Lmu2,iflag,kininv) result(res)
    real(dp) :: z,mu2(:),Lmu2(:),kininv(8)
    real(dp) :: res(2,2,size(Lmu2))
    integer  :: iflag
    real(dp) :: pqq(size(Lmu2)),restmp_Qq2(size(Lmu2)),restmp_QqQl(size(Lmu2)),restmp_Ql2(size(Lmu2)),lz
    real(dp) :: spart,ec,e3,e4,eta13,eta14,eta23,eta24,eta34

    if(iflag == 0) then
       lz = log(z)
    elseif (iflag == 1) then
       lz = zero
    else
       lz = -99
       print *, 'error'
       stop
    endif

    spart = kininv(1)
    ec = half*sqrt(spart)
    e3 = kininv(2)
    e4 = kininv(3)
    eta13 = kininv(4)
    eta14 = kininv(5)
    eta23 = kininv(6)
    eta24 = kininv(7)
    eta34 = kininv(8)
       
    pqq = PqqTNLO(z,Lmu2,iflag)
    
    restmp_Qq2  = PqqNNLO(z,Lmu2,iflag)

    restmp_QqQl = pqq*(&
         Geq(e3,e4,mu2,eta13,eta14,eta23,eta24) + &
         log(ec**2/mu2)*log(eta13*eta24/eta14/eta23) &
         -2*log(e3*eta13/e4/eta14)*lz)
    
    restmp_Ql2 = pqq * Ge2(e3,e4,ec,eta34)

    !--

    res(1,1,:) = Cf*(Qdn2*restmp_Qq2(:)+Qdn*Q_lep*restmp_QqQl(:)+Q_lep2*restmp_Ql2(:))
    res(1,2,:) = Cf*(Qup2*restmp_Qq2(:)+Qup*Q_lep*restmp_QqQl(:)+Q_lep2*restmp_Ql2(:))
    res(2,1,:) = Cf*(Qdn2*restmp_Qq2(:)-Qdn*Q_lep*restmp_QqQl(:)+Q_lep2*restmp_Ql2(:))
    res(2,2,:) = Cf*(Qup2*restmp_Qq2(:)-Qup*Q_lep*restmp_QqQl(:)+Q_lep2*restmp_Ql2(:))
        
  end function intsub_qqb_ga_sub12_z2

  function intsub_qqb_ga_sub12_1zb(zb,mu2,Lmu2,iflag,kininv) result(res)
    real(dp) :: zb,mu2(:),Lmu2(:),kininv(8)
    real(dp) :: res(2,2,size(Lmu2))
    integer  :: iflag
    real(dp) :: pqq(size(Lmu2)),restmp_Qq2(size(Lmu2)),restmp_QqQl(size(Lmu2)),restmp_Ql2(size(Lmu2)),lzb
    real(dp) :: spart,ec,e3,e4,eta13,eta14,eta23,eta24,eta34

    if(iflag == 0) then
       lzb = log(zb)
    elseif (iflag == 1) then
       lzb = zero
    else
       lzb = -99
       print *, 'error'
       stop
    endif

    spart = kininv(1)
    ec = half*sqrt(spart)
    e3 = kininv(2)
    e4 = kininv(3)
    eta13 = kininv(4)
    eta14 = kininv(5)
    eta23 = kininv(6)
    eta24 = kininv(7)
    eta34 = kininv(8)
       
    pqq = PqqTNLO(zb,Lmu2,iflag)
    
    restmp_Qq2  = PqqNNLO(zb,Lmu2,iflag)

    restmp_QqQl = pqq*(&
         Geq(e3,e4,mu2,eta13,eta14,eta23,eta24) + &
         log(ec**2/mu2)*log(eta13*eta24/eta14/eta23) &
         +2*log(e3*eta23/e4/eta24)*lzb)
    
    restmp_Ql2 = pqq * Ge2(e3,e4,ec,eta34)

    !--

    res(1,1,:) = Cf*(Qdn2*restmp_Qq2(:)+Qdn*Q_lep*restmp_QqQl(:)+Q_lep2*restmp_Ql2(:))
    res(1,2,:) = Cf*(Qup2*restmp_Qq2(:)+Qup*Q_lep*restmp_QqQl(:)+Q_lep2*restmp_Ql2(:))
    res(2,1,:) = Cf*(Qdn2*restmp_Qq2(:)-Qdn*Q_lep*restmp_QqQl(:)+Q_lep2*restmp_Ql2(:))
    res(2,2,:) = Cf*(Qup2*restmp_Qq2(:)-Qup*Q_lep*restmp_QqQl(:)+Q_lep2*restmp_Ql2(:))
    
  end function intsub_qqb_ga_sub12_1zb

  function intsub_qqb_ga_sub12_12(Lmu2,kininv) result(res)
    real(dp) :: Lmu2(:),kininv(8)
    real(dp) :: res(2,2,size(Lmu2))
    real(dp) :: restmp_Qq2(size(Lmu2)),restmp_QqQl(size(Lmu2)),restmp_Ql2(size(Lmu2))
    real(dp) :: spart,ec,e3,e4,eta13,eta14,eta23,eta24,eta34,tmp(1)

    spart = kininv(1)
    ec = half*sqrt(spart)
    e3 = kininv(2)
    e4 = kininv(3)
    eta13 = kininv(4)
    eta14 = kininv(5)
    eta23 = kininv(6)
    eta24 = kininv(7)
    eta34 = kininv(8)

    restmp_Qq2 = 16._dp/45._dp*pisq**2 + four*(pisq+8*zeta3)*Lmu2 &
         +(9-four/three*pisq)*Lmu2**2

    tmp = Geq(e3,e4,[ec**2],eta13,eta14,eta23,eta24)
    restmp_QqQl = tmp(1)*(two/three*pisq + 3*Lmu2)

    restmp_Ql2 = Ge2(e3,e4,ec,eta34)*(two/three*pisq+3*Lmu2)

    !--

    res(1,1,:) = Cf*(Qdn2*restmp_Qq2(:)+Qdn*Q_lep*restmp_QqQl(:)+Q_lep2*restmp_Ql2(:))
    res(1,2,:) = Cf*(Qup2*restmp_Qq2(:)+Qup*Q_lep*restmp_QqQl(:)+Q_lep2*restmp_Ql2(:))
    res(2,1,:) = Cf*(Qdn2*restmp_Qq2(:)-Qdn*Q_lep*restmp_QqQl(:)+Q_lep2*restmp_Ql2(:))
    res(2,2,:) = Cf*(Qup2*restmp_Qq2(:)-Qup*Q_lep*restmp_QqQl(:)+Q_lep2*restmp_Ql2(:))
    
  end function intsub_qqb_ga_sub12_12

  !-- aq channel below

  function intsub_aq_onloqcd_z2(z,Lmu2,w11,kininv) result(res)
    real(dp) :: z,Lmu2(:),w11,kininv(7)
    real(dp) :: res(2,2,size(Lmu2))
    real(dp) :: eta51L,l51,restmp(size(Lmu2))

    eta51L = kininv(6)
    l51    = log(half*eta51L)
    
    restmp = PqgTNLO(z,Lmu2) + w11*l51*PqgAPR(z)

    res(1,1,:) = xn*Qdn2*restmp(:)
    res(1,2,:) = xn*Qup2*restmp(:)
    res(2,1,:) = xn*Qdn2*restmp(:)
    res(2,2,:) = xn*Qup2*restmp(:)
    
  end function intsub_aq_onloqcd_z2
  
  function intsub_aq_onloewk_1zb(zb,Lmu2,w22,kininv,iflag) result(res)
    real(dp) :: zb,Lmu2(:),w22,kininv(7)
    real(dp) :: res(size(Lmu2))
    integer  :: iflag
    real(dp) :: eta52L,l52

    eta52L = kininv(7)
    l52    = log(eta52L*half)

    res = Cf*(PqqTNLO(zb,Lmu2,iflag) + w22*l52*PqqAPR(zb,iflag))
    
  end function intsub_aq_onloewk_1zb

  function intsub_aq_onloewk_12(mu2,w56_11,w56_22,proc) result(res)
    real(dp) :: mu2(:),w56_11,w56_22
    type(KinConfig), intent(in) :: proc
    real(dp) :: res(size(mu2))
    real(dp) :: spart,ec,e5
    real(dp) :: eta51L,eta52L,eta51,eta52
    real(dp) :: Lecoe5,Li2ometa

    spart = proc%Lim_kininv2(1)
    ec = half*sqrt(spart)
    e5     = proc%Lim_kininv(5)

    eta51L = proc%Lim_kininv(6)
    eta52L = proc%Lim_kininv(7)

    eta51 = proc%Lim_etaij(1,5)
    eta52 = proc%Lim_etaij(2,5)

    Lecoe5 = log(ec/e5)
    Li2ometa = real(dilog2(one-eta52),kind=dp)

    res = 6.5_dp + Lecoe5**2 - pisq + 2*Li2ometa + (3 + 2*Lecoe5)*log(4*eta52L) + &
         (3*Log(ec**2/mu2))/2._dp

    if (w56_11.ne.zero) then
       res = res + (-1.5_dp - 2*Lecoe5)*log(eta51L/(1 - eta51))*w56_11
    endif

    if (w56_22.ne.zero) then
       res = res + (-1.5_dp - 2*Lecoe5)*log(eta52L/(1 - eta52))*w56_22
    endif

    res = res*Cf

  end function intsub_aq_onloewk_12

  function intsub_aq_sub12_zzb(z,zb,Lmu2,iflagzb) result(res)
    real(dp) :: z,zb,Lmu2(:)
    real(dp) :: res(2,2,size(Lmu2))
    real(dp) :: restmp(size(Lmu2))
    integer  :: iflagzb

    restmp = PqgTNLO(z,Lmu2)*PqqTNLO(zb,Lmu2,iflagzb)

    res(1,1,:) = xn*Cf*Qdn2 * restmp !-- d db
    res(1,2,:) = xn*Cf*Qup2 * restmp !-- u ub
    res(2,1,:) = xn*Cf*Qdn2 * restmp !-- db d
    res(2,2,:) = xn*Cf*Qup2 * restmp !-- ub u

  end function intsub_aq_sub12_zzb

  function intsub_aq_sub12_z2(z,Lmu2) result(res)
    real(dp) :: z,Lmu2(:)
    real(dp) :: res(2,2,size(Lmu2))
    real(dp) :: restmp(size(Lmu2))
   
    restmp  = PqgNNLO(z,Lmu2)

    res(1,1,:) = Cf*xn*Qdn2*restmp(:)
    res(1,2,:) = Cf*xn*Qup2*restmp(:)
    res(2,1,:) = Cf*xn*Qdn2*restmp(:)
    res(2,2,:) = Cf*xn*Qup2*restmp(:)
        
  end function intsub_aq_sub12_z2

  function intsub_aq_sub12_1zb(zb,Lmu2) result(res)
    real(dp) :: zb,Lmu2(:)
    real(dp) :: res(2,2,size(Lmu2))
    real(dp) :: restmp(size(Lmu2))

    restmp  = PaqNNLO(zb,Lmu2)

    res(1,1,:) = Cf*Qdn2*restmp(:)
    res(1,2,:) = Cf*Qup2*restmp(:)
    res(2,1,:) = Cf*Qdn2*restmp(:)
    res(2,2,:) = Cf*Qup2*restmp(:)
    
  end function intsub_aq_sub12_1zb

  !-- ag channel bellow

  function intsub_ag_onloqcd_z2_split(z,Lmu2,w11,kininv) result(res)
    real(dp) :: z,Lmu2(:),w11,kininv(7)
    real(dp) :: res(2,2,size(Lmu2))
    real(dp) :: eta51L,l51
    real(dp) :: restmp(size(Lmu2))

    eta51L = kininv(6)
    l51    = log(eta51L)
    
    restmp = (PqgTNLO(z,Lmu2) + w11*l51*PqgAPR(z))

    res(1,1,:) = xn*Qdn2*restmp(:)
    res(1,2,:) = xn*Qup2*restmp(:)
    res(2,1,:) = xn*Qdn2*restmp(:)
    res(2,2,:) = xn*Qup2*restmp(:)

  end function intsub_ag_onloqcd_z2_split

  function intsub_ag_onloqcd_z2(z,Lmu2) result(res)
    real(dp) :: z,Lmu2(:)
    real(dp) :: res(2,2,size(Lmu2))
    real(dp) :: restmp(size(Lmu2))
    
    restmp = PqgTNLO(z,Lmu2)

    res(1,1,:) = xn*Qdn2*restmp(:)
    res(1,2,:) = xn*Qup2*restmp(:)
    res(2,1,:) = xn*Qdn2*restmp(:)
    res(2,2,:) = xn*Qup2*restmp(:)
    
  end function intsub_ag_onloqcd_z2

  function intsub_ag_sub12_zzb(z,zb,Lmu2) result(res)
    real(dp) :: z,zb,Lmu2(:)
    real(dp) :: res(2,2,size(Lmu2))
    real(dp) :: restmp(size(Lmu2))

    restmp = PqgTNLO(z,Lmu2)*PqgTNLO(zb,Lmu2)

    res(1,1,:) = tr*xn*Qdn2*restmp(:)
    res(1,2,:) = tr*xn*Qup2*restmp(:)
    res(2,1,:) = tr*xn*Qdn2*restmp(:)
    res(2,2,:) = tr*xn*Qup2*restmp(:)

  end function intsub_ag_sub12_zzb

  function intsub_ag_sub12_1zb(zb,Lmu2) result(res)
    real(dp) :: zb,Lmu2(:)
    real(dp) :: res(2,2,size(Lmu2))
    real(dp) :: restmp(size(Lmu2))

    restmp = Paq2NNLO(zb,Lmu2)

    res(1,1,:) = tr*Qdn2*restmp(:)
    res(1,2,:) = tr*Qup2*restmp(:)
    res(2,1,:) = tr*Qdn2*restmp(:)
    res(2,1,:) = tr*Qup2*restmp(:)

  end function intsub_ag_sub12_1zb
  
  !-- gq channel below
  function intsub_gq_onloqcd_12(mu2,proc) result(res)
    real(dp) :: mu2(:)
    type(KinConfig) :: proc
    real(dp) :: res(2,2,size(mu2))
    real(dp) :: restmp_Qq2(size(mu2)),restmp_QqQl,restmp_Ql2,tmp(1)
    real(dp) :: spart,ec,e3,e4,e5,eta15,eta25,eta35,eta45,eta23,eta24,eta34

    spart = proc%Lim_KinInv2(1)
    ec = half*sqrt(spart)
    e3 = proc%Lim_KinInv2(2)
    e4 = proc%Lim_KinInv2(3)
    e5 = proc%Lim_KinInv(5)

    eta15 = proc%Lim_KinInv(6)
    eta25 = proc%Lim_KinInv(7)

    eta35 = proc%Lim_etaij(3,5)
    eta45 = proc%Lim_etaij(4,5)
    
    eta23 = proc%Lim_KinInv2(6)
    eta24 = proc%Lim_KinInv2(7)
    eta34 = proc%Lim_KinInv2(8)

    restmp_Qq2  = Gq2b(e5,ec,mu2,eta15,eta25)
    
    tmp = Geq(e3,e4,[ec**2],eta35,eta45,eta23,eta24)
    restmp_QqQl = tmp(1)-2*log(e5/ec)*log(e3*eta35/e4/eta45)
    
    restmp_Ql2  = Ge2(e3,e4,ec,eta34)

    !--
    res(1,1,:) = (Qdn2*restmp_Qq2(:)+Qdn*Q_lep*restmp_QqQl+Q_lep2*restmp_Ql2)
    res(1,2,:) = (Qup2*restmp_Qq2(:)+Qup*Q_lep*restmp_QqQl+Q_lep2*restmp_Ql2)
    res(2,1,:) = (Qdn2*restmp_Qq2(:)-Qdn*Q_lep*restmp_QqQl+Q_lep2*restmp_Ql2)
    res(2,2,:) = (Qup2*restmp_Qq2(:)-Qup*Q_lep*restmp_QqQl+Q_lep2*restmp_Ql2)

  end function intsub_gq_onloqcd_12

  function intsub_gq_onloqcd_1zb(z,Lmu2,iflag) result(res)
    real(dp) :: z,Lmu2(:)
    real(dp) :: res(2,2,size(Lmu2)),restmp(size(Lmu2))
    integer  :: iflag

    restmp = PqqTNLO(z,Lmu2,iflag)

    res(1,1,:) = Qdn2*restmp
    res(1,2,:) = Qup2*restmp
    res(2,1,:) = Qdn2*restmp
    res(2,2,:) = Qup2*restmp

  end function intsub_gq_onloqcd_1zb
  
  !--- Functions below

  function PqqAPR(z,iflag) result(res)
    real(dp) :: res
    real(dp) :: z
    integer  :: iflag

    if (iflag == 0) then
       res = two/(one-z)-(one+z)
    elseif (iflag == 1) then
       res = -two/(one-z)
    else
       res = -99
       print *, 'error'
       stop
    endif

  end function PqqAPR

  function PqgAPR(z) result(res)
    real(dp) :: res
    real(dp) :: z

    res = one-two*z+two*z**2

  end function PqgAPR
  
  function PqqTNLO(z,Lmu2,iflag) result(res)
    real(dp) :: z,Lmu2(:)
    real(dp) :: res(size(Lmu2))
    integer  :: iflag
    real(dp) :: omz,lomz

    res = zero
    
    omz  = one-z
    lomz = log(omz)
    
    if (iflag == 0) then !-- reg + plus
       res = four*lomz/omz - 2*(one+z)*lomz + omz + Lmu2 * (2/omz-(one+z))
    elseif (iflag == 1) then !-- -plus
       res = -four*lomz/omz - 2*Lmu2/omz
    endif
    
  end function PqqTNLO
  
  function PqgTNLO(z,Lmu2) result(res)
    real(dp) :: z,Lmu2(:)
    real(dp) :: res(size(Lmu2))
    real(dp) :: lomz,pqg

    res = zero

    pqg  = one-two*z+two*z**2
    lomz = log(one-z)

    res = 2*(1 - z)*z + pqg*Lmu2 + 2*pqg*lomz
    
  end function PqgTNLO

  function PqqNNLO(z,Lmu2,iflag) result(res)
    real(dp) :: z,Lmu2(:)
    real(dp) :: res(size(Lmu2))
    integer  :: iflag
    real(dp) :: omz,lomz
    real(dp) :: lz,Li2omz,Li3omz,Li3mZ3,pqq

    res = zero
    
    omz  = one-z
    lomz = log(omz)
    
    if (iflag == 0) then !-- reg + plus

       lz = log(z)
       Li2omz = real(dilog2(one-z),kind=dp)
       Li3omz = trilog(one-z)
       Li3mZ3 = trilog(z)-zeta3
       pqq = (one+z**2)/omz
       
       res = (16*lomz**3)/omz + (8*omz*pisq)/3._dp +  &
            lz*(5 - 4*Li2omz*pqq + (10*pisq*pqq)/3._dp - 9*z) + Li2omz*(6 - 6*z) +  &
            lz**2*(2 - z) - 8*lomz**3*(1 + z) + 2*(-9 + 8*z) - (Li3omz*(10 - 6*z**2))/omz -  &
            (4*Li3mZ3*(5 + 3*z**2))/omz + (lz**3*(15 + 13*z**2))/(6 - 6*z) +  &
            lomz**2*(-4*omz - (lz*(-5 + 3*z**2))/omz) +  &
            Lmu2**2*(12/omz + (8*lomz)/omz - 4*lomz*(1 + z) - 4*(2 + z) -  &
            (lz*(1 + 3*z**2))/omz) +  &
            lomz*(12 - 8*lz**2*pqq - 9*z - (2*Li2omz*(-7 + z**2))/omz +  &
            (2*lz*(7 - 2*z + 7*z**2))/omz) +  &
            Lmu2*(12 + (8*Li2omz)/omz + (24*lomz)/omz + (24*lomz**2)/omz - 10*z -  &
            12*lomz**2*(1 + z) + lz**2*(1 + z) + (4*lz*(1 + z + z**2))/omz +  &
            lomz*((-8*lz*z**2)/omz - 8*(2 + z))) + (32*zeta3)/omz - 16*(1 + z)*zeta3
       
    elseif (iflag == 1) then

       res =  (-12*Lmu2**2)/omz - (24*Lmu2*lomz)/omz - (8*Lmu2**2*lomz)/omz - &
            (24*Lmu2*lomz**2)/omz - (16*lomz**3)/omz - (32*zeta3)/omz

    end if
       
  end function PqqNNLO

  function PqgNNLO(z,Lmu2) result(res)
    real(dp) :: z,Lmu2(:)
    real(dp) :: res(size(Lmu2))
    real(dp) :: omz,lomz
    real(dp) :: lz,Li2omz,Li3omz,Li3z

    res = zero
    
    omz  = one-z
    lomz = log(omz)
    
    lz = log(z)
    Li2omz = real(dilog2(one-z),kind=dp)
    Li3omz = trilog(one-z)
    Li3z = trilog(z)

    res = -17.25_dp - 3*Li2omz + (255*z)/4._dp - 49*z**2 + Li3z*(-9 + 18*z - 14*z**2) +  &
         Li3omz*(-9 + 18*z - 10*z**2) + lz**2*(-0.375_dp - z/2._dp + z**2) +  &
         (11*lomz**3*(1 - 2*z + 2*z**2))/6._dp + lz**3*(1.25_dp - (5*z)/2._dp + (7*z**2)/3._dp) +  &
         (pisq*(17 - 26*z + 18*z**2))/12._dp +  &
         lomz**2*(-7 + 21*z - 17*z**2 + lz*(2.5_dp - 5*z + z**2)) +  &
         Lmu2**2*(1.25_dp - 2*z + 3*z**2 + lz*(-0.5_dp + z - 2*z**2) +  &
         lomz*(1 - 2*z + 2*z**2)) +  &
         lz*(1 + (35*z)/4._dp - 8*z**2 + Li2omz*(-1 + 2*z - 2*z**2) +  &
         pisq*(1.5_dp - 3*z + 3*z**2)) +  &
         lomz*(16 - (109*z)/2._dp + 44*z**2 + lz**2*(-3.5_dp + 7*z - 7*z**2) +  &
         Li2omz*(5 - 10*z + 2*z**2) + (pisq*(1 - 2*z + 2*z**2))/3._dp +  &
         lz*(3 - 14*z + 14*z**2)) +  &
         Lmu2*(8 + Li2omz*(2 - 4*z) + lz**2*(0.5_dp - z) - (41*z)/2._dp + 15*z**2 +  &
         (pisq*(1 - 2*z + 2*z**2))/3._dp + lomz**2*(3 - 6*z + 6*z**2) +  &
         lz*(1.5_dp - 10*z + 10*z**2) + lomz*(-1 + 8*z - 4*z**2 - 4*lz*z**2)) +  &
         (18 - 36*z + 32*z**2)*zeta3
       
  end function PqgNNLO

  function PaqNNLO(z,Lmu2) result(res)
    real(dp) :: z,Lmu2(:)
    real(dp) :: res(size(Lmu2))
    real(dp) :: omz,lomz
    real(dp) :: lz,Li2omz,Li3omz,Li3z

    res = zero
    
    omz  = one-z
    lomz = log(omz)
    
    lz = log(z)
    Li2omz = real(dilog2(one-z),kind=dp)
    Li3omz = trilog(one-z)
    Li3z = trilog(z)

    res =  36.5_dp + lomz**2*((58 - 42/z - 27*z)/4._dp + lz*(4 - 2*z)) + (lz**3*(2 - z))/12._dp +  &
         4*Li3omz*(-2 + z) + 2*Li3z*(-2 + z) - 27/z - (29*z)/4._dp +  &
         lz**2*(2.5_dp + (17*z)/8._dp) + (lz*(-3 + 5*z))/4._dp + Li2omz*(8 + 5*z) +  &
         (pisq*(-38 + 30/z + 17*z))/12._dp + (13*lomz**3*(2 - 2*z + z**2))/(6._dp*z) +  &
         Lmu2**2*(1 + lz*(1 - z/2._dp) - z/4._dp + lomz*(-2 + 2/z + z)) +  &
         lomz*(-23 + Li2omz*(8 - 4*z) + 18/z - (3*z)/2._dp + lz*(8 + 5*z) -  &
         (5*pisq*(2 - 2*z + z**2))/(3._dp*z)) +  &
         Lmu2*(-7.5_dp + lomz*(10 + lz*(4 - 2*z) - 6/z - 4*z) + Li2omz*(4 - 2*z) +  &
         lz**2*(1 - z/2._dp) + 5/z - 2*z + lz*(4 + (5*z)/2._dp) +  &
         lomz**2*(-6 + 6/z + 3*z) - (2*pisq*(2 - 2*z + z**2))/(3._dp*z)) +  &
         (-12 + 16/z + 6*z)*zeta3
    
  end function PaqNNLO

  !-- for the ga channel
  function Paq2NNLO(z,Lmu2) result(res)
    real(dp) :: z,Lmu2(:)
    real(dp) :: res(size(Lmu2))
    real(dp) :: omz,lomz
    real(dp) :: lz,Li2z,Li3omz,Li3z

    res = zero
    
    omz  = one-z
    lomz = log(omz)
    
    lz = log(z)
    Li2z = real(dilog2(z),kind=dp)

    Li3omz = trilog(one-z)
    Li3z = trilog(z)

    res = (lz**2*(-5 - 11*z))/2._dp - 16*Li3omz*(1 + z) - 8*Li3z*(1 + z) +  &
         (lz**3*(1 + z))/3._dp + lz*(2 + 6*z) +  &
         (206 + 1056*z - 840*z**2 - 422*z**3)/(27._dp*z) +  &
         (4*pisq*(-2 - 3*z - 3*z**2 + 2*z**3))/(9._dp*z) +  &
         Li2z*(4 + 12*z - 16*lomz*(1 + z)) +  &
         Lmu2**2*(1 + 4/(3._dp*z) - z - (4*z**2)/3._dp + 2*lz*(1 + z)) +  &
         lomz**2*(-8*lz*(1 + z) + (4*(3 + 4/z - 3*z - 4*z**2))/3._dp) +  &
         Lmu2*(lz*(-2 - 6*z) - 8*Li2z*(1 + z) + 2*lz**2*(1 + z) + (4*pisq*(1 + z))/3._dp +  &
         (4*lomz*(3 + 4/z - 3*z - 4*z**2))/3._dp +  &
         (2*(-20 - 57*z + 39*z**2 + 38*z**3))/(9._dp*z)) +  &
         lomz*((8*pisq*(1 + z))/3._dp + (4*(-20 - 57*z + 39*z**2 + 38*z**3))/(9._dp*z)) +  &
         8*(1 + z)*zeta3 
     
  end function Paq2NNLO

  !--
  
  function Geq(e3,e4,sc2,eta13,eta14,eta23,eta24) result(res)
    real(dp) :: e3,e4,sc2(:),eta13,eta14,eta23,eta24
    real(dp) :: res(size(sc2))
    real(dp) :: Li2_13,Li2_14,Li2_23,Li2_24

    Li2_13 = real(dilog2(one-eta13),kind=dp)
    Li2_14 = real(dilog2(one-eta14),kind=dp)
    Li2_23 = real(dilog2(one-eta23),kind=dp)
    Li2_24 = real(dilog2(one-eta24),kind=dp)

    res = 3*log((eta13*eta24)/(eta14*eta23)) + log(eta23/eta13)*log(e3**2/sc2) + &
         log(eta14/eta24)*log(e4**2/sc2) + 2*Li2_13 - &
         2*Li2_14 - 2*Li2_23 + 2*Li2_24

  end function Geq

  function Ge2(e3,e4,ec,eta34) result(res)
    real(dp) :: res
    real(dp) :: e3,e4,ec,eta34
    real(dp) :: Li2_34

    Li2_34 = real(dilog2(one-eta34),kind=dp)

    res =  13 - (2*pisq)/3._dp + log(e3/e4)**2 + (3 - 2*log((e3*e4)/ec**2))*Log(eta34) + &
         2*Li2_34

  end function Ge2

  !-- function for gqb sector
  function Gq2b(e5,ec,mu2,eta15,eta25) result(res)
    real(dp) :: e5,ec,mu2(:),eta15,eta25
    real(dp) :: res(size(mu2))
    real(dp) :: Li2_25,L5

    Li2_25 = real(dilog2(one-eta25),kind=dp)
    L5 = log(e5/ec)

    res = 6.5_dp - pisq + L5**2 + (-1.5_dp + 2*L5)*log(eta15/eta25) + &
         2*(1.5_dp - L5)*log(4*eta25) + (3*Log(ec**2/mu2))/2._dp + &
         2*Li2_25

  end function Gq2b
  
end module mod_intsub_fc

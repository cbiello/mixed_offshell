module mod_xsects_nnlo_subonloqcd_ag_fc
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_process
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_splittings_bare
  use mod_hoppet_tools
  use mod_hoppet_nnlo_fc
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_subonloqcd_ag(2)
#endif

  private

  public :: xsect_nnlo_subonloqcd_ag_fc
  public :: xsect_nnlo_subonloqcd_ag_fc_nohoppet       !-- as in the paper, i.e. no a/c splitting
  public :: xsect_nnlo_subonloqcd_ag_fc_nohoppet_split !-- with a/c splitting

contains

  function xsect_nnlo_subonloqcd_ag_fc(yRnd,ff,vegasweight)
    use mod_kinematics_nlo
    use mod_aux_sectors
    integer :: xsect_nnlo_subonloqcd_ag_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_qg(2),kin(2)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(2,2),res_lo(2,2)
    real(dp)    :: z,s5i

    xsect_nnlo_subonloqcd_ag_fc = 0

    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       print *, 'overriding input'
    endif
#endif
    
    if (xx(xE)*(one-xx(xRHO)).lt.buff_r) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    call kinematics_nlo_is(yr=xx,HardProc=HardProc,C2Lim=C2Lim,compute_etas=.false.)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_qg(1) = zero

    else    

       call res_tree_g_qg(HardProc%AmpMom,res_nlo)
       res_nlo(:,1) = xn*Qdn2*res_nlo(:,1)
       res_nlo(:,2) = xn*Qup2*res_nlo(:,2)
       call get_respdf_hoppet(xPij,PDFs,qg_lumi,1,1,HardProc,res_nlo,respdf,myPDFs1_Lmu=[xPij_Lmu])

       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintNLO_qg(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(2) = zero
       FintNLO_qg(2) = zero

    else

       call res_tree_qqb(C2Lim%AmpMom,res_lo)
       res_lo(:,1) = xn*Qdn2*res_lo(:,1)
       res_lo(:,2) = xn*Qup2*res_lo(:,2)
       call get_respdf_hoppet(xPij,PDFs,qg_lumi,1,1,C2Lim,res_lo,respdf,myPDFs1_Lmu=[xPij_Lmu])       
       z   = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Tr*Pgq_spav(z)/(one-z))&
            * C2Lim%wgt

       FintNLO_qg(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()
    
    call check_ff(ff,xx,FintNLO_qg)
    
#if(_withchecks == 1)
    FintNNLO_subonloqcd_ag = FintNLO_qg
#endif

  end function xsect_nnlo_subonloqcd_ag_fc

  !--

  function xsect_nnlo_subonloqcd_ag_fc_nohoppet(yRnd,ff,vegasweight)
    use mod_kinematics_nlo_z_fc
    use mod_aux_sectors_fc
    use mod_intsub_fc
    integer :: xsect_nnlo_subonloqcd_ag_fc_nohoppet
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_qg(2),kin(2)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(2,2),res_lo(2,2)
    real(dp)    :: zz,s5i
    real(dp)    :: z,Lmu2(ipdf),int_sub(2,2,ipdf),restmp(1,1,ipdf),factor(1,1,ipdf)
    integer     :: i,j
    integer     :: nupdown(2)
    
    xsect_nnlo_subonloqcd_ag_fc_nohoppet = 0

    nupdown = [ndn,nup]
    
    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z = buff+onet*real(yRnd(kNLO_max+1),dp)

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       z = yRnd(kNLO_max+1)
       print *, 'overriding input'
    endif
#endif
    
    if (xx(xE)*(one-xx(xRHO)).lt.buff_r) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()

    !--

    !-- FLM[z,2,3,4,5]
    
    call kinematics_nlo_is_z(yr=xx,z1=z,z2=one,HardProc=HardProc,C2Lim=C2Lim,compute_etas=.false.)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_qg(1) = zero

    else    

       call get_Lmu2(HardProc,Lmu2)

       call res_tree_g_qg(HardProc%AmpMom,res_nlo)

       int_sub = intsub_ag_onloqcd_z2(z,Lmu2)

       restmp = zero
       do i = 1,2
          do j = 1,2
             restmp(1,1,:) = restmp(1,1,:) + int_sub(j,i,:)*res_nlo(j,i)*nupdown(i)
          enddo
       enddo
       factor(1,1,:) = restmp(1,1,:)/restmp(1,1,1)

       call get_respdf(ag_lumi,1,1,HardProc,restmp(:,:,1),respdf,factor)

       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintNLO_qg(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(2) = zero
       FintNLO_qg(2) = zero

    else

       call get_Lmu2(C2Lim,Lmu2)
       
       call res_tree_qqb(C2Lim%AmpMom,res_lo)

       int_sub = intsub_ag_onloqcd_z2(z,Lmu2)

       restmp = zero
       do i = 1,2
          do j = 1,2
             restmp(1,1,:) = restmp(1,1,:) + int_sub(j,i,:)*res_lo(j,i)*nupdown(i)
          enddo
       enddo
       factor(1,1,:) = restmp(1,1,:)/restmp(1,1,1)

       call get_respdf(ag_lumi,1,1,C2Lim,restmp(:,:,1),respdf,factor)
       
       zz  = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Tr*Pgq_spav(zz)/(one-zz))&
            * C2Lim%wgt

       FintNLO_qg(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()
    
    call check_ff(ff,xx,FintNLO_qg)
    
#if(_withchecks == 1)
    FintNNLO_subonloqcd_ag = FintNLO_qg
#endif

  end function xsect_nnlo_subonloqcd_ag_fc_nohoppet

  !-- RR split into a,c
  function xsect_nnlo_subonloqcd_ag_fc_nohoppet_split(yRnd,ff,vegasweight)
    use mod_kinematics_nlo_z_fc
    use mod_aux_sectors_fc
    use mod_intsub_fc
    integer :: xsect_nnlo_subonloqcd_ag_fc_nohoppet_split
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: HardProc,C2Lim
    !--
    real(dp)    :: xx(kNLO_max_full)
    real(dp)    :: FintNLO_qg(2),kin(2)
    real(dp)    :: respdf(ipdf)
    real(dp)    :: res_nlo(2,2),res_lo(2,2)
    real(dp)    :: zz,s5i
    real(dp)    :: z,Lmu2(ipdf),int_sub(2,2,ipdf),restmp(1,1,ipdf),factor(1,1,ipdf),w11
    integer     :: i,j
    integer     :: nupdown(2)
    
    xsect_nnlo_subonloqcd_ag_fc_nohoppet_split = 0

    nupdown = [ndn,nup]
    
    ff(1) = zero

    xx(1:kNLO_max)=buff+onet*real(yRnd(1:kNLO_max),dp)
    call random_number(xx(kNLO_max_full))
    z = buff+onet*real(yRnd(kNLO_max+1),dp)

#if (_withchecks == 1)
    if (override) then
       xx(1:kNLO_max_full) = yRnd(1:kNLO_max_full)
       z = yRnd(kNLO_max+1)
       print *, 'overriding input'
    endif
#endif
    
    if (xx(xE)*(one-xx(xRHO)).lt.buff_r) then
       failed_points = failed_points + 1
       return
    endif
    
    call open_histo()
    
    !--

    !-- FLM[z,2,3,4,5]
    
    call kinematics_nlo_is_z(yr=xx,z1=z,z2=one,HardProc=HardProc,C2Lim=C2Lim,compute_etas=.false.)
    HardProc%ids(1:5) = [0,0,id_el,-id_el,id_q]
    C2Lim%ids(1:4) = [0,0,id_el,-id_el]

    !-- Hard
    call cut_histo(HardProc)
    if (HardProc%makecut.or.HardProc%flag) then

       kin(1) = zero
       FintNLO_qg(1) = zero

    else    

       call get_Lmu2(HardProc,Lmu2)

       call res_tree_g_qg(HardProc%AmpMom,res_nlo)
       call get_w11(HardProc,w11)
       int_sub = intsub_ag_onloqcd_z2_split(z,Lmu2,w11,HardProc%Lim_KinInv)

       restmp = zero
       do i = 1,2
          do j = 1,2
             restmp(1,1,:) = restmp(1,1,:) + int_sub(j,i,:)*res_nlo(j,i)*nupdown(i)
          enddo
       enddo
       factor(1,1,:) = restmp(1,1,:)/restmp(1,1,1)

       call get_respdf(ag_lumi,1,1,HardProc,restmp(:,:,1),respdf,factor)

       respdf = respdf*HardProc%wgt

       kin(1) = respdf(1)
       FintNLO_qg(1) = kin(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- C2
    call cut_histo(C2Lim)
    
    if (C2Lim%makecut.or.C2Lim%flag) then

       kin(2) = zero
       FintNLO_qg(2) = zero

    else

       call get_Lmu2(C2Lim,Lmu2)
       
       call res_tree_qqb(C2Lim%AmpMom,res_lo)
       call get_w11(C2Lim,w11)
       int_sub = intsub_ag_onloqcd_z2_split(z,Lmu2,w11,C2Lim%Lim_KinInv)

       restmp = zero
       do i = 1,2
          do j = 1,2
             restmp(1,1,:) = restmp(1,1,:) + int_sub(j,i,:)*res_lo(j,i)*nupdown(i)
          enddo
       enddo
       factor(1,1,:) = restmp(1,1,:)/restmp(1,1,1)

       call get_respdf(ag_lumi,1,1,C2Lim,restmp(:,:,1),respdf,factor)
       
       zz  = C2Lim%Lim_KinInv(1)
       s5i = C2Lim%Lim_KinInv(2)

       !-- use Pqg so that in the soft limit we do not have 1/[1-(1-x)]
       respdf = (-one)*respdf*(two/s5i)*(Tr*Pgq_spav(zz)/(one-zz))&
            * C2Lim%wgt

       FintNLO_qg(2) = respdf(1)
       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = sum(kin)
    call close_histo()
    
    call check_ff(ff,xx,FintNLO_qg)
    
#if(_withchecks == 1)
    FintNNLO_subonloqcd_ag = FintNLO_qg
#endif

  end function xsect_nnlo_subonloqcd_ag_fc_nohoppet_split
  
  !-- TC damping factor with 6||1
  subroutine get_w11(proc,damp)
    type(KinConfig), intent(in) :: proc
    real(dp), intent(out) :: damp
    real(dp) :: eta15,eta25
    
    eta15 = proc%Lim_etaij(1,5)
    eta25 = proc%Lim_etaij(2,5)
    
    damp = eta25*(one+eta15)
    
  end subroutine get_w11
  
end module mod_xsects_nnlo_subonloqcd_ag_fc
  

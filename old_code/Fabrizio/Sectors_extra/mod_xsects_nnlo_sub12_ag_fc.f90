module mod_xsects_nnlo_sub12_ag_fc
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_kinematics_gen
  use mod_kinematics_lo_z_fc
  use mod_process
  use mod_auxfunctions
  use mod_aux_sectors_fc
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_intsub_fc
  use mod_hoppet_tools
  use mod_hoppet_nnlo_fc
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_sub12_ag_fc(2)
#endif

  private

  public :: xsect_nnlo_sub12_ag_fc

contains

  function xsect_nnlo_sub12_ag_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_sub12_ag_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(2),respdf(ipdf)
    real(dp)    :: res_lo(2,2),int_sub(2,2,ipdf)
    real(dp)    :: z,zb
    real(dp)    :: Lmu2(ipdf)
    real(dp)    :: factor(1,1,ipdf),restmp(1,1,ipdf)
    integer     :: i,j
    integer     :: nupdown(2)

    xsect_nnlo_sub12_ag_fc = 0

    nupdown = [ndn,nup]
    
    ff(1) = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))
    z  = buff+onet*real(yRnd(kLO_max+1),dp)
    zb = buff+onet*real(yRnd(kLO_max+2),dp)
    
#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       z  = yRnd(kLO_max+1)
       zb = yRnd(kLO_max+2)
       print *, 'overriding input'
    endif
#endif

    call open_histo()
    
    !-- [z,zb]
    call kinematics_lo_z_fc(xx,z,zb,LOProc)
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero
       
    else    

       call get_Lmu2(LOProc,Lmu2)

       call res_tree_qqb(LOProc%AmpMom,res_lo)
       int_sub = intsub_ag_sub12_zzb(z,zb,Lmu2)

       restmp = zero
       do i = 1,2
          do j = 1,2
             restmp(1,1,:) = restmp(1,1,:) + int_sub(j,i,:)*res_lo(j,i)*nupdown(i)
          enddo
       enddo
       factor(1,1,:) = restmp(1,1,:)/restmp(1,1,1)
             
       call get_respdf(ag_lumi,1,1,LOProc,restmp(:,:,1),respdf,factor)

       respdf = respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[z,zb] 1', respdf
      
    endif

    !-- [1,zb]
    call kinematics_lo_z_fc(xx,one,zb,LOProc)
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(2) = zero
       
    else    

       call get_Lmu2(LOProc,Lmu2)
       
       call res_treeAA_aa(LOProc%AmpMom,res_lo(1,1))
       int_sub = intsub_ag_sub12_1zb(zb,Lmu2)

       restmp = zero
       !-- qb channel already accounted for
       do i = 1,2
          restmp(1,1,:) = restmp(1,1,:) + int_sub(1,i,:)*res_lo(1,1)*nupdown(i)
       enddo
       factor(1,1,:) = restmp(1,1,:)/restmp(1,1,1)
       
       call get_respdf(ag_lumi,1,1,LOProc,restmp(:,:,1),respdf,factor)

       respdf = respdf*LOProc%wgt

       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 2', respdf
       !pause
       
    endif
    
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNNLO_sub12_ag_fc = kin
#endif

  end function xsect_nnlo_sub12_ag_fc
  
end module mod_xsects_nnlo_sub12_ag_fc
  

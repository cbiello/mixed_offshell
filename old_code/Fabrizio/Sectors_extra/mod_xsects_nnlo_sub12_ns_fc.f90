module mod_xsects_nnlo_sub12_ns_fc
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
  real(dp), public, save :: FintNNLO_sub12_fc(4)
#endif

  private

  public :: xsect_nnlo_sub12_test !-- nloqcd subtraction, as a test
  public :: xsect_nnlo_sub12_ns_ga_fc
  public :: xsect_nnlo_sub12_ns_qqbORqq_fc

contains

  !--
  
  function xsect_nnlo_sub12_test(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_sub12_test
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(3),respdf(ipdf)
    real(dp)    :: res_lo(2,2),int_sub(2,2,ipdf)
    real(dp)    :: z,zb
    real(dp)    :: Lmu2(ipdf)

    xsect_nnlo_sub12_test = 0

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

    !-- [z,2]
    call kinematics_lo_z_fc(xx,z,one,LOProc)
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero
       
    else    

       call get_Lmu2(LOProc,Lmu2)
       
       call res_tree_qqb(LOProc%AmpMom,res_lo)
       int_sub = intsub_test_z(z,Lmu2,0)
       call get_respdf(ns_lumi,1,0,LOProc,res_lo,respdf,int_sub)

       respdf = respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- [1,zb]
    call kinematics_lo_z_fc(xx,one,zb,LOProc)
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(2) = zero
       
    else    

       call get_Lmu2(LOProc,Lmu2)
       
       call res_tree_qqb(LOProc%AmpMom,res_lo)
       int_sub = intsub_test_z(zb,Lmu2,0)
       call get_respdf(ns_lumi,1,0,LOProc,res_lo,respdf,int_sub)

       respdf = respdf*LOProc%wgt

       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif

    !-- [1,2]
    call kinematics_lo_z_fc(xx,one,one,LOProc)
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(3) = zero
       
    else    

       call get_Lmu2(LOProc,Lmu2)
       
       call res_tree_qqb(LOProc%AmpMom,res_lo)
       int_sub = intsub_test_z(z,Lmu2,1) + &
            intsub_test_z(zb,Lmu2,1) + &
            intsub_test_12(Lmu2)
       call get_respdf(ns_lumi,1,0,LOProc,res_lo,respdf,int_sub)

       respdf = respdf*LOProc%wgt

       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
        
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNNLO_sub12_fc(1:3) = kin
#endif

  end function xsect_nnlo_sub12_test

  !--

  
  function xsect_nnlo_sub12_ns_ga_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_sub12_ns_ga_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(4),respdf(ipdf)
    real(dp)    :: res_lo(2,2),int_sub(2,2,ipdf)
    real(dp)    :: z,zb
    real(dp)    :: Lmu2(ipdf),my_mu2(ipdf)

    xsect_nnlo_sub12_ns_ga_fc = 0

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
       int_sub = intsub_qqb_ga_sub12_zzb(z,zb,Lmu2,0,0)
       call get_respdf(ns_lumi,1,1,LOProc,res_lo,respdf,int_sub)

       respdf = respdf*LOProc%wgt

       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[z,zb] 1', respdf
      
    endif

    !-- [z,2]
    call kinematics_lo_z_fc(xx,z,one,LOProc)
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(2) = zero
       
    else    

       call get_Lmu2(LOProc,Lmu2,my_mu2)
       
       call res_tree_qqb(LOProc%AmpMom,res_lo)
       int_sub = intsub_qqb_ga_sub12_zzb(z,zb,Lmu2,0,1) + &
            intsub_qqb_ga_sub12_z2(z,my_mu2,Lmu2,0,LOProc%Lim_KinInv2)
       call get_respdf(ns_lumi,1,1,LOProc,res_lo,respdf,int_sub)

       respdf = respdf*LOProc%wgt

       kin(2) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[z,2 ] 2', respdf
       
    endif

    !-- [1,zb]
    call kinematics_lo_z_fc(xx,one,zb,LOProc)
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(3) = zero
       
    else    

       call get_Lmu2(LOProc,Lmu2,my_mu2)
       
       call res_tree_qqb(LOProc%AmpMom,res_lo)
       int_sub = intsub_qqb_ga_sub12_zzb(z,zb,Lmu2,1,0) + &
            intsub_qqb_ga_sub12_1zb(zb,my_mu2,Lmu2,0,LOProc%Lim_KinInv2)
       call get_respdf(ns_lumi,1,1,LOProc,res_lo,respdf,int_sub)

       respdf = respdf*LOProc%wgt

       kin(3) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,zb] 3', respdf
       
    endif

    !-- [1,2]
    call kinematics_lo_z_fc(xx,one,one,LOProc)
    LOProc%ids(1:4) = [0,0,id_el,-id_el]
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(4) = zero
       
    else    

       call get_Lmu2(LOProc,Lmu2,my_mu2)
       
       call res_tree_qqb(LOProc%AmpMom,res_lo)
       int_sub = intsub_qqb_ga_sub12_zzb(z,zb,Lmu2,1,1) + &
            intsub_qqb_ga_sub12_1zb(zb,my_mu2,Lmu2,1,LOProc%Lim_KinInv2) + &
            intsub_qqb_ga_sub12_z2(z,my_mu2,Lmu2,1,LOProc%Lim_KinInv2) + &
            intsub_qqb_ga_sub12_12(Lmu2,LOProc%Lim_KinInv2)

       call get_respdf(ns_lumi,1,1,LOProc,res_lo,respdf,int_sub)

       respdf = respdf*LOProc%wgt

       kin(4) = respdf(1)

       call fill_histo(respdf,vegasweight)

       !print *, '[1,2 ] 4', respdf
       !pause
       
    endif
    
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNNLO_sub12_fc = kin
#endif

  end function xsect_nnlo_sub12_ns_ga_fc

  !-- qqb and qq identical, only different splitting function
  function xsect_nnlo_sub12_ns_qqbORqq_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_sub12_ns_qqbORqq_fc
    real(dp15) :: yRnd(30),ff(1),vegasweight
    type(KinConfig) :: LOProc
    !--
    real(dp)    :: xx(kLO_max_full)
    real(dp)    :: kin(1),respdf(ipdf),respdf_1(ipdf),respdf_2(ipdf)
    real(dp)    :: res_lo(2,2)

    xsect_nnlo_sub12_ns_qqbORqq_fc = 0

    ff(1) = zero

    xx(1:kLO_max)=buff+onet*real(yRnd(1:kLO_max),dp)
    call random_number(xx(kLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kLO_max_full) = yRnd(1:kLO_max_full)
       print *, 'overriding input'
    endif
#endif

    call open_histo()

    call kinematics_lo(xx,LOProc)
    call cut_histo(LOProc)

    if (LOProc%makecut.or.LOProc%flag) then

       kin(1) = zero
       
    else    

       call res_tree_qqb(LOProc%AmpMom,res_lo)
       res_lo(:,1) = res_lo(:,1) * Qdn2
       res_lo(:,2) = res_lo(:,2) * Qup2
       
       call get_respdf_hoppet(xPij,PDFs,ns_lumi,1,1,LOProc,res_lo,respdf_1,myPDFs1_Lmu=[xPij_Lmu])
       respdf_1 = Cf*respdf_1*LOProc%wgt

       call get_respdf_hoppet(PDFs,xPij,ns_lumi,1,1,LOProc,res_lo,respdf_2,myPDFs2_Lmu=[xPij_Lmu])
       respdf_2 = Cf*respdf_2*LOProc%wgt

       respdf = respdf_1 + respdf_2
       
       kin(1) = respdf(1)

       call fill_histo(respdf,vegasweight)

    endif
    
    ff(1) = kin(1)
    call close_histo()

    call check_ff(ff,xx,kin)
    
#if(_withchecks == 1)
    FintNNLO_sub12_fc(1:1) = kin
#endif

  end function xsect_nnlo_sub12_ns_qqbORqq_fc

  
end module mod_xsects_nnlo_sub12_ns_fc
  

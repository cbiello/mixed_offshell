module mod_xsects_nnlo_rr_ag_fc
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_aux_sectors
  use mod_process
  use mod_histo
  use mod_cut_histo
  use mod_amplitudes_tree_ppll
  use mod_splittings_bare
  use mod_eikonals
  use mod_partitions
  use mod_auxfunctions
  implicit none
#if(_withchecks == 1)
  real(dp), public, save :: FintNNLO_rr_tc_ag_fc(4)
#endif

  public :: xsect_nnlo_rr_5262a_ag_fc,xsect_nnlo_rr_5262c_ag_fc
  public :: xsect_nnlo_rr_5262b_ag_fc,xsect_nnlo_rr_5262d_ag_fc

  public :: xsect_nnlo_rr_5161_ag_fc !-- as in the paper, no splitting
  
  private

  !! -- a(p1) g(p2) -> e-(p3) e+(p4) q(p5) qb(p6) 

contains

  function xsect_nnlo_rr_5262b_ag_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5262b_ag_fc
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5262b_ag_fc = xsect_nnlo_rr_5262bd_ag_fc(&
      yRnd,ff,vegasweight,2)

  end function xsect_nnlo_rr_5262b_ag_fc

  function xsect_nnlo_rr_5262d_ag_fc(yRnd,ff,vegasweight)
    integer :: xsect_nnlo_rr_5262d_ag_fc
    real(dp15) :: yRnd(30), ff(1), vegasweight

    xsect_nnlo_rr_5262d_ag_fc = xsect_nnlo_rr_5262bd_ag_fc(&
      yRnd,ff,vegasweight,4)

  end function xsect_nnlo_rr_5262d_ag_fc

  !!*************************************************************************!!

  !-- sector as in the paper right now, i.e. no splitting
  function xsect_nnlo_rr_5161_ag_fc(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_dc_ii !-- only for parameters etc, hard copied below
    implicit none
    integer :: xsect_nnlo_rr_5161_ag_fc
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: C5Lim,C6Lim
    !--
    real(dp) :: xx(kNNLO_max_full)
    real(dp) :: FintNNLO_ag_fc(3),kin(3)
    real(dp) :: respdf(ipdf)
    real(dp) :: res_tmp1(1,1),res_tmp2(2,2)
    real(dp) :: z5,z6,s15,s16
    real(dp) :: damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5161_ag_fc = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    !! Use this even if no soft singularities
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    !xx(xRHO6) = 1E-10_dp
    
    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_ii_5i6i_fc(xx(1:kNNLO_max_full), 1, 2, &
         HardProc=HardProc,C5Lim=C5Lim,C6Lim=C6Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_q,-id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

       kin(1) = zero
       FintNNLO_ag_fc(1) = zero
       
    else
       
       call res_tree_ga_ag(HardProc%AmpMom,res_tmp1(1,1))
       call get_respdf(ag_lumi,1,1,HardProc,res_tmp1,respdf)
       
       call partition_nnlo_qcd(HardProc,damp,iconf,1,1)
       
       respdf = respdf * HardProc%wgt * damp
       
       FintNNLO_ag_fc(1) = respdf(1)
       kin(1) = respdf(1)
       call fill_histo(respdf,vegasweight)

       !print *, '1 ',FintNNLO_ag_fc(1)

    endif

    !!-----------------------------------------------------------------------!!
    !!                               C5                                      !!
    !!-----------------------------------------------------------------------!!
    
    C5Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C5Lim%npart = 5
    
    call cut_histo(C5Lim)
    
    if (C5Lim%makecut) then
       
       kin(2) = zero
       FintNNLO_ag_fc(2) = zero
       
    else
       
       call res_tree_g_qg(C5Lim%AmpMom,res_tmp2)
       res_tmp1(1,1) = res_tmp2(2,1)*ndn*Qdn2 + res_tmp2(2,2)*nup*Qup2
       call get_respdf(ag_lumi,1,1,C5Lim,res_tmp1,respdf)
       
       call partition_nnlo_qcd(C5Lim,damp,iconf,1,1)
       
       s15 = C5Lim%Lim_sij(1,5)
       z5  = C5Lim%Lim_z(1)

       respdf = -respdf*xn*(two/s15*Pqg(z5))*C5Lim%wgt*damp
       
       FintNNLO_ag_fc(2) = respdf(1)
       kin(2) = respdf(1)
       
       call fill_histo(respdf,vegasweight)
       
       !print *, '2 ',FintNNLO_ag_fc(2),FintNNLO_ag_fc(2)/FintNNLO_ag_fc(1),one+FintNNLO_ag_fc(2)/FintNNLO_ag_fc(1)
       
    endif

    
    !!-----------------------------------------------------------------------!!
    !!                               C6                                      !!
    !!-----------------------------------------------------------------------!!
    
    C6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C6Lim%npart = 5
    
    call cut_histo(C6Lim)
    
    if (C6Lim%makecut) then
       
       kin(3) = zero
       FintNNLO_ag_fc(3) = zero
       
    else

       call res_tree_g_qg(C6Lim%AmpMom,res_tmp2)
       res_tmp1(1,1) = res_tmp2(1,1)*ndn*Qdn2 + res_tmp2(1,2)*nup*Qup2
       call get_respdf(ag_lumi,1,1,C6Lim,res_tmp1,respdf)
       
       call partition_nnlo_qcd(C6Lim,damp,iconf,1,1)
       
       s16 = C6Lim%Lim_sij(1,6)
       z6  = C6Lim%Lim_z(2)
       
       respdf = -respdf*xn*(two/s16*Pqg(z6))*C6Lim%wgt*damp
       
       FintNNLO_ag_fc(3) = respdf(1)
       kin(3) = respdf(1)
       
       call fill_histo(respdf,vegasweight)

       !print *, '3 ',FintNNLO_ag_fc(3),FintNNLO_ag_fc(3)/FintNNLO_ag_fc(1),one+FintNNLO_ag_fc(3)/FintNNLO_ag_fc(1)
       
    endif
    
    ff(1) = sum(kin)
    call close_histo()
    
    call check_ff(ff,xx,FintNNLO_ag_fc)
    
#if(_withchecks == 1)
    FintNNLO_rr_tc_ag_fc(1:3) = FintNNLO_ag_fc
#endif

    !if (ff(1).ne.zero) then
    !   print *, 'ff', ff(1)
    !   pause
    !endif

  contains

    !-- j is the other direction
    subroutine kinematics_nnlo_ii_5i6i_fc(yr, i, j,  &
         HardProc, C5Lim, C6Lim)
      use mod_kinematics_gen
      use mod_kinematics_lept
      use mod_aux_kinematics
      integer, parameter :: NNLOdim_vegas = tauy_max + 2 + 3 + 3
      !
      integer, parameter :: kNNLO_min = tauy_max + 1
      integer, parameter :: kNNLO_max = tauy_max + 2 + 3 + 3
      integer, parameter :: kNNLO_max_full = kNNLO_max + 2
      !
      integer, parameter :: xE5 = kNNLO_min
      integer, parameter :: xE6 = kNNLO_min + 1
      integer, parameter :: xRHO5 = kNNLO_min + 2
      integer, parameter :: xRHO6 = kNNLO_min + 4
      !!
      integer, parameter :: xPHI5   = kNNLO_min + 3
      integer, parameter :: xPHI6   = kNNLO_min + 5
      !
      integer, parameter :: xLMIN  = kNNLO_min + 6
      integer, parameter :: xLMAX  = kNNLO_min + 8
      !
      integer, parameter :: npart_lo = 4, npart_nlo = 5, npart_nnlo = 6
      !
      real(dp),        intent(in)    :: yr(:)
      integer,         intent(in)    :: i, j
      type(KinConfig), intent(inout) :: HardProc
      type(KinConfig), optional, intent(inout) :: C5Lim, C6Lim
      !---
      real(dp) :: mv, mv2, tau, ylab, jac
      real(dp) :: p1(4), p2(4), p5(4), p6(4), pv(4)
      real(dp) :: sgn, kallenF
      real(dp) :: x1, x2, x3, x4, x5, x6
      real(dp) :: cos5i, sin5i, cos6i, sin6i
      real(dp) :: phi5, phi6, cos_phi5, cos_phi6, sin_phi5, sin_phi6
      real(dp) :: xi1, xi2, spart, sqrts
      real(dp) :: nlept(3,2), ns(4,6)
      real(dp) :: eta5i, eta6i, eta5j, eta6j
      real(dp) :: eta51, eta61, eta52, eta62, eta65
      
      !-- process-dependent part
      call get_tauy(yr(1:tauy_max),tau,ylab,jac)
      mv2 = tau*sh
      mv  = sqrt(mv2)
      
      !-- Initialization of kinematic configurations
      call initialize_config(HardProc)                             !  1
      if(present(C5Lim))       call initialize_config(C5Lim)       !  2
      if(present(C6Lim))       call initialize_config(C6Lim)       !  3

      x1 = yr(xE5)   ! E_5
      x2 = yr(xE6)   ! E_6
      x3 = yr(xRHO5) ! \theta_{5i}
      x4 = yr(xRHO6) ! \theta_{6i}
      x5 = yr(xPHI5) ! \phi_{5}
      x6 = yr(xPHI6) ! related to \lambda
      
      !! Choose sign of \cos\theta depending on emitter 
      !! {i,j} = {1,2} or {i,j} = {2,1}
      if(i .eq. 1) then
         sgn = +one
      elseif(i .eq. 2) then
         sgn = -one
      else
         print *, 'wrong emitter'
         sgn = zero
         stop
      endif
      
      !!-----------------------------------------------------------------------!!
      !!                      STRUCTURE 1: Hard Process.                       !!
      !!-----------------------------------------------------------------------!!
      
      !! Azimuthal variables
      phi5 = twopi*x5 
      phi6 = twopi*x6
      cos_phi5 = cos(phi5)
      sin_phi5 = sin(phi5)
      cos_phi6 = cos(phi6)
      sin_phi6 = sin(phi6)
      
      !! Polar variables
      eta5i = x3
      eta5j = one - x3
      eta6i = x4
      eta6j = one - x4
      
      cos5i = sgn*( one - 2*x3)
      sin5i = two*sqrt(x3*(one-x3))
      
      cos6i = sgn*( one - 2*x4)
      sin6i = two*sqrt(x4*(one-x4))
      
      eta65 = half*(one - sin5i*sin6i*cos(phi5-phi6) - cos5i*cos6i)
      
      !! set eta variables in the limit
      HardProc%Lim_etaij(1,2) = one
      HardProc%Lim_etaij(i,5) = eta5i
      HardProc%Lim_etaij(j,5) = eta5j
      HardProc%Lim_etaij(i,6) = eta6i
      HardProc%Lim_etaij(j,6) = eta6j
      HardProc%Lim_etaij(5,6) = eta65
      
      eta51 = HardProc%Lim_etaij(1,5)
      eta52 = HardProc%Lim_etaij(2,5)
      eta61 = HardProc%Lim_etaij(1,6)
      eta62 = HardProc%Lim_etaij(2,6)
      
      ns(:,1) = [one,zero,zero, one]
      ns(:,2) = [one,zero,zero,-one]
      
      ns(:,5) = [one,sin5i*cos_phi5,sin5i*sin_phi5,cos5i]
      ns(:,6) = [one,sin6i*cos_phi6,sin6i*sin_phi6,cos6i]
      
      !! Compute partonic fractions {xi1, xi2}, partonic s and its square root
      !! and check if point is physical
      call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
           spart,sqrts,xi1,xi2,HardProc%flag)
      
      if (.not. HardProc%flag) then

         HardProc%npart = npart_nnlo
         
         HardProc%PartFrac = (/xi1,xi2/)
         
         p1 = half*sqrts*ns(:,1)
         p2 = half*sqrts*ns(:,2)
         p5 = half*sqrts*x1*ns(:,5)
         p6 = half*sqrts*x2*ns(:,6)
         
         pv = p1 + p2 - p5 - p6
         
         call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, HardProc)
         ns(:,3) = (/one,nlept(:,1)/)
         ns(:,4) = (/one,nlept(:,2)/)
         
         HardProc%AmpMom(:,1) = p1
         HardProc%AmpMom(:,2) = p2
         HardProc%AmpMom(:,5) = p5
         HardProc%AmpMom(:,6) = p6
         
         HardProc%wgt = one/8.0_dp/pi * kallenF &
              * (spart/two)**2 & 
              * x1*x2          & 
              * one/two/spart  &
              * spart/mv2      &
              * jac
         
      endif

      !!-----------------------------------------------------------------------!!
      !!                        STRUCTURE 2: C5i FLM.                          !!
      !!-----------------------------------------------------------------------!!
      
      if(present(C5Lim)) then
         
         !! x3 -> 0
         C5Lim%Lim_etaij(1,2) = one
         C5Lim%Lim_etaij(i,5) = zero
         C5Lim%Lim_etaij(j,5) = one
         C5Lim%Lim_etaij(i,6) = x4
         C5Lim%Lim_etaij(j,6) = one - x4
         C5Lim%Lim_etaij(5,6) = x4
         
         eta51 = C5Lim%Lim_etaij(1,5)
         eta52 = C5Lim%Lim_etaij(2,5)
         eta61 = C5Lim%Lim_etaij(1,6)
         eta62 = C5Lim%Lim_etaij(2,6)
         eta65 = x4
         
         call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
              spart,sqrts,xi1,xi2,C5Lim%flag)
         
         if (.not. C5Lim%flag) then 
            
            C5Lim%npart = npart_nlo         
            
            C5Lim%PartFrac = (/xi1,xi2/)
            
            ns(:,5) = [one,zero,zero,sgn]
            ns(:,6) = [one,sin6i*cos_phi6,sin6i*sin_phi6,cos6i]
            
            p1 = half*sqrts*ns(:,1)
            p2 = half*sqrts*ns(:,2)
            p5 = half*sqrts*x1*ns(:,5)
            p6 = half*sqrts*x2*ns(:,6)
            
            if(i .eq. 1) then
               p1 = p1 - p5
            elseif(i .eq. 2) then
               p2 = p2 - p5
            endif
            
            pv = p1 + p2 - p6
            call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C5Lim)
            ns(:,3) = (/one,nlept(:,1)/)
            ns(:,4) = (/one,nlept(:,2)/)
            
            C5Lim%AmpMom(:,1) = p1
            C5Lim%AmpMom(:,2) = p2
            C5Lim%AmpMom(:,5) = p6
            
            C5Lim%Lim_z(1) = one/(one - x1)
            C5Lim%Lim_sij(i,5) = spart*x1*HardProc%Lim_etaij(i,5)
            
            C5Lim%wgt = one/8.0_dp/pi * kallenF &
                 * (spart/two)**2 & 
                 * x1*x2          & 
                 * one/two/spart  & 
                 * spart/MV2      &
                 * jac
            
         endif
      endif
      
      !!-----------------------------------------------------------------------!!
      !!                        STRUCTURE 5: C6i FLM.                          !!
      !!-----------------------------------------------------------------------!!
      
      if(present(C6Lim)) then
         
         !! x4 -> 0
         C6Lim%Lim_etaij(1,2) = one
         C6Lim%Lim_etaij(i,5) = x3
         C6Lim%Lim_etaij(j,5) = one - x3
         C6Lim%Lim_etaij(i,6) = zero
         C6Lim%Lim_etaij(j,6) = one
         C6Lim%Lim_etaij(5,6) = x3

         eta51 = C6Lim%Lim_etaij(1,5)
         eta52 = C6Lim%Lim_etaij(2,5)
         eta61 = C6Lim%Lim_etaij(1,6)
         eta62 = C6Lim%Lim_etaij(2,6)
         eta65 = x3
         
         call get_part_s_xi_nnlo(mv2,ylab,x1,x2,eta65,eta51,eta61,eta52,eta62, &
              spart,sqrts,xi1,xi2,C6Lim%flag)
         
         if (.not. C6Lim%flag) then
            
            C6Lim%npart = npart_nlo
            
            C6Lim%PartFrac = (/xi1,xi2/)
            
            ns(:,5) = [one,sin5i*cos_phi5,sin5i*sin_phi5,cos5i]
            ns(:,6) = [one,zero,zero,sgn]
            
            p1 = half*sqrts*ns(:,1)
            p2 = half*sqrts*ns(:,2)
            p5 = half*sqrts*x1*ns(:,5)
            p6 = half*sqrts*x2*ns(:,6)
            
            if(i .eq. 1) then
               p1 = p1 - p6
            elseif(i .eq. 2) then
               p2 = p2 - p6
            endif

            pv = p1 + p2 - p5
            call get_lept_mom(mv, pv, yr(xLMIN:xLMAX), nlept, kallenF, C6Lim)
            ns(:,3) = (/one,nlept(:,1)/)
            ns(:,4) = (/one,nlept(:,2)/)
            
            C6Lim%AmpMom(:,1) = p1
            C6Lim%AmpMom(:,2) = p2
            C6Lim%AmpMom(:,5) = p5
            
            C6Lim%Lim_z(2) = one/(one - x2)
            C6Lim%Lim_sij(i,6) = spart*x2*HardProc%Lim_etaij(i,6)
            
            C6Lim%wgt = one/8.0_dp/pi * kallenF &
                 * (spart/two)**2 & 
                 * x1*x2          & 
                 * one/two/spart  & 
                 * spart/MV2      &
                 * jac
            
         endif
      end if

    end subroutine kinematics_nnlo_ii_5i6i_fc
      
  end function xsect_nnlo_rr_5161_ag_fc
  
  function xsect_nnlo_rr_5262a_ag_fc(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5262a_ag_fc
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: TCLim, TCC6Lim,C6Lim
    !--
    real(dp) :: xx(kNNLO_max_full)
    real(dp) :: FintNNLO_ag_fc(4),kin(3)
    real(dp) :: respdf(ipdf),respdf_vect(2,ipdf)
    real(dp) :: res_tmp1(1,1),res_tmp2(2,2)
    real(dp) :: z1,z5,z6,s25,s26,s56,si2
    real(dp) :: damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5262a_ag_fc = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_abcd_5i6iac(&
         xx(1:kNNLO_max_full), 2, 1, 1, &
         HardProc,TCLim=TCLim,TCC6Lim=TCC6Lim,C6Lim=C6Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_q,-id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_ag_fc(1) = zero

    else

      call res_tree_ga_ag(HardProc%AmpMom,res_tmp1(1,1))
      call get_respdf(ag_lumi,1,1,HardProc,res_tmp1,respdf)

      call partition_nnlo_qcd(HardProc,damp,iconf,2,2)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_ag_fc(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

      !print *, '1 ',FintNNLO_ag_fc(1)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          TC + TCC6                                    !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(2) = zero
      FintNNLO_ag_fc(2:3) = zero

    else
       
      call res_treeAA_aa(TCLim%AmpMom,res_tmp1(1,1))
      res_tmp1 = res_tmp1 * (Qup2*nup + Qdn2*ndn) !-- qqb/qbq taken care in Pgqqb

      call get_respdf(ag_lumi,1,1,TCLim,res_tmp1,respdf)

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s25 = TCLim%Lim_sij(2,5)
      s26 = TCLim%Lim_sij(2,6)

      !-- damp = one
      respdf_vect(1,:) = -Tr*(Pgqqb_spav_ab(-s25,-s26,s56,z1,z5,z6))*TCLim%wgt*respdf
      FintNNLO_ag_fc(2) = respdf_vect(1,1)

      !-- TC + C6
      !respdf_vect_part(2,:) = respdf*(Cf*two*Pqg(z5)/si5)*(two*Pqg(z6)/si6)*TCC6Lim%wgt
      z5 = TCC6Lim%Lim_z(1)
      z6 = TCC6Lim%Lim_z(2)
      si2 = TCC6Lim%Lim_KinInv(1)
      respdf_vect(2,:) = +Tr*(Pgq_spav(z5)*two/si2)*(Pqq(z6)*two/s26)*TCC6Lim%wgt*respdf
      FintNNLO_ag_fc(3) = respdf_vect(2,1)

      respdf = sum(respdf_vect(1:2,:),1)
      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

      !print *, '2 ',FintNNLO_ag_fc(2)
      !print *, '3 ',FintNNLO_ag_fc(3)

    endif

    !!-----------------------------------------------------------------------!!
    !!                               C6                                      !!
    !!-----------------------------------------------------------------------!!

    C6Lim%ids = [0,0,id_el,-id_el,id_q,0]
    C6Lim%npart = 5

    call cut_histo(C6Lim)

    if (C6Lim%makecut) then

      kin(3) = zero
      FintNNLO_ag_fc(4) = zero

    else
       
      call res_tree_a_aq(C6Lim%AmpMom,res_tmp2)
      res_tmp1(1,1) = res_tmp2(2,1)*ndn + res_tmp2(2,2)*nup
      call get_respdf(ag_lumi,1,1,C6Lim,res_tmp1,respdf)

      call partition_nnlo_qcd(C6Lim,damp,iconf,2,2)

      s26 = C6Lim%Lim_sij(2,6)
      z6  = C6Lim%Lim_z(2)
      
      respdf = -respdf*Tr*(two/s26*Pqg(z6))*C6Lim%wgt*damp
            
      FintNNLO_ag_fc(4) = respdf(1)
      kin(3) = respdf(1)

      call fill_histo(respdf,vegasweight)
      
      !print *, '4 ',FintNNLO_ag_fc(4)
      
    endif

    !!-----------------------------------------------------------------------!!
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ag_fc)

#if(_withchecks == 1)
    FintNNLO_rr_tc_ag_fc = FintNNLO_ag_fc
#endif

    !if (ff(1).ne.zero) then
    !   print *, 'ff', ff(1)
    !   pause
    !endif

  end function xsect_nnlo_rr_5262a_ag_fc

  function xsect_nnlo_rr_5262c_ag_fc(yRnd,ff,vegasweight)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5262c_ag_fc
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: TCLim, TCC5Lim,C5Lim
    !--
    real(dp) :: xx(kNNLO_max_full)
    real(dp) :: FintNNLO_ag_fc(4),kin(3)
    real(dp) :: respdf(ipdf),respdf_vect(2,ipdf)
    real(dp) :: res_tmp1(1,1),res_tmp2(2,2)
    real(dp) :: z1,z5,z6,s25,s26,s56,si2
    real(dp) :: damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5262c_ag_fc = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_abcd_5i6iac(&
         xx(1:kNNLO_max_full), 2, 1, 3, &
         HardProc,TCLim=TCLim,TCC5Lim=TCC5Lim,C5Lim=C5Lim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_q,-id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_ag_fc(1) = zero

    else

      call res_tree_ga_ag(HardProc%AmpMom,res_tmp1(1,1))
      call get_respdf(ag_lumi,1,1,HardProc,res_tmp1,respdf)

      call partition_nnlo_qcd(HardProc,damp,iconf,2,2)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_ag_fc(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

      !print *, '1 ',FintNNLO_ag_fc(1)

    endif

    !!-----------------------------------------------------------------------!!
    !!                          TC + TCC5                                    !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(2) = zero
      FintNNLO_ag_fc(2:3) = zero

    else
       
      call res_treeAA_aa(TCLim%AmpMom,res_tmp1(1,1))
      res_tmp1 = res_tmp1 * (Qup2*nup + Qdn2*ndn) !-- qqb/qbq taken care in Pgqqb

      call get_respdf(ag_lumi,1,1,TCLim,res_tmp1,respdf)

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s25 = TCLim%Lim_sij(2,5)
      s26 = TCLim%Lim_sij(2,6)

      !-- damp = one
      respdf_vect(1,:) = -Tr*(Pgqqb_spav_ab(-s25,-s26,s56,z1,z5,z6))*TCLim%wgt*respdf
      FintNNLO_ag_fc(2) = respdf_vect(1,1)

      !-- TC + C5
      !respdf_vect_part(2,:) = respdf*(Cf*two*Pqg(z5)/si5)*(two*Pqg(z6)/si6)*TCC6Lim%wgt
      z5 = TCC5Lim%Lim_z(1)
      z6 = TCC5Lim%Lim_z(2)
      si2 = TCC5Lim%Lim_KinInv(1)
      respdf_vect(2,:) = +Tr*(Pgq_spav(z6)*two/si2)*(Pqq(z5)*two/s25)*TCC5Lim%wgt*respdf
      FintNNLO_ag_fc(3) = respdf_vect(2,1)

      respdf = sum(respdf_vect(1:2,:),1)
      kin(2) = respdf(1)
      call fill_histo(respdf,vegasweight)

      !print *, '2 ',FintNNLO_ag_fc(2)
      !print *, '3 ',FintNNLO_ag_fc(3)

    endif

    !!-----------------------------------------------------------------------!!
    !!                               C5                                      !!
    !!-----------------------------------------------------------------------!!

    C5Lim%ids = [0,0,id_el,-id_el,-id_q,0]
    C5Lim%npart = 5

    call cut_histo(C5Lim)

    if (C5Lim%makecut) then

      kin(3) = zero
      FintNNLO_ag_fc(4) = zero

    else
       
      call res_tree_a_aq(C5Lim%AmpMom,res_tmp2)
      res_tmp1(1,1) = res_tmp2(1,1)*ndn + res_tmp2(1,2)*nup
      call get_respdf(ag_lumi,1,1,C5Lim,res_tmp1,respdf)

      call partition_nnlo_qcd(C5Lim,damp,iconf,2,2)

      s25 = C5Lim%Lim_sij(2,5)
      z5  = C5Lim%Lim_z(1)
      
      respdf = -respdf*Tr*(two/s25*Pqg(z5))*C5Lim%wgt*damp
            
      FintNNLO_ag_fc(4) = respdf(1)
      kin(3) = respdf(1)

      call fill_histo(respdf,vegasweight)
      
      !print *, '4 ',FintNNLO_ag_fc(4)
      
    endif

    !!-----------------------------------------------------------------------!!
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ag_fc)

#if(_withchecks == 1)
    FintNNLO_rr_tc_ag_fc = FintNNLO_ag_fc
#endif

    !if (ff(1).ne.zero) then
    !   print *, 'ff', ff(1)
    !   pause
    !endif

  end function xsect_nnlo_rr_5262c_ag_fc

  function xsect_nnlo_rr_5262bd_ag_fc(yRnd,ff,vegasweight,sec)
    use mod_kinematics_nnlo_tc_is_eu_abcd
    implicit none
    integer :: xsect_nnlo_rr_5262bd_ag_fc,sec
    real(dp15) :: yRnd(30), ff(1), vegasweight
    type(KinConfig) :: HardProc
    type(KinConfig) :: TCLim
    !--
    real(dp) :: xx(kNNLO_max_full)
    real(dp) :: FintNNLO_ag_fc(2),kin(2)
    real(dp) :: respdf(ipdf)
    real(dp) :: res_tmp1(1,1)
    real(dp) :: z1,z5,z6,s25,s26,s56
    real(dp) :: damp
    integer, parameter :: iconf(2,2) = reshape([1,2, 5,6],[2,2])

    xsect_nnlo_rr_5262bd_ag_fc = 0

    ff(1) = zero

    !! Technical buffers
    xx(1:kNNLO_max)=buff+onet*real(yRnd(1:kNNLO_max),dp)
    call random_number(xx(kNNLO_max_full-1))
    call random_number(xx(kNNLO_max_full))

#if (_withchecks == 1)
    if (override) then
       xx(1:kNNLO_max_full) = yRnd(1:kNNLO_max_full)
       print *, 'overriding input'
    endif
#endif

    !! Checks on buffers and if point should be failed immediately
    if (xx(xE5)*xx(xE6)*xx(xRHO5)*xx(xRHO6).lt.buff_rr) then
       failed_points = failed_points + 1
       return
    endif

    call open_histo()

    !! Generate kinematics
    call kinematics_nnlo_abcd_5i6ibd(&
         xx(1:kNNLO_max_full), 2, 1, sec, &
         HardProc,TCLim=TCLim)

    !!-----------------------------------------------------------------------!!
    !!                             Hard Process                              !!
    !!-----------------------------------------------------------------------!!

    HardProc%ids = [0,0,id_el,-id_el,id_q,-id_q]
    HardProc%npart = 6

    call cut_histo(HardProc)

    if (HardProc%makecut) then

      kin(1) = zero
      FintNNLO_ag_fc(1) = zero

    else

      call res_tree_ga_ag(HardProc%AmpMom,res_tmp1(1,1))
      call get_respdf(ag_lumi,1,1,HardProc,res_tmp1,respdf)

      call partition_nnlo_qcd(HardProc,damp,iconf,2,2)

      respdf = respdf * HardProc%wgt * damp

      FintNNLO_ag_fc(1) = respdf(1)
      kin(1) = respdf(1)
      call fill_histo(respdf,vegasweight)

      !print *, '1 ',FintNNLO_ag_fc(1)

    endif

    !!-----------------------------------------------------------------------!!
    !!                                TC                                     !!
    !!-----------------------------------------------------------------------!!

    TCLim%ids = [0,0,id_el,-id_el,0,0]
    TCLim%npart = 4

    call cut_histo(TCLim)

    if (TCLim%makecut) then

      kin(2) = zero
      FintNNLO_ag_fc(2) = zero

    else
       
      call res_treeAA_aa(TCLim%AmpMom,res_tmp1(1,1))
      res_tmp1 = res_tmp1 * (Qup2*nup + Qdn2*ndn) !-- qqb/qbq taken care in Pgqqb

      call get_respdf(ag_lumi,1,1,TCLim,res_tmp1,respdf)

      z5 = TCLim%Lim_z(1)
      z6 = TCLim%Lim_z(2)
      z1 = TCLim%Lim_z(3)
      s56 = TCLim%Lim_sij(5,6)
      s25 = TCLim%Lim_sij(2,5)
      s26 = TCLim%Lim_sij(2,6)

      !-- damp = one
      respdf = -Tr*(Pgqqb_spav_ab(-s25,-s26,s56,z1,z5,z6))*TCLim%wgt*respdf

      FintNNLO_ag_fc(2) = respdf(1)
      kin(2) = respdf(1)
      
      call fill_histo(respdf,vegasweight)

      !print *, '2 ',FintNNLO_ag_fc(2)
      
    endif

    !!-----------------------------------------------------------------------!!
    
    ff(1) = sum(kin)
    call close_histo()

    call check_ff(ff,xx,FintNNLO_ag_fc)

#if(_withchecks == 1)
    FintNNLO_rr_tc_ag_fc(1:2) = FintNNLO_ag_fc
    FintNNLO_rr_tc_ag_fc(3:4) = zero
#endif

    !if (ff(1).ne.zero) then
    !   print *, 'ff', ff(1)
    !   pause
    !endif

  end function xsect_nnlo_rr_5262bd_ag_fc  

end module mod_xsects_nnlo_rr_ag_fc


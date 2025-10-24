module mod_check_amp
  use mod_types
  use mod_consts_dp
  use mod_proc_parms
  use mod_auxfunctions
  use mod_amplitudes_tree_ppll

  public :: check_amp_lo
!  public :: check_amp_nloqcd_r
!  public :: check_amp_nloewk_r

contains

  subroutine check_amp()
    integer            :: ipoint, i,j
    real(dp)           :: p(1:4,6,20),p_lo(4,4),p_nlo(4,5)
    real(dp)           :: mgamps(1:20,-2:2,-2:2),ouramps(-2:2,-2:2)
    logical            :: amps_agree,all_amps_agree
    real(dp)           :: tol=1E-6_dp
    

    mgamps = zero   ! initialize MG amps to zero
    
    if (corr .eq. 'lo') then
       call import_MGresults_lo(p,mgamps)
    elseif (corr .eq. 'nloqcd') then
       call import_MGresults_nloqcd(p,mgamps)
    endif
    
     amps_agree = .true.
     all_amps_agree = .true.

     do ipoint = 1,20
        print *, "Checking point #",ipoint

        if (corr .eq. 'lo') then
           p_lo(:,1:4) = p(:,1:4,ipoint)
           call get_our_amps_lo(p_lo,ouramps)
        elseif (corr .eq. 'nloqcd') then
           p_nlo(:,1:5) = p(:,1:5,ipoint)
           call get_our_amps_nloqcd(p_nlo,ouramps)
        endif

        
        do i = -2,2
           do j = -2,2
!             print *, i,j!,mgamps(ipoint,i,j),ouramps(i,j)
              call compare_amps(mgamps(ipoint,i,j),ouramps(i,j),tol,amps_agree)
              if (.not. amps_agree) all_amps_agree = .false.
           enddo
        enddo
     enddo

     if (all_amps_agree) then
        print *, " "
        print *, "All amplitudes agree within a relative difference of ", tol, " for all phase space points"
	print *, " "
     else
        print *, " "
        print *, "Ampitudes do not all agree"
     endif

     stop


   end subroutine check_amp




   ! *************** functions to import MG results *******************

   subroutine import_MGresults_lo(p,mgamps)
     real(dp), intent(out)           :: p(1:4,6,20),mgamps(1:20,-2:2,-2:2)
     real(dp)          :: uux_emep(20), ddx_emep(20), uxu_emep(20), dxd_emep(20)
     real(dp)           :: dux_emvex(20), uxd_emvex(20)
     real(dp)           :: udx_veep(20), dxu_veep(20)

     
#if (_Vcharge == 0) 
#include "./MG_output/MG_qqb_emep.out"
#include "./MG_output/MG_qbq_emep.out"
    mgamps(:,+2,-2) = uux_emep(:)
    mgamps(:,+1,-1) = ddx_emep(:)
    mgamps(:,-2,+2) = uxu_emep(:)
    mgamps(:,-1,+1) = dxd_emep(:)

!
#elif (_Vcharge == 1)
#include "./MG_output/MG_qbq_veep.out"
#include "./MG_output/MG_qqb_veep.out"
    mgamps(:,-1,+2) = dxu_veep(:)
    mgamps(:,+2,-1) = udx_veep(:)
     !
#elif (_Vcharge == -1)
#include "./MG_output/MG_qqb_emvx.out"
#include "./MG_output/MG_qbq_emvx.out"
    mgamps(:,+1,-2) = dux_emvex(:)
    mgamps(:,-2,+1) = uxd_emvex(:)
#endif

  end subroutine import_MGresults_lo


  
  subroutine import_MGresults_nloqcd(p,mgamps)
     real(dp), intent(out)           :: p(1:4,6,20),mgamps(1:20,-2:2,-2:2)
     real(dp)           :: uux_emepg(20), ddx_emepg(20), uxu_emepg(20), dxd_emepg(20)
     real(dp)           :: gux_emepux(20), gdx_emepdx(20), gu_emepu(20), gd_emepd(20)
     real(dp)           :: uxg_emepux(20), dxg_emepdx(20), ug_emepu(20), dg_emepd(20)
     real(dp)           :: dux_emvexg(20), uxd_emvexg(20)
     real(dp)           :: dg_emvexu(20), gd_emvexu(20), gux_emvexdx(20), uxg_emvexdx(20)
     real(dp)           :: udx_veepg(20), dxu_veepg(20)
     real(dp)           :: ug_veepd(20), gu_veepd(20), gdx_veepux(20), dxg_veepux(20)

     
#if (_Vcharge == 0) 
#include "./MG_output/MG_qqb_emep_g.out"
#include "./MG_output/MG_qbq_emep_g.out"
#include "./MG_output/MG_qg_emep_q.out"
#include "./MG_output/MG_qbg_emep_qb.out"
#include "./MG_output/MG_gqb_emep_qb.out"
#include "./MG_output/MG_gq_emep_q.out"
    mgamps(:,+2,-2) = uux_emepg(:)
    mgamps(:,+1,-1) = ddx_emepg(:)
    mgamps(:,-2,+2) = uxu_emepg(:)
    mgamps(:,-1,+1) = dxd_emepg(:)
    mgamps(:,+2,0) = ug_emepu(:)
    mgamps(:,+1,0) = dg_emepd(:)
    mgamps(:,-2,0) = uxg_emepux(:)
    mgamps(:,-1,0) = dxg_emepdx(:)
    mgamps(:,0,-2) = gux_emepux(:)
    mgamps(:,0,-1) = gdx_emepdx(:)
    mgamps(:,0,+2) = gu_emepu(:)
    mgamps(:,0,+1) = gd_emepd(:)


#elif (_Vcharge == 1)
#include "./MG_output/MG_qqb_veep_g.out"
#include "./MG_output/MG_qbq_veep_g.out"
#include "./MG_output/MG_qg_veep_q.out"
#include "./MG_output/MG_qbg_veep_qb.out"
#include "./MG_output/MG_gq_veep_q.out"
#include "./MG_output/MG_gqb_veep_qb.out"

    mgamps(:,+2,-1) = udx_veepg(:)
    mgamps(:,-1,+2) = dxu_veepg(:)
    mgamps(:,0,-1) = gdx_veepux(:)
    mgamps(:,-1,0) = dxg_veepux(:)
    mgamps(:,0,+2) = gu_veepd(:)
    mgamps(:,+2,0) = ug_veepd(:)


#elif (_Vcharge == -1)
#include "./MG_output/MG_qqb_emvx_g.out"
#include "./MG_output/MG_qbq_emvx_g.out"
#include "./MG_output/MG_qg_emvx_q.out"
#include "./MG_output/MG_qbg_emvx_qb.out"
#include "./MG_output/MG_gq_emvx_q.out"
#include "./MG_output/MG_gqb_emvx_qb.out"
    
    mgamps(:,+1,-2) = dux_emvexg(:)
    mgamps(:,-2,+1) = uxd_emvexg(:)
    mgamps(:,0,-2) = gux_emvexdx(:)
    mgamps(:,-2,0) = uxg_emvexdx(:)
    mgamps(:,0,+1) = gd_emvexu(:)
    mgamps(:,+1,0) = dg_emvexu(:)
  
    
#endif

  end subroutine import_MGresults_nloqcd



   ! ************* functions to get our amplitudes *************

   subroutine get_our_amps_lo(p,amp)
     real(dp), intent(in)    :: p(4,4)
     real(dp), intent(out)   :: amp(-2:2,-2:2)
     real(dp)                :: res(-5:7,-5:7),res1(2,2)


     call res_tree_qqb_gen(p,res)
     amp(-2:2,-2:2) = res(-2:2,-2:2)

   end subroutine get_our_amps_lo

   
   subroutine get_our_amps_nloqcd(p,amp)
     real(dp), intent(in)    :: p(4,5)
     real(dp), intent(out)   :: amp(-2:2,-2:2)
     real(dp)                :: res_qqb(-5:7,-5:7), res_gq(-5:7,-5:7), res_qg(-5:7,-5:7)

     amp = zero
     call res_tree_g_qqb_gen(p,res_qqb)
     call res_tree_g_gq_gen(p,res_gq)
     call res_tree_g_qg_gen(p,res_qg)
     amp(-2:2,-2:2) = res_qqb(-2:2,-2:2) + res_gq(-2:2,-2:2) + res_qg(-2:2,-2:2)

     amp = amp * (0.118_dp)*four*pi   ! gs^2

   end subroutine get_our_amps_nloqcd

     
! ************* FUNCTIONS "COMPARE_AMPS" *************     
   subroutine compare_amps(amp1,amp2,tol,agree)
     real(dp), intent(in)     :: amp1,amp2,tol
     logical, intent(out)     :: agree

     if (amp1 .ne. amp1 .or. amp2 .ne. amp2) then
        print *, "amplitude is nan" ,amp1,amp2
	agree=.false.
        return
     endif

     if (abs(amp1) .lt. 1E-30_dp ) then    ! its zero                                                                                                                                                                                                                           
	if (abs(amp2) .gt. 1E-30_dp) then
           agree=.false.
           print *, "Amplitudes do not agree"
           print *, "mgamps=",amp1
           print *, "ouramps=",amp2
	else
           agree=.true.
	endif
     else
        if (abs((amp1-amp2)/amp1) .gt. tol) then
           agree=.false.
           print *, "Amplitudes do not agree"
           print *, "mgamp=",amp1
           print *, "ouramp=",amp2
           print *, "ratio=",amp1/amp2
	else
           agree=.true.
           print *, "Amplitudes agree!", amp1,amp2
	endif
     endif

     
   end subroutine compare_amps


 end module mod_check_amp



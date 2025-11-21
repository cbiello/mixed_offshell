module mod_check_amp
  use mod_types
  use mod_consts_dp
  use mod_proc_parms
  use mod_parms
  use mod_auxfunctions
  use mod_amplitudes_tree_ppll

  public :: check_amp

contains

  subroutine check_amp()
    integer            :: ipoint, i,j
    real(dp)           :: p(1:4,6,20),p_lo(4,4),p_nlo(4,5),p_qcdew(4,6)
    real(dp)           :: mgamps(1:20,-2:3,-2:3),ouramps(-2:3,-2:3)
    logical            :: amps_agree,all_amps_agree
    real(dp)           :: tol=1E-6_dp
    

    mgamps = zero   ! initialize MG amps to zero
    
    if (corr .eq. 'lo') then
       call import_MGresults_lo(p,mgamps)
    elseif (corr .eq. 'nloqcd') then
       call import_MGresults_nloqcd(p,mgamps)
    elseif (corr .eq. 'nloew') then
       call import_MGresults_nloew(p,mgamps)
    elseif (corr .eq. 'qcdew') then
       call import_MGresults_qcdew(p,mgamps)
       
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
        elseif (corr .eq. 'nloew') then
           p_nlo(:,1:5) = p(:,1:5,ipoint)
           call get_our_amps_nloew(p_nlo,ouramps)
        elseif (corr .eq. 'qcdew') then
           p_qcdew(:,1:6) = p(:,1:6,ipoint)
           call get_our_amps_qcdew(p_qcdew,ouramps)
        endif

        
        do i = -2,3
           do j = -2,3
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
     real(dp), intent(out)           :: p(1:4,6,20),mgamps(1:20,-2:3,-2:3)
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
     real(dp), intent(out)           :: p(1:4,6,20),mgamps(1:20,-2:3,-2:3)
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



  subroutine import_MGresults_nloew(p,mgamps)
     real(dp), intent(out)           :: p(1:4,6,20),mgamps(1:20,-2:3,-2:3)
     real(dp)           :: uux_emepa(20), ddx_emepa(20), uxu_emepa(20), dxd_emepa(20)
     real(dp)           :: aux_emepux(20), adx_emepdx(20), au_emepu(20), ad_emepd(20)
     real(dp)           :: uxa_emepux(20), dxa_emepdx(20), ua_emepu(20), da_emepd(20)
     real(dp)           :: dux_emvexa(20), uxd_emvexa(20)
     real(dp)           :: da_emvexu(20), ad_emvexu(20), aux_emvexdx(20), uxa_emvexdx(20)
     real(dp)           :: udx_veepa(20), dxu_veepa(20)
     real(dp)           :: ua_veepd(20), au_veepd(20), adx_veepux(20), dxa_veepux(20)

     
#if (_Vcharge == 0) 
#include "./MG_output/MG_qqb_emep_a.out"
#include "./MG_output/MG_qbq_emep_a.out"
#include "./MG_output/MG_qa_emep_q.out"
#include "./MG_output/MG_aq_emep_q.out"
#include "./MG_output/MG_qba_emep_qb.out"
#include "./MG_output/MG_aqb_emep_qb.out"
     
    mgamps(:,+2,-2) = uux_emepa(:)
    mgamps(:,+1,-1) = ddx_emepa(:)
    mgamps(:,-2,+2) = uxu_emepa(:)
    mgamps(:,-1,+1) = dxd_emepa(:)
    mgamps(:,+2,3) = ua_emepu(:)   ! for checks, put photon in position 3
    mgamps(:,+1,3) = da_emepd(:)
    mgamps(:,-2,3) = uxa_emepux(:)  
    mgamps(:,-1,3) = dxa_emepdx(:)
    
    mgamps(:,3,+2) = au_emepu(:)   
    mgamps(:,3,+1) = ad_emepd(:)
    mgamps(:,3,-2) = aux_emepux(:) 
    mgamps(:,3,-1) = adx_emepdx(:)


#elif (_Vcharge == 1)
#include "./MG_output/MG_qqb_veep_a.out"
#include "./MG_output/MG_qbq_veep_a.out"
#include "./MG_output/MG_qa_veep_q.out"
#include "./MG_output/MG_qba_veep_qb.out"
#include "./MG_output/MG_aq_veep_q.out"
#include "./MG_output/MG_aqb_veep_qb.out"

    mgamps(:,+2,-1) = udx_veepa(:)
    mgamps(:,-1,+2) = dxu_veepa(:)
    mgamps(:,3,-1) = adx_veepux(:)
    mgamps(:,-1,3) = dxa_veepux(:)
    mgamps(:,3,+2) = au_veepd(:)
    mgamps(:,+2,3) = ua_veepd(:)


#elif (_Vcharge == -1)
#include "./MG_output/MG_qqb_emvx_a.out"
#include "./MG_output/MG_qbq_emvx_a.out"
#include "./MG_output/MG_qa_emvx_q.out"
#include "./MG_output/MG_qba_emvx_qb.out"
#include "./MG_output/MG_aq_emvx_q.out"
#include "./MG_output/MG_aqb_emvx_qb.out"
    
    mgamps(:,+1,-2) = dux_emvexa(:)
    mgamps(:,-2,+1) = uxd_emvexa(:)
    mgamps(:,3,-2) = aux_emvexdx(:)
    mgamps(:,-2,3) = uxa_emvexdx(:)
    mgamps(:,3,+1) = ad_emvexu(:)
    mgamps(:,+1,3) = da_emvexu(:)
  
    
#endif

  end subroutine import_MGresults_nloew






  subroutine import_MGresults_qcdew(p,mgamps)
     real(dp), intent(out)           :: p(1:4,6,20),mgamps(1:20,-2:3,-2:3)
     real(dp)           :: uux_emepga(20), ddx_emepga(20), uxu_emepga(20), dxd_emepga(20)
     real(dp)           :: gux_emepuxa(20), gdx_emepdxa(20), gu_emepua(20), gd_emepda(20)
     real(dp)           :: uxg_emepuxa(20), dxg_emepdxa(20), ug_emepua(20), dg_emepda(20)
     real(dp)           :: aux_emepgux(20), adx_emepgdx(20), au_emepgu(20), ad_emepgd(20)
     real(dp)           :: uxa_emepgux(20), dxa_emepgdx(20), ua_emepgu(20), da_emepgd(20)
     real(dp)           :: ag_emepuux(20),ag_emepddx(20),ga_emepuux(20),ga_emepddx(20)
     !
     real(dp)           :: udx_veepga(20), dxu_veepga(20),adx_veepgux(20),dxa_veepgux(20),au_veepgd(20),ua_veepgd(20)
     real(dp)           :: gdx_veepuxa(20),dxg_veepuxa(20), gu_veepda(20),ug_veepda(20),ga_veepdux(20),ag_veepdux(20)
     !
     real(dp)           :: dux_emvexga(20), uxd_emvexga(20),aux_emvexgdx(20),uxa_emvexgdx(20),ad_emvexgu(20),da_emvexgu(20)
     real(dp)           :: gux_emvexdxa(20),uxg_emvexdxa(20),gd_emvexua(20),dg_emvexua(20),ga_emvexudx(20),ag_emvexudx(20)
     

     
#if (_Vcharge == 0) 
#include "./MG_output/MG_qqb_emep_ga.out"
#include "./MG_output/MG_qbq_emep_ga.out"

#include "./MG_output/MG_gq_emep_qa.out"
#include "./MG_output/MG_qg_emep_qa.out"
#include "./MG_output/MG_gqb_emep_qba.out"
#include "./MG_output/MG_qbg_emep_qba.out"

#include "./MG_output/MG_aq_emep_gq.out"
#include "./MG_output/MG_qba_emep_gqb.out"
#include "./MG_output/MG_qa_emep_gq.out"
#include "./MG_output/MG_aqb_emep_gqb.out"

#include "./MG_output/MG_ag_emep_qqb.out"
#include "./MG_output/MG_ga_emep_qqb.out"
     
    mgamps(:,+2,-2) = uux_emepga(:)
    mgamps(:,+1,-1) = ddx_emepga(:)
    mgamps(:,-2,+2) = uxu_emepga(:)
    mgamps(:,-1,+1) = dxd_emepga(:)
    mgamps(:,+2,0) = ug_emepua(:)
    mgamps(:,+1,0) = dg_emepda(:)
    mgamps(:,0,+2) = gu_emepua(:)
    mgamps(:,0,+1) = gd_emepda(:)
    mgamps(:,-2,0) = uxg_emepuxa(:)
    mgamps(:,-1,0) = dxg_emepdxa(:)
    mgamps(:,0,-2) = gux_emepuxa(:)
    mgamps(:,0,-1) = gdx_emepdxa(:)

    mgamps(:,-2,3) = uxa_emepgux(:)  
    mgamps(:,-1,3) = dxa_emepgdx(:)
    mgamps(:,+2,3) = ua_emepgu(:)  
    mgamps(:,+1,3) = da_emepgd(:)
!    
    mgamps(:,3,+2) = au_emepgu(:)   
    mgamps(:,3,+1) = ad_emepgd(:)
    mgamps(:,3,-2) = aux_emepgux(:) 
    mgamps(:,3,-1) = adx_emepgdx(:)

    mgamps(:,0,3) = ga_emepuux(:)*nup + ga_emepddx(:)*ndn
    mgamps(:,3,0) = ag_emepuux(:)*nup + ag_emepddx(:)*ndn

    !
!
#elif (_Vcharge == 1)
#include "./MG_output/MG_qqb_veep_ga.out"
#include "./MG_output/MG_qbq_veep_ga.out"

#include "./MG_output/MG_gq_veep_qa.out"
#include "./MG_output/MG_qg_veep_qa.out"
#include "./MG_output/MG_gqb_veep_qba.out"
#include "./MG_output/MG_qbg_veep_qba.out"

#include "./MG_output/MG_aq_veep_gq.out"
#include "./MG_output/MG_qba_veep_gqb.out"
#include "./MG_output/MG_qa_veep_gq.out"
#include "./MG_output/MG_aqb_veep_gqb.out"

#include "./MG_output/MG_ag_veep_qqb.out"
#include "./MG_output/MG_ga_veep_qqb.out"

!
    mgamps(:,+2,-1) = udx_veepga(:)
    mgamps(:,-1,+2) = dxu_veepga(:)
    mgamps(:,3,-1) =  adx_veepgux(:)
    mgamps(:,-1,3) =  dxa_veepgux(:)
    mgamps(:,3,+2) =  au_veepgd(:)
    mgamps(:,+2,3) =  ua_veepgd(:)

    mgamps(:,0,-1) =  gdx_veepuxa(:)
    mgamps(:,-1,0) =  dxg_veepuxa(:)
    mgamps(:,0,+2) =  gu_veepda(:)
    mgamps(:,+2,0) =  ug_veepda(:)

    mgamps(:,0,3) = ga_veepdux(:)*nup
    mgamps(:,3,0) = ag_veepdux(:)*nup
!
!
#elif (_Vcharge == -1)

#include "./MG_output/MG_qqb_emvx_ga.out"
#include "./MG_output/MG_qbq_emvx_ga.out"

#include "./MG_output/MG_gq_emvx_qa.out"
#include "./MG_output/MG_qg_emvx_qa.out"
#include "./MG_output/MG_gqb_emvx_qba.out"
#include "./MG_output/MG_qbg_emvx_qba.out"

#include "./MG_output/MG_aq_emvx_gq.out"
#include "./MG_output/MG_qba_emvx_gqb.out"
#include "./MG_output/MG_qa_emvx_gq.out"
#include "./MG_output/MG_aqb_emvx_gqb.out"

#include "./MG_output/MG_ag_emvx_qqb.out"
#include "./MG_output/MG_ga_emvx_qqb.out"
!    
    mgamps(:,+1,-2) = dux_emvexga(:)
    mgamps(:,-2,+1) = uxd_emvexga(:)
    mgamps(:,3,-2)  = aux_emvexgdx(:)
    mgamps(:,-2,3)  = uxa_emvexgdx(:)
    mgamps(:,3,+1)  = ad_emvexgu(:)
    mgamps(:,+1,3)  = da_emvexgu(:)

    mgamps(:,0,-2)  = gux_emvexdxa(:)
    mgamps(:,-2,0)  = uxg_emvexdxa(:)
    mgamps(:,0,+1)  = gd_emvexua(:)
    mgamps(:,+1,0)  = dg_emvexua(:)

    mgamps(:,0,3) = ga_emvexudx(:)*nup
    mgamps(:,3,0) = ag_emvexudx(:)*nup
!  
!    
#endif

  end subroutine import_MGresults_qcdew
  



   ! ************* functions to get our amplitudes *************

   subroutine get_our_amps_lo(p,amp)
     real(dp), intent(in)    :: p(4,4)
     real(dp), intent(out)   :: amp(-2:3,-2:3)
     real(dp)                :: res(-5:7,-5:7),res1(2,2)


     call res_tree_qqb_gen(p,res)
     amp(-2:2,-2:2) = res(-2:2,-2:2)

   end subroutine get_our_amps_lo

   
   subroutine get_our_amps_nloqcd(p,amp)
     real(dp), intent(in)    :: p(4,5)
     real(dp), intent(out)   :: amp(-2:3,-2:3)
     real(dp)                :: res_qqb(-5:7,-5:7), res_gq(-5:7,-5:7), res_qg(-5:7,-5:7)

     amp = zero
     call res_tree_g_qqb_gen(p,res_qqb)
     call res_tree_g_gq_gen(p,res_gq)
     call res_tree_g_qg_gen(p,res_qg)
     amp(-2:2,-2:2) = res_qqb(-2:2,-2:2) + res_gq(-2:2,-2:2) + res_qg(-2:2,-2:2)

     amp = amp * (0.118_dp)*four*pi   ! gs^2

   end subroutine get_our_amps_nloqcd



   subroutine get_our_amps_nloew(p,amp)
     real(dp), intent(in)    :: p(4,5)
     real(dp), intent(out)   :: amp(-2:3,-2:3)
     real(dp)                :: res_qqb(-5:7,-5:7), res_aq(-5:7,-5:7), res_qa(-5:7,-5:7)

     amp = zero
     call res_tree_a_qqb_gen(p,res_qqb)
     call res_tree_a_aq_gen(p,res_aq)
     call res_tree_a_qa_gen(p,res_qa)
     amp(-2:2,-2:2) = res_qqb(-2:2,-2:2)
     amp(-2:2,3) = res_qa(-2:2,7)      ! for the checks, put the photon in position 3
     amp(3,-2:2) = res_aq(7,-2:2)


     amp = amp * eesq

   end subroutine get_our_amps_nloew



   subroutine get_our_amps_qcdew(p,amp)
     real(dp), intent(in)    :: p(4,6)
     real(dp), intent(out)   :: amp(-2:3,-2:3)
     real(dp)                :: res_qqb(-5:7,-5:7), res_aq(-5:7,-5:7), res_qa(-5:7,-5:7), res_gq(-5:7,-5:7), res_qg(-5:7,-5:7), res_ga(-5:7,-5:7),  res_ag(-5:7,-5:7)
     

!     amp = zero
     call res_tree_ga_qqb_gen(p,res_qqb)
     call res_tree_ga_aq_gen(p,res_aq)
     call res_tree_ga_qa_gen(p,res_qa)
     call res_tree_ga_gq_gen(p,res_gq)
     call res_tree_ga_qg_gen(p,res_qg)
     call res_tree_ga_ga_gen(p,res_ga)
     call res_tree_ga_ag_gen(p,res_ag)
     amp(-2:2,-2:2) = res_qqb(-2:2,-2:2) + res_gq(-2:2,-2:2) + res_qg(-2:2,-2:2)
     amp(-2:2,3) = res_qa(-2:2,7) + res_ga(-2:2,7)      ! for the checks, put the photon in position 3
     amp(3,-2:2) = res_aq(7,-2:2) + res_ag(7,-2:2)      ! for the checks, put the photon in position 3

!

     amp = amp * eesq * (0.118_dp)*four*pi   ! gs^2

   end subroutine get_our_amps_qcdew


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



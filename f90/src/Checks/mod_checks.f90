module mod_checks
  !-- container for checking routines
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  use mod_auxfunctions
  implicit none
  private

  public :: do_checks
  
contains

  subroutine do_checks()
    integer :: icheck

    print *, 'Doing checks'

    icheck = 6
    
    if (icheck == 1) then 
       call check_rv_amps_all()
    elseif (icheck == 2) then
       call check_rr_ag_amp()
    elseif (icheck == 3) then
       call check_quark_amp()
    elseif (icheck == 4) then
       call check_aa_amp()
    elseif (icheck == 5) then
       call check_my_rv_amps()
    elseif (icheck == 6) then
       call check_my_nf_amps()       
    endif
    
    stop
    
  end subroutine do_checks

  !--

  subroutine check_my_nf_amps()
    use mod_amplitudes_tree_ppll
    use mod_amplitudes_twol_nf_ppll
    real(dp) :: p4(4,4)
    real(dp) :: res0(2,2),res0check(2,2),res2(2,2)

    p4(:,1) = [157.63994801951441_dp,  0.0000000000000000_dp, 0.0000000000000000_dp, 157.63994801951441_dp]
    p4(:,2) = [157.63994801951441_dp,  0.0000000000000000_dp, 0.0000000000000000_dp,-157.63994801951441_dp]
    p4(:,3) = [157.63994801951441_dp,  63.055694515696473_dp, 60.751183318642511_dp, 131.08633157626718_dp]
    p4(:,4) = [157.63994801951441_dp, -63.055694515696473_dp,-60.751183318642511_dp,-131.08633157626718_dp]
    
    call res_tree_qqb(p4,res0)
    call res_twol_nf_qqb(p4,res0check,res2)
    print *, 'res0     ',res0
    print *, 'res0check',res0check
    
  end subroutine check_my_nf_amps
    
  subroutine check_my_rv_amps()
    use mod_ol_interface
    use mod_amplitudes_tree_ppll
    use mod_amplitudes_loop_ppll
    real(dp) :: p5(4,5),p4(4,4)
    real(dp) :: mu
    real(dp) :: res0(2,2),res1(2,2)
    real(dp) :: myres0(2,2),myres1(2,2)

    mu = 200.0_dp

    p5(:,1) = [157.63994801951441_dp,  0.0000000000000000_dp, 0.0000000000000000_dp, 157.63994801951441_dp]
    p5(:,2) = [157.63994801951441_dp,  0.0000000000000000_dp, 0.0000000000000000_dp,-157.63994801951441_dp]
    p5(:,3) = [152.53965460656173_dp,  61.015586361462930_dp, 58.785635473683840_dp, 126.84515564423882_dp]
    p5(:,4) = [156.36208303788214_dp, -60.266881434355547_dp,-64.767855301256915_dp,-128.92683558436670_dp]
    p5(:,5) = [6.3781583945849265_dp,-0.74870492710738068_dp, 5.9822198275730800_dp, 2.0816799401278798_dp]

    p4(:,1) = [157.63994801951441_dp,  0.0000000000000000_dp, 0.0000000000000000_dp, 157.63994801951441_dp]
    p4(:,2) = [157.63994801951441_dp,  0.0000000000000000_dp, 0.0000000000000000_dp,-157.63994801951441_dp]
    p4(:,3) = [157.63994801951441_dp,  63.055694515696473_dp, 60.751183318642511_dp, 131.08633157626718_dp]
    p4(:,4) = [157.63994801951441_dp, -63.055694515696473_dp,-60.751183318642511_dp,-131.08633157626718_dp]

    corr = 'nnlo'
   
    !-- QCD loops
    ch   = 'ns'
    sec = 'rvqcd'
    call initialise_ol(0)
    call ol_res_qcdloop_qqb(p4,res0,res1)
    call res_qcdloop_qqb(p4,myres0,myres1)
    call printout('qqb',res0,res1,myres0,myres1)

    !-- check 5 point ns
    call res_tree_a_qqb(p5,myres0)
    call ol_res_qcdloop_a_qqb(p5,res0,res1)
    call res_qcdloop_a_qqb(p5,myres0,myres1)
    call printout('a_qqb',res0,res1,myres0,myres1)

    !-- check 5 point aq
    ch   = 'aq'
    call initialise_ol(0)
    call ol_res_qcdloop_a_aq(p5,res0,res1)
    call res_qcdloop_a_aq(p5,myres0,myres1)
    call printout('aq',res0,res1,myres0,myres1)

    !-- check 5 point qa
    ch   = 'qa'
    call initialise_ol(0)
    call ol_res_qcdloop_a_qa(p5,res0,res1)
    call res_qcdloop_a_qa(p5,myres0,myres1)
    call printout('qa',res0,res1,myres0,myres1)

  contains

    subroutine printout(string,res0,res1,myres0,myres1)
      character(*), intent(in) :: string
      real(dp), intent(in) :: res0(:,:),res1(:,:),myres0(:,:),myres1(:,:)

      print *, string
      print *, 'OL 0', res0
      print *, 'me 0', myres0
      print *, '    ', one-res0/myres0
      print *, ''
      print *, 'OL 1', res1
      print *, 'me 1', myres1
      print *, '    ', one-res1/myres1
      print *, ''
    end subroutine printout
    
  end subroutine check_my_rv_amps
    
  !--

  subroutine check_aa_amp()
    use mod_amplitudes_tree_ppll
    real(dp) :: p5(4,5)
    real(dp) :: res

    !-- for now, madgraph point
    p5(:,1) = (/0.5000000E+03_dp,  0.0000000E+00_dp,  0.0000000E+00_dp,  0.5000000E+03_dp/)
    p5(:,2) = (/0.5000000E+03_dp,  0.0000000E+00_dp,  0.0000000E+00_dp, -0.5000000E+03_dp/)
    p5(:,3) = (/0.4585788E+03_dp,  0.1694532E+03_dp,  0.3796537E+03_dp, -0.1935025E+03_dp/)
    p5(:,4) = (/0.3640666E+03_dp, -0.1832987E+02_dp, -0.3477043E+03_dp,  0.1063496E+03_dp/)
    p5(:,5) = (/0.1773546E+03_dp, -0.1511234E+03_dp, -0.3194936E+02_dp,  0.8715287E+02_dp/)

    call res_treeAA_a_aa(p5,res)
    print *, 'res', res * eesq
    stop
  end subroutine check_aa_amp
    
  subroutine check_rv_amps_all()
    use mod_ol_interface
    use mod_amplitudes_tree_ppll
    use mod_amplitudes_loop_ppll
    real(dp) :: p5(4,5),p4(4,4)
    real(dp) :: mu
    real(dp) :: res0(2,2),res1(2,2),res0b(2,3),res1b(2,3)
    real(dp) :: myres0(2,2)

    mu = 200.0_dp

    p5(:,1) = [157.63994801951441_dp,  0.0000000000000000_dp, 0.0000000000000000_dp, 157.63994801951441_dp]
    p5(:,2) = [157.63994801951441_dp,  0.0000000000000000_dp, 0.0000000000000000_dp,-157.63994801951441_dp]
    p5(:,3) = [152.53965460656173_dp,  61.015586361462930_dp, 58.785635473683840_dp, 126.84515564423882_dp]
    p5(:,4) = [156.36208303788214_dp, -60.266881434355547_dp,-64.767855301256915_dp,-128.92683558436670_dp]
    p5(:,5) = [6.3781583945849265_dp,-0.74870492710738068_dp, 5.9822198275730800_dp, 2.0816799401278798_dp]

    p4(:,1) = [157.63994801951441_dp,  0.0000000000000000_dp, 0.0000000000000000_dp, 157.63994801951441_dp]
    p4(:,2) = [157.63994801951441_dp,  0.0000000000000000_dp, 0.0000000000000000_dp,-157.63994801951441_dp]
    p4(:,3) = [157.63994801951441_dp,  63.055694515696473_dp, 60.751183318642511_dp, 131.08633157626718_dp]
    p4(:,4) = [157.63994801951441_dp, -63.055694515696473_dp,-60.751183318642511_dp,-131.08633157626718_dp]

    corr = 'nnlo'
    ch   = 'ns'

    !-- EW loops
    sec = 'rv_g'
    call initialise_ol(0)

    call res_ewkloop_qqb(p4,res0b,res1b)
    call res_tree_qqb(p4,myres0)
    call printout('ew-ns',res0b,res1b,myres0)

    call res_ewkloop_g_qqb(p5,res0b,res1b)
    call res_tree_g_qqb(p5,myres0)
    call printout('ew-ns-g',res0b,res1b,myres0)
    
    !-- QCD loops
    sec = 'rv53_a'
    call initialise_ol(0)
    call res_qcdloop_qqb(p4,res0,res1)
    call res_tree_qqb(p4,myres0)
    call printout('qcd-ns',res0,res1,myres0)

    call res_qcdloop_a_qqb(p5,res0,res1)
    call res_tree_a_qqb(p5,myres0)
    call printout('qcd-ns-a',res0,res1,myres0)

    !-- gq and qg
    corr = 'nnlo'
    ch   = 'gq'
    call initialise_ol(0)
    call res_ewkloop_g_gq(p5,res0b,res1b)
    call res_tree_g_gq(p5,myres0)
    call printout('ew-gq',res0b,res1b,myres0)

    ch   = 'qg'
    call initialise_ol(0)
    call res_ewkloop_g_qg(p5,res0b,res1b)
    call res_tree_g_qg(p5,myres0)
    call printout('ew-qg',res0b,res1b,myres0)

    !-- aq and qa
    corr = 'nnlo'
    ch   = 'aq'
    call initialise_ol(0)
    call res_qcdloop_a_aq(p5,res0,res1)
    call res_tree_a_aq(p5,myres0)
    call printout('qcd-aq',res0,res1,myres0)

    ch   = 'qa'
    call initialise_ol(0)
    call res_qcdloop_a_qa(p5,res0,res1)
    call res_tree_a_qa(p5,myres0)
    call printout('qcd-qa',res0,res1,myres0)

  contains
    
    subroutine printout(string,res0_ol,res1_ol,res0_us)
      character(*), intent(in) :: string
      real(dp), intent(in) :: res0_ol(:,:),res1_ol(:,:),res0_us(:,:)
      character(30) :: proc(2,3)
      integer  :: len,i,j
      real(dp) :: res0us

      len = size(res0_ol,2)
      
      if (string.eq.'qcd-ns') then
         proc(1,1) = 'd db -> e+ e- qcd';proc(1,2) = 'u ub -> e+ e- qcd'
         proc(2,1) = 'db d -> e+ e- qcd';proc(2,2) = 'ub u -> e+ e- qcd'
      elseif (string.eq.'qcd-ns-a') then
         proc(1,1) = 'd db -> e+ e- γ qcd';proc(1,2) = 'u ub -> e+ e- γ qcd'
         proc(2,1) = 'db d -> e+ e- γ qcd';proc(2,2) = 'ub u -> e+ e- γ qcd'
      elseif (string.eq.'ew-ns') then
         proc(1,1) = 'd db -> e+ e- ew';proc(1,2) = 'u ub -> e+ e- ew'
         proc(2,1) = 'db d -> e+ e- ew';proc(2,2) = 'ub u -> e+ e- ew'
         proc(1,3) = 'b bb -> e+ e- ew';proc(2,3) = 'bb b -> e+ e- ew'
      elseif (string.eq.'ew-ns-g') then
         proc(1,1) = 'd db -> e+ e- g ew';proc(1,2) = 'u ub -> e+ e- g ew'
         proc(2,1) = 'db d -> e+ e- g ew';proc(2,2) = 'ub u -> e+ e- g ew'
         proc(1,3) = 'b bb -> e+ e- g ew';proc(2,3) = 'bb b -> e+ e- g ew'
      elseif (string.eq.'ew-gq') then
         proc(1,1) = 'g db -> e+ e- db ew';proc(1,2) = 'g ub -> e+ e- ub ew'
         proc(2,1) = 'g d  -> e+ e- d  ew';proc(2,2) = 'g u  -> e+ e- u  ew'
         proc(1,3) = 'g bb -> e+ e- bb ew';proc(2,3) = 'g b  -> e+ e- b ew'
      elseif (string.eq.'ew-qg') then
         proc(1,1) = 'd g  -> e+ e- d  ew';proc(1,2) = 'u g  -> e+ e- u  ew'
         proc(2,1) = 'db g -> e+ e- db ew';proc(2,2) = 'ub g -> e+ e- ub ew'
         proc(1,3) = 'b g  -> e+ e- b  ew';proc(2,3) = 'bb g -> e+ e- bb ew'
      elseif (string.eq.'qcd-aq') then
         proc(1,1) = 'γ db -> e+ e- db qcd';proc(1,2) = 'γ ub -> e+ e- ub qcd'
         proc(2,1) = 'γ d  -> e+ e- d  qcd';proc(2,2) = 'γ u  -> e+ e- u  qcd'
         proc(1,3) = 'γ bb -> e+ e- bb qcd';proc(2,3) = 'γ b  -> e+ e- b qcd'
      elseif (string.eq.'qcd-qa') then
         proc(1,1) = 'd γ  -> e+ e- d  qcd';proc(1,2) = 'u γ  -> e+ e- u  qcd'
         proc(2,1) = 'db γ -> e+ e- db qcd';proc(2,2) = 'ub γ -> e+ e- ub qcd'
         proc(1,3) = 'b γ  -> e+ e- b  qcd';proc(2,3) = 'bb γ -> e+ e- bb qcd'
      endif

      print *, 'channel, tree_us, tree_ol, 1-ratio, fin_ol'
      
      do i = 1,2
         do j = 1,len
            if (j.eq.3) then
               res0us = res0_us(i,1)
            else
               res0us = res0_us(i,j)
            endif
            print *, trim(proc(i,j)), res0us,res0_ol(i,j),one-res0us/res0_ol(i,j),res1_ol(i,j)
         enddo
      enddo

      print *, ''
      
    end subroutine printout
    
  end subroutine check_rv_amps_all

  subroutine check_rr_ag_amp()
    use mod_amplitudes_tree_ppll
    real(dp) :: p6(4,6),res_jj(2,2),res_jj_ol(2,2)

    p6(:,1) = [500.00000000000000_dp, 0.0000000000000000_dp, 0.0000000000000000_dp, 500.00000000000000_dp]     
    p6(:,2) = [500.00000000000000_dp, 0.0000000000000000_dp, 0.0000000000000000_dp,-500.00000000000000_dp]     
    p6(:,3) = [254.11173663453982_dp,-20.788191863948018_dp, 130.41070543048184_dp, 217.10291035261022_dp]     
    p6(:,4) = [308.27642719975739_dp,-304.46346817638226_dp, 26.163062947104724_dp,-40.642911439539645_dp]     
    p6(:,5) = [27.039412347527474_dp, 6.7685428263938441_dp, 5.6881141966812629_dp, 25.553121237774590_dp]     
    p6(:,6) = [410.57242381817514_dp, 318.48311721393640_dp,-162.26188257426779_dp,-202.01312015084520_dp]

    !---------------------------------------------------------------------------------------------!
    print *,'=================================================================================================='
    print *,'                 Checks of tree-level six-point amplitudes with g and a emissions                 '
    print *,'=================================================================================================='
    print *,''
    print *, '                               Result                  Benchmark                   Check         '
    print *,'--------------------------------------------------------------------------------------------------'
    !----------------------------------------------------------------------------------------------------------!
    call res_tree_ga_qqb(p6,res_jj)
    call get_ol_res_tree_ga_qqb(res_jj_ol)
    print *, 'd  db -> e- e+ g a  ', res_jj(1,1), res_jj_ol(1,1), abs(one-res_jj(1,1)/res_jj_ol(1,1))
    print *, 'db d  -> e- e+ g a  ', res_jj(2,1), res_jj_ol(2,1), abs(one-res_jj(2,1)/res_jj_ol(2,1))
    print *, 'u  ub -> e- e+ g a  ', res_jj(1,2), res_jj_ol(1,2), abs(one-res_jj(1,2)/res_jj_ol(1,2))
    print *, 'ub u  -> e- e+ g a  ', res_jj(2,2), res_jj_ol(2,2), abs(one-res_jj(2,2)/res_jj_ol(2,2))
    print *,''
    !----------------------------------------------------------------------------------------------------------!
    call res_tree_ga_gq(p6,res_jj)
    call get_ol_res_tree_ga_gq(res_jj_ol)
    print *, 'g db -> e- e+ db a  ', res_jj(1,1), res_jj_ol(1,1), abs(one-res_jj(1,1)/res_jj_ol(1,1))
    print *, 'g d  -> e- e+ d  a  ', res_jj(2,1), res_jj_ol(2,1), abs(one-res_jj(2,1)/res_jj_ol(2,1))
    print *, 'g ub -> e- e+ ub a  ', res_jj(1,2), res_jj_ol(1,2), abs(one-res_jj(1,2)/res_jj_ol(1,2))
    print *, 'g u  -> e- e+ u  a  ', res_jj(2,2), res_jj_ol(2,2), abs(one-res_jj(2,2)/res_jj_ol(2,2))
    print *,''
    !----------------------------------------------------------------------------------------------------------!
    call res_tree_ga_qg(p6,res_jj)
    call get_ol_res_tree_ga_qg(res_jj_ol)
    print *, 'd  g -> e- e+ db a  ', res_jj(1,1), res_jj_ol(1,1), abs(one-res_jj(1,1)/res_jj_ol(1,1))
    print *, 'db g -> e- e+ d  a  ', res_jj(2,1), res_jj_ol(2,1), abs(one-res_jj(2,1)/res_jj_ol(2,1))
    print *, 'u  g -> e- e+ ub a  ', res_jj(1,2), res_jj_ol(1,2), abs(one-res_jj(1,2)/res_jj_ol(1,2))
    print *, 'ub g -> e- e+ u  a  ', res_jj(2,2), res_jj_ol(2,2), abs(one-res_jj(2,2)/res_jj_ol(2,2))
    print *,''
    !----------------------------------------------------------------------------------------------------------!
    call res_tree_ga_aq(p6,res_jj)
    call get_ol_res_tree_ga_aq(res_jj_ol)
    print *, 'a db -> e- e+ db g  ', res_jj(1,1), res_jj_ol(1,1), abs(one-res_jj(1,1)/res_jj_ol(1,1))
    print *, 'a d  -> e- e+ d  g  ', res_jj(2,1), res_jj_ol(2,1), abs(one-res_jj(2,1)/res_jj_ol(2,1))
    print *, 'a ub -> e- e+ ub g  ', res_jj(1,2), res_jj_ol(1,2), abs(one-res_jj(1,2)/res_jj_ol(1,2))
    print *, 'a u  -> e- e+ u  g  ', res_jj(2,2), res_jj_ol(2,2), abs(one-res_jj(2,2)/res_jj_ol(2,2))
    print *,''
    !----------------------------------------------------------------------------------------------------------!
    call res_tree_ga_qa(p6,res_jj)
    call get_ol_res_tree_ga_qa(res_jj_ol)
    print *, 'd  a -> e- e+ g db  ', res_jj(1,1), res_jj_ol(1,1), abs(one-res_jj(1,1)/res_jj_ol(1,1))
    print *, 'db a -> e- e+ g d   ', res_jj(2,1), res_jj_ol(2,1), abs(one-res_jj(2,1)/res_jj_ol(2,1))
    print *, 'u  a -> e- e+ g ub  ', res_jj(1,2), res_jj_ol(1,2), abs(one-res_jj(1,2)/res_jj_ol(1,2))
    print *, 'ub a -> e- e+ g u   ', res_jj(2,2), res_jj_ol(2,2), abs(one-res_jj(2,2)/res_jj_ol(2,2))
    print *,''
    !----------------------------------------------------------------------------------------------------------!
    call res_tree_ga_ag(p6,res_jj(1,1))
    call get_ol_res_tree_ga_ag(res_jj_ol(1,1))
    call res_tree_ga_ga(p6,res_jj(2,1))
    call get_ol_res_tree_ga_ga(res_jj_ol(2,1))
    print *, 'a  g -> e- e+ q qb  ', res_jj(1,1), res_jj_ol(1,1), abs(one-res_jj(1,1)/res_jj_ol(1,1))
    print *, 'g  a -> e- e+ q qb  ', res_jj(2,1), res_jj_ol(2,1), abs(one-res_jj(2,1)/res_jj_ol(2,1))
    print *,''
    !----------------------------------------------------------------------------------------------------------!
    print *,'=================================================================================================='
    print *,''

  end subroutine check_rr_ag_amp

  subroutine check_quark_amp()
    use mod_amplitudes_tree_ppll
    real(dp) :: p_ex(4,6),res_jj(4,2)
    real(dp) :: rescheck_jj(4,2)

    print *, 'amplitude removed'
    
    !! madgraph point
    print *, 'point 0'
    p_ex(:,1) = [0.5000000E+03_dp,  0.0000000E+00_dp,  0.0000000E+00_dp,  0.5000000E+03_dp]
    p_ex(:,2) = [0.5000000E+03_dp,  0.0000000E+00_dp,  0.0000000E+00_dp, -0.5000000E+03_dp]
    p_ex(:,3) = [0.8855133E+02_dp, -0.2210069E+02_dp,  0.4008035E+02_dp, -0.7580543E+02_dp]
    p_ex(:,4) = [0.3283294E+03_dp, -0.1038496E+03_dp, -0.3019338E+03_dp,  0.7649492E+02_dp]
    p_ex(:,5) = [0.1523581E+03_dp, -0.1058810E+03_dp, -0.9770964E+02_dp,  0.4954839E+02_dp]
    p_ex(:,6) = [0.4307611E+03_dp,  0.2318313E+03_dp,  0.3595630E+03_dp, -0.5023788E+02_dp]
    !
    call get_ol('point 0 full',rescheck_jj)
    !call restree_jj_qqll(p_ex,res_jj)
    call printout(res_jj,rescheck_jj)
    
    !! point-1
    print *, 'point 1'
    p_ex(:,1) = [262.77102682298408_dp, 0.0000000000000000_dp, 0.0000000000000000_dp, 262.77102682298408_dp]
    p_ex(:,2) = [262.77102682298408_dp, 0.0000000000000000_dp, 0.0000000000000000_dp,-262.77102682298408_dp]
    p_ex(:,3) = [174.05103057901002_dp, 139.41349232905276_dp, 102.04828739298388_dp, 21.066239398421189_dp]
    p_ex(:,4) = [221.35331169751839_dp,-143.46339406934567_dp,-152.43150832971710_dp, 71.973456425891897_dp]
    p_ex(:,5) = [16.660985559150575_dp,-10.894554696680853_dp, 6.7542261396301644_dp, 10.643192520009668_dp]
    p_ex(:,6) = [113.47672581028918_dp, 14.944456436973752_dp, 43.628994797103068_dp,-103.68288834432275_dp]
    call get_ol('point 1 full',rescheck_jj)
    !call restree_jj_qqll(p_ex,res_jj)
    call printout(res_jj,rescheck_jj)
    
    !! point-2
    print *, 'point 2'
    p_ex(:,1) = [388.97716089524152_dp, 0.0000000000000000_dp, 0.0000000000000000_dp, 388.97716089524152_dp]     
    p_ex(:,2) = [388.97716089524152_dp, 0.0000000000000000_dp, 0.0000000000000000_dp,-388.97716089524152_dp]     
    p_ex(:,3) = [245.10973460211017_dp,-216.51599299646449_dp, 23.190803334289470_dp, 112.52463469917210_dp]     
    p_ex(:,4) = [346.05707316175665_dp, 202.24658698584167_dp, 20.246350722765982_dp,-280.07481361284624_dp]     
    p_ex(:,5) = [23.436287516016595_dp,-2.5655495234645538_dp, 16.440442793195587_dp, 16.504222760848492_dp]     
    p_ex(:,6) = [163.35122651059964_dp, 16.834955534087353_dp,-59.877596850251038_dp, 151.04595615282568_dp]
    call get_ol('point 2 full',rescheck_jj)
    !call restree_jj_qqll(p_ex,res_jj)
    call printout(res_jj,rescheck_jj)
    
    !! point-3
    print *, 'point 3'
    p_ex(:,1) = [515.14919741575125_dp, 0.0000000000000000_dp, 0.0000000000000000_dp, 515.14919741575125_dp]
    p_ex(:,2) = [515.14919741575125_dp, 0.0000000000000000_dp, 0.0000000000000000_dp,-515.14919741575125_dp]
    p_ex(:,3) = [279.87732720788176_dp,-148.20751868121872_dp, 57.092632220451264_dp, 230.44800072364308_dp]
    p_ex(:,4) = [332.85866313699307_dp,-28.274654875998941_dp, 18.486637019755392_dp,-331.13996703612980_dp]
    p_ex(:,5) = [194.03907727526635_dp, 144.44758583658182_dp, 128.61793784073748_dp, 15.603990560632166_dp]
    p_ex(:,6) = [223.52332721136131_dp, 32.034587720635834_dp,-204.19720708094414_dp, 85.087975751854572_dp]
    call get_ol('point 3 full',rescheck_jj)
    !call restree_jj_qqll(p_ex,res_jj)
    call printout(res_jj,rescheck_jj)
    
  contains

    subroutine printout(res,res_ol)
      real(dp), intent(in) :: res(4,2),res_ol(4,2)

      !-- to compare against madgraph: multiply by eesq * (0.118_dp*4*pi)
      print *, 'u ub -> u ub em ep   ', res(1,2), res_ol(1,2), one-res(1,2)/res_ol(1,2)
      print *, 'd db -> d db em ep   ', res(1,1), res_ol(1,1), one-res(1,1)/res_ol(1,1)
      print *, ''
      print *, 'ub u -> u ub em ep   ', res(2,2), res_ol(2,2), one-res(2,2)/res_ol(2,2)
      print *, 'db d -> d db em ep   ', res(2,1), res_ol(2,1), one-res(2,1)/res_ol(2,1)
      print *, ''
      print *, 'u u -> u u em ep     ', res(3,2), res_ol(3,2), one-res(3,2)/res_ol(3,2)
      print *, 'd d -> d d em ep     ', res(3,1), res_ol(3,1), one-res(3,1)/res_ol(3,1)
      print *, ''
      print *, 'ub ub -> ub ub em ep ', res(4,2), res_ol(4,2), one-res(4,2)/res_ol(4,2)
      print *, 'db db -> db db em ep ', res(4,1), res_ol(4,1), one-res(4,1)/res_ol(4,1)
      print *, ''

    end subroutine printout

    subroutine get_ol(str,res)
      character(*), intent(in) :: str
      real(dp), intent(out) :: res(4,2)

      if (str.eq.'point 0/z') then
         res(1,:) = [-2.1229959413148309E-013_dp,-7.6763541788693725E-012_dp]
         res(2,:) = [-1.2393174055680087E-013_dp,-2.1175950570899554E-012_dp]
         res(3,:) = [-1.7030297653586373E-013_dp, 3.1219572019801917E-013_dp]
         res(4,:) = [ 8.2783968861789399E-014_dp,-1.7124998429832077E-012_dp]
      elseif(str.eq.'point 1/z') then
         res(1,:) = [ 7.5165213125212843E-011_dp,-9.8569678485796261E-010_dp]
         res(2,:) = [-1.5959037272301080E-012_dp,-4.3471367187512824E-011_dp]
         res(3,:) = [-4.2014375065389229E-013_dp, 2.3911926628542874E-011_dp]
         res(4,:) = [ 2.1327084692632060E-012_dp, 3.4891088692060724E-012_dp]
      elseif(str.eq.'point 2/z') then
         res(1,:) = [-2.2053989832747040E-011_dp,-3.7977198649773765E-010_dp]
         res(2,:) = [-5.6151545919891327E-011_dp,-4.9886841111369055E-010_dp]
         res(3,:) = [ 1.2901443143620188E-012_dp, 2.7326972880573484E-011_dp]
         res(4,:) = [ 1.8471996352604184E-012_dp, 2.2870530313385747E-011_dp]
      elseif(str.eq.'point 3/z') then
         res(1,:) = [-2.9311702066329598E-013_dp,-1.3966308344494352E-012_dp]
         res(2,:) = [ 1.4149724702416870E-013_dp,-6.6437509783551799E-012_dp]
         res(3,:) = [ 1.5493826711033617E-016_dp, 1.3367158188499441E-013_dp]
         res(4,:) = [ 1.1087652401379337E-014_dp, 4.6209868810841963E-014_dp]
      endif

      if (str.eq.'point 0/a') then
         res(1,2) = -2.8353760575051047E-012_dp; res(1,1) = -6.4860039464737552E-012_dp
         res(2,2) = -1.1678857620677020E-012_dp; res(2,1) = -2.5207570407848939E-012_dp
         res(3,2) = -4.1680956489422371E-013_dp; res(3,1) = -2.2106972462314633E-013_dp
         res(4,2) = -7.8080538943201123E-013_dp; res(4,1) = -4.6428280061545814E-013_dp
      elseif (str.eq.'point 1/a') then;
         res(1,2) =  1.3395142150429449E-010_dp; res(1,1) =  1.6183727944875318E-010_dp
         res(2,2) =  1.2106572206969817E-011_dp; res(2,1) =  4.2783127201081166E-011_dp
         res(3,2) =  2.2715209459496166E-012_dp; res(3,1) =  7.7121909495367383E-012_dp
         res(4,2) =  3.6500996216381944E-012_dp; res(4,1) =  9.7733565973782522E-012_dp
      elseif (str.eq.'point 2/a') then
         res(1,2) =  1.2781911455789292E-010_dp; res(1,1) =  1.7009014819054887E-010_dp
         res(2,2) =  1.1673281313562431E-010_dp; res(2,1) =  1.7049643291574140E-010_dp
         res(3,2) =  5.0378013956362630E-012_dp; res(3,1) =  1.0012833123860131E-011_dp
         res(4,2) =  7.1070468131645155E-012_dp; res(4,1) =  1.4503436995170582E-011_dp
      elseif (str.eq.'point 3/a') then
         res(1,2) = -1.2588339794919218E-012_dp; res(1,1) = -3.7925768287241910E-012_dp
         res(2,2) = -1.0602025463890169E-012_dp; res(2,1) = -3.8382652395413529E-012_dp
         res(3,2) = -3.2729007406980446E-014_dp; res(3,1) =  1.7538052209544243E-013_dp
         res(4,2) = -1.4619242928171652E-014_dp; res(4,1) =  2.3908230384599513E-013_dp
      endif

      if (str.eq.'point 0 full') then
         res(1,2) = -5.3060092210267509E-012_dp; res(1,1) = -2.6231798905905379E-012_dp
         res(2,2) = -5.6181570865486473E-012_dp; res(2,1) = -3.4890769486160454E-012_dp
         res(3,2) =  9.4881339370047883E-013_dp; res(3,1) =  3.8824871950140589E-013_dp
         res(4,2) = -9.0180309808961239E-012_dp; res(4,1) = -8.2497927275456559E-013_dp
      elseif (str.eq.'point 1 full') then
         res(1,2) = -1.9011110644048614E-010_dp; res(1,1) = -3.6069685112667100E-013_dp
         res(2,2) = -5.8367524614137758E-011_dp; res(2,1) =  7.6967821128870871E-011_dp
         res(3,2) = -5.3592398628793560E-012_dp; res(3,1) =  1.0524938158215907E-011_dp
         res(4,2) =  1.4146102262577612E-011_dp; res(4,1) =  3.1784317507109282E-011_dp
      elseif (str.eq.'point 2 full') then
         res(1,2) = -3.2459460905983662E-010_dp; res(1,1) =  2.3912023765571005E-010_dp
         res(2,2) = -1.1537439591391660E-010_dp; res(2,1) =  2.5922326442277396E-011_dp
         res(3,2) =  1.9267561186288120E-011_dp; res(3,1) =  3.8415937766313461E-012_dp
         res(4,2) =  8.9961066845565064E-011_dp; res(4,1) =  3.9504639883337337E-011_dp
      elseif (str.eq.'point 3 full') then
         res(1,2) = -9.6904541926565827E-012_dp; res(1,1) = -9.9437012326132266E-012_dp
         res(2,2) = -2.1088446406297348E-012_dp; res(2,1) = -6.6562199103618892E-013_dp
         res(3,2) = -5.3009612309998847E-013_dp; res(3,1) =  2.0322632721081374E-013_dp
         res(4,2) = -1.0422059837307901E-013_dp; res(4,1) =  6.3192362078973648E-013_dp
      endif
           
    end subroutine get_ol

  end subroutine check_quark_amp

  !---------------------------------------------------------------------------!
  !--         OpenLoops benchmarks for tree-level 6-point amplitudes        --!
  !---------------------------------------------------------------------------!
  subroutine get_ol_res_tree_ga_qqb(res)
    real(dp), intent(out) :: res(2,2)
    !-- OpenLoops result for point p6
    if(cm_scheme) then
      res(1,1) = 3.8449097258451839E-010_dp ! d dx -> e- e+ g a
      res(2,1) = 1.1737018297516170E-009_dp ! dx d -> e- e+ g a
      res(1,2) = 1.6513063853585264E-008_dp ! u ux -> e- e+ g a
      res(2,2) = 9.1386058733605980E-010_dp ! ux u -> e- e+ g a
    else
      res(1,1) = 3.8424737982616583E-010_dp ! d dx -> e- e+ g a
      res(2,1) = 1.1728052948247910E-009_dp ! dx d -> e- e+ g a
      res(1,2) = 1.6499672181720152E-008_dp ! u ux -> e- e+ g a
      res(2,2) = 9.1300566167208474E-010_dp ! ux u -> e- e+ g a
    endif
  end subroutine get_ol_res_tree_ga_qqb

  subroutine get_ol_res_tree_ga_gq(res)
    real(dp), intent(out) :: res(2,2)
    !-- OpenLoops result for point p6
    if(cm_scheme) then
      res(1,1) = 3.1384543723625038E-012_dp ! g dx -> e- e+ dx a
      res(2,1) = 1.1873979020107822E-011_dp ! g d  -> e- e+ d  a
      res(1,2) = 1.6405468414351902E-010_dp ! g ux -> e- e+ ux a
      res(2,2) = 9.0699575329169744E-012_dp ! g u  -> e- e+ u  a
    else
      res(1,1) = 3.1364628807550253E-012_dp ! g dx -> e- e+ dx a
      res(2,1) = 1.1864903900424446E-011_dp ! g d  -> e- e+ d  a
      res(1,2) = 1.6392184660869063E-010_dp ! g ux -> e- e+ ux a
      res(2,2) = 9.0615956825852931E-012_dp ! g u  -> e- e+ u  a
    endif
  end subroutine get_ol_res_tree_ga_gq

  subroutine get_ol_res_tree_ga_qg(res)
    real(dp), intent(out) :: res(2,2)
    !-- OpenLoops result for point p6
    if(cm_scheme) then
      res(1,1) = 3.3088502727303343E-013_dp ! dx g -> e- e+ dx a
      res(2,1) = 4.0022950285405131E-013_dp ! d  g -> e- e+ d  a
      res(1,2) = 3.7022795728663668E-012_dp ! ux g -> e- e+ ux a
      res(2,2) = 1.3601383356465315E-012_dp ! u  g -> e- e+ u  a
    else
      res(1,1) = 3.3066182036842907E-013_dp ! dx g -> e- e+ dx a
      res(2,1) = 3.9993345961429756E-013_dp ! d  g -> e- e+ d  a
      res(1,2) = 3.6992279276553808E-012_dp ! ux g -> e- e+ ux a
      res(2,2) = 1.3590047815671916E-012_dp ! u  g -> e- e+ u  a
    endif
  end subroutine get_ol_res_tree_ga_qg

  subroutine get_ol_res_tree_ga_aq(res)
    real(dp), intent(out) :: res(2,2)
    !-- OpenLoops result for point p6
    if(cm_scheme) then
      res(1,1) = 2.5237106557134734E-009_dp ! a dx -> e- e+ g dx
      res(2,1) = 1.2022053806456836E-009_dp ! a d  -> e- e+ g d 
      res(1,2) = 4.2796546399072683E-009_dp ! a ux -> e- e+ g ux
      res(2,2) = 2.7656040245207594E-009_dp ! a u  -> e- e+ g u 
    else
      res(1,1) = 2.5219010759438396E-009_dp ! a dx -> e- e+ g dx
      res(2,1) = 1.2013111014668857E-009_dp ! a d  -> e- e+ g d 
      res(1,2) = 4.2761783032598248E-009_dp ! a ux -> e- e+ g ux
      res(2,2) = 2.7632721079044771E-009_dp ! a u  -> e- e+ g u 
    endif
  end subroutine get_ol_res_tree_ga_aq

  subroutine get_ol_res_tree_ga_qa(res)
    real(dp), intent(out) :: res(2,2)
    !-- OpenLoops result for point p6
    if(cm_scheme) then
      res(1,1) = 1.1480242834354026E-008_dp ! dx a -> e- e+ g dx
      res(2,1) = 1.9009331168370175E-009_dp ! d  a -> e- e+ g d 
      res(1,2) = 3.9920243263975607E-009_dp ! ux a -> e- e+ g ux
      res(2,2) = 3.3624232384793256E-008_dp ! u  a -> e- e+ g u 
    else
      res(1,1) = 1.1471889743953151E-008_dp ! dx a -> e- e+ g dx
      res(2,1) = 1.8992848612165219E-009_dp ! d  a -> e- e+ g d 
      res(1,2) = 3.9889082869132259E-009_dp ! ux a -> e- e+ g ux
      res(2,2) = 3.3597746937612605E-008_dp ! u  a -> e- e+ g u 
    endif
  end subroutine get_ol_res_tree_ga_qa

  subroutine get_ol_res_tree_ga_ga(res)
    real(dp), intent(out) :: res
    real(dp) :: res_up, res_dn
    !-- OpenLoops result for point p6
    if(cm_scheme) then
      res_dn = 1.9707883295247697E-011_dp ! a g -> e- e+ d dx
      res_up = 2.7222722241845016E-010_dp ! a g -> e- e+ u ux
    else
      res_dn = 1.9691035776726455E-011_dp ! a g -> e- e+ d dx
      res_up = 2.7201417400670489E-010_dp ! a g -> e- e+ u ux
    endif
    res = ndn*res_dn + nup*res_up
  end subroutine get_ol_res_tree_ga_ga

  subroutine get_ol_res_tree_ga_ag(res)
    real(dp), intent(out) :: res
    real(dp) :: res_up, res_dn
    !-- OpenLoops result for point p6
    if(cm_scheme) then
      res_dn = 9.9716650945148547E-011_dp ! a g -> e- e+ d dx
      res_up = 8.6549199402927041E-011_dp ! a g -> e- e+ u ux
    else
      res_dn = 9.9641682621505481E-011_dp ! a g -> e- e+ d dx
      res_up = 8.6491212766078802E-011_dp ! a g -> e- e+ u ux
    endif
    res = ndn*res_dn + nup*res_up
  end subroutine get_ol_res_tree_ga_ag

end module mod_checks

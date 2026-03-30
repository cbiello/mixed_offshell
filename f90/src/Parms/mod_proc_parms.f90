module mod_proc_parms
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_parser
  implicit none
  private

  real(dp), public, save :: sh

  real(dp), public, save :: mu,xmuR,xmuF
  logical, public, save :: dynscale

  real(dp), public, save :: buff
  real(dp), public, save :: onet

  real(dp), public, save :: buff_rr,buff_rv,buff_z,buff_r
  
  real(dp), public, save :: units

  character(30), public, save :: pdfset,sec,corr,ch,histodir
  integer, public, save :: pdfmem

  character(strlen), public :: outputfile,pref

  integer, public, save :: failed_points,icount_good,icount_nan,icount_acc
  integer, public, save :: acc_failed_points,acc_icount_good,acc_icount_nan,acc_icount_acc

  public :: get_proc_parms,help_proc_parms
  public :: set_outfile
  public :: set_proc_parms,print_proc_parms
  public :: print_counters,reset_counters
  public :: set_buff

#if (_Vcharge == 0)
  real(dp), public, parameter :: current_charge = zero
#elif (_Vcharge == -1)
  real(dp), public, parameter :: current_charge = -one
#elif (_Vcharge == +1)
  real(dp), public, parameter :: current_charge = +one
#endif
  
  real(dp), public, parameter :: Q_lep = Qel
  real(dp), public, parameter :: Q_lep2 = Q_lep**2
  real(dp), public, save :: cL_lep
  real(dp), public, save :: cR_lep
  complex(dp), public, save :: cms_cL_lep
  complex(dp), public, save :: cms_cR_lep

  logical, public, save  :: override   = .false.
  logical, public, save  :: checkpoles = .false.

  !-- for checking limits
  logical, public, save  :: checklim = .false.
  real(dp), public, save :: lim

  !-- cuts, etc

  real(dp), public, save :: qmin,qmax,q2min,q2max,q2min_tech,q2max_tech
  real(dp), public, save :: ptlep_cut,ylep_cut,ya_beam
  real(dp), public, save :: ptlm_min,ptlm_max,ptlp_min,ptlp_max
  real(dp), public, save :: qmin_rescaling
  real(dp), public, save :: qmax_rescaling
  
  real(dp), public, save :: ptmiss_cut,mt_min


  real(dp), public, save :: R_phot
  logical, public, save :: rec_phot_beam !-- whether collinear photons should be recombined in beam

  integer, public, save :: taumode = 1

contains
  
  subroutine get_proc_parms()
    character(30) :: str_units
    real(dp) :: ecoll

    !-- default directory for output histograms
    histodir = trim(string_val_opt('-histodir','./histo'))
    histodir = trim(histodir)//'/'

    !-- collider energy
    ecoll = real_val_opt('-ecoll',14000._dp)
    sh = ecoll**2

    !-- PDFs
    pdfset = trim(string_val_opt('-pdfset','NNPDF31_nnlo_as_0118_luxqed'))
    pdfmem = int_val_opt('-pdfmem',0)
    
    !-- scale
    mu = real_val_opt('-mu',91.1876_dp)
    dynscale = log_val_opt_witharg('-dynscale',.true.)
    !
    xmuR = real_val_opt('-xmuR',one)
    xmuF = real_val_opt('-xmuF',one)

    !-- tau=mll^2/S generation mode
    taumode = int_val_opt('-taumode',1)

    !-- cuts
    qmin = real_val_opt('-qmin',200._dp)
    qmax = real_val_opt('-qmax',ecoll)
    !
    q2min = qmin**2
    q2max = qmax**2
    qmin_rescaling = real_val_opt('-qmin_rsc',20._dp)
    qmax_rescaling = real_val_opt('-qmax_rsc',20._dp)
    !
    ptlep_cut = real_val_opt('-ptlep_cut',20._dp)
    ylep_cut  = real_val_opt('-ylep_cut',2.5_dp)
    !
    ptmiss_cut = real_val_opt('-ptmiss_cut',20._dp)
    mt_min = real_val_opt('-mt_min',60._dp)
    !
    ptlm_min = real_val_opt('-ptlm_min',ptlep_cut)
    ptlm_max = real_val_opt('-ptlm_max',ecoll)
    ptlp_min = real_val_opt('-ptlp_min',ptlep_cut)
    ptlp_max = real_val_opt('-ptlp_max',ecoll)
    ya_beam = real_val_opt('-ya_beam',3.0_dp)
    !
    R_phot = real_val_opt('-R_phot',0.1_dp)
    !
    rec_phot_beam = log_val_opt_witharg('-phot_beam',.false.)

    !-- units
    str_units = trim(string_val_opt('-units','fb'))
    if (str_units=='fb') then
       units = GeVtoFb
    elseif (str_units=='pb') then
       units = GeVtoPb
    elseif (str_units=='nb') then
       units = GeVtoNb
    else
       print *, 'units not implemeted'
       stop
    endif
    
    !-- buffers
    buff = real_val_opt('-buff',1E-12_dp)
    !buff_rr = real_val_opt('-buff_rr',1E-12_dp)
    buff_rr = real_val_opt('-buff_rr',1E-10_dp)
    buff_rv = real_val_opt('-buff_rv',1E-12_dp)
    buff_z = real_val_opt('-buff_z',1E-12_dp)
    buff_r = real_val_opt('-buff_r',1E-12_dp)

    !-- output file and grid
    outputfile = trim(string_val_opt('-outfile','auto'))
    pref       = trim(string_val_opt('-pref','no'))

    !-- check poles
    checkpoles = log_val_opt_witharg('-checkpoles',.false.)

    !-- checklimits
    checklim = log_val_opt('-checklim',.false.)
    lim = real_val_opt('-lim',1E-8_dp)

  end subroutine get_proc_parms

  subroutine help_proc_parms(idev)
    integer, intent(in) :: idev

    write(idev,*) ''
    write(idev,*) ' -coll 14000'
    write(idev,*) ''
    write(idev,*) ' -pdfset NNPDF31_nnlo_as_0118_luxqed'
    write(idev,*) ' -pdfmem 0'
    write(idev,*) ''
    write(idev,*) ' -dynscale     --> runs with mu = mll/2'
    write(idev,*) ' -mu 91.1876   --> if dynscale = .false.'
    write(idev,*) ' -xmuR 1'
    write(idev,*) ' -xmuF 1'
    write(idev,*) ''
    write(idev,*) ' -taumode 1'
    write(idev,*) ''
    write(idev,*) ' -qmin_rsc 20'
    write(idev,*) ' -qmax_rsc 20'
    write(idev,*) ' -qmin 300'
    write(idev,*) ' -qmax 14000'
    write(idev,*) ' -ptlep_cut 20'
    write(idev,*) ' -ylep_cut 2.5'
    write(idev,*) ' -R_phot 0.1_dp'
    write(idev,*) ''
    write(idev,*) ' -units fb/pb/nb'
    write(idev,*) ''
    write(idev,*) ' -buff 1D-12'
    write(idev,*) ' -buff_rr 1D-12'
    write(idev,*) ' -buff_rv 1D-12'
    write(idev,*) ' -buff_z 1D-12'
    write(idev,*) ' -buff_r 1D-12'
    write(idev,*) ''
    write(idev,*) ' -outfile auto --> override automatic file naming'
    write(idev,*) ' -pref no  --> add a prefix to automatic file naming'
    write(idev,*) ''
    write(idev,*) ' -checkpoles --> check 1L poles, if config debug is active'
    write(idev,*) ''
    write(idev,*) ' -checklim --> check limits'
    write(idev,*) ' -lim      --> how deep we want to go with the limit'
    
  end subroutine help_proc_parms
  
  !-- without the ``histo'' prefix and suffix
  subroutine set_outfile(pref,myoutfile,mygridfile,seednr)
    character(*), intent(in) :: pref
    character(strlen) :: tmpfile,myoutfile,mygridfile
    integer :: seednr
    character(4) :: ch_xmur,ch_xmuf
    character(len=5) :: seedstr

    write(seedstr, '(I5.5)') seednr

    if (pref.eq.'no') then
#if defined(ANALYSIS)
       tmpfile = ANALYSIS
#else
       tmpfile = 'test'
#endif    
    else
       tmpfile = trim(pref)
    endif
    tmpfile = trim(tmpfile)//"_"//trim(corr)//'_'//trim(ch)//'_'//trim(sec)
    if (dynscale) then
       tmpfile=trim(tmpfile)//"_dynscale"
    else
       tmpfile=trim(tmpfile)//"_fixscale"
    endif
    write(ch_xmur,'(F4.2)') xmur
    write(ch_xmuf,'(F4.2)') xmuf
    tmpfile = trim(tmpfile)//"_xMuR_"//trim(ch_xmur)//"_xMuF_"//trim(ch_xmuf)//'_s'//trim(seedstr)

    if (myoutfile  .eq. 'auto') myoutfile  = tmpfile
    if (mygridfile .eq. 'auto') mygridfile = tmpfile

    myoutfile  = trim(histodir)//trim(myoutfile)//'.histo'
    mygridfile = './grid/'//trim(mygridfile)//'.grid'
    
  end subroutine set_outfile

  subroutine set_proc_parms()

    failed_points = 0
    call set_buff()

    cL_lep = cLel
    cR_lep = cRel

    cms_cL_lep = cms_cLel
    cms_cR_lep = cms_cRel

    call set_tech_q2minmax()
    
  end subroutine set_proc_parms

  subroutine set_buff()
    onet = 1._dp - 2._dp * buff
  end subroutine set_buff

  subroutine set_tech_q2minmax()
    
    q2min_tech = q2min
    q2max_tech = q2max
    
    !-- set q2min technical
    if (corr.eq.'nloewk') then
       if (ch.eq.'ns' .and. sec.eq.'r_is') call rescale_q2min()
    endif

    if (corr.eq.'nnlo') then
       if (ch.eq.'ns_ga') then
          if (sec(1:7).eq.'rr_5161')  call rescale_q2min()
          if (sec(1:7).eq.'rr_5262')  call rescale_q2min()
          if (sec(1:7).eq.'rr_5162')  call rescale_q2min()
          if (sec(1:7).eq.'rr_5261')  call rescale_q2min()
       elseif (ch.eq.'ns') then
          if (sec(1:8).eq.'rvqcd_is') call rescale_q2min()
          if (sec(1:9).eq.'s_oewk_is') call rescale_q2min()
       elseif (ch.eq.'gq') then
          if (sec(1:7).eq.'rr_5161') call rescale_q2min()
          if (sec(1:7).eq.'rr_5162') call rescale_q2min()
          if (sec(1:9).eq.'s_oewk_is') call rescale_q2min()
       elseif (ch.eq.'qg') then
          if (sec(1:7).eq.'rr_5262') call rescale_q2min()
          if (sec(1:7).eq.'rr_5261') call rescale_q2min()
          if (sec(1:9).eq.'s_oewk_is') call rescale_q2min()
       endif
    endif

    !-- set q2max technical
    if (corr.eq.'nloewk') then
       if (ch.eq.'ns' .and. sec(1:4).eq.'r_fs') call rescale_q2max()
       if (ch.eq.'aa' .and. sec(1:4).eq.'r_fs') call rescale_q2max()
    endif

    if (corr.eq.'nnlo') then
       if (ch.eq.'ns_ga') then
          if (sec(1:7).eq.'rr_5163')  call rescale_q2max()
          if (sec(1:7).eq.'rr_5164')  call rescale_q2max()
          if (sec(1:7).eq.'rr_5263')  call rescale_q2max()
          if (sec(1:7).eq.'rr_5264')  call rescale_q2max()
       elseif (ch.eq.'ns') then
          if (sec(1:8).eq.'rvqcd_fs') call rescale_q2max()
          if (sec(1:9).eq.'s_oewk_fs') call rescale_q2max()
       elseif (ch.eq.'gq') then
          if (sec(1:7).eq.'rr_5163') call rescale_q2max()
          if (sec(1:7).eq.'rr_5164') call rescale_q2max()
          if (sec(1:9).eq.'s_oewk_fs') call rescale_q2max()
       elseif (ch.eq.'qg') then
          if (sec(1:7).eq.'rr_5263') call rescale_q2max()
          if (sec(1:7).eq.'rr_5264') call rescale_q2max()
          if (sec(1:9).eq.'s_oewk_fs') call rescale_q2max()
       endif
    endif
    
    !-- set q2max technical
    if (abs(q2max-sh).lt.1E-14_dp) then
       return
    else
       write (6,*) '# warning, q2max cut-off smaller than collider energy'
    endif

    if (q2max_tech.gt.sh) then
       write (6,*) '# q2max technical cut-off larger than collider energy -> set it to sh'
       q2max_tech = sh
    endif
    
  contains
    
    subroutine rescale_q2min()
      q2min_tech = q2min/qmin_rescaling**2
    end subroutine rescale_q2min

    subroutine rescale_q2max()
      q2max_tech = q2max*qmax_rescaling**2
    end subroutine rescale_q2max

  end subroutine set_tech_q2minmax

  !---------------------

  subroutine print_counters(outdev,iter,acc)
    integer, intent(in) :: outdev
    logical, intent(in) :: iter,acc

    write(outdev,*) '#'
    
    if(iter) then
       write(outdev,*) '# Iteration, NaN points:                         ', icount_nan
       write(outdev,*) '# Iteration, Points failing the technical cuts:  ', failed_points
       write(outdev,*) '# Iteration, Good points passing technical cuts: ', icount_good
       write(outdev,*) '# Iteration, Points passing physical cuts:       ', icount_acc
       write(outdev,*) '# Iteration, Acceptance:                         ', real(icount_acc,kind=dp)/real(icount_good,kind=dp)
       write(outdev,*) '#'
    endif

    if(acc) then
       write(outdev,*) '# NaN points:                         ', acc_icount_nan
       write(outdev,*) '# Points failing the technical cuts:  ', acc_failed_points
       write(outdev,*) '# Good points passing technical cuts: ', acc_icount_good
       write(outdev,*) '# Points passing physical cuts:       ', acc_icount_acc
       write(outdev,*) '# Acceptance:                         ', real(acc_icount_acc,kind=dp)/real(acc_icount_good,kind=dp)   
    endif
    
    write(outdev,*) '#'
        
  end subroutine print_counters

  subroutine reset_counters(all)
    logical, intent(in) :: all 

    icount_nan    = 0
    failed_points = 0
    icount_good   = 0
    icount_acc    = 0

    if (all) then
       acc_icount_nan    = 0
       acc_failed_points = 0
       acc_icount_good   = 0
       acc_icount_acc    = 0
    endif
    
  end subroutine reset_counters
    
  subroutine print_proc_parms(outdev)
    integer, intent(in) :: outdev
    real(dp15) :: alphasPDF

    write(outdev,*) '# sqrt(sh) = ', sqrt(sh)
    write(outdev,*) '#'
    write(outdev,*) '# corr = ', trim(adjustl(corr))
    write(outdev,*) '# sec  = ', trim(adjustl(sec))
    write(outdev,*) '# ch   = ', trim(adjustl(ch))
    write(outdev,*) '#'
    write(outdev,*) '# pdfset = ', trim(adjustl(pdfset))
    write(outdev,*) '# pdfmem = ', pdfmem
    write(outdev,*) '#'
    if (dynscale) then
       write(outdev,*) '# using dynamical scale, see cut_histo'
    else
       write(outdev,*) '# using fixed scale'
       write(outdev,*) '# mu = ', mu
    endif
    write(outdev,*) '# xmuR =        ', xmuR
    write(outdev,*) '# xmuF =        ', xmuF
    write(outdev,*) '#'
    write(outdev,*) '# alphasPDF(mu) = ', alphasPDF(real(mu,kind=dp15))
    write(outdev,*) '#'
    write(outdev,*) '# qmin = ',qmin
    write(outdev,*) '# qmax = ',qmax
    write(outdev,*) '# taumode = ', taumode
    write(outdev,*) '#'
    write(outdev,*) '# technical cut on qmin = ',sqrt(abs(q2min_tech))
    write(outdev,*) '# technical cut on qmax = ',sqrt(abs(q2max_tech))
    write(outdev,*) '#'
    write(outdev,*) '# ptlep_cut = ',ptlep_cut
    write(outdev,*) '# ylep_cut  = ',ylep_cut
    write(outdev,*) '#'
    write(outdev,*) '# R_phot =', R_phot, '# for isolation'
    write(outdev,*) '#'
    write(outdev,*) '# buff =      ', buff
    write(outdev,*) '# buff_rr =      ', buff_rr
    write(outdev,*) '# buff_rv =      ', buff_rv
    write(outdev,*) '# buff_z  =      ', buff_z
    write(outdev,*) '# buff_r  =      ', buff_r
    write(outdev,*) '#'
    write(outdev,*) '# lepton charge = ', Q_lep
    write(outdev,*) '# lepton cL = ', cL_lep
    write(outdev,*) '# lepton cR = ', cR_lep
    write(outdev,*) '#'
    write(outdev,*) '# out = ', outputfile
    write(outdev,*) '#'
    if (checkpoles) then
       write(outdev,*) '# checking poles of 1L amplitudes'
    endif
    if (checklim) then
       write(outdev,*) '# checking limits with lim depth = ', lim
    endif
        
  end subroutine print_proc_parms
  
end module mod_proc_parms


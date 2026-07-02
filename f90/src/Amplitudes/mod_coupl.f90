module mod_coupl
  use mod_types
  use mod_consts_dp
  use mod_parms
  use mod_proc_parms
  implicit none
  private

  public :: get_coupl, get_coupl_w, quark_type, get_coupl_gen

  interface get_coupl
     module procedure get_coupl_real, get_coupl_cmplx
  end interface get_coupl

    interface get_coupl_gen
     module procedure get_coupl_real_gen, get_coupl_cmplx_gen
  end interface get_coupl_gen


  !interface get_coupl_w
  !   module procedure get_coupl_w_real, get_coupl_w_cmplx
  !end interface get_coupl_w

contains





  subroutine get_coupl_real_gen(mll,Q1,Q2,cL1,cL2,cR1,cR2,coupl,cin)
    integer,optional :: cin !if charged current mediated, cin=1
    real(dp), intent(in)     :: mll
    real(dp), intent(in)     :: Q1(:),Q2(:),cL1(:),cL2(:),cR1(:),cR2(:)
    complex(dp), intent(out) :: coupl(1:3,-1:1,-1:1)  ! 1=up-type, 2-dn type, 3=fl changing
    complex(dp) :: propZ, propW

    propZ = mll/(mll-mzsq_prop)
    propW = mll/(mll-mwsq_prop)

#if (_Vcharge == 0)
    ! Z boson
    coupl = zero
    coupl(1:2,-1,-1) = Q1(1:2)*Q2(1:2) + propZ * cL1(1:2) * cL2(1:2)
    coupl(1:2,-1,+1) = Q1(1:2)*Q2(1:2) + propZ * cL1(1:2) * cR2(1:2)
    coupl(1:2,+1,-1) = Q1(1:2)*Q2(1:2) + propZ * cR1(1:2) * cL2(1:2)
    coupl(1:2,+1,+1) = Q1(1:2)*Q2(1:2) + propZ * cR1(1:2) * cR2(1:2)
#else
    !W+ or W-
    coupl = zero
    coupl(3,-1,-1) = propW * cLW * cLWnue
    print *, clw,clwnue
    stop
#endif

  end subroutine get_coupl_real_gen

  subroutine get_coupl_cmplx_gen(mll,Q1,Q2,cL1,cL2,cR1,cR2,coupl,cin)
    integer,optional :: cin !if charged current mediated, cin=1
    real(dp), intent(in)     :: mll
    real(dp), intent(in)     :: Q1(:),Q2(:)
    complex(dp), intent(in)  :: cL1(:),cL2(:),cR1(:),cR2(:)
    complex(dp), intent(out) :: coupl(1:3,-1:1,-1:1)
    complex(dp) :: propZ, propW

    propZ = mll/(mll-mzsq_prop)
    propW = mll/(mll-mwsq_prop)


#if (_Vcharge == 0)
    ! Z boson
    coupl = zero
    coupl(1:2,-1,-1) = Q1(1:2)*Q2(1:2) + propZ * cL1(1:2) * cL2(1:2)
    coupl(1:2,-1,+1) = Q1(1:2)*Q2(1:2) + propZ * cL1(1:2) * cR2(1:2)
    coupl(1:2,+1,-1) = Q1(1:2)*Q2(1:2) + propZ * cR1(1:2) * cL2(1:2)
    coupl(1:2,+1,+1) = Q1(1:2)*Q2(1:2) + propZ * cR1(1:2) * cR2(1:2)       
#else
    ! W+ or W-        
    coupl = czero
    coupl(3,-1,-1) = propW * cms_cLW* cms_cLWnue
#endif
    
  end subroutine get_coupl_cmplx_gen
  


  
  subroutine get_coupl_real(mll,Q1,Q2,cL1,cL2,cR1,cR2,coupl,cin)
    integer,optional :: cin !if charged current mediated, cin=1
    real(dp), intent(in)     :: mll
    real(dp), intent(in)     :: Q1(:),Q2(:),cL1(:),cL2(:),cR1(:),cR2(:)
    complex(dp), intent(out) :: coupl(1:2,-1:1,-1:1)  ! 1=up-type, 2-dn type, 3=fl changing
    complex(dp) :: propZ, propW

    propZ = mll/(mll-mzsq_prop)
    propW = mll/(mll-mwsq_prop)

#if (_Vcharge == 0)
    ! Z boson
    coupl = zero
    coupl(1:2,-1,-1) = Q1(1:2)*Q2(1:2) + propZ * cL1(1:2) * cL2(1:2)
    coupl(1:2,-1,+1) = Q1(1:2)*Q2(1:2) + propZ * cL1(1:2) * cR2(1:2)
    coupl(1:2,+1,-1) = Q1(1:2)*Q2(1:2) + propZ * cR1(1:2) * cL2(1:2)
    coupl(1:2,+1,+1) = Q1(1:2)*Q2(1:2) + propZ * cR1(1:2) * cR2(1:2)
#else
    !W+ or W-
    coupl = zero
    coupl(:,-1,-1) = propW * cLW * cLWnue
#endif

  end subroutine get_coupl_real

  subroutine get_coupl_cmplx(mll,Q1,Q2,cL1,cL2,cR1,cR2,coupl,cin)
    integer,optional :: cin !if charged current mediated, cin=1
    real(dp), intent(in)     :: mll
    real(dp), intent(in)     :: Q1(:),Q2(:)
    complex(dp), intent(in)  :: cL1(:),cL2(:),cR1(:),cR2(:)
    complex(dp), intent(out) :: coupl(1:2,-1:1,-1:1)
    complex(dp) :: propZ, propW

    propZ = mll/(mll-mzsq_prop)
    propW = mll/(mll-mwsq_prop)

#if (_Vcharge == 0)
    ! Z boson
        coupl = zero
        coupl(1:2,-1,-1) = Q1(1:2)*Q2(1:2) + propZ * cL1(1:2) * cL2(1:2)
        coupl(1:2,-1,+1) = Q1(1:2)*Q2(1:2) + propZ * cL1(1:2) * cR2(1:2)
        coupl(1:2,+1,-1) = Q1(1:2)*Q2(1:2) + propZ * cR1(1:2) * cL2(1:2)
        coupl(1:2,+1,+1) = Q1(1:2)*Q2(1:2) + propZ * cR1(1:2) * cR2(1:2)       
#else
    ! W+ or W-        
        coupl = czero
        coupl(:,-1,-1) = propW * cms_cLW* cms_cLWnue
#endif

      end subroutine get_coupl_cmplx

      
  !-- coupling of w to quark and lepton line
  !-- index: down/up, L/R_quark, L/R_lepton
  subroutine get_coupl_w(mll,coupl)
    real(dp), intent(in) :: mll
    complex(dp), intent(out) :: coupl(1:2,-1:1,-1:1)
    complex(dp) :: propW

    !chiara 
    !real(dp) :: mwMG, gawMG
    !mwMG = 80.419002445756163_dp
    !gawMG = 2.0476000000000001_dp
    !propW = mll/(mll -  mwMG**2 + ci * mwMG * gawMG )

    propW = mll/(mll-mwsq_prop)

    !u + dbar > W+
    !d + ubar > W-
    coupl(:,-1,-1) = propW * cms_cLWud * cms_cLWnue
    coupl(:,-1,+1) = czero
    coupl(:,+1,-1) = czero
    coupl(:,+1,+1) = czero

    return

  end subroutine get_coupl_w

  
  ! !-- coupling of a/z to quark and lepton line
  ! !-- index: down/up, L/R_quark, L/R_lepton
  ! subroutine get_coupl_az(mll,coupl)
  !   real(dp), intent(in) :: mll
  !   complex(dp), intent(out) :: coupl(-1:1,-1:1,1:2)
  !   complex(dp) :: propZ
    
  !   propZ = mll/(mll-mzsq_prop)
    
  !   coupl(-1,-1,1) = eesq * (Qdn*Q_lep + propZ * cLdn * cL_lep)
  !   coupl(-1,+1,1) = eesq * (Qdn*Q_lep + propZ * cLdn * cR_lep)
  !   coupl(+1,-1,1) = eesq * (Qdn*Q_lep + propZ * cRdn * cL_lep)
  !   coupl(+1,+1,1) = eesq * (Qdn*Q_lep + propZ * cRdn * cR_lep)
                                      
  !   coupl(-1,-1,2) = eesq * (Qup*Q_lep + propZ * cLup * cL_lep)
  !   coupl(-1,+1,2) = eesq * (Qup*Q_lep + propZ * cLup * cR_lep)
  !   coupl(+1,-1,2) = eesq * (Qup*Q_lep + propZ * cRup * cL_lep)
  !   coupl(+1,+1,2) = eesq * (Qup*Q_lep + propZ * cRup * cR_lep)

  !   return

  ! end subroutine get_coupl_az

  ! subroutine get_coupl_az_cms(mll,coupl)
  !   real(dp), intent(in) :: mll
  !   complex(dp), intent(out) :: coupl(-1:1,-1:1,1:2)
  !   complex(dp) :: propZ
    
  !   propZ = mll/(mll-mzsq_prop)
    
  !   coupl(-1,-1,1) = eesq * (Qdn*Q_lep + propZ * cms_cLdn * cms_cL_lep)
  !   coupl(-1,+1,1) = eesq * (Qdn*Q_lep + propZ * cms_cLdn * cms_cR_lep)
  !   coupl(+1,-1,1) = eesq * (Qdn*Q_lep + propZ * cms_cRdn * cms_cL_lep)
  !   coupl(+1,+1,1) = eesq * (Qdn*Q_lep + propZ * cms_cRdn * cms_cR_lep)
                                      
  !   coupl(-1,-1,2) = eesq * (Qup*Q_lep + propZ * cms_cLup * cms_cL_lep)
  !   coupl(-1,+1,2) = eesq * (Qup*Q_lep + propZ * cms_cLup * cms_cR_lep)
  !   coupl(+1,-1,2) = eesq * (Qup*Q_lep + propZ * cms_cRup * cms_cL_lep)
  !   coupl(+1,+1,2) = eesq * (Qup*Q_lep + propZ * cms_cRup * cms_cR_lep)

  !   return

  ! end subroutine get_coupl_az_cms

  ! subroutine get_coupl_az_noeesq(mll,coupl)
  !   real(dp), intent(in) :: mll
  !   complex(dp), intent(out) :: coupl(-1:1,-1:1,1:2)
  !   complex(dp) :: propZ
    
  !   propZ = mll/(mll-mzsq + ci*mz*gaz)

  !   coupl(-1,-1,1) = (Qdn*Q_lep + propZ * cLdn * cL_lep)
  !   coupl(-1,+1,1) = (Qdn*Q_lep + propZ * cLdn * cR_lep)
  !   coupl(+1,-1,1) = (Qdn*Q_lep + propZ * cRdn * cL_lep)
  !   coupl(+1,+1,1) = (Qdn*Q_lep + propZ * cRdn * cR_lep)
                               
  !   coupl(-1,-1,2) = (Qup*Q_lep + propZ * cLup * cL_lep)
  !   coupl(-1,+1,2) = (Qup*Q_lep + propZ * cLup * cR_lep)
  !   coupl(+1,-1,2) = (Qup*Q_lep + propZ * cRup * cL_lep)
  !   coupl(+1,+1,2) = (Qup*Q_lep + propZ * cRup * cR_lep)

  !   return

  ! end subroutine get_coupl_az_noeesq

  ! !-- coupling of a/z to quark lines, of the same flavour
  ! !-- index: down/up, L/R_quark, L/R_quark
  ! !-- careful: this has eesq factored out
  ! subroutine get_coupl_az_qq_noeesq(mll,coupl)
  !   real(dp), intent(in) :: mll
  !   complex(dp), intent(out) :: coupl(-1:1,-1:1,1:2)
  !   complex(dp) :: propZ
    
  !   propZ = mll/(mll-mzsq + ci*mz*gaz)

  !   coupl(-1,-1,1) = (Qdn*Qdn + propZ * cLdn * cLdn) ! * eesq
  !   coupl(-1,+1,1) = (Qdn*Qdn + propZ * cLdn * cRdn) ! * eesq
  !   coupl(+1,-1,1) = (Qdn*Qdn + propZ * cRdn * cLdn) ! * eesq
  !   coupl(+1,+1,1) = (Qdn*Qdn + propZ * cRdn * cRdn) ! * eesq

  !   coupl(-1,-1,2) = (Qup*Qup + propZ * cLup * cLup) ! * eesq
  !   coupl(-1,+1,2) = (Qup*Qup + propZ * cLup * cRup) ! * eesq
  !   coupl(+1,-1,2) = (Qup*Qup + propZ * cRup * cLup) ! * eesq
  !   coupl(+1,+1,2) = (Qup*Qup + propZ * cRup * cRup) ! * eesq

  !   return
  ! end subroutine get_coupl_az_qq_noeesq

  
  function quark_type(i)
    ! maps PDG-like quark number (u=2, db=-1, etc) to 1=dn-type, 2=up-type, used in get_coupl
    integer, intent(in)   :: i
    integer               :: quark_type

    
    if (mod(abs(i),2) .eq. 0) then       ! up-type quarks
       quark_type = 2
    elseif (mod(abs(i),2) .eq. 1) then    ! down-type quarks
       quark_type = 1
    endif


  end function quark_type

end module mod_coupl

  !-- coupling of a/z to quark and lepton line
  !-- index: down/up, L/R_quark, L/R_lepton
  subroutine get_coupl_az(mll,coupl)
    real(dp), intent(in) :: mll
    complex(dp), intent(out) :: coupl(-1:1,-1:1,1:2)
    complex(dp) :: propZ
    
    propZ = mll/(mll-mzsq + ci*mz*gaz)

#if (_zOnly == 0)
    coupl(-1,-1,1) = eesq * (Qdn*Q_lep + propZ * cLdn * cL_lep)
    coupl(-1,+1,1) = eesq * (Qdn*Q_lep + propZ * cLdn * cR_lep)
    coupl(+1,-1,1) = eesq * (Qdn*Q_lep + propZ * cRdn * cL_lep)
    coupl(+1,+1,1) = eesq * (Qdn*Q_lep + propZ * cRdn * cR_lep)
                                      
    coupl(-1,-1,2) = eesq * (Qup*Q_lep + propZ * cLup * cL_lep)
    coupl(-1,+1,2) = eesq * (Qup*Q_lep + propZ * cLup * cR_lep)
    coupl(+1,-1,2) = eesq * (Qup*Q_lep + propZ * cRup * cL_lep)
    coupl(+1,+1,2) = eesq * (Qup*Q_lep + propZ * cRup * cR_lep)
#elif (_zOnly == 1)
    coupl(-1,-1,1) = eesq * (propZ * cLdn * cL_lep)
    coupl(-1,+1,1) = eesq * (propZ * cLdn * cR_lep)
    coupl(+1,-1,1) = eesq * (propZ * cRdn * cL_lep)
    coupl(+1,+1,1) = eesq * (propZ * cRdn * cR_lep)
                                      
    coupl(-1,-1,2) = eesq * (propZ * cLup * cL_lep)
    coupl(-1,+1,2) = eesq * (propZ * cLup * cR_lep)
    coupl(+1,-1,2) = eesq * (propZ * cRup * cL_lep)
    coupl(+1,+1,2) = eesq * (propZ * cRup * cR_lep)
#endif

    return

  end subroutine get_coupl_az

  !-- coupling of w to quark and lepton line
  !-- index: down/up, L/R_quark, L/R_lepton
  subroutine get_coupl_w(mll,coupl)
    real(dp), intent(in) :: mll
    complex(dp), intent(out) :: coupl(-1:1,-1:1,1:2)
    complex(dp) :: propW
    
    propW = mll/(mll-mwsq + ci*mw*gaw)
    
    !u + dbar > W+
    coupl(-1,-1,1) = eesq * (propW * cLWud * cLWnue)
    coupl(-1,+1,1) = czero
    coupl(+1,-1,1) = czero
    coupl(+1,+1,1) = czero
    
    !d + ubar > W-
    coupl(-1,-1,2) = eesq * (propW * cLWud * cLWnue)
    coupl(-1,+1,2) = czero
    coupl(+1,-1,2) = czero
    coupl(+1,+1,2) = czero

    return

  end subroutine get_coupl_w

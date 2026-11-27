  !-- coupling of a/z to quark and lepton line
  !-- index: down/up, L/R_quark, L/R_lepton
  subroutine get_coupl_az(mll,coupl)
    real(dp), intent(in) :: mll
    complex(dp), intent(out) :: coupl(-1:1,-1:1,1:2)
    complex(dp) :: propZ
    
    propZ = mll/(mll-mzsq + ci*mz*gaz)


    coupl(-1,-1,1) = eesq * (Qdn*Q_lep*zero + propZ * cLdn * cL_lep)
    coupl(-1,+1,1) = eesq * (Qdn*Q_lep*zero + propZ * cLdn * cR_lep)
    coupl(+1,-1,1) = eesq * (Qdn*Q_lep*zero + propZ * cRdn * cL_lep)
    coupl(+1,+1,1) = eesq * (Qdn*Q_lep*zero + propZ * cRdn * cR_lep)                                      
    coupl(-1,-1,2) = eesq * (Qup*Q_lep*zero + propZ * cLup * cL_lep)
    coupl(-1,+1,2) = eesq * (Qup*Q_lep*zero + propZ * cLup * cR_lep)
    coupl(+1,-1,2) = eesq * (Qup*Q_lep*zero + propZ * cRup * cL_lep)
    coupl(+1,+1,2) = eesq * (Qup*Q_lep*zero + propZ * cRup * cR_lep)

    return

  end subroutine get_coupl_az

    write(iunit,*) '#Run parameters '
    write(iunit,*) '#Process ', 'az'
    write(iunit,*) '#E_coll        ', sqrt(sh)
    write(iunit,*) '#'
    write(iunit,*) '#Phase space   ', _PhaseSpaceType
    write(iunit,*) '#'
    write(iunit,*) '#vegasNc0      ', vegasNc0
    write(iunit,*) '#vegasNc1      ', vegasNc1
    write(iunit,*) '#vegasIt0      ', vegasIt0
    write(iunit,*) '#vegasIt1      ', vegasIt1
    write(iunit,*) '#seed          ', seed
    write(iunit,*) '#'
    write(iunit,*) '#qmin          ', sqrt(myq2min)
    write(iunit,*) '#qmax          ', sqrt(myq2max)
    write(iunit,*) '#'
    write(iunit,*) '#mu0           ', sqrt(musq_in)
    write(iunit,*) '#mu0 could be overwritten in cut_histo by user'
    write(iunit,*) '#'
    write(iunit,*) '#xMuR          ', xMuR
    write(iunit,*) '#xMuF          ', xMuF
    write(iunit,*) '#'
    write(iunit,*) '#PDF set       ', pdffile
    write(iunit,*) '#nf            ', nf
    write(iunit,*) '#nf_pdf        ', nf_pdf
    write(iunit,*) '#'
    write(iunit,*) '#bufftc        ', bufftc
    write(iunit,*) '#buffrr        ', buffrr
    write(iunit,*) '#buffrr_pow    ', buffrr_pow
    write(iunit,*) '#buffrv        ', buffrv
    write(iunit,*) '#buff          ', buff
    write(iunit,*) '#'
    write(iunit,*) '#withmax       ', withmax




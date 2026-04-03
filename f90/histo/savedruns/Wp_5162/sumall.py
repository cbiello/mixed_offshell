runname="WpH_anom_fid"
chlist="qg gq gg ns s"
#chlist="qg"
scalelist="0.50 1.00 2.00"
#scalelist="1.00"
#coupllist="cdW_1.0_noSM cuW_1.0_noSM cHq3_1.0_noSM cdW_1.0_withSM cHud_1.0_noSM SM test1 test2"
coupllist="cdW_1.0_noSM cuW_1.0_noSM cHq3_1.0_noSM cdW_1.0_withSM cHud_1.0_noSM SM test1 test2"
#coupllist="SM"

## r0 = GamH / GammaHbb * BR(H->bb)
## with BR(H->bb) = 0.5824, dGam(H->bb)=1.92637679909396175582 for massive b quarks, GamH=4.165 MeV
r0=1.2592012  


rm *avg*histo
rm *summed*histo
rm *dNNLO*histo

for coupl in ${coupllist};do
    couplout=$coupl
    if [ $coupl == "cdW_1.0_withSM" ]; then
	couplout="cHq3_1.0_withSM"
    fi
    if [ $coupl == "test1" ]; then
	scalelist="1.00"
    elif [ $coupl == "test2" ]; then
	scalelist="1.00"
    else
	scalelist="0.50 1.00 2.00"
	#scalelist="1.00"
    fi
    for scale in ${scalelist};do
	
	for ch in ${chlist}; do
	    if [ $ch = "ns" ]; then
		seclist="rr4151a rr4151b rr4151c rr4151d rr4252a rr4252b rr4252c rr4252d rr4152 rr4251 vv rv sub124 sub12 subv schs schr schv"
		#		seclist="rr4151b rr4151c rr4151d rr4252a rr4252b rr4252c rr4252d rr4152 rr4251 rv sub124 sub12 subv schs schr schv vv"
	    elif [ $ch = "qg" ]; then
		seclist="rr4151a rr4151b rr4151c rr4151d rr4252a rr4252b rr4252c rr4252d rr4152 rr4251 subv sub124 sub12 rv schs schr"
		#		seclist="subv"
	    elif [ $ch = "gq" ]; then
		seclist="rr4151a rr4151b rr4151c rr4151d rr4252a rr4252b rr4252c rr4252d rr4152 rr4251 subv sub124 sub12 rv schs schr" 
	    elif [ $ch = "s" ]; then
		seclist="rr4151a rr4151b rr4151c rr4151d rr4252a rr4252b rr4252c rr4252d rr4152 rr4251"
	    elif [ $ch = "gg" ]; then
		seclist="sub124 rr sub12"
	    fi
	    
	    for sec in ${seclist}; do
		myfiles=${runname}_*_${coupl}_s*_nnlo_${sec}_${ch}_xMuR_${scale}_xMuF_${scale}.histo
		avgfile=${runname}_avg_${couplout}_nnlo_${sec}_${ch}_xMuR_${scale}_xMuF_${scale}.histo
		r0file=${runname}_avg_${couplout}_nnlo_r0_${sec}_${ch}_xMuR_${scale}_xMuF_${scale}.histo
		rrsummer avg ${myfiles} ${avgfile}
		rrsummer mul $r0 ${avgfile} ${r0file}
	    done
	    
	    myfiles=${runname}_avg_${couplout}_nnlo_r0_*_${ch}_xMuR_${scale}_xMuF_${scale}.histo
	    sumfiles=${runname}_summed_${couplout}_nnlo_${ch}_xMuR_${scale}_xMuF_${scale}.histo
	    rrsummer add ${myfiles} ${sumfiles}
	done
	#	rrsummer add ${runname}_summed_${couplout}_nnlo_gq_xMuR_${scale}_xMuF_${scale}.histo ${runname}_summed_${couplout}_nnlo_qg_xMuR_${scale}_xMuF_${scale}.histo ${runname}_summed_${couplout}_nnlo_gqPqg_xMuR_${scale}_xMuF_${scale}.histo
	#	rrsummer add ${runname}_summed_${couplout}_nnlo_ns_xMuR_${scale}_xMuF_${scale}.histo ${runname}_summed_${couplout}_nnlo_s_xMuR_${scale}_xMuF_${scale}.histo ${runname}_summed_${couplout}_nnlo_nsPs_xMuR_${scale}_xMuF_${scale}.histo
	
	myfiles=${runname}_avg_${couplout}_nnlo_r0_*_*_xMuR_${scale}_xMuF_${scale}.histo
	sumfiles=${runname}_${couplout}_dNNLO_xMuR_${scale}_xMuF_${scale}.histo
	rrsummer add ${myfiles} ${sumfiles}
	nlofile=../NLO/${runname}_${couplout}_NLO_xMuR_${scale}_xMuF_${scale}.histo
	nnlofile=${runname}_${couplout}_NNLO_xMuR_${scale}_xMuF_${scale}.histo
	rrsummer add ${sumfiles} ${nlofile} ${nnlofile}
    done
done

for coupl in ${coupllist};do
    couplout=$coupl
    if [ $coupl == "cdW_1.0_withSM" ]; then
	couplout="cHq3_1.0_withSM"
    fi
    for scale in ${scalelist};do
	sumfiles=${runname}_${couplout}_dNNLO_xMuR_${scale}_xMuF_${scale}.histo
	rrsummer sbn 1 ${sumfiles}
#	for ch in ${chlist}; do
#	    rrsummer sbn 1 ${runname}_summed_${couplout}_nnlo_${ch}_xMuR_${scale}_xMuF_${scale}.histo
#	done
    done
done


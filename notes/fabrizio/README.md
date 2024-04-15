- results_for_fabrizio_complete.nb: simplified results from Chiara and Federica
- subtraction_f.wl: Fabrizio's notebook
- comparison_against_chiara_and_federica.nb: compare the two notebook. Chiara's and Federica's results hard-typed, Fabrizio's notebook needs to be run

- to get the results:
grep 'Accumulated result' * | sed 's/PC_ATLAS_160601736_/res[/;s/_dynscale_xMuR_1.00_xMuF_1.00.histo:/]=/;s/_/,/g;s/# Accumulated result://;s/E/*10^/' > res-200-300.mat
for i in *.histo; do res=`sed -n '/^ 1    1/p' $i | awk '{print $4}'`; echo $i $res | sed 's/PC_ATLAS_160601736_/res[/;s/_dynscale_xMuR_1.00_xMuF_1.00.histo/]=/;s/_/,/g;s/# Accumulated result://;s/E/*10^/' >> res2.mat; done

- to get Raoul's QCD results:
for i in *.histo; do res=`sed -n '/^ 1    1/p' $i | awk '{print $4}'`; echo $i $res | sed 's/DY_QCD_az_gavin_/res[/;s/_avg_xMuR_1.00_xMuF_1.00.histo/]=/;s/_/,/g;s/# Accumulated result://;s/E/*10^/;s/nnlo/nnloqcd/' >> resQCD.mat; done

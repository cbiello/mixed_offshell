# gnuplot file, template generated automatically with
# /Users/fabrizio/tools/bin/gptemplate.py test_plots

set term pdfcairo enhanced font 'Helvetica,14pt' lw 1.5 size 10cm,8cm
set  style data lines
set lt 1 pt 4 ps 0.6 lw 1.5 lc rgb '#f00000'
set lt 2 pt 5 ps 0.6 lw 1.5 lc rgb '#00c000'
set lt 3 pt 6 ps 0.6 lw 1.5 lc rgb '#0000e0'
set lt 4 pt 7 ps 0.6 lw 1.5 lc rgb '#ff8000'
set lt 5 pt 8 ps 0.6 lw 1.5 lc rgb '#5070ff'
set key spacing 1.3

set output 'test_plots.pdf'



# put your code here

file = '../f90/histo/test_nnlo_ns_rr5161a_fixscale_xMuR_1.00_xMuF_1.00.histo'
plt  = '< mergeidx_mine.pl -fc -f '.file.' '

observables = "mll yll ptlm ptlp ylm ylp"
names = "m_{ll} y_{ll} p_{⊥,l^{-}} p_{⊥,l^{+}} y_{l^{-}} y_{l^{+}}"

do for [i=1:words(observables)] {

   plot plt.word(observables,i) u 3:4 w st lt 3 t word(names,i)

}	
   


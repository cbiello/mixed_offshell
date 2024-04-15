# gnuplot file, template generated automatically with
# /Users/fabrizio/tools/bin/gptemplate.py test_plots

set term pdfcairo enhanced font 'Helvetica,16pt' lw 1.5 size 10cm,8cm
set  style data lines
set lt 1 pt 4 ps 0.6 lw 1.5 lc rgb '#f00000'
set lt 2 pt 5 ps 0.6 lw 1.5 lc rgb '#00c000'
set lt 3 pt 6 ps 0.6 lw 1.5 lc rgb '#0000e0'
set lt 4 pt 7 ps 0.6 lw 1.5 lc rgb '#ff8000'
set lt 5 pt 8 ps 0.6 lw 1.5 lc rgb '#5070ff'
set key spacing 1.3

set output 'test_plots_lo.pdf'



# put your code here

myunit = 'fb'

file = '../f90/histo/test_lo_ns_lo_dynscale_xMuR_1.00_xMuF_1.00.histo'
plt  = '< mergeidx_mine.pl -fc -f '.file.' '

observables = "mll mll yll ptlm ptlm ptlp ptlm ylm ylp"
names =    "m_{ll} m_{ll} y_{ll} p_{⊥,l^{-}} p_{⊥,l^{-}} p_{⊥,l^{+}} p_{⊥,l^{+}} y_{l^{-}} y_{l^{+}}"
logscale = "0 1 0 0 1 0 1 0 0"
units = "[GeV] [GeV] 0 [GeV] [GeV] [GeV] [GeV] 0 0"

do for [i=1:words(observables)] {

   obs = word(names,i)

   if(word(logscale,i) eq "1") {set logscale y} else {unset logscale y}
   if(word(units,i) eq "0") {set xlabel obs} else {set xlabel obs.' '.word(units,i)}

   ytitle = "dσ/d".obs
   set ylabel ytitle

   plot plt.word(observables,i) u 3:4 w st lt 3 t 'LO'

}	
   


reset
set term post eps enhanced color "Helvetica" 18


xsecA=94.603114210973331;
xsecB=111.17682659786097;
xsecC=169.16502871531844;
xsecD=47.625489850834896;
fileA="WH_basiccuts_m20220609_s25258_couplA_lo_xMuR_1.00_xMuF_1.00.histo"
fileB="WH_basiccuts_m20220609_s25258_couplB_lo_xMuR_1.00_xMuF_1.00.histo"
fileC="WH_basiccuts_m20220609_s25258_couplC_lo_xMuR_1.00_xMuF_1.00.histo"
fileD="WH_basiccuts_m20220609_s25258_couplD_lo_xMuR_1.00_xMuF_1.00.histo"
unset yrange;	
#### m(VH)

xsecA=52.386858895490327;
xsecB=63.994883530776747;
xsecC=116.86668193409521;
xsecD=35.504369031529748;

set output "mVH.eps"
NHi=2
set xrange[0:950]
set xlabel "m(VH)"
set ylabel "d{/Symbol s}/dm(VH) [fb/GeV]"

plot fileA u ($3):($1==NHi ? ($4/xsecA):1/0) w steps ti "A", \
     fileB u ($3):($1==NHi ? ($4/xsecB):1/0) w steps ti "B", \
     fileC u ($3):($1==NHi ? ($4/xsecC):1/0) w steps ti "C", \
     fileD u ($3):($1==NHi ? ($4/xsecD):1/0) w steps ti "D"


#### pT(H)

xsecA=58.393999701317448;
xsecB=71.286370868924735;
xsecC=125.62422218471278;
xsecD=37.456820263814315;



set output "pTH.eps"
NHi=3
set xrange[0:450]
set xlabel "pT(H)"
set ylabel "d{/Symbol s}/dpT(H) [fb/GeV]"

plot fileA u ($3):($1==NHi ? ($4/xsecA):1/0) w steps ti "A", \
     fileB u ($3):($1==NHi ? ($4/xsecB):1/0) w steps ti "B", \
     fileC u ($3):($1==NHi ? ($4/xsecC):1/0) w steps ti "C", \
     fileD u ($3):($1==NHi ? ($4/xsecD):1/0) w steps ti "D"


#### y(H)

xsecA=94.603114210973331;
xsecB=111.17682659786097;
xsecC=169.16502871531844;
xsecD=47.625489850834896;

set output "yH.eps"
NHi=4
set xrange[-5:5]
set xlabel "y(H)"
set ylabel "d{/Symbol s}/dy(H) [fb]"

plot fileA u ($3):($1==NHi ? ($4/xsecA):1/0) w steps ti "A", \
     fileB u ($3):($1==NHi ? ($4/xsecB):1/0) w steps ti "B", \
     fileC u ($3):($1==NHi ? ($4/xsecC):1/0) w steps ti "C", \
     fileD u ($3):($1==NHi ? ($4/xsecD):1/0) w steps ti "D"

#### pT(l)

xsecA=90.777723623600110;
xsecB=84.564354142725222;
xsecC=154.97017780790537;
xsecD=44.226402827972905;

set output "pTl.eps"
NHi=6
set xrange[0:450]
set xlabel "pT(l)"
set ylabel "d{/Symbol s}/dpT(l) [fb/GeV]"

plot fileA u ($3):($1==NHi ? ($4/xsecA):1/0) w steps ti "A", \
     fileB u ($3):($1==NHi ? ($4/xsecB):1/0) w steps ti "B", \
     fileC u ($3):($1==NHi ? ($4/xsecC):1/0) w steps ti "C", \
     fileD u ($3):($1==NHi ? ($4/xsecD):1/0) w steps ti "D"

#### y(l)

xsecA=93.919769881857832;
xsecB=110.39120122172012;
xsecC=167.56832848200335;
xsecD=47.128462874288289;

set output "yl.eps"
NHi=7
set xrange[-2.5:2.5]
set xlabel "y(l)"
set ylabel "d{/Symbol s}/dy(l) [fb]"

plot fileA u ($3):($1==NHi ? ($4/xsecA):1/0) w steps ti "A", \
     fileB u ($3):($1==NHi ? ($4/xsecB):1/0) w steps ti "B", \
     fileC u ($3):($1==NHi ? ($4/xsecC):1/0) w steps ti "C", \
     fileD u ($3):($1==NHi ? ($4/xsecD):1/0) w steps ti "D"

#### pT(miss)

xsecA=71.542660679116210;
xsecB=107.35471364280585;
xsecC=155.12924374969214;
xsecD=44.263001441738034;

set output "pTmiss.eps"
NHi=8
set xrange[0:450]
set xlabel "pT(l)"
set ylabel "d{/Symbol s}/dpT(miss) [fb/GeV]"

plot fileA u ($3):($1==NHi ? ($4/xsecA):1/0) w steps ti "A", \
     fileB u ($3):($1==NHi ? ($4/xsecB):1/0) w steps ti "B", \
     fileC u ($3):($1==NHi ? ($4/xsecC):1/0) w steps ti "C", \
     fileD u ($3):($1==NHi ? ($4/xsecD):1/0) w steps ti "D"

#### m(VH_rec)


set output "mVH_rec.eps"
NHi=9
set xrange[0:950]
set xlabel "m(VH rec)"
set ylabel "d{/Symbol s}/dm(VH rec) [fb/GeV]"

plot fileA u ($3):($1==NHi ? ($4/xsecA):1/0) w steps ti "A", \
     fileB u ($3):($1==NHi ? ($4/xsecB):1/0) w steps ti "B", \
     fileC u ($3):($1==NHi ? ($4/xsecC):1/0) w steps ti "C", \
     fileD u ($3):($1==NHi ? ($4/xsecD):1/0) w steps ti "D"


#### pT(H_rec)

xsecA=59.451659987524259;
xsecB=71.190980599864929;
xsecC=113.82975857170008;
xsecD=32.053142755873196;

set output "pTH_rec.eps"
NHi=11
set xrange[0:450]
set xlabel "pT(H rec)"
set ylabel "d{/Symbol s}/dpT(H rec) [fb/GeV]"

plot fileA u ($3):($1==NHi ? ($4/xsecA):1/0) w steps ti "A", \
     fileB u ($3):($1==NHi ? ($4/xsecB):1/0) w steps ti "B", \
     fileC u ($3):($1==NHi ? ($4/xsecC):1/0) w steps ti "C", \
     fileD u ($3):($1==NHi ? ($4/xsecD):1/0) w steps ti "D"


#### y(H_rec)

xsecA=60.480747586602909;
xsecB=72.284421090650724
xsecC=114.97940824754329;
xsecD=32.340651595081745;

set output "yH_rec.eps"
NHi=12
set xrange[-2.5:2.5]
set xlabel "y(H rec)"
set ylabel "d{/Symbol s}/dy(H rec) [fb]"

plot fileA u ($3):($1==NHi ? ($4/xsecA):1/0) w steps ti "A", \
     fileB u ($3):($1==NHi ? ($4/xsecB):1/0) w steps ti "B", \
     fileC u ($3):($1==NHi ? ($4/xsecC):1/0) w steps ti "C", \
     fileD u ($3):($1==NHi ? ($4/xsecD):1/0) w steps ti "D"

#### pT(b-jet)

xsecA=52.487832437686393;
xsecB=63.419019615123943;
xsecC=105.24577862647050;
xsecD=30.078264489954677;

set output "pTbjet.eps"
NHi=14
set xrange[0:450]
set xlabel "pT(b-jet)"
set ylabel "d{/Symbol s}/dpT(b-jet) [fb/GeV]"

plot fileA u ($3):($1==NHi ? ($4/xsecA):1/0) w steps ti "A", \
     fileB u ($3):($1==NHi ? ($4/xsecB):1/0) w steps ti "B", \
     fileC u ($3):($1==NHi ? ($4/xsecC):1/0) w steps ti "C", \
     fileD u ($3):($1==NHi ? ($4/xsecD):1/0) w steps ti "D"



#### y(b-jet)

xsecA=60.480881379986435;
xsecB=72.284380681124389;
xsecC=114.97970485482213;
xsecD=32.340681832958580;

set output "ybjet.eps"
NHi=15
set xrange[-2.5:2.5]
set xlabel "y(b-jet)"
set ylabel "d{/Symbol s}/dy(b-jet) [fb]"

plot fileA u ($3):($1==NHi ? ($4/xsecA):1/0) w steps ti "A", \
     fileB u ($3):($1==NHi ? ($4/xsecB):1/0) w steps ti "B", \
     fileC u ($3):($1==NHi ? ($4/xsecC):1/0) w steps ti "C", \
     fileD u ($3):($1==NHi ? ($4/xsecD):1/0) w steps ti "D"

#### DeltaR(bjet-l)

xsecA=48.217903789746572;
xsecB=53.301651376460143;
xsecC=90.481211535543139;
xsecD=26.299030335892240;

set output "DeltaR_bl.eps"
NHi=19
set xrange[0:3.1]
set xlabel "DeltaR(bl)"
set ylabel "d{/Symbol s}/dDeltaR(b-l) [fb]"

plot fileA u ($3):($1==NHi ? ($4/xsecA):1/0) w steps ti "A", \
     fileB u ($3):($1==NHi ? ($4/xsecB):1/0) w steps ti "B", \
     fileC u ($3):($1==NHi ? ($4/xsecC):1/0) w steps ti "C", \
     fileD u ($3):($1==NHi ? ($4/xsecD):1/0) w steps ti "D"

#### DeltaR(bb)

xsecA=59.183396035568911;
xsecB=70.945315121513730;
xsecC=112.37550068562962;
xsecD=31.059655342015230;

set output "DeltaR_bb.eps"
NHi=20
set xrange[0:3.1]
set xlabel "DeltaR(bb)"
set ylabel "d{/Symbol s}/dDeltaR(bb) [fb]"

plot fileA u ($3):($1==NHi ? ($4/xsecA):1/0) w steps ti "A", \
     fileB u ($3):($1==NHi ? ($4/xsecB):1/0) w steps ti "B", \
     fileC u ($3):($1==NHi ? ($4/xsecC):1/0) w steps ti "C", \
     fileD u ($3):($1==NHi ? ($4/xsecD):1/0) w steps ti "D"

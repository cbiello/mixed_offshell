from plotter_settings import *

################################################################################
initializetex()

################################################################################
#>> output file
pltname='ptl'

#>> index
NHi=6


#>> legend
lab0=r"${\rm d}\sigma_{\rm SM}$"
lab1=r"${\rm d}\sigma_{UW}$"
lab2=r"${\rm d}\sigma_{DW}$"
lab3=r"${\rm d}\sigma_{HQ3}$"
lab4=r"${\rm d}\sigma_{I}$"
lab5=r"${\rm d}\sigma_{HUD}$"


#>> colors:
col0=colazure
col1=colmagenta
col2=colorange
col3=colplum
col4=colgoldenrod
col5=colgreen
col6=colBlack

#>> axes labels
XLABEL0=r"$p_{T,\ell}  {\rm [ GeV ]} $"
YLABEL1=r"$1/\sigma_{\rm NNLO} \,\, {\rm d}\sigma_{\rm NNLO}/{\rm d}p_{T,\ell}$ [ 1/GeV ]"
YLABEL2=r"NLO/LO"
YLABEL3=r"NNLO/NLO"
#>> xlimit
XMIN0=20
XMAX0=200

YMIN0=3e-5
YMAX0=0.02

YMIN1=0.8
YMAX1=1.7

YMIN2=0.9
YMAX2=1.13

my_minor_xticks=numpy.arange(0,1000,20)
my_minor_yticks1=numpy.arange(0.0,5.0,0.25)

################################################################################
#output='figs/'+pltname+'.pdf'
output='./'+pltname+'.pdf'

#>> load data
### SM
xx,val_SM_LO,err_SM_LO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_SM_LO_xMuR_1.00_xMuF_1.00_ptl.histo')
xx,val_SM_NLO,err_SM_NLO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_SM_NLO_xMuR_1.00_xMuF_1.00_ptl.histo')
xx,val_SM_NNLO,err_SM_NNLO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_SM_NNLO_xMuR_1.00_xMuF_1.00_ptl.histo')
### cuWsq
xx,val_cuWsq_LO,err0_cuWsq_LO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cuW_1.0_noSM_LO_xMuR_1.00_xMuF_1.00_ptl.histo')
xx,val_cuWsq_NLO,err0_cuWsq_NLO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cuW_1.0_noSM_NLO_xMuR_1.00_xMuF_1.00_ptl.histo')
xx,val_cuWsq_NNLO,err0_cuWsq_NNLO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cuW_1.0_noSM_NNLO_xMuR_1.00_xMuF_1.00_ptl.histo')
### cdWsq
xx,val_cdWsq_LO,err0_cdWsq_LO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cdW_1.0_noSM_LO_xMuR_1.00_xMuF_1.00_ptl.histo')
xx,val_cdWsq_NLO,err0_cdWsq_NLO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cdW_1.0_noSM_NLO_xMuR_1.00_xMuF_1.00_ptl.histo')
xx,val_cdWsq_NNLO,err0_cdWsq_NNLO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cdW_1.0_noSM_NNLO_xMuR_1.00_xMuF_1.00_ptl.histo')

### cHq3^2
xx,val_cHq3sq_LO,err0_cHq3sq_LO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cHq3_1.0_noSM_LO_xMuR_1.00_xMuF_1.00_ptl.histo')
xx,val_cHq3sq_NLO,err0_cHq3sq_NLO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cHq3_1.0_noSM_NLO_xMuR_1.00_xMuF_1.00_ptl.histo')
xx,val_cHq3sq_NNLO,err0_cHq3sq_NNLO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cHq3_1.0_noSM_NNLO_xMuR_1.00_xMuF_1.00_ptl.histo')

### cHq3
xx,val_cHq3_LO,err0_cHq3_LO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cHq3Intf_LO_xMuR_1.00_xMuF_1.00_ptl.histo')
xx,val_cHq3_NLO,err0_cHq3_NLO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cHq3Intf_NLO_xMuR_1.00_xMuF_1.00_ptl.histo')
xx,val_cHq3_NNLO,err0_cHq3_NNLO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cHq3Intf_NNLO_xMuR_1.00_xMuF_1.00_ptl.histo')

### cHud^2
xx,val_cHudsq_LO,err0_cHudsq_LO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cHud_1.0_noSM_LO_xMuR_1.00_xMuF_1.00_ptl.histo')
xx,val_cHudsq_NLO,err0_cHudsq_NLO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cHud_1.0_noSM_NLO_xMuR_1.00_xMuF_1.00_ptl.histo')
xx,val_cHudsq_NNLO,err0_cHudsq_NNLO,bins=load_histo(NHi,'../data/nf_5/fid/prepare_for_plotting/WpH_Hundk_fid_cHud_1.0_noSM_NNLO_xMuR_1.00_xMuF_1.00_ptl.histo')



## normalize all NNLO results to 1
xsecSM=58.60
xseccuWsq=83.17
xseccdWsq=55.17
xseccHq3sq=70.15
xseccHq3=106.80
xseccHudsq=19.20


val_SM_NNLO_resc = val_SM_NNLO/xsecSM
val_cuWsq_NNLO_resc=val_cuWsq_NNLO/xseccuWsq
val_cdWsq_NNLO_resc=val_cdWsq_NNLO/xseccdWsq
val_cHq3sq_NNLO_resc=val_cHq3sq_NNLO/xseccHq3sq
val_cHq3_NNLO_resc=val_cHq3_NNLO/xseccHq3
val_cHudsq_NNLO_resc=val_cHudsq_NNLO/xseccHudsq

### central value
#val0=val0a
#val1=val1a
#val2=val2a
#val3=val3a
#val4=val4a

binsize=(xx[2]-xx[1]) ###(bins[1]-bins[0])
print(binsize)
val_SM_NNLO_resc=val_SM_NNLO_resc/binsize
val_cuWsq_NNLO_resc=val_cuWsq_NNLO_resc/binsize
val_cdWsq_NNLO_resc=val_cdWsq_NNLO_resc/binsize
val_cHq3sq_NNLO_resc=val_cHq3sq_NNLO_resc/binsize
val_cHq3_NNLO_resc=val_cHq3_NNLO_resc/binsize
val_cHudsq_NNLO_resc=val_cHudsq_NNLO_resc/binsize


################################################################################
fig_dims=[]
fig_dims=figsize(fig_dims)
fig = plt.figure(figsize=fig_dims)

pdf_pages = PdfPages(output)

################################################################################
#>> panels
#old ax1=plt.subplot2grid((5,1),(0,0),rowspan=3)
#old ax2=plt.subplot2grid((5,1),(3,0),rowspan=2)
gs = fig.add_gridspec(3, hspace=0,height_ratios=[3,1,1])
axs = gs.subplots(sharex=True, sharey=False)
ax1=axs[0]
ax2=axs[1]
ax3=axs[2]


#>> some settings
plt.setp(ax1.get_xticklabels(),visible=False)
plt.setp(ax2.get_xticklabels(),visible=True)

ax1.set_ylabel(YLABEL1)
ax1.set_xticks(my_minor_xticks,minor=True)
ax1.grid(dashes=[10,10],linewidth=0.1,which='both')
#ax1.set_yscale('log')
ax1.tick_params(direction='in')

ax2.set_ylabel(YLABEL2)
ax2.set_xticks(my_minor_xticks,minor=True)
ax2.set_yticks(my_minor_yticks1,minor=True)
ax2.grid(dashes=[10,10],linewidth=0.1,which='both')
ax2.tick_params(direction='in')

ax3.set_ylabel(YLABEL3)
ax3.set_xticks(my_minor_xticks,minor=True)
ax3.set_yticks(my_minor_yticks1,minor=True)
ax3.grid(dashes=[10,10],linewidth=0.1,which='both')
ax3.tick_params(direction='in')


#>> plot things (panel [1],linewidth=1.0)
p0a,= ax1.step(xx,val_SM_NNLO_resc,where='post',color=col0[1],linewidth=1.0)
p1a,= ax1.step(xx,val_cuWsq_NNLO_resc,where='post',color=col1[1],linewidth=1.0)
p2a,= ax1.step(xx,val_cdWsq_NNLO_resc,where='post',color=col2[1],linewidth=1.0)
p3a,= ax1.step(xx,val_cHq3sq_NNLO_resc,where='post',color=col3[1],linewidth=1.0)
p4a,= ax1.step(xx,val_cHq3_NNLO_resc,where='post',color=col4[1],linewidth=1.0)
p5a,= ax1.step(xx,val_cHudsq_NNLO_resc,where='post',color=col5[1],linewidth=1.0)

p0b,= [plt.Rectangle((0,0),0,0,facecolor='none',linewidth=0.2,edgecolor=col0[0],alpha=None)]
p1b,= [plt.Rectangle((0,0),0,0,facecolor='none',linewidth=0.2,edgecolor=col1[0],alpha=None)]
p2b,= [plt.Rectangle((0,0),0,0,facecolor='none',linewidth=0.2,edgecolor=col2[0],alpha=None)]
p3b,= [plt.Rectangle((0,0),0,0,facecolor='none',linewidth=0.2,edgecolor=col3[0],alpha=None)]
p4b,= [plt.Rectangle((0,0),0,0,facecolor='none',linewidth=0.2,edgecolor=col4[0],alpha=None)]
p5b,= [plt.Rectangle((0,0),0,0,facecolor='none',linewidth=0.2,edgecolor=col5[0],alpha=None)]
#lines = ((p0a,p0b))
lines = ((p0a,p0b),(p1a,p1b),(p2a,p2b),(p3a,p3b),(p4a,p4b),(p5a,p5b))
legend=(lab0,lab1,lab2,lab3,lab4,lab5)

ax1.set_xlim(xmin=XMIN0, xmax=XMAX0)
ax1.set_ylim(ymin=YMIN0, ymax=YMAX0)

#>> plot things (panel 2)
ax2.step(xx, myratio(val_SM_NLO,val_SM_LO), where='post',color=col0[1],linewidth=1.0)
ax2.step(xx, myratio(val_cuWsq_NLO,val_cuWsq_LO), where='post',color=col1[1],linewidth=1.0)
ax2.step(xx, myratio(val_cdWsq_NLO,val_cdWsq_LO), where='post',color=col2[1],linewidth=1.0)
ax2.step(xx, myratio(val_cHq3sq_NLO,val_cHq3sq_LO), where='post',color=col3[1],linewidth=1.0)
ax2.step(xx, myratio(val_cHq3_NLO,val_cHq3_LO), where='post',color=col4[1],linewidth=1.0)
ax2.step(xx, myratio(val_cHudsq_NLO,val_cHudsq_LO), where='post',color=col5[1],linewidth=1.0)
ax2.step(xx, myratio(val_cHudsq_NLO,val_cHudsq_NLO), where='post',color=col6[1],linewidth=1.0,linestyle='dashed')  # this is just to plot 1
ax2.set_xlim(xmin=XMIN0, xmax=XMAX0)
ax2.set_ylim(ymin=YMIN1, ymax=YMAX1)


#>> plot things (panel 2)
ax3.step(xx, myratio(val_SM_NNLO,val_SM_NLO), where='post',color=col0[1],linewidth=1.0)
ax3.step(xx, myratio(val_cuWsq_NNLO,val_cuWsq_NLO), where='post',color=col1[1],linewidth=1.0)
ax3.step(xx, myratio(val_cdWsq_NNLO,val_cdWsq_NLO), where='post',color=col2[1],linewidth=1.0)
ax3.step(xx, myratio(val_cHq3sq_NNLO,val_cHq3sq_NLO), where='post',color=col3[1],linewidth=1.0)
ax3.step(xx, myratio(val_cHq3_NNLO,val_cHq3_NLO), where='post',color=col4[1],linewidth=1.0)
ax3.step(xx, myratio(val_cHudsq_NNLO,val_cHudsq_NLO), where='post',color=col5[1],linewidth=1.0)
ax3.step(xx, myratio(val_cHudsq_NNLO,val_cHudsq_NNLO), where='post',color=col6[1],linewidth=1.0,linestyle='dashed')  # this is just to plot 1
ax3.set_xlim(xmin=XMIN0, xmax=XMAX0)
ax3.set_ylim(ymin=YMIN2, ymax=YMAX2)

#>> xlabel
plt.xlabel(XLABEL0)

#>> legend
leg=ax1.legend(lines,legend,loc='upper right')
##
rect = leg.get_frame()
rect.set_linewidth(0.0)

################################################################################
#>> put any labels in here
#ax1.annotate(r'150 GeV $< p_{t,W} <$ 250~GeV', xy=(0.490,0.90), xycoords="axes fraction", color="#000033")
ax1.annotate(r'LHC 13 TeV, NNLO QCD',                       xy=(0.330,0.9), xycoords="axes fraction", color="#ff0000")
#ax1.annotate(r'\bf $\rightarrow$ jet radius: R=0.4',xy=(0.55,0.82),xycoords="axes fraction")

################################################################################
#>> save figure
pdf_pages.savefig(fig,bbox_inches='tight')
pdf_pages.close()


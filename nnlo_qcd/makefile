Here = $(PWD)
Exec = ./nnlo_wp_qp

ObjectDir = $(Here)/objects
ModuleDir = $(Here)/modules
VegasDir=$(Here)/vegas_lepage

include identity.mk

Opt = Yes

ConfigFile= $(Here)/config_wp.cfg

##################################################################################################

ifeq ($(Iam),fabrizio)
   #
   F95compiler = gfortran -cpp @$(ConfigFile)
   ifeq ($(Opt),Yes)
     fcompOpts = -O2 -I $(ModuleDir) -ffree-line-length-1000
   else
     fcompOpts = -O0 -I $(ModuleDir) -ffree-line-length-1000
   endif
   #
   SEDcmd=sed
   LHAPDF=/Users/fabrizio/tools/lib

else ifeq ($(Iam),kirill)
   #
   F95compiler = ifort @$(ConfigFile)
   ifeq ($(Opt),Yes)
     fcompOpts  = -O2 -g -fpp -vec-report0 -module $(ModuleDir)
   else
     fcompOpts  = -O0 -fpp -nozero -debug all -check uninit -check bounds -check pointer -ftrapuv -module $(ModuleDir)
   endif
   #
   SEDcmd=gsed
   LHAPDF=/Users/kirill/u1/LHAPDF/lib

else ifeq ($(Iam),raoul)
   #
   F95compiler = ifort  @$(ConfigFile)
   ifeq ($(Opt),Yes)
     fcompOpts  = -O2 -g -fpp -vec-report0 -module $(ModuleDir)
   else
     fcompOpts  = -O0 -fpp -nozero -debug all -check uninit -check bounds -check pointer -ftrapuv -module $(ModuleDir)
   endif
   #
   SEDcmd=sed
   LHAPDF=/x2/roentsch/Lib/Install/LHAPDF-6.1.5/lib64

else ifeq ($(Iam),raoul-laptop)
   #
   F95compiler = gfortran -cpp @$(ConfigFile)
   ifeq ($(Opt),Yes)
     fcompOpts = -O2 -I $(ModuleDir) -ffree-line-length-1000
   else
     fcompOpts = -O0 -I $(ModuleDir) -ffree-line-length-1000
   endif
   #
   SEDcmd=sed
   LHAPDF=/home/raoul/Lib/Install/LHAPDF-6.1.5/lib

else ifeq ($(Iam),lxplus-fabrizio)
   #
   F95compiler = ifort @$(ConfigFile)
   ifeq ($(Opt),Yes)
     fcompOpts  = -O2 -g -fpp -vec-report0 -module $(ModuleDir)
   else
     fcompOpts  = -O0 -fpp -nozero -debug all -check uninit -check bounds -check pointer -ftrapuv -module $(ModuleDir)
   endif
   #
   SEDcmd=sed
   LHAPDF=/afs/cern.ch/user/f/facaola/tools/lib

else ifeq ($(Iam),lxplus)

   F95compiler = ifort @$(ConfigFile)
   ifeq ($(Opt),Yes)
     fcompOpts  = -O2 -g -fpp -vec-report0 -module $(ModuleDir)
   else
     fcompOpts  = -O0 -fpp -nozero -debug all -check uninit -check bounds -check pointer -ftrapuv -module $(ModuleDir)
   endif	

   F95compiler = gfortran -cpp @$(ConfigFile)
   ifeq ($(Opt),Yes)
     fcompOpts = -O2 -I $(ModuleDir) -ffree-line-length-1000
   else
     fcompOpts = -O0 -I $(ModuleDir) -ffree-line-length-1000
   endif
   #
   SEDcmd=sed
#   LHAPDF=/afs/cern.ch/user/f/facaola/tools/lib
   LHAPDF=/afs/cern.ch/user/r/rarontsc/Lib/Install/Lxplus/LHAPDF-6.2.0/lib
#   LHAPDF=/afs/cern.ch/user/r/rarontsc/Lib/Install/Lxplus/LHAPDF-6.2.0/lib

else ifeq ($(Iam),ttpmarvin)
   #-- KIT TTP (ttpmarvin and the SLURM nodes): make Iam=ttpmarvin
   F95compiler = gfortran -cpp @$(ConfigFile)
   ifeq ($(Opt),Yes)
     fcompOpts = -O2 -I $(ModuleDir) -ffree-line-length-1000
   else
     fcompOpts = -O0 -I $(ModuleDir) -ffree-line-length-1000
   endif
   #
   SEDcmd=sed
   LHAPDF=/users/ttp/signorile/Desktop/tools/lhapdf-install/lib64
   fcompOpts += -Wl,-rpath=$(LHAPDF)

endif

##################################################################################################

fcomp = $(F95compiler) $(fcompOpts) -I$(VegasDir) -L$(LHAPDF) -lLHAPDF 

fcomp4main = $(fcomp) -o $(Exec) $(ALLCODE)

##################################################################################################

BasicDef = \
      ./BasicDef/mod_types.f90 \
      ./BasicDef/common_def.f90 \
      ./BasicDef/mod_consts_dp.f90 \
      ./BasicDef/mod_consts_debug.f90 \

BasicDefObj = \
      $(ObjectDir)/mod_types.o \
      $(ObjectDir)/common_def.o \
      $(ObjectDir)/mod_consts_dp.o \
      $(ObjectDir)/mod_consts_debug.o \

BasicDefInc = \
      ./BasicDef/inc_procconsts_dp.f90 \

AuxFunctions = \
      ./AuxFunctions/cpolylog_main.F \
      ./AuxFunctions/mod_auxfunctions.f90 \
      ./AuxFunctions/mod_procfunctions.f90 \
      ./AuxFunctions/mod_kinfunctions.f90 \
      ./AuxFunctions/mod_lumi.f90 \
      ./AuxFunctions/mod_auxlumi_wp.f90 \

AuxFunctionsObj = \
      $(ObjectDir)/cpolylog_main.o \
      $(ObjectDir)/mod_auxfunctions.o \
      $(ObjectDir)/mod_procfunctions.o \
      $(ObjectDir)/mod_kinfunctions.o \
      $(ObjectDir)/mod_lumi.o \
      $(ObjectDir)/mod_auxlumi_wp.o \

GenProc = \
      ./GenProc/mod_process.f90 \
      ./GenProc/mod_allocate_arrays.f90 \

GenProcObj = \
      $(ObjectDir)/mod_process.o \
      $(ObjectDir)/mod_allocate_arrays.o \

Histograms = \
      ./Histograms/mod_histograms.f90 \

HistogramsObj = \
      $(ObjectDir)/mod_histograms.o \

HistogramsInc = \
      ./Histograms/cut_histo_w.f90 \
      ./Histograms/check_sumbin.f90 \
      ./Histograms/prepare_and_write.f90 \

KinematicsGen = \
      ./Kinematics/mod_kinematics_gen.f90 \
      ./Kinematics/mod_kinematics_lept.f90 \

KinematicsGenObj = \
      $(ObjectDir)/mod_kinematics_gen.o \
      $(ObjectDir)/mod_kinematics_lept.o \

KinematicsLO = \
      ./Kinematics/mod_kinematicsLO.f90 \

KinematicsLOObj = \
      $(ObjectDir)/mod_kinematicsLO.o \

KinematicsNLO = \
      ./Kinematics/mod_kinematicsNLO.f90 \
       ./Kinematics/mod_kinematicsNLOz.f90 \

KinematicsNLOObj = \
      $(ObjectDir)/mod_kinematicsNLO.o \
       $(ObjectDir)/mod_kinematicsNLOz.o \

# KinematicsNNLO2q = \
#       ./Kinematics/mod_kinematicsNNLO_2q.f90 \

# KinematicsNNLO2qObj = \
#       $(ObjectDir)/mod_kinematicsNNLO_2q.o \

# KinematicsNNLOgg = \
#       ./Kinematics/mod_kinematicsNNLO_gg.f90 \

# KinematicsNNLOggObj = \
#       $(ObjectDir)/mod_kinematicsNNLO_gg.o \

KinematicsNNLOEordAord = \
      ./Kinematics/mod_kinematicsNNLO_EordAord.f90 \

KinematicsNNLOEordAordObj = \
      $(ObjectDir)/mod_kinematicsNNLO_EordAord.o \

KinematicsNNLOEunordAord = \
      ./Kinematics/mod_kinematicsNNLO_EunordAord.f90 \

KinematicsNNLOEunordAordObj = \
      $(ObjectDir)/mod_kinematicsNNLO_EunordAord.o \

KinematicsNNLOEunordAunord = \
      ./Kinematics/mod_kinematicsNNLO_EunordAunord.f90 \

KinematicsNNLOEunordAunordObj = \
      $(ObjectDir)/mod_kinematicsNNLO_EunordAunord.o \

Amplitudes = \
      ./Amplitudes/mod_amplitudes_tree_wp.f90 \
      ./Amplitudes/mod_amplitudes_loop_wp.f90 \

AmplitudesObj = \
      $(ObjectDir)/mod_amplitudes_tree_wp.o \
      $(ObjectDir)/mod_amplitudes_loop_wp.o \

AmplitudesInc = \

SingLimits = \
	./SingLimits/mod_splitting.f90 \
	./SingLimits/mod_splitting_aux.f90 \
	./SingLimits/mod_limits_dy_ns.f90 \
	./SingLimits/mod_limits_dy_qg.f90 \
	./SingLimits/mod_limits_dy_gq.f90 \
	./SingLimits/mod_limits_rv_dy_ns.f90 \
	./SingLimits/mod_limits_rv_dy_qg.f90 \

SingLimitsObj = \
	$(ObjectDir)/mod_splitting.o \
	$(ObjectDir)/mod_splitting_aux.o \
	$(ObjectDir)/mod_limits_dy_ns.o \
	$(ObjectDir)/mod_limits_dy_qg.o \
	$(ObjectDir)/mod_limits_dy_gq.o \
	$(ObjectDir)/mod_limits_rv_dy_ns.o \
	$(ObjectDir)/mod_limits_rv_dy_qg.o \

 TCSplit = \
       ./TCSplit/mod_kinematics_tc_ord.f90 \
       ./TCSplit/mod_kinematics_tc_unord.f90 \
       ./TCSplit/mod_tripleclimits_q_gg_qstar.f90 \
       ./TCSplit/mod_my_tripleclimits_q_qqb_qstar.f90 \
       ./TCSplit/mod_tripleclimits_g_qg_qstar.f90 \

 TCSplitObj = \
       $(ObjectDir)/mod_kinematics_tc_ord.o \
       $(ObjectDir)/mod_kinematics_tc_unord.o \
       $(ObjectDir)/mod_tripleclimits_q_gg_qstar.o \
       $(ObjectDir)/mod_my_tripleclimits_q_qqb_qstar.o \
       $(ObjectDir)/mod_tripleclimits_g_qg_qstar.o \

 IntSub = \
       ./IntSub/mod_subtrfn_nlo_dy.f90 \
       ./IntSub/mod_subtrfn_nnlo_dy_ns.f90 \
       ./IntSub/mod_subtrfn_nnlo_dy_s.f90 \
       ./IntSub/mod_subtrfn_nnlo_dy_gq.f90 \
       ./IntSub/mod_subtrfn_nnlo_dy_gg.f90 \
#
 IntSubObj = \
       $(ObjectDir)/mod_subtrfn_nlo_dy.o \
       $(ObjectDir)/mod_subtrfn_nnlo_dy_ns.o \
       $(ObjectDir)/mod_subtrfn_nnlo_dy_s.o \
       $(ObjectDir)/mod_subtrfn_nnlo_dy_gq.o \
       $(ObjectDir)/mod_subtrfn_nnlo_dy_gg.o \

SectorsLO = \
      ./Sectors/mod_xsects_lo.f90 \

SectorsLOObj = \
      $(ObjectDir)/mod_xsects_lo.o \

SectorsNLO = \
      ./Sectors/mod_xsects_nlo_wp.f90 \

SectorsNLOObj = \
      $(ObjectDir)/mod_xsects_nlo_wp.o \

SectorsNNLOsub = \
       ./Sectors/mod_xsects_nnlo_sub_wp_qqb.f90 \
       ./Sectors/mod_xsects_nnlo_sub_wp_gq.f90 \
       ./Sectors/mod_xsects_nnlo_sub_wp_qg.f90 \
       ./Sectors/mod_xsects_nnlo_sub_wp_gg.f90 \


 SectorsNNLOsubObj = \
       $(ObjectDir)/mod_xsects_nnlo_sub_wp_qqb.o \
       $(ObjectDir)/mod_xsects_nnlo_sub_wp_gq.o \
       $(ObjectDir)/mod_xsects_nnlo_sub_wp_qg.o \
       $(ObjectDir)/mod_xsects_nnlo_sub_wp_gg.o \

SectorsNNLOrrNS = \
       ./Sectors/mod_xsects_nnlo_rr_wp_ns.f90 \

SectorsNNLOrrNSObj = \
       $(ObjectDir)/mod_xsects_nnlo_rr_wp_ns.o \

SectorsNNLOrrQG = \
       ./Sectors/mod_xsects_nnlo_rr_wp_qg.f90 \
       ./Sectors/mod_xsects_nnlo_rr_wp_gq.f90 \

SectorsNNLOrrQGObj = \
       $(ObjectDir)/mod_xsects_nnlo_rr_wp_qg.o \
       $(ObjectDir)/mod_xsects_nnlo_rr_wp_gq.o \

SectorsNNLOrrS = \
       ./Sectors/mod_xsects_nnlo_rr_wp_s.f90 \

SectorsNNLOrrSObj = \
       $(ObjectDir)/mod_xsects_nnlo_rr_wp_s.o \

SectorsNNLOrrGG = \
       ./Sectors/mod_xsects_nnlo_rr_wp_gg.f90 \

SectorsNNLOrrGGObj = \
       $(ObjectDir)/mod_xsects_nnlo_rr_wp_gg.o \

SectorsNNLOrv = \
      ./Sectors/mod_xsects_nnlo_rv_wp_ns.f90 \
      ./Sectors/mod_xsects_nnlo_rv_wp_qg.f90 \
      ./Sectors/mod_xsects_nnlo_rv_wp_gq.f90 \

SectorsNNLOrvObj = \
      $(ObjectDir)/mod_xsects_nnlo_rv_wp_ns.o \
      $(ObjectDir)/mod_xsects_nnlo_rv_wp_qg.o \
      $(ObjectDir)/mod_xsects_nnlo_rv_wp_gq.o \

SectorsNNLOvv = \
       ./Sectors/mod_xsects_nnlo_vv_wp_ns.f90 \

SectorsNNLOvvObj = \
       $(ObjectDir)/mod_xsects_nnlo_vv_wp_ns.o \

Checks = \
      ./Checks/mod_check_lim_nlo.f90 \
      ./Checks/mod_check_lim_rr_wp_ns.f90 \
      ./Checks/mod_check_lim_rr_wp_qg.f90 \
      ./Checks/mod_check_lim_rr_wp_gq.f90 \
      ./Checks/mod_check_lim_rr_wp_s.f90 \
      ./Checks/mod_check_lim_rr_wp_gg.f90 \
      ./Checks/mod_check_lim_rv.f90 \
      ./Checks/mod_check_lim_sub.f90 \
#      ./Checks/mod_check_lim_rr_h_gq.f90 \
#      ./Checks/mod_check_lim_rr_h_qq.f90 \


ChecksObj = \
      $(ObjectDir)/mod_check_lim_nlo.o \
      $(ObjectDir)/mod_check_lim_rr_wp_ns.o \
      $(ObjectDir)/mod_check_lim_rr_wp_qg.o \
      $(ObjectDir)/mod_check_lim_rr_wp_gq.o \
      $(ObjectDir)/mod_check_lim_rr_wp_s.o \
      $(ObjectDir)/mod_check_lim_rr_wp_gg.o \
      $(ObjectDir)/mod_check_lim_rv.o \
      $(ObjectDir)/mod_check_lim_sub.o \
#      $(ObjectDir)/mod_check_lim_rr_h_gq.o \
#      $(ObjectDir)/mod_check_lim_rr_h_qq.o \


Vegas = \
      vegas_lepage/vegas.f \
      mod_vegas_lepage.f90 \

VegasObj = \
      $(ObjectDir)/vegas.o \
      $(ObjectDir)/mod_vegas_lepage.o \

VegasInc = \
      vegas_lepage/maxwt.f \
      vegas_lepage/ran1.f \
      vegas_lepage/vegas_common.f \

Extra = \
      Chooser/mod_for_chooser.f90 \
      Chooser/mod_for_chooser_proc.f90 \
      Chooser/mod_tc_proc.f90 \
      Chooser/mod_lim_proc.f90 \
      Chooser/mod_run_proc.f90 \
      mod_chooser.f90 \

ExtraObj = \
      $(ObjectDir)/mod_for_chooser.o \
      $(ObjectDir)/mod_for_chooser_proc.o \
      $(ObjectDir)/mod_tc_proc.o \
      $(ObjectDir)/mod_lim_proc.o \
      $(ObjectDir)/mod_run_proc.o \
      $(ObjectDir)/mod_chooser.o \

Main = \
      main.f90 \

MainObj = \
      $(ObjectDir)/main.o \

AdditionalInc = \
      printout.f90 \
      gen_dump.script \


##################################################################################################

OURCODEBASIC = $(BasicDefObj) $(AuxFunctionsObj) $(GenProcObj) $(HistogramsObj) $(KinematicsGenObj) $(KinematicsLOObj) $(KinematicsNLOObj) $(KinematicsNNLO2qObj) $(KinematicsNNLOggObj) $(KinematicsNNLOqgObj) $(KinematicsNNLOEordAordObj) $(KinematicsNNLOEunordAordObj) $(KinematicsNNLOEunordAunordObj)

OURCODEAMPS = $(AmplitudesObj)

OURCODELIM = $(SingLimitsObj)

OURCODEINTSUB = $(TCSplitObj) $(IntSubObj)

OURCODERR = $(SectorsNNLOrrNSObj) $(SectorsNNLOrrQGObj) $(SectorsNNLOrrSObj) $(SectorsNNLOrrGGObj) 

OURCODESEC = $(SectorsLOObj) $(SectorsNLOObj) $(SectorsNNLOsubObj) $(OURCODERR) $(SectorsNNLOrvObj) $(SectorsNNLOvvObj)

OURCODECHECK = $(ChecksObj)

OURCODEMAIN = $(VegasObj) $(ExtraObj) $(MainObj)


OURCODE = $(OURCODEBASIC) $(OURCODEAMPS) $(OURCODELIM) $(OURCODEINTSUB) $(OURCODESEC) $(OURCODECHECK) $(OURCODEMAIN)

ALLCODE = $(OURCODE)

##################################################################################################

nnlo_wp: $(ALLCODE) 
	@echo " linking"
	@echo " executable file is " $(Exec)
	@echo " "
	$(fcomp4main)
	if [ -f *.o ]; then mv *.o $(ObjectDir)/.; fi
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi


$(BasicDefObj): $(BasicDef) $(config) $(BasicDecInc)
	@echo " compiling BasicDef "
	$(fcomp) -c $(BasicDef) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(AuxFunctionsObj): $(AuxFunctions) $(config) 
	@echo " compiling AuxFunctions "
	$(fcomp) -c $(AuxFunctions) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(GenProcObj): $(GenProc) $(config) 
	@echo " compiling GenProc "
	$(fcomp) -c $(GenProc) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(HistogramsObj): $(Histograms) $(config) $(HistogramsInc)
	@echo " compiling Histograms "
	sh gen_dump.script Histograms/cut_histo_w.f90 $(SEDcmd)
	$(fcomp) -c $(Histograms) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(KinematicsGenObj): $(KinematicsGen) $(config) 
	@echo " compiling KinematicsGen "
	$(fcomp) -c $(KinematicsGen) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(KinematicsLOObj): $(KinematicsLO) $(config) 
	@echo " compiling KinematicsLO "
	$(fcomp) -c $(KinematicsLO) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(KinematicsNLOObj): $(KinematicsNLO) $(config) 
	@echo " compiling KinematicsNLO "
	$(fcomp) -c $(KinematicsNLO) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(KinematicsNNLO2qObj): $(KinematicsNNLO2q) $(config) 
	@echo " compiling KinematicsNNLO2q "
	$(fcomp) -c $(KinematicsNNLO2q) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(KinematicsNNLOggObj): $(KinematicsNNLOgg) $(config) 
	@echo " compiling KinematicsNNLOgg "
	$(fcomp) -c $(KinematicsNNLOgg) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(KinematicsNNLOqgObj): $(KinematicsNNLOqg) $(config) 
	@echo " compiling KinematicsNNLOqg "
	$(fcomp) -c $(KinematicsNNLOqg) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(KinematicsNNLOEordAordObj): $(KinematicsNNLOEordAord) $(config) 
	@echo " compiling KinematicsNNLOEordAord "
	$(fcomp) -c $(KinematicsNNLOEordAord) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(KinematicsNNLOEunordAordObj): $(KinematicsNNLOEunordAord) $(config) 
	@echo " compiling KinematicsNNLOEunordAord "
	$(fcomp) -c $(KinematicsNNLOEunordAord) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(KinematicsNNLOEunordAunordObj): $(KinematicsNNLOEunordAunord) $(config) 
	@echo " compiling KinematicsNNLOEunordAunord "
	$(fcomp) -c $(KinematicsNNLOEunordAunord) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(AmplitudesObj): $(Amplitudes) $(config) $(AmplitudesInc)
	@echo " compiling Amplitudes "
	$(fcomp) -c $(Amplitudes) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(SingLimitsObj): $(SingLimits) $(config) 
	@echo " compiling SingLimits "
	$(fcomp) -c $(SingLimits) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(TCSplitObj): $(TCSplit) $(config) 
	@echo " compiling TCSplit "
	$(fcomp) -c $(TCSplit) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(IntSubObj): $(IntSub) $(config) 
	@echo " compiling IntSub "
	$(fcomp) -c $(IntSub) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(SectorsLOObj): $(SectorsLO) $(config) 
	@echo " compiling SectorsLO "
	$(fcomp) -c $(SectorsLO) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(SectorsNLOObj): $(SectorsNLO) $(config) 
	@echo " compiling SectorsNLO "
	$(fcomp) -c $(SectorsNLO) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(SectorsNNLOsubObj): $(SectorsNNLOsub) $(config) 
	@echo " compiling SectorsNNLOsub "
	$(fcomp) -c $(SectorsNNLOsub) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(SectorsNNLOrrNSObj): $(SectorsNNLOrrNS) $(config) 
	@echo " compiling SectorsNNLOrrNS "
	$(fcomp) -c $(SectorsNNLOrrNS) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(SectorsNNLOrrQGObj): $(SectorsNNLOrrQG) $(config) 
	@echo " compiling SectorsNNLOrrQG "
	$(fcomp) -c $(SectorsNNLOrrQG) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(SectorsNNLOrrSObj): $(SectorsNNLOrrS) $(config) 
	@echo " compiling SectorsNNLOrrS "
	$(fcomp) -c $(SectorsNNLOrrS) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(SectorsNNLOrrGGObj): $(SectorsNNLOrrGG) $(config) 
	@echo " compiling SectorsNNLOrrGG "
	$(fcomp) -c $(SectorsNNLOrrGG) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(SectorsNNLOrvObj): $(SectorsNNLOrv) $(config) 
	@echo " compiling SectorsNNLOrv "
	$(fcomp) -c $(SectorsNNLOrv) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(SectorsNNLOvvObj): $(SectorsNNLOvv) $(config) 
	@echo " compiling SectorsNNLOvv "
	$(fcomp) -c $(SectorsNNLOvv) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(ChecksObj): $(Checks) $(config) 
	@echo " compiling Checks "
	$(fcomp) -c $(Checks) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(VegasObj): $(Vegas) $(config) $(VegasInc)
	@echo " compiling Vegas "
	$(fcomp) -c $(Vegas) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(ExtraObj): $(Extra) $(config) 
	@echo " compiling Extra "
	$(fcomp) -c $(Extra) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi

$(MainObj): $(Main) $(config) $(AdditionalInc)
	@echo " compiling Main "
	$(fcomp) -c $(Main) 
	mv *.o $(ObjectDir)/.
	if [ "`echo *.mod`" != "*.mod" ]; then mv *.mod $(ModuleDir)/.; fi


######################################################################

clean:
	rm -f *~
	rm -f ./BasicDef/*~
	rm -f ./Kinematics/*~
	rm -f ./Histograms/*~
	rm -f ./AuxFunctions/*~
	rm -f ./Amplitudes/*~
	rm -f ./Amplitudes/CheckAmps/*~
	rm -f ./SingLimits/*~
	rm -f ./TCSplit/*~
	rm -f ./IntSub/*~
	rm -f ./Sectors/*~
	rm -f ./Sectors/include/*~
	rm -f ./GenProc/*~
	rm -f ./GenProc/include/*~
	rm -f ./vegas_lepage/*~
	rm -f ./Checks/*~
	rm -f ./Checks/CheckAmps/*~
	rm -f ./Chooser/*~

remove:
	rm -f $(Exec)
	rm -f $(ModuleDir)/*
	rm -f $(ObjectDir)/*

wiperes:
	rm -f histo/*.dat*
	rm -f histo/*.histo*
	rm -f grid/*.vg

madremove:
	rm -f P0_*.dat
	rm -f P1_*.dat
	rm -f P2_*.dat

madwipe:
	rm -f Checks/CheckAmps/open_*.f90
	rm -f Checks/CheckAmps/close_*.f90
	rm -f Checks/CheckAmps/write4mad_*.f90
	rm -f Amplitudes/CheckAmps/checkmad_*.f90
	rm -f MadScripts/tomad_*.sh 


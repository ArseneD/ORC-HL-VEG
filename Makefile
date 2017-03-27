#- $Id: AA_make.ldef 12 2010-11-05 15:42:13Z mmaipsl $
#---------------------------------------------------------------------
#-
#-
#- $Id: AA_make.gdef 2312 2014-08-04 14:02:17Z mafoipsl $
#-
#- Validate the correlation between the target and the environment
#-
UTIL_DIR = ../../util
WW_h_t = $(shell cat $(UTIL_DIR)/.host_target)
WW_h_w = $(shell $(UTIL_DIR)/w_i_h)
WW_t_e = $(shell $(UTIL_DIR)/w_i_e $(WW_h_t) $(WW_h_w))
ifeq "$(WW_t_e)" "NO"
 $(error )
endif
#-
######-Q- ada      F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
####-Q- curie  F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
######-Q- cur_mono  F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
#- Global definitions for gfortran, generic system
M_K = make
P_C = cpp
FCM_ARCH = gfortran 
P_O = -fpreprocessed -P -C -traditional $(P_P)
F_C = gfortran -c -cpp
F_D =
F_P = -fdefault-real-8
w_w = -O3 -funroll-all-loops $(F_D) $(F_P) -I$(MODDIR)
F_O = $(w_w) -J$(MODDIR)
F_L = gfortran
M_M = 0
L_X = 0
L_O =
A_C = ar -rs
A_G = ar -x
C_C = cc -c
C_O =
C_L = cc
#-
NCDF_INC = /opt/local/include
NCDF_LIB = -L/opt/local/lib -lnetcdff -lnetcdf
##-Q- gfortran  NCDF_LIB = -L/usr/local/lib -lnetcdf
#-
#####-Q- lxiv8    F_O = -DCPP_PARA -O3 -check bounds -traceback $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR) -fp-model precise
#####-Q- lxiv8    F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
######-Q- ciment F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
RM      = rm -f
STRIP   = strip
SIZE    = size

#- $Id: AA_make 1249 2013-04-04 12:00:37Z josefine.ghattas $
all_lib : libioipsl libparallel libparameters liborglob libstomate libsechiba

libioipsl:
	(cd ../IOIPSL/src; $(M_K) -f Makefile)

libparameters :
	(cd src_parameters ; $(M_K) -f Makefile)

libparallel :
	(cd ../../modeles/ORCHIDEE/src_parallel ; $(M_K) -f Makefile)

liborglob :
	(cd src_global ; $(M_K) -f Makefile)

libstomate :
	(cd src_stomate ; $(M_K) -f Makefile)

libsechiba :
	(cd src_sechiba ; $(M_K) -f Makefile)

driver : all_lib
	(cd src_driver ; $(M_K) -f Makefile)

config : 
	(cd src_parameters; $(M_K) -f Makefile config)
	(cd src_sechiba; $(M_K) -f Makefile config)
	(cd src_stomate; $(M_K) -f Makefile config)

clean : 
	(cd src_parameters; $(M_K) -f Makefile clean)
	(cd src_parallel; $(M_K) -f Makefile clean)
	(cd src_global; $(M_K) -f Makefile clean)
	(cd src_sechiba; $(M_K) -f Makefile clean)
	(cd src_stomate; $(M_K) -f Makefile clean)
	(cd src_driver; $(M_K) -f Makefile clean)
	(rm -f ../../lib/*)

doc:
	doxygen Doxyfile_ORCHIDEE
	export TEXINPUTS="${TEXINPUTS}:${PWD}/DOC//"; export BIBINPUTS="${BIBINPUTS}:${PWD}/DOC/BIB//"; cd docs/latex ; latex refman.tex 

bib:
	export TEXINPUTS="${TEXINPUTS}:${PWD}/DOC//"; export BIBINPUTS="${BIBINPUTS}:${PWD}/DOC/BIB//"; cd docs/latex ; bibtex refman 

index:
	export TEXINPUTS="${TEXINPUTS}:${PWD}/DOC//"; export BIBINPUTS="${BIBINPUTS}:${PWD}/DOC/BIB//"; cd docs/latex ; makeindex refman 

toc:
	export TEXINPUTS="${TEXINPUTS}:${PWD}/DOC//"; export BIBINPUTS="${BIBINPUTS}:${PWD}/DOC/BIB//"; cd docs/latex ; latex refman.tex 


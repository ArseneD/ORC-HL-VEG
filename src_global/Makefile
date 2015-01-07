#-
#- $Id: AA_make.ldef 1311 2013-06-04 08:13:54Z josefine.ghattas $
#-
#---------------------------------------------------------------------
#- Creation des elements relatifs a ORGLOB
#---------------------------------------------------------------------
SHELL = /bin/sh
#---------------------------------------------------------------------
LIBDIR = ../../../lib
BINDIR = ../../../bin
MODDIR = $(LIBDIR)
#---------------------------------------------------------------------
ORCHIDEE_LIB = $(LIBDIR)/liborchidee.a
MODEL_LIB = $(LIBDIR)/liborglob.a
SXMODEL_LIB = $(MODEL_LIB)
#-
#- $Id: AA_make.gdef 2312 2014-08-04 14:02:17Z mafoipsl $
#-
#- Validate the correlation between the target and the environment
#-
UTIL_DIR = ../../../util
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
RM      = rm -f
STRIP   = strip
SIZE    = size

#-
#- $Id: AA_make 1311 2013-06-04 08:13:54Z josefine.ghattas $
#-
PARAM_LIB = $(LIBDIR)/libparameters.a
SXPARAM_LIB = $(PARAM_LIB)
#-
MODS1 = grid.f90           \
        solar.f90          \
	interpol_help.f90  \
	gauss_jordan_method.f90  

OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)
#-
all:
	$(M_K) libparallel
	$(M_K) libparameters
	$(M_K) m_all
	@echo orglob is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

libparallel:
	(cd ../src_parallel; $(M_K) -f Makefile)

libparameters:
	(cd ../src_parameters; $(M_K) -f Makefile)

$(MODEL_LIB)(%.o) : %.f90
	$(F_C) $(F_O) -I$(NCDF_INC) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	$(A_C) $(ORCHIDEE_LIB) $*.o
	$(RM) $*.o

config : 
	$(BINDIR)/Fparser -name ORGLOB $(MODS1)
	echo 'Configuration of ORGLOB done'

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(grid.o): \
  $(PARAM_LIB)(constantes.o)

$(MODEL_LIB)(interpol_help.o): \
  $(PARAM_LIB)(constantes.o)

$(MODEL_LIB)(gauss_jordan_method.o): \
  $(PARAM_LIB)(constantes.o)

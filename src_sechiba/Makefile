#-
#- $Id: AA_make.ldef 1311 2013-06-04 08:13:54Z josefine.ghattas $
#-
#---------------------------------------------------------------------
#- Creation des elements relatifs a SECHIBA
#---------------------------------------------------------------------
SHELL = /bin/sh
#---------------------------------------------------------------------
LIBDIR = ../../../lib
BINDIR = ../../../bin
MODDIR = $(LIBDIR)
#---------------------------------------------------------------------
ORCHIDEE_LIB = $(LIBDIR)/liborchidee.a
MODEL_LIB = $(LIBDIR)/libsechiba.a
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
######-Q- ciment F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
RM      = rm -f
STRIP   = strip
SIZE    = size

#-
#- $Id: AA_make 2222 2014-07-08 05:51:19Z josefine.ghattas $
#-
PARALLEL_LIB = $(LIBDIR)/libparallel.a
SXPARALLEL_LIB = $(PARALLEL_LIB)
#-
PARAM_LIB = $(LIBDIR)/libparameters.a
SXPARAM_LIB = $(PARAM_LIB)
#-
ORGLOB_LIB = $(LIBDIR)/liborglob.a
SXORGLOB_LIB = $(ORGLOB_LIB)
#-
STOMATE_LIB = $(LIBDIR)/libstomate.a
SXSTOMATE_LIB = $(STOMATE_LIB)
#-
MODS1 =              \
  qsat_moisture.f90 \
  sechiba_io_p.f90   \
  sechiba_io.f90     \
  watchout.f90       \
  slowproc.f90       \
  diffuco.f90        \
  condveg.f90        \
  enerbil.f90        \
  explicitsnow.f90       \
  hydrol.f90         \
  hydrolc.f90        \
  thermosoil.f90     \
  gammad_inc.f90     \
  hydro_subgrid.f90  \
  init_top.f90       \
  routing.f90        \
  sechiba.f90        \
  intersurf.f90 
OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)
#-
all:
	$(M_K) libparallel
	$(M_K) libparameters
	$(M_K) liborglob
	$(M_K) libstomate
	$(M_K) m_all
	@echo sechiba is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

libparallel:
	(cd ../src_parallel; $(M_K) -f Makefile)

libparameters:
	(cd ../src_parameters; $(M_K) -f Makefile)

liborglob:
	(cd ../src_global; $(M_K) -f Makefile)

libstomate:
	(cd ../src_stomate; $(M_K) -f Makefile)

$(MODEL_LIB)(%.o) : %.f90
	$(F_C) $(F_O) -I$(NCDF_INC) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	$(A_C) $(ORCHIDEE_LIB) $*.o
	$(RM) $*.o

config : 
	$(BINDIR)/Fparser -name SECHIBA $(MODS1)
	echo 'Configuration of SECHIBA done'

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(sechiba.o):          \
  $(PARAM_LIB)(pft_parameters.o)  \
  $(MODEL_LIB)(sechiba_io_p.o)    \
  $(MODEL_LIB)(sechiba_io.o)      \
  $(MODEL_LIB)(diffuco.o)         \
  $(MODEL_LIB)(condveg.o)         \
  $(MODEL_LIB)(enerbil.o)         \
  $(MODEL_LIB)(hydrolc.o)         \
  $(MODEL_LIB)(hydrol.o)          \
  $(MODEL_LIB)(thermosoil.o)      \
  $(MODEL_LIB)(slowproc.o)        \
  $(MODEL_LIB)(routing.o)
 
$(MODEL_LIB)(sechiba_io_p.o): \
  $(PARAM_LIB)(constantes.o)  \
  $(PARALLEL_LIB)(mod_orchidee_para.o)

$(MODEL_LIB)(sechiba_io.o): \
  $(MODEL_LIB)(sechiba_io_p.o) \
  $(PARAM_LIB)(constantes.o)

$(MODEL_LIB)(watchout.o): \
  $(ORGLOB_LIB)(grid.o)

$(MODEL_LIB)(hydrol.o):           \
  $(PARAM_LIB)(constantes.o)  \
  $(MODEL_LIB)(explicitsnow.o)  \
  $(MODEL_LIB)(sechiba_io.o)      \
  $(MODEL_LIB)(init_top.o)        \
  $(MODEL_LIB)(hydro_subgrid.o)   \
  $(ORGLOB_LIB)(grid.o)

$(MODEL_LIB)(explicitsnow.o):         \
  $(PARAM_LIB)(pft_parameters.o)  \
  $(PARAM_LIB)(constantes.o)      \
  $(PARAM_LIB)(constantes_soil.o) \
  $(MODEL_LIB)(qsat_moisture.o)   \
  $(MODEL_LIB)(sechiba_io.o)

$(MODEL_LIB)(hydro_subgrid.o):    \
  $(PARAM_LIB)(constantes.o)      \
  $(PARAM_LIB)(constantes_soil.o) \

$(MODEL_LIB)(init_top.o):         \
  $(PARAM_LIB)(constantes.o)      \
  $(PARAM_LIB)(pft_parameters.o)  \
  $(MODEL_LIB)(gammad_inc.o)

$(MODEL_LIB)(hydrolc.o):          \
  $(PARAM_LIB)(constantes.o)  \
  $(MODEL_LIB)(explicitsnow.o)  \
  $(MODEL_LIB)(sechiba_io.o)

$(MODEL_LIB)(slowproc.o):         \
  $(PARAM_LIB)(pft_parameters.o)  \
  $(STOMATE_LIB)(stomate.o)       \
  $(MODEL_LIB)(sechiba_io.o)      \
  $(ORGLOB_LIB)(interpol_help.o)

$(MODEL_LIB)(diffuco.o):          \
  $(PARAM_LIB)(pft_parameters.o)  \
  $(MODEL_LIB)(qsat_moisture.o)  \
  $(MODEL_LIB)(sechiba_io.o)

$(MODEL_LIB)(enerbil.o):          \
  $(PARAM_LIB)(constantes.o)      \
  $(MODEL_LIB)(qsat_moisture.o)   \
  $(MODEL_LIB)(explicitsnow.o)        \
  $(MODEL_LIB)(sechiba_io.o)

$(MODEL_LIB)(condveg.o):          \
  $(PARAM_LIB)(pft_parameters.o)  \
  $(PARAM_LIB)(constantes.o)  \
  $(PARAM_LIB)(constantes_soil.o)  \
  $(PARAM_LIB)(pft_parameters.o)  \
  $(MODEL_LIB)(sechiba_io.o) \
  $(ORGLOB_LIB)(interpol_help.o)

$(MODEL_LIB)(thermosoil.o): \
  $(PARAM_LIB)(constantes.o)  \
  $(PARAM_LIB)(constantes_soil.o)  \
  $(PARAM_LIB)(pft_parameters.o)  \
  $(MODEL_LIB)(sechiba_io.o)

$(MODEL_LIB)(routing.o):          \
  $(PARAM_LIB)(constantes.o)  \
  $(MODEL_LIB)(sechiba_io.o)

$(MODEL_LIB)(intersurf.o):   \
  $(PARAM_LIB)(constantes.o)  \
  $(ORGLOB_LIB)(grid.o)      \
  $(MODEL_LIB)(watchout.o)   \
  $(MODEL_LIB)(thermosoil.o) \
  $(MODEL_LIB)(sechiba.o)

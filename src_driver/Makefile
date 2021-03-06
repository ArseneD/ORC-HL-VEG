#-
#- $Id: AA_make.ldef 1042 2012-10-26 14:32:10Z josefine.ghattas $
#-
SHELL = /bin/sh
#---------------------------------------------------------------------
LIBDIR = ../../../lib
MODDIR = $(LIBDIR)
BINDIR = ../../../bin
SCPDIR = ../../scripts
ORDIR  = ..
#---------------------------------------------------------------------
MODEL_LIB = $(LIBDIR)/liborchidee.a
SXMODEL_LIB = $(MODEL_LIB)
#-
IOIPSL_LIB = $(LIBDIR)/libioipsl.a
SXIOIPSL_LIB = $(IOIPSL_LIB)
#-
PARAM_LIB = $(LIBDIR)/libparameters.a
SXPARAM_LIB = $(PARAM_LIB)
#-
PARALLEL_LIB = $(LIBDIR)/libparallel.a
SXPARALLEL_LIB = $(PARALLEL_LIB)
#-
ORGLOB_LIB = $(LIBDIR)/liborglob.a
SXORGLOB_LIB = $(ORGLOB_LIB)
#-
SECHIBA_LIB = $(LIBDIR)/libsechiba.a
SXSECHIBA_LIB = $(SECHIBA_LIB)
#-
STOMATE_LIB = $(LIBDIR)/libstomate.a
SXSTOMATE_LIB = $(STOMATE_LIB)
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
#- $Id: AA_make 1042 2012-10-26 14:32:10Z josefine.ghattas $
#-
MODS1 = dim2_driver.f90 \
        weather.f90 \
        readdim2.f90 \
        teststomate.f90 \
        forcesoil.f90 \
        getprec.f90

all : orchidee_ol teststomate forcesoil

orchidee_ol:
	$(M_K) m_orchidee_ol
	@echo orchidee_ol is OK

m_orchidee_ol: $(BINDIR)/orchidee_ol
$(BINDIR)/orchidee_ol: \
  $(SXIOIPSL_LIB) $(SXPARAM_LIB) $(SXPARALLEL_LIB) $(SXORGLOB_LIB) $(SXSTOMATE_LIB) $(SXSECHIBA_LIB) \
  $(MODEL_LIB)(weather.o) $(MODEL_LIB)(readdim2.o) \
  $(MODEL_LIB)(dim2_driver.o)
	$(A_G) $(MODEL_LIB) dim2_driver.o
	$(F_L) -o $(BINDIR)/orchidee_ol dim2_driver.o $(L_O) \
          $(SXMODEL_LIB) $(SXSECHIBA_LIB) $(SXSTOMATE_LIB) \
          $(SXPARALLEL_LIB) $(SXORGLOB_LIB) $(SXPARAM_LIB) $(SXIOIPSL_LIB) $(NCDF_LIB) 
	$(RM) dim2_driver.o

teststomate:
	$(M_K) m_teststomate
	@echo teststomate is OK

m_teststomate: $(BINDIR)/teststomate
$(BINDIR)/teststomate: \
  $(SXIOIPSL_LIB) $(SXPARALLEL_LIB) $(SXPARAM_LIB) $(SXORGLOB_LIB) $(SXSTOMATE_LIB) $(SXSECHIBA_LIB) \
  $(MODEL_LIB)(teststomate.o)
	$(A_G) $(MODEL_LIB) teststomate.o
	$(F_L) -o $(BINDIR)/teststomate teststomate.o $(L_O) \
          $(SXMODEL_LIB) $(SXSECHIBA_LIB) $(SXSTOMATE_LIB) \
          $(SXPARALLEL_LIB) $(SXORGLOB_LIB) $(SXPARAM_LIB) $(SXIOIPSL_LIB) $(NCDF_LIB)
	$(RM) teststomate.o

forcesoil:
	$(M_K) m_forcesoil
	@echo m_forcesoil is OK

m_forcesoil: $(BINDIR)/forcesoil
$(BINDIR)/forcesoil: \
  $(SXIOIPSL_LIB) $(SXPARALLEL_LIB) $(SXPARAM_LIB) $(SXORGLOB_LIB) $(SXSTOMATE_LIB) $(SXSECHIBA_LIB) \
  $(MODEL_LIB)(forcesoil.o)
	$(A_G) $(MODEL_LIB) forcesoil.o
	$(F_L) -o $(BINDIR)/forcesoil forcesoil.o $(L_O) \
          $(SXMODEL_LIB) $(SXSECHIBA_LIB) $(SXSTOMATE_LIB) \
          $(SXPARALLEL_LIB) $(SXORGLOB_LIB) $(SXPARAM_LIB) $(SXIOIPSL_LIB) $(NCDF_LIB)
	$(RM) forcesoil.o

libioipsl:
	(cd ../../IOIPSL/src; $(M_K) -f Makefile)

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) -I$(NCDF_INC) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	$(RM) $*.o


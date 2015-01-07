#-
#- $Id: AA_make.ldef 1311 2013-06-04 08:13:54Z josefine.ghattas $
#-
#---------------------------------------------------------------------
#- Creation des elements relatifs a PARALLEL
#---------------------------------------------------------------------
SHELL = /bin/sh
#---------------------------------------------------------------------
LIBDIR = ../../../lib
BINDIR = ../../../bin
MODDIR = $(LIBDIR)
ORDIR  = ..
#---------------------------------------------------------------------
ORCHIDEE_LIB = $(LIBDIR)/liborchidee.a
MODEL_LIB = $(LIBDIR)/libparallel.a
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
#- $Id: AA_make 1931 2014-02-26 13:12:56Z josefine.ghattas $
#-
MODS1 = mod_orchidee_para_var.F90 \
        xios_orchidee.f90 \
        mod_orchidee_mpi_data.F90 \
	mod_orchidee_omp_data.F90 \
	timer.f90 \
        mod_orchidee_mpi_transfert.F90 \
        mod_orchidee_omp_transfert.F90 \
        mod_orchidee_transfert_para.F90 \
        mod_orchidee_para.F90 \
	tools_para.f90 \
	ioipsl_para.f90 \
	orch_write_field.f90 \
	write_field.f90 \
	orch_write_field_p.f90 \
	write_field_p.f90


PARAM_LIB=$(LIBDIR)/libparameters.a
FILTERf90=$(filter %.f90,$(MODS1))
FILTERF90=$(filter %.F90,$(MODS1))
OBJSMODS1=$(FILTERf90:.f90=.o) $(FILTERF90:.F90=.o)
#-
.PRECIOUS : $(MODEL_LIB)
#-
all:
	$(M_K) m_all
	@echo parallel is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

$(MODEL_LIB)(%.o) : %.f90
	$(F_C) $(F_O) -I$(NCDF_INC) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	$(A_C) $(ORCHIDEE_LIB) $*.o
	$(RM) $*.o

$(MODEL_LIB)(%.o) : %.F90
	$(F_C) $(F_O) -I$(NCDF_INC) $*.F90
	$(A_C) $(MODEL_LIB) $*.o
	$(A_C) $(ORCHIDEE_LIB) $*.o
	$(RM) $*.o


$(PARAM_LIB)(pft_parameters_var.o) : ../src_parameters/pft_parameters_var.f90
	$(F_C) $(F_O) -I$(NCDF_INC) ../src_parameters/pft_parameters_var.f90
	$(A_C) $(PARAM_LIB) pft_parameters_var.o
	$(A_C) $(ORCHIDEE_LIB) pft_parameters_var.o
	$(RM) ../src_parallel/pft_parameters_var.o

$(PARAM_LIB)(constantes_var.o) : ../src_parameters/constantes_var.f90
	$(F_C) $(F_O) -I$(NCDF_INC) ../src_parameters/constantes_var.f90
	$(A_C) $(PARAM_LIB) constantes_var.o
	$(A_C) $(ORCHIDEE_LIB) constantes_var.o
	$(RM) ../src_parallel/constantes_var.o

$(PARAM_LIB)(constantes_soil_var.o) : ../src_parameters/constantes_soil_var.f90
	$(F_C) $(F_O) -I$(NCDF_INC) ../src_parameters/constantes_soil_var.f90
	$(A_C) $(PARAM_LIB) constantes_soil_var.o
	$(A_C) $(ORCHIDEE_LIB) constantes_soil_var.o
	$(RM) ../src_parallel/constantes_soil_var.o

config : 
	$(BINDIR)/Fparser -name PARALLEL $(MODS1)
	echo 'Configuration of PARALLEL done'

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(mpi_dummy.o):

$(MODEL_LIB)(mod_orchidee_para_var.o):

$(MODEL_LIB)(xios_orchidee.o): \
       $(MODEL_LIB)(mod_orchidee_para_var.o) \
       $(MODEL_LIB)(mod_orchidee_transfert_para.o) \
       $(MODEL_LIB)(ioipsl_para.o) \
       $(PARAM_LIB)(pft_parameters_var.o) \
       $(PARAM_LIB)(constantes_var.o) \
       $(PARAM_LIB)(constantes_soil_var.o)

$(MODEL_LIB)(mod_orchidee_mpi_data.o): \
    $(MODEL_LIB)(mod_orchidee_para_var.o) \
    $(MODEL_LIB)(xios_orchidee.o) 

$(MODEL_LIB)(mod_orchidee_omp_data.o): \
    $(MODEL_LIB)(mod_orchidee_para_var.o)

$(MODEL_LIB)(timer.o): \
    $(MODEL_LIB)(mod_orchidee_para_var.o)

$(MODEL_LIB)(orch_write_field.o): \
    $(MODEL_LIB)(mod_orchidee_para.o)

$(MODEL_LIB)(write_field.o): \
  $(MODEL_LIB)(orch_write_field.o)

$(MODEL_LIB)(mod_orchidee_mpi_transfert.o): \
  $(MODEL_LIB)(mod_orchidee_para_var.o) \
  $(MODEL_LIB)(timer.o)

$(MODEL_LIB)(mod_orchidee_omp_transfert.o): \
  $(MODEL_LIB)(mod_orchidee_omp_data.o)

$(MODEL_LIB)(mod_orchidee_transfert_para.o): \
  $(MODEL_LIB)(mod_orchidee_para_var.o) \
  $(MODEL_LIB)(mod_orchidee_mpi_transfert.o) \
  $(MODEL_LIB)(mod_orchidee_omp_transfert.o)

$(MODEL_LIB)(mod_orchidee_para.o): \
  $(MODEL_LIB)(mod_orchidee_mpi_data.o) \
  $(MODEL_LIB)(mod_orchidee_omp_data.o) \
  $(MODEL_LIB)(mod_orchidee_transfert_para.o)

$(MODEL_LIB)(ioipsl_para.o): \
  $(MODEL_LIB)(mod_orchidee_transfert_para.o) \
  $(MODEL_LIB)(mod_orchidee_para_var.o)

$(MODEL_LIB)(tools_para.o): \
  $(MODEL_LIB)(mod_orchidee_para_var.o)

$(MODEL_LIB)(orch_write_field_p.o): \
  $(MODEL_LIB)(mod_orchidee_para.o) \
  $(MODEL_LIB)(write_field.o)

$(MODEL_LIB)(write_field_p.o): \
  $(MODEL_LIB)(orch_write_field_p.o)


.SUFFIXES: .F .o

all: dummy aerosols_monan

dummy:
	echo "****** compiling aerosols_physics_monan ******"

OBJS = \
	module_aer_dust.o 

aerosols_monan: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)

# DEPENDENCIES:
# module_aer_dust.o:      

clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../physics_mmm -I../../../framework -I../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../physics_mmm -I../../../framework -I../../../external/esmf_time_f90
endif

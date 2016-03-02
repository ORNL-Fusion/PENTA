#-----------------------------------------------------------------------
# This is the makefile for PENTA3
#	7/14/2010-7/20/2010 JL
#-----------------------------------------------------------------------

SHELL = /bin/sh

#FCOMP = pgf90
#FLDR = pgf90
#FFLAGS = -fast
##FFLAGS = -fast -Msmartalloc
#FFLAGS_DBG = -g -Mbounds 
#LDRFLAGS = 

FCOMP = gfortran
FLDR = gfortran
FFLAGS = -O3
FFLAGS_DBG = -g -Wall -frange-check -fbounds-check 
LDRFLAGS =

SOURCEPATH = ../SRC/
EXECUTABLE = PENTA3
DEBUG_EX = PENTA3_dbg

# External library paths and libraries

LIBDIR = -L.
LIBS = 

#-----------------------------------------------------------------------

# The following macro lists all of the .90 object files that are
# needed to build the executable.  The "\" signifies that the line 
# is continued

OBJS90 = penta.o read_input_file_mod.o penta_functions_mod.o penta_subroutines.o\
 penta_math_routines_mod.o bspline.o pprof_pass.o vmec_var_pass.o coeff_var_pass.o

# This macro lists the modules to be compiled which have identically named .f90 files

MODS90 = penta_kind_mod.mod io_unit_spec.mod read_input_file_mod.mod phys_const.mod penta_functions_mod.mod\
 penta_math_routines_mod.mod vmec_var_pass.mod pprof_pass.mod coeff_var_pass.mod bspline.mod

#-----------------------------------------------------------------------

# The first dependency list is the default, so it should be
# used for the final executable.  Other dependencies may get
# invoked in the process of creating the final executable.
# Note that $(macro) evaluates a macro or an environment variable,
# and the list of commands or rules follow the dependency.

$(EXECUTABLE) : $(OBJS90)
	@echo "Creating "$@" in directory "$(PWD)"."
	$(FLDR) -o $@ $(LDRFLAGS) $(OBJS90) $(LIBDIR) $(LIBS)

# The following dependency is similar, but it conditionally
# replaces the value of FFLAGS with FFLAGS_DBG when
# processing the rules for the targets

debug: FFLAGS = $(FFLAGS_DBG)
debug: $(OBJS90)
	@echo "Creating "$(DEBUG_EX)" in directory "$(PWD)"."
	$(FLDR) -o $(DEBUG_EX) $(LDRFLAGS) $(OBJS90) $(LIBDIR) $(LIBS)

#-----------------------------------------------------------------------

# The following are static pattern rules that provide
# instructions for creating object files from the fortran
# source files.  The older suffix rule would have been ".f.o:"

# The $< at the end is a macro for the current source file.

$(MODS90) : %.mod :  $(SOURCEPATH)%.f90
	$(FCOMP) $(FFLAGS) -c $<

$(OBJS90) : %.o:  $(SOURCEPATH)%.f90
	$(FCOMP) $(FFLAGS) -c $<

#-----------------------------------------------------------------------

# A list of module dependencies ensures that module information
# specific to a particular source is available.


coeff_var_pass.mod : penta_kind_mod.mod
io_unit_spec.mod : penta_kind_mod.mod
penta_math_routines_mod.mod : penta_kind_mod.mod phys_const.mod
phys_const.mod : penta_kind_mod.mod
pprof_pass.mod : penta_kind_mod.mod
vmec_var_pass.mod : penta_kind_mod.mod
penta_functions_mod.mod : penta_kind_mod.mod phys_const.mod\
 penta_math_routines_mod.mod bspline.mod
read_input_file_mod.mod : penta_kind_mod.mod vmec_var_pass.mod\
 pprof_pass.mod bspline.mod coeff_var_pass.mod io_unit_spec.mod\
 bspline.mod
$(SOURCEPATH)penta_subroutines.f90: penta_kind_mod.mod phys_const.mod\
 penta_functions_mod.mod pprof_pass.mod vmec_var_pass.mod bspline.mod
$(SOURCEPATH)penta.f90: penta_kind_mod.mod io_unit_spec.mod\
 read_input_file_mod.mod vmec_var_pass.mod phys_const.mod\
 penta_functions_mod.mod pprof_pass.mod coeff_var_pass.mod\
 penta_math_routines_mod.mod

#-----------------------------------------------------------------------

# Specifying 'clean' dependencies is also handy.  With no
# listed dependencies, clean is considered a 'phony' target.
# The dash in front of a command means to continue even if that
# command returns a nonzero error code.

clean:
	-rm *.o *.mod
	@echo "Cleaned object and mod files from "$(PWD)"."




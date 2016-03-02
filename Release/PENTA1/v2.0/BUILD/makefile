#-----------------------------------------------------------------------
# This is the makefile for PENTA1
#	7/24/2009 JL
#-----------------------------------------------------------------------

SHELL = /bin/sh


FCOMP = gfortran
FLDR = gfortran
FFLAGS = -O3
FFLAGS_DBG = -g -Wall -frange-check -fbounds-check 
LDRFLAGS = 

SOURCEPATH = ../SRC/
EXECUTABLE = PENTA1
DEBUG_EX = PENTA1_dbg

# External library paths and libraries

LIBDIR = -L.
LIBS = 


#-----------------------------------------------------------------------

OBJS = penta_imp.o penta_modules.o penta_subroutines.o \
 penta_functions_mod.o quanc8.o read_input_file_mod.o 

OBJS90 = bspline.o

# This is a module-list macro.

MODS = penta_modules.mod read_input_file_mod.mod penta_subroutines.mod \
penta_functions_mod.mod penta_kind_mod.mod
MODS90 = bspline.mod numeric.mod

#-----------------------------------------------------------------------

# The first dependency list is the default, so it should be
# used for the final executable.  Other dependencies may get
# invoked in the process of creating the final executable.
# Note that $(macro) evaluates a macro or an environment variable,
# and the list of commands or rules follow the dependency.

$(EXECUTABLE) : $(OBJS) $(OBJS90)
	@echo "Creating "$@" in directory "$(PWD)"."
	$(FLDR) -o $@ $(LDRFLAGS) $(OBJS) $(OBJS90) $(LIBDIR) $(LIBS)

# The following dependency is similar, but it conditionally
# replaces the value of FFLAGS with FFLAGS_DBG when
# processing the rules for the target aph_dbg.

debug : FFLAGS = $(FFLAGS_DBG)
debug :  $(OBJS) $(OBJS90)
	@echo "Creating "$(DEBUG_EX)" in directory "$(PWD)"."
	$(FLDR) -o $(DEBUG_EX) $(LDRFLAGS) $(OBJS) $(OBJS90) $(LIBDIR) $(LIBS)

#-----------------------------------------------------------------------

# The following is a static pattern rule that provides
# instructions for creating object files from the fortran
# source files.  The older suffix rule would have been ".f.o:"

# The symbol at the end is a macro for the current source file.

$(OBJS) : %.o : $(SOURCEPATH)%.f
	$(FCOMP) $(FFLAGS) -c $<
$(OBJS90) : %.o : $(SOURCEPATH)%.f90
	$(FCOMP) $(FFLAGS) -c $<

$(MODS) : %.mod : $(SOURCEPATH)%.f
	$(FCOMP) $(FFLAGS) -c $<

$(MODS90) : %.mod : $(SOURCEPATH)%.f90
	$(FCOMP) $(FFLAGS) -c $<

# A list of module dependencies ensures that module information
# specific to a particular source is available.

$(SOURCEPATH)penta_imp.f: read_input_file_mod.mod\
 penta_modules.mod penta_subroutines.mod bspline.mod penta_functions_mod.mod
$(SOURCEPATH)penta_subroutines.f: penta_functions_mod.mod
read_input_file_mod.mod : bspline.mod penta_modules.mod
penta_subroutines.mod : bspline.mod penta_modules.mod penta_functions_mod.mod
penta_functions_mod.mod : bspline.mod

#-----------------------------------------------------------------------

# Specifying 'clean' dependencies is also handy.  With no
# listed dependencies, clean is considered a 'phony' target.
# The dash in front of a command means to continue even if that
# command returns a nonzero error code.

clean:
	-rm *.o *.mod
	@echo "Cleaned object and mod files from "$(PWD)"."




##################################################################
#                        Makefile for f90                         #
###################################################################

#BGP_SYS = /bgsys/drivers/ppcfloor/comm

########################     compiler     #########################

#F90 = ifort
F90 = mpif90
F77 = mpif90
#F90 = $(BGP_SYS)/bin/mpixlf90
#F77 = $(BGP_SYS)/bin/mpixlf77

########################  compiler flags  #########################

#F90FLAGS= -c
#F90FLAGS= -c 
F90FLAGS= -c #-check bounds -g -traceback#-c #-fbacktrace -fbounds-check 
#F90FLAGS = -c -warn -CB -debug extended
F77FLAGS= -c #-g -traceback
LFLAGS =
#PREP = scalasca -inst
#PREP = skin

########################  objects alpha   #########################

INIT    = .
SRCDIR  = $(INIT)
OBJ     = $(INIT)
OBJDIR  = $(OBJ)
CALCDIR = $(INIT)

OBJECTS = $(OBJ)/declaration.o\
          $(OBJ)/start.o\
          $(OBJ)/y_grid_substrate.o\
          $(OBJ)/geometry_sym.o\
          $(OBJ)/stats.o\
          $(OBJ)/spectra.o\
          $(OBJ)/tridLU_3D.o\
          $(OBJ)/rft_buff.o\
          $(OBJ)/cft_buff.o\
          $(OBJ)/FOU3D.o\
          $(OBJ)/littleharsh.o

########################      build       #########################

littleharsh :printmsgA $(OBJECTS)
	@echo
	@echo Linking... 
	$(PREP) $(F90) -o $@ $(OBJECTS) $(LFLAGS)
	@echo 
	@echo littleharsh built, congratulations.
	@echo 

########################     compile      #########################

$(OBJDIR)/declaration.o : $(SRCDIR)/declaration.f90 $(SRCDIR)/makefile
	@echo compiling declaration.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/declaration.f90 

$(OBJDIR)/start.o : $(SRCDIR)/start.f90 $(SRCDIR)/makefile
	@echo compiling start.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/start.f90

$(OBJDIR)/y_grid_substrate.o : $(SRCDIR)/y_grid_substrate.f90 $(SRCDIR)/makefile
	@echo compiling y_grid_substrate.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/y_grid_substrate.f90

$(OBJDIR)/geometry_sym.o : $(SRCDIR)/geometry_sym.f90 $(SRCDIR)/makefile
	@echo compiling geometry_sym.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/geometry_sym.f90

$(OBJDIR)/stats.o : $(SRCDIR)/stats.f90 $(SRCDIR)/makefile
	@echo compiling stats.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/stats.f90

$(OBJDIR)/spectra.o : $(SRCDIR)/spectra.f90 $(SRCDIR)/makefile
	@echo compiling spectra.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/spectra.f90

$(OBJDIR)/tridLU_3D.o : $(SRCDIR)/tridLU_3D.f90 $(SRCDIR)/makefile
	@echo compiling tridLU_3D.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/tridLU_3D.f90

$(OBJDIR)/FOU3D.o : $(SRCDIR)/FOU3D.f90 $(SRCDIR)/makefile
	@echo compiling FOU3D.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/FOU3D.f90

$(OBJDIR)/rft_buff.o : $(SRCDIR)/rft_buff.f $(SRCDIR)/makefile
	@echo compiling rft_buff.f
	@cd $(OBJDIR); $(PREP) $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/rft_buff.f

$(OBJDIR)/cft_buff.o : $(SRCDIR)/cft_buff.f $(SRCDIR)/makefile
	@echo compiling cft_buff.f
	@cd $(OBJDIR); $(PREP) $(F77) $(F77FLAGS) -I$(SRCDIR) $(SRCDIR)/cft_buff.f

$(OBJDIR)/littleharsh.o : $(SRCDIR)/littleharsh.f90  $(SRCDIR)/makefile
	@echo compiling littleharsh.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/littleharsh.f90

########################      message     #########################

printmsgA :
	@echo
	@echo Building ...
	@echo F77 Compiler flags : $(F77FLAGS)
	@echo F90 Compiler flags : $(F90FLAGS)
	@echo Linker   flags : $(LFLAGS)
	@echo Prepend  flags : $(PREP)

########################      clean       #########################

clean:
	@find . \( -name '*.o' \) -exec rm {} \;

########################   end of file    #########################

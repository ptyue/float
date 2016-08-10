#
# partflow executable makefile
#
# Need to set environment variables PFL_DIR (partflow directory)
# and files.
#
#
# set machine dependent FLAGS
#
include ../flags/copts.ifc
#
PROG       = partflow
SRC1_F     = partflow.f
SRC2_F     = rmesh.f tecpost.f rmshgrmp.f flow129bdry.f rmshgmsh.f
SRC_H      = ../include/parameter2d.h
FLDIR      = -L$(PFL_DIR)/lib -L/usr/lib64
LIBS       = -lflows -lblas -lgrummp -lgmsh
ALL_LIBS   = $(FLDIR) $(LIBS)  
OBJ1_F     = $(SRC1_F:.f=.o)
OBJ2_F     = $(SRC2_F:.f=.o)
#
#	set GRUMMP varialbes
#
GRUMLIBDIR = $(HOME)/mesh/GRUMMP-0.4.0alpha3/lib/libO/linux-gnu
GRUMINCDIR = $(HOME)/mesh/GRUMMP-0.4.0alpha3/include
GRUMINC = -I$(GRUMINCDIR)
GRUMLIBS = -L$(GRUMLIBDIR) -lGR_3D -lGR_surf -lGR_2D -lGR_base -lOptMS -lSUMAAlog_lite -lm
#
#	set Gmsh variables
#
GMSHDIR = /home/ptyue/mesh/gmsh2.8.4/gmsh-2.8.4-source
GMSHLIBS = -L$(GMSHDIR)/lib -lGmsh
GMSHINC     = -I$(GMSHDIR)/Common -I$(GMSHDIR)/Fltk -I$(GMSHDIR)/Geo -I$(GMSHDIR)/Graphics -I$(GMSHDIR)/Mesh -I$(GMSHDIR)/Numeric -I$(GMSHDIR)/Parser -I$(GMSHDIR)/Plugin -I$(GMSHDIR)/Post -I$(GMSHDIR)/Qt -I$(GMSHDIR)/Solver -I$(GMSHDIR)/lib/Common


all : $(PROG)
#
.KEEP_STATE:
#
# define dependencies of include files
#
$(OBJ1_F) : $(SRC_H)
#
# compile programs
#
$(PROG): $(OBJ1_F) $(OBJ2_F) $(FULL_LIBS)
	$(FC) $(LNKOPT) -o $(PROG) $(OBJ1_F) $(OBJ2_F) $(ALL_LIBS) $(GRUMINC) $(GRUMLIBS) $(GMSHINC) $(GMSHLIBS)
	mv $(PROG) $(PFL_DIR)/bin/$(PROG)
#
tags: $(SRC1_F) $(SRC2_F)
	ctags -x $< > tags
clean :
	/bin/rm -f *.o *.out
#cleanall :
#	/bin/rm -f *.o *.a *.out $(PFL_DIR)/bin/*.out $(PFL_DIR)/lib/*.a
# DO NOT DELETE

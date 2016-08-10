include ../flags/copts.ifc
LIB = libgrummp.a
#LIBDIR = $(HOME)/GRUMMP/GRUMMP-0.2.2alpha3/lib/libO/linux
#INCDIR = $(HOME)/GRUMMP/GRUMMP-0.2.2alpha3/include
LIBDIR = $(HOME)/mesh/GRUMMP-0.4.0alpha3/lib/libO/linux
INCDIR = $(HOME)/mesh/GRUMMP-0.4.0alpha3/include
INC = -I$(INCDIR)
LIB_3D = -L$(LIBDIR) -lGR_3D -lGR_surf -lGR_2D -lGR_base -lOptMS -lSUMAAlog_lite -lm
LIBS = libGR_2D.a  libGR_base.a  libOptMS.a  libSUMAAlog_lite.a

#intrgrummp:GRUMMP_IF_2D.o intrgrummp.o
#	$(FC) -o main Fortran.o GRUMMP_IF_2D.o intrgrummp.o $(INC) $(LIB_3D)

all: $(LIB)
GRUMMP_IF_2D.o:GRUMMP_IF_2D.cxx
	$(CC) -c $(CFLAGS)  GRUMMP_IF_2D.cxx $(INC) 
intrgrummp.o:intrgrummp.f 
	$(FC) -c $(FFLAGS) intrgrummp.f
$(LIB):GRUMMP_IF_2D.o intrgrummp.o
	ar $(ARARGS) $@ $?
	mv $@ $(PFL_DIR)/lib/.
	(cd $(PFL_DIR)/lib; $(RANLIB) $@)
clean:
	rm -f *.o $(PFL_DIR)/lib/$(LIB)

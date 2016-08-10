include ../flags/copts.ifc

GMSHDIR = /home/ptyue/mesh/gmsh2.8.4/gmsh-2.8.4-source
GMSHLIBS = -L$(GMSHDIR)/lib -lGmsh
GMSHINC     = -I$(GMSHDIR)/Common -I$(GMSHDIR)/Fltk -I$(GMSHDIR)/Geo -I$(GMSHDIR)/Graphics -I$(GMSHDIR)/Mesh -I$(GMSHDIR)/Numeric -I$(GMSHDIR)/Parser -I$(GMSHDIR)/Plugin -I$(GMSHDIR)/Post -I$(GMSHDIR)/Qt -I$(GMSHDIR)/Solver -I$(GMSHDIR)/lib/Common
SRCC     = tecoutput.cpp extractGMSH.cpp GMSH_interface.cpp
SRCF     = writebgm.f tecpostgmsh.f epc.f mshreorder.f 

OBJF = $(SRCF:.f=.o)
OBJC = $(SRCC:.cpp=.o)
LIB = libgmsh.a

all: $(OBJF) $(OBJC)  $(LIB)

$(LIB): $(OBJF) $(OBJC)
	ar $(ARARGS) $@ $? 
	mv $@ $(PFL_DIR)/lib/.
	(cd $(PFL_DIR)/lib ; $(RANLIB) $@)


.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@ $(GMSHINC) $(GMSHLIBS)
.f.o:
	$(FC) -c $(FFLAGS) $< -o $@

clean:
	rm -f *.o $(PROG)


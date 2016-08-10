#
# partflow library makefile
#
# Need to set environment variables PFL_DIR (partflow directory) and 
# and files.
#
#
# set machine dependent FLAGS
#
include ../flags/copts.ifc
#
#PFL_DIR    = $(HOME)/results/contact_gmsh
LIB   = libflows.a
SRC_F = partbc.f beginc.f shap1d.f shap2d.f forces.f front.f thin.f \
	matnew.f pflow.f  prefro.f chkdat.f stream.f ndloc.f mshmov.f   \
	cvnust.f mesh2d.f bdnrfn.f prtmsh.f solver.f wreslt2.f \
	mshold.f preflw.f pstflw.f randist.f wrstrt.f \
	update.f vortct.f xperid.f intplt.f inibd.f masbal.f \
	collision.f matlen.f poisson.f matsld.f DropDf.f timestep.f \
        matamphi.f rstextra.f matcontact.f boundpre.f dropspread.f
SRC_H =
SRC_C = wreslt.c
OBJ_F = $(SRC_F:.f=.o)
OBJ_C = $(SRC_C:.c=.o)
all : $(LIB)
#
.KEEP_STATE:
#
# define dependencies of include files
#
$(OBJ_F) : $(SRC_H)
#
# create library and put a copy in $(PFL_DIR)/lib
#
.PRECIOUS: $(LIB)
$(LIB): $(OBJ_F) $(OBJ_C)
	ar $(ARARGS) $@ $? 
	mv $@ $(PFL_DIR)/lib/.
	(cd $(PFL_DIR)/lib ; $(RANLIB) $@)
#
clean :
	/bin/rm -f *.o *.a
#cleanall :
#	/bin/rm -f *.o *.a $(PFL_DIR)/lib/$(LIB)

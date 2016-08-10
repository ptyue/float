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
LIB   = libblas.a
SRC_F=    chpr.f    dcabs1.f  dsymm.f   izamax.f  ssbmv.f   strsv.f   zherk.f \
caxpy.f   crotg.f   dcopy.f   dsymv.f   lsame.f   sscal.f   xerbla.f  zhpmv.f \
ccopy.f   cscal.f   ddot.f    dsyr2.f   Makefile  sspmv.f   zaxpy.f   zhpr2.f \
cdotc.f   csrot.f   dgbmv.f   dsyr2k.f  sasum.f   sspr2.f   zcopy.f   zhpr.f \
cdotu.f   csscal.f  dgemm.f   dsyr.f    saxpy.f   sspr.f    zdotc.f   zrotg.f \
cgbmv.f   cswap.f   dgemv.f   dsyrk.f   scasum.f  sswap.f   zdotu.f   zscal.f \
cgemm.f   csymm.f   dger.f    dtbmv.f   scnrm2.f  ssymm.f   zdrot.f   zswap.f \
cgemv.f   csyr2k.f  dnrm2.f   dtbsv.f   scopy.f   ssymv.f   zdscal.f  zsymm.f \
cgerc.f   csyrk.f   drot.f    dtpmv.f   sdot.f    ssyr2.f   zgbmv.f   zsyr2k.f \
cgeru.f   ctbmv.f   drotg.f   dtpsv.f   sdsdot.f  ssyr2k.f  zgemm.f   zsyrk.f \
chbmv.f   ctbsv.f   drotm.f   dtrmm.f   sgbmv.f   ssyr.f    zgemv.f   ztbmv.f \
chemm.f   ctpmv.f   drotmg.f  dtrmv.f   sgemm.f   ssyrk.f   zgerc.f   ztbsv.f \
chemv.f   ctpsv.f   dsbmv.f   dtrsm.f   sgemv.f   stbmv.f   zgeru.f   ztpmv.f \
cher2.f   ctrmm.f   dscal.f   dtrsv.f   sger.f    stbsv.f   zhbmv.f   ztpsv.f \
cher2k.f  ctrmv.f   dsdot.f   dzasum.f  snrm2.f   stpmv.f   zhemm.f   ztrmm.f \
cher.f    ctrsm.f   dspmv.f   dznrm2.f  srot.f    stpsv.f   zhemv.f   ztrmv.f \
cherk.f   ctrsv.f   dspr2.f   icamax.f  srotg.f   strmm.f   zher2.f   ztrsm.f \
chpmv.f   dasum.f   dspr.f    idamax.f  srotm.f   strmv.f   zher2k.f  ztrsv.f \
chpr2.f   daxpy.f   dswap.f   isamax.f  srotmg.f  strsm.f   zher.f

SRC_H =
#SRC_C = wreslt.c
OBJ_F = $(SRC_F:.f=.o)
OBJ_C = $(SRC_C:.c=.o)
#
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

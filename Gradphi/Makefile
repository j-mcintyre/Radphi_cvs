# gnumake

CPP_Linux = cc -E
CPP_SunOS = /usr/ucb/cc -E

OS := $(shell uname)
CPP = $(CPP_$(OS))

MAIN = gradphi.o

#-------------------------------------------------------------------
# The symbol SPECIAL is used to active any of the following switches
# to enable custom versions of the simulation program.
#SW_RATES = -DSAVE_RATES_NTUPLE
#SW_BREMS = -DDIAGNOSE_BREMSSTRAHLUNG -DWRITE_ALL
#SW_BREMS = -DDIAGNOSE_BREMSSTRAHLUNGC=10
#SW_CEREN = -DFOLLOW_CERENKOV -DEARLY_DETECTION
#SW_CEREN = -DFOLLOW_CERENKOV=1
#SW_CEREN = -DOLD_CERENKOV_N0
#SW_SPECIAL = -DWRITE_ALL -DLGD_FOR_GLUEX
SW_SPECIAL = -DWRITE_ALL
#-------------------------------------------------------------------

SPECIAL = $(SW_RATES) $(SW_BREMS) $(SW_CEREN) $(SW_SPECIAL)

AUXS =

OBJS = gustep.o gukine.o uginit.o uglast.o guout.o gsvertx.o\
       major.o matdef.o meddef.o uhinit.o gutrev.o \
       detector.o lgd.o lgpmt.o rpd.o cpv.o upv.o trgt.o placelgd.o \
       placerpd.o placecpv.o placeupv.o placetrgt.o sensitive.o angles.o \
       placebrad.o placeclas.o placecols.o placeepip.o placegtub.o \
       placehole.o placewall.o placeweep.o beampipe.o beamtube.o \
       clasmgnt.o clastgt.o collim.o herccover.o radiator.o leadwall.o \
       shield.o sweepmgnt.o dnut.o placednut.o placecover.o \
       lgdframe.o placeframe.o platform.o placeplatform.o \
       savehits.o upvhits.o rpdhits.o cpvhits.o lgdhits.o triggered.o\
       inputEvents.o outputEvents.o LGDsave.o RPDsave.o CPVsave.o UPVsave.o \
       sort.o timeAnalysis.o bsd.o bgv.o bsdhits.o bgvhits.o BSDsave.o \
       timeAnalysisBSD.o placebsd.o placebgv.o BGVsave.o timeAnalysisCPV.o \
       ratestudy.o cerenkov.o background.o

FIXES = gxphys.o timel.o

INCS = geometry.inc hitconstants.inc ratestudy.inc cerenkov.inc

LDOPTS = -Wl,-uuginit_,-uuglast_,-ugutrev_,-ugustep_,-uguout_,-ugukine_

LDXOPTS = -Wl,-uuginit_,-uuglast_,-ugutrev_,-ugustep_,-uguout_,-ugukine_,-ugxphys_,-ugxcs_,-export-dynamic

LIBFILE = $(LIBDIR)/libgradphi.a

include $(RADPHI_HOME)/makefiles/rules.make

#CFLAGS = -O $(SPECIAL)
CFLAGS = -g $(SPECIAL)
#FFLAGS = -O $(SPECIAL) -DCERNLIB_MOTIF -D_FILE_OFFSET_BITS=64
FFLAGS = -g $(SPECIAL) -DCERNLIB_MOTIF -D_FILE_OFFSET_BITS=64 -DCERNLIB_CG

F77 = gfortran

all: gradphi++.x 
lib: $(LIBFILE)
exe:
exe-install:

gradphi++.x:	gxint.o $(AUXS) $(INCS) $(LIBFILE)
	gxint -v $(CERN_LEVEL) -g 321 -m -d Motif -o $@ -L$(LIBDIR) \
	 	-ldata -lpthread -- \
		$(LDXOPTS) $(FFLAGS) gxint.o $(AUXS) $(LIBFILE) \

gradphi.x:	$(MAIN) $(AUXS) $(INCS) $(LIBFILE)
	$(F77)	$(FFLAGS) -o $@ $(LDOPTS) $(MAIN) $(AUXS) \
		$(LIBFILE) -L$(LIBDIR) -ldata -lpthread \
		`$(CERNLIB) geant321 pawlib mathlib graflib grafX11` \
		$(NETLIB) $(XLIBS)

timel.o:	timel.c
	$(CC) $(CFLAGS) -c $(INCLUDE) -I${CERN_ROOT}/include/ $^

gxint.o:	gxint.F
	$(F77) $(FFLAGS) -c $(FINCLUDE) -I${CERN_ROOT}/include/ $^

$(LIBFILE): $(OBJS) $(FIXES)
	$(AR) rv $@ $?

.F.o:	$(INCS)
	$(F77) $(FFLAGS) -c $(FINCLUDE) -I${CERN_ROOT}/include/ $^

.c.o:	$(INCS) outputEvents.h

dumpEvent:	dumpEvent.c 
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ -L$(LIBDIR) -ldata -lpthread

cleaner:
	find . -name paw.metafile -exec rm {} \;
	find . -name last.kumac -exec rm {} \;
	find . -name "fort.*" -exec rm {} \;
	find . -name geant.hbook -exec rm {} \;
	find . -name "*.o" -exec rm {} \;

work:
	touch ${WORK}/events.in

save:
	@rm ${WORK}/simData.itape ${WORK}/control.in ${WORK}/events.in \
	    ${WORK}/*.x


#! gnumake


BBOOK = -lwabbit
GEOM = -lLGD
MAP = -lparam -lmapmanager
DATA = -lUtil -ldata -lpthread
MATH = -lm
CHEES = -lChees
OBJS_libChees = Chees.o  Input.o Particle.o Four.o Kinmaster.o NosPicker.o Box.o Keeper.o SqMatrix3.o Three.o util.o gen.o collide.o loadParam.o ReadUserIntensity.o CheesCloth.o
INCLUDE_FILES = Box.h Input.h Keeper.h Kinmaster.h NosPicker.h Particle.h collide.h loadParam.h targets.h util.h ReadUserIntensity.h 

OBJS_chees = main.o
LIBS_chees = $(CHEES)  $(BBOOK)  $(MAP)
LIBFILE_local = ./libChees.a
LIBS_chees_local = $(LIBFILE_local) $(GEOM) $(BBOOK) $(DATA) $(MAP) $(FLIB) $(MATH)

OBJS_YLM = createYLM.o
LIBS_YLM = $(BBOOK) $(MATH)

OBJS_CDF = createCDF.o
LIBS_CDF = $(BBOOK)

LIBFILE = $(LIBDIR)/libChees.a

#OBJS_CHECKBEAMPICKER= checkBeamPicker.o
#LIBS_CHECKBEAMPICKER= E852beamPicker.o  $(CHEES) $(BBOOK) $(DATA) $(MAP) $(FLIB) $(MATH)

OBJS_BEAMCHECK= beamCheck.o
LIBS_BEAMCHECK=  $(BBOOK) $(DATA) $(MAP) $(FLIB) $(MATH)  	 

all: $(INCLUDE_FILES) lib Chees createYLM createCDF

include $(RADPHI_HOME)/makefiles/rules.make

exe:
exe-install:


lib: $(INCLUDE_FILES) $(LIBFILE)
libLocal: $(INCLUDE_FILES) $(LIBFILE_local)  

$(LIBFILE): $(OBJS_libChees)
	$(AR) rv $@ $?

$(LIBFILE_local): $(OBJS_libChees)
	ar rv $@ $?

mc: libChees.a Cebaf

cebaf.o: targetRegion.h

Cebaf: cebaf.o util.o $(CHEES) -lLGD $(MAP) $(BBOOK) $(GEOM) $(DATA)
	time $(CXX) $(ARCHFLAGS) -o $@ $^ -lm -lX11

Chees: $(OBJS_chees) $(LIBS_chees) 
	time $(CXX) $(ARCHFLAGS) $(CCFLAGS) -o $@ $^ -lm -lX11

CheesLoCal: $(OBJS_chees) $(LIBS_chees_local) 
	time $(CXX) $(CCFLAGS) -o $@ $^ 

createYLM: $(OBJS_YLM) $(LIBS_YLM)
	time cc $(CFLAGS) -o $@ $^

createCDF: $(OBJS_CDF) $(LIBS_CDF)
	time cc $(CFLAGS) -o $@ $^

checkBeamPicker: $(OBJS_CHECKBEAMPICKER) $(LIBS_CHECKBEAMPICKER)
	time $(CXX) $(CCFLAGS) -o $@ $^ 

beamCheck: $(OBJS_BEAMCHECK) $(LIBS_BEAMCHECK)
	time cc $(CFLAGS) -o $@ $^

mcwrap:		mcwrap.o fffnmc.o genFermi.o
	#$(F77) -g -ftrap=%all -static -o $@ $^ -I../include \
	$(F77) -g -o $@ $^ -I../include \
	-L$(LIBDIR) $(DATA) \
	`$(CERNLIB) mathlib packlib` \
	$(NETLIB)

mcwrap.o:	mcwrap.c
	$(CC) -g -c $(ARCHFLAGS) $(CFLAGS) $^ -I../include

fffnmc.o:	fffnmc.F
	$(F77) -c -g -fno-automatic $(FINCLUDE) $^

genFermi.o:	genFermi.F
	$(F77) -c -g -fno-automatic $(FINCLUDE) $^

-include *.depend

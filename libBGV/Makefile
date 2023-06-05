#! gnumake

OBJS = bgvClusterize.o makePhotonBGV.o MCsmearBGV.o
LIBFILE = $(LIBDIR)/libBGV.a


include $(RADPHI_HOME)/makefiles/rules.make

all: lib

exe:
exe-install:

lib: $(LIBFILE)

$(LIBFILE): $(OBJS)
	$(AR) rv $@ $?

-include *.depend

makePhotonBGV.o : -lm

MCsmearBGV.o : -lm

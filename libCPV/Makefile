#! gnumake

OBJS = cpvGeom.o
LIBFILE = $(LIBDIR)/libCPV.a


include $(RADPHI_HOME)/makefiles/rules.make

all: lib

exe:
exe-install:

lib: $(LIBFILE)

$(LIBFILE): $(OBJS)
	$(AR) rv $@ $?

-include *.depend

#! gnumake

OBJS = lgdGeom.o lgdClusterIU.o lgdClusterUC.o lgdfGeom.o lgdUtil.o\
       esum.o lgdHole.o
LIBFILE = $(LIBDIR)/libLGD.a


include $(RADPHI_HOME)/makefiles/rules.make

all: lib

exe:
exe-install:

lib: $(LIBFILE)

$(LIBFILE): $(OBJS)
	$(AR) rv $@ $?

-include *.depend

#! gnumake

LIB_OBJS = average.o linearfit.o cern_X.o cern_fit.o cern_noX.o cern_minuit.o particles.o mtime.o mwait.o encoderHandle.o

all: lib exe exe-install

include $(RADPHI_HOME)/makefiles/rules.make

LIBFILE = $(LIBDIR)/libUtil.a


lib: $(LIBFILE)

$(LIBFILE): $(LIB_OBJS)
	$(AR) rv $@ $?

exe:
exe-install:




-include *.depend





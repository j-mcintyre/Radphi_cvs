#! gnumake

all: Event

include $(RADPHI_HOME)/makefiles/rules.make

OBJS = main.o Event.o init.o feedback.o setup.o draw.o input.o

Event:  $(OBJS) -lHv -lmakehits  -lLGD -lBGV -ldata -lUtil -lparam -lmapmanager 
	$(F77) $(ARCHFLAGS) -o $@ $^ $(NETLIB) $(XLIBS) -lm -lpthread \
        -L$(CERN_ROOT)/lib -lmathlib -lkernlib -lstdc++


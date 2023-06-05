#! gnumake


all: rodd

include $(RADPHI_HOME)/makefiles/rules.make

rodd : main.o commands.o runs.o errors.o parse.o process.o -lmakehits -lLGD -lBGV -ldata -lUtil -lparam -lmapmanager -lclib
	$(CC) -g $(ARCHFLAGS) -o $@ $^ \
	`$(CERNLIB) packlib mathlib` \
	$(XLIBS) $(FLIB) $(NETLIB) -lpthread -lstdc++


-include *.o.depend

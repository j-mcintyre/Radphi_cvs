#! gnumake

all: lib exe exe-install

include $(RADPHI_HOME)/makefiles/rules.make

LDOPTS :=

lib:
exe: readRaw dumpEvent
exe-install: ntuplet mctuplet


readRaw: readRaw.o -ldata
	$(CC) $(ARCHFLAGS) $(CFLAGS) -o $@ $^ \
	-lpthread

readTest: readTest.o -ldata
	$(CC) $(ARCHFLAGS) $(CFLAGS) -o $@ $^ \
	-lpthread

dumpEvent: dumpEvent.o -ldata
	$(CC) $(ARCHFLAGS) $(CFLAGS) -o $@ $^ \
	-lpthread

sampleAnalysis: sampleAnalysis.o -lLGD -lBGV -ldata -lUtil -lclib -lparam -lmapmanager -lCPV
	$(CC) $(ARCHFLAGS) -o $@ $^ \
	`$(CERNLIB) mathlib packlib` \
	 $(XLIBS) $(NETLIB) $(FLIB) \
	-lpthread

sampleMC: sampleMC.o -lLGD -lBGV -ldata -lUtil -lclib -lparam -lmapmanager -lCPV
	$(CC) $(ARCHFLAGS) -o $@ $^ \
	`$(CERNLIB) mathlib packlib` \
	 $(XLIBS) $(NETLIB) $(FLIB) \
	-lpthread

plotMC: plotMC.o -lLGD -lBGV -ldata -lUtil -lclib -lparam -lmapmanager -lCPV
	$(CC) $(ARCHFLAGS) -o $@ $^ \
	`$(CERNLIB) mathlib packlib` \
	 $(XLIBS) $(NETLIB) $(FLIB) \
	-lpthread

scalers: scalers.c scdef.c -lUtil -ldata -lparam -lmapmanager
	$(CC) $(ARCHFLAGS) $(CFLAGS) $(DLLFLAG) \
	-I../include -o scalers.sl $^ -lpthread

ntuple: ntuple.o ntdef.o hidef.o merge_clusters.o -lmakehits -lBGV -lLGD -lUtil -ldata -lparam -lmapmanager
	$(CC) $(ARCHFLAGS) $(CFLAGS) $(DLLFLAG) -o ntuple.sl $^ \
	-lpthread -lstdc++

ntuplet: ntuplet.F ntuple.o ntdef.o hidef.o merge_clusters.o -lmakehits -lBGV -lLGD -lUtil -ldata -lparam -lmapmanager -lclib
	$(F77) $(ARCHFLAGS) $(FFLAGS) $(LDOPTS) -g -o ntuplet $^ \
	`$(CERNLIB) pawlib graflib grafX11 packlib mathlib` \
	$(NETLIB) $(XLIBS) \
	-lpthread -lstdc++

mctuple: mctuple.o mcdef.o hidef.o merge_clusters.o -lmakehits -lBGV -lLGD -lUtil -ldata -lparam -lmapmanager -lclib
	$(CC) $(ARCHFLAGS) $(CFLAGS) $(DLLFLAG) -o mctuple.sl $^ \
	-lpthread -lstdc++

mctuplet: mctuplet.F mctuple.o mcdef.o hidef.o merge_clusters.o -lmakehits -lBGV -lLGD -lUtil -ldata -lparam -lmapmanager -lclib
	$(F77) $(ARCHFLAGS) $(FFLAGS) $(LDOPTS) -g -o mctuplet $^ \
	`$(CERNLIB) pawlib graflib grafX11 packlib mathlib` \
	 $(NETLIB) $(XLIBS) \
	 -lpthread -lstdc++

psfReduce: psfReduce.o -lmakehits -lLGD -lBGV -ldata -lUtil -lparam -lmapmanager -lCPV -lclib
	$(CC) $(ARCHFLAGS) -o $@ $^ \
	`$(CERNLIB) mathlib packlib`  \
	$(XLIBS) $(NETLIB) $(FLIB) \
	-lpthread -lstdc++

work:

save:

-include *.o.depend


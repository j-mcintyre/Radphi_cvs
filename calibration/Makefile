#! gnumake

all:

include $(RADPHI_HOME)/makefiles/rules.make


lgdtune: lgdtune.o -lmakehits -lLGD -lBGV -ldata -lUtil -lparam -lmapmanager -lclib -llapack3 -lblas
	$(CC) $(ARCHFLAGS) $(CDBG) -o $@ $^ \
	`$(CERNLIB) packlib mathlib` \
	$(XLIBS) $(FLIB) -lpthread

extractPedestals: extractPedestals.o -lLGD -ldata -lUtil -lparam -lmapmanager -lpthread
	$(CC) $(ARCHFLAGS) $(COPT) -o $@ $^ $(XLIBS) -lm

lgdCalibrate: lgdCalibrate.o  cal784.o calWrapper.o  -lLGD -lmakehits -lBGV -ldata -lUtil -lparam -lmapmanager -lclib -lpacklib
	f77 $(ARCHFLAGS) -o  $@ $^ $(XLIBS) $(FLIB) -lpthread

setNewVoltages: setNewVoltages.o -lm -lclib
	f77 $(ARCHFLAGS) -o  $@ $^ $(XLIBS) $(FLIB)

ccToMap: ccToMap.o -ldata -lmapmanager -lpthread
	f77 $(ARCHFLAGS) -o  $@ $^ $(XLIBS) $(FLIB)

dan3: dan3.o  -lLGD -ldata -lUtil -lparam -lmapmanager -lclib -lpacklib 
	f77 $(ARCHFLAGS) -o  $@ $^ $(XLIBS) $(FLIB) -lpthread

pedestalsToMap: pedestalsToMap.o -ldata -lmapmanager
	$(CC) $(ARCHFLAGS) $(COPT) -o $@ $^ -lpthread

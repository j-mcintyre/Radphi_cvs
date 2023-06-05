#! gnumake

all: extractPedestals

include $(RADPHI_HOME)/makefiles/rules.make

set_adc_trims: set_adc_trims.o 
	$(CC) $(ARCHFLAGS) $(COPT) -o $@ $^ -lX11 -lm

extractPedestals: extractPedestals.o -lLGD -ldata -lUtil -lparam -lmapmanager
	$(CC) $(ARCHFLAGS) $(COPT) -o $@ $^ -lX11 -lm

#lgdCalibrate: lgdCalibrate.o cal784.o calWrapper.o ../tapeWriter/libwrapper.o -lLGD -ldata -lUtil -lparam -lmapmanager -lclib -lpacklib
lgdCalibrate: lgdCalibrate.o cal784.o calWrapper.o -lLGD -ldata -lUtil -lparam -lmapmanager -lclib -lpacklib
	f77 $(ARCHFLAGS) -o  $@ $^ -lX11 $(FLIB) 

dan3: dan3.o  -lLGD -ldata -lUtil -lparam -lmapmanager -lclib -lpacklib
	f77 $(ARCHFLAGS) -o  $@ $^ -lX11 $(FLIB)

pedestalsToMap: pedestalsToMap.o -ldata -lmapmanager
	$(CC) $(ARCHFLAGS) $(COPT) -o $@ $^

#! gnumake

#all: gap_2k

GLOBAL := 1
include $(RADPHI_HOME)/makefiles/rules.make

data_put: data_put.o -lLGD -lUtil -lmakehits -ldata -lparam -lmapmanager -lBGV -lclib -lpacklib
	gfortran $(ARCHFLAGS) -o  $@ $^ -lX11 $(FLIB)


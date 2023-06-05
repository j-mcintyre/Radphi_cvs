all:
	cd clib && make clean && make lib
	cd libdata && make clean && make lib
	cd mapmanager && make clean && make && make exe-install
	cd libParam && make clean && make lib && make exe-install
	cd libUtil && make clean && make lib
	cd libBGV && make clean && make lib
	cd libCPV && make clean && make lib
	cd libHv && make clean && make lib
	cd libLGD && make clean && make lib
	cd taggerdata && make clean && make
	cd makehits && make clean && make lib
	cd tapeWriter && make clean && make
	cd Event && make clean && make && make
	cd Gradphi && make clean && make all
	#cd MonteCarlo && make clean && make
	cd calibration && make clean && make
	cd RODD && make clean && make
	cd Examples && make clean && make

clean:
	find $(RADPHI_HOME) -name "*.o" -exec rm {} \;
	find $(RADPHI_HOME) -name "*.a" -exec rm {} \;
	find $(RADPHI_HOME) -name "*.depend" -exec rm {} \;

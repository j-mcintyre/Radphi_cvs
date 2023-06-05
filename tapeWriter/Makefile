#! gnumake


all: lib exe exe-install


include $(RADPHI_HOME)/makefiles/rules.make

LOPT := 

lib:
exe: wrapFile labelTape countEvents compact partition
exe-install: wrapFile.install compact.install partition.install

wrapFile: wrapFile.o -ldata
	$(CC) $(ARCHFLAGS) $(COPT) $(LOPT) -o $@ $^ -lpthread

labelTape: labelTape.o -ldata -lUtil
	$(CC) $(ARCHFLAGS) $(COPT) $(LOPT) -o $@ $^ -lpthread

countEvents: countEvents.o -ldata
	$(CC) $(ARCHFLAGS) $(COPT) $(LOPT) -o $@ $^ -lpthread

seekItape: seekItape.o -ldata
	$(CC) $(ARCHFLAGS) $(COPT) $(LOPT) -o $@ $^ -lpthread

dumpDirToTapes: dumpDirToTapes.o -ldata
	$(CC) $(ARCHFLAGS) $(COPT) $(LOPT) -o $@ $^ -lpthread

feedFromDD: feedFromDD.o -ldd -ldata -lUtil -lmsql
	$(CC) $(ARCHFLAGS) $(LOPT) -o $@ $^ $(NETLIB) $(FLIB) -lpthread

compact: compact.o -lmakehits -lLGD -lBGV -ldata -lUtil -lclib -lparam -lmapmanager
	$(CC) $(ARCHFLAGS) $(COPT) $(LOPT) -o $@ $^ \
	`$(CERNLIB) packlib` \
	 $(XLIBS) $(NETLIB) $(FLIB) -lpthread -lstdc++

partition: partition.o -ldata -lUtil -lclib -lparam -lmapmanager
	$(CC) $(ARCHFLAGS) $(COPT) $(LOPT) -o $@ $^ \
        `$(CERNLIB) packlib` \
	 $(XLIBS) $(NETLIB) $(FLIB) -lpthread -lstdc++

-include *.o.depend

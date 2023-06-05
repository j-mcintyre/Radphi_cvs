#! gnumake

F77 := gfortran

RADPHI_LIBS = $(RADPHI_HOME)/lib.$(RADPHI_MAKETARGET)
CERN_LIBS = $(CERN_ROOT)/lib

OBJS = fileIO.c.o eventIO.c.o fitImage.c.o psf.c.o setupHbook.o

.F.o:
	$(F77) -g -c $^

psfit.sl: psfit.F $(OBJS)
	$(F77) -g -shared -o $@ $^ -L$(RADPHI_LIBS) -lmakehits -lLGD -lBGV -lUtil -ldata -lparam -lmapmanager -lpthread \
	`cernlib packlib`

fitImage.sl:	fitImage.c.o eventIO.c.o fileIO.c.o psf.c.o
	$(F77) -g -shared -o $@ $^ -L$(RADPHI_LIBS) -lmakehits -lLGD -lBGV -lUtil -ldata -lparam -lmapmanager -lpthread 

psf.c.o:	psf.c psf.h
	cc -g -c -I. -I$(RADPHI_HOME)/include -o $@ $<

fitImage.c.o:	fitImage.c fitImage.h
	cc -g -c -I. -I$(RADPHI_HOME)/include -o $@ $<

eventIO.sl:	eventIO.c.o fileIO.c.o
	$(F77) -g -shared -o $@ $^ -L$(RADPHI_LIBS) -lmakehits -lLGD -lBGV -lUtil -ldata -lparam -lmapmanager -lpthread 

fileIO.sl:	fileIO.c.o
	$(F77) -g -shared -o $@ $^ -L$(RADPHI_LIBS) -ldata -lpthread

fileIO.c.o:	fileIO.c fileIO.h
	cc -g -c -I. -I$(RADPHI_HOME)/include -o $@ $<

eventIO.c.o:	eventIO.c  eventIO.h eventIOfriends.c lgdStrings.c
	cc -g -c -I. -I$(RADPHI_HOME)/include -o $@ $<

pict5: 
	latex pict5.tex
	dvips -o pict5.ps pict5.dvi


5deg:
	rm file.in
	ln -s file_t05.in file.in

10deg:
	rm file.in
	ln -s file_t10.in file.in

15deg:
	rm file.in
	ln -s file_t15.in file.in

20deg:
	rm file.in
	ln -s file_t20.in file.in


cleanex:
	rm *.dvi *.log *.aux

clean:
	rm *.o *.sl 

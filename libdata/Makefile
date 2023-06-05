#! gnumake
all: lib
GLOBAL:=1
include $(RADPHI_HOME)/makefiles/rules.make

clean: 
	rm -f *.o *.oo *.depend
	rm -f endianConverterAuto.c
	@echo " "
	@echo "NOTE: the warnings about 'over-riding commands'"
	@echo "and 'ignoring old commands' are a feature"
	@echo "see Makefile comments for details -CPS "
	@echo " "
# endianConverterAuto.c is a perl script generated c file.  It 
# takes care of all the endian conversions within an itape file,
# by parsing through the include files to figure out what the 
# structure of each type of group is (So Richard tells me).  It is
# apparently not in the depencies in the make file correctly 
# (possibly it cannot be for technical reasons) but it has to be 
# removed before the make, so that it is freshly generated.

LIBFILE = $(LIBDIR)/libdata.a

OBJS = data.o dataIO.o evt_dataIO.o tapeIO.o crc.oo tapeIoctl.o endianConverterAuto.o endianConverter.o evio.o swap_util.o unpackMap.o disIO.o dataSeekTape.o unpackData.o libdecode.o thawData.o pedestals.o tapedb.o calibration.o libwrapper.o detector_params.o fitting.o

lib: $(LIBFILE)
exe:
exe-install:

$(LIBFILE): $(OBJS)
	$(AR) rv $@ $?

test_put: test_put.o -ldata -lUtil -lparam -lmapmanager
	$(CC) $(ARCHFLAGS) $(COPT) -o $@ $^ -lX11 -lm

endianConverterAuto.o: endianConverterAuto.c

endianConverterAuto.c: endianStubGen.perl
	$(PERL) endianStubGen.perl

endianConverter.o: endianConverterAuto.c



#
# Include all the .depend files
#
-include *.depend
#
# end file
#

push: 
	cd ../Examples;\
	gmake;\
        ftp jlabs3.jlab.org;\
        cd ../libdata

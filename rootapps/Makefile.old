#include $(RADPHI_HOME)/rootclasses/rules.make
CC := g++
F77 := gfortran

.SUFFIXES:	.so

.f.o:;	@echo making $@
	$(F77) -c $<

.C.o:;	@echo making $@
	@CC -c $(CDBFLAGS) -I$(ROOTINC) -I$(RADPHI_HOME)/include $<

.c.o:;	@echo making $@
	@cc -c $(CDBFLAGS) $<

.o.so:;	@echo "Building" $@
	@$(LD) $(SOFLAGS) $< -o $@
	@echo "done"

h2root: h2root.cxx hbook.o cernlib.o kernlib.o zebra.o
	$(CC) -fPIC -I$(shell root-config --incdir) -o $@ $^ \
	`root-config --libs` -l crypt -lgfortran

myroot: myroot.o
	@echo "Linking $@ ..."
	@$(LD) -o $@ $^ $(LDFLAGS) \
	-L$(RADPHI_HOME)/rootclasses -lrootclasses \
	-R$(RADPHI_HOME)/rootclasses \
	$(GLIBS)
	@echo "done"

dumpEvent: dumpEvent.o
	@echo "Linking $@ ..."
	@$(LD) -o $@ $^ $(LDFLAGS) $(GLIBS)
	@echo "done"

itape2tree: itape2tree.o
	@echo "Linking $@ ..."
	@$(LD) -o $@ $^ $(LDFLAGS) \
	-L$(RADPHI_HOME)/rootclasses -lrootclasses \
	-R$(RADPHI_HOME)/rootclasses \
        $(GLIBS)
	@echo "done"

clean:;	@rm -f *.o core

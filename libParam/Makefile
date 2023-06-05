#! gnumake



OBJS = param.o paramMM.o
LIBFILE = $(LIBDIR)/libparam.a

all:  lib

include $(RADPHI_HOME)/makefiles/rules.make

lib: $(LIBFILE)
exe: xrdbToMap xrdbFromMap
exe-install: xrdbToMap.install xrdbFromMap.install


$(LIBFILE): $(OBJS)
	$(AR) rv $@ $?

xrdbToMap: xrdbToMap.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS) $(CFLAGS) -o $@ $^ -lpthread

xrdbFromMap: xrdbFromMap.o -lmapmanager -ldata
	$(CC) $(ARCHFLAGS) $(CFLAGS) -o $@ $^ -lpthread

#
# Include all the .depend files
#
-include *.depend
#
# end file
#

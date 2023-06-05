#! gnumake

MYFLAGS = -DLGD_CLUSTERIZER_UC
#MYFLAGS = -DLGD_CLUSTERIZER_IU

include $(RADPHI_HOME)/makefiles/rules.make

OBJS = makeHits.o setup_makeHits.o make_lgd_hits.o make_bsd_hits.o \
       make_photons.o  make_lgd_clusters.o make_bsd_pixels.o \
       make_bgv_hits.o make_cpv_hits.o make_upv_hits.o make_mesons.o\
       make_tagger_hits.o make_time_list.o makePhotonPatch.o \
       setup_global_utility_structures.o make_recoil.o make_tagger_clusters.o \
       make_tagger_photons.o make_bgv_clusters.o

LIBFILE = $(LIBDIR)/libmakehits.a

CFLAGS := $(CFLAGS) $(MYFLAGS)
ACPP   := $(ACPP) $(MYFLAGS)

all: lib

exe:
exe-install:

lib: $(LIBFILE)
	echo $(LIBDIR)

$(LIBFILE): $(OBJS)
	$(AR) rv $@ $?

-include *.depend


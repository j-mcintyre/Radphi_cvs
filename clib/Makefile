#! gnumake
#
# Makefile for clib
#

OBJS= flip.o flip_integer.o flip_short.o flipint.o flipshort.o \
      gauss.o \
      BetaGamma.o \
      DotProduct3.o CrossProduct3.o \
      EffMass.o \
      add3.o add4.o \
      v3sub.o v4sub.o \
      nanCheck.o \
      lorentz.o make4vector.o \
      negatime.o p3vector.o p4vector.o p4vectors.o sum4vectors.o \
      sum3vectors.o v3diff.o v3mag.o v3magsq.o v3add.o \
      v3sum.o v4add.o v4diff.o v4mag.o v4magsq.o v4sum.o \
      nPartMass.o scalerMult4V.o n4VectSum.o Ctime.o effMass.o \
      ScalerMult.o ctheta.o

all: lib exe exe-install

include $(RADPHI_HOME)/makefiles/rules.make


CLIBFILE= $(LIBDIR)/libclib.a


lib: $(CLIBFILE)
exe:
exe-install:

$(CLIBFILE):  $(OBJS)
	$(AR) rv $(CLIBFILE) $?

-include *.depend

# end file








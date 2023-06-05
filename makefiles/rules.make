#
# Rules for Makefiles
#


#
# List of supported targets-
# irix-mips4	IRIX 6.2, n32 mode, mips4, debugging
# irix-mips4-o	IRIX 6.2, n32 mode, mips4, optimized
# irix-mips3    IRIX 6.2, n32 mode, mips3, debugging
# irix-mips3-o  IRIX 6.2, n32 mode, mips3, optimized
# irix-mips2    IRIX 6.2, mips2, debugging
# irix-mips2-o  IRIX 6.2, mips2, optimized
# hpux		HPUX 10.x, debugging
# hpux-o	HPUX 10.x, optimized
# solaris	Solaris 2.5 (== SunOS 5.5), debugging
# solaris-o	Solaris 2.5,(== SunOS 5.5), optimized
# dec		DEC OSF/1 4.x, debugging
# dec-o		DEC OSF/1 4.x, optimized
# aix		AIX 4.1, debugging
# aix-o		AIX 4.1, optimized
# linux5	Redhat Linux 5.x, debugging
# linux5-o	Redhat Linux 5.x, optimized
# linux		Redhat Linux 6.x, debugging
# linux-o	Redhat Linux 6.x, optimized
# vxworks-ppc   VxWorks 5.3 

PWD := $(shell pwd)

#
# GLIB- by default, global libraries (~radphi/lib.xxx) are not modified
#

GLIB := 0

#
# XLIBS defines the location of the X11 libraries
# Do NOT put on the dependencies list
#

XLIBS_irix-mips4	:= -lFm_c -lXm_s -lXt_s -lXext -lX11_s -lmalloc -lPW
XLIBS_irix-mips4-o 	:= $(XLIBS_irix-mips4)
XLIBS_irix-mips3 	:= $(XLIBS_irix-mips4)
XLIBS_irix-mips3-o 	:= $(XLIBS_irix-mips4)
XLIBS_irix-mips2       	:= $(XLIBS_irix-mips4)
XLIBS_irix-mips2-o     	:= $(XLIBS_irix-mips4)
XLIBS_dec		:= -lFm_c -lXm -lXt -lX11 -lPW
XLIBS_dec-o		:= $(XLIBS_dec)
XLIBS_aix 		:= -lFm_c -lXm -lXt -lX11 -lPW
XLIBS_aix-o		:= $(XLIBS_rs6000)
XLIBS_hpux		:= -lsomeone_figure_it_out
XLIBS_hpux-o		:= $(XLIBS_hpux)
XLIBS_solaris	       	:= -lXm -lXt -lX11 -lw
XLIBS_solaris-o		:= $(XLIBS_solaris)
XLIBS_linux5	       	:= -lXp -lXpm -lSM -lXm -lXt -lICE -lXext -lX11
XLIBS_linux5-o		:= $(XLIBS_linux5)
XLIBS_linux	       	:= -lXpm -lXm -lXmu -lXt -lSM -lICE -lXext -lX11 -lXp -ldl
ifdef IU
XLIBS_linux             := -lXpm -lXt -lSM -lICE -lXext -lX11
endif
XLIBS_linux-o		:= $(XLIBS_linux)


XLIBS:= $(XLIBS_$(RADPHI_MAKETARGET))

PEXLIBS= -lXmu -lPEX5

#
# INCLUDE defines the location of the C and C++ include files
#

BASE_INCLUDE 		:= -I. -I../include -I$(RADPHI_HOME)/include -I$(RADPHI_HOME)/include/CC -DFUNCPROTO=15 
INCLUDE_irix-mips4 	:= $(BASE_INCLUDE) -DIRIX
INCLUDE_irix-mips4-o 	:= $(INCLUDE_irix-mips4)
INCLUDE_irix-mips3 	:= $(INCLUDE_irix-mips4)
INCLUDE_irix-mips3-o 	:= $(INCLUDE_irix-mips4)
INCLUDE_irix-mips2     	:= $(INCLUDE_irix-mips4)
INCLUDE_irix-mips2-o 	:= $(INCLUDE_irix-mips4)
INCLUDE_hpux		:= $(BASE_INCLUDE) -DHPUX 
INCLUDE_hpux-o		:= $(INCLUDE_hpux)
INCLUDE_dec		:= $(BASE_INCLUDE) -DDEC
INCLUDE_dec-o		:= $(INCLUDE_dec)
INCLUDE_aix		:= $(BASE_INCLUDE) -DAIX
INCLUDE_aix-o		:= $(INCLUDE_aix)
INCLUDE_solaris		:= $(BASE_INCLUDE) -DSOLARIS -I/usr/openwin/include -I/usr/dt/include -I${CODA}/common/include
INCLUDE_solaris-o	:= $(INCLUDE_solaris)
INCLUDE_linux5		:= $(BASE_INCLUDE) -DLINUX
INCLUDE_linux5-o	:= $(INCLUDE_linux5)
INCLUDE_linux		:= $(BASE_INCLUDE) -DLINUX
INCLUDE_linux-o		:= $(INCLUDE_linux)
VXWORKS_ROOT            := /usr/local/tornado/ppc2/target
INCLUDE_vxworks-ppc	:= $(BASE_INCLUDE) -I$(VXWORKS_ROOT)/h \
                           -I$(VXWORKS_ROOT)/h/net \
                           -mcpu=604 -DCPU=PPC604 -DVXWORKS -D_GNU_TOOL\
                            -nostdinc -fno-for-scope -fno-builtin -fvolatile \
                           -DVXWORKSPPC

INCLUDE:= $(INCLUDE_$(RADPHI_MAKETARGET))

#
# FINCLUDE defines the location of F77 include files
#    *and* any special f77 switches 
#

FINCLUDE_aix 		:=-qextname -I. -I../include -I$(RADPHI_HOME)/include
FINCLUDE_aix-o		:= $(FINCLUDE_rs6000)
FINCLUDE_irix-mips4	:= -static -Wf,-I. -Wf,-I../include -Wf,-I$(RADPHI_HOME)/include $(INCLUDE)
FINCLUDE_irix-mips4-o	:= $(FINCLUDE_irix-mips4)
FINCLUDE_irix-mips3	:= $(FINCLUDE_irix-mips4)
FINCLUDE_irix-mips3-o	:= $(FINCLUDE_irix-mips4)
FINCLUDE_irix-mips2	:= $(FINCLUDE_irix-mips4)
FINCLUDE_irix-mips2-o	:= $(FINCLUDE_irix-mips4)
FINCLUDE_hpux		:= +ppu -C -K -I. -I../include -I$(RADPHI_HOME)/include
FINCLUDE_hpux-o		:= $(FINCLUDE_hpux)
FINCLUDE_solaris	:= -e -I. -I../include -I$(RADPHI_HOME)/include
FINCLUDE_solaris-o	:= $(FINCLUDE_solaris)
FINCLUDE_dec		:= -fpe1 $(FINCLUDE_irix-mips4)
FINCLUDE_dec-o		:= $(FINCLUDE_dec)
FINCLUDE_linux5		:= -fno-second-underscore -I. -I../include -I$(RADPHI_HOME)/include -DLINUX
FINCLUDE_linux5-o	:= $(FINCLUDE_linux5)
FINCLUDE_linux		:= -fno-second-underscore -I. -I../include -I$(RADPHI_HOME)/include -DLINUX
FINCLUDE_linux-o	:= $(FINCLUDE_linux)

FINCLUDE:= $(FINCLUDE_$(RADPHI_MAKETARGET))

#
# Here: COPT- defines the flags used to build the %.oo objects (both C and f77)
#       CDBG- defines the flags used to build the %.o objects,
#			  note that for the optimised targets (such as irix5o)
#			  CDBG is explicitely set equal to COPT
#

COPT_default		:= -O
COPT_irix-mips4		:= -O3
COPT_irix-mips4-o	:= $(COPT_irix-mips4)
COPT_irix-mips3		:= $(COPT_irix-mips4)
COPT_irix-mips3-o	:= $(COPT_irix-mips4)
COPT_irix-mips2		:= -O2
COPT_irix-mips2-o	:= $(COPT_irix-mips2)
COPT_aix		:= -O3 -qstrict
COPT_aix-o		:= $(COPT_aix)
COPT_solaris		:= -O
COPT_solaris-o		:= $(COPT_solaris)
COPT_dec=		:= -g3 -O2
COPT_dec-o		:= $(COPT_dec)
COPT_hpux		:= -O
COPT_hpux-o		:= $(COPT_hpux)
COPT_linux5		:= -O
COPT_linux5-o		:= $(COPT_linux5)
COPT_linux		:= -O3
#ifdef IU
COPT_linux              := -O3
#endif
COPT_linux-o		:= $(COPT_linux)
COPT_vxworks-ppc	:= -g
COPT_vxworks-ppc-o	:= -O

CDBG_default		:= -g
CDBG_irix-mips4-o	:= $(COPT_irix-mips4-o)
CDBG_irix-mips3-o	:= $(COPT_irix-mips4-o)
CDBG_irix-mips2-o	:= $(COPT_irix-mips2-o)
CDBG_aix-o		:= $(COPT_aix-o)
CDBG_solaris-o		:= $(COPT_solaris-o)
CDBG_hpux-o		:= $(COPT_hpux-o)
CDBG_dec-o		:= $(COPT_dec-o)
CDBG_linux5-o		:= $(COPT_linux5-o)
CDBG_linux-o		:= $(COPT_linux-o)

COPT:=$(COPT_$(RADPHI_MAKETARGET))
CDBG:=$(CDBG_$(RADPHI_MAKETARGET))

COPT := -O1

ifndef CDBG
CDBG:=$(CDBG_default)
endif
ifndef COPT
COPT:=$(COPT_default)
endif
ifdef Optimizer
CDBG:=$(COPT)
endif

#
# Various architecture-dependant flags, common between C, f77
#

ARCHFLAGS_irix-mips4		:= -n32 -mips4
ARCHFLAGS_irix-mips4-o		:= -n32 -mips4 -WK,-inline -OPT:IEEE_arithmetic=3:roundoff=3 
ARCHFLAGS_irix-mips3		:=  -n32 -mips3
ARCHFLAGS_irix-mips3-o		:= -n32 -mips3 -WK,-inline -OPT:IEEE_arithmetic=3:roundoff=3 
ARCHFLAGS_irix-mips2		:= -32 -mips2
ARCHFLAGS_irix-mips2-o		:= -32 -mips2

ARCHFLAGS:= $(ARCHFLAGS_$(RADPHI_MAKETARGET))

#
# Various target-dependant C compiler flags
#

CFLAGS_irix-mips4	:= -fullwarn -woff 1174,1209,1233
CFLAGS_irix-mips4-o     := $(CFLAGS_irix-mips4)
CFLAGS_irix-mips3	:= $(CFLAGS_irix-mips4)
CFLAGS_irix-mips3-o	:= $(CFLAGS_irix-mips4)
CFLAGS_irix-mips2	:= -fullwarn
CFLAGS_irix-mips2-o	:= -fullwarn
CFLAGS_dec		:= -ieee_with_inexact
CFLAGS_dec-o		:= $(CFLAGS_dec)
CFLAGS_hpux		:= -Ae
CFLAGS_hpux-o		:= $(CFLAGS_hpux)
CFLAGS_linux		:= -D_FILE_OFFSET_BITS=64
CFLAGS_linux-o		:= $(CFLAGS_linux)

CFLAGS:= $(CFLAGS_$(RADPHI_MAKETARGET))

DLLFLAG_solaris		:= -G
DLLFLAG_linux5		:= -shared
DLLFLAG_linux		:= -shared

DLLFLAG:= $(DLLFLAG_$(RADPHI_MAKETARGET))

#
# C++ stuff
#

CCFLAGS_irix-mips4	:= -ansi
CCFLAGS_irix-mips4-o	:= $(CCFLAGS_irix-mips4)
CCFLAGS_irix-mips3	:= $(CCFLAGS_irix-mips4)
CCFLAGS_irix-mips3-o	:= $(CCFLAGS_irix-mips4)
CCFLAGS_irix-mips2	:= $(CCFLAGS_irix-mips4)
CCFLAGS_irix-mips2-o	:= $(CCFLAGS_irix-mips4)

CCFLAGS:= $(CCFLAGS_$(RADPHI_MAKETARGET))

#
# The irix6 C++ preprocessor had to be kludged (???)
#


CCPPFLAGS = $(CCPPFLAGS_$(RADPHI_MAKETARGET))

#
# Location of the Fortran run-time libraries. 
# Do NOT put on the dependencies list
#

FLIB_aix		:= -lxlf90 -lxlf -lm
FLIB_aix-o		:= $(FLIB_rs6000)
FLIB_irix-mips4		:= -lftn -lm
FLIB_irix-mips4-o	:= $(FLIB_irix-mips4)
FLIB_irix-mips3		:= $(FLIB_irix-mips4)
FLIB_irix-mips3-o	:= $(FLIB_irix-mips4)
FLIB_irix-mips2		:= $(FLIB_irix-mips4)
FLIB_irix-mips2-o	:= $(FLIB_irix-mips4)
FLIB_solaris		:= -lM77 -lF77 -lnsl -ldl -lsocket -lresolv -lsunmath -lm
FLIB_solaris-o		:= $(FLIB_solaris)
FLIB_hpux		:= -lf /usr/lib/libisamstub.sl
FLIB_hpux-o		:= $(FLIB_hpux)
FLIB_dec		:= -lfor -lFutil -lUfor
FLIB_dec-o:=		:= $(FLIB_dec)
FLIB_linux5		:= -lgfortran -lnsl -ldl -lresolv -lm
FLIB_linux5-o		:= $(FLIB_linux5)
FLIB_linux		:= -lgfortran -lnsl -ldl -lresolv -lm
FLIB_linux-o		:= $(FLIB_linux)

FLIB:= $(FLIB_$(RADPHI_MAKETARGET))

#
# Additional stuff put directly on the link command line
#	(not on the dependancies list)
#

SYSLIB:= $(SYSLIB_$(RADPHI_MAKETARGET))

#
# Location of the CERN library is given by the CERNLIBDIR
# environment variable. If it is not set, the default
# is in $RADPHI_HOME/lib.$RADPHI_MAKETARGET
#

CERNLIB:= $(CERN_ROOT)/bin/cernlib


#
# Additional C libraries
#

CLIB:= $(CLIB_$(RADPHI_MAKETARGET))/$(TARGET_ARCH)


#
# Networking libraries (for Solaris)
#
NETLIB_solaris := -lsocket -lnsl
NETLIB_solaris-o := $(NETLIB_solaris)
NETLIB_linux5 := -lnsl
NETLIB_linux5-o := $(NETLIB_linux5)
NETLIB_linux := -lnsl
NETLIB_linux-o := $(NETLIB_linux)

NETLIB := $(NETLIB_$(RADPHI_MAKETARGET))

#
# Flags to use when compiling C code that references F77 common blocks
#


COMMON_irix-mips4	:= -common
COMMON_irix-mips4-o	:= $(COMMON_irix-mips4)
COMMON_irix-mips3	:= $(COMMON_irix-mips4)
COMMON_irix-mips3-o	:= $(COMMON_irix-mips4)
COMMON_irix-mips2	:= $(COMMON_irix-mips4)
COMMON_irix-mips2-o	:= $(COMMON_irix-mips4)

COMMON:=$(COMMON_$(RADPHI_MAKETARGET))


#
# ACPP: the C preprocessor to create dependacies list
#

ACPP_default:=		cc -M
ACPP_hpux:=		acpp -M
ACPP_hpux-o:=		acpp -M
ACPP_decunix:=		cc -M
ACPP_decunixo:=		cc -M
ACPP_solaris:=		/usr/ccs/lib/cpp -M
ACPP_solaris-o:=	/usr/ccs/lib/cpp -M
ACPP_irix5:=		cc -M
ACPP_irix5o:=		cc -M
ACPP_irix5mips1o:=	cc -M
ACPP_irix5mips2o:=	cc -M
ACPP_irix6n32mips3o:=	cc -M
ACPP_linux5:=		gcc -M
ACPP_linux5o:=		gcc -M
ACPP_linux:=		gcc -M
ACPP_linuxo:=		gcc -M
ACPP_aix:=		cpp -M
ACPP_aix-o:=		cpp -M
ACPP_vxworks-ppc:=	ccppc -E $(COPT)

ACPP:=$(ACPP_$(RADPHI_MAKETARGET))

ifndef ACPP
ACPP=$(ACPP_default)
endif

CC := gcc
ifndef CC
CC := cc
endif

#F77      := f77
F77      := gfortran

#
# LIBDIR: all the library files are created/updated there
#		(instead of ~radphi/lib)
#

ifndef LIBDIR
ifeq ($(GLIB),0)
LIBDIR := ../lib.$(RADPHI_MAKETARGET)/$(TARGET_ARCH)
else
LIBDIR := $(RADPHI_HOME)/lib.$(RADPHI_MAKETARGET)/$(TARGET_ARCH)
endif
endif

#
# VPATH: library search path
# This should NOT include system paths 
#

ifeq ($(GLIB),0)
VPATH := ../lib.$(RADPHI_MAKETARGET)/$(TARGET_ARCH)/:$(RADPHI_HOME)/lib.$(RADPHI_MAKETARGET)/$(TARGET_ARCH)
else
VPATH := $(RADPHI_HOME)/lib.$(RADPHI_MAKETARGET)/$(TARGET_ARCH)
endif

VPATH := $(CERNLIBDIR):$(VPATH)

#
# AR: the 'ar' command- to set group-writable permissions
#			on newly created libraries
#

ifeq ($(TARGET),solaris)
AR= umask 2; /usr/xpg4/bin/ar
endif

ifndef AR
AR= umask 2; ar
endif

ifndef RADPHI_PERL
PERL := perl
else
PERL := $(RADPHI_PERL)
endif

.PHONY : all

%: %.o
	-mv -f $@ $@.sav
	cc -o $@ $^

%: %.c
	$(MAKE) $*.o
	$(MAKE) $*

#
#
# These are the main rules to build object files-
#
# %.o --- build "normal" object files (usually with -g)
# %.oo -- build "optimised" object files (always with -O)
#
#


%.oc: %.c
	$(ACPP) $< $(INCLUDE) > $*.oc.depend
	$(CC) -c $< $(CDBG) $(COMMON) $(ARCHFLAGS) $(CFLAGS) $(INCLUDE)
	mv $*.o $*.oc

%.oo: %.c
	$(ACPP) $< $(INCLUDE) > $*.oo.depend
	$(CC) -c $< $(COPT) $(ARCHFLAGS) $(CFLAGS) $(INCLUDE)
	mv $*.o $*.oo

%.oo: %.f
	$(FC) -c $< $(COPT) $(ARCHFLAGS) $(FFLAGS) $(FINCLUDE)
	mv $*.o $*.oo

%.o: %.f
	$(FC) -c $< $(CDBG) $(ARCHFLAGS) $(FFLAGS) $(FINCLUDE)

%.o: %.c
	$(ACPP) $< $(INCLUDE) > $@.depend
	$(CC) -c $< $(CDBG) $(ARCHFLAGS) $(CFLAGS) $(INCLUDE)

#
#
#    Here we define the C++ stuff- which compiler we use
#    and which libraries we link with. This may or may
#    not be working right now
# 

CXX_irix-mips4 		:= CC
CXX_irix-mips4-o 	:= $(CXX_irix-mips4)
CXX_irix-mips3  	:= $(CXX_irix-mips4)
CXX_irix-mips3-o 	:= $(CXX_irix-mips4)
CXX_irix-mips2	 	:= $(CXX_irix-mips4)
CXX_irix-mips2-o 	:= $(CXX_irix-mips4)
CXX_solaris		:= g++
CXX_solaris-o		:= $(CXX_solaris)
CXX_aix			:= g++
CXX_aix-o		:= $(CXX_solaris)
CXX_hpux		:= g++
CXX_hpux-o		:= $(CXX_solaris)
CXX_linux5		:= g++
CXX_linux5-o		:= $(CXX_linux)
CXX_linux		:= g++
CXX_linux-o		:= $(CXX_linux)

CXX= $(CXX_$(RADPHI_MAKETARGET))

CCOPT= $(COPT)
CCDBG= $(CDBG)
ifeq ($(RADPHI_MAKETARGET),hpux)
CCDBG= -O
endif
ifeq ($(RADPHI_MAKETARGET),hpux-o)
CCDBG= -O
endif
CCARCHFLAGS= $(ARCHFLAGS)

ifeq ($(GCC),1)
CXX= g++
CCOPT= -O
CCDBG= $(CDBG)
ifeq ($(RADPHI_MAKETARGET),hpux)
CCDBG= -O
endif
ifeq ($(RADPHI_MAKETARGET),hpux-o)
CCDBG= -O
endif
CCARCHFLAGS= $(GARCHFLAGS)
CCFLAGS:= -Wall -I$(RADPHI_HOME)/include/CC/ISO
endif

%.o: %.C
#	$(ACPP) $< $(INCLUDE) > $@.depend
	$(CXX) $(CCARCHFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCDBG) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)

%.oo: %.C
	$(CXX) $(CCPPFLAGS) -c  -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCOPT) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)
	mv $*.o $*.oo

%.o: %.c++
	$(CXX) $(CCPPFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCDBG) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)

%.oo: %.c++
	$(CXX) $(CCPPFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCOPT) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)
	mv $*.o $*.oo

%.o: %.cxx
	$(CXX) $(CCPPFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCDBG) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)

%.oo: %.cxx
	$(CXX) $(CCPPFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCOPT) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)
	mv $*.o $*.oo

%.o: %.cc
	$(CXX) $(CCPPFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCDBG) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)

%.oo: %.cc
	$(CXX) $(CCPPFLAGS) -c -M $<  $(INCLUDE) > $@.depend
	$(CXX) -c $< $(CCOPT) $(CCARCHFLAGS) $(CCFLAGS) $(INCLUDE)
	mv $*.o $*.oo


#
# end of C++ stuff
# 

#
# Helper rules to install stuff in ../bin.xxx
#

INSTALL := $(PERL) $(RADPHI_HOME)/scripts/install.perl

%.install: %
	@$(INSTALL) -s $^ ../bin.$(RADPHI_MAKETARGET)/$(TARGET_ARCH)

%.install-save: %
	@$(INSTALL) -s $^ ../bin.$(RADPHI_MAKETARGET)/$(TARGET_ARCH)


#
# Make everything precious- 
# to avoid removal of source files after targets are built.
#

.PRECIOUS: %

Debug:
	@echo RADPHI_HOME\\t        $(RADPHI_HOME)
	@echo LIBDIR\\t\\t          $(LIBDIR)
	@echo MAKETARGET\\t         $(RADPHI_MAKETARGET)
	@echo VPATH\\t\\t           $(VPATH)
	@echo ARCHFLAGS\\t          $(ARCHFLAGS)
	@echo CC\\t\\t              $(CC)
	@echo CFLAGS\\t\\t          $(CFLAGS)
	@echo CDBG\\t\\t            $(CDBG)
	@echo CCFLAGS\\t\\t         $(CCFLAGS)
	@echo INCLUDE\\t\\t         $(INCLUDE)


clean:
	rm -f *.o *.oo *.depend

#
# end file
#


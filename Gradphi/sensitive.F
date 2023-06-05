c February 13, 2003
c   Added hits declaration for LGD photocathodes in the case where
c   we want to do explicit Cerenkov generation and following, and
c   count the photoelectrons at the phototube [rtj].
c - Richard Jones
c
c--------------------------------------------------------------------
c Defines sensitive volumes
c--------------------------------------------------------------------
c Author: Andriy Kurylov
c
      subroutine sensitive
      implicit none
#define CERNLIB_TYPE 1
#include "geant321/gcsets.inc"      
      character*4 CHNMSV(7), CHNAMH(12)
      integer NBITSV(7), NBITSH(12)
      real ORIG(12), FACT(12)
c some parameters for sensitive volume definitions
      character*4 CHSET, CHDET
      integer NV, NH, IDTYP
c -----------------------------------------------------------------------
c      
c hits arrays for LGD Cerenkov model response
      CHSET='CALO'
      CHDET='BLOK'
      NV=7
      NH=6
      IDTYP=0
      CHNMSV(1)='WAL1'
      CHNMSV(2)='GBX1'
      CHNMSV(3)='LGB1'
      CHNMSV(4)='WAL2'
      CHNMSV(5)='GBX2'
      CHNMSV(6)='LGB2'
      CHNMSV(7)='BLKC'
      NBITSV(1)=8
      NBITSV(2)=8
      NBITSV(3)=8
      NBITSV(4)=8
      NBITSV(5)=8
      NBITSV(6)=8
      NBITSV(7)=8
      CHNAMH(1)='X   '
      CHNAMH(2)='Y   '
      CHNAMH(3)='Z   '
      CHNAMH(4)='Ein '
      CHNAMH(5)='Eout'   
      CHNAMH(6)='NoPE'
      NBITSH(1)=32
      NBITSH(2)=32
      NBITSH(3)=32
      NBITSH(4)=32
      NBITSH(5)=32
      NBITSH(6)=32
      ORIG(1)=1000.
      ORIG(2)=1000.
      ORIG(3)=1000.
      ORIG(4)=0.
      ORIG(5)=0.
      ORIG(6)=0.
      FACT(1)=1.0e+4
      FACT(2)=1.0e+4
      FACT(3)=1.0e+4
      FACT(4)=1.0e+6
      FACT(5)=1.0e+6
      FACT(6)=1.0e+3
c declaring BLOK as a sensitive volume
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)

c hits arrays for LGD Cerenkov track-following response
      CHSET='CALO'
      CHDET='PMTC'
      NV=1
      NH=6
      IDTYP=1
      CHNMSV(1)='PMTL'
c the rest of the parameters have not changed since the definition of BLOK 
c now PMTC is a sensitive volume
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)
c-----------------------------------------------------------------------------     
c now specify the rest of the arrays for other volumes
c
      NBITSV(1)=8
      CHNAMH(1)='X   '
      CHNAMH(2)='Y   '
      CHNAMH(3)='Z   '
      CHNAMH(4)='Etot'
      CHNAMH(5)='Elos'   
      CHNAMH(6)='ToF '
      CHNAMH(7)='TDC '
      NBITSH(1)=32
      NBITSH(2)=32
      NBITSH(3)=32
      NBITSH(4)=32
      NBITSH(5)=32
      NBITSH(6)=32
      NBITSH(7)=32 
      ORIG(1)=1000.
      ORIG(2)=1000.
      ORIG(3)=1000.
      ORIG(4)=0.
      ORIG(5)=0.
      ORIG(6)=100e-9
      ORIG(7)=100e-9
      FACT(1)=1.0e+4
      FACT(2)=1.0e+4
      FACT(3)=1.0e+4
      FACT(4)=1.0e+6
      FACT(5)=1.0e+6
      FACT(6)=1.0e+12
      FACT(7)=1.0e+12      
c next goes RPD. Here DSEN and THSN are the volumes
      CHSET='RPD '
      CHDET='DSEN'
      NV=1
      NH=7
      IDTYP=10
      CHNMSV(1)='CEL1'
      NBITSV(1)=32
c the rest of the parameters have not changed since the definition of BLOK
c now DSEN is a sensitive volume
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)
c
      CHSET='RPD '
      CHDET='THSN'
      NV=2
      NH=7
      IDTYP=11
      CHNMSV(1)='CEL1'
      CHNMSV(2)='THSN'
      NBITSV(1)=32
      NBITSV(2)=32
c the rest of the parameters have not changed since the definition of BLOK 
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)
c------------------------------------------------------------------------------
c next goes CPV
      CHSET='CPV '
      CHDET='CPSA'
      NV=1
      NH=7
      IDTYP=20
      CHNMSV(1)='CPVA'
      NBITSV(1)=32
c the rest of the parameters have not changed since the definition of BLOK 
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)
      CHDET='CPSB'
      IDTYP=21
      CHNMSV(1)='CPVB'
      NBITSV(1)=32
c the rest of the parameters have not changed since the definition of BLOK 
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)
      CHDET='CPSC'
      IDTYP=22
      CHNMSV(1)='CPVC'
      NBITSV(1)=32
c the rest of the parameters have not changed since the definition of BLOK 
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)
      CHDET='CPSD'
      IDTYP=23
      CHNMSV(1)='CPVD'
      NBITSV(1)=32
c the rest of the parameters have not changed since the definition of BLOK 
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)
      CHDET='CPSE'
      IDTYP=24
      CHNMSV(1)='CPVE'
      NBITSV(1)=32
c the rest of the parameters have not changed since the definition of BLOK 
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)
c---------------------------------------------------------------------------
c UPV
      CHSET='UPV '
      CHDET='UPVS'
      NV=1
      NH=7
      IDTYP=30
      CHNMSV(1)='UPVA'
      NBITSV(1)=32
c the rest of the parameters have not changed since the definition of BLOK
c now DSEN is a sensitive volume
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)
c---------------------------------------------------------------------------
c  BSD
c must define BSR1, BSDR2, and BSR3 to be in this set
      CHSET='BSD '
      CHDET='BSR1'
      NV=1
      NH=7
      IDTYP=40
      CHNMSV(1)='BSR1'
      NBITSV(1)=8
      NBITSH(6)=32
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)
c now for BSR2
      CHDET='BSR2'
      CHNMSV(1)='BSR2'
      IDTYP=41
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)
c now for BSR3
      CHDET='BSR3'
      CHNMSV(1)='BSR3'
      IDTYP=42
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)
c---------------------------------------------------------------------------
c  This section by Tom Bogue 6/5/99
c  BGV
c  see BGVsave.c for a description of the hit params.
      CHSET='BGV '
      CHDET='BGVS'
      NV=1
      NH=12
      CHNAMH(1)='Xin '
      CHNAMH(2)='Yin '
      CHNAMH(3)='Zin '
      CHNAMH(4)='tmin'
      CHNAMH(5)='EtUp'
      CHNAMH(6)='EtDn'
      CHNAMH(7)='ECHA'
      CHNAMH(8)='ENEU'
      CHNAMH(9)='EIn '
      CHNAMH(10)='EOut'
      CHNAMH(11)='Eup '
      CHNAMH(12)='Edn '
      ORIG(4)=100e-9
      ORIG(5)=0
      ORIG(6)=0
      ORIG(7)=0
      ORIG(8)=0
      ORIG(9)=0
      ORIG(10)=0
      ORIG(11)=0
      ORIG(12)=0
      FACT(4)=1.0e+12
      FACT(5)=1.0e+6
      FACT(6)=1.0e+6
      FACT(7)=1.0e+6
      FACT(8)=1.0e+6
      FACT(9)=1.0e+6
      FACT(10)=1.0e+6
      FACT(11)=1.0e+6
      FACT(12)=1.0e+6
c record energy in keV, time in ps and distance in um, just like everywhere else
      NBITSH(8)=32
      NBITSH(9)=32
      NBITSH(10)=32
      NBITSH(11)=32
      NBITSH(12)=32
      IDTYP=50
      CHNMSV(1)='BGVS'
      NBITSV(1)=8
      call GSDET(CHSET,CHDET,NV,CHNMSV,
     +           NBITSV,IDTYP,100,100,ISET,IDET)
      call GSDETH(CHSET,CHDET,NH,CHNAMH,NBITSH,ORIG,FACT)
      end

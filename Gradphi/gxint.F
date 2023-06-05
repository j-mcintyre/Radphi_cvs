*
* $Id: gxint.F,v 1.4 2005/09/24 18:10:09 radphi Exp $
*
* $Log: gxint.F,v $
* Revision 1.4  2005/09/24 18:10:09  radphi
* BGVsave.c, CPVsave.c
*    - remove unused parameters tres,convert from argument list [rtj]
* LGDsave.c
*    - remove arguments lgdnx,lgdny and replace with globals LGD_ROWS and
*      LGD_COLUMNS [rtj]
*    - remove argument convert and replace with global LGD_GEV_PER_PE
* Makefile
*    - include fixes for gxcs.F
* RPDsave.c, UPVsave.c
*    - remove unused parameters tres from argument lists [rtj]
* angles.F
*    - introduce implicit none [rtj]
* background.F
*    - introduce implicit none [rtj]
*    - express bggate constants in s [rtj]
* beampipe.F, bgv.F, bsd.F, clasmgnt.F, clastgt.F, cerenkov.F, collim.F, cpv.F
*    - introduce implicit none [rtj]
* bgvhits.F
*    - introduce implicit none [rtj]
*    - improve the treatment of multi-hit pulse timing in the BGV [rtj]
* bsdhits.F
*    - introduce implicit none [rtj]
*    - improve the treatment of pulse-height attenuation in the BSD [rtj]
* control.in
*    - comment out the addition of EM background to events [rtj]
* cpvhits.F, dnut.F, detector.F, gxvertx.F, gufld.F, gustep.F, gutrev.F, gxint.F
*    - introduce implicit none [rtj]
* gukine.F
*    - introduce implicit none [rtj]
*    - remove arguments Nx,Ny,convert,tres from save_next_event() [rtj]
* geometry.inc
*    - adopt a uniform naming convention for geometric constants [rtj]
*    - increase the material in the barrel support structure [rtj]
* hitconstants.inc
*    - adopt a uniform naming convention for response model parameters [rtj]
*    - moved digitization parameters to c header file digitization.h [rtj]
* hercover.F, lgd.F, lgdframe.F, lgdhits.F, lgpmt.F, major.F, matdef.F, meddef.F
*    - introduce implicit none [rtj]
* outputEvents.c
*    - removed unused arguments from argument lists [rtj]
* place*.F, platform.F, radiator.F, ratestudy.F, rpd.F, rpdhits.F
*    - introduce implicit none [rtj]
* savehits.F, sensitive.F, shield.F, sweepmgnt.F
*    - introduce implicit none [rtj]
* timeAnalysis*.c
*    - exchanged argument tres for global TDC_PULSE_RESOLUTION [rtj]
*    - introduced discriminator thresholds XXX_TDC_THRESHOLD to get more
*      realistic counting rates in detectors [rtj]
* trgt.F, uginit.F, uglast.F, uhinit.F, upv.F, upvhits.F
*    - introduce implicit none [rtj]
*
* Revision 1.3  2001/06/25 16:47:23  radphi
* Changes made by jonesrt@hector
* checking in a new version of Makefile that works with rules.make -rtj
*
* Revision 1.2  1997/01/07 10:25:42  cernlib
* Remove #ifdef CERNLIB_MAIN; this shall be done via Imakefile.
*
* Revision 1.1.1.1  1995/10/24 10:21:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      PROGRAM GXINT
      implicit none
*
*     GEANT main program. To link with the MOTIF user interface
*     the routine GPAWPP(NWGEAN,NWPAW) should be called, whereas
*     the routine GPAW(NWGEAN,NWPAW) gives access to the basic
*     graphics version.
*
      INTEGER NWGEAN,NWPAW
#if !defined(CERNLIB_IBM)
      PARAMETER (NWGEAN=3000000,NWPAW=1000000)
#endif
#if defined(CERNLIB_IBM)
      PARAMETER (NWGEAN=1000000,NWPAW=500000)
#endif
      REAL GEANT,PAW
      COMMON/GCBANK/GEANT(NWGEAN)
      COMMON/PAWC/PAW(NWPAW)
#if defined(CERNLIB_HPUX)
      ON REAL UNDERFLOW IGNORE
#endif
*
#if defined(CERNLIB_IBM)
      CALL INITC
      CALL ERRSET(151,999,-1)
#endif
#if (defined(CERNLIB_MOTIF))&&(!defined(CERNLIB_IBM))
      CALL GPAWPP(NWGEAN,NWPAW)
#endif
#if !defined(CERNLIB_MOTIF)||defined(CERNLIB_IBM)
      CALL GPAW(NWGEAN,NWPAW)
#endif
*
      END
      SUBROUTINE QNEXT
      END
#if !defined(CERNLIB_CZ)
      SUBROUTINE CZOPEN
      END
      SUBROUTINE CZTCP
      END
      SUBROUTINE CZCLOS
      END
      SUBROUTINE CZPUTA
      END
#endif
#if defined(CERNLIB_IBM)
      FUNCTION IOSCLR()
      IOSCLR=0
      END
#endif

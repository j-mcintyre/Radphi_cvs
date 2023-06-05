      subroutine gsvertx(vert,ntbeam,nttarg,tofv,ubuf,nubuf,nvtx)
      implicit none

C.    * This is a stub interface to the real GSVERT geant routine
C.    * that includes the TOFG in the argument list.  From c it is
C.    * difficult to load a fortran common block element like TOFG
C.    * before the call to GSVERT so this routine stands in between.

#define CERNLIB_TYPE 1
#include "geant321/gctrak.inc"

      real vert(3)
      integer ntbeam,nttarg
      real tofv,ubuf
      integer nubuf,nvtx

      TOFG=tofv
      call GSVERT(vert,ntbeam,nttarg,ubuf,nubuf,nvtx)
      end

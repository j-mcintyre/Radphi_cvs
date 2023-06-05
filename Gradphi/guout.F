c November 18, 2002
c   Merged the standard subroutine with special code for measuring
c   charged particle rates in the CPV and BSD layers.  The special
c   code is selected with the switch SAVE_RATES_NTUPLE.  To see the
c   main code see source file ratestudy.F.
c - Richard Jones
c
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GUOUT
      implicit none
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine called at the end of each event             *
C.    *                                                                *
C.    *    ==>Called by : GTRIG                                        *
C.    *                                                                *
C.    ******************************************************************
C.
#define CERNLIB_TYPE 1
#include "geant321/pilot.h"
#include "geant321/gckine.inc"
#include "geant321/gcsets.inc"
#include "geant321/gcflag.inc"

#     include "geometry.inc"
#     include "hitconstants.inc"

      external save_next_event

      real mantissa
#if defined SAVE_RATES_NTUPLE
# include "ratestudy.inc"
      character*4 chset,chdet
      integer nvdim,nhdim,nhmax
      integer itrs,numvs(1)
      integer itra(1)
      integer numv(1,1)
      integer nhitu,nhitv,nhitb
      real hits(10)
      logical sawsomething
      logical hexist
      external hexist
#endif
C
c if this event is a power of 10, print the run number
      mantissa=IDEVT
      do while (mantissa.gt.10)
        mantissa=mantissa/10
      enddo
      if (int(mantissa).eq.mantissa) 
     +  print *,'event #',IDEVT,' simulated'
      call save_next_event(IDRUN)
C       
C.    ------------------------------------------------------------------
C.
#if defined SAVE_RATES_NTUPLE
      if (.not.hexist(ntID)) then
         call hcdir('//RZfile',' ')
         call hbnt(ntID,'hits on scintillator planes',' ')
         call hbname(ntID,'hits',ntvars,ntnames)
      endif

      do i=1,50
        ntvars(i)=0
      enddo
      sawsomething=.false.

      chset='BRAD'
      chdet='UPVS'
      nvdim=1
      nhdim=10
      nhmax=1
      itrs=0
      numvs(1)=0
      call gfhits(chset,chdet,nvdim,nhdim,nhmax,itrs,numvs,
     +            itra,numv,hits,nhitu)
      if (nhitu.gt.0) then
        sawsomething=.true.
        Xupv=hits(1)
        Yupv=hits(2)
        Zupv=hits(3)
        Tupv=hits(4)
        Nnupv=hits(5)
        Ncupv=hits(6)
        Enupv=hits(7)
        Ecupv=hits(8)
        Xnupv=hits(9)
        Xcupv=hits(10)
      endif

      Nncpv=0
      Nccpv=0
      Encpv=0
      Eccpv=0
      Xncpv=0
      Xccpv=0
      chset='BRAD'
      chdet='CPSA'
      nvdim=1
      nhdim=10
      nhmax=1
      itrs=0
      numvs(1)=0
      call gfhits(chset,chdet,nvdim,nhdim,nhmax,itrs,numvs,
     +            itra,numv,hits,nhitv)
      if (nhitv.gt.0) then
        sawsomething=.true.
        Xcpv=hits(1)
        Ycpv=hits(2)
        Zcpv=hits(3)
        Tcpv=hits(4)
        Nncpv=Nncpv+hits(5)
        Nccpv=Nccpv+hits(6)
        Encpv=Encpv+hits(7)
        Eccpv=Eccpv+hits(8)
        Xncpv=Xncpv+hits(9)
        Xccpv=Xccpv+hits(10)
      endif
      chdet='CPSB'
      call gfhits(chset,chdet,nvdim,nhdim,nhmax,itrs,numvs,
     +            itra,numv,hits,nhitv)
      if (nhitv.gt.0) then
        sawsomething=.true.
        Xcpv=hits(1)
        Ycpv=hits(2)
        Zcpv=hits(3)
        Tcpv=hits(4)
        Nncpv=Nncpv+hits(5)
        Nccpv=Nccpv+hits(6)
        Encpv=Encpv+hits(7)
        Eccpv=Eccpv+hits(8)
        Xncpv=Xncpv+hits(9)
        Xccpv=Xccpv+hits(10)
      endif
      chdet='CPSC'
      call gfhits(chset,chdet,nvdim,nhdim,nhmax,itrs,numvs,
     +            itra,numv,hits,nhitv)
      if (nhitv.gt.0) then
        sawsomething=.true.
        Xcpv=hits(1)
        Ycpv=hits(2)
        Zcpv=hits(3)
        Tcpv=hits(4)
        Nncpv=Nncpv+hits(5)
        Nccpv=Nccpv+hits(6)
        Encpv=Encpv+hits(7)
        Eccpv=Eccpv+hits(8)
        Xncpv=Xncpv+hits(9)
        Xccpv=Xccpv+hits(10)
      endif
      chdet='CPSD'
      call gfhits(chset,chdet,nvdim,nhdim,nhmax,itrs,numvs,
     +            itra,numv,hits,nhitv)
      if (nhitv.gt.0) then
        sawsomething=.true.
        Xcpv=hits(1)
        Ycpv=hits(2)
        Zcpv=hits(3)
        Tcpv=hits(4)
        Nncpv=Nncpv+hits(5)
        Nccpv=Nccpv+hits(6)
        Encpv=Encpv+hits(7)
        Eccpv=Eccpv+hits(8)
        Xncpv=Xncpv+hits(9)
        Xccpv=Xccpv+hits(10)
      endif
      chdet='CPSE'
      call gfhits(chset,chdet,nvdim,nhdim,nhmax,itrs,numvs,
     +            itra,numv,hits,nhitv)
      if (nhitv.gt.0) then
        sawsomething=.true.
        Xcpv=hits(1)
        Ycpv=hits(2)
        Zcpv=hits(3)
        Tcpv=hits(4)
        Nncpv=Nncpv+hits(5)
        Nccpv=Nccpv+hits(6)
        Encpv=Encpv+hits(7)
        Eccpv=Eccpv+hits(8)
        Xncpv=Xncpv+hits(9)
        Xccpv=Xccpv+hits(10)
      endif

      chset='BRAD'
      chdet='BSR1'
      nvdim=1
      nhdim=10
      nhmax=1
      itrs=0
      numvs(1)=0
      call gfhits(chset,chdet,nvdim,nhdim,nhmax,itrs,numvs,
     +            itra,numv,hits,nhitb)
      if (nhitb.gt.0) then
        sawsomething=.true.
        Xbsd(1)=hits(1)
        Ybsd(1)=hits(2)
        Zbsd(1)=hits(3)
        Tbsd(1)=hits(4)
        Nnbsd(1)=hits(5)
        Ncbsd(1)=hits(6)
        Enbsd(1)=hits(7)
        Ecbsd(1)=hits(8)
        Xnbsd(1)=hits(9)
        Xcbsd(1)=hits(10)
      endif
      chdet='BSR2'
      call gfhits(chset,chdet,nvdim,nhdim,nhmax,itrs,numvs,
     +            itra,numv,hits,nhitb)
      if (nhitb.gt.0) then
        sawsomething=.true.
        Xbsd(2)=hits(1)
        Ybsd(2)=hits(2)
        Zbsd(2)=hits(3)
        Tbsd(2)=hits(4)
        Nnbsd(2)=hits(5)
        Ncbsd(2)=hits(6)
        Enbsd(2)=hits(7)
        Ecbsd(2)=hits(8)
        Xnbsd(2)=hits(9)
        Xcbsd(2)=hits(10)
      endif
      chdet='BSR3'
      call gfhits(chset,chdet,nvdim,nhdim,nhmax,itrs,numvs,
     +            itra,numv,hits,nhitb)
      if (nhitb.gt.0) then
        sawsomething=.true.
        Xbsd(3)=hits(1)
        Ybsd(3)=hits(2)
        Zbsd(3)=hits(3)
        Tbsd(3)=hits(4)
        Nnbsd(3)=hits(5)
        Ncbsd(3)=hits(6)
        Enbsd(3)=hits(7)
        Ecbsd(3)=hits(8)
        Xnbsd(3)=hits(9)
        Xcbsd(3)=hits(10)
      endif
      
      if (sawsomething) then
        call hfnt(ntID)
      endif
#endif
      END

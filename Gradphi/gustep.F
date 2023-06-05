c
c February 13, 2003
c   Introduced a special value for the FOLLOW_CERENKOV switch to
c   allow for generation of a detailed ntuple of shower results.
c   To use this feature, set FOLLOW_CERENKOV to an integer greater
c   than zero.
c - Richard Jones
c
c November 18, 2002
c   Merged the standard subroutine with special code for measuring
c   the bremsstrahlung beam energy spectrum at the radiator and on
c   the target.  This code is selected by compiling with the switch
c   DIAGNOSE_BREMSSTRAHLUNG turned on.
c - Richard Jones
c
c November 18, 2002
c   Merged the standard subroutine with special code for measuring
c   charged particle rates in the CPV and BSD layers.  The special
c   code is selected with the switch SAVE_RATES_NTUPLE.  To see the
c   main code see source file ratestudy.F.
c - Richard Jones
c
c November 18, 2002
c   Introduced switch FOLLOW_CERENKOV to enable a simulation mode
c   where each Cerenkov photon is produced individually in the LGD
c   and followed using ray optics to the point of detection or
c   absorption.  The switch EARLY_DETECTION can be used to request
c   that the detection take place at the point of generation instead
c   of at the tube.
c - Richard Jones
c------------------------------------------------------------------
*CMZ :          30/08/94  14.16.48  by  S.Ravndal
*-- Author :
      SUBROUTINE GUSTEP
      implicit none
C.
C.    *
C.    *       User routine called at the end of each tracking step
C.    *       INWVOL is different from 0 when the track has reached
C.    *              a volume boundary
C.    *       ISTOP is different from 0 if the track has stopped
C.    *
C.
#define CERNLIB_TYPE 1
#include "geant321/gctmed.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcphys.inc"

      character*4 CNAMES(99)
      equivalence (NAMES(1),CNAMES(1))
      logical hexist
      external hexist

      real disp(3)
      integer slipsteps
      save slipsteps
      integer ip

#if defined SAVE_RATES_NTUPLE
      logical firstseen
      data firstseen/.true./
      save firstseen
#endif

#if defined FOLLOW_CERENKOV
# include "geant321/gcbank.inc"
# include "geant321/gcjloc.inc"
#     include "cerenkov.inc"
      integer istat
      real rnd(MXPHOT)
      real rat,effic
      integer maxstore
      parameter (maxstore=10000)
      real vlast,enu,E0,x0,cost0,store
      common /follower/vlast(6),enu,E0,x0(3),cost0,store(7,maxstore)
      integer ibin,ilo,ihi
      logical firstcall
      data firstcall/.true./
#endif

#if defined DIAGNOSE_BREMSSTRAHLUNG
      common /diagBrems/ptype,x0(3),E0
      integer ptype
      real x0,E0
      save /diagBrems/
      if (.not.hexist(8)) then
        call hbook1(8,'spectrum of beam in radiator',600,0.,6.0,0.)
        call hbook1(9,'spectrum of beam on target',600,0.,6.0,0.)
# if DIAGNOSE_BREMSSTRAHLUNG == 10
        call hbnt(10,'all hits in LGD','//RZfile')
        call hbname(10,'/lgdhits',ptype,'ptype:i')
        call hbname(10,'/lgdhits',x0,'x0(3):r')
        call hbname(10,'/lgdhits',E0,'E0:r')
# endif
      endif
#endif

# if defined SAVE_RATES_NTUPLE
      if (firstseen) then
        call declare_hits
        firstseen=.false.
      endif
# endif

#if defined FOLLOW_CERENKOV && FOLLOW_CERENKOV
      if (firstcall) then
        call hbnt(1,'showers',' ')
        call hbname(1,'/prof',IEVENT,'ievent:i')
        call hbname(1,'/prof',x0,'x0(3):r')
        call hbname(1,'/prof',cost0,'cost0:r')
        call hbname(1,'/prof',enu,'Enu:r')
        call hbname(1,'/prof',E0,'E0:r')
        firstcall=.false.
      endif
#endif

      if ((INWVOL.gt.0).and.(STEP.le.EPSIL)) then
        slipsteps=slipsteps+1
        if (slipsteps.gt.10) then
          call GRNDM(disp,3)
          VECT(1)=VECT(1)+(0.5-disp(1))*PREC
          VECT(2)=VECT(2)+(0.5-disp(2))*PREC
          VECT(3)=VECT(3)+(0.5-disp(3))*PREC
        endif
      else
        slipsteps=0
      endif
C
C Ngkin is > 0 if secondary particle have been generated,
C store them to the particle stak
C
      if (idebug.ne.0) then
          call gdebug
      endif

      do ip=1,NGKINE
        if (GKIN(5,ip).ne.4) then
          call GSKING(ip)
        endif
      enddo

      if (CNAMES(nlevel).eq.'BRAD') then
        if (ipart.eq.3) then
          if (inwvol.eq.1) then
            zintbr=1e-6
          else
            istop=1
          endif
#if defined DIAGNOSE_BREMSSTRAHLUNG
        elseif (inwvol.eq.1) then
          if (ipart.eq.1) then
            if (GETOT.gt.0) then
              call hfill(8,GETOT,0.,1.)
            endif
          endif
#endif
        endif
#if DIAGNOSE_BREMSSTRAHLUNG == 10
      elseif ((inwvol.eq.1.and.nstep.gt.0).and.
     +        (CNAMES(nlevel).eq.'BLOK'.or.
     +         CNAMES(nlevel).eq.'BGVS')) then
        ptype=IPART
        x0(1)=VECT(1)
        x0(2)=VECT(2)
        x0(3)=VECT(3)
        E0=GEKIN
        call hfnt(10)
#endif
#if defined DIAGNOSE_BREMSSTRAHLUNG
      elseif (CNAMES(nlevel).eq.'TRGT') then
        if (inwvol.eq.1) then
          if (ipart.eq.1) then
            if (GETOT.gt.0) then
              call hfill(9,GETOT,0.,1.)
            endif
          endif
        endif
#endif
      endif

#if defined FOLLOW_CERENKOV
C
C---- Also store any generated C photons
      if (CNAMES(NLEVEL).eq.'BLOK') then
        call GRNDM(rnd,NGPHOT)
        do i=1,NGPHOT
          if (XPHOT(7,i).ge.Q(JTCKOV+NPCKOV+1)) then
            rat=1.
            ibin=NPCKOV-1
          elseif (VECT(7).lt.Q(JTCKOV+2)) then
            rat=0.
            ibin=1
          else
            ilo = 1
            ihi = NPCKOV
   10       ibin = (ilo+ihi)/2
            if (Q(JTCKOV+ibin+1).lt.XPHOT(7,i)) then
              ilo = ibin
            else
              ihi = ibin
            endif
            if (ihi-ilo.gt.1) go to 10
            ibin = ilo
            rat = (XPHOT(7,i)-Q(JTCKOV+ibin+1))/
     +            (Q(JTCKOV+ibin+2)-Q(JTCKOV+ibin+1))
          endif
          effic=(1-rat)*effPMT(ibin)+rat*effPMT(ibin+1)
          if (rnd(i).lt.effic) then
            call GSKPHO(i)
# if FOLLOW_CERENKOV
            n = IQ(JSTAK+1)
            if (n.gt.maxstore) stop 'gustep - stack overflow'
            store(1,n) = (VECT(1)+vlast(1))/2
            store(2,n) = (VECT(2)+vlast(2))/2
            store(3,n) = (VECT(3)+vlast(3))/2
            store(4,n) = (VECT(1)-vlast(1))/STEP
            store(5,n) = (VECT(2)-vlast(2))/STEP
            store(6,n) = (VECT(3)-vlast(3))/STEP
            store(7,n) = GEKIN+DESTEP/2
          else if (rnd(i).gt.0.99) then
            enu=XPHOT(7,i)
            cost0=vlast(6)
            x0(1)=(VECT(1)+vlast(1))/2
            x0(2)=(VECT(2)+vlast(2))/2
            x0(3)=-(VECT(3)+vlast(3))/2
            E0=GEKIN+DESTEP/2
            call HFNT(1)
# endif
          endif
        enddo
# if FOLLOW_CERENKOV
        vlast(1)=VECT(1)
        vlast(2)=VECT(2)
        vlast(3)=VECT(3)
        vlast(4)=VECT(4)
        vlast(5)=VECT(5)
        vlast(6)=VECT(6)
# endif
      endif
#endif

#if defined SAVE_RATES_NTUPLE
      if (isvol.ne.0) then
        if (sleng.gt.0) then
          if (inwvol.ne.0) then
            call register_hit
          endif
        endif
      endif
#else
C
C       Saving hits; main subroutine
C
      if (ISVOL.ne.0) then
       call savehits
      endif 
#endif

#if defined FOLLOW_CERENKOV && FOLLOW_CERENKOV
# if defined EARLY_DETECTION
C      
C---- Generate hits for Cerenkov photons from GEANT
      if (IPART.eq.50) then
        x0(1)=store(1,ISTAK)
        x0(2)=store(2,ISTAK)
        x0(3)=store(3,ISTAK)
        cost0=store(6,ISTAK)
        E0=store(7,ISTAK)
        enu=GETOT
        call HFNT(1)
        ISTOP=1
      endif
# else
      if (ISTOP.gt.0) then
        if ((IPART.eq.50).and.(DESTEP.gt.0).and.
     +      (CNAMES(NLEVEL).eq.'PMTC')) then
          x0(1)=store(1,ISTAK)
          x0(2)=store(2,ISTAK)
          x0(3)=store(3,ISTAK)
          cost0=store(6,ISTAK)
          E0=store(7,ISTAK)
          enu=GETOT
          call HFNT(1)
        endif
      endif
# endif
#endif
      END

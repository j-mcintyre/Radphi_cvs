*cmz :          30/08/94  16.21.00  by  S.Ravndal
*-- Author :
      SUBROUTINE GUKINE
      implicit none
C
C Generates Kinematics for primary track
C
#define CERNLIB_TYPE
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcscan.inc"
C
      REAL VERTEX(3),PLAB(3),RNDM(3)
      SAVE VERTEX,PLAB
      DATA VERTEX/3*0./
      DATA PLAB  /3*0./
      REAL PABS,ETOT
      CHARACTER*20 CHNPAR
      INTEGER NUBUF,IADR
      REAL UBUF(999)
      INTEGER IK
      REAL THETA,PHI
      INTEGER NVERT,NT
      INTEGER store_mc_event, store_mc_vertex
      EXTERNAL store_mc_event, store_mc_vertex

      character inspec*9
      parameter (inspec='events.in')
      character infile*120
      logical firstCall,oneTrack
      data firstCall/.true./
      data oneTrack/.false./
      save firstCall,oneTrack
      save infile
      integer iseed1,iseed2

      logical load_next_event
      external load_next_event

      if (firstCall) then
        firstCall=.false.
        call GRNDMQ(iseed1,iseed2,0,'G')
        if ((iseed1.eq.0).and.(iseed2.eq.0)) then
           write(LOUT,980) 'GUKINE error: '
980        format(1x,a14,'initial random number generator seeds are 0!')
           stop 'cannot continue without good random numbers'
        endif
        open(unit=48,file=inspec,status='old',err=99)
        read(48,'(a120)',end=20) infile  ! end of control file => STOP
        if (infile(1:1).eq.' ') then
          oneTrack=.true.
        endif
        goto 90
99        write(LOUT,990) inspec
          write(LOUT,992)
990       format(' GUKINE - open failed on input control file ',a10)
992       format(' Failing over to one-track generation mode.')
          oneTrack=.true.
90      continue
      endif

      if (oneTrack) goto 10

c--------------------------------------------------------------------
c Now follows the code for the case that we are reading event kinematics
c from a file.  The interface to event generation was designed in c so
c here we just get the next input filename from the input control file
c and pass it to a c routine which does the real work of specifying the
c kinematics for the event.  The c routine returns a logical value of
c true each time it successfully loads the next event from the file.
c When it returns false, it is time to advance to the next input file
c (if any) in the input control file list.  When finished with the input
c control file list, the job is done.

      do while (.not.load_next_event(infile))
        read(48,'(a120)',end=20) infile
        if (infile(1:1).eq.' ') then
          goto 10
        endif
      enddo
      call superimposeBg
      return
 
20    IEORUN=1
      IEOTRI=1
      oneTrack=.true.
      return

10    oneTrack=.true.
c This is the default gukine code for the GEANT examples.  Here it is
c called one-track mode, programmed by the data card KINE or from the
c interactive command line by the command KINE, as in:
c      kine <particle type> <energy> <theta> <phi> [ <dE> <dtheta> <dphi> ]
c--------------------------------------------------------------------------
C
C If the particle number is equal to the particle ID the direction
C cosines are distributed randomly, otherwise the angles theta and phi
C given by the KINE data card are taken.  If a range of energies and
C angles is desired then the sizes of the respective ranges are given
C in arguments 4, 5 and 6 of the KINE data card.
C
      VERTEX(1)=VSCAN(1)
      VERTEX(2)=VSCAN(2)
      VERTEX(3)=VSCAN(3)
      IF (IKINE.GT.100)THEN
         IK=IKINE-100
         CALL GRNDM(RNDM,3)
         PABS=PKINE(1)+PKINE(4)*(RNDM(1)-0.5)
         THETA=(PKINE(2)+PKINE(5)*(RNDM(2)-0.5))*DEGRAD
         PHI=(PKINE(3)+PKINE(6)*(RNDM(3)-0.5))*DEGRAD
      ELSE
         IK=IKINE
         CALL GRNDM(RNDM,2)
         THETA=PI*RNDM(1)
         PHI=TWOPI*RNDM(2)
      ENDIF
C
C Calculating the momementum in the LAB-frame
C
      PLAB(1) = PABS*SIN(THETA)*COS(PHI)
      PLAB(2) = PABS*SIN(THETA)*SIN(PHI)
      PLAB(3) = PABS*COS(THETA)
C
C Storing the first vertex and retrieving the actual vertex number NVERT
C in the JVERTX data structure
C
      IADR = 0
      NUBUF = store_mc_event(UBUF)/4
      CALL GSRUNG(NUBUF,UBUF,IADR)
      CALL GFPART(IK,CHNPAR,ITRTYP,AMASS,CHARGE,TLIFE,UBUF,NUBUF)
      ETOT = SQRT(PABS**2 + AMASS**2)
      NUBUF = store_mc_vertex(UBUF,IK,VERTEX,PLAB,ETOT)/4
      CALL GSVERT(VERTEX,0,0,UBUF,NUBUF,NVERT)
C
C Store the particle IK with its given momentum in the particle stack
C and attach the primary vertex NVERT to it
C
      CALL GSKINE(PLAB,IK,NVERT,0,0,NT)
C
C Kinematics debug (controlled by ISWIT(1) )
C
      IF (IDEBUG.EQ.1) THEN
          IF (ISWIT(1).EQ.1) THEN
              CALL GPRINT('VERT',0)
              CALL GPRINT('KINE',0)
          END IF
      END IF
C
      END

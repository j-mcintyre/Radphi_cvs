*
* -R.T. Jones  March 11, 2004
*  Added support for a user-defined weights that determine the partition
*  of the ntuple or chain among the slaves.  Default weights of unity are
*  assigned to each slave at piaf startup.  The overall normalization of
*  the weights is arbitrary.  The command to assign new weights is:
*       paw> piaf/message reload <s> <w>
*  where <s> is a slave index 1..NSLAVE and <w> is an integer weight that
*  represents the relative number of rows to be processed by that slave.
*  To disable, comment out the next line.
#define WEIGHTED_NTUPLE_PARTITIONING 1
*
* -R.T. Jones  July 30, 2003
*  Added support for ntuple chains longer than 2^31 events.
*  To do that we had to change a number of variables and arguments
*  from int to long long int, and hoped not to break anything.  So
*  far we believe that this should work on piaf without any changes
*  to paw itself.  To disable, comment out the following line.
#define LONG_LONG_CHAIN 1
*
* -R.T.Jones - May 30, 2003
*  Provide a mechanism to disable the automatic purging of histograms and
*  vectors before every nt/xxxx command.  This autopurge feature makes the
*  piaf service stateless, and is normally desirable.  However there are 
*  when it would be very handy to be able to save state information in
*  histograms and vectors that have a lifetime beyond the scope of a single
*  command.  In this state manual purging is still possible through an
*  explicit piaf/message to Master and/or Slaves.  The default behavior is
*  still to autopurge, but users can change this by the following commands.
*             piaf/message keep vectors
*             piaf/message keep histos
*             piaf/message keep none		[default]
*  The following line enables the above commands.
#define NTCMD_AUTOPURGE_SELECTABLE
*
* -R.T. Jones  May 21, 2003
*  When a .f77 or .c file is uploaded to piaf via the CSEXEC command, any
*  former instance of the corresponding executable is dissociated from the
*  master and slaves and deleted, then the new code is compiled on the
*  master into a new .sl executable.  However when a .sl or .csl file is
*  uploaded, it automatically overwrites any former instance of the
*  executable (which may be currently active) and usually crashes the piaf
*  master.  The fix is to rename the .sl or .csl code to a unique name on
*  piaf before running it.  Then when a new .sl or .csl file is uploaded,
*  the running executable can be dissociated before the new one is renamed
*  over top of the old one.  This also solves the problem of two people
*  uploading .sl or .csl files with the same name.  Select with the option
#define RENAME_SHARED_EXECS 1
*
* -R.T. Jones  May 21, 2003
*  When a .sl or .csl file is uploaded to piaf via the CSEXEC command, it
*  is immediately loaded into piafserv on the master and slaves.  However
*  when a .f77 or .c file is uploaded, it is compiled and loaded only on
*  the piaf master.  This gives an inconsistent behavior for compiled and
*  uncompiled programs.  The .sl behavior is better because it is often
*  convenient to be able to load a suite of functions via a single "call"
*  command on the client, and then access any of the functions in subsequent
*  ntuple commands.  As it is, the slaves will not recognize the functions
*  contained in the .f or .c file until after they are loaded, say by a
*  dummy ntuple pass such as nt/plot 1.suiteloader.f77 ! 1.  The behavior
*  has been changed to load all functions on piaf master and slaves in
*  response to a call command, such as call suiteloader.f77 if you have set
#define LOAD_SLAVES_ON_CSEXEC 1
*
* -R.T. Jones  April 30, 2003
*  Add a new option to HINPF command to allow client code to push a histogram
*  to the slaves without having it reset on arrival.  The option format is
*     HINPF     <id>  noreset
*  where the extra noreset option is recognized if the following is defined. 
#define HINFP_OPTION_NORESET 1
*
* -R.T. Jones  April 30, 2003
*  Add a new facility to allow a client to execute arbitrary commands on
*  piaf master or slave servers.  This is done through a specially formatted
*  piaf/message command on the client that is of the form
*     PAW> piaf/message [Master | Slaves] <kuip command>
*  This command allows any arbitrary kuip command to be executed on the piaf
*  hosts, including RSHELL, so this is a huge security hole.  There are good
*  reasons NOT to document this feature.  It is enabled by the following line.
#define GENERAL_COMMAND_MESSAGES 1
*
* -R.T. Jones  October 3, 2001
*  In standard design, when a dll is no longer needed, a unlink (rm) command
*  is issued both on the master and on each slave.  In the case where the
*  master and slaves all share a common area (eg. nfs) this produces a fault,
*  which is fixed by having only the master delete stale dll's.
*  To suppress this change, comment out the following line.
#define PIAF_SHARED_WORKDIR 1
*
* $Id: psexec.F,v 1.4 2005/11/01 16:31:11 radphi Exp $
*
* $Log: psexec.F,v $
* Revision 1.4  2005/11/01 16:31:11  radphi
* psexec.F, pslvio.F
*    - changed the field lengths on piaf/status messages to prevent
*      the 80-character message from overflowing [rtj]
*
* Revision 1.3  2005/11/01 15:34:30  radphi
* psexec.F
*    - just a spelling error in the comments [rtj]
* pslvio.F
*    - make the Format statement more robust [rtj]
* qp_execute.c, qp_hbook_if.c qp_hbook_if.h
*    - fix a bug in qp_exec_init (qp_execute.c) where it was passing a short
*      row count (PCNEXT_OPTIMIZATION only, in the case where the slave
*      computes a non-contiguous slice of the ntuple chain) and so reserving
*      too little space in the ntuple data cache to hold the data [rtj]
*    - introduced a new function h_nrows_piaf() to return the TOTAL number
*      of rows in a slave's ntuple processing slice, as opposed to the size
*      of the first contiguous segment which is returned by h_range_piaf()
*      and h_next_nt(). [rtj]
*
* Revision 1.2  2005/04/19 21:58:50  radphi
* psexec.F
*    - tell slaves to set their local comis work area to /tmp to avoid collisions
*      in filenames of work files used during comis compile step, in the case
*      where two piaf processes running on different cluster nodes happen to
*      have the same PID [rtj]
* README
*    - describe the above change [rtj]
*
* Revision 1.1.1.1  2004/03/20 20:05:12  radphi
* cernlib piaf package with enhancements by Richard.T.Jones@uconn.edu
*
* Revision 1.7  1998/05/05 08:30:15  couet
* - The last CVS comment produced some errors on AIX
*
* Revision 1.6  1997/08/04 12:29:54  dinofm
* The 'stagealloc' comment has been modified to 'Fortran' syntax 
* to allow compilation on AIX (reported by M. Dahlinger/GSI).
*
* Revision 1.5  1997/02/14 16:18:46  dinofm
* PSSYSTEM function is used to execute SHELL commands on all  platforms.
*
* Revision 1.4  1996/12/04 16:56:21  dinofm
* Take care of the new SHIFT STAGER.
*
* Revision 1.3  1996/09/11 15:01:34  couet
* - Uses hldnt
*
* Revision 1.2  1996/04/24 08:47:37  dinofm
* The PIAF/MESSAGE 'back door' handling allows to retrieve a PQRY message
* (used, e.g., to ask information to the Piaf server).
*
* Revision 1.1.1.1  1996/03/01 11:39:27  mclareni
* Paw
*
*
#include "paw/pilot.h"
*CMZ :          31/01/96  11.30.32  by  Timo Hakulinen
*-- Author :    Alfred Nathaniel   28/04/93
      SUBROUTINE PSEXEC(HANGUP)
*
* Receive a client message and execute it
*
#include "hbook/czbuff.inc"
#include "hbook/hcopt.inc"
#include "hbook/hcpiaf.inc"
#include "paw/pawchn.inc"
#include "paw/rzcxio.inc"
#include "paw/quest.inc"
*

* When compiled on SUN we assume is the CS2, thus stagealloc is used 

#if defined(CERNLIB_SUN)
#  define CERNLIB_STAGEALLOC
#endif

      INTEGER PSSYSTEM

      INTEGER HANGUP

      PARAMETER    (MAXDIR = 100)
      CHARACTER*16  CHDIRR(MAXDIR)
      CHARACTER*1   CHTYPE
      LOGICAL       PROPAG

      CHARACTER CHMAIL*80, CHMESS*6, CHLINE*74, CHOPT*14, CHDIR*60
      CHARACTER OBNAME*40, OBCLAS*20, STEXT*20, LTEXT*80, CHPID*12
      CHARACTER CHLCL*80

      SAVE      NF

#ifdef NTCMD_AUTOPURGE_SELECTABLE
      integer hkeep,vkeep
      save hkeep,vkeep
      data hkeep,vkeep/0,0/
#endif

*
* Listen to client (in case of interrupt this is necessary)
*
      CALL PFSOCK(0)
*
* In case of a SIGHUP signal try to close the master and slave servers in
* a clean way so that PALOGS will be called
*
      IF (HANGUP .NE. 0) THEN
         CHMAIL = '*CLOSE'
         GOTO 2
      ENDIF

      NSLBADO = NSLBAD
*
* Loop over client messages
*
 1    CONTINUE
*
* In case a slave crashes send a new group view to the remaining slaves
*
      IF (MASTPF .AND. NSLBAD.NE.NSLBADO) THEN
         CALL PSMGRP(NSLAVE)
         NSLBADO = NSLBAD
      ENDIF
*
* After 100 hours of inactivity Piaf will shut itself down (an orphaned
* slave will terminate after 101 hours)
*
      IF (MASTPF) THEN
         CALL PSALARM(100)
      ELSE
         CALL PSALARM(101)
      ENDIF

      CALL CZGETA(CHMAIL,ISTAT)
      IF (ISTAT.NE.0) CHMAIL = '*CLOSE'

 2    CALL DEBUG

      CHMESS=CHMAIL(:6)
      CHLINE=CHMAIL(7:)

      CALL PFLOG(8,'PSEXEC',CHMAIL)

      IF (CHMESS(:5).EQ.'KUIP ') THEN
*
* Execute a KUIP command
*
         IF (CHMESS(6:6).EQ.'+') THEN
*--- long command line
            READ(CHLINE,'(I6)') L
            CALL CZGETC(L,ISTAT)
         ELSE
            L=LENOCC(CHLINE)
            CHBUF(:L)=CHLINE
         ENDIF
#ifdef NTCMD_AUTOPURGE_SELECTABLE
         if (CHBUF(1:15).eq.'histo/delete 0;') then
            if ((hkeep.ne.0).and.(vkeep.ne.0)) then
              CHBUF = 'nt/cuts $0 -'
              L=LENOCC(CHBUF)
            else if (hkeep.ne.0) then
              CHBUF = 'nt/cuts $0 -; vec/delete *'
              L=LENOCC(CHBUF)
            else if (vkeep.ne.0) then
              CHBUF = 'nt/cuts $0 -; histo/delete 0'
              L=LENOCC(CHBUF)
            endif
         endif
#endif
         CALL PSKUIP(CHBUF(1:L),ISTAT)

      ELSEIF (CHMESS.EQ.'PACUTS') THEN
*
* Define graphical cuts
*
         CALL PSPAW(CHMESS, CHLINE, ISTAT)

      ELSEIF (CHMESS.EQ.'PCHROP') THEN
*
* Open chain member
*
         CALL PSPAW(CHMESS, CHLINE, ISTAT)

      ELSEIF (CHMESS.EQ.'GROUP:') THEN
*
* Get the current slave group view: # of active slaves and own unique id.
* This message will only be send by the master to all active slave servers.
* If the group view changes free all dynamic Ntuple buffers, since the
* range of events to be handled by each slave depends on the size of the
* group.
*
         NGOLD = NGSIZE
         READ(CHLINE,'(I5,I5)') NGSIZE, MYSID
         IF (NGSIZE .NE. NGOLD) CALL HBFREE(0)

      ELSEIF (CHMESS.EQ.'HCDIR:') THEN
*
* Change directory
*
         CHOPT=CHLINE(:14)
         CHDIR=CHLINE(15:)
         IQUEST(1)=0
         IF (INDEX(CHOPT,'R').NE.0) THEN
            L=LEN(CHDIR)
         ELSE
            L=LENOCC(CHDIR)
            IF (L.EQ.0) L=1
         ENDIF
         CALL HCDIR(CHDIR(:L),CHOPT)
         IQ1=IQUEST(1)
*
* If directory change was succesful propagate to slave servers
*
         IF (MASTPF .AND. IQ1.EQ.0) THEN
            CALL PSCAST(CHMAIL, NSLAVE, ISTAT)
            CALL PSLOOP(NSLAVE, ISTAT)
         ENDIF
*
         CALL PSPLOG(IQUEST(1),ISTAT)
         IF (MASTPF .AND. IQ1.EQ.0 .AND. INDEX(CHOPT,'R').NE.0) THEN
            CALL CZPUTA(CHDIR,ISTAT)
         ENDIF

      ELSEIF (CHMESS.EQ.'HLDIR:') THEN
*
* List directory
*
         CHOPT=CHLINE(:14)
         CHDIR=CHLINE(15:)
         L=LENOCC(CHDIR)
         IF (L.EQ.0) L=1
         CALL HLDIR(CHDIR(:L),CHOPT)
         CALL PSPLOG(IQUEST(1),ISTAT)

      ELSEIF (CHMESS .EQ. 'HRDIR:') THEN
*
* Retrieve list of subdirectories
*
         READ(CHLINE, '(I4)') MXDIR
         CALL HRDIR(MAXDIR, CHDIRR, NDIR)
         CALL PSPLOG(IQUEST(1),ISTAT)
         NDIR = MIN(MXDIR, NDIR)
         WRITE(CHBUF, '(I4)') NDIR
         DO 12 I = 1, NDIR
            CHBUF(5+16*(I-1):5+16*I-1) = CHDIRR(I)
  12     CONTINUE
         CALL CZPUTC(5+16*NDIR-1, ISTAT)

      ELSEIF (CHMESS .EQ. 'HLNEXT') THEN
*
* Return element by element the contents of a file
*
         READ(CHLINE(2:),'(I8,A)') IDH, CHOPT
         CALL HLNEXT(IDH, CHTYPE, CHDIR, CHOPT)
         IQUEST(1)  = IDH
         CHMAIL     = CHTYPE
         CHMAIL(2:) = CHDIR
         CALL PSPLOG(IQUEST(1),ISTAT)
         CALL CZPUTA(CHMAIL,ISTAT)

      ELSEIF (CHMESS .EQ. 'HGETNT') THEN
*
* Execute HNTLD in case of a chain and return the number of rows in the chain
*
* Strip-off the variable from the Ntuple ID. In case it is a function and
* it contains vectors the Comis translator tries to get the vector while
* the client does not expect that (dead-lock).
*
         I = INDEX(CHLINE, '.')
         IF (I .GT. 2) THEN
            CHDIR = CHLINE(2:I-1)
         ELSE
            CHDIR = CHLINE(2:)
         ENDIF
         CALL HNTLD(CHDIR)
         IQUEST(1) = NCHEVT
#ifdef LONG_LONG_CHAIN
*
* If the row count overflows integer*4 then lie about the length,
* because piaf will process the full chain regardless.
*
         IF (IQUEST(1).LT.0) THEN
            IQUEST(1) = 99999999
         ENDIF
#endif
         CALL PSPLOG(IQUEST(1),ISTAT)

      ELSEIF (CHMESS.EQ.'CZIN  ') THEN
*
* Read request from client's HRZIN
*
         CALL PSRZIN(CHLINE)

      ELSEIF (CHMESS.EQ.'HINPF') THEN
*
* Client wants to push a histogram to the server
*
         READ(CHLINE,'(I12)') IDH
#ifdef HINFP_OPTION_NORESET
         if (SLAVPF.and.(index(chline,'noreset').gt.0)) then
            SLAVPF = .false.
            call HINPF(idh,0)
            SLAVPF = .true.
         else
            call HINPF(idh,0)
         endif
#else
         CALL HINPF(IDH,0)
#endif
*
* Propagate to active slave servers
*
         IF (MASTPF) THEN
            DO 26 I = 1, NSLAVE
               CALL PFSOCK(I)
#ifdef HINFP_OPTION_NORESET
               call CZPUTA(chmail,istat)
               call HOUTPF(idh)
#else
               CALL PFHOUT(IDH, ISTAT)
#endif
26          CONTINUE
            CALL PFSOCK(0)
         ENDIF

      ELSEIF (CHMESS.EQ.'INQRZ:') THEN
*
* Client wants to know the record length of an RZ file
*
         LRECL=0
         CALL RZOPEN(10,'INQRZ',CHLINE,'X',LRECL,ISTAT)
         IF (ISTAT.EQ.0) THEN
            WRITE(CHMAIL,'(I8)') LRECL*4
            CALL RZEND('INQRZ')
            CLOSE(10)
         ELSE
            WRITE(CHMAIL,'(I8,A)') 0,' Cannot determine RZ record size'
         ENDIF
         CALL CZPUTA(CHMAIL,ISTAT)

      ELSEIF (CHMESS.EQ.'GETA :') THEN
*
* Client wants to receive a file
*
         CALL PSGETA(CLINPF,LOGLPF)
         CALL PSPLOG(IQUEST(1),ISTAT)

      ELSEIF (CHMESS.EQ.'PUTA :') THEN
*
* Client wants to push a file
*
         LF = 0
         IF (CHLINE .NE. ' ') THEN
            READ(CHLINE,'(I12)') LF
         ENDIF
         CALL PSPUTA(CLINPF,LOGLPF,LF)
         CALL PSPLOG(IQUEST(1),ISTAT)

      ELSEIF (CHMESS.EQ.'CSEXEC') THEN
*
* The pushed Fortran file has to be COMIS compiled
*
         CALL PAWCS

#ifdef PIAF_SHARED_WORKDIR
         IF (.NOT.MASTPF) CALL CSHLOPT('/tmp', 'path')
#endif

         PROPAG = .TRUE.
         LF = LENOCC(CHLINE)
         IF (CHLINE(LF-1:LF) .EQ. '77') PROPAG = .FALSE.
         IF (CHLINE(LF-1:LF) .EQ. '.c') PROPAG = .FALSE.

         IF (MASTPF .AND. .NOT.PROPAG) THEN
*
* Tell slaves and master to unload shared library
*
            CALL HITOC(MASPID, CHPID, LP, IERROR)
            I = INDEXB(CHLINE(:LF), '.')
            WRITE(CHSMPF,'(A)') 'CSRMSL'//CHLINE(:I-1)//'_'//CHPID(1:LP)
            CALL PSCAST(CHSMPF,  NSLVPF, ISTAT)

            CALL CSRMSL(CHSMPF(7:))
            CALL UNLINKF(CHSMPF(7:LENOCC(CHSMPF))//'.sl')
#ifdef RENAME_SHARED_EXECS
         ELSEIF (MASTPF) THEN
            CALL HITOC(MASPID, CHPID, LP, IERROR)
            I = INDEXB(CHLINE(:LF), '.')
            WRITE(CHSMPF,'(A)') 'CSRMSL'//CHLINE(:I-1)//'_'//CHPID(1:LP)
            CALL PSCAST(CHSMPF,  NSLVPF, ISTAT)

            CALL CSRMSL(CHSMPF(7:))
            CHLMPF = CHSMPF(7:LENOCC(CHSMPF))//CHLINE(I:LF)
            CALL UNLINKF(CHLMPF)
            CALL PSCP(CHLINE(:LF)//' '//CHLMPF)
            IF (CHLINE(I:LF).EQ.'.sl') THEN
              CALL PSCP(CHLINE(:I-1)//'.f '//
     +                  CHSMPF(7:LENOCC(CHSMPF))//'.f')
            ELSEIF (CHLINE(I:LF).EQ.'.csl') THEN
              CALL PSCP(CHLINE(:I-1)//'.c '//
     +                  CHSMPF(7:LENOCC(CHSMPF))//'.c')
            ENDIF
            CHLINE = CHLMPF
#endif
         ENDIF

         CHLMPF = CHLINE
         CALL KUHOME(CHLMPF, LF)
         CALL CSEXEC('!FILE '//CHLMPF,ISTAT)
         IQUEST(1)=ISTAT

*
* If Comis compilation was successful compile also on slave servers
*
         IF (MASTPF .AND. PROPAG .AND. ISTAT.EQ.0) THEN
#ifdef RENAME_SHARED_EXECS
            CALL PSCAST('CSEXEC'//CHLMPF, NSLAVE, ISTAT)
#else
            CALL PSCAST(CHMAIL, NSLAVE, ISTAT)
#endif
            CALL PSLOOP(NSLAVE, ISTAT)
#ifdef LOAD_SLAVES_ON_CSEXEC
         ELSEIF (MASTPF .AND. ISTAT.EQ.0) THEN
            CALL PSCAST('CSEXEC'//CHLMPF, NSLAVE, ISTAT)
            CALL PSLOOP(NSLAVE, ISTAT)
#endif
         ENDIF
*
         CALL PSPLOG(IQUEST(1),ISTAT)

      ELSEIF (CHMESS.EQ.'CSRMSL') THEN
*
* Unload shared library
* This message will only be send by the master to all slave servers
*
         CHLMPF = CHLINE
         CALL KUHOME(CHLMPF, LF)
         CALL CSRMSL(CHLMPF)
#ifndef PIAF_SHARED_WORKDIR
         CALL UNLINKF(CHLINE(1:LENOCC(CHLINE))//'.sl')
#endif

      ELSEIF (CHMESS.EQ.'SHELL:') THEN
*
* Execute a local shell command
*
*-- avoid security leak, e.g 'ls; rm -r /'
         L=LENOCC(CHLINE)
         IF (INDEX(CHLINE(:L),';').EQ.0 .AND.
     +       INDEX(CHLINE(:L),'|').EQ.0 .AND.
     +       INDEX(CHLINE(:L),'&').EQ.0 .AND.
     +       INDEX(CHLINE(:L),'>').EQ.0 .AND.
     +       INDEX(CHLINE(:L),'<').EQ.0 .AND.
     +       INDEX(CHLINE(:L),'`').EQ.0 .AND.
     +      (CHLINE(:4).EQ.'cat' .OR. CHLINE(:3).EQ.'cp' .OR.
     +       CHLINE(:3).EQ.'ls'  .OR. CHLINE(:3).EQ.'mv' .OR.
     +       CHLINE(:3).EQ.'rm'  .OR. CHLINE(:4).EQ.'pwd' .OR.
     +       CHLINE(:6).EQ.'mkdir' .OR. CHLINE(:6).EQ.'rmdir')) THEN
#if defined(CERNLIB_CORE)||defined(CERNLIB_STAGEALLOC)
*
* Special case for rm and DPM.
*
            IF (CHLINE(:3) .EQ. 'rm') THEN
               CALL PSRM(CHLINE(3:L))
*
* Special case for cp and DPM.
*
            ELSEIF (CHLINE(:3) .EQ. 'cp') THEN
               CALL PSCP(CHLINE(3:L))
*
* Special case for mv and DPM.
*
            ELSEIF (CHLINE(:3) .EQ. 'mv') THEN
               CALL PSMV(CHLINE(3:L))
            ELSE
               IQUEST(1)=PSSYSTEM(CHLINE(:L))
            ENDIF
#else
            IQUEST(1)=PSSYSTEM(CHLINE(:L))
#endif
         ELSE
            PRINT *,' Invalid command: ',CHLINE(:L)
            IQUEST(1)=1
         ENDIF
         CALL PSPLOG(IQUEST(1),ISTAT)

      ELSEIF (CHMESS.EQ.'MESS :') THEN
*
* Back door for new commands
*
         L = LENOCC(CHLINE)
         CALL DECODEMESS (CHLINE,IRET)
         IF (IRET.EQ.0) THEN
#ifdef GENERAL_COMMAND_MESSAGES
            if (chline(1:6).eq.'Master') then
               call kuexec(chline(8:L))
            else if (chline(1:6).eq.'Slaves') then
               if (mastpf) then
                  call pscast(chmail,nslave,istat)
                  call psloop(nslave,istat)
               else
                  call kuexec(chline(8:L))
               endif
# ifdef NTCMD_AUTOPURGE_SELECTABLE
            else if (chline(1:9).eq.'keep none') then
               vkeep = 0
               hkeep = 0
            else if (chline(1:9).eq.'keep vect') then
               vkeep = 1
            else if (chline(1:9).eq.'keep hist') then
               hkeep = 1
# endif
# ifdef WEIGHTED_NTUPLE_PARTITIONING
            else if (chline(1:6).eq.'reload') then
               read(chline(7:),*,err=6) islave,iweight
               call setNTweight(islave,iweight)
               if (mastpf) then
                  call pscast(chmail,nslave,istat)
                  call psloop(nslave,istat)
               endif
    6          continue
# endif
            else
               PRINT *,' Invalid message: ',CHLINE(:L)
            endif
#else
            PRINT *,' Invalid message: ',CHLINE(:L)
#endif
         ELSE
            PRINT *,CHLINE
         ENDIF
         IQUEST(1) = 0
         CALL PSPLOG(IQUEST(1),ISTAT)

      ELSEIF (CHMESS.EQ.'PING :') THEN
*
* Noop, just to check if server is still alive
*

      ELSEIF (CHMESS.EQ.'HCOPT:') THEN
*
* Set HCOPT words
*
         READ(CHLINE,'(10I5)') (ICOPT(I),I=1,10)
*
* Propagate to slave servers
*
         IF (MASTPF) THEN
            CALL PSCAST(CHMAIL, NSLAVE, ISTAT)
         ENDIF

      ELSEIF (CHMESS.EQ.'LOGLEV') THEN
*
* Set log level
*
         READ(CHLINE,'(I12)') LOGLPF
*
* Propagate log level to slave servers
*
         IF (MASTPF) THEN
            CALL PSCAST(CHMAIL, NSLAVE, ISTAT)
         ENDIF

      ELSEIF (CHMESS.EQ.'STATUS') THEN
*
* Report current server status, including IO and CPU statistics of slaves
*
         CALL PSLVIO(NSLVPF, ISTAT)
         TOTIO = RZXIN + RZXOUT + SRIN + SROUT
         CALL PSSTAT(NSLAVE,TOTIO,SVMIO,STUSER,STCPU)
         CALL PSPLOG(IQUEST(1),ISTAT)

      ELSEIF (CHMESS.EQ.'IOSTAT') THEN
*
* Collect the I/O (real and virtual) and CPU statistics from all slaves
* This message will only be sent by the master to all slave servers
*
            CALL HBVM(VIO)
            CALL PSCPU(TUSER, TCPU)
            WRITE(CHSMPF,'(3F15.0,2F15.3)') RZXIN, RZXOUT, VIO, TUSER,
     +                                      TCPU
            CALL CZPUTA(CHSMPF, ISTAT)

      ELSEIF (CHMESS.EQ.'GETCWD') THEN
*
* Return current working directory
*
         CALL GETWDF(CHMAIL)
         DO 10 I = LENOCC(CHMAIL), 1, -1
            IF (CHMAIL(I:I) .EQ. '/') THEN
               CALL CZPUTA('~'//CHMAIL(I:),ISTAT)
               GOTO 1
            ENDIF
  10     CONTINUE
         CALL CZPUTA(CHMAIL,ISTAT)

      ELSEIF (CHMESS.EQ.'SCAND:') THEN
*
* Return one by one the files in the Piaf directory
*
         IF (CHLINE .NE. ' ') THEN
            NF   = 0
            LP   = LENOCC(CHLINE)
            CALL PSSCAN(CHLINE, LP, NF, OBNAME, OBCLAS, STEXT, LTEXT)
         ELSE
            NF = NF + 1
            CALL PSSCAN(' ', 0, NF, OBNAME, OBCLAS, STEXT, LTEXT)
         ENDIF
         CALL CZPUTA(OBNAME//OBCLAS//STEXT, ISTAT)
         IF (OBNAME .NE. ' ') THEN
            CALL CZPUTA(LTEXT, ISTAT)
         ENDIF

      ELSEIF (CHMESS.EQ.'MODE :') THEN
*
* Switch server mode
*
         IF (CHLINE(1:1).EQ.'?') THEN
            IF (NSLAVE .EQ. 0) THEN
               PRINT *,' Piaf server running in sequential mode'
            ELSE
               PRINT 10000, NSLAVE
            ENDIF
         ELSEIF (CHLINE(1:3).EQ.'SEQ') THEN
            NSLAVE = 0
            PRINT *,' Piaf server set to sequential mode'
         ELSEIF (CHLINE(1:3).EQ.'PAR') THEN
            READ(CHLINE(4:),'(I12)') NSLAV
            CALL PSPAR(NSLAV, 1)
         ELSE
            PRINT *,' Invalid Piaf server mode ', CHLINE(1:3)
         ENDIF
*
* Tell client that mode change has finished
*
         CALL CZPUTA('MODE*',ISTAT)
         CALL PSPLOG(IQUEST(1),ISTAT)

      ELSEIF (CHMESS.EQ.'STAGE:') THEN
*
* Submit a staging request
*
         CALL CZGETA(CHMAIL,ISTAT)
         IF (INDEX(CHLINE(:8),'L').NE.0) THEN
            CALL CZGETA(CHLCL,ISTAT)
         ELSE
            CHLCL=' '
         ENDIF
         CALL PSTAGE(CHMAIL,CHLINE(9:),CHLINE(:8),CHLCL)
         CALL PSPLOG(IQUEST(1),ISTAT)

      ELSEIF (CHMESS.EQ.'*CLOSE') THEN
*
* Close all slaves
*
         IF (MASTPF) THEN
            CALL PSLVIO(NSLVPF, ISTAT)
            CALL PSCAST(CHMESS, NSLVPF, ISTAT)
         ENDIF
*
* Close connection
*
         GOTO 9

      ELSE
*
* Unknown message: close connection
*
         L=LENOCC(CHMAIL)
         IF (L.EQ.0) L=1
         PRINT '(2A)','  Unknown message: ',CHMAIL(:L)
         CHMAIL = '*CLOSE'
         GOTO 2
      ENDIF
*
      GOTO 1
*
 9    CONTINUE
*
10000 FORMAT('  Piaf server running in parallel mode (',I1,' slaves)')
*
      END

# ifdef WEIGHTED_NTUPLE_PARTITIONING
      subroutine setNTweight(islave,iweight)
      integer islave,iweight
      integer max_slave_weights
      parameter (max_slave_weights=256)
      common /slave_weights/slave_weight(max_slave_weights)
      integer slave_weight
      save /slave_weights/
      data slave_weight/max_slave_weights*1/
      if (islave.lt.0.or.islave.gt.max_slave_weights) then
        print *, '  Invalid slave index',islave
      elseif (iweight.lt.1) then
        print *, '  Invalid weight factor',iweight
      elseif (islave.gt.0) then
        slave_weight(islave) = iweight
      else
        do is=1,max_slave_weights
          slave_weight(is) = iweight
        enddo
      endif
      end

      integer function getNTweight(islave)
      integer islave
      integer max_slave_weights
      parameter (max_slave_weights=256)
      common /slave_weights/slave_weight(max_slave_weights)
      integer slave_weight
      if (islave.le.0.or.islave.gt.max_slave_weights) then
        getNTweight = 0
      else
        getNTweight = slave_weight(islave)
      endif
      end
#endif

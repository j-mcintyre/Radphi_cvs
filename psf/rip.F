*---------------------------------------------
* rip.F - a stand-alone host for paw packages
* ---------
* Richard Jones
* September 2, 2002
*---------------------------------------------
      program rip
      implicit none
      integer hspace
      parameter (hspace=20000000)
      common /pawc/ hq(hspace)
      common /quest/ iquest(100)
      real hq
      integer iquest
      integer istat

C---- initialization of HBOOK
      call HLIMIT(hspace)

C---- do user processing here
      call lgdStrings('rip%d.out')

C---- Termination of HBOOK
#ifdef RZ_FILE
      iquest(10) = 256000  ! extend RZ quota to 2^32
      call HROPEN(50,'RZfile',RZ_FILE,'NQE',RZ_BSIZE,istat)
      call HROUT(0,icycle,' ')
      call HREND('RZfile')
#endif
      call exit(0)
      END

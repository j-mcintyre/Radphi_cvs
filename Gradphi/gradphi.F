      implicit none
      integer ispace,hspace
      parameter (ispace=2000000)
      parameter (hspace=50000)
      common /gcbank/ q(ispace)
      common /pawc/ hq(hspace)
      common /quest/ iquest(100)
      real q,hq
      integer iquest
      real secmax
      parameter (secmax=300000.)
      integer istat,icycle

C---- Initialization of HBOOK, ZEBRA, clock
      call GZEBRA(ispace)
      call HLIMIT(-hspace)
      call TIMEST(secmax)

C---- Initialization of GEANT
      call HPLINT(0)
      call UGINIT

C---- Simulation
      iquest(10) = 256000  ! extend RZ quota to 2^32
      call HROPEN(50,'RZfile','geant.hbook','NQE',1024,istat)
        call GRUN
        call HROUT(0,icycle,' ')
      call HREND('RZfile')

C---- Termination ----
      CALL UGLAST
      WRITE(0,999)
  999 FORMAT(//,' ***** NORMAL EXIT FROM RADPHI MONTECARLO *****',//)
      END

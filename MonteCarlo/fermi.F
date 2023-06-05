      real pFermiRMS(2)
      real pFermiPart(2)
c     data pFermiRMS/0.090,0.400/
      data pFermiRMS/0.093,0.320/
c     data pFermiPart/0.450,1.0/
      data pFermiPart/0.900,1.0/
      real pi
      parameter (pi=3.141592654)
      integer comp, part
      real p
      real pRL(5)
      integer i
      do i=1,1000000
      p=rand(0)
      do part=1,2
        if (p.le.pFermiPart(part)) goto 5
      enddo
    5 continue
      do comp=1,3
        rho=sqrt(-2.*alog(rand(0)))
        phi=2*pi*rand(0)
        pRL(comp)=pFermiRMS(part)*rho*cos(phi)
      enddo
      pRL(5)=sqrt(pRL(1)**2+pRL(2)**2+pRL(3)**2)
      pRL(4)=sqrt(pRL(5)**2+mR**2)
      print *, pRL
      end do
      end


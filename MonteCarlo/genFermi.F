C-------------------------------------------------------------------------
C genFermi - generates a diffractive photoproduction event from a nucleus.
C	The reaction is written as  0 + A -> 1 + N + R  in the lab where A
C       where A is a nuclear target, N is an ejected nucleon and R is the
C	recoiling residual (a-1) nucleus.  It generates the scattering angle
C	according to the diffractive cross section exp(bt), and the azimuthal
C       angle of the scattering plane is generated uniformly [0,2pi].
C
C	The model assumes quasi-free scattering from a single nucleon, and
C	uses a Gaussian approximation for the nucleon momentum distribution
C	in the nucleus.  There is no final-state interaction between the
C	recoil nucleon and the residual	nucleus.
C
C Modification history:
C - October 14, 2004 -rtj
C       modified the Fermi distribution to agree with a double-Gaussian
C       fit to the real data transverse momentum distribution.
C - March 29, 2001 -rtj
C       repaired a bug in the use of pFermiRMS. The value provided was a
C       total r.m.s. whereas it was being used as the r.m.s. of a projection
C       along one of the three axes.  Now it is divided by sqrt(3).
C-------------------------------------------------------------------------
      subroutine genFermi(pbeam,b,m1,m2,p1L,p2L)
      implicit none
      real pbeam(4),b
      real p1L(5),p2L(5)
      real m1,m2
      real s,t,t0
      real rho,phi
      real L(4,4)
      real beta(3),beta2
      real gamma,gfact
      real phi0L,phi1L
      real p0L(5),pRL(5)
      real p0cm(5),p1cm(5)
      real sintcm,costcm,phicm
      real m0,mN,mR
      data m0,mR,mN/8.39479, 7.480, 0.93828/
      real pFermiRMS(2)  ! Gaussian RMS for Fermi motion, 2 Gaussian fit
      real pFermiPart(2) ! respective weights for 2-Gaussian fit
c     data pFermiRMS/0.00110,0.00360/  ! uncomment to suppress Fermi motion
c      data pFermiRMS/0.110,0.360/     ! original 2-Gaussian setup
c      data pFermiPart/0.42,1.0/
c      data pFermiRMS/0.110,0.290/      ! make contribution from wider
c      data pFermiPart/0.70,1.0/        ! Gaussian smaller (mcwrap1.3)
c      data pFermiRMS/0.110,0.440/       ! very wide Gaussian 
c      data pFermiPart/0.70,1.0/         ! (mcwrap1.4)
c      data pFermiRMS/0.090,0.400/       ! wide and slim Gaussian 
c      data pFermiPart/0.450,1.0/         ! (mcwrap1.5) to run with b=6
c      data pFermiRMS/0.090,0.290/       ! wide and slim Gaussian 
c      data pFermiPart/0.400,1.0/         ! (mcwrap1.6) to run with b=2
c      data pFermiRMS/0.090,0.400/       ! wide and slim Gaussian 
c      data pFermiPart/0.500,1.0/         ! (mcwrap1.7) to run with b=4
c      data pFermiRMS/0.100,0.400/       ! wide and slim Gaussian 
c      data pFermiPart/0.500,1.0/         ! (mcwrap1.8) to run with b=6
      data pFermiRMS/0.093,0.320/       ! wide and slim Gaussian 
      data pFermiPart/0.900,1.0/         ! (mcwrap1.8) to run with b=6
      real pi
      parameter (pi=3.141592654)
      real ranf
      external ranf
      real Ekin,pperp
      integer part,comp
      real p
      real uhat(3),vhat(3),what(3),vXu(3)
      real vhatXuhat,uhatXzhat

C start again if kinematics fall below threshold
    1 continue

C generate the momentum of the recoil nucleus; this part does not scatter
      p=ranf()
      do part=1,2
        if (p.le.pFermiPart(part)) goto 5
      enddo
    5 continue
      do comp=1,3
        rho=sqrt(-2.*alog(ranf()))
        phi=2*pi*ranf()
        pRL(comp)=pFermiRMS(part)*rho*cos(phi)
      enddo
      pRL(5)=sqrt(pRL(1)**2+pRL(2)**2+pRL(3)**2)
      pRL(4)=sqrt(pRL(5)**2+mR**2)

C find the the c.m. frame of the (gamma,N) reaction
      p0L(1)=pbeam(1)-pRL(1) 
      p0L(2)=pbeam(2)-pRL(2) 
      p0L(3)=pbeam(3)-pRL(3)
c     p0L(4)=pbeam(4)+sqrt(mN**2+pRL(5)**2)	!on-shell approximation
c     p0L(4)=pbeam(4)+sqrt(0.9**2+pRL(5)**2)	!off-shell approximation
      p0L(4)=pbeam(4)+m0-pRL(4)			!energy conservation
      p0L(5)=sqrt(p0L(1)**2+p0L(2)**2+p0L(3)**2)
      s=p0L(4)**2-p0L(5)**2
      beta(1)=p0L(1)/p0L(4)
      beta(2)=p0L(2)/p0L(4)
      beta(3)=p0L(3)/p0L(4)
      beta2=beta(1)**2+beta(2)**2+beta(3)**2
      gamma=1/Sqrt(1-beta2)
      gfact=(gamma-1)/beta2
      L(1,1)=1+beta(1)*beta(1)*gfact
      L(1,2)=beta(1)*beta(2)*gfact
      L(1,3)=beta(1)*beta(3)*gfact
      L(2,2)=1+beta(2)*beta(2)*gfact
      L(2,3)=beta(2)*beta(3)*gfact
      L(3,3)=1+beta(3)*beta(3)*gfact
      L(1,4)=-gamma*beta(1)
      L(2,4)=-gamma*beta(2)
      L(3,4)=-gamma*beta(3)
      L(4,4)=gamma

c boost the initial beam particle into the (gamma,p) c.m. frame
      p0cm(1)=L(1,1)*pbeam(1)+L(1,2)*pbeam(2)+L(1,3)*pbeam(3)
     +       +L(1,4)*pbeam(4)
      p0cm(2)=L(1,2)*pbeam(1)+L(2,2)*pbeam(2)+L(2,3)*pbeam(3)
     +       +L(2,4)*pbeam(4)
      p0cm(3)=L(1,3)*pbeam(1)+L(2,3)*pbeam(2)+L(3,3)*pbeam(3)
     +       +L(3,4)*pbeam(4)
      p0cm(4)=L(1,4)*pbeam(1)+L(2,4)*pbeam(2)+L(3,4)*pbeam(3)
     +       +L(4,4)*pbeam(4)
      p0cm(5)=p0cm(4)

c generate t and the angles in the c.m. frame
      if (s.lt.0) then				!this is rare, but possible
        go to 1
      endif
      p1cm(4)=(s+m1**2-m2**2)/(2*sqrt(s))
      if (p1cm(4).lt.m1) then			!below threshold, try again
        go to 1
      endif
      p1cm(5)=sqrt(p1cm(4)**2-m1**2)
      t0=(p0cm(4)-p1cm(4))**2-(p0cm(5)-p1cm(5))**2
      t=t0+alog(ranf()+1e-30)/b
      costcm=1+(t-t0)/(2*P1cm(5)*P0cm(5))
      if (costcm.lt.-1) then			!beyond tmax, try again
        go to 1
      endif
      sintcm=sqrt(1-costcm**2)
      phicm=2*pi*ranf()

c generate the scattered direction of the forward system in cm frame
      uhat(1)=p0cm(1)/p0cm(5)
      uhat(2)=p0cm(2)/p0cm(5)
      uhat(3)=p0cm(3)/p0cm(5)
      uhatXzhat=sqrt(uhat(1)**2+uhat(2)**2)
      vhat(1)=uhat(2)/uhatXzhat
      vhat(2)=-uhat(1)/uhatXzhat
      vhat(3)=0
      vXu(1)=vhat(2)*uhat(3)-vhat(3)*uhat(2)
      vXu(2)=vhat(3)*uhat(1)-vhat(1)*uhat(3)
      vXu(3)=vhat(1)*uhat(2)-vhat(2)*uhat(1)
      vhatXuhat=sqrt(vXu(1)**2+vXu(2)**2+vXu(3)**2)
      what(1)=vXu(1)/vhatXuhat
      what(2)=vXu(2)/vhatXuhat
      what(3)=vXu(3)/vhatXuhat
      p1cm(1)=p1cm(5)*(sintcm*(cos(phicm)*vhat(1)+sin(phicm)*what(1))
     +                +costcm*uhat(1))
      p1cm(2)=p1cm(5)*(sintcm*(cos(phicm)*vhat(2)+sin(phicm)*what(2))
     +                +costcm*uhat(2))
      p1cm(3)=p1cm(5)*(sintcm*(cos(phicm)*vhat(3)+sin(phicm)*what(3))
     +                +costcm*uhat(3))

c boost back to the lab frame
      p1L(1)=L(1,1)*p1cm(1)+L(1,2)*p1cm(2)+L(1,3)*p1cm(3)-L(1,4)*p1cm(4)
      p1L(2)=L(1,2)*p1cm(1)+L(2,2)*p1cm(2)+L(2,3)*p1cm(3)-L(2,4)*p1cm(4)
      p1L(3)=L(1,3)*p1cm(1)+L(2,3)*p1cm(2)+L(3,3)*p1cm(3)-L(3,4)*p1cm(4)
      p1L(4)=L(4,4)*p1cm(4)-L(1,4)*p1cm(1)-L(2,4)*p1cm(2)-L(3,4)*p1cm(3)

c calculate the recoil nucleon momentum by subtraction
      p2L(1)=p0L(1)-p1L(1)
      p2L(2)=p0L(2)-p1L(2)
      p2L(3)=p0L(3)-p1L(3)
      p2L(4)=p0L(4)-p1L(4)
      p2L(5)=sqrt(p2L(1)**2+p2L(2)**2+p2L(3)**2)
      pperp=sqrt(p2L(1)**2+p2L(2)**2)
      Ekin=sqrt(p2L(5)**2+m2**2)-m2
      if (Ekin.ne.Ekin) then
        write(6,*) 'genFermi Fatal error 1 - this should never happen!'
      endif

c     call ipm(1,p2L(3),pperp)
c     if ((atan2(pperp,p2L(3))*180/3.1416.gt.32).and.
c    +    (atan2(pperp,p2L(3))*180/3.1416.lt.90)) then
c       call hfill(9,Ekin*exp((72-35*p2L(3)/pperp)/120),0.,1.)
c     endif
      end

      real function ranf()
      real rndm(1)
      call ranlux(rndm, 1)
      ranf=rndm(1)
      end

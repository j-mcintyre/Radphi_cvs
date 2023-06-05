      Subroutine c
      implicit none
c     include "geometry.inc"
c     include "cerenkov.inc"
      Integer nbins
      Parameter (nbins=51)
      Real Ephot(nbins)
      Real eff0(nbins),eff1(nbins),effPMT(nbins)
      Real rindexLGD(nbins),alenLGD(nbins)
      Real rindexWrap(nbins),alenWrap(nbins)
      Real rindexGap(nbins),alenGap(nbins)
      Real rindexAIR(nbins),alenAIR(nbins)
      Real rindexPMT(nbins),alenPMT(nbins)
      Common /Cerenkov/Ephot,eff0,eff1,effPMT,rindexLGD,alenLGD,
     +                 rindexWrap,alenWrap,rindexGap,alenGap,
     +                 rindexAIR,alenAIR,rindexPMT,alenPMT
c     save /Cerenkov/
      Data eff0/nbins*0.0/
      Data eff1/nbins*1.0/
      Data effPMT/
     + 0.0001,0.0001,0.0002,0.0002,0.0002,0.0002,0.0002,0.0001,
     + 0.0001,0.0001,0.0002,0.0004,0.0008,0.0015,0.0027,0.0045,
     + 0.0070,0.0104,0.0148,0.0203,0.0268,0.0346,0.0435,0.0535,
     + 0.0645,0.0763,0.0888,0.1020,0.1156,0.1296,0.1438,0.1582,
     + 0.1727,0.1873,0.2016,0.2156,0.2290,0.2412,0.2519,0.2604,
     + 0.2662,0.2691,0.2689,0.2657,0.2598,0.2512,0.2389,0.2208,
     + 0.1943,0.1583,0.1135/
      Data alenLGD/
     + 275.83,258.83,258.83,251.08,251.08,275.83,275.83,275.83,
     + 275.83,275.83,275.83,267.07,275.83,251.08,251.08,258.83,
     + 236.88,236.88,251.08,267.07,275.83,275.83,275.83,275.83,
     + 275.83,275.83,258.83,251.08,243.77,236.88,230.35,224.16,
     + 212.72,207.42,202.37,188.57,184.37,176.49,162.55,153.43,
     + 127.07, 91.92, 61.65, 43.27, 24.58, 11.82,  4.75,  1.69,
     +   1.60,  1.60,  1.60/
      Data rindexLGD/nbins*1.62/
c     Data alenGap/nbins*0.001/		! black paper wrapping (air gap)
c     Data rindexGap/nbins*1.00/	! is a very absorptive dielectric
c     Data alenGap/nbins*0.10/		! aluminum wrapping (no air gap)
c     Data rindexGap/nbins*0.00/	! is a reflective metal surface
      Data alenGap/nbins*30000./	! transparent air gap
      Data rindexGap/nbins*1.00/	! for total-internal-reflection
      Data alenWrap/nbins*0.10/		! aluminum wrapping (outside gap)
      Data rindexWrap/nbins*0.00/	! is a reflective metal surface
      Data alenAIR/nbins*30000./
      Data rindexAIR/nbins*1.00/
      Data alenPMT/nbins*0.001/
      Data rindexPMT/nbins*1.52/
      Data Ephot/
     + 1.54e-9,1.56e-9,1.58e-9,1.60e-9,1.62e-9,1.65e-9,1.67e-9,1.69e-9,
     + 1.72e-9,1.74e-9,1.76e-9,1.79e-9,1.82e-9,1.84e-9,1.87e-9,1.90e-9,
     + 1.93e-9,1.96e-9,1.99e-9,2.02e-9,2.06e-9,2.09e-9,2.13e-9,2.17e-9,
     + 2.21e-9,2.25e-9,2.29e-9,2.33e-9,2.37e-9,2.42e-9,2.47e-9,2.52e-9,
     + 2.57e-9,2.63e-9,2.68e-9,2.74e-9,2.81e-9,2.87e-9,2.94e-9,3.01e-9,
     + 3.09e-9,3.17e-9,3.25e-9,3.34e-9,3.43e-9,3.53e-9,3.63e-9,3.74e-9,
     + 3.86e-9,3.98e-9,4.12e-9/
      real cereN0(180)
      data cereN0/
     +     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     +     0.0000,0.0000,0.0000,0.0370,0.0740,0.1110,0.2220,
     +     1.0360,1.8130,2.4050,2.8490,3.2190,3.5150,3.8110,
     +     4.0330,4.1810,4.4030,4.5140,4.6250,4.7360,4.8470,
     +     4.9580,5.0320,5.0690,5.1060,5.1800,5.2170,5.1800,
     +     5.2170,5.2170,5.1800,5.1800,5.2170,5.1800,5.1060,
     +     5.1060,4.9950,4.9210,4.8470,4.7730,4.6990,4.6250,
     +     4.5510,4.5140,4.4770,4.3660,4.3660,4.3660,4.3290,
     +     4.2920,4.2920,4.2920,4.2180,4.1810,4.1810,4.1440,
     +     4.0700,4.0330,3.8850,3.8110,3.7370,3.6630,3.5520,
     +     3.4410,3.4410,3.2560,3.1450,3.0710,2.9230,2.8120,
     +     2.6640,2.5160,2.3680,2.1460,2.0350,1.8500,1.6650,
     +     1.4430,1.1840,0.9250,0.6290,0.3700,0.0740,90*0.0/

      integer i
      vector eff(51)
      vector enu(51)
      vector att(51)
      vector cn0(180)
      do i=1,51
         eff(i) = effPMT(i)
         enu(i) = Ephot(i)
         att(i) = alenLGD(i)
      enddo
      do i=1,180
         cn0(i) = cereN0(i)
      enddo
      end

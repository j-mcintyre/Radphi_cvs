c------------------------------------------------------------------------
c Defines materials for detector description
c----------------------------------------------------------------------
c Author: Anmdriy Kurylov
c
c Appended by Tom Bogue for the summer 99 update of the monte carlo
c
      subroutine matdef
      implicit none
#     include "geometry.inc"
      character*20 namate
      integer imate,nwbuf,nlmat
      real dens,radl,absl,ubuf
      real A(10),Z(10),wmat(10)
      nwbuf=0   
c plast
      imate=mat_PLAST
      namate='plast'
      A(1)=12.0111
      A(2)=1.00797
      Z(1)=6.0
      Z(2)=1.0
      dens=1.049
      nlmat=2
      wmat(1)=0.91549
      wmat(2)=0.08451
      call GSMIXT(imate,namate,A,Z,dens,nlmat,wmat)
c CsI
      imate=mat_CsI
      namate='CsI'
      A(1)=132.905
      A(2)=126.9044 
      Z(1)=55.0
      Z(2)=53.0
      dens=4.526
      nlmat=2
      wmat(1)=0.511548
      wmat(2)=0.488452
      call GSMIXT(imate,namate,A,Z,dens,nlmat,wmat)
      
c acrylic
      imate=mat_ACRYLIC
      namate='acrylic'
      A(1)=12.0111         
      A(2)=1.00797
      A(3)=15.9994
      Z(1)=6.0
      Z(2)=1.0
      Z(3)=8.0
      dens=1.181
      nlmat=3
      wmat(1)=0.59985
      wmat(2)=0.08054
      wmat(3)=0.31961
      call GSMIXT(imate,namate,A,Z,dens,nlmat,wmat)
      
c mylar = C5 H4 02 (also stands in for Kaptan)
      imate=mat_mylar
      namate='mylar'
      A(1)=12.011
      A(2)=1.008
      A(3)=15.999
      Z(1)=6.
      Z(2)=1.
      Z(3)=8.
      dens=1.39
      nlmat=-3
      wmat(1)=5.
      wmat(2)=4.
      wmat(3)=2.
      call GSMIXT(imate,namate,A,Z,dens,nlmat,wmat)
      
c ch2
      imate=mat_CH2
      namate='ch2'
      A(1)=12.0111
      A(2)=1.00797
      Z(1)=6.0
      Z(2)=1.0
      dens=0.935
      nlmat=-2
      wmat(1)=1.
      wmat(2)=2.
      call GSMIXT(imate,namate,A,Z,dens,nlmat,wmat)
      
C pyrex
      imate=mat_PYREX
      namate='pyrex'
      A(1)=28.09
      A(2)=10.82
      A(3)=23.00
      A(4)=16.00
      Z(1)=14.
      Z(2)=5.
      Z(3)=11.
      Z(4)=8.
      dens=2.23
      nlmat=+4
      wmat(1)=0.385
      wmat(2)=0.039
      wmat(3)=0.038
      wmat(4)=0.538
      call GSMIXT(imate,namate,A,Z,dens,nlmat,wmat)
      
C shield - looks like concrete, mostly dirt
      imate=mat_SHIELD
      namate='shield'
      A(1)=16.00
      A(2)=28.09
      A(3)=40.08
      A(4)=22.99
      A(5)=55.85
      A(6)=26.98
      Z(1)=8.
      Z(2)=14.
      Z(3)=20.
      Z(4)=11.
      Z(5)=26.
      Z(6)=13.
      dens=2.50
      nlmat=+6
      wmat(1)=0.52
      wmat(2)=0.325
      wmat(3)=0.06
      wmat(4)=0.015
      wmat(5)=0.04
      wmat(6)=0.04
      call GSMIXT(imate,namate,A,Z,dens,nlmat,wmat)
      
C lead-glass: SF3
      imate=mat_LEAD_GLASS_SF3
      namate='lead-glass'
      A(1)=207.19
      A(2)=28.09
      A(3)=39.10
      A(4)=22.99
      A(5)=16.00
      Z(1)=82.
      Z(2)=14.
      Z(3)=19.
      Z(4)=11.
      Z(5)=8.
      dens=3.6
      nlmat=-5
      wmat(1)=0.068
      wmat(2)=0.241
      wmat(3)=0.075
      wmat(4)=0.019
      wmat(5)=0.597
      call GSMIXT(imate,namate,A,Z,dens,nlmat,wmat)
      
c helium 
      imate=mat_HELIUM
      namate='helium-gas'
      A(1)=4.00
      Z(1)=2.
      dens=0.178e-3
      radl=5.30e5
      nwbuf=0
      call GSMATE(imate,namate,A,Z,dens,radl,absl,ubuf,nwbuf)
      
c silicon (for silicon trackers)
      imate=mat_silicon
      namate='silicon'
      A(1)=28.09
      Z(1)=14.
      dens=2.33
      radl=9.36
      nwbuf=0
      call GSMATE(imate,namate,A,Z,dens,radl,absl,ubuf,nwbuf)
      
C G10 circuit board: 60% SiO2, 40% epoxy; guess at epoxy = CH
      imate=mat_G10
      namate='G10-board'
      A(1)=28.09
      A(2)=16.00
      A(3)=12.01
      A(4)=1.01
      Z(1)=14.
      Z(2)=8.
      Z(3)=6.
      Z(4)=1.
      dens=1.39
      nlmat=+4
      wmat(1)=0.28
      wmat(2)=0.32
      wmat(3)=0.37
      wmat(4)=0.03
      call GSMIXT(imate,namate,A,Z,dens,nlmat,wmat)
      
c nickel (for collimators)
      imate=mat_nickel
      namate='nickel'
      A(1)=58.69
      Z(1)=28.
      dens=8.90
      nlmat=1
      wmat(1)=1.
      call GSMIXT(imate,namate,A,Z,dens,nlmat,wmat)      
c
c begin section for summer 99 update
c
c lead_scifi in the BGV
      imate=mat_lead_scifi
      A(1)=207.19
      A(2)=121.76
      A(3)=12.01
      A(4)=1.01
      Z(1)=82.0
      Z(2)=51.0
      Z(3)=6.0
      Z(4)=1.0
      dens=4.58
      wmat(1)=0.806
      wmat(2)=0.051
      wmat(3)=0.131
      wmat(4)=0.012
      call GSMIXT(imate,'lead_scifi',A,Z,dens,4,wmat)
c
c begin section for summer 2020 update
c
c carbon fiber epoxy in the BSD
      imate=mat_carbon_fiber_epoxy
      A(1)=12
      A(2)=16
      A(3)=1
      A(4)=14
      Z(1)=6
      Z(1)=8
      Z(1)=1
      Z(1)=7
c      dens=2.0
      dens=10.0
      wmat(1)=0.85
      wmat(2)=0.06
      wmat(3)=0.03
      wmat(4)=0.06
      call GSMIXT(imate,'carbon_fiber_epoxy',A,Z,dens,4,wmat)
c
c scintillator wrapping in the BSD
      imate=mat_scint_wrap
      A(1)=12
      A(2)=16
      A(3)=27
      A(4)=1
      Z(1)=6
      Z(2)=8
      Z(3)=13
      Z(4)=1
      dens=1.5
      wmat(1)=0.4
      wmat(2)=0.3
      wmat(3)=0.25
      wmat(4)=0.05
      call GSMIXT(imate,'scint_wrap',A,Z,dens,4,wmat)
      end

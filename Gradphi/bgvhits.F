c This routine saves hits into JHITS from the BGV detector
c Written 6/7/99 by Tom Bogue
c this rountine assumes that it is only called when in the lead part of BGVs
c
#define EXTRA_DOWNSTREAM_DELAY 0.5e-9

      subroutine bgvhits
      implicit none
#define CERNLIB_TYPE 1
c
c this file defines common variables that describe the particle, specifically
c INWVOL, which is a flag to tell if the particle is entering or exiting the
c current volume
c TOFG, the time of flight of the particle
c ISTOP, a flag that tells if the particle will stop being tracked
c DESTEP, the energy lost this step
c and VECT, the postion/momentum of the particle
#include"geant321/gctrak.inc"
c these two files define common variables that are used in 
c the call to GSCHIT, specifically
c ISET, IDET, ITRA, and IHDET,IHSET
#include"geant321/gcsets.inc"
#include"geant321/gckine.inc"
c
c this file contains the hitparameters for the BGV, 
c namely cEffective_BGV, bgvAttenuationUp, bgvAttenuationDown
#     include "hitconstants.inc"
c this file contains geometrical parameters for the BGV, like
c  LeadLength_BGV
#     include "geometry.inc"
c
c the hit parameters, which should be saved from call to call
      real hit(12)
      save hit
c here is a return variable, not used
      integer IHIT
c hit = (x,y,z,tmin,EtUp,EtDn,ECHA,ENEU,Ein,Elost,Eup,Edn)
c tmin is the time of the first hit in this section
c EtUp is the LaPlace transform of the signal at the upstream end
c  of this section, sampled at some characteristic rise-time value.
c EtDn is the LaPlace transform of the signal at the downstream end
c  of this section, sampled at some characteristic rise-time value.
c ECHA is the energy lost to charged particles
c ENEU is the energy lost to neutral particles
c Eup is the signal after attenuation received at the upstream end
c Edn is the signal after attenuation received at the downstream end
      real TrackNums
      real VolumeNums
      real NumHitsReturned
c these parameters are used to recall the old hit parameters
      real dtUp,dtDn
c these are the time required for the signal to reach the 
c  upstream/downstream lightguides.
      real dEUp,dEDn
      real x(3)
c this is the vector of the hit in the local coordinate sys.
c
c if it's in a new volume,
      if (INWVOL .eq. 1) then
c       initialize the system
c first load the old values, to restore tUp and tDn
        CALL GFHITS(IHSET,IHDET,NVNAME,4,1,ITRA,NUMBV,
     +              TrackNums,VolumeNums,hit,NumHitsReturned)
c IHSET is the name of this set of detectors,
c IHDET is the name of this detector,
c NVNAME is the number of volume numbers needed,
c 6 is the number of hit parameters to restore
c 1 is the number of hits to be returned
c ITRA is the current track number
c NUMBV is the set of volume numbers needed to identify this detector
c TrackNums is the track number of the hits loaded
c VolumeNums is just a buffer to put the volumes
c hit is where to put the hit
c NumHitsReturned is the number of hits returned.  since this is a calorimeter,
c  and uses GSCHIT, there should only be 1 hit.
c   
        hit(1)=VECT(1)
        hit(2)=VECT(2)
        hit(3)=VECT(3)
c xin=x;yin=y;zin=z
c note: VECT, and now hit, is in global coordinates
        if (NumHitsReturned .eq. 0) then
          hit(4)=TOFG
        endif
c tin=t
        hit(5)=0
c EtUp=0, initialize the LaPlace transform of the waveform
        hit(6)=0
c EtDn=0, initialize the LaPlace transform of the waveform
        hit(7)=0
c ECHA=0, since no energy has yet been lost to a charged particle
        hit(8)=0
c ENEU=0, since no energy has been lost to neutrals either
        hit (9)=GETOT
c Ein=ETOT, since the total energy of the particle has been added to the system
        hit(10)=0
c Eout=0 since no energy from this particle has left yet.
        hit(11)=0
        hit(12)=0
c no energy received at either end yet, so Eup=Edn=0
      endif
c if this is a charged particle, 
      if (CHARGE .ne. 0) then
c update tUp and tDn
c the first step to this is to find dtUp and dtDn, 
c  the time required for the signal to reach the upstream and downstream ends
c Now to convert from the MRS to the local cood. sys.
        CALL GMTOD(VECT,x,1)
c we want to add the effective energy at either end
c we assume the signal will attenuate exponentially
        dEup=DESTEP*exp(-x(3)/bgvAttenuationUp)
        dEdn=DESTEP*exp(-(LeadLength_BGV-x(3))/
     +                       bgvAttenuationDown)
c  Now we just calculate the time for the signal to propagate
c  remember, the origin of this sys is the upstream end.
        dtUp=x(3)/cEffectiveBGV
        dtDn=(LeadLength_BGV-x(3))/cEffectiveBGV
c  We keep track of two energy sums: one unweighted and one weighted
c  by the factor exp(-t/tau) where tau is some time scale of order of
c  the signal rise time.  These two numbers approximate the shape of
c  the signal by measuring its LaPlace transform at the origin and
c  at one other point determined by tau=bgvRisetime.
        hit(11)=hit(11)+dEup
        hit( 5)=hit( 5)+dEup*exp(-(TOFG+dtUp)/bgvRisetime)
        hit(12)=hit(12)+dEdn
        hit( 6)=hit( 6)+dEdn*exp(-(TOFG+dtDn)/bgvRisetime)
c now we test to see if these arrived sooner than the earlier hits
c and are over discriminator threshold
        if (TOFG+dtUp .lt. hit(4)) then
          hit(4)=TOFG+dtUp
        endif
      endif
c now to deal with the case of an old track
      if (INWVOL .ne. 1) then
c i.e. this isn't a new track
c if its a charged particle, sum its energy loss to ECHA
        if (CHARGE .ne. 0) hit(7)=hit(7)+DESTEP
c and check to see if this is the last step for this particle 
        if ((INWVOL .eq. 2) .or. (ISTOP .ne. 0)) then 
c if this is a neutral particle which is stopped by the BGV,
c then add its energy into the E from neutrals, ENEU
          if ((CHARGE .eq. 0) .AND. (ISTOP .eq. 2)) 
     +       hit(8)=DESTEP
c now add the rest of the energy into Eout
          hit(10)=GETOT
c now save the hit into JHITS
          if (hit(11)+hit(12) .gt. 0.0) then
            call GSCHIT(ISET,IDET,ITRA,NUMBV,hit,8,IHIT)
          endif
        endif
      endif
      end

      function extract_bgv_t0(hits,end)
      implicit none
      real extract_bgv_t0
      real hits(12)
      integer end
#     include "hitconstants.inc"
c
      extract_bgv_t0=1.0e-6
      if (end .eq. 0) then	!upstream
        if ((hits(5) .gt. 0) .and.
     +      (hits(11) .ge. BGVupstreamThreshold)) then
          extract_bgv_t0=bgvRisetime*log(hits(11)/hits(5))
        endif
      elseif (end .eq. 1) then	!downstream
        if ((hits(6) .gt. 0) .and.
     +      (hits(12) .ge. BGVdownstreamThreshold)) then
          extract_bgv_t0=bgvRisetime*log(hits(12)/hits(6))
     +                   +EXTRA_DOWNSTREAM_DELAY
        endif
      endif
      end

c--------------------------------------------------------------
c This file contains some constants such as effective
c speed of light inside the scintillator, index of refraction in lg block
c etc.
c----------------------------------------------------------------
c
c Mar 26, 2005 -rtj
c   Moved parameters that are only used at digitization time into
c   c header file digitizer.h
c Nov 18, 2002 -rtj
c   Added new array cereN0 that contains the N0 parameter as a
c    function of track angle with respect to the block axis, required
c    for new functionality in lgdhits.F to work.
c Jan 17, 2001 -rtj
c   Modified the LGD refractive index to its correct value of 1.62
c   Modified the lgdAttenuation length to about 170 cm based on
c    spectral transmission data from S. Teige.
c   Modified the "N0" and "convert" parameter values to restore the
c    approximately correct values for mean and r.m.s. LGD responses.
c
c Updated 6/99 by Tom Bogue
c Author: Andriy Kurylov
c
      real lgdAttenuation
      parameter (lgdAttenuation=0.07) !same as above for lgd blocks
c
c The above lgdAttenuation factor is obsolete, being replaced by
c a double exponential formula governed by the following constants
c
      real lgdtaillen,lgdtail,lgdriselen,lgdrise,lgdzmiddle
      parameter (lgdtaillen=88,lgdriselen=3.1)
      parameter (lgdtail=1.03,lgdrise=0.00031)
      parameter (lgdzmiddle=-10.3)
c
      real refindex,cerewind
      parameter (refindex=1.62) ! index of refraction for lgd
      parameter (cerewind=1.74) ! Cerenkov spectral window (eV)
c
c The N0 parameter is now obsolete,
c replaced by a angle-dependent function - see below
      real N0
      parameter (N0=2.0) ! number of Cerenkov photons per cm
c
c The number of photoelectrons collected on the photocathode per unit
c charged track length in the shower depends in an important way on the
c angle of the radiating track with respect to the block axis.  In
c particular, there is a dead spot in the collection efficiency where
c the charged track is nearly parallel to the z-axis.  Note that the
c cereN0 factor tabulated below replaces the former N0 parameter above,
c which reflects an average over track angles in the shower.
c
      real cereN0        ! number of Cerenkov photons per cm
      common /Ceregen/cereN0(180)  ! as a function of angle (degrees)
c
c BGV parameters
c effective speed of light in lightguide in cm/s
c note: cEffective is < c/n because the light bounces around
      real cEffectiveBGV
      parameter (cEffectiveBGV=1.80e+10)
c The attenuation length of light in the fibers of the BGV
      real bgvAttenuationUp,bgvAttenuationDown
      parameter (bgvAttenuationUp=130.0,bgvAttenuationDown=75.0)
c The bgvRisetime parameter is used to sample the LaPlace transform
c of the signal in the two ends of the BGV, later used during packing
c to extract an estimate for the leading edge of the pulse.
      real bgvRisetime
      parameter (bgvRisetime=10.0e-9)
c the threshold for the TDC
      real bgvUpstreamThreshold
      parameter (bgvUpstreamThreshold = 0.015)
      real bgvDownstreamThreshold
      parameter (bgvDownstreamThreshold = 0.015)
c
c BSD parmeters
c the effective speed of light due to n and bouncing
      real cEffectiveBSD
      parameter (cEffectiveBSD= 1.4e+10)
c z-component of the effective speed of light
c cSubzEffective = cEffective of BSD/sqrt(1+(R*(phiMax/Zmax))^2)
      real cSubzEffectiveRing1BSD
      parameter (cSubzEffectiveRing1BSD=cEffectiveBSD*0.7)
      real cSubzEffectiveRing2BSD
      parameter (cSubzEffectiveRing2BSD=cEffectiveBSD*0.7)
      real cSubzEffectiveRing3BSD
      parameter (cSubzEffectiveRing3BSD=cEffectiveBSD*1.0)
      real bsdattenuation
      parameter (bsdattenuation=0.027)  !inverse attenuation length for bsd
      real bsdattenuationRing1
      parameter (bsdattenuationRing1=bsdattenuation)
      real bsdattenuationRing2
      parameter (bsdattenuationRing2=bsdattenuation)
      real bsdattenuationRing3
      parameter (bsdattenuationRing3=bsdattenuation)
c 
c CPV hit parameters
c half length's of the scintillating strips
      real halflengthCPSA
      parameter (halflengthCPSA=28.496)
      real halflengthCPSB
      parameter (halflengthCPSB=22.146)
      real halflengthCPSC
      parameter (halflengthCPSC=28.258)
      real halflengthCPSD
      parameter (halflengthCPSD=26.194)
      real halflengthCPSE
      parameter (halflengthCPSE=26.194)
      real cEffectiveCPV,cpvattenuation
      parameter (cEffectiveCPV=2e+10) !effective speed of light in cpv
      parameter (cpvattenuation=0.01)  !inverse attenuation length for cpv

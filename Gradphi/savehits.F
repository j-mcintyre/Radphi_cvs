*CMZ :          30/08/94  14.16.48  by  S.Ravndal
*-- Author : Andriy Kurylov
      subroutine savehits
      implicit none
C.
c This subroutine saves hits; what it saves depends 
c on the sensitive volume
C.
#define CERNLIB_TYPE 1
#include "geant321/gctmed.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcsets.inc"
C.
      CHARACTER*4 CNAMES(20)
      EQUIVALENCE (NAMES(1),CNAMES(1))
C.
c      
c subroutine rpdhits controls rpd; for rpd two things are 
c saved for each hit: adc (amplitude digital convertor) and
c tdc (time digital convertor). RPD works as a scintillator
c
c subroutine upvhits does the same thing for upv as rpdhits
c does for rpd. The same things are saved
c
c subroutine cpvhits does the same thing for cpv as rpdhits 
c does for rpd. The same things are saved
c
c subroutine lgdhits saves hits for lgd. LGD is a calorimeter
c and things that are saved are: tdc and number of Cerenkov 
c photons emitted by the particle. Attenuation is taken 
c into account.
c
c now we just decide in which volume we are and then call 
c the corresponding subroutine
c
c---------------------------------------------------------------------
      if (CNAMES(NLEVEL).eq.'DSEN'.or.CNAMES(NLEVEL).eq.'THSN') then
       call rpdhits
      elseif ((CNAMES(NLEVEL).eq.'CPSA').or.
     +        (CNAMES(NLEVEL).eq.'CPSB').or.
     +        (CNAMES(NLEVEL).eq.'CPSC').or.
     +        (CNAMES(NLEVEL).eq.'CPSD').or.
     +        (CNAMES(NLEVEL).eq.'CPSE')) then
       call cpvhits
      elseif (CNAMES(NLEVEL).eq.'UPVS') then
       call upvhits
      elseif ((CNAMES(NLEVEL).eq.'BLOK').or.
     +        (CNAMES(NLEVEL).eq.'PMTC')) then
       call lgdhits
      elseif (CNAMES(NLEVEL).eq.'BGVS') then 
       call bgvhits
      elseif (CNAMES(NLEVEL).eq.'BSR1' .or. CNAMES(NLEVEL).eq.'BSR2'
     +    .or.CNAMES(NLEVEL).eq.'BSR3') then
       call bsdhits
      endif
      end

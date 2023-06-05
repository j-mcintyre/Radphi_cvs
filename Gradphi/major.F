c----------------------------------------------------------------------------
c calls subroutines for that define materials, media and geometry
c Author: Andriy Kurylov
c         
      subroutine major   
      implicit none
      call GMATE
      call matdef
      call meddef
      call detector
      call lgd
      call lgpmt
      call rpd
      call cpv
      call upv
      call beampipe
      call beamtube
      call clasmgnt
      call clastgt
      call collim
      call herccover
      call radiator
      call leadwall
      call shield
      call sweepmgnt
      call dnut
      call trgt
      call lgdframe
      call platform
      call bsd
      call bgv
c now placing things
      call placelgd
c     call placerpd
      call placecpv
      call placeupv
      call placetrgt      
      call placebrad
      call placeepip
      call placeclas
      call placehole
      call placegtub
      call placewall
      call placecols
c     call placeweep
      call placednut
      call placecover
      call placeframe
      call placeplatform
      call placebsd
      call placebgv
      call sensitive
      call GGCLOS
      call GPART
      call GPHYSI
      end

c November 18, 2002
c   Special code to store and save an ntuple containing information
c   about hits in the BSD and CPV scintillator layers.  This code
c   was used to estimate the charged particle rates in the scintillators
c   for the Y2000 run.  To enable it, Gradphi must be compiled with the
c   SAVE_RATES_NTUPLE switch turned on.
c - Richard Jones
c
      subroutine declare_hits
      implicit none
      character*4 chset,chdet,chnmsv(1)
      integer nv,nbitsv(1)
      integer idtyp,nwhi,nwdi
      integer iset,idet
      integer nh,nbitsh(10)
      character*4 chnamh(10)
      real orig(10),fact(10)
      data orig/3*1e6,7*0/
      data fact/3*1e3,1e12,2*1,4*1e6/
      data nbitsh/10*32/
      common /hitsetID/iset
      save /hitsetID/

      chset='DIAG'
      chdet='UPVS'
      nv=1
      chnmsv(1)='UPVS'
      nbitsv(1)=32
      idtyp=z'face'
      nwhi=1000
      nwdi=10
      call gsdet(chset,chdet,nv,chnmsv,nbitsv,idtyp,nwhi,nwdi,iset,idet)
      nh=10
      chnamh(1)='x   '
      chnamh(2)='y   '
      chnamh(3)='z   '
      chnamh(4)='t   '
      chnamh(5)='Nnup'
      chnamh(6)='Ncup'
      chnamh(7)='Enup'
      chnamh(8)='Ecup'
      chnamh(9)='Xnup'
      chnamh(10)='Xcup'
      call gsdeth(chset,chdet,nh,chnamh,nbitsh,orig,fact)

      chset='DIAG'
      chdet='CPSA'
      nv=1
      chnmsv(1)='CPVA'
      nbitsv(1)=32
      idtyp=z'face'
      nwhi=1000
      nwdi=10
      call gsdet(chset,chdet,nv,chnmsv,nbitsv,idtyp,nwhi,nwdi,iset,idet)
      nh=10
      chnamh(1)='x   '
      chnamh(2)='y   '
      chnamh(3)='z   '
      chnamh(4)='t   '
      chnamh(5)='Nncv'
      chnamh(6)='Nccv'
      chnamh(7)='Encv'
      chnamh(8)='Eccv'
      chnamh(9)='Xncv'
      chnamh(10)='Xccv'
      call gsdeth(chset,chdet,nh,chnamh,nbitsh,orig,fact)
      chset='DIAG'
      chdet='CPSB'
      nv=1
      chnmsv(1)='CPVB'
      nbitsv(1)=32
      idtyp=z'face'
      nwhi=1000
      nwdi=10
      call gsdet(chset,chdet,nv,chnmsv,nbitsv,idtyp,nwhi,nwdi,iset,idet)
      nh=10
      chnamh(1)='x   '
      chnamh(2)='y   '
      chnamh(3)='z   '
      chnamh(4)='t   '
      chnamh(5)='Nncv'
      chnamh(6)='Nccv'
      chnamh(7)='Encv'
      chnamh(8)='Eccv'
      chnamh(9)='Xncv'
      chnamh(10)='Xccv'
      call gsdeth(chset,chdet,nh,chnamh,nbitsh,orig,fact)
      chset='DIAG'
      chdet='CPSC'
      nv=1
      chnmsv(1)='CPVC'
      nbitsv(1)=32
      idtyp=z'face'
      nwhi=1000
      nwdi=10
      call gsdet(chset,chdet,nv,chnmsv,nbitsv,idtyp,nwhi,nwdi,iset,idet)
      nh=10
      chnamh(1)='x   '
      chnamh(2)='y   '
      chnamh(3)='z   '
      chnamh(4)='t   '
      chnamh(5)='Nncv'
      chnamh(6)='Nccv'
      chnamh(7)='Encv'
      chnamh(8)='Eccv'
      chnamh(9)='Xncv'
      chnamh(10)='Xccv'
      call gsdeth(chset,chdet,nh,chnamh,nbitsh,orig,fact)
      chset='DIAG'
      chdet='CPSD'
      nv=1
      chnmsv(1)='CPVD'
      nbitsv(1)=32
      idtyp=z'face'
      nwhi=1000
      nwdi=10
      call gsdet(chset,chdet,nv,chnmsv,nbitsv,idtyp,nwhi,nwdi,iset,idet)
      nh=10
      chnamh(1)='x   '
      chnamh(2)='y   '
      chnamh(3)='z   '
      chnamh(4)='t   '
      chnamh(5)='Nncv'
      chnamh(6)='Nccv'
      chnamh(7)='Encv'
      chnamh(8)='Eccv'
      chnamh(9)='Xncv'
      chnamh(10)='Xccv'
      call gsdeth(chset,chdet,nh,chnamh,nbitsh,orig,fact)
      chset='DIAG'
      chdet='CPSE'
      nv=1
      chnmsv(1)='CPVE'
      nbitsv(1)=32
      idtyp=z'face'
      nwhi=1000
      nwdi=10
      call gsdet(chset,chdet,nv,chnmsv,nbitsv,idtyp,nwhi,nwdi,iset,idet)
      nh=10
      chnamh(1)='x   '
      chnamh(2)='y   '
      chnamh(3)='z   '
      chnamh(4)='t   '
      chnamh(5)='Nncv'
      chnamh(6)='Nccv'
      chnamh(7)='Encv'
      chnamh(8)='Eccv'
      chnamh(9)='Xncv'
      chnamh(10)='Xccv'
      call gsdeth(chset,chdet,nh,chnamh,nbitsh,orig,fact)

      chset='DIAG'
      chdet='BSR1'
      nv=1
      chnmsv(1)='BSR1'
      nbitsv(1)=32
      idtyp=z'face'
      nwhi=1000
      nwdi=10
      call gsdet(chset,chdet,nv,chnmsv,nbitsv,idtyp,nwhi,nwdi,iset,idet)
      nh=10
      chnamh(1)='x   '
      chnamh(2)='y   '
      chnamh(3)='z   '
      chnamh(4)='t   '
      chnamh(5)='Nnbs'
      chnamh(6)='Ncbs'
      chnamh(7)='Enbs'
      chnamh(8)='Ecbs'
      chnamh(9)='Xnbs'
      chnamh(10)='Xcbs'
      call gsdeth(chset,chdet,nh,chnamh,nbitsh,orig,fact)
      chset='DIAG'
      chdet='BSR2'
      nv=1
      chnmsv(1)='BSR2'
      nbitsv(1)=32
      idtyp=z'face'
      nwhi=1000
      nwdi=10
      call gsdet(chset,chdet,nv,chnmsv,nbitsv,idtyp,nwhi,nwdi,iset,idet)
      nh=10
      chnamh(1)='x   '
      chnamh(2)='y   '
      chnamh(3)='z   '
      chnamh(4)='t   '
      chnamh(5)='Nnbs'
      chnamh(6)='Ncbs'
      chnamh(7)='Enbs'
      chnamh(8)='Ecbs'
      chnamh(9)='Xnbs'
      chnamh(10)='Xcbs'
      call gsdeth(chset,chdet,nh,chnamh,nbitsh,orig,fact)
      chset='DIAG'
      chdet='BSR3'
      nv=1
      chnmsv(1)='BSR3'
      nbitsv(1)=32
      idtyp=z'face'
      nwhi=1000
      nwdi=10
      chnamh(1)='x   '
      chnamh(2)='y   '
      chnamh(3)='z   '
      chnamh(4)='t   '
      chnamh(5)='Nnbs'
      chnamh(6)='Ncbs'
      chnamh(7)='Enbs'
      chnamh(8)='Ecbs'
      chnamh(9)='Xnbs'
      chnamh(10)='Xcbs'
      call gsdet(chset,chdet,nv,chnmsv,nbitsv,idtyp,nwhi,nwdi,iset,idet)
      nh=10
      call gsdeth(chset,chdet,nh,chnamh,nbitsh,orig,fact)
      end


      subroutine register_hit

#include "geant321/pilot.h"
#include "geant321/gcsets.inc"
#include "geant321/gctrak.inc"
#include "geant321/gckine.inc"

      integer ihit
      integer nhsum
      real hits(10)
      integer i
      common /hitsetID/myset
      integer myset
  
      if (iset.ne.myset) return

      do i=5,10
        hits(i)=0
      enddo
      hits(1)=vect(1)
      hits(2)=vect(2)
      hits(3)=vect(3)
      hits(4)=tofg
      if (charge.eq.0) then
        if (inwvol.eq.1) then
          hits(5)=1
          hits(7)=gekin
      write(6,*) 'neutral entry',gekin
        elseif (inwvol.eq.2) then
          hits(9)=gekin
      write(6,*) 'neutral exit',gekin
        endif
      else
        if (inwvol.eq.1) then
          hits(6)=1
          hits(8)=gekin
      write(6,*) 'charged entry',gekin
        elseif (inwvol.eq.2) then
          hits(10)=gekin
      write(6,*) 'charged exit',gekin
        endif
      endif
      nhsum=6
      write(6,*) '*',hits
      call gschit(iset,idet,itra,numbv,hits,nhsum,ihit)
      write(6,*) '+',hits,ihit
      end

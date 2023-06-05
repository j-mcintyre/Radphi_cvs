      subroutine wc
      end

      subroutine wcprint(n,x,y,z)
      integer n
      real x(*),y(*),z(*)
c gcvolu.inc
      COMMON/GCVOLU/NLEVEL,NAMES(15),NUMBER(15),
     +LVOLUM(15),LINDEX(15),INFROM,NLEVMX,NLDEV(15),LINMX(15),
     +GTRAN(3,15),GRMAT(10,15),GONLY(15),GLX(3)
c gtvolu.inc
      INTEGER NLEVEL,NAMES,NUMBER,LVOLUM,LINDEX,INFROM,NLEVMX,
     +        NLDEV,LINMX
      REAL GTRAN,GRMAT,GONLY,GLX
c      
      real xc(3),ubuf(99)
      character*20 natmed
      do i=1,n
        xc(1)=x(i)
        xc(2)=y(i)
        xc(3)=z(i)
        call gmedia(xc,numed)
        call gftmed(numed,natmed,nmat,isvol,ifield,fieldm,tmaxfd,
     +              stemax,deemax,epsil,stmin,ubuf,nwbuf)
        write(6,1000) i,xc,natmed,(names(l),number(l),l=1,nlevel)
      enddo
 1000 format(' point',i3,':',3g12.5,1x,a20,99('/',a4,i3))
      end

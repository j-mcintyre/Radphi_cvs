c this subroutine converts Euler angles into 
c angle convention used by GEANT
      subroutine ANGLES(thetae,phie,psie)
       implicit none
       real thetae,phie,psie
       real thetax,thetay,thetaz,phix,phiy,phiz     
       common thetax,thetay,thetaz,phix,phiy,phiz     
       real rm(3,3), rm0(3,3), rm1(3,3), rm2(3,3), rm3(3,3)
       real tn, sum
       integer i,j,k
c phi rotation      
#ifdef SOLARIS
       rm1(1,1)=COSD(phie)
       rm1(1,2)=SIND(phie)
#else
       rm1(1,1)=DCOS(dble(phie))
       rm1(1,2)=DSIN(dble(phie))
#endif
       rm1(1,3)=0.
       rm1(2,1)=-rm1(1,2)
       rm1(2,2)=rm1(1,1)
       rm1(2,3)=0.
       rm1(3,1)=0.
       rm1(3,2)=0.
       rm1(3,3)=1.
c theta rotation
       rm2(1,1)=1.
       rm2(1,2)=0.
       rm2(1,3)=0.
       rm2(2,1)=0.
#ifdef SOLARIS
       rm2(2,2)=COSD(thetae)
       rm2(2,3)=SIND(thetae)
#else
       rm2(2,2)=DCOS(dble(thetae))
       rm2(2,3)=DSIN(dble(thetae))
#endif
       rm2(3,1)=0.
       rm2(3,2)=-rm2(2,3)
       rm2(3,3)=rm2(2,2)
c psi rotation
#ifdef SOLARIS
       rm3(1,1)=COSD(psie)
       rm3(1,2)=SIND(psie)
#else
       rm3(1,1)=DCOS(dble(psie))
       rm3(1,2)=DSIN(dble(psie))
#endif
       rm3(1,3)=0.
       rm3(2,1)=-rm3(1,2)
       rm3(2,2)=rm3(1,1)
       rm3(2,3)=0.
       rm3(3,1)=0.
       rm3(3,2)=0.
       rm3(3,3)=1.
c compute rotational matrix
       do i=1,3
        do j=1,3
         sum=0.
         do k=1,3
          sum=sum+rm2(i,k)*rm1(k,j)
         enddo
         rm0(i,j)=sum
        enddo
       enddo
       do i=1,3
        do j=1,3
         sum=0.
         do k=1,3
          sum=sum+rm3(i,k)*rm0(k,j)
         enddo
         rm(i,j)=sum
        enddo
       enddo
c m is the rotational matrix
c now compute GEANT angles
c thetax and phix
       thetax=360*ACOS(rm(1,3))/(2*3.1415926)
       if (rm(1,1).ne.0.) then
        tn=rm(1,2)/rm(1,1)
        phix=360*ATAN(tn)/(2*3.1415926)         
       elseif(rm(1,2).ne.0.) then
        tn=rm(1,1)/rm(1,2)
        phix=360*ATAN(tn)/(2*3.1415926)
        phix=90.-phix
       else
        phix=0.
       endif  
c thetay and phiy
       thetay=360*ACOS(rm(2,3))/(2*3.1415926)
       if (rm(2,1).ne.0.) then
        tn=rm(2,2)/rm(2,1)
        phiy=360*ATAN(tn)/(2*3.1415926)         
       elseif(rm(2,2).ne.0.) then
        tn=rm(2,1)/rm(2,2)
        phiy=360*ATAN(tn)/(2*3.1415926)
        phiy=90.-phiy
       else
        phiy=0.
       endif  
c thetax and phiz
       thetaz=360*ACOS(rm(3,3))/(2*3.1415926)
       if (rm(3,1).ne.0.) then
        tn=rm(3,2)/rm(3,1)
        phiz=360*ATAN(tn)/(2*3.1415926)         
       elseif(rm(3,2).ne.0.) then
        tn=rm(3,1)/rm(3,2)
        phiz=360*ATAN(tn)/(2*3.1415926)
        phiz=90.-phiz
       else
        phiz=0.
       endif
       return
      end

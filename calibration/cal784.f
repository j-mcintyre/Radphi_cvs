c
c     routines to calibrate the 71x43 LGD
c
c     calling routine declares storage for the arrays defining
c     the matrix and vector in the system of linear equations
c     defining the calibration constants. (The *'s are 3053)
c
c     There are two sets of routines: gcal0,gcal and gsolv
c     And                           : ecal0,ecal and esolv
c
c     the ecal routines take single particles (electrons or photons)
c     as input. Very much like the calibrator for P2.
c
c     The gcal routines take reconstructed photon events as input.
c     The arguments are similar to the output from the constrained
c     fitting routines.
c
c     DON'T use them at the same time, they will frag each others
c     variables.
c
c     call xcal0 to initialize
c     call xcal  to accumulate data for the linear fit.
c     call xsolv to get the new calibration constants.
c
c     S. Teige 16/Feb/94
c
      subroutine gcal0(ncc,store,x,u)
      implicit none
      real        store(-203:203,*) !packed matrix
      real        x(*)              !energy sums vector
      real        u(*)              !scratch for back substitution
      integer     ncc(*)            !number of events contributing to a block
      integer i,j                   !loop indicies
c
      do i=1,784
         x(i)=0.
         u(i)=0.
         ncc(i)=0
         do j=-203,203
            store(j,i)=0.
         enddo
      enddo
      return
      end
      subroutine gcal(eprime,nbclus,ebclus,mbclus,ipoint,
     $  npair,cc,ncc,store,x,u)
      implicit none
c
c     accumulate the matrix and vector in the equation system for the cc.
c     Changed variables to integrate with constrained fitting routines.
c
c
c     input arguments
c
      real        store(-203:203,*) !packed matrix
      real        x(*)              !energy sums vector
      real        u(*)              !scratch for back substitution
      integer     ncc(*)            !number of events contributing to a block
      real        cc(*)             !current calibration constants
      real        eprime(*)         !fitted cluster energies
      integer     nbclus(8,*)       !pointers to the blocks in the clusters
      real        ebclus(8,*)       !the energies in the blocks
      integer     mbclus(*)         !the number of blocks in the clusters
      integer     ipoint(4,2)       !pointers to photon pairs
      integer     npair             !the number of pairs in this event
c
c     loop indicies etc.
c
      integer     i,j,k,l,m,n,ip
      real        va           !scratch real (Anyone know Bob Handler?)
      real        efit         !  
      real        sigma         !the expected energy resolution
                                !(acts like a convergance accellerator)
c*********************************************************************
      do l=1,npair
         do m=1,2
            ip=ipoint(l,m)
            efit=eprime(ip)
            sigma=(0.06*sqrt(efit)+0.02*efit)**2
            do j=1,mbclus(ip)
               k=nbclus(ip,j)
               ncc(k)=ncc(k)+1
               va=ebclus(ip,j)/cc(k)
               x(k)=x(k)+va*efit/sigma
               do i=1,mbclus(ip)
                  n=nbclus(ip,i)
                  store(k-n,n)=store(k-n,n)
     $                        +(ebclus(ip,i)*va)/(sigma*cc(n))
               enddo
            enddo
         enddo
      enddo
      return
      end
      subroutine gsolv(cc,ccnew,ncc,store,x,u)
      implicit none
c
c     solve the linear system of equations for the
c     new calibration constants.
c
c
c     input arguments
c
      real        store(-203:203,*) !packed matrix
      real        x(*)              !energy sums vector
      real        u(*)              !scratch for back substitution
      integer     ncc(*)            !number of events contributing to a block
      real        cc(*)             !initial values of calibration constants
      real        ccnew(*)          !"improved" values of the cc
      integer     i,j,k             !loop indicies
      real        va,sum            !scratch real
      integer jlim
      real tmp
c
c     force blocks with little or no data
c
      do i=1,784
         if(ncc(i).lt.10) then
            store(0,i)=1.e+20
            x(i)=cc(i)*1.e+20
         else
            print*, 'Channel ',i,' has enough data'
         endif
         if(store(0,i).le.0.1) then
             store(0,i)=1.e+20
             x(i)=cc(i)*1.e+20
         endif
      enddo
c
c     solve the matrix equation
c

      do i=1,783

         jlim = i+203
         if (jlim.gt.784) jlim = 784

         do j=i+1,jlim

C         do j=i+1,i+203
C            if(j.le.784) then

               va=store(j-i,i)/store(0,i)
               if(va.ne.0.) then
                  x(j)=x(j)-va*x(i)      

                  do k=i+1,jlim

C                  do k=i+1,i+203
C                     if(k.le.784) then
                     tmp = store(i-k,k)
                     if(tmp.ne.0) then
                        store(j-k,k)=store(j-k,k)-VA*tmp
C     endif
                     endif
                  enddo
               endif
C            endif
         enddo
      enddo
c
c     solve for u (cc) here
c

      print*,"Negative CCs are listed below:"

      u(784)=x(784)/store(0,784)
      ccnew(784)=u(784)
      do i=783,1,-1
         sum=0.
         do j=1,203
            if((j+i).le.784) then
               sum=sum+store(j,i)*u(j+i)
            endif
         enddo
         u(i)=(x(i)-sum)/store(0,i)

         
c
c        fix for negative cc
c
         if(u(i).lt.0.) then
            print*,"   Channel:",i,", old cc:",cc(i),", new cc:",u(i)
            u(i)=cc(i)
         endif
         ccnew(i)=u(i)
      enddo
      return
      end
c
c for single particle calibration, i.e. electrons.
c
      subroutine ecal0(ncc,store,x,u)
      implicit none
c
c     input arguments
c
      real        store(-203:203,*) !packed matrix
      real        x(*)              !energy sums vector
      real        u(*)              !scratch for back substitution
      integer     ncc(*)            !number of events contributing to a block
c
      integer i,j                   !loop indicies
c
      do i=1,784
         x(i)=0.
         u(i)=0.
         ncc(i)=0
         do j=-203,203
            store(j,i)=0.
         enddo
      enddo
      return
      end
      subroutine ecal(efit,nbclus,ebclus,mbclus,ip,
     $  cc,ncc,store,x,u)
      implicit none
c
c     accumulate the matrix and vector in the equation system for the cc.
c
c
c     input arguments
c
      real        store(-203:203,*) !packed matrix
      real        x(*)              !energy sums vector
      real        u(*)              !scratch for back substitution
      integer     ncc(*)            !number of events contributing to a block
      real        cc(*)             !current calibration constants
      real        efit              !the right answer for the energy
      integer     nbclus(8,*)       !pointers to the blocks in the clusters
      real        ebclus(8,*)       !the energies in the blocks
      integer     mbclus(*)         !the number of blocks in the clusters
      integer     ip                !pointer to the cluster number
                                    !to be used (=1 for electrons)
c
c     loop indicies etc.
c
      integer     i,j,k,n
      real        va                !scratch real (Anyone know Bob Handler?)
      real        sigma             !expected energy resolution
c*********************************************************************
      sigma=(0.06*sqrt(efit)+0.02*efit)**2
      do j=1,mbclus(ip)
         k=nbclus(ip,j)
         ncc(k)=ncc(k)+1
         va=ebclus(ip,j)/cc(k)
         x(k)=x(k)+va*efit/sigma
         do i=1,mbclus(ip)
            n=nbclus(ip,i)
            store(k-n,n)=store(k-n,n)
     $                        +(ebclus(ip,i)*va)/(sigma*cc(n))
         enddo
      enddo
      return
      end
      subroutine esolv(cc,ccnew,ncc,store,x,u)
      implicit none
c
c     solve the linear system of equations for the
c     new calibration constants. Provided for consistancy.
c     (Actually just calls gsolv)
c
      real        store(-203:203,*)     !packed matrix
      real        x(*)           
      real        u(*)
      real        cc(*),ccnew(*)
      integer     ncc(*)
      call gsolv(cc,ccnew,ncc,store,x,u)
      return
      end

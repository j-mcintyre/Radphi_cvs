      program balance
c
c     gain balance the tubes using calibration constants
c
c     SWT 9/Jun/98
c
      implicit none
c
      real cc(1000)
      integer volt(1000)
c
      real ccgoal /0.00126/   ! aim for this value for all cc's
c
      real gain_coef /9.362/ ! exponent in gain vs voltage function
c
      integer i
      integer k,tmpv
      real tmpcc
      integer imax,itenp,itenm
      integer irow,icol
      real delta
c
c     get the current voltages
c
      open(unit=77,type='old',name='new_voltages.dat',
     $     form='formatted')
c
c     get the current calibration constants
c
      open(unit=78,type='old',name='cc.dat',
     $     form='formatted')
c
c     initialize
c
      do i=1,1000
         cc(i)=0.
         volt(i)=0
      enddo
c
      do i=1,1000
         read(77,*,end=100),k,tmpv
         volt(k)=tmpv
      enddo
 100  continue
c
      do i=1,1000
         read(78,*,end=200),tmpcc,k
         cc(k)=tmpcc
      enddo
 200  continue
c
c     start balancing gains.
c
      imax=0
      itenp=0
      itenm=0
      do i=1,800
         if((volt(i).ne.0).and.(cc(i).ne.0.)) then
            tmpv=volt(i)*exp(alog(cc(i)/ccgoal)/gain_coef)
c
c           have calculated the new voltage, check to see if it 
c           is reasonable.
c
            delta=1. - float(tmpv)/float(volt(i))
c
c           increase more than 10% ?
c
            if(delta.ge.0.10) then
               tmpv=1.1*float(volt(i))
               itenp=itenp+1
            endif
c
c           maximum voltage is 1950
c
            if(tmpv.gt.1950) then
               tmpv=1950
               if(volt(i).ne.1950) then
                  imax=imax+1
                  print*,' pushing maximum ',i,volt(i),tmpv
               endif
            endif
c
c           decrease more than 10%
c
            if(delta.le.-0.10) then
               tmpv=0.90*float(volt(i))
               itenm=itenm+1
            endif
c
c           write down the answer
c
            irow=i/28
            icol=i-28*irow
            write(80,22) i,tmpv,volt(i),tmpv-volt(i),irow,icol
 22         format(6(1x,i5))
            volt(i)=tmpv
         endif
      enddo
      print*,' number tending over maximum    ',imax
      print*,' number wanting to increase 10% ',itenp
      print*,' number wanting to decrease 10% ',itenm
      do i=1,800
         if(volt(i).ne.0) then
            write(79,*) i,volt(i)
         endif
      enddo
      stop
      end

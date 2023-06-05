* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
	subroutine psfit() 	! package loader (dummy)
	end

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
*  
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
	subroutine fitter(count)
        implicit none
        integer count
        character*1000 instr
        include 'psfit.inc'
 
        integer ret,icond,space,run,step
        integer eventnext,eventrun

        real ufunc
        double precision chi
        external eventnext,eventrun
        external FCN,ufunc
       
        instr = "par.in"
        step = report_step
 
        call MNSETI('MINIUT usage for psf clusterizer')
c        call MNCOMD(FCN,'show title',icond,ufunc)

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
*  Fit single event 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
        space = index(instr,' ')
        call setuphbook()
        if ( count.eq.0 ) then
          lunout=6
          errout=6
          fiterror=0
          event = eventnext()
          run = eventrun()
          call setpar0(instr(:space-1),lunout)
          call fit1('d',lunout)
          Write(6,*) 'Report for event', event,fiterror,' RunNo',run
          call report(lunout)
          return
        elseif ( count.gt.0 ) then
          lunout=98
          errout=97
          open(unit=lunout,file='fit.out',status='unknown')
          open(unit=errout,file='err.out',status='unknown')
	  do i=1,count 
            event = eventnext()
	    if (event.le.0) return
            run = eventrun()
            fiterror=0
            call setpar0(instr(:space-1),lunout)
            call fit1('d',lunout)
            Write(lunout,*) 'Report for EVENT', event,' COUNT',i,
     *                       ' RUN',run
	    call report(lunout)
            if (fiterror.ne.0) then
              Write(errout,*) 'EVENT', event,' COUNT',i,
     *                       ' RUN',run
              Write(errout,*) '* * * * * * * *'
            endif
            if (i/step*step.eq.i) then
              write(6,*) i
              if (i/step.ge.10) step = step*10 
            endif
	  enddo
          close(lunout)
          close(errout)
          return
        else
          write(6,*) 'The <count> has to be non-negative!'
	endif

c        call MNCOMD(FCN,'SET printout -1',icond,ufunc)
c        call MNCOMD(FCN,'call FCN 0',icond,ufunc)

   	end


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
* User method to skip $count of events in the input file list
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
	 integer function fwd(count)
	 implicit none
         integer count
         include 'psfit.inc'
  
         integer eventnext
         external eventnext

         do i=1,count
           event = eventnext()
           if (event.le.0) then 
             fwd = i
             return 
           endif
         enddo
         write(6,*) 'Skipped ',count, ' events!'
         fwd = count
         return 
         end

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
* User method to plot moments of X and Y hit coordinates 
* within a shower image.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
         integer function moms(count)
         implicit none
         integer count
         include 'psfit.inc'
 
         integer mom,step
         integer eventnext,fitgetmom
         external eventnext,fitgetmom
     
         step=report_step
         call setuphbook()
         do i=1,count
           event = eventnext()
           if (event.le.0) return
           mom = fitgetmom() 
           if (mom.lt.0) then 
             write(6,*) 'Error in moments...',i
             return
           endif

           if (i/step*step.eq.i) then
             write(6,*) i
             if (i/step.ge.10) step = step*10 
           endif
         enddo
         moms = count
         return
         end

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
* User method to initialize parameters for MINUIT
* using valaues from the current event. 
* Parameter names, steps ... are defined in psfit.inc
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
	subroutine setpar0(carfile,out)
	implicit none
        integer out
        character*(*) carfile
        include 'psfit.inc'

        character*20 command(MAX_PAR)
        integer icond
        
        integer  np,nc
        real ufunc
        external FCN,ufunc
       
C  Read additional pars and commands from input parameter file
C  The format is the same as for MINIUT command SET PAR
C  execept that the first stand-alone number determines the number
C  of parameters, while the second stand-alone number determines
C  the number of MINIUT commands.
        lunin=99
        open(unit=lunin,file=carfile,status='old',err=9) 
        read(lunin,*) np
        do i=1,np
          read(lunin,*) j, pname(j),pstart(j),pstep(j),pmin(j),pmax(j) 
        enddo     
        read(lunin,*) nc
        do i=1,nc
          read(lunin,*) command(i)
        enddo     
        close(lunin)

        call MNINIT(5,out,7)

C get init values
        call fitgetpar0(nparf,pstart) 
        npars=np+nparf
        do i=1,npars
          ipar(i) = i
	enddo
       
	do i=1,npars
C For the first parameter (energy) set MINUIT fitting step to 
C the pre-deffined percentage of energy.
          if(i.eq.1) then
            call MNPARM(ipar(i),pname(i),pstart(i),pstep(i)*pstart(i),
     *                  pmin(i),pmax(i),icond)
          else
            call MNPARM(ipar(i),pname(i),pstart(i),pstep(i),
     *                  pmin(i),pmax(i),icond)
          endif
          if (icond.ne.0) then
	    write(6,*) 'Unable to define parameter No.',ipar(i)
            return
	  endif
	enddo

C execute additional commands 
        do i=1,nc
          call MNCOMD(FCN,command(i),icond,ufunc)
          if (icond.ne.0) then
	    write(6,*) 'Unable to execute command',command(i)
            return
          endif
	enddo
        return
9	write(6,*) 'ERROR - could not open file ',carfile
	end

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
* MINUIT user function
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
        subroutine FCN(npar,grad,fval,xval,iflag,func)
        implicit none
        integer npar
        double precision fval,grad(*)
        double precision xval(*)
        integer iflag
        external func
        integer Initialize,FindGradient,Terminate,FindFunction
        parameter (Initialize=1,FindGradient=2,Terminate=3,
     *             FindFunction=4)
        integer Printout,Summary
	parameter (Printout=0,Summary=-1)
* _____________________________________________________________________
* |    MINUIT USER FUNCTION ARGS                                      |
* |                                                                   |
* |     npar    number of variable parameters                         |
* |     grad    vector of derivatives                                 |
* |     fval    function value calculated                             |
* |     xval    external parameter values (must not be changed)       |
* |     iflag   mode flag:                                            |
* |                         1 = initialization call                   |
* |                         2 = gradient call => calculate derivatives|
* |                         3 = termination call                      |
* |                         4 = normal entry => calculate f           |
* |                         5 = reserved                              |
* |                        >5 => user options                         |
* +-------------------------------------------------------------------"
*	Fit parameters 
*		1	Energy	Total shower energy 
*	        2	X0	Shower center coord. in the  LGD plane
*  		3	Y0		-||-
* 		4	Sigma	Sigma for dummy one gauss-fit
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        include 'psfit.inc'

        double precision fmin,fedm,errdef
        integer npari,nparx,istat
       
        double precision xval1(MAX_PAR)
        double precision xstep(MAX_PAR)
        double precision fitgetchi2
        external fitgetchi2

        if (iflag.eq.Initialize ) then
          return        	
	endif

C NOTE: Npars is generaly different then npar
C       npar - number of fitting pars 
C	npars - number of all pars
C       One needs npars to set proper number of parameters for Chi2, 
C       and npar for ndof!

        call MNSTAT(fmin,fedm,errdef,npari,nparx,istat)

	fval = fitgetchi2(nparx,xval,ndof)
        ndof = ndof - npar

        if (iflag.eq.FindFunction) then
          return
        endif

        if (iflag.eq.FindGradient) then
          write(6,*) 'Better ask MINUIT to FindGradient,'
          write(6,*) 'return ...'
          return
          xstep(1) = pstep(1)*xval(1)
          do i=2,npar
            xstep(i) = pstep(i)
          enddo
          do i=1,npar
            do j=1,npar
              xval1(j) = xval(j)
              if (i.eq.j) xval1(j)=xval(j)+xstep(j)
	    enddo
            grad(i) = (fitgetchi2(npar,xval1,ndof) - 
     *               fitgetchi2(npar,xval,ndof))/xstep(i)
          enddo
        endif
  
        if(iflag.eq.Printout) then
        write(6,*) 'Chi2:',fval,'  Ndof:',ndof,'  Chi2/Ndof:',fval/ndof
          write(6,*) 'Parameters:',(xval(i),i=1,nparx)
C          write(6,*) 'Gradient:',(grad(i),i=1,npar)
        endif
        end

* * * * * * * * * * * * * * * * * * * * * * * * * * *
* User-callable method to launch an interactive     *
* MINUIT session                                    *
* * * * * * * * * * * * * * * * * * * * * * * * * * *
	subroutine fit()
        implicit none
        include 'psfit.inc'
	external FCN
        external ufunc
         call MINUIT(FCN,ufunc)
        end

* * * * * * * * * * * * * * * * * * * * * * * * * * *
* User-callable method to launch one fit using      *
* MINUIT 'migrad' procedure			    *
* NOTE - currently carfile is dumy variable since   *
*        fir parameters are initialized from data   *
* * * * * * * * * * * * * * * * * * * * * * * * * * *
        Subroutine fit1(carfile,out)
        Implicit none
        integer out
        Character*(*) carfile
        include 'psfit.inc'
        external FCN
        external ufunc

        double precision fmin,fedm,errdef
        integer npari,nparx,istat
        integer icondn

        call MNCOMD(FCN,'SET printout -1',icondn,ufunc)
        if (out.ne.6) then
          call MNCOMD(FCN,'SET OUT 98',icondn,ufunc)
        endif
        write(out,*) 'first call to migrad...'
        call MNCOMD(FCN,'migrad',icondn,ufunc)
        if (icondn.ne.0) then
          write(out,*) 'migrad failed, retrying...'
          call MNCOMD(FCN,'migrad',icondn,ufunc)
          if (icondn.ne.0) then
            write(errout,*) 'HELP - cannot find minimum, going out'
            write(out,*) 'HELP - cannot find minimum, going out!'
            fiterror = 1
c            write(6,*) 'try "set output 6" to see output on screen,'
c            write(6,*) 'and "set interactive" to see command prompt'
c            write(6,*) 'When satisfied with fit, return with control-D.'
c            call MNINPU(5,istat)
c            call MNINTR(FCN,ufunc)
          endif
        endif

        call MNSTAT(fmin,fedm,errdef,npari,nparx,istat)

        write(out,1000) fmin,ndof,npari,nparx
1000    format('Fit1 found minimum',f15.5,
     +       ' and',i3,' degrees of freedom,',
     +       ' varing',i3,' of' ,i3, 'pars') 

        end

* * * * * * * * * * * * * * * * * * * * * * * * * * *
* User-callable method to collect fitting results 
*
* I tried to fill all histograms within one module (fithf). However, 
* it seems profile histograming does not work in C (might be because hf2 
* is a variant of hfill on which profile cannot be appliied).
* Profile histograms are defined in kumac (outside of fitting module)
* This requires manual histogram reset, but provides filling over
* several fitting loops.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        Subroutine report(out)
        Implicit none
        integer out
        include 'psfit.inc'
        external FCN, ufunc

        integer icondn
        double precision pval(MAX_PAR)
        double precision perr(MAX_PAR)
        double precision plow(MAX_PAR)
        double precision pupr(MAX_PAR)

        real E,X,Y,R1,chi2,chi2ndof
        real Etot, Emax, Edif
        real theta,phi,Z0
        real thetaFit,phiFit
        real sinx,siny
        real sinxf,sinyf

        double precision fmin,fedm,errdef
        integer npari,nparx
        integer id6,id7
        integer fithf, eventmc
        real Ephot,Px,Py,Pz
        real eventphoton,mcphoton
        
        external fithf
        external eventphoton, mcphoton, eventmc
        
        if (out.ne.6) then
          call MNCOMD(FCN,'SET OUT 98',icondn,ufunc)
        endif

        call MNSTAT(fmin,fedm,errdef,npari,nparx,icondn)
        do i=1,nparx
          call MNPOUT(i,pname(i),pval(i),perr(i),
     *                plow(i),pupr(i),icondn)
        enddo  

        Z0=LGD_PLANE
        E=pval(1)
        X=pval(2)
	Y=pval(3)
        R1=pval(4)
        thetaFit = atan2(sqrt(X*X+Y*Y),Z0)        
        phiFit = atan2(Y,X)
        sinxf = sin(thetaFit)*cos(phiFit)
        sinyf = sin(thetaFit)*sin(phiFit)

        chi2=fmin
        chi2ndof = fmin/ndof

        isMC = eventmc()
        if (isMC.gt.0) then
          Ephot = mcphoton(Px,Py,Pz)
          theta = atan2(sqrt(Px*Px+Py*Py),Pz)
          phi = atan2(Py,Px) 
        else
          Ephot = eventphoton(Px,Py,Pz)
          theta = atan2(sqrt(Px*Px+Py*Py),Pz)
          phi = atan2(Py,Px) 
        endif
        sinx = Px/Ephot
        siny = Py/Ephot
    
        Edif = abs((Ephot-E)/Ephot)
        write(out,*) "In: E",Ephot," theta",theta," phi",phi
        write(out,*) "FFit: E",E," theta",thetaFit," phi",phiFit,
     *                 " R",R1,fiterror
        
        if ( abs(sinx-sinxf).gt.0.1.or.Edif.gt.2.or.
     *       abs(siny-sinyf).gt.0.1) then 
          write(out,*) "WARNING - Fitter missed The Photon!"
          write(errout,*) "In: E",Ephot," theta",theta," phi",phi
          write(errout,*) "FFit: E",E," theta",thetaFit," phi",phiFit,
     *                 " R",R1
          fiterror = 1
        endif 

        if (out.eq.6.or.fiterror.eq.0) then
          icondn = fithf(nparx,pval)
          if (icondn.lt.0) then
            write(out,*) 'No histograms from fitImage!'
            return
          endif
          call  hfill(30, chi2, 0., 1.)
          call  hfill(31, float(ndof), 0., 1.)
          call  hfill(32, chi2ndof, 0., 1.)
C -- energy slices
          nbins=9
          do i=1,nbins
            id6 = 6000+i-1
            id7 = 7000+i-1
            if (Ebonds(i).lt.Ephot.and.Ephot.le.Ebonds(i+1)) then
              call  hfill(id6, chi2, 0., 1.)
              call  hfill(id7, chi2ndof, 0., 1.)
            endif 
          enddo
C -- profiles ---

        endif

	end

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* User-callable method to execute a single MINUIT command
* through the function-call interface. The procedure assumes 
* that MINUIT has already been initialized via a call setpar0  
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      Integer Function exec1(command)
      Implicit none
      Character*(*) command
      include 'psfit.inc'
      external FCN
      external ufunc
      integer icondn

C      open(unit=lunout,file='cpvfit.out',status='unknown')
C      open(unit=lunsav,file='cpvfit.sav',status='unknown')
       call MNCOMD(FCN,command,icondn,ufunc)
       if (icondn.ne.0) then
         write(6,*) 'Error from MNCOMD: icondn=',icondn
         write(6,*) ' command was: ',command
       endif
C       close(unit=lunout)
C       close(unit=lunsav)
       exec1=icondn
      end



* * * * * * * * * * * * * * * * * * * * * * * * * * *
* Dummyfroutine needed by MINUIT                  *
* * * * * * * * * * * * * * * * * * * * * * * * * * *
        real function ufunc(x)
        implicit none
        real x
         ufunc = 0.
        end


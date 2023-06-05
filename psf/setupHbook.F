       subroutine setuphbook()
       implicit none
       include 'psfit.inc'
       
       integer id
       logical HEXIST
       external HEXIST

* * * * * * * * * * * * * * * * * * * * 
* Reserved HISTOGRAM ID numbers
* ! BEFORE BOOKING - make sure ID is defined in psfit.inc
* - 20   --- LGD, photons, MC  kinematics ..
* - 30   --- Chi2, ndof ..
* - 200  --- LGD hits 2D
* - 300  --- parameters from fitter
* - 350  --- moments m1,m2,V,S,K(x,y) 
* - 400  --- variables from fitter
* - 500  --- profiles from fitter (currently booked from psf.kumac)!
* ......
* - 10xx --- LGD yield in the xx E-slice 
* - 11xx --- LGD hits in the xx E-slice 
* - 20yy --- PSF        
* - 30xx --- Errors
* - 60xx --- Chi2 
* - 70xx --- Chi2/nfod 
* - 80xx --- LGD yield bin by bin 
* - 81xx --- LGD hits  -- || --
* - 90xx --- PSF       -- || -- 
* * * * * * * * * * * * * * * * * * * * *

       write(6,*) 'Number of histograms', hnum
       do i=1,hnum
         id = hid(i)
         if (id.ne.0) then 
         if (.not.HEXIST(id)) then					! hid
	   if (i.eq.1) call hbook1(id,'Etot lgd',150,0.,6.,0)		! 20
	   if (i.eq.2) call hbook1(id,'Emax lgd',100,0.,4.,0)		! 21	
	   if (i.eq.3) call hbook1(id,'Etot img',150,0.,6.,0)		! 22
	   if (i.eq.4) call hbook1(id,'Emax img',100,0.,4.,0)		! 23
	   if (i.eq.5) call hbook1(id,'Ephot',150,0.,6.,0)		! 24
	   if (i.eq.6) call hbook1(id,'Emc',150,0.,6.,0)		! 25
	   if (i.eq.7) call hbook1(id,'Theta Phot',100,0.,1.,0)		! 26
	   if (i.eq.8) call hbook1(id,'Phi Phot',200,-1.,1.,0)		! 27
	   if (i.eq.9) call hbook1(id,'Theta MC',100,0.,1.,0)		! 28
	   if (i.eq.10) call hbook1(id,'Phi MC',200,-1.,1.,0)		! 29


	   if (i.eq.11) call hbook1(id,'Chi2',500,0.,500.,0)		! 30
	   if (i.eq.12) call hbook1(id,'Ndof',100,-0.5,99.5,0)		! 31
	   if (i.eq.13) call hbook1(id,'Chi2/Ndof ',500,0.,50.,0)	! 32
	   if (i.eq.14) call hbook1(id,'E bins',300,0.,6.,0)		! 32

	   if (i.eq.21) call hbook2(id,'LGD hits',28,-0.5,27.5,		! 200
     * 					          28,-0.5,27.5,0)
	   if (i.eq.22) call hbook2(id,'psf hits',28,-0.5,27.5,		! 201
     * 					          28,-0.5,27.5,0)
	   if (i.eq.23) call hbook2(id,'Errors',28,-0.5,27.5,		! 202
     * 					          28,-0.5,27.5,0)
	   if (i.eq.24) call hbook1(id,'Data bins',49,-0.5,48.5,0)	! 208
	   if (i.eq.25) call hbook1(id,'PSF  bins',49,-0.5,48.5,0)	! 209


	   if (i.eq.31) call hbook1(id,'E fit',150,0.,6.,0)		! 300
	   if (i.eq.32) call hbook1(id,'X0 fit',120,-60.,60.,0)		! 301
	   if (i.eq.33) call hbook1(id,'Y0 fit',120,-60.,60.,0)		! 302
	   if (i.eq.34) call hbook1(id,'R fit',100,0.,10.,0)		! 303
	   if (i.eq.37) call hbook1(id,'Theta fit',100,0.,1.,0)		! 306
	   if (i.eq.38) call hbook1(id,'Phi fit ',200,-1.,1.,0)		! 307

	   if (i.eq.41) call hbook1(id,'Avg X',120,-60.,60.,0)		! 350
	   if (i.eq.42) call hbook1(id,'Rms X',100,0.,10.,0)	! 351
	   if (i.eq.43) call hbook1(id,'Skew X',200,-10.,10.,0)		! 352
	   if (i.eq.44) call hbook1(id,'Kurtosis Y',400,-10.,30.,0)	! 353
	   if (i.eq.45) call hbook1(id,'Avg Y',120,-60.,60.,0)		! 360
	   if (i.eq.46) call hbook1(id,'Rms Y',100,0.,10.,0)	! 361
	   if (i.eq.47) call hbook1(id,'Skew Y',200,-10.,10.,0)		! 362
	   if (i.eq.48) call hbook1(id,'Kurtosis Y',400,-10.,30.,0)	! 363

	   if (i.eq.51) call hbook2(id,'Efit vs Ephot',			! 400
     *				   60,0.,6.,60,0.,6.,0)		
	   if (i.eq.52) call hbook2(id,'R vs Ephot',			! 401
     *				   60,0.,6.,50,0.,10.,0)		
	   if (i.eq.53) call hbook2(id,'TheFit vs Ephot',		! 402
     *				   60,0.,6.,20,0.,1.,0)		
	   if (i.eq.54) call hbook2(id,'PhiFit vs Ephot',		! 403
     *				   60,0.,6.,40,-1.,1.,0)		
	   if (i.eq.55) call hbook2(id,'Efit vs Emc',			! 404
     *				   60,0.,6.,60,0.,6.,0)		
	   if (i.eq.56) call hbook2(id,'R vs Emc',			! 405
     *				   60,0.,6.,50,0.,10.,0)		
	   if (i.eq.57) call hbook2(id,'TheFit vs Emc',			! 406
     *				   60,0.,6.,20,0.,1.,0)		
	   if (i.eq.58) call hbook2(id,'PhiFit vs Emc',			! 407
     *				   60,0.,6.,40,-1.,1.,0)		
	   if (i.eq.59) call hbook2(id,'Etot vs Emc',			! 408
     *				   60,0.,6.,60,0.,6.,0)		

	   if (i.eq.61) call hbprof(id,'Efit vs Ephot',			! 500
     *				   60,0.,6.,0.,10.,'e')		
	   if (i.eq.62) call hbprof(id,'R vs Ephot',			! 501
     *				   60,0.,6.,0.,10.,'e')		
	   if (i.eq.63) call hbprof(id,'TheFit vs Ephot',		! 502
     *				   60,0.,6.,0.,2.,'e')		
	   if (i.eq.64) call hbprof(id,'PhiFit vs Ephot',		! 503
     *				   60,0.,6.,-2.,2.,'e')		
	   if (i.eq.65) call hbprof(id,'Efit vs Emc',			! 504
     *				   60,0.,6.,0.,10.,'e')		
	   if (i.eq.66) call hbprof(id,'R vs Emc',			! 505
     *				   60,0.,6.,0.,10.,'e')		
	   if (i.eq.67) call hbprof(id,'TheFit vs Emc',			! 506
     *				   60,0.,6.,0.,2.,'e')		
	   if (i.eq.68) call hbprof(id,'PhiFit vs Emc',			! 507
     *				   60,0.,6.,-2.,2.,'e')		
	   if (i.eq.69) call hbprof(id,'Etot vs Emc',			! 508
     *				   60,0.,6.,0.,10.,'e')		
	   if (i.eq.70) call hbprof(id,'Ephot vs Emc',			! 509
     *				   60,0.,6.,0.,10.,'e')		

 
           if (id.ge.1000.and.id.lt.2000) call hbook2(id,'LGD hits',	!1[01]xx
     *                                    28,-0.5,27.5,28,-0.5,27.5,0)
           if (id.ge.2000.and.id.lt.3000) call hbook2(id,'PSF hits',	! 20xx
     *                                    28,-0.5,27.5,28,-0.5,27.5,0)
           if (id.ge.3000.and.id.lt.4000) call hbook2(id,'Errors', 	! 30xx
     *                                    28,-0.5,27.5,28,-0.5,27.5,0)

           if (id.ge.6000.and.id.lt.7000) call hbook1(id,'Chi2 ', 	! 60xx
     *                                      500,0.,500.,0)
           if (id.ge.7000.and.id.lt.8000) call hbook1(id,'Chi2/ndof',   ! 70xx
     *                                      500,0.,50.,0)
           if (id.ge.8000.and.id.lt.9000) call hbook1(id,'LGD bins',    !8[01]xx
     *                                      49,-0.5,48.5,0)
           if (id.ge.9000.and.id.lt.10000) call hbook1(id,'PSF bins',   ! 90xx
     *                                      49,-0.5,48.5,0)
	 else
	   call hreset(id,'')
	 endif
         endif
       enddo
       end


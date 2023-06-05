	integer MAX_HIST, hnum  	 	! Maximum and current number
        parameter (MAX_HIST=10,hnum=3)		! of global histograms 
        integer hid(MAX_HIST)
        COMMON /idhist/ hid


	integer MAX_FILES 
        parameter (MAX_FILES=100)     
   
	integer NDIMC,SHOWER_SIZE		! global parameters defined in 
	parameter (NDIMC=784,SHOWER_SIZE=9)	! psf.h also

C --  Event data ---
	integer nhits,nhitsch(NDIMC)
	real etotal,etotalch(NDIMC)
	COMMON /cevent/ nhits,etotal,nhitsch,etotalch

C -- Shower data ---
	integer crow,ccol
	real xc,yc,emax,etot,yield(SHOWER_SIZE,SHOWER_SIZE)
	COMMON /cshower/ crow,ccol,xc,yc,emax,etot,yield


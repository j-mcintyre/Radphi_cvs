* * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Histogram setup
* files: psfit.inc, setupHbook.F
* variables: 	- hnum		curent number of histograms
* 
* * * * * * * * * * * * * * * * * * * * * * * * * * * *
	integer MAX_HIST, hnum  	 	! Maximum and current number
        parameter (MAX_HIST=200,hnum=170)	! of global histograms 
        integer hid(MAX_HIST)

        integer i,j,isMC

        integer report_step
        parameter (report_step=100)	! progres report step  


* * * * * * * * * * * * * * * * * * * * 
* Reserved HISTOGRAM ID numbers
* ! BEFORE BOOKING - make shure ID is defined in psfit.inc
* - 20   --- LGD, photons, MC kinematics ..
* - 30   --- Chi2, ndof ...
* - 200  --- LGD hits 2D
* - 400  --- 2D plots from fitter
* - 300  --- parameters from fitter
* - 350  --- moments
* - 500  --- profiles from fitter
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

      DATA (hid(i),i=1,hnum)/20, 21, 22, 23, 24, 25, 26, 27, 28, 29, !  1 - 10 
     *                       30, 31, 32, 33,  0,  0,  0,  0,  0,  0, ! 11 - 20
     *                      200,201,202,208,209,  0,  0,  0,  0,  0, ! 21 - 30
     *                      300,301,302,303,304,305,306,307,  0,  0, ! 31 - 40
     *                      350,351,352,353,360,361,362,363,  0,  0, ! 41 - 50
     *                      400,401,402,403,404,405,406,407,408,  0, ! 51 - 60
     *                      500,501,502,503,504,505,506,507,508,509, ! 61 - 70
     *                      510,511,  0,  0,514,515,  0,  0,  0,  0, ! 71 - 80
     *             1000,1001,1002,1003,1004,1005,1006,1007,1008,  0, ! 81 - 90
     *             1100,1101,1102,1103,1104,1105,1106,1107,1108,  0, ! 91 -100
     *             2000,2001,2002,2003,2004,2005,2006,2007,2008,  0, !101 -110
     *             3000,3001,3002,3003,3004,3005,3006,3007,3008,  0, !111 -120
     *             6000,6001,6002,6003,6004,6005,6006,6007,6008,  0, !121 -130
     *             7000,7001,7002,7003,7004,7005,7006,7007,7008,  0, !131 -140
     *             8000,8001,8002,8003,8004,8005,8006,8007,8008,  0, !141 -150
     *             8100,8101,8102,8103,8104,8105,8106,8107,8108,  0, !151 -160
     *             9000,9001,9002,9003,9004,9005,9006,9007,9008,  0/ !161 -170

      DATA (hid(i),i=hnum+1,MAX_HIST) /30*0/

       integer lunin,lunout,lunsav,fiterror,errout
       COMMON /fiterr/ fiterror,errout

       integer MAX_E_BIN, nbins
       parameter (MAX_E_BIN=20)
       real Ebonds(MAX_E_BIN)

       DATA Ebonds/  0.0, 0.2, 0.4, 0.8, 1.2,
     *               1.6, 2.0, 2.8, 3.6, 6.0,
     *               0.0, 0.0, 0.0, 0.0, 0.0,
     *               0.0, 0.0, 0.0, 0.0, 0.0/ 

c       COMMON /eslices/ nbins,Ebonds


        real LGD_PLANE
        parameter (LGD_PLANE=120.)

* * * * * * * * * * * * * * * * * * * * * * * * * * * *
* nparf = internal, npars = external pars
* - nparf is the number of final fit parameters (E,X,Y)
* - npars is total number of params supplied to MINIUT
* * * * * * * * * * * * * * * * * * * * * * * * * * * *
        integer MAX_PAR,lgdrows,lgdcols
        integer nparf,npars 	

        parameter (MAX_PAR=50,nparf=3,lgdrows=28,lgdcols=28)

        integer ipar(MAX_PAR)
        CHARACTER*8 pname(MAX_PAR)
        double precision pstart(MAX_PAR)  
        double precision pstep(MAX_PAR)   
        double precision pmin(MAX_PAR)  
        double precision pmax(MAX_PAR)

* * * * * * * * * * * * * * * * * * * * * * * * * *
* Init values for the first 3 parameters: names, limits ...
* are going to be taken from data. 
* The step value in the case of Energy refes to percentage of Etot.
* * * * * * * * * * * * * * * * * * * * * * * * * *
        DATA (pname(i),i=1,3) /' Energy', ' X0', ' Y0'/
        DATA (pstep(i),i=1,3) /0.01d0,0.2d0,0.2d0/
        DATA (pmin(i),i=1,3) /0.d0,2*0.d0/
        DATA (pmax(i),i=1,3) /100.d0,2*0.d0/

	integer ndof,event,events
	COMMON /fitdata/ ndof


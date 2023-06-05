* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
* This is a fortran skeleton - dummy package loader
* 			     - acctual cod is in *.F file 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
	subroutine psfit() 	! package loader (dummy)
	implicit none
         write(6,*) 'This should not be called'
	end

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
* 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
	subroutine fitter(count)
        implicit none
        integer count
         write(6,*) 'This should not be called'
   	end

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
* 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
        integer function fwd(count)
	implicit none
        integer count
         write(6,*) 'This should not be called'
	end

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        integer function moms(count)
        implicit none
        integer count
         write(6,*) 'This should not be called'
        end

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
* 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
	subroutine setpar0(carfile,out)
	implicit none
        integer out
        character*(*) carfile
         write(6,*) 'This should not be called'
   	end
	 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
* 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
        Subroutine report(out)
	implicit none
        integer out
         write(6,*) 'This should not be called'
	end

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
* 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
        subroutine FCN(npar,grad,fval,xval,iflag,func)
        implicit none
        integer npar
        double precision fval,grad(0:*)
        double precision xval(0:*)
        real func
        integer iflag
        external func
        integer Initialize,FindGradient,Terminate,FindFunction
        parameter (Initialize=1,FindGradient=2,Terminate=3,
     *            FindFunction=4)
        integer Printout,Summary
	parameter (Printout=0,Summary=-1)
* _____________________________________________________________________
* |    MINUIT USER FUNCTION                                           |
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
*
         write(6,*) 'This should not be called'
   	end


* * * * * * * * * * * * * * * * * * * * * * * * * * *
* User-callable method to launch an interactive     *
* MINUIT session                                    *
* * * * * * * * * * * * * * * * * * * * * * * * * * *
	subroutine fit()
        implicit none
         write(6,*) 'This should not be called'
	end

* * * * * * * * * * * * * * * * * * * * * * * * * * *
* User-callable method to execute a single fit      *
* step, starting from parameters saved in a file    *
* * * * * * * * * * * * * * * * * * * * * * * * * * *
        subroutine fit1(carfile,out)
        implicit none
        integer out
        character*(*) carfile
         write(6,*) 'This should not be called'
        end

        subroutine fitN(carfile)
        implicit none
        character*(*) carfile
         write(6,*) 'This should not be called'
        end

* * * * * * * * * * * * * * * * * * * * * * * * * * *
* User-callable method to execute a single MINUIT   *
* card deck, usually to load a set of parameters    *
* * * * * * * * * * * * * * * * * * * * * * * * * * *
*        subroutine load1(carfile)
*        implicit none
*        character*(*) carfile
*         write(6,*) 'This should not be called'
*        end


* * * * * * * * * * * * * * * * * * * * * * * * * * *
* User-callable method to execute a single MINUIT   *
* command, through the function-call interface      *
* The procedure assumes that MINUIT has already     *
* been initialized via a call to load1 or fit1.     *
* * * * * * * * * * * * * * * * * * * * * * * * * * *
        integer function exec1(command)
        implicit none
        character*(*) command
         write(6,*) 'This should not be called'
        end
* * * * * * * * * * * * * * * * * * * * * * * * * * *
* User-callable method to execute a single MINUIT   *
* command file through the function-call interface  *
* The procedure assumes that MINUIT has already     *
* been initialized via a call to load1 or fit1.     *
* * * * * * * * * * * * * * * * * * * * * * * * * * *
*        integer function execF(carfile)
*        implicit none
*        character*(*) command
*         write(6,*) 'This should not be called'
*        end
      
* * * * * * * * * * * * * * * * * * * * * * * * * * *
* Dummyfroutine needed by MINUIT                  *
* * * * * * * * * * * * * * * * * * * * * * * * * * *
        real function ufunc(x)
        implicit none
        real x
         ufunc = 0.
        end


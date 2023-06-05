// CHEES header file
////////////////////////////////////////////////////
//                                                //
// NosPicker class header file                    //
//                                                //
// MSM 6 July 1995                                //
//                                                //
////////////////////////////////////////////////////

extern "C" 
{
  #include <wabbit.h>
}


#ifndef _NOSPICKER_CLASS_INCLUDED_
#define _NOSPICKER_CLASS_INCLUDED_

class NosPicker
{
 
 public:
  
  wabbitData_t *mData;           // stores data from histo
  /* public functions */
  
  NosPicker();                      // constructor
  ~NosPicker();                     // destructor

  double Pick();                    // actually picks a number from the CDF
  double Pick_real_phi(Four overall_CM_vec,
		       Four CMvec,
		       double *product_phi);
  int InitNosPicker(int, char *);   // initializes histogram
	

  private:
    
   /* private data */
      
   wabbitHeader_t mHeader[1];     // stores header from histo

   /* private functions */

   double doInterpolation(int, double);    // does linear interpolation
   int doBinarySearch(int, int, double);   // performs binary search

};

#endif
    

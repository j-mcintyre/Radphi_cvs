////////////////////////////////////////////////////////////////////////////
//                                                                        //
// NosPicker class code                                                   //
//                                                                        //
// MSM 6 July 1995                                                        //
//                                                                        //
// NosPicker stands for "Non-Object Code Specified" Picker, and it is and //
// object responsible for getting data from a CDF                         //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ntypes.h>
#include <Four.h>
#include "util.h"

#include "NosPicker.h"

extern Four into_rest_frame_of(Four vec1,Four vec2);
double costheta_between(Four vec1, Four vec2);
extern float angle_4v(Four vec1,Four vec2);

////////////////////////////////////////////////////////////////////////////

NosPicker::NosPicker()   // constructor
{
  mData = NULL;
}

////////////////////////////////////////////////////////////////////////////

int NosPicker::InitNosPicker(int histID, char *bbookFile)  // inits histo
{
  // InitNosPicker returns a 0 if the histogram file is 
  
  int error;

  fprintf(stderr, "Initializing histogram %d from file %s...", histID, 
	  bbookFile);

  /* now store header and data */
  if((error = bbGetHisto(histID, bbookFile, mHeader, &mData)))
    {
      switch(error)
	{
	case 1:
	  fprintf(stderr, "InitNosPicker: Problems were encountered with ");
	  fprintf(stderr, "the opening or format of file %s\n", bbookFile);
	  return error;
	case 2:
	  return 0;
	case 3:
	  fprintf(stderr, "InitNosPicker: Requested histogram id %d ", histID);
	  fprintf(stderr, "was missing from file %s\n", bbookFile);
	  return error;
	default:
	  return 0;
	}
    }  
  
  fprintf(stderr, "done\n");

  return error;
}

////////////////////////////////////////////////////////////////////////////

NosPicker::~NosPicker()
{
  if(mData)
    free(mData);
}

////////////////////////////////////////////////////////////////////////////

double NosPicker::Pick()
{
  /* first get a random number between 0 and 1 */
  double number = unfRand();
  
  
  /* now interpolate linearly between the upper and lower bins */
  if(number != (0.0 || 1.0))
    return doInterpolation(doBinarySearch(0,mHeader->nBins1-1,number), number);
  else if(number == 0.0)
    return 0.0;
  
  return 1.0;
}

////////////////////////////////////////////////
double NosPicker::Pick_real_phi(Four in_recoil,
				Four overall_cm_vec,
				double *product_phi)
{
  int hit=0;
  Four daut,recoil;
  double cos_theta_daut,sin_theta_daut,phi_daut;
  double cos_prnt_daut;
  double test_val;
  /* first isolate incoming data */
  recoil=in_recoil;
  /* renormalize recoil (recoil is now no longer a valid 4-vector)*/
  double recoil_total=sqrt((recoil.vec.space.x*recoil.vec.space.x)+
			   (recoil.vec.space.y*recoil.vec.space.y)+
			   (recoil.vec.space.z*recoil.vec.space.z));
  recoil.vec.space.x/=recoil_total;
  recoil.vec.space.y/=recoil_total;
  recoil.vec.space.z/=recoil_total;
  while(hit==0){
    cos_theta_daut=unfRand();
    cos_theta_daut=(cos_theta_daut*2.0)-1.0;
    /*fprintf(stdout,"cos_theta(daughter wrt recoil)=%f\t",cos_theta_daut);*/
    phi_daut=unfRand()*2.0*3.14159;
    daut.vec.space.z=cos_theta_daut;
    sin_theta_daut=sin(acos(cos_theta_daut));
    daut.vec.space.x=sin_theta_daut*cos(phi_daut);
    daut.vec.space.y=sin_theta_daut*sin(phi_daut);
    /*    fprintf(stdout,"x=%f y=%f z=%f\t",daut.vec.space.x,daut.vec.space.y,daut.vec.space.z);*/
    /*    fprintf(stdout,"px=%f py=%f pz=%f\n",recoil.vec.space.x,recoil.vec.space.y,recoil.vec.space.z); */
    /* construct cos(theta_(particle*daughter)) */
    cos_prnt_daut=(daut.vec.space.x*recoil.vec.space.x)+(daut.vec.space.y*recoil.vec.space.y)+(daut.vec.space.z*recoil.vec.space.z);
    /* cos_prnt_daut*=1.0; */
    test_val=unfRand();
    /* scale if needed, according to magic function */
    test_val=test_val*2.0;

    /* Here is where the magic function is stored */
    /*    fprintf(stdout,"%f %f %f\n",test_val,cos_prnt_daut,(1.0+(cos_prnt_daut*cos_prnt_daut))); */
    /*    if(test_val<(1.0+(cos_prnt_daut*cos_prnt_daut))) */
    if(cos_prnt_daut>.3&&cos_prnt_daut<.5)
      hit=1;
  }

  /*  fprintf(stdout,"helicity cos_theta=%f\n",cos_prnt_daut);*/
  /* move theta and phi into GJ frame */
  /*  fprintf(stdout,"recoil, daughter,angle in Picker:\n");
  recoil.print();
  daut.print(); 
  double angle_between=angle_4v(recoil,daut);
  fprintf(stdout,"angle between is=%f\n",angle_between);*/
  Four transrec=into_rest_frame_of(in_recoil,overall_cm_vec);
  Four dautGJ=into_rest_frame_of(daut,overall_cm_vec);
  /*  fprintf(stdout,"recoil should be the same:\n");
  transrec.print(); */
  
  *product_phi=phi_daut;
  /*  fprintf(stdout,"return value = %f\n",cos_theta_daut);*/
  /*  return cos_theta_daut; */
  double dautmag=sqrt((dautGJ.vec.space.x*dautGJ.vec.space.x)+
		      (dautGJ.vec.space.y*dautGJ.vec.space.y)+
		      (dautGJ.vec.space.z*dautGJ.vec.space.z));
  if(dautGJ.vec.space.x>0.0){
    (*product_phi)=atan(dautGJ.vec.space.y/dautGJ.vec.space.x);
    if(*product_phi<0.0)
      *product_phi+=6.283185;
  }
  else{
         (*product_phi)=atan(dautGJ.vec.space.y/dautGJ.vec.space.x)+3.14159;
  }

  /*fprintf(stdout,"costheta=%f,phi=%f.\n",dautGJ.vec.space.z/dautmag,*product_phi); */
  return dautGJ.vec.space.z/dautmag;
  /*  return .314159; */
}


////////////////////////////////////////////////////////////////////////////

double NosPicker::doInterpolation(int binNo, double randomNo)
{
  float binWidth = 0.0;
  float slope = 0.0;

  /* first find the bin width */
  binWidth= (mHeader->upperLimit1-mHeader->lowerLimit1) /
            ((float)(mHeader->nBins1));

  /* now find the slope of the liner interpolation */
  if(binNo)
    slope = (mData->data[binNo] - mData->data[binNo-1])/binWidth;
  else
    slope = (mData->data[binNo])/(binWidth);
    

  /* finally solve for x */
  return (randomNo - (double)mData->data[binNo])/(double)slope +
  ((double)mHeader->lowerLimit1 + binWidth*(double)binNo + (double)binWidth);
}

////////////////////////////////////////////////////////////////////////////

int NosPicker::doBinarySearch(int firstBin, int lastBin, double randomNo)
{
  int middleBin = (firstBin+lastBin)/2;

  if(firstBin == lastBin)
    return middleBin;

  if(randomNo < mData->data[middleBin])
    return doBinarySearch(firstBin, middleBin, randomNo);
  else if(randomNo > mData->data[middleBin])
    return doBinarySearch(middleBin+1, lastBin, randomNo);
  else 
    return middleBin;
}

double costheta_between(Four vec1, Four vec2)
{
  double one_dot_two,angle,arg;
  one_dot_two=(vec1.vec.space.x*vec2.vec.space.x)+(vec1.vec.space.y*vec2.vec.space.y)+(vec1.vec.space.z*vec2.vec.space.z);
  arg = one_dot_two/(vec1.mag3()*vec2.mag3());
  if(fabs(arg)>1){
    if (arg<0)
      return(PI);
    else
      return(0);
  }
  angle=one_dot_two/(vec1.mag3()*vec2.mag3());
  return angle;
}

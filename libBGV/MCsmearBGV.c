/*
 * MCsmearBGV : apply channel-dependent fluctuations to the simulation hit
 *              time and energy values.  This is done so that the standard
 *              E and t extraction that is done for real data in the function
 *              bgvClusterize() also works for simulated data.
 *
 * richard.t.jones@uconn.edu, October 4, 2005
 */


#include <stdio.h>
#include <makeHits.h>
#include <math.h>

// If you want to use the low-bias nonlinearity correction
// in makePhotonBGV() then comment out the following line.
#define SIMPLE_LEAKAGE_CORRECTION 1

#define MINIMUM_UPSTREAM_ENERGY 0.005
#define MINIMUM_DOWNSTREAM_ENERGY 0.005
#define PE_PER_GEV_UPSTREAM 120
#define PE_PER_GEV_DOWNSTREAM 35

/*
 * The following values must be consistent with those used 
 * in the simulation, in the hitconstants.inc include file.
 */

#define LEAD_LENGTH_BGV 86.5
#define C_EFFECTIVE_BGV 16.7
#define BGV_ATTENUATION_UP 130.0
#define BGV_ATTENUATION_DOWN 75.0
#define EXTRA_DOWNSTREAM_DELAY 0.5

/* cernlib functions */
void rannor_(float *a1,float *a2);   /* unit normal random generator */
float rngama_(float *p);             /* continuum Poisson generator */

void MCsmearBGV(bgv_hits_t *hits, time_list_t *tlist)
{
  int ih;
  float energy;
  float a1,a2;
  for (ih=0;ih<hits->nhits;ih++) {
    if (hits->hit[ih].times[0] && hits->hit[ih].times[1]) {
      int channel=hits->hit[ih].channel;
      float tu=tlist->element[hits->hit[ih].tindex[0]].le;
      float td=tlist->element[hits->hit[ih].tindex[1]].le;
      float Eu=hits->hit[ih].energy[0];
      float Ed=hits->hit[ih].energy[1];
      float sigma_tu=sqrt(0.3+0.060/Eu);
      float sigma_td=sqrt(0.3+0.140/Ed);
      float pu=Eu*PE_PER_GEV_UPSTREAM;
      float pd=Ed*PE_PER_GEV_DOWNSTREAM;
      float un0, un1;
      rannor_(&un0,&un1);
      tlist->element[hits->hit[ih].tindex[0]].le=tu+=sigma_tu*un0;
      tlist->element[hits->hit[ih].tindex[1]].le=td+=sigma_td*un1;
      Eu*=rngama_(&pu)/pu;
      Ed*=rngama_(&pd)/pd;
#ifdef SIMPLE_LEAKAGE_CORRECTION
      Eu*=1.45;
      Ed*=1.45;
#endif
      if (fabs(tu) < 20. && fabs(td) < 20. && fabs(tu-td) < 10.) {
        float z=(tu-td+EXTRA_DOWNSTREAM_DELAY)*C_EFFECTIVE_BGV/2;
        Eu*=exp((LEAD_LENGTH_BGV/2+z)/BGV_ATTENUATION_UP);
        Ed*=exp((LEAD_LENGTH_BGV/2-z)/BGV_ATTENUATION_DOWN);
        hits->hit[ih].energy[0]=Eu*exp(-z/bgv_atten_length[channel+24]);
        hits->hit[ih].energy[1]=Ed*exp(z/bgv_atten_length[channel]);
      }
    }
  }
}

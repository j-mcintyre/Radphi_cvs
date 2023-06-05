/*     make_bsd_pixels.c 
 *     R.T. Jones    June 9 2000
 *
 * May 17, 2005  R.T. Jones
 *  1) fixed the definition of "energy" member of pixel structure to
 *     estimate the dE/dx seen by one thickness of barrel scintillator.
 *     To accomplish this, I divide each measured energy by a gain factor
 *     that has been adjusted to take out the combination of path-length
 *     and attenuation factors, so that the most-probable "energy" value
 *     for a given particle is independent of direction in the barrel.
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <itypes.h> 
#include <iitypes.h>
#include <eventType.h>
#include <ntypes.h> 

#include <disData.h>
#include <pedestals.h>
#include <calibration.h> 
#include <makeHits.h>

#define  S_PITCH   (2.0*M_PI)/24.0
/* #define DEBUGGING 1 */

#define VLIGHT_BSD_RIGHT        9.8
#define VLIGHT_BSD_LEFT         9.8
#define VLIGHT_BSD_STRAIGHT     14.
#define Z0_BSD_RIGHT            12
#define Z0_BSD_LEFT             12.
#define Z0_BSD_STRAIGHT         12.

#define PIXEL_TIME_WIN		5.

int make_bsd_pixels(itape_header_t *event, bsd_pixels_t *p_pixels, int maxpix)
{

  bsd_hits_t *hits=NULL;
  bsd_pixels_t *pixels=NULL;
  time_list_t *timelist=NULL;
  int i,j,k;

  static char pixelTable[12][12][24];
  static int initOnce=0;
  if (!initOnce) {
    initOnce++;
    for (i=0;i<12;i++) { /* i is the right paddle index here */
      for (j=0;j<12;j++) { /* j is the left paddle index */
        int ipixel=(j-i+12)%12; /* this is the raw z position of the pixel */
        int iazim0=(2*i+ipixel)%24;
        int iazim1=(iazim0+1)%24;
        if(ipixel>=(BSD_N_RINGS))
          ipixel=-1; /* this corrects for pixel past the end of the barrel */
        for (k=0;k<24;k++) {
          pixelTable[i][j][k] = -1;
        }
        pixelTable[i][j][iazim0] = ipixel;
        pixelTable[i][j][iazim1] = ipixel;
      }
    }
  }

  if (p_pixels == NULL) {
    pixels = malloc(sizeof_bsd_pixels_t(maxpix));
  }
  else{
    pixels = p_pixels;
  }

  pixels->npixels = 0;

  if ((hits = data_getGroup(event,GROUP_BSD_HITS,0)) == NULL) {
    if (p_pixels == NULL)
      free(pixels);
    return(1);
  }
  if ((timelist = data_getGroup(event,GROUP_TIME_LIST,0)) == NULL) {
    fprintf(stderr,"Error in make_bsd_pixels: Failed to get time_list\n");
  }
  
  pixels->npixels = 0;
#ifdef DEBUGGING 
  fprintf(stderr,"make_bsd_pixels: hits structure has %d hits.\n",hits->nhits);
#endif
  if (hits->nhits > 2) {
    int hitlist[3][24];
    /* list of hits for each layer.  The first index is the layer, 
     * the second is the index within that layer.  The code below
     * only allows one hit per paddle per layer to be listed, so 
     * the 24 cannot by definition be exceeded.
     */
    int nhits[] = {0,0,0}; /* counts up hits in layers */
    for (i=0;i<hits->nhits;i++) {
      int layer=hits->hit[i].channel/12;
      layer=layer>2 ? 2:layer;
      /* the preceeding two statements make sure that the following 
       * convention is followed:
       * layer 0 (right twisting paddles): channels 0 to 11
       * layer 1 (left): channels 12 to 23
       * layer 2 (straight): channels 24 to 47
       */
      if (hits->hit[i].channel > 47)
        continue;
      /* The "channels" higher than 47 are the TDCed OR signals; 
       * thus the pixel reconstruction needs to ignore them.
       */
      if ((hits->hit[i].times > 0) ||
          (hits->hit[i].energy > AMIN_PIXEL)) {
#ifdef DEBUGGING 
        fprintf(stderr,"Hit number %d\n",i);
        if(hits->hit[i].times > 0)
          fprintf(stderr,"BSD hit by TDC.\n");
        if(hits->hit[i].energy > AMIN_PIXEL)
          fprintf(stderr,"BSD hit by ADC.\n");
#endif
        for (j=0;j<nhits[layer];j++) {
          if (hits->hit[i].channel ==
              hits->hit[hitlist[layer][j]].channel)
            break;
        }
        hitlist[layer][nhits[layer]++] = i;
      }
    }
    for (i=0;i<nhits[0];i++) { /* the i-th hit in the R layer */
      for (j=0;j<nhits[1];j++) { /* the j-th hit in the L layer */
        for (k=0;k<nhits[2];k++) { /* the k-th hit in the S layer */
          int ipix=pixelTable [hits->hit[hitlist[0][i]].channel]
                              [hits->hit[hitlist[1][j]].channel-12]
                              [hits->hit[hitlist[2][k]].channel-24];
          if (ipix >= 0) {
            /* if it's a valid pixel.  The pixel table is filled out with 
             * the z index if it's a valid pixel (from 0 to 7 inclusive),
             * or -1 if it's not a valid pixel
             */
            if (pixels->npixels == maxpix) {
#ifdef VERBOSE
              fprintf(stderr,"Warning in make_bsd_pixels: ");
              fprintf(stderr,"buffer overflow with %d pixels, truncating.\n",
                      maxpix+1);
#endif
              goto breakout;
            }
            pixels->pixel[pixels->npixels].ring = ipix;
            pixels->pixel[pixels->npixels].right = 
                                    hits->hit[hitlist[0][i]].channel;
            pixels->pixel[pixels->npixels].left =
                                    hits->hit[hitlist[1][j]].channel;
            pixels->pixel[pixels->npixels].straight =
                                    hits->hit[hitlist[2][k]].channel;
            pixels->pixel[pixels->npixels].z = bsd_ring_midz[ipix];
            pixels->pixel[pixels->npixels].phi =
                (pixels->pixel[pixels->npixels].straight - 22.5)*S_PITCH;
            if (pixels->pixel[pixels->npixels].phi >= 2.0*M_PI)
              pixels->pixel[pixels->npixels].phi -= 2.0*M_PI;
            pixels->pixel[pixels->npixels].energy[0] =
                		     hits->hit[hitlist[0][i]].energy /
                bsd_pixel_gain[ipix][hits->hit[hitlist[0][i]].channel];
            pixels->pixel[pixels->npixels].energy[1] =
                                     hits->hit[hitlist[1][j]].energy /
                bsd_pixel_gain[ipix][hits->hit[hitlist[1][j]].channel];
            pixels->pixel[pixels->npixels].energy[2] =
                                     hits->hit[hitlist[2][k]].energy /
                bsd_pixel_gain[ipix][hits->hit[hitlist[2][k]].channel];
            pixels->pixel[pixels->npixels].time[0] = FLT_MAX;
            pixels->pixel[pixels->npixels].time[1] = FLT_MAX;
            pixels->pixel[pixels->npixels].time[2] = FLT_MAX;
#ifdef DEBUGGING
            fprintf(stderr,"valid pixel: timelist=%d\n",timelist);
            fprintf(stderr,"hits in l0:%d, l1:%d,  l2:%d \n",
                    hits->hit[hitlist[0][i]].times,
                    hits->hit[hitlist[1][i]].times,
                    hits->hit[hitlist[2][i]].times);
#endif
            if (timelist) {
              int l1t0=hits->hit[hitlist[0][i]].tindex;
              int l2t0=hits->hit[hitlist[1][j]].tindex;
              int l3t0=hits->hit[hitlist[2][k]].tindex;
              int l1t1=hits->hit[hitlist[0][i]].times+l1t0;
              int l2t1=hits->hit[hitlist[1][j]].times+l2t0;
              int l3t1=hits->hit[hitlist[2][k]].times+l3t0;
              int i1,i2,i3;
              int goodtimes=0;
              for (i1=l1t0;i1<l1t1;i1++) {
                double t1=timelist->element[i1].le;
                for (i2=l2t0;i2<l2t1;i2++) {
                  double t2=timelist->element[i2].le;
                  for (i3=l3t0;i3<l3t1;i3++) {
                    double t3=timelist->element[i3].le;
                    if ((fabs(t1-t2) < PIXEL_TIME_WIN) &&
                        (fabs(t1-t3) < PIXEL_TIME_WIN) &&
                        (fabs(t2-t3) < PIXEL_TIME_WIN)) {
                      pixels->pixel[pixels->npixels].time[0] = t1;
                      pixels->pixel[pixels->npixels].time[1] = t2;
                      pixels->pixel[pixels->npixels].time[2] = t3;
                      ++goodtimes;
                    }
                  }
                }
              }
              /* Apply a time correction to take out the propagation
               * delay of the light in the scintillator from the point
               * of production to the pmt - courtesy of E. Smith
               */
              pixels->pixel[pixels->npixels].time[0] -= 
                 (pixels->pixel[pixels->npixels].z - Z0_BSD_RIGHT) /
	         VLIGHT_BSD_RIGHT;
              pixels->pixel[pixels->npixels].time[1] -=
                 (pixels->pixel[pixels->npixels].z - Z0_BSD_LEFT) /
	         VLIGHT_BSD_LEFT;
              pixels->pixel[pixels->npixels].time[2] -=
                 (pixels->pixel[pixels->npixels].z - Z0_BSD_STRAIGHT) /
	         VLIGHT_BSD_STRAIGHT;
              if (goodtimes) {
	        ++pixels->npixels;
              }
            }
            else {
	      ++pixels->npixels;
            }
          }
        }
      }
    }
  }

breakout:
  if (p_pixels == NULL) {
    int size = sizeof_bsd_pixels_t(pixels->npixels);
    bsd_pixels_t *tmp = data_addGroup(event,BUFSIZE,GROUP_BSD_PIXELS,0,size);
    memcpy(tmp,pixels,size);
    free(pixels); 
  }

  return(0);
}

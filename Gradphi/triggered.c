/* This function returns true (1) if the current event structures contain
 * a hit pattern that passes the online trigger cuts.  The cuts are:
 *  1) DOR condition: at least one block in LGD is over HIGH_THRESHOLD
 *  2) MAM condition: total energy in LGD is over MAM_THRESHOLD
 *  3) BSD condition: at least one TDC hit in each of three BSD layers
*/

#define HIGH_THRESHOLD 80
#define LOW_THRESHOLD 150
#define MAM_THRESHOLD 120
#define MAGIC_MAM_CONSTANT 16
#define GEV_PER_CHANNEL 0.0018
#define AVERAGE_PEDESTAL 35

#include <stdlib.h>
#include <ntypes.h>
#include <disData.h>
#include <lgdCluster.h>
#include <iitypes.h>

int triggered (itape_header_t* event)
{
  lgd_hits_t *lgdHits;
  int NoverThresh = 0;
  float CoverThresh = 0.;
  bsd_hits_t *bsdHits;
  int NwithTDC[] = {0,0,0,0};
  int n;

  if ((lgdHits = data_getGroup(event,GROUP_LGD_HITS,0)) == 0) {
    return 0;
  }
  for (n = 0; n < lgdHits->nhits; n++) {
    float adc = (lgdHits->hit[n].energy / GEV_PER_CHANNEL) + AVERAGE_PEDESTAL;
    if (adc > LOW_THRESHOLD) {
      ++NoverThresh;
    }
    if (adc > HIGH_THRESHOLD) {
      CoverThresh += adc;
    }
  }
  if ((NoverThresh == 0) ||
      (CoverThresh < MAM_THRESHOLD * MAGIC_MAM_CONSTANT)) {
    return 0;
  }
  
  if ((bsdHits = data_getGroup(event,GROUP_BSD_HITS,0)) == 0) {
    return 0;
  }
  for (n = 0; n < bsdHits->nhits; n++) {
    if (bsdHits->hit[n].times) {
      int layer = bsdHits->hit[n].channel / 12;
      ++NwithTDC[layer];
    }
  }
  if ((NwithTDC[0] == 0) ||
      (NwithTDC[1] == 0) ||
      (NwithTDC[2] + NwithTDC[3] == 0)) {
    return 0;
  }
  
  return 1;
}

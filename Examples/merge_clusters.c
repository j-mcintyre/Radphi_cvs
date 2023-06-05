#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iitypes.h>

//#define MAX_MASS_FOR_MERGING 0.1
#define MAX_MASS_FOR_MERGING 0.
#define SQR(X) ((X)*(X))

void merge_lowmass_cluster_pairs(photons_t *photons)
{
   int i_lightest, j_lightest;
   int i,j;
   double m_lightest = MAX_MASS_FOR_MERGING;
   for (i=0; i<photons->nPhotons; ++i) {
      if (photons->photon[i].cluster < 0)
         continue;
      for (j=i+1; j<photons->nPhotons; ++j) {
         if (photons->photon[j].cluster < 0)
            continue;
         double mass = sqrt(SQR(photons->photon[i].energy +
                                photons->photon[j].energy) -
                            SQR(photons->photon[i].momentum.x +
                                photons->photon[j].momentum.x) -
                            SQR(photons->photon[i].momentum.y +
                                photons->photon[j].momentum.y) -
                            SQR(photons->photon[i].momentum.z +
                                photons->photon[j].momentum.z)
                           );
         if (mass < m_lightest) {
            m_lightest = mass;
            i_lightest = i;
            j_lightest = j;
         }
      }
   }
   if (m_lightest < MAX_MASS_FOR_MERGING) {
      i = i_lightest;
      j = j_lightest;
      if (photons->photon[j].energy > photons->photon[i].energy) {
         photon_t copy = photons->photon[j];
         photons->photon[j] = photons->photon[i];
         photons->photon[i] = copy;
      }
      photons->photon[i].energy += photons->photon[j].energy;
      photons->photon[i].momentum.x += photons->photon[j].momentum.x;
      photons->photon[i].momentum.y += photons->photon[j].momentum.y;
      photons->photon[i].momentum.z += photons->photon[j].momentum.z;
      for (j++; j<photons->nPhotons; j++) {
         photons->photon[j-1] = photons->photon[j];
      }
      --photons->nPhotons;
      //printf("found a light pair with mass %f, pair was %d,%d, photon count is now %d\n",
      //       i_lightest,j_lightest,m_lightest,photons->nPhotons);
   }
}

extern "C" {
#include <math.h>
#include <stdlib.h>
#include <lgdGeom.h>
#include <lgdCluster.h>
}

#include <iostream>

#include <lgdClusterUC.h>

#if ! defined MIN_CLUSTER_BLOCK_COUNT
# define MIN_CLUSTER_BLOCK_COUNT 2
#endif

#if ! defined MIN_CLUSTER_SEED_ENERGY
# define MIN_CLUSTER_SEED_ENERGY 0.04
#endif

#define LOG2_WEIGHTING
#define SQR(x) (x)*(x)

/* 
 * Enable one of the following switches to choose the shape of the core
 * of the expectation function, in terms of shower energy.
 */
//#define HUNGRY_JACK        //High energy limit at sigma0=3.0
//#define SLIM_JIM           //High energy limit at sigma0=2.0
#define DEFAULT_PROFILE      //Standard shower parameterization

extern "C" {

void sortzv_(float* const a, int* index, const int* const n,
             const int* const mode, const int* const nway,
             const int* const nsort);

int ambiguous_events = 0;

/* Clusterize the hits given in <hits>. Returns cluster info
 * int <userClusterHits> and <clusters>.  These structures are
 * allocated below if the pointers are null on entry.
 */

int lgdClusterizeUC(lgd_hits_t *hits,lgd_hits_t **userClusterHits,
		    lgd_clusters_t **userClusters)
{
   float a[999];
   for (int ih = 0; ih < hits->nhits; ih++) {
      a[ih] = hits->hit[ih].energy;
   }
   int iord[999];
   const int one=1, zero=0;
   sortzv_(a,iord,&hits->nhits,&one,&one,&zero);

   LGDcluster::setHitlist(hits);
   LGDcluster* clusterList[999];
   int clusterCount = 0;
   int iter;
   for (iter=0; iter < 99; iter++) {

      // 1. At beginning of iteration, recompute info for all clusters.
      //    If something changed, return all hits to the pool and repeat.

      bool something_changed = false;
      for (int c = 0; c < clusterCount; c++) {
         something_changed |= clusterList[c]->update();
      }
      if (something_changed) {
         for (int c = 0; c < clusterCount; c++) {
            clusterList[c]->resetHits();
         }
      }
      else if (iter > 0) {
         break;
      }

      // 2. Look for blocks with energy large enough to require formation
      //    of a new cluster, and assign them as cluster seeds.

      for (int oh = 0; oh < hits->nhits; oh++) {
         int ih = iord[oh]-1;
         double energy = hits->hit[ih].energy;
         if (energy < MIN_CLUSTER_SEED_ENERGY)
            break;
         double totalAllowed = 0;
         for (int c = 0; c < clusterCount; c++) {
            totalAllowed += clusterList[c]->getEallowed(ih);
         }
         if (energy > totalAllowed) {
            clusterList[clusterCount] = new LGDcluster();
            clusterList[clusterCount]->addHit(ih,1.);
            clusterList[clusterCount]->update();
            ++clusterCount;
         }
         else if (iter > 0) {
            for (int c = 0; c < clusterCount; c++) {
               if (clusterList[c]->getHits())
                  continue;
               totalAllowed -= clusterList[c]->getEallowed(ih);
               if (energy > totalAllowed) {
                  clusterList[c]->addHit(ih,1.);
                  break;
               }
            }
         }
      }

      // 3. Share all non-seed blocks among seeded clusters, where
      //    any cluster shares a block if it expects at least 1 KeV in it.

      for (int oh = 0; oh < hits->nhits; oh++) {
         int ih = iord[oh]-1;
         if (LGDcluster::getUsed(ih) < 0)
            continue;
         double totalExpected = 0;
         for (int c = 0; c < clusterCount; c++) {
            if (clusterList[c]->getHits() > 0) {
               totalExpected += clusterList[c]->getEexpected(ih);
            }
         }
         for (int c = 0; c < clusterCount; c++) {
            if (clusterList[c]->getHits() > 0) {
               double expected = clusterList[c]->getEexpected(ih);
               if (expected > 1e-6) {
                  clusterList[c]->addHit(ih,expected/totalExpected);
               }
            }
         }
      }
   }
   
   if (iter == 99) {
      ++ambiguous_events;
   }
   
   if (*userClusters == 0) {
      *userClusters = 
          (lgd_clusters_t*) malloc(sizeof(lgd_clusters_t)*clusterCount);
   }
   (*userClusters)->nClusters = 0;
   int totalBlockCount = 0;
   for (int c = 0; c < clusterCount; c++) {
      int blockCount = clusterList[c]->getHits();
      (*userClusters)->nClusters += (blockCount < MIN_CLUSTER_BLOCK_COUNT)? 0:1;
      totalBlockCount += blockCount;
   }
   if (*userClusterHits == 0) {
      *userClusterHits = 
          (lgd_hits_t*) malloc(sizeof(lgd_hits_t)*totalBlockCount);
   }

   int nch = 0;
   int ncl = 0;
   for (int c = 0; (c < clusterCount) && (ncl < 8); c++) {
      int hitlist[999];
      int blockCount = clusterList[c]->getHits(hitlist,hits->nhits);
      if (blockCount < MIN_CLUSTER_BLOCK_COUNT) {
         continue;
      }
      else {
         for (int ih = 0; ih < blockCount; ih++) {
           (*userClusterHits)->hit[nch+ih].channel = 
                      hits->hit[hitlist[ih]].channel;
           (*userClusterHits)->hit[nch+ih].energy =
                      hits->hit[hitlist[ih]].energy;
         }
         (*userClusters)->cluster[ncl].nBlocks = blockCount;
         (*userClusters)->cluster[ncl].firstBlock = nch;
         (*userClusters)->cluster[ncl].flags = 0;
         (*userClusters)->cluster[ncl].energy = clusterList[c]->getEnergy();
         (*userClusters)->cluster[ncl].space = clusterList[c]->getCentroid();
         (*userClusters)->cluster[ncl].width = clusterList[c]->getRMS();
         nch += blockCount;
         ++ncl;
      }
   }
   (*userClusterHits)->nhits = nch;
   (*userClusters)->nClusters = ncl;

   for (int c = 0; c < clusterCount; c++) {
      delete clusterList[c];
   }
}
}


const lgd_hits_t* LGDcluster::fHitlist = 0;
int* LGDcluster::fHitused = 0;

LGDcluster::LGDcluster()
{
   fEnergy = 0;
   fEmax = 0;
   fCentroid.x = 0;
   fCentroid.y = 0;
   fCentroid.z = 0;
   fRMS = 0;
   fRMS_u = 0;
   fRMS_v = 0;
   fNhits = 0;

   if (fHitlist && fHitlist->nhits > 0) {
      fHit = new int[fHitlist->nhits];
      fHitf = new double[fHitlist->nhits];
      fEallowed = new double[fHitlist->nhits];
      fEexpected = new double[fHitlist->nhits];
   }
   else {
      fHit = 0;
      fHitf = 0;
      fEallowed = 0;
      fEexpected = 0;
   }
}

LGDcluster::~LGDcluster()
{
   if (fHit) delete [] fHit;
   if (fHitf) delete [] fHitf;
   if (fEallowed) delete [] fEallowed;
   if (fEexpected) delete [] fEexpected;
}

void LGDcluster::setHitlist(const lgd_hits_t* const hits)
{
   fHitlist = hits;
   if (fHitused) {
      delete [] fHitused;
      fHitused = 0;
   }
   if (hits->nhits > 0) {
      fHitused = new int[hits->nhits];
      for (int h = 0; h < hits->nhits; h++) {
         fHitused[h] = 0;
      }
   }
}

int LGDcluster::addHit(const int ihit, const double frac)
{
   if (ihit >= 0 && ihit < fHitlist->nhits && fNhits < fHitlist->nhits) {
      fHit[fNhits] = ihit;
      fHitf[fNhits] = frac;
      if (fNhits == 0) {
        fHitused[ihit] = -1;	// special used code for cluster seeds
      }
      else {
        ++fHitused[ihit];	// otherwise, just count owning clusters
      }
      ++fNhits;
      return 0;
   }
   else {
      return 1;
   }
}

void LGDcluster::resetHits()
{
   if (fNhits) {
      fHitused[fHit[0]] = 0;
      for (int h = 1; h < fNhits; h++) {
         --fHitused[fHit[h]];
      }
      fNhits = 0;
   }
}

bool LGDcluster::update()
{
   double energy = 0;
   for (int h = 0; h < fNhits; h++) {
      int ih = fHit[h];
      double frac = fHitf[h];
      energy += fHitlist->hit[ih].energy*frac;
   }
   double Emax=0;
   if (fNhits > 0) Emax = fHitlist->hit[fHit[0]].energy;

#ifdef LOG2_WEIGHTING
  /* This complicated centroid algorithm was copied from
   * Scott Teige's lgdClusterIU.c code -- don't ask [rtj]
   */
   vector3_t centroid;
   centroid.x = 0;
   centroid.y = 0;
   centroid.z = 0;
   double weight = 0.0;
   double weightSum = 0.0;
   double centerWeight = 0.0;
   double neighborMaxWeight = 0.0;
   double logFraction = exp(-0.23*(energy));
   double currentOffset = log(energy)/7.0 + 3.7;
   vector3_t centerBlock;
   for (int h = 0; h < fNhits; h++) {
      int ih = fHit[h];
      double frac = fHitf[h];
      vector3_t block;
      lgdLocalCoord(fHitlist->hit[ih].channel,&block);
   /*
    * Find the weight - offset determines minimum energy of block
    * to be used in weighting, since negative weights are thrown out
    */
      weight = currentOffset + log(fHitlist->hit[ih].energy*frac/energy);
      if (h == 0) {
         centerWeight = weight;
         lgdLocalCoord(fHitlist->hit[ih].channel,&centerBlock);
      }
      else {
         neighborMaxWeight =
               (neighborMaxWeight < weight)? weight:neighborMaxWeight;
      }
      if (weight > 0) {
         centroid.x += block.x*weight;
         centroid.y += block.y*weight;
         weightSum += weight;
      }
   }
   /*
    * Now patch up the center block's weight if it's got a neighbor
    * in the cluster that had positive weight
    */
   if (neighborMaxWeight > 0) {
      centroid.x += (logFraction-1)*(centerWeight-neighborMaxWeight)
                           *centerBlock.x;
      centroid.y += (logFraction-1)*(centerWeight-neighborMaxWeight)
                           *centerBlock.y;
      weightSum += (logFraction-1)*(centerWeight-neighborMaxWeight);
   }
   centroid.x /= weightSum;
   centroid.y /= weightSum;
#else
   vector3_t centroid;
   centroid.x = 0;
   centroid.y = 0;
   centroid.z = 0;
   for (int h = 0; h < fNhits; h++) {
      int ih = fHit[h];
      double frac = fHitf[h];
      vector3_t block;
      lgdLocalCoord(fHitlist->hit[ih].channel,&block);
      centroid.x += block.x*(fHitlist->hit[ih].energy*frac);
      centroid.y += block.y*(fHitlist->hit[ih].energy*frac);
   }
   centroid.x /= energy;
   centroid.y /= energy;
#endif

   double RMS = 0;
   double RMS_u = 0;
   double RMS_v = 0;
   for (int h = 0; h < fNhits; h++) {
      int ih = fHit[h];
      double frac = fHitf[h];
      vector3_t block;
      lgdLocalCoord(fHitlist->hit[ih].channel,&block);
      RMS += fHitlist->hit[ih].energy*frac
             *(SQR(block.x-centroid.x)+SQR(block.y-centroid.y));
      double u0 = sqrt(SQR(centroid.x)+SQR(centroid.y));
      double v0 = 0;
      double phi = atan2(centroid.y,centroid.x);
      double u = block.x*cos(phi)+block.y*sin(phi);
      double v =-block.x*sin(phi)+block.y*cos(phi);
      RMS_u += fHitlist->hit[ih].energy*frac*SQR(u-u0);
      RMS_v += fHitlist->hit[ih].energy*frac*SQR(v-v0);
   }
   fRMS = sqrt(RMS)/energy;
   fRMS_u = sqrt(RMS_u)/energy;
   fRMS_v = sqrt(RMS_v)/energy;

   bool something_changed = false;
   if (fabs(energy-fEnergy) > 0.001) {
      fEnergy = energy;
      something_changed = true;
   }
   if (fabs(Emax-fEmax) > 0.001) {
      fEmax = Emax;
      something_changed = true;
   }
   if (fabs(centroid.x-fCentroid.x) > 0.1 ||
       fabs(centroid.y-fCentroid.y) > 0.1) {
      fCentroid = centroid;
      something_changed = true;
   }
   if (something_changed) {
      for (int ih = 0; ih < fHitlist->nhits; ih++) {
         shower_profile(ih,fEallowed[ih],fEexpected[ih]);
      }
   }
   return something_changed;
}

#define MOLIER_RADIUS 8.5
#define MAX_SHOWER_RADIUS 20

void LGDcluster::shower_profile(const int ihit,
                                double& Eallowed, double& Eexpected) const
{
   Eallowed = Eexpected = 0;
   if (fEnergy == 0) return;
   vector3_t r;
   lgdLocalCoord(fHitlist->hit[ihit].channel,&r);
   double dist = sqrt(SQR(r.x-fCentroid.x)+SQR(r.y-fCentroid.y));
   if (dist > MAX_SHOWER_RADIUS) return;
   double theta = atan2(sqrt(SQR(fCentroid.x)+SQR(fCentroid.y)),120.);
   double phi = atan2(fCentroid.y,fCentroid.x);
   double u0 = sqrt(SQR(fCentroid.x)+SQR(fCentroid.y));   
   double v0 = 0;
   double u = r.x*cos(phi)+r.y*sin(phi);
   double v =-r.x*sin(phi)+r.y*cos(phi);
#if defined HUNGRY_JACK
   double vVar = SQR(2.0+2.5/sqrt(1.0+SQR(fEnergy-1.)));
   double uVar = vVar+SQR(20*theta);
   double vTail= 2.0 + .3*fEnergy + .7*sqrt(fEnergy) +10*theta;
   double uTail= vTail+20.*theta;
#elif defined SLIM_JIM
   double vVar = SQR(3.2+2.2/sqrt(2.5+SQR(fEnergy-1.)));
   double uVar = vVar+SQR(SQR(8*theta));
   double vTail= 4.2 + 0.5*log(0.9+fEnergy);
   double uTail= vTail+0.*theta+SQR(10*theta);
#elif defined DEFAULT_PROFILE
   double vVar = SQR(MOLIER_RADIUS/2.3);
//   double vVar = SQR(MOLIER_RADIUS/2.3  +0.8*log(fEnergy+1.));
   double uVar = vVar+SQR(SQR(8*theta));
   double vTail = 4.5+0.9*log(fEnergy+0.05);
   double uTail = vTail+SQR(10*theta);
#endif
   double core = exp(-0.5*SQR(SQR(u-u0)/uVar + SQR(v-v0)/vVar));
   double tail = exp(-sqrt(SQR((u-u0)/uTail)+SQR((v-v0)/vTail)));
   Eexpected = fEnergy*core;
   Eallowed = 2*fEmax*core + (0.2+0.5*log(fEmax+1.))*tail;

   if ((dist <= 4.) && (Eallowed < fEmax) ) {
      std::cerr << "Warning: lgdClusterUC Eallowed value out of range!\n";
      Eallowed = fEmax;
   }
}


#if ! defined CLASS_LGD_CLUSTER_
# define CLASS_LGD_CLUSTER_

#include <ntypes.h>

class LGDcluster
{
public :
   LGDcluster();
   ~LGDcluster();

   static void setHitlist(const lgd_hits_t* const hits);

   double getEexpected(const int ihit) const;
   double getEallowed(const int ihit) const;
   double getEnergy() const;
   double getEmax() const;
   vector3_t getCentroid() const;
   double getRMS() const;
   double getRMS_u() const;
   double getRMS_v() const;
   int getHits() const;
   int getHits(int hitlist[], int nhits) const;
   static int getUsed(const int ihit);
   int addHit(const int ihit, const double frac);
   void resetHits();
   bool update();

private :
   void shower_profile(const int ihit,
                       double& Eallowed, double& Eexpected) const ;

   double fEnergy;		// total cluster energy (GeV) or 0 if stale
   double fEmax;                // energy in the first block of the cluster
   vector3_t fCentroid;		// cluster centroid position (cm)
   double fRMS;			// cluster r.m.s. size (cm)
   double fRMS_u;		// cluster r.m.s. size in radial direction (cm)
   double fRMS_v;		// cluster r.m.s. size in azimuth direction (cm)
   int fNhits;			// number of hits owned by this cluster
   int *fHit;			// index list of hits owned by this cluster
   double *fHitf;		// list of hit fractions owned by this cluster
   static const lgd_hits_t* fHitlist;	// pointer to user's hit list
   static int *fHitused;	// number of clusters that use hit,
                                // or -1 if it is a cluster seed
   double *fEexpected;		// expected energy of hit by cluster (GeV)
   double *fEallowed;		// allowed energy of hit by cluster (GeV)
};

inline double LGDcluster::getEexpected(const int ihit) const
{
   if (ihit >= 0 && fHitlist && ihit < fHitlist->nhits)
      return fEexpected[ihit];
   else
      return 0;
}

inline double LGDcluster::getEallowed(const int ihit) const
{
   if (ihit >= 0 && fHitlist && ihit < fHitlist->nhits)
      return fEallowed[ihit];
   else
      return 0;
}

inline double LGDcluster::getEnergy() const
{
   return fEnergy;
}

inline double LGDcluster::getEmax() const
{
   return fEmax;
}

inline vector3_t LGDcluster::getCentroid() const
{
   return fCentroid;
}

inline double LGDcluster::getRMS() const
{
   return fRMS;
}

inline double LGDcluster::getRMS_u() const
{
   return fRMS_u;
}

inline double LGDcluster::getRMS_v() const
{
   return fRMS_v;
}

inline int LGDcluster::getHits() const
{
   return fNhits;
}

inline int LGDcluster::getHits(int hitlist[], const int nhits) const
{
   int ih;
   for (ih = 0; (ih < fNhits) && (ih < nhits); ih++) {
      hitlist[ih] = fHit[ih];
   }
   return ih;
}

inline int LGDcluster::getUsed(const int ihit)
{
   if (ihit >= 0 && fHitlist && ihit < fHitlist->nhits)
      return fHitused[ihit];
   else
      return 0;
}

#endif

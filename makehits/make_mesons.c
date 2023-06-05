/*     make_mesons.c
 *
 *     R.T. Jones	 March 2 2001
 *     - original version
 *
 *     R.T. Jones        April 8, 2003
 *     - removed the rule that suppressed the re-use of daughters that were
 *       used to form compound mesons at one level from being considered as
 *       constituents at subsequent levels.  This rule had been implemented
 *       before because of fears that re-use of photons from pi0's might
 *       create excessive combinatoric background in eta->2gamma and other
 *       channels.  However what it did was effectively to make it impossible
 *       to find eta->3pi0 through eta->sigma,pi0 after sigma->2pi0 because
 *       the sigma search allocated all of the pi0's to sigmas and left none
 *       free for the eta search.  It is now clear that the original search
 *       should cover the full combinatoric space of possibilities, and leave
 *       the resolution of ambiguous reconstructions to a later stage in the
 *       processing.  The rule can be recovered by defining the conditional
 *       HIDE_USED_CHILDREN.
 *
 * This function creates and fills the mesons_t structure in a similar fashion
 * to how make_photons creates and fills the photons_t structure.  This
 * algorithm is very simple and flexible, owing to the simple n-gamma final
 * state.  It allows the user to modify its behavior to make it find the
 * correct pattern in just about any event.
 *
 * Method:
 * The procedure is driven by an ordered list of decays for which to search.
 * A decay is specified by a mother type and width, and two daughter types.
 * The width is used to form an acceptance window around the true mass of the
 * mother.  A meson found at one level can be used to form a compound meson
 * at a subsequent level.  All acceptable combinations are kept, except that
 * it is forbidden that a compound decay use any given photon more than once
 * in its decay chain.
 *
 * Example:
 * The meson search path is specified in the call to setup_makeMesons.  For
 * example, suppose the call to setup_makeMesons gives the following list.
 *
 *          mother       width         daughter1    daughter2
 *          ------       -----         ---------    ---------
 *      1)   pi0         0.050          photon        photon
 *      2)   eta         0.050          photon        photon
 *      3)   a0_980      0.300          pi0           eta
 *      4)   phi         0.050          a0            gamma
 * 
 * In the above example, a search is first made for all photon pairs that
 * match a pi0.  All combinations are tried, so that any given photon may
 * appear as a constituent of more than one pi0.  The second step is an eta
 * search with a window close to the rms experimental mass resolution.  The
 * third step keeps all pi0,eta pairs in the mass range 700-1300 MeV, making
 * sure they are disjoint in their photon membership.  Finally a strict
 * condition is imposed on the phi->a0,gamma.  Note that the bachelor photon
 * in a phi->a0,gamma solution may also have belonged to a pi0 pair that
 * does not appear in the a0 decay tree; this does not contradict the above
 * reconstruction provided that another pi0 was found to satisfy a0->pi0,eta
 * decay kinematics.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>  
#include <string.h>
#include <math.h>

#include <iitypes.h> 

#include <disData.h>
#include <lgdCluster.h>
#include <lgdUtil.h>
#include <makeHits.h>

typedef struct {
   Particle_t type;
   vector3_t momentum;
   float energy;
   float mass;
   int ndecay;
   Particle_t dtype[2];
   int daughter[2];
   unsigned int used;
} picking_list_t;

picking_list_t *prod;
int nprod;

int setup_makeMesons(itape_header_t *event, meson_search_list_t *list)
{

   if (meson_list) {
      free(meson_list);
   }
   if (list) {
      if (list->stages > 0) {
         int mi,li;
         meson_list = malloc(sizeof_meson_search_list_t(list->stages));
         meson_list->stages = list->stages;
         for (mi=0;mi<list->stages;mi++) {
            meson_list->try[mi] = list->try[mi];
         }
      }
   }
   else {
      int stage=0;
      meson_list = malloc(sizeof_meson_search_list_t(9));

/* Define the default meson search algorithm here.  */

      meson_list->try[stage].meson = Pi0;
      meson_list->try[stage].daughter[0] = Gamma;
      meson_list->try[stage].daughter[1] = Gamma;
      meson_list->try[stage].width = 0.100;
      ++stage;
     
      meson_list->try[stage].meson = Eta;
      meson_list->try[stage].daughter[0] = Gamma;
      meson_list->try[stage].daughter[1] = Gamma;
      meson_list->try[stage].width = 0.150;
      ++stage;
    
      meson_list->try[stage].meson = omega;
      meson_list->try[stage].daughter[0] = Pi0;
      meson_list->try[stage].daughter[1] = Gamma;
      meson_list->try[stage].width = 0.200;
      ++stage;
    
      meson_list->try[stage].meson = PhiMeson;
      meson_list->try[stage].daughter[0] = Eta;
      meson_list->try[stage].daughter[1] = Gamma;
      meson_list->try[stage].width = 0.200;
      ++stage;
    
      meson_list->try[stage].meson = EtaPrime;
      meson_list->try[stage].daughter[0] = Gamma;
      meson_list->try[stage].daughter[1] = Gamma;
      meson_list->try[stage].width = 0.200;
      ++stage;
    
      meson_list->try[stage].meson = f0_980;
      meson_list->try[stage].daughter[0] = Pi0;
      meson_list->try[stage].daughter[1] = Pi0;
      meson_list->try[stage].width = 0.300;
      ++stage;
    
      meson_list->try[stage].meson = a0_980;
      meson_list->try[stage].daughter[0] = Eta;
      meson_list->try[stage].daughter[1] = Pi0;
      meson_list->try[stage].width = 0.300;
      ++stage;

      meson_list->stages = stage;
   }
   return 0;
}

int make_mesons(itape_header_t *event, mesons_t *p_mesons, int max)
{
  mesons_t *mesons;
  photons_t *photons;
  vector3_t vertex={0.0, 0.0, 0.0};
  int stage;
  int i,j;

  if (meson_list == NULL) {
     return 1;
  }
  photons = data_getGroup(event,GROUP_PHOTONS,0);
  if (photons == NULL) {
#ifdef VERBOSE
     fprintf(stderr,"make_mesons: no photons group, cannot make mesons!\n");
#endif
     return 1;
  }
  if (photons->nPhotons > 32) {
     fprintf(stderr,"make_mesons: cannot deal with %d photons in one event!\n",
             photons->nPhotons);
     return 1;
  }

/* Here we form an internal particle table for the clustering algorithm
 * to use.  It starts off with the final-state photons and as new mesons
 * are found they are appended to the list. Never mind that a photon is
 * not a meson, the photons will not be copied into the output table.
 */

  prod = malloc(sizeof(picking_list_t)*9999);
  for (nprod=0; nprod<photons->nPhotons; nprod++) {
     prod[nprod].type = Gamma;
     prod[nprod].momentum = photons->photon[nprod].momentum;
     prod[nprod].energy = photons->photon[nprod].energy;
     prod[nprod].mass = 0;
     prod[nprod].ndecay = 0;
     prod[nprod].used = (1 << nprod);
  }

/* Pick the mesons according to the list */

  for (stage=0; stage<meson_list->stages; stage++) {
     Particle_t mother = meson_list->try[stage].meson;
     float mass_min = ParticleMass(mother) - meson_list->try[stage].width;
     float mass_max = ParticleMass(mother) + meson_list->try[stage].width;
     Particle_t type1 = meson_list->try[stage].daughter[0];
     Particle_t type2 = meson_list->try[stage].daughter[1];
     int npicking = nprod;
     int c1,c2;
     for (c1=0; c1<npicking; c1++) {
        if (prod[c1].type == type1) {
           if (type1 == type2)
              c2 = c1+1;
           else
              c2=0;
           for ( ; c2<npicking; c2++) {
              if (prod[c2].type == type2) {
                 if ((prod[c1].used & prod[c2].used) == 0) {
                    vector3_t p;
                    float p_0;
                    float mass;
                    p.x = prod[c1].momentum.x + prod[c2].momentum.x;
                    p.y = prod[c1].momentum.y + prod[c2].momentum.y;
                    p.z = prod[c1].momentum.z + prod[c2].momentum.z;
                    p_0 = prod[c1].energy + prod[c2].energy;
                    mass = sqrt(p_0*p_0 - p.x*p.x - p.y*p.y -p.z*p.z);
                    if ((mass > mass_min) && (mass < mass_max)) {
                       prod[nprod].type = mother;
                       prod[nprod].momentum.x = p.x;
                       prod[nprod].momentum.y = p.y;
                       prod[nprod].momentum.z = p.z;
                       prod[nprod].energy = p_0;
                       prod[nprod].mass = mass;
                       prod[nprod].ndecay = 2;
                       prod[nprod].dtype[0] = type1;
                       prod[nprod].dtype[1] = type2;
                       prod[nprod].daughter[0] = c1;
                       prod[nprod].daughter[1] = c2;
                       prod[nprod].used = (prod[c1].used | prod[c2].used);
                       ++nprod;
                    }
                 }
              }
           }
        }
     }
#if defined HIDE_USED_CHILDREN
     for ( ; npicking<nprod; npicking++) {
        int child;
        for (child=0; child<prod[npicking].ndecay; child++) {
           int kid = prod[npicking].daughter[child];
           prod[kid].used = 0xffffffff;
        }
     }
#endif
  }

/* Save the mesons to a permanent structure */

  if (p_mesons == NULL){
    mesons = malloc(sizeof_mesons_t(max));
  }
  else{
    mesons = p_mesons;
  }

  mesons->nMesons = nprod - photons->nPhotons;
  if (mesons->nMesons > max){
#ifdef VERBOSE
    fprintf(stderr,"Warning in make_mesons: ");
    fprintf(stderr,"buffer overflow at %d entries, truncating.\n",max+1);
#endif
    mesons->nMesons = max;
  }

  nprod = photons->nPhotons;
  for (i=0; i<mesons->nMesons; i++) {
     mesons->meson[i].type = prod[nprod].type;
     mesons->meson[i].momentum = prod[nprod].momentum;
     mesons->meson[i].vertex = vertex;
     mesons->meson[i].energy = prod[nprod].energy;
     mesons->meson[i].mass = prod[nprod].mass;
     mesons->meson[i].ndecay = prod[nprod].ndecay;
     for (j=0; j<mesons->meson[i].ndecay; j++) {
        mesons->meson[i].dtype[j] = prod[nprod].dtype[j];
        mesons->meson[i].daughter[j] = 
           (mesons->meson[i].dtype[j] == Gamma) ?
             prod[nprod].daughter[j]
           : prod[nprod].daughter[j] - photons->nPhotons;
     }
     ++nprod;
  }

  free(prod);

  if (p_mesons == NULL){
    int size = sizeof_mesons_t(mesons->nMesons);
    mesons_t *tmp = data_addGroup(event,BUFSIZE,GROUP_MESONS,0,size);
    memcpy(tmp,mesons,size);
    free(mesons);
  }

  return 0;
}

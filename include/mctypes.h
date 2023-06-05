/*
 * mctypes.h - data structures for the MonteCarlo data group
 *
 * The Monte Carlo data group, like the raw and unpacked data groups,
 * contains descriptive information about an event, in this case from
 * simulation.  A file containing only the Monte Carlo group for each
 * event is produced by an event generator for input to the simulation
 * program (Gradphi for GEANT simulation of the Radphi experiment).
 * Normally the output from the simulation appends the MC group to
 * the output event for reference later, in case the parameters of
 * the original process would be of interest in the analysis.
*/

#ifndef MCTYPES_INCLUDED
#define MCTYPES_INCLUDED

#include <ntypes.h>      /* for int32, etc... definitions */

#include <particleType.h> /* define the basic particle types */

typedef struct
{
  uint32 kind;		/* particle types are defined in particleType.h */
  vector4_t momentum;	/* in (GeV,GeV/c) */
} mc_particle_t;

typedef struct
{
  uint32 npart;		/* number of particles produced at this vertex */
  uint32 kind;		/* particle type in case of decay vertex, ow 0 */
  vector4_t origin;	/* relative to nominal target center in cm     */
  mc_particle_t part[1]; /* list of particles produced at this vertex  */
} mc_vertex_t;

typedef struct
{
  uint32 npart;		/* number of initial-state particles           */
  int32 type;		/* user-defined reaction type                  */
  float32 weight;	/* for weighted simulations, usually 1         */
  uint32 nfinalVert;    /* number of vertices in the final state       */
  mc_particle_t part[1]; /* list of particles in the initial state     */
} mc_event_t;

#define sizeof_mc_event_t(N)  \
        (sizeof(mc_event_t) + (((N)>0)?(N)-1:0)*sizeof(mc_particle_t))

#define sizeof_mc_vertex_t(N) \
        (sizeof(mc_vertex_t) + (((N)>0)?(N)-1:0)*sizeof(mc_particle_t))

#endif

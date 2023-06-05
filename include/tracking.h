/*
 * tracking.h --- MPS track finding and fitting stuff
 *
*/

#ifndef TRACKING_INCLUDED
#define TRACKING_INCLUDED

#include <ntypes.h>
#include <disData.h>

#ifndef BIT
#define BIT(n) (1<<(n))
#endif
/*
BIT     HEX
0       0x1 
1       0x2 
2       0x4 
3       0x8 
4       0x10 
5       0x20 
6       0x40 
7       0x80 
8       0x100 
9       0x200 
10      0x400 
11      0x800 
12      0x1000 
13      0x2000 
14      0x4000 
15      0x8000 
16      0x10000 
17      0x20000 
18      0x40000 
19      0x80000 
20      0x100000 
*/


/*
 * The TRK_F_xxx flags defined here are used in the 'flags' field of the tracks and vertices
 * data structures.
*/

#define TRK_F_BEAM          BIT(16+0)     /* this bit is set for the beam track */
#define TRK_F_RECOIL        BIT(16+1)     /* this bit is set for the missing momentum and charge track */
#define TRK_F_LGDCLUSTER    BIT(16+2)     /* this bit is set for tracks created from LGD clusters      */
#define TRK_F_SECVTX        BIT(16+5)     /* this bit is set for tracks connecting vertices            */
#define TRK_F_FAKE          BIT(16+6)     /* this bit is set for the "fake" beam track */

#define VTX_F_SINGLE        0x0001        /* single-vertex fit (no multiple vertices allowed) */
#define VTX_F_PRIMARY       0x0002        /* primary vertex (with the beam track)             */
#define VTX_F_SECONDARY     0x0004        /* secondary vertex */


typedef struct
{
  int32      index;         /* 
		             * for TRK_F_MPSTRACK tracks   - index into the tracks_t array
		             * for TRK_F_LGDCLUSTER tracks - index into the array of LGD clusters.
		             * for TRK_F_TCYLTRACK tracks  - index into the array of TCYL tracks
			     * for TRK_F_SECVTX tracks     - index into the array of vertices: index
			     *                                                of the secondary vertex
			     *                                                this track connects to.
		            */
  int32      vertex;        /* number of the vertex (counting from 0) */
  vector3_t  p;             /* momentum 3-vector in GeV/c  */
  int32      charge;        /* charge (usually -1, 0 or 1) */
  float32    p_cov[3*3];    /* covariance matrix (usually TRK_F_COV2_xxx or TRK_F_COV3_xxx, as set of flags) */
  Particle_t particleType;
  float32    xi2pdf;        /* This xi2pdf is for the fit with the vertex */
  uint32     flags;         /* see TRK_F_xxx flags above */
} geo_track_t;

typedef struct
{
  vector3_t  v;              /* vertex position */
  float32    v_cov[3*3];     /* "covariance" matrix. Contents is defined by the TRK_F_COVxxx bits set in flags */
  int32      ntracks;        /* number of tracks on this vertex */
  int32      firstTrack;     /* offset of the first track into the geo_track_t array */
  int32      vertexType;
  uint32     flags;          /* see VTX_F_xxx flags above */
  float32    xi2pdf;         /* "quality of the fit"      */
} geo_vertex_t;

/* GROUP_GEO_VERTICES */

typedef struct
{
  int32        nvertices;
  geo_vertex_t vertex[1];
} geo_vertices_rec_t;

/* GROUP_GEO_TRACKS */

typedef struct
{
  int32       ntracks;
  geo_track_t track[1];
} geo_tracks_rec_t;

#endif
/* end file */

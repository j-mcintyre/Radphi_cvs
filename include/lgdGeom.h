/* lgdGeom.h 4.1 
   10 Dec 1995 18:48:23
   /auto/home/lemond/e852/source/include/SCCS/s.lgdGeom.h
*/

#include <ntypes.h>
#include <encoderHandle.h>

#ifndef LGDGEOM_INCLUDED
#define LGDGEOM_INCLUDED


typedef struct{
  int row; /* row of the block */
  int col; /* col of the block */
  int channel; /* logical channel number */
  int hole; /* True if this is a hole (corner/beam hole) */
  vector3_t blockSize;
  vector3_t space; /* X/Y/Z in LGD coordinates */
} lgdGeom_t;


#define LGDGEOM_OK 0
#define LGDGEOM_READERR 1

#define LGDGEOM_MAXHOLES 3053

#define LGDGEOM_BYRAWADC 1
#define LGDGEOM_BYROWCOL 2
#define LGDGEOM_BYCHANNEL 3
#define LGDGEOM_RAWERROR -1
#define LGDGEOM_ROWCOLERROR -2
#define LGDGEOM_CHANNELERROR -3
#define LGDGEOM_METHODERR -4
#define LGDGEOM_ISAHOLE -5

#define LGDGEOM_COORDERR -99999.0

#define LGDGEOM_TRUE 1
#define LGDGEOM_FALSE 0


int lgdGeomSetup(int runNo);
int lgdGeomIsAnEdge(int channel);

/*
 * in the 'lgdGeomGetNeighbors' call, the 'neighbor' parameter should have at least 8 elements
 * or memory would be overwritten
*/

int lgdGeomGetNeighbors(int channel,int neighborChannels[8],int *numberOfNeighbors);

/*
 * in the 'lgdGetGeom' call, the 'method' argument is one of the LGDGEOM_BYxxxx defined above
*/

int lgdGetGeom(unsigned long adcChannel,lgdGeom_t *block,int method);

int lgdGetChannelRaw(unsigned long ADCWord, int *channel);

/*
 * the lgdLocal...(...) functions return positions of the CENTER of the block
*/

float lgdLocalX(int channel);
float lgdLocalY(int channel);
float lgdLocalZ(int channel);
int   lgdLocalCoord(int channel,vector3_t *local);
void  lgdLocalToGlobalCoord(vector3_t local,vector3_t *mps);

void  lgdLocateXporter(encoderData_t *encoders);

void lgdGlobalToLocalCoord(vector3_t *local,vector3_t mps);
int  lgdGlobalCoord(int channel,vector3_t *mps);

void lgdGetGeomParam(int *nchannels,int *nrows, int *ncolumns);

/*
 * lgdBlockSize(...) returns one LGD block size: x, y, and z
*/

void lgdModuleSize(vector3_t *block);
void lgdBlockSize(int channel, vector3_t *block);
int lgdLocalPositionToGeomType(vector3_t local, lgdGeom_t *item);
int lgdGlobalPositionToGeomType(vector3_t global,lgdGeom_t *item);
#endif

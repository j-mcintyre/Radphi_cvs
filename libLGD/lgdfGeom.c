/* 
        lgdfGeom.c
         Created           :  28-MAR-1994 by Thom Sulanke
*/

#include "lgdGeom.h"

/***********************************************************************/
/*                                                                     */
/*   lgdfGeom.c                                                        */
/*   ----------                                                        */
/*                                                                     */
/*         Created     :  28-MAR-1994    Author : Thom Sulanke         */
/*         Purpose     :  FORTRAN interface routines for lgdGeom       */
/*                                                                     */
/***********************************************************************/
 
int lgdfgeomsetup_( int *runNum )
 
{
  return lgdGeomSetup( *runNum );
}


void lgdfgetgeomparam_( int *nChannel, int *nRows, int *nCols )

{
  lgdGetGeomParam( nChannel, nRows, nCols );
}


int lgdfgeomisanedge_( int *channel )

{
  return lgdGeomIsAnEdge( *channel );
}


int lgdfgeomgetneighbors_( int *channel, int neighbors[], int *nNeighbors )

{
  return lgdGeomGetNeighbors( *channel, neighbors, nNeighbors );
}


int lgdfgetgeom_( int *adc, lgdGeom_t *block, int *method )

{
  return lgdGetGeom( *adc, block, *method );
}


int lgdfgetchannelraw_( int *adc, int *channel )

{
  return lgdGetChannelRaw( *adc, channel );
}


float lgdflocalx_( int *channel )

{
  return lgdLocalX( *channel );
}


float lgdflocaly_( int *channel )

{
  return lgdLocalY( *channel );
}


float lgdflocalz_( int *channel )

{
  return lgdLocalZ( *channel );
}


int lgdflocalcoord_( int *channel, vector3_t *position )

{
  return lgdLocalCoord( *channel, position );
}


/*void lgdflocatexporter_( int *encoderx, int *encodery )

{
  lgdLocateXporter( *encoderx, *encodery);
}

*/
void lgdflocaltoglobalcoord_( vector3_t *lgdcoord, vector3_t *mpscoord )

{
  lgdLocalToGlobalCoord( *lgdcoord, mpscoord );
}


void lgdfglobaltolocalcoord_( vector3_t *lgdcoord, vector3_t *mpscoord )

{
  lgdGlobalToLocalCoord( lgdcoord, *mpscoord );
}


int lgdfglobalcoord_( int *channel, vector3_t *position )

{
  return lgdGlobalCoord( *channel, position );
}

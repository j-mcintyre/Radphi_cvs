/*
   $Log: lgdCal.h,v $
   Revision 1.1  1997/05/24 00:49:36  radphi
   Initial revision by lfcrob@jlabs2
   Lgd Calibration routines

*/

#ifndef LGDCAL_H_INCLUDED
#define LGDCAL_H_INCLUDED

#include <lgdCluster.h>
#define CAL_PRIM_SIZE 784
#define CAL_SEC_SIZE 407
void ecal0(int ncc[CAL_PRIM_SIZE],float store[CAL_PRIM_SIZE][CAL_SEC_SIZE],float x[CAL_PRIM_SIZE],
           float u[CAL_PRIM_SIZE]);
void ecal(float efit,lgd_cluster_t *cluster,lgd_hits_t *hits,
          float cc[CAL_PRIM_SIZE],int ncc[CAL_PRIM_SIZE],float store[CAL_PRIM_SIZE][CAL_SEC_SIZE],
          float x[CAL_PRIM_SIZE],float u[CAL_PRIM_SIZE]);
void esolv(float cc[CAL_PRIM_SIZE],float newcc[CAL_PRIM_SIZE],int ncc[CAL_PRIM_SIZE],
           float store[CAL_PRIM_SIZE][CAL_SEC_SIZE],float x[CAL_PRIM_SIZE],float u[CAL_PRIM_SIZE]);

#endif

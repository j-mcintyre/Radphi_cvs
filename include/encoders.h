/*
  $Log: encoders.h,v $
  Revision 1.1  1997/04/08 21:09:49  radphi
  Initial revision

  */

#ifndef ENCODERS_H_INCLUDED
#define ENCODERS_H_INCLUDED

int getNormEncodeChannel(int *norm);
int getLGDEncodeChannels(int *xChannel, int *yChannel);
int getRPDEncodeChannels(int *xChannel, int *yChannel, int *zChannel);
int readNormEncoder(int *norm);
int readLGDEncoders(int *x, int *y);
int readRPDEncoders(int *x, int *y, int *z);

int readNormGain(int *norm);
int readLGDGain(int *x, int *y);
int readRPDGain(int *x, int *y, int *z);

int initEncoders(void);

#endif

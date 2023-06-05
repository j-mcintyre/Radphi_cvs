/*
 *  pcntwk.h  --
 *	Map the /PCNTWK/ and /PCNTW1/ commons
 *
 *  Author:   Richard Jones <richard.t.jones@uconn.edu>
 *
 */

#ifndef CERN_PCNTWK
#define CERN_PCNTWK

#include	<cfortran/cfortran.h>


typedef struct {
	int	npass;
	int	ipass;
	int	iminev;
	int	imaxev;
	int	misbyt;
	int	nchain;
	int	ntcont;
	int	iedone;
} pcntwk_def;

#define PCNTWK COMMON_BLOCK(PCNTWK,pcntwk)
COMMON_BLOCK_DEF(pcntwk_def,PCNTWK);


typedef struct {
	char 	crhost[16];
} pcntw1_def;

#define PCNTW1 COMMON_BLOCK(PCNTW1,pcntw1)
COMMON_BLOCK_DEF(pcntw1_def,PCNTW1);

#endif

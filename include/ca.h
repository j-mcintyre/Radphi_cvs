/*-----------------------------------------------------------------------------
 * Copyright (c) 1991,1992 Southeastern Universities Research Association,
 *                         Continuous Electron Beam Accelerator Facility
 *
 * This software was developed under a United States Government license
 * described in the NOTICE file included as part of this distribution.
 *
 * CEBAF Data Acquisition Group, 12000 Jefferson Ave., Newport News, VA 23606
 * Email: coda@cebaf.gov  Tel: (804) 249-7101  Fax: (804) 249-7363
 *-----------------------------------------------------------------------------
 * 
 * Description:
 *	Function prototypes for CAMAC standard routines
 *	
 * Author:  Chip Watson, CEBAF Data Acquisition Group
 *
 * Revision History:
 *   $Log: ca.h,v $
 *   Revision 1.1  1997/01/30 19:16:18  lfcrob
 *   Initial revision
 *
 *   Revision 1.1.1.1  1996/08/21 19:18:46  heyes
 *   Imported sources
 *
*	  Revision 1.1  94/03/16  07:53:28  07:53:28  heyes (Graham Heyes)
*	  Initial revision
*	  
*	  Revision 1.1  94/03/15  11:47:38  11:47:38  heyes (Graham Heyes)
*	  Initial revision
*	  
 *	  Revision 1.3  93/05/13  14:31:52  14:31:52  watson (Chip Watson)
 *	  fix ccci macro expansion
 *	  
 *	  Revision 1.2  1992/12/08  18:52:54  watson
 *	  add block I/O routines, conditional inline expansion
 *
 *	  Revision 1.1  1992/06/05  19:25:13  watson
 *	  Initial revision
 *
 */

void caopen(char *server,int *success);
void ccinit(int b);
void cfga(int fa[],int exta[],int intc[],int qa[],int cb[]);
void cfmad(int f,int extb[],int intc[],int cb[]);
void cfubc(int f,int extb[],int intc[],int cb[]);
void cfubl(int f,int extb[],int intc[],int cb[]);
void cfubr(int f,int extb[],int intc[],int cb[]);
void cglam(int lam,int *b,int *c,int *n,int *m,int inta[]);
void cgreg(int lam,int *b,int *c,int *n,int *a);
void csga(int fa[],int exta[],unsigned short intc[],int qa[],int cb[]);
void csmad(int f,int extb[],unsigned short intc[],int cb[]);
void csubc(int f,int extb[],unsigned short intc[],int cb[]);
void csubl(int f,int extb[],unsigned short intc[],int cb[]);
void csubr(int f,int extb[],unsigned short intc[],int cb[]);

#ifndef CAMAC_BUILD_INLINE

void cdreg(int *ext,int b,int c,int n,int a);
void cfsa(int f,int ext,int *data,int *q);
void cccz(int ext);
void cccc(int ext);
void ccci(int ext,int l);
void ctci(int ext,int *l);
void cccd(int ext,int l);
void ctcd(int ext,int *l);
void ccgl(int ext,int l);
void ctgl(int ext,int *l);
void cdlam(int *lam,int b,int c,int n,int a,int *inta);
void cclm(int lam,int l);
void cclc(int lam);
void ctlm(int lam,int *l);
void cclwt(int lam);
void ctstat(int *istat);
void cssa(int f,int ext,unsigned short *data,int *q);

#else

#include "ks2917.h"

#define cdreg(ext,b,c,n,a) CDREG(ext,b,c,n,a)
#define cfsa(f,ext,dat,q)  CFSA(f,ext,dat,q)
#define cccz(ext)	   CCCZ(ext)
#define cccc(ext)	   CCCC(ext)
#define ccci(ext,l)	   CCCI(ext,l)
#define ctci(ext,l)	   CTCI(ext,l)
#define cccd(ext,l)	   CCCD(ext,l)
#define ctcd(ext,l)	   CTCD(ext,l)
#define ccgl(lam,l)	   CCGL(lam,l)
#define ctgl(ext,l)	   CTGL(ext,l)
#define cdlam(lam,b,c,n,a,inta) CDLAM(lam,b,c,n,a,inta)
#define cclc(lam)	   CCLC(lam)
#define cclm(lam,l)	   CCLM(lam,l)
#define ccrgl(lam)	   CCRGL(lam)
#define ctstat(istat)	   CTSTAT(istat)
#define cssa(f,ext,dat,q)  CSSA(f,ext,dat,q)

#endif

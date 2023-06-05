
/* 
 * $Log: camac.h,v $
 * Revision 1.9  2000/05/02 15:10:45  radphi
 * Changes made by radphi@urs3
 * adding multi-crate, and turning off controller inhibit line
 *
 * Revision 1.8  2000/05/02 14:34:54  radphi
 * Changes made by radphi@urs3
 * updated to allow more than 1 camac crate
 *
 * Revision 1.7  1997/04/08 21:08:34  radphi
 * Removed #define's for module slots
 *
 * Revision 1.6  1997/03/28 22:04:17  lfcrob
 * *** empty log message ***
 *
 * Revision 1.5  1997/02/12 18:53:40  lfcrob
 * Fixed LGDXP slot
 *
 * Revision 1.4  1997/02/10  01:52:10  lfcrob
 * fixed BRANCH to be 0
 *
 * Revision 1.3  1997/02/10  01:34:53  lfcrob
 * Started adding some $defines for module locations
 *
 * Revision 1.2  1997/01/31  17:23:07  lfcrob
 * Added the rest of the protos
 *
 * Revision 1.1  1997/01/30  15:42:26  lfcrob
 * Initial revision
 *
 */

/* Note - protottypes are being added as they are needed ;) */

#ifndef CAMAC_H_INCLUDED
#define CAMAC_H_INCLUDED

#define BRANCH 0
#define CRATE 0

#define MAX_N_CRATES 2

void caopen(char *server, int *success);
void ctstat(int *istat);
void ccinit(int b);

void cdreg(int *ext, int b, int c, int n, int a);
void cfsa(int f, int ext, int *data, int *q);
void cccz(int ext);
void cccc(int ext);
void ccci(int ext, int logic);
void ctci(int ext, int *logic);
void cccd(int ext, int logic);
void ctcd(int ext, int *logic);
void ctgl(int ext, int *logic);
void cdlam(int *lam, int b, int c, int n, int a, int *inta);
void cclm(int lam, int logic);
void ctlm(int lam, int *logic);


void cfga(int fa[],int exta[], int intc[], int qa[], int cb[]);
void cfmad(int f, int extb[],int intc[], int cb[]);
void cfubc(int f, int ext, int intc[], int cb[]);
void cfubl(int f, int ext, int intc[], int cb[]);
void cfubr(int f, int ext, int intc[], int cb[]);
void csga(int fa[],int exta[], int intc[], int qa[], int cb[]);
void csmad(int f, int extb[], int intt[], int cb[]);
void csubc(int f, int ext, int intt[], int cb[]);
void csubl(int f, int ext, int intt[], int cb[]);
void csubr(int f, int ext, int intt[], int cb[]);


#endif

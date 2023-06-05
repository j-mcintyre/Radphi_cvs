/*
 * This header file is the definition of the paw ntuple blocks for the
 * Radphi hits information.  The hits are one level above raw data in
 * that basic calibration corrections have been applied and clustering
 * of signals that belong to a single hit has been carried out.  The
 * the list of blocks currently defined are:
 *      /bsd_hits/	BSD hit data
 *      /bsd_pixs/	BSD pixel data
 *	/bgv_hits/	BGV hit data
 *	/cpv_hits/	CPV hit data
 *	/upv_hits/	UPV hit data
 *	/tag_hits/	TAG hit data
 *	/time_list/	multihit time data
 * Timelist is used for multihit timing information, and is shared
 * between all of the blocks except lgd_hits.  New groups are to be
 * added as the corresponding hits structures are defined and filled.
 *	/lgd_hits/	LGD hit data
 *	/photons/       reconstructed photons
 *	/mesons/	reconstructed mesons
 *	/tag_coin/	tagger coincidence data
 *	/recoils/	barrel scintillator clusters
 *	/bclusters/	barrel gamma clusters
 *	/tag_align/	tagger time table for alignment studies
 */

#include <clib.h>

#define NT_FORM_BSD_HIT  "nhbsd[0,400]:I,chbsd(nhbsd):I,Ebsd(nhbsd):R,\
                          ntbsd(nhbsd):I,t1bsd(nhbsd):I"
#define NT_MAX_BSD_HIT 400

typedef struct {
    int nhbsd;
    int chbsd[NT_MAX_BSD_HIT];
    float Ebsd[NT_MAX_BSD_HIT];
    int ntbsd[NT_MAX_BSD_HIT];
    int t1bsd[NT_MAX_BSD_HIT];
} nt_bsd_hits_t;

#define NT_FORM_BSD_PIX  "npix[0,200]:I,ipix(npix):I,\
                          rpix(npix):I,lpix(npix):I,spix(npix):I,\
                          zpix(npix):R,phipix(npix):R,\
                          rpixt(npix):R,lpixt(npix):R,spixt(npix):R,\
                          rpixE(npix):R,lpixE(npix):R,spixE(npix):R"
#define NT_MAX_BSD_PIX 200

typedef struct {
    int npix;
    int ipix[NT_MAX_BSD_PIX];
    int rpix[NT_MAX_BSD_PIX];
    int lpix[NT_MAX_BSD_PIX];
    int spix[NT_MAX_BSD_PIX];
    float zpix[NT_MAX_BSD_PIX];
    float phipix[NT_MAX_BSD_PIX];
    float rpixt[NT_MAX_BSD_PIX];
    float lpixt[NT_MAX_BSD_PIX];
    float spixt[NT_MAX_BSD_PIX];
    float rpixE[NT_MAX_BSD_PIX];
    float lpixE[NT_MAX_BSD_PIX];
    float spixE[NT_MAX_BSD_PIX];
} nt_bsd_pixs_t;

#define NT_FORM_BGV_HIT  "nhbgv[0,400]:I,chbgv(nhbgv):I,EbgvDwn(nhbgv):R,\
                          EbgvUp(nhbgv):R,ntbgvDwn(nhbgv):I,t1bgvDwn(nhbgv):I,\
                          ntbgvUp(nhbgv):I,t1bgvUp(nhbgv):I"
#define NT_MAX_BGV_HIT 400

typedef struct {
    int nhbgv;
    int chbgv[NT_MAX_BGV_HIT];
    float EbgvDwn[NT_MAX_BGV_HIT];
    float EbgvUp[NT_MAX_BGV_HIT];
    int ntbgvDwn[NT_MAX_BGV_HIT];
    int t1bgvDwn[NT_MAX_BGV_HIT];
    int ntbgvUp[NT_MAX_BGV_HIT];
    int t1bgvUp[NT_MAX_BGV_HIT];
} nt_bgv_hits_t;

#define NT_FORM_CPV_HIT  "nhcpv[0,250]:I,chcpv(nhcpv):I,Ecpv(nhcpv):R,\
                          ntcpv(nhcpv):I,t1cpv(nhcpv):I"
#define NT_MAX_CPV_HIT 250

typedef struct {
    int nhcpv;
    int chcpv[NT_MAX_CPV_HIT];
    float Ecpv[NT_MAX_CPV_HIT];
    int ntcpv[NT_MAX_CPV_HIT];
    int t1cpv[NT_MAX_CPV_HIT];
} nt_cpv_hits_t;

#define NT_FORM_UPV_HIT  "nhupv[0,50]:I,chupv(nhupv):I,Eupv(nhupv):R,\
                          ntupv(nhupv):I,t1upv(nhupv):I"
#define NT_MAX_UPV_HIT 50

typedef struct {
    int nhupv;
    int chupv[NT_MAX_UPV_HIT];
    float Eupv[NT_MAX_UPV_HIT];
    int ntupv[NT_MAX_UPV_HIT];
    int t1upv[NT_MAX_UPV_HIT];
} nt_upv_hits_t;

#define NT_FORM_TAG_HIT  "nhtag[0,320]:I,chtag(nhtag):I,Etag0(nhtag):R,\
                          Etag1(nhtag):R,nttagL(nhtag):I,nttagR(nhtag):I,\
                          t1tagL(nhtag):I,t1tagR(nhtag):I"
#define NT_MAX_TAG_HIT 320

typedef struct {
    int nhtag;
    int chtag[NT_MAX_TAG_HIT];
    float Etag0[NT_MAX_TAG_HIT];
    float Etag1[NT_MAX_TAG_HIT];
    int nttagL[NT_MAX_TAG_HIT];
    int nttagR[NT_MAX_TAG_HIT];
    int t1tagL[NT_MAX_TAG_HIT];
    int t1tagR[NT_MAX_TAG_HIT];
} nt_tag_hits_t;

#define NT_FORM_LGD_HIT  "nhlgd[0,640]:I,chlgd(nhlgd):I,Elgd(nhlgd):R,\
                          clust(nhlgd):I"
#define NT_MAX_LGD_HIT 640

typedef struct {
    int nhlgd;
    int chlgd[NT_MAX_LGD_HIT];
    float Elgd[NT_MAX_LGD_HIT];
    int clust[NT_MAX_LGD_HIT];
} nt_lgd_hits_t;

#define NT_FORM_TIME_LIST "ntimes[0,500]:I,le(ntimes):R"
#define NT_MAX_TIME_LIST 500

typedef struct {
    int ntimes;
    float le[NT_MAX_TIME_LIST];
} nt_time_list_t;

#define NT_FORM_PHOTONS "nphot[0,100]:I,nfrwd:I,Efrwd:R,pvect(4,nphot):R"
#define NT_MAX_PHOTONS 100

typedef struct {
    int nphot;
    int nfrwd;
    float Efrwd;
    vector4_t pvect[NT_MAX_PHOTONS];
} nt_photons_t;

#define NT_FORM_MESONS "nmes[0,500]:I,mtype(nmes):I,\
                        ptot(4,nmes):R,amass(nmes):R,\
                        idtype(2,nmes):I,ichild(2,nmes):I"
#define NT_MAX_MESONS 500

typedef struct {
    int nmes;
    int mtype[NT_MAX_MESONS];
    vector4_t ptot[NT_MAX_MESONS];
    float amass[NT_MAX_MESONS];
    int idtype[NT_MAX_MESONS][2];
    int ichild[NT_MAX_MESONS][2];
} nt_mesons_t;

#define NT_FORM_TAG_COIN "ncoin[0,30]:I,cochan(ncoin):I,\
                          cotime(ncoin):R,coenergy(ncoin):R,\
                          tagweight(ncoin):R"
#define NT_MAX_TAG_COIN 30

typedef struct {
    int ncoin;
    int cochan[NT_MAX_TAG_COIN];
    float cotime[NT_MAX_TAG_COIN];
    float coenergy[NT_MAX_TAG_COIN];
    float tagweight[NT_MAX_TAG_COIN];
} nt_tag_coin_t;

#define NT_FORM_RECOILS "nrec[0,100]:I,trec0:R,\
                         therec(nrec):R,phirec(nrec):R,dErec(nrec):R,\
                         Erec(nrec):R,trec(nrec):R,mrec(nrec):I"
#define NT_MAX_RECOILS 100

typedef struct {
    int nrec;
    float trec0;
    float therec[NT_MAX_RECOILS];
    float phirec[NT_MAX_RECOILS];
    float dErec[NT_MAX_RECOILS];
    float Erec[NT_MAX_RECOILS];
    float trec[NT_MAX_RECOILS];
    int mrec[NT_MAX_RECOILS];
} nt_recoils_t;

#define NT_FORM_BCLUSTERS "nbcl[0,24]:I,bcE(nbcl):R,\
                           bcphi(nbcl):R,bcz(nbcl):R,bct(nbcl):R,\
                           bcsE(nbcl):R"
#define NT_MAX_BCLUSTERS 24


typedef struct {
    int nbcl;
    float bcE[NT_MAX_BCLUSTERS];
    float bcphi[NT_MAX_BCLUSTERS];
    float bcz[NT_MAX_BCLUSTERS];
    float bct[NT_MAX_BCLUSTERS];
    float bcsE[NT_MAX_BCLUSTERS];
} nt_bclusters_t;

#define NT_FORM_TAG_ALIGN "taref:R,ntal[0,250]:I,chtal(ntal):I,tal(ntal):R"
#define NT_MAX_TAG_ALIGN 250

typedef struct {
    float taref;
    int ntal;
    int chtal[NT_MAX_TAG_ALIGN];
    float tal[NT_MAX_TAG_ALIGN];
} nt_tag_align_t;

void appendnt_bsd_hits(int id, nt_bsd_hits_t *p_bsd);
void appendnt_bsd_pixs(int id, nt_bsd_pixs_t *p_pix);
void appendnt_bgv_hits(int id, nt_bgv_hits_t *p_bgv);
void appendnt_cpv_hits(int id, nt_cpv_hits_t *p_cpv);
void appendnt_upv_hits(int id, nt_upv_hits_t *p_upv);
void appendnt_lgd_hits(int id, nt_lgd_hits_t *p_lgd);
void appendnt_tag_hits(int id, nt_tag_hits_t *p_tag);
void appendnt_time_list(int id, nt_time_list_t *p_tlist);
void appendnt_photons(int id, nt_photons_t *p_clust);
void appendnt_mesons(int id, nt_mesons_t *p_mes);
void appendnt_recoil(int id, nt_recoils_t *p_rec);
void appendnt_bclusters(int id, nt_bclusters_t *p_bclust);
void appendnt_tag_coin(int id, nt_tag_coin_t *p_coin);
void appendnt_tag_align(int id, nt_tag_align_t *p_talign);

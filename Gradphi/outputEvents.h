#define GRADPHI_VERSION_RELEASE 200
#define BUFSIZE 100000

/* GEANT external function for fetching hits;
 * no use to call this function from guout itself
 * because it is necessary to know all the arrays 
 * for saving hits anyway, and call of gfhits from 
 * here will prevent passing huge array into a
 * subroutine
*/

void gfhits_(char *chset, char *chdet, int *nvdim, int *nhdim,
             int *nhmax, int *itrs, int *numvs, 
             int *itra, int *numbv, float *hits,
             int *nhits, int len1,int len2);
             
/* This function returns a variable distributed according to Poisson
 * distribution with the mean value mean
*/
void gpoiss_(float *mean, int *value, int *number);

/* This function will return the pointer to the array of user words saved 
 * with the corresponding vertex; in this case this array is mc_event 
 * structure
*/
void gfvert_(int *nvtx,float *vert,int *ntbeam,int *nttarg,float *tofg,
             float *ubuf,int *nubuf);

int LGDsave(itape_header_t *event,
            int nvdimlgd,int nhdimlgd,int nhmax);
int CPVsave(itape_header_t *event,int nvdimcpv,int nhdimcpv,
            int nhmax,int *index,time_list_t *time_list);
int RPDsave(itape_header_t *event,int nvdimrpd,int nhdimrpd,
            int nhmax,int *index,time_list_t *time_list);
int UPVsave(itape_header_t *event,int nvdimupv,int nhdimupv,
            int nhmax,int *index,time_list_t *time_list);
int BSDsave(itape_header_t *event,int nvdimupv,int nhdimupv,
            int nhmax,int *index,time_list_t *time_list);
int BGVsave(itape_header_t *event,
            int nvdimbgv,int nhdimbgv,int nhmaxbgv,
            int *index,time_list_t *time_list);
void sort(float *hits,float *sorthits,int dimension,int nhits);
void timeAnalysis (int *channel,int chdim,int *numbv,float *energy, 
                   int *times,float *timelist,int listdim,float *sorthits,
                   int nhits,int key,int nvdim,int nhdim);
void timeAnalysisBSD (int *channel,int chdim,int *numbv,float *energy,
                   int *times,float *timelist,int listdim,float *sorthits,
                   int nhits,int key,int nvdim,int nhdim);
void timeAnalysisCPV (int *numbv,float *energy,
                   int *times,float *timelist,int listdim,float *sorthits,
                   int nhits,int key,int nvdim,int nhdim);
float extract_bgv_t0_ (float *hits, int *end);

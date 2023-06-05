/*
 * disData.h
 *
*/

#ifndef disDataH_INCLUDED
#define disDataH_INCLUDED


#include <rectypes.h>

/********************************************

An event will look like this:

---------------------------------------------
event      header              itape_header_t
           ----------------------------------
           group1              group_header_t
                               data
           ----------------------------------
           group2              group_header_t
                               data
           ----------------------------------
	   ...
           ----------------------------------
           groupN              group_header_t
                               data
---------------------------------------------

Assigned group numbers:

    RANGE           SUBSYSTEM                    

   0 -   99 - Raw data, internal use
 100 -  199 - RPD,BSD,BGV data structures                
 200 -  299 - LGD data structures                
 300 -  399 - Tagger data structures             
 400 -  499 - CPV and UPV data structures
 500 -  599 - Trigger data structures
 600 -  699 - Misc.

**********************************************************************/

/*
 * Note: This file is used to automatically generate little <--> big
 *       endian data conversion routines. For each data group,
 *       it needs to know the data type and the include file
 *       where the data type is defined.
 *
 * Therefore each data group should be defined using the following
 *       format:
 *
 * #define GROUP_XXX        NUM  ** header.h:type_t - user comments **
 *
 *
 * Some data groups come in more than one structure shape, distinguished
 * by sub-group index 0,1...  This allows groups to have an internal
 * hierarchical structure, or to have their structures evolve over time
 * and keep the software backward-compatible.  For data groups with more
 * than one valid sub-group index, a comma-separated list of <type_t>
 * entries should appear in the comment section shown above, as in
 *
 * #define GROUP_XXX        NUM  ** header.h:type0_t,type1_t[,...] **
 *
 * This requires that all of the sub-group structure definitions are kept
 * in the same header file, which is the logical thing anyway.
 *
 * There exist "special" predefined types and header files that can be used
 *       when there is no "real" data type or header file corresponding
 *       to the data group:
 *
 *   NOCONVERT_t --- a data type that does not need endian conversion,
 *                        for example an array of char.
 *   INTARRAY_t  --- a data type that should be converted as an array
 *                        of 32-bit integers, for example the raw data group.
 *
*/

#define GROUP_NOGROUP        -1  /* IGNORE */
#define GROUP_DELETED         0  /* groups deleted with data_removeGroup()    disData.h:NOCONVERT_t */
#define GROUP_RAW             1  /* raw data group, see iitypes.h INTARRAY_t */
#define GROUP_FROZEN         11  /* itypes.h:frozen_groups_t */

#define GROUP_RPD_ADCS      110  /* itypes.h:adc_values_t */
#define GROUP_RPD_TDCS      111  /* itypes.h:tdc_values_t */
#define GROUP_RPD_HITS      112  /* iitypes.h:rpd_hits_t */
#define GROUP_RPD_COINS     113  /* iitypes.h:rpd_coins_t */

#define GROUP_BSD_ADCS      120  /* itypes.h:adc_values_t */
#define GROUP_BSD_TDCS_LONG     121  /* itypes.h:tdc_values_t */
#define GROUP_BSD_TDCS_SHORT    122  /* itypes.h:tdc_values_t */
#define GROUP_BSD_HITS       123  /* iitypes.h:bsd_hits_t */
#define GROUP_BSD_PIXELS     124  /* iitypes.h:bsd_pixels_t */
#define GROUP_BSD_TDCS GROUP_BSD_TDCS_LONG

#define GROUP_BGV_ADCS      130  /* itypes.h:adc_values_t */
#define GROUP_BGV_TDCS      131  /* itypes.h:tdc_values_t */
#define GROUP_BGV_HITS      132  /* iitypes.h:bgv_hits_t */
#define GROUP_BGV_COINS     133  /* iitypes.h:rpd_coins_t */
#define GROUP_BGV_CLUSTERS  134  /* libBGV.h:bgv_clusters_t */

#define GROUP_LGD_ADCS      210  /* itypes.h:adc_values_t */
#define GROUP_LGD_HITS      211  /* lgdCluster.h:lgd_hits_t */
#define GROUP_ESUM_ADCS     212  /* itypes.h:adc_values_t */
#define GROUP_MAM_REGS      213  /* itypes.h:reg_values_t */
#define GROUP_LGD_CLUSTERS      222  /* lgdCluster.h:lgd_clusters_t */
#define GROUP_LGD_CLUSTER_HITS  223  /* lgdCluster.h:lgd_hits_t  */
#define GROUP_PHOTONS           230  /* iitypes.h:photons_t  */
#define GROUP_MESONS            231  /* iitypes.h:mesons_t  */

#define GROUP_LGD_RECONSTRUCT 250 /* lgdCluster.h:lgd_reconstruct_t */

#define GROUP_TAGGER_TDCS         310 /* itypes.h:tdc_values_t */
#define GROUP_TAGGER_HITS         311 /* iitypes.h:tagger_hits_t */
#define GROUP_TAGGER_CLUSTERS     312 /* iitypes.h:tagger_clusters_t */
#define GROUP_TAGGER_PHOTONS      313 /* iitypes.h:tagger_photons_t */

#define GROUP_CPV_ADCS      410   /* itypes.h:adc_values_t */
#define GROUP_CPV_TDCS      411   /* itypes.h:tdc_values_t */
#define GROUP_CPV_HITS      412   /* iitypes.h:cpv_hits_t */

#define GROUP_UPV_ADCS      450   /* itypes.h:adc_values_t */
#define GROUP_UPV_TDCS      451   /* itypes.h:tdc_values_t */
#define GROUP_UPV_HITS      452   /* iitypes.h:upv_hits_t */

#define GROUP_SCALERS       610   /* itypes.h:scaler_values_t */
#define GROUP_TIME_LIST     611   /* iitypes.h:time_list_t */
#define GROUP_MC_EVENT      612   /* mctypes.h:mc_event_t,mc_vertex_t */

#define GROUP_RECOIL        710   /* iitypes.h:recoil_t */
#define GROUP_3_PHOTON_FITS 723   /* iitypes.h:three_photon_fits_t */
#define GROUP_3_PHOTON_ERR_MATRICES 733 /* IGNORE */
#define GROUP_ESR_NPARTICLE 802 /*esr.h:esr_nparticle_t */
#ifdef AIX
#define TRANSCOMPUTERCODE 1
#endif
#ifdef IRIX
#define TRANSCOMPUTERCODE 1
#endif
#ifdef DEC
#define TRANSCOMPUTERCODE 2
#endif
#ifdef SOLARIS
#define TRANSCOMPUTERCODE 1
#endif
#ifdef HPUX
#define TRANSCOMPUTERCODE 1
#endif
#ifdef INTEL_LINUX
#define TRANSCOMPUTERCODE 2
#endif
#ifdef PPC_LINUX
#define TRANSCOMPUTERCODE 1
#endif 
/*  note on ppc-linux ... ppc-linux runs on powerPC chips which are big-endian */
#ifndef TRANSCOMPUTERCODE
#define TRANSCOMPUTERCODE -1
#endif

#include <ntypes.h>

typedef struct
{
  uint32   length;       /* record length in bytes      */
  uint32   type;         /* should always be TYPE_ITAPE */
  uint32   transComputerCode;
  uint32   ngroups;
   int32   runNo;
  uint32   spillNo;
   int32   eventNo;      /* 0-BOS, (-1)-EOS, others- events number */
  uint16   eventType;    /* event type (decoded latch word)        */
   int16   trigger;      /* trigger type */
   int32   time;         /* time as returned by time()             */
  uint32   latch;        /* virtual latch word                     */
} itape_header_t;

typedef struct
{
  uint32   length;
  uint32   type;
  uint32   subtype;
} group_header_t;


/*
 * 'data_newItape' will fill initialize an empty itape
*/

int data_newItape(itape_header_t *event);

/*
 * 'data_getGroup' returns the pointer to the group data or NULL
*/

void* data_getGroup(const itape_header_t *event, uint32 group, uint32 subgroup );

size_t data_getGroupSize(const itape_header_t *event, uint32 group, uint32 subgroup);

/*
 * 'data_getGroupHeader' returns the pointer to the group header or NULL
*/

group_header_t* data_getGroupHeader(const itape_header_t *event,uint32 group, uint32 subgroup);

void* data_addGroup(itape_header_t *event, size_t bufsize,
		    uint32 group, uint32 subgroup, size_t length);

group_header_t* data_addGroup_header(itape_header_t *event,size_t bufsize,
				     uint32  group, uint32 subgroup, size_t length);

/*
 * 'data_removeGroup' returns 0 if group found and removed, 1 otherwise
 *    the space occupied by the removed group is not deallocated and is lost
 *    until 'data_clean' is called to compact the itape
*/

int data_removeGroup(itape_header_t *event, uint32 group, uint32 subgroup);


int data_removeAllgroups(itape_header_t *event);  /* returns 1 */
int data_saveGroups(itape_header_t *event,int nsave,int *isave, int *sisave);  /* returns 1 */

/*
 * returns the address of the end of used space inside the buffer
*/

void* data_getAdr(const itape_header_t *event);

/*
 * 'data_clean' removes the unused parts of the itape created by the 'removeGroup' functions
 *
 *    NOTICE: this call moves data inside the itape and therefore all pointers
 *            into the itape (for example, pointers returned by 'getGroup')
 *            become invalid.
*/

int data_clean(itape_header_t *event);

/*
 * Calculates a new CRC for the itape
*/

int data_addCRC(itape_header_t *event);

/*
 * Check the CRC. If it is broken, return 1, otherwise, return 0
*/

int data_checkCRC(const itape_header_t *event);
int data_checkCRC1(const itape_header_t *event,size_t bufferLength);

/*
 * perform big <--> little endian conversion
*/

int endian_convertItape(itape_header_t *event,int itapeLength);

/*
 * The 'data_listGroups' function will return-
 *    in 'ngroups' - number of groups in the event
 *    in 'groupsList' - a pointer to the array of groups
 *    in 'groupsSize' - a pointer to the array of their lengths
 *
 * All the three arguments can be NULL if the returned information is not needed
 * The arrays should be free'ed after use
*/

int data_listGroups(const itape_header_t *event,int*ngroups,uint32  **groupsList, uint32 **subgroupsList, size_t **groupsSize);

#endif
/* end file */

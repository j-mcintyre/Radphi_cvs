/*
 * tapeData.h --- format of data written on tapes
 *
*/

#include <rectypes.h>
#include <ntypes.h>

/*
   our type will be TYPE_TAPEHEADER and TYPE_TAPEFOOTER
*/

typedef struct
{
  uint32 length;   /* in bytes */
  uint32 type;     /* should be TYPE_TAPEHEADER */
  uint32 transComputerCode; /* here just in case..., should be 1 */
  uint32 ngroups;  /* should be 0 (for compatibility with itape records) */
  uint32 tapeSeqNumber;
  uint32 startTime;
  char devname[256];        /* written on this device */
  char from[256];           /* data from ... */
  char crc[4];
} tapeHeader_t;

/* end file */

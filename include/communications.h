/*
 * communications.h
 *
*/

#include <disIO.h>

#define TransComputerCode 1 /* defines order of bytes for words--
	Silicon Graphics is defined to be 1, Vax is 0 */

#define DATA_RECORD 1
#define COMMAND     2
#define ITAPE       3
#define ITAPE_RECORD 4
#define TIMEOUT	    5

#define MAP_REQUEST 6

#define EIGHT_MM	9	

typedef struct {
    int transcomputer;
    int nbytes;
    int type;
} DatagramHeader_t;


typedef struct /* NOTE: all integers are 'network' endian (see 'man ntohl') */
{
  int requestCode; /* from the MAP_REQUEST_XXX set  */
  int dataType;    /* from the MAP_DATATYPE_XXX set */
  int dataCount;   /* number of 'dataType' items (array length) */

  int accessKey;   /* map manager access key */

  int    systemName_offset;
  int subsystemName_offset;
  int      itemName_offset;
  int          data_offset;
  
  char data[1];
} map_request_t;

/* end file */

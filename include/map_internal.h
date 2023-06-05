/* 
   map_internal.h        declarations for internal map routines
   Created           : 19-OCT-1992 by Thom Sulanke
*/
 
#define VERSION 2      /* data base format version */
#define ROOT  0        /* offset in file of root of tree */

#ifdef UNDEF

/* The IBM XL (eXtra Large) C compiler does not
   understand this
*/

#ifndef _SIZE_T_
#define _SIZE_T_
typedef unsigned int size_t;
#endif /* !_SIZE_T_ */

#endif

#ifndef vaxc
#include <sys/types.h> /* to get the correct size_t */
#else
#include <types.h> /* to get the correct size_t */
#endif

typedef long pointer_t;

#define NULL_LOC 0    /* null pointer within file */

typedef struct subtableheader {
  int negversion;      /* cause version 1 software to fail */
  int nullname;
  int table_length;
  int table_used;
  pointer_t table;
} subtableheader_t;

typedef struct nametable {
  char name[MAP_NAMELENGTH];
  pointer_t loc;
} nametable_t;

#define INIT_NAMETABLE_SIZE 10

typedef struct subsystem {
  pointer_t next;
  char name[MAP_NAMELENGTH];
  pointer_t list;
  pointer_t table;
  int table_length;
  int table_used;
} subsystem_t;

#define SUBSYSTEM_LIST_OFFSET 20

typedef struct item {
  pointer_t next;
  char name[MAP_NAMELENGTH];
  pointer_t list;
  int length;
  int type;
  pointer_t table;
  int table_length;
  int table_used;
  long spare[1];
} item_t;

#define ITEM_LIST_OFFSET 20
#define INIT_TABLE_SIZE 50
#define TABLE_FACTOR 5

typedef struct arrayheader {
  pointer_t next;
  int time;
  int length;
  long spare[1];
} arrayheader_t;

#ifdef __LP64__
#  define ARRAY_VALUES_OFFSET 24
#else
#  define ARRAY_VALUES_OFFSET 16
#endif

typedef struct name {
  pointer_t next;
  char name[MAP_NAMELENGTH];
} name_t;

typedef struct table {
  int time;
  pointer_t loc;
} table_t;

/* function prototypes */
int 	    map_add_item1(const char filename[], const char subsystemname[], 
			  const char itemname[], int arraylength, 
			  int arraytype, int table_size);
int         map_add_sub(int fd, const char subsystemname[], int nitems);
int         map_close_map( int fd );
int         map_dump( const char filename[] );
void        map_error(int error_code, char *fmt, ... );
void        map_fillname(char to[], const char from[]); 
int   map_find_item(int fd, const char subsystemname[], 
			  const char itemname[],pointer_t *adr );
int   map_find_pred(int fd, pointer_t start, const char name[], 
			  int table_used, 
			  pointer_t table_loc, int *tindex, 
			  pointer_t *succ_loc,pointer_t *adr);
int    map_find_same(int fd, pointer_t start, const char name[], 
			  int table_used, 
			  pointer_t table_loc, int *tindex, 
			  pointer_t *succ_loc,pointer_t *adr);
int   map_find_sub(int fd, const char subsystemname[],pointer_t *adr);
pointer_t   map_find_time_pred(int fd, pointer_t start, int atime);
int 	    map_fix_head(const char filename[], int nsubsystems); 
int 	    map_fix_item(const char filename[], const char subsystemname[], 
			 const char itemname[], int nvals); 
int 	    map_fix_sub(const char filename[], const char subsystemname[], 
			int nitems); 
int         map_open_ro(const char filename[]);
int         map_open_ro_lock(const char filename[]);
int         map_open_rw(const char filename[]);
int         map_overwrite(const void *chunck, size_t nbytes, int fd, 
			  pointer_t offset);
int         map_read(void *chunck, size_t nbytes, int fd, pointer_t offset);
int         map_read_char(char s[], size_t nbytes, int fd, pointer_t offset);
int         map_read_float(float f[], size_t nbytes, int fd, pointer_t offset);
int         map_read_int(int i[], size_t nbytes, int fd, pointer_t offset);
int         map_read_pointer(pointer_t *pointer,
                             size_t nbytes, int fd, pointer_t offset);
int         map_read_name(name_t *name,
                             size_t nbytes, int fd, pointer_t offset);
int         map_read_table(table_t *table,
                             size_t nbytes, int fd, pointer_t offset);
int         map_read_nametable(nametable_t *nametable,
                             size_t nbytes, int fd, pointer_t offset);
int         map_read_arrayheader(arrayheader_t *arrayheader,
                             size_t nbytes, int fd, pointer_t offset);
int         map_read_subsystem(subsystem_t *subsystem,
                             size_t nbytes, int fd, pointer_t offset);
int         map_read_subtableheader(subtableheader_t *subtableheader,
                             size_t nbytes, int fd, pointer_t offset);
int         map_read_item(item_t *item,
                             size_t nbytes, int fd, pointer_t offset);
int         map_write_char(const char s[],
                             size_t nbytes, int fd, pointer_t *adr);
int         map_write_float(const float f[],
                             size_t nbytes, int fd, pointer_t *adr);
int         map_write_int(const int i[], 
                             size_t nbytes, int fd, pointer_t *adr);
int         map_write_pointer(const pointer_t *pointer,
                             size_t nbytes, int fd, pointer_t *adr);
int         map_write_name(const name_t *name,
                             size_t nbytes, int fd, pointer_t *adr);
int         map_write_table(const table_t *table,
                             size_t nbytes, int fd, pointer_t *adr);
int         map_write_nametable(const nametable_t *nametable,
                             size_t nbytes, int fd, pointer_t *adr);
int         map_write_arrayheader(const arrayheader_t *arrayheader,
                             size_t nbytes, int fd, pointer_t *adr);
int         map_write_subsystem(const subsystem_t *subsystem,
                             size_t nbytes, int fd, pointer_t *adr);
int         map_write_subtableheader(const subtableheader_t *subtableheader,
                             size_t nbytes, int fd, pointer_t *adr);
int         map_write_item(const item_t *item,
                             size_t nbytes, int fd, pointer_t *adr);
int         map_ovewrite_char(const char s[],
                             size_t nbytes, int fd, pointer_t offset);
int         map_ovewrite_float(const float f[],
                             size_t nbytes, int fd, pointer_t offset);
int         map_ovewrite_int(const int i[],
                             size_t nbytes, int fd, pointer_t offset);
int         map_ovewrite_pointer(const pointer_t *pointer,
                             size_t nbytes, int fd, pointer_t offset);
int         map_ovewrite_name(const name_t *name,
                             size_t nbytes, int fd, pointer_t offset);
int         map_ovewrite_table(const table_t *table,
                             size_t nbytes, int fd, pointer_t offset);
int         map_ovewrite_nametable(const nametable_t *nametable,
                             size_t nbytes, int fd, pointer_t offset);
int         map_ovewrite_arrayheader(const arrayheader_t *arrayheader,
                             size_t nbytes, int fd, pointer_t offset);
int         map_ovewrite_subsystem(const subsystem_t *subsystem,
                             size_t nbytes, int fd, pointer_t offset);
int         map_ovewrite_subtableheader(const subtableheader_t *subtableheader,
                             size_t nbytes, int fd, pointer_t offset);
int         map_ovewrite_item(const item_t *item,
                             size_t nbytes, int fd, pointer_t offset);
int 	    map_stat_head(const char filename[], int *nsubsystems); 
void        map_syserror(int error_code, char *fmt, ... );
int	    map_time_pred(int fd, pointer_t pointer, item_t item, int atime,
			  int *tindex, pointer_t *suc_loc, int *suc_time,pointer_t *adr);
void        map_warn(int error_code, char *fmt, ... );
int         map_write(const void *chunck, size_t nbytes, int fd,pointer_t *adr);

/* libdata function for little/big-endian conversions */

#include <ntypes.h>
extern int endian_convertTypeN_uint32(uint32*dptr,int dcount);

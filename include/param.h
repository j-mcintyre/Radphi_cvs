/*
 * param.h - parameters database
*/

#ifndef include_paramH
#define include_paramH

#ifndef __CINT__
static const char sccsid_paramH[] = "@(#)param.h\t4.1\tCreated 10 Dec 1995 18:48:48, \tcompiled "__DATE__;
#endif

#include <ntypes.h>

/* set the verbose flag (default = 1) */
void SetMapVerbose(int);

/*
 * clear the database
*/


int param_clean(void);

/*
 * load a '.db' file
*/

int param_loadFile(const char*file);

/*
 * Load the database for given run:
 *
 * Database files are looked-up in the following order:
 *   1) E852_CONFIG/run.0/NAME.db
 *   2) E852_CONFIG/run.RUNNO/NAME.db  (where RUNNO is the run number, i.e. run.6311/space.db)
 *   3) Map manager file E852_CONFIG/maps/db-files.Map, subsystem 'db-files',
 *                       item 'NAME'.
*/

int param_loadDatabase(const char*name, int run);

/* The current database is saved somwhere in the E852_CONFIG area */

int param_saveDatabase(int run);

int param_putValue(const char*name,const char*value);

/*
 * Check if the requested parameter has a value,
 * return 1 if yes, 0 if parameter not found.
 *
 * Note: There should be no ***silent*** param_getValue(...)
 *       A silent param_getValue(...) can be simulated
 *       by a param_valueExists(...)/param_getValue(...) pair.
*/

int param_valueExists(const char *name);


const char* param_getValue(const char*name);
const char* param_getValueI(const char*name,int i);

const char* param_getString(const char*name,const char*deflt);
const char* param_getStringI(const char*name,int i,const char*deflt);

unsigned long param_getUlong(const char*name,unsigned long deflt);
unsigned long param_getUlongI(const char*name,int i,unsigned long deflt);

uint32 param_get_uint32(const char*name,uint32 deflt);
uint32 param_getI_uint32(const char*name,int i,uint32 deflt);

int32 param_get_int32(const char*name,int32 deflt);
int32 param_getI_int32(const char*name,int i,int32 deflt);

int param_getInt(const char*name,int deflt);
int param_getIntI(const char*name,int i,int deflt);

double param_getDouble(const char*name,double deflt);
double param_getDoubleI(const char*name,int i,double deflt);

const vector3_t* param_get3vector(const char*name,const vector3_t* deflt);
const vector3_t* param_get3vectorI(const char*name,int i,const vector3_t* deflt);

const vector4_t* param_get4vector(const char*name,const vector4_t* deflt);
const vector4_t* param_get4vectorI(const char*name,int i,const vector4_t* deflt);

int param_getIntArray(const char*name,int ival[],int nval);
int param_getDblArray(const char*name,double dval[],int nval);

const char* param_FSII(const char*fmt,const char*s,int v1,int v2);
const char* param_FII(const char*fmt,int v1,int v2);
const char* param_FI(const char*fmt,int v);
const char* param_FS(const char*fmt,const char*s);

/*
 * These are the Map Manager interface routines
*/

int param_createMMfile(const char* sysName);

int   param_getMM_int  (int runNo,const char*sysName,const char*subsysName,const char*itemName);
float param_getMM_float(int runNo,const char*sysName,const char*subsysName,const char*itemName);

int   param_getMM1_int  (int runNo,const char*sysName,const char*subsysName,const char*itemName,int deflt,int*erunNo);
float param_getMM1_float(int runNo,const char*sysName,const char*subsysName,const char*itemName,float deflt,int*erunNo);

int param_putMMarray_float(int runNo,const char*sysName,const char*subsysName,const char*itemName,int size,const float array[]);
int param_getMMarray_float(int runNo,const char*sysName,const char*subsysName,const char*itemName,int size,float array[]);
int param_getMM1array_float(int runNo,const char*sysName,const char*subsysName,const char*itemName,int size,float array[],int* erunNo);

int param_putMMarray_int(int runNo,const char*sysName,const char*subsysName,const char*itemName,int size,const int array[]);
int param_getMMarray_int(int runNo,const char*sysName,const char*subsysName,const char*itemName,int size,int array[]);
int param_getMM1array_int(int runNo,const char*sysName,const char*subsysName,const char*itemName,int size,int array[],int* erunNo);

/* NOTE: *size returns size of the string,
         memory region pointed to by *array should be free()ed.
*/

int param_getMMarray_string(int runNo,const char*sysName,const char*subsysName,const char*itemName,
			    int *size,char *array[]);
int param_getMM1array_string(int runNo,const char*sysName,const char*subsysName,const char*itemName,
			     int *size,char *array[],int* erunNo);

/* get the item array length (in bytes for strings) */
int param_getMMarray_size(int runNo,const char*sysName,const char*subsysName,const char*itemName);

#endif
/* end file */

/*
 * paramMM.c --- Map Manager interface routines
 *
*/

static const char rcsid[]="$Id";
/* 
$Log: paramMM.c,v $
Revision 1.2  1997/03/06 15:35:11  lfcrob
Fixed for RADPHI

 * Revision 1.1  1996/11/26  18:20:09  lfcrob
 * Initial revision
 *
*/

#include <stdio.h>
#include <stdlib.h>

#include <param.h>
#include <map_manager.h>


int paramMM_doPrintTime = 0;   /* toggle this to enable/disable map manager timing */

#define MAPMANAGER_INT    0
#define MAPMANAGER_FLOAT  1
#define MAPMANAGER_STRING 2

static int iabs(int i) { if (i>=0) return i; else return -i; }

static const char *filename(const char* sysName)
{
  static const char *env = NULL;
  static char fname[1024];

  if (env==NULL) env = getenv("RADPHI_CONFIG");
  if (env==NULL)
    {
      fprintf(stderr,"param_getMM_filename: Environment variable RADPHI_CONFIG is not set. Exiting...\n");
      exit(1);
    }

  sprintf(fname,"%s/maps/%s.map",env,sysName);

  return fname;
}

int param_getMM1array_float(int runNo,
			    const char* sysName,const char *subsysName,const char* itemName,
			    int size,float array[],int *erunNo)
{
  int err = 0;
  int firsttime;


  err = map_get_float(filename(sysName),subsysName,itemName,size,array,iabs(runNo),&firsttime);
  if (err != MAP_OK) return -1;

  if (firsttime < 0) err = (-1);  /* map manager error */

  if (erunNo)
    *erunNo = firsttime;

  return err;
}

int param_getMMarray_float(int runNo,
			   const char* sysName,const char *subsysName,const char* itemName,
			   int size,float array[])
{
  return param_getMM1array_float(runNo,sysName,subsysName,itemName,size,array,NULL);
}

float param_getMM1_float(int runNo,
			 const char* sysName,const char *subsysName,const char* itemName,
			 float deflt,
			 int *erunNo)
{
  int err = 0;
  int firsttime;
  float val = deflt;


  err = map_get_float(filename(sysName),subsysName,itemName,1,&val,iabs(runNo),&firsttime);
  if (err != MAP_OK)
    {
      if (erunNo) *erunNo = 0;
      return deflt;
    }


  if (firsttime < 0) val = deflt;  /* map manager error */

  if (erunNo) *erunNo = firsttime;

  return val;
}

float param_getMM_float(int runNo,const char* sysName,const char *subsysName,const char* itemName)
{
  return param_getMM1_float(runNo,sysName,subsysName,itemName,0,NULL);
}

int param_getMM_int(int runNo,const char* sysName,const char *subsysName,const char* itemName)
{
  return param_getMM1_int(runNo,sysName,subsysName,itemName,0,NULL);
}

int param_getMM1_int(int runNo,const char* sysName,const char *subsysName,const char* itemName,int deflt,int*erunNo)
{
  int err = 0;
  int firsttime;
  int val = 0;

  runNo = iabs(runNo);


  err = map_get_int(filename(sysName),subsysName,itemName,1,&val,runNo,&firsttime);
  if (err != MAP_OK) return 0;


  if (firsttime < 0) val = deflt;  /* map manager error */

  if (erunNo)
    *erunNo = firsttime;

  return val;
}

int param_getMMarray_string(int runNo,
			    const char* sysName,const char *subsysName,const char* itemName,
			    int *size,char *array[])
{
  return param_getMM1array_string(runNo,sysName,subsysName,itemName,size,array,NULL);
}

int param_getMM1array_string(int runNo,
			    const char* sysName,const char *subsysName,const char* itemName,
			    int *size,char *array[],int *erunNo)
{
  int err = 0;
  int firsttime;
  int maxsize = 10*1024;

  runNo = iabs(runNo);


  if (erunNo) *erunNo = 0;

  while (1) /* the stupid map manager cannot tell us the size of the array,
	       but it tells us (indirectly), when there is not enough space
	       by setting *size equal to maxsize.

	       Other routines also expect that the buffer is larger than *size
	       by at least one character.
	     */
    {
      *array = malloc(maxsize);

      err = map_get_char(filename(sysName),subsysName,itemName,maxsize,*array,runNo,&firsttime,size);
      if (err != MAP_OK)
	{
	  free(*array);
	  *array = NULL;
	  return -1;
	}

      if (*size < maxsize) break; /* here we make sure that the buffer is
				     bigger than data
				   */

      free(*array);
      *array = NULL;

      maxsize *= 2;

/*
      fprintf(stderr,"param.getMMarray_string: Notice: Increasing size to %d for %s/%s/%s\n",
	      maxsize,sysName,subsysName,itemName);
*/
    }


  if (firsttime < 0) err = (-1);  /* map manager error */

  if (erunNo) *erunNo = firsttime;

  return err;
}

int param_getMM1array_int(int runNo,
			  const char* sysName,const char *subsysName,const char* itemName,
			  int size,int array[],int *erunNo)
{
  int err = 0;
  int firsttime;

  runNo = iabs(runNo);


  err = map_get_int(filename(sysName),subsysName,itemName,size,array,runNo,&firsttime);
  if (err != MAP_OK) return -1;


  if (firsttime < 0) err = (-1);  /* map manager error */

  if (erunNo)
    *erunNo = firsttime;

  return err;
}

int param_getMMarray_int(int runNo,
			 const char* sysName,const char *subsysName,const char* itemName,
			 int size,int array[])
{
  return param_getMM1array_int(runNo,sysName,subsysName,itemName,size,array,NULL);
}

int param_getMMarray_size(int runNo,
			 const char* sysName,const char *subsysName,const char* itemName)
{
  int length,type;
  int err = 0;

  err = map_stat_item(filename(sysName),subsysName,itemName,&length,&type);
  if (err != MAP_OK) return -1;

  return length;
}

int param_putMMarray_float(int runNo,
			   const char* sysName,const char *subsysName,const char* itemName,
			   int size,const float array[])
{
  int err = 0;
  int firsttime;

  runNo = iabs(runNo);

  firsttime = -runNo;

  err = map_put_float(filename(sysName),subsysName,itemName,size,array,firsttime);
  if (err != MAP_OK) 
    {
      if (err == MAP_USER_WARN_NOMATCH)
	{
	  err = map_add_item(filename(sysName),subsysName,itemName,size,MAPMANAGER_FLOAT);
	  if (err != MAP_OK) return -1;

	  err = map_put_float(filename(sysName),subsysName,itemName,size,array,firsttime);
	  if (err != MAP_OK) return -1;
	}
      else 
	return -1;
    }

  return 0;
}

int param_putMMarray_int(int runNo,
			 const char* sysName,const char *subsysName,const char* itemName,
			 int size,const int array[])
{
  int err = 0;
  int firsttime;

  runNo = iabs(runNo);

  firsttime = -runNo;

/* First, try to create the item (just in case it does not exist yet...) */

  err = map_put_int(filename(sysName),subsysName,itemName,size,array,firsttime);
  if (err != MAP_OK) 
    {
      if (err == MAP_USER_WARN_NOMATCH)
	{
	  err = map_add_item(filename(sysName),subsysName,itemName,size,MAPMANAGER_INT);
	  if (err != MAP_OK) return -1;

	  err = map_put_int(filename(sysName),subsysName,itemName,size,array,firsttime);
	  if (err != MAP_OK) return -1;
	}
      else 
	return -1;
    }

  return 0;
}

int param_createMMfile(const char* name)
{
  int ret;

  ret = map_create(filename(name));
  if (ret != MAP_OK)
    {
      fprintf(stderr,"param_createMMfile: map_create(%s) failed\n",name);
      return -1;
    }

  return 0;
}

/* end file */

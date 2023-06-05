/*
 * param.c - parameters database
*/

static const char rcsid[]="$Id";
/* 
$Log: param.c,v $
Revision 1.2  1997/03/06 15:34:58  lfcrob
Fixed for RADPHI

 * Revision 1.1  1996/11/26  18:20:04  lfcrob
 * Initial revision
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <X11/Xfuncproto.h>
#include <X11/Xlib.h>
#include <X11/Xresource.h>


#include <param.h>

int paramVerbose = 1;

static XrmDatabase paramdb = NULL;

static int iabs(int i) { if (i>=0) return i; else return -i; }

void SetMapVerbose(int i)
{
  paramVerbose = i;
}

int param_loadFile(const char*file)
{
  XrmDatabase new = XrmGetFileDatabase(file);
  if (new==NULL)
    {
      if (errno==ENOENT) return -1;
      fprintf(stderr,"param.loadFile: Error loading file %s: %s",file,strerror(errno));
      /*esend(ALARM_WARNING,"param_loadFile: Cannot load file","%s",file);*/
      return -1;
    }
  XrmMergeDatabases(new,&paramdb);
  return 0;
}

int param_clean(void)
{
  if (paramdb!=NULL) XrmDestroyDatabase(paramdb);
  paramdb = NULL;
  return 0;
}

int param_loadIncludes(int run)
{
  char flist[256];
  const char *ptr,*ptr1;
  const char*iptr[20];
  int n,i;
  int ret, ret1 = 0;

  run = iabs(run);

  if (!param_valueExists("include")) /* nothing to include */
    return 0;

  ptr = param_getValue("include");
  if (!ptr) return 0;

  strncpy(flist,ptr,sizeof(flist)-2);

  param_putValue("include"," ");

  ptr = strtok(flist," \t");
  n = 0;
  while (ptr)
    {
      iptr[n] = ptr;
      n++;
      if (n>19) break;
      ptr = strtok(NULL," \t");
    }

  for (i=0; i<n; i++)
    {
      ret = param_loadDatabase(iptr[i],run);
      if (ret) ret1 = (-1);
    }
  
  return ret1;
}

static void padPrintf(int len,const char*str)
{
  int lll = len - strlen(str);

  fprintf(stderr,"'%s'",str);

  while (lll-- > 0) fputc(' ',stderr);
}

int param_loadDatabase(const char*name,int runNo)
{
  char fname[256];
  const char* dir = getenv("RADPHI_CONFIG");
  int ret;

  runNo = iabs(runNo);

  if (dir==NULL)
    {
      fprintf(stderr,"param_load: Error: Environment variable RADPHI_CONFIG is not set.\n");
      return -1;
    }
  if (paramVerbose) {
    fprintf(stderr,"param.loadDatabase: Loading database ");
    padPrintf(16,name);
    fprintf(stderr,"\tfor run %d...",runNo);
  }

  while (1)
    {
  
      /* Try to load from $RADPHI_CONFIG/run.0 directory */

      sprintf(fname,"%s/run.%d/%s.xrdb",dir,0,name);
      
      if (param_loadFile(fname)==0)
	{
	  if (paramVerbose)
	    fprintf(stderr,"from run.%d/%s.xrdb\n",0,name);
	  ret = 0;
	  break;
	}
  
      /* Try to load from the Map Manager file */

      {
	int t = runNo;
	int size = 0;
	char *file = NULL;
	int erunNo = 0;
	
	if (t==0) t = 99999; /*
			       "default" or "current" run 0 corresponds to time 99999
			       in the map manager file
			       */
	
	ret = param_getMM1array_string(t,"xrdbFiles","xrdbFiles",name,
				      &size,&file,&erunNo);
	if ((ret!=0)||(file==NULL))
	  {
	    fprintf(stderr,"param.loadDatabase: Error reading the Map Manager file\n");
	    /*esend(ALARM_ERROR,"param_load: Cannot load Map Manager database","(%s) for run %d",name,runNo);*/
	    return -1;
	  }
	
	file[size] = 0;
	
	{
	  XrmDatabase new = XrmGetStringDatabase(file);
	  XrmMergeDatabases(new,&paramdb);
	}

	free(file);
	if (paramVerbose)
	  fprintf(stderr,"from the Map Manager, key %d\n",erunNo);
      }

      ret = 0;
      break;
    }

  ret = param_loadIncludes(runNo);

  sprintf(fname,"%d",runNo); param_putValue("runNo",fname);
  return ret;
}

int param_saveDatabase(int run)
{
  char fname[256];
  FILE* fp;
  const char* dir = getenv("RADPHI_CONFIG");

  run = iabs(run);
  
  sprintf(fname,"%s/used/Run%d.db",dir,run);

  if (paramdb!=NULL)
    {
      XrmPutFileDatabase(paramdb,fname);
    }
  return 0;
}

const char* param_getValueI(const char*name,int i)
{
  return param_getValue(param_FI(name,i));
}

int param_valueExists(const char*name)
{
  static XrmValue value;
  Bool xret;
  char*type = NULL;
  char b[256];
  strcpy(b,name);

  while (1)
    {
      xret = XrmGetResource(paramdb,b,NULL,&type,&value);
      if (xret!=True)
	{
	  /*fprintf(stderr,"param_get: Cannot get value for %s\n",b);*/
	  /*esend(ALARM_WARNING,"param_get: Cannot get value for","%s",b);*/
	  return 0;
	}
      if (*value.addr!='@') break;
      sprintf(b,"%s",(char*)value.addr+1);
    }

  return 1;
}

const char* param_getValue(const char*name)
{
  static XrmValue value;
  Bool xret;
  char*type = NULL;
  char b[256];
  strcpy(b,name);

  while (1)
    {
      xret = XrmGetResource(paramdb,b,NULL,&type,&value);
      if (xret!=True)
	{
	  fprintf(stderr,"param_get: Cannot get value for %s\n",b);
	  /*esend(ALARM_WARNING,"param_get: Cannot get value for","%s",b);*/
	  return NULL;
	}
      if (*value.addr!='@') break;
      sprintf(b,"%s",(char*)value.addr+1);
    }

  return value.addr;
}

int param_putValue(const char*name,const char*value)
{
  XrmPutStringResource(&paramdb,name,value);
  return 0;
}

const char* param_getStringI(const char*name,int i,const char* deflt)
{
  return param_getString(param_FI(name,i),deflt);
}

const char* param_getString(const char*name,const char* deflt)
{
  const char* val = param_getValue(name);
  if (val==NULL) return deflt;
  else return val;
}

unsigned long param_getUlongI(const char*name,int i,unsigned long deflt)
{
  return param_getUlong(param_FI(name,i),deflt);
}

unsigned long param_getUlong(const char*name,unsigned long deflt)
{
  const char* val = param_getValue(name);
  if (val==NULL) return deflt;
  else return strtoul(val,NULL,0);
}

uint32 param_getI_uint32(const char*name,int i,uint32 deflt)
{
  return param_get_uint32(param_FI(name,i),deflt);
}

uint32 param_get_uint32(const char*name,uint32 deflt)
{
  const char* val = param_getValue(name);
  if (val==NULL) return deflt;
  else return (uint32)strtoul(val,NULL,0);
}

int32 param_getI_int32(const char*name,int i,int32 deflt)
{
  return param_get_int32(param_FI(name,i),deflt);
}

int32 param_get_int32(const char*name,int32 deflt)
{
  const char* val = param_getValue(name);
  if (val==NULL) return deflt;
  else return (int32)strtol(val,NULL,0);
}


int param_getIntI(const char*name,int i,int deflt)
{
  return param_getInt(param_FI(name,i),deflt);
}

int param_getInt(const char*name,int deflt)
{
  const char* val = param_getValue(name);
  if (val==NULL) return deflt;
  return (int)strtol(val,NULL,0);
}

double param_getDoubleI(const char*name,int i,double deflt)
{
  return param_getDouble(param_FI(name,i),deflt);
}

double param_getDouble(const char*name,double deflt)
{
  const char* val = param_getValue(name);
  if (val==NULL) return deflt;
  return strtod(val,NULL);
}

const vector3_t* param_get3vectorI(const char*name,int i,const vector3_t* deflt)
{
  return param_get3vector(param_FI(name,i),deflt);
}

const vector3_t* param_get3vector(const char*name,const vector3_t* deflt)
{
  static vector3_t v;
  const char* val = param_getValue(name);
  if (val==NULL) return deflt;
  v.x = strtod(val,(void*)&val);
  v.y = strtod(val,(void*)&val);
  v.z = strtod(val,(void*)&val);
  return &v;
}

const vector4_t* param_get4vectorI(const char*name,int i,const vector4_t* deflt)
{
  return param_get4vector(param_FI(name,i),deflt);
}

const vector4_t* param_get4vector(const char*name,const vector4_t* deflt)
{
  static vector4_t v;
  const char* val = param_getValue(name);
  if (val==NULL) return deflt;
  v.space.x = strtod(val,(void*)&val);
  v.space.y = strtod(val,(void*)&val);
  v.space.z = strtod(val,(void*)&val);
  v.t = strtod(val,(void*)&val);
  return &v;
}

int param_getIntArray(const char*name,int ival[],int nval)
{
  int i,v;
  const char* val = param_getValue(name);
  char* ptr;
  if (val==NULL) return 0;

  i = 0;
  while (i<nval)
    {
      v = (int)strtol(val,&ptr,0); 
      if (ptr == val) break;
      ival[i++] = v;
      val = ptr;
    }
  return i;
}

int param_getDblArray(const char*name,double dval[],int nval)
{
  int i;
  double v;
  const char* val = param_getValue(name);
  char* ptr;
  if (val==NULL) return 0;

  i = 0;
  while (i<nval)
    {
      v = strtod(val,&ptr); 
      if (ptr == val) break;
      dval[i++] = v;
      val = ptr;
    }
  return i;
}

const char* param_FS(const char*fmt,const char*s)
{
  static char buf[256];
  sprintf(buf,fmt,s);
  return buf;
}

const char* param_FSII(const char*fmt,const char*s,int v1,int v2)
{
  static char buf[256];
  sprintf(buf,fmt,s,v1,v2);
  return buf;
}

const char* param_FII(const char*fmt,int v1,int v2)
{
  static char buf[256];
  sprintf(buf,fmt,v1,v2);
  return buf;
}

const char* param_FI(const char*fmt,int v)
{
  static char buf[256];
  sprintf(buf,fmt,v);
  return buf;
}

/* end file */

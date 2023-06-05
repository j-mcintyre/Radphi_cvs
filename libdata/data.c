/*
 * data.c
 *
*/


#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <time.h>
#include <limits.h>

#include <disData.h>
#include <crc.h>
#ifdef USE_ALARM
#include <alarm.h>
#endif

int data_check_crc = 1;

int data_newItape(itape_header_t *ptr)
{
  itape_header_t *ih = ptr;
  ih->length = sizeof(itape_header_t);
  ih->type = TYPE_ITAPE;
  ih->transComputerCode = TRANSCOMPUTERCODE;
  ih->ngroups = 0;
  ih->runNo = 0;
  ih->spillNo = 0;
  ih->eventNo = 0;
  ih->trigger   = 0;
  ih->eventType = 0;
  ih->time = time(NULL);
  ih->latch = 0;
  return 0;
}

group_header_t* data_getGroupHeader(const itape_header_t *buffer,uint32 group,uint32 subgroup)
{
  const itape_header_t *ih = buffer;
  const group_header_t *gh;
  int sz;
  int i,n;

  if (ih==NULL) return NULL;

  sz = ih->length;
  n  = ih->ngroups;
  gh = (group_header_t*)(((char*)ih) + sizeof(itape_header_t));
  sz -= sizeof(itape_header_t);
  for (i=0; (sz>0) && (i<n); i++)
    {
      size_t length = gh->length;

      if (length == 0)
	{
	  fprintf(stderr,"data_getGroup: Itape broken: group id(%d),size(%d) has length zero\n",gh->type,gh->length);
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"data_getGroup: Itape broken","group id(%d),size(%d) has length zero\n",gh->type,gh->length);
#endif
	  return NULL;
	}

      if (length > sz)
	{
	  fprintf(stderr,"data_getGroup: Itape broken: group id(%d),size(%d) bigger than itape (size=%d)\n",gh->type,gh->length,sz);
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"data_getGroup: Itape broken","group id(%d),size(%d) bigger than itape (size=%d)\n",gh->type,gh->length,sz);
#endif
	  return NULL;
	}

      if (length&3) /* check that we stay aligned at 4 bytes */
	{
	  fprintf(stderr,"data_getGroup: Itape broken: group id(%d),size(%d) length is not multiple of 4 bytes.\n",gh->type,gh->length);
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"data_getGroup: Itape broken","group id(%d),size(%d) length is not multiple of 4 bytes.\n",gh->type,gh->length);
#endif
	  return NULL;
	}

      if ((gh->type == group) && (gh->subtype == subgroup)) return (group_header_t*)gh;

      gh = (group_header_t*)(((char*)gh) + length);
      sz -= length;
    }

  if ((i != n) || (sz != 0))
    {
      fprintf(stderr,"data_getGroup: Itape broken: group id(%d),size(%d) bigger than itape (size=%d)\n",gh->type,gh->length,sz);
#ifdef USE_ALARM
      esend(ALARM_ERROR,"data_getGroup: Itape broken","group id(%d),size(%d) bigger than itape (size=%d)\n",gh->type,gh->length,sz);
#endif
    }

  return NULL;
}

void* data_getGroup(const itape_header_t *buffer, uint32 group, uint32 subgroup)
{
  char*p = (char*) data_getGroupHeader(buffer,group,subgroup);
  if (p==NULL) return NULL;
  p += sizeof(group_header_t);
  return p;
}

size_t data_getGroupSize(const itape_header_t *buffer,uint32 group, uint32 subgroup)
{
  group_header_t *p = data_getGroupHeader(buffer,group,subgroup);
  if (p==NULL) return ULONG_MAX;
  return p->length - sizeof(group_header_t);
}

group_header_t *data_addGroup_header(itape_header_t *buffer,size_t bufsize,uint32 type, uint32 subtype,size_t length)
{
  itape_header_t *ih = buffer;
  group_header_t *gh;
  int i;

  if (ih==NULL) return NULL;

  data_removeGroup(ih,type,subtype);

  if (ih->length + sizeof(group_header_t) + length>=bufsize)
    {
      fprintf(stderr,"data_addGroup: Buffer overflow: buffer length is %d, used %d, create group type %d length %d\n",bufsize,ih->length,type,length);
#ifdef USE_ALARM
      esend(ALARM_ERROR,"data_addGroup: Buffer overflow","buffer length is %d, used %d, create group type %d length %d\n",bufsize,ih->length,type,length);
#endif
      return NULL;
    }

  gh = (void*)((char*)ih + sizeof(*ih));
  for (i=0; i<ih->ngroups; i++)
    {
      gh = (void*)((char*)gh + gh->length);
    }
/*
  gh = (group_header_t*) (p + ih->length);
*/
  gh->length = (uint32)length;
  gh->type = type;
  gh->subtype = subtype;

  ih->length += (int)sizeof(group_header_t) + length;
  ih->ngroups ++;

  return gh;
}

void *data_addGroup(itape_header_t *buffer,size_t bufsize,uint32 type,uint32 subtype,size_t length)
{
  itape_header_t *ih = buffer;
  char *p = (void *)buffer;
  group_header_t *gh;

  if (ih==NULL) return NULL;

  data_removeGroup(buffer,type,subtype);

  if ((length%4)!=0)
    {
      /* align everything to a 4 byte boundary */

      length = ((length/4) + 1)*4;
    }

  if (ih->length + sizeof(group_header_t) + length>=bufsize)
    {
      fprintf(stderr,"data_addGroup: Buffer overflow: buffer length is %d, used %d, create group type %d length %d\n",bufsize,ih->length,type,length);
#ifdef USE_ALARM
      esend(ALARM_ERROR,"data_addGroup: Buffer overflow","buffer length is %d, used %d, create group type %d length %d\n",bufsize,ih->length,type,length);
#endif
      return NULL;
    }

  gh = (group_header_t*) (p + ih->length);
  gh->length = (uint32)sizeof(group_header_t) + (uint32)length;
  gh->type = type;
  gh->subtype = subtype;

  ih->length += (int)sizeof(group_header_t) + length;
  ih->ngroups ++;

  return ((char*)gh)+sizeof(group_header_t);
}

void *data_getAdr(const itape_header_t *buffer)
{
  const itape_header_t *ih = buffer;
  return ((char*)buffer + ih->length + sizeof(group_header_t));
}

int data_addCRC(itape_header_t *buffer)
{
  unsigned char *p = (void *)buffer;
  itape_header_t *ih = buffer;
  uint32 crc = data_crc(ih,ih->length);
  p[ih->length] = (unsigned char)(crc & 0xff);
  p[ih->length+1] = (unsigned char)((crc>>8) & 0xff);
  p[ih->length+2] = (unsigned char)((crc>>16) & 0xff);
  p[ih->length+3] = (unsigned char)((crc>>24) & 0xff);
  return 0;
}     

int data_checkCRC1(const itape_header_t *buffer,size_t bufferLength)
{
  if (data_check_crc)
    {
      const unsigned char *p = (void *) buffer;
      uint32 crc = data_crc(buffer,bufferLength);
      uint32 oldcrc = p[bufferLength] | p[bufferLength+1]<<8 |
	p[bufferLength+2]<<16 | p[bufferLength+3]<<24;
/*
fprintf(stderr,"crc: new: 0x%x data: 0x%x\n",crc,oldcrc);
*/
      return crc!=oldcrc;
    }
  else return 0;
}     

int data_checkCRC(const itape_header_t *buffer)
{
  if (data_check_crc)
    {
      const itape_header_t *ih = buffer;

      return data_checkCRC1(ih,ih->length);
    }
  else return 0;
}     

int data_clean(itape_header_t * buffer)
{
  char* wptr = (void *)buffer;
  const char* rptr = (void *)buffer;
  itape_header_t* ih = buffer;
  group_header_t* gh = (void *)buffer;
  uint32 sz;
  uint32 count;
  uint32 size;

  count = 0;
  size = ih->length;

  if (size < sizeof(itape_header_t))
    {
      fprintf(stderr,"data_clean: Itape broken (0): Itape is too short (size=%d)\n",size);
#ifdef USE_ALARM
      esend(ALARM_ERROR,"data_clean: Itape broken","Itape is too short (size=%d)\n",size);
#endif
      return -1;
    }

  rptr += sizeof(itape_header_t);
  size -= sizeof(itape_header_t);

  gh = (void*) rptr;
  while (gh->type!=0)
    {
      uint32 length = gh->length;

      if ((length == 0)||(length > size))
	{
	  fprintf(stderr,"data_clean: Itape broken (1): group id(%d),size(%d) bigger than itape (size=%d)\n",
		  gh->type,length,size);
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"data_clean: Itape broken","group id(%d),size(%d) bigger than itape (size=%d)\n",
		gh->type,length,size);
#endif
	  return -1;
	}

      if (length&3)
	{
	  fprintf(stderr,"data_clean: Itape broken (1): group id(%d),size(%d) length is not a multiple of 4 (size=%d)\n",
		  gh->type,length,size);
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"data_clean: Itape broken","group id(%d),size(%d) length is not a multiple of 4 (size=%d)\n",
		gh->type,length,size);
#endif
	  return -1;
	}

      rptr += length;
      size -= length;
      gh = (void*) rptr;
      count ++;

      if (size==0)
	return 0;
    }

  if (size==0)
    return 0;

  wptr = (void*)rptr;

  while (size>0)
    {
      while (gh->type == 0)
	{
	  uint32 length = gh->length;

	  if ((length == 0)||(length > size))
	    {
	      fprintf(stderr,"data_clean: Itape broken (2): group id(%d),size(%d) bigger than itape (size=%d)\n",
		      gh->type,gh->length,size);
#ifdef USE_ALARM
	      esend(ALARM_ERROR,"data_clean: Itape broken","group id(%d),size(%d) bigger than itape (size=%d)\n",
		    gh->type,gh->length,size);
#endif
	      return -1;
	    }

	  if (length&3)
	    {
	      fprintf(stderr,"data_clean: Itape broken (2): group id(%d),size(%d) length is not a multiple of 4 (size=%d)\n",
		      gh->type,length,size);
#ifdef USE_ALARM
	      esend(ALARM_ERROR,"data_clean: Itape broken","group id(%d),size(%d) length is not a multiple of 4 (size=%d)\n",
		    gh->type,length,size);
#endif
	      return -1;
	    }

	  rptr += length;
	  size -= length;
	  gh = (void*) rptr;
	  if (size==0)
	    break;
	}

      if (size==0) break;
      sz = gh->length;

      if ((sz == 0)||(sz > size))
	{
	  fprintf(stderr,"data_clean: Itape broken (3): group id(%d),size(%d) bigger than itape (size=%d)\n",
		  gh->type,gh->length,size);
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"data_clean: Itape broken","group id(%d),size(%d) bigger than itape (size=%d)\n",
		gh->type,gh->length,size);
#endif
	  return -1;
	}

      if (sz&3)
	{
	  fprintf(stderr,"data_clean: Itape broken (3): group id(%d),size(%d) length is not a multiple of 4 (size=%d)\n",
		  gh->type,gh->length,size);
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"data_clean: Itape broken","group id(%d),size(%d) length is not a multiple of 4 (size=%d)\n",
		gh->type,gh->length,size);
#endif
	  return -1;
	}

#ifdef UNDEF /* this was a crashing memmove() scare... */

      if ((rptr >= wptr)&&(rptr <= wptr + sz))
	memmove(wptr,rptr,sz);
      else
	memcpy(wptr,rptr,sz);
#endif

      memmove(wptr,rptr,sz);

      wptr += sz;
      rptr += sz;
      size -= sz;
      count ++;
      gh = (void*) rptr;
    }

  ih->ngroups = count;
  ih->length = (uint32)((char*)wptr - (char*)buffer);

  if (ih->length&3)
    {
      fprintf(stderr,"data_clean: Itape broken (4): Itape length is not a multiple of 4 (size=%d)\n",
	      ih->length);
#ifdef USE_ALARM
      esend(ALARM_ERROR,"data_clean: Itape broken","Itape length is not a multiple of 4 (size=%d)\n",
	    ih->length);
#endif
      return -1;
    }

  return 0;
}

int data_removeGroup(itape_header_t *buffer,uint32 group,uint32 subgroup)
{
  while (1)
    {
      group_header_t* gh = data_getGroupHeader(buffer,group,subgroup);
      if (gh==NULL) return 1;
      gh->type = 0;
    }
}

int data_saveGroups(itape_header_t *buffer,int nsave,int *isave,int *sisave)
{
  char*p = (void *)buffer;
  itape_header_t *h = buffer;
  group_header_t *g;
  int sz;
  int i,n;
  int isv,save;

  if (h==NULL) return 0;

  sz = h->length;
  n = h->ngroups;
  p += sizeof(itape_header_t);
  sz -= sizeof(itape_header_t);
  for (i=0; i<n; i++) {
      g = (group_header_t*) p;
      /* if g is in the list to save - don't remove it */
      for (save = 0,isv = 0; isv < nsave; ++isv) {
	if ((g->type == isave[isv]) && (g->subtype == sisave[isv]))
	  save = 1;
      }
      if (!save) 
	g->type = 0;
      if ((g->length>sz)||(g->length==0)) {
	  fprintf(stderr,"data_saveGroups: Itape broken: group id(%d),size(%d) bigger than itape (size=%d)\n",g->type,g->length,sz);
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"data_saveGroups: Itape broken","group id(%d),size(%d) bigger than itape (size=%d)\n",g->type,g->length,sz);
#endif
	  return -1;
	}
    
      p += g->length;
      sz -= g->length;
    }

  return 0;
}

int data_removeAllgroups(itape_header_t *buffer)
{
  char*p = (void *)buffer;
  itape_header_t *h = buffer;
  group_header_t *g;
  int sz;
  int i,n;

  if (h==NULL) return 0;

  sz = h->length;
  n = h->ngroups;
  p += sizeof(itape_header_t);
  sz -= sizeof(itape_header_t);
  for (i=0; i<n; i++) {
      g = (group_header_t*) p;
      g->type = 0;
      if ((g->length>sz)||(g->length==0)) {
	  fprintf(stderr,"data_removeAllgroups: Itape broken: group id(%d),size(%d) bigger than itape (size=%d)\n",g->type,g->length,sz);
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"data_removeAllgroups: Itape broken","group id(%d),size(%d) bigger than itape (size=%d)\n",g->type,g->length,sz);
#endif
	  return -1;
	}
      p += g->length;
      sz -= g->length;
    }
  return 0;
}

int data_listGroups(const itape_header_t *event,int *ngroups_ret,uint32 **list_ret, uint32 **slist_ret, size_t **size_ret)
{
  const char*p = (void *)event;
  const itape_header_t *h = event;
  const group_header_t *g;
  int sz;
  int i,n;
  int ngroups = 0;

  if (ngroups_ret) (*ngroups_ret) = 0;

  if (h==NULL)
	return -1;

  if (h->type != TYPE_ITAPE)
	return -1;

  sz = h->length;
  n = h->ngroups;

  if (ngroups_ret) (*ngroups_ret) = 0;
  if (list_ret) (*list_ret) = (uint32 *)malloc(sizeof(list_ret[0]) * n);
  if (slist_ret) (*slist_ret) = (uint32 *)malloc(sizeof(slist_ret[0]) * n);
  if (size_ret) (*size_ret) = (size_t *)malloc(sizeof(size_ret[0]) * n);

  p += sizeof(itape_header_t);
  sz -= sizeof(itape_header_t);
  for (i=0; i<n; i++)
    {
      g = (group_header_t*) p;
      if ((g->length > sz)||(g->length <= 0)||(g->length&3))
	{
	  fprintf(stderr,"data_listGroups: Itape broken: group id(%d),size(%d) bigger than itape (size=%d)\n",g->type,g->length,sz);
#ifdef USE_ALARM
	  esend(ALARM_ERROR,"data_listGroups: Itape broken","group id(%d),size(%d) bigger than itape (size=%d)\n",g->type,g->length,sz);
#endif
	  if (ngroups_ret) (*ngroups_ret) = ngroups;
	  return -1;
	}

      if (list_ret) (*list_ret)[ngroups] = g->type;
      if (slist_ret) (*slist_ret)[ngroups] = g->subtype;
      if (size_ret) (*size_ret)[ngroups] = g->length;
      ngroups++;

      p += g->length;
      sz -= g->length;
    }

  if (ngroups_ret) (*ngroups_ret) = ngroups;

  return 0;
}


/* end file */

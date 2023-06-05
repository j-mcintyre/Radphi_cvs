/*
 * endianConverter.c
*/

/*
  $Log: endianConverter.c,v $
  Revision 1.3  2000/01/18 04:03:03  radphi
  Changes made by jonesrt@golem
  fixed endian-conversion package to allow different structures with the
  same group number but different subgroups -- to support MonteCarlo group. -rtj

 * Revision 1.2  1996/11/15  15:08:58  lfcrob
 * Fixed funxtion dec's to use itape_header_t * instead of void *
 * Added RCS logging stuff
 *
*/

static char rcsid[] = "$Id: endianConverter.c,v 1.3 2000/01/18 04:03:03 radphi Exp $";

#include <stdio.h>
#include <stdlib.h>

#include <disData.h>
#include <ntypes.h>

#include <endianConverterAuto.h>

int endian_convertType1_NOCONVERT_t(void *w,int len)
{
  return 0;
}

int endian_convertType1_uint32(uint32 *w,int len)
{
  *w = (((*w)&0xff)<<24)|(((*w)&0xff00)<<8)|(((*w)&0xff0000)>>8)|(((*w)&0xff000000)>>24);
  
  return 0;
}

int endian_convertTypeN_uint32(uint32 ww[],int count)
{
  int i;

  for (i=0; i<count; i++)
    ww[i] = ((ww[i]&0xff)<<24)|((ww[i]&0xff00)<<8)|((ww[i]&0xff0000)>>8)|((ww[i]&0xff000000)>>24);

  return 0;
}

int endian_convertType1_uint16(uint16 *w,int len)
{
  *w = (((*w)&0xff)<<8)|(((*w)&0xff00)>>8);
  
  return 0;
}

int endian_convertTypeN_uint16(uint16 ww[],int count)
{
  int i;

  for (i=0; i<count; i++)
    ww[i] = ((ww[i]&0xff)<<8)|((ww[i]&0xff00)>>8);

  return 0;
}

int endian_convertTypeN_vector3_t(vector3_t ww[],int count)
{
  return endian_convertTypeN_uint32((void*)ww,count*3);
}

int endian_convertTypeN_vector4_t(vector4_t ww[],int count)
{
  return endian_convertTypeN_uint32((void*)ww,count*4);
}

int endian_convertType1_vector3_t(vector3_t *ww,int len)
{
  return endian_convertTypeN_uint32((void*)ww,3);
}

int endian_convertType1_vector4_t(vector4_t *ww,int len)
{
  return endian_convertTypeN_uint32((void*)ww,4);
}

int endian_convertType1_double(double *ww,int len)
{
  uint32 *w = (uint32*)ww;
  uint32 tmp;

  endian_convertTypeN_uint32(w,2);

  tmp = w[0]; w[0] = w[1]; w[1] = tmp;

  return 0;
}

int endian_convertTypeN_double(double ww[],int count)
{
  int i;

  for (i=0; i<count; i++)
    endian_convertType1_double(&ww[i],sizeof(ww[i]));

  return 0;
}

int endian_convertTypeN_uint8(uint8 ww[],int count)
{
  return 0;
}

int endian_convertType1_uint8(uint8 *ww,int len)
{
  return 0;
}

int endian_convertType1_INTARRAY_t(void *ww,int count)
{
  return endian_convertTypeN_uint32(ww,count);
}

int endian_convertItape(itape_header_t *itape,int itapeLength)
{
  int ret = 0;
  itape_header_t *ih = itape;
  uint32 type;

  type = ih->type;

  endian_convertType1_uint32(&type,sizeof(uint32));

  switch (type)
    {
    case TYPE_TAPEHEADER:
      endian_convertType1_tapeHeader_t((void *)itape,itapeLength);
      break;

    case TYPE_ITAPE:
      {
	int err = 0;
	int ngroups,igroup;
	char *ptr;
	group_header_t *hdr;
	ptr = (void *)itape;

	endian_convertType1_itape_header_t(itape,itapeLength);

	ngroups = ih->ngroups;

	ptr = (ptr + sizeof(itape_header_t));

	for (igroup=0; igroup<ngroups; igroup++)
	  {
	    uint32 length;
	    group_header_t *hdr = (void*)ptr;
	    
	    /* fprintf(stderr,"itape: 0x%x, hdr: 0x%x, length: 0x%x, group: 0x%x\n",itape,hdr,hdr->length,hdr->type); */

	    endian_convertType1_group_header_t(hdr,sizeof(group_header_t));
	    
	    length = hdr->length;

	    if (length >= itapeLength)
	      {
		fprintf(stderr,"endian_convertItape: Error: Broken data group: len: 0x%x, group: 0x%x.\n",
			hdr->length,hdr->type);
		err |= 1;

		break;  /* this event is broken- do not try to convert it. give up now. */
	      }

	    err |= endian_convertGroup(hdr,ptr + sizeof(group_header_t),hdr->length - sizeof(group_header_t));

	    ptr += hdr->length;
	  }

	if (err)
	  {
	    fprintf(stderr,"endian_convertItape: Error: The event was only partially converted because of previous errors.\n");
	    ret = 1;
	  }

	data_addCRC(itape);
      }
      break;
    }

  return ret;
}

/* end file */

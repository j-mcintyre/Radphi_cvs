/*
  $Log: nPartMass.c,v $
  Revision 1.2  1997/05/29 20:26:44  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: nPartMass.c,v 1.2 1997/05/29 20:26:44 radphi Exp $";

#include <stdio.h>
#include <stdarg.h>

#include <itypes.h>
#include <clib.h>

float nPartMass(int nPart,...)
{
  va_list args;
  vector4_t *pArray[20];
  float mass;
  int ptr_no=0;

  va_start(args,nPart);
  
  while(ptr_no < nPart)
    (pArray[ptr_no++] = va_arg(args,void *));
  mass = EffMass(nPart,(const vector4_t **)pArray);
  va_end(args);
  return(mass);
}

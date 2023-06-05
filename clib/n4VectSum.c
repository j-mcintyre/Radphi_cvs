/*
  $Log: n4VectSum.c,v $
  Revision 1.2  1997/05/29 20:26:43  radphi
  Changes made by lfcrob@dustbunny
  ANSI-fied and Radphi-fied :)

*/

static const char rcsid[] = "$Id: n4VectSum.c,v 1.2 1997/05/29 20:26:43 radphi Exp $";

#include <stdio.h>
#include <stdarg.h>
#include <itypes.h>
#include <clib.h>



vector4_t n4VectSum(int nVect,...)
{
  va_list args;
  vector4_t *pArray[20];
  vector4_t answer={0,0,0,0};
  int ptr_no=0;
  int index;

  va_start(args,nVect);
  while(ptr_no<nVect)
    (pArray[ptr_no++] = va_arg(args,void *));
  for(index=0;index<nVect;index++){
    answer.t += pArray[index]->t;
    answer.space.x += pArray[index]->space.x;
    answer.space.y += pArray[index]->space.y;
    answer.space.z += pArray[index]->space.z;
  }
  return(answer);
}

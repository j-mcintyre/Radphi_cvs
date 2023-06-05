#include "Box.h"


Box::Box(double a, double b, double c, double d)   // Box constructor.
{                                              
  left = a;
  right = b;
  top = c;
  bottom = d;
}


int Box::RectCut(point *pos)          // returns 1,2,3 or 4 
                                      // if location is within
{                                     // range defined by Box.
				       
  int t;
  if(pos->y > top)
    return(1);				     
  if(pos->y < bottom)
    return(2);				     
  if(pos->x < left)
    return(3);				     
  if(pos->x > right)
    return(4);
  return(0);				     
//  t = ((pos->x > left) && (pos->x < right));
//  return (t && (pos->y > bottom) && (pos->y < top));
}


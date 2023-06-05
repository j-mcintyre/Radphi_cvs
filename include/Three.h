// CHEES header file

/*
    Three.h

    header for the Three class
     (a three vector class)

*/

#include <ntypes.h>
#include <stdio.h>

#ifndef _THREE_CLASS_INCLUDED_
#define _THREE_CLASS_INCLUDED_ 

class Four;

class Three
{
 public:

  dvector3_t vec;

  Three(double = 0.0, double = 0.0, double = 0.0);
  Three(dvector3_t);
  Three(vector3_t);
  Three(vector4_t);

  Three operator+(Three &);
  Three operator-(Three &);
  void operator+=(Three &);
  friend Three operator*(Three &, float);
  friend Three operator*(float,Three &);
  friend Three operator/(Three &, float);
  double operator*(Three &);
  Three operator/(Three &);

  void Print(void);
  double Mag(void);  /* magnitude */
  
};

/* misc (non-member) functions */

  double Det(Three &,Three &,Three &);  /* determinant */

#endif

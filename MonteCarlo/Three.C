/* 
    Three.C

    basic three vector class +
    3x3 matrix determinant (uses three Threes for the matrix)

*/

#include "Three.h"
#include <math.h>

/* --------------------------------------------------------------------- */


Three::Three(double a, double b, double c)
{
  vec.x = a;
  vec.y = b;
  vec.z = c;
}

Three::Three(dvector3_t one)
{
  vec.x = one.x;
  vec.y = one.y;
  vec.z = one.z;
}

Three::Three(vector3_t one)
{
  vec.x = one.x;
  vec.y = one.y;
  vec.z = one.z;
}

Three::Three(vector4_t one)
{
  vec.x = one.space.x;
  vec.y = one.space.y;
  vec.z = one.space.z;
}

Three Three::operator+(Three &one)
{
  Three temp(vec.x + one.vec.x, vec.y + one.vec.y, vec.z + one.vec.z);
  return temp;
}

Three Three::operator-(Three &one)
{
  Three temp(vec.x - one.vec.x, vec.y - one.vec.y, vec.z - one.vec.z);
  return temp;
}

void Three::operator+=(Three &one)
{
  vec.x += one.vec.x;
  vec.y += one.vec.y;
  vec.z += one.vec.z;
}


Three operator*(Three &one, float t)
{
  Three temp(one.vec.x * t, one.vec.y * t, one.vec.z * t);
  return temp;
}

Three operator*(float t,Three &one)
{
  Three temp(one.vec.x * t, one.vec.y * t, one.vec.z * t);
  return temp;
}

Three operator/(Three &one, float t)
{
  Three temp(one.vec.x / t, one.vec.y / t, one.vec.z / t);
  return temp;
}

double Three::operator*(Three &one)
{
  return (one.vec.x * vec.x + one.vec.y * vec.y + one.vec.z * vec.z);
}


Three Three::operator/(Three &one)  /* cross product */
{
  Three temp( (vec.y * one.vec.z - vec.z * one.vec.y), 
              (vec.z * one.vec.x - vec.x * one.vec.z), 
              (vec.x * one.vec.y - vec.y * one.vec.x)      );

  return temp;

}

double Det(Three &row1, Three &row2, Three &row3 )
{
  /* 
     find determinant of a 3x3 matrix defined by 

     matrix =  |-          -|
               |    row1    |
               |    row2    |
               |    row3    |
               |-          -|

  */

  return   row1.vec.x * (row2.vec.y * row3.vec.z - row3.vec.y * row2.vec.z )
         - row1.vec.y * (row2.vec.x * row3.vec.z - row3.vec.x * row2.vec.z )
         + row1.vec.z * (row2.vec.x * row3.vec.y - row3.vec.x * row2.vec.y ) ; 

}

void Three::Print(void)
{
  printf ("\n x: %f y:%f z:%f\n",vec.x,vec.y,vec.z);
}

double Three::Mag(void)
{
 return sqrt( vec.x  * vec.x + vec.y * vec.y + vec.z * vec.z );
}



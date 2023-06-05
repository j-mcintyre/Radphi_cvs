/*
    SqMatrix3 object.

    A simple 3x3 matrix class to help with Cramer's Rule calculations.

    PCF 23 Aug 1995

*/

#include <SqMatrix3.h>
#include "Three.h"

//////////////////////////////////////////////////////////////

SqMatrix3::SqMatrix3()
{

}

//////////////////////////////////////////////////////////////

SqMatrix3::SqMatrix3(Three row1,Three row2,Three row3)
{
  n[1][1] = row1.vec.x;
  n[1][2] = row1.vec.y;
  n[1][3] = row1.vec.z;

  n[2][1] = row2.vec.x;
  n[2][2] = row2.vec.y;
  n[2][3] = row2.vec.z;

  n[3][1] = row3.vec.x;
  n[3][2] = row3.vec.y;
  n[3][3] = row3.vec.z;
}

//////////////////////////////////////////////////////////////

void SqMatrix3::Set(Three row1,Three row2,Three row3)
{
  n[1][1] = row1.vec.x;
  n[1][2] = row1.vec.y;
  n[1][3] = row1.vec.z;

  n[2][1] = row2.vec.x;
  n[2][2] = row2.vec.y;
  n[2][3] = row2.vec.z;

  n[3][1] = row3.vec.x;
  n[3][2] = row3.vec.y;
  n[3][3] = row3.vec.z;
}

//////////////////////////////////////////////////////////////

SqMatrix3 SqMatrix3::operator*(SqMatrix3 &one)  /* matrix times matrix */
{
SqMatrix3 temp;

   temp.n[1][1] = n[1][1]*one.n[1][1]  +  n[1][2]*one.n[2][1]  +  n[1][3]*one.n[3][1];
   temp.n[1][2] = n[1][1]*one.n[1][2]  +  n[1][2]*one.n[2][2]  +  n[1][3]*one.n[3][2];
   temp.n[1][3] = n[1][1]*one.n[1][3]  +  n[1][2]*one.n[2][3]  +  n[1][3]*one.n[3][3];
   
   temp.n[2][1] = n[2][1]*one.n[1][1]  +  n[2][2]*one.n[2][1]  +  n[2][3]*one.n[3][1];
   temp.n[2][2] = n[2][1]*one.n[1][2]  +  n[2][2]*one.n[2][2]  +  n[2][3]*one.n[3][2];
   temp.n[2][3] = n[2][1]*one.n[1][3]  +  n[2][2]*one.n[2][3]  +  n[2][3]*one.n[3][3];
 
   temp.n[3][1] = n[3][1]*one.n[1][1]  +  n[3][2]*one.n[2][1]  +  n[3][3]*one.n[3][1];
   temp.n[3][2] = n[3][1]*one.n[1][2]  +  n[3][2]*one.n[2][2]  +  n[3][3]*one.n[3][2];
   temp.n[3][3] = n[3][1]*one.n[1][3]  +  n[3][2]*one.n[2][3]  +  n[3][3]*one.n[3][3];


	return temp;
}

//////////////////////////////////////////////////////////////

Four SqMatrix3::operator*(Four &one) /* 3x3 matrix times spatial 
                                          component of vector */
{
Four temp;
 
 temp.vec.space.x = n[1][1]*one.vec.space.x + n[1][2]*one.vec.space.y + 
                    n[1][3]*one.vec.space.z;
                    
 temp.vec.space.y = n[2][1]*one.vec.space.x + n[2][2]*one.vec.space.y + 
                    n[2][3]*one.vec.space.z;

 temp.vec.space.z = n[3][1]*one.vec.space.x + n[3][2]*one.vec.space.y + 
                    n[3][3]*one.vec.space.z;
                    
                    
 temp.vec.t = one.vec.t;

   return temp;
   
}

Three SqMatrix3::operator*(Three &one) /* 3x3 matrix times spatial 
                                          component of vector */
{
Three temp;
 
 temp.vec.x = n[1][1]*one.vec.x + n[1][2]*one.vec.y + 
                    n[1][3]*one.vec.z;
                    
 temp.vec.y = n[2][1]*one.vec.x + n[2][2]*one.vec.y + 
                    n[2][3]*one.vec.z;

 temp.vec.z = n[3][1]*one.vec.x + n[3][2]*one.vec.y + 
                    n[3][3]*one.vec.z;
                    
                    

   return temp;
   
}
//////////////////////////////////////////////////////////////

double SqMatrix3::Determinant(void)
{

  return - n[1][3]*n[2][2]*n[3][1] + n[1][2]*n[2][3]*n[3][1]
         + n[1][3]*n[2][1]*n[3][2] - n[1][1]*n[2][3]*n[3][2]
         - n[1][2]*n[2][1]*n[3][3] + n[1][1]*n[2][2]*n[3][3] ;
}

//////////////////////////////////////////////////////////////

SqMatrix3 SqMatrix3::ReplaceColumn(int colNo, Three newCol)
{
SqMatrix3 answer;

  answer.n[1][1] = n[1][1];
  answer.n[2][1] = n[2][1];
  answer.n[3][1] = n[3][1]; 

  answer.n[1][2] = n[1][2]; 
  answer.n[2][2] = n[2][2]; 
  answer.n[3][2] = n[3][2]; 

  answer.n[1][3] = n[1][3]; 
  answer.n[2][3] = n[2][3]; 
  answer.n[3][3] = n[3][3]; 

  answer.n[1][colNo] = newCol.vec.x;
  answer.n[2][colNo] = newCol.vec.y;
  answer.n[3][colNo] = newCol.vec.z;

  return answer;

}

//////////////////////////////////////////////////////////////

void SqMatrix3::print(void)
{
 printf("\n");
 printf("%f %f %f\n",n[1][1],n[1][2],n[1][3]);
 printf("%f %f %f\n",n[2][1],n[2][2],n[2][3]);
 printf("%f %f %f\n",n[3][1],n[3][2],n[3][3]);
 printf("\n");

}

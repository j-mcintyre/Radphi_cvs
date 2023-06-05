/////////////////////////////////////////////////////////
// Four class code                                     //
//                                                     //
// PCF 19 June 1995                                    //
/////////////////////////////////////////////////////////

#include "Four.h"


Four::Four(double x, double y, double z, double t)
{
  vec.t = t;                 // Four class constructor
  vec.space.x = x;
  vec.space.y = y;
  vec.space.z = z;
}


Four::Four(const Four &one)
{                                    // Copy constructor
  vec.space.x = one.vec.space.x;
  vec.space.y = one.vec.space.y;
  vec.space.z = one.vec.space.z;
  vec.t = one.vec.t;
}


Four::Four(double one[4])
{                              // Constructor based on array
  vec.space.x = one[0];
  vec.space.y = one[1];
  vec.space.z = one[2];
  vec.t = one[3];
}


void Four::print(void)        // Output vector to screen
{
  cout << vec.space.x << ',' << vec.space.y << ',' << vec.space.z;
  cout << ',' << vec.t << endl;
}

  
Four::Four(dvector3_t one, double t)
{                                    // Constructor based on 3 vector
  vec.space.x = one.x;
  vec.space.y = one.y;
  vec.space.z = one.z;
  vec.t = t;
}

Four::Four(dvector4_t one)           // Constructor based on 4 vector
{
  vec.space.x = one.space.x;
  vec.space.y = one.space.y;
  vec.space.z = one.space.z;
  vec.t = one.t;
}

Four::Four(vector4_t one)           // Constructor based on 4 vector
{
  vec.space.x = one.space.x;
  vec.space.y = one.space.y;
  vec.space.z = one.space.z;
  vec.t = one.t;
}

void Four::Set(double x, double y, double z, double t)
{                                    // Set values in a 4 vector
  vec.t = t;
  vec.space.x = x;
  vec.space.y = y;
  vec.space.z = z;
}

Four operator*(Four &one, double c)
{                                   // Scalar multiplication
  Four temp;

  temp.vec.space.x = one.vec.space.x * c;
  temp.vec.space.y = one.vec.space.y * c;
  temp.vec.space.z = one.vec.space.z * c;
  temp.vec.t = one.vec.t * c;

  return temp;
}

double Four::dot3(Four* one)
{                                 // Does a dot product for the space vectors
 double answer;

 answer = (one->vec.space.x * vec.space.x) +
          (one->vec.space.y * vec.space.y) +  
          (one->vec.space.z * vec.space.z);  

 return answer;
}

double Four::operator*(Four &one)
{

// dot product of four vectors

  double answer;


  answer = one.vec.t * vec.t - (
           
           (one.vec.space.x * vec.space.x) +
           (one.vec.space.y * vec.space.y) +  
           (one.vec.space.z * vec.space.z)

                               );  

  return answer;
}


dvector3_t Four::operator/(Four &one)
{
 /* three vector cross product */

  dvector3_t temp;

  temp.x =  (vec.space.y * one.vec.space.z) - (one.vec.space.y * vec.space.z);
  temp.y =  (one.vec.space.x * vec.space.z) - (vec.space.x * one.vec.space.z);
  temp.z =  (vec.space.x * one.vec.space.y) - (one.vec.space.x * vec.space.y);

  return temp;
}


Four Four::operator+(Four &one)
{                              // Addition of two four vectors
  Four temp;

  temp.vec.space.x = one.vec.space.x + vec.space.x;
  temp.vec.space.y = one.vec.space.y + vec.space.y;
  temp.vec.space.z = one.vec.space.z + vec.space.z;
  temp.vec.t =  one.vec.t + vec.t;

  return temp;
}


void Four::operator=(dvector4_t &one)     // Set one 4 vector equal to another
{
  vec.space.x = one.space.x;
  vec.space.y = one.space.y;
  vec.space.z = one.space.z;
  vec.t = one.t;
}

void Four::operator=(vector4_t &one)     // Set one 4 vector equal to another
{
  vec.space.x = one.space.x;
  vec.space.y = one.space.y;
  vec.space.z = one.space.z;
  vec.t = one.t;
}

int Four::operator>(Four &one)          // Determines if the magnitude of the
{                                   // space vectors of the first is greater
  return (mag3() > one.mag3());     // than those of the second 
}


int Four::operator<(Four &one)      // Magnitude of space vectors less than
{
  return (mag3() < one.mag3());
}


void Four::operator+=(Four &one)      // Adds corresponding vectors of second
{                                   // object to that of first
  vec.space.x += one.vec.space.x;
  vec.space.y += one.vec.space.y;
  vec.space.z += one.vec.space.z;

  vec.t += one.vec.t;
}


void Four::operator*=(double c)     // Sets vector equal to itself after 
{                                // scalar multiplication
  vec.space.x *= c;
  vec.space.y *= c;
  vec.space.z *= c;

  vec.t *= c;
}


double Four::mag3(void)                 // Returns the magnitude of the 
{                                     // space vectors
  double x = vec.space.x*vec.space.x +
            vec.space.y*vec.space.y +
            vec.space.z*vec.space.z;

  return sqrt(x);
}


double Four::mag4(void)               // Returns the magnitude of all the
{                                   // vectors
  double x = mag3();
  x = (vec.t * vec.t) - (x * x);
  return sqrt(x);
}


double Four::mag4Sq(void)             // Returns the square of the magnitude
{                                  // of all the vectors
  double x = mag3();
  x = (vec.t * vec.t) - (x * x);
  return x;
}



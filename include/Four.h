// CHEES header file
// Four class header file

#include <iostream>
#include <math.h>
#include <ntypes.h>
using namespace std;

#ifndef _FOUR_CLASS_INCLUDED_
#define _FOUR_CLASS_INCLUDED_

class Particle;
class Keeper;

class Four{
 public:

  friend ostream &operator<<(ostream &, Particle &);
  friend class Keeper;
  friend class Particle;
// Member Functions ////////////////////////////////////////

  Four(double = 0.0, double = 0.0, double = 0.0, double = 0.0);
  Four(const Four &);
  Four(double[4]);
  Four(dvector3_t, double);
  Four(dvector4_t);
  Four(vector4_t);

 void Set(double,double,double,double); // sets the four vec indicies 

  friend Four operator*(Four &, double);
  double operator*(Four &);
  dvector3_t operator/(Four &);   /* 3-vector cross product */
  Four operator+(Four &);
  void operator=(dvector4_t &);
  void operator=(vector4_t &);
  int operator>(Four &);
  int operator<(Four &);
  void operator+=(Four &);
  void operator*=(double);

  double mag3(void);
  double mag4(void);
  double mag4Sq(void);

  double dot3(Four*);  // dots the three vector parts

  void print(void);

// Data /////////////////////////////////////////////////////

  dvector4_t vec;

 private:

};

#endif







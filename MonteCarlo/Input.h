// CHEES header file
// Input.h
// reads input file and sets up parameters for CHEES 
//
// PCF 21 June 1995
//

// input file struct defs + link to input reading routine
extern "C" {
#include <mc_param.h>
}


#ifndef _INPUT_CLASS_INCLUDED_
#define _INPUT_CLASS_INCLUDED_

class Particle;

class Input
{

public:

// data members

mc_param_t mInParams;

// member functions

      Input();
int   ReadInput(char*);    // returns true if everything went okay 
void  WriteHeader(void);   // writes output file header
void  SetUpParticles(void);   // creates and sets up Particles using input data 

private:

void  FillTree(Particle*,int); // sets up daughter particles 

};

#endif

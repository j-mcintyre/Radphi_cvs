// CHEES header file

/* 
    3x3 matrix class

    Just a solution that does what I need it to do; not meant to be some great
    matrix class.

    PCF 23 Aug 1995

*/

#ifndef _3X3MATRIX_CLASS_INCLUDED_
#define _3X3MATRIX_CLASS_INCLUDED_

#include <Three.h>
#include <Four.h>

class SqMatrix3
{
public:

// Data  ////////////////////////////////////////////////////////////////

	double n[4][4];   /* the actual 3x3 array of data */
                          /* I don't use 0 here -- I start from 1 */

// Member Functions /////////////////////////////////////////////////////

  	SqMatrix3();
	SqMatrix3(Three,Three,Three);  /* make matrix using three vectors */
		
void      Set(Three,Three,Three); /* set matrix entries */

double 	  Determinant(void);             /* find determinant of the matrix  */
SqMatrix3 ReplaceColumn(int,Three);      /* returns a matrix obtained by
                                            replacing a column of the current
                                            matrix.  Useful for Cramer's rule.
                                         */
SqMatrix3 operator*(SqMatrix3 &);  /* matrix times matrix */
Four 	  operator*(Four &); /* matrix times 4-vector. Ignores t component */
Three 	  operator*(Three &); /* matrix times vector */

void print(void);                  /* print out matrix */

};

#endif

    

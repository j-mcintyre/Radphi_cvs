/* CHEES header file */

/* 
   fileStructs.h

   Various output file formats that Chees uses (see Keeper.C)

*/

#ifndef _CHEES_FILE_STRUCTS_
#define _CHEES_FILE_STRUCTS_

typedef struct
{
         int            numVectors;
unsigned int            globalAcceptBits;
         vector3_t      principleVertex;

} c_file_header;

typedef struct
{
         dvector4_t     vec;
unsigned int            PersonalAcceptBits;

} c_file_data;

#endif

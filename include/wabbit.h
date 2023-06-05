#ifndef WABBIT_H_INCLUDED
#define WABBIT_H_INCLUDED

#define MAX_TUPLE 100

#define TYPE_HISTO 0x1
#define TYPE_TUPLE 0x2
#define TYPE_HEADER 0x4
#define TYPE_DATA 0x8

#define TYPE_READ 0x1
#define TYPE_WRITE 0x2

#define INT 0x1
#define FLOAT 0x2
#define VECTOR_3 0x4
#define VECTOR_4 0x8
#define ARRAY 0x1000
#define ARRAY_VAR 0x2000
#define ARRAY_FIXED 0x4000
#define RANGED 0x8000

#define FIT_NO_LIMIT 0x1
#define GAUSSIAN 0x100000
#define USER_GUESS 0x100

typedef struct{
  int nBins;
  int nEntries;
  float nOver[8];
  float data[1];
}wabbitData_t;

/*

  Overflow regions :

     0 | 1 | 2
     _________
     7 |   | 3
     _________
     6 | 5 | 4


*/

typedef struct{
  int  id;
  char name[80];
  int nBins1;
  float lowerLimit1;
  float upperLimit1;
  int nBins2;
  float lowerLimit2;
  float upperLimit2;
  float max;
} wabbitHeader_t;

typedef struct{
  int id;
  char name[80];
  int type;
  int size;
  int index;
  int lowRange;
  int highRange;
  int nChunks;
  int chunkSz;
  void *address;
} tupleElement_t;

typedef struct{
  int id;
  int type;
  int nElements;
  int size;
  char name[80];  
  char fname[80];
} tupleHeader_t;

typedef struct{
  double value;
  double paramError;
  double minusError;
  double plusError;
  double globalCorCo;
} bunnyFit_t;

int bbook1(int id,char *name,int nBins,float lowerLimit,float upperLimit,float max);

int bbook2(int id,char *name,int nBins1,float lowerLimit1,float upperLimit1,
	   int nBins2, float lowerLimit2,float upperLimit2,float max);

int bz1(int id); /* zero the histogram */
int bf1(int id,float value,float weight);
int bf2(int id,float value1,float value2,float weight);
int bbsave(int id, char *fname);


int bbGetHisto(int id,char *filename,wabbitHeader_t *header,wabbitData_t **data);
/* Note : User must free(*data) when done */


extern int bbfit_mninit1;  /* minuit output channels */
extern int bbfit_mninit2;
extern int bbfit_mninit3;

int bbfit(int id, bunnyFit_t *result, double **covar, int flag, double low, double high);


int bbnt(int id, char *name, char *file,int bufsize);
int bbname(int id,void *address,char *name,int type);
int bfnt(int id);
int btout(int id);
int entrySize(tupleElement_t *elementList,int entry);

int bbGetTuple(int id);
int bbGetFormat(int id,tupleHeader_t *header,tupleElement_t **list);
int bbntr(int id,char *fname);
int bbntc(int id);
#endif


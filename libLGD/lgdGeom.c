#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map_manager.h>
#include <lgdGeom.h>
#include <encoderHandle.h>

static int nChannels;
static int nRows;
static int nCols;
static int nHoles;
static int hole[LGDGEOM_MAXHOLES];
static int encoderZero[2];
static float encoderPitch[2];
static float xOffset,yOffset,zOffset;
static float CountsPerRot = 512.0;
static float zPos;
static float x0;
static float y0;
static int centerChannel;
static float centerOffset[2]={0,0};

static encoderData_t encoders;
static int defaultXporter[2] = {0,0};

static int *isAnEdge=NULL;
static int **neighborsSave=NULL;
static int *nNeighborsSave=NULL;
static int *rowSave=NULL;
static int *colSave=NULL;
static int *holeSave=NULL;
static vector3_t *posSave=NULL;
static int **byRowCol= NULL;
static float *blockSize=NULL;


static void lgdGeomFindEdges(void);
static void lgdGeomFindNeighbors(void);
static void lgdBuildRowCol(void);
static void lgdBuildLocalCoord(void);
static void lgdGeomBuildByRowCol(void);
static void lgdGeomBuildHoles(void);

int lgdGeomSetup(int runNum)

     /* Sets up the geometry package. Reads data in from the map and
	builds some usefull arrays */

{
  int ret;
  int index;
  int time;
  char *dir;
  char *fileName = "maps/lgdGeom.map";
  char *mapFile;
  
  float *xtmp= NULL;
  float *ytmp = NULL;


  /* For MC users */

  runNum = abs(runNum);

  dir = getenv("RADPHI_CONFIG");

  if (dir == NULL)
    {
      fprintf(stderr,"lgdGeomSetup: Error: RADPHI_CONFIG is not defined. Exiting.\n");
      exit(1);
    }

  mapFile = malloc(strlen(dir) + strlen(fileName) +2);
  strcpy(mapFile,dir);
  strcat(mapFile,"/");
  strcat(mapFile,fileName);

  /* Free up storage if it has already been used */

  if(isAnEdge)
    free(isAnEdge);
  if(holeSave)
    free(holeSave);
  if(nNeighborsSave)
    free(nNeighborsSave);
  if(rowSave)
    free(rowSave);
  if(colSave)
    free(colSave);
  if(posSave)
    free(posSave);
  if(neighborsSave){
    for(index=0;index<nChannels;index++)
      free(neighborsSave[index]);
    free(neighborsSave);
  }
  if(byRowCol){
    for(index=0;index<nRows;index++)
      free(byRowCol[index]);
    free(byRowCol);
  }
  
  /* Read the Map */

  if(map_get_int(mapFile,"geometry","nChannels",1,&nChannels,runNum,&time))
    return(LGDGEOM_READERR);
  if(map_get_int(mapFile,"geometry","nRows",1,&nRows,runNum,&time))
     return(LGDGEOM_READERR);
  if(map_get_int(mapFile,"geometry","nCols",1,&nCols,runNum,&time))
    return(LGDGEOM_READERR);
  if(map_get_int(mapFile,"geometry","nHoles",1,&nHoles,runNum,&time))
    return(LGDGEOM_READERR);
  if(map_get_int(mapFile,"geometry","hole",LGDGEOM_MAXHOLES,hole,runNum,&time))
    return(LGDGEOM_READERR);
  if(map_get_int(mapFile,"geometry","encoderZero",2,encoderZero,runNum,&time))
    return(LGDGEOM_READERR);
  if(map_get_float(mapFile,"geometry","encoderPitch",2,encoderPitch,runNum,&time))
    return(LGDGEOM_READERR);
  if(map_get_float(mapFile,"geometry","zPos",1,&zPos,runNum,&time))
    return(LGDGEOM_READERR);
  if(map_get_float(mapFile,"geometry","x0",1,&x0,runNum,&time))
    return(LGDGEOM_READERR);
  if(map_get_float(mapFile,"geometry","y0",1,&y0,runNum,&time))
    return(LGDGEOM_READERR);
  if(map_get_int(mapFile,"geometry","centerChannel",1,&centerChannel,runNum,&time))
    return(LGDGEOM_READERR);
  if(map_get_float(mapFile,"geometry","centerOffset",2,centerOffset,runNum,&time))
    return(LGDGEOM_READERR);
  if(getEncoderData(&encoders,abs(runNum)))
    return(LGDGEOM_READERR); 
  

  /* Allocate storage */

  isAnEdge = malloc(nChannels*sizeof(int));
  nNeighborsSave = malloc(nChannels*sizeof(int));
  neighborsSave = malloc(nChannels*sizeof(int *));
  for(index=0;index<nChannels;index++)
    neighborsSave[index] = malloc(8*sizeof(int));
  rowSave = malloc(nChannels*sizeof(int));
  colSave = malloc(nChannels*sizeof(int));
  posSave = malloc(nChannels*sizeof(vector3_t));
  holeSave = malloc(nChannels*sizeof(int));
  memset(holeSave,0,sizeof(int)*nChannels);
  
  lgdGeomBuildHoles();

  lgdLocateXporter(&encoders);

  /* Fill some handy arrays */

  if(nChannels != 324){
    lgdBuildRowCol();
    lgdBuildLocalCoord();
    lgdGeomFindNeighbors();
    lgdGeomFindEdges();
    return(LGDGEOM_OK);
  }
  else{
    int array[324][8];
    int channel =0;
    
    if(map_get_int(mapFile,"geometry","p2row",nChannels,rowSave,runNum,&time))
      return(LGDGEOM_READERR);
    if(map_get_int(mapFile,"geometry","p2col",nChannels,colSave,runNum,&time))
      return(LGDGEOM_READERR);
    xtmp = malloc(nChannels*sizeof(float));
    ytmp = malloc(nChannels*sizeof(float));
    if(map_get_float(mapFile,"geometry","p2xChannel",nChannels,xtmp,runNum,&time))
      return(LGDGEOM_READERR);
    if(map_get_float(mapFile,"geometry","p2yChannel",nChannels,ytmp,runNum,&time))
      return(LGDGEOM_READERR);
    for(index=0;index<nChannels;index++){
      posSave[index].x = xtmp[index];
      posSave[index].y = ytmp[index];
      posSave[index].z = 0.0;      
    }
    free(xtmp);
    free(ytmp);
    xtmp = NULL;
    ytmp = NULL;
    if(map_get_int(mapFile,"geometry","p2nNeighbors",nChannels,nNeighborsSave,runNum,&time))
      return(LGDGEOM_READERR);
    if(map_get_int(mapFile,"geometry","p2neighbor",nChannels*8,&array[0][0],runNum,&time))
      return(LGDGEOM_READERR);
    for(channel=0;channel<nChannels;channel++){
      for(index=0;index<8;index++)
	neighborsSave[channel][index] = array[channel][index];
    }    
    lgdGeomFindEdges();
    lgdGeomBuildByRowCol();
    return(LGDGEOM_OK);
  }
    
}
  
static void lgdGeomBuildHoles(void)
{
  int index;

  for(index=0;index<nHoles;index++)
    holeSave[hole[index]] = LGDGEOM_TRUE;
    
}
      
static void lgdGeomBuildByRowCol(void)
{
  int row;
  int col;
  int channel;
  int index;

  byRowCol = malloc(nRows*sizeof(int *));
  for(row=0;row<nRows;row++)
    byRowCol[row] = malloc(nCols * sizeof(int));
  for(row=0;row<nRows;row++){
    for(col=0;col<nCols;col++){
      for(channel=0;channel<nChannels;channel++)
	if( rowSave[channel] == row && colSave[channel] == col)
	  byRowCol[row][col] = channel;
    }
  }
}

void lgdGetGeomParam(int *nChannelsU,int *nRowsU, int *nColsU)
{
  *nChannelsU = nChannels;
  *nRowsU = nRows;
  *nColsU = nCols;
  return;
}

static void lgdBuildRowCol(void)
{
  int channel;
  int index;

  for(channel=0;channel<nChannels;channel++){
    rowSave[channel] = channel/nCols;
    colSave[channel] = channel%nCols;
  }
}

static void lgdBuildLocalCoord(void)

{
  int channel;

  for(channel=0;channel<nChannels;channel++){
    posSave[channel].x = (colSave[channel] - colSave[centerChannel])*4.0 + centerOffset[0];
    posSave[channel].y = (rowSave[channel] - rowSave[centerChannel])*4.0 + centerOffset[1];
    posSave[channel].z = 0.0;
  }
}

static void lgdGeomFindNeighbors(void)
     
     /* fills the neighbor arrays. nNeighborsSave[channel] is the number of
	neighbors. neighborsSave[channel][index] contains nNeighbor items,
	the channles of the valid neighbors */

{
  lgdGeom_t current;
  lgdGeom_t neighbor;
  
  for(current.channel=0;current.channel<nChannels;current.channel++){
    lgdGetGeom(0,&current,LGDGEOM_BYCHANNEL);
    nNeighborsSave[current.channel] = 0;
    neighbor.row = current.row+1;
    neighbor.col = current.col-1;
    if(!lgdGetGeom(0,&neighbor,LGDGEOM_BYROWCOL)){
      neighborsSave[current.channel][nNeighborsSave[current.channel]] = neighbor.channel;
      nNeighborsSave[current.channel]++;
    }
    neighbor.row = current.row+1;
    neighbor.col = current.col;
    if(!lgdGetGeom(0,&neighbor,LGDGEOM_BYROWCOL)){
      neighborsSave[current.channel][nNeighborsSave[current.channel]] = 
	neighbor.channel;
      nNeighborsSave[current.channel]++;
    }
    neighbor.row = current.row+1;
    neighbor.col = current.col+1;
    if(!lgdGetGeom(0,&neighbor,LGDGEOM_BYROWCOL)){
      neighborsSave[current.channel][nNeighborsSave[current.channel]] = 
	neighbor.channel;
      nNeighborsSave[current.channel]++;
    }
    neighbor.row = current.row;
    neighbor.col = current.col-1;
    if(!lgdGetGeom(0,&neighbor,LGDGEOM_BYROWCOL)){
      neighborsSave[current.channel][nNeighborsSave[current.channel]] = 
	neighbor.channel;
      nNeighborsSave[current.channel]++;
    }
    neighbor.row = current.row;
    neighbor.col = current.col+1;
    if(!lgdGetGeom(0,&neighbor,LGDGEOM_BYROWCOL)){
      neighborsSave[current.channel][nNeighborsSave[current.channel]] = 
	neighbor.channel;
      nNeighborsSave[current.channel]++;
    }
    neighbor.row = current.row-1;
    neighbor.col = current.col-1;
    if(!lgdGetGeom(0,&neighbor,LGDGEOM_BYROWCOL)){
      neighborsSave[current.channel][nNeighborsSave[current.channel]] = 
	neighbor.channel;
      nNeighborsSave[current.channel]++;
    }
    neighbor.row = current.row-1;
    neighbor.col = current.col;
    if(!lgdGetGeom(0,&neighbor,LGDGEOM_BYROWCOL)){
      neighborsSave[current.channel][nNeighborsSave[current.channel]] = 
	neighbor.channel;
      nNeighborsSave[current.channel]++;
    }
    neighbor.row = current.row-1;
    neighbor.col = current.col+1;
    if(!lgdGetGeom(0,&neighbor,LGDGEOM_BYROWCOL)){
      neighborsSave[current.channel][nNeighborsSave[current.channel]] = 
	neighbor.channel;
      nNeighborsSave[current.channel]++;
    }
  }
}

static void lgdGeomFindEdges(void)

     /* fills the isAnEdge array: LGDGEOM_TRUE if the module is an edge
	or a hole. */

{
  lgdGeom_t block;
  int neighbors[8];
  int nNeighbors;
  
  for(block.channel=0;block.channel<nChannels;block.channel++){
    lgdGetGeom(0,&block,LGDGEOM_BYCHANNEL);
    lgdGeomGetNeighbors(block.channel,neighbors,&nNeighbors);
    if(nNeighbors != 8){
      isAnEdge[block.channel] = LGDGEOM_TRUE;
    }
  }
}

int lgdGeomIsAnEdge(int channel)
     
     /* Returns LGDGEOM_TRUE if the channel is an edge or a hole */

{
  if(channel >=0 && channel < nChannels)
    return(isAnEdge[channel]);
  else
    return(LGDGEOM_CHANNELERROR);
}


int lgdGeomGetNeighbors(int channel, int *neighbor,int *nNeighbors)

     /* Fills <neighbor> and sets nNeighbors */

{
  if(channel >= 0 && channel < nChannels){
    *nNeighbors = nNeighborsSave[channel];
    memcpy(neighbor,neighborsSave[channel],8*sizeof(int));
    return(LGDGEOM_OK);
  }
  else
    return(LGDGEOM_CHANNELERROR);
}

int lgdGetChannelRaw(unsigned long ADCWord, int *channel)

{
  int row,col;
  
  col = ((unsigned int)ADCWord & 0xff000000) >> 24;
  row = ((unsigned int)ADCWord & 0xff0000) >> 16;
  if(row < nRows && col < nCols){
    *channel = (row * nCols) + col;
    return(LGDGEOM_OK);
  }
  else{
    return(LGDGEOM_RAWERROR);
  }
}

int lgdGetGeom(unsigned long ADCChannel,lgdGeom_t *item,int index)

     /* Fills <item> based on <index>, the indexing scheme */

{
  int loop;

  switch(index){
  case LGDGEOM_BYRAWADC:
    if(nChannels == 324)
      return(LGDGEOM_RAWERROR);
    item->col = ((unsigned int)ADCChannel & 0xff00) >> 8;
    item->row = ((unsigned int)ADCChannel & 0xff);
    if(item->row > nRows-1 || item->col > nCols-1){
      item->row = -1;
      item->col = -1;
      item->channel = -1;
      item->hole = LGDGEOM_TRUE;
      return(LGDGEOM_RAWERROR);
    }
    item->channel = (item->row * nCols) + item->col;
    break;
  case LGDGEOM_BYROWCOL:
    if(item->row > nRows-1 || item->col > nCols-1 ||
       item->row < 0 || item->col < 0){
      item->row = -1;
      item->col = -1;
      item->channel = -1;
      item->hole = LGDGEOM_TRUE;
      return(LGDGEOM_ROWCOLERROR);    
    }
    else
      if(nChannels == 324)
	item->channel = byRowCol[item->row][item->col];
      else
	item->channel = (item->row * nCols) + item->col;
    break;
  case LGDGEOM_BYCHANNEL:
    if( (item->channel < 0) || (item->channel > nChannels-1) ){
      item->row = -1;
      item->col = -1;
      item->channel = -1;
      item->hole = LGDGEOM_TRUE;
      return(LGDGEOM_CHANNELERROR);
    }
    item->col = colSave[item->channel];
    item->row = rowSave[item->channel];
    break;   
  default:
    return(LGDGEOM_METHODERR);
  }
  lgdLocalCoord(item->channel,&item->space);
  if(holeSave[item->channel]){
    item->hole = LGDGEOM_TRUE;
    return(LGDGEOM_ISAHOLE);
  }
  item->hole = LGDGEOM_FALSE;
  return(LGDGEOM_OK);
}

float lgdLocalX(int index)

{
  if(index >=0 && index < nChannels)
    return(posSave[index].x);
  else
    return(LGDGEOM_COORDERR);
}

float lgdLocalY(int index)

{
  if(index >=0 && index < nChannels)
    return(posSave[index].y);
  else
    return(LGDGEOM_COORDERR);
}

float lgdLocalZ(int index)

{
  if(index >=0 && index < nChannels)
    return(posSave[index].z);
  else
    return(LGDGEOM_COORDERR);
}

int lgdLocalCoord(int index, vector3_t *local)
{
  if((index<0) || (index >= nChannels))
    return(LGDGEOM_CHANNELERROR);

  (*local) = posSave[index];
  return(LGDGEOM_OK);
}

void lgdLocateXporter(encoderData_t *encoders)

{

  /* x0(y0) is the offset between local and Global coord when 
     encoder is at encoderZero 
     Note that increasing y encoder => increasing y
               increasing x encoder => decreasing x
  */

  xOffset = x0-encoderToLength(encoders,ENC_LGD_X);
  yOffset = encoderToLength(encoders,ENC_LGD_Y)-y0;
  zOffset = zPos;
}

void lgdLocalToGlobalCoord(vector3_t local, vector3_t *mps)
{
  mps->x = local.x + 0; /* xOffset not correctly set right now */
  mps->y = local.y + 0; /* yOffset not correctly set right now */
  mps->z = local.z + zOffset;
}

void lgdGlobalToLocalCoord(vector3_t *local, vector3_t mps)
{
  local->x = mps.x - 0; /* xOffset not correctly set right now */
  local->y = mps.y - 0; /* yOffset not correctly set right now */
  local->z = mps.z - zOffset;
}

int lgdGlobalCoord(int index, vector3_t *mps)

{  
  if(index >=0 && index < nChannels){
    lgdLocalToGlobalCoord(posSave[index],mps);
    return(LGDGEOM_OK);
  }
  else
    return(LGDGEOM_CHANNELERROR);
}

void lgdModuleSize(vector3_t *block)
{
  block->x = 4.0;
  block->y = 4.0;
  block->z = 40.0;
}
void lgdBlockSize(int channel,vector3_t *block)
{
  block->x = 4.0;
  block->y = 4.0;
  block->z = 45.0;
}

int lgdLocalPositionToRC(vector3_t local,int *row, int*col)
{
  static vector3_t blockSize;
  static init=0;
  /*  if(!init){
    lgdBlockSize(&blockSize);
    init=1;
  }*/
  
  *col=(int)((local.x)/4.0) + colSave[centerChannel];
  if(local.x < 0) (*col)--;
  *row=(int)((local.y)/4.0) + rowSave[centerChannel];
  if(local.y < 0) (*row)--;
  return(1);
}

int lgdLocalPositionToGeomType(vector3_t local, lgdGeom_t *item)
{
  static vector3_t moduleSize;
  static init=0;
  if(!init){
    lgdModuleSize(&moduleSize);
    init=1;
  }
  
  item->col=(int)((local.x)/4.0) + colSave[centerChannel];
  if(local.x < 0) (item->col)--;
  item->row=(int)((local.y)/4.0) + rowSave[centerChannel];
  if(local.y < 0) (item->row)--;
  return(lgdGetGeom(0,item,LGDGEOM_BYROWCOL));
}
  
int lgdGlobalPositionToGeomType(vector3_t global,lgdGeom_t *item)
{
  vector3_t local;
  lgdGlobalToLocalCoord( &local, global);
  return(lgdLocalPositionToGeomType(local,item));
}

/* end file */

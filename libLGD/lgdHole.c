#include <lgdHole.h>
#include <ntypes.h>
#include <math.h>
#include <string.h>

#define SQR(x) (x)*(x) 

/*
 * Returns the distance that path between a block and cluster 
*  passes through the lgd hole.
 */
double compressHole(const vector3_t *block, const vector3_t *cluster) {

  vector3_t point[2];
  memset(point,0,sizeof(float)*2*3);
/* Do nothing if cluster is inside the hole */
  if ( ((cluster->x > HOLE_X_LOW) && (cluster->x < HOLE_X_HIGH)) &&
       ((cluster->y > HOLE_Y_LOW) && (cluster->y < HOLE_Y_HIGH)) ) return 0;
/*
* First, get the intersection of the line between the cluster and the block
* with the hole's wall paralel to y-axis (x fixed). If the y-coordinate
* is within hole boundaries take it.
*/
  int count=0;
  float y,Cx;
  for (Cx = HOLE_X_LOW; Cx <= HOLE_X_HIGH; Cx += (HOLE_X_HIGH - HOLE_X_LOW) ) {
    y = getY_Xfix(Cx,block,cluster);
    if ((y >= HOLE_Y_LOW) && (y <= HOLE_Y_HIGH)) {
       point[count].x = Cx;
       point[count].y = y;
       count++;
    }
  }
/*
* If the cluster-block line goes through the hole return the distance.
If not, check for crossing between walls pararle to x-axis.
*/
  if (count == 2)  
     return sqrt(SQR(point[0].x-point[1].x)+SQR(point[0].y-point[1].y));
  
  float x,Cy;
  for (Cy = HOLE_Y_LOW; Cy <= HOLE_Y_HIGH; Cy += (HOLE_Y_HIGH - HOLE_Y_LOW) ) {
    x = getX_Yfix(Cy,block,cluster); 
    if ((x >= HOLE_X_LOW) && (x <= HOLE_X_HIGH)) {
      point[count].x = x;
      point[count].y = Cy;
      count++;
    }    
  }
  if (count == 2) 
     return sqrt(SQR(point[0].x-point[1].x)+SQR(point[0].y-point[1].y));
  
  return 0;
}

/* 
 *  Return X coordinate of the intersection between the line defined by
 *  two ponts in XY space (say block and cluster) and fixed Y line. 
 */
float getX_Yfix(float const yFix, const vector3_t *block, const vector3_t *clust) {

  if (block->y == clust->y) {
    return -10000.; 
  }
  else {
    return (block->x - clust->x)*(yFix - clust->y)/(block->y - clust->y) + clust->x;   
  }
}

/* 
 *  Return Y coordinate of the intersection between the line defined by
 *  two ponts in XY space (say block and cluster) and fixed X line. 
 */
float getY_Xfix(float const xFix, const vector3_t *block, const vector3_t *clust) {

  if (block->x == clust->x) {
    return -10000.; /* Some unreasonable value for the line that is paralel */
  }                 /* with fixed X line*/
  else {
    return (block->y - clust->y)*(xFix - clust->x)/(block->x - clust->x) + clust->y;
  }
}

int testhole() {
 float x0,y0,x,y;

 return 0;
}

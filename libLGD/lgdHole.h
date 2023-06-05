#include <ntypes.h>

#define HOLE_X_LOW -4.0
#define HOLE_X_HIGH 4.0
#define HOLE_Y_LOW -4.0
#define HOLE_Y_HIGH 4.0

float getY_Xfix(float const xFix, const vector3_t *block, const vector3_t *clust);
float getX_Yfix(float const yFix, const vector3_t *block, const vector3_t *clust);
double compressHole(const vector3_t *block, const vector3_t *center);                                                           



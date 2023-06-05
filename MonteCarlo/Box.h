// CHEES header file 
// header file for Box class

#ifndef _BOX_CLASS_INCLUDED_
#define _BOX_CLASS_INCLUDED_

typedef struct {double x, y; } point;
typedef struct {point lgd, dea, csi, beam; } position;


class Box{
 public:
  Box(double = 0, double = 0, double = 0, double = 0);
  int RectCut(point *);
 private:
  double left;
  double right;
  double top;
  double bottom;
};


#endif

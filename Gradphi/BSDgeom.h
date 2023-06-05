/* this is information on the BSD for use by c routines */
#define BSDscintHeight 70.0
/* here is the z coordinate of the downstream end of the BSD Scints.*/
#define BSDz0 -35.0
/* outer ring */
#define NumScintsRing3 24
#define NumScintsRing2 12
#define NumScintsRing1 12
/*these are how much each scint twists from the upstream end to the
  downstream end */
/*the angle is measured as if looking downstream. I.e.
 phi = phi0+twist*(z-z0)/zmax, where z+ is downstream */
#define BSDtwistRing3 0.0
#define BSDtwistRing2 -120.0
#define BSDtwistRing1 120.0
/* here is the angles of the start of the first strip in every ring */
//#define BSDphi0Ring1  0.0
#define BSDphi0Ring1 15.0
//#define BSDphi0Ring2 30.0
#define BSDphi0Ring2 15.0
#define BSDphi0Ring3 15.0

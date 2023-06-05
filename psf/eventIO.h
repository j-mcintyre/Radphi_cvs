/* public data members */

/* public function members (c interface) */
int eventGetNext(void);
/* get characteristic LGD kinematics */
float eventGetLgdEtot(void);
float eventGetLgdEmax(int *ch);
float eventGetLgdE(int ch);
float eventGetLgdErc(int row, int col);
/* get shower characteristic variables */
float eventGetEtot(void);
float eventGetEmax(int *ch);
float eventGetE(int ch);
float eventGetErc(int row, int col);
/* get X,Y and R,C from LGD block (channel) */
void eventGetLgdXY(int row, int col, float *x, float *y);
void eventGetLgdRC(int ch, int *row, int *col);

/*  public function members (fortran interface) */
int eventtest_(int *count);
void eventclose_(void);
int eventnext_(void);
int eventmc_(void);
void eventimage_(void);
float eventphoton_(float *Px, float *Py, float *Pz);
float mcphoton_(float *Px, float *Py, float *Pz);
int  eventrun_();
void eventdef_(int *defin);

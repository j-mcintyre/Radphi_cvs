#include <Chees.h>


extern "C"{
void CutCheesWrap(char *file,mc_param_t *param,particles_db_t **part);
int GenerateOneWrap(cheesInfo_t *info);
void configureLgdGeometryWrap(int,int,int);
void setCheesRunNoWrap(int);
}

void CutCheesWrap(char *file,mc_param_t *param,particles_db_t **part)
{
  CutChees(file,param,part);
}

int GenerateOneWrap(cheesInfo_t *info)
{
  return(GenerateOne(info));
}


void setCheesRunNoWrap(int runNo)
{
  setCheesRunNo(runNo);
}

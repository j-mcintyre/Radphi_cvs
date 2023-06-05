
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>

#include <param.h>
#include <mc_param.h>
#include <assert.h>

static char rcsid[] = "$Id: loadParam.c,v 1.5 1999/02/20 21:59:07 radphi Exp $";


/*
   $Log: loadParam.c,v $
   Revision 1.5  1999/02/20 21:59:07  radphi
   Changes made by crsteffe@dustbunny
   updating for >3 photon final states

 * Revision 1.4  1997/02/12  20:46:41  lfcrob
 * Fixed up CEBAF geom
 *
 * Revision 1.3  1997/01/29  15:47:13  lfcrob
 * Fixed loading of BREM beam parameters
 *
   */
#define NAN sqrt(-1)

int loadBeam(mc_beam_t *beam,int *);
int loadTarget(int *target,int *fermi);
int loadProduct(int *product);
int loadRecoil(int *recoil);
int loadTransfer(mc_transfer_t *transfer);
int loadCodes(mc_codes_t *codes);
int loadParticle(char *name);
void setParticleID(char *name,Particle_t *);
void setForceReload(int force);
float param_getFloatV(char *base,char *name, float def); 
float param_getFloatVI(char *base,char *name,int i,float def);
int param_getIntV(char *base,char *name, int def); 
int param_getIntVI(char *base,char *name, int i,int  def);
char *param_getValueV(char *base,char *name);
char *param_getValueVI(char *base,char *name,int i);
int param_getBoolean(char *name);
int param_getBooleanV(char *base,char *name);
void cleanString(char *);
void printParticle(particle_db_t *part,particles_db_t *particleDb);
int loadDecayAngles(particle_db_t *part);
int readHistogramDef(char *root,histogramDef_t *histo);

int checkBoolean(int value);
int checkFloat(float value1);
int checkInt(int value1);

particles_db_t *particleDb=NULL;

static int maxParticles=0;
static int forceReload=0;

static void *overrideFile=NULL;

void flushParam(void)
{
  if(particleDb)
    particleDb->nParticles=0;
  param_clean();
}

int cheesParamAddFile(char *filename)

{
  size_t len;

  len = strlen(filename);
  if(len){
    overrideFile = calloc(sizeof(char),len+1);
    strcpy(overrideFile,filename);
    return(0);
  }
  return(1);
}

int initParam(char *userFile)

{
  char *libDir=NULL;
  char filename[1024];

  if((libDir = getenv("RADPHI_CONFIG")) == NULL){
    fprintf(stderr,"Environment variable RADPHI_CONFIG undefined.\n");
    return(PARAM_ERROR);
 }
  
  /* load the master database */

  strcpy(filename,libDir);
  strcat(filename,"/MonteCarlo/particles.xrdb");
/*  sprintf(filename,"%s/%s",libDir,PARTICLE_DB);*/

  if(param_loadFile(filename)){
    fprintf(stderr,"Unable to load master database %s.\n",filename);
    return(PARAM_ERROR);
  }
  
  /* load the user database */

  
  if(param_loadFile(userFile)){
    fprintf(stderr,"Unable to load user database %s.\n",userFile);
    return(PARAM_ERROR);
  }

  /* Load the override database, if given */

  if(overrideFile)
    if(param_loadFile(overrideFile)){
      fprintf(stderr,"Unable to load override database %s.\n",overrideFile);
      return(PARAM_ERROR);
    }
  
  return(PARAM_OK);
}

int loadParam(mc_param_t *param,particles_db_t **userParticleDb)
{
  char *name="pi+";
  int index;
  int ret=0;

  if(!param){
    fprintf(stderr,"loadParam : bad pointer to mc_param_t\n");
    return(PARAM_ERROR);
  }

  if(!userParticleDb){
    fprintf(stderr,"loadParam : bad pointer to particles_db_t\n");
    return(PARAM_ERROR);
  }

  ret |= loadBeam(&param->beam,&param->reaction.beam);
  ret |= loadTarget(&param->reaction.target,&param->codes.fermiTarget);
  ret |= loadProduct(&param->reaction.product);
  ret |= loadRecoil(&param->reaction.recoil);
  ret |= loadTransfer(&param->transfer);
  ret |= loadCodes(&param->codes);
  *userParticleDb = particleDb;
  return(ret);
}

int loadCodes(mc_codes_t *codes)

{
  char *tmp;
  int ret=0;

  codes->nIterations = param_getInt("setup.nIterations",0);
  ret |= checkInt(codes->nIterations);
  codes->doHisto = param_getBoolean("setup.doHistos");
  ret |= checkBoolean(codes->doHisto);
  codes->returnAll = param_getBoolean("setup.returnAll");
  ret |= checkBoolean(codes->returnAll);
  codes->separationCut = (float)param_getDouble("setup.separationCut",NAN);
  ret |= checkFloat(codes->separationCut);
  codes->energyCut = (float)param_getDouble("setup.energyCut",NAN);
  ret |= checkFloat(codes->energyCut);
  codes->doEnergySmear=param_getBoolean("setup.doEnergySmear");
  ret |= checkBoolean(codes->doEnergySmear);
  codes->doPositionSmear=param_getBoolean("setup.doPositionSmear");
  ret |= checkBoolean(codes->doPositionSmear);
  
  if(codes->doEnergySmear){
    codes->energySmearFloor = (float)param_getDouble("setup.energySmear.floor",NAN);
    ret |= checkFloat(codes->energySmearFloor);
    codes->energySmearStat = (float)param_getDouble("setup.energySmear.stat",NAN);
    ret |= checkFloat(codes->energySmearStat);
  }
  if(codes->doPositionSmear){
    codes->positionSmear[0] = (float)param_getDouble("setup.x.positionSmear",NAN);
    ret |= checkFloat(codes->positionSmear[0]);
    codes->positionSmear[1] = (float)param_getDouble("setup.y.positionSmear",NAN);
    ret |= checkFloat(codes->positionSmear[1]);
  }


  tmp = (char*)param_getValue("setup.geometry");
  if(!tmp)
    return(PARAM_ERROR);
  cleanString(tmp);
  codes->geometry = 0;
  if(!strcmp(tmp,"E852"))
    codes->geometry = GEOM_E852;
  if(!strcmp(tmp,"CEBAF_328"))
    codes->geometry = GEOM_CEBAF_328;
  if(!strcmp(tmp,"CEBAF_620"))
    codes->geometry = GEOM_CEBAF_620;
  if(!strcmp(tmp,"GAMS2000"))
    codes->geometry = GEOM_GAMS2000;
  if(!strcmp(tmp,"GAMS4000"))
    codes->geometry = GEOM_GAMS4000;
  if(!codes->geometry){
    fprintf(stderr,"loadCodes : Geometry undefined\n");
    return(PARAM_ERROR);    
  }

  tmp = (char*)param_getValue("setup.output");
  if(!tmp)
    return(PARAM_ERROR);
  cleanString(tmp);
  codes->outFileType = 0;
  if(!strcmp(tmp,"NONE"))
    codes->outFileType = OUT_NONE;
  if(!strcmp(tmp,"C"))
    codes->outFileType = OUT_C;
  if(!strcmp(tmp,"C++"))
    codes->outFileType = OUT_CPLUSPLUS;
  if(!strcmp(tmp,"FORTRAN"))
    codes->outFileType = OUT_FORTRAN;
  if(!strcmp(tmp,"ASCII"))
    codes->outFileType = OUT_ASCII;
  if(!codes->outFileType){
    fprintf(stderr,"loadCodes : Output type undefined\n");
    return(PARAM_ERROR);    
  }
  if(codes->outFileType != OUT_NONE){
    tmp = (char *)param_getValue("setup.outputFile");
    if(!tmp)
      return(PARAM_ERROR);
    cleanString(tmp);
    strcpy(codes->outFile,tmp);
  }
  return(ret);
}

int loadTransfer(mc_transfer_t *transfer)
{
  char *tmp;
  char *tmpName;
  int ret=0;

  tmp = (char *)param_getValue("transfer.type");
  if(!tmp)
    return(PARAM_ERROR);
  transfer->transferCode=0;
  cleanString(tmp);
  if(!strcmp(tmp,"TDIST_1"))
    transfer->transferCode = TDIST_1;
  if(!strcmp(tmp,"TDIST_2"))
    transfer->transferCode = TDIST_2;
  if(!strcmp(tmp,"TDIST_3"))
    transfer->transferCode = TDIST_3;
  if(!strcmp(tmp,"USER"))
    transfer->transferCode = TDIST_USER;
  if(!strcmp(tmp,"DELTA"))
    transfer->transferCode = TDIST_DELTA;
  if(!strcmp(tmp,"FLAT"))
    transfer->transferCode = TDIST_FLAT;
  switch(transfer->transferCode){
  case TDIST_1:
    transfer->exponent = (float)param_getDouble("transfer.exponent",NAN);
    ret |= checkFloat(transfer->exponent);
    break;
  case TDIST_2:
    transfer->exponent = (float)param_getDouble("transfer.exponent",NAN);
    ret |= checkFloat(transfer->exponent);
    break;
  case TDIST_3:
    transfer->exponent = (float)param_getDouble("transfer.exponent",NAN);
    ret |= checkFloat(transfer->exponent);
    transfer->exchangeMass = (float)param_getDouble("transfer.exchangeMass",NAN);
    ret |= checkFloat(transfer->exchangeMass);
    transfer->spinPower = (float)param_getDouble("transfer.spinPower",NAN);
    ret |= checkFloat(transfer->spinPower);
    transfer->tPower = (float)param_getDouble("transfer.tPower",NAN);
    ret |= checkFloat(transfer->tPower);
    break;
  case TDIST_USER:
    ret |= readHistogramDef("transfer.histo",&(transfer->histo));
    break;
  case TDIST_FLAT:
    transfer->tWidth = (float)param_getDouble("transfer.t.width",NAN);
    ret |= checkFloat(transfer->tWidth);
  case TDIST_DELTA:
    transfer->t = (float)param_getDouble("transfer.t.value",NAN);
    ret |= checkFloat(transfer->t);
    break;
  default:
    fprintf(stderr,"Unknown transfer type : %s\n",tmp);
    return(PARAM_ERROR);
  }
  return(ret);
}

int loadTarget(int *target, int *fermi)

{
  char *tmp;

  tmp = (char *)param_getValue("target.type");
  if(!tmp)
    return(PARAM_ERROR);
  cleanString(tmp);
  *target = loadParticle(tmp);
  if(*target == -1)
    return(PARAM_ERROR);
  *fermi = param_getBoolean("target.fermi");
  return(checkBoolean(*fermi));  
}


int loadProduct(int *product)

{
  char *tmp;

  tmp = (char *)param_getValue("product.type");
  if(!tmp)
    return(PARAM_ERROR);
  cleanString(tmp);
  *product = loadParticle(tmp);
  if(*product == -1)
    return(PARAM_ERROR);
  return(PARAM_OK);
}

int loadRecoil(int *recoil)

{
  char *tmp;

  tmp = (char *)param_getValue("recoil.type");
  if(!tmp)
    return(PARAM_ERROR);
  cleanString(tmp);
  *recoil = loadParticle(tmp);
  if(*recoil == -1)
    return(PARAM_ERROR);
  return(PARAM_OK);
}

int loadBeam(mc_beam_t *beam, int *beamParticle)

{
  char *tmp;
  int ret=0;

  tmp = (char *)param_getValue("beam.type");
  if(!tmp)
    return(PARAM_ERROR);
  cleanString(tmp);
  *beamParticle = loadParticle(tmp);
  if(*beamParticle == -1)
    ret |= PARAM_ERROR;
  
  tmp = (char *)param_getValue("beam.type");
  if(!tmp){
    return(PARAM_ERROR);
  }
  cleanString(tmp);
  tmp = (char *)param_getValue("beam.momentum.type");
  if(!tmp)
    return(PARAM_ERROR);
  cleanString(tmp);
  if(!strcmp(tmp,"FIXED"))
    beam->pCode = BEAM_MOMENTUM_FIXED;
  if(!strcmp(tmp,"GAUSSIAN"))
    beam->pCode = BEAM_MOMENTUM_GAUSSIAN;
  if(!strcmp(tmp,"FLAT"))
    beam->pCode = BEAM_MOMENTUM_FLAT;
  if(!strcmp(tmp,"BREM"))
    beam->pCode = BEAM_MOMENTUM_BREM;
  if(!strcmp(tmp,"BNL_94"))
    beam->pCode = BEAM_MOMENTUM_BNL_94;
  switch(beam->pCode){
  case BEAM_MOMENTUM_BREM:
    beam->bremThetaMax = (float)param_getDouble("beam.momentum.bremThetaMax",NAN);
    ret |= checkFloat(beam->bremThetaMax);
  case BEAM_MOMENTUM_FIXED:
    beam->p = (float)param_getDouble("beam.momentum.value",NAN);
    ret |= checkFloat(beam->p);
    break;
  case BEAM_MOMENTUM_FLAT:
  case BEAM_MOMENTUM_GAUSSIAN:
    beam->p = (float)param_getDouble("beam.momentum.value",NAN);
    ret |= checkFloat(beam->p);
    beam->pWidth = (float)param_getDouble("beam.momentum.sigma",NAN);
    ret |= checkFloat(beam->pWidth);    
    break;
  case BEAM_MOMENTUM_BNL_94:
    return(ret);
  default:
    fprintf(stderr,"loadBeam : unknown momentum distribution - %s\n",tmp);
    return(PARAM_ERROR);
  }

   
  tmp = (char *)param_getValue("beam.spatial.type");
  if(!tmp)
    return(PARAM_ERROR);
  cleanString(tmp);
  if(!strcmp(tmp,"NONE"))
    beam->spreadCode=BEAM_SPATIAL_NONE;
  if(!strcmp(tmp,"CIRCLE"))
    beam->spreadCode=BEAM_SPATIAL_CIRCLE;
  if(!strcmp(tmp,"GAUSSIAN"))
    beam->spreadCode=BEAM_SPATIAL_GAUSSIAN;
  switch(beam->spreadCode){
  case BEAM_SPATIAL_NONE:
    break;
  case BEAM_SPATIAL_GAUSSIAN:
    beam->sigma[0] = (float)param_getDouble("beam.spatial.x.sigma",NAN);
    ret |= checkFloat(beam->sigma[0]);
    beam->sigma[1] = (float)param_getDouble("beam.spatial.y.sigma",NAN);
    ret |= checkFloat(beam->sigma[1]);
    break;
  case BEAM_SPATIAL_CIRCLE:
    beam->radius = (float)param_getDouble("beam.spatial.radius",NAN);
    ret |= checkFloat(beam->radius);
    break;
  default:
    fprintf(stderr,"loadBeam : unknown spatial distribution - %s\n",tmp);
    return(PARAM_ERROR);
  }

  if(beam->pCode != BEAM_MOMENTUM_BREM){
    tmp = (char *)param_getValue("beam.angular.type");
    if(!tmp)
      return(PARAM_ERROR);
    cleanString(tmp);
    if(!strcmp(tmp,"NONE"))
      beam->angularCode=BEAM_ANGULAR_NONE;
    if(!strcmp(tmp,"DIVERG"))
      beam->angularCode=BEAM_ANGULAR_DIVERG;
    switch(beam->angularCode){
    case BEAM_ANGULAR_NONE:
      break;
    case BEAM_ANGULAR_DIVERG:
      beam->diverg[0] = (float)param_getDouble("beam.angular.x.sigma",NAN);
      ret |= checkFloat(beam->diverg[0]);
      beam->diverg[1] = (float)param_getDouble("beam.angular.y.sigma",NAN);
      ret |= checkFloat(beam->diverg[1]);
      break;
    default:
      fprintf(stderr,"loadBeam : unknown angular distribution - %s\n",tmp);
      return(PARAM_ERROR);
    }
  }
  return(ret);
}

void setForceReload(int force)
{
  forceReload = force;
}

int loadParticle(char *name)

{
  float mass;
  int index,index2;
  int length;
  void *tmp = NULL;
  unsigned long size;
  particle_db_t *part=NULL;
  int nPart=-1;
  char *code;
  char tmpName[1024];
  char *products;
  char product[24];
  int nProd=0;
  int ret=0;
  float branchSum=0;

  /* Search for the particle in question */
  
  if(!name)
    return(LOAD_PARTICLE_FAIL);
  if(particleDb)
    for(index=0;index<particleDb->nParticles;index++)
      if(!strcmp(name,particleDb->particle[index].name)){
	if(!forceReload)
	  return(index);
	else{
	  nPart = index;
	}
      }

  
  if(!particleDb){
    if(!maxParticles){
      maxParticles = param_getInt("setup.maxParticles",0);
    }
    size = sizeof(particles_db_t) + maxParticles*sizeof(particle_db_t);
    particleDb = malloc(size);
    if(!particleDb){
      fprintf(stderr,"loadParticle : malloc() failed.\n");
      return(LOAD_PARTICLE_FAIL);
    }
    particleDb->nParticles=0;
  }
  if(!(forceReload && (nPart != -1))){
    nPart = particleDb->nParticles;
    strcpy(particleDb->particle[nPart].name,name);
    particleDb->nParticles++;
  }
  part = &particleDb->particle[nPart];
  if(particleDb->nParticles > maxParticles){
    fprintf(stderr,"loadParticle : too many particles, increase value of setup.maxParticles\n");
    return(LOAD_PARTICLE_FAIL);
  }
  
  part->isPhoton = param_getBooleanV(name,".isPhoton");
  ret != checkBoolean(part->isPhoton);
  code = param_getValueV(name,".massCode");
  if(!code)
    return(LOAD_PARTICLE_FAIL);
  cleanString(code);
  part->massCode = 0;
  if(!strcmp(code,"DELTA"))
    part->massCode = MASS_DELTA;
  if(!strcmp(code,"GAUSSIAN"))
    part->massCode = MASS_GAUSSIAN;
  if(!strcmp(code,"BW"))
    part->massCode = MASS_BW;
  if(!strcmp(code,"FLAT"))
    part->massCode = MASS_FLAT;
  if(!strcmp(code,"USER"))
    part->massCode = MASS_USER;
  switch(part->massCode){
  case MASS_DELTA:    
    part->mass = param_getFloatV(name,".mass",NAN);
    ret |= checkFloat(part->mass);
    break;
  case MASS_BW:
  case MASS_GAUSSIAN:
  case MASS_FLAT:
    part->mass = param_getFloatV(name,".mass",NAN);
    ret |= checkFloat(part->mass);
    part->sigma = param_getFloatV(name,".width",NAN);
    ret |= checkFloat(part->sigma);
    break;
  case MASS_USER:
    part->massCode = MASS_USER;
    sprintf(tmpName,"%s.mass",name);
    ret |= readHistogramDef(tmpName,&(part->massHisto));
    break;
  default:
    fprintf(stderr,"loadParticle : unknown mass code %s.\n",code);
  }      

  part->charge = param_getIntV(name,".charge",INT_MAX);
  ret |= (part->charge == INT_MAX)?1:0;

  /* Load the spin crap */

  part->spin = param_getIntV(name,".spin",INT_MAX);
  if(part->spin == INT_MAX)
    ret |= PARAM_ERROR;
  
  /* Get the decay angle info */
  ret |= loadDecayAngles(part);
  part->lifetime = param_getFloatV(name,".lifetime",NAN);
  ret |= checkFloat(part->lifetime);
  part->returnData = param_getBooleanV(name,".returnData");
  ret |= checkBoolean(part->returnData);
  part->nDecays = param_getIntV(name,".nDecays",-1);
  if(part->nDecays < 0)
    return(LOAD_PARTICLE_FAIL);
  

  /* Now read in the decays */

  for(index=0;index<part->nDecays;index++){
    part->decay[index].branchingRatio = param_getFloatVI(name,".decay.%d.branch",index+1,NAN);
    ret |=checkFloat(part->decay[index].branchingRatio);
    products = param_getValueVI(name,".decay.%d.products",index+1);
    
    if(!products)
      return(LOAD_PARTICLE_FAIL);
    
    /* Now parse the products */

    length=0;
    for(index2=0;index2<strlen(products);index2++){
      if(products[index2] != ' '){
	if(products[index2] == ','){	  
	  product[length] = '\0';
	  part->decay[index].product[nProd++] = loadParticle(product);
	  if(part->decay[index].product[nProd-1] == -1)
	    return(LOAD_PARTICLE_FAIL);
	  length=0;
	}
	else{
	  product[length++] = products[index2];
	}
      }      
    }
    if(length){
      product[length] = '\0';
      part->decay[index].product[nProd++] = loadParticle(product);
      if(part->decay[index].product[nProd-1] == -1)
	return(LOAD_PARTICLE_FAIL);
    }
    part->decay[index].nProducts=nProd;
    nProd=0;
  }
  /* Verify the branching ratios */

  if(part->nDecays == 1){
    if(part->decay[0].branchingRatio != 1.0){
      part->decay[0].branchingRatio=1;
      fprintf(stderr,"loadParticle : %s has only one decay mode, BR set to 1.0\n",part->name);
    }
  }
  if(part->nDecays > 1){
    for(index=0;index<part->nDecays;index++)
      branchSum += part->decay[index].branchingRatio;
    if(branchSum != 1){
      fprintf(stderr,"loadParticle : %s BR's sum to %7.4f instead of 1.0\n",part->name,branchSum);
      return(LOAD_PARTICLE_FAIL);
    }
  }

  if(ret)
    return(LOAD_PARTICLE_FAIL);
  setParticleID(part->name,&part->id);
  return(nPart);
}

float param_getFloatVI(char *base,char *name,int i,float def)

{
  char tmp[1024];

  sprintf(tmp,"%s%s",base,name);
  return((float)param_getDoubleI(tmp,i,def));
}		 
float param_getFloatV(char *base,char *name, float def)

{
  char tmp[1024];

  sprintf(tmp,"%s%s",base,name);
  return((float)param_getDouble(tmp,def));
}
		 
int param_getIntVI(char *base,char *name, int i,int  def)

{
  char tmp[1024];

  sprintf(tmp,"%s%s",base,name);
  return((int)param_getIntI(tmp,i,def));
}
		 
int param_getIntV(char *base,char *name, int  def)

{
  char tmp[1024];

  sprintf(tmp,"%s%s",base,name);
  return((int)param_getInt(tmp,def));
}
	 
char *param_getValueV(char *base,char *name)

{
  char tmp[1024];

  sprintf(tmp,"%s%s",base,name);
  return((char *)param_getValue(tmp));
}	
 
char *param_getValueVI(char *base,char *name,int i)

{
  char tmp[1024];

  sprintf(tmp,"%s%s",base,name);
  return((char *)param_getValueI(tmp,i));
}
 
int param_getBoolean(char *name)

{
  char *value;

  value = (char *)param_getValue(name);
  if(!value){
    return(BOOLEAN_FAIL);
  }
  cleanString(value);
  if(toupper(value[0]) == 'T')
    return(PARAM_ERROR);
  if(toupper(value[0]) == 'F')
    return(PARAM_OK);
  return(BOOLEAN_FAIL);
}

int param_getBooleanV(char *base,char *name)

{
  char tmp[1024];
  char *value;

  sprintf(tmp,"%s%s",base,name);
  value = (char *)param_getValue(tmp);
  if(!value){
    return(BOOLEAN_FAIL);
  }
  cleanString(value);  
  if(toupper(value[0]) == 'T')
    return(PARAM_ERROR);
  if(toupper(value[0]) == 'F')
    return(PARAM_OK);
  return(BOOLEAN_FAIL);
}

void cleanString(char *str)

{
  int i;

  for(i=0;i<strlen(str);i++){
    if(str[i] == ' '){
      str[i]='\0';
      return;
    }
  }
}


void printParticle(particle_db_t *part,particles_db_t *particleDb)


{
  int index,index2;

  printf("\n\nName : %s\n",part->name);
  printf("Mass code : 0x%x\n",part->massCode);
  switch(part->massCode){
  case MASS_DELTA:
    printf("Mass : %f\n",part->mass);
    break;
  case MASS_BW:
  case MASS_GAUSSIAN:
  case MASS_FLAT:
    printf("Mass : %f, Width : %f\n",part->mass,part->sigma);
    break;
  case MASS_USER:
    break;
  default:
    break;
  }
  printf("Charge : %d\n",part->charge);
  printf("Lifetime : %f\n",part->lifetime);
  printf("Return : %d\n",part->returnData);
  printf("nDecays : %d\n",part->nDecays);
  for(index=0;index<part->nDecays;index++){
    printf("Decay %2d :\n",index);
    printf("\tBranch : %f\n",part->decay[index].branchingRatio);
    printf("\tProducts : %d\n",part->decay[index].nProducts);    
    for(index2=0;index2<part->decay[index].nProducts;index2++){
      printf("\t\t%d : %s\n",index2,particleDb->particle[part->decay[index].product[index2]].name);
    }
  }
}

int checkBoolean(int value)

{
  if(value==-1)
    return(PARAM_ERROR);
  return(PARAM_OK);
}

int checkFloat(float value1)
{
  return(isnan((double)value1));
}

int checkInt(int value1)
{
  if(value1==0)
    return(PARAM_ERROR);
  return(PARAM_OK);
}


int loadDecayAngles(particle_db_t *part)

{
  char *code;
  int tmpCode;
  char *name = part->name;
  int ret=0;
  char tmpName[1024];

  /* Set the 2 and 3 body codes to 0 */

  part->twoBody.code = 0;
  part->threeBody.code = 0;

  /* Do theta first */

  code = param_getValueV(name,".2Body.theta.code");
  fprintf(stdout,"Debug: .2Body.theta.code for %s = %s.\n",name,code);
  if(!code)
    return(PARAM_ERROR);
  cleanString(code);
  tmpCode=0;
  if(!strcmp(code,"YLM"))
    tmpCode = DECAY_THETA_YLM;
  if(!strcmp(code,"UNIFORM"))
    tmpCode = DECAY_THETA_UNIFORM;
  if(!strcmp(code,"FLAT"))
    tmpCode = DECAY_THETA_FLAT;
  if(!strcmp(code,"USER"))
    tmpCode = DECAY_THETA_USER;
  if(!strcmp(code,"CORRELATED"))
    tmpCode = DECAY_THETA_CORRELATED;
  if(!strcmp(code,"REAL_PHI"))
    tmpCode = DECAY_THETA_REAL_PHI;
  if(!tmpCode){
    fprintf(stderr,"loadParticle: Unknown 2 body theta distribution code : %s\n",code);
  }
  part->twoBody.code |= tmpCode;
  
  /* Handle User histogram */

  if(part->twoBody.code &  DECAY_THETA_CORRELATED){
    char *tmpFile;
    sprintf(tmpName,"%s.2Body.thetaHisto",name);
    tmpFile = param_getValueV(tmpName,".filename");
    if(!tmpFile)
      ret |= PARAM_ERROR;
    else{
      cleanString(tmpFile);
      strcpy(part->twoBody.thetaHisto.filename,tmpFile);
    }
  }

  if(part->twoBody.code &  DECAY_THETA_USER){
    sprintf(tmpName,"%s.2Body.thetaHisto",name);
    ret |= readHistogramDef(tmpName,&part->twoBody.thetaHisto);
  }

  /* Handle "flat" distibution */

  if(part->twoBody.code &  DECAY_THETA_FLAT){
    part->twoBody.theta = param_getFloatV(name,".2Body.theta.value",NAN);
    ret |= checkFloat(part->twoBody.theta);
    part->twoBody.thetaWidth = param_getFloatV(name,".2Body.theta.width",NAN);
    ret |= checkFloat(part->twoBody.thetaWidth);
  }
  
  /* Now do phi */

  
  if(!(part->twoBody.code & DECAY_THETA_CORRELATED)){
    code = param_getValueV(name,".2Body.phi.code");
    if(!code)
      return(PARAM_ERROR);
    cleanString(code);
    tmpCode=0;
    if(!strcmp(code,"UNIFORM"))
      tmpCode = DECAY_PHI_UNIFORM;
    if(!strcmp(code,"FLAT"))
      tmpCode = DECAY_PHI_FLAT;
    if(!strcmp(code,"USER"))
      tmpCode = DECAY_PHI_USER;
    if(!tmpCode){
      fprintf(stderr,"loadParticle: Unknown 2 body phi distribution code : %s\n",code);
    }
    part->twoBody.code |= tmpCode;
    
    /* Handle User histogram */
    
    if(part->twoBody.code &  DECAY_PHI_USER){
      sprintf(tmpName,"%s.2Body.phiHisto",name);
      ret |= readHistogramDef(tmpName,&part->twoBody.phiHisto);
    }
    
    /* Handle "flat" distibution */
    
    if(part->twoBody.code &  DECAY_PHI_FLAT){
      part->twoBody.phi = param_getFloatV(name,".2Body.phi.value",NAN);
      ret |= checkFloat(part->twoBody.phi);
      part->twoBody.phiWidth = param_getFloatV(name,".2Body.phi.width",NAN);
      ret |= checkFloat(part->twoBody.phiWidth);
    }
  }

  /* On to the exciting world of 3-Body */


  /* Do theta first */

  code = param_getValueV(name,".3Body.theta.code");
  if(!code)
    return(PARAM_ERROR);
  cleanString(code);
  tmpCode=0;
  if(!strcmp(code,"DALITZ"))
    tmpCode = DECAY_THETA_DALITZ;
  if(!strcmp(code,"USER"))
    tmpCode = DECAY_THETA_USER;
  if(!tmpCode){
    fprintf(stderr,"loadParticle: Unknown 3 body theta distribution code : %s\n",code);
  }
  part->threeBody.code |= tmpCode;

  /* Histogram is the only option here */

  if(part->threeBody.code &  DECAY_THETA_USER){
    sprintf(tmpName,"%s.3Body.thetaHisto",name);
    ret |= readHistogramDef(tmpName,&part->threeBody.thetaHisto);
  }

  /* Last but not least, 3-body phi */

  code = param_getValueV(name,".3Body.phi.code");
  if(!code)
    return(PARAM_ERROR);
  cleanString(code);
  tmpCode=0;
  if(!strcmp(code,"UNIFORM"))
    tmpCode = DECAY_PHI_UNIFORM;
  if(!strcmp(code,"FLAT"))
    tmpCode = DECAY_PHI_FLAT;
  if(!strcmp(code,"USER"))
    tmpCode = DECAY_PHI_USER;
  if(!tmpCode){
    fprintf(stderr,"loadParticle: Unknown 3 body phi distribution code : %s\n",code);
  }
  part->threeBody.code |= tmpCode;

   /* Handle User histogram */

  if(part->threeBody.code &  DECAY_PHI_USER){
    sprintf(tmpName,"%s.3Body.phiHisto",name);
    ret |= readHistogramDef(tmpName,&part->threeBody.phiHisto);
  }

  /* Handle "flat" distibution */

  if(part->threeBody.code &  DECAY_PHI_FLAT){
    part->threeBody.phi = param_getFloatV(name,".3Body.phi.value",NAN);
    ret |= checkFloat(part->threeBody.phi);
    part->threeBody.phiWidth = param_getFloatV(name,".3Body..phi.width",NAN);
    ret |= checkFloat(part->threeBody.phiWidth);
  }
  return(ret);
}


int readHistogramDef(char *root,histogramDef_t *histo)
{
  char *tmpName;
  int ret=0;

  tmpName = param_getValueV(root,".filename");
  if(!tmpName)
    return(PARAM_ERROR);
  cleanString(tmpName);
  strcpy(histo->filename,tmpName);

  histo->id = param_getIntV(root,".id",0);
  ret |= checkInt(histo->id);
  return(ret);
}

void setParticleID(char *name,Particle_t *id)
{
  if(!strcmp(name,"pi0")){
    *id = Pi0;
    return;
  }
  if(!strcmp(name,"pi+")){
    *id = PiPlus;
    return;
  }
  if(!strcmp(name,"pi-")){
    *id = PiMinus;
    return;
  }
  if(!strcmp(name,"gamma")){
    *id = Gamma;
    return;
  }
  if(!strcmp(name,"proton")){
    *id = Proton;
    return;
  }
  if(!strcmp(name,"neutron")){
    *id = Neutron;
    return;
  }
  if(!strcmp(name,"eta")){
    *id = Eta;
    return;
  }
  if(!strcmp(name,"omega")){
    *id = omega;
    return;
  }
  if(!strcmp(name,"rho0")){
    *id = Rho0;
    return;
  }
  if(!strcmp(name,"rho-")){
    *id = RhoMinus;
    return;
  }
  if(!strcmp(name,"rho+")){
    *id = RhoPlus;
    return;
  }
  if(!strcmp(name,"etaprime")){
    *id = EtaPrime;
    return;
  }
  *id=Unknown;
  return;
}

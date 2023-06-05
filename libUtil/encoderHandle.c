/*
  $Log#
  */


#include <stdio.h>
#include <stdlib.h>

#include <encoderHandle.h>
#include <map_manager.h>

static char rcsid[] = "$Id";

static char map[1024];
static int initEncoderHandle(void);

static int initEncoderHandle(void)
{ 
  static char *dir=NULL;
  
  if(dir != NULL)
    return(0);
  if((dir = getenv("RADPHI_CONFIG")) == NULL){
    fprintf(stderr,"encodersToMap: RADPHI_CONFIG not defined\n");
    return(1);
  }
  sprintf(map,"%s/maps/encoders.map",dir);
  return(0);
}

int normalizeEncoder(int value,int gain)
{
  int index;
  int realGain=1;
  int realValue;
  
  for(index=0;index<gain;index++)
    realGain=realGain*2;
  realValue = value/realGain;
  return(realValue);
}

int getEncoderData(encoderData_t *encoders, int time)
{
  int actualTime;

  if(initEncoderHandle())
    return(1);
  map_get_int(map,"lgd","x",1,&encoders->lgd.x,time,&actualTime);
  map_get_int(map,"lgd","y",1,&encoders->lgd.y,time,&actualTime);
  map_get_int(map,"lgd","xGain",1,&encoders->lgdGain.x,time,&actualTime);
  map_get_int(map,"lgd","yGain",1,&encoders->lgdGain.y,time,&actualTime);
  map_get_int(map,"rpd","x",1,&encoders->rpd.x,time,&actualTime);
  map_get_int(map,"rpd","y",1,&encoders->rpd.y,time,&actualTime);
  map_get_int(map,"rpd","z",1,&encoders->rpd.z,time,&actualTime);
  map_get_int(map,"rpd","xGain",1,&encoders->rpdGain.x,time,&actualTime);
  map_get_int(map,"rpd","yGain",1,&encoders->rpdGain.y,time,&actualTime);
  map_get_int(map,"rpd","zGain",1,&encoders->rpdGain.z,time,&actualTime);
  map_get_int(map,"norm","value",1,&encoders->norm,time,&actualTime);
  map_get_int(map,"norm","gain",1,&encoders->normGain,time,&actualTime);
  return(0);
  
}

int encodersToMap(encoderData_t *encoders, int time)
   
{
  encoderData_t lastData;
  if(initEncoderHandle())
    return(1);
  getEncoderData(&lastData,0);
  
  /* Do LGD */
  
  if(encoders->lgdGain.x != lastData.lgdGain.x){
    map_put_int(map,"lgd","xGain",1,&encoders->lgdGain.x,time);
  }
  if(encoders->lgdGain.y != lastData.lgdGain.y){
    map_put_int(map,"lgd","yGain",1,&encoders->lgdGain.y,time);
  }
  if(abs(encoders->lgd.x - lastData.lgd.x) > ENCODER_THRESH){
    map_put_int(map,"lgd","x",1,&encoders->lgd.x,time);
  }
  if(abs(encoders->lgd.y - lastData.lgd.y) > ENCODER_THRESH){
    map_put_int(map,"lgd","y",1,&encoders->lgd.y,time);
  }


  /* RPD */

  if(encoders->rpdGain.x != lastData.rpdGain.x){
    map_put_int(map,"rpd","xGain",1,&encoders->rpdGain.x,time);
  }
  if(encoders->rpdGain.y != lastData.rpdGain.y){
    map_put_int(map,"rpd","yGain",1,&encoders->rpdGain.y,time);
  }
  if(encoders->rpdGain.z != lastData.rpdGain.z){
    map_put_int(map,"rpd","zGain",1,&encoders->rpdGain.z,time);
  }
  if(abs(encoders->rpd.x - lastData.rpd.x) > ENCODER_THRESH){
    map_put_int(map,"rpd","x",1,&encoders->rpd.x,time);
  }
  if(abs(encoders->rpd.y - lastData.rpd.y) > ENCODER_THRESH){
    map_put_int(map,"rpd","y",1,&encoders->rpd.y,time);
  }
  if(abs(encoders->rpd.z - lastData.rpd.z) > ENCODER_THRESH){
    map_put_int(map,"rpd","z",1,&encoders->rpd.z,time);
  }

  /* Norm */


  if(encoders->normGain != lastData.normGain){
    map_put_int(map,"norm","gain",1,&encoders->normGain,time);
  }
  if(abs(encoders->norm - lastData.norm) > ENCODER_THRESH){
    map_put_int(map,"norm","value",1,&encoders->norm,time);
  }
  return(0);
}
float encoderToLength(encoderData_t *encoders,int which)
{
  float inchesPerOhm=0;
  float voltage;
  float sourceV;
  float length;

  sourceV = encoders->norm*VOLTS_PER_COUNT;
  switch(which){
  case ENC_LGD_X:
    inchesPerOhm = IPO_LGX_X;
    voltage = encoders->lgd.x*VOLTS_PER_COUNT;
    break;
  case ENC_LGD_Y:
    inchesPerOhm = IPO_LGX_Y;
    voltage = encoders->lgd.y*VOLTS_PER_COUNT;
    break;
  case ENC_RPD_X:
    inchesPerOhm = IPO_RPD_X;
    voltage = encoders->rpd.x*VOLTS_PER_COUNT;
    break;
  case ENC_RPD_Y:
    inchesPerOhm = IPO_RPD_Y;
    voltage = encoders->rpd.y*VOLTS_PER_COUNT;
    break;
  case ENC_RPD_Z:
    inchesPerOhm = IPO_RPD_Z;
    voltage = encoders->rpd.z*VOLTS_PER_COUNT;
    break;
  }
  switch(which){
  case ENC_RPD_Z:
    length= 20.0*voltage/sourceV;
    break;
  default:
    length = inchesPerOhm*500.0*((1.0/(1.0-(voltage/sourceV)))-1.0);
    break;
  }
  /* Convert to cm */

  return(length*2.54);
}

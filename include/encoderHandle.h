/*
  $Log: encoderHandle.h,v $
  Revision 1.1  1997/04/08 21:14:42  radphi
  Initial revision

  */
#ifndef ENCODER_HANDLE_H_INCLUDED
#define ENCODER_HANDLE_H_INCLUDED

typedef struct{
  int x;
  int y;
  int z;
}ivector3_t;

typedef struct{
  ivector3_t lgd;
  ivector3_t rpd;
  ivector3_t lgdGain;
  ivector3_t rpdGain;
  int norm;
  int normGain;
}encoderData_t;

#define VOLTS_PER_COUNT 0.000305175

#define IPO_LGX_X 0.05
#define IPO_LGX_Y 0.05
#define IPO_RPD_X 0.05
#define IPO_RPD_Y 0.02
#define IPO_RPD_Z 0.02

#define ENC_LGD_X 1
#define ENC_LGD_Y 2
#define ENC_RPD_X 3
#define ENC_RPD_Y 4
#define ENC_RPD_Z 5

int normalizeEncoder(int value,int gain);
int getEncoderData(encoderData_t *encoders, int time);
int encodersToMap(encoderData_t *encoders, int time);
float encoderToLength(encoderData_t *encoders,int which);

#define ENCODER_THRESH 5

#endif

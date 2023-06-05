/* $Log: rpdXport.h,v $
/* Revision 1.1  1997/03/03 18:15:09  lfcrob
/* Initial revision
/*
 */

#ifndef RPDXPORT_H_INCLUDED
#define RPDXPORT_H_INCLUDED

#define RPDXP_OK 0
#define RPDXP_ERROR 1
#define RPDXP_CAMAC_ERROR 2

#define RPDXP_MOTOR_0 0
#define RPDXP_MOTOR_1 1

#define RPDXP_CW 0x1
#define RPDXP_CCW 0x2

#define RPDXP_LEFT RPDXP_CW
#define RPDXP_RIGHT RPDXP_CCW
#define RPDXP_UP RPDXP_CW
#define RPDXP_DOWN RPDXP_CCW


#define RPDXP_HORIZONTAL RPDXP_MOTOR_0
#define RPDXP_VERTICAL RPDXP_MOTOR_1

#define RPDXP_MAX_STEPS 0x7fff
#define RPDXP_CCW_MASK 0x8000

#define RPDXP_CW_LIMIT 0x1
#define RPDXP_CCW_LIMIT 0x2

#define RPDXP_HOR_CW_LIMIT 0x1
#define RPDXP_HOR_CCW_LIMIT 0x2
#define RPDXP_HOR_DONE 0x4
#define RPDXP_VER_CW_LIMIT 0x8
#define RPDXP_VER_CCW_LIMIT 0x10
#define RPDXP_VER_DONE 0x20
#define RPDXP_NO_POWER 0x40
#define RPDXP_POWER_STUCK 0x80

int rpdXP_init(void);
int rpdXP_status(int *status);
int rpdXP_move(int motor,int direction,int nSteps);
#endif

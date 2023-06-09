/*
  $Log: rpdXP.h,v $
  Revision 1.1  1997/04/08 21:09:24  radphi
  Initial revision

  */

#include <camacControl.h>


#ifndef RPDXP_H_INCLUDED
#define RPD_XP_H_INCLUDED

#define RPDXP_VERTICAL_ID 0
#define RPDXP_HORIZONTAL_ID 1
#define RPDXP_Z_ID 2

#define RPDXP_VERTICAL 1
#define RPDXP_HORIZONTAL 0
#define RPDXP_Z 0

#define RPDXP_PZ JOERGER_SMCR_CCW
#define RPDXP_MZ JOERGER_SMCR_CW
#define RPDXP_LEFT JOERGER_SMCR_CCW
#define RPDXP_RIGHT JOERGER_SMCR_CW
#define RPDXP_UP JOERGER_SMCR_CW
#define RPDXP_DOWN JOERGER_SMCR_CCW


#define RPDXP_Z_CW_LIMIT JOERGER_SMCR_M0_CW_LIMIT
#define RPDXP_Z_CCW_LIMIT JOERGER_SMCR_M0_CCW_LIMIT
#define RPDXP_Z_DONE JOERGER_SMCR_M0_DONE

#define RPDXP_HOR_CW_LIMIT JOERGER_SMCR_M0_CW_LIMIT
#define RPDXP_HOR_CCW_LIMIT JOERGER_SMCR_M0_CCW_LIMIT
#define RPDXP_HOR_DONE JOERGER_SMCR_M0_DONE
#define RPDXP_VER_CW_LIMIT JOERGER_SMCR_M1_CW_LIMIT
#define RPDXP_VER_CCW_LIMIT JOERGER_SMCR_M1_CCW_LIMIT
#define RPDXP_VER_DONE JOERGER_SMCR_M1_DONE
#define RPDXP_NO_POWER JOERGER_SMCR_NO_POWER
#define RPDXP_POWER_STUCK JOERGER_SMCR_POWER_STUCK

#define RPDXP_UP_LIMIT RPDXP_VER_CW_LIMIT
#define RPDXP_DOWN_LIMIT RPDXP_VER_CCW_LIMIT
#define RPDXP_RIGHT_LIMIT RPDXP_HOR_CW_LIMIT
#define RPDXP_LEFT_LIMIT RPDXP_HOR_CCW_LIMIT
#define RPDXP_PZ_LIMIT RPDXP_Z_CCW_LIMIT
#define RPDXP_MZ_LIMIT RPDXP_Z_CW_LIMIT

#define STEPS_PER_CM_HOR 203.5
/*#define STEPS_PER_CM_HOR 193.55*/
#define STEPS_PER_CM_VER 30715.2231
#define STEPS_PER_CM_Z 1.0

int rpdXP_getXYModule(int *myCrate, int *mySlot);
int rpdXP_getZModule(int *myCrate, int *mySlot);

#endif

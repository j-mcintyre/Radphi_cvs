/* $Log: lgdXP.h,v $
/* Revision 1.2  1997/04/15 16:14:51  radphi
/* Converted to use libcamacControl
/* Added cm/sec defines
/*
 * Revision 1.1  1997/02/10 01:35:52  lfcrob
 * Initial revision
 *
 */

#ifndef LGDXP_H_INCLUDED
#define  LGDXP_H_INCLUDED 

#define LGDXP_1_CW_ON      0x1
#define LGDXP_1_CW_LIMIT   0x2
#define LGDXP_1_CCW_ON     0x4
#define LGDXP_1_CCW_LIMIT  0x8
#define LGDXP_2_CW_ON      0x10
#define LGDXP_2_CW_LIMIT   0x20
#define LGDXP_2_CCW_ON     0x40
#define LGDXP_2_CCW_LIMIT  0x80

#define LGDXP_MOVE_UP 0x1
#define LGDXP_MOVE_DOWN 0x2
#define LGDXP_MOVE_BLEFT 0x4
#define LGDXP_MOVE_BRIGHT 0x8

#define LGDXP_MOVE_UP_ON  LGDXP_2_CW_ON
#define LGDXP_MOVE_UP_LIMIT  LGDXP_2_CW_LIMIT
#define LGDXP_MOVE_DOWN_ON  LGDXP_2_CCW_ON
#define LGDXP_MOVE_DOWN_LIMIT  LGDXP_2_CCW_LIMIT
#define LGDXP_MOVE_BLEFT_ON  LGDXP_1_CW_ON
#define LGDXP_MOVE_BLEFT_LIMIT  LGDXP_1_CW_LIMIT
#define LGDXP_MOVE_BRIGHT_ON  LGDXP_1_CCW_ON
#define LGDXP_MOVE_BRIGHT_LIMIT  LGDXP_1_CCW_LIMIT

#define LGDXP_CM_PER_SEC_HOR 0.269875
#define LGDXP_SEC_PER_CM_HOR 3.71
/*#define LGDXP_CM_PER_SEC_VER 0.245237*/
/*#define LGDXP_CM_PER_SEC_VER 0.091963*/
#define LGDXP_SEC_PER_CM_UP   10.1981
#define LGDXP_SEC_PER_CM_DOWN 10.2161

#endif

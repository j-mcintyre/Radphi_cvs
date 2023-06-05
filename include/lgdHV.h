

#define LGDHV_OK 0
#define LGDHV_READY 0
#define LGDHV_ERROR 1
#define LGDHV_BUSY 2
#define LGDHV_CAMAC_ERROR 3


int lgdHV_init(void);
int lgdHV_ready(void);
int lgdHV_setChannel(int row, int col, int voltage);
int lgdHV_setRow(int row, int voltage);
int lgdHV_setCol(int col, int voltage);
int lgdHV_setAll(int voltage);
int lgdHV_pulseChannel(int row, int col);
int lgdHV_pulseRow(int row);
int lgdHV_pulseCol(int col);
int lgdHV_pulseAll(void);


#ifndef LGDCONTROLDB_H_INCLUDED
#define LGDCONTROLDB_H_INCLUDED

#define NUM_CHANNELS 784

#define LGDDB_TRUE 1
#define LGDDB_FALSE 0
#define LGDDB_OK 0
#define LGDDB_MAP_ERROR -1
#define LGDDB_BADCHANNEL -2

#define LGDDB_INSTRUMENTED 1
#define LGDDB_NOT_INSTRUMENTED 0
#define LGDDB_UNKNOWN_VOLTAGE -1
#define LGDDB_NO_CABLE -1
#define LGDDB_NOT_FOUND -1
#define LGDDB_BASE_SETPOINT 1000
#define LGDDB_MAX_SETPOINT 2022


typedef struct{
  int channel;
  int cable;
  int instrumented;
  int voltage;
  int setpoint;    
} lgdChannel_t;

typedef struct{
  int nChannels;
  lgdChannel_t channel[1];
} lgdChannels_t;

int lgdDB_load(int time,char *userDB,char *dumpFile);
int lgdDB_putCable(int channel,int newCable);
int lgdDB_putInstrumented(int channel,int newStatus);
int lgdDB_putVoltage(int channel,int newVoltage);
int lgdDB_putSetpoint(int channel,int newSetpoint);
int lgdDB_changed(void);
int lgdDB_update(int newTime);
int lgdDB_getCable(int channel);
int lgdDB_getInstrumented(int channel);
int lgdDB_getVoltage(int channel);
int lgdDB_getSetpoint(int channel);
int lgdDB_getChannel(int channel,lgdChannel_t *item);
int lgdDB_getChannelForCable(int cable);

int lgdDB_readDump(char *filename);
int lgdDB_writeDump(char *filename);
#endif

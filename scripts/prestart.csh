#! /bin/csh

#
#
# To be run at DAQ prestart
#
#

source /home/radphi/setup

/usr/bin/date >>&! /radphi/logfiles/triggerLoad.log
master_trigger_control.perl `get_run_number` >>&! /radphi/logfiles/triggerLoad.log
/usr/local/bin/perl $RADPHI_HOME/scripts/killTriggerd.perl
~radphi/swdev/dev/triggerControl/triggerd -q >>&! /radphi/logfiles/triggerd.log &
#$RADPHI_BIN/runLog  >>&! /radphi/logfiles/runs.log 

# end file


\#! /bin/csh


#
#
#
#  To be run at DAQ end
#
#
#

source /home/radphi/setup

$RADPHI_BIN/beep
xclock -geom 1024x666 -bg red -hl yellow -update 1 -chime &
/usr/local/bin/perl $RADPHI_HOME/scripts/killTriggerd.perl
/usr/local/bin/perl $RADPHI_HOME/scripts/kill_generic.perl prnt_trg
#end file

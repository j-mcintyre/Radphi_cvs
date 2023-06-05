#! /bin/csh


#
# Start the Dispatcher
#

source /home/radphi/setup

cd $RADPHI_LOGDIR
touch ONLINE_RUNNING
$RADPHI_BIN/log -nDispatcher -e $RADPHI_BIN/Dispatcher -H 10357 -socketBuffer 150000 -sizeWRITEqueue 100 -sizeRAWqueue 100 
rm ONLINE_RUNNING


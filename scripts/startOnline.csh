#! /bin/csh -f
#
#  Starts the radphi online Dispatcher system
#    miserably hacked into by D. Armstrong
#   this should only be run from urs2 !!

source /home/radphi/setup

echo "Starting Dispatcher"

$RADPHI_BIN/Dispatcher -H 10357 -socketBuffer 150000 -sizeWRITEqueue 100 -sizeRAWqueue 100 >& /dev/null &

sleep 3

echo "Starting OLA"

$RADPHI_BIN/OLA -d urs2.cebaf.gov:10357 >& /dev/null &

echo "Starting scalerd"

$RADPHI_BIN/scalerd >& /dev/null & 

echo "Starting onlinehisto"

$RADPHI_BIN/onlinehisto >& /dev/null & 

echo "Starting feedFromDD"

rsh urs1.jlab.org -l radphi "source /home/radphi/setup; source /home/radphi/dosetupcoda; /home/radphi/bin.solaris/feedFromDD -Durs2.jlab.org:10357 >& /dev/null "  >& /dev/null &
 




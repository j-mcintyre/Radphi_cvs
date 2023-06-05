#!/usr/bin/csh
# script to execute dat_put
#
# elton 7/11/2000
#
echo "Run Number is "$1
rm tagdata 
set filename = "~mammei/radphi/paw/tagger_mean$1.dat"
echo "File Name is "$filename 
ln -s $filename tagdata
data_put $1
#
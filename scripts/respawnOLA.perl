#! /usr/local/bin/perl


$file = "ONLINE_RUNNING";
while(-e $file){
    $cmd = "$ENV{'RADPHI_BIN'}/log -nOLA -e $ENV{'RADPHI_BIN'}/OLA -d $ENV{'RADPHI_DISPATCHER'}";
    system($cmd);
}

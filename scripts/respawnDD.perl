#! /usr/local/bin/perl

$file = "ONLINE_RUNNING";
while(-e $file){
    $cmd = "$ENV{'RADPHI_BIN'}/log -nfeedFromDD -e $ENV{'RADPHI_BIN'}/feedFromDD -D$ENV{'RADPHI_DISPATCHER'}";
#    print "$cmd\n";
    system($cmd);
}

#! /usr/local/bin/perl


$file = "ONLINE_RUNNING";
while(-e $file){
    $cmd = "$ENV{'RADPHI_BIN'}/log -nscalerd -e $ENV{'RADPHI_BIN'}/scalerd";
    system($cmd);
}

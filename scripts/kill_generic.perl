#! /usr/local/bin/perl

if(!@ARGV[0]){
    die("kill_generic.perl invoked with no arguments!\n");
}

$cmd = "ps -u radphi | grep @ARGV[0]";

open(PS,"$cmd |");
while (<PS>){
    chop;
    $PID = (split())[0];
    $cmd = "kill -9 $PID";
    system($cmd);
}

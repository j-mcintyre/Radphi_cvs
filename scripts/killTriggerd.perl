#! /usr/local/bin/perl



$cmd = "ps -u radphi | grep triggerd";

open(PS,"$cmd |");
while (<PS>){
    chop;
    $PID = (split())[0];
    $cmd = "kill -9 $PID";
    system($cmd);
}

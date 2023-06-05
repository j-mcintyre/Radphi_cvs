#!/apps/bin/perl
$cmd="trig_setup_1.perl @ARGV[0]";
system($cmd);

$cmd="/home/radphi/swdev/dev/triggerControl/triggerSetup2 /tmp/default.camac.xrdb /tmp/camac.xrdb >/dev/null";
$customtrigger=system($cmd);
if($customtrigger){
    $cmd="/usr/bin/rm /tmp/camac.xrdb";
    system($cmd);
    $cmd="xterm -e /home/radphi/swdev/dev/triggerControl/triggerSetup /tmp/default.camac.xrdb /tmp/camac.xrdb";
    system($cmd);
}
$modfile="/tmp/camac.xrdb";
if(open(MODFILE,$modfile)){
    $cmd="xrdbToMap -r /tmp/camac.xrdb @ARGV[0]";
    system($cmd);
    $cmd="triggerLoad /tmp/camac.xrdb";
    system($cmd);
    $cmd="xterm -e /home/radphi/swdev/dev/triggerControl/prnt_trg /tmp/camac.xrdb &";
    system($cmd);
}
else{
    print "perl: loading default trigger\n";
    $cmd="triggerLoad /tmp/default.camac.xrdb";
    system($cmd);
    $cmd="xterm -e /home/radphi/swdev/dev/triggerControl/prnt_trg /tmp/default.camac.xrdb &";
    system($cmd);
}



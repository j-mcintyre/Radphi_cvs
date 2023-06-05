#! /usr/local/bin/perl


require "getopts.pl";

&Getopts('sfhv');

if($opt_h){
    &Usage();
    exit(0);
}

if($#ARGV < 1){
    &Usage();
    exit(0);
}

$dir = $ARGV[$#ARGV];
$file = $ARGV[$#ARGV-1];
$installFile = $dir."/".$file;
if($opt_v){
    print "Install $file in $dir\n";
}
if(!-e $dir){
    print STDERR "ERROR : $dir does not exists!\n";
    exit(1);
}
if(!-e $file){
    print STDERR "ERROR : $file does not exists!\n";
    exit(1);
}

if(-e $installFile){
    $ctime1 = (stat($installFile))[10];
    $ctime2 = (stat($file))[10];
    if(($ctime2 < $ctime1) && !($opt_f)){
	if($opt_v){
	    print "$installFile is newer, no action taken.\n";
	}
	    exit(0);
    }
}
if(($opt_s) && (-e $installFile)){
    $cmd = "mv $installFile $installFile.O";
    !system($cmd) || die "Failed to make backup of $installFile\n";
}

# Do this install

$cmd = "cp $file $installFile";
!system($cmd) || die "Failed to install $file in $installFile\n";


# END


sub Usage{
    print STDERR "Usage:\n";
    print STDERR "  install.perl <Options> file dir\n";
    print STDERR "Options are:\n";
    print STDERR "  -v         verbose output\n";
    print STDERR "  -s         save old files in <file>.O\n";
    print STDERR "  -f         force installation over newer file.\n";
    print STDERR "  -h         print this message.\n";

}

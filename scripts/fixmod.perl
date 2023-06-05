#! /usr1/local/bin/perl

$home = "/group/radphi/";
@systems =("source","config","bin.*","lib.*","include","makefiles");
$dirmode = 0775;

($UID,$GID) = (getpwnam("radphi"))[2,3];

foreach $system (@systems){
    $rootdir = $home.$system;
    if($rootdir =~ /\*/){
	$rootdir =~ s/\*//;
	opendir(ROOT,$home);
	@allfiles = readdir(ROOT);
	closedir(ROOT);    
	foreach $file (@allfiles){
	    if($file =~ /$system/){
		doDir($rootdir.$file);
	    }
	}
    }
    else{
	doDir($rootdir);
    }
}
sub doDir{
    local($rootdir) = $_[0];
    opendir(DIR,$rootdir);
    local(@allfiles) = readdir(DIR);
    closedir(DIR);
    foreach $file (@allfiles){

# We ignore . and .. entries

	if(($file ne ".") && ($file ne "..")){
	    local($fullfile) = $rootdir."/".$file;
	    if(-d $rootdir."/".$file){
		chown($UID,$GID,$fullfile);
		chmod($dirmode,$fullfile);
		system("chmod g+s $fullfile");
		doDir($fullfile);
	    }
	    else{ 
		if (-f  $fullfile){
#		    chown($UID,$GID,$fullfile);
		}
	    }
	}
    }
}

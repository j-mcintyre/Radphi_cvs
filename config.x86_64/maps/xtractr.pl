#!/usr/bin/perl

our $extracting = 0;

while (<>) {
   my $line = $_;
   if ($line =~ /^ *Item:  *([^\.]*)\./) {
      $item = $1;
   }
   elsif ($line =~ /^ *Time:  *([0-9]*),/) {
      $run = $1;
      print "extracting ${item}_$run.xrdb\n";
      if ($extracting) {
         system("../../libParam/xrdbFromMap $item $run >${item}_$run.xrdb") || die;
      }
      else {
         unlink "$item.xrdb";
         symlink "${item}_$run.xrdb", "$item.xrdb" || die;
         system("../../libParam/xrdbToMap $item.xrdb $run") || die;
     }
   }
}

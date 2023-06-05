#!/usr/bin/perl

use FileHandle;

while (<>) {
   my $file = $_;
   chop $file;
   my $in = new FileHandle($file);
   my $out = new FileHandle(">$file.out");
   my @out = ();
   while (<$in>) {
      my $line = $_;
      if ($line =~ /^$/) {
         push @out, $line;
      }
      else {
         print $out @out if (@out);
         print $out $line;
         @out = ();
      }
   }
   print "truncated ",scalar @out," lines from file $file\n";
   $out->close();
   $in->close();
}

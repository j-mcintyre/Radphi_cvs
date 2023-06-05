#!/usr/bin/perl
#
# lgd_status.pl :  reads run-dependent channel status information
#		   from text files and formats the output in the
#		   ascii map format.  To install the results into
#		   the map, use maptoascii to dump the current map
#		   contents into a text file, then use a text editor
#		   to insert the output from lgd_status.pl at the
#		   appropriate place in the map.  The final step is
#		   to convert the map back to binary using mapfrascii.
#
# Richard Jones, University of Connecticut
# September 18, 2002
#-------------------

open(FD,"lgd_hole_blocks.txt") || die "cannot open lgd_hole_blocks.txt";
fill_hash(1);
close(FD);

open(FD,"<lgd_cold_blocks.txt") || die "cannot open lgd_cold_blocks.txt";
fill_hash(2);
close(FD);

open(FD,"<lgd_hot_blocks.txt") || die "cannot open lgd_hot_blocks.txt";
fill_hash(4);
close(FD);

if ($ARGV[0]) {
   my $chan = $ARGV[0];
   my $code = 0;
   foreach $run (sort keys %mask) {
      if (${$mask{$run}}[$chan] != $code) {
         $code = ${$mask{$run}}[$chan];
         print "channel $chan changed to $code at run $run\n";
      }
   }
   exit;
}

foreach $run (sort {$b <=> $a} keys %mask) {
   printf("C %d %d\n",$run,784);
   for ($chan=0; $chan<784;) {
      for ($i=0; $chan<784 && $i<20; $i++) {
         $code = ${mask{$run}}[$chan++];
         $code = 4 if ($code == 6);
         printf(" %d",$code);
      }
      print "\n";
   }
}

sub fill_hash
{
   my $code = "@_";
   while ($line = <FD>) {
      chop $line;
      @line = split(' ',$line);
      $chan = shift @line;
      while ($range = shift @line) {
        set_flag($range,$chan,$code);
      }
   }
}

sub set_flag
{
   my ($range,$chan,$code) = @_;
   if ($range =~ /-/) {
      ($run0, $run1) = split(/-/,$range);
   }
   else {
      $run0 = $range;
      $run1 = $range;
   }
   $run1++;
   $mask{$run0} = allocate($run0);
   foreach $run (sort keys %mask) {
      next if ($run < $run0);
      last if ($run > $run1);
      ${$mask{$run}}[$chan] |= $code;
   }
   $mask{$run1} = allocate($run1);
   ${$mask{$run1}}[$chan] &= ~($code);
}

sub allocate
{
   my $run = "@_";
   if (! defined $mask{0}) {
      for ($i=0;$i<784;$i++) {$default_mask[$i] = 0;}
      $mask{0} = \@default_mask;
   }
   foreach $trun (sort keys %mask) {
      last if ($trun > $run);
      $ref = $trun;
      last if ($trun == $run);
   }
   return [ @{$mask{$ref}} ];
}

#!/usr/bin/perl
#
# rundb.cgi - Searchable index to the Radphi run logbook database
#	      for use as a CGI script on a web server
#
# Richard Jones
# July 18, 1999 -RTJ
#	v1.0  -  original version
#
# Oct. 12, 2000 -RTJ
#	v1.1  -  switched to use Berkeley DB database
#
# Oct. 23, 2000 -RTJ
#	v1.2  -	 added hooks for selecting staged disks and for "not"
#                patterns to select on trigger and comment fields.
#
     
use CGI;
use DB_File;
use Fcntl;
use POSIX;

@path_info = split('/',$ENV{'PATH_INFO'});

$DBfilename = "/home/radphi/runbook.db";
tie (%dbase, 'DB_File', $DBfilename, &O_RDONLY, 0644);

my (%month,%year);
%month = qw(Jan 0 Feb 1 Mar 2 Apr 3 May 4 Jun 5
            Jul 6 Aug 7 Sep 8 Oct 9 Nov 10 Dec 11);
%year  = qw(1996 96 1997 97 1998 98 1999 99 2000 100
            2001 101 2002 102 2003 103 2004 104 2005 105);

$query = new CGI;
$query->autoEscape(0);

print $query->header('text/html');

print
   $query->start_html(-title=>'Radphi Run Index',
                      -author=>'jonesrt@uconnvm.uconn.edu',
                      -BGCOLOR=>'80e0f0'),"\n",
   "<H1 ALIGN=CENTER> Radphi Data Logbook</H1>\n";

# Check if a form is being submitted
$action = $query->param('Submit');
if ("$action" eq 'Search') 
{
   process_submitted_form($query);
}
elsif ("$path_info[1]")
{
   print_detailed_info(@path_info);
   exit;
}

print
   $query->startform,
   "<TABLE WIDTH=600 ALIGN=CENTER>\n",
   "<TR><TD><NOBR>",
   $query->submit(-name=>'Submit',
                  -value=>'Search'),
   " for up to ",
   $query->textfield(-name=>"runCount",
                     -default=>"10",
                     -size=>4,
                     -maxlength=>4),
   " runs starting at run ",
   $query->textfield(-name=>"runStart",
                     -default=>"",
                     -size=>5,
                     -maxlength=>5),
   " with specified beam energy ",
   $query->textfield(-name=>"energySpec",
                     -default=>"",
                     -size=>4,
                     -maxlength=>4),
   "GeV",
   "</NOBR></TD></TR>\n",
   "<TR><TD><NOBR>",
   "between ",
   $query->popup_menu(-name=>"monthStart",
                      -values=>[0,1,2,3,4,5,6,7,8,9,10,11],
                      -labels=>{0=>"Jan", 1=>"Feb", 2=>"Mar",
                                3=>"Apr", 4=>"May", 5=>"Jun",
                                6=>"Jul", 7=>"Aug", 8=>"Sep",
                                9=>"Oct",10=>"Nov",11=>"Dec"},
                      -default=>"0"),
   " ",
   $query->textfield(-name=>"dayStart",
                     -default=>" 1",
                     -size=>2,
                     -maxlength=>3),
   ", ",
   $query->popup_menu(-name=>"yearStart",
                      -values=>[96,97,98,99,100,101,102,103,104,105],
                      -labels=>{96=>"1996", 97=>"1997", 98=>"1998",
                                99=>"1999",100=>"2000",101=>"2001",
                               102=>"2002",103=>"2003",104=>"2004",
                               105=>"2005"},
                      -default=>"96"),
   " at ",
   $query->textfield(-name=>"timeStart",
                     -default=>"00:00",
                     -size=>5,
                     -maxlength=>5),
   " and ",
   $query->popup_menu(-name=>"monthStop",
                      -values=>[0,1,2,3,4,5,6,7,8,9,10,11],
                      -labels=>{0=>"Jan", 1=>"Feb", 2=>"Mar",
                                3=>"Apr", 4=>"May", 5=>"Jun",
                                6=>"Jul", 7=>"Aug", 8=>"Sep",
                                9=>"Oct",10=>"Nov",11=>"Dec"},
                      -default=>"11"),
   " ",
   $query->textfield(-name=>"dayStop",
                     -default=>"31",
                     -size=>2,
                     -maxlength=>3),
   ", ",
   $query->popup_menu(-name=>"yearStop",
                      -values=>[96,97,98,99,100,101,102,103,104,105],
                      -labels=>{96=>"1996", 97=>"1997", 98=>"1998",
                                99=>"1999",100=>"2000",101=>"2001",
                               102=>"2002",103=>"2003",104=>"2004",
                               105=>"2005"},
                      -default=>"105"),
   " at ",
   $query->textfield(-name=>"timeStop",
                     -default=>"23:59",
                     -size=>5,
                     -maxlength=>5),
   "</NOBR></TD></TR>\n",
   "<TR><TD><NOBR>",
   "archived on ",
   $query->textfield(-name=>"archivedOn",
                     -default=>"",
                     -size=>5,
                     -maxlength=>5),
   " but not on ",
   $query->textfield(-name=>"archivedNotOn",
                     -default=>"",
                     -size=>5,
                     -maxlength=>5),
   " and on disk at ",
   $query->popup_menu(-name=>"onDiskAt",
                      -values=>["","jlab","IUfarm","uconn"],
                      -labels=>{""=>"", "jlab"=>"Jlab",
                                "IUfarm"=>"IU farm", "uconn"=>"UConn"},
                      -default=>""),
   " but not at ",
   $query->popup_menu(-name=>"onDiskNotAt",
                      -values=>["","jlab","IUfarm","uconn"],
                      -labels=>{""=>"", "jlab"=>"Jlab",
                                "IUfarm"=>"IU farm", "uconn"=>"UConn"},
                      -default=>""),
   "</NOBR></TD></TR>\n",
   "<TR><TD><NOBR>",
   "whose trigger specification contains the expression ",
   $query->textfield(-name=>"triggerExpr",
                     -default=>"",
                     -size=>12,
                     -maxlength=>80),
   " not ",
   $query->textfield(-name=>"triggerNotExpr",
                     -default=>"",
                     -size=>12,
                     -maxlength=>80),
   "</NOBR></TD></TR>\n",
   "<TR><TD><NOBR>",
   "and whose comment contains the expression ",
   $query->textfield(-name=>"commentExpr",
                     -default=>"",
                     -size=>14,
                     -maxlength=>80),
   " not ",
   $query->textfield(-name=>"commentNotExpr",
                     -default=>"",
                     -size=>14,
                     -maxlength=>80),
   "</NOBR></TD></TR>\n",
   "<TR><TD>",
   "</TD></TR>",
   "</TABLE>\n",
   $query->end_form,"\n",
   $query->end_html,"\n";

# End of main

sub process_submitted_form
{
   my($query) = @_;

   $entry = $query->param('runCount');
   if ($entry =~ /^ *\d\d* *$/)
   {
      $runCount = $entry;
   }
   else
   {
      $runCount = 0;
   }
   $entry = $query->param('runStart');
   if ($entry =~ /^ *\d\d* *$/)
   {
      $runStart = $entry;
   }
   else
   {
      $runStart = 0;
   }
   $entry = $query->param('dayStart');
   if ($entry =~ /^ *\d\d* *$/)
   {
      $dayStart = $entry;
   }
   else
   {
      $dayStart = 0;
   }
   $entry = $query->param('timeStart');
   if ($entry =~ /^ *\d\d*:\d\d* *$/)
   {
      ($hour, $min) = split(":",$entry);
   }
   elsif ($entry =~ /^ *\d\d* *$/)
   {
      ($hour, $min) = ($entry, 0);
      
   }
   else
   {
      ($hour, $min) = (0, 0);
   }
   $wday = 0;
   $yday = 0;
   $time0 = mktime(0, $min, $hour, $dayStart,
                   $query->param('monthStart'),
                   $query->param('yearStart'),
                   $wday,$yday,-1);
   $entry = $query->param('dayStop');
   if ($entry =~ /^ *\d\d* *$/)
   {
      $dayStop = $entry;
   }
   else
   {
      $dayStop = 0;
   }
   $entry = $query->param('timeStop');
   if ($entry =~ /^ *\d\d*:\d\d* *$/)
   {
      ($hour, $min) = split(":",$entry);
   }
   elsif ($entry =~ /^ *\d\d* *$/)
   {
      ($hour, $min) = ($entry, 0);
      
   }
   else
   {
      ($hour, $min) = (0, 0);
   }
   $time1 = mktime(0, $min, $hour, $dayStop,
                   $query->param('monthStop'),
                   $query->param('yearStop'),
                   $wday,$yday,-1);
   $energySpec = $query->param('energySpec');
   $archivedOn = $query->param('archivedOn');
   $archivedNotOn = $query->param('archivedNotOn');
   $onDiskAt = $query->param('onDiskAt');
   $onDiskNotAt = $query->param('onDiskNotAt');
   $triggerExpr = $query->param('triggerExpr');
   $triggerNotExpr = $query->param('triggerNotExpr');
   $commentExpr = $query->param('commentExpr');
   $commentNotExpr = $query->param('commentNotExpr');

   @fullList = sort(keys %dbase);
   for ($i=0;$i<$#fullList;$i++)
   {
      ($run, $junk) = split('\.',$fullList[$i]);
      if ($run =~ /^\d+$/)
      {
         $run0 = $run;
         last;
      }
   }
   for ($i=$#fullList;$i>0;$i--)
   {
      ($run, $junk) = split('\.',$fullList[$i]);
      if ($run =~ /^\d+$/)
      {
         $run1 = $run;
         last;
      }
   }

   $#runList = -1;
   $run0 = ($runStart > $run0) ? $runStart : $run0;
   RUN: for ($run=$run0;$run<=$run1;$run++)
   {
      next RUN if (! $dbase{"$run.date"});
      next RUN if ($dbase{"$run.date"} < $time0);
      last RUN if ($dbase{"$run.date"} > $time1);
      if ("$energySpec")
      {
         next RUN if ($dbase{"$run.energy"} !~ /$energySpec/i)
      }
      if ("$archivedOn")
      {
         next RUN if ($dbase{"$run.archive"} !~ /$archivedOn/i)
      }
      if ("$archivedNotOn")
      {
         next RUN if ($dbase{"$run.archive"} =~ /$archivedNotOn/i)
      }
      if ("$triggerExpr")
      {
         next RUN if ($dbase{"$run.trigger"} !~ /$triggerExpr/i)
      }
      if ("$triggerNotExpr")
      {
         next RUN if ($dbase{"$run.trigger"} =~ /$triggerNotExpr/i)
      }
      if ("$commentExpr")
      {
         next RUN if ($dbase{"$run.comment"} !~ /$commentExpr/i)
      }
      if ("$commentNotExpr")
      {
         next RUN if ($dbase{"$run.comment"} =~ /$commentNotExpr/i)
      }
      if ("$onDiskAt")
      {
         next RUN if (grep(/$run\.$onDiskAt/,@fullList) == 0);
      }
      if ("$onDiskNotAt")
      {
         next RUN if (grep(/$run\.$onDiskNotAt/,@fullList) != 0);
      }
      $runList[++$#runList] = $run;
      $query->param("runStart",$run+1);
      last RUN if ($#runList == $runCount-1);
   }
   if ($#runList < 0)
   {
      print "<P ALIGN=CENTER>\n",
            "&nbsp;<BR>&nbsp;<BR>&nbsp;<BR>\n",
            "This search returned no entries<BR>\n";
      return;
   }

   print "<TABLE ALIGN=CENTER BGCOLOR=ffffff ",
         "CELLPADDING=3 BORDER=4>\n",
         "<TR><TH>Run</TH><TH>Energy</TH>\n",
         "<TH>Trigger</TH><TH>Date</TH>\n",
         "<TH>File</TH><TH>Archive</TH></TR>\n";
   foreach $run (@runList)
   {
      my $name = $ENV{"REQUEST_URI"};
      $name =~ s/\/$//;
      print "<TR>\n",
            "<TD ALIGN=CENTER VALIGN=CENTER ROWSPAN=2>\n",
            "<A HREF=\"$name/$run\">$run</A></TD>\n",
            "<TD ALIGN=CENTER>",$dbase{"$run.energy"},"</TD>\n",
            "<TD>",$dbase{"$run.trigger"},"</TD>\n",
            "<TD>",ctime($dbase{"$run.date"}),"</TD>\n",
            "<TD>",$dbase{"$run.filename"},
            " (",$dbase{"$run.size"},") </TD>\n",
            "<TD ALIGN=CENTER>",$dbase{"$run.archive"},"</TD>\n",
            "</TR>\n",
            "<TR>\n",
            "<TD COLSPAN=5>&nbsp;",$dbase{"$run.comment"},"</TD>\n",
            "</TR>\n";
   }
   print "</TABLE>\n";
}

sub print_detailed_info
{
   my $run = $_[1];
   print "<FONT SIZE=+2><U>Summary for run $run:</U></FONT><BR><BR>\n",
         "<UL>\n",
         "<LI><B>Beam Energy:</B> ",$dbase{"$run.energy"},"\n",
         "<LI><B>Trigger:</B> ",$dbase{"$run.trigger"},"\n",
         "<LI><B>Comment:</B> ",$dbase{"$run.comment"},"\n",
         "<LI><B>Archive:</B> ",$dbase{"$run.archive"},"\n",
         "<BR>&nbsp;<BR>\n";
   if ($dbase{"$run.archive"} =~ /mss/i) {
      print
         "<LI><B>Status on mss:</B>\n<BR><UL>\n",
         " <LI><B>created:</B> ",ctime($dbase{"$run.date"}),"\n",
         " <LI><B>data file:</B> ",$dbase{"$run.filename"},
         " (",$dbase{"$run.size"},")\n",
         "</UL>&nbsp;<BR>\n";
   }
   if ($dbase{"$run.archive"} =~ /hpss/i) {
      print
         "<LI><B>Status on hpss:</B>\n<BR><UL>\n",
         " <LI><B>created:</B> ",ctime($dbase{"$run.hpss.date"}),"\n",
         " <LI><B>data file:</B> ",$dbase{"$run.hpss.filename"},
         " (",size_in_units($dbase{"$run.hpss.size"}),")\n",
         "</UL>&nbsp;<BR>\n";
   }
   %sites = ("uconn","UConn","IUfarm","IU farm","jlab","Jlab");
   foreach $site (keys %sites) {
      local %areas = qw(pvfs0 /pvfs0 pvfs1 /pvfs1
                        raid /radphi raw /radphi/raw_data
                        scratch /work/hallb/e94016
                        jcache /cache/hallb/e94016/raw);
      local @list = grep(/$run\.$site.*.size/, keys %dbase);
      if (@list > 0) {
         print "<LI><B>Status on disk at $sites{$site}:</B>\n<BR><UL>\n";
         foreach $key (@list) {
            my ($run, $site, $area, $text);
            ($run, $site, $area, $text) = split(/\./,$key);
            print
         " <LI><B>created:</B> ",ctime($dbase{"$run.$site.$area.date"}),"\n",
         " <LI><B>found at:</B> ",$areas{$area},
         "/",$dbase{"$run.$site.$area.filename"},
         " (",size_in_units($dbase{"$run.$site.$area.size"}),")\n";
         }
         print "</UL>&nbsp;<BR>\n";
      }
   }
   if (exists $dbase{"$run.epics.TimeStamp"}) {
      print "<LI><B>Epics information from run header:</B>\n",
            "time stamp ",$dbase{"$run.epics.TimeStamp"},"\n<BR><UL>\n";
      for $name (split ':',$dbase{"epics.names"}) {
         if (exists $dbase{"$run.epics.$name"}) {
            print " <LI>$name = ",$dbase{"$run.epics.$name"},"\n";
         }
      }
      print "</UL>&nbsp;<BR>\n";
   }
   print $query->end_form,"\n",
         $query->end_html,"\n";
}

sub size_in_units
{
   my $size = "@_";
   if ($size > 1024*1024)
   {
      $size = sprintf("%.1f MB",$size/1024./1024.);
   }
   elsif ($size > 1024)
   {
      $size = sprintf("%.1f KB",$size/1024.);
   }
   else
   {
      $size = "$size B";
   }
}

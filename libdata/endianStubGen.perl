#!perl5

$sccsid = "@(#)endianStubGen.perl\t1.5\tCreated 10/1/96 22:47:08";

$verbose = grep ("-v", @ARGV);

$outName = "endianConverterAuto.c";
$hdrName = "endianConverterAuto.h";
$depName = "$outName.depend";

$home = $ENV{"RADPHI_HOME"};

die "Error: env. RADPHI_HOME is not defined.\n" if ! defined $home;

$includeDir = "$home/include";
$includeMain = "$includeDir/disData.h";

die "Error: $includeDir is not a directory.\n" if ! -d $includeDir;
die "Error: cannot read $includeMain: $!\n" if ! -r $includeMain;

print STDERR "Reading $includeMain...\n" if $verbose;

open(IN,"<$includeMain") || die "Error: cannot read $includeMain: $!\n";

$headers{"disData.h"} = "itape_header_t";
$types{"itape_header_t"} = "disData.h";

$headers{"disData.h"} .= "group_header_t";
$types{"group_header_t"} = "disData.h";

$headers{"tapeData.h"} = "tapeHeader_t";
$types{"tapeHeader_t"} = "tapeData.h";

while (<IN>)
{
    chop; # remove the \n at the end.

    next if ! /^#define/;
    next if ! /GROUP_/;

    next if /IGNORE/;

#    print $_,"\n";

    s/\t/ /g;

    (s/ +/ /g);

    /^#define (GROUP_.*) ([0-9]+)(.*)/;

#    print "Parse: <$1> <$2> <$3>\n";

    $group = int($2);

    if (($group > 0)||($2 eq "0"))
    {
	$name = $1;
	$comm = $3;

	$comm =~ s/\)/ /g;
	$comm =~ s/\(/ /g;
	$comm =~ s/:/ /g;
	$comm =~ s/;/ /g;
	$comm =~ s/\*/ /g;

	@comments = split(" ",$comm);

	$type = "none";
	$header = "none";

	foreach $comment (sort @comments)
	{
	    $type = $comment if $comment =~ /_t$/;
	    $header = $comment if $comment =~ /\.h$/;

#	    print "Comment: $comment\n";
	}

#	print "Group: $group, name: $name, type: $type, header: $header\n";

	if ($type eq "none")
	{
	    print STDERR "Error: Group: $name, No type information\n";
	    next;
	}

	if ($header eq "none")
	{
	    print STDERR "Error: Group: $name, No header information\n";
	    next;
	}

	if (defined $numbers{$group})
	{
	    print STDERR "Error: Duplicated group: $group, new name: $name, previous name: $numbers{$group}\n";
	    next;
	}

	$gnumber{$name} = $group;
	$gtype{$name}   = $type;
	$gheader{$name} = $header;

	$numbers{$group} = $name;
	$headers{$header} .= " $type" if ! defined $types{$type};
        foreach $item (split(",",$type))
        {
	    $types{$item} = $header if ! defined $types{$item};
        }

	print STDERR "Group: $group, name: $name, type: $type, header: $header\n" if $verbose;
    }
    else
    {
	/^#define (GROUP_.*) (GROUP_.*)/;

	#print "<$1> <$2>\n";

	if (defined $gnumber{$2})
	{
	    #print STDERR "Warning: Duplicated group: $1, same as $2.\n";
	    next;
	}

	print STDERR "Error parsing line: <$_>: Parse: <$1> <$2> <$3>\n";
    }
}

close(IN);

print STDERR "Reading header files...\n";

foreach $header (sort keys %headers)
{
    $r = "$includeDir/$header";

    print STDERR "Header: $header, types: $headers{$header}\n" if $verbose;

    if (! -r $r)
    {
	print STDERR "Error: Cannot read $r: $!\n";
	next;
    }

    open (IN, "< $r") || die "Error: Cannot read $r: $!\n";

    $file = "";

    while (<IN>)
    {
	s/\t/ /g;      # replace TABs with spaces

	$file .= $_;
    }

    close IN;

#    print STDERR "File before: [$file]\n";

    $file =~ s/\/\*(.*?)\*\///gs; # remove comments
    $file =~ s/ +/ /gs;           # remove extra spaces

#    print STDERR "File after: [$file]\n";

    @lines = split(/\n/,$file);

    $str = "";
    $inStr = 0;

    foreach $_ (@lines)
    {
	$_ .= "\n";

	next if /^\n/;  # kill empty lines
	next if /^ \n/; # kill empty lines

	chop; # remove trailing \n

	next if /^int/;
	next if /^void/;
	next if /^#include/;
	next if /^#define/;
	next if /^#endif/;
	next if /^#if/;

	#print "Read: ",$_,"\n";

	if (/^typedef (.*) (.*);$/)
	{
	    print STDERR "typedef: ($1) ($2)\n" if $verbose;

	    $typedef{$2} = $typedef{$1};
	}

	$inStr = 1 if /typedef/;
	$inStr = 1 if /{/;

	$str .= $_ if $inStr;

	$inStr = 2 if /}/;

	if ($inStr == 2)
	{
#	    $str =~ s/\/\*(.*)\*\///g;  # remove comments
#	    $str =~ s/ +/ /g;           # remove all extra spaces

	    $str =~ /(.*){(.*)}(.*);/;      # match: "typedef xxx { content } name;"

	    #print STDERR "Str: ($str) ($1) ($2) ($3)\n";

	    $union   = $1;
	    $content = $2;
	    $name    = $3;      $name =~ s/ +//;

	    if ($union =~ /struct/)
	    {
		$typedef{$name} = $content;

		print STDERR "Header: $header, Struct: $name, Content: [$content]\n" if $verbose;
	    }

	    if ($union =~ /union/)
	    {
		$typedef{$name} = "undefined";

		$typeeq{$name} = "INTARRAY_t";

		print STDERR "Header: $header, Union: $name.\n" if $verbose;
		print STDERR "Warning: Type $name is a union: don't know how to convert unions, will use INTARRAY_t and hope for the best.\n";
	    }

	    $inStr = 0;
	    $str = "";
	}
    }
}

print STDERR "Creating fake types...\n";

$external{"NOCONVERT_t"} = "void";
$external{"INTARRAY_t"} = "void";
$external{"uint8"}      = "uint8";
$external{"uint16"}     = "uint16";
$external{"uint32"}     = "uint32";
$external{"double"}     = "double";
$external{"vector3_t"}  = "vector3_t";
$external{"vector4_t"}  = "vector4_t";

$typeeq{"int"} = "uint32";
$typeeq{"int16"} = "uint16";
$typeeq{"int32"} = "uint32";
$typeeq{"Particle_t"} = "uint32";
$typeeq{"eventClass_t"} = "uint32";
$typeeq{"float"} = "uint32";
$typeeq{"float32"} = "uint32";
$typeeq{"char"} = "uint8";
#$typeeq{"lgd_mam_t"} = "mam_list_t";

print STDERR "Generating endian-converter routines:\n";

open (OUT,">$outName") || die "Error: Cannot write to $outName: $!\n";

open (HDR,">$hdrName") || die "Error: Cannot write to $hdrName: $!\n";

open (DEP,">$depName") || die "Error: Cannot write to $depName: $!\n";

print DEP "$outName: \\\n";

print OUT "/*\n";
print OUT "   $outName --- endian converter\n";
print OUT "\n";
print OUT "   this file is automatically generated by $sccsid\n";
print OUT "*/\n";
print OUT "\n";
print OUT "\n";
print OUT "#include <stdio.h>\n";
print OUT "\n";
print OUT "#include <$hdrName>\n";
print OUT "\n";

print DEP "\t$hdrName\\\n";

print HDR "/*\n";
print HDR "   $hdrName --- endian converter\n";
print HDR "\n";
print HDR "   this file is automatically generated by 'endianStubGen.perl'\n";
print HDR "*/\n";
print HDR "\n";
print HDR "\n";
print HDR "#include <disData.h>\n";
print HDR "\n";

print DEP "\t$includeMain\\\n";

print STDERR "Generating type converters...\n";

print HDR "\n";

foreach $type (sort keys %external)
{
    $typedef{$type} = "external";

    print HDR "extern int endian_convertType1_$type($external{$type}*dptr,int dlength);\n";
    print HDR "extern int endian_convertTypeN_$type($external{$type}*dptr,int dcount);\n";
}

print HDR "\n";

@alltypes = sort keys %types;

for ($itype=0; $itype<=$#alltypes; $itype++)
{
    $type = $alltypes[$itype];

    print STDERR "Type: $type\n" if $verbose;

    $def = $typedef{$type};

    $def = $typedef{$typeeq{$type}} if defined $typeeq{$type};

    if (! defined $def)
    {
	print STDERR "Error: No definition for type: $type.\n";
	next;
    }

    next if $def eq "external";

    #print STDERR "Def: $def\n";

    @def = split(";",$def);

    $i = 0;

    foreach $d (@def)
    {
	$d =~ s/^ +//g;

	($t,$n) = split(" ",$d);

	next if $t eq "";

	@n1 = ($n);

	if ($n =~ /,/)
	{
	    @n1 = split(",",$n);
	}

	foreach $n (@n1)
	{
	    $c = 1;
	    $a = 0;

	    if ($n =~ /(.*)\[(.*)\]/)
	    {
		$n = $1;
		$c = $2;
		$a = 1;
	    }

	    $d[$i] = $d;
	    $t[$i] = $t;
	    $n[$i] = $n;
	    $c[$i] = $c;
	    $a[$i] = $a;

	    $i++;
	}
    }

    $count = $i;

    $include = $types{$type};
    $include = "" if $include eq "predefined";
    $include = "" if $includeonce{$include};

    print HDR "\n";

    if ($include ne "")
    {
	$includeonce{$include} = 1;
	print HDR "#include <$include>\n" unless $include eq "";
	print HDR "\n";

	print DEP "\t$includeDir/$include\\\n";
    }

    $out = "";
    $out1 = "";
    $out2 = "";

    $isNotAnArray = 1;

    for ($i=0; $i<$count; $i++)
    {
	$t = $typeeq{$t[$i]};
	$t = $t[$i] if !defined $t;

	if (defined $typedef{$t})
	{
	    if ((!defined($types{$t}))&&($typedef{$t} ne "external"))
	    {
		$types{$t} = "";
		$alltypes[++$#alltypes] = $t;

		$out1 .=  "int endian_convertType1_$t($t*dptr,int dlength);\n";
		$out1 .=  "int endian_convertTypeN_$t($t dptr[],int dcount);\n";
	    }
	}
	else
	{
	    $missing{$t} .= " $type";
	}

	$cast = "";
	$cast = "($external{$t}"."*)" if defined $external{$t};

	if ($a[$i] == 0)
	{
	    $out2 .= "   /* ($d[$i])\t --> ($t[$i]|$t)\t($n[$i]) */\n";
	
	    $out .=  "   endian_convertType1_$t($cast"."&dptr->$n[$i],sizeof(dptr->$n[$i]));\n";
	}
	else
	{
	    $c = $c[$i];

	    if (($i==$count-1)&&($c == 1))
	    {
		# suspect it is an array-type group

		if ($count == 2)
		{
		    $c = "dptr->$n[0]";
		    $isNotAnArray = 0;
		}
		elsif ($n[$i-1] =~ /^n/)
		{
		    $c = "dptr->$n[$i-1]";
		    $isNotAnArray = 0;
		}
		else
		{
		    $jj = 0;
		    $jjindex = (-1);

		    for ($j=0; $j < $i; $j++)
		    {
			if ($n[$j] =~ /^n/)
			{
			    $jj ++;
			    $jjindex = $j if $jj == 1;
			}
		    }

		    $c = "dptr->$n[$jjindex]";
		    $isNotAnArray = 0;
		}
	    }

	    if (($isNotAnArray == 0)&&($count>2))
	    {
		print STDERR "Notice: Found an array type: $type: $n[$i]\[$c\].\n";
	    }

	    $out2 .= "   /* ($d[$i])\t --> ($t[$i]|$t)\t($n[$i])[$c] */\n";
	
	    $out .=  "   endian_convertTypeN_$t($cast"."dptr->$n[$i],$c);\n";
	}
    }

    print HDR $out1;
    print HDR "\n";
    print HDR "int endian_convertType1_$type($type *dptr,int dlength);\n";
    print HDR "int endian_convertTypeN_$type($type *dptr,int dcount);\n" if $isNotAnArray;

    print OUT "\n";
    print OUT "int endian_convertType1_$type($type *dptr,int dlength)\n";
    print OUT "{\n";
#    print OUT "\n";
#    print OUT "   /* definition: $def */\n";
#    print OUT "\n";

    print OUT "\n";
    print OUT $out2;
    print OUT "\n";
    print OUT $out;
    print OUT "\n";
    print OUT "   return 0;\n";
    print OUT "}\n";

    if ($isNotAnArray)
    {
	print OUT "\n";
	print OUT "\n";
	print OUT "int endian_convertTypeN_$type($type dptr[],int dcount)\n";
	print OUT "{\n";
	print OUT "   int i;\n";
	print OUT "\n";
	print OUT "   for (i=0; i<dcount; i++)\n";
	print OUT "      endian_convertType1_$type(&dptr[i],sizeof(dptr[i]));\n";
	print OUT "\n";
	print OUT "   return 0;\n";
	print OUT "}\n";
	print OUT "\n";
    }

}

foreach $m (sort keys %missing)
{
    print STDERR "Missing type definition: $m, needed for: $missing{$m}\n";
}

print STDERR "Generating group converters...\n";

print OUT "\n";
print OUT "\n";
print OUT "\n";
print OUT "\n";

print HDR "\n";
print HDR "int endian_convertGroup(group_header_t *hdr,void *dptr,int dlength);\n";
print HDR "\n";


print OUT "int endian_convertGroup(group_header_t *hdr,void *dptr,int dlength)\n";
print OUT "{\n";
print OUT "  int group=hdr->type;\n";
print OUT "  int subgroup=hdr->subtype;\n";
print OUT "\n";
print OUT "  switch (group)\n";
print OUT "    {\n";

foreach $name (sort keys %gnumber)
{
    print STDERR "Group: $name, type: $type\n" if $verbose;

    print OUT "      case $name:\n";
    print OUT "         switch (subgroup)\n";
    print OUT "           {\n";

    @types = split(",",$gtype{$name});
    for ($subtype = 0; $subtype < @types; $subtype++)
    {
        $type = $types[$subtype];
        $type = $typeeq{$type} if defined $typeeq{$type};
        print OUT "          case $subtype:\n";
        print OUT "            return endian_convertType1_$type(dptr,dlength);\n";
    }
    print OUT "          default:\n";
    print OUT "            fprintf(stderr,\"$outName:endian_convertGroup: Error: Unknown subgroup %d found in group $name.\\n\",subgroup);\n";
    print OUT "            return -1;\n";
    print OUT "          }\n";
}

print OUT "    }\n";

print OUT "\n";
print OUT "  fprintf(stderr,\"$outName:endian_convertGroup: Error: Don't know how to convert group %d.\\n\",group);\n";
print OUT "\n";
print OUT "  return -1;\n";
print OUT "}\n";
print OUT "\n";
print OUT "/* end file */\n";
close OUT;

print HDR "\n";
print HDR "/* end file */\n";
close HDR;

print DEP "\n";
close DEP;

print STDERR "Done.\n";

#system "cat $outName";

#end file


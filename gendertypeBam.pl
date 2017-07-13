#!/usr/bin/perl

use strict;
use warnings;

my $bamFile = $ARGV[0];

my $chromCount = `/oicr/local/analysis/sw//samtools/samtools-0.1.18/samtools view $bamFile | cut -f 3 | uniq -c`;

my $fRatioCut = 40;		# x/y ratio greater than this is typed female
my $mRatioCut = 10;		# x/y ratio less than this is male


my $xCount = 0;
my $yCount = 0;

for my $l (split(/\n/, $chromCount))
{
	if ($l =~ /(.*) chrX/)
	{
		$xCount = $1;
		$xCount =~ s/^ *//;
	}
	if ($l =~ /(.*) chrY/)
	{
		$yCount = $1;
		$yCount =~ s/^ *//;
	}
}

$bamFile =~ s/.*\///;


my $ratio = 0;
if ($yCount > 0)
{
	$ratio = $xCount / $yCount;
}

my $gender = "Unknown";

if ($ratio > $fRatioCut)
{
	$gender = "Female";
}
elsif (($ratio < $mRatioCut) and ($ratio > 0))		# ratio of zero or less is not a callable outcome
{
	$gender = "Male";
}


print "$bamFile\t$xCount\t$yCount\t$ratio\t$gender\n";


#!/usr/bin/perl

use strict;
use warnings;

my $bamFile = $ARGV[0];
my $breakpointA = $ARGV[1];		# as chr12:25398281+
my $breakpointB = $ARGV[2];		# as chr17:928830200-

my $window = 500;
my $samtools = "/oicr/local/analysis/sw//samtools/samtools-0.1.18/samtools";

unless (-e $bamFile)
{
	die "$bamFile doesn't exist\n";
}

my ($chrA, $chrB, $startA, $startB, $endA, $endB, $strandA, $strandB);

if ($breakpointA =~ /(chr.*):(.*)(.)/)
{
	$chrA = $1;
	$strandA = $3;

	if ($strandA eq "-")
	{
		$startA = $2;
		$endA = $startA + $window;
	}
		elsif ($strandA eq "+")
	{
		$endA = $2;
		$startA = $endA - $window;
	}
	else
	{
		die "Couldn't parse $breakpointA, expected chr12:25398281+\n";
	}
}
else
{
	die "Couldn't parse $breakpointA, expected chr12:25398281+\n";
}

if ($breakpointB =~ /(chr.*):(.*)(.)/)
{
	$chrB = $1;
	$strandB = $3;

	if ($strandB eq "+")
	{
		$startB = $2;
		$endB = $startB + $window;
	}
	elsif ($strandB eq "-")
	{
		$endB = $2;
		$startB = $endB - $window;
	}
	else
	{
		die "Couldn't parse $breakpointB, expected chr12:25398281+\n";
	}
}
else
{
	die "Couldn't parse $breakpointB, expected chr12:25398281+\n";
}




my $bpAreads = `$samtools view -F 4 $bamFile $chrA:$startA-$endA`;
my $bpBreads = `$samtools view -F 4 $bamFile $chrB:$startB-$endB`;
#print "$bpAreads\n";

if ($chrA eq $chrB)
{
	$chrA = "=";
	$chrB = "=";
}

my @f;
for my $read (split/\n/, $bpAreads)
{
	@f = split(/\t/, $read);

	if (($f[6] eq $chrB) and ($f[7] >= $startB) and ($f[7] <= $endB))
	{
		print $read . "\n";
	}
}

for my $read (split/\n/, $bpBreads)
{
	@f = split(/\t/, $read);
	if (($f[6] eq $chrA) and ($f[7] >= $startA) and ($f[7] <= $endA))
	{
		print $read . "\n";
	}
}


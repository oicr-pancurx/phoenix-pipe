#!/usr/bin/perl

use strict;
use warnings;

my $minSoftClips = 3;
my $minSupportingPairs = 3;
my $minNormalPairs = 1;

my $minSize = 125;
my $pairEvidenceMinSize = 1000;		# if event is too small (and an INS or DEL) normal pairs may span both breakpoints as well as discordant ones

# read CREST file from STDIN
my $normalBam = $ARGV[0];
my $tumourBam = $ARGV[1];


# CREST Filtering			Discordant Pairs			
# 	Size	Softclip	Tumour	Normal	CREST Orientation	Expected Pair Orientation
# 	Small DEL	<1000	>=3 on both	>=0	>=0	+ +	F1R2/R1F2
# 	DEL	>1000	>=0	>=3	<1	+ +	F1R2/R1F2
# 	ITX	>125	>=0	>=3	<1	+ - or - +	F1F2/F2F1 (+ -) or R1R2/R2R1  (- +)
# 	INV	>125	>=0	>=3	<1	+ - or - +	F1F2/F2F1 (+ -) or R1R2/R2R1  (- +)
# 	Small INS	<1000 and > 125	>=3 on both	>=0	>=0	+ +	F1R2/R1F2
# 	INS	>1000	>=0	>=3	<1	+ +	F1R2/R1F2
# 	CTX		>=0	>=3	<1	anything	anything


my $l;
my @f;

my $pass;

my ($chr1,$pos1,$strand1,$softClips1,$chr2,$pos2,$strand2,$softClips2,$type,$size);
my ($tumourSupportingPairs,$normalSupportingPairs);

while ($l = <STDIN>)
{
	$pass = 1;

	chomp $l;
	@f = split(/\t/, $l);

	$chr1 = $f[0];
	$pos1 = $f[1];
	$strand1 = $f[2];
	$softClips1 = $f[3];

	$chr2 = $f[4];
	$pos2 = $f[5];
	$strand2 = $f[6];
	$softClips2 = $f[7];

	$type = $f[8];

	$size = abs($pos2 - $pos1);

	$tumourSupportingPairs = countSupportingPairs($tumourBam,$chr1,$pos1,$strand1,$chr2,$pos2,$strand2);
	$normalSupportingPairs = countSupportingPairs($normalBam,$chr1,$pos1,$strand1,$chr2,$pos2,$strand2);

	if (($chr1 =~ /chrUn/) or ($chr2 =~ /chrUn/))
	{
		$pass = 0;		# don't know how to handle unknown chromosomes downstream
	}

	if (($type ne "DEL") and ($size < $minSize))
	{
		$pass = 0;		# all non-deletion events have a minimum size
	}

	if (($type eq "DEL") or ($type eq "INS"))
	{
		if ($size < $pairEvidenceMinSize)
		{
			if (($softClips1 < $minSoftClips) or ($softClips2 < $minSoftClips))
			{
				$pass = 0;		# small insertions or deletions must have soft clip support for both breakpoints
			}
		}
		else		# large INS or DEL
		{
			if (($tumourSupportingPairs < $minSupportingPairs) or ($normalSupportingPairs > $minNormalPairs))
			{
				$pass = 0;		# large insertions or deletions must have supporting discordant pairs in the tumour and no supporting discordant pairs in the normal
			}
		}
	}
	else		# ITX, INV or CTX
	{
		if (($tumourSupportingPairs < $minSupportingPairs) or ($normalSupportingPairs > $minNormalPairs))
		{
			$pass = 0;		# translocations and inversions must have supporting discordant pairs in the tumour and no supporting discordant pairs in the normal
		}
	}


	if ($pass == 1)
	{
		print "$chr1\t$pos1\t$strand1\t$softClips1\t$chr2\t$pos2\t$strand2\t$softClips2\t$type\t$size\t$tumourSupportingPairs\t$normalSupportingPairs\tPASS\n";
	}
	else
	{
		print "$chr1\t$pos1\t$strand1\t$softClips1\t$chr2\t$pos2\t$strand2\t$softClips2\t$type\t$size\t$tumourSupportingPairs\t$normalSupportingPairs\tFAIL\n";
	}
}


sub countSupportingPairs
{
	my $bamFile = shift;
	my $chr1 = shift;
	my $pos1 = shift;
	my $strand1 = shift;
	my $chr2 = shift;
	my $pos2 = shift;
	my $strand2 = shift;

	my @reads = listReadsAtBPs($bamFile,"$chr1:$pos1$strand1","$chr2:$pos2$strand2");
	my %supportingPairs;

	my @f;
	my $readName;
	my $readReverse;
	my $mateReverse;

	for my $r (@reads)
	{
		@f = split(/\t/, $r);
		$readName = $f[0];

		if ($f[1] & 16)
		{
			$readReverse = 1;
		}
		else
		{
			$readReverse = 0;
		}

		if ($f[1] & 32)
		{
			$mateReverse = 1;
		}
		else
		{
			$mateReverse = 0;
		}

#		print "$readReverse $mateReverse $r\n";

		if ("$strand1 $strand2" eq "+ +")
		{
			if ($readReverse != $mateReverse)
			{
				$supportingPairs{$readName}++;
			}
		}
		elsif ("$strand1 $strand2" eq "+ -")
		{
			if (($readReverse == 0) and ($mateReverse == 0))
			{
				$supportingPairs{$readName}++;
			}
		}
		elsif ("$strand1 $strand2" eq "- +")
		{
			if (($readReverse == 1) and ($mateReverse == 1))
			{
				$supportingPairs{$readName}++;
			}
		}
	}

	return scalar(keys %supportingPairs);
}



sub listReadsAtBPs
{
	my $bamFile = $_[0];
	my $breakpointA = $_[1];		# as chr12:25398281+
	my $breakpointB = $_[2];		# as chr17:928830200-
	
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
	
	my @resultReads;

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
			push(@resultReads, $read);
		}
	}
	
	for my $read (split/\n/, $bpBreads)
	{
		@f = split(/\t/, $read);
		if (($f[6] eq $chrA) and ($f[7] >= $startA) and ($f[7] <= $endA))
		{
			push(@resultReads, $read);
		}
	}

	return @resultReads;
}


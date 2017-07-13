#!/usr/bin/perl

use strict;
use warnings;
use Tabix;

my @f;
my $l;
my @db;

my $nearSize = 5000;
my $refSeq = $ARGV[0];

my ($chr1, $pos1, $chr2, $pos2, $type);
my $dbSNPstring;

my $tabix = Tabix->new(-data => $refSeq, -index => "$refSeq.tbi");
my $iter;

my $leftSplit;
my $rightSplit;
my $leftNear;
my $rightNear;
my $between;

while ($l = <STDIN>)
{
	chomp $l;
	@f = split(/\t/, $l);

	if ($l =~ /^#/)
	{
		print $l . "\n";
	}
	elsif ($l =~ /^chrom1/)
	{
		print "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t";
		print "split_left_bp\tsplit_right_bp\tnear_left_bp\tnear_right_bp\tbetween_bps";
		for (my $i = 5; $i < scalar @f; $i++)
		{
			print "\t$f[$i]";
		}
		print "\n";
	}
	else
	{
		$chr1 = $f[0];
		$pos1 = $f[1];
		$chr2 = $f[2];
		$pos2 = $f[3];
		$type = $f[4];

		# test for genes split by the left or right bp

		$leftSplit = genesInRange($chr1, $pos1 - 1, $pos1, $tabix);
		$rightSplit = genesInRange($chr2, $pos2 - 1, $pos2, $tabix);


		# test for genes within $nearSize bp of the left or right bp
		
		$leftNear = genesInRange($chr1, $pos1 - $nearSize, $pos1 + $nearSize, $tabix);
		$rightNear = genesInRange($chr2, $pos2 - $nearSize, $pos2 + $nearSize, $tabix);


		# test for genes that are within DUP or DEL type events
		$between = ".";
		if (($type eq "DEL") or ($type eq "DUP"))
		{
			if ($pos1 < $pos2)
			{
				$between = genesInRange($chr1, $pos1, $pos2, $tabix);
			}
			else
			{
				$between = genesInRange($chr1, $pos2, $pos1, $tabix);
			}
		}


		print "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t";
		print "$leftSplit\t$rightSplit\t$leftNear\t$rightNear\t$between";
		for (my $i = 5; $i < scalar @f; $i++)
		{
			print "\t$f[$i]";
		}
		print "\n";
	}

}


sub genesInRange
{
	my $chr = shift;
	my $start = shift;
	my $end = shift;

	$iter = $tabix->query($chr, $start, $end);

	my @fields;
	my @geneInfo;

	my $genes = ".";

	if (defined $iter->{"_"})		# happens if the contig isn't in the bed file
	{
		while ($l = $tabix->read($iter))
		{
			chomp $l;
			@fields = split(/\t/, $l);
			@geneInfo = split(/,/, $fields[3]);

			if ($genes eq ".")
			{
				$genes = "$geneInfo[0]";
			}
			else
			{
				$genes = "$genes,$geneInfo[0]";
			}
		}
	}

	return $genes;
}




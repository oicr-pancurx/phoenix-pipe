#!/usr/bin/perl

use strict;
use warnings;
use Tabix;

my $bed = $ARGV[0];
#my $bed = "/oicr/data/reference/genomes/homo_sapiens_mc/Agilent/SureSelectHumanAllExonV4/S03723314_Regions.merged.bed.gz";

my $l;
my @f;

my $bedReturn;

my $tabix = Tabix->new(-data => $bed, -index => "$bed.tbi");
my $iter;


while ($l = <STDIN>)
{
	chomp $l;
	@f = split(/\t/, $l);

	if ($l =~ /^#/)
	{
		print $l . "\n";
	}
	else
	{
		$iter = $tabix->query($f[0],$f[1] - 1,$f[1]);

		if (defined $iter->{"_"})		# happens if the contig isn't in the bed file
		{
			if ($bedReturn = $tabix->read($iter))
			{
				print $l . "\n";
			}
		}
	}
}




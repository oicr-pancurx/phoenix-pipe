#!/usr/bin/perl

use strict;
use warnings;
use Tabix;
use Data::Dumper;

my $bed = $ARGV[0];
my $trackName = $ARGV[1];

my $l;
my @f;

my $chr;
my $start;
my $end;

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
		$chr = $f[0];
		$start = $f[1] - 1;
		$end = $f[1];


		$iter = $tabix->query($chr,$start,$end);

		if (defined $iter->{"_"})		# happens if the contig isn't in the bed file
		{
			if ($bedReturn = $tabix->read($iter))
			{
				$f[7] .= ";TRACK=$trackName";
				
				print $f[0];
				for (my $i = 1; $i < scalar(@f); $i++)
				{
					print "\t$f[$i]";
				}
				print "\n";
			}
			else
			{
				print $l . "\n";
			}
		}
		else
		{
			print $l . "\n";
		}
	}
}





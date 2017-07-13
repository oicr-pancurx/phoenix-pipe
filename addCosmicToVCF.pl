#!/usr/bin/perl

use strict;
use warnings;
use Tabix;

my @f;
my $l;
my @db;

my $cosmicLine;

my $cosmic = $ARGV[0];

my ($chr, $pos, $ref, $alt, $dbRef, $dbAlt);
my $dbSNPstring;

my $tabix = Tabix->new(-data => $cosmic, -index => "$cosmic.tbi");
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
		$pos = $f[1];
		$ref = $f[3];
		$alt = $f[4];

		$iter = $tabix->query($chr,$pos - 1,$pos);

		if (defined $iter->{"_"})		# happens if the contig isn't in the bed file
		{
			while ($cosmicLine = $tabix->read($iter))
			{
				@db = split(/\t/, $cosmicLine);
				for my $dbAlt (split(/,/, $db[4]))
				{
					if (($db[3] eq $f[3]) and ($dbAlt eq $f[4]))
					{
						$f[7] .= ";COSMIC=$db[2]";
					}
				}
			}
		}

		print $f[0];
		for (my $i = 1; $i < scalar (@f); $i++)
		{
			print "\t$f[$i]";
		}
		print "\n";
	}

}





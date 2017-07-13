#!/usr/bin/perl

use strict;
use warnings;


my $l;
my @f;

my $nSamp;
my $nCol = -1;


while ($l = <STDIN>)
{
	chomp $l;

	if ($l =~ /matched_sample_id="(.*)"/)
	{
		$nSamp = $1;
		print "$l\n";
	}
	elsif ($l =~ /^#CHROM/)
	{
		@f = split(/\t/, $l);

		for (my $i = 0; $i < scalar @f; $i++)
		{
			if ($f[$i] eq $nSamp)
			{
				$nCol = $i;
			}
		}

		print "$l\n";
	}
	elsif ($l =~ /^#/)
	{
		print "$l\n";
	}
	else
	{
		if ($nCol == -1)
		{
			die "Couldn't detect normal column\n";
		}

		@f = split(/\t/, $l);

		if (($f[$nCol] =~ /^0\/1/) or ($f[$nCol] =~ /^1\/1/))
		{
			print "$l\n";
		}
	}
}




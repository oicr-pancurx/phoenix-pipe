#!/usr/bin/perl

use strict;
use warnings;

my $l;
my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$gt1,$gt2);

my ($nor,$tum);
my %gt;

my ($depth,$freq);

my $minDepth = 10;
my $minQual = 499;
my $minFreq = 0.33;
my $maxFreq = 0.66;

print "CHR	POS	REF_COUNT	VAR_COUNT\n";

while ($l = <STDIN>)
{
	chomp $l;

	if ($l =~ /##DCC=<analyzed_sample_id="(.*)">/)
	{
		$tum = $1;
	}
	elsif ($l =~ /##DCC=<matched_sample_id="(.*)">/)
	{
		$nor = $1;
	}
	elsif ($l =~ /^#CHROM/)
	{
		($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$gt1,$gt2) = split(/\t/, $l);
	}
	elsif ($l =~ /^#/)
	{
	}
	else
	{
		($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$gt{$gt1},$gt{$gt2}) = split(/\t/, $l);

		if ($gt{$nor} =~ /^0\/1:(.*?),(.*?):/)
		{
			$depth = $1 + $2;
			$freq = 0;
			if ($depth > 0)
			{
				$freq = $2 / $depth;
			}
			if ($gt{$tum} =~ /^.\/.:(.*?),(.*?):/)
			{
				if (($depth >= $minDepth) and ($qual >= $minQual) and ($freq >= $minFreq) and ($freq <= $maxFreq))
				{
					print "$chr\t$pos\t$1\t$2\n";
				}
			}
		}
	}


}

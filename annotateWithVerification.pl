#!/usr/bin/perl

use strict;
use warnings;

my $valFile = $ARGV[0];


my $l;
my @f;

my %val;

open (VAL, "<$valFile") or die "Couldn't open [$valFile]\n";

while ($l = <VAL>)
{
	chomp $l;
	@f = split(/\t/, $l);

	unless (exists $val{"$f[0]\t$f[1]\t$f[3]\t$f[4]"})
	{
		$val{"$f[0]\t$f[1]\t$f[3]\t$f[4]"} = $f[6];
	}
	else
	{
		unless (($f[6] eq "INSUFFICIENT_COVERAGE") or ($val{"$f[0]\t$f[1]\t$f[3]\t$f[4]"} eq "TRUE-SOMATIC"))
		{
			$val{"$f[0]\t$f[1]\t$f[3]\t$f[4]"} = $f[6];
		}
	}
}

close VAL;


my $verified = 0;
my $notVerified = 0;



while ($l = <STDIN>)
{
	unless ($l =~ /^#/)
	{
		chomp $l;
		@f = split(/\t/, $l);

		if (exists $val{"$f[0]\t$f[1]\t$f[3]\t$f[4]"})
		{
			$verified++;

			if ($f[7] =~ /.*VERIFICATION=/)
			{
				if ($f[7] =~ /.*VERIFICATION=(.*?);.*/)
				{
					if ($1 ne $val{"$f[0]\t$f[1]\t$f[3]\t$f[4]"})
					{
						warn "Verification results inconsistent with what is already in the vcf ($1 ne " . $val{"$f[0]\t$f[1]\t$f[3]\t$f[4]"} . "), not updating\n";
					}
				}
				elsif ($f[7] =~ /.*VERIFICATION=(.*)/)
				{
					if ($1 ne $val{"$f[0]\t$f[1]\t$f[3]\t$f[4]"})
					{
						warn "Verification results inconsistent with what is already in the vcf ($1 ne " . $val{"$f[0]\t$f[1]\t$f[3]\t$f[4]"} . "), not updating\n";
					}
				}
				print $l . "\n";
			}
			else
			{
				$f[7] .= ";VERIFICATION=" . $val{"$f[0]\t$f[1]\t$f[3]\t$f[4]"};
				
				print "$f[0]";
				for (my $i = 1; $i < scalar(@f); $i++)
				{
					print "\t$f[$i]";
				}
				print "\n";
			}

		}
		else
		{
			$notVerified++;
			print $l . "\n";
		}
	}
	else
	{
		print $l;
	}
}

warn "\t$verified verified, $notVerified not verified\n\n";




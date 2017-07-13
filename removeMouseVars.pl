#!/usr/bin/perl

use strict;
use warnings;
use Tabix;

my $mouseList = $ARGV[0];
my $mouseOutput = $ARGV[1];

my $l;
my @f;

my $mouseLine;
my @m;

my $tabix = Tabix->new(-data => $mouseList, -index => "$mouseList.tbi");
my $iter;

my $foundMouse;

open (MOUSEOUT, ">$mouseOutput") or die "Couldn't open >$mouseOutput\n";


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
			$foundMouse = 0;
			while ($mouseLine = $tabix->read($iter))
			{
				@m = split(/\t/, $mouseLine);
				if (($m[0] eq $f[0]) and ($m[1] eq $f[1]) and ($m[3] eq $f[3]) and ($m[4] eq $f[4]))
				{
					$foundMouse = 1;
				}
			}
			if ($foundMouse == 0)
			{
				print $l . "\n";
			}
			else
			{
				print MOUSEOUT $l . "\n";
			}
		}
		else
		{
			print $l . "\n";
		}
	}
}

close MOUSEOUT;



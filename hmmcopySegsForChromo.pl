#!/usr/bin/perl

use strict;
use warnings;


my $l;
my @f;

my $count = 0;

while ($l = <STDIN>)
{
	chomp $l;
	@f = split(/\t/, $l);

	if ($count == 0)
	{
		print "\t$f[0]\t$f[1]\t$f[2]\t$f[4]\n";
	}
	else
	{
		print "$count\t$f[0]\t$f[1]\t$f[2]\t$f[4]\n";
	}
	$count++;

}


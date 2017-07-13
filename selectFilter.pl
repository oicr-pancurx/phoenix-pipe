#!/usr/bin/perl

use strict;
use warnings;

my $l;
my @f;

my $filterTerm = $ARGV[0];

while ($l = <STDIN>)
{
	chomp $l;
	@f = split(/\t/, $l);

	if ($l =~ /^#/)
	{
		print $l . "\n";
	}
	elsif ($f[6] eq $filterTerm)
	{
		print $l . "\n";
	}
}



#!/usr/bin/perl

use strict;
use warnings;

my $l;
my @f;

while ($l = <STDIN>)
{
	chomp $l;

	if ($l =~ /^#/)
	{
		print $l . "\n";
	}
	else
	{
		@f = split(/\t/, $l)
		{

		}
	}
}

#!/usr/bin/perl

use strict;
use warnings;

my $l;

while ($l = <STDIN>)
{
	print "@" . $l;

	$l = <STDIN>;
	print $l;

	$l = <STDIN>;
	print "+" . $l;

	$l = <STDIN>;
	print $l;

}

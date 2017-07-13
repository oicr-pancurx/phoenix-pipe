#!/usr/bin/perl

use strict;
use warnings;

until (-e "gatk_done")
{
	sleep 300;		# check every 5 minutes
}


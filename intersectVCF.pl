#!/usr/bin/perl

use strict;
use warnings;

my $fileA = $ARGV[0];
my $fileB = $ARGV[1];	# prints from this one

my $alignTool = $ARGV[2];
my $callTool = $ARGV[3];

my $alignString;
my $callString;

open (FILEA, "$fileA") or die;
open (FILEB, "$fileB") or die;

my %fileAhash;
my $line;
my @fields;

while ($line = <FILEA>)
{
	unless ($line =~ /^#/)
	{
		chomp $line;
		@fields = split(/\t/, $line);

		unless ($fields[0] =~ /_/)		# don't want no random chromosomes!
		{
			$fileAhash{"$fields[0]\t$fields[1]\t$fields[3]\t$fields[4]"} = $line;
		}
	}
}

my $altDepth;
my $totalDepth;
my $freq;

if ($alignTool eq "bwa/0.6.2")
{
	$alignString = "BWA 0.6.2 - http://bio-bwa.sourceforge.net/";
}
elsif ($alignTool eq "novocraft/3.00.05")
{
	$alignString = "Novoalign 3.00.05 - http://novocraft.com/main/index.php";
}

if ($callTool eq "mutect/1.1.4,strelka/v1.0.7")
{
	$callString = "Strelka 1.0.7 - https://sites.google.com/site/strelkasomaticvariantcaller/,MuTect 1.1.4 - http://www.broadinstitute.org/cancer/cga/mutect";
}
else
{
	die "Not handling $callTool";
}


while ($line = <FILEB>)
{
	chomp $line;
	@fields = split(/\t/, $line);

	if ($line =~ /^#/)
	{
		if ($line =~ /##DCC=<alignment_algorithm/)
		{
			print "##DCC=<alignment_algorithm=\"$alignString\">\n";
		}
		elsif ($line =~ /##DCC=<variation_calling_algorithm/)
		{
			print "##DCC=<variation_calling_algorithm=\"$callString\">\n";
		}
		else
		{
			print $line . "\n";
		}
	}
	elsif (length($fields[3]) ne length($fields[4]))	# print indels
	{
		unless ($fields[0] =~ /_/)		# don't want no random chromosomes!
		{
			print $line . "\n";
		}
	}
	elsif (exists $fileAhash{"$fields[0]\t$fields[1]\t$fields[3]\t$fields[4]"})
	{
		print $line . "\n";
	}
}

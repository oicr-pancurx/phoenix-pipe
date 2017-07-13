#!/usr/bin/perl

use strict;
use warnings;

my $funseqFile = $ARGV[0];

my %funseq;
my @headerLines;

my $l;
my @f;


# pull variant function file into hash
open (FILE, "<$funseqFile") or die "Couldn't open $funseqFile\n";
while ($l = <FILE>)
{
	chomp $l;
	@f = split(/\t/, $l);

	if ($l =~ /^#/)
	{
		unless (($l =~ /^#CHROM/) or ($l =~/^##fileformat/))
		{
			push(@headerLines, $l);
		}
	}
	else
	{
		$f[0] = "chr$f[0]";

		$f[7] =~ s/.*?CDS=/;CDS=/;

		$funseq{"$f[0]\t$f[1]\t$f[3]\t$f[4]"} = $f[7];;
	}
}
close FILE;





# iterate over vcf to check for annotation and add to info column
# add note to header?

while ($l = <STDIN>)
{
	chomp $l;
	@f = split(/\t/, $l);

	if ($l =~ /^#/)
	{
		if ($l =~ /^#CHROM/)
		{
			for my $funseqLine (@headerLines)
			{
				print $funseqLine . "\n";
			}
		}
		print $l . "\n";
	}
	else
	{
		if (exists $funseq{"$f[0]\t$f[1]\t$f[3]\t$f[4]"})
		{
			$f[7] .= $funseq{"$f[0]\t$f[1]\t$f[3]\t$f[4]"};
		}

		print $f[0];
		for (my $i = 1; $i < scalar(@f); $i++)
		{
			print "\t$f[$i]";
		}
		print "\n";
	}


}






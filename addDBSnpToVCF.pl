#!/usr/bin/perl

use strict;
use warnings;
use Tabix;

my @f;
my $l;
my @db;

my $dbSNPline;

my $dbSNP = $ARGV[0];
#my $dbSNP = "/oicr/data/reference/genomes/homo_sapiens_mc/dbSNP/hg19_random/Genomic/dbSNP137/dbSNP137_chr.vcf.gz";

my ($chr, $pos, $ref, $alt, $dbRef, $dbAlt);
my $dbSNPstring;

my $tabix = Tabix->new(-data => $dbSNP, -index => "$dbSNP.tbi");
my $iter;
my $start;
my $end;

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
		$chr = $f[0];
		$pos = $f[1];
		$ref = $f[3];
		$alt = $f[4];

		$start = $pos - 1;
		$end = $pos;

		$iter = $tabix->query($chr,$start,$end);


		if (defined $iter->{"_"})		# happens if the contig isn't in the bed file
		{
			while ($dbSNPline = $tabix->read($iter))
			{
				@db = split(/\t/, $dbSNPline);
				for my $dbAlt (split(/,/, $db[4]))
				{
					if (($db[3] eq $f[3]) and ($dbAlt eq $f[4]))
					{
						$f[2] = $db[2];
						
                        if (length($f[3]) == length($f[4]))             # snv
                        {
                                $f[7] .= ";DBSNP_ALLELE=$f[3]/$f[4]";
                                $f[7] .= ";DBSNP_STRAND=1";
                        }
                        elsif (length($f[3]) > length($f[4]))           # deletion
                        {
                                $dbRef = substr($f[3],1);         # remove the anchor base
                                $dbAlt = "-" x length($dbRef);
                                $f[7] .= ";DBSNP_ALLELE=$dbRef/$dbAlt";
								$f[7] .= ";DBSNP_STRAND=1";
                        }
                        elsif (length($f[3]) < length($f[4]))           # insertion
                        {
                                $dbAlt = substr($f[4],1);         # remove the anchor base
                                $dbRef = "-" x length($dbAlt);
								$f[7] .= ";DBSNP_ALLELE=$dbRef/$dbAlt";
								$f[7] .= ";DBSNP_STRAND=1";
                        }

						for my $info (split(/;/, $db[7]))
						{
							if ($info =~ /GMAF=.*/)
							{
								$f[7] .= ";DBSNP_$info";
							}
							elsif ($info =~ /G5A/)
							{
								$f[7] .= ";DBSNP_$info";
							}
							elsif ($info =~ /G5/)
							{
								$f[7] .= ";DBSNP_$info";
							}
						}
					}
				}
			}
		}
	
		print $f[0];
		for (my $i = 1; $i < scalar (@f); $i++)
		{
			print "\t$f[$i]";
		}
		print "\n";
	}

	}

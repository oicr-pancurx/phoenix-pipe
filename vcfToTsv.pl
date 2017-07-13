#!/usr/bin/perl

use strict;
use warnings;

my $qualCut = 50;

my $l;
my @f;
my @info;

my $sample = "";
my $seqType = "";
my $varType = "";
my $prep = "";

my $tSample;
my $nSample;

my %sampleCount;
my %callCount;
my %metadata;

$sampleCount{"Somatic exome"} = 0;
$sampleCount{"Somatic wgs"} = 0;
$sampleCount{"Germline exome"} = 0;
$sampleCount{"Germline wgs"} = 0;

my $freq;
my $depth;
my $prediction;
my $reason;
my $dbSNP;
my $common;
my $cosmic;
my $tracks;
my $gene;
my $consequence;
my $annovar;
my ($refBases,$altBases,$totBases);

my %allCalls;
my @sampOrder;

my $outputPrefix = $ARGV[0];

my $delOutput = "$outputPrefix-deleterious.tsv";
my $allOutput = "$outputPrefix-all.tsv";
my $matrixOutput = "$outputPrefix-matrix.tsv";

open (DEL, ">$delOutput") or die "Couldn't open $delOutput\n";
open (ALL, ">$allOutput") or die "Couldn't open $allOutput\n";
open (MATRIX, ">$matrixOutput") or die "Couldn't open $matrixOutput\n";


print DEL "Sample\tPrep\tSeq Type\tVar Type\tGene\tChromosome\tPosition\tRef\tAlt\tFrequency\tDepth\tdbSNP\tCOSMIC\tConsequence\tAnnovar Annotation(s)\n";
print ALL "Sample\tPrep\tSeq Type\tVar Type\tGene\tChromosome\tPosition\tRef\tAlt\tFrequency\tDepth\tdbSNP\tCOSMIC\tConsequence\tAnnovar Annotation(s)\n";

my $tCol;

while ($l = <STDIN>)
{
	if ($l =~ /analyzed_sample_id="(.*)"/)
	{
		$tSample = $1;
	}
	elsif ($l =~ /matched_sample_id="(.*)"/)
	{
		$nSample = $1;
		print DEL "\n";
		print ALL "\n";
	}
	elsif ($l =~ /variation_calling_algorithm="(.*?) /)
	{
		if ($1 eq "Strelka")
		{
			$varType = "Somatic";
			$sample = $tSample;
		}
		elsif ($1 eq "GATK")
		{
			$varType = "Germline";
			$sample = $nSample;
		}
		else
		{
			$varType = "";
			$sample = "";
		}
		$sampleCount{$varType}++;
	}
	elsif ($l =~ /sequencing_strategy=.,term="(.*)"/)
	{
		if ($1 eq "WXS")
		{
			$seqType = "exome";
		}
		elsif ($1 eq "WGS")
		{
			$seqType = "wgs";
		}
		else
		{
			$seqType = "";
		}
		$sampleCount{all}++;
		$sampleCount{$seqType}++;

		$sampleCount{"$varType $seqType"}++;
	}
	elsif ($l =~ /^#CHROM/)
	{
		chomp $l;
		@f = split(/\t/, $l);
		
		if ($f[9] eq "$sample")
		{
			$tCol = 9;
		}
		else
		{
			$tCol = 10;		# either the sample, or strelka's TUMOR column
		}


		if ($sample =~ /.*526/)
		{
			$prep = "lcm";
		}
		elsif ($sample =~ /ASHPC.*/)
		{
			$prep = "flow";
		}
		else
		{
			$prep = "bulk";
		}

		push (@sampOrder, "$sample $seqType");

		warn "Processing $sample $seqType $prep $varType\n";
	}
	elsif ($l =~ /^#/)
	{
	}
	else
	{
		chomp $l;
		@f = split(/\t/, $l);

		if ($f[5] eq ".")
		{
			$f[5] = $qualCut;		# if no quals, take em all
		}

		if ($f[5] >= $qualCut)
		{
	
			$freq = "";
			$depth = "";
			$dbSNP = $f[2];
			$common = ".";
			$cosmic = ".";
			$tracks = ".";
			$gene = ".";
			$consequence = ".";
			$prediction = ".";

	
			@info = split(/;/, $f[7]);
	
			for my $i (@info)
			{
				if ($i =~ /COSMIC=(.*)/)
				{
					$cosmic = $1;
				}
	
				elsif ($i =~ /TRACK=(.*)/)
				{
					if ($tracks eq ".")
					{
						$tracks = "$1";
					}
					else
					{
						$tracks .= ",$1";
					}
				}
	
				elsif ($i =~ /DBSNP_G5$/)
				{
					if ($common eq ".")
					{
						$common = "G5";
					}
					else
					{
						$common .= ",G5";
					}
				}
	
				elsif ($i =~ /DBSNP_G5A/)
				{
					if ($common eq ".")
					{
						$common = "G5A";
					}
					else
					{
						$common .= ",G5A";
					}
				}
	
				elsif ($i =~ /ANNOVAR_EXONIC=(.*?),(.*)/)
				{
					$consequence = $1;
					$annovar = $2;
				}
	
				elsif ($i =~ /ANNOVAR=(.*?),(.*)/)
				{
					if ($consequence eq ".")
					{
						$consequence = $1;
						$annovar = $2;
					}
					$gene = $2;
					$gene =~ s/,.*//;
					if ($gene =~ /dist/)
					{
						$gene = ".";
					}
					$gene =~ s/\(.*//;
	
				}
			}
	
	
	        if ($f[8] eq "GT:AD:DP:GQ:PL")      # GATK format
	        {
				if ($f[$tCol] =~ /\.\/\./)
				{
					$tracks = "skip";
				}
				elsif ($f[$tCol] =~ /^0\/0/)
				{
					$tracks = "skip";
				}
	            elsif ($f[$tCol] =~ /^.*?:(.*?),(.*?):(.*?):.*/)
	            {
	                $refBases = $1;
	                $altBases = $2;
	                $totBases = $3;
	            }
	            else
	            {
	                warn "Assumed GATK format, couldn't parse frequencies from $f[$tCol]\n";

					$refBases = 0;
					$altBases = 0;
					$totBases = 0;
	            }
	        }
	        elsif ($f[8] eq "DP:FDP:SDP:SUBDP:AU:CU:GU:TU")     # strelka snv format
	        {
	            if ($f[10] =~ /^(.*?):.*:(.*?),.*?:(.*?),.*?:(.*?),.*?:(.*?),.*?$/)
	            {
	                $totBases = $1;
	
	                if ($f[4] eq "A")
	                {
	                    $altBases = $2;
	                }
	                elsif ($f[4] eq "C")
	                {
	                    $altBases = $3;
	                }
	                elsif ($f[4] eq "G")
	                {
	                    $altBases = $4;
	                }
	                elsif ($f[4] eq "T")
	                {
	                    $altBases = $5;
	                }
	            }
	            else
	            {
	                warn "Assumed Strelka format, couldn't parse frequencies from $f[10]\n";
					$refBases = 0;
					$altBases = 0;
					$totBases = 0;
	            }
	        }
			elsif ($f[8] eq "DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50")		# strelka indel format
			{
				if ($f[10] =~ /^(.*?):.*?:.*?,.*?:(.*?),.*/)
				{
					$totBases = $1;
					$altBases = $2;
				}
				else
				{
					warn "Assumed Strelka indel format, couldn't parse frequencies from $f[10]\n";
					$refBases = 0;
					$altBases = 0;
					$totBases = 0;
				}
			}
			
			unless ($tracks eq "skip")
			{
				$depth = $totBases;
				$freq = 0;
				if ($totBases > 0)
				{
					$freq = sprintf("%.2f", $altBases / $totBases * 100);
				}
			}
			else
			{
				$depth = 0;
				$freq = 0;
			}

			$reason = "";

			if ($freq > 95)
			{
				$reason .= ",freq";
			}
			if ($common ne ".")
			{
				$reason .= ",$common";
			}
			$reason =~ s/^,//;
	
#			if (($freq > 95) or ($common ne ".") or (exists $panelOfNormals{"$f[0]\t$f[1]\t$f[3]\t$f[4]"}))
#			{
#				$prediction = "Likely Germline ($reason)";
#				$prediction = "Likely Germline";
#			}
#			elsif ($cosmic ne ".")
#			{
#				$prediction = "Likely Somatic";
#			}
#			else
#			{
#				$prediction = "Potentially Somatic";
#			}
			

			if (($consequence =~ /nonsynonymous/) or ($consequence =~ /stopgain/) or ($consequence =~ /splicing/) or ($consequence =~ /^frameshift/))
			{
				print DEL "$sample\t$prep\t$seqType\t$varType\t$gene\t$f[0]\t$f[1]\t$f[3]\t$f[4]\t";
				print DEL "$freq%\t$depth\t$dbSNP\t$cosmic\t$consequence\t$annovar";
				print DEL "\n";

				$allCalls{"$gene\t$f[0]\t$f[1]\t$f[3]\t$f[4]"}{"$sample $seqType"} = "$freq% ($depth)";
			
				unless (exists $metadata{"$gene\t$f[0]\t$f[1]\t$f[3]\t$f[4]"})
				{
					$metadata{"$gene\t$f[0]\t$f[1]\t$f[3]\t$f[4]"} = "\t$dbSNP\t$cosmic\t$consequence\t$annovar";
				}
				unless ($metadata{"$gene\t$f[0]\t$f[1]\t$f[3]\t$f[4]"} =~ /COSM/)		# don't want to overwrite cosmic calls - should really annotate the base calls
				{
					$metadata{"$gene\t$f[0]\t$f[1]\t$f[3]\t$f[4]"} = "\t$dbSNP\t$cosmic\t$consequence\t$annovar";
				}

				$callCount{"$gene\t$f[0]\t$f[1]\t$f[3]\t$f[4]"}{"$varType $seqType"}{$sample} = 1;		# don't want to count references more than once
			}

			print ALL "$sample\t$prep\t$seqType\t$varType\t$gene\t$f[0]\t$f[1]\t$f[3]\t$f[4]\t";
			print ALL "$freq%\t$depth\t$dbSNP\t$cosmic\t$consequence\t$annovar";
			print ALL "\n";
	
		}
	}
}

if (1 == 1)
{
	print MATRIX "Gene\tChromosome\tPosition\tRef\tAlt\tdbSNP\tCOSMIC\tConsequence\tAnnovar Annotation\tSomatic Exome Calls ($sampleCount{'Somatic exome'} total)\tSomatic WGS Calls ($sampleCount{'Somatic wgs'} total)\tGermline Exome Calls ($sampleCount{'Germline exome'} total)\tGermline WGS Calls ($sampleCount{'Germline wgs'} total)";
	for my $samp (sort @sampOrder)
	{
		print MATRIX "\t$samp";
	}
	print MATRIX "\n";
	my $timesCalled;
	for my $call (sort keys %allCalls)
	{
		print MATRIX $call;
		print MATRIX $metadata{$call};


		print MATRIX "\t" . scalar(keys %{ $callCount{$call}{"Somatic exome"} });
		print MATRIX "\t" . scalar(keys %{ $callCount{$call}{"Somatic wgs"} });
		print MATRIX "\t" . scalar(keys %{ $callCount{$call}{"Germline exome"} });
		print MATRIX "\t" . scalar(keys %{ $callCount{$call}{"Germline wgs"} });

		for my $samp (sort @sampOrder)
		{
			if (exists $allCalls{$call}{$samp})
			{
				print MATRIX "\t$allCalls{$call}{$samp}";
			}
			else
			{
				print MATRIX "\t";
			}
		}

		print MATRIX "\n";
	}
}


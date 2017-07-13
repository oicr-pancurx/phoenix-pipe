#!/usr/bin/perl

use strict;
use warnings;

use Bio::DB::Fasta;

my $samplePath = $ARGV[0];
my $sample = "";
if (($samplePath =~ /^.*\/(.*?)\/exome/) or ($samplePath =~ /^.*\/(.*?)\/wgs/))
{
	$sample = $1;
}
else
{
	die "Couldn't parse sample from $samplePath\n";
}

my %dataPath = (
	"somatic" => "/bwa/0.6.2/final_strelka-mutect/$sample.final.vcf",
	"germline" => "/bwa/0.6.2/gatk/1.3.16/final/$sample.filtered.final.vcf",
	"cn_points" => "/bwa/0.6.2/HMMcopy/0.1.1/*.segments_somatic.seg",
	"cn_states" => "/bwa/0.6.2/HMMcopy/0.1.1/*.cnv_somatic_segments",
	"structural" => "/bwa/0.6.2/crest/alpha/$sample.predSV.txt"
);

my $outputFile = "$sample.onePage.png";





#my %colours = (
#	"C>A" => "#a6cee3",
#	"C>G" => "#1f78b4",
#	"C>T" => "#b2df8a",	# transversion
#	"T>A" => "#33a02c",
#	"T>C" => "#fb9a99",	# transversion
#	"T>G" => "#e31a1c",
#	"del" => "#984ea3",
#	"ins" => "#ff7f00",
#);

my %colours = (
	"C>A" => "#332288",
	"C>G" => "#88CCEE",
	"C>T" => "#117733",	# transversion
	"T>A" => "#DDCC77",
	"T>C" => "#CC6677",	# transversion
	"T>G" => "#AA4499",
	"del" => "#44AA99",
	"ins" => "#882255",
);

my %cnColours = (
	"1" => "#3A89C9",
	"2" => "#99C7EC",
	"3" => "#E6F5FE",
	"4" => "#FFE3AA",
	"5" => "#F5A275",
	"6" => "#D24D3E"
);

my %bpColours = (
	"CTX" => "#4477AA",
	"DEL" => "#117733",
	"INS" => "#DDCC77",
	"INV" => "#CC6677",
	"ITX" => "#AA4499",
);


my %genesOfInterest = ("AKT3",1,"PARP1",1,"NRAS",1,"NOTCH2",1,"NOTCH2NL",1,"PIK3C2B",1,"PTCH2",1,"PIK3CD",1,"PIK3R3",1,"DPYD",1,"GSTM1",1,"EPHA2",1,"MTHFR",1,"PLK3",1,"CDC7",1,"WNT2B",1,"WNT3A",1,"WNT4",1,"WNT9A",1,"ARID1A",1,"MCL1",1,"PTEN",1,"RET",1,"MGMT",1,"FGFR2",1,"CYP2C19",1,"CYP2C9",1,"CYP2C8",1,"CYP2E1",1,"ABCC2",1,"WNT8B",1,"TET1",1,"CYP17A1",1,"CHEK1",1,"HRAS",1,"PIK3C2A",1,"CCND1",1,"ATM",1,"WEE1",1,"RRM1",1,"GSTP1",1,"WNT11",1,"CBL",1,"PGR",1,"SLC22A6",1,"SLCO2B1",1,"WT1",1,"KRAS",1,"ErbB3",1,"MDM2",1,"CDK4",1,"PIK3C2G",1,"WNT1",1,"PTPN11",1,"PXN",1,"GLI1",1,"WNT10B",1,"WNT5B",1,"SLCO1B1",1,"SLCO1B3",1,"SOCS2",1,"BRCA2",1,"FLT1",1,"FLT3",1,"RB1",1,"ERCC5",1,"AKT1",1,"PARP2",1,"ESR2",1,"TEP1",1,"PGF",1,"IGF1R",1,"MAP2K1",1,"IDH2",1,"CYP1A1",1,"DLL4",1,"CYP1A2",1,"CDH1",1,"TSC2",1,"PLK1",1,"TUBB3",1,"ERCC4",1,"SULT1A1",1,"SOCS1",1,"SH2B",1,"BRCA1",1,"ErbB2",1,"TP53",1,"RARA",1,"TOP2A",1,"AURKB",1,"PIK3R5",1,"NF1",1,"WNT3",1,"WNT9B",1,"SOCS3",1,"STAT3",1,"BCL2",1,"DCC",1,"PIK3C3",1,"TYMS",1,"SMAD4",1,"ROCK1",1,"AKT2",1,"MAP2K2",1,"NOTCH3",1,"XRCC1",1,"AURKC",1,"ERCC1",1,"PIK3R2",1,"ERCC2",1,"STK11",1,"CYP2A6",1,"CYP2B6",1,"AXL",1,"ALK",1,"MSH2",1,"ErbB4",1,"IDH1",1,"UGT1A1",1,"GLI2",1,"ERCC3",1,"WNT10A",1,"WNT6",1,"SRC",1,"AURKA",1,"TOP1",1,"GART",1,"BCR",1,"CHEK2",1,"CYP2D6",1,"GSTT1",1,"NF2",1,"EWSR1",1,"WNT7B",1,"PIK3CA",1,"MLH1",1,"RAF1",1,"CTNNB1",1,"VHL",1,"PIK3CB",1,"FANCD2",1,"TOP2B",1,"ATR",1,"MST1R",1,"TERC",1,"RASSF1",1,"WNT5A",1,"WNT7A",1,"SLC15A2",1,"KIT",1,"KDR",1,"PDGFRA",1,"FGFR3",1,"EIF4E",1,"PLK4",1,"ABCG2",1,"UGT2B15",1,"UGT2B17",1,"UGT2B7",1,"TET2",1,"PDGFRB",1,"DHFR",1,"FLT4",1,"PIK3R1",1,"TERT",1,"PLK2",1,"WNT8A",1,"NPM1",1,"ESR1",1,"NOTCH4",1,"CCND3",1,"TUBB",1,"SLC29A1",1,"TPMT",1,"SLC22A1",1,"SLC22A2",1,"EGFR",1,"BRAF",1,"MET",1,"SHH",1,"SMO",1,"PIK3CG",1,"CDK5",1,"CYP3A4",1,"EPHB4",1,"GLI3",1,"CYP3A5",1,"WNT16",1,"WNT2",1,"ABCB1",1,"FGFR1",1,"PTK2",1,"PTK2B",1,"ANGPT1",1,"ANGPT2",1,"NAT1",1,"NAT2",1,"ABL1",1,"JAK2",1,"CDKN2A",1,"NOTCH1",1,"PTCH1",1,"TSC1",1,"ROR2",1,"AR",1,"ARAF",1);

my %chrTopOffset = (
	"chr1" => 0,
	"chr2" => 249250621,
	"chr3" => 492449994,
	"chr4" => 690472424,
	"chr5" => 881626700,
	"chr6" => 1062541960,
	"chr7" => 1233657027,
	"chr8" => 1392795690,
	"end" => 1539159712,
);
my %chrBotOffset =(
	"chr9" => 0,
	"chr10" => 141213431,
	"chr11" => 276748178,
	"chr12" => 411754694,
	"chr13" => 545606589,
	"chr14" => 660776467,
	"chr15" => 768126007,
	"chr16" => 870657399,
	"chr17" => 961012152,
	"chr18" => 1042207362,
	"chr19" => 1120284610,
	"chr20" => 1179413593,
	"chr21" => 1242439113,
	"chr22" => 1290569008,
	"chrX" => 1341873574,
	"chrY" => 1497144134,
	"end" => 1556517700,
);

my $l;
my $lastL;
my @f;

my %count = (
	"snv" => 0,
	"ins" => 0,
	"del" => 0,
	"nonsyn" => 0,
	"frameins" => 0,
	"framedel" => 0,
	"tv" => 0,
	"kataegis" => 0
);
my %vcfHash;
my @vcfOrder;

my $tSample;
my $nSample;

my $tCoverage;
my $nCoverage;

my $tCol;
my $nCol;

my $date = `date`;
chomp $date;

my ($chr, $pos, $id, $ref, $alt, $freq, $depth, $consequence, $gene);
my ($refBases, $altBases, $totBases);
my $info;
my $type;
my $change;
my $context;

my $lastChrSNV = "null";
my $lastPosSNV = "null";
my $prevDist;

my %comp = (
	"A>C" => "T>G",
	"A>G" => "T>C",
	"A>T" => "T>A",
	"G>A" => "C>T",
	"G>C" => "C>G",
	"G>T" => "C>A",
);

# make gene of interest ticks while reading somatic vcf
my (@genePosTop, @genePosBot, @geneNameTop, @geneNameBot);

open (FILE, "<$samplePath/$dataPath{somatic}") or die "Couldn't open $samplePath/$dataPath{somatic}\n";
while ($l = <FILE>)
{
	if ($l =~ /##DCC=<analyzed_sample_id="(.*?)">/)
	{
		$tSample = $1;
	}
	elsif ($l =~ /##DCC=<matched_sample_id="(.*?)">/)
	{
		$nSample = $1;
	}
	elsif ($l =~ /##DCC=<analyzed_seq_coverage=(.*?)>/)
	{
		$tCoverage = sprintf("%0.2f", $1);
	}
	elsif ($l =~ /##DCC=<matched_seq_coverage=(.*?)>/)
	{
		$nCoverage = sprintf("%0.2f", $1);
	}
	elsif ($l =~ /^#CHROM/)
	{
		chomp $l;
		@f = split(/\t/, $l);
		if ($f[9] eq $tSample)
		{
			$tCol = 9;
			$nCol = 10;
		}
		elsif ($f[10] eq $tSample)
		{
			$tCol = 10;
			$nCol = 9;
		}
		elsif ($f[10] eq "TUMOR")
		{
			$tCol = 10;
			$nCol = 9;
		}
		else
		{
			die "Couldn't detect genotype column!\n";
		}
	}

	unless ($l =~ /^#/)
	{
		chomp $l;
		@f = split(/\t/, $l);

		$chr = $f[0];
		$pos = $f[1];
		$id = $f[2];
		$ref = $f[3];
		$info = $f[7];
		for $alt (split(/,/, $f[4]))
		{
			push(@vcfOrder, "$chr:$pos");
			if (length($ref) == length($alt))
			{
				$type = "snv";
			
				$change = "$ref>$alt";
				if (exists $comp{$change})
				{
					$change = $comp{$change};
				}

				$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{change} = $change;

				$count{snv}++;


				if (($change eq "C>T") or ($change eq "T>C"))
				{
					$count{tv}++;
				}

				if ($chr eq $lastChrSNV)
				{
					$prevDist = $pos - $lastPosSNV;		# assuming sorted input!
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{prevdist} = $prevDist;
				}


				$lastChrSNV = $chr;
				$lastPosSNV = $pos;

			}
			elsif (length($ref) > length($alt))
			{
				$type = "del";

				$count{del}++;

				$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{size} = length($ref) - length($alt);
			}
			else
			{
				$type = "ins";

				$count{ins}++;

				$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{size} = length($alt) - length($ref);
			}

			$gene = "";
			$consequence = "";

			for $info (split(/;/, $f[7]))
			{
				if ($info =~ /COSMIC=(.*)/)
				{
					if ($id eq ".")
					{
						$id = $1;
					}
					else
					{
						$id .= ",$1";
					}
				}
				elsif ($info =~ /ANNOVAR=exonic,(.*)/)
				{
					$gene = $1;
				}
				elsif ($info =~ /ANNOVAR=splicing,(.*)\((.*)\)/)
				{
					$gene = $1;
					$consequence =~ "splicing,$1($2)";
				}
				elsif ($info =~ /ANNOVAR_EXONIC=(.*)/)
				{
					$consequence = $1;
					$consequence =~ s/$gene/ $gene/g;
				}
			}

			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{id} = $id;
			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{gene} = $gene;
			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} = $consequence;

			if ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^frameshift/)
			{
				$count{"frame$type"}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";

					if (exists $chrTopOffset{$chr})
					{
						push(@genePosTop, $pos + $chrTopOffset{$chr});
						push(@geneNameTop, "fs:$gene");
					}
					elsif (exists $chrBotOffset{$chr})
					{
						push(@genePosBot, $pos + $chrBotOffset{$chr});
						push(@geneNameBot, "fs:$gene");
					}
				}
			}
			elsif ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^nonsynonymous/)
			{
				$count{nonsyn}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";

					if (exists $chrTopOffset{$chr})
					{
						push(@genePosTop, $pos + $chrTopOffset{$chr});
						push(@geneNameTop, "ns:$gene");
					}
					elsif (exists $chrBotOffset{$chr})
					{
						push(@genePosBot, $pos + $chrBotOffset{$chr});
						push(@geneNameBot, "ns:$gene");
					}
				}
			}
			elsif ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^stopgain/)
			{
				$count{nonsyn}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";

					if (exists $chrTopOffset{$chr})
					{
						push(@genePosTop, $pos + $chrTopOffset{$chr});
						push(@geneNameTop, "sg:$gene");
					}
					elsif (exists $chrBotOffset{$chr})
					{
						push(@genePosBot, $pos + $chrBotOffset{$chr});
						push(@geneNameBot, "sg:$gene");
					}
				}
			}
			elsif ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^stoploss/)
			{
				$count{nonsyn}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";

					if (exists $chrTopOffset{$chr})
					{
						push(@genePosTop, $pos + $chrTopOffset{$chr});
						push(@geneNameTop, "sl:$gene");
					}
					elsif (exists $chrBotOffset{$chr})
					{
						push(@genePosBot, $pos + $chrBotOffset{$chr});
						push(@geneNameBot, "sl:$gene");
					}
				}
			}
			elsif ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^splicing/)
			{
				$count{"splice$type"}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";

					if (exists $chrTopOffset{$chr})
					{
						push(@genePosTop, $pos + $chrTopOffset{$chr});
						push(@geneNameTop, "sp:$gene");
					}
					elsif (exists $chrBotOffset{$chr})
					{
						push(@genePosBot, $pos + $chrBotOffset{$chr});
						push(@geneNameBot, "sp:$gene");
					}
				}
			}

			$freq = "";
			$depth = "";

	        if ($f[8] eq "GT:AD:DP:GQ:PL")      # GATK format
	        {
	            if ($f[$tCol] =~ /^.*?:(.*?),(.*?):(.*?):.*/)
	            {
	                $refBases = $1;
	                $altBases = $2;
	                $totBases = $3;
	
					$depth = $totBases;
					$freq = $altBases / $depth;
	            }
	            else
	            {
	                die "Assumed GATK format, couldn't parse frequencies from $f[$tCol]\n";
	            }
	        }
	        elsif ($f[8] eq "DP:FDP:SDP:SUBDP:AU:CU:GU:TU")     # strelka format
	        {
	            if ($f[$tCol] =~ /^(.*?):.*:(.*?),.*?:(.*?),.*?:(.*?),.*?:(.*?),.*?$/)
	            {
	                $totBases = $1;
	
	                if ($alt eq "A")
	                {
	                    $altBases = $2;
	                }
	                elsif ($alt eq "C")
	                {
	                    $altBases = $3;
	                }
	                elsif ($alt eq "G")
	                {
	                    $altBases = $4;
	                }
	                elsif ($alt eq "T")
	                {
	                    $altBases = $5;
	                }
					
					$depth = $totBases;
					$freq = $altBases / $depth;
	            }
	            else
	            {
	                die "Assumed Strelka format, couldn't parse frequencies from $f[$tCol]\n";
	            }
	        }
			elsif ($f[8] eq "DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50")
			{
				if ($f[$tCol] =~ /^(.*?):.*?:.*?,.*?:(.*?),.*/)
				{
					$totBases = $1;
					$altBases = $2;
	
					$depth = $totBases;
					$freq = $altBases / $depth;
				}
	            else
	            {
	                die "Assumed Strelka format, couldn't parse frequencies from $f[$tCol]\n";
	            }
			}
	
	
			unless ($freq eq "")
			{
				$freq = sprintf("%.4f", $freq);
			}

			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{tfreq} = $freq;
			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{tdepth} = $depth;



			$freq = "";
			$depth = "";

	        if ($f[8] eq "GT:AD:DP:GQ:PL")      # GATK format
	        {
	            if ($f[$nCol] =~ /^.*?:(.*?),(.*?):(.*?):.*/)
	            {
	                $refBases = $1;
	                $altBases = $2;
	                $totBases = $3;
	
					$depth = $totBases;
					$freq = $altBases / $depth;
	            }
	            else
	            {
	                die "Assumed GATK format, couldn't parse frequencies from $f[$nCol]\n";
	            }
	        }
	        elsif ($f[8] eq "DP:FDP:SDP:SUBDP:AU:CU:GU:TU")     # strelka format
	        {
	            if ($f[$nCol] =~ /^(.*?):.*:(.*?),.*?:(.*?),.*?:(.*?),.*?:(.*?),.*?$/)
	            {
	                $totBases = $1;
	
	                if ($alt eq "A")
	                {
	                    $altBases = $2;
	                }
	                elsif ($alt eq "C")
	                {
	                    $altBases = $3;
	                }
	                elsif ($alt eq "G")
	                {
	                    $altBases = $4;
	                }
	                elsif ($alt eq "T")
	                {
	                    $altBases = $5;
	                }
					
					$depth = $totBases;
					$freq = $altBases / $depth;
	            }
	            else
	            {
	                die "Assumed Strelka format, couldn't parse frequencies from $f[$nCol]\n";
	            }
	        }
			elsif ($f[8] eq "DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50")
			{
				if ($f[$nCol] =~ /^(.*?):.*?:.*?,.*?:(.*?),.*/)
				{
					$totBases = $1;
					$altBases = $2;
	
					$depth = $totBases;
					$freq = 0;
					if ($depth > 0)
					{
						$freq = $altBases / $depth;
					}
				}
	            else
	            {
	                die "Assumed Strelka format, couldn't parse frequencies from $f[$nCol]\n";
	            }
			}
	
	
			unless ($freq eq "")
			{
				$freq = sprintf("%.4f", $freq);
			}

			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{nfreq} = $freq;
			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{ndepth} = $depth;

		}
	}
}
close FILE;


my $titvRatio = "NA";
if ($count{tv} > 0)
{
	$titvRatio = sprintf("%0.4f", ($count{snv} - $count{tv}) / $count{tv});
}







# grab het freqs from GATK

my @hetFreqsTop;
my @hetFreqsBot;
my @hetPosTop;
my @hetPosBot;

my %hetLowAvg;
my %hetHighAvg;
my $hetWindow = 10000;

my $germN;
my $germT;

my $freq1;
my $freq2;

open (FILE, "<$samplePath/$dataPath{germline}") or die "Couldn't open $samplePath/$dataPath{germline}\n";

while ($l = <FILE>)
{
	if ($l =~ /##DCC=<analyzed_sample_id="(.*)">/)
	{
		$germT = $1;
	}
	elsif ($l =~ /##DCC=<matched_sample_id="(.*)">/)
	{
		$germN = $1;
	}
	elsif ($l =~ /#CHROM/)
	{
		chomp $l;
		@f = split(/\t/, $l);

		if ($f[9] eq $germT)
		{
			$germT = 9;
			$germN = 10;
		}
		elsif ($f[9] eq $germN)
		{
			$germN = 9;
			$germT = 10;
		}
		else
		{
			die "Couldn't determine tumour/normal genotype columns in $samplePath/$dataPath{germline}\n";
		}
	}
	elsif ($l =~ /^#/)
	{
	}
	else
	{
		chomp $l;

		@f = split(/\t/, $l);

		if (length($f[3]) == length($f[4]))
		{
			if ($f[$germN] =~ /^0\/1:/)
			{
				unless ($f[$germT] =~ /^\.\/\.:/)
				{
					if ($f[$germT] =~ m/^.*?:(.*?),(.*?):(.*?):/)
					{
						$freq1 = $1/$3;
						$freq2 = $2/$3;

						if ($freq1 < $freq2)
						{
							$hetLowAvg{$f[0]}{int($f[1] / $hetWindow) * $hetWindow + ($hetWindow / 2)} = $freq1;
							$hetHighAvg{$f[0]}{int($f[1] / $hetWindow) * $hetWindow + ($hetWindow / 2)}  = $freq2;
#							push(@{ $hetLowAvg{$f[0]}{int($f[1] / $hetWindow) * $hetWindow + ($hetWindow / 2)} }, $freq1);
#							push(@{ $hetHighAvg{$f[0]}{int($f[1] / $hetWindow) * $hetWindow + ($hetWindow / 2)} }, $freq2);
						}
						else
						{
							$hetLowAvg{$f[0]}{int($f[1] / $hetWindow) * $hetWindow + ($hetWindow / 2)} = $freq2;
							$hetHighAvg{$f[0]}{int($f[1] / $hetWindow) * $hetWindow + ($hetWindow / 2)}  = $freq1;
#							push(@{ $hetLowAvg{$f[0]}{int($f[1] / $hetWindow) * $hetWindow + ($hetWindow / 2)} }, $freq2);
#							push(@{ $hetHighAvg{$f[0]}{int($f[1] / $hetWindow)* $hetWindow + ($hetWindow / 2)} }, $freq1);
						}
					}
				}
			}
		}
	}
}

close FILE;

for $chr (keys %hetLowAvg)
{
	if (exists $chrTopOffset{$chr})
	{
		for $pos (keys %{ $hetLowAvg{$chr} })
		{
			push(@hetFreqsTop, $hetLowAvg{$chr}{$pos});
#			push(@hetFreqsTop, mean($hetLowAvg{$chr}{$pos}));
			push(@hetPosTop, $pos + $chrTopOffset{$chr});
		}
		for $pos (keys %{ $hetHighAvg{$chr} })
		{
			push(@hetFreqsTop, $hetHighAvg{$chr}{$pos});
#			push(@hetFreqsTop, mean($hetHighAvg{$chr}{$pos}));
			push(@hetPosTop, $pos + $chrTopOffset{$chr});
		}
	}
	elsif (exists $chrBotOffset{$chr})
	{
		for $pos (keys %{ $hetLowAvg{$chr} })
		{
			push(@hetFreqsBot, $hetLowAvg{$chr}{$pos});
#			push(@hetFreqsBot, mean($hetLowAvg{$chr}{$pos}));
			push(@hetPosBot, $pos + $chrBotOffset{$chr});
		}
		for $pos (keys %{ $hetHighAvg{$chr} })
		{
			push(@hetFreqsBot, $hetHighAvg{$chr}{$pos});
#			push(@hetFreqsBot, mean($hetHighAvg{$chr}{$pos}));
			push(@hetPosBot, $pos + $chrBotOffset{$chr});
		}
	}
}





# pull cn from HMMcopy
my (@cnStateTop, @cnStartTop, @cnEndTop, @cnStateBot, @cnStartBot, @cnEndBot);
my (@cnPosTop, @cnValTop, @cnPosBot, @cnValBot, @cnColTop, @cnColBot);
my %cnPoints;

my %cnStateType;

my $file = `ls $samplePath/$dataPath{cn_states}`;
chomp $file;

open (FILE, $file) or die "Couldn't open $file\n";

while ($l = <FILE>)
{
	unless ($l =~ /^chr\t/)
	{
		chomp $l;
		@f = split(/\t/, $l);

		if (exists $chrTopOffset{$f[0]})
		{
			push (@cnStartTop, $f[1] + $chrTopOffset{$f[0]});
			push (@cnEndTop, $f[2] + $chrTopOffset{$f[0]});
			push (@cnStateTop, $f[4]);
		}
		elsif (exists $chrBotOffset{$f[0]})
		{
			push (@cnStartBot, $f[1] + $chrBotOffset{$f[0]});
			push (@cnEndBot, $f[2] + $chrBotOffset{$f[0]});
			push (@cnStateBot, $f[4]);
		}

		push(@{ $cnStateType{$f[0]}{pos} }, $f[2]);
		push(@{ $cnStateType{$f[0]}{state} }, $f[3]);

	}
}
close FILE;

$file = `ls $samplePath/$dataPath{cn_points}`;
chomp $file;

open (FILE, $file) or die "Couldn't open $file\n";


while ($l = <FILE>)
{
	unless ($l =~ /^sample\t/)
	{
		chomp $l;
		@f = split(/\t/, $l);

		$chr = $f[1];
		$pos = ($f[2] + $f[3]) / 2;

		push(@{ $cnPoints{$chr}{int($pos / $hetWindow) * $hetWindow + ($hetWindow / 2)} }, $f[4]);
	}
}
close FILE;

my $i;

for $chr (keys %cnPoints)
{
	$i = 0;
	for $pos (sort {$a <=> $b} keys %{ $cnPoints{$chr} })
	{
		if (exists $chrTopOffset{$chr})
		{
			push(@cnPosTop, $pos + $chrTopOffset{$chr});
			push(@cnValTop, mean($cnPoints{$chr}{$pos}));

			while (($pos > $cnStateType{$chr}{pos}[$i]) and ($i < (scalar @{ $cnStateType{$chr}{pos} } - 1)))
			{
				$i++;
			}
			push(@cnColTop, $cnColours{$cnStateType{$chr}{state}[$i]});
		}
		elsif (exists $chrBotOffset{$chr})
		{
			push(@cnPosBot, $pos + $chrBotOffset{$chr});
			push(@cnValBot, mean($cnPoints{$chr}{$pos}));
			while (($pos > $cnStateType{$chr}{pos}[$i]) and ($i < (scalar @{ $cnStateType{$chr}{pos} } - 1)))
			{
				$i++;
			}
			push(@cnColBot, $cnColours{$cnStateType{$chr}{state}[$i]});
		}

	}
}



# get SV from crest
my (@breakX1, @breakX2, @breakY1, @breakY2, @breakCol);
my ($leftChr,$leftPos,$rightChr,$rightPos);
my $bpWidth = 1500000000;

open (FILE, "$samplePath/$dataPath{structural}") or die "Couldn't open $samplePath/$dataPath{structural}\n";

while ($l = <FILE>)
{
	chomp $l;
	@f = split(/\t/, $l);

	$leftChr = $f[0];
	$leftPos = $f[1];
	$rightChr = $f[4];
	$rightPos = $f[5];
	$type = $f[8];

	if ((exists $chrTopOffset{$leftChr}) and (exists $chrTopOffset{$rightChr}))
	{
		$leftPos = (($leftPos + $chrTopOffset{$leftChr}) / $chrTopOffset{end}) * $bpWidth;
		$rightPos = (($rightPos + $chrTopOffset{$rightChr}) / $chrTopOffset{end}) * $bpWidth;

		push(@breakX1,$leftPos);
		push(@breakY1, 1);

		push(@breakX2, ($leftPos + $rightPos) / 2);
		push(@breakY2, 0.5);
		
		push(@breakX1, ($leftPos + $rightPos) / 2);
		push(@breakY1, 0.5);

		push(@breakX2, $rightPos);
		push(@breakY2, 1);

		push(@breakCol,$bpColours{$type});
		push(@breakCol,$bpColours{$type});
	}
	elsif ((exists $chrBotOffset{$leftChr}) and (exists $chrBotOffset{$rightChr}))
	{
		$leftPos = (($leftPos + $chrBotOffset{$leftChr}) / $chrBotOffset{end}) * $bpWidth;
		$rightPos = (($rightPos + $chrBotOffset{$rightChr}) / $chrBotOffset{end}) * $bpWidth;

		push(@breakX1,$leftPos);
		push(@breakY1, 0);

		push(@breakX2, ($leftPos + $rightPos) / 2);
		push(@breakY2, 0.5);
		
		push(@breakX1, ($leftPos + $rightPos) / 2);
		push(@breakY1, 0.5);

		push(@breakX2, $rightPos);
		push(@breakY2, 0);

		push(@breakCol,$bpColours{$type});
		push(@breakCol,$bpColours{$type});
	}
	elsif ((exists $chrTopOffset{$leftChr}) and (exists $chrBotOffset{$rightChr}))
	{
		$leftPos = (($leftPos + $chrTopOffset{$leftChr}) / $chrTopOffset{end}) * $bpWidth;
		$rightPos = (($rightPos + $chrBotOffset{$rightChr}) / $chrBotOffset{end}) * $bpWidth;

		push(@breakX1,$leftPos);
		push(@breakY1, 1);

		push(@breakX2, ($leftPos + $rightPos) / 2);
		push(@breakY2, 0.5);
		
		push(@breakX1, ($leftPos + $rightPos) / 2);
		push(@breakY1, 0.5);

		push(@breakX2, $rightPos);
		push(@breakY2, 0);

		push(@breakCol,$bpColours{$type});
		push(@breakCol,$bpColours{$type});
	}
	elsif ((exists $chrBotOffset{$leftChr}) and (exists $chrTopOffset{$rightChr}))
	{
		$leftPos = (($leftPos + $chrBotOffset{$leftChr}) / $chrBotOffset{end}) * $bpWidth;
		$rightPos = (($rightPos + $chrTopOffset{$rightChr}) / $chrTopOffset{end}) * $bpWidth;

		push(@breakX1,$leftPos);
		push(@breakY1, 0);

		push(@breakX2, ($leftPos + $rightPos) / 2);
		push(@breakY2, 0.5);
		
		push(@breakX1, ($leftPos + $rightPos) / 2);
		push(@breakY1, 0.5);
		
		push(@breakX2, $rightPos);
		push(@breakY2, 1);

		push(@breakCol,$bpColours{$type});
		push(@breakCol,$bpColours{$type});
	}


}
close FILE;


# make rainfall plot
my @positionsTop;
my @positionsBot;

my @distancesTop;
my @distancesBot;

my @coloursTop;
my @coloursBot;

my $rainMaxTop = 0;
my $rainMaxBot = 0;

for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{snv} })
		{
			if (exists $chrTopOffset{$chr})
			{
				if (exists $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist})
				{
					push(@positionsTop, $pos + $chrTopOffset{$chr});
					push(@distancesTop, $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist});
					push(@coloursTop, $colours{$vcfHash{$chr}{$pos}{snv}{$alt}{change}});

					if ($vcfHash{$chr}{$pos}{snv}{$alt}{prevdist} > $rainMaxTop)
					{
						$rainMaxTop = $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist};
					}
				}
			}
			elsif (exists $chrBotOffset{$chr})
			{
				if (exists $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist})
				{
					push(@positionsBot, $pos + $chrBotOffset{$chr});
					push(@distancesBot, $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist});
					push(@coloursBot, $colours{$vcfHash{$chr}{$pos}{snv}{$alt}{change}});

					if ($vcfHash{$chr}{$pos}{snv}{$alt}{prevdist} > $rainMaxBot)
					{
						$rainMaxBot = $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist};
					}
				}
			}
		}
	}
}

my @verticalsTop;
my @verticalsBot;
for $chr (keys %chrTopOffset)
{
	push(@verticalsTop, $chrTopOffset{$chr});
}
for $chr (keys %chrBotOffset)
{
	push(@verticalsBot, $chrBotOffset{$chr});
}

my @chrListTop = qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 end/;
my @atLabelsTop;
my @chrLabelsTop;
for (my $i = 0; $i < (scalar @chrListTop - 1); $i++)
{
	push(@atLabelsTop, ($chrTopOffset{$chrListTop[$i]} + $chrTopOffset{$chrListTop[$i + 1]}) / 2);

	push(@chrLabelsTop, $chrListTop[$i]);
	$chrLabelsTop[$i] =~ s/chr//;
}

my @chrListBot = qw/chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY end/;
my @atLabelsBot;
my @chrLabelsBot;
for (my $i = 0; $i < (scalar @chrListBot - 1); $i++)
{
	push(@atLabelsBot, ($chrBotOffset{$chrListBot[$i]} + $chrBotOffset{$chrListBot[$i + 1]}) / 2);

	push(@chrLabelsBot, $chrListBot[$i]);
	$chrLabelsBot[$i] =~ s/chr//;
}



# read karotypes
$file = "/u/rdenroche/hg19_bands.txt";
open (FILE, $file) or die "Couldn't open $file\n";

my (@karoStartTop, @karoStartBot, @karoEndTop, @karoEndBot, @karoColTop, @karoColBot);

my %karoColour = (
	"acen" => "#CC6677",
	"gpos25" => "grey75",
	"gpos50" => "grey50",
	"gpos75" => "grey25",
	"gpos100" => "grey0"
);

while ($l = <FILE>)
{
	unless ($l =~ /^#/)
	{
		chomp $l;
		@f = split(/\t/, $l);

		if (exists $chrTopOffset{$f[0]})
		{
			push(@karoStartTop, $f[1] + $chrTopOffset{$f[0]});
			push(@karoEndTop, $f[2] + $chrTopOffset{$f[0]});
			if (exists $karoColour{$f[4]})
			{
				push(@karoColTop, $karoColour{$f[4]});
			}
			else
			{
				push(@karoColTop, "white");
			}
		}
		elsif (exists $chrBotOffset{$f[0]})
		{
			push(@karoStartBot, $f[1] + $chrBotOffset{$f[0]});
			push(@karoEndBot, $f[2] + $chrBotOffset{$f[0]});
			if (exists $karoColour{$f[4]})
			{
				push(@karoColBot, $karoColour{$f[4]});
			}
			else
			{
				push(@karoColBot, "white");
			}
		}


	}
}
close FILE;










open (RFILE, ">$outputFile.R") or die "Couldn't open >$outputFile.R\n";


print RFILE "xvalsTop <- c($positionsTop[0]";
for (my $i = 1; $i < scalar @positionsTop; $i++)
{
	print RFILE ", \"$positionsTop[$i]\"";
}
print RFILE ")\n";

print RFILE "yvalsTop <- c($distancesTop[0]";
for (my $i = 1; $i < scalar @distancesTop; $i++)
{
	print RFILE ", \"$distancesTop[$i]\"";
}
print RFILE ")\n";

print RFILE "colsTop <- c(\"$coloursTop[0]\"";
for (my $i = 1; $i < scalar @coloursTop; $i++)
{
	print RFILE ", \"$coloursTop[$i]\"";
}
print RFILE ")\n";

print RFILE "vertsTop <- c($verticalsTop[0]";
for (my $i = 1; $i < scalar @verticalsTop; $i++)
{
	print RFILE ", $verticalsTop[$i]";
}
print RFILE ")\n";

print RFILE "atLabelsTop <- c($atLabelsTop[0]";
for (my $i = 1; $i < scalar @atLabelsTop; $i++)
{
	print RFILE ", $atLabelsTop[$i]";
}
print RFILE ")\n";

print RFILE "chrLabelsTop <- c(\"$chrLabelsTop[0]\"";
for (my $i = 1; $i < scalar @chrLabelsTop; $i++)
{
	print RFILE ", \"$chrLabelsTop[$i]\"";
}
print RFILE ")\n";



print RFILE "xvalsBot <- c($positionsBot[0]";
for (my $i = 1; $i < scalar @positionsBot; $i++)
{
	print RFILE ", \"$positionsBot[$i]\"";
}
print RFILE ")\n";

print RFILE "yvalsBot <- c($distancesBot[0]";
for (my $i = 1; $i < scalar @distancesBot; $i++)
{
	print RFILE ", \"$distancesBot[$i]\"";
}
print RFILE ")\n";

print RFILE "colsBot <- c(\"$coloursBot[0]\"";
for (my $i = 1; $i < scalar @coloursBot; $i++)
{
	print RFILE ", \"$coloursBot[$i]\"";
}
print RFILE ")\n";

print RFILE "vertsBot <- c($verticalsBot[0]";
for (my $i = 1; $i < scalar @verticalsBot; $i++)
{
	print RFILE ", $verticalsBot[$i]";
}
print RFILE ")\n";

print RFILE "atLabelsBot <- c($atLabelsBot[0]";
for (my $i = 1; $i < scalar @atLabelsBot; $i++)
{
	print RFILE ", $atLabelsBot[$i]";
}
print RFILE ")\n";

print RFILE "chrLabelsBot <- c(\"$chrLabelsBot[0]\"";
for (my $i = 1; $i < scalar @chrLabelsBot; $i++)
{
	print RFILE ", \"$chrLabelsBot[$i]\"";
}
print RFILE ")\n";




# het freq plot R values
print RFILE "hetFreqsTop <- c($hetFreqsTop[0]";
for (my $i = 1; $i < scalar @hetFreqsTop; $i++)
{
	print RFILE ", $hetFreqsTop[$i]";
}
print RFILE ")\n";

print RFILE "hetFreqsBot <- c($hetFreqsBot[0]";
for (my $i = 1; $i < scalar @hetFreqsBot; $i++)
{
	print RFILE ", $hetFreqsBot[$i]";
}
print RFILE ")\n";

print RFILE "hetPosTop <- c($hetPosTop[0]";
for (my $i = 1; $i < scalar @hetPosTop; $i++)
{
	print RFILE ", $hetPosTop[$i]";
}
print RFILE ")\n";

print RFILE "hetPosBot <- c($hetPosBot[0]";
for (my $i = 1; $i < scalar @hetPosBot; $i++)
{
	print RFILE ", $hetPosBot[$i]";
}
print RFILE ")\n";



# cn state plot R values

print RFILE "cnStateTop <- c($cnStateTop[0]";
for (my $i = 1; $i < scalar @cnStateTop; $i++)
{
	print RFILE ", $cnStateTop[$i]";
}
print RFILE ")\n";

print RFILE "cnStartTop <- c($cnStartTop[0]";
for (my $i = 1; $i < scalar @cnStartTop; $i++)
{
	print RFILE ", $cnStartTop[$i]";
}
print RFILE ")\n";

print RFILE "cnEndTop <- c($cnEndTop[0]";
for (my $i = 1; $i < scalar @cnEndTop; $i++)
{
	print RFILE ", $cnEndTop[$i]";
}
print RFILE ")\n";

print RFILE "cnStateBot <- c($cnStateBot[0]";
for (my $i = 1; $i < scalar @cnStateBot; $i++)
{
	print RFILE ", $cnStateBot[$i]";
}
print RFILE ")\n";

print RFILE "cnStartBot <- c($cnStartBot[0]";
for (my $i = 1; $i < scalar @cnStartBot; $i++)
{
	print RFILE ", $cnStartBot[$i]";
}
print RFILE ")\n";

print RFILE "cnEndBot <- c($cnEndBot[0]";
for (my $i = 1; $i < scalar @cnEndBot; $i++)
{
	print RFILE ", $cnEndBot[$i]";
}
print RFILE ")\n";


# cn point plot R values

print RFILE "cnPosTop <- c($cnPosTop[0]";
for (my $i = 1; $i < scalar @cnPosTop; $i++)
{
	print RFILE ", $cnPosTop[$i]";
}
print RFILE ")\n";

print RFILE "cnValTop <- c($cnValTop[0]";
for (my $i = 1; $i < scalar @cnValTop; $i++)
{
	print RFILE ", $cnValTop[$i]";
}
print RFILE ")\n";

print RFILE "cnPosBot <- c($cnPosBot[0]";
for (my $i = 1; $i < scalar @cnPosBot; $i++)
{
	print RFILE ", $cnPosBot[$i]";
}
print RFILE ")\n";

print RFILE "cnValBot <- c($cnValBot[0]";
for (my $i = 1; $i < scalar @cnValBot; $i++)
{
	print RFILE ", $cnValBot[$i]";
}
print RFILE ")\n";

print RFILE "cnColTop <- c(\"$cnColTop[0]\"";
for (my $i = 1; $i < scalar @cnColTop; $i++)
{
	print RFILE ", \"$cnColTop[$i]\"";
}
print RFILE ")\n";

print RFILE "cnColBot <- c(\"$cnColBot[0]\"";
for (my $i = 1; $i < scalar @cnColBot; $i++)
{
	print RFILE ", \"$cnColBot[$i]\"";
}
print RFILE ")\n";


# karo bands
#my (@karoStartTop, @karoStartBot, @karoEndTop, @karoEndBot, @karoColTop, @karoColBot);

print RFILE "karoStartTop <- c($karoStartTop[0]";
for (my $i = 1; $i < scalar @karoStartTop; $i++)
{
	print RFILE ", $karoStartTop[$i]";
}
print RFILE ")\n";

print RFILE "karoEndTop <- c($karoEndTop[0]";
for (my $i = 1; $i < scalar @karoEndTop; $i++)
{
	print RFILE ", $karoEndTop[$i]";
}
print RFILE ")\n";

print RFILE "karoColTop <- c(\"$karoColTop[0]\"";
for (my $i = 1; $i < scalar @karoColTop; $i++)
{
	print RFILE ", \"$karoColTop[$i]\"";
}
print RFILE ")\n";

print RFILE "karoStartBot <- c($karoStartBot[0]";
for (my $i = 1; $i < scalar @karoStartBot; $i++)
{
	print RFILE ", $karoStartBot[$i]";
}
print RFILE ")\n";

print RFILE "karoEndBot <- c($karoEndBot[0]";
for (my $i = 1; $i < scalar @karoEndBot; $i++)
{
	print RFILE ", $karoEndBot[$i]";
}
print RFILE ")\n";

print RFILE "karoColBot <- c(\"$karoColBot[0]\"";
for (my $i = 1; $i < scalar @karoColBot; $i++)
{
	print RFILE ", \"$karoColBot[$i]\"";
}
print RFILE ")\n";


my $breakWidth = 12000000;
my @karoBreakTop;
my @karoBreakBot;

for $chr (qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 end/)
{
	push (@karoBreakTop, $chrTopOffset{$chr});
}
for $chr (qw/chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY end/)
{
	push (@karoBreakBot, $chrBotOffset{$chr});
}

# triangles
print RFILE "karoBreakTop <- c($karoBreakTop[0]," . ($karoBreakTop[0] - $breakWidth) . "," . ($karoBreakTop[0] + $breakWidth) . ",NA,$karoBreakTop[0]," . ($karoBreakTop[0] - $breakWidth) . "," . ($karoBreakTop[0] + $breakWidth) . ",NA";
for (my $i = 1; $i < scalar @karoBreakTop; $i++)
{
	print RFILE ",$karoBreakTop[$i]," . ($karoBreakTop[$i] - $breakWidth) . "," . ($karoBreakTop[$i] + $breakWidth) . ",NA,$karoBreakTop[$i]," . ($karoBreakTop[$i] - $breakWidth) . "," . ($karoBreakTop[$i] + $breakWidth) . ",NA";
}
print RFILE ")\n";

print RFILE "karoBreakTopY <- c(0.5,1,1,NA,0.5,0,0,NA";
for (my $i = 1; $i < scalar @karoBreakTop; $i++)
{
	print RFILE ",0.5,1,1,NA,0.5,0,0,NA";
}
print RFILE ")\n";

print RFILE "karoBreakBot <- c($karoBreakBot[0]," . ($karoBreakBot[0] - $breakWidth) . "," . ($karoBreakBot[0] + $breakWidth) . ",NA,$karoBreakBot[0]," . ($karoBreakBot[0] - $breakWidth) . "," . ($karoBreakBot[0] + $breakWidth) . ",NA";
for (my $i = 1; $i < scalar @karoBreakBot; $i++)
{
	print RFILE ",$karoBreakBot[$i]," . ($karoBreakBot[$i] - $breakWidth) . "," . ($karoBreakBot[$i] + $breakWidth) . ",NA,$karoBreakBot[$i]," . ($karoBreakBot[$i] - $breakWidth) . "," . ($karoBreakBot[$i] + $breakWidth) . ",NA";
}
print RFILE ")\n";

print RFILE "karoBreakBotY <- c(0.5,1,1,NA,0.5,0,0,NA";
for (my $i = 1; $i < scalar @karoBreakBot; $i++)
{
	print RFILE ",0.5,1,1,NA,0.5,0,0,NA";
}
print RFILE ")\n";


# border
print RFILE "karoBorderTop <- c($karoBreakTop[0]," . ($karoBreakTop[0] + $breakWidth);
for (my $i = 1; $i < (scalar @karoBreakTop - 1); $i++)
{
	print RFILE "," . ($karoBreakTop[$i] - $breakWidth) . ",$karoBreakTop[$i]," . ($karoBreakTop[$i] + $breakWidth);
}
print RFILE "," . ($karoBreakTop[-1] - $breakWidth) . ",$karoBreakTop[-1])\n";

print RFILE "karoBorderTopOver <- c(0.5,1";
for (my $i = 1; $i < (scalar @karoBreakTop - 1); $i++)
{
	print RFILE ",1,0.5,1";
}
print RFILE ",1,0.5)\n";

print RFILE "karoBorderTopUnder <- c(0.5,0";
for (my $i = 1; $i < (scalar @karoBreakTop - 1); $i++)
{
	print RFILE ",0,0.5,0";
}
print RFILE ",0,0.5)\n";

print RFILE "karoBorderBot <- c($karoBreakBot[0]," . ($karoBreakBot[0] + $breakWidth);
for (my $i = 1; $i < (scalar @karoBreakBot - 1); $i++)
{
	print RFILE "," . ($karoBreakBot[$i] - $breakWidth) . ",$karoBreakBot[$i]," . ($karoBreakBot[$i] + $breakWidth);
}
print RFILE "," . ($karoBreakBot[-1] - $breakWidth) . ",$karoBreakBot[-1])\n";

print RFILE "karoBorderBotOver <- c(0.5,1";
for (my $i = 1; $i < (scalar @karoBreakBot - 1); $i++)
{
	print RFILE ",1,0.5,1";
}
print RFILE ",1,0.5)\n";

print RFILE "karoBorderBotUnder <- c(0.5,0";
for (my $i = 1; $i < (scalar @karoBreakBot - 1); $i++)
{
	print RFILE ",0,0.5,0";
}
print RFILE ",0,0.5)\n";

# print sv
if (scalar (@breakX1 > 0))
{
	print RFILE "breakX1 <- c($breakX1[0]";
	for (my $i = 1; $i < scalar @breakX1; $i++)
	{
		print RFILE ", $breakX1[$i]";
	}
	print RFILE ")\n";
	
	print RFILE "breakY1 <- c($breakY1[0]";
	for (my $i = 1; $i < scalar @breakY1; $i++)
	{
		print RFILE ", $breakY1[$i]";
	}
	print RFILE ")\n";
	print RFILE "breakX2 <- c($breakX2[0]";
	for (my $i = 1; $i < scalar @breakX2; $i++)
	{
	print RFILE ", $breakX2[$i]";
	}
	print RFILE ")\n";
	
	print RFILE "breakY2 <- c($breakY2[0]";
	for (my $i = 1; $i < scalar @breakY2; $i++)
	{
		print RFILE ", $breakY2[$i]";
	}
	print RFILE ")\n";

	print RFILE "breakCol <- c(\"$breakCol[0]\"";
	for (my $i = 1; $i < scalar @breakCol; $i++)
	{
		print RFILE ", \"$breakCol[$i]\"";
	}
	print RFILE ")\n";
}


# genes of interest ticks
my $geneBump = 18000000;
if (scalar @geneNameTop > 0)
{
	print RFILE "geneNameTop <- c(\"$geneNameTop[0]\"";
	for (my $i = 1; $i < scalar @geneNameTop; $i++)
	{
		print RFILE ", \"$geneNameTop[$i]\"";
	}
	print RFILE ")\n";
	
	print RFILE "genePosTop <- c($genePosTop[0]";
	for (my $i = 1; $i < scalar @genePosTop; $i++)
	{
		if ($genePosTop[$i - 1] + $geneBump > $genePosTop[$i])
		{
			$genePosTop[$i] += $geneBump - ($genePosTop[$i] - $genePosTop[$i - 1]);
		}
		print RFILE ", \"$genePosTop[$i]\"";
	}
	print RFILE ")\n";
}

if (scalar @geneNameBot > 0)
{
	print RFILE "geneNameBot <- c(\"$geneNameBot[0]\"";
	for (my $i = 1; $i < scalar @geneNameBot; $i++)
	{
		print RFILE ", \"$geneNameBot[$i]\"";
	}
	print RFILE ")\n";
	
	print RFILE "genePosBot <- c($genePosBot[0]";
	for (my $i = 1; $i < scalar @genePosBot; $i++)
	{
		if ($genePosBot[$i - 1] + $geneBump > $genePosBot[$i])
		{
			$genePosBot[$i] += $geneBump - ($genePosBot[$i] - $genePosBot[$i - 1]);
		}
		print RFILE ", \"$genePosBot[$i]\"";
	}
	print RFILE ")\n";
}


my $headerString = "   $tSample (${tCoverage}x) // $nSample (${nCoverage}x)     $count{snv} somatic SNVs, " . ($count{ins} + $count{del}) . " somatic indels";
my $footerString = "   $outputFile $date";



print RFILE "png(filename = \"$outputFile\", width = 2200, height = 1700)\n";

print RFILE "par(oma=c(2.5,0,2.5,0),mgp=c(1,0,-1), cex=2)\n";
print RFILE "par(fig=c(0,1,0.86,1), mar=c(0, 4, 0, 0) + 0.1)\n";
print RFILE "plot(c(0,$chrTopOffset{end}), c(1,$rainMaxTop), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", log = \"y\", main=\"\", xlab=\"\", ylab=\"Somatic SNV\\nRainfall\",)\n";
print RFILE "abline(v=vertsTop, col=\"darkgrey\")\n";
print RFILE "points(xvalsTop, yvalsTop, col=colsTop, pch = 20)\n";
print RFILE "axis(side=3, at=atLabelsTop, labels=chrLabelsTop, tick=FALSE)\n";
if (scalar @geneNameTop > 0)
{
	print RFILE "text(genePosTop,1, labels=geneNameTop, pos=4, srt=45, cex=0.75)\n";
}

print RFILE "par(fig=c(0,1,0.72,0.86), mar=c(0, 4, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrTopOffset{end}), c(0, 1), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", ylim=c(0, 1), main=\"\", xlab=\"\", ylab=\"Germline Het\\nFreq in Tumour\",)\n";
print RFILE "abline(v=vertsTop, col=\"darkgrey\")\n";
print RFILE "points(hetPosTop, hetFreqsTop, col=\"black\", pch=\".\")\n";

print RFILE "par(fig=c(0,1,0.58,0.72), mar=c(0, 4, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrTopOffset{end}), c(-4,4), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", main=\"\", xlab=\"\", ylab=\"Copy Number\\nVariation\",)\n";
print RFILE "abline(v=vertsTop, col=\"darkgrey\")\n";
print RFILE "points(cnPosTop, cnValTop, col=cnColTop, pch=20)\n";
print RFILE "segments(cnStartTop, cnStateTop, cnEndTop, cnStateTop, col=\"black\", lwd=2)\n";

print RFILE "par(fig=c(0,1,0.55,0.58), mar=c(0, 4, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrTopOffset{end}), c(0,1),type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\")\n";
print RFILE "rect(karoStartTop, 0, karoEndTop, 1, col=karoColTop, border=NA)\n";
print RFILE "polygon(karoBreakTop,karoBreakTopY,col=\"white\",border=NA)\n";
print RFILE "lines(karoBorderTop,karoBorderTopOver,col=\"black\")\n";
print RFILE "lines(karoBorderTop,karoBorderTopUnder,col=\"black\")\n";


print RFILE "par(fig=c(0,1,0.45,0.55), mar=c(0, 4, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$bpWidth), c(0,1), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"Structural\\nBreak Points\")\n";

if (scalar (@breakX1 > 0))
{
	print RFILE "segments(breakX1,breakY1,breakX2,breakY2,col=breakCol)\n";
}


print RFILE "par(fig=c(0,1,0.42,0.45), mar=c(0, 4, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrBotOffset{end}), c(0,1), type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\")\n";
print RFILE "rect(karoStartBot, 0, karoEndBot, 1, col=karoColBot, border=NA)\n";
print RFILE "polygon(karoBreakBot,karoBreakBotY,col=\"white\",border=NA)\n";
print RFILE "lines(karoBorderBot,karoBorderBotOver,col=\"black\")\n";
print RFILE "lines(karoBorderBot,karoBorderBotUnder,col=\"black\")\n";

print RFILE "par(fig=c(0,1,0,0.14), mar=c(0, 4, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrBotOffset{end}), c(1,$rainMaxBot), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", log = \"y\", main=\"\", xlab=\"\", ylab=\"Somatic SNV\\nRainfall\",)\n";
print RFILE "abline(v=vertsBot, col=\"darkgrey\")\n";
print RFILE "points(xvalsBot, yvalsBot, col=colsBot, pch=20)\n";
print RFILE "axis(side=1, at=atLabelsBot, labels=chrLabelsBot, tick=FALSE)\n";
if (scalar @geneNameBot > 0)
{
	print RFILE "text(genePosBot,1, labels=geneNameBot, pos=4, srt=45, cex=0.75)\n";
}

print RFILE "par(fig=c(0,1,0.14,0.28), mar=c(0, 4, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrBotOffset{end}), c(0, 1), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", ylim=c(0, 1), main=\"\", xlab=\"\", ylab=\"Germline Het\\nFreq in Tumour\",)\n";
print RFILE "abline(v=vertsBot, col=\"darkgrey\")\n";
print RFILE "points(hetPosBot, hetFreqsBot, col=\"black\", pch=\".\")\n";

print RFILE "par(fig=c(0,1,0.28,0.42), mar=c(0, 4, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrBotOffset{end}), c(-4,4), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", main=\"\", xlab=\"\", ylab=\"Copy Number\\nVariation\",)\n";
print RFILE "abline(v=vertsBot, col=\"darkgrey\")\n";
print RFILE "points(cnPosBot, cnValBot, col=cnColBot, pch=20)\n";
print RFILE "segments(cnStartBot, cnStateBot, cnEndBot, cnStateBot, col=\"black\", lwd=2)\n";

print RFILE "mtext(\"$headerString\", NORTH<-3, line=1, adj=0, outer=TRUE, cex=2)\n";
print RFILE "mtext(\"$footerString\", SOUTH<-1, line=1, adj=0, outer=TRUE, cex=1.5)\n";

print RFILE "dev.off()\n";


close RFILE;
`Rscript $outputFile.R`;





sub mean
{
	my $sum = 0;
	my $count = 0;
	my $array = $_[0];
	for my $val (@{$array})
	{
		$sum += $val;
		$count++;
	}
	return $sum/$count;
}



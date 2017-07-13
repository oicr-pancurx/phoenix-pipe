#!/usr/bin/perl

use strict;
use warnings;

use Bio::DB::Fasta;

my $vcfFile = $ARGV[0];
my $seqType = $ARGV[1];
my $outputFile = "$vcfFile.html";
my $outputDir = "${vcfFile}_plots";

my $relativeDir = $outputDir;
$relativeDir =~ s/.*\///;



my %fastaHandles = ("path" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/fasta/");
my %referenceHash = ("chunkSize" => 10000);



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


my %genesOfInterest = ("AKT3",1,"PARP1",1,"NRAS",1,"NOTCH2",1,"NOTCH2NL",1,"PIK3C2B",1,"PTCH2",1,"PIK3CD",1,"PIK3R3",1,"DPYD",1,"GSTM1",1,"EPHA2",1,"MTHFR",1,"PLK3",1,"CDC7",1,"WNT2B",1,"WNT3A",1,"WNT4",1,"WNT9A",1,"ARID1A",1,"MCL1",1,"PTEN",1,"RET",1,"MGMT",1,"FGFR2",1,"CYP2C19",1,"CYP2C9",1,"CYP2C8",1,"CYP2E1",1,"ABCC2",1,"WNT8B",1,"TET1",1,"CYP17A1",1,"CHEK1",1,"HRAS",1,"PIK3C2A",1,"CCND1",1,"ATM",1,"WEE1",1,"RRM1",1,"GSTP1",1,"WNT11",1,"CBL",1,"PGR",1,"SLC22A6",1,"SLCO2B1",1,"WT1",1,"KRAS",1,"ErbB3",1,"MDM2",1,"CDK4",1,"PIK3C2G",1,"WNT1",1,"PTPN11",1,"PXN",1,"GLI1",1,"WNT10B",1,"WNT5B",1,"SLCO1B1",1,"SLCO1B3",1,"SOCS2",1,"BRCA2",1,"FLT1",1,"FLT3",1,"RB1",1,"ERCC5",1,"AKT1",1,"PARP2",1,"ESR2",1,"TEP1",1,"PGF",1,"IGF1R",1,"MAP2K1",1,"IDH2",1,"CYP1A1",1,"DLL4",1,"CYP1A2",1,"CDH1",1,"TSC2",1,"PLK1",1,"TUBB3",1,"ERCC4",1,"SULT1A1",1,"SOCS1",1,"SH2B",1,"BRCA1",1,"ErbB2",1,"TP53",1,"RARA",1,"TOP2A",1,"AURKB",1,"PIK3R5",1,"NF1",1,"WNT3",1,"WNT9B",1,"SOCS3",1,"STAT3",1,"BCL2",1,"DCC",1,"PIK3C3",1,"TYMS",1,"SMAD4",1,"ROCK1",1,"AKT2",1,"MAP2K2",1,"NOTCH3",1,"XRCC1",1,"AURKC",1,"ERCC1",1,"PIK3R2",1,"ERCC2",1,"STK11",1,"CYP2A6",1,"CYP2B6",1,"AXL",1,"ALK",1,"MSH2",1,"ErbB4",1,"IDH1",1,"UGT1A1",1,"GLI2",1,"ERCC3",1,"WNT10A",1,"WNT6",1,"SRC",1,"AURKA",1,"TOP1",1,"GART",1,"BCR",1,"CHEK2",1,"CYP2D6",1,"GSTT1",1,"NF2",1,"EWSR1",1,"WNT7B",1,"PIK3CA",1,"MLH1",1,"RAF1",1,"CTNNB1",1,"VHL",1,"PIK3CB",1,"FANCD2",1,"TOP2B",1,"ATR",1,"MST1R",1,"TERC",1,"RASSF1",1,"WNT5A",1,"WNT7A",1,"SLC15A2",1,"KIT",1,"KDR",1,"PDGFRA",1,"FGFR3",1,"EIF4E",1,"PLK4",1,"ABCG2",1,"UGT2B15",1,"UGT2B17",1,"UGT2B7",1,"TET2",1,"PDGFRB",1,"DHFR",1,"FLT4",1,"PIK3R1",1,"TERT",1,"PLK2",1,"WNT8A",1,"NPM1",1,"ESR1",1,"NOTCH4",1,"CCND3",1,"TUBB",1,"SLC29A1",1,"TPMT",1,"SLC22A1",1,"SLC22A2",1,"EGFR",1,"BRAF",1,"MET",1,"SHH",1,"SMO",1,"PIK3CG",1,"CDK5",1,"CYP3A4",1,"EPHB4",1,"GLI3",1,"CYP3A5",1,"WNT16",1,"WNT2",1,"ABCB1",1,"FGFR1",1,"PTK2",1,"PTK2B",1,"ANGPT1",1,"ANGPT2",1,"NAT1",1,"NAT2",1,"ABL1",1,"JAK2",1,"CDKN2A",1,"NOTCH1",1,"PTCH1",1,"TSC1",1,"ROR2",1,"AR",1,"ARAF",1,"RNF43",1,"MAP2K4",1,"TGFBR2",1,"KDM6A",1);

my %chrOffset = (
	"chr1" => 0,
	"chr2" => 249250621,
	"chr3" => 492449994,
	"chr4" => 690472424,
	"chr5" => 881626700,
	"chr6" => 1062541960,
	"chr7" => 1233657027,
	"chr8" => 1392795690,
	"chr9" => 1539159712,
	"chr10" => 1680373143,
	"chr11" => 1815907890,
	"chr12" => 1950914406,
	"chr13" => 2084766301,
	"chr14" => 2199936179,
	"chr15" => 2307285719,
	"chr16" => 2409817111,
	"chr17" => 2500171864,
	"chr18" => 2581367074,
	"chr19" => 2659444322,
	"chr20" => 2718573305,
	"chr21" => 2781598825,
	"chr22" => 2829728720,
	"chrX" => 2881033286,
	"chrY" => 3036303846,
	"end" => 3095677412
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

my %contextComp = (
	"AA" => "TT",
	"AC" => "GT",
	"AG" => "CT",
	"AT" => "AT",
	"CA" => "TG",
	"CC" => "GG",
	"CG" => "CG",
	"CT" => "AG",
	"GA" => "TC",
	"GC" => "GC",
	"GG" => "CC",
	"GT" => "AC",
	"TA" => "TA",
	"TC" => "GA",
	"TG" => "CA",
	"TT" => "AA",
);

open (FILE, "<$vcfFile") or die "Couldn't open $vcfFile\n";
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
				$context = getBase($chr, $pos - 1, \%referenceHash, \%fastaHandles) . getBase($chr, $pos + 1, \%referenceHash, \%fastaHandles);
				if (exists $comp{$change})
				{
					$change = $comp{$change};
					$context = $contextComp{$context};
				}

				$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{change} = $change;
				$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{context} = $context;

				$count{snv}++;



				if (($change eq "C>T") or ($change eq "T>C"))
				{
					$count{ts}++;
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

			if ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /frameshift/)		# matches frameshift and nonframeshift
			{
				$count{"frame$type"}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";
				}
			}
			elsif ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^nonsynonymous/)
			{
				$count{nonsyn}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";
				}
			}
			elsif ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^stopgain/)
			{
				$count{nonsyn}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";
				}
			}
			elsif ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^stoploss/)
			{
				$count{nonsyn}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";
				}
			}
			elsif ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^splicing/)
			{
				$count{"splice$type"}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";
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


			freeMemChunks($chr, $pos, \%referenceHash);

		}
	}
}
close FILE;


my $titvRatio = "NA";
if (($count{snv} - $count{ts}) > 0)
{
	$titvRatio = sprintf("%0.4f", $count{ts} / ($count{snv} - $count{ts}));
}


# print html document
open (HTML, ">$outputFile") or die "Couldn't open $outputFile for write\n";

print HTML "<html>\n<head>\n";
print HTML "<title>VCF Report for $tSample</title>\n";
print HTML "<script src=\"./sorttable.js\"></script>\n";
print HTML "<style type=\"text/css\">\nth, td {\n  padding: 3px !important;\n}\ntable\n{\nborder-collapse:collapse;\n}\n/* Sortable tables */\ntable.sortable thead {\n\tbackground-color:#eee;\n\tcolor:#000000;\n\tfont-weight: bold;\n\tcursor: default;\n}\n</style>\n";
print HTML "</head>\n<body>\n";

print HTML "<p>VCF Report generated $date.</p>\n";

print HTML "<h2>Analyzed sample: $tSample</h2>\n";
print HTML "<h2>Matched sample: $nSample</h2>\n";
print HTML "<p>$tSample covered at ${tCoverage}x</p>\n";
print HTML "<p>$nSample covered at ${nCoverage}x</p>\n";
print HTML "<hr>\n";
print HTML "<p><a href=\"../../../../../$tSample/wgs/bwa/0.6.2/collapsed/$tSample.bam\">$tSample.bam</a><br>\n";
print HTML "<a href=\"../../../../../$nSample/wgs/bwa/0.6.2/collapsed/$nSample.bam\">$nSample.bam</a></p>\n";
if ($seqType eq "wgs")
{
	print HTML "<p><a href=\"../HMMcopy/0.1.1/\">HMMcopy results directory</a></p>\n";
	print HTML "<p><a href=\"../integration/$tSample.onePage.png\"><img src=\"../integration/$tSample.onePage.png\" width=\"220\" height=\"170\"></a></p>\n";
}
print HTML "<hr>\n";

print HTML "<h2>Total somatic SNVs called: $count{snv}</h2>\n";
print HTML "<p>Nonsynonymous SNVs called: $count{nonsyn}</p>\n";
print HTML "<p>Ti/Tv ratio: $count{ts}/" . ($count{snv} - $count{ts}) . " ($titvRatio)</p>\n";
print HTML "<p>Kataegis clusters detected: $count{kataegis}</p>\n";

print HTML "<img border=\"0\" src=\"$relativeDir/singleNucPlot.png\"/>\n";
print HTML "<img border=\"0\" src=\"$relativeDir/triNucPlot.png\"/>\n<br>\n";
print HTML "<img border=\"0\" src=\"$relativeDir/rainfallPlot.png\"/>\n<br>\n";
print HTML "<img border=\"0\" src=\"$relativeDir/snvNormalDepthPlot.png\"/>\n";
print HTML "<img border=\"0\" src=\"$relativeDir/snvTumourDepthPlot.png\"/>\n";
print HTML "<img border=\"0\" src=\"$relativeDir/snvTumourFreqPlot.png\"/>\n<br>\n";
print HTML "<img border=\"0\" src=\"$relativeDir/snvTumourFreqPlot-del.png\"/>\n<br>\n";

print HTML "<hr>\n";

print HTML "<h2>Total somatic indels called: " . ($count{ins} + $count{del}) . "</h2>\n";
print HTML "<p>Insertions called: $count{ins}, Frameshift insertions called: $count{frameins}</p>\n";
print HTML "<p>Deletions called: $count{del}, Frameshift deletions called: $count{framedel}</p>\n";

print HTML "<img border=\"0\" src=\"$relativeDir/insLength.png\"/>\n";
print HTML "<img border=\"0\" src=\"$relativeDir/delLength.png\"/>\n<br>\n";
print HTML "<img border=\"0\" src=\"$relativeDir/indelNormalDepthPlot.png\"/>\n";
print HTML "<img border=\"0\" src=\"$relativeDir/indelTumourDepthPlot.png\"/>\n";
print HTML "<img border=\"0\" src=\"$relativeDir/indelTumourFreqPlot.png\"/>\n<br>\n";

print HTML "<hr>\n";

print HTML "<h2>Hits in genes of interest</h2>\n";

print HTML "<table border=\"1\" class=\"sortable\">\n<thead>\n<tr>\n";
print HTML "<th>Gene</th><th>Chr</th><th>Position</th><th>Ref</th><th>Alt</th><th>Frequency</th><th>Depth</th><th>External ID</th><th>Consequence</th>\n</tr>\n</thead>\n<tbody>\n";


# $subCounts{$vcfHash{$chr}{$pos}{snv}{$alt}
for $l (@vcfOrder)
{
	($chr, $pos) = split(/:/, $l);

	for $type (keys %{ $vcfHash{$chr}{$pos} })
	{
		for $change (keys %{ $vcfHash{$chr}{$pos}{$type} })
		{
			if (exists $vcfHash{$chr}{$pos}{$type}{$change}{print})
			{
				$id = $vcfHash{$chr}{$pos}{$type}{$change}{id};
				($ref, $alt) = split(/>/, $change);

				$gene = $vcfHash{$chr}{$pos}{$type}{$change}{gene};
				$freq = $vcfHash{$chr}{$pos}{$type}{$change}{tfreq};
				$depth = $vcfHash{$chr}{$pos}{$type}{$change}{tdepth};
				$consequence = $vcfHash{$chr}{$pos}{$type}{$change}{consequence};
				$consequence =~ s/$gene/ $gene/g;

				print HTML "<tr><td>$gene</td><td>$chr</td><td>$pos</td><td>$ref</td><td>$alt</td><td>$freq</td><td>$depth</td><td>$id</td><td>$consequence</td></tr>\n";
			}
		}
	}

}
print HTML "</tbody>\n</table>\n";

print HTML "</html>\n";

close HTML;



# make plots
unless (-d $outputDir)
{
	`mkdir $outputDir`;
}

# make snv substitution type plot
my @subTypes = qw/C>A C>G C>T T>A T>C T>G/;
my %subCounts;

for $change (@subTypes)
{
	$subCounts{$change} = 0;
}

# $vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{change}
for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{snv} })
		{
			$subCounts{$vcfHash{$chr}{$pos}{snv}{$alt}{change}}++;
		}
	}
}


open (RFILE, ">$outputDir/singleNucPlot.R") or die "Couldn't open $outputDir/singleNucPlot.R\n";

print RFILE "values <- c($subCounts{$subTypes[0]}";
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	print RFILE ", $subCounts{$subTypes[$i]}";
}
print RFILE ")\n";

print RFILE "labels <- c(\"$subTypes[0]\"";
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	print RFILE ", \"$subTypes[$i]\"";
}
print RFILE ")\n";

print RFILE "cols <- c(\"$colours{$subTypes[0]}\"";
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	print RFILE ", \"$colours{$subTypes[$i]}\"";
}
print RFILE ")\n";

print RFILE "png(filename = \"$outputDir/singleNucPlot.png\", width = 525, height = 350)\n";
print RFILE "barplot(values, col=cols, names.arg=labels, main=\"Substitution Types\", xlab=\"\", ylab=\"Count\", border=\"black\")\n";
#print RFILE "barplot(values)\n";
print RFILE "dev.off()\n";


close RFILE;
`Rscript $outputDir/singleNucPlot.R`;





# make trinucleotide context plot
my @contextTypes = qw/AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT/;
%subCounts = ();

for my $change (@subTypes)
{
	for my $context (@contextTypes)
	{
		$subCounts{$change}{$context} = 0;
	}
}

for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{snv} })
		{
			$subCounts{$vcfHash{$chr}{$pos}{snv}{$alt}{change}}{$vcfHash{$chr}{$pos}{snv}{$alt}{context}}++;
		}
	}
}

my $contextList = "(";
for my $context (@contextTypes)
{
	$contextList .= substr($context, 0, 1) . "N" . substr($context, 1) . " ";
}
$contextList =~ s/ $/)/;

open (RFILE, ">$outputDir/triNucPlot.R") or die "Couldn't open $outputDir/triNucPlot.R\n";

print RFILE "values <- c($subCounts{$subTypes[0]}{$contextTypes[0]}";
for (my $j = 1; $j < scalar @contextTypes; $j++)
{
	print RFILE ", $subCounts{$subTypes[0]}{$contextTypes[$j]}";
}
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	for (my $j = 0; $j < scalar @contextTypes; $j++)
	{
		print RFILE ", $subCounts{$subTypes[$i]}{$contextTypes[$j]}";
	}
}
print RFILE ")\n";

print RFILE "labels <- c(\"$subTypes[0]\"";
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	print RFILE ", \"$subTypes[$i]\"";
}
print RFILE ")\n";

print RFILE "cols <- c(\"$colours{$subTypes[0]}\"";
for (my $j = 1; $j < scalar @contextTypes; $j++)
{
	print RFILE ", \"$colours{$subTypes[0]}\"";
}
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	for (my $j = 0; $j < scalar @contextTypes; $j++)
	{
		print RFILE ", \"$colours{$subTypes[$i]}\"";
	}
}
print RFILE ")\n";

print RFILE "png(filename = \"$outputDir/triNucPlot.png\", width = 525, height = 350)\n";
print RFILE "barplot(matrix(values,16,6), col=cols, names.arg=labels, main=\"Trinucleotide Context\", xlab=\"$contextList\", ylab=\"Count\", border=\"black\", beside=TRUE)\n";
print RFILE "dev.off()\n";


close RFILE;
`Rscript $outputDir/triNucPlot.R`;





# make rainfall plot
my @positions;
my @distances;
my @colours;


for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{snv} })
		{
			if (exists $chrOffset{$chr})
			{
				if (exists $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist})
				{
					push(@positions, $pos + $chrOffset{$chr});
					push(@distances, $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist});
					push(@colours, $colours{$vcfHash{$chr}{$pos}{snv}{$alt}{change}});
				}
			}
		}
	}
}

my @verticals;
for $chr (keys %chrOffset)
{
	push(@verticals, $chrOffset{$chr});
}

my @chrList = qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY end/;
my @atLabels;
my @chrLabels;
for (my $i = 0; $i < (scalar @chrList - 1); $i++)
{
	push(@atLabels, ($chrOffset{$chrList[$i]} + $chrOffset{$chrList[$i + 1]}) / 2);

	push(@chrLabels, $chrList[$i]);
	$chrLabels[$i] =~ s/chr//;
}


open (RFILE, ">$outputDir/rainfallPlot.R") or die "Couldn't open >$outputDir/rainfallPlot.R\n";


print RFILE "xvals <- c($positions[0]";
for (my $i = 1; $i < scalar @positions; $i++)
{
	print RFILE ", \"$positions[$i]\"";
}
print RFILE ")\n";

print RFILE "yvals <- c($distances[0]";
for (my $i = 1; $i < scalar @distances; $i++)
{
	print RFILE ", \"$distances[$i]\"";
}
print RFILE ")\n";

print RFILE "cols <- c(\"$colours[0]\"";
for (my $i = 1; $i < scalar @colours; $i++)
{
	print RFILE ", \"$colours[$i]\"";
}
print RFILE ")\n";

print RFILE "verts <- c($verticals[0]";
for (my $i = 1; $i < scalar @verticals; $i++)
{
	print RFILE ", $verticals[$i]";
}
print RFILE ")\n";

print RFILE "atLabels <- c($atLabels[0]";
for (my $i = 1; $i < scalar @atLabels; $i++)
{
	print RFILE ", $atLabels[$i]";
}
print RFILE ")\n";

print RFILE "chrLabels <- c(\"$chrLabels[0]\"";
for (my $i = 1; $i < scalar @chrLabels; $i++)
{
	print RFILE ", \"$chrLabels[$i]\"";
}
print RFILE ")\n";

print RFILE "png(filename = \"$outputDir/rainfallPlot.png\", width = 1100, height = 400)\n";
print RFILE "plot(xvals, yvals, type = \"n\", bty=\"n\", xaxt=\"n\", log = \"y\", main=\"Rainfall Plot\", xlab=\"Genomic Position\", ylab=\"Genomic Distance Between SNVs\",)\n";
print RFILE "abline(v=verts, col=\"darkgrey\")\n";
print RFILE "points(xvals, yvals, col=cols, pch = 19)\n";
print RFILE "axis(side=1, at=atLabels, labels=chrLabels, tick=FALSE)\n";
print RFILE "dev.off()\n";


close RFILE;
`Rscript $outputDir/rainfallPlot.R`;


# make indel plots

my %insCounts;
my @insLabels = qw/1 2 3 4 5 6 7 8 9 10+/;
my $indelSize;

for $indelSize (@insLabels)
{
	$insCounts{$indelSize} = 0;
}

for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{ins} })
		{
			$indelSize = $vcfHash{$chr}{$pos}{ins}{$alt}{size};
			if ($indelSize >= 10)
			{
				$indelSize = "10+";
			}

			$insCounts{$indelSize}++;
		}
	}
}

open (RFILE, ">$outputDir/insLength.R") or die "Couldn't open >$outputDir/insLength.R\n";

print RFILE "labs <- c(\"$insLabels[0]\"";
for (my $i = 1; $i < scalar @insLabels; $i++)
{
	print RFILE ", \"$insLabels[$i]\"";
}
print RFILE ")\n";

print RFILE "vals <- c($insCounts{$insLabels[0]}";
for (my $i = 1; $i < scalar @insLabels; $i++)
{
	print RFILE ", $insCounts{$insLabels[$i]}";
}
print RFILE ")\n";

print RFILE "png(filename = \"$outputDir/insLength.png\", width = 525, height = 350)\n";
print RFILE "barplot(vals, col=\"$colours{ins}\", names.arg=labs, main=\"Insertion Size Distribution\", xlab=\"Size (bp)\", ylab=\"Count\", border=\"black\")\n";
print RFILE "dev.off()\n";


close RFILE;
`Rscript $outputDir/insLength.R`;



my %delCounts;
my @delLabels = qw/1 2 3 4 5 6 7 8 9 10+/;

for $indelSize (@delLabels)
{
	$delCounts{$indelSize} = 0;
}

for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{del} })
		{
			$indelSize = $vcfHash{$chr}{$pos}{del}{$alt}{size};
			if ($indelSize >= 10)
			{
				$indelSize = "10+";
			}

			$delCounts{$indelSize}++;
		}
	}
}

open (RFILE, ">$outputDir/delLength.R") or die "Couldn't open >$outputDir/delLength.R\n";

print RFILE "labs <- c(\"$delLabels[0]\"";
for (my $i = 1; $i < scalar @delLabels; $i++)
{
	print RFILE ", \"$delLabels[$i]\"";
}
print RFILE ")\n";

print RFILE "vals <- c($delCounts{$delLabels[0]}";
for (my $i = 1; $i < scalar @delLabels; $i++)
{
	print RFILE ", $delCounts{$delLabels[$i]}";
}
print RFILE ")\n";

print RFILE "png(filename = \"$outputDir/delLength.png\", width = 525, height = 350)\n";
print RFILE "barplot(vals, col=\"$colours{del}\", names.arg=labs, main=\"Deletion Size Distribution\", xlab=\"Size (bp)\", ylab=\"Count\", border=\"black\")\n";
print RFILE "dev.off()\n";


close RFILE;
`Rscript $outputDir/delLength.R`;



# make depth plots

# snv normal DOC distribution

my $bin;
my $maxBin = 0;
my $maxMaxBin = 250;
my %depthBins;

for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{snv} })
		{
			$bin = int($vcfHash{$chr}{$pos}{snv}{$alt}{ndepth} / 5) * 5;
			$depthBins{$vcfHash{$chr}{$pos}{snv}{$alt}{change}}{$bin}++;
			if ($bin > $maxBin)
			{
				$maxBin = $bin;
			}
		}
	}
}

if ($maxBin > $maxMaxBin)
{
	$maxBin = $maxMaxBin;
}

my $numBins;
for $change (@subTypes)
{
	$numBins = 0;
	for ($bin = 0; $bin <= $maxBin; $bin += 5)
	{
		$numBins++;
		unless (exists $depthBins{$change}{$bin})
		{
			$depthBins{$change}{$bin} = 0;
		}
	}
}

open (RFILE, ">$outputDir/snvNormalDepthPlot.R") or die "Couldn't open >$outputDir/snvNormalDepthPlot.R\n";

print RFILE "vals <- c($depthBins{$subTypes[0]}{0}";
for ($bin = 5; $bin <= $maxBin; $bin += 5)
{
	print RFILE ", $depthBins{$subTypes[0]}{$bin}";
}
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	for ($bin = 0; $bin <= $maxBin; $bin += 5)
	{
		print RFILE ", $depthBins{$subTypes[$i]}{$bin}";
	}
}
print RFILE ")\n";

print RFILE "cols <- c(\"$colours{$subTypes[0]}\"";
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	print RFILE ", \"$colours{$subTypes[$i]}\"";
}
print RFILE ")\n";

print RFILE "labs <- c(\"0\"";
for ($bin = 20; $bin <= $maxBin; $bin += 20)
{
	print RFILE ", \"$bin\"";
}
print RFILE ")\n";

print RFILE "ticks <-c(0";
for (my $i = 4; $i <= ($maxBin/5); $i += 4)
{
	print RFILE ", $i";
}
print RFILE ")\n";

print RFILE "png(filename = \"$outputDir/snvNormalDepthPlot.png\", width = 525, height = 350)\n";
print RFILE "bp <- barplot(matrix(vals,6,$numBins,byrow=TRUE), col=cols, main=\"Normal DOC of Reported SNVs\", xlab=\"Depth of Coverage\", ylab=\"Count\", border=\"black\")\n";
print RFILE "axis(side=1, at=bp[1 + ticks] - 0.6, lab=labs)\n";
print RFILE "dev.off()\n";


close RFILE;
`Rscript $outputDir/snvNormalDepthPlot.R`;



# snv tumour DOC distribution

$maxBin = 0;
%depthBins = ();

for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{snv} })
		{
			$bin = int($vcfHash{$chr}{$pos}{snv}{$alt}{tdepth} / 5) * 5;
			$depthBins{$vcfHash{$chr}{$pos}{snv}{$alt}{change}}{$bin}++;
			if ($bin > $maxBin)
			{
				$maxBin = $bin;
			}
		}
	}
}

if ($maxBin > $maxMaxBin)
{
	$maxBin = $maxMaxBin;
}
for $change (@subTypes)
{
	$numBins = 0;
	for ($bin = 0; $bin <= $maxBin; $bin += 5)
	{
		$numBins++;
		unless (exists $depthBins{$change}{$bin})
		{
			$depthBins{$change}{$bin} = 0;
		}
	}
}

open (RFILE, ">$outputDir/snvTumourDepthPlot.R") or die "Couldn't open >$outputDir/snvTumourDepthPlot.R\n";

print RFILE "vals <- c($depthBins{$subTypes[0]}{0}";
for ($bin = 5; $bin <= $maxBin; $bin += 5)
{
	print RFILE ", $depthBins{$subTypes[0]}{$bin}";
}
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	for ($bin = 0; $bin <= $maxBin; $bin += 5)
	{
		print RFILE ", $depthBins{$subTypes[$i]}{$bin}";
	}
}
print RFILE ")\n";

print RFILE "cols <- c(\"$colours{$subTypes[0]}\"";
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	print RFILE ", \"$colours{$subTypes[$i]}\"";
}
print RFILE ")\n";

print RFILE "labs <- c(\"0\"";
for ($bin = 20; $bin <= $maxBin; $bin += 20)
{
	print RFILE ", \"$bin\"";
}
print RFILE ")\n";

print RFILE "ticks <-c(0";
for (my $i = 4; $i <= ($maxBin/5); $i += 4)
{
	print RFILE ", $i";
}
print RFILE ")\n";

print RFILE "png(filename = \"$outputDir/snvTumourDepthPlot.png\", width = 525, height = 350)\n";
print RFILE "bp <- barplot(matrix(vals,6,$numBins,byrow=TRUE), col=cols, main=\"Tumour DOC of Reported SNVs\", xlab=\"Depth of Coverage\", ylab=\"Count\", border=\"black\")\n";
print RFILE "axis(side=1, at=bp[1 + ticks] - 0.6, lab=labs)\n";
print RFILE "dev.off()\n";


close RFILE;
`Rscript $outputDir/snvTumourDepthPlot.R`;




# snv freq distribution

$maxBin = 0;
%depthBins = ();

my $binSize = 2;
$maxBin = 100;

for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{snv} })
		{
			$bin = int(($vcfHash{$chr}{$pos}{snv}{$alt}{tfreq} * 100) / $binSize) * $binSize;
			if ($bin == $maxBin)
			{
				$bin = $maxBin - $binSize;
			}
			$depthBins{$vcfHash{$chr}{$pos}{snv}{$alt}{change}}{$bin}++;
		}
	}
}

for $change (@subTypes)
{
	$numBins = 0;
	for ($bin = 0; $bin < $maxBin; $bin += $binSize)
	{
		$numBins++;
		unless (exists $depthBins{$change}{$bin})
		{
			$depthBins{$change}{$bin} = 0;
		}
	}
}

open (RFILE, ">$outputDir/snvTumourFreqPlot.R") or die "Couldn't open >$outputDir/snvTumourFreqPlot.R\n";

print RFILE "vals <- c($depthBins{$subTypes[0]}{0}";
for ($bin = $binSize; $bin < $maxBin; $bin += $binSize)
{
	print RFILE ", $depthBins{$subTypes[0]}{$bin}";
}
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	for ($bin = 0; $bin < $maxBin; $bin += $binSize)
	{
		print RFILE ", $depthBins{$subTypes[$i]}{$bin}";
	}
}
print RFILE ")\n";

print RFILE "cols <- c(\"$colours{$subTypes[0]}\"";
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	print RFILE ", \"$colours{$subTypes[$i]}\"";
}
print RFILE ")\n";

print RFILE "labs <- c(\"0%\"";
for ($bin = $binSize * 5; $bin <= $maxBin; $bin += $binSize * 5)
{
	print RFILE ", \"$bin%\"";
}
print RFILE ")\n";

print RFILE "ticks <-c(0.1";
for (my $i = 5; $i <= ($maxBin/$binSize); $i += 5)
{
	print RFILE ", " . ($i + 0.2*$i + 0.1);
}
print RFILE ")\n";

print RFILE "png(filename = \"$outputDir/snvTumourFreqPlot.png\", width = 1100, height = 400)\n";
print RFILE "bp <- barplot(matrix(vals,6,$numBins,byrow=TRUE), col=cols, main=\"Tumour Frequency of Reported SNVs\", xlab=\"Frequency\", ylab=\"Count\", border=\"black\")\n";
print RFILE "axis(side=1, at=ticks, lab=labs)\n";
print RFILE "dev.off()\n";


close RFILE;
`Rscript $outputDir/snvTumourFreqPlot.R`;






# non-synonymous snv freq distribution

$maxBin = 0;
%depthBins = ();

$binSize = 2;
$maxBin = 100;

for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{snv} })
		{
			$consequence = $vcfHash{$chr}{$pos}{snv}{$alt}{consequence};
			if (($consequence =~ /^nonsynonymous/) or ($consequence =~ /^stopgain/) or ($consequence =~ /^frameshift/) or ($consequence =~ /^splicing/))
			{
				$bin = int(($vcfHash{$chr}{$pos}{snv}{$alt}{tfreq} * 100) / $binSize) * $binSize;
				if ($bin == $maxBin)
				{
					$bin = $maxBin - $binSize;
				}
				$depthBins{$vcfHash{$chr}{$pos}{snv}{$alt}{change}}{$bin}++;
			}
		}
	}
}

for $change (@subTypes)
{
	$numBins = 0;
	for ($bin = 0; $bin < $maxBin; $bin += $binSize)
	{
		$numBins++;
		unless (exists $depthBins{$change}{$bin})
		{
			$depthBins{$change}{$bin} = 0;
		}
	}
}

open (RFILE, ">$outputDir/snvTumourFreqPlot-del.R") or die "Couldn't open >$outputDir/snvTumourFreqPlot-del.R\n";

print RFILE "vals <- c($depthBins{$subTypes[0]}{0}";
for ($bin = $binSize; $bin < $maxBin; $bin += $binSize)
{
	print RFILE ", $depthBins{$subTypes[0]}{$bin}";
}
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	for ($bin = 0; $bin < $maxBin; $bin += $binSize)
	{
		print RFILE ", $depthBins{$subTypes[$i]}{$bin}";
	}
}
print RFILE ")\n";

print RFILE "cols <- c(\"$colours{$subTypes[0]}\"";
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	print RFILE ", \"$colours{$subTypes[$i]}\"";
}
print RFILE ")\n";

print RFILE "labs <- c(\"0%\"";
for ($bin = $binSize * 5; $bin <= $maxBin; $bin += $binSize * 5)
{
	print RFILE ", \"$bin%\"";
}
print RFILE ")\n";

print RFILE "ticks <-c(0.1";
for (my $i = 5; $i <= ($maxBin/$binSize); $i += 5)
{
	print RFILE ", " . ($i + 0.2*$i + 0.1);
}
print RFILE ")\n";

print RFILE "png(filename = \"$outputDir/snvTumourFreqPlot-del.png\", width = 1100, height = 400)\n";
print RFILE "bp <- barplot(matrix(vals,6,$numBins,byrow=TRUE), col=cols, main=\"Tumour Frequency of Deleterious SNVs\", xlab=\"Frequency\", ylab=\"Count\", border=\"black\")\n";
print RFILE "axis(side=1, at=ticks, lab=labs)\n";
print RFILE "dev.off()\n";


close RFILE;
`Rscript $outputDir/snvTumourFreqPlot-del.R`;






# indel normal DOC distribution

$maxBin = 0;
%depthBins = ();
@subTypes = qw/ins del/;

for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{ins} })
		{
			$bin = int($vcfHash{$chr}{$pos}{ins}{$alt}{ndepth} / 5) * 5;
			$depthBins{ins}{$bin}++;
			if ($bin > $maxBin)
			{
				$maxBin = $bin;
			}
		}
		for $alt (keys %{ $vcfHash{$chr}{$pos}{del} })
		{
			$bin = int($vcfHash{$chr}{$pos}{del}{$alt}{ndepth} / 5) * 5;
			$depthBins{del}{$bin}++;
			if ($bin > $maxBin)
			{
				$maxBin = $bin;
			}
		}
	}
}
if ($maxBin > $maxMaxBin)
{
	$maxBin = $maxMaxBin;
}
for $change (@subTypes)
{
	$numBins = 0;
	for ($bin = 0; $bin <= $maxBin; $bin += 5)
	{
		$numBins++;
		unless (exists $depthBins{$change}{$bin})
		{
			$depthBins{$change}{$bin} = 0;
		}
	}
}

open (RFILE, ">$outputDir/indelNormalDepthPlot.R") or die "Couldn't open >$outputDir/indelNormalDepthPlot.R\n";

print RFILE "vals <- c($depthBins{$subTypes[0]}{0}";
for ($bin = 5; $bin <= $maxBin; $bin += 5)
{
	print RFILE ", $depthBins{$subTypes[0]}{$bin}";
}
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	for ($bin = 0; $bin <= $maxBin; $bin += 5)
	{
		print RFILE ", $depthBins{$subTypes[$i]}{$bin}";
	}
}
print RFILE ")\n";

print RFILE "cols <- c(\"$colours{$subTypes[0]}\"";
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	print RFILE ", \"$colours{$subTypes[$i]}\"";
}
print RFILE ")\n";

print RFILE "labs <- c(\"0\"";
for ($bin = 20; $bin <= $maxBin; $bin += 20)
{
	print RFILE ", \"$bin\"";
}
print RFILE ")\n";

print RFILE "ticks <-c(0";
for (my $i = 4; $i <= ($maxBin/5); $i += 4)
{
	print RFILE ", $i";
}
print RFILE ")\n";

print RFILE "png(filename = \"$outputDir/indelNormalDepthPlot.png\", width = 525, height = 350)\n";
print RFILE "bp <- barplot(matrix(vals,2,$numBins,byrow=TRUE), col=cols, main=\"Normal DOC of Reported Indels\", xlab=\"Depth of Coverage\", ylab=\"Count\", border=\"black\")\n";
print RFILE "axis(side=1, at=bp[1 + ticks] - 0.6, lab=labs)\n";
print RFILE "dev.off()\n";


close RFILE;
`Rscript $outputDir/indelNormalDepthPlot.R`;



# indel tumour DOC distribution

$maxBin = 0;
%depthBins = ();

for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{ins} })
		{
			$bin = int($vcfHash{$chr}{$pos}{ins}{$alt}{tdepth} / 5) * 5;
			$depthBins{ins}{$bin}++;
			if ($bin > $maxBin)
			{
				$maxBin = $bin;
			}
		}
		for $alt (keys %{ $vcfHash{$chr}{$pos}{del} })
		{
			$bin = int($vcfHash{$chr}{$pos}{del}{$alt}{tdepth} / 5) * 5;
			$depthBins{del}{$bin}++;
			if ($bin > $maxBin)
			{
				$maxBin = $bin;
			}
		}
	}
}

if ($maxBin > $maxMaxBin)
{
	$maxBin = $maxMaxBin;
}
for $change (@subTypes)
{
	$numBins = 0;
	for ($bin = 0; $bin <= $maxBin; $bin += 5)
	{
		$numBins++;
		unless (exists $depthBins{$change}{$bin})
		{
			$depthBins{$change}{$bin} = 0;
		}
	}
}

open (RFILE, ">$outputDir/indelTumourDepthPlot.R") or die "Couldn't open >$outputDir/indelTumourDepthPlot.R\n";

print RFILE "vals <- c($depthBins{$subTypes[0]}{0}";
for ($bin = 5; $bin <= $maxBin; $bin += 5)
{
	print RFILE ", $depthBins{$subTypes[0]}{$bin}";
}
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	for ($bin = 0; $bin <= $maxBin; $bin += 5)
	{
		print RFILE ", $depthBins{$subTypes[$i]}{$bin}";
	}
}
print RFILE ")\n";

print RFILE "cols <- c(\"$colours{$subTypes[0]}\"";
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	print RFILE ", \"$colours{$subTypes[$i]}\"";
}
print RFILE ")\n";

print RFILE "labs <- c(\"0\"";
for ($bin = 20; $bin <= $maxBin; $bin += 20)
{
	print RFILE ", \"$bin\"";
}
print RFILE ")\n";

print RFILE "ticks <-c(0";
for (my $i = 4; $i <= ($maxBin/5); $i += 4)
{
	print RFILE ", $i";
}
print RFILE ")\n";

print RFILE "png(filename = \"$outputDir/indelTumourDepthPlot.png\", width = 525, height = 350)\n";
print RFILE "bp <- barplot(matrix(vals,2,$numBins,byrow=TRUE), col=cols, main=\"Tumour DOC of Reported Indels\", xlab=\"Depth of Coverage\", ylab=\"Count\", border=\"black\")\n";
print RFILE "axis(side=1, at=bp[1 + ticks] - 0.6, lab=labs)\n";
print RFILE "dev.off()\n";


close RFILE;
`Rscript $outputDir/indelTumourDepthPlot.R`;




# indel freq distribution

$maxBin = 0;
%depthBins = ();

$binSize = 2;
$maxBin = 100;

for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{ins} })
		{
			$bin = int(($vcfHash{$chr}{$pos}{ins}{$alt}{tfreq} * 100) / $binSize) * $binSize;
			if ($bin == $maxBin)
			{
				$bin = $maxBin - $binSize;
			}
			$depthBins{ins}{$bin}++;
		}
		for $alt (keys %{ $vcfHash{$chr}{$pos}{del} })
		{
			$bin = int(($vcfHash{$chr}{$pos}{del}{$alt}{tfreq} * 100) / $binSize) * $binSize;
			if ($bin == $maxBin)
			{
				$bin = $maxBin - $binSize;
			}
			$depthBins{del}{$bin}++;
		}
	}
}

for $change (@subTypes)
{
	$numBins = 0;
	for ($bin = 0; $bin < $maxBin; $bin += $binSize)
	{
		$numBins++;
		unless (exists $depthBins{$change}{$bin})
		{
			$depthBins{$change}{$bin} = 0;
		}
	}
}

open (RFILE, ">$outputDir/indelTumourFreqPlot.R") or die "Couldn't open >$outputDir/indelTumourFreqPlot.R\n";

print RFILE "vals <- c($depthBins{$subTypes[0]}{0}";
for ($bin = $binSize; $bin < $maxBin; $bin += $binSize)
{
	print RFILE ", $depthBins{$subTypes[0]}{$bin}";
}
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	for ($bin = 0; $bin < $maxBin; $bin += $binSize)
	{
		print RFILE ", $depthBins{$subTypes[$i]}{$bin}";
	}
}
print RFILE ")\n";

print RFILE "cols <- c(\"$colours{$subTypes[0]}\"";
for (my $i = 1; $i < scalar @subTypes; $i++)
{
	print RFILE ", \"$colours{$subTypes[$i]}\"";
}
print RFILE ")\n";

print RFILE "labs <- c(\"0%\"";
for ($bin = $binSize * 5; $bin <= $maxBin; $bin += $binSize * 5)
{
	print RFILE ", \"$bin%\"";
}
print RFILE ")\n";

print RFILE "ticks <-c(0.1";
for (my $i = 5; $i <= ($maxBin/$binSize); $i += 5)
{
	print RFILE ", " . ($i + 0.2*$i + 0.1);
}
print RFILE ")\n";

print RFILE "png(filename = \"$outputDir/indelTumourFreqPlot.png\", width = 1100, height = 400)\n";
print RFILE "bp <- barplot(matrix(vals,2,$numBins,byrow=TRUE), col=cols, main=\"Tumour Frequency of Reported Indels\", xlab=\"Frequency\", ylab=\"Count\", border=\"black\")\n";
print RFILE "axis(side=1, at=ticks, lab=labs)\n";
print RFILE "dev.off()\n";


close RFILE;
`Rscript $outputDir/indelTumourFreqPlot.R`;









# getBase returns the base at a specific position in the reference.  It will initialize fasta handles and pull new chunks of reference if necessary
# input is chromosome, position, reference hash and fasta handles
# output is a single base (and the reference hash and fasta handles may be modified)
sub getBase
{
    my $chr = $_[0];
    my $pos = $_[1];
    my $reference = $_[2];
    my $fastaHandles = $_[3];

    my $chunkStart = int(($pos - 1) / $reference->{"chunkSize"}) * $reference->{"chunkSize"} + 1;       # +1 because the first base in the reference is 1, $pos - 1 so that multiples of chunk size resolve to the correct chunk
    my $chunkEnd = $chunkStart + $reference->{"chunkSize"} - 1;

    unless (exists $reference->{$chr}{$chunkStart}{$pos})       # if the position isn't in our hash, we need to get a new chunk from the reference
    {
        unless (exists ($fastaHandles->{$chr}))     # create a handle for the chromosome fasta, if it doesn't exist
        {
#           warn "Creating fasta handle for $chr\n";
            $fastaHandles->{$chr} = Bio::DB::Fasta->new("$fastaHandles{path}/$chr.fa");
        }

#       warn "Pulling $chr:$chunkStart-$chunkEnd from fasta\n";
        my $newChunk = uc($fastaHandles{$chr}->seq($chr, $chunkStart, $chunkEnd));
        my $i = $chunkStart;
        for my $base (split("", $newChunk))
        {
            $reference->{$chr}{$chunkStart}{$i} = $base;
            $i++;
        }
    }
#   warn "returning $reference->{$chr}{$chunkStart}{$pos}\n";
    return $reference->{$chr}{$chunkStart}{$pos};
}


# getRange returns a string of bases from the reference in the specified range by calling getBase
# input is chromosome, start pos, end pos, reference hash and fasta handles
# output is a string of bases
sub getRange
{
    my $chr = $_[0];
    my $start = $_[1];
    my $end = $_[2];
    my $reference = $_[3];
    my $fastaHandles = $_[4];

    my $seq = "";

    for (my $p = $start; $p <= $end; $p++)
    {
        $seq .= getBase($chr, $p, $reference, $fastaHandles);
#       warn "Got base: $chr:$p\t$seq\n";
    }

    return $seq;
}


# freeMemChunks tests if the next indel to be processed is on a different chromosome or more than a chunk away from the reference sequences currently in memory
#   if there exist chunks that we won't need again (assuming the input is sorted) the chunks will be removed from the reference hash
# input is the chromosome and position of the current indel, and the reference hash
# there is no output
sub freeMemChunks
{
    my $chr = $_[0];
    my $pos = $_[1];
    my $reference = $_[2];

    # delete chunks from non-current chromosomes
    for my $refChr (keys %$reference)
    {
        if (($refChr ne $chr) and ($refChr ne "chunkSize"))
        {
#           warn "deleting all chunks for $refChr.\n";
            delete $reference->{$refChr};
        }
    }

    # delete chunks if they are more than 1.5 chunks away from the current indel
    # 1.5 so that we are at least in the middle of the current chunk before dropping the previous one
    for my $chunkPos (keys %{ $reference->{$chr} })
    {
        if ($chunkPos < ($pos - (1.5 * $reference->{"chunkSize"})))
        {
#           warn "deleting $chr:$chunkPos chunk.\n";
            delete $reference->{$chr}{$chunkPos};
        }
    }

    return;
}



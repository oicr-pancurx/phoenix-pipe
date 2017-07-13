#!/usr/bin/perl

use strict;
use warnings;

use JSON;

my $oldTimestamp = $ARGV[0];		# old timestamp to pull comments from

my $timestamp = `date +%y%m%d`;
chomp $timestamp;

unless (-d "./statusReport")
{
	`mkdir ./statusReport`;
}

my $wgsSeqFile = "./statusReport/PCSI-wgs-seq-report-$timestamp.tsv";
my $exomeSeqFile = "./statusReport/PCSI-exome-seq-report-$timestamp.tsv";

my $wgsAnalysisFile = "./statusReport/PCSI-wgs-analysis-report-$timestamp.tsv";
my $exomeAnalysisFile = "./statusReport/PCSI-exome-analysis-report-$timestamp.tsv";

my $allFile = "./statusReport/PCSI-all-report-$timestamp.tsv";

my $onlyMatched = 0;

my @sequencingHeaders = ("Institute","External ID","Donor","Prep Type","Sample","Deep KRAS","Mouse","Estimated Cellularity","","Lanes","Coverage","Status","Notes","Oldest Lane","Newest Lane");
my @analysisHeaders = ("Donor","Sample","Prep Type","Status","Notes","Action Items","","Deep KRAS","Deep KRAS Variant","Mouse","Coverage","Mismatch Error","Indel Error","","SNVs","Non-Silent SNVs","Ts/Tv","Ts/Tv Ratio","Indels","","CTX","DEL","INS","INV","ITX","","Copy State: 1","2","3","4","5","6");
my @allHeaders = ("Institute","External ID","Donor","Prep Type","Seq Type","Sample","Deep KRAS","Deep KRAS Variant","Mouse","","Lanes","Coverage","Mismatch Error","Indel Error","Status","Notes","Action Items","Oldest Lane","Newest Lane","SNVs","Non-Silent SNVs","Ts/Tv","Ts/Tv Ratio","Indels","","CTX","DEL","INS","INV","ITX","","Copy State: 1","2","3","4","5","6");

# these columns will keep the values from the old spreadsheet
my %seqKeepHeaders = ("Institute" => 1, "External ID" => 1, "Prep Type" => 1, "Estimated Cellularity" => 1, "Status" => 1, "Notes" => 1);
my %analysisKeepHeaders = ("Prep Type" => 1, "Status" => 1, "Notes" => 1, "Action Items" => 1);

# these columns will have an entry in the Notes column if they are updated
my %seqWarnHeaders = ("Lanes" => 1,"Coverage" => 1);
my %analysisWarnHeaders = ("Coverage" => 1, "SNVs" => 1, "CTX" => 1, "Copy State: 1" => 1);

my $xenomePath = "xenome/1.0.1-r/";
my $deepKRASpath = "deepKRAS";
my @aligners = qw(bwa/0.6.2);
my @strelkaVersions = qw(strelka/v1.0.7/ gatk/2.4.9/strelka/v1.0.7/);
my @mutectVersions = qw(mutect/1.1.4/ gatk/2.4.9/mutect/1.1.4/);

my $hmmcopyDir = "HMMcopy/0.1.1";
my $crestDir = "final_crest-delly";


my $date;
my $donor;
my $samp;
my $seqType;

my %sampleToDonor;

my $l;
my @f; 


my $oldWgsSeqFile;
my $oldExomeSeqFile;
my $oldWgsAnalysisFile;
my $oldExomeAnalysisFile;

my %oldSeqHash;
my %oldAnalysisHash;

# read old files
if (defined $oldTimestamp)
{
	$oldWgsSeqFile = "./statusReport/PCSI-wgs-seq-report-$oldTimestamp.tsv";
	$oldExomeSeqFile = "./statusReport/PCSI-exome-seq-report-$oldTimestamp.tsv";

	$oldWgsAnalysisFile = "./statusReport/PCSI-wgs-analysis-report-$oldTimestamp.tsv";
	$oldExomeAnalysisFile = "./statusReport/PCSI-exome-analysis-report-$oldTimestamp.tsv";

	readOldFile($oldWgsSeqFile, "wgs", \%oldSeqHash);
	readOldFile($oldExomeSeqFile, "exome", \%oldSeqHash);

	readOldFile($oldWgsAnalysisFile, "wgs", \%oldAnalysisHash);
	readOldFile($oldExomeAnalysisFile, "exome", \%oldAnalysisHash);
}
else
{
	warn "No comments imported from old reports!\n";
}

sub readOldFile
{
	my $file = $_[0];
	my $seqType = $_[1];
	my $oldSampleHash = $_[2];

	my $l;
	my @f;
	my $gotHeader = 0;

	my @headers;
	my %row;

	open (FILE, "<$file") or die "Couldn't open $file\n";

	while ($l = <FILE>)
	{
		chomp $l;
		if ($gotHeader == 0)
		{
			if ($l =~ /\tStatus\tNotes\t/)		# ugly way to detect header!!
			{
				@headers = split(/\t/, $l);
				$gotHeader = 1;

			}
		}
		else
		{
			@f = split(/\t/, $l);
			%row = ();
			for (my $i = 0; $i < scalar(@f); $i++)
			{
				$row{$headers[$i]} = $f[$i];
			}
			for my $i (@headers)
			{
				if (exists $row{$i})
				{
					$oldSampleHash->{$row{"Donor"}}{$row{"Sample"}}{$seqType}{$i} = $row{$i};
				}
			}
		}
	}
}








my $samples = `ls PCSI*/*/wgs/fastq/*gz PCSI*/*/exome/fastq/*gz PCSI*/*/wgs/fastq/*fastq PCSI*/*/exome/fastq/*fastq`;
chomp $samples;
my @samples = split(/\n/, $samples);

my %sampleHash;

for my $sample (@samples)
{
	$sample =~ /(.*?)\/(.*?)\/(.*?)\//;
	$donor = $1;
	$samp = $2;
	$seqType = $3;
	$sampleHash{$donor}{$samp}{$seqType}{fastqCount}++;


	$sample =~ /.*\/.*?_.*?_.*?_.*?_.*?_.*?_.*?_.*?_.*?_(.*?)_.*/;
	$date = $1;

	unless ($date =~ /^\d+$/)
	{
		warn $sample . "\n";
	}

	if (exists $sampleHash{$donor}{$samp}{$seqType}{oldestLane})
	{
		if ($date < $sampleHash{$donor}{$samp}{$seqType}{oldestLane})
		{
			$sampleHash{$donor}{$samp}{$seqType}{oldestLane} = $date;
		}
	}
	else
	{
		$sampleHash{$donor}{$samp}{$seqType}{oldestLane} = $date;
	}

	if (exists $sampleHash{$donor}{$samp}{$seqType}{newestLane})
	{
		if ($date > $sampleHash{$donor}{$samp}{$seqType}{newestLane})
		{
			$sampleHash{$donor}{$samp}{$seqType}{newestLane} = $date;
		}
	}
	else
	{
		$sampleHash{$donor}{$samp}{$seqType}{newestLane} = $date;
	}

}
# finalize
for my $donor (sort keys %sampleHash)
{
	for my $sample (sort keys %{ $sampleHash{$donor} })
	{
		for $seqType (sort keys %{ $sampleHash{$donor}{$sample} })
		{
			$sampleHash{$donor}{$sample}{$seqType}{"Donor"} = $donor;
			$sampleHash{$donor}{$sample}{$seqType}{"Sample"} = $sample;
			$sampleHash{$donor}{$sample}{$seqType}{"Seq Type"} = $seqType;

			if ($sample =~ /^ASHPC/)
			{
				$sampleHash{$donor}{$sample}{$seqType}{"Prep Type"} = "flow";
			}
			elsif ($sample =~ /526$/)
			{
				$sampleHash{$donor}{$sample}{$seqType}{"Prep Type"} = "lcm";
			}
			else
			{
				$sampleHash{$donor}{$sample}{$seqType}{"Prep Type"} = "bulk";
			}

			$sampleHash{$donor}{$sample}{$seqType}{"Lanes"} = $sampleHash{$donor}{$sample}{$seqType}{fastqCount} / 2;

			$sampleHash{$donor}{$sample}{$seqType}{oldestLane} =~ s/(..)(..)(..)/20$1-$2-$3/;
			$sampleHash{$donor}{$sample}{$seqType}{newestLane} =~ s/(..)(..)(..)/20$1-$2-$3/;

			$sampleHash{$donor}{$sample}{$seqType}{"Oldest Lane"} = $sampleHash{$donor}{$sample}{$seqType}{oldestLane};
			$sampleHash{$donor}{$sample}{$seqType}{"Newest Lane"} = $sampleHash{$donor}{$sample}{$seqType}{newestLane};


			$sample =~ s/(.*?_.*?_.*?_.*?)_.*/$1/;		# trim group IDs
			$sampleToDonor{$sample} = $donor;
		}
	}
}


# get external names and institutes from LIMS
my %externalIDs;
my %institutes;
my $name;

warn "\nGetting external IDs from LIMS\n";
my $externalList = `export MODULEPATH=/oicr/local/Modules/versions:/oicr/local/Modules/3.2.7/modulefiles:/oicr/local/Modules/modulefiles:/oicr/local/analysis/Modules/modulefiles; eval \`/oicr/local/Modules/3.2.7/bin/modulecmd bash load spb-lims\`; spb-lims-library-attributes --lib PCSI --attr "External Name"; spb-lims-library-attributes --lib ASHPC --attr "External Name"; spb-lims-library-attributes --lib RAMP --attr "External Name"`;

warn "\nGetting external institutess from LIMS\n";
my $instituteList = `export MODULEPATH=/oicr/local/Modules/versions:/oicr/local/Modules/3.2.7/modulefiles:/oicr/local/Modules/modulefiles:/oicr/local/analysis/Modules/modulefiles; eval \`/oicr/local/Modules/3.2.7/bin/modulecmd bash load spb-lims\`; spb-lims-library-attributes --lib PCSI --attr "Institute"; spb-lims-library-attributes --lib ASHPC --attr "Institute"; spb-lims-library-attributes --lib RAMP --attr "Institute"`;

for $l (split(/\n/, $externalList))
{
	if ($l =~ /(.*),"(.*)"/)
	{
		$samp = $1;
		$name = $2;

		$samp =~ s/(.*?_.*?_.*?_.*?)_.*/$1/;

		if (exists $sampleToDonor{$samp})
		{
			unless (exists $externalIDs{$sampleToDonor{$samp}})
			{
				$externalIDs{$sampleToDonor{$samp}} = "";
			}
			$externalIDs{$sampleToDonor{$samp}} .= "$name,";
		}
	}
}

for $l (split(/\n/, $instituteList))
{
	if ($l =~ /(.*),"(.*)"/)
	{
		$samp = $1;
		$name = $2;
		$samp =~ s/(.*?_.*?_.*?_.*?)_.*/$1/;
		if (exists $sampleToDonor{$samp})
		{
			unless (exists $institutes{$sampleToDonor{$samp}})
			{
				$institutes{$sampleToDonor{$samp}} = "";
			}
			$institutes{$sampleToDonor{$samp}} .= "$name,";
		}
	}
}

for my $donor (keys %sampleHash)
{
	for my $sample (sort keys %{ $sampleHash{$donor} })
	{
		for $seqType (sort keys %{ $sampleHash{$donor}{$sample} })
		{
			if (exists $externalIDs{$donor})
			{
				$externalIDs{$donor} =~ s/,$//;
				$sampleHash{$donor}{$sample}{$seqType}{"External ID"} = $externalIDs{$donor};
			}

			if (exists $institutes{$donor})
			{
				$institutes{$donor} =~ s/,$//;
				$sampleHash{$donor}{$sample}{$seqType}{"Institute"} = $institutes{$donor};
			}
		}
	}
}


my $hasR;
my $numSamps;

if ($onlyMatched == 1)
{
	# prune R or PX sets
	for my $donor (keys %sampleHash)
	{
		$hasR = 0;
		$numSamps = 0;
	
		for my $sample (sort keys %{ $sampleHash{$donor} })
		{
			$numSamps++;
			if ($sample =~ /.*R$/)
			{
				$hasR = 1;
			}
		}
	
		if ($numSamps < 2)
		{
			delete $sampleHash{$donor};
		}
		elsif ($hasR == 0)
		{
			delete $sampleHash{$donor};
		}
	
	}
}

# get mouse contamination rate from xenome
warn "\nGetting mouse contamination from xenome log\n";

my $xenomeLogs;
my @xenomeLogs;

my $sum;
my $count;

for my $donor (sort keys %sampleHash)
{
	for my $sample (sort keys %{ $sampleHash{$donor} })
	{
		for $seqType (keys %{ $sampleHash{$donor}{$sample} })
		{
			if (($sample =~ /.*X$/) or ($sample =~ /.*X_526$/))
			{
					$sum = 0;
					$count = 0;
		
					if (-d "$donor/$sample/$seqType/$xenomePath/")
					{
						$xenomeLogs = `tail -n 5 $donor/$sample/$seqType/$xenomePath/*.log | grep mouse`;
						chomp $xenomeLogs;
						@xenomeLogs = split(/\n/, $xenomeLogs);
			
							for my $xenome (@xenomeLogs)
						{
						if ($xenome =~ /.*\t(.*)\tmouse/)
							{
								$sum += $1;
								$count++;
							}
							else
							{
								warn "Couldn't parse $xenome from $donor $sample\n";
							}
						}
						if ($count > 0)
						{
							$sampleHash{$donor}{$sample}{$seqType}{"Mouse"} = ($sum/$count) / 100;
						}
						else
						{
							$sampleHash{$donor}{$sample}{$seqType}{"Mouse"} = "NA";
						}
					}
					else
					{
						$sampleHash{$donor}{$sample}{$seqType}{"Mouse"} = "NA";
					}
				}
			else
			{
				$sampleHash{$donor}{$sample}{$seqType}{"Mouse"} = "";
			}
		}
	}
}


# get deep KRAS
warn "\nGetting deep KRAS cellularity estimate\n";


my $krasBams;
my @krasBams;
my $newestDate;
my $newestBam;
my $bestFreq;
my $bestCall;

for my $donor (sort keys %sampleHash)
{
	for my $sample (sort keys %{ $sampleHash{$donor} })
	{
		for $seqType (sort keys %{ $sampleHash{$donor}{$sample} })
		{
			if (-d "$donor/$sample/$deepKRASpath/")
			{
				$krasBams = `ls $donor/$sample/$deepKRASpath/*.calls`;
				@krasBams = split(/\n/, $krasBams);
				$newestDate = 0;
				for my $bam (@krasBams)
				{
					$bam =~ /.*\/R_(20..)_(..)_(..)_.*/;
					$date = "$1$2$3";
					
					if ($date > $newestDate)
					{
						$newestBam = $bam;
						$newestDate = $date;
					}
				}
	
				$bestFreq = 0;
				$bestCall = "";
				open (BAM, "<$newestBam") or die "Couldn't open $newestBam\n";
				while ($l = <BAM>)
				{
					chomp $l;
					@f = split(/\t/, $l);
	
					if (($f[1] == 25398284) or ($f[1] == 25398285) or ($f[1] == 25380275))
					{
						if ($f[5] > $bestFreq)
						{
							$bestFreq = $f[5];
							$bestCall = "$f[0]:$f[1] $f[2]>$f[3]";
						}
					}
				}
				close BAM;
	
	
				$sampleHash{$donor}{$sample}{$seqType}{"Deep KRAS Variant"} = $bestCall;
				$sampleHash{$donor}{$sample}{$seqType}{"Deep KRAS"} = $bestFreq;
	
	
			}
		}
	}
}

# get doc from json 
warn "\nGetting coverage metric from JSON files\n";

my $jsonHash;

my ($percentOnTarget, $estimatedYield, $estimatedCoverage);

for my $aligner (@aligners)
{
	for my $donor (sort keys %sampleHash)
	{
		for my $sample (sort keys %{ $sampleHash{$donor} })
		{
			for $seqType (sort keys %{ $sampleHash{$donor}{$sample} })
			{
				if (-e "$donor/$sample/$seqType/$aligner/json/$sample.json")
				{
					open (JSON, "<$donor/$sample/$seqType/$aligner/json/$sample.json") or die "Couldn't open $donor/$sample/$seqType/$aligner/json/$sample.json\n";
					if ($l = <JSON>)
					{
						unless ($l =~ /\n/)
						{
							$jsonHash = decode_json($l);
	
							if ($jsonHash->{"mapped reads"} == 0)
							{
								$sampleHash{$donor}{$sample}{$seqType}{"Coverage"} = 0;
							}
							else
							{
								$percentOnTarget = ($jsonHash->{"reads on target"} / $jsonHash->{"mapped reads"}) * 100;
								$estimatedYield = int($jsonHash->{"aligned bases"} * ($percentOnTarget / 100));
								$estimatedCoverage = $estimatedYield / $jsonHash->{"target size"};
								$estimatedCoverage = sprintf "%.2f", $estimatedCoverage;
	
								$sampleHash{$donor}{$sample}{$seqType}{"Coverage"} = $estimatedCoverage;
	
	
								$sampleHash{$donor}{$sample}{$seqType}{"Mismatch Error"} = sprintf("%.2f", $jsonHash->{"mismatch bases"} / $jsonHash->{"aligned bases"} * 100);
								$sampleHash{$donor}{$sample}{$seqType}{"Indel Error"} = sprintf("%.2f", ($jsonHash->{"inserted bases"} + $jsonHash->{"deleted bases"}) / $jsonHash->{"aligned bases"} * 100);
							}
	
						}
						else
						{
							warn "$donor/$sample/$seqType/$aligner/json/$sample.json contains usage\n";
							$sampleHash{$donor}{$sample}{$seqType}{"Coverage"} = "NA";
							$sampleHash{$donor}{$sample}{$seqType}{"Mismatch Error"} = "NA";
							$sampleHash{$donor}{$sample}{$seqType}{"Indel Error"} = "NA";
						}
					}
					else
					{
						warn "$donor/$sample/$seqType/$aligner/json/$sample.json is empty\n";
						$sampleHash{$donor}{$sample}{$seqType}{"Coverage"} = "NA";
						$sampleHash{$donor}{$sample}{$seqType}{"Mismatch Error"} = "NA";
						$sampleHash{$donor}{$sample}{$seqType}{"Indel Error"} = "NA";
					}
					close JSON;
				}
				else
				{
					$sampleHash{$donor}{$sample}{$seqType}{"Coverage"} = "NA";
					$sampleHash{$donor}{$sample}{$seqType}{"Mismatch Error"} = "NA";
					$sampleHash{$donor}{$sample}{$seqType}{"Indel Error"} = "NA";
				}
			}
		}
	}

}



my ($file, $ref, $alt, $tv, $tstv);

# get variants
warn "\nCounting Variants\n";
for my $donor (sort keys %sampleHash)
{
	for my $sample (sort keys %{ $sampleHash{$donor} })
	{
		for $seqType (sort keys %{ $sampleHash{$donor}{$sample} })
		{
			unless ($sample =~ /.*R/)
			{
				#$file = "$donor/$sample/$seqType/$aligner/$mutect/mutect_$sample.vcf";
				$file = "$donor/$sample/$seqType/bwa/0.6.2/final_strelka-mutect/$sample.final.vcf";
				if (-e $file)
				{
					open (VCF, "<$file") or die "Couldn't open $file\n";
	
					$sampleHash{$donor}{$sample}{$seqType}{totalSNV} = 0;
					$sampleHash{$donor}{$sample}{$seqType}{nonsilentSNV} = 0;
					$sampleHash{$donor}{$sample}{$seqType}{transSNV} = 0;
					$sampleHash{$donor}{$sample}{$seqType}{travSNV} = 0;
					$sampleHash{$donor}{$sample}{$seqType}{tstvSNV} = 0;
					$sampleHash{$donor}{$sample}{$seqType}{totalIndel} = 0;
	
					while ($l = <VCF>)
					{
						chomp $l;
						unless ($l =~ /^#/)
						{
							@f = split(/\t/, $l);
							
							if (length($f[3]) == length($f[4]))
							{
								$sampleHash{$donor}{$sample}{$seqType}{totalSNV}++;
								if ($f[7] =~ /ANNOVAR_EXONIC=(.*),/)
								{
									if ($1 ne "synonymous-SNV")
									{
										$sampleHash{$donor}{$sample}{$seqType}{nonsilentSNV}++
									}
								}
	
								$ref = $f[3];
								$alt = $f[4];
	
								if ((($ref eq "A") and ($alt eq "G")) or (($ref eq "G") and ($alt eq "A")) or (($ref eq "C") and ($alt eq "T")) or (($ref eq "T") and ($alt eq "C")))
								{
									$sampleHash{$donor}{$sample}{$seqType}{transSNV}++;
								}
	
							}
							else
							{
								$sampleHash{$donor}{$sample}{$seqType}{totalIndel}++;
							}
						}
					}
	
					$sampleHash{$donor}{$sample}{$seqType}{travSNV} = $sampleHash{$donor}{$sample}{$seqType}{totalSNV} - $sampleHash{$donor}{$sample}{$seqType}{transSNV};
	
					if ($sampleHash{$donor}{$sample}{$seqType}{travSNV} > 0)
					{
						$sampleHash{$donor}{$sample}{$seqType}{tstvSNV} = sprintf("%.2f", $sampleHash{$donor}{$sample}{$seqType}{transSNV} / $sampleHash{$donor}{$sample}{$seqType}{travSNV});
					}
					else
					{
						$sampleHash{$donor}{$sample}{$seqType}{tstvSNV} = "nan";
					}
	
					close VCF;
				}
				else
				{
					$sampleHash{$donor}{$sample}{$seqType}{totalSNV} = "NA";
					$sampleHash{$donor}{$sample}{$seqType}{nonsilentSNV} = "NA";
					$sampleHash{$donor}{$sample}{$seqType}{transSNV} = "NA";
					$sampleHash{$donor}{$sample}{$seqType}{travSNV} = "NA";
					$sampleHash{$donor}{$sample}{$seqType}{tstvSNV} = "NA";
					$sampleHash{$donor}{$sample}{$seqType}{totalIndel} = "NA";
				}
			}
			else
			{
				$sampleHash{$donor}{$sample}{$seqType}{totalSNV} = "";
				$sampleHash{$donor}{$sample}{$seqType}{nonsilentSNV} = "";
				$sampleHash{$donor}{$sample}{$seqType}{transSNV} = "";
				$sampleHash{$donor}{$sample}{$seqType}{travSNV} = "";
				$sampleHash{$donor}{$sample}{$seqType}{tstvSNV} = "";
				$sampleHash{$donor}{$sample}{$seqType}{totalIndel} = "";
			}
		}
	}
}
# finalize
for my $donor (sort keys %sampleHash)
{
	for my $sample (sort keys %{ $sampleHash{$donor} })
	{
		for $seqType (sort keys %{ $sampleHash{$donor}{$sample} })
		{
			if ($sampleHash{$donor}{$sample}{$seqType}{transSNV} ne "")
			{
				$sampleHash{$donor}{$sample}{$seqType}{"SNVs"} = $sampleHash{$donor}{$sample}{$seqType}{totalSNV};
				$sampleHash{$donor}{$sample}{$seqType}{"Non-Silent SNVs"} = $sampleHash{$donor}{$sample}{$seqType}{nonsilentSNV};
				$sampleHash{$donor}{$sample}{$seqType}{"Ts/Tv"} = "$sampleHash{$donor}{$sample}{$seqType}{transSNV}/$sampleHash{$donor}{$sample}{$seqType}{travSNV}";
				$sampleHash{$donor}{$sample}{$seqType}{"Ts/Tv Ratio"} = $sampleHash{$donor}{$sample}{$seqType}{tstvSNV};
				$sampleHash{$donor}{$sample}{$seqType}{"Indels"} = $sampleHash{$donor}{$sample}{$seqType}{totalIndel};
			}
			else
			{
				$sampleHash{$donor}{$sample}{$seqType}{"SNVs"} = $sampleHash{$donor}{$sample}{$seqType}{totalSNV};
				$sampleHash{$donor}{$sample}{$seqType}{"Non-Silent SNVs"} = $sampleHash{$donor}{$sample}{$seqType}{nonsilentSNV};
				$sampleHash{$donor}{$sample}{$seqType}{"Ts/Tv"} = "";
				$sampleHash{$donor}{$sample}{$seqType}{"Ts/Tv Ratio"} = $sampleHash{$donor}{$sample}{$seqType}{tstvSNV};
				$sampleHash{$donor}{$sample}{$seqType}{"Indels"} = $sampleHash{$donor}{$sample}{$seqType}{totalIndel};
			}
		}
	}
}

if ($seqType eq "wgs")
{
	# count structural variants
	warn "\nCounting Structural Variants\n";
	for my $donor (sort keys %sampleHash)
	{
		for my $sample (sort keys %{ $sampleHash{$donor} })
		{
			for $seqType (sort keys %{ $sampleHash{$donor}{$sample} })
			{
				unless ($sample =~ /.*R/)
				{
					$file = "$donor/$sample/$seqType/bwa/0.6.2/$crestDir/$sample.annotatedSV.tsv";
					if (-e $file)
					{
						open (FILE, "<$file") or die "Couldn't open file\n";
		
						while ($l = <FILE>)
						{
							chomp $l;
							@f = split(/\t/,$l);
							$sampleHash{$donor}{$sample}{$seqType}{$f[4]}++;
						}
					}
					else
					{
					}
				}
			}
		}
	}
	
	
	# copy number states
	warn "\nCounting Copy Number States\n";
	
	for my $donor (sort keys %sampleHash)
	{
		for my $sample (sort keys %{ $sampleHash{$donor} })
		{
			for $seqType (sort keys %{ $sampleHash{$donor}{$sample} })
			{
				unless ($sample =~ /.*R/)
				{
					$file = `ls $donor/$sample/$seqType/bwa/0.6.2/$hmmcopyDir/*.cnv_somatic_segments`;
					chomp $file;
					if (-e $file)
					{
						open (FILE, "<$file") or die "Couldn't open file\n";
		
						while ($l = <FILE>)
						{
							chomp $l;
							@f = split(/\t/,$l);
							unless ($f[3] eq "state")
							{
								if ($f[3] == 1)
								{
									$sampleHash{$donor}{$sample}{$seqType}{"Copy State: $f[3]"}++;
								}
								else
								{
									$sampleHash{$donor}{$sample}{$seqType}{"$f[3]"}++;
								}
							}
						}
					}
					else
					{
					}
				}
			}
		}
	}
}





warn "\nPrinting reports!\n\n";

printReport(\%sampleHash, \@sequencingHeaders, "wgs", $wgsSeqFile, \%oldSeqHash, \%seqKeepHeaders, \%seqWarnHeaders);
printReport(\%sampleHash, \@analysisHeaders, "wgs", $wgsAnalysisFile, \%oldAnalysisHash, \%analysisKeepHeaders, \%analysisWarnHeaders);

printReport(\%sampleHash, \@sequencingHeaders, "exome", $exomeSeqFile, \%oldSeqHash, \%seqKeepHeaders, \%seqWarnHeaders);
printReport(\%sampleHash, \@analysisHeaders, "exome", $exomeAnalysisFile, \%oldAnalysisHash, \%analysisKeepHeaders, \%analysisWarnHeaders);

printAllReport(\%sampleHash, \@allHeaders, $allFile);

sub printReport
{
	my $sampleHash = $_[0];
	my $headerRef = $_[1];
	my $type = $_[2];
	my $fileName = $_[3];
	my $oldHash = $_[4];
	my $keepHeaders = $_[5];
	my $warnHeaders = $_[6];

	# merge old and new sample hierarchies
	my %newSampleHash;

	for my $donor (sort keys %{ $oldHash })
	{
		for my $sample (sort keys %{ $oldHash->{$donor} })
		{
			for my $seqType (sort keys %{ $oldHash->{$donor}{$sample} })
			{
				if ($seqType eq "$type")
				{
					for my $header (@{ $headerRef })
					{
						if (exists $oldHash->{$donor}{$sample}{$type}{$header})
						{
							$newSampleHash{$donor}{$sample}{$type}{$header} = $oldHash->{$donor}{$sample}{$type}{$header};
						}
					}
				}
			}
		}
	}

	for my $donor (sort keys %{ $sampleHash })
	{
		for my $sample (sort keys %{ $sampleHash->{$donor} })
		{
			for $seqType (sort keys %{ $sampleHash->{$donor}{$sample} })
			{
				if ($seqType eq "$type")
				{
					for my $header (@{ $headerRef })
					{
						if (exists $sampleHash->{$donor}{$sample}{$type}{$header})
						{
							if (exists $keepHeaders->{$header})
							{
								if ($oldHash->{$donor}{$sample}{$type}{$header} ne "")
								{
									$newSampleHash{$donor}{$sample}{$type}{$header} = $oldHash->{$donor}{$sample}{$type}{$header};
								}
								else
								{
									$newSampleHash{$donor}{$sample}{$type}{$header} = $sampleHash->{$donor}{$sample}{$type}{$header};
								}
							}
							elsif (exists $warnHeaders->{$header})
							{
								if (exists $oldHash->{$donor}{$sample}{$type}{$header})
								{
									unless ("$sampleHash->{$donor}{$sample}{$type}{$header}" eq "$oldHash->{$donor}{$sample}{$type}{$header}")
									{
										$newSampleHash{$donor}{$sample}{$type}{$header} = $sampleHash->{$donor}{$sample}{$type}{$header};
										$newSampleHash{$donor}{$sample}{$type}{"Notes"} .= "; New $header!";
									}
								}
								else
								{
									$newSampleHash{$donor}{$sample}{$type}{$header} = $sampleHash->{$donor}{$sample}{$type}{$header};
								}

							}
							else
							{
								$newSampleHash{$donor}{$sample}{$type}{$header} = $sampleHash->{$donor}{$sample}{$type}{$header};
							}
						}
					}
				}
			}
		}
	}


	my $timestamp = `date`;
	chomp $timestamp;

	open (FILE, ">$fileName") or die "Couldn't open $fileName\n";

	print FILE "This report was generated $timestamp, contents may have been modified since then.\n\n";


	if ($fileName =~ /seq/)		# not an idea way to control this
	{
		# count complete and sequenced donor pairs
		my %complete;
		my %sequenced;
	
		my %tissueStatus;
	
		my $tissue;
		my $paired;

		my $institute;
	
		my @preps = ("bulk","lcm","flow");
	
		for my $prep (@preps)
		{
			for my $donor (sort keys %newSampleHash )
			{
				$institute = "";
				$paired = 0;
				%tissueStatus = ();
				for my $sample (sort keys %{ $newSampleHash{$donor} })
				{
					if (exists $newSampleHash{$donor}{$sample}{$seqType}{"Institute"})
					{
						$institute = $newSampleHash{$donor}{$sample}{$seqType}{"Institute"};
					}
					for $seqType (sort keys %{ $newSampleHash{$donor}{$sample} })
					{
						if ($seqType eq "$type")
						{
							if ($newSampleHash{$donor}{$sample}{$seqType}{"Prep Type"} eq $prep)
							{
								if ($newSampleHash{$donor}{$sample}{$type}{"Sample"} =~ /^.*?_.*?_.*?_(.)/)
								{
									$tissue = $1;
									$tissueStatus{$tissue} = $newSampleHash{$donor}{$sample}{$seqType}{"Status"};
								}
							}
						}
					}
				}
	
				for $tissue (qw/P X C/)
				{
					if ((exists $tissueStatus{R}) and (exists $tissueStatus{$tissue}))
					{
						$paired = 1;
						$sequenced{$prep}{"TN Pair"}{$donor} = 1;
						$sequenced{$prep}{"R$tissue"}{$donor} = 1;
	
						$sequenced{total}{"TN Pair"}{$donor} = 1;
						$sequenced{total}{"R$tissue"}{$donor} = 1;
		
						$sequenced{$prep}{$institute}{$donor} = 1;
						$sequenced{total}{$institute}{$donor} = 1;

						if (($tissueStatus{R} eq "complete") and ($tissueStatus{$tissue} eq "complete"))
						{
							$complete{$prep}{"TN Pair"}{$donor} = 1;
							$complete{$prep}{"R$tissue"}{$donor} = 1;
	
							$complete{total}{"TN Pair"}{$donor} = 1;
							$complete{total}{"R$tissue"}{$donor} = 1;

							$complete{$prep}{$institute}{$donor} = 1;
							$complete{total}{$institute}{$donor} = 1;

#							warn "$donor ($prep) from $institute was counted\n";

						}
					}
				}
				if ($paired == 0)
				{
					if ((exists $tissueStatus{R}) or (exists $tissueStatus{P}) or (exists $tissueStatus{X}) or (exists $tissueStatus{C}))
					{
						$sequenced{$prep}{"Unpaired"}{$donor} = 1;
						$sequenced{total}{"Unpaired"}{$donor} = 1;


					}
				}
			}
		}

		print FILE "\t\tCompleted Donors:\t\t\t\t\t\t\t\t\t\tSequenced Donors:\n";
		print FILE "\t\tTN Pair\tRP\tRX\tRC\t\tUHN\tMayo\tMGH\tSHSC\tUNMC\t\t\tTN Pair\tRP\tRX\tRC\tUnpaired\t\tUHN\tMayo\tMGH\tSHSC\tUNMC\n";

		warn "\t\tCompleted Donors:\t\t\t\t\t\t\t\t\t\tSequenced Donors:\n";
		warn "\t\tTN Pair\tRP\tRX\tRC\t\tUHN\tMayo\tMGH\tSHSC\tUNMC\t\t\tTN Pair\tRP\tRX\tRC\tUnpaired\t\tUHN\tMayo\tMGH\tSHSC\tUNMC\n";

		@preps = ("bulk","lcm","flow","total");

		for my $prep (@preps)
		{
			print FILE "\t$prep";
			print FILE "\t" . scalar keys %{ $complete{$prep}{"TN Pair"} };
			print FILE "\t" . scalar keys %{ $complete{$prep}{"RP"} };
			print FILE "\t" . scalar keys %{ $complete{$prep}{"RX"} };
			print FILE "\t" . scalar keys %{ $complete{$prep}{"RC"} };
			print FILE "\t";
			print FILE "\t" . scalar keys %{ $complete{$prep}{"UHN"} };
			print FILE "\t" . scalar keys %{ $complete{$prep}{"Mayo"} };
			print FILE "\t" . scalar keys %{ $complete{$prep}{"MGH"} };
			print FILE "\t" . scalar keys %{ $complete{$prep}{"SHSC"} };
			print FILE "\t" . scalar keys %{ $complete{$prep}{"UNMC"} };
			print FILE "\t\t";

			print FILE "\t$prep";
			print FILE "\t" . scalar keys %{ $sequenced{$prep}{"TN Pair"} };
			print FILE "\t" . scalar keys %{ $sequenced{$prep}{"RP"} };
			print FILE "\t" . scalar keys %{ $sequenced{$prep}{"RX"} };
			print FILE "\t" . scalar keys %{ $sequenced{$prep}{"RC"} };
			print FILE "\t" . scalar keys %{ $sequenced{$prep}{"Unpaired"} };
			print FILE "\t";
			print FILE "\t" . scalar keys %{ $sequenced{$prep}{"UHN"} };
			print FILE "\t" . scalar keys %{ $sequenced{$prep}{"Mayo"} };
			print FILE "\t" . scalar keys %{ $sequenced{$prep}{"MGH"} };
			print FILE "\t" . scalar keys %{ $sequenced{$prep}{"SHSC"} };
			print FILE "\t" . scalar keys %{ $sequenced{$prep}{"UNMC"} };
			print FILE "\n";

			warn "\t$prep" . "\t" . scalar keys(%{ $complete{$prep}{"TN Pair"} }) . "\t" . scalar keys( %{ $complete{$prep}{"RP"} }) . "\t" . scalar keys( %{ $complete{$prep}{"RX"} }) . "\t" . scalar keys( %{ $complete{$prep}{"RC"} }) . "\t\t" . scalar keys( %{ $complete{$prep}{"UHN"} }) . "\t" . scalar keys( %{ $complete{$prep}{"Mayo"} }) . "\t" . scalar keys( %{ $complete{$prep}{"MGH"} }) . "\t" . scalar keys( %{ $complete{$prep}{"SHSC"} }) . "\t" . scalar keys( %{ $complete{$prep}{"UNMC"} }) . "\t\t" . "\t$prep" . "\t" . scalar keys( %{ $sequenced{$prep}{"TN Pair"} }) . "\t" . scalar keys( %{ $sequenced{$prep}{"RP"} }) . "\t" . scalar keys( %{ $sequenced{$prep}{"RX"} }) . "\t" . scalar keys( %{ $sequenced{$prep}{"RC"} }) . "\t" . scalar keys( %{ $sequenced{$prep}{"Unpaired"} }) . "\t\t" . scalar keys( %{ $sequenced{$prep}{"UHN"} }) . "\t" . scalar keys( %{ $sequenced{$prep}{"Mayo"} }) . "\t" . scalar keys( %{ $sequenced{$prep}{"MGH"} }) . "\t" . scalar keys( %{ $sequenced{$prep}{"SHSC"} }) . "\t" . scalar keys( %{ $sequenced{$prep}{"UNMC"} }) . "\n";
		}

		print FILE "\n\n";

	}






	for my $header (@{ $headerRef })
	{
		print FILE "$header\t";
	}
	print FILE "\n\n";
	
	
	my $printed;
	for my $donor (sort keys %newSampleHash )
	{
		$printed = 0;
		for my $sample (sort keys %{ $newSampleHash{$donor} })
		{
			for $seqType (sort keys %{ $newSampleHash{$donor}{$sample} })
			{
				if ($seqType eq "$type")
				{
					$printed = 1;
					for my $header (@{ $headerRef })
					{
						if (exists $newSampleHash{$donor}{$sample}{$seqType}{$header})
						{
							print FILE "$newSampleHash{$donor}{$sample}{$seqType}{$header}\t";
						}
						else
						{
							print FILE "\t";
						}
					}
					print FILE "\n";
				}
			}
	
		}
		if ($printed == 1)
		{
			print FILE "\n";
		}
	}
}


sub printAllReport
{
	my $sampleHash = $_[0];
	my $headerRef = $_[1];
	my $fileName = $_[2];

	my $timestamp = `date`;
	chomp $timestamp;

	open (FILE, ">$fileName") or die "Couldn't open $fileName\n";

	print FILE "This report was generated $timestamp, contents may have been modified since then.\n\n";

	for my $header (@{ $headerRef })
	{
		print FILE "$header\t";
	}
	print FILE "\n";
	
	
	my $printed;
	for my $donor (sort keys %{ $sampleHash })
	{
		$printed = 0;
		for my $sample (sort keys %{ $sampleHash->{$donor} })
		{
			for $seqType (sort keys %{ $sampleHash->{$donor}{$sample} })
			{
				if ($seqType eq "exome")
				{
					$printed = 1;
					for my $header (@{ $headerRef })
					{
						if (exists $sampleHash->{$donor}{$sample}{$seqType}{$header})
						{
							print FILE "$sampleHash->{$donor}{$sample}{$seqType}{$header}\t";
						}
						else
						{
							print FILE "\t";
						}
					}
					print FILE "\n";
				}
			}
	
		}
		for my $sample (sort keys %{ $sampleHash->{$donor} })
		{
			for $seqType (sort keys %{ $sampleHash->{$donor}{$sample} })
			{
				if ($seqType eq "wgs")
				{
					$printed = 1;
					for my $header (@{ $headerRef })
					{
						if (exists $sampleHash->{$donor}{$sample}{$seqType}{$header})
						{
							print FILE "$sampleHash->{$donor}{$sample}{$seqType}{$header}\t";
						}
						else
						{
							print FILE "\t";
						}
					}
					print FILE "\n";
				}
			}
		}

		if ($printed == 1)
		{
			print FILE "\n";
		}
	}
}



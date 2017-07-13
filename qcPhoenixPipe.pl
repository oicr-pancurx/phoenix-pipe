#!/usr/bin/perl

# takes path up to seq type as input
# eg: PCSI0666/PCSI_0666_Lv_M_526/wgs

use strict;
use warnings;

my $path = shift;
my $type = shift;		# either tumour or normal

my $wc;
my $ls;
my $warnFound;
my $failFound;

my %doTest = (
	"lanes" => 1,
	"collapsed" => 1,
	"json" => 1,
	"strelka/v1.0.7" => 1,
	"mutect/1.1.4" => 1,
	"final_strelka-mutect" => 1,
	"gatk/1.3.16" => 1,
	"final_gatk-germline" => 1,
	"crest/alpha" => 1,
	"delly/0.5.5" => 1,
	"final_crest-delly" => 1,
	"HMMcopy/0.1.1" => 1,
	"celluloid/v11.2" => 1,
	"polysolver/1.0" => 1,
	"netMHC/pan-2.8/polysolver/1.0" => 1,
	"integration" => 1,
);


my $dir;

########## test lane bams ##########
$dir = "lanes";
if ($doTest{$dir} == 1)
{
	$warnFound = "";
	$failFound = "";

	# directory should exist
	$failFound .= checkDir("$path/$dir/", $dir);

	# bwa log files should have reached completion
	$failFound .= checkLogs("$path/$dir/", $dir);
	
	# there should be a bam file with a header for every fastq file

	# the number of reads in the bam file should match the number of reads in the fastq files
	
	
	print $warnFound;
	print $failFound;
}

########## test merged bam ##########
$dir = "collapsed";
if ($doTest{$dir} == 1)
{
	$warnFound = "";
	$failFound = "";

	# directory should exist
	$failFound .= checkDir("$path/$dir/", $dir);

	# log files should have reached completion
	$failFound .= checkLogs("$path/$dir/", $dir);

	# bam file should exist and have non-zero size
	$ls = `ls $path/$dir/*.bam`;
	$failFound .= checkFile($ls, $dir);
	
	# the number of reads in the lane bam files (and the fastqs) should match the number of reads in the merged/collapsed bam
	
	print $warnFound;
	print $failFound;
}



########## test json files ##########
$dir = "json";
if ($doTest{$dir} == 1)
{
	$warnFound = "";
	$failFound = "";

	# directory should exist
	$failFound .= checkDir("$path/$dir/", $dir);

	# log files should have reached completion
	$failFound .= checkLogs("$path/$dir/", $dir);

	# json files should exist for each lane and be non-empty
	$ls = `ls $path/$dir/*.json`;
	$failFound .= checkFile($ls, $dir);
	
	print $warnFound;
	print $failFound;
}


if ($type eq "tumour")
{

	########## test stelka ##########
	$dir = "strelka/v1.0.7";
	if ($doTest{$dir} == 1)
	{
		$warnFound = "";
		$failFound = "";
	
		# directory should exist
		$failFound .= checkDir("$path/$dir/", $dir);
	
		# log files should have reached completion
		$failFound .= checkLogs("$path/$dir/", $dir);		# there are other .log files in the annovar diretory that are not from phoenixPipe
		$failFound .= checkLogs("$path/$dir/final/", $dir);
	
		# annovar log file should have reached completion
		$failFound .= checkLogs("$path/$dir/annovar/2013-06-21/*phoenix", $dir);
	
		# annovar output should exist and not be empty
		$ls = `ls $path/$dir/annovar/2013-06-21/*.annovar.variant_function $path/$dir/annovar/2013-06-21/*.annovar.exonic_variant_function`;
		$failFound .= checkFile($ls, $dir);
	
		# funseq log file should have reached completion
		$failFound .= checkLogs("$path/$dir/funseq/0.1/", $dir);
	
		# funseq vcf should exist and not be empty
		$ls = `ls $path/$dir/funseq/0.1/out/*.FunSEQ.vcf`;
		$failFound .= checkFile($ls, $dir);
		
		# final log file should have reached completion
		$failFound .= checkLogs("$path/$dir/final/", $dir);
	
		# vcf files should not be empty
		$ls = `ls $path/$dir/*.vcf $path/$dir/final/*final.vcf`;
		$failFound .= checkFile($ls, $dir);
	
		# vcf files should have variants from each chromosome
		$ls = `ls $path/$dir/*.vcf $path/$dir/final/*final.vcf`;
		$warnFound .= checkVcfForChroms($ls, $dir);
		
		print $warnFound;
		print $failFound;
	}
	
	########## test mutect ##########
	$dir = "mutect/1.1.4";
	if ($doTest{$dir} == 1)
	{
		$warnFound = "";
		$failFound = "";
	
		# directory should exist
		$failFound .= checkDir("$path/$dir/", $dir);
	
		# log files should have reached completion
		$failFound .= checkLogs("$path/$dir/", $dir);
		$failFound .= checkLogs("$path/$dir/tempsplit/", $dir);
	
		# annovar log file should have reached completion
		$failFound .= checkLogs("$path/$dir/annovar/2013-06-21/*phoenix", $dir);
	
		# annovar output should exist and not be empty
		$ls = `ls $path/$dir/annovar/2013-06-21/*.annovar.variant_function $path/$dir/annovar/2013-06-21/*.annovar.exonic_variant_function`;
		$failFound .= checkFile($ls, $dir);
	
		# funseq log file should have reached completion
		$failFound .= checkLogs("$path/$dir/funseq/0.1/", $dir);
	
		# funseq vcf should exist and not be empty
		$ls = `ls $path/$dir/funseq/0.1/out/*.FunSEQ.vcf`;
		$failFound .= checkFile($ls, $dir);
	
		# final vcf file should have snvs from each chromosome
		$ls = `ls $path/$dir/*.vcf $path/$dir/final/*.vcf`;
		$warnFound .= checkVcfForChroms($ls, $dir);
		
		print $warnFound;
		print $failFound;
	}
	
	########## test final_strelka-mutect ##########
	$dir = "final_strelka-mutect";
	if ($doTest{$dir} == 1)
	{
		$warnFound = "";
		$failFound = "";
	
		# directory should exist
		$failFound .= checkDir("$path/$dir/", $dir);
	
		# log files should have reached completion
		$failFound .= checkLogs("$path/$dir/", $dir);
	
		# final vcf file should have snvs and indels from each chromosome
		$ls = `ls $path/$dir/*.vcf`;
		$warnFound .= checkVcfForChroms($ls, $dir);
		
		print $warnFound;
		print $failFound;
	}
	
	########## test gatk ##########
	$dir = "gatk/1.3.16";
	if ($doTest{$dir} == 1)
	{
		$warnFound = "";
		$failFound = "";
	
		# directory should exist
		$failFound .= checkDir("$path/$dir/", $dir);
	
		# log files should have reached completion
		$failFound .= checkLogs("$path/$dir/", $dir);
	
		# gatk logs should have no error messages (grep -i "error message" */*/wgs/bwa/*/gatk/*/gatkLog/*)
		$failFound .= checkFilesForMessage("$path/$dir/gatkLog/*", $dir, "-i \"error message\"");
	
		
		# annovar log file should have reached completion
		$failFound .= checkLogs("$path/$dir/annovar/2013-06-21/*phoenix", $dir);
	
		# annovar output should exist and not be empty
		$ls = `ls $path/$dir/annovar/2013-06-21/*.annovar.variant_function $path/$dir/annovar/2013-06-21/*.annovar.exonic_variant_function`;
		$failFound .= checkFile($ls, $dir);
		
		# funseq log file should have reached completion
		$failFound .= checkLogs("$path/$dir/funseq/0.1/", $dir);
	
		# funseq vcf should exist and not be empty
		$ls = `ls $path/$dir/funseq/0.1/out/*.FunSEQ.vcf`;
		$failFound .= checkFile($ls, $dir);
		
		# final vcf file should have snvs and indels from each chromosome
		$ls = `ls $path/$dir/*.vcf $path/$dir/final/*.vcf`;
		$warnFound .= checkVcfForChroms($ls, $dir);
		
		print $warnFound;
		print $failFound;
	}
	
	########## test final_gatk-germline ##########
	$dir = "final_gatk-germline";
	if ($doTest{$dir} == 1)
	{
		$warnFound = "";
		$failFound = "";
	
		# directory should exist
		$failFound .= checkDir("$path/$dir/", $dir);
	
		# log files should have reached completion
		$failFound .= checkLogs("$path/$dir/", $dir);
	
		# final vcf file should have snvs and indels from each chromosome
		$ls = `ls $path/$dir/*.vcf`;
		$warnFound .= checkVcfForChroms($ls, $dir);
		
		print $warnFound;
		print $failFound;
	}
	
	########## test crest ##########
	$dir = "crest/alpha";
	if ($doTest{$dir} == 1)
	{
		$warnFound = "";
		$failFound = "";
	
		# directory should exist
		$failFound .= checkDir("$path/$dir/", $dir);
	
		# log files should have reached completion
		$failFound .= checkLogs("$path/$dir/extractSClip/", $dir);
		$failFound .= checkLogs("$path/$dir/tempSplit/", $dir);
		$failFound .= checkLogs("$path/$dir/filter/", $dir);
	
		# log files should not have blat server down messages
		$failFound .= checkFilesForMessage("$path/$dir/tempSplit/*.log", $dir, "Sorry");
		
		# output shouldn't be empty
		$ls = `ls $path/$dir/*.predSV.txt $path/$dir/filter/*.predSV.txt`;
		$failFound .= checkFile($ls, $dir);
	
		
		print $warnFound;
		print $failFound;
	}
	
	########## test delly ##########
	$dir = "delly/0.5.5";
	if ($doTest{$dir} == 1)
	{
		$warnFound = "";
		$failFound = "";
	
		# directory should exist
		$failFound .= checkDir("$path/$dir/", $dir);
	
		# log files should have reached completion
		$failFound .= checkLogs("$path/$dir/", $dir);
	
		# output should contain SVs from each chr
		$ls = `ls $path/$dir/*.vcf`;
		$failFound .= checkVcfForChroms($ls, $dir);
	
		# filtered output shouldn't be empty
		$ls = `ls $path//$dir/filter/*.vcf`;
		$failFound .= checkFile($ls, $dir);
		
		print $warnFound;
		print $failFound;
	}
	
	########## test final_crest-delly ##########
	$dir = "final_crest-delly";
	if ($doTest{$dir} == 1)
	{
		$warnFound = "";
		$failFound = "";
	
		# directory should exist
		$failFound .= checkDir("$path/$dir/", $dir);
	
		# log files should have reached completion
		$failFound .= checkLogs("$path/$dir/", $dir);
	
		#output shouldn't be empty
		$ls = `ls $path/$dir/*.tsv`;
		$failFound .= checkFile($ls, $dir);
		
		print $warnFound;
		print $failFound;
	}
	
	########## test HMMcopy ##########
	$dir = "HMMcopy/0.1.1";
	if ($doTest{$dir} == 1)
	{
		$warnFound = "";
		$failFound = "";
	
		# directory should exist
		$failFound .= checkDir("$path/$dir/", $dir);
	
		# log files should have reached completion
		$failFound .= checkLogs("$path/$dir/", $dir);
	
		# wig files for tumour and normal should be non-empty
		$ls = `ls $path/$dir/*.wig`;
		$failFound .= checkFile($ls, $dir);
		
		# tumour segments file should have a segment for each chromosome
		$ls = `ls $path/$dir/*.cnv_somatic_segments`;
		$failFound .= checkVcfForChroms($ls, $dir);
		
		print $warnFound;
		print $failFound;
	}
	
	########## test celluloid ##########
	$dir = "celluloid/v11.2";
	if ($doTest{$dir} == 1)
	{
		$warnFound = "";
		$failFound = "";
	
		# directory should exist
		$failFound .= checkDir("$path/$dir/", $dir);
	
		# log files should have reached completion
		$failFound .= checkLogs("$path/$dir/", $dir);
	
		# solution directory should exists
		$failFound .= checkDir("$path/$dir/solution/", $dir);
		
		# param file should be non-empty
		
		# segment file should have a segment for each chromosome
		
		print $warnFound;
		print $failFound;
	}
	
	########## test polysolver ##########
	$dir = "polysolver/1.0";
	if ($doTest{$dir} == 1)
	{
		$warnFound = "";
		$failFound = "";
	
		# directory should exist
		$failFound .= checkDir("$path/$dir/", $dir);
	
		# log files should have reached completion
		$failFound .= checkLogs("$path/$dir/", $dir);
	
		# log files should not have issues locating the perl module
		$failFound .= checkFilesForMessage("$path/$dir/*.log", $dir, "Can't locate Bio/DB/Sam.pm");

		# hla type file should be non-empty and not contain the default results
		$warnFound .= checkFilesDoNotContain("$path/$dir/*.HLA.txt", $dir, "A*01:01\nB*07:02\nC*01:02");

		
		print $warnFound;
		print $failFound;
	}
	
	########## test netMHC ##########
	$dir = "netMHC/pan-2.8/polysolver/1.0";
	if ($doTest{$dir} == 1)
	{
		$warnFound = "";
		$failFound = "";
	
		# directory should exist
		$failFound .= checkDir("$path/$dir/", $dir);
	
		# log files should have reached completion
		$failFound .= checkLogs("$path/$dir/", $dir);
	
		# txt for each hla type should contain results
		
		print $warnFound;
		print $failFound;
	}
	
	
}


sub checkDir
{
	my $dir = shift;
	my $id = shift;

	if (-d $dir)
	{
		return "";
	}
	else
	{
		return "$id FAIL: the directory $path/lanes/ does not exist\n";
	}
}

sub checkLogs
{
	my $dir = shift;
	my $id = shift;

	my $ls = `ls $dir*.log`;
	chomp $ls;

	my $tail;
	my $message = "";

	my $doneMessageDate = 1472702400;	# Sept 1st, 2016

	for my $log (split(/\n/, $ls))
	{
		$tail = `tail -n 1 $log`;
		chomp $tail;

		if ((stat($log))[9] > $doneMessageDate)
		{
			unless ($tail =~ /phoenixPipe\/.*done/)
			{
				$message .= "$id FAIL: log file $log did not end with phoenixPipe/.*-done\n";
			}
		}
	}

	return $message;

}

sub checkFilesForMessage
{
	my $dir = shift;
	my $id = shift;
	my $phrase = shift;

	my $grep = `grep $phrase $dir | wc -l`;
	chomp $grep;

	my $message = "";

	if ($grep > 0)
	{
		$message .= "$id FAIL: log files in $dir contained the phrase $phrase\n";
	}

	return $message;

}

sub checkFile
{
	my $ls = shift;
	my $id = shift;

	chomp $ls;

	my $message = "";

	for my $file (split(/\n/, $ls))
	{
		unless (-e $file)
		{
			$message .= "$id FAIL: file $file does not exist\n";
		}
		else
		{
			unless (-s $file)
			{
				$message .= "$id FAIL: file $file is empty\n";
			}
		}
	}

	return $message;
}


sub checkVcfForChroms
{
	my $ls = shift;
	my $id = shift;

	chomp $ls;

	my $uniq;
	my @chroms = qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY/;
	my %inVcf;
	my $message = "";

	for my $file (split(/\n/, $ls))
	{
		$uniq = `grep -v ^# $file | cut -f 1 | uniq`;
		chomp $uniq;

		for my $l (split(/\n/, $uniq))
		{
			$inVcf{$l}++;
		}

		for my $chr (@chroms)
		{
			unless (exists $inVcf{$chr})
			{
				$message .= "$id WARN: no variants on $chr in $file\n";
			}
		}
	}

	return $message;
}




sub checkFilesDoNotContain
{
	my $fileList = shift;
	my $id = shift;
	my $badContent = shift;

	warn "ls $fileList\n";

	my $files = `ls $fileList`;
	chomp $files;

	my $message = "";
	my $content;

	my $sanitizedContent = $badContent;
	$sanitizedContent =~ s/\n/;/g;

	for my $file (split/\n/, $files)
	{
		$content = `cat $file`;
		chomp $content;

		if ($content eq $badContent)
		{
			$message .= "$id WARN: txt file in $fileList matched the undesired content: $sanitizedContent\n";
		}
	}

	return $message;

}


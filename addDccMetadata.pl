#!/usr/bin/perl

use strict;
use warnings;

use JSON;

my $bamModule = $ARGV[0];
my $callModule = $ARGV[1];
my $seqType = $ARGV[2];
my $pathToNormalJSON = $ARGV[3];
my $pathToTumourJSON = $ARGV[4];

my $normalSample = $ARGV[5];
my $tumourSample = $ARGV[6];
my $donor = $ARGV[7];

my $l;
my @f;


my $analyzed_sample_id = $tumourSample;
my $matched_sample_id = $normalSample;
my $donor_id = $donor;
$donor_id =~ s/(\d+)/_$1/;		# putting the underscore back in


my $assembly_version = "1,term=\"GRCh37\"";
my $platform = "60,term=\"Illumina HiSeq\"";

my $experimental_protocol = "";
if ($seqType eq "exome")
{
	$experimental_protocol = "Exome capture with Agilent SureSelectHumanAllExonV4 - protocol pending http://oicr.on.ca/";
}
elsif ($seqType eq "wgs")
{
	if ($tumourSample =~ /526$/)
	{
		if ($callModule =~ /gatk/)
		{
			# germline is from reference, which is never LCM
			$experimental_protocol = "Whole genome shotgun sequencing on bulk material - protocol pending http://oicr.on.ca/";
		}
		else
		{
			$experimental_protocol = "Whole genome shotgun sequencing on laser capture microdissected material - protocol pending http://oicr.on.ca/";
		}
	}
	elsif ($tumourSample =~ /ASHPC/)
	{
		$experimental_protocol = "Whole genome shotgun sequencing on flow sorted material - protocol pending http://oicr.on.ca/";
	}
	else
	{
		$experimental_protocol = "Whole genome shotgun sequencing on bulk material - protocol pending http://oicr.on.ca/";
	}
}

my $base_calling_algorithm = "CASAVA 1.8.2 http://support.illumina.com/sequencing/sequencing_software/casava.ilmn";

my $alignment_algorithm = "";
if ($bamModule =~ /bwa\/0.6.2/)
{
	$alignment_algorithm = "BWA 0.6.2 http://bio-bwa.sourceforge.net/";
}
elsif ($bamModule =~ /novocraft\/3.00.05/)
{
	$alignment_algorithm = "Novoalign 3.00.05 http://novocraft.com/main/index.php";
}

my $variation_calling_algorithm = "";
if ($callModule =~ /strelka\/v1.0.7/)
{
	$variation_calling_algorithm = "Strelka 1.0.7 https://sites.google.com/site/strelkasomaticvariantcaller/";
}
elsif ($callModule =~ /mutect\/1.1.4/)
{
	$variation_calling_algorithm = "MuTect 1.1.4 http://www.broadinstitute.org/cancer/cga/mutect";
}
elsif ($callModule =~ /gatk\/1\.3\.16/)
{
	$variation_calling_algorithm = "GATK 1.3.16 http://www.broadinstitute.org/gatk/";
}

my $other_analysis_algorithm = "";
if ($bamModule =~ /gatk\/2.4.9/)
{
	$other_analysis_algorithm = "Picard 1.90 http://picard.sourceforge.net/,GATK 2.4.9 http://www.broadinstitute.org/gatk/,ANNOVAR http://www.openbioinformatics.org/annovar/";
}
else
{
	$other_analysis_algorithm = "Picard 1.90 http://picard.sourceforge.net/,ANNOVAR http://www.openbioinformatics.org/annovar/";
}



my $sequencing_strategy = "";
if ($seqType eq "exome")
{
	$sequencing_strategy = "3,term=\"WXS\"";
}
elsif ($seqType eq "wgs")
{
	$sequencing_strategy = "1,term=\"WGS\"";
}


my $analyzed_seq_coverage = "";
my $matched_seq_coverage = "";
my %jsonHash;
if (-e $pathToNormalJSON)
{
	open (FILE, $pathToNormalJSON) or warn "Couldn't open $pathToNormalJSON.\n";
	if ($l = <FILE>)
	{
		$jsonHash{j} = decode_json($l);

		$matched_seq_coverage = ($jsonHash{j}{"aligned bases"} * ($jsonHash{j}{"reads on target"} / $jsonHash{j}{"mapped reads"}) ) / $jsonHash{j}{"target size"};
	}
	else
	{
		warn "Empty JSON file!\n";
	}
}
if (-e $pathToTumourJSON)
{
	open (FILE, $pathToTumourJSON) or warn "Couldn't open $pathToTumourJSON.\n";
	if ($l = <FILE>)
	{
		$jsonHash{j} = decode_json($l);

		$analyzed_seq_coverage = ($jsonHash{j}{"aligned bases"} * ($jsonHash{j}{"reads on target"} / $jsonHash{j}{"mapped reads"}) ) / $jsonHash{j}{"target size"};
	}
	else
	{
		warn "Empty JSON file!\n";
	}
}







while ($l = <STDIN>)
{
	if ($l =~ /^#CHROM/)
	{
		print "##DCC=<donor_id=\"$donor_id\">\n";
		print "##DCC=<analyzed_sample_id=\"$analyzed_sample_id\">\n";
		print "##DCC=<matched_sample_id=\"$matched_sample_id\">\n";
		print "##DCC=<assembly_version=$assembly_version>\n";
		print "##DCC=<platform=$platform>\n";
		print "##DCC=<experimental_protocol=\"$experimental_protocol\">\n";
		print "##DCC=<base_calling_algorithm=\"$base_calling_algorithm\">\n";
		print "##DCC=<alignment_algorithm=\"$alignment_algorithm\">\n";
		print "##DCC=<variation_calling_algorithm=\"$variation_calling_algorithm\">\n";
		print "##DCC=<other_analysis_algorithm=\"$other_analysis_algorithm\">\n";
		print "##DCC=<sequencing_strategy=$sequencing_strategy>\n";
		print "##DCC=<analyzed_seq_coverage=$analyzed_seq_coverage>\n";
		print "##DCC=<matched_seq_coverage=$matched_seq_coverage>\n";
		


		print $l;
	}
	else
	{
		print $l;
	}

}


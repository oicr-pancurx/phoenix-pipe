#!/usr/bin/perl

use strict;
use warnings;

my $date = `date`;
chomp $date;

# ( echo "<html><head></head><body>"; for i in `ls */*/wgs/*/*/final*/*final.vcf.html; ls */*/exome/*/*/final*/*final.vcf.html`; do DONOR=`echo $i | cut -f 1 -d "/"`; SAMPLE=`echo $i | cut -f 2 -d "/"`; TYPE=`echo $i | cut -f 3 -d "/"`; echo "<a href=\"$i\">$DONOR $SAMPLE $TYPE</a><br>"; done; echo "</body></html>") > vcfReport.html


open (HTML, ">vcfReport.html") or die "Couldn't open >vcfReport.html\n";


print HTML "<html><head>\n";
print HTML "<title>VCF Report</title>\n";
print HTML "<script src=\"./sorttable.js\"></script>\n";
print HTML "<style type=\"text/css\">\nth, td {\n  padding: 3px !important;\n}\ntable\n{\nborder-collapse:collapse;\n}\n/* Sortable tables */\ntable.sortable thead {\n\tbackground-color:#eee;\n\tcolor:#000000;\n\tfont-weight: bold;\n\tcursor: default;\n}\n</style>\n";
print HTML "</head><body>\n";

print HTML "<p>VCF Report generated $date.</p>\n";
print HTML "<p><a href=\"./jsonReport/wgs/PCSI_report.html\">WGS Coverage QC Report</a></p>\n";
print HTML "<p><a href=\"./jsonReport/exome/PCSI_report.html\">Exome Coverage QC Report</a></p>\n";

# Donor Sample SeqType LastModified


print HTML "<table border=\"1\" class=\"sortable\">\n<thead>\n<tr>\n";
print HTML "<th>Donor</th><th>Sample</th><th>Sequencing Type</th><th>Normal Coverage (x)</th><th>Tumour Coverage (x)</th><th>Somatic SNVs</th><th>Somatic Indels</th><th>Last Modified</th>\n</tr>\n</thead>\n<tbody>\n";


my ($sample, $donor, $seqType, $timestamp, $lines, $variantLines, $nCoverage, $tCoverage, $numSNVs, $numIndels);
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst);

my @covs;

my $fileString = `ls */*/*/*/*/final*/*final.vcf.html`;
chomp $fileString;

for my $file (split(/\n/, $fileString))
{
	if ($file =~ /^(.*?)\/(.*?)\/(.*?)\/.*/)
	{
		$donor = $1;
		$sample = $2;
		$seqType = $3;
	}
	else
	{
		die "Couldn't parse $file\n";
	}

	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime((stat($file))[9]);
	$year += 1900;
	$mon += 1;
	if ($mon < 10)
	{
		$mon = "0$mon";
	}
	if ($mday < 10)
	{
		$mday = "0$mday";
	}
	$timestamp = "$year$mon$mday - " . localtime((stat($file))[9]) . "\n";

	$lines = `grep "covered at" $file`;
	chomp $lines;
	@covs = split(/\n/, $lines);

	$tCoverage = $covs[0];
	$nCoverage = $covs[1];

	$tCoverage =~ s/.*covered at (.*?)x<\/p>/$1/;
	$nCoverage =~ s/.*covered at (.*?)x<\/p>/$1/;


	$lines = `grep "Total somatic SNV" $file`;
	chomp $lines;
	$numSNVs = $lines;
	$numSNVs =~ s/<h2>Total somatic SNVs called: (.*?)<\/h2>/$1/;


	$lines = `grep "Total somatic indels" $file`;
	chomp $lines;
	$numIndels = $lines;
	$numIndels =~ s/<h2>Total somatic indels called: (.*?)<\/h2>/$1/;


	unless ($donor =~ /MPCC/)
	{
		print HTML "<tr><td><a href=\"$file\">$donor</a></td><td><a href=\"$file\">$sample</a></td><td><a href=\"$file\">$seqType</a></td><td>$nCoverage</td><td>$tCoverage</td><td>$numSNVs</td><td>$numIndels</td><td>$timestamp</td></tr>\n";
	}

}

print HTML "</tbody></table>\n";



print HTML "</body></html>\n";


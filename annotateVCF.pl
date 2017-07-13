#!/usr/bin/perl

use strict;
use warnings;

my $vcf = $ARGV[0];
my $vcfBase = `basename $vcf`;
chomp $vcfBase;

my %refHash = (
	"annovar" => "annovar/2013-06-21",

    "add-annovar-script" => "/u/rdenroche/git/spb-analysis-tools/phoenixPipe/addAnnovarToVCF.pl",
    "add-dbsnp-script" => "/u/rdenroche/git/spb-analysis-tools/phoenixPipe/addDBSnpToVCF.pl",
	"add-cosmic-script" => "/u/rdenroche/git/spb-analysis-tools/phoenixPipe/addCosmicToVCF.pl",
    "add-bed-track-script" => "/u/rdenroche/git/spb-analysis-tools/phoenixPipe/addBedTrackToVCF.pl",
    "remove-mouse-script" => "/u/rdenroche/git/spb-analysis-tools/phoenixPipe/removeMouseVars.pl",
    "select-bed-script" => "/u/rdenroche/git/spb-analysis-tools/phoenixPipe/selectBed.pl",
    "select-filter-script" => "/u/rdenroche/git/spb-analysis-tools/phoenixPipe/selectFilter.pl",
    "annotate-verification-script" => "/u/rdenroche/git/spb-analysis-tools/phoenixPipe/annotateWithVerification.pl",
    "add-dcc-metadata" => "/u/rdenroche/git/spb-analysis-tools/phoenixPipe/addDccMetadata.pl",

    "bwa/0.6.2-mouse-blacklist" => "/isilon4/seqprodbio/phoenix/PCSI/bwamouse/mouseall-bwa_0.6.2.vcf.gz",
    "novocraft/3.00.05-mouse-blacklist" => "/isilon4/seqprodbio/phoenix/PCSI/oldmouse/mouseall-novocraft_2.07.05.vcf.gz",
    "exome-target-tabix" => "/oicr/data/reference/genomes/homo_sapiens_mc/Agilent/SureSelectHumanAllExonV4/S03723314_Regions.merged.bed.gz",
    "dbsnp-tabix" => "/oicr/data/reference/genomes/homo_sapiens_mc/dbSNP/hg19_random/Genomic/dbSNP137/dbSNP137_chr.vcf.gz",
	"cosmic-tabix" => "/oicr/data/reference/genomes/homo_sapiens_mc/Cosmic/v64/CosmicCodingMuts_v64_02042013_noLimit_modified_contig_names.vcf.gz",
    "rmsk" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/tracks/UCSC-rmsk-130723.bed.gz",
    "segdup" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/tracks/UCSC-genomicSuperDups-130723.bed.gz",
    "simpleRepeat" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/tracks/UCSC-simpleRepeat-130723.bed.gz"
);


unless (-e $vcf)
{
	die "$vcf does not exist\n";
}
if (-d "$vcf-annovar")
{
	die "Annovar output already exists.  Please rm -r $vcf-annovar\n";
}

`mkdir $vcf-annovar`;

open (CMD, ">$vcf-annovar/annotate.cmd") or die;

print CMD "module load $refHash{annovar}; module load tabix; convert2annovar.pl --format vcf4 --includeinfo $vcf > $vcf-annovar/$vcfBase.annovar; annotate_variation.pl -geneanno $vcf-annovar/$vcfBase.annovar -buildver hg19 \$HUMANDB\n";

print CMD "cat $vcf | $refHash{'select-filter-script'} PASS | $refHash{'add-annovar-script'} $vcf-annovar/$vcfBase.annovar.variant_function $vcf-annovar/$vcfBase.annovar.exonic_variant_function | $refHash{'add-dbsnp-script'} $refHash{'dbsnp-tabix'} | $refHash{'add-bed-track-script'} $refHash{rmsk} rmsk | $refHash{'add-bed-track-script'} $refHash{segdup} segdup | $refHash{'add-bed-track-script'} $refHash{simpleRepeat} simpleRepeat | $refHash{'add-cosmic-script'} $refHash{'cosmic-tabix'}  > $vcf.annotated.vcf\n";

`qsub -cwd -b y -N anno$vcfBase -e $vcf-annovar/$vcfBase.log -o $vcf-annovar/$vcfBase.log -l h_vmem=4g \"bash $vcf-annovar/annotate.cmd\" > $vcf-annovar/$vcfBase.log`;









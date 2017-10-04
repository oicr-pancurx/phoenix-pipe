#!/usr/bin/perl

use strict;
use warnings;

my $commandLine = join " ", $0, @ARGV;

my $nFastqDir = $ARGV[0];
my $tFastqDir = $ARGV[1];



my $prodCode = "/.mounts/labs/PCSI/production/";
my $pipeCode = "$prodCode/phoenix-pipe";
my $callCode = "$prodCode/simple-caller";
my $parseCode = "$prodCode/phoenix-parse";
my $reportCode = "$prodCode/phoenix-report";
my $jsonCode = "$prodCode/json-report";

my %refHash = (
#	"sge_queue" => "clinical,transient",	# not supported by all subroutines yet
	"sge_queue" => "default,mid_mem,high_mem,marathon,transient,production",	# not supported by all subroutines yet
	"hg19_random" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/fasta/hg19_random.fa",
	"dbSNP" => "/oicr/data/reference/genomes/homo_sapiens_mc/dbSNP/hg19_random/Genomic/dbSNP142/dbsnp142_chr_noINS_noDup.vcf",
	"cosmic" => "/oicr/data/genomes/homo_sapiens_mc/MuTect/hg19_cosmic_OICR_v54_120711.vcf",
	"novocraft/2.07.14" => "/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/novocraft/hg19_random.nix",
	"novocraft/2.08.02" => "/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/novocraft/hg19_random.nix",
	"novocraft/3.00.05" => "/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/novocraft/3.00.05/hg19_random.nix",
	"bwa/0.7.4" => "/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/bwa/0.7.4/hg19_random.fa",
	"bwa/0.6.2" => "/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/bwa/0.6.2/hg19_random.fa",
	"bwa/0.7.12" => "/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/bwa/0.7.12/hg19_random.fa",
	"picard" => "picard/1.90",
	"samtools" => "samtools/0.1.19",
	"strelka/v1.0.7-exome" => "/oicr/local/analysis/sw/strelka/strelka_workflow-1.0.7/etc/strelka_config_bwa_exome_default.ini",
	"strelka/v1.0.7-wgs" => "/oicr/local/analysis/sw/strelka/strelka_workflow-1.0.7/etc/strelka_config_bwa_wgs_default.ini",
	"strelka/v0.4.7-exome" => "/oicr/local/analysis/sw/strelka/strelka_workflow-0.4.7/strelka/etc/strelka_config_bwa_exome.ini",
	"strelka/v0.4.7-wgs" => "/oicr/local/analysis/sw/strelka/strelka_workflow-0.4.7/strelka/etc/strelka_config_bwa_wgs.ini",
	"xenome/1.0.1-r" => "/oicr/data/reference/genomes/xenome/2013Mar/idx",
	"xenome-fix" => "$pipeCode/fixXenome.pl",
	"samStats.pl" => "$jsonCode/samStats.pl",
#	"exome-target" => "/oicr/data/genomes/homo_sapiens/Agilent/SureSelect_Whole_Exome_ICGC_Sanger/GRCh37hg19/sanger.exons.bed.hg19",		# assuming agilent capture!
	"exome-target" => "/oicr/data/reference/genomes/homo_sapiens_mc/Agilent/SureSelectHumanAllExonV4/S03723314_Regions.bed",		# assuming agilent capture!
	"wgs-target" => "/oicr/data/genomes/homo_sapiens/UCSC/Genomic/UCSC_hg19_random/hg19_random.genome.sizes.bed",

	"add-annovar-script" => "$pipeCode/addAnnovarToVCF.pl",
	"add-funseq-script" => "$pipeCode/addFunseqToVCF.pl",
	"add-dbsnp-script" => "$pipeCode/addDBSnpToVCF.pl",
	"add-cosmic-script" => "$pipeCode/addCosmicToVCF.pl",
	"add-bed-track-script" => "$pipeCode/addBedTrackToVCF.pl",
	"remove-mouse-script" => "$pipeCode/removeMouseVars.pl",
	"select-bed-script" => "$pipeCode/selectBed.pl",
	"select-filter-script" => "$pipeCode/selectFilter.pl",
	"annotate-verification-script" => "$pipeCode/annotateWithVerification.pl",
	"add-dcc-metadata" => "$pipeCode/addDccMetadata.pl",
	"intersect-vcf-script" => "$pipeCode/intersectVCF.pl",

	"gatk/1.3.16 helper" => "$pipeCode/oicr_gatk_helper_1_3.sh",
	"gatk/1.3.16 waiter" => "$pipeCode/gatkWaiter.pl",
	"gatk germline filter" => "$pipeCode/pullGermlineFromGATK.pl",

	"bwa/0.6.2-mouse-blacklist" => "/oicr/data/reference/genomes/homo_sapiens_mc/BWA-0.6.2_Mouse_Blacklist/mouseall-bwa_0.6.2.vcf.gz",
	"novocraft/3.00.05-mouse-blacklist" => "/.mounts/labs/prod/phoenix/PCSI/oldmouse/mouseall-novocraft_2.07.05.vcf.gz",
	"exome-target-tabix" => "/oicr/data/reference/genomes/homo_sapiens_mc/Agilent/SureSelectHumanAllExonV4/S03723314_Regions.merged.bed.gz",
	"dbsnp-tabix" => "/oicr/data/reference/genomes/homo_sapiens_mc/dbSNP/hg19_random/Genomic/dbSNP137/dbSNP137_chr.vcf.gz",
	"cosmic-tabix" => "/oicr/data/reference/genomes/homo_sapiens_mc/Cosmic/v64/CosmicCodingMuts_v64_02042013_noLimit_modified_contig_names.vcf.gz",
	"rmsk" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/tracks/UCSC-rmsk-130723.bed.gz",
	"segdup" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/tracks/UCSC-genomicSuperDups-130723.bed.gz",
	"simpleRepeat" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/tracks/UCSC-simpleRepeat-130723.bed.gz",

	"R module" => "R/3.3.0",
	"HMM R script" => "$pipeCode/HMM_pipe_R_v4.sh",

	"blat ref" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/blat/hg19_random.fa.2bit",
	"blat host" => "cn5-80",
	"blat port" => "9998",
	# restart blat server (after qrsh'ing to specific node): nohup gfServer start cn5-83 9998 /oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/blat/hg19_random.fa.2bit -canStop &
	"crest filter script" => "$pipeCode/crestFilter_v2.pl",

	"gavin python" => "/.mounts/labs/steinlab/public/gavin/local/bin/python",		# currently required for delly filtering script
	"delly filter script" => "$pipeCode/dellyFilter.py",
	"sv merge script" => "$pipeCode/svMerge.py",
	"sv annotate script" => "$pipeCode/addGenesToSV.pl",
	"refseq genes tabix" => "/oicr/data/reference/genomes/homo_sapiens_mc/refSeq/refSeq_genes.bed.gz",

	"somatic vcf report" => "$pipeCode/vcfReport.pl",

	"sid bed file" => "$callCode/SID_both.bed",
	"rob caller script" => "$callCode/robCaller_v2.0.pl",
	"rob genotyper script" => "$callCode/robGenotyper.pl",

	"one page script" => "$pipeCode/onePageReport_v2.pl",

	"ATHLATES_novocraft/2.07.06" => "/oicr/local/analysis/sw/ATHLATES/Athlates_2014_04_26/db/ref/hla.clean.fasta.ndx",
	"ATHLATES_bed_path" => "/oicr/local/analysis/sw/ATHLATES/Athlates_2014_04_26/db/bed",
	"ATHLATES_msa_path" => "/oicr/local/analysis/sw/ATHLATES/Athlates_2014_04_26/db/msa",
	"vcf_to_peptide" => "$pipeCode/vcfToPeptide.pl",

	"gc_wig" => "/oicr/data/reference/genomes/homo_sapiens_mc/celluloid/gc.wig",
	"map_wig" => "/oicr/data/reference/genomes/homo_sapiens_mc/celluloid/map.wig",
	"celluloid_mchan" => "$pipeCode/run_celluloid_v0.11_mchan.R",
	"celluloid_gatk_to_af" => "$pipeCode/vcf2het.csh",

	"celluloid_v11" => "$pipeCode/celluloid_v0.11.0_pipeline.r",
	"tabix_bed_script" => "$pipeCode/tabixIndexCelluloid.sh",

	"gatk_to_af" => "$pipeCode/gatkToHetAR.pl",
	"cna_segment_to_fixed" => "$pipeCode/hmmcopySegsForChromo.pl",
	"chromo_script" => "$pipeCode/chromothripsis.R",

	"germline pathogen script" => "$pipeCode/germlinePathogenicFilter.sh",

	"perl lib path" => "export PERL5LIB=/.mounts/labs/PCSI/perl/lib/x86_64-linux-gnu/perl/5.20.2/:/.mounts/labs/PCSI/perl/share/perl/5.20.2:/.mounts/labs/PCSI/perl/lib/perl5/:/.mounts/labs/PCSI/perl/lib/perl5/x86_64-linux-gnu-thread-multi/",

	"qc script" => "$pipeCode/qcPhoenixPipe.pl",
	"phoenix parser" => "$parseCode/phoenixParser.pl",
	"phoenix reporter" => "$reportCode/phoenixReporter.pl",
	"short reporter" => "$reportCode/shortReport.pl",
	"report plots" => "$reportCode/runPlots.sh",

	"cosmicSigNNLS script" => "$pipeCode/vcfToSigs.pl",
);


my $doSub;
my $doClean;

my $sgeUniquifier;

if (exists $ARGV[2])
{
	$doSub = $ARGV[2];
}
else
{
	$doSub = "nope";
}

if (exists $ARGV[3])
{
	$doClean = $ARGV[3];
}
else
{
	$doClean = "nope";
}

if (exists $ARGV[4])
{
	$sgeUniquifier = $ARGV[4];
}
else
{
	$sgeUniquifier = "";
}


my %nMeta;
my %tMeta;

$nMeta{sge_uniq} = $sgeUniquifier;
$tMeta{sge_uniq} = $sgeUniquifier;

getMeta($nFastqDir, \%nMeta);

unless ($tFastqDir eq "single")
{
	getMeta($tFastqDir, \%tMeta);

	if ($doSub eq "go")
	{
		doXenome("xenome/1.0.1-r", \%tMeta, \%refHash);		# only acts if sample type is X
	
	
		doBWAsai("bwa/0.6.2", \%nMeta, \%refHash);
		doBWA("bwa/0.6.2", \%nMeta, \%refHash);
		doBWAsai("bwa/0.6.2", \%tMeta, \%refHash);
		doBWA("bwa/0.6.2", \%tMeta, \%refHash);
		
		doSamStatsLane("bwa/0.6.2", \%nMeta, \%refHash);
		doSamStatsLane("bwa/0.6.2", \%tMeta, \%refHash);
	
		doFilter("bwa/0.6.2", \%nMeta, \%refHash);
		doFilter("bwa/0.6.2", \%tMeta, \%refHash);
	
		doMergeCollapse("bwa/0.6.2", "false", \%nMeta, \%refHash);
		doMergeCollapse("bwa/0.6.2", "false", \%tMeta, \%refHash);
	
		doSamStatsMerge("bwa/0.6.2", \%nMeta, \%refHash);
		doSamStatsMerge("bwa/0.6.2", \%tMeta, \%refHash);
	
	
#		doStrelka("strelka/v0.4.7", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);
		doStrelka("strelka/v1.0.7", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);
		doAnnovar("annovar/2013-06-21", "bwa/0.6.2", "strelka/v1.0.7", \%tMeta, \%refHash);
		doFunSeq("funseq/0.1", "bwa/0.6.2", "strelka/v1.0.7", "annovar/2013-06-21", \%tMeta, \%refHash);
		doAnnotation("bwa/0.6.2", "strelka/v1.0.7", "annovar/2013-06-21", "funseq/0.1", \%nMeta, \%tMeta, \%refHash);

#		doMutect("mutect/1.1.4", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);
		doMutect_split("mutect/1.1.4", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);
		doAnnovar("annovar/2013-06-21", "bwa/0.6.2", "mutect/1.1.4", \%tMeta, \%refHash);
		doFunSeq("funseq/0.1", "bwa/0.6.2", "mutect/1.1.4", "annovar/2013-06-21", \%tMeta, \%refHash);
		doAnnotation("bwa/0.6.2", "mutect/1.1.4", "annovar/2013-06-21", "funseq/0.1", \%nMeta, \%tMeta, \%refHash);

		doIntersect("bwa/0.6.2", "mutect/1.1.4", "strelka/v1.0.7", "bwa/0.6.2/final_strelka-mutect", \%tMeta, \%refHash);
		doMoreAnnovar("bwa/0.6.2", "final_strelka-mutect", \%tMeta, \%refHash);
	
		doGATK("gatk/1.3.16", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);
		doAnnovar("annovar/2013-06-21", "bwa/0.6.2", "gatk/1.3.16", \%tMeta, \%refHash);
		doFunSeq("funseq/0.1", "bwa/0.6.2", "gatk/1.3.16", "annovar/2013-06-21", \%tMeta, \%refHash);
		doAnnotation("bwa/0.6.2", "gatk/1.3.16", "annovar/2013-06-21", "funseq/0.1", \%nMeta, \%tMeta, \%refHash);
		doFilterGATK("bwa/0.6.2", "gatk/1.3.16", \%tMeta, \%refHash);
		doMoreAnnovar("bwa/0.6.2", "final_gatk-germline", \%tMeta, \%refHash);


		doGenotype("bwa/0.6.2", \%nMeta, \%refHash);
		doGenotype("bwa/0.6.2", \%tMeta, \%refHash);

		doGendertype("bwa/0.6.2", \%nMeta, \%refHash);
		doGendertype("bwa/0.6.2", \%tMeta, \%refHash);
	
		doSomatotype("bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);
		doSomatotype("bwa/0.6.2", \%tMeta, \%tMeta, \%refHash);

		#doATHLATES("athlates/1.0", \%nMeta, \%tMeta, \%refHash);
		#doNetMHC("netMHC/3.4", "bwa/0.6.2", "athlates/1.0", \%nMeta, \%tMeta, \%refHash);
		#doNetMHC("netMHC/pan-2.8", "bwa/0.6.2", "athlates/1.0", \%nMeta, \%tMeta, \%refHash);

		doPOLYSOLVER("polysolver/1.0", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);
		doNetMHC("netMHC/pan-2.8", "bwa/0.6.2", "polysolver/1.0", \%nMeta, \%tMeta, \%refHash);

		doCosmicSigNNLS("R/3.3.0", "bwa/0.6.2", \%tMeta, \%refHash);


	#	doRealignGATK("gatk/2.4.9", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);
	#	doRecalGATK("gatk/2.4.9", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);
	
	#	doStrelka("strelka/v0.4.7", "bwa/0.6.2/gatk/2.4.9", \%nMeta, \%tMeta, \%refHash);
	#	doStrelka("strelka/v1.0.7", "bwa/0.6.2/gatk/2.4.9", \%nMeta, \%tMeta, \%refHash);
	
	#	doMutect("mutect/1.1.4", "bwa/0.6.2/gatk/2.4.9", \%nMeta, \%tMeta, \%refHash);
	
	
	
	#	doAnnovar("annovar/2013-06-21", "bwa/0.6.2", "strelka/v0.4.7", \%tMeta, \%refHash);
	#	doAnnotation("bwa/0.6.2", "strelka/v0.4.7", "annovar/2013-06-21", \%nMeta, \%tMeta, \%refHash);
	
	
	
	#	doAnnovar("annovar/2013-06-21", "bwa/0.6.2/gatk/2.4.9", "mutect/1.1.4", \%tMeta, \%refHash);
	#	doAnnovar("annovar/2013-06-21", "bwa/0.6.2/gatk/2.4.9", "strelka/v1.0.7", \%tMeta, \%refHash);
	#	doAnnovar("annovar/2013-06-21", "bwa/0.6.2/gatk/2.4.9", "strelka/v0.4.7", \%tMeta, \%refHash);
	
	#	doAnnotation("bwa/0.6.2/gatk/2.4.9", "mutect/1.1.4", "annovar/2013-06-21", \%nMeta, \%tMeta, \%refHash);
	#	doAnnotation("bwa/0.6.2/gatk/2.4.9", "strelka/v1.0.7", "annovar/2013-06-21", \%nMeta, \%tMeta, \%refHash);
	#	doAnnotation("bwa/0.6.2/gatk/2.4.9", "strelka/v0.4.7", "annovar/2013-06-21", \%nMeta, \%tMeta, \%refHash);
	
	
		if ($tMeta{"sample_type"} eq "wgs")
		{
			doHMMcopy("HMMcopy/0.1.1", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);

#			doCREST("crest/1.0.1", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);
#			doFilterCREST("crest/1.0.1", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);

			doCREST("crest/alpha", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);
			doFilterCREST("crest/alpha", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);

			doDelly("delly/0.5.5", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);
			doFilterDelly("delly/0.5.5", "bwa/0.6.2", \%tMeta, \%refHash);

			doMergeSV("crest/alpha", "delly/0.5.5", "bwa/0.6.2", \%tMeta, \%refHash);

			doCelluloidMC("R/3.3.0","bwa/0.6.2",\%tMeta, \%refHash);
			doChromothripsis("R/3.3.0","bwa/0.6.2",\%tMeta, \%refHash);

			doCelluloid("R/3.3.0","bwa/0.6.2",\%tMeta, \%refHash);

			doIntegration("bwa/0.6.2", \%tMeta, \%refHash);
		}
	
	
	
	#	doNovoalign("novocraft/3.00.05", \%nMeta, \%refHash);
	#	doNovoalign("novocraft/3.00.05", \%tMeta, \%refHash);
	
	#	doSamStatsLane("novocraft/3.00.05", \%nMeta, \%refHash);
	#	doSamStatsLane("novocraft/3.00.05", \%tMeta, \%refHash);
	
	#	doFilter("novocraft/3.00.05", \%nMeta, \%refHash);
	#	doFilter("novocraft/3.00.05", \%tMeta, \%refHash);
		
	#	doMergeCollapse("novocraft/3.00.05", \%nMeta, \%refHash);
	#	doMergeCollapse("novocraft/3.00.05", \%tMeta, \%refHash);
	
	#	doSamStatsMerge("novocraft/3.00.05", \%nMeta, \%refHash);
	#	doSamStatsMerge("novocraft/3.00.05", \%tMeta, \%refHash);
	
	
	#	doStrelka("strelka/v0.4.7", "novocraft/3.00.05", \%nMeta, \%tMeta, \%refHash);
	#	doStrelka("strelka/v1.0.7", "novocraft/3.00.05", \%nMeta, \%tMeta, \%refHash);
		
	#	doMutect("mutect/1.1.4", "novocraft/3.00.05", \%nMeta, \%tMeta, \%refHash);
	
	#	doRealignGATK("gatk/2.4.9", "novocraft/3.00.05", \%nMeta, \%tMeta, \%refHash);
	#	doRecalGATK("gatk/2.4.9", "novocraft/3.00.05", \%nMeta, \%tMeta, \%refHash);
	
	#	doStrelka("strelka/v0.4.7", "novocraft/3.00.05/gatk/2.4.9", \%nMeta, \%tMeta, \%refHash);
	#	doStrelka("strelka/v1.0.7", "novocraft/3.00.05/gatk/2.4.9", \%nMeta, \%tMeta, \%refHash);
	
	#	doMutect("mutect/1.1.4", "novocraft/3.00.05/gatk/2.4.9", \%nMeta, \%tMeta, \%refHash);
	
	
	#	doAnnovar("annovar/2013-06-21", "novocraft/3.00.05", "mutect/1.1.4", \%tMeta, \%refHash);
	#	doAnnovar("annovar/2013-06-21", "novocraft/3.00.05", "strelka/v1.0.7", \%tMeta, \%refHash);
	#	doAnnovar("annovar/2013-06-21", "novocraft/3.00.05", "strelka/v0.4.7", \%tMeta, \%refHash);
	
	#	doAnnovar("annovar/2013-06-21", "novocraft/3.00.05/gatk/2.4.9", "mutect/1.1.4", \%tMeta, \%refHash);
	#	doAnnovar("annovar/2013-06-21", "novocraft/3.00.05/gatk/2.4.9", "strelka/v1.0.7", \%tMeta, \%refHash);
	#	doAnnovar("annovar/2013-06-21", "novocraft/3.00.05/gatk/2.4.9", "strelka/v0.4.7", \%tMeta, \%refHash);
	
	
	#	doAnnotation("novocraft/3.00.05", "mutect/1.1.4", "annovar/2013-06-21", \%nMeta, \%tMeta, \%refHash);
	#	doAnnotation("novocraft/3.00.05", "strelka/v1.0.7", "annovar/2013-06-21", \%nMeta, \%tMeta, \%refHash);
	#	doAnnotation("novocraft/3.00.05", "strelka/v0.4.7", "annovar/2013-06-21", \%nMeta, \%tMeta, \%refHash);
	
	#	doAnnotation("novocraft/3.00.05/gatk/2.4.9", "mutect/1.1.4", "annovar/2013-06-21", \%nMeta, \%tMeta, \%refHash);
	#	doAnnotation("novocraft/3.00.05/gatk/2.4.9", "strelka/v1.0.7", "annovar/2013-06-21", \%nMeta, \%tMeta, \%refHash);
	#	doAnnotation("novocraft/3.00.05/gatk/2.4.9", "strelka/v0.4.7", "annovar/2013-06-21", \%nMeta, \%tMeta, \%refHash);
	
	
	
		doQC("bwa/0.6.2", \%tMeta, \%nMeta, \%refHash, $commandLine);
		doPhoenixParser("bwa/0.6.2", \%tMeta, \%refHash);
	
		doProvenance(\%nMeta, \%tMeta);
	
	#	doVarScan2("varscan/2.3.2", "novocraft/3.00.05", \%nMeta, \%tMeta);
	#	doVarScan2("varscan/2.3.2", "bwa/0.6.2", \%nMeta, \%tMeta);
	
	
	#	doGATK("gatk/1.3.16", "novocraft/3.00.05", \%nMeta, \%tMeta);
	#	doGATK("gatk/1.3.16", "bwa/0.6.2", \%nMeta, \%tMeta);
	
	#	doTimSomatic("novocraft/3.00.05/gatk/1.3.16", \%nMeta, \%tMeta);
	#	doTimSomatic("bwa/0.6.2/gatk/1.3.16", \%nMeta, \%tMeta);
	
	#	doSplitGATK("novocraft/3.00.03/gatk/1.3.16", \%nMeta, \%tMeta);
	#	doSplitGATK("bwa/0.7.0/gatk/1.3.16", \%nMeta, \%tMeta);
	
	#	doStrelka("strelka/v0.4.7", "novocraft/3.00.03/gatk/1.3.16", \%nMeta, \%tMeta, $seqType);
	#	doStrelka("strelka/v0.4.7", "bwa/0.7.0/gatk/1.3.16", \%nMeta, \%tMeta, $seqType);
	
	#	doMutect("mutect/1.1.4", "novocraft/3.00.03/gatk/1.3.16", \%nMeta, \%tMeta, $seqType);
	#	doMutect("mutect/1.1.4", "bwa/0.7.0/gatk/1.3.16", \%nMeta, \%tMeta, $seqType);
	
	#	doVarScan2("varscan/2.3.2", "novocraft/3.00.03/gatk/1.3.16", \%nMeta, \%tMeta);
	#	doVarScan2("varscan/2.3.2", "bwa/0.7.0/gatk/1.3.16", \%nMeta, \%tMeta);
	
	
	#	doNewGATK("gatk/2.4.9", "novocraft/3.00.05", \%nMeta, \%tMeta, \%refHash);
	#	doNewGATK("gatk/2.4.9", "bwa/0.6.2", \%nMeta, \%tMeta, \%refHash);
	
	
	#	doVarScan2("varscan/2.3.2", "novocraft/3.00.03/gatk/2.4.9", \%nMeta, \%tMeta);
	#	doVarScan2("varscan/2.3.2", "bwa/0.7.0/gatk/2.4.9", \%nMeta, \%tMeta);
	
	#	doClean(\%nMeta, $doClean);
	#	doClean(\%tMeta, $doClean);
	}
	elsif ($doSub eq "mem")
	{
		doXenome("xenome/1.0.1-r", \%nMeta, \%refHash);
	
		doBWAmem("bwa/0.7.12", \%nMeta, \%refHash);
		doBWAmem("bwa/0.7.12", \%tMeta, \%refHash);
	
		doSamStatsLane("bwa/0.7.12", \%nMeta, \%refHash);
		doFilter("bwa/0.7.12", \%nMeta, \%refHash);
		doMergeCollapse("bwa/0.7.12", "false", \%nMeta, \%refHash);
		doSamStatsMerge("bwa/0.7.12", \%nMeta, \%refHash);
	
		doSamStatsLane("bwa/0.7.12", \%tMeta, \%refHash);
		doFilter("bwa/0.7.12", \%tMeta, \%refHash);
		doMergeCollapse("bwa/0.7.12", "false", \%tMeta, \%refHash);
		doSamStatsMerge("bwa/0.7.12", \%tMeta, \%refHash);

		doStrelka("strelka/v1.0.7", "bwa/0.7.12", \%nMeta, \%tMeta, \%refHash);
		doAnnovar("annovar/2013-06-21", "bwa/0.7.12", "strelka/v1.0.7", \%tMeta, \%refHash);
		doFunSeq("funseq/0.1", "bwa/0.7.12", "strelka/v1.0.7", "annovar/2013-06-21", \%tMeta, \%refHash);
		doAnnotation("bwa/0.7.12", "strelka/v1.0.7", "annovar/2013-06-21", "funseq/0.1", \%nMeta, \%tMeta, \%refHash);
	
		doMutect_split("mutect/1.1.4", "bwa/0.7.12", \%nMeta, \%tMeta, \%refHash);
		doAnnovar("annovar/2013-06-21", "bwa/0.7.12", "mutect/1.1.4", \%tMeta, \%refHash);
		doFunSeq("funseq/0.1", "bwa/0.7.12", "mutect/1.1.4", "annovar/2013-06-21", \%tMeta, \%refHash);
		doAnnotation("bwa/0.7.12", "mutect/1.1.4", "annovar/2013-06-21", "funseq/0.1", \%nMeta, \%tMeta, \%refHash);
	
		doIntersect("bwa/0.7.12", "mutect/1.1.4", "strelka/v1.0.7", "bwa/0.7.12/final_strelka-mutect", \%tMeta, \%refHash);
	
		doHMMcopy("HMMcopy/0.1.1", "bwa/0.7.12", \%nMeta, \%tMeta, \%refHash);
	
		doCREST("crest/alpha", "bwa/0.7.12", \%nMeta, \%tMeta, \%refHash);
		doFilterCREST("crest/alpha", "bwa/0.7.12", \%nMeta, \%tMeta, \%refHash);
	
		doDelly("delly/0.5.5", "bwa/0.7.12", \%nMeta, \%tMeta, \%refHash);
		doFilterDelly("delly/0.5.5", "bwa/0.7.12", \%tMeta, \%refHash);
	
		doMergeSV("crest/alpha", "delly/0.5.5", "bwa/0.7.12", \%tMeta, \%refHash);
	
		doGATK("gatk/1.3.16", "bwa/0.7.12", \%nMeta, \%tMeta, \%refHash);
		doAnnovar("annovar/2013-06-21", "bwa/0.7.12", "gatk/1.3.16", \%tMeta, \%refHash);
		doFunSeq("funseq/0.1", "bwa/0.7.12", "gatk/1.3.16", "annovar/2013-06-21", \%tMeta, \%refHash);
		doAnnotation("bwa/0.7.12", "gatk/1.3.16", "annovar/2013-06-21", "funseq/0.1", \%nMeta, \%tMeta, \%refHash);
		doFilterGATK("bwa/0.7.12", "gatk/1.3.16", \%tMeta, \%refHash);
	
		doIntegration("bwa/0.7.12", \%tMeta, \%refHash);
	
		doGenotype("bwa/0.7.12", \%nMeta, \%refHash);
		doGenotype("bwa/0.7.12", \%tMeta, \%refHash);
	
		doGendertype("bwa/0.7.12", \%nMeta, \%refHash);
		doGendertype("bwa/0.7.12", \%tMeta, \%refHash);
	
		doSomatotype("bwa/0.7.12", \%nMeta, \%tMeta, \%refHash);
	
		doGenotype("bwa/0.7.12", \%nMeta, \%refHash);
		doGendertype("bwa/0.7.12", \%tMeta, \%refHash);
	}
}
elsif ($doSub eq "go")
{
	doXenome("xenome/1.0.1-r", \%nMeta, \%refHash);		# only acts if sample type is X

	doBWAsai("bwa/0.6.2", \%nMeta, \%refHash);
	doBWA("bwa/0.6.2", \%nMeta, \%refHash);
	
	doSamStatsLane("bwa/0.6.2", \%nMeta, \%refHash);

	doFilter("bwa/0.6.2", \%nMeta, \%refHash);

	doMergeCollapse("bwa/0.6.2", "false", \%nMeta, \%refHash);

	doSamStatsMerge("bwa/0.6.2", \%nMeta, \%refHash);

	doQCsingle("bwa/0.6.2", \%nMeta, \%refHash, $commandLine);

	#doProvenance(\%nMeta, \%tMeta); needs to be refactored to handle single samples
}
elsif ($doSub eq "gatk")
{


	doXenome("xenome/1.0.1-r", \%nMeta, \%refHash);		# only acts if sample type is X

	doBWAsai("bwa/0.6.2", \%nMeta, \%refHash);
	doBWA("bwa/0.6.2", \%nMeta, \%refHash);
	
	doSamStatsLane("bwa/0.6.2", \%nMeta, \%refHash);

	doFilter("bwa/0.6.2", \%nMeta, \%refHash);

	doMergeCollapse("bwa/0.6.2", "false",\%nMeta, \%refHash);

	doSamStatsMerge("bwa/0.6.2", \%nMeta, \%refHash);

	doSingleGATK("gatk/1.3.16", "bwa/0.6.2", \%nMeta, \%refHash);
	doAnnovar("annovar/2013-06-21", "bwa/0.6.2", "gatk/1.3.16", \%nMeta, \%refHash);
	doFunSeq("funseq/0.1", "bwa/0.6.2", "gatk/1.3.16", "annovar/2013-06-21", \%nMeta, \%refHash);
	doSingleAnnotation("bwa/0.6.2", "gatk/1.3.16", "annovar/2013-06-21", "funseq/0.1", \%nMeta, \%refHash);
	doFilterGATK("bwa/0.6.2", "gatk/1.3.16", \%nMeta, \%refHash);

}
elsif ($doSub eq "mega")
{
	doXenome("xenome/1.0.1-r", \%nMeta, \%refHash);

	doBWAmem("bwa/0.7.12", \%nMeta, \%refHash);

	doSamStatsLane("bwa/0.7.12", \%nMeta, \%refHash);
	doFilter("bwa/0.7.12", \%nMeta, \%refHash);
	doMergeCollapse("bwa/0.7.12", "false", \%nMeta, \%refHash);
	doSamStatsMerge("bwa/0.7.12", \%nMeta, \%refHash);

	doGenotype("bwa/0.7.12", \%nMeta, \%refHash);
	doGendertype("bwa/0.7.12", \%nMeta, \%refHash);

	doRealignGATK("gatk/3.5.0", "bwa/0.7.12", "wgs", \%nMeta, "single", \%refHash);
	doRecalGATK("gatk/3.5.0", "bwa/0.7.12", "wgs", \%nMeta, "single", \%refHash);
	doHaplotypeGATK("gatk/3.5.0", "bwa/0.7.12", "wgs", \%nMeta, \%refHash);
}


sub getMeta
{
	my $dir = $_[0];
	my $metaHash = $_[1];

	my $cwd = `pwd`;
	chomp $cwd;

	$metaHash->{"working_dir"} = $cwd;

	my $absPath;

	my $fileList = `ls $dir`;
	chomp $fileList;
	my @files = split(/\n/, $fileList);

	my ($library, $sample, $sampleGroup, $sampleType, $tissueType, $run, $lane, $barcode, $instrument, $flowcell, $swidIUS, $qualFormat);
	my ($file1, $file2);


	$metaHash->{username} = `whoami`;
	chomp $metaHash->{username};
	
	$metaHash->{time} = `date`;
	chomp $metaHash->{time};

	$metaHash->{date} = `date +%Y%m%d-%H%M%S`;
	chomp $metaHash->{date};


	for $file1 (@files)
	{
		$library = "nullynull";
		# SWID_9555452_PCSI_0713_Ab_M_PE_463_WG_526_170517_D00343_0169_BCAYV8ANXX_AGTTCC_L006_R1_001.fastq.gz
		if ($file1 =~ /^SWID_(.*?)_(.*?_.*?_.*?_.*?_.*?_.*?_.*?)_.*?_(.*?_.*?_.*?_.*?)_(.*?)_L00(.)_R1_001\.fastq\.gz$/)
		{
			$swidIUS = $1;
			$library = $2;
			$run = $3;
			$barcode = $4;
			$lane = $5;
			$qualFormat = "sanger";		# probably
		}
		# SWID_298986_CME_0003_St_R_PE_406_EX_121023_SN1068_0102_AC1CEHACXX_NoIndex_L005_R1_001.fastq.gz
		elsif ($file1 =~ /^SWID_(.*?)_(.*?_.*?_.*?_.*?_.*?_.*?_.*?)_(.*?_.*?_.*?_.*?)_(.*?)_L00(.)_R1_001\.fastq\.gz$/)
		{
			$swidIUS = $1;
			$library = $2;
			$run = $3;
			$barcode = $4;
			$lane = $5;
			$qualFormat = "sanger";		# probably
		}
		# SWID_7633_PCSI_0039_Pa_X_PE_356_EX_110805_SN803_0063_AB01E5ACXX_CGATGT_L002_R1_001.fastq
		elsif ($file1 =~ /^SWID_(.*?)_(.*?_.*?_.*?_.*?_.*?_.*?_.*?)_(.*?_.*?_.*?_.*?)_(.*?)_L00(.)_R1_001\.fastq$/)
		{
			$swidIUS = $1;
			$library = $2;
			$run = $3;
			$barcode = $4;
			$lane = $5;
			$qualFormat = "sanger";		# probably
		}
		# $sampleGroup/$sample/$seqType/fastq/SWID_${swidIUS}_${library}_${run}_${barcode}_L00${lane}_R${readNum}_sequence.txt.gz
		elsif ($file1 =~ /^SWID_(.*?)_(.*?_.*?_.*?_.*?_.*?_.*?_.*?)_(.*?_.*?_.*?_.*)_(.*?)_L00(.)_R1_sequence\.txt\.gz$/)
		{
			$swidIUS = $1;
			$library = $2;
			$run = $3;
			$barcode = $4;
			$lane = $5;
			$qualFormat = "illumina";		# probably :/
		}
		# PCSI_0001_Ly_R_PE_400_EX_NoIndex_6_110401_SN393_0121_A81CGEABXX.bam_R1.fastq.gz
		elsif ($file1 =~ /^(PCSI_...._.._._.._..._..)_(.*?)_(.)_(.*).bam_R1.fastq.gz/)
		{
			$swidIUS = "PCSIfromBam";
			$library = $1;
			$barcode = $2;
			$lane = $3;
			$run = $4;
			$qualFormat = "sanger";
		}

		# POPS_150306_SN1068_0179_AC6CC0ACXX_RD-43277_ACGTATCA_L006_R1.fastq.gz
		elsif ($file1 =~ /POPS_(.*?)_(RD-.*?)_(.*?)_L00(.)_R1.fastq.gz/)
		{
			$swidIUS = "POPS";
			$library = "$2-$3";
			$run = $1;
			$barcode = $3;
			$lane = $4;
			$qualFormat = "sanger";
		}
		# HI.0655.004.Index_1.99-001_pool_R1.fastq.gz
		elsif ($file1 =~ /^(HI\.....)\.00(.)\.(Index_.*?)\.(.*?)_R1\.fastq\.gz$/)
		{
			$swidIUS = "QCPS";
			$run = $1;
			$lane = $2;
			$barcode = $3;
			$library = $4;
			$qualFormat = "sanger";
		}
		# AC1GPPACXX_8_Library_BIA4882A1AR21_R1.fastq.gz - aussie format
		elsif ($file1 =~ /^(.*?)_(.*?)_(.*)_R1\.fastq\.gz$/)
		{
			$swidIUS = "QMCG";
			$run = $1;
			$lane = $2;
			$library = $3;
			$barcode = "NoIndex";
			$qualFormat = "sanger";		# probably...
		}
		# 130215_0424_H06Y9ADXX_1_IL-TP-016_1.sanfastq.gz
		elsif ($file1 =~ /^(.*?_.*?_.*?)_(.)_(.*?)_1\.sanfastq\.gz/)
		{
			$swidIUS = "SCOTLAND";
			$run = $1;
			$lane = $2;
			$library = $3;
			$barcode = "NoIndex";
			$qualFormat = "sanger";
		}
		# SSR#####_1.fastq.gz
		elsif ($file1 =~ /^(SRR.*?)_1\.fastq\.gz$/)
		{
			$swidIUS = "SRA";
			$run = "SRA";
			$lane = "SRA";
			$library = $1;
			$barcode = "NoIndex";
			$qualFormat = "sanger";		# I hope
		}
		# FDCG_0001_1_UN_Whole_C1_KHSCR_L12002_C7BH9ACXX_ATGCCTAA_L001_R1_001.fastq.gz
		elsif ($file1 =~ /^(FDCG_.*)_(.*?)_(.*?)_L00(.)_R1_001\.fastq\.gz/)
		{
			$swidIUS = "FDGC";
			$library = $1;
			$run = $2;
			$barcode = $3;
			$lane = $4;
			$qualFormat = "sanger";		# hopefully!
		}
        elsif ($file1 =~ /^(LP.*?-DNA.*?)\.bam_R1\.fastq\.gz/)
        {
            $swidIUS = "JHMI";
            $library = $1;
            $run = "JHMI";
            $barcode = "NoIndex";
            $lane = "JHMI";
            $qualFormat = "sanger";         # fingers crossed
        }

		unless ($library eq "nullynull")
		{
			# FDCG for Ashton
			if ($library =~ /^FDGC.*/)
			{
				$instrument = "FDGC";
				$flowcell = $run;
				$tissueType = "R";
			}
			elsif ($library =~ /^(.*?_.*?_.*?_.*?)_.*_(.*?)$/)
			{
				$sample = $1;
				$sampleType = $2;

				$sampleGroup = $sample;
				$sampleGroup =~ s/_.._.//;
				$sampleGroup =~ s/_//;


				if ($sampleType eq "EX")
				{
					$sampleType = "exome";
				}
				elsif ($sampleType eq "WG")
				{
					$sampleType = "wgs";
				}
				else
				{
					$sampleType = "other";
				}
			}
			elsif ($library =~ /^Library.*/)
			{
				# Aussie bam - sample names will come from path
				$instrument = "QCMG";
				$flowcell = "QCMG";
			}
			elsif ($library =~ /^IL-TP.*/)
			{
				# Scotland fastq - sample names will come from path
				$instrument = "SCOTLAND";
				$flowcell = "SCOTLAND";
				$tissueType = "R";
			}
			elsif ($library =~ /^RD-.*/)
			{
				# POPS fastq - sample names will come from path
				$tissueType = "R";
			}
			elsif ($swidIUS eq "QCPS")
			{
				# QPCS fastq - sample names from path
				$instrument = "QPCS";
				$flowcell = "QPCS";
				$tissueType = "R";
			}
			elsif ($library =~ /SRR/)
			{
				#SRA file
				$instrument = "SRA";
				$flowcell = "SRA";
				$tissueType = "R";
			}
            elsif ($library =~ /^LP/)
            {
                $instrument = "JHMI";
                $flowcell = "JHMI";
                $tissueType = "R";
            }
			else
			{
				die "Couldn't parse library $library\n";
			}

			if ($run =~ /^.*?_(.*?)_.*?_(.*?)$/)
			{
				$instrument = $1;
				$flowcell = $2;
			}

			$file2 = $file1;
			$file2 =~ s/R1_001\.fastq\.gz/R2_001.fastq.gz/;
			$file2 =~ s/R1_001\.fastq/R2_001.fastq/;
			$file2 =~ s/R1_sequence\.txt\.gz/R2_sequence.txt.gz/;
			$file2 =~ s/R1\.fastq\.gz/R2.fastq.gz/;
			$file2 =~ s/1\.sanfastq\.gz/2.sanfastq.gz/;
			$file2 =~ s/_1\.fastq\.gz/_2.fastq.gz/;

			unless ((-e "$dir/$file2") or (-l "$dir/$file2"))
			{
				$file2 = "Not found";
			}

			if ($dir =~ /unaligned_bam/)
			{
				$file2 = "Not found";
			}


			# respect dir donor, sample_groupID and type if it fits the format
			if ($dir =~ /(.*)\/(.*)\/(.*)\/fastq/)
			{
				$sampleGroup = $1;
				$sample = $2;
				$sampleType = $3;
				$sampleGroup =~ s/^.*\///g;  #pmk 20141110
				
			}
			if ($dir =~ /(.*)\/(.*)\/(.*)\/unaligned_bam/)
			{
				$sampleGroup = $1;
				$sample = $2;
				$sampleType = $3;
				$sampleGroup =~ s/^.*\///g;  # pmk 20141110
			}

			unless (($instrument eq "SCOTLAND") or ($library =~ /^RD-.*/) or ($instrument eq "QPCS") or ($instrument eq "SRA") or ($instrument eq "JHMI"))
			{
				$tissueType = (split(/_/, $sample))[3];
			}


			print "Read 1:       $dir/$file1\n";
			print "Read 2:       $dir/$file2\n";
			print "Sample Group: $sampleGroup\n";
			print "Sample:       $sample\n";
			print "Tissue Type:  $tissueType\n";
			print "Sample Type:  $sampleType\n";
			print "Library:      $library\n";
			print "Barcode:      $barcode\n";
			print "Lane:         $lane\n";
			print "Run:          $run\n";
			print "Instrument:   $instrument\n";
			print "Flowcell:     $flowcell\n";
			print "\n";

			$metaHash->{"fastq"}{$file1}{"read1"} = "$cwd/$dir/$file1";
			$metaHash->{"fastq"}{$file1}{"read2"} = "$cwd/$dir/$file2";

			$metaHash->{"fastq"}{$file1}{"library"} = $library;
			$metaHash->{"fastq"}{$file1}{"barcode"} = $barcode;
			$metaHash->{"fastq"}{$file1}{"lane"} = $lane;
			$metaHash->{"fastq"}{$file1}{"run"} = $run;
			$metaHash->{"fastq"}{$file1}{"instrument"} = $instrument;
			$metaHash->{"fastq"}{$file1}{"flowcell"} = $flowcell;
			$metaHash->{"fastq"}{$file1}{"swid_ius"} = $swidIUS;
			$metaHash->{"fastq"}{$file1}{"quality_format"} = $qualFormat;


			if (exists $metaHash->{"sample"})
			{
				unless (($sample eq $metaHash->{"sample"}) and ($sampleType eq $metaHash->{"sample_type"}) and ($sampleGroup eq $metaHash->{"sample_group"}) and ($tissueType eq $metaHash->{"tissue_type"}))
				{
					die "Multiple samples ($sampleGroup, $sample, $sampleType, $tissueType ne $metaHash->{sample_group}, $metaHash->{sample}, $metaHash->{sample_type}, $metaHash->{tissue_type})!!\n";
				}
			}
			else
			{
				$metaHash->{"sample"} = $sample;
				$metaHash->{"sample_type"} = $sampleType;
				$metaHash->{"sample_group"} = $sampleGroup;
				$metaHash->{"tissue_type"} = $tissueType;
			}

			$absPath = `readlink -f $cwd/$dir/$file1`;
			chomp $absPath;
			$metaHash->{provenence}{"$cwd/$dir/$file1"}{fastq} = "ln -s $absPath $cwd/$dir/$file1\n";

			unless ($file2 eq "Not found")
			{
				$absPath = `readlink -f $cwd/$dir/$file2`;
				chomp $absPath;
				$metaHash->{provenence}{"$cwd/$dir/$file2"}{fastq} = "ln -s $absPath $cwd/$dir/$file2\n";
			}
		}
	}
}


sub doXenome
{
	my $module = $_[0];
	my $metaHash = $_[1];
	my $refHash = $_[2];

	my $pathRoot = buildPath($metaHash);
	my $dir = "$pathRoot/$module/";

	my $sgePre = "x$metaHash->{sge_uniq}";
	my $outName;

	my $files;
	my $command;

	my @pieces = qw(human ambiguous both neither);
	my $catCommand;

	if (($metaHash->{tissue_type} eq "X") or ($metaHash->{tissue_type} eq "O"))
	{
		`mkdir -p $dir`;

		for my $lane (sort keys %{ $metaHash->{fastq} })
		{
	
			$outName = "$metaHash->{fastq}{$lane}{library}_$metaHash->{fastq}{$lane}{barcode}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{run}";
			if ($metaHash->{fastq}{$lane}{read2} =~ /.*Not found/)
			{
				$files = "-i $metaHash->{fastq}{$lane}{read1}";
	
				$catCommand = "cat";
				for my $piece (@pieces)
				{
					$catCommand .= " $dir/${outName}_${piece}.fastq";
				}
				$catCommand .= " | $refHash->{'xenome-fix'} > $dir/$outName.read1.fastq;";
			}
			else
			{
				$files = "--pairs -i $metaHash->{fastq}{$lane}{read1} -i $metaHash->{fastq}{$lane}{read2}";
	
				$catCommand = "";
				for my $read (qw(1 2))
				{
					$catCommand .= " cat";
					for my $piece (@pieces)
					{
						$catCommand .= " $dir/${outName}_${piece}_$read.fastq";
					}
					$catCommand .= " | $refHash->{'xenome-fix'} > $dir/$outName.read$read.fastq;";
				}
			}
			
			$command = "module load $module; mkdir $dir/$outName.tmp; xenome classify -M 16 --tmp-dir $dir/$outName.tmp -P $refHash->{$module} --graft-name human --host-name mouse --output-filename-prefix $dir/$outName $files; $catCommand\n";
	
			unless (-e "$dir/$outName.touch")
			{
				`touch $dir/$outName.touch`;		# so that future runs won't cause concurrent processing
	
				open (SUBFILE, ">$dir/$outName.sub") or die;
				print SUBFILE "qsub -cwd -b y -q $refHash->{sge_queue} -N $sgePre$outName -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=24g -q $refHash->{sge_queue} \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
				close SUBFILE;
	
				open (COMMAND, ">$dir/$outName.cmd") or die;
				print COMMAND $command;
				print COMMAND "\necho phoenixPipe/$module-done\n";
				close COMMAND;
	
				`bash $dir/$outName.sub`;
			}
			else
			{
				warn " [XENOME] Skipping: $dir/$outName.touch\n";
			}
	
			$metaHash->{fastq}{$lane}{xenome}{read1} = "$dir/$outName.read1.fastq";
			push (@{ $metaHash->{provenence}{"$dir/$outName.read1.fastq"}{input}}, $metaHash->{fastq}{$lane}{read1});
			$metaHash->{provenence}{"$dir/$outName.read1.fastq"}{command} = "$dir/$outName.cmd";

			if ($metaHash->{fastq}{$lane}{read2} =~ /.*Not found/)
			{
				$metaHash->{fastq}{$lane}{xenome}{read2} = "Not found";
			}
			else
			{
				$metaHash->{fastq}{$lane}{xenome}{read2} = "$dir/$outName.read2.fastq";
				push (@{ $metaHash->{provenence}{"$dir/$outName.read2.fastq"}{input}}, $metaHash->{fastq}{$lane}{read2});
				$metaHash->{provenence}{"$dir/$outName.read2.fastq"}{command} = "$dir/$outName.cmd";
			}
	
			$metaHash->{fastq}{$lane}{xenome}{hold_jid} = "$sgePre$outName";
	
			push (@{ $metaHash->{clean}{full} }, "$dir/$outName.read1.fastq");
			push (@{ $metaHash->{clean}{full} }, "$dir/$outName.read2.fastq");
	
			push (@{ $metaHash->{clean}{temp}{files} }, "$dir/$outName.tmp");
			push (@{ $metaHash->{clean}{temp}{hold_jids} }, "$sgePre$outName");
	
			push (@{ $metaHash->{sub_list} }, "$dir/$outName.sub");
			push (@{ $metaHash->{cmd_list} }, "$dir/$outName.cmd");



		}
	}
	else
	{
		warn " [XENOME] Skipping $metaHash->{sample} (not a xenograft)\n";
	}


}


sub doNovoalign
{
	my $module = $_[0];
	my $metaHash = $_[1];
	my $refHash = $_[2];

	my $pathRoot = buildPath($metaHash);
	my $dir = "$pathRoot/$module/lanes/";
	`mkdir -p $dir`;

	my $sgePre = "n$metaHash->{sge_uniq}";
	my $outName;

	my $files;
	my $command;
	my $hold_jid;

	my $fastqFormat = "";

	for my $lane (sort keys %{ $metaHash->{fastq} })
	{

		$outName = "$metaHash->{fastq}{$lane}{library}_$metaHash->{fastq}{$lane}{barcode}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{run}";

		if ($metaHash->{fastq}{$lane}{read2} =~ /.*Not found/)
		{
			if ($metaHash->{tissue_type} eq "X")
			{
				$files = "$metaHash->{fastq}{$lane}{xenome}{read1}";
				$hold_jid = "-hold_jid $metaHash->{fastq}{$lane}{xenome}{hold_jid}";

				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $metaHash->{fastq}{$lane}{xenome}{read1});
			}
			else
			{
				$files = "$metaHash->{fastq}{$lane}{read1}";
				$hold_jid = "";

				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $metaHash->{fastq}{$lane}{read1});
			}
		}
		else
		{
			if ($metaHash->{tissue_type} eq "X")
			{
				$files = "$metaHash->{fastq}{$lane}{xenome}{read1} $metaHash->{fastq}{$lane}{xenome}{read2}";
				$hold_jid = "-hold_jid $metaHash->{fastq}{$lane}{xenome}{hold_jid}";

				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $metaHash->{fastq}{$lane}{xenome}{read1});
				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $metaHash->{fastq}{$lane}{xenome}{read2});
			}
			else
			{
				$files = "$metaHash->{fastq}{$lane}{read1} $metaHash->{fastq}{$lane}{read2}";
				$hold_jid = "";

				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $metaHash->{fastq}{$lane}{read1});
				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $metaHash->{fastq}{$lane}{read2});
			}
		}

		if ($metaHash->{fastq}{$lane}{"quality_format"} eq "illumina")
		{
			$fastqFormat = "-F ILMFQ";
		}
		elsif ($metaHash->{fastq}{$lane}{"quality_format"} eq "sanger")
		{
			$fastqFormat = "-F ILM1.8";
		}

		unless (-e "$dir/$outName.touch")
		{
			`touch $dir/$outName.touch`;		# so that future runs won't cause concurrent processing

			open (SUBFILE, ">$dir/$outName.sub") or die;
			print SUBFILE "qsub -cwd -b y -q $refHash->{sge_queue} -N $sgePre$outName $hold_jid -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=24g,exclusive=1 -q $refHash->{sge_queue} \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
			close SUBFILE;

			open (COMMAND, ">$dir/$outName.cmd") or die;
			print COMMAND "module load $module; module load $refHash->{picard}; novoalign $fastqFormat -d $refHash->{$module} -f $files -r ALL 1 -R 0 -oSAM \$\'\@RG\\tID:$metaHash->{fastq}{$lane}{run}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{barcode}\\tPU:$metaHash->{fastq}{$lane}{run}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{barcode}\\tLB:$metaHash->{fastq}{$lane}{library}\\tSM:$metaHash->{sample}\\tPL:Illumina\' > $dir/$outName.sam; java -Xmx10g -jar \$PICARDROOT/SortSam.jar I=$dir/$outName.sam O=$dir/$outName.bam SO=coordinate CREATE_INDEX=true TMP_DIR=$dir/$outName.tmp; rm $dir/$outName.sam\n";
			print COMMAND "\necho phoenixPipe/$module-done\n";
			close COMMAND;

			`bash $dir/$outName.sub`;
		}
		else
		{
			warn " [NOVOALIGN] Skipping alignment: $dir/$outName.touch\n";
		}

		$metaHash->{bam}{$module}{lanes}{$outName}{file} = "$dir/$outName.bam";
		$metaHash->{bam}{$module}{lanes}{$outName}{index} = "$dir/$outName.bai";
		$metaHash->{bam}{$module}{lanes}{$outName}{hold_jid} = "$sgePre$outName";

		push (@{ $metaHash->{clean}{full} }, "$dir/$outName.bam");
		push (@{ $metaHash->{clean}{full} }, "$dir/$outName.bai");

		push (@{ $metaHash->{clean}{temp}{files} }, "$dir/$outName.tmp");
		push (@{ $metaHash->{clean}{temp}{hold_jids} }, "$sgePre$outName");

		$metaHash->{provenence}{"$dir/$outName.bam"}{command} = "$dir/$outName.cmd";
	}

}




sub doBWAmem
{
	my $module = $_[0];
	my $metaHash = $_[1];
	my $refHash = $_[2];

	my $pathRoot = buildPath($metaHash);
	my $dir = "$pathRoot/$module/lanes/";
	`mkdir -p $dir`;

	my $sgePre = "n$metaHash->{sge_uniq}";
	my $outName;

	my $files;
	my $command;
	my $hold_jid;
	my $qualFlag;

	my $convertQuals;
	my $cleanQuals;
	my $newFiles;

	for my $lane (sort keys %{ $metaHash->{fastq} })
	{

		$outName = "$metaHash->{fastq}{$lane}{library}_$metaHash->{fastq}{$lane}{barcode}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{run}";


		if ($metaHash->{fastq}{$lane}{read2} =~ /.*Not found/)
		{
			if ($metaHash->{tissue_type} eq "X")
			{
				$files = "$metaHash->{fastq}{$lane}{xenome}{read1}";
				$hold_jid = "-hold_jid $metaHash->{fastq}{$lane}{xenome}{hold_jid}";

				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $metaHash->{fastq}{$lane}{xenome}{read1});
			}
			else
			{
				$files = "$metaHash->{fastq}{$lane}{read1}";
				$hold_jid = "";

				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $metaHash->{fastq}{$lane}{read1});
			}
		}
		else
		{
			if ($metaHash->{tissue_type} eq "X")
			{
				$files = "$metaHash->{fastq}{$lane}{xenome}{read1} $metaHash->{fastq}{$lane}{xenome}{read2}";
				$hold_jid = "-hold_jid $metaHash->{fastq}{$lane}{xenome}{hold_jid}";

				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $metaHash->{fastq}{$lane}{xenome}{read1});
				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $metaHash->{fastq}{$lane}{xenome}{read2});
			}
			else
			{
				$files = "$metaHash->{fastq}{$lane}{read1} $metaHash->{fastq}{$lane}{read2}";
				$hold_jid = "";

				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $metaHash->{fastq}{$lane}{read1});
				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $metaHash->{fastq}{$lane}{read2});
			}
		}
		
		$convertQuals = "";
		$cleanQuals = "";
		$newFiles = "";
		if ($metaHash->{fastq}{$lane}{"quality_format"} eq "illumina")
		{
			$convertQuals = "module load seqtk/150618;";
			for my $f (split(/ /, $files))
			{
				$convertQuals .= " seqtk seq -VQ64 $f > ${f}.qualFix.fastq;";
				$cleanQuals .= "; rm ${f}.qualFix.fastq";
				$newFiles .= " ${f}.qualFix.fastq";
			}
			$files = $newFiles;
		}

		unless (-e "$dir/$outName.touch")
		{
			`touch $dir/$outName.touch`;		# so that future runs won't cause concurrent processing


			open (SUBFILE, ">$dir/$outName.sub") or die;
#			print SUBFILE "qsub -cwd -b y -N $sgePre$outName $hold_jid -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=24g,exclusive=1 -q $refHash->{sge_queue} \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
			print SUBFILE "qsub -cwd -b y -N $sgePre$outName $hold_jid -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=6g -pe smp 4 -q $refHash->{sge_queue} \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
			close SUBFILE;

			open (COMMAND, ">$dir/$outName.cmd") or die;
#			print COMMAND "$convertQuals module load $module; module load $refHash->{picard}; bwa mem -t `grep processor /proc/cpuinfo | wc -l` -M $refHash->{$module} $files -R '\@RG\\tID:$metaHash->{fastq}{$lane}{run}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{barcode}\\tPU:$metaHash->{fastq}{$lane}{run}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{barcode}\\tLB:$metaHash->{fastq}{$lane}{library}\\tSM:$metaHash->{sample}\\tPL:Illumina\' | java -Xmx1g -jar \$PICARDROOT/SortSam.jar I=/dev/stdin O=$dir/$outName.bam SO=coordinate CREATE_INDEX=true TMP_DIR=$dir/$outName.tmp $cleanQuals\n";
			print COMMAND "$convertQuals module load $module; module load $refHash->{picard}; bwa mem -t 4 -M $refHash->{$module} $files -R '\@RG\\tID:$metaHash->{fastq}{$lane}{run}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{barcode}\\tPU:$metaHash->{fastq}{$lane}{run}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{barcode}\\tLB:$metaHash->{fastq}{$lane}{library}\\tSM:$metaHash->{sample}\\tPL:Illumina\' | java -Xmx2g -jar \$PICARDROOT/SortSam.jar I=/dev/stdin O=$dir/$outName.bam SO=coordinate CREATE_INDEX=true TMP_DIR=$dir/$outName.tmp $cleanQuals\nmodule load samtools; samtools flagstat $dir/$outName.bam > $dir/$outName.bam.flagstat; zcat $files | wc -l > $dir/$outName.fastq_line_count\n";
			print COMMAND "\necho phoenixPipe/$module-done\n";
			close COMMAND;

			`bash $dir/$outName.sub`;
		}
		else
		{
			warn " [BWAMEM] Skipping alignment: $dir/$outName.touch\n";
		}

		$metaHash->{bam}{$module}{lanes}{$outName}{file} = "$dir/$outName.bam";
		$metaHash->{bam}{$module}{lanes}{$outName}{index} = "$dir/$outName.bai";
		$metaHash->{bam}{$module}{lanes}{$outName}{hold_jid} = "$sgePre$outName";

		push (@{ $metaHash->{clean}{full} }, "$dir/$outName.bam");
		push (@{ $metaHash->{clean}{full} }, "$dir/$outName.bai");

		push (@{ $metaHash->{clean}{temp}{files} }, "$dir/$outName.tmp");
		push (@{ $metaHash->{clean}{temp}{hold_jids} }, "$sgePre$outName");

		$metaHash->{provenence}{"$dir/$outName.bam"}{command} = "$dir/$outName.cmd";
	}

}




sub doBWAsai
{
	my $module = $_[0];
	my $metaHash = $_[1];
	my $refHash = $_[2];

	my $bwaReference = $refHash->{$module};

	my $pathRoot = buildPath($metaHash);
	my $dir = "$pathRoot/$module/lanes/";
	`mkdir -p $dir`;

	my $sgePre = "b$metaHash->{sge_uniq}";
	my $outName;

	my @reads;
	my $command;
	my $hold_jid;

	my $qualFlag;

	for my $lane (sort keys %{ $metaHash->{fastq} })
	{
		if ($metaHash->{fastq}{$lane}{read2} =~ /.*Not found/)
		{
			@reads = qw(read1);
		}
		else
		{
			@reads = qw(read1 read2);
		}

		if ($metaHash->{fastq}{$lane}{"quality_format"} eq "illumina")
		{
			$qualFlag = "-I";
		}
		else
		{
			$qualFlag = "";
		}

		if ($metaHash->{fastq}{$lane}{read1} =~ /fastq/)	# a little sloppy...
		{
			for my $read (@reads)
			{
				$outName = "$metaHash->{fastq}{$lane}{library}_$metaHash->{fastq}{$lane}{barcode}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{run}-$read";
		
				if (($metaHash->{tissue_type} eq "X") or ($metaHash->{tissue_type} eq "O"))
				{
#					$command = "module load $module; bwa aln $qualFlag -t `grep processor /proc/cpuinfo | wc -l` $bwaReference $metaHash->{fastq}{$lane}{xenome}{$read} > $dir/${outName}.sai\n";
					$command = "module load $module; bwa aln $qualFlag -t 4 $bwaReference $metaHash->{fastq}{$lane}{xenome}{$read} > $dir/${outName}.sai\nwc -l $metaHash->{fastq}{$lane}{xenome}{$read} > $dir/${outName}.fastq_line_count";
					$hold_jid = "-hold_jid $metaHash->{fastq}{$lane}{xenome}{hold_jid}";
	
					push (@{ $metaHash->{provenence}{"$dir/$outName.sai"}{input}}, $metaHash->{fastq}{$lane}{xenome}{$read});
				}
				else
				{
#					$command = "module load $module; bwa aln $qualFlag -t `grep processor /proc/cpuinfo | wc -l` $bwaReference $metaHash->{fastq}{$lane}{$read} > $dir/${outName}.sai\n";
					$command = "module load $module; bwa aln $qualFlag -t 4 $bwaReference $metaHash->{fastq}{$lane}{$read} > $dir/${outName}.sai\nzcat $metaHash->{fastq}{$lane}{$read} | wc -l > $dir/${outName}.fastq_line_count";
#					$command = "module load $module; bwa aln $qualFlag -t 24 $bwaReference $metaHash->{fastq}{$lane}{$read} > $dir/${outName}.sai\nzcat $metaHash->{fastq}{$lane}{$read} | wc -l > $dir/${outName}.fastq_line_count";
					$hold_jid = "";
	
					push (@{ $metaHash->{provenence}{"$dir/$outName.sai"}{input}}, $metaHash->{fastq}{$lane}{$read});
				}
	
				unless ((-e "$dir/${outName}.touch") or (-e "$dir/../collapsed/$metaHash->{sample}.touch"))
				{
					`touch $dir/${outName}.touch`;		# so that future runs won't cause concurrent processing
		
					open (SUBFILE, ">$dir/$outName.sub") or die;
#					print SUBFILE "qsub -cwd -b y -N $sgePre$outName $hold_jid -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=24g,exclusive=1 -q $refHash->{sge_queue} \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
					print SUBFILE "qsub -cwd -b y -N $sgePre$outName $hold_jid -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=6g -pe smp 4 -q $refHash->{sge_queue} \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
#					print SUBFILE "qsub -cwd -b y -N $sgePre$outName $hold_jid -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=120g -q $refHash->{sge_queue} \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
					close SUBFILE;
		
					open (COMMAND, ">$dir/$outName.cmd") or die;
					print COMMAND $command;
					print COMMAND "\necho phoenixPipe/$module-done\n";
					close COMMAND;
	
					`bash $dir/$outName.sub`;
				}
				else
				{
					warn " [BWA] Skipping sai generation: $dir/${outName}.touch\n";
				}
				
				$metaHash->{sai}{$module}{lanes}{$outName}{file} = "$dir/${outName}.sai";
				$metaHash->{sai}{$module}{lanes}{$outName}{hold_jid} = "$sgePre$outName";
	
				push (@{ $metaHash->{clean}{temp}{files} }, "$dir/$outName.sai");
				push (@{ $metaHash->{clean}{temp}{hold_jids} }, "$sgePre$outName");
	
				push (@{ $metaHash->{sub_list} }, "$dir/$outName.sub");
				push (@{ $metaHash->{cmd_list} }, "$dir/$outName.cmd");
	
				$metaHash->{provenence}{"$dir/$outName.sai"}{command} = "$dir/$outName.cmd";
			}
		}
		elsif ($metaHash->{fastq}{$lane}{read1} =~ /unaligned_bam/)
		{

			# no xenome support yet
			$hold_jid = "";

			for my $r (qw/1 2/)
			{
				$outName = "$metaHash->{fastq}{$lane}{library}_$metaHash->{fastq}{$lane}{barcode}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{run}-$r";

				unless ((-e "$dir/${outName}.touch") or (-e "$dir/../collapsed/$metaHash->{sample}.touch"))
				{
					`touch $dir/${outName}.touch`;
	
#					$command = "module load $module; bwa aln $qualFlag -t `grep processor /proc/cpuinfo | wc -l` $bwaReference -b$r $metaHash->{fastq}{$lane}{read1} > $dir/${outName}.sai\n";
					$command = "module load $module; bwa aln $qualFlag -t 4 $bwaReference -b$r $metaHash->{fastq}{$lane}{read1} > $dir/${outName}.sai\n";
	
					open (SUBFILE, ">$dir/$outName.sub") or die;
#					print SUBFILE "qsub -cwd -b y -N $sgePre$outName $hold_jid -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=24g,exclusive=1 -q $refHash->{sge_queue} \" bash $dir/$outName.cmd\" > $dir/$outName.log\n";
					print SUBFILE "qsub -cwd -b y -N $sgePre$outName $hold_jid -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=6g -pe smp 4 -q $refHash->{sge_queue} \" bash $dir/$outName.cmd\" > $dir/$outName.log\n";
					close SUBFILE;
	
					open (COMMAND, ">$dir/$outName.cmd") or die;
					print COMMAND $command;
					print COMMAND "\necho phoenixPipe/$module-done\n";
					close COMMAND;
	
					`bash $dir/$outName.sub`;
				}
				else
				{
					warn " [BWA] Skipping sai generation: $dir/${outName}.touch\n";
				}

				$metaHash->{sai}{$module}{lanes}{$outName}{file} = "$dir/${outName}.sai";
				$metaHash->{sai}{$module}{lanes}{$outName}{hold_jid} = "$sgePre$outName";

				push (@{ $metaHash->{clean}{temp}{files} }, "$dir/$outName.sai");
				push (@{ $metaHash->{clean}{temp}{hold_jids} }, "$sgePre$outName");

				push (@{ $metaHash->{sub_list} }, "$dir/$outName.sub");
				push (@{ $metaHash->{cmd_list} }, "$dir/$outName.cmd");

				$metaHash->{provenence}{"$dir/$outName.sai"}{command} = "$dir/$outName.cmd";
			}
		}
	

	}



}

sub doBWA
{
	my $module = $_[0];
	my $metaHash = $_[1];
	my $refHash = $_[2];

	my $bwaReference = $refHash->{$module};

	my $pathRoot = buildPath($metaHash);
	my $dir = "$pathRoot/$module/lanes/";
	`mkdir -p $dir`;

	my $sgePre = "b$metaHash->{sge_uniq}";
	my $outName;
	my $command;
	my $provenence = "";

	my $fastqs;
	my $sais;
	my $hold_jids;
	my $bwaMode;

	for my $lane (sort keys %{ $metaHash->{fastq} })
	{
		if ($metaHash->{fastq}{$lane}{read1} =~ /fastq/)
		{
			$outName = "$metaHash->{fastq}{$lane}{library}_$metaHash->{fastq}{$lane}{barcode}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{run}";
	
			if ($metaHash->{fastq}{$lane}{read2} =~ /.*Not found/)
			{
				$bwaMode = "samse";
				if (($metaHash->{tissue_type} eq "X") or ($metaHash->{tissue_type} eq "O"))
				{
					$fastqs = "$metaHash->{fastq}{$lane}{xenome}{read1}";
				}
				else
				{
					$fastqs = "$metaHash->{fastq}{$lane}{read1}";
				}
				$sais = "$metaHash->{sai}{$module}{lanes}{\"$outName-read1\"}{file}";
				$hold_jids = "$metaHash->{sai}{$module}{lanes}{\"$outName-read1\"}{hold_jid}";
	
				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, "$metaHash->{sai}{$module}{lanes}{\"$outName-read1\"}{file}");
			}
			else
			{
				$bwaMode = "sampe";
	
				if (($metaHash->{tissue_type} eq "X") or ($metaHash->{tissue_type} eq "O"))
				{
					$fastqs = "$metaHash->{fastq}{$lane}{xenome}{read1} $metaHash->{fastq}{$lane}{xenome}{read2}";
				}
				else
				{
					$fastqs = "$metaHash->{fastq}{$lane}{read1} $metaHash->{fastq}{$lane}{read2}";
				}
				$sais = "$metaHash->{sai}{$module}{lanes}{\"$outName-read1\"}{file} $metaHash->{sai}{$module}{lanes}{\"$outName-read2\"}{file}";
				$hold_jids = "$metaHash->{sai}{$module}{lanes}{\"$outName-read1\"}{hold_jid},$metaHash->{sai}{$module}{lanes}{\"$outName-read2\"}{hold_jid}";
				
				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, "$metaHash->{sai}{$module}{lanes}{\"$outName-read1\"}{file}");
				push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, "$metaHash->{sai}{$module}{lanes}{\"$outName-read2\"}{file}");
			}
	
			$command = "module load $module; module load $refHash->{picard}; module load $refHash->{samtools}; bwa $bwaMode $bwaReference $sais $fastqs | java -Xmx2g -jar \$PICARDROOT/AddOrReplaceReadGroups.jar RGID=$metaHash->{fastq}{$lane}{run}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{barcode} RGLB=$metaHash->{fastq}{$lane}{library} RGPL=Illumina RGPU=$metaHash->{fastq}{$lane}{run}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{barcode} RGSM=$metaHash->{sample} VALIDATION_STRINGENCY=SILENT I=/dev/stdin O=$dir/$outName.bam SO=coordinate CREATE_INDEX=true TMP_DIR=$dir/$outName.tmp; rm $sais; module load samtools; samtools flagstat $dir/$outName.bam > $dir/$outName.bam.flagstat\n";
	
			unless ((-e "$dir/$outName.touch") or (-e "$dir/../collapsed/$metaHash->{sample}.touch"))
			{
				`touch $dir/$outName.touch`;		# so that future runs won't cause concurrent processing
	
	
				open (SUBFILE, ">$dir/$outName.sub") or die;
				print SUBFILE "qsub -cwd -b y -N $sgePre$outName -hold_jid $hold_jids -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=16g -q $refHash->{sge_queue} \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
				close SUBFILE;
	
				open (COMMAND, ">$dir/$outName.cmd") or die;
				print COMMAND $command;
				print COMMAND "\necho phoenixPipe/$module-done\n";
				close COMMAND;
	
				`bash $dir/$outName.sub`;
			}
			else
			{
				warn " [BWA] Skipping alignment: $dir/$outName.touch\n";
			}
	
	
			$metaHash->{bam}{$module}{lanes}{$outName}{file} = "$dir/$outName.bam";
			$metaHash->{bam}{$module}{lanes}{$outName}{index} = "$dir/$outName.bai";
			$metaHash->{bam}{$module}{lanes}{$outName}{hold_jid} = "$sgePre$outName";
	
			push (@{ $metaHash->{clean}{full} }, "$dir/$outName.bam");
			push (@{ $metaHash->{clean}{full} }, "$dir/$outName.bai");
	
			push (@{ $metaHash->{clean}{temp}{files} }, "$dir/$outName.tmp");
			push (@{ $metaHash->{clean}{temp}{hold_jids} }, "$sgePre$outName");
	
			$metaHash->{provenence}{"$dir/$outName.bam"}{command} = "$dir/$outName.cmd";
		}
		elsif ($metaHash->{fastq}{$lane}{read1} =~ /unaligned_bam/)
		{
			$outName = "$metaHash->{fastq}{$lane}{library}_$metaHash->{fastq}{$lane}{barcode}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{run}";

			# no xenome
			
			$fastqs = "$metaHash->{fastq}{$lane}{read1} $metaHash->{fastq}{$lane}{read1}";		# bam file twice, apparently
			$sais = "$metaHash->{sai}{$module}{lanes}{\"$outName-1\"}{file} $metaHash->{sai}{$module}{lanes}{\"$outName-2\"}{file}";
			$hold_jids = "$metaHash->{sai}{$module}{lanes}{\"$outName-1\"}{hold_jid},$metaHash->{sai}{$module}{lanes}{\"$outName-2\"}{hold_jid}";

			push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, "$metaHash->{sai}{$module}{lanes}{\"$outName-1\"}{file}");
			push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, "$metaHash->{sai}{$module}{lanes}{\"$outName-2\"}{file}");

			$bwaMode = "sampe";
			
			$command = "module load $module; module load $refHash->{picard}; module load $refHash->{samtools}; bwa $bwaMode $bwaReference $sais $fastqs | java -Xmx2g -jar \$PICARDROOT/AddOrReplaceReadGroups.jar RGID=$metaHash->{fastq}{$lane}{run}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{barcode} RGLB=$metaHash->{fastq}{$lane}{library} RGPL=Illumina RGPU=$metaHash->{fastq}{$lane}{run}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{barcode} RGSM=$metaHash->{sample} VALIDATION_STRINGENCY=SILENT I=/dev/stdin O=$dir/$outName.bam SO=coordinate CREATE_INDEX=true TMP_DIR=$dir/$outName.tmp; rm $sais\n";
	
			unless ((-e "$dir/$outName.touch") or (-e "$dir/../collapsed/$metaHash->{sample}.touch"))
			{
				`touch $dir/$outName.touch`;		# so that future runs won't cause concurrent processing
	
	
				open (SUBFILE, ">$dir/$outName.sub") or die;
				print SUBFILE "qsub -cwd -b y -N $sgePre$outName -hold_jid $hold_jids -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=16g -q $refHash->{sge_queue} \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
				close SUBFILE;
	
				open (COMMAND, ">$dir/$outName.cmd") or die;
				print COMMAND $command;
				print COMMAND "\necho phoenixPipe/$module-done\n";
				close COMMAND;
	
				`bash $dir/$outName.sub`;

				$metaHash->{bam}{$module}{lanes}{$outName}{file} = "$dir/$outName.bam";
				$metaHash->{bam}{$module}{lanes}{$outName}{index} = "$dir/$outName.bai";
				$metaHash->{bam}{$module}{lanes}{$outName}{hold_jid} = "$sgePre$outName";

				push (@{ $metaHash->{clean}{full} }, "$dir/$outName.bam");
				push (@{ $metaHash->{clean}{full} }, "$dir/$outName.bai");
				push (@{ $metaHash->{clean}{temp}{files} }, "$dir/$outName.tmp");
				push (@{ $metaHash->{clean}{temp}{hold_jids} }, "$sgePre$outName");
				$metaHash->{provenence}{"$dir/$outName.bam"}{command} = "$dir/$outName.cmd";
			}
			else
			{
				warn " [BWA] Skipping alignment: $dir/$outName.touch\n";
			}
		}
	}
}



sub doFilter
{
	my $module = $_[0];
	my $metaHash = $_[1];
	my $refHash = $_[2];

	# filter each lane
	my $pathRoot = buildPath($metaHash);
	my $dir = "$pathRoot/$module/filtered";
	`mkdir -p $dir`;

	my $sgePre = "f$metaHash->{sge_uniq}";

	my $outName;
	my $inFile;


	for $outName (sort keys %{ $metaHash->{bam}{$module}{lanes} })
	{
		$inFile = $metaHash->{bam}{$module}{lanes}{$outName}{file};
		unless (-e "$dir/$outName.touch")
		{
			`touch $dir/$outName.touch`;

			open (SUBFILE, ">$dir/$outName.sub") or die;
			print SUBFILE "qsub -cwd -b y -N $sgePre$outName -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=2g -q $refHash->{sge_queue} -hold_jid $metaHash->{bam}{$module}{lanes}{$outName}{hold_jid} \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
			close SUBFILE;

			open (COMMAND, ">$dir/$outName.cmd") or die;
			if ($module =~ /novocraft/)
			{
				print COMMAND "ln -s $inFile $dir/$outName.bam\n";		# not doing any filtering!
#				print COMMAND "module load $refHash->{samtools}; samtools view -u -F 4 $inFile > $dir/$outName.bam\n";
#				print COMMAND "module load $refHash->{samtools}; samtools view -u -F 4 $inFile | samtools view -u -F 256 - | samtools view -b -q 30 - > $dir/$outName.bam\n";		# probably don't need the -F 256
			}
			elsif ($module =~ /bwa/)
			{
				print COMMAND "ln -s $inFile $dir/$outName.bam\n";		# not doing any filtering!
#				print COMMAND "module load $refHash->{samtools}; samtools view -u -F 4 $inFile > $dir/$outName.bam\n";
#				print COMMAND "module load $refHash->{samtools}; samtools view -u -F 4 $inFile | samtools view -b -q 5 - > $dir/$outName.bam\n";
			}
			print COMMAND "\necho phoenixPipe/filter-done\n";
			close COMMAND;

			`bash $dir/$outName.sub`;

		}
		else
		{
			warn " [FILTER] Skipping: $dir/$outName.touch\n";
		}

		$metaHash->{bam}{$module}{filtered}{$outName}{file} = "$dir/$outName.bam";
		$metaHash->{bam}{$module}{filtered}{$outName}{hold_jid} = "$sgePre$outName";

		push (@{ $metaHash->{clean}{temp}{files} }, "$dir/$outName.bam");
		push (@{ $metaHash->{clean}{temp}{hold_jids} }, "$sgePre$outName");

		$metaHash->{provenence}{"$dir/$outName.bam"}{command} = "$dir/$outName.cmd";
		push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $inFile);
	}
	close SUBFILE;
}




sub doGenotype
{
	my $module = $_[0];
	my $metaHash = $_[1];
	my $refHash = $_[2];

	my $pathRoot = buildPath($metaHash);
	my $dir = "$pathRoot/$module/sampleIdentity/genotype";
	`mkdir -p $dir`;

	my $sgePre = "sid$metaHash->{sge_uniq}";

	my $sidBed = $refHash->{"sid bed file"};

	my $robCaller = $refHash->{"rob caller script"};
	my $robGenotyper = $refHash->{"rob genotyper script"};

	my $perlLib = $refHash{'perl lib path'};

	my $outName;
	my $inFile;


	for $outName (sort keys %{ $metaHash->{bam}{$module}{lanes} })
	{
		$inFile = $metaHash->{bam}{$module}{lanes}{$outName}{file};
		unless (-e "$dir/$outName.touch")
		{
			`touch $dir/$outName.touch`;

			open (SUBFILE, ">$dir/$outName.sub") or die;
			print SUBFILE "qsub -cwd -b y -N $sgePre$outName -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=8g -q $refHash->{sge_queue} -hold_jid $metaHash->{bam}{$module}{lanes}{$outName}{hold_jid} \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
			close SUBFILE;

			open (COMMAND, ">$dir/$outName.cmd") or die;
			print COMMAND "module load samtools; $perlLib; samtools view -F 4 $inFile -L $sidBed | $robCaller -q 0 -d 0 -f 0 | $robGenotyper $sidBed > $dir/$outName.genotype\n";
			print COMMAND "\necho phoenixPipe/genotype-done\n";
			close COMMAND;

			`bash $dir/$outName.sub`;

		}
		else
		{
			warn " [GENOTYPE] Skipping: $dir/$outName.touch\n";
		}

		$metaHash->{bam}{$module}{genotype}{$outName}{file} = "$dir/$outName.genotype";
		$metaHash->{bam}{$module}{genotype}{$outName}{hold_jid} = "$sgePre$outName";

		$metaHash->{provenence}{"$dir/$outName.genotype"}{command} = "$dir/$outName.cmd";
		push (@{ $metaHash->{provenence}{"$dir/$outName.genotype"}{input}}, $inFile);
	}
	close SUBFILE;
}




sub doGendertype
{
	my $module = $_[0];
	my $metaHash = $_[1];
	my $refHash = $_[2];

	my $pathRoot = buildPath($metaHash);
	my $dir = "$pathRoot/$module/sampleIdentity/gendertype";
	`mkdir -p $dir`;

	my $sgePre = "sid$metaHash->{sge_uniq}";


	my $outName;
	my $inFile;


	for $outName (sort keys %{ $metaHash->{bam}{$module}{lanes} })
	{
		$inFile = $metaHash->{bam}{$module}{lanes}{$outName}{file};
		unless (-e "$dir/$outName.touch")
		{
			`touch $dir/$outName.touch`;

			open (SUBFILE, ">$dir/$outName.sub") or die;
			print SUBFILE "qsub -cwd -b y -N $sgePre$outName -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=2g -q $refHash->{sge_queue} -hold_jid $metaHash->{bam}{$module}{lanes}{$outName}{hold_jid} \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
			close SUBFILE;

			open (COMMAND, ">$dir/$outName.cmd") or die;
			print COMMAND "module load samtools; samtools view -F 4 $inFile chrX chrY | cut -f 3 | uniq -c > $dir/$outName.gendertype\n";
			print COMMAND "\necho phoenixPipe/gendertype-done\n";
			close COMMAND;

			`bash $dir/$outName.sub`;

		}
		else
		{
			warn " [GENDERTYPE] Skipping: $dir/$outName.touch\n";
		}

		$metaHash->{bam}{$module}{gendertype}{$outName}{file} = "$dir/$outName.gendertype";
		$metaHash->{bam}{$module}{gendertype}{$outName}{hold_jid} = "$sgePre$outName";

		$metaHash->{provenence}{"$dir/$outName.gendertype"}{command} = "$dir/$outName.cmd";
		push (@{ $metaHash->{provenence}{"$dir/$outName.gendertype"}{input}}, $inFile);
	}
	close SUBFILE;
}




sub doSomatotype
{
	my $module = $_[0];
	my $metaHash = $_[1];
	my $tumourHash = $_[2];
	my $refHash = $_[3];

	my $pathRoot = buildPath($metaHash);
	my $dir = "$pathRoot/$module/sampleIdentity/somatotype";
	`mkdir -p $dir`;

	my $sgePre = "sid$metaHash->{sge_uniq}";

	my $sidBed = $tumourHash->{vcf}{final}{file};
	my $holdJid = $tumourHash->{vcf}{final}{hold_jid};

	my $robCaller = $refHash->{"rob caller script"};
	my $robGenotyper = $refHash->{"rob genotyper script"};


	my $outName;
	my $inFile;


	for $outName (sort keys %{ $metaHash->{bam}{$module}{lanes} })
	{
		$inFile = $metaHash->{bam}{$module}{lanes}{$outName}{file};
		unless (-e "$dir/$outName.touch")
		{
			`touch $dir/$outName.touch`;

			open (SUBFILE, ">$dir/$outName.sub") or die;
			print SUBFILE "qsub -cwd -b y -N $sgePre$outName -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=8g -q $refHash->{sge_queue} -hold_jid $metaHash->{bam}{$module}{lanes}{$outName}{hold_jid},$holdJid \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
			close SUBFILE;

			open (COMMAND, ">$dir/$outName.cmd") or die;
			print COMMAND "module load samtools; samtools view -F 4 $inFile -L $sidBed | $robCaller -q 0 -d 0 -f 0 | $robGenotyper $sidBed > $dir/$outName.somatotype\n";
			print COMMAND "\necho phoenixPipe/somatotype-done\n";
			close COMMAND;

			`bash $dir/$outName.sub`;

		}
		else
		{
			warn " [SOMATOTYPE] Skipping: $dir/$outName.touch\n";
		}

		$metaHash->{bam}{$module}{somatotype}{$outName}{file} = "$dir/$outName.somatotype";
		$metaHash->{bam}{$module}{somatotype}{$outName}{hold_jid} = "$sgePre$outName";

		$metaHash->{provenence}{"$dir/$outName.somatotype"}{command} = "$dir/$outName.cmd";
		push (@{ $metaHash->{provenence}{"$dir/$outName.somatotype"}{input}}, $inFile);
	}
	close SUBFILE;
}







sub doMergeCollapse
{
	my $module = $_[0];
	my $rmDups = $_[1];
	my $metaHash = $_[2];
	my $refHash = $_[3];

	# merge and collapse lanes together
	my $pathRoot = buildPath($metaHash);
	my $dir = "$pathRoot/$module/collapsed";
	`mkdir -p $dir`;

	my $sgePre = "c$metaHash->{sge_uniq}";
	my $inputString = "";
	my $holdString = "";
	my $filterName;
	my $inFile;

	my $outName = $metaHash->{sample};

	for $filterName (sort keys %{ $metaHash->{bam}{$module}{filtered} })
	{
		$inFile = $metaHash->{bam}{$module}{filtered}{$filterName}{file};
		$inputString = "$inputString I=$inFile";
		$holdString = "$metaHash->{bam}{$module}{filtered}{$filterName}{hold_jid},$holdString";

		push (@{ $metaHash->{provenence}{"$dir/$outName.bam"}{input}}, $inFile);
	}

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;


		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -N $sgePre$outName -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=16g -q $refHash->{sge_queue} -hold_jid $holdString \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "module load $refHash->{picard}; java -Xmx4g -jar \$PICARDROOT/MarkDuplicates.jar TMP_DIR=$dir/picardTmp $inputString O=$dir/$outName.bam METRICS_FILE=$dir/$outName.mmm ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=$rmDups CREATE_INDEX=true\nmodule load samtools; samtools flagstat $dir/$outName.bam > $dir/$outName.bam.flagstat";
		print COMMAND "\necho phoenixPipe/collapse-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [COLLAPSE] Skipping: $dir/$outName.touch\n";
	}

	$metaHash->{bam}{$module}{final}{file} = "$dir/$outName.bam";
	$metaHash->{bam}{$module}{final}{index} = "$dir/$outName.bai";
	$metaHash->{bam}{$module}{final}{hold_jid} = "$sgePre$outName";
	
	push (@{ $metaHash->{clean}{full} }, "$dir/$outName.bam");

	$metaHash->{provenence}{"$dir/$outName.bam"}{command} = "$dir/$outName.cmd";
}




sub doSamStatsLane
{
	my $bamModule = $_[0];
	my $metaHash = $_[1];
	my $refHash = $_[2];

	my $pathRoot = buildPath($metaHash);
	my $dir = "$pathRoot/$bamModule/json/";
	`mkdir -p $dir`;

	my $sgePre = "j$metaHash->{sge_uniq}";
	my $outName;
	my $bam;
	my $hold_jid;
	my $target = $refHash->{"$metaHash->{sample_type}-target"};

	my ($sample_group, $sample, $library, $barcode, $instrument, $run_name, $lane_num, $sequencing_type);

	my $perlLib = $refHash{'perl lib path'};

	for my $lane (sort keys %{ $metaHash->{fastq} })
	{
		$outName = "$metaHash->{fastq}{$lane}{library}_$metaHash->{fastq}{$lane}{barcode}_$metaHash->{fastq}{$lane}{lane}_$metaHash->{fastq}{$lane}{run}";
		$bam = "$metaHash->{bam}{$bamModule}{lanes}{$outName}{file}";
		$hold_jid = "-hold_jid $metaHash->{bam}{$bamModule}{lanes}{$outName}{hold_jid}";

		$sample_group = $metaHash->{sample_group};
		$sample = $metaHash->{sample};
		$library = $metaHash->{fastq}{$lane}{library};
		$barcode = $metaHash->{fastq}{$lane}{barcode};
		$instrument = $metaHash->{fastq}{$lane}{instrument};
		$run_name = $metaHash->{fastq}{$lane}{run};
		$lane_num = $metaHash->{fastq}{$lane}{lane};
		$sequencing_type = $metaHash->{sample_type};

		unless (-e "$dir/$outName.touch")
		{
			`touch $dir/$outName.touch`;

			open (SUBFILE, ">$dir/$outName.sub") or die;
			print SUBFILE "qsub -cwd -b y  -N $sgePre$outName -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=2g -q $refHash->{sge_queue} $hold_jid \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
			close SUBFILE;

			open (COMMAND, ">$dir/$outName.cmd") or die;
			print COMMAND "$perlLib; module load $refHash->{samtools}; samtools view $bam | $refHash->{'samStats.pl'} -q 20 -b $bam -r $target -j '{\"sample\":\"$sample\",\"library\":\"$library\",\"barcode\":\"$barcode\",\"instrument\":\"$instrument\",\"run name\":\"$run_name\",\"lane\":\"$lane_num\",\"sample group\":\"$sample_group\",\"sequencing type\":\"$sequencing_type\"}' > $dir/$outName.json\n";
			print COMMAND "\necho phoenixPipe/samStats-done\n";
			close COMMAND;

			`bash $dir/$outName.sub`;
		}
		else
		{
			warn " [SAMSTATS] Skipping: $dir/$outName.touch\n";
		}

		$metaHash->{bam}{$bamModule}{lanes}{$outName}{json}{file} = "$dir/$outName.json";
		$metaHash->{bam}{$bamModule}{lanes}{$outName}{json}{hold_jid} = "$sgePre$outName";
		
		$metaHash->{provenence}{"$dir/$outName.json"}{command} = "$dir/$outName.cmd";
		push (@{ $metaHash->{provenence}{"$dir/$outName.json"}{input}}, $bam);

	}
}



sub doSamStatsMerge
{
	my $bamModule = $_[0];
	my $metaHash = $_[1];
	my $refHash = $_[2];

	my $pathRoot = buildPath($metaHash);
	my $dir = "$pathRoot/$bamModule/json/";
	`mkdir -p $dir`;

	my $sgePre = "j$metaHash->{sge_uniq}";
	my $outName;
	my $bam;
	my $hold_jid;
	my $target = $refHash->{"$metaHash->{sample_type}-target"};

	my ($sample_group, $sample, $library, $sequencing_type);

	$outName = "$metaHash->{sample}";
	$bam = "$metaHash->{bam}{$bamModule}{final}{file}";
	$hold_jid = "-hold_jid $metaHash->{bam}{$bamModule}{final}{hold_jid}";

	$sample_group = $metaHash->{sample_group};
	$sample = $metaHash->{sample};
	$library = "merged";
	$sequencing_type = $metaHash->{sample_type};

	my $perlLib = $refHash{'perl lib path'};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;

		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -N $sgePre$outName -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=2g -q $refHash->{sge_queue} $hold_jid \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "$perlLib; module load $refHash->{samtools}; samtools view $bam | $refHash->{'samStats.pl'} -q 20 -b $bam -r $target -j '{\"sample\":\"$sample\",\"library\":\"$library\",\"sample group\":\"$sample_group\",\"sequencing type\":\"$sequencing_type\"}' > $dir/$outName.json\n";
		print COMMAND "\necho phoenixPipe/samStats-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [SAMSTATS] Skipping: $dir/$outName.touch\n";
	}

	$metaHash->{bam}{$bamModule}{final}{json}{file} = "$dir/$outName.json";
	$metaHash->{bam}{$bamModule}{final}{json}{hold_jid} = "$sgePre$outName";
	
	$metaHash->{provenence}{"$dir/$outName.json"}{command} = "$dir/$outName.cmd";
	push (@{ $metaHash->{provenence}{"$dir/$outName.json"}{input}}, $bam);

}



sub doStrelka
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $normalHash = $_[2];
	my $tumourHash = $_[3];
	my $refHash = $_[4];

	my $seqType = $tumourHash->{sample_type};

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/$module";
	`mkdir -p $dir`;

	my $sgePre = "s$tumourHash->{sge_uniq}";

	my $normalBam = $normalHash->{bam}{$bamModule}{final}{file};
	my $tumourBam = $tumourHash->{bam}{$bamModule}{final}{file};
	my $normalBai;
	my $tumourBai;
	my $holdJid = "$normalHash->{bam}{$bamModule}{final}{hold_jid},$tumourHash->{bam}{$bamModule}{final}{hold_jid}";

	my $fastaRef = $refHash->{hg19_random};

	my $config;

	my $outName = $tumourHash->{sample};


	if ($seqType eq "exome")
	{
		$config = $refHash->{"$module-exome"};
	}
	elsif ($seqType eq "wgs")
	{
		$config = $refHash->{"$module-wgs"};
	}
	else
	{
		die " [STRELKA] Aborting! The sequence type must be either exome or wgs!\n";
	}

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=16g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;
	
		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "module load $module; cp $config $dir/config.ini; configureStrelkaWorkflow.pl --normal $normalBam --tumor $tumourBam --ref $fastaRef --config $dir/config.ini --output-dir $dir/analysis; cd $dir/analysis; qmake -inherit -- -j 16; qmake -inherit -- -j 16\n";		# double qmake to catch rare failures
		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;


		# fix index names dagnabbit
		unless (-e "$normalBam.bai")
		{
			$normalBai = $normalBam;
			$normalBai =~ s/\.bam$/\.bai/;
			$normalBai =~ s/.*\///;
			`ln -s $normalBai $normalBam.bai`;
		}

		unless (-e "$tumourBam.bai")
		{
			$tumourBai = $tumourBam;
			$tumourBai =~ s/\.bam$/\.bai/;
			$tumourBai =~ s/.*\///;
			`ln -s $tumourBai $tumourBam.bai`;
		}

		# links to results, why not?
		`ln -s ./analysis/results/passed.somatic.snvs.vcf $dir/strelka_$outName-passed.somatic.snvs.vcf`;
		`ln -s ./analysis/results/passed.somatic.indels.vcf $dir/strelka_$outName-passed.somatic.indels.vcf`;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [STRELKA] Skipping: $dir/$outName.touch\n";
	}

	$tumourHash->{vcf}{$bamModule}{$module}{snv}{file} = "$dir/strelka_$outName-passed.somatic.snvs.vcf";
	$tumourHash->{vcf}{$bamModule}{$module}{snv}{hold_jid} = "$sgePre$outName";

	$tumourHash->{vcf}{$bamModule}{$module}{indel}{file} = "$dir/strelka_$outName-passed.somatic.indels.vcf";

	$tumourHash->{provenence}{"$dir/strelka_$outName-passed.somatic.snvs.vcf"}{command} = "$dir/$outName.cmd";
	push (@{ $tumourHash->{provenence}{"$dir/strelka_$outName-passed.somatic.snvs.vcf"}{input}}, $normalBam);
	push (@{ $tumourHash->{provenence}{"$dir/strelka_$outName-passed.somatic.snvs.vcf"}{input}}, $tumourBam);

	$tumourHash->{provenence}{"$dir/strelka_$outName-passed.somatic.indels.vcf"}{command} = "$dir/$outName.cmd";
	push (@{ $tumourHash->{provenence}{"$dir/strelka_$outName-passed.somatic.indels.vcf"}{input}}, $normalBam);
	push (@{ $tumourHash->{provenence}{"$dir/strelka_$outName-passed.somatic.indels.vcf"}{input}}, $tumourBam);

}




sub doMutect
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $normalHash = $_[2];
	my $tumourHash = $_[3];
	my $refHash = $_[4];

	my $seqType = $tumourHash->{sample_type};

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/$module";
	`mkdir -p $dir`;

	my $sgePre = "m$tumourHash->{sge_uniq}";

	my $normalBam = $normalHash->{bam}{$bamModule}{final}{file};
	my $tumourBam = $tumourHash->{bam}{$bamModule}{final}{file};
	my $holdJid = "$normalHash->{bam}{$bamModule}{final}{hold_jid},$tumourHash->{bam}{$bamModule}{final}{hold_jid}";

	my $fastaRef = $refHash->{"hg19_random"};
	my $dbsnpRef = $refHash->{"dbSNP"};
	my $cosmicRef = $refHash->{"cosmic"};


	my $outName = $tumourHash->{sample};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=16g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "module load $module; java -Xmx10g -Djava.io.tmpdir=$dir/tmp -jar \$MUTECTROOT/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence $fastaRef --dbsnp $dbsnpRef --cosmic $cosmicRef --input_file:tumor $tumourBam --input_file:normal $normalBam --out $dir/mutect.call_stats.txt --coverage_file $dir/mutect.coverage.wig.txt --vcf $dir/mutect_$outName.vcf\n";
		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [MUTECT] Skipping: $dir/$outName.touch\n";
	}

	$tumourHash->{vcf}{$bamModule}{$module}{snv}{file} = "$dir/mutect_$outName.vcf";
	$tumourHash->{vcf}{$bamModule}{$module}{snv}{hold_jid} = "$sgePre$outName";

	$tumourHash->{provenence}{"$dir/mutect_$outName.vcf"}{command} = "$dir/$outName.cmd";
	push (@{ $tumourHash->{provenence}{"$dir/mutect_$outName.vcf"}{input}}, $normalBam);
	push (@{ $tumourHash->{provenence}{"$dir/mutect_$outName.vcf"}{input}}, $tumourBam);

}



sub doMutect_split
{
    my $module = $_[0];
    my $bamModule = $_[1];
    my $normalHash = $_[2];
    my $tumourHash = $_[3];
    my $refHash = $_[4];

    my $seqType = $tumourHash->{sample_type};

    # put somatic calls in the tumour dir
    my $pathRoot = buildPath($tumourHash);
    my $dir = "$pathRoot/$bamModule/$module";
    `mkdir -p $dir`;

    my $splitdir = "$dir/tempsplit";
    `mkdir -p $splitdir`;

    my $sgePre = "m$tumourHash->{sge_uniq}";

    my $normalBam = $normalHash->{bam}{$bamModule}{final}{file};
    my $tumourBam = $tumourHash->{bam}{$bamModule}{final}{file};
    my $holdJid = "$normalHash->{bam}{$bamModule}{final}{hold_jid},$tumourHash->{bam}{$bamModule}{final}{hold_jid}";

    my $fastaRef = $refHash->{"hg19_random"};
    my $dbsnpRef = $refHash->{"dbSNP"};
    my $cosmicRef = $refHash->{"cosmic"};

    my $outName = $tumourHash->{sample};

    my %chrom=map{split /\t/, $_} split /\n/,`cat $fastaRef.fai | cut -f 1,2`;

	my $commandFileList = "";

    unless(-e "$splitdir/$outName.touch")
    {
        `touch $splitdir/$outName.touch`;
        my @vcf;
        open (SUBFILE, ">$splitdir/$outName.sub") or die;
        while(my($chrom,$chrlen) = each %chrom)
        {
            print SUBFILE "qsub -cwd -b y -l h_vmem=16g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $splitdir/$outName.$chrom.log -o $splitdir/$outName.$chrom.log \"bash $splitdir/$outName.$chrom.cmd\" > $splitdir/$outName.$chrom.log\n";

            open (COMMAND, ">$splitdir/$outName.$chrom.cmd") or die;
            print COMMAND "module load $module; ".
                        "java -Xmx10g -Djava.io.tmpdir=$splitdir/tmp_$chrom -jar \$MUTECTROOT/muTect-1.1.4.jar ".
                        "--analysis_type MuTect ".
                        "--reference_sequence $fastaRef ".
                        "--dbsnp $dbsnpRef ".
                        "--cosmic $cosmicRef ".
                        "--input_file:tumor $tumourBam ".
                        "--input_file:normal $normalBam ".
                        "--intervals $chrom:1-$chrlen ".
                        "--out $splitdir/mutect.call_stats.$chrom.txt ".
                        "--coverage_file $splitdir/mutect.coverage.wig.$chrom.txt ".
                        "--vcf $splitdir/mutect_$outName.$chrom.vcf\n";
			print COMMAND "\necho phoenixPipe/$module-done\n";
            close COMMAND;

            push(@vcf,"$splitdir/mutect_$outName.$chrom.vcf");

        }
        close SUBFILE;
        `bash $splitdir/$outName.sub`;


        $holdJid="$sgePre$outName";  ### hold on the name used for mutect analysis
        open (SUBFILE, ">$dir/$outName.merge.sub") or die;
        print SUBFILE "qsub -cwd -b y -l h_vmem=16g -q $refHash->{sge_queue} -N $sgePre$outName.merge -hold_jid $holdJid -e $dir/$outName.merge.log -o $dir/$outName.merge.log \"bash $dir/$outName.merge.cmd\" > $dir/$outName.merge.log\n";
        close SUBFILE;



        open (COMMAND, ">$dir/$outName.merge.cmd") or die;
        print COMMAND "module load vcftools/0.1.14; mkdir $splitdir/tmp; TMPDIR=$splitdir/tmp; vcf=`ls $splitdir/*vcf`;vcf-concat \$vcf | vcf-sort > $dir/mutect_$outName.merged.vcf\n";		# hardcoded module - thanks Larry!
		print COMMAND "\necho phoenixPipe/$module-done\n";
        close COMMAND;

        `bash $dir/$outName.merge.sub`;
    }
    else
    {
        warn " [MUTECT_SPLIT] Skipping: $splitdir/$outName.touch\n";
    }


    $tumourHash->{vcf}{$bamModule}{$module}{snv}{file} = "$dir/mutect_$outName.merged.vcf";
    $tumourHash->{vcf}{$bamModule}{$module}{snv}{hold_jid} = "$sgePre$outName.merge";




	while(my($chrom,$chrlen) = each %chrom)
	{
		$commandFileList .= "$splitdir/$outName.$chrom.cmd,";
	}
	$commandFileList .= "$dir/$outName.merge.cmd";	# last command, so no comma

    $tumourHash->{provenence}{"$dir/mutect_$outName.merged.vcf"}{command} = $commandFileList;
    push (@{ $tumourHash->{provenence}{"$dir/mutect_$outName.merged.vcf"}{input}}, $normalBam);
    push (@{ $tumourHash->{provenence}{"$dir/mutect_$outName.merged.vcf"}{input}}, $tumourBam);

}





sub doAnnovar
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $callModule = $_[2];
	my $tumourHash = $_[3];
	my $refHash = $_[4];

	# somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/$callModule/$module";
	`mkdir -p $dir`;

	my $sgePre = "a$tumourHash->{sge_uniq}";

	my $snvVcf = $tumourHash->{vcf}{$bamModule}{$callModule}{snv}{file};
	my $indelVcf = "none";
	if (($callModule =~ /strelka/) or ($callModule =~ /gatk/))
	{
		$indelVcf = $tumourHash->{vcf}{$bamModule}{$callModule}{indel}{file};
	}

	my $holdJid = "$tumourHash->{vcf}{$bamModule}{$callModule}{snv}{hold_jid}";


	my $outName = $tumourHash->{sample};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=4g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.phoenix.log -o $dir/$outName.phoenix.log \"bash $dir/$outName.cmd\" > $dir/$outName.phoenix.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "module load $module; convert2annovar.pl --format vcf4 --includeinfo $snvVcf > $dir/$outName.snvs.annovar; annotate_variation.pl -geneanno $dir/$outName.snvs.annovar -buildver hg19 \$HUMANDB\n";

		unless ($indelVcf eq "none")
		{
			print COMMAND "convert2annovar.pl --format vcf4 --includeinfo $indelVcf > $dir/$outName.indels.annovar; annotate_variation.pl -geneanno $dir/$outName.indels.annovar -buildver hg19 \$HUMANDB\n";
		}
		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [ANNOVAR] Skipping: $dir/$outName.touch\n";
	}


	$tumourHash->{vcf}{$bamModule}{$callModule}{$module}{snv}{exonic_variant_function}{file} = "$dir/$outName.snvs.annovar.exonic_variant_function";
	$tumourHash->{vcf}{$bamModule}{$callModule}{$module}{snv}{variant_function}{file} = "$dir/$outName.snvs.annovar.variant_function";
	$tumourHash->{vcf}{$bamModule}{$callModule}{$module}{snv}{hold_jid} = "$sgePre$outName";
	
	unless ($indelVcf eq "none")
	{
		$tumourHash->{vcf}{$bamModule}{$callModule}{$module}{indel}{exonic_variant_function}{file} = "$dir/$outName.indels.annovar.exonic_variant_function";
		$tumourHash->{vcf}{$bamModule}{$callModule}{$module}{indel}{variant_function}{file} = "$dir/$outName.indels.annovar.variant_function";
		$tumourHash->{vcf}{$bamModule}{$callModule}{$module}{indel}{hold_jid} = "$sgePre$outName";
	}


	$tumourHash->{provenence}{"$dir/$outName.snvs.annovar"}{command} = "$dir/$outName.cmd";
	push (@{ $tumourHash->{provenence}{"$dir/$outName.snvs.annovar"}{input}}, $snvVcf);

	unless ($indelVcf eq "none")
	{
		$tumourHash->{provenence}{"$dir/$outName.indels.annovar"}{command} = "$dir/$outName.cmd";
		push (@{ $tumourHash->{provenence}{"$dir/$outName.indels.annovar"}{input}}, $indelVcf);
	}

}




sub doFunSeq
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $callModule = $_[2];
	my $annovarModule = $_[3];		# to hold for
	my $tumourHash = $_[4];
	my $refHash = $_[5];

	# somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/$callModule/$module";
	`mkdir -p $dir`;

	my $sgePre = "u$tumourHash->{sge_uniq}";

	my $snvVcf = $tumourHash->{vcf}{$bamModule}{$callModule}{snv}{file};

	my $funseqMode = 1;		# somatic
	if ($callModule =~ /gatk/)
	{
		$funseqMode = 2;		# germline
	}

	my $holdJid = "$tumourHash->{vcf}{$bamModule}{$callModule}{snv}{hold_jid},$tumourHash->{vcf}{$bamModule}{$callModule}{$annovarModule}{snv}{hold_jid}";


	my $outName = $tumourHash->{sample};

	my $funseqOut = `basename $snvVcf`;
	chomp $funseqOut;
	$funseqOut =~ s/(.*?)\..*/$1/;

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		`mkdir $dir/out`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=16g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "module load $module; cd $dir; funseq.sh -f $snvVcf -maf 0 -m $funseqMode -inf vcf -outf vcf\n";

		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [FUNSEQ] Skipping: $dir/$outName.touch\n";
	}


	$tumourHash->{vcf}{$bamModule}{$callModule}{$module}{snv}{file} = "$dir/out/$funseqOut.FunSEQ.vcf";
	$tumourHash->{vcf}{$bamModule}{$callModule}{$module}{snv}{hold_jid} = "$sgePre$outName";

	# temporary to avoid die on provenence disconnect
#	$tumourHash->{provenence}{"$dir/out/$outName.FunSEQ.vcf"}{command} = "$dir/$outName.cmd";
#	push (@{ $tumourHash->{provenence}{"$dir/out/$funseqOut.FunSEQ.vcf"}{input}}, $snvVcf);

}







sub doAnnotation
{
	my $bamModule = $_[0];
	my $callModule = $_[1];
	my $annovarModule = $_[2];
	my $funseqModule = $_[3];
	my $normalHash = $_[4];
	my $tumourHash = $_[5];
	my $refHash = $_[6];

	my $seqType = $tumourHash->{sample_type};

	# somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/$callModule/final";
	`mkdir -p $dir`;

	my $sgePre = "o$tumourHash->{sge_uniq}";

	my $selectFilter = $refHash->{"select-filter-script"};
	my $selectBed = $refHash->{"select-bed-script"};
	my $addAnnovar = $refHash->{"add-annovar-script"};
	my $addFunseq = $refHash->{"add-funseq-script"};
	my $addDBSNP = $refHash->{"add-dbsnp-script"};
	my $addBed = $refHash->{"add-bed-track-script"};
	my $removeMouse = $refHash->{"remove-mouse-script"};
	my $addCosmic = $refHash{'add-cosmic-script'};

	my $perlLib = $refHash{'perl lib path'};

	my $addVerification = $refHash->{"annotate-verification-script"};
	my $verifiedFile = "$pathRoot/../verification/$tumourHash->{sample}.verified.vcf";

	my $addDCC = $refHash->{"add-dcc-metadata"};
	my $pathToTumourJSON = "$tumourHash->{bam}{$bamModule}{final}{json}{file}";
	my $pathToNormalJSON = "$normalHash->{bam}{$bamModule}{final}{json}{file}";

	my $mouseBlacklist = $refHash->{"bwa/0.6.2-mouse-blacklist"};	# assuming bwa
	if ($bamModule =~ /novocraft/)
	{
		$mouseBlacklist = $refHash->{"novocraft/3.00.05-mouse-blacklist"};
	}

	my $exomePath = $refHash->{"exome-target-tabix"};
	my $dbSNPpath = $refHash->{"dbsnp-tabix"};
	my $cosmicPath = $refHash->{"cosmic-tabix"};

	my $trackScript = "";
	for my $track (qw(rmsk segdup simpleRepeat))
	{
		$trackScript .= " | $addBed $refHash->{$track} $track ";
	}
	$trackScript .= " | $addBed $exomePath exomeTarget ";


	my $snvVcf = $tumourHash->{vcf}{$bamModule}{$callModule}{snv}{file};
	my $snvAnnoVF = $tumourHash->{vcf}{$bamModule}{$callModule}{$annovarModule}{snv}{variant_function}{file};
	my $snvAnnoEVF = $tumourHash->{vcf}{$bamModule}{$callModule}{$annovarModule}{snv}{exonic_variant_function}{file};
	my $snvFunseq = $tumourHash->{vcf}{$bamModule}{$callModule}{$funseqModule}{snv}{file};

	my $indelVcf = "none";
	my $indelAnnoVF = "none";
	my $indelAnnoEVF = "none";

	if (($callModule =~ /strelka/) or ($callModule =~ /gatk/))
	{
		$indelVcf = $tumourHash->{vcf}{$bamModule}{$callModule}{indel}{file};
		$indelAnnoVF = $tumourHash->{vcf}{$bamModule}{$callModule}{$annovarModule}{indel}{variant_function}{file};
		$indelAnnoEVF = $tumourHash->{vcf}{$bamModule}{$callModule}{$annovarModule}{indel}{exonic_variant_function}{file};
	}

	my $holdJid = "$tumourHash->{vcf}{$bamModule}{$callModule}{snv}{hold_jid},$tumourHash->{vcf}{$bamModule}{$callModule}{$annovarModule}{snv}{hold_jid},$tumourHash->{vcf}{$bamModule}{$callModule}{$funseqModule}{snv}{hold_jid},$tumourHash->{bam}{$bamModule}{final}{json}{hold_jid},$normalHash->{bam}{$bamModule}{final}{json}{hold_jid}";


	my $outName = $tumourHash->{sample};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=16g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "$perlLib; cat $snvVcf | $addDCC $bamModule $callModule $seqType $pathToNormalJSON $pathToTumourJSON $normalHash->{sample} $tumourHash->{sample} $tumourHash->{sample_group} | $selectFilter PASS | $addAnnovar $snvAnnoVF $snvAnnoEVF | $addFunseq $snvFunseq | $addDBSNP $dbSNPpath | $addCosmic $cosmicPath  $trackScript";		# trackScript includes leading pipe
		
		if ($seqType eq "exome")
		{
			print COMMAND " | $selectBed $exomePath ";
		}

		if ($tumourHash->{tissue_type} eq "X")
		{
			print COMMAND " | $removeMouse $mouseBlacklist $dir/$outName.mouse-snv.vcf ";
		}

		if (-e $verifiedFile)
		{
			print COMMAND " | $addVerification $verifiedFile ";
		}

		print COMMAND " > $dir/$outName.final.vcf\n";


		unless ($indelVcf eq "none")
		{
			print COMMAND "cat $indelVcf | $selectFilter PASS | $addAnnovar $indelAnnoVF $indelAnnoEVF | $addDBSNP $dbSNPpath | $addCosmic $cosmicPath  $trackScript";		# trackScript includes leading pipe
			
			if ($seqType eq "exome")
			{
				print COMMAND " | $selectBed $exomePath ";
			}

			if ($tumourHash->{tissue_type} eq "X")
			{
				print COMMAND " | $removeMouse $mouseBlacklist $dir/$outName.mouse-indel.vcf ";
			}

			if (-e $verifiedFile)
			{
				print COMMAND " | $addVerification $verifiedFile ";
			}

			print COMMAND " | grep -v \"^#\" >> $dir/$outName.final.vcf\n";
		}

		print COMMAND "(grep '^#' $dir/$outName.final.vcf; grep -v '^#' $dir/$outName.final.vcf | sed 's/chrX/23/' | sed 's/chrY/24/' | sed 's/chrM/25/' | sed 's/chr//' | sort -n -k1,1 -k2,2n | sed -e 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' | sed 's/chr25/chrM/') > $dir/randomtempfile.txt; mv $dir/randomtempfile.txt $dir/$outName.final.vcf\n";

		print COMMAND "\necho phoenixPipe/annotation-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [ANNOTATION] Skipping: $dir/$outName.touch\n";
	}


	$tumourHash->{vcf}{$bamModule}{$callModule}{final}{file} = "$dir/$outName.final.vcf";
	$tumourHash->{vcf}{$bamModule}{$callModule}{final}{hold_jid} = "$sgePre$outName";
	


	$tumourHash->{provenence}{"$dir/$outName.final.vcf"}{command} = "$dir/$outName.cmd";
	push (@{ $tumourHash->{provenence}{"$dir/$outName.final.vcf"}{input}}, $snvVcf);

	unless ($indelVcf eq "none")
	{
		$tumourHash->{provenence}{"$dir/$outName.final.vcf"}{command} = "$dir/$outName.cmd";
		push (@{ $tumourHash->{provenence}{"$dir/$outName.final.vcf"}{input}}, $indelVcf);
	}

}



sub doSingleAnnotation
{
	my $bamModule = $_[0];
	my $callModule = $_[1];
	my $annovarModule = $_[2];
	my $funseqModule = $_[3];
	my $normalHash = $_[4];
	my $refHash = $_[5];

	my $seqType = $normalHash->{sample_type};

	my $pathRoot = buildPath($normalHash);
	my $dir = "$pathRoot/$bamModule/$callModule/final";
	`mkdir -p $dir`;

	my $sgePre = "o$normalHash->{sge_uniq}";

	my $selectFilter = $refHash->{"select-filter-script"};
	my $selectBed = $refHash->{"select-bed-script"};
	my $addAnnovar = $refHash->{"add-annovar-script"};
	my $addFunseq = $refHash->{"add-funseq-script"};
	my $addDBSNP = $refHash->{"add-dbsnp-script"};
	my $addBed = $refHash->{"add-bed-track-script"};
	my $removeMouse = $refHash->{"remove-mouse-script"};
	my $addCosmic = $refHash{'add-cosmic-script'};

	my $addVerification = $refHash->{"annotate-verification-script"};
	my $verifiedFile = "$pathRoot/../verification/$normalHash->{sample}.verified.vcf";

	my $addDCC = $refHash->{"add-dcc-metadata"};
	my $pathToNormalJSON = "$normalHash->{bam}{$bamModule}{final}{json}{file}";

	my $mouseBlacklist = $refHash->{"bwa/0.6.2-mouse-blacklist"};	# assuming bwa
	if ($bamModule =~ /novocraft/)
	{
		$mouseBlacklist = $refHash->{"novocraft/3.00.05-mouse-blacklist"};
	}

	my $exomePath = $refHash->{"exome-target-tabix"};
	my $dbSNPpath = $refHash->{"dbsnp-tabix"};
	my $cosmicPath = $refHash->{"cosmic-tabix"};

	my $trackScript = "";
	for my $track (qw(rmsk segdup simpleRepeat))
	{
		$trackScript .= " | $addBed $refHash->{$track} $track ";
	}
	$trackScript .= " | $addBed $exomePath exomeTarget ";


	my $snvVcf = $normalHash->{vcf}{$bamModule}{$callModule}{snv}{file};
	my $snvAnnoVF = $normalHash->{vcf}{$bamModule}{$callModule}{$annovarModule}{snv}{variant_function}{file};
	my $snvAnnoEVF = $normalHash->{vcf}{$bamModule}{$callModule}{$annovarModule}{snv}{exonic_variant_function}{file};
	my $snvFunseq = $normalHash->{vcf}{$bamModule}{$callModule}{$funseqModule}{snv}{file};

	my $indelVcf = "none";
	my $indelAnnoVF = "none";
	my $indelAnnoEVF = "none";

	if (($callModule =~ /strelka/) or ($callModule =~ /gatk/))
	{
		$indelVcf = $normalHash->{vcf}{$bamModule}{$callModule}{indel}{file};
		$indelAnnoVF = $normalHash->{vcf}{$bamModule}{$callModule}{$annovarModule}{indel}{variant_function}{file};
		$indelAnnoEVF = $normalHash->{vcf}{$bamModule}{$callModule}{$annovarModule}{indel}{exonic_variant_function}{file};
	}

	my $holdJid = "$normalHash->{vcf}{$bamModule}{$callModule}{snv}{hold_jid},$normalHash->{vcf}{$bamModule}{$callModule}{$annovarModule}{snv}{hold_jid},$normalHash->{vcf}{$bamModule}{$callModule}{$funseqModule}{snv}{hold_jid},$normalHash->{bam}{$bamModule}{final}{json}{hold_jid}";


	my $outName = $normalHash->{sample};
	my $perlLib = $refHash{'perl lib path'};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=4g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "$perlLib; cat $snvVcf | $selectFilter PASS | $addAnnovar $snvAnnoVF $snvAnnoEVF | $addFunseq $snvFunseq | $addDBSNP $dbSNPpath | $addCosmic $cosmicPath  $trackScript";		# trackScript includes leading pipe
		
		if ($seqType eq "exome")
		{
			print COMMAND " | $selectBed $exomePath ";
		}

		if ($normalHash->{tissue_type} eq "X")
		{
			print COMMAND " | $removeMouse $mouseBlacklist $dir/$outName.mouse-snv.vcf ";
		}

		if (-e $verifiedFile)
		{
			print COMMAND " | $addVerification $verifiedFile ";
		}

		print COMMAND " > $dir/$outName.final.vcf\n";


		unless ($indelVcf eq "none")
		{
			print COMMAND "cat $indelVcf | $selectFilter PASS | $addAnnovar $indelAnnoVF $indelAnnoEVF | $addDBSNP $dbSNPpath | $addCosmic $cosmicPath  $trackScript";		# trackScript includes leading pipe
			
			if ($seqType eq "exome")
			{
				print COMMAND " | $selectBed $exomePath ";
			}

			if ($normalHash->{tissue_type} eq "X")
			{
				print COMMAND " | $removeMouse $mouseBlacklist $dir/$outName.mouse-indel.vcf ";
			}

			if (-e $verifiedFile)
			{
				print COMMAND " | $addVerification $verifiedFile ";
			}

			print COMMAND " | grep -v \"^#\" >> $dir/$outName.final.vcf\n";
		}

		print COMMAND "(grep '^#' $dir/$outName.final.vcf; grep -v '^#' $dir/$outName.final.vcf | sed 's/chrX/23/' | sed 's/chrY/24/' | sed 's/chrM/25/' | sed 's/chr//' | sort -n -k1,1 -k2,2n | sed -e 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' | sed 's/chr25/chrM/') > $dir/randomtempfile.txt; mv $dir/randomtempfile.txt $dir/$outName.final.vcf\n";

		print COMMAND "\necho phoenixPipe/annotation-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [ANNOTATION] Skipping: $dir/$outName.touch\n";
	}


	$normalHash->{vcf}{$bamModule}{$callModule}{final}{file} = "$dir/$outName.final.vcf";
	$normalHash->{vcf}{$bamModule}{$callModule}{final}{hold_jid} = "$sgePre$outName";
	


	$normalHash->{provenence}{"$dir/$outName.final.vcf"}{command} = "$dir/$outName.cmd";
	push (@{ $normalHash->{provenence}{"$dir/$outName.final.vcf"}{input}}, $snvVcf);

	unless ($indelVcf eq "none")
	{
		$normalHash->{provenence}{"$dir/$outName.final.vcf"}{command} = "$dir/$outName.cmd";
		push (@{ $normalHash->{provenence}{"$dir/$outName.final.vcf"}{input}}, $indelVcf);
	}

}



sub doIntersect
{
	my $bamModule = $_[0];
	my $callModule1 = $_[1];
	my $callModule2 = $_[2];
	my $dirName = $_[3];
	my $tumourHash = $_[4];
	my $refHash = $_[5];

	my $seqType = $tumourHash->{sample_type};

	# somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);

	my $dir = "$pathRoot/$dirName";
	`mkdir -p $dir`;

	my $sgePre = "i$tumourHash->{sge_uniq}";

	my $intersectScript = $refHash->{"intersect-vcf-script"};

	my $alignString = "$bamModule";
	my $callString = "$callModule1,$callModule2";

	my $reportScript = $refHash->{"somatic vcf report"};

	my $file1 = $tumourHash->{vcf}{$bamModule}{$callModule1}{final}{file};
	my $file2 = $tumourHash->{vcf}{$bamModule}{$callModule2}{final}{file};

	my $holdJid = "$tumourHash->{vcf}{$bamModule}{$callModule1}{final}{hold_jid},$tumourHash->{vcf}{$bamModule}{$callModule2}{final}{hold_jid}";


	my $outName = $tumourHash->{sample};

	my $perlLib = $refHash{'perl lib path'};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=2g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;		#rmsk segdup simpleRepeat
#		print COMMAND "$perlLib; $intersectScript $file1 $file2 $alignString $callString | grep -v 'TRACK=rmsk' | grep -v 'TRACK=segdup' | grep -v 'TRACK=simpleRepeat' | grep -v 'DBSNP_G5' > $dir/$outName.final.vcf\n";
		print COMMAND "$perlLib; $intersectScript $file1 $file2 $alignString $callString > $dir/$outName.final.vcf\n";
		print COMMAND "module load R/3.3.0; $reportScript $dir/$outName.final.vcf $seqType\n";
		print COMMAND "\necho phoenixPipe/intersect-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [INTERSECT] Skipping: $dir/$outName.touch\n";
	}


	$tumourHash->{vcf}{final}{file} = "$dir/$outName.final.vcf";
	$tumourHash->{vcf}{final}{hold_jid} = "$sgePre$outName";

	$tumourHash->{provenence}{"$dir/$outName.final.vcf"}{command} = "$dir/$outName.cmd";
	push (@{ $tumourHash->{provenence}{"$dir/$outName.final.vcf"}{input}}, $file1);
	push (@{ $tumourHash->{provenence}{"$dir/$outName.final.vcf"}{input}}, $file2);
}




sub doIntegration
{
	my $bamModule = $_[0];
	my $tumourHash = $_[1];
	my $refHash = $_[2];

	# somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);

	my $dir = "$pathRoot/$bamModule/integration";
	`mkdir -p $dir`;

	my $sgePre = "z$tumourHash->{sge_uniq}";

	my $onePageScript = $refHash->{"one page script"};

	# sv file doesnt exist yet (waiting for the merge) so this'll be a warning for now
	my $holdJid = "$tumourHash->{vcf}{final}{hold_jid},$tumourHash->{vcf}{germline}{hold_jid},$tumourHash->{sv}{final}{hold_jid},$tumourHash->{cnv}{final}{hold_jid}";

	my $outName = $tumourHash->{sample};
	my $perlLib = $refHash{'perl lib path'};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=16g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "cd $dir; $perlLib; module load R/3.3.0; $onePageScript $pathRoot\n";
		print COMMAND "\necho phoenixPipe/integration-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [INTEGRATION] Skipping: $dir/$outName.touch\n";
	}

	# don't need these yet
#	$tumourHash->{vcf}{final}{file} = "$dir/$outName.final.vcf";
#	$tumourHash->{vcf}{final}{hold_jid} = "$sgePre$outName";

#	$tumourHash->{provenence}{"$dir/$outName.final.vcf"}{command} = "$dir/$outName.cmd";
#	push (@{ $tumourHash->{provenence}{"$dir/$outName.final.vcf"}{input}}, $file1);
#	push (@{ $tumourHash->{provenence}{"$dir/$outName.final.vcf"}{input}}, $file2);
}







sub doFilterGATK
{
	my $bamModule = $_[0];
	my $callModule = $_[1];
	my $tumourHash = $_[2];
	my $refHash = $_[3];

	my $gatkQualCut = 50;

	my $seqType = $tumourHash->{sample_type};

	# somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);

	my $dir = "$pathRoot/$bamModule/final_gatk-germline";
	`mkdir -p $dir`;

	my $sgePre = "i$tumourHash->{sge_uniq}";

	my $file = $tumourHash->{vcf}{$bamModule}{$callModule}{final}{file};

	my $holdJid = "$tumourHash->{vcf}{$bamModule}{$callModule}{final}{hold_jid}";

	my $germlineScript = $refHash->{"gatk germline filter"};
	my $pathoGL = $refHash->{"germline pathogen script"};

	my $outName = "$tumourHash->{sample}";

	my $perlLib = $refHash{'perl lib path'};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=16g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;		#rmsk segdup simpleRepeat
#		print COMMAND "$perlLib; (grep \"^#\" $file; grep -v \"^#\" $file | awk \'{if (\$6 >= $gatkQualCut) print \$0}\') | grep -v 'TRACK=rmsk' $file | grep -v 'TRACK=segdup' | grep -v 'TRACK=simpleRepeat' > $dir/$outName.final.vcf\n";
		print COMMAND "$perlLib; (grep \"^#\" $file; grep -v \"^#\" $file | awk \'{if (\$6 >= $gatkQualCut) print \$0}\') | $germlineScript > $dir/$outName.germline.final.vcf\n";
#		print COMMAND "cd $dir; bash $pathoGL $dir/$outName.germline.final.vcf\n";
		print COMMAND "\necho phoenixPipe/gatkFilter-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [GATK_FILTER] Skipping: $dir/$outName.touch\n";
	}


	$tumourHash->{vcf}{germline}{file} = "$dir/$outName.germline.final.vcf";
	$tumourHash->{vcf}{germline}{hold_jid} = "$sgePre$outName";

	$tumourHash->{provenence}{"$dir/$outName.germline.final.vcf"}{command} = "$dir/$outName.cmd";
	push (@{ $tumourHash->{provenence}{"$dir/$outName.germline.final.vcf"}{input}}, $file);
}





sub doGATK
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $normalHash = $_[2];
	my $tumourHash = $_[3];
	my $refHash = $_[4];

	my $seqType = $tumourHash->{sample_type};

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/$module";
	`mkdir -p $dir`;

	my $sgePre = "g$tumourHash->{sge_uniq}";

	my $normalBam = $normalHash->{bam}{$bamModule}{final}{file};
	my $tumourBam = $tumourHash->{bam}{$bamModule}{final}{file};
	my $holdJid = "$normalHash->{bam}{$bamModule}{final}{hold_jid},$tumourHash->{bam}{$bamModule}{final}{hold_jid}";

	my $gatkHelper;
	my $gatkWaiter;
	if (exists $refHash->{"$module helper"})
	{
		$gatkHelper = $refHash->{"$module helper"};
		$gatkWaiter = $refHash->{"$module waiter"};
	}
	else
	{
		die " [GATK] Aborting: no helper script found for $module\n\n";
	}

	my $outName = $tumourHash->{sample};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "cd $dir; mkdir bam; ln -s $normalBam bam/.; ln -s $tumourBam bam/.; $gatkHelper -b bam -t $seqType -n $outName > gatk; bash gatk\n";
		print COMMAND "$gatkWaiter\n";
		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [GATK] Skipping: $dir/$outName.touch\n";
	}

#	depricated
#	$normalHash->{bam}{"$bamModule/$module"}{realnrecal}{files_pre} = "$dir/bam_files/$outName.bam.realigned.recal.chr";
#	$normalHash->{bam}{"$bamModule/$module"}{realnrecal}{hold_jid} = "$sgePre$outName,gatk_$outName\\*";
#	$tumourHash->{bam}{"$bamModule/$module"}{realnrecal}{files_pre} = "$dir/bam_files/$outName.bam.realigned.recal.chr";
#	$tumourHash->{bam}{"$bamModule/$module"}{realnrecal}{hold_jid} = "$sgePre$outName,gatk_$outName\\*";

	$tumourHash->{vcf}{$bamModule}{$module}{snv}{file} = "$dir/$outName.bam.realigned.recal.bam.snps.raw.filtered.annotated.vcf";
	$tumourHash->{vcf}{$bamModule}{$module}{snv}{hold_jid} = "$sgePre$outName";
	$tumourHash->{vcf}{$bamModule}{$module}{indel}{file} = "$dir/$outName.bam.realigned.recal.bam.indels.raw.filtered.annotated.vcf";
	$tumourHash->{vcf}{$bamModule}{$module}{indel}{hold_jid} = "$sgePre$outName";

	$tumourHash->{provenence}{"$dir/$outName.bam.realigned.recal.bam.snps.raw.filtered.annotated.vcf"}{command} = "$dir/$outName.cmd";
	$tumourHash->{provenence}{"$dir/$outName.bam.realigned.recal.bam.indels.raw.filtered.annotated.vcf"}{command} = "$dir/$outName.cmd";

	push (@{ $tumourHash->{provenence}{"$dir/$outName.bam.realigned.recal.bam.snps.raw.filtered.annotated.vcf"}{input}}, $normalBam);
	push (@{ $tumourHash->{provenence}{"$dir/$outName.bam.realigned.recal.bam.snps.raw.filtered.annotated.vcf"}{input}}, $tumourBam);
	push (@{ $tumourHash->{provenence}{"$dir/$outName.bam.realigned.recal.bam.indels.raw.filtered.annotated.vcf"}{input}}, $normalBam);
	push (@{ $tumourHash->{provenence}{"$dir/$outName.bam.realigned.recal.bam.indels.raw.filtered.annotated.vcf"}{input}}, $tumourBam);

}





sub doSingleGATK
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $normalHash = $_[2];
	my $refHash = $_[3];

	my $seqType = $normalHash->{sample_type};

	my $pathRoot = buildPath($normalHash);
	my $dir = "$pathRoot/$bamModule/$module";
	`mkdir -p $dir`;

	my $sgePre = "g$normalHash->{sge_uniq}";

	my $normalBam = $normalHash->{bam}{$bamModule}{final}{file};
	my $holdJid = "$normalHash->{bam}{$bamModule}{final}{hold_jid}";

	my $gatkHelper;
	my $gatkWaiter;
	if (exists $refHash->{"$module helper"})
	{
		$gatkHelper = $refHash->{"$module helper"};
		$gatkWaiter = $refHash->{"$module waiter"};
	}
	else
	{
		die " [GATK] Aborting: no helper script found for $module\n\n";
	}

	my $outName = $normalHash->{sample};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "cd $dir; mkdir bam; ln -s $normalBam bam/.; $gatkHelper -b bam -t $seqType -n $outName > gatk; bash gatk\n";
		print COMMAND "$gatkWaiter\n";
		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [GATK] Skipping: $dir/$outName.touch\n";
	}

	$normalHash->{vcf}{$bamModule}{$module}{snv}{file} = "$dir/$outName.bam.realigned.recal.bam.snps.raw.filtered.annotated.vcf";
	$normalHash->{vcf}{$bamModule}{$module}{snv}{hold_jid} = "$sgePre$outName";
	$normalHash->{vcf}{$bamModule}{$module}{indel}{file} = "$dir/$outName.bam.realigned.recal.bam.indels.raw.filtered.annotated.vcf";
	$normalHash->{vcf}{$bamModule}{$module}{indel}{hold_jid} = "$sgePre$outName";

	$normalHash->{provenence}{"$dir/$outName.bam.realigned.recal.bam.snps.raw.filtered.annotated.vcf"}{command} = "$dir/$outName.cmd";
	$normalHash->{provenence}{"$dir/$outName.bam.realigned.recal.bam.indels.raw.filtered.annotated.vcf"}{command} = "$dir/$outName.cmd";

	push (@{ $normalHash->{provenence}{"$dir/$outName.bam.realigned.recal.bam.snps.raw.filtered.annotated.vcf"}{input}}, $normalBam);
	push (@{ $normalHash->{provenence}{"$dir/$outName.bam.realigned.recal.bam.indels.raw.filtered.annotated.vcf"}{input}}, $normalBam);

}













# depricated!
sub doTimSomatic
{
	my $module = $_[0];
	my $normalHash = $_[1];
	my $tumourHash = $_[2];

	my $seqType = $tumourHash->{sample_type};

	# somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$module";

	my $sgePre = "t$tumourHash->{sge_uniq}";

	my $holdJid = $tumourHash->{vcf}{$module}{gatkout}{hold_jid};

	my $helper = "/u/rdenroche/svn/seq_prod_bio/postProcessing/somaticClassifier/trunk/doPairwiseClassify";


	my $outName = $tumourHash->{sample};

	if (-e $dir)
	{
		unless (-e "$dir/somatic.touch")
		{
			`touch $dir/somatic.touch`;

			open (SUBFILE, ">$dir/TSsub_$outName") or die;
			print SUBFILE "qsub -cwd -b y -l h_vmem=2g -N $sgePre$outName -hold_jid $holdJid -e $dir/TS-$outName.log -o $dir/TS-$outName.log \"bash $dir/TScommand_$outName\"\n";
			close SUBFILE;
	
			open (COMMAND, ">$dir/TScommand_$outName") or die;
			print COMMAND "cd $dir; $helper $outName\n";
			close COMMAND;

			`bash $dir/TSsub_$outName`;
		}
		else
		{
			warn " [TIMSOMATIC] Skipping analysis: $dir/somatic.touch\n";
		}
	}
	else
	{
		die " [TIMSOMATIC] Aborting! Input does not exist: $dir\n";
	}


}


# depricated!
# could do the merge and sample split separately
sub doSplitGATK
{
	my $module = $_[0];
	my $normalHash = $_[1];
	my $tumourHash = $_[2];

	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$module/split";
	`mkdir -p $dir`;

	my $touchFile = "splitgatk.touch";

	my $sgePre = "p$tumourHash->{sge_uniq}";

	my $realignedPrefix = $tumourHash->{bam}{$module}{realnrecal}{files_pre};
	my $holdJid = $tumourHash->{bam}{$module}{realnrecal}{hold_jid};

	my $fileList;
	my @files;

	my $inputList = "";


	# create readgroup lists for each sample
	open (NORMALRG, ">$dir/$normalHash->{sample}.rgList.txt") or die;
	for my $ius (keys %{ $normalHash->{fastq} })
	{
		print NORMALRG "$normalHash->{fastq}{$ius}{run}_$normalHash->{fastq}{$ius}{lane}_$normalHash->{fastq}{$ius}{barcode}\n";
	}
	close NORMALRG;

	open (TUMOURRG, ">$dir/$tumourHash->{sample}.rgList.txt") or die;
	for my $ius (keys %{ $tumourHash->{fastq} })
	{
		print TUMOURRG "$tumourHash->{fastq}{$ius}{run}_$tumourHash->{fastq}{$ius}{lane}_$tumourHash->{fastq}{$ius}{barcode}\n";
	}
	close TUMOURRG;

	my $outName;

	if ($fileList = `ls $realignedPrefix*bam`)
	{
		unless (-e "$dir/$touchFile")
		{
			`touch $dir/$touchFile`;
			`touch $dir/$normalHash->{sample}.realn.recal.bam`;
			`touch $dir/$tumourHash->{sample}.realn.recal.bam`;

			open (SUBFILE, ">$dir/sub_splitGATK") or die;

			print SUBFILE "qsub -cwd -b y -l h_vmem=16g -N ${sgePre}$tumourHash->{sample}-mergeGATK -hold_jid $holdJid -e $dir/mergeGATK.log -o $dir/mergeGATK.log \"bash $dir/command_mergeGATK\"\n";

			$outName = $normalHash->{sample};
			print SUBFILE "qsub -cwd -b y -l h_vmem=8g -N $sgePre$outName -hold_jid $holdJid,${sgePre}$tumourHash->{sample}-mergeGATK -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/command_$outName\"\n";

			$outName = $tumourHash->{sample};
			print SUBFILE "qsub -cwd -b y -l h_vmem=8g -N $sgePre$outName -hold_jid $holdJid,${sgePre}$tumourHash->{sample}-mergeGATK -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/command_$outName\"\n";

			close SUBFILE;

			chomp $fileList;
			@files = split(/\n/, $fileList);

			# create list of chr files
			for my $file (@files)
			{
				$inputList = "$inputList I=$file";
			}

			open (COMMAND, ">$dir/command_mergeGATK") or die;
			print COMMAND "module load picard; java -Xmx6g -jar \$PICARDROOT/MergeSamFiles.jar $inputList O=$dir/mergedGATK.realigned.recal.bam CREATE_INDEX=true ASSUME_SORTED=true SO=coordinate TMP_DIR=$dir/picardTmp\n";
			close COMMAND;


			$outName = $normalHash->{sample};
			open (COMMAND, ">$dir/command_$outName") or die;
			print COMMAND "module load samtools; samtools view -hbR $dir/$outName.rgList.txt $dir/mergedGATK.realigned.recal.bam > $dir/$outName.realn.recal.bam; samtools index $dir/$outName.realn.recal.bam\n";
			close COMMAND;

			$outName = $tumourHash->{sample};
			open (COMMAND, ">$dir/command_$outName") or die;
			print COMMAND "module load samtools; samtools view -hbR $dir/$outName.rgList.txt $dir/mergedGATK.realigned.recal.bam > $dir/$outName.realn.recal.bam; samtools index $dir/$outName.realn.recal.bam\n";
			close COMMAND;

			`bash $dir/sub_splitGATK`;
		}
		else
		{
			warn " [SPLITGATK] Skipping analysis: $dir/$touchFile\n";
		}
	}
	else
	{
		die " [SPLITGATK] Aborting! Input does not exist: $realignedPrefix*\n";
	}

	$outName = $normalHash->{sample};
	$normalHash->{bam}{$module}{final}{file} = "$dir/$outName.realn.recal.bam";
	$normalHash->{bam}{$module}{final}{hold_jid} = "$sgePre$outName";

	$outName = $tumourHash->{sample};
	$tumourHash->{bam}{$module}{final}{file} = "$dir/$outName.realn.recal.bam";
	$tumourHash->{bam}{$module}{final}{hold_jid} = "$sgePre$outName";

}

# depricated!
sub doNewGATK
{

	my $module = $_[0];
	my $bamModule = $_[1];
	my $normalHash = $_[2];
	my $tumourHash = $_[3];

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/$module";
	`mkdir -p $dir`;

	my $sgePre = "g$tumourHash->{sge_uniq}";

	my $normalBam = $normalHash->{bam}{$bamModule}{final}{file};
	my $tumourBam = $tumourHash->{bam}{$bamModule}{final}{file};
	my $holdJid = "$normalHash->{bam}{$bamModule}{final}{hold_jid},$tumourHash->{bam}{$bamModule}{final}{hold_jid}";

	my $gatkHelper;
	if ($module eq "gatk/2.4.9")
	{
		$gatkHelper = "/u/rdenroche/GATKv4_realign_recal.pl";
	}
	else
	{
		die " [GATK] Aborting! Not a supported module.\n";
	}


	my $outName = "g249";

	if ((-e $normalBam) and (-e $tumourBam))
	{
		unless (-e "$dir/gatk.touch")
		{
			`touch $dir/gatk.touch`;
			`touch $dir/$normalHash->{sample}.realigned.recal.bam`;
			`touch $dir/$tumourHash->{sample}.realigned.recal.bam`;

			open (SUBFILE, ">$dir/sub_$outName") or die;
			print SUBFILE "qsub -cwd -b y -l h_vmem=2g -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/command_$outName\"\n";
			close SUBFILE;
	
			open (COMMAND, ">$dir/command_$outName") or die;
			print COMMAND "cd $dir; $gatkHelper $outName $normalBam $tumourBam; sleep 300\n";
			close COMMAND;

			`bash $dir/sub_$outName`;
		}
		else
		{
			warn " [GATK] Skipping analysis: $dir/gatk.touch\n";
		}
	}
	else
	{
		die " [GATK] Aborting! Input does not exist: $normalBam or $tumourBam\n";
	}

	$normalHash->{bam}{"$bamModule/$module"}{final}{file} = "$dir/$normalHash->{sample}.realigned.recal.bam";
	$normalHash->{bam}{"$bamModule/$module"}{final}{hold_jid} = "$sgePre$outName,pqr_$outName-$normalHash->{sample}.realigned.bam";

	$tumourHash->{bam}{"$bamModule/$module"}{final}{file} = "$dir/$tumourHash->{sample}.realigned.recal.bam";
	$tumourHash->{bam}{"$bamModule/$module"}{final}{hold_jid} = "$sgePre$outName,pqr_$outName-$tumourHash->{sample}.realigned.bam";


}



sub doRealignGATK
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $target = $_[2];
	my $normalHash = $_[3];
	my $tumourHash = $_[4];
	my $refHash = $_[5];


	my $pathRoot = buildPath($normalHash);
	my $sgePre = "r$normalHash->{sge_uniq}";
	my $normalBam = $normalHash->{bam}{$bamModule}{final}{file};
	my $tumourBam;
	my $bamInput = "-I $normalBam";
	my $holdJid = "$normalHash->{bam}{$bamModule}{final}{hold_jid}";
	my $outName = $normalHash->{sample};

	my $targetParam = "";
	if ($target eq "exome")
	{
		$targetParam = " -L " . $refHash->{"exome-target"} . " ";
	}

	unless ($tumourHash eq "single")
	{
		$pathRoot = buildPath($tumourHash);
		$sgePre = "r$tumourHash->{sge_uniq}";
		$tumourBam = $tumourHash->{bam}{$bamModule}{final}{file};
		$bamInput .= " -I $tumourBam";
		$holdJid .= ",$tumourHash->{bam}{$bamModule}{final}{hold_jid}";
		$outName = $tumourHash->{sample};
	}

	my $dir = "$pathRoot/$bamModule/$module/realign";
	`mkdir -p $dir`;

	my $fastaRef = $refHash->{"hg19_random"};
	my $dbsnpRef = $refHash->{"dbSNP"};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
#		print SUBFILE "qsub -cwd -b y -l h_vmem=24g,exclusive=1 -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		print SUBFILE "qsub -cwd -b y -l h_vmem=4g -pe smp 6 -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "cd $dir\n";		# so that -nWayOut works right
#		print COMMAND "module load $module; java -Xmx16g -Djava.io.tmpdir=$dir/$outName.tmp -jar \$GATKROOT/GenomeAnalysisTK.jar -T RealignerTargetCreator $bamInput -R $fastaRef --known $dbsnpRef $targetParam -o $dir/$outName.realignerTarget.intervals -nt `grep processor /proc/cpuinfo | wc -l`\n";
		print COMMAND "module load $module; java -Xmx16g -Djava.io.tmpdir=$dir/$outName.tmp -jar \$GATKROOT/GenomeAnalysisTK.jar -T RealignerTargetCreator $bamInput -R $fastaRef --known $dbsnpRef $targetParam -o $dir/$outName.realignerTarget.intervals -nt 8\n";

		print COMMAND "java -Xmx16g -Djava.io.tmpdir=$dir/$outName.tmp -jar \$GATKROOT/GenomeAnalysisTK.jar -T IndelRealigner -targetIntervals $dir/$outName.realignerTarget.intervals $bamInput -R $fastaRef -known $dbsnpRef $targetParam -nWayOut .realigned.bam\n";
		print COMMAND "\necho phoenixPipe/$module/realign-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [GATK-REALIGN] Skipping: $dir/$outName.touch\n";
	}


	$normalHash->{bam}{"$bamModule/$module"}{realigned}{file} = "$dir/$normalHash->{sample}.realigned.bam";
	$normalHash->{bam}{"$bamModule/$module"}{realigned}{hold_jid} = "$sgePre$outName";

	$normalHash->{provenence}{"$dir/$normalHash->{sample}.realigned.bam"}{command} = "$dir/$outName.cmd";
	push (@{ $normalHash->{provenence}{"$dir/$normalHash->{sample}.realigned.bam"}{input}}, $normalBam);
	
	unless ($tumourHash eq "single")
	{
		$tumourHash->{bam}{"$bamModule/$module"}{realigned}{file} = "$dir/$tumourHash->{sample}.realigned.bam";
		$tumourHash->{bam}{"$bamModule/$module"}{realigned}{hold_jid} = "$sgePre$outName";

		$tumourHash->{provenence}{"$dir/$tumourHash->{sample}.realigned.bam"}{command} = "$dir/$outName.cmd";
		push (@{ $tumourHash->{provenence}{"$dir/$tumourHash->{sample}.realigned.bam"}{input}}, $tumourBam);
	}
}




sub doRecalGATK
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $target = $_[2];
	my $normalHash = $_[3];
	my $tumourHash = $_[4];
	my $refHash = $_[5];

	my $pathRoot = buildPath($normalHash);
	my $sgePre = "q$normalHash->{sge_uniq}";
	my $normalBam = $normalHash->{bam}{"$bamModule/$module"}{realigned}{file};
	my $tumourBam;
	my $holdJid = "$normalHash->{bam}{\"$bamModule/$module\"}{realigned}{hold_jid}";

	my $targetParam = "";
	if ($target eq "exome")
	{
		$targetParam = " -L " . $refHash->{"exome-target"} . " ";
	}

	unless ($tumourHash eq "single")
	{
		$pathRoot = buildPath($tumourHash);
		$sgePre = "q$tumourHash->{sge_uniq}";
		$tumourBam = $tumourHash->{bam}{"$bamModule/$module"}{realigned}{file};
		$holdJid .= ",$normalHash->{bam}{\"$bamModule/$module\"}{realigned}{hold_jid}";
	}

	my $dir = "$pathRoot/$bamModule/$module/recal";
	`mkdir -p $dir`;

	my $fastaRef = $refHash->{"hg19_random"};
	my $dbsnpRef = $refHash->{"dbSNP"};


	my $outName = $normalHash->{sample};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=24g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "module load $module; java -Xmx16g -Djava.io.tmpdir=$dir/$outName.tmp -jar \$GATKROOT/GenomeAnalysisTK.jar -T BaseRecalibrator -I $normalBam -R $fastaRef -knownSites $dbsnpRef $targetParam -o $dir/$outName.recalibrationData.grp\n";

		print COMMAND "module load $module; java -Xmx16g -Djava.io.tmpdir=$dir/$outName.tmp -jar \$GATKROOT/GenomeAnalysisTK.jar -T PrintReads -I $normalBam -R $fastaRef -BQSR $dir/$outName.recalibrationData.grp $targetParam -o $dir/$outName.realn.recal.bam\n";
		print COMMAND "\necho phoenixPipe/$module/recal-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [GATK-RECAL] Skipping: $dir/$outName.touch\n";
	}


	$normalHash->{bam}{"$bamModule/$module"}{final}{file} = "$dir/$outName.realn.recal.bam";
	$normalHash->{bam}{"$bamModule/$module"}{final}{hold_jid} = "$sgePre$outName";
	# needs cleaning
	
	$normalHash->{provenence}{"$dir/$outName.realn.recal.bam"}{command} = "$dir/$outName.cmd";
	push (@{ $normalHash->{provenence}{"$dir/$outName.realn.recal.bam"}{input}}, $normalBam);


	unless ($tumourHash eq "single")
	{
		$outName = $tumourHash->{sample};
	
		unless (-e "$dir/$outName.touch")
		{
			`touch $dir/$outName.touch`;
			open (SUBFILE, ">$dir/$outName.sub") or die;
			print SUBFILE "qsub -cwd -b y -l h_vmem=24g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
			close SUBFILE;
	
			open (COMMAND, ">$dir/$outName.cmd") or die;
			print COMMAND "module load $module; java -Xmx20g -Djava.io.tmpdir=$dir/$outName.tmp -jar \$GATKROOT/GenomeAnalysisTK.jar -T BaseRecalibrator -I $tumourBam -R $fastaRef -knownSites $dbsnpRef -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o $dir/$outName.recalibrationData.grp\n";
	
			print COMMAND "module load $module; java -Xmx20g -Djava.io.tmpdir=$dir/$outName.tmp -jar \$GATKROOT/GenomeAnalysisTK.jar -T PrintReads -I $tumourBam -R $fastaRef -BQSR $dir/$outName.recalibrationData.grp -o $dir/$outName.realn.recal.bam\n";
			print COMMAND "\necho phoenixPipe/$module/recal-done\n";
			close COMMAND;
	
			`bash $dir/$outName.sub`;
		}
		else
		{
			warn " [GATK-RECAL] Skipping: $dir/$outName.touch\n";
		}
	
	
		$tumourHash->{bam}{"$bamModule/$module"}{final}{file} = "$dir/$outName.realn.recal.bam";
		$tumourHash->{bam}{"$bamModule/$module"}{final}{hold_jid} = "$sgePre$outName";
		# needs cleaning
	
		$tumourHash->{provenence}{"$dir/$outName.realn.recal.bam"}{command} = "$dir/$outName.cmd";
		push (@{ $tumourHash->{provenence}{"$dir/$outName.realn.recal.bam"}{input}}, $tumourBam);
	}
}




sub doHaplotypeGATK
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $target = $_[2];
	my $normalHash = $_[3];
	my $refHash = $_[4];

	my $pathRoot = buildPath($normalHash);
	my $sgePre = "h$normalHash->{sge_uniq}";
	my $normalBam = $normalHash->{bam}{"$bamModule/$module"}{final}{file};
	my $tumourBam;
	my $holdJid = "$normalHash->{bam}{\"$bamModule/$module\"}{final}{hold_jid}";

	my $dir = "$pathRoot/$bamModule/$module/";
	`mkdir -p $dir`;

	my $fastaRef = $refHash->{"hg19_random"};
	my $dbsnpRef = $refHash->{"dbSNP"};

	my $targetParam = "";
	if ($target eq "exome")
	{
		$targetParam = " -L " . $refHash->{"exome-target"} . " ";
	}

	my $outName = $normalHash->{sample};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
#		print SUBFILE "qsub -cwd -b y -l h_vmem=24g,exclusive=1 -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		print SUBFILE "qsub -cwd -b y -l h_vmem=4g -pe smp 6 -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
#		print SUBFILE "qsub -cwd -b y -l h_vmem=8g -pe smp 6 -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
#		print COMMAND "module load $module; java -Xmx10g -Djava.io.tmpdir=$dir/$outName.tmp -jar \$GATKROOT/GenomeAnalysisTK.jar -T HaplotypeCaller -I $normalBam -R $fastaRef --dbsnp $dbsnpRef $targetParam --emitRefConfidence GVCF -nct `grep processor /proc/cpuinfo | wc -l` -o $dir/$outName.raw.g.vcf\n";
		print COMMAND "module load $module; java -Xmx16g -Djava.io.tmpdir=$dir/$outName.tmp -jar \$GATKROOT/GenomeAnalysisTK.jar -T HaplotypeCaller -I $normalBam -R $fastaRef --dbsnp $dbsnpRef $targetParam --emitRefConfidence GVCF -nct 6 -o $dir/$outName.raw.g.vcf\n";
#		print COMMAND "module load $module; java -Xmx32g -Djava.io.tmpdir=$dir/$outName.tmp -jar \$GATKROOT/GenomeAnalysisTK.jar -T HaplotypeCaller -I $normalBam -R $fastaRef --dbsnp $dbsnpRef $targetParam --emitRefConfidence GVCF -nct 6 -o $dir/$outName.raw.g.vcf\n";
		print COMMAND "\necho phoenixPipe/$module/haplo-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [GATK-HAPLO] Skipping: $dir/$outName.touch\n";
	}


	$normalHash->{vcf}{"$bamModule/$module"}{raw_gvcf}{file} = "$dir/$outName.realn.recal.bam";
	$normalHash->{vcf}{"$bamModule/$module"}{raw_gvcf}{hold_jid} = "$sgePre$outName";
	
	$normalHash->{provenence}{"$dir/$outName.raw.g.vcf"}{command} = "$dir/$outName.cmd";
	push (@{ $normalHash->{provenence}{"$dir/$outName.raw.g.vcf"}{input}}, $normalBam);


}





# depricated!
sub doVarScan2
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $normalHash = $_[2];
	my $tumourHash = $_[3];
	
	my $seqType = $tumourHash->{sample_type};

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/$module";
	`mkdir -p $dir`;

	my $sgePre = "v$tumourHash->{sge_uniq}";

	my $normalBam = $normalHash->{bam}{$bamModule}{final}{file};
	my $tumourBam = $tumourHash->{bam}{$bamModule}{final}{file};
	my $holdJid = "$normalHash->{bam}{$bamModule}{final}{hold_jid},$tumourHash->{bam}{$bamModule}{final}{hold_jid}";

	my $samtools = "samtools/0.1.19";
	my $fastaRef = "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/samtools/hg19_random.fa";

	my $varscan = "/u/rdenroche/varscan/VarScan.v2.3.2.jar";

	my $outName = $tumourHash->{sample};

	if ((-e $normalBam) and (-e $tumourBam))
	{
		unless (-e "$dir/varscan.touch")
		{
			`touch $dir/varscan.touch`;
			open (SUBFILE, ">$dir/sub_$outName") or die;
			print SUBFILE "qsub -cwd -b y -l h_vmem=16g -N $sgePre$outName-pileupN -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/command_pileup-$normalHash->{sample}\"\n";
			print SUBFILE "qsub -cwd -b y -l h_vmem=16g -N $sgePre$outName-pileupT -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/command_pileup-$tumourHash->{sample}\"\n";
			print SUBFILE "qsub -cwd -b y -l h_vmem=16g -N $sgePre$outName -hold_jid $holdJid,$sgePre$outName-pileupN,$sgePre$outName-pileupT -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/command_$outName\"\n";
			close SUBFILE;
	
			open (COMMAND, ">$dir/command_pileup-$normalHash->{sample}") or die;
			print COMMAND "module load $samtools; samtools mpileup -f $fastaRef -q 10 $normalBam > $dir/$normalHash->{sample}.mpileup\n";
			close COMMAND;

			open (COMMAND, ">$dir/command_pileup-$tumourHash->{sample}") or die;
			print COMMAND "module load $samtools; samtools mpileup -f $fastaRef -q 10 $tumourBam > $dir/$tumourHash->{sample}.mpileup\n";
			close COMMAND;

			open (COMMAND, ">$dir/command_$outName") or die;
			print COMMAND "java -Xmx10g -jar $varscan somatic $dir/$normalHash->{sample}.mpileup $dir/$tumourHash->{sample}.mpileup $dir/$outName --min-coverage 4 --min-var-freq 0.05 --p-value 0.05 --somatic-p-value 0.05 --output-vcf 1\n";
			close COMMAND;

			`bash $dir/sub_$outName`;
		}
		else
		{
			warn " [VARSCAN2] Skipping analysis: $dir/varscan.touch\n";
		}
	}
	else
	{
		die " [VARSCAN2] Aborting! Input does not exist: $normalBam or $tumourBam\n";
	}

}


sub doHMMcopy
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $normalHash = $_[2];
	my $tumourHash = $_[3];
	my $refHash = $_[4];

	my $seqType = $tumourHash->{sample_type};

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/$module";
	`mkdir -p $dir`;

	my $sgePre = "h$tumourHash->{sge_uniq}";

	my $normalBam = $normalHash->{bam}{$bamModule}{final}{file};
	my $tumourBam = $tumourHash->{bam}{$bamModule}{final}{file};
	my $holdJid = "$normalHash->{bam}{$bamModule}{final}{hold_jid},$tumourHash->{bam}{$bamModule}{final}{hold_jid}";

	my $Rmodule = $refHash->{"R module"};
	my $HMMscript = $refHash->{"HMM R script"};

	my $outName = $tumourHash->{sample};

	my $tumourName = $tumourHash->{sample};
	my $normalName = $normalHash->{sample};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=16g -q $refHash->{sge_queue} -N $sgePre$outName.normal -hold_jid $holdJid -e $dir/$outName.normal.log -o $dir/$outName.normal.log \"bash $dir/$outName.normal.cmd\" > $dir/$outName.normal.log\n";
		print SUBFILE "qsub -cwd -b y -l h_vmem=16g -q $refHash->{sge_queue} -N $sgePre$outName.tumour -hold_jid $holdJid -e $dir/$outName.tumour.log -o $dir/$outName.tumour.log \"bash $dir/$outName.tumour.cmd\" > $dir/$outName.tumour.log\n";
		print SUBFILE "qsub -cwd -b y -l h_vmem=32g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $sgePre$outName.normal,$sgePre$outName.tumour -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.normal.cmd") or die;
		print COMMAND "module load $module; readCounter $normalBam > $dir/$normalName.wig\n";
		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;

		open (COMMAND, ">$dir/$outName.tumour.cmd") or die;
		print COMMAND "module load $module; readCounter $tumourBam > $dir/$tumourName.wig\n";
		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "module load $module; module load $Rmodule; $HMMscript $dir/$tumourName.wig $dir/$normalName.wig $dir > $dir/$outName.R; R < $dir/$outName.R --vanilla\n";
		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [HMMcopy] Skipping: $dir/$outName.touch\n";
	}

	my $hmmPre = $outName;
	$hmmPre =~ s/^(.*?_.*?)_.*/$1/;
	$tumourHash->{cnv}{final}{png_file} = "$dir/$hmmPre.genomeCNV.png";		# not exactly useful
	$tumourHash->{cnv}{final}{seg_file} = "$dir/$hmmPre*.segments_somatic.seg";
	$tumourHash->{cnv}{final}{segment_file} = "$dir/$hmmPre*.cnv_somatic_segments";

	$tumourHash->{cnv}{final}{tumour_wig} = "$dir/$tumourName.wig";
	$tumourHash->{cnv}{final}{normal_wig} = "$dir/$normalName.wig";

	$tumourHash->{cnv}{final}{hold_jid} = "$sgePre$outName";

	$tumourHash->{provenence}{"$dir/$outName.HMMcopy"}{command} = "$dir/$outName.normal.cmd,$dir/$outName.tumour.cmd,$dir/$outName.cmd";
	push (@{ $tumourHash->{provenence}{"$dir/$outName.HMMcopy"}{input}}, $normalBam);
	push (@{ $tumourHash->{provenence}{"$dir/$outName.HMMcopy"}{input}}, $tumourBam);

}




sub doCREST
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $normalHash = $_[2];
	my $tumourHash = $_[3];
	my $refHash = $_[4];

	my $seqType = $tumourHash->{sample_type};

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/$module";
	`mkdir -p $dir`;
	`mkdir -p $dir/extractSClip`;
	`mkdir -p $dir/tempSplit`;

	my $sgePre = "e$tumourHash->{sge_uniq}";

	my $normalBam = $normalHash->{bam}{$bamModule}{final}{file};
	my $tumourBam = $tumourHash->{bam}{$bamModule}{final}{file};
	my $holdJid = "$normalHash->{bam}{$bamModule}{final}{hold_jid},$tumourHash->{bam}{$bamModule}{final}{hold_jid}";

	my $fasta = $refHash->{"hg19_random"};
	my $blatRef = $refHash->{"blat ref"};
	my $blatHost = $refHash->{"blat host"};
	my $blatPort = $refHash->{"blat port"};

	my $outName = $tumourHash->{sample};

	my $tumourName = $tumourHash->{sample};
	my $normalName = $normalHash->{sample};

	my $catOrder = "";
	my $commandFileList = "";

	my $normalBai;
	my $tumourBai;

	my $perlLib = $refHash{'perl lib path'};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;

		# fix index names dagnabbit
		unless (-e "$normalBam.bai")
		{
			$normalBai = $normalBam;
			$normalBai =~ s/\.bam$/\.bai/;
			$normalBai =~ s/.*\///;
			`ln -s $normalBai $normalBam.bai`;
		}

		unless (-e "$tumourBam.bai")
		{
			$tumourBai = $tumourBam;
			$tumourBai =~ s/\.bam$/\.bai/;
			$tumourBai =~ s/.*\///;
			`ln -s $tumourBai $tumourBam.bai`;
		}

		for my $chr (qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y/)
		{
			print SUBFILE "qsub -cwd -b y -l h_vmem=8g -q $refHash->{sge_queue} -N $sgePre$outName.extract -hold_jid $holdJid -e $dir/extractSClip/$outName.chr$chr.log -o $dir/extractSClip/$outName.chr$chr.log \"bash $dir/extractSClip/$outName.chr$chr.cmd\" > $dir/extractSClip/$outName.chr$chr.log\n";

			open (COMMAND, ">$dir/extractSClip/$outName.chr$chr.cmd") or die;
			print COMMAND "$perlLib; module load $module; extractSClip.pl -i $tumourBam --ref_genome $fasta -r chr$chr -o $dir/extractSClip/\n";
			print COMMAND "\necho phoenixPipe/$module-done\n";
			close COMMAND;
		}

		print SUBFILE "qsub -cwd -b y -l h_vmem=4g -q $refHash->{sge_queue} -N $sgePre$outName.merge -hold_jid $sgePre$outName.extract -e $dir/extractSClip/$outName.log -o $dir/extractSClip/$outName.log \"bash $dir/extractSClip/$outName.cmd\" > $dir/extractSClip/$outName.log\n";

		open (COMMAND, ">$dir/extractSClip/$outName.cmd") or die;
		print COMMAND "cat $dir/extractSClip/*.cover > $dir/extractSClip/$tumourName.cover\n";
		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;

		for my $chr (qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y/)
		{
			print SUBFILE "qsub -q $refHash->{sge_queue} -cwd -b y -l h_vmem=96g -N $sgePre$outName -hold_jid $sgePre$outName.merge -e $dir/tempSplit/$outName.chr$chr.log -o $dir/tempSplit/$outName.chr$chr.log \"bash $dir/tempSplit/$outName.chr$chr.cmd\" > $dir/tempSplit/$outName.chr$chr.log\n";

			open (COMMAND, ">$dir/tempSplit/$outName.chr$chr.cmd") or die;
			print COMMAND "$perlLib; module load $module; CREST.pl -f $dir/extractSClip/$tumourName.cover -d $tumourBam -g $normalBam --ref_genome $fasta -t $blatRef --2bitdir / --blatserver $blatHost --blatport $blatPort -r chr$chr -p SV_$chr -o $dir/tempSplit/\n";
			print COMMAND "\necho phoenixPipe/$module-done\n";
			close COMMAND;

			$catOrder .= " $dir/tempSplit/SV_$chr.predSV.txt";
		}


		print SUBFILE "qsub -cwd -b y -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $sgePre$outName -e $dir/tempSplit/$outName.log -o $dir/tempSplit/$outName.log \"bash $dir/tempSplit/$outName.cmd\" > $dir/tempSplit/$outName.log\n";

		open (COMMAND, ">$dir/tempSplit/$outName.cmd") or die;
		print COMMAND "cat $catOrder > $dir/$outName.predSV.txt\n";
		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;
		

		close SUBFILE;
		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [CREST] Skipping: $dir/$outName.touch\n";
	}

	$tumourHash->{sv}{$module}{file} = "$dir/$outName.predSV.txt";
	$tumourHash->{sv}{$module}{hold_jid} = "$sgePre$outName";



	for my $chr (qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y/)
	{
		$commandFileList .= "$dir/extractSClip/$outName.chr$chr.cmd,";
	}
	$commandFileList .= "$dir/extractSClip/$outName.cmd,";
	for my $chr (qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y/)
	{
		$commandFileList .= "$dir/tempSplit/$outName.chr$chr.cmd,";
	}
	$commandFileList .= "$dir/extractSClip/$outName.cmd";		# last command, so no comma


	$tumourHash->{provenence}{"$dir/$outName.CREST"}{command} = $commandFileList;
	push (@{ $tumourHash->{provenence}{"$dir/$outName.CREST"}{input}}, $normalBam);
	push (@{ $tumourHash->{provenence}{"$dir/$outName.CREST"}{input}}, $tumourBam);

}




sub doFilterCREST
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $normalHash = $_[2];
	my $tumourHash = $_[3];
	my $refHash = $_[4];

	my $seqType = $tumourHash->{sample_type};

	my $filterScript = $refHash{"crest filter script"};

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/$module/filter";
	`mkdir -p $dir`;

	my $sgePre = "ef$tumourHash->{sge_uniq}";

	my $tumourBam = $tumourHash->{bam}{$bamModule}{final}{file};
	my $normalBam = $normalHash->{bam}{$bamModule}{final}{file};
	my $crestFile = $tumourHash->{sv}{$module}{file};
	my $holdJid = $tumourHash->{sv}{$module}{hold_jid};

	my $outName = $tumourHash->{sample};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;
		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=4g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "cat $crestFile | $filterScript $normalBam $tumourBam | grep PASS > $dir/$outName.filtered.predSV.txt; cat $crestFile | $filterScript $normalBam $tumourBam | grep FAIL > $dir/$outName.failed.predSV.txt\n";
		print COMMAND "\necho phoenixPipe/$module/filter-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [CREST_FITLER] Skipping: $dir/$outName.touch\n";
	}

	$tumourHash->{sv}{$module}{filter_file} = "$dir/$outName.filtered.predSV.txt";
	$tumourHash->{sv}{$module}{filter_hold_jid} = "$sgePre$outName";


#	$tumourHash->{provenence}{"$dir/$outName.filtered.predSV.txt"}{command} = "$dir/$outName.cmd";
#	push (@{ $tumourHash->{provenence}{"$dir/$outName.filtered.predSV.txt"}{input}}, "$dir/$outName.CREST");

}









sub doDelly
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $normalHash = $_[2];
	my $tumourHash = $_[3];
	my $refHash = $_[4];

	my $seqType = $tumourHash->{sample_type};

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/$module";
	`mkdir -p $dir`;

	my $sgePre = "d$tumourHash->{sge_uniq}";

	my $normalBam = $normalHash->{bam}{$bamModule}{final}{file};
	my $tumourBam = $tumourHash->{bam}{$bamModule}{final}{file};
	my $holdJid = "$normalHash->{bam}{$bamModule}{final}{hold_jid},$tumourHash->{bam}{$bamModule}{final}{hold_jid}";

	my $fasta = $refHash->{"hg19_random"};

	my $outName = $tumourHash->{sample};

	my $tumourName = $tumourHash->{sample};
	my $normalName = $normalHash->{sample};

	my %memReqs = (
		"DEL" => "32g",
		"DUP" => "16g",
		"INV" => "16g",
		"TRA" => "16g",
	);

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;

		open (SUBFILE, ">$dir/$outName.sub") or die;

		for my $svType (qw/DEL DUP INV TRA/)
		{

			print SUBFILE "qsub -cwd -b y -l h_vmem=$memReqs{$svType} -q $refHash->{sge_queue} -N $sgePre$outName.$svType -hold_jid $holdJid -e $dir/$outName.$svType.log -o $dir/$outName.$svType.log \"bash $dir/$outName.$svType.cmd\" > $dir/$outName.$svType.log\n";

			open (COMMAND, ">$dir/$outName.$svType.cmd") or die;
			print COMMAND "module load $module; delly -t $svType -o $dir/$outName.$svType.vcf -g $fasta $tumourBam $normalBam\n";
			print COMMAND "\necho phoenixPipe/$module/$svType-done\n";
			close COMMAND;
		}

		close SUBFILE;
		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [DELLY] Skipping: $dir/$outName.touch\n";
	}


	$tumourHash->{sv}{$module}{file} = "$dir/$outName";
	$tumourHash->{sv}{$module}{hold_jid} = "$sgePre$outName.DEL,$sgePre$outName.DUP,$sgePre$outName.INV,$sgePre$outName.TRA";

#	$tumourHash->{provenence}{"$dir/$outName.CREST"}{command} = $commandFileList;

#	push (@{ $tumourHash->{provenence}{"$dir/$outName.CREST"}{input}}, $normalBam);
#	push (@{ $tumourHash->{provenence}{"$dir/$outName.CREST"}{input}}, $tumourBam);

}





sub doFilterDelly
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $tumourHash = $_[2];
	my $refHash = $_[3];

	my $seqType = $tumourHash->{sample_type};

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/$module/filter";
	`mkdir -p $dir`;

	my $sgePre = "df$tumourHash->{sge_uniq}";

	my $dellyDir = $tumourHash->{sv}{$module}{file};
	my $holdJid = $tumourHash->{sv}{$module}{hold_jid};

	my $outName = $tumourHash->{sample};

	my $python = $refHash->{"gavin python"};		# required for filter script dependencies
	my $filterScript = $refHash->{"delly filter script"};

	my $tumourName = $tumourHash->{sample};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;

		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=2g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "cd $dir; ln -s ../*.vcf .; $python $filterScript; rm $outName.DEL.vcf $outName.DUP.vcf $outName.INV.vcf $outName.TRA.vcf; (grep ^# $outName.DEL.filtered.vcf; cat *.vcf | grep -v ^#) | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/chrM/chr25/' | sed 's/^chr//' | sort -n | sed 's/^/chr/' | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' | sed 's/chr25/chrM/' > $dir/$outName.merged.vcf\n";
		print COMMAND "\necho phoenixPipe/$module/filter-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [DELLY_FILTER] Skipping: $dir/$outName.touch\n";
	}


	$tumourHash->{sv}{$module}{filter_file} = "$dir/$outName.merged.vcf";
	$tumourHash->{sv}{$module}{filter_hold_jid} = "$sgePre$outName";

#	$tumourHash->{provenence}{"$dir/$outName.CREST"}{command} = $commandFileList;

#	push (@{ $tumourHash->{provenence}{"$dir/$outName.CREST"}{input}}, $normalBam);
#	push (@{ $tumourHash->{provenence}{"$dir/$outName.CREST"}{input}}, $tumourBam);

}


sub doMergeSV
{
	my $module1 = $_[0];
	my $module2 = $_[1];
	my $bamModule = $_[2];
	my $tumourHash = $_[3];
	my $refHash = $_[4];

	my $seqType = $tumourHash->{sample_type};

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tumourHash);
	my $dir = "$pathRoot/$bamModule/final_crest-delly/";
	`mkdir -p $dir`;

	my $sgePre = "sv$tumourHash->{sge_uniq}";

	my $svFile1 = $tumourHash->{sv}{$module1}{filter_file};
	my $svFile2 = $tumourHash->{sv}{$module2}{filter_file};

	if ($module2 =~ /delly/)
	{
		$svFile2 =~ s/\.merged\.vcf$//;
	}

	my $holdJid = "$tumourHash->{sv}{$module1}{filter_hold_jid},$tumourHash->{sv}{$module2}{filter_hold_jid}";

	my $outName = $tumourHash->{sample};

	my $python = $refHash->{"gavin python"};		# required for filter script dependencies
	my $mergeScript = $refHash->{"sv merge script"};
	my $annotateScript = $refHash->{"sv annotate script"};
	my $refSeqGenes = $refHash->{"refseq genes tabix"};

	my $perlLib = $refHash->{'perl lib path'};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;

		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=2g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "$perlLib; $python $mergeScript $svFile1 $svFile2 $dir/$outName.unannotatedSV.tsv; cat $dir/$outName.unannotatedSV.tsv | $annotateScript $refSeqGenes > $dir/$outName.annotatedSV.tsv\n";
		
		print COMMAND "\necho phoenixPipe/svMerge-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [SV_MERGE] Skipping: $dir/$outName.touch\n";
	}


	$tumourHash->{sv}{final}{file} = "$dir/$outName.annotatedSV.tsv";
	$tumourHash->{sv}{final}{hold_jid} = "$sgePre$outName";

#	$tumourHash->{provenence}{"$dir/$outName.CREST"}{command} = $commandFileList;

#	push (@{ $tumourHash->{provenence}{"$dir/$outName.CREST"}{input}}, $normalBam);
#	push (@{ $tumourHash->{provenence}{"$dir/$outName.CREST"}{input}}, $tumourBam);

}



sub doATHLATES
{
	my $module = $_[0];
	my $nHash = $_[1];
	my $tHash = $_[2];
	my $refHash = $_[3];

	my $novoModule = "novocraft/2.07.06";		# may depend on these specific versions...
	my $samtoolsModule = "samtools/0.1.18";

	my $novoRef = $refHash->{"ATHLATES_$novoModule"};
	my $bedPath = $refHash->{ATHLATES_bed_path};
	my $msaPath = $refHash->{ATHLATES_msa_path};

	my $pathRoot = buildPath($nHash);
	unless ($tHash eq "single")
	{
		$pathRoot = buildPath($tHash);
	}
	my $dir = "$pathRoot/$module/";
	`mkdir -p $dir`;

	`mkdir -p $dir/bam`;
	`mkdir -p $dir/bed`;

	my $sgePre = "at$nHash->{sge_uniq}";
	my $outName;
	my $outPath;

	my $files;
	my $command;
	my $hold_jid;

	my $fastqFormat = "";

	my $laneCount = 0;
	my $laneFileList = "";
	my $laneJobList = "";

	for my $lane (sort keys %{ $nHash->{fastq} })
	{
		$outPath = "bam/";
		$outName = "$nHash->{fastq}{$lane}{library}_$nHash->{fastq}{$lane}{barcode}_$nHash->{fastq}{$lane}{lane}_$nHash->{fastq}{$lane}{run}";
		
		$laneCount++;
		$laneFileList .= " $dir/$outPath/$outName.sort.bam";
		$laneJobList .= ",$sgePre$outName";

		if ($nHash->{fastq}{$lane}{read2} =~ /.*Not found/)
		{
			if ($nHash->{tissue_type} eq "X")
			{
				$files = "$nHash->{fastq}{$lane}{xenome}{read1}";
				$hold_jid = "-hold_jid $nHash->{fastq}{$lane}{xenome}{hold_jid}";
			}
			else
			{
				$files = "$nHash->{fastq}{$lane}{read1}";
				$hold_jid = "";
			}
		}
		else
		{
			if ($nHash->{tissue_type} eq "X")
			{
				$files = "$nHash->{fastq}{$lane}{xenome}{read1} $nHash->{fastq}{$lane}{xenome}{read2}";
				$hold_jid = "-hold_jid $nHash->{fastq}{$lane}{xenome}{hold_jid}";
			}
			else
			{
				$files = "$nHash->{fastq}{$lane}{read1} $nHash->{fastq}{$lane}{read2}";
				$hold_jid = "";
			}
		}

		if ($nHash->{fastq}{$lane}{"quality_format"} eq "illumina")
		{
			$fastqFormat = "-F ILMFQ";
		}
		elsif ($nHash->{fastq}{$lane}{"quality_format"} eq "sanger")
		{
			$fastqFormat = "";
		}

		unless (-e "$dir/$outPath/$outName.touch")
		{
			`touch $dir/$outPath/$outName.touch`;		# so that future runs won't cause concurrent processing

			open (SUBFILE, ">$dir/$outPath/$outName.sub") or die;
			print SUBFILE "qsub -cwd -b y -q $refHash->{sge_queue} -N $sgePre$outName $hold_jid -e $dir/$outPath/$outName.log -o $dir/$outPath/$outName.log -l h_vmem=24g -l exclusive=1 -q $refHash->{sge_queue} \"bash $dir/$outPath/$outName.cmd\" > $dir/$outPath/$outName.log\n";
			close SUBFILE;

			open (COMMAND, ">$dir/$outPath/$outName.cmd") or die;
			print COMMAND "module load $novoModule; module load $samtoolsModule; novoalign $fastqFormat -d $novoRef -f $files -t 30 -o SAM -r ALL -l 80 -e 100 | samtools view -bS -h -F 4 - > $dir/$outPath/$outName.bam\n";
			print COMMAND "samtools sort $dir/$outPath/$outName.bam $dir/$outPath/$outName.sort\n";
			print COMMAND "\necho phoenixPipe/$module-done\n";
			close COMMAND;

			`bash $dir/$outPath/$outName.sub`;
		}
		else
		{
			warn " [ATHLATES] Skipping alignment: $dir/$outPath/$outName.touch\n";
		}
	}

	$laneFileList =~ s/^ //;
	$laneJobList =~ s/^,//;

	$outName = $nHash->{sample};
	unless ($tHash eq "single")
	{
		$outName = $tHash->{sample};
	}
	$outPath = "$dir/bed";

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;

		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $laneJobList -e $dir/$outName.log -o $dir/$outName.log -l h_vmem=2g \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;

		if ($laneCount > 1)
		{
			print COMMAND "module load $module; module load $samtoolsModule; samtools merge $dir/$outName.bam $laneFileList; samtools sort $dir/$outName.bam $dir/$outName.sort\n";
		}
		else
		{
			print COMMAND "module load $module; module load $samtoolsModule; ln -s $laneFileList $dir/$outName.sort.bam\n";
		}

		for my $i (qw/A non-A B non-B C non-C/)
		{
			print COMMAND "samtools view -b -L $bedPath/hla.$i.bed $dir/$outName.sort.bam > $outPath/$i.bam; samtools view -h -o $outPath/$i.sam $outPath/$i.bam; (grep \"^@\" $outPath/$i.sam; grep -v \"^@\" $outPath/$i.sam | sort -k1,1 -k3,3) > $outPath/$i.sort.sam; samtools view -bS $outPath/$i.sort.sam > $outPath/$i.sort.bam\n";
		}

		for my $i (qw/A B C/)
		{
			`mkdir $dir/HLA-$i`;
			print COMMAND "typing -bam $outPath/$i.sort.bam -exlbam $outPath/non-$i.sort.bam -msa $msaPath/${i}_nuc.txt -o $dir/HLA-$i/output.$i\n";
		}
		print COMMAND "cat $dir/HLA*/*.typing.txt | athlatesToTypes.pl > $dir/$outName.HLA.txt\n";
		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`
	}
	else
	{
		warn " [ATHLATES] Skipping: $dir/$outName.touch\n";
	}


	$nHash->{hla}{$module}{file} = "$dir/$outName.HLA.txt";
	$nHash->{hla}{$module}{hold_jid} = "$sgePre$outName";

}



sub doNetMHC
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $hlaModule = $_[2];
	my $nHash = $_[3];
	my $tHash = $_[4];
	my $refHash = $_[5];

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tHash);
	my $dir = "$pathRoot/$bamModule/$module/$hlaModule/";
	`mkdir -p $dir`;

	my $sgePre = "na$tHash->{sge_uniq}";

	my $hlaFile = $nHash->{hla}{$hlaModule}{file};
	my $vcfFile = $tHash->{vcf}{final}{file};

	my $holdJid = "$nHash->{hla}{$hlaModule}{hold_jid},$tHash->{vcf}{final}{hold_jid}";

	my $outName = $tHash->{sample};

	my $vcfToPeptide = $refHash->{vcf_to_peptide};
	my $perlLib = $refHash->{"perl lib path"};

	my $program = "netMHC";
	if ($module =~ /pan/)
	{
		$program = "netMHCpan";
	}

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;

		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=2g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "$perlLib; module load $module; cat $vcfFile | $vcfToPeptide > $dir/$outName.peptideMap; cut -f 5 $dir/$outName.peptideMap | sort | uniq > $dir/$outName.peptides; for i in `cat $hlaFile | sed 's/\\*//' | sed 's/^/HLA-/'`; do $program -a \$i -p $dir/$outName.peptides > $dir/$outName.\$i.txt; done\n";
		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [netMHC] Skipping: $dir/$outName.touch\n";
	}


}



sub doPOLYSOLVER
{
	my $module = $_[0];
	my $bamModule = $_[1];
	my $nHash = $_[2];
	my $tHash = $_[3];
	my $refHash = $_[4];

	my $pathRoot = buildPath($nHash);
	unless ($tHash eq "single")
	{
		$pathRoot = buildPath($tHash);
	}
	my $dir = "$pathRoot/$bamModule/$module/";
	`mkdir -p $dir`;

	my $sgePre = "ps$nHash->{sge_uniq}";

	my $file = $nHash->{bam}{$bamModule}{final}{file};
	my $hold_jid = $nHash->{bam}{$bamModule}{final}{hold_jid};

	my $outName = $nHash->{sample};
	unless ($tHash eq "single")
	{
		$outName = $tHash->{sample};
	}

	my $perlLib = $refHash{'perl lib path'};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;

		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=24g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $hold_jid -e $dir/$outName.log -o $dir/$outName.log  \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "$perlLib; module load $module; bash \$POLYSOLVERROOT/scripts/shell_call_hla_type $file Unknown 1 hg19 STDFQ 0 $dir/out; cat $dir/out/winners.hla.txt | sed 's/hla_/\\n/g' | sed 's/a_/A*/' | sed 's/b_/B*/' | sed 's/c_/C*/' | sed 's/_/:/g' | sed 's/  //' | grep -v \"-\" | cut -f 1,2 -d \":\" | uniq > $dir/$outName.HLA.txt; rm $dir/out/*.sam $dir/out/*.bam $dir/out/*.lik1 $dir/out/*.lik2\n";
		print COMMAND "\necho phoenixPipe/$module-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`
	}
	else
	{
		warn " [POLYSOLVER] Skipping: $dir/$outName.touch\n";
	}


	$nHash->{hla}{$module}{file} = "$dir/$outName.HLA.txt";
	$nHash->{hla}{$module}{hold_jid} = "$sgePre$outName";

}

sub doCelluloidMC
{
	my $rModule = $_[0];
	my $bamModule = $_[1];
	my $tHash = $_[2];
	my $refHash = $_[3];

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tHash);
	my $dir = "$pathRoot/$bamModule/celluloid/mchan/";
	`mkdir -p $dir`;

	my $sgePre = "cm$tHash->{sge_uniq}";

	my $germlineFile = $tHash->{vcf}{germline}{file};
	my $tumourWig = $tHash->{cnv}{final}{tumour_wig};
	my $normalWig = $tHash->{cnv}{final}{normal_wig};

	my $gcWig = $refHash->{gc_wig};
	my $mapWig = $refHash->{map_wig};

	my $outName = $tHash->{sample};

	my $germlineAF = "$outName.af.tsv";

	my $holdJid = "$tHash->{vcf}{germline}{hold_jid},$tHash->{cnv}{final}{hold_jid}";

	my $gatkToAF = $refHash->{celluloid_gatk_to_af};
	my $celluloidScript = $refHash->{celluloid_mchan};

	my $perlLib = $refHash{'perl lib path'};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;

		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=16g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "$perlLib; module load $rModule; cd $dir; $gatkToAF $germlineFile $tHash->{sample} $germlineAF 30; Rscript $celluloidScript -T $tumourWig -N $normalWig -G $gcWig -M $mapWig -A $germlineAF -O $outName\n";
		print COMMAND "\necho phoenixPipe/celluloidMC-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [CELLULOID_MCHAN] Skipping: $dir/$outName.touch\n";
	}

	$tHash->{celluloid_mchan}{af_file} = "$dir/$germlineAF";
	$tHash->{celluloid_mchan}{seg_file} = "$dir/*$outName.seg";
	$tHash->{celluloid_mchan}{segments_file} = "$dir/*$outName.segments";
	$tHash->{celluloid_mchan}{hold_jid} = "$sgePre$outName";
	

}


sub doCelluloid
{
	my $rModule = $_[0];
	my $bamModule = $_[1];
	my $tHash = $_[2];
	my $refHash = $_[3];

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tHash);
	my $dir = "$pathRoot/$bamModule/celluloid/v11.2/";
	`mkdir -p $dir`;

	my $sgePre = "cl$tHash->{sge_uniq}";

	my $germlineFile = $tHash->{vcf}{germline}{file};
	my $tumourWig = $tHash->{cnv}{final}{tumour_wig};
	my $normalWig = $tHash->{cnv}{final}{normal_wig};

	my $gcWig = $refHash->{gc_wig};
	my $mapWig = $refHash->{map_wig};

	my $outName = $tHash->{sample};

	my $germlineAF = "$outName.af.tsv";

	my $holdJid = "$tHash->{vcf}{germline}{hold_jid},$tHash->{cnv}{final}{hold_jid}";

	my $gatkToAF = $refHash->{celluloid_gatk_to_af};
	my $celluloidScript = $refHash->{celluloid_v11};
	my $tabixScript = $refHash->{tabix_bed_script};

	my $perlLib = $refHash{'perl lib path'};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;

		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=16g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "$perlLib; module load $rModule; cd $dir; $gatkToAF $germlineFile $tHash->{sample} $germlineAF 30; Rscript $celluloidScript $tumourWig $normalWig $gcWig $mapWig $germlineAF $outName; for i in `ls solution*/segments_*.txt`; do bash $tabixScript \$i; done; ln -s solution1 solution; \n";
		print COMMAND "\necho phoenixPipe/celluloid-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [CELLULOID] Skipping: $dir/$outName.touch\n";
	}

#	$tHash->{celluloid_mchan}{af_file} = "$dir/$germlineAF";
#	$tHash->{celluloid_mchan}{seg_file} = "$dir/*$outName.seg";
#	$tHash->{celluloid_mchan}{segments_file} = "$dir/*$outName.segments";
#	$tHash->{celluloid_mchan}{hold_jid} = "$sgePre$outName";
	

}


sub doChromothripsis
{
	my $rModule = $_[0];
	my $bamModule = $_[1];
	my $tHash = $_[2];
	my $refHash = $_[3];

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tHash);
	my $dir = "$pathRoot/$bamModule/celluloid/mchan/chromothripsis";
	`mkdir -p $dir`;

	my $sgePre = "cr$tHash->{sge_uniq}";

	my $germlineFile = $tHash->{celluloid_mchan}{af_file};
	my $cnaSegFile = $tHash->{celluloid_mchan}{seg_file};
	my $cnaSegmentFile = $tHash->{celluloid_mchan}{segments_file};
	my $svFile = $tHash->{sv}{final}{file};

	my $outName = $tHash->{sample};

	my $germlineAF = "$outName.af.tsv";

	my $holdJid = "$tHash->{celluloid_mchan}{hold_jid},$tHash->{sv}{final}{hold_jid}";

	my $gatkToAF = "sed 's/ /	/g'";

	my $chromoScript = $refHash->{chromo_script};

	my $perlLib = $refHash{'perl lib path'};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;

		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=16g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "$perlLib; module load $rModule; cd $dir; cat $germlineFile | $gatkToAF > $germlineAF; Rscript $chromoScript -t $svFile -n $cnaSegmentFile -s $cnaSegFile -a $germlineAF -o $outName\n";
		print COMMAND "\necho phoenixPipe/chromothripsis-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [CHROMOTHRIPSIS] Skipping: $dir/$outName.touch\n";
	}


}


sub doQC
{
	my $bamModule = $_[0];
	my $tHash = $_[1];
	my $nHash = $_[2];
	my $refHash = $_[3];
	my $commandLine = $_[4];

	# put qc in the sample root
	my $path = buildPath($tHash);

	`mkdir -p $path`;

	my $dir = $path;

	my $sgePre = "qc$tHash->{sge_uniq}";
	my $outName = $tHash->{sample};
	my $holdJid = "\\*$outName\\*";

	my $qcScript = $refHash->{"qc script"};
	my $perlLib = $refHash{'perl lib path'};

	open (SUBFILE, ">$dir/$outName.qc.sub") or die;
	print SUBFILE "qsub -cwd -b y -l h_vmem=2g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.qc.log -o $dir/$outName.qc.log \"bash $dir/$outName.qc.cmd\" > $dir/$outName.qc.log\n";
	close SUBFILE;

	open (COMMAND, ">$dir/$outName.qc.cmd") or die;
	print COMMAND "$perlLib; $qcScript $path/$bamModule tumour > $dir/$outName.qc.txt\n";
	print COMMAND "\necho phoenixPipe/qc-done\n";
	close COMMAND;

	open (COMMANDLINE, ">$dir/$outName.commandLine") or die;
	print COMMANDLINE "$commandLine\n";
	close COMMANDLINE;

	`bash $dir/$outName.qc.sub`;



	$path = buildPath($nHash);
	`mkdir -p $path`;

	$dir = $path;

	$sgePre = "qc$nHash->{sge_uniq}";
	$outName = $nHash->{sample};

	$holdJid = "\*$outName\*";

	open (SUBFILE, ">$dir/$outName.qc.sub") or die;
	print SUBFILE "qsub -cwd -b y -l h_vmem=2g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.qc.log -o $dir/$outName.qc.log \"bash $dir/$outName.qc.cmd\" > $dir/$outName.qc.log\n";
	close SUBFILE;

	open (COMMAND, ">$dir/$outName.qc.cmd") or die;
	print COMMAND "$perlLib; $qcScript $path/$bamModule normal > $dir/$outName.qc.txt\n";
	print COMMAND "\necho phoenixPipe/qc-done\n";
	close COMMAND;

	`bash $dir/$outName.qc.sub`;
}


sub doQCsingle
{
	my $bamModule = $_[0];
	my $nHash = $_[1];
	my $refHash = $_[2];
	my $commandLine = $_[3];

	# put qc in the sample root

	my $path = buildPath($nHash);
	`mkdir -p $path`;

	my $dir = $path;

	my $perlLib = $refHash{'perl lib path'};

	my $sgePre = "qc$nHash->{sge_uniq}";
	my $outName = $nHash->{sample};

	my $holdJid = "\*$outName\*";

	my $qcScript = $refHash->{"qc script"};

	open (SUBFILE, ">$dir/$outName.qc.sub") or die;
	print SUBFILE "qsub -cwd -b y -l h_vmem=2g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.qc.log -o $dir/$outName.qc.log \"bash $dir/$outName.qc.cmd\" > $dir/$outName.qc.log\n";
	close SUBFILE;

	open (COMMAND, ">$dir/$outName.qc.cmd") or die;
	print COMMAND "$perlLib; $qcScript $path/$bamModule normal > $dir/$outName.qc.txt\n";
	print COMMAND "\necho phoenixPipe/qc-done\n";
	close COMMAND;

	open (COMMANDLINE, ">$dir/$outName.commandLine") or die;
	print COMMANDLINE "$commandLine\n";
	close COMMANDLINE;

	`bash $dir/$outName.qc.sub`;


}


sub doPhoenixParser
{
	my $bamModule = $_[0];
	my $tHash = $_[1];
	my $refHash = $_[2];

	# put qc in the sample root

	my $path = buildPath($tHash);
	my $dir = "$path/$bamModule/results";
	`mkdir -p $dir`;

	my $perlLib = $refHash{'perl lib path'};

	my $sgePre = "pp$tHash->{sge_uniq}";
	my $outName = $tHash->{sample};

	my $holdJid = "\\\*$outName\\\*";

	my $ppScript = $refHash->{"phoenix parser"};

	my $rModule = $refHash->{"R module"};
	my $prScript = $refHash->{"phoenix reporter"};
	my $srScript = $refHash->{"short reporter"};
	my $plotScript = $refHash->{"report plots"};

	open (SUBFILE, ">$dir/$outName.sub") or die;
	print SUBFILE "qsub -cwd -b y -l h_vmem=2g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
	close SUBFILE;

	open (COMMAND, ">$dir/$outName.cmd") or die;
	print COMMAND "$perlLib; $ppScript $tHash->{working_dir} $tHash->{sample_group} $tHash->{sample} $tHash->{sample_type} $dir\n";
	print COMMAND "mkdir $dir/plots; module load $rModule; $plotScript $dir/$outName.summary.csv $dir/$outName.variants.csv $dir/plots/$outName; $prScript $dir/$outName.summary.csv $dir/$outName.variants.csv $dir; $srScript $dir/$outName.summary.csv $dir/$outName.variants.csv $dir\n";
	print COMMAND "\necho phoenixPipe/phoenixParser-done\n";
	close COMMAND;

	`bash $dir/$outName.sub`;


}


sub doMoreAnnovar
{
	my $bamModule = shift;
	my $module = shift;
	my $tHash = shift;
	my $refHash = shift;

	my $path = buildPath($tHash);
	my $dir = "$path/$bamModule/$module/annovar";
	`mkdir -p $dir`;

	my $sgePre = "ma$tHash->{sge_uniq}";
	my $outName = $tHash->{sample};

	my $vcfFile = $tHash->{vcf}{final}{file};
	my $holdJid = $tHash->{vcf}{final}{hold_jid};
	my $vmem = "16g";
	my $perlLib = $refHash{'perl lib path'};

	if ($module =~ /germline/)
	{
		$vcfFile = $tHash->{vcf}{germline}{file};
		$holdJid = $tHash->{vcf}{germline}{hold_jid};
		$vmem = "16g";
	}

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;

		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=$vmem -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "$perlLib; module load annovar/2014-07-15\n";
		for my $chr (qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y/)
		{
			print COMMAND "(grep ^# $vcfFile; grep -w chr$chr $vcfFile) > $dir/$outName.chr$chr.vcf; convert2annovar.pl -format vcf4old $dir/$outName.chr$chr.vcf -outfile $dir/$outName.chr$chr.input.vcf -includeinfo; table_annovar.pl $dir/$outName.chr$chr.input.vcf /oicr/local/analysis/sw/annovar/humandb/ -buildver hg19 -remove -protocol refGene,avsnp142,popfreq_all_20150413,popfreq_max_20150413,ljb26_all,clinvar_20150330 -operation g,f,f,f,f,f -otherinfo -nastring NA -out $dir/$outName.chr$chr.output.vcf\n";
		}
		print COMMAND "echo \\# `head -n 1 $dir/$outName.chr1.output.vcf.hg19_multianno.txt` > $dir/$outName.output.vcf.hg19_multianno.txt; cat";
		for my $chr (qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y/)
		{
			print COMMAND " $dir/$outName.chr$chr.output.vcf.hg19_multianno.txt";
		}
		print COMMAND " | grep -v Func.refGene >> $dir/$outName.output.vcf.hg19_multianno.txt\n";

		print COMMAND "module load tabix; cat $dir/$outName.output.vcf.hg19_multianno.txt | bgzip > $dir/$outName.output.vcf.hg19_multianno.txt.gz; tabix -p bed $dir/$outName.output.vcf.hg19_multianno.txt.gz\n";

		print COMMAND "rm $dir/$outName.*.txt $dir/$outName.*.vcf $dir/$outName.*.invalid_input\n";

		print COMMAND "\necho phoenixPipe/more_annovar-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [MORE_ANNOVAR] Skipping $dir/$outName.touch\n"
	}
}




sub doCosmicSigNNLS
{
	my $rModule = $_[0];
	my $bamModule = $_[1];
	my $tHash = $_[2];
	my $refHash = $_[3];

	# put somatic calls in the tumour dir
	my $pathRoot = buildPath($tHash);
	my $dir = "$pathRoot/$bamModule/cosmicSigNNLS";
	`mkdir -p $dir`;

	my $sgePre = "cs$tHash->{sge_uniq}";

	my $vcfFile = $tHash->{vcf}{final}{file};

	my $outName = $tHash->{sample};


	my $holdJid = "$tHash->{vcf}{final}{hold_jid}";

	my $perlLib = $refHash{'perl lib path'};

	my $script = $refHash->{"cosmicSigNNLS script"};

	unless (-e "$dir/$outName.touch")
	{
		`touch $dir/$outName.touch`;

		open (SUBFILE, ">$dir/$outName.sub") or die;
		print SUBFILE "qsub -cwd -b y -l h_vmem=4g -q $refHash->{sge_queue} -N $sgePre$outName -hold_jid $holdJid -e $dir/$outName.log -o $dir/$outName.log \"bash $dir/$outName.cmd\" > $dir/$outName.log\n";
		close SUBFILE;

		open (COMMAND, ">$dir/$outName.cmd") or die;
		print COMMAND "$perlLib; module load $rModule; cd $dir; $script $vcfFile $dir $outName\n";
		print COMMAND "\necho phoenixPipe/cosmicSigNNLS-done\n";
		close COMMAND;

		`bash $dir/$outName.sub`;
	}
	else
	{
		warn " [COSMICSIGNNLS] Skipping: $dir/$outName.touch\n";
	}


}







sub buildPath
{
	my $metaHash = $_[0];

	my $path = "$metaHash->{working_dir}/$metaHash->{sample_group}/$metaHash->{sample}/$metaHash->{sample_type}";
	return $path;
}


sub doClean
{
	my $metaHash = $_[0];
	my $doClean = $_[1];

	my $pathRoot = buildPath($metaHash);


	open (COMMAND, ">$pathRoot/$metaHash->{sample}-$metaHash->{sample_type}_cleanTempFiles.cmd") or die "Couldn't open >$pathRoot/$metaHash->{sample}-$metaHash->{sample_type}_cleanTempFiles.cmd\n";
	for my $file (@{ $metaHash->{clean}{temp}{files} })
	{
		print COMMAND "rm -rf $file\n";
	}
	close COMMAND;

	my $hold_jids = "";
	for my $jid (@{ $metaHash->{clean}{temp}{hold_jids} })
	{
		$hold_jids .= "$jid,";
	}
	open (SUBFILE, ">$pathRoot/$metaHash->{sample}-$metaHash->{sample_type}_cleanTempFiles.sub") or die;
	print SUBFILE "qsub -cwd -b y -N z$metaHash->{sample}-$metaHash->{sample_type}_cleanTempFiles -hold_jid $hold_jids \"bash $pathRoot/$metaHash->{sample}-$metaHash->{sample_type}_cleanTempFiles.cmd\"\n";
	close SUBFILE;
	
	if ($doClean eq "clean")
	{
		`bash $pathRoot/$metaHash->{sample}/$metaHash->{sample_type}_cleanTempFiles.sub`;
	}

	open (COMMAND, ">$pathRoot/$metaHash->{sample}-$metaHash->{sample_type}_cleanAllFiles.cmd") or die;
	for my $file (@{ $metaHash->{clean}{temp}{files} })
	{
		print COMMAND "rm -rf $file\n";
	}
	for my $file (@{ $metaHash->{clean}{full} })
	{
		print COMMAND "rm -rf $file\n";
	}
	close COMMAND;
}




sub doProvenance
{
	my $normalHash = $_[0];
	my $tumourHash = $_[1];

	my $l;
	my $commands;

	for my $output (keys %{ $normalHash->{provenence} })
	{
		$commands = getCommands($output, $normalHash, $tumourHash, 0);		# 0 is the depth, part of the recursive formatting

		open (FILE, ">$output.provenence.txt") or die "Couldn't open >$output.provenence.txt\n";

		print FILE "\n";
		print FILE "Provenence For:  $output\n\n";
		print FILE "Normal Sample:   $normalHash->{sample}\n";
		print FILE "Sequencing Type: $normalHash->{sample_type}\n";
		print FILE "Submitting User: $normalHash->{username}\n\n";

		print FILE "Command History:\n";
		print FILE "$commands\n";

		close FILE;
	}

	for my $output (keys %{ $tumourHash->{provenence} })
	{
		$commands = getCommands($output, $normalHash, $tumourHash, 0);		# 0 is the depth, part of the recursive formatting

		open (FILE, ">$output.provenence.txt") or die "Couldn't open >$output.provenence.txt\n";

		print FILE "\n";
		print FILE "Provenence For:  $output\n\n";
		print FILE "Normal Sample:   $normalHash->{sample}\n";
		print FILE "Tumour Sample:   $tumourHash->{sample}\n";
		print FILE "Sequencing Type: $tumourHash->{sample_type}\n";
		print FILE "Submitting User: $tumourHash->{username}\n\n";

		print FILE "Command History:\n";
		print FILE "$commands\n";

		close FILE;
	}

}


sub getCommands
{
	my $file = $_[0];
	my $normalHash = $_[1];
	my $tumourHash = $_[2];
	my $depth = $_[3];

	my $command = "";
	my $l;

	my $pad = "";
	for (my $i = 0; $i < $depth; $i++)
	{
		$pad .= "\t";
	}


	# open file and get commands
	if (exists $normalHash->{provenence}{$file}{command})
	{
		for my $cmd (split(/,/, $normalHash->{provenence}{$file}{command}))
		{
			open (COMMAND, "<$cmd") or die "Couldn't open $cmd for $file\n";
			while ($l = <COMMAND>)
			{
				$command .= "$pad$l";
			}
		}
	}
	elsif (exists $tumourHash->{provenence}{$file}{command})
	{
		for my $cmd (split(/,/, $tumourHash->{provenence}{$file}{command}))
		{
			open (COMMAND, "<$cmd") or die "Couldn't open $cmd for $file\n";
			while ($l = <COMMAND>)
			{
				$command .= "$pad$l";
			}
		}
	}
	elsif (exists $normalHash->{provenence}{$file}{fastq})
	{
		$command .= "$pad$normalHash->{provenence}{$file}{fastq}";
	}
	elsif (exists $tumourHash->{provenence}{$file}{fastq})
	{
		$command .= "$pad$tumourHash->{provenence}{$file}{fastq}";
	}
	else
	{
		die "Provenence disconnect - couldn't find command for $file\n";
	}
	
	close COMMAND;

	# if input, getCommands
	if (exists $normalHash->{provenence}{$file}{input})
	{
		for my $input (@{ $normalHash->{provenence}{$file}{input} })
		{
			$command .= getCommands($input, $normalHash, $tumourHash, $depth + 1);
		}
	}
	if (exists $tumourHash->{provenence}{$file}{input})
	{
		for my $input (@{ $tumourHash->{provenence}{$file}{input} })
		{
			$command .= getCommands($input, $normalHash, $tumourHash, $depth + 1);
		}
	}
	
	
	return $command;
	
}




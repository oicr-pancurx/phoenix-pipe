#!/bin/bash 
#set -x
#
# Shell script written for GATK 
#
VERSION="$Revision$"
bamname=""

function version ()
{
        printf "$PROGNAME $VERSION\n"
        exit 2
}


# Find out how I was called and with what arguments
PROGNAME=`basename $0`
PROGARGS=$@

#usage
function usage ()
{
   printf "\n\n"
   printf "usage: $PROGNAME \n\n"
	printf "\t-b directory with bam files\n" 
   printf "\t-f fasta reference (default /oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/gatk/resources/hg19_random.fa) \n"
	printf "\t-d dbSNP file (default /oicr/data/reference/genomes/homo_sapiens_mc/dbSNP/hg19_random/Genomic/dbSNP142/dbsnp142_chr.vcf)\n"
   printf "\t-t sequencing type (ie exome or wgs -- default is exome)\n"
	printf "\t-v gatk version (default to 1.3.16)\n"
	printf "\t-p picard version (default to 1.40)\n"
	printf "\t-j java version (default to 1.6.0_21)\n"
   printf "\t-n specify the base output filename rather than derive it from the input .bam files (default derive from input .bams)\n"
   printf "\n\n"
   exit 1
}

# Process the command line arguments
function cmdargs ()
{
        args=`getopt b:f:d:v:t:n:p:j:T: $*`

        set -- $args
        for optarg
        do
        case $optarg in
        	-b)
				INPUT_DIR=$2; shift; shift;
				;;
			-f)
				FASTA=$2; shift; shift
				;;
			-d)
				DBSNP=$2; shift; shift 
				;;
			-t)
				TYPE=$2; shift; shift
				;;
			-v)    
				VER=$2; shift; shift
				;;
			-p)	
				PICVER=$2; shift; shift
				;;
			-j)
				JAVAVER=$2; shift; shift
				;;
         -n)
            bamname="$2.bam"; shift; shift
            ;;
			-T)
				#Better if this is set in your environment!
				TEMP_DIR=$2; shift; shift
				;;
			--)
				shift; break
				;;
			esac
		done
}

function default (){

	if [ -z $INPUT_DIR ] ; then
		usage
	fi

	if [ -z $FASTA ] ; then
		FASTA="/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/gatk/resources/hg19_random.fa"
	fi

	if [[ $DBSNP =~ vcf$ ]]; then
		DBSNP="$DBSNP"
	elif [ -z $DBSNP ]; then 
		DBSNP="/oicr/data/reference/genomes/homo_sapiens_mc/dbSNP/hg19_random/Genomic/dbSNP142/dbsnp142_chr_noINS_noDup.vcf"
	else 
		echo "The dbSNP input file \"$DBSNP\" does not seem a valid .rod or .vcf file (by extension inspection)"
		usage
	fi

	if [ -z $TYPE ] ; then
		TYPE="exome"
	fi
	if [ -z $VER ] ; then 
		VER="1.3.16"
	fi

	if [ -z $PICVER ] ; then
		PICVER="1.90"
	fi

	if [ -z $JAVAVER ] ; then
		JAVAVER="1.6.0_21"
	fi

	# Each cluster node MUST have a /scratch dir mounted to host specific NFS scratch space 
	if [ -z $TEMP_DIR ] ; then
		TMPDIR="/scratch"
	fi

	#load modules
	module load picard/$PICVER
	module load gatk/$VER
	module load java/$JAVAVER
	module load samtools
}

function directories() {
	
	echo "#!/bin/bash"
	echo "# OICR SEQPRODBIO Team oicr_gatk_helper2.sh script - supports GATK 1.3 and above"
	echo ""
	echo ""

	LOG_DIR="gatkLog"
	echo "if [ ! -d $LOG_DIR ] ; then"
	echo "  mkdir $LOG_DIR"
	echo "fi"
	echo 

	# The /scratch dir is not mounted on hn1. No point in doing this check really
	TMP_DIR="./tmp"
	echo "TMP_DIR=$TMP_DIR"
	#echo "if [[ ! -d \$TMP_DIR ]] ; then"
	#echo "  mkdir \$TMP_DIR"
	#echo "fi"
	#echo "if [[ ! -w \$TMP_DIR ]]; then"
	#echo "	echo \"The TMP dir '\$TMP_DIR' is not writable!\""
	#echo "   exit;"
	#echo "fi"
	#echo 
	
	# rob - no slow bas020 scratch for us!
	TMP_DIR="./tmp"
	echo "if [[ ! -d $TMP_DIR ]] ; then"
	echo "  mkdir $TMP_DIR"
	echo "fi"
	echo

	INTERVAL_DIR="interval"
	echo "if [ ! -d $INTERVAL_DIR ] ; then"
	echo "  mkdir $INTERVAL_DIR"
	echo "fi"
	echo

	BAM_DIR="bam_files"
	echo "if [ ! -d $BAM_DIR ] ; then"
	echo "  mkdir $BAM_DIR"
	echo "fi"
	echo

	VCF_DIR="vcf"
	echo "if [ ! -d $VCF_DIR ] ; then"
	echo " mkdir $VCF_DIR"
	echo "fi"
	echo

   METRICS_DIR="metrics"
   echo "if [ ! -d $METRICS_DIR ] ; then"
   echo " mkdir $METRICS_DIR"
   echo "fi"
   echo

}


function GATK(){
#to get all chromosomes from header of bam file -assumes all files in dir are the same reference and order --> 1 per job
name=`ls $INPUT_DIR/*bam | cut -f1 | head -1`
chr=( `samtools view -H $name | grep "^@SQ" | cut -f 2 | grep SN | sed s/SN:// | grep -v "_"` )
SGE_PREFIX="gatk_"
bname=""

BAMINPUT=""
FIRSTTIME="1";

for i in `ls $INPUT_DIR/*.bam`
do
	TRIM=`basename $i | cut -f 1 -d .`
	LENGTH=${#TRIM}
	
	if [ "${FIRSTTIME}" == "1" ]
	then
		bname=${TRIM:0:$(($LENGTH - 1))}
		FIRSTTIME="0";
	fi

	BAMINPUT="${BAMINPUT} -I ${i}"

	bname="${bname}${TRIM:$((${LENGTH} - 1))}"
done

# Use the bamname if it was specified, otherwise use generated name
if [ -z $bamname ];then
   bname="${bname}.bam"
else
   bname=$bamname
fi




for i in `ls $INPUT_DIR/*.bam`
do
#index input bam
	echo
	echo "# index the input bam file"
	echo "qsub -cwd -N ${SGE_PREFIX}${bname}_index -e $LOG_DIR -o $LOG_DIR -b y -l h_vmem=16G \"module load java/$JAVAVER; module load picard/$PICVER; time java -Xmx4g -Djava.io.tmpdir=\${TMP_DIR} -jar \\\${PICARDROOT}/BuildBamIndex.jar INPUT=${i} VALIDATION_STRINGENCY=LENIENT TMP_DIR=./tmp\""
done

#realigner target creator
	echo
	echo "# realigner target creator "	
	for j in ${chr[@]} ; do
		echo -n "module load java/$JAVAVER; module load gatk/${VER}; time java -Xmx16g -Djava.io.tmpdir=./tmp/ -jar \${GATKROOT}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $FASTA -o ./${INTERVAL_DIR}/${bname}.${j}.intervals $BAMINPUT -L ${j} --known $DBSNP -et NO_ET" > ./$SCRIPT_DIR/RealignerTargetCreator_${j}.sh;
		echo "qsub -cwd -N  ${SGE_PREFIX}${bname}_realignerTargetCreator  -hold_jid  ${SGE_PREFIX}${bname}_index -e $LOG_DIR -o $LOG_DIR -l mf=20G -l h_vmem=96G -b y \"bash ./$SCRIPT_DIR/RealignerTargetCreator_${j}.sh \" ";		
	done
	
#indel realigner
	echo "# Indel realigner"
	for j in ${chr[@]} ; do
		echo -n  "module load java/$JAVAVER; module load gatk/${VER}; time java -Xmx20g -Djava.io.tmpdir=./tmp/ -jar \${GATKROOT}/GenomeAnalysisTK.jar $BAMINPUT -R $FASTA -T IndelRealigner -L ${j} -targetIntervals ./${INTERVAL_DIR}/${bname}.${j}.intervals -o ./${BAM_DIR}/${bname}.realigned.${j}.bam -known $DBSNP --maxReadsInMemory 10000000 -et NO_ET" > ./$SCRIPT_DIR/IndelRealigner_${j}.sh
	
		echo "qsub -cwd  -N  ${SGE_PREFIX}${bname}_indelRealigner -hold_jid  ${SGE_PREFIX}${bname}_realignerTargetCreator -e $LOG_DIR -o $LOG_DIR -l h_vmem=96G -b y \"bash ./$SCRIPT_DIR/IndelRealigner_${j}.sh \" ";
		done

#fix mates
	echo
	echo "# fix mate "
	for j in ${chr[@]} ; do
		echo "qsub -cwd  -N ${SGE_PREFIX}${bname}_fixmate -hold_jid ${SGE_PREFIX}${bname}_indelRealigner -e $LOG_DIR -o $LOG_DIR -b y -l h_vmem=16G \" module load java/$JAVAVER; module load picard/$PICVER; time java -jar -Xmx4g -Djava.io.tmpdir=\${TMP_DIR} \\\${PICARDROOT}/FixMateInformation.jar INPUT=./${BAM_DIR}/${bname}.realigned.${j}.bam OUTPUT=./${BAM_DIR}/${bname}.realigned.${j}.fixmate.bam SO=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=./tmp\"";
	done	

# rob - no need to index2 when CREATE_INDEX=true!

#index bam
#TAB 20111118 - workaround to avoid error 'Mate negative strand flag should not be set for unpaired read' or 'MRNM should not be set for unpaired read'
#The fix is to prevent this - either in the aligner - or at worse fix the FixMate step...
#Also, should the .bai files actually be .bam.bai? I guess it depends on which tool is using the bam
#	echo
#	echo "#create index "
#	for j in ${chr[@]} ; do
#		echo "qsub -cwd -N ${SGE_PREFIX}${bname}_index2 -hold_jid ${SGE_PREFIX}${bname}_fixmate -e $LOG_DIR -o $LOG_DIR -b y -l h_vmem=6G \"module load java/$JAVAVER; module load picard/$PICVER; time java -Xmx4g -Djava.io.tmpdir=\${TMP_DIR} -jar \\\${PICARDROOT}/BuildBamIndex.jar INPUT=./${BAM_DIR}/${bname}.realigned.${j}.fixmate.bam VALIDATION_STRINGENCY=SILENT\""
#	done

#count covariates for recalibration
	echo 
	echo "# count covariates "
	echo "R=\"$FASTA\""  > ./$SCRIPT_DIR/CountCovariates.sh
	echo "dbSNP=\"${DBSNP}\""  >> ./$SCRIPT_DIR/CountCovariates.sh
	echo "TMP_DIR=\"${TMP_DIR}\"" >> ./$SCRIPT_DIR/CountCovariates.sh
	echo -n "module load java/$JAVAVER; module load gatk/${VER}; time java -Xmx16g -Djava.io.tmpdir=\${TMP_DIR} -jar \${GATKROOT}/GenomeAnalysisTK.jar -l INFO -R $FASTA -knownSites ${DBSNP} " >> ./$SCRIPT_DIR/CountCovariates.sh
	for j in ${chr[@]} ; do
		echo -n "-I ./${BAM_DIR}/${bname}.realigned.${j}.fixmate.bam " >> ./$SCRIPT_DIR/CountCovariates.sh
	done
	echo "  -T CountCovariates  -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile ./${BAM_DIR}/${bname}.realigned.recal_data.csv -nt 8 -et NO_ET" >> ./$SCRIPT_DIR/CountCovariates.sh
	#TAB 20111130 - noted removed '-pe smp 8' since this is no longer implemented and will block the job forever...
	echo "qsub -cwd -b y -N ${SGE_PREFIX}${bname}_countCovariates -hold_jid ${SGE_PREFIX}${bname}_fixmate -e $LOG_DIR -o $LOG_DIR -l h_vmem=96G  \"bash ./$SCRIPT_DIR/CountCovariates.sh \" ";

#table recalibration 
	echo 
	echo "# table recal"
	for j in ${chr[@]} ; do
		echo "R=\"$FASTA\""  > ./$SCRIPT_DIR/TableRecal_${j}.sh
		echo "dbSNP=\"${DBSNP}\""  >> ./$SCRIPT_DIR/TableRecal_${j}.sh
	   echo "TMP_DIR=\"${TMP_DIR}\"" >> ./$SCRIPT_DIR/TableRecal_${j}.sh
		echo -n "module load java/$JAVAVER; module load gatk/${VER}; time java -Xmx16g -Djava.io.tmpdir=\${TMP_DIR} -jar \${GATKROOT}/GenomeAnalysisTK.jar -T TableRecalibration --preserve_qscores_less_than 5 -l INFO -R $FASTA  ">> ./$SCRIPT_DIR/TableRecal_${j}.sh
		for k in ${chr[@]} ; do 
			echo -n "-I ./${BAM_DIR}/${bname}.realigned.${k}.fixmate.bam " >>./$SCRIPT_DIR/TableRecal_${j}.sh
		done
		echo "--out ./${BAM_DIR}/${bname}.realigned.recal.${j}.bam -recalFile ./${BAM_DIR}/${bname}.realigned.recal_data.csv -L ${j} -et NO_ET" >> ./$SCRIPT_DIR/TableRecal_${j}.sh
	done

	for j in ${chr[@]} ; do 
		echo "qsub -cwd  -N  ${SGE_PREFIX}${bname}_table_recal -hold_jid ${SGE_PREFIX}${bname}_countCovariates -e $LOG_DIR -o $LOG_DIR -b y -l h_vmem=96G \"bash ./$SCRIPT_DIR/TableRecal_${j}.sh \" ";
	done

#index bam
#TAB 20111118 - workaround to avoid error 'Mate negative strand flag should not be set for unpaired read' or 'MRNM should not be set for unpaired read'
#The fix is to prevent this - either in the aligner - or at worse fix the FixMate step...
#Also, should the .bai files actually be .bam.bai? I guess it depends on which tool is using the bam
	echo 
	echo "# index"
	for j in ${chr[@]} ; do
		echo "qsub -cwd -b y -N  ${SGE_PREFIX}${bname}_index3 -hold_jid ${SGE_PREFIX}${bname}_table_recal -e $LOG_DIR -o $LOG_DIR -l h_vmem=16G \" module load java/$JAVAVER; module load picard/$PICVER; time java -Xmx4g -Djava.io.tmpdir=\${TMP_DIR} -jar \\\${PICARDROOT}/BuildBamIndex.jar I=./${BAM_DIR}/${bname}.realigned.recal.${j}.bam VALIDATION_STRINGENCY=SILENT TMP_DIR=./tmp\"";
	done

# TAB 20111129 - Add AnalyzeCovariates step here

#snp calling
# TAB - removed -pe smp 8 from qsub...
	echo
	echo "# unifiedGenotyper"
	for j in ${chr[@]} ; do
        	echo "qsub -cwd -b y -N ${SGE_PREFIX}${bname}_unifiedGenotyper -hold_jid ${SGE_PREFIX}${bname}_index3 -e $LOG_DIR -o $LOG_DIR -l h_vmem=16G \"module load java/$JAVAVER; module load gatk/${VER}; time java -Xmx4g -Djava.io.tmpdir=\${TMP_DIR} -jar \\\${GATKROOT}/GenomeAnalysisTK.jar -R $FASTA -T UnifiedGenotyper -l INFO -I ./${BAM_DIR}/${bname}.realigned.recal.${j}.bam --dbsnp $DBSNP -o $VCF_DIR/${bname}.realigned.recal.bam.snps.raw.${j}.vcf -stand_call_conf 30 -stand_emit_conf 1.0 -metrics $METRICS_DIR/${bname}.${j}.snp.metrics -nt 8 -L ${j} --computeSLOD -et NO_ET\"";
         done

#merge snp output
	echo 
	echo "#merge snps"
	echo "grep \"^#\" ./$VCF_DIR/${bname}.realigned.recal.bam.snps.raw.${chr[0]}.vcf  | grep -v \"##UnifiedGenotyper\" > ${bname}.realigned.recal.bam.snps.raw.vcf "> ./$SCRIPT_DIR/snps_merge.sh 
	for j in ${chr[@]} ; do
		echo "grep \"^chr\" ./$VCF_DIR/${bname}.realigned.recal.bam.snps.raw.${j}.vcf >> ${bname}.realigned.recal.bam.snps.raw.vcf  " >> ./$SCRIPT_DIR/snps_merge.sh 
	done
	echo "qsub -cwd -b y -N ${SGE_PREFIX}${bname}_snpMerge -hold_jid ${SGE_PREFIX}${bname}_unifiedGenotyper -e $LOG_DIR -o $LOG_DIR  \"bash ./$SCRIPT_DIR/snps_merge.sh \"  ";

#indel calling
	echo 
   echo "# unifiedGenotyper indel"
   for j in ${chr[@]} ; do
   echo "qsub -cwd -N ${SGE_PREFIX}${bname}_unifiedGenotyperIndel -hold_jid ${SGE_PREFIX}${bname}_index3 -e $LOG_DIR -o $LOG_DIR -b y -l h_vmem=16G \" module load java/$JAVAVER; module load gatk/${VER}; time java -Xmx4g -Djava.io.tmpdir=\${TMP_DIR}  -jar \\\${GATKROOT}/GenomeAnalysisTK.jar -R $FASTA -T UnifiedGenotyper -glm INDEL -l INFO -I ./${BAM_DIR}/${bname}.realigned.recal.${j}.bam --dbsnp $DBSNP -o ./$VCF_DIR/${bname}.realigned.recal.bam.indels.raw.${j}.vcf -stand_call_conf 30 -stand_emit_conf 1.0 -metrics $METRICS_DIR/${bname}.${j}.indel.metrics -nt 8 -L ${j} --computeSLOD -et NO_ET \"";
   done

#merge indel output
	echo
	echo "#merge indel"
	echo "grep \"^#\" ./$VCF_DIR/${bname}.realigned.recal.bam.indels.raw.${chr[0]}.vcf | grep -v \"##UnifiedGenotyper\" > ${bname}.realigned.recal.bam.indels.raw.vcf " > ./$SCRIPT_DIR/indel_merge.sh
	for j in ${chr[@]} ; do
		echo "grep \"^chr\" ./$VCF_DIR/${bname}.realigned.recal.bam.indels.raw.${j}.vcf >> ${bname}.realigned.recal.bam.indels.raw.vcf " >> ./$SCRIPT_DIR/indel_merge.sh
	done
	echo "qsub -cwd -b y -N ${SGE_PREFIX}${bname}_indelMerge -hold_jid  ${SGE_PREFIX}${bname}_unifiedGenotyperIndel -e $LOG_DIR -o $LOG_DIR  \"bash ./$SCRIPT_DIR/indel_merge.sh \" ";


 #filter indel	 
 # TAB Checked to comply with http://www.broadinstitute.org/gsa/wiki/index.php/Best_Practice_Variant_Detection_with_the_GATK_v2
 # TAB This requires updating, and implementation of variant quality score recalibration
	echo
	echo "# filter indels "
	for j in ${chr[@]} ; do
		echo "module load java/$JAVAVER; module load gatk/${VER}; time java -Xmx16g -Djava.io.tmpdir=${TMP_DIR}  -jar \${GATKROOT}/GenomeAnalysisTK.jar -T VariantFiltration -R $FASTA --variant ./$VCF_DIR/${bname}.realigned.recal.bam.indels.raw.${j}.vcf -o ./$VCF_DIR/${bname}.realigned.recal.bam.indels.raw.filtered.${j}.vcf -L ${j} --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\"  --filterExpression \"FS > 60.0\" --filterName \"StrandBiasFishers\" --filterExpression \"SB > -10.0\" --filterName \"StrandBias\" --filterExpression \"DP < 5 \" --filterName \"LowCoverage\" --filterExpression \"QUAL < 30.0\"  --filterName \"VeryLowQual\" --filterExpression \"QUAL > 30.0 && QUAL < 50.0 \"  --filterName \"LowQual\" -et NO_ET " > ./$SCRIPT_DIR/indelFilter_${j}.sh;
	done
	for j in ${chr[@]} ; do
		echo "qsub -cwd  -N ${SGE_PREFIX}${bname}_indelFiltering -hold_jid ${SGE_PREFIX}${bname}_unifiedGenotyperIndel -e $LOG_DIR -o $LOG_DIR -l h_vmem=96G -b y \"bash ./$SCRIPT_DIR/indelFilter_${j}.sh \""
	done

#merge filtered indel
	echo
	echo "#merge filtered indel"
	echo "grep \"^#\" ./$VCF_DIR/${bname}.realigned.recal.bam.indels.raw.filtered.${chr[0]}.vcf  | grep -v \"##UnifiedGenotyper\" > ${bname}.realigned.recal.bam.indels.raw.filtered.vcf " > ./$SCRIPT_DIR/indel_filter_merge.sh
	for j in ${chr[@]} ; do
		echo "grep \"^chr\" ./$VCF_DIR/${bname}.realigned.recal.bam.indels.raw.filtered.${j}.vcf >> ${bname}.realigned.recal.bam.indels.raw.filtered.vcf " >> ./$SCRIPT_DIR/indel_filter_merge.sh
	done
	echo "qsub -cwd -b y -N ${SGE_PREFIX}${bname}_indelFilterMerge -hold_jid ${SGE_PREFIX}${bname}_indelFiltering -e $LOG_DIR -o $LOG_DIR  \"bash ./$SCRIPT_DIR/indel_filter_merge.sh\""

#filter snps
# TAB Checked to comply with http://www.broadinstitute.org/gsa/wiki/index.php/Best_Practice_Variant_Detection_with_the_GATK_v2
# For wgs data, increased the depth cutoff from DP > 100 to DP > 500 - we have some deep wgs's
# TAB This requires updating, and implementation of variant quality score recalibration
	echo
	echo "#filter snps "
	if [ $TYPE = "exome" ]; then
		for j in ${chr[@]} ; do
		echo "module load gatk/${VER}; time java -Xmx16g -Djava.io.tmpdir=${TMP_DIR} -jar \${GATKROOT}/GenomeAnalysisTK.jar -T VariantFiltration -R $FASTA  --variant ./$VCF_DIR/${bname}.realigned.recal.bam.snps.raw.${j}.vcf -o ./$VCF_DIR/${bname}.realigned.recal.bam.snps.raw.filtered.${j}.vcf -L ${j} --clusterWindowSize 10 --clusterSize 3 --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\" --mask ./$VCF_DIR/${bname}.realigned.recal.bam.indels.raw.${j}.vcf --maskName InDel --filterExpression \"DP < 5 \" --filterName \"LowCoverage\" --filterExpression \"QUAL < 30.0 \" --filterName \"VeryLowQual\" --filterExpression \"QUAL > 30.0 && QUAL < 50.0 \" --filterName \"LowQual\" --filterExpression \"QD < 2.0 \" --filterName \"LowQD\" --filterExpression \"FS > 60.0 \" --filterName \"StrandBiasFishers\" --filterExpression \"SB > -10.0\" --filterName \"StrandBias\" -et NO_ET " > ./$SCRIPT_DIR/snpFilter_${j}.sh
	done
	fi
	if [ $TYPE = "wgs" ]; then
		for j in ${chr[@]} ; do
		echo "module load java/$JAVAVER; module load gatk/${VER}; time java -Xmx16g -Djava.io.tmpdir=${TMP_DIR} -jar \${GATKROOT}/GenomeAnalysisTK.jar -T VariantFiltration -R $FASTA  --variant ./$VCF_DIR/${bname}.realigned.recal.bam.snps.raw.${j}.vcf -o ./$VCF_DIR/${bname}.realigned.recal.bam.snps.raw.filtered.${j}.vcf -L ${j} --clusterWindowSize 10 --clusterSize 3 --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\" --mask ./$VCF_DIR/${bname}.realigned.recal.bam.indels.raw.${j}.vcf --maskName InDel --filterExpression \"DP < 5 \" --filterName \"LowCoverage\" --filterExpression \"QUAL < 30.0 \" --filterName \"VeryLowQual\" --filterExpression \"QUAL > 30.0 && QUAL < 50.0 \" --filterName \"LowQual\" --filterExpression \"QD < 2.0 \" --filterName \"LowQD\" --filterExpression \"FS > 60.0 \" --filterName \"StrandBiasFishers\" --filterExpression \"SB > -0.1\" --filterName \"StrandBias\" --filterExpression \"DP > 500 \" --filterName \"ExcessiveDepth\" -et NO_ET " > ./$SCRIPT_DIR/snpFilter_${j}.sh;
	done
	fi
	for j in ${chr[@]} ; do
		echo "qsub -cwd  -N ${SGE_PREFIX}${bname}_snpFiltering -hold_jid ${SGE_PREFIX}${bname}_unifiedGenotyper,${SGE_PREFIX}${bname}_unifiedGenotyperIndel -e $LOG_DIR -o $LOG_DIR -l h_vmem=96G -b y \"bash ./$SCRIPT_DIR/snpFilter_${j}.sh \" ";
	done
	
#merge filtered snps
        echo
	echo "#merge filtered snps"
	echo "grep \"^#\" ./$VCF_DIR/${bname}.realigned.recal.bam.snps.raw.filtered.${chr[0]}.vcf  | grep -v \"##UnifiedGenotyper\" > ${bname}.realigned.recal.bam.snps.raw.filtered.vcf " > ./$SCRIPT_DIR/snps_filter_merge.sh
	for j in ${chr[@]} ; do
		echo "grep \"^chr\" ./$VCF_DIR/${bname}.realigned.recal.bam.snps.raw.filtered.${j}.vcf >> ${bname}.realigned.recal.bam.snps.raw.filtered.vcf " >> ./$SCRIPT_DIR/snps_filter_merge.sh
	done
	echo "qsub -cwd -b y -N ${SGE_PREFIX}${bname}_snpsFilterMerge -hold_jid ${SGE_PREFIX}${bname}_snpFiltering -e $LOG_DIR -o $LOG_DIR  \"bash ./$SCRIPT_DIR/snps_filter_merge.sh\" ";


#annotate vcf
echo
echo "#annotate vcf"
rm -rf resources
ln -s /oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/gatk/resources

R="hg19_random.fa"
dbSNP="dbSNP142_chr.vcf"
g1000="1000g_20100804_chr.vcf"
HapMap="HapMap_r27_nr_hg19_chr.vcf"
nmblExomeTargetVCF="/.mounts/labs/spbprojects/PCSI_2012/exomeTargets/nimblegen_hg19_2.1M_Human_Exome.gatk.vcf.gz"
agilentExomeTargetVCF="/.mounts/labs/spbprojects/PCSI_2012/exomeTargets/agilent_icgc_sanger.exons.hg19.sorted.fixed.gatk.vcf.gz"
inputSNPVCF="${bname}.realigned.recal.bam.snps.raw.filtered.vcf"
inputINDELVCF="${bname}.realigned.recal.bam.indels.raw.filtered.vcf"
PREFIX=${inputVCF%%.*}
VCFSNPOUTPREFIX=${inputSNPVCF%%.vcf}
VCFINDELOUTPREFIX=${inputINDELVCF%%.vcf}

# Generate Variant Annotator.
echo "module load gatk/${VER}; \${JAVA_HOME}/bin/java -Xmx4g -Djava.io.tmpdir=${TMP_DIR} -jar \${GATKROOT}/GenomeAnalysisTK.jar -et NO_ET -l INFO -R resources/$R -o ${VCFSNPOUTPREFIX}.annotated.vcf --dbsnp resources/$dbSNP -resource:snp resources/$dbSNP -E snp.dbSNPBuildID -E snp.G5 -E snp.G5A -E snp.GMAF -comp:1KG_CEU resources/${g1000} -comp:HapMap resources/${HapMap} -comp:Nmblhg19 $nmblExomeTargetVCF -comp:AgilentICGChg19 $agilentExomeTargetVCF --variant ${inputSNPVCF} -T VariantAnnotator" > ./$SCRIPT_DIR/annotateSNP.sh
echo "module load gatk/${VER}; \${JAVA_HOME}/bin/java -Xmx4g -Djava.io.tmpdir=${TMP_DIR} -jar \${GATKROOT}/GenomeAnalysisTK.jar -et NO_ET -l INFO -R resources/$R -o ${VCFINDELOUTPREFIX}.annotated.vcf --dbsnp resources/$dbSNP -resource:snp resources/$dbSNP -E snp.dbSNPBuildID -E snp.G5 -E snp.G5A -E snp.GMAF -comp:1KG_CEU resources/${g1000} -comp:HapMap resources/${HapMap} -comp:Nmblhg19 $nmblExomeTargetVCF -comp:AgilentICGChg19 $agilentExomeTargetVCF --variant ${inputINDELVCF} -T VariantAnnotator" > ./$SCRIPT_DIR/annotateINDEL.sh

echo "qsub -cwd -N ${SGE_PREFIX}${bname}_VariantAnnotatorSNP -hold_jid ${SGE_PREFIX}${bname}_snpsFilterMerge -e gatkLog -o gatkLog -l h_vmem=16G -b y \"bash ./$SCRIPT_DIR/annotateSNP.sh\"";
echo "qsub -cwd -N ${SGE_PREFIX}${bname}_VariantAnnotatorINDEL -hold_jid ${SGE_PREFIX}${bname}_indelFilterMerge -e gatkLog -o gatkLog -l h_vmem=16G -b y \"bash ./$SCRIPT_DIR/annotateINDEL.sh\"";

# delete temporary bam files

echo
echo "#delete temporary bam files"

#echo "qsub -cwd -b y -N ${SGE_PREFIX}${bname}_clean_tmp_files -hold_jid ${SGE_PREFIX}${bname}_table_recal -e gatkLog -o gatkLog \"rm bam_files/*.bam.realigned.chr*.ba*\""
echo "qsub -cwd -b y -N ${SGE_PREFIX}${bname}_clean_tmp_files -hold_jid ${SGE_PREFIX}${bname}_VariantAnnotatorSNP,${SGE_PREFIX}${bname}_VariantAnnotatorINDEL -e gatkLog -o gatkLog \"rm bam_files/*.ba?\""


# touch a file to indicate gatk is done
echo
echo "# touch a file to indicate gatk is done"
echo "qsub -cwd -b y -N ${SGE_PREFIX}${bname}_done -hold_jid ${SGE_PREFIX}${bname}_VariantAnnotatorSNP,${SGE_PREFIX}${bname}_VariantAnnotatorINDEL -e gatkLog -o gatkLog \"touch gatk_done\""
}

function callableLociWalker
{
# run callable loci walker
	echo
	echo "# run callable loci walker"
	# build recal bam files list
	BAMFILES=""
	for j in  ${chr[@]}
	do
		BAMFILES="${BAMFILES} -I ./${BAM_DIR}/${bname}.realigned.recal.${j}.bam"
	done

	if [ "${TYPE}" == "exome" ]
	then
		MINDEPTH=8
	elif [ "${TYPE}" == "wgs" ]
	then
		MINDEPTH=5
	fi

	echo "module load java/${JAVAVER}; module load gatk/${VER}; java -Xmx4g -Djava.io.tmpdir=${TMP_DIR} -jar \${GATKROOT}/GenomeAnalysisTK.jar -l INFO -R ${FASTA} -o ${bname}.realigned.recal.callable -format BED -maxDepth 5000 -minDepth ${MINDEPTH} -summary ${bname}.realigned.recal.callable.summary -T CallableLoci ${BAMFILES}" > ./${SCRIPT_DIR}/callable_loci.sh
	chmod 764 ./${SCRIPT_DIR}/callable_loci.sh

	echo "qsub -cwd -b y -N ${SGE_PREFIX}${bname}_callableLociWalker -hold_jid ${SGE_PREFIX}${bname}_index3 -e $LOG_DIR -o $LOG_DIR -l mf=4G ./${SCRIPT_DIR}/callable_loci.sh"

	echo "qsub -cwd -b y -N ${SGE_PREFIX}${bname}_grepForCallable -hold_jid ${SGE_PREFIX}${bname}_callableLociWalker -e $LOG_DIR -o $LOG_DIR \"grep CALLABLE ${bname}.realigned.recal.callable | sed 's/\\ /\\t/g' > ${bname}.realigned.recal.callable.bed\""

	echo "qsub -cwd -b y -N ${SGE_PREFIX}${bname}_callableBedToVcf -hold_jid ${SGE_PREFIX}${bname}_grepForCallable -e $LOG_DIR -o $LOG_DIR -l mf=48G /home/rdenroche/GATK/bedToVcf ${bname}.realigned.recal.callable.bed"

	echo "qsub -cwd -b y -N ${SGE_PREFIX}${bname}_gzipAndIndexCallableVcf -hold_jid ${SGE_PREFIX}${bname}_callableBedToVcf -e $LOG_DIR -o $LOG_DIR \"module load tabix; bgzip ${bname}.realigned.recal.callable.vcf; tabix -f -p vcf ${bname}.realigned.recal.callable.vcf.gz\""

}	

	
SCRIPT_DIR="scripts"
	if [ ! -d $SCRIPT_DIR ] ; then
 	       mkdir $SCRIPT_DIR
        fi

cmdargs $PROGARGS
default
directories
GATK
#callableLociWalker 	# not sure if this still works as a standalone function - some variables may be out of scope -rob

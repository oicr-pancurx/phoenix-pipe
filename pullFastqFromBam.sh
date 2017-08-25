#!/bin/bash

BAM=$1;
R1="_R1.fastq";
R2="_R2.fastq";

OUTNAME=`echo $BAM | sed 's/\.bam$//'`;

qsub -cwd -b y -N b2f`basename $BAM` -l h_vmem=8g -e $OUTNAME.b2f.log -o $OUTNAME.b2f.log "module load samtools; samtools view $BAM | /u/rdenroche/pullFastqFromBam.pl $BAM$R1 $BAM$R2"


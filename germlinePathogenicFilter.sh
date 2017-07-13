#!/bin/sh

VCF=$1;
name=`basename $VCF | cut -f 1 -d "."`;

mkdir annovar
(cut -f 1-5,9-11 $VCF | grep "^#CHR"; grep -wf /u/aconnor/lists/hpcs_genes $VCF | cut -f 1-5,9-11) > annovar/$name.truncated.vcf

module load annovar/2014-07-15
convert2annovar.pl -format vcf4old annovar/$name.truncated.vcf -outfile annovar/$name.input.vcf -includeinfo

table_annovar.pl annovar/$name.input.vcf /oicr/local/analysis/sw/annovar/humandb/ -buildver hg19 -remove -protocol refGene,avsnp142,popfreq_all_20150413,popfreq_max_20150413,ljb26_all,clinvar_20150330 -operation g,f,f,f,f,f -otherinfo -nastring NA -out annovar/$name.output.vcf

awk '{ if (($13 < 0.01) && ($19 < 0.01) && ($27 < 0.01)) print $0}' annovar/$name.output.vcf.hg19_multianno.txt > annovar/$name.rare.vcf
grep -wf /u/aconnor/lists/hpcs_genes annovar/$name.output.vcf.hg19_multianno.txt > annovar/$name.genes.vcf
grep -wf /u/aconnor/lists/non_silent_variants annovar/$name.output.vcf.hg19_multianno.txt > annovar/$name.nonsil.vcf

grep -wf /u/aconnor/lists/hpcs_genes annovar/$name.rare.vcf | grep -wf /u/aconnor/lists/non_silent_variants > annovar/$name.genes_nonsil_rare.vcf

grep pathogenic annovar/$name.genes_nonsil_rare.vcf | grep -v non-pathogenic > $name.pathogenic.vcf
grep -wf /u/aconnor/lists/inactivating_variants annovar/$name.genes_nonsil_rare.vcf | grep -v pathogenic > $name.likely_pathogenic.vcf




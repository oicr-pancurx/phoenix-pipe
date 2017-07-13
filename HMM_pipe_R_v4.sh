

PRI=`basename $1| cut -f 1 -d .` 
REF=`basename $2| cut -f 1 -d .` 

OUTDIR=$3

PTYPE=`basename $1| cut -f 3 -d _ |cut -f 1 -d .`
RTYPE=`basename $2| cut -f 3 -d _ |cut -f 1 -d .`

SAM=`basename $1| cut -f 1-2 -d _`

#echo "PRI----$PRI"
#echo "REF----$REF"
#echo "PTYPE----$PTYPE"
#echo "RTYPE----$RTYPE"
#echo "SAM----$SAM"




echo ".libPaths(new=\"/.mounts/labs/PCSI/R/3.3/\")" 
echo "library(HMMcopy)" 



echo "tum_uncorrected_reads <- wigsToRangedData(\"$1\",\"/oicr/data/reference/genomes/homo_sapiens_mc/HMMcopy/gc_hg19.wig\",\"/oicr/data/reference/genomes/homo_sapiens_mc/HMMcopy/hg19_map.wig\")"
echo "norm_uncorrected_reads <- wigsToRangedData(\"$2\",\"/oicr/data/reference/genomes/homo_sapiens_mc/HMMcopy/gc_hg19.wig\",\"/oicr/data/reference/genomes/homo_sapiens_mc/HMMcopy/hg19_map.wig\")"



echo "tum_corrected_copy <- correctReadcount(tum_uncorrected_reads)"
echo "somatic_corrected_copy <- correctReadcount(tum_uncorrected_reads)"
echo "norm_corrected_copy <- correctReadcount(norm_uncorrected_reads)"

echo "somatic_corrected_copy\$copy <- tum_corrected_copy\$copy - norm_corrected_copy\$copy"

echo "param <- HMMsegment(tum_corrected_copy, getparam = TRUE)"
#echo "param\$mu <- log(c(1, 1.4, 2, 2.7, 3, 4.5) / 2, 2)"
#echo "param\$m <- param\$mu"
echo "param\$e <- 0.999999999999999"
echo "param\$strength <- 1e+30"





#echo "segmented_copy_tum <- HMMsegment(tum_corrected_copy)"
echo "segmented_copy_tum <- HMMsegment(tum_corrected_copy, param)"
#echo "segmented_copy_norm <- HMMsegment(norm_corrected_copy)"
echo "segmented_copy_somatic <- HMMsegment(somatic_corrected_copy, param)"
echo "segmented_copy_norm <- HMMsegment(norm_corrected_copy, param)"
#echo "segmented_copy_somatic <- HMMsegment(somatic_corrected_copy)"

####### Output seg files
echo "rangedDataToSeg(tum_corrected_copy, file = \"$OUTDIR/$SAM.$PTYPE.segments_tumor.seg\")"
echo "rangedDataToSeg(norm_corrected_copy, file = \"$OUTDIR/$SAM.$RTYPE.segments_norm.seg\")"
echo "rangedDataToSeg(somatic_corrected_copy, file = \"$OUTDIR/$SAM.$PTYPE$RTYPE.segments_somatic.seg\")"
########################
####### Output segments
echo "chr_segments_tum <- segmented_copy_tum\$segs[grep(\"_\", segmented_copy_tum\$segs\$chr,invert=TRUE),]"
echo "write.table(chr_segments_tum,file=\"$OUTDIR/$SAM$PTYPE.cnv_tum_segments\",quote=FALSE,sep=\"\t\",row.names=FALSE)"


echo "chr_segments_norm <- segmented_copy_norm\$segs[grep(\"_\", segmented_copy_norm\$segs\$chr,invert=TRUE),]"
echo "write.table(chr_segments_norm,file=\"$OUTDIR/$SAM$RTYPE.cnv_norm_segments\",quote=FALSE,sep=\"\t\",row.names=FALSE)"


echo "chr_segments_somatic <- segmented_copy_somatic\$segs[grep(\"_\", segmented_copy_somatic\$segs\$chr,invert=TRUE),]"
echo "write.table(chr_segments_somatic,file=\"$OUTDIR/$SAM$PTYPE$RTYPE.cnv_somatic_segments\",quote=FALSE,sep=\"\t\",row.names=FALSE)"
#########################

echo "plotchr <- function (correctOutput, segmentOutput, chr = space(correctOutput)[1]," 
echo "...)" 
echo "{"
echo "if (is.null(segmentOutput\$segs)) {"
echo "         warning(\"Processed segments now found, automatically processing\")"
echo "         segmentOutput\$segs <- processSegments(segments\$segs," 
echo "             space(correctOutput), start(correctOutput), end(correctOutput)," 
echo "             correctOutput\$copy)"
echo "     }"
echo "     segs <- segmentOutput\$segs"
echo "     correctOutput\$state <- segmentOutput\$state"
echo "     cols <- stateCols()"
echo "     range <- quantile(correctOutput\$copy, na.rm = TRUE, prob = c(0.01," 
echo "         0.99))"
echo "     a <- correctOutput[as.character(chr)]"
echo "     b <- segs[segs\$chr == chr, ]"
echo "     plot(start(a), a\$copy, col = cols[as.numeric(as.character(a\$state))]," 
echo "         ylim = c(-4,4), ...)"
echo "     for (k in 1:nrow(b)) {"
echo "         lines(c(b\$start[k], b\$end[k]), rep(b\$median[k], 2), lwd = 3," 
echo "             col = \"green\")"
echo "     }"
echo " }"


echo "pngplotchr <- function (correctOutput, segmentOutput, chr = space(correctOutput)[1]," 
echo "...)" 
echo "{"
echo "if (is.null(segmentOutput\$segs)) {"
echo "         warning(\"Processed segments now found, automatically processing\")"
echo "         segmentOutput\$segs <- processSegments(segments\$segs," 
echo "             space(correctOutput), start(correctOutput), end(correctOutput)," 
echo "             correctOutput\$copy)"
echo "     }"
echo "     segs <- segmentOutput\$segs"
echo "     correctOutput\$state <- segmentOutput\$state"
echo "     cols <- stateCols()"
echo "     range <- quantile(correctOutput\$copy, na.rm = TRUE, prob = c(0.01," 
echo "         0.99))"
echo "     a <- correctOutput[as.character(chr)]"
echo "     b <- segs[segs\$chr == chr, ]"
echo "     plot(start(a), a\$copy, col = cols[as.numeric(as.character(a\$state))]," 
echo "         xaxt='n',yaxt='n',ylim = c(-3,3), ...)"
#echo "     for (k in 1:nrow(b)) {"
#echo "         lines(c(b\$start[k], b\$end[k]), rep(b\$median[k], 2), lwd = 3," 
#echo "             col = \"green\")"
#echo "     }"
echo " }"









####Output#####
#PDF of CNV plots#
echo "pdf(\"$OUTDIR/$SAM.$PTYPE.CNV_Tumor.pdf\",10,75)"
echo " par(mfrow=c(24,1))"

echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr1\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 1 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr2\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 2 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr3\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 3 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr4\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 4 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr5\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 5 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr6\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 6 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr7\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 7 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr8\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 8 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr9\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 9 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr10\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 10 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr11\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 11 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr12\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 12 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr13\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 13 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr14\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 14 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr15\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 15 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr16\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 16 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr17\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 17 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr18\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 18 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr19\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 19 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr20\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 20 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr21\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 21 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chr22\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 22 Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chrX\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome X Position\")"
echo "plotchr(tum_corrected_copy,segmented_copy_tum,chr=\"chrY\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome Y Position\")"

echo "dev.off()"


echo "pdf(\"$OUTDIR/$SAM.$RTYPE.CNV_Normal.pdf\",10,75)"
echo " par(mfrow=c(24,1))"

echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr1\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 1 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr2\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 2 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr3\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 3 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr4\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 4 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr5\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 5 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr6\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 6 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr7\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 7 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr8\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 8 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr9\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 9 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr10\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 10 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr11\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 11 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr12\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 12 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr13\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 13 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr14\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 14 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr15\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 15 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr16\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 16 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr17\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 17 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr18\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 18 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr19\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 19 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr20\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 20 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr21\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 21 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chr22\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 22 Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chrX\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome X Position\")"
echo "plotchr(norm_corrected_copy,segmented_copy_norm,chr=\"chrY\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome Y Position\")"

echo "dev.off()"


echo "pdf(\"$OUTDIR/$SAM.$PTYPE$RTYPE.CNV_Somatic.pdf\",10,75)"
echo " par(mfrow=c(24,1))"

echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr1\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 1 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr2\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 2 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr3\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 3 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr4\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 4 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr5\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 5 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr6\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 6 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr7\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 7 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr8\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 8 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr9\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome 9 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr10\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 10 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr11\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 11 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr12\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 12 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr13\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 13 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr14\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 14 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr15\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 15 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr16\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 16 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr17\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 17 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr18\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 18 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr19\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 19 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr20\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 20 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr21\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 21 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr22\", pch= \".\",ylab = \"Tumour Copy Number\", xlab = \"Chromosome 22 Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chrX\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome X Position\")"
echo "plotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chrY\", pch= \".\",ylab =  \"Tumour Copy Number\", xlab = \"Chromosome Y Position\")"

echo "dev.off()"

echo "png(\"$OUTDIR/$SAM.genomeCNV.png\", width=1200,height=600)"
echo "par(mfrow=c(1,24),mar=c(4,0,4,0),mgp=c(0,1,0),oma=c(1,1,1,1))"

echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr1\", pch= \".\",xlab = \"chr1\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr2\", pch= \".\",xlab = \"chr2\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr3\", pch= \".\",xlab = \"chr3\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr4\", pch= \".\",xlab = \"chr4\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr5\", pch= \".\",xlab = \"chr5\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr6\", pch= \".\",xlab = \"chr6\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr7\", pch= \".\",xlab = \"chr7\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr8\", pch= \".\",xlab = \"chr8\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr9\", pch= \".\",xlab = \"chr9\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr10\", pch= \".\",xlab = \"chr10\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr11\", pch= \".\",xlab = \"chr11\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr12\", pch= \".\",xlab = \"chr12\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr13\", pch= \".\",xlab = \"chr13\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr14\", pch= \".\",xlab = \"chr14\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr15\", pch= \".\",xlab = \"chr15\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr16\", pch= \".\",xlab = \"chr16\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr17\", pch= \".\",xlab = \"chr17\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr18\", pch= \".\",xlab = \"chr18\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr19\", pch= \".\",xlab = \"chr19\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr20\", pch= \".\",xlab = \"chr20\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr21\", pch= \".\",xlab = \"chr21\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chr22\", pch= \".\",xlab = \"chr22\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chrX\", pch= \".\",xlab = \"chrX\")"
echo "pngplotchr(somatic_corrected_copy,segmented_copy_somatic,chr=\"chrY\", pch= \".\",xlab = \"chrY\")"

echo "dev.off()"






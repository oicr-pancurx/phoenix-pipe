#### chromothripsis.R ##############################################################################
# script used check from evidence of chromothripsis.

### PREAMBLE #######################################################################################
# load library
library(getopt);
options(stringsAsFactors=FALSE);

### FUNCTIONS ######################################################################################
# usage
# function to return usage
usage <- function(){
	usage.text <-
	'\nchromothripsis.R
		-h help
		-t sv (required)
		-n copy number (required)
		-s seg file (required)
		-a AR file (required)
		-b bed file of genome size (defaults to /oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/bedfiles/hg19_random.genome.sizes.bed)
		-o output file prefix (ie sample ID is recommended so name is appended to each file). defaults to ""
	EXAMPLES:
		Rscript chromothripsis.R -t ASHPC_0008_Pa_P.annotatedSV.tsv -n ASHPC_0008_Pa_P_corrected.segments -o ASHPC0008\n\n';
	return(usage.text);
	}

### OBTAIN COMMAND LINE ARGUMENTS ##################################################################
params <- matrix(
	c(
		'help', 'h', 0, "logical",
		'sv', 't', 1, "character",
		'copynumber', 'n', 1, "character",
		'seg', 's', 1, "character",
		'AR', 'a', 1, "character",
		'bed', 'b', 2, "character",
		'output', 'o', 2, "character"
		),
	ncol = 4,
	byrow = TRUE
	);

opt = getopt(params);

# verified required argments and sets defaults
if (!is.null(opt$help) || is.null(opt$sv) || is.null(opt$copynumber) || is.null(opt$seg) || is.null(opt$AR)){
	cat(usage());
	q(status=1)
	}

if(is.null(opt$bed)){
	opt$bed <- "/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/bedfiles/hg19_random.genome.sizes.bed";
	}

### PREAMBLE #######################################################################################
# load library for stats
library(EMT);
library(adehabitatLT);
library(tseries);
library(lattice);
library(latticeExtra);
library(ReorderCluster);
library(gridExtra);
library(grid);
library(png);
library(stringr);
library(chromAL);

# general parameters
date.value <- Sys.Date();

# get file prefix
file.prefix <- return.file.prefix(
	stored.date = date.value,
	output.prefix = opt$output
	);

# store chr info
chr <- paste("chr",c(1:22, 'X','Y'), sep="");
output.cutoff <- 240;

### READ IN DATA ###################################################################################
print(paste("Reading in size file: ", opt$bed, sep = ""));
size <- read.table(
	file = opt$bed,
	as.is = TRUE,
	sep = "\t",
	header = FALSE
	);
colnames(size) <- c('chr', 'start', 'end');

print(paste("Reading in structural variant file: ", opt$sv, sep = ""));
sv.all <- read.table(
	file = opt$sv,
	sep = "\t",
	as.is = TRUE,
	header = TRUE,
	comment = "",
	flush = TRUE
	);
sv.all <- check.duplicate.sv(sv.all);

print(paste("Reading in copy number file: ", opt$copynumber, sep = ""));
cnv <- read.table(
	file = opt$copynumber,
	as.is = TRUE,
	header = TRUE
	);
cnv <- cnv[!is.na(cnv$median),];

### JOIN CLASSIFICATION ############################################################################
print("Classifying joins");
sv <- sv.all;
sv$chrom1 <- factor(sv$chrom1, levels = size$chr, ordered = TRUE);
sv$chrom2 <- factor(sv$chrom2, levels = size$chr, ordered = TRUE);

sv <- classify.joins(sv.data = sv);

write.table(
	x = sv,
	file = paste(file.prefix, '_sv_classified.txt', sep = ""),
	sep = "\t",
	quote = FALSE,
	col.names = TRUE,
	row.names = FALSE
	);

# removing unnecessary columns
sv <- sv[,c('chrom1', 'pos1', 'chrom2', 'pos2', 'type', 'delly_CT', 'crest_strand1', 'crest_strand2', 'class')];
sv <- unique(sv);

### MERGE SV AND CNV DATASETS ######################################################################
# will remove closest cnv within 1kb of sv
print("Merging sv and cnv data");
cnv <- cnv[cnv$chr %in% chr,];

all.breaks <- combine.sv.cnv.breakspoints(
	sv.data = sv,
	cnv.data = cnv,
	chromosomes = size[size$chr %in% chr,]
	);

all.breaks <- filter.duplicates(
	merged.breakpoints = all.breaks,
	distance = 1000
	);

breaks <- all.breaks[!all.breaks$dup,];

write.table(
	x = breaks[c('chr', 'pos', 'type')],
	file = paste(file.prefix, "_cnv_sv_breakpoints.txt", sep=""),
	sep = "\t",
	row.names = FALSE,
	quote = FALSE
	);

### SUMMARIZE SV AND CNV DATA ######################################################################
counts.all <- summarize.sv.cnv.counts(sv.all, cnv, size[size$chr %in% chr,]);
counts <- counts.all[c('SV', 'CNV')];
rownames(counts) <- counts.all$chr;

### SIGNIFICANT SV AND CNV COUNTS ##################################################################
cutoff.cnv.sv <- 8;

print("Plotting cnv vs sv");
significant.values.cnv.sv <- counts[counts$CNV >= cutoff.cnv.sv | counts$SV >= cutoff.cnv.sv,];
significant.values.cnv.sv$sv.pos <- significant.values.cnv.sv$SV + 1.5;
if(length(duplicated(significant.values.cnv.sv$SV)) > 0){
	significant.values.cnv.sv[duplicated(significant.values.cnv.sv$SV),]$sv.pos <- significant.values.cnv.sv[duplicated(significant.values.cnv.sv$SV),]$SV - 1;
	}
significant.values.cnv.sv$cnv.pos <- significant.values.cnv.sv$CNV + 1.5;
if(length(duplicated(significant.values.cnv.sv$CNV)) > 0){
	significant.values.cnv.sv[duplicated(significant.values.cnv.sv$CNV),]$cnv.pos <- significant.values.cnv.sv[duplicated(significant.values.cnv.sv$CNV),]$CNV - 1;
	}

scatterplot.sv.cnv.object <- create.sv.cnv.scatterplot(
	equation.to.use = counts$SV ~ counts$CNV,
	data = counts,
	sig.hits.x.coord = as.vector(significant.values.cnv.sv$cnv.pos),
	sig.hits.y.coord = as.vector(significant.values.cnv.sv$sv.pos),
	sig.hits.labels = as.vector(substr(rownames(significant.values.cnv.sv), 4, 5)),
	xlim = c(-1, round(max(c(counts$SV, counts$CNV)))+5),
	ylim = c(-1, round(max(c(counts$SV, counts$CNV)))+5),
	main = opt$output
	);

output.plot(
	trellis.objects = list(scatterplot.sv.cnv.object),
	filename = paste(file.prefix, "_cnv_vs_sv.tiff", sep="")
	);

significant.values.cnv.sv <- counts[counts$CNV >= cutoff.cnv.sv & counts$SV >= cutoff.cnv.sv,];
write.table(
	counts,
	file = paste(file.prefix, "_cnv_vs_sv.txt", sep=""),
	sep = "\t",
	col.names = NA,
	quote = FALSE
	);

### RUN HYPERGENOMETIC TEST ########################################################################
# set adjusted pvalue threshold
adjusted.p.threshold <- 0.1;

summary.all.break.counts <- as.data.frame(table(breaks$chr));
summary.all.break.counts <- summary.all.break.counts[summary.all.break.counts$Var1 %in% chr,];
summary.all.break.counts <- merge(
	x = size[size$chr %in% summary.all.break.counts$Var1,],
	y = summary.all.break.counts,
	by.x = 'chr',
	by.y = 'Var1'
	);

summary.all.break.counts$chr <- factor(x = summary.all.break.counts$chr, levels = chr, ordered = TRUE);
summary.all.break.counts <- summary.all.break.counts[order(summary.all.break.counts$chr),];

hyper.test.results.SV.CNV <- run.hyper.test(
	count.data = summary.all.break.counts[1:4],
	type = 'full'
	);

write.table(
	hyper.test.results.SV.CNV,
	file = paste(file.prefix, "hyper_test_results_SV_CNV.txt", sep="_"),
	quote = FALSE,
	col.names = NA,
	sep = "\t"
	);

hyper.test.results.SV.CNV.object <- create.clustering.plot(
	data = hyper.test.results.SV.CNV['neg.log10.adjusted.pvalue'],
	cutoff = adjusted.p.threshold,
	main = "SV & CNV Count Distribution"
	);

output.plot(
	trellis.objects = list(hyper.test.results.SV.CNV.object),
	filename = paste(file.prefix, "hyper_test_results_SV_CNV.tiff", sep="_"),
	width = 12
	);

### DISTRIBUTION PLOT ##############################################################################
print("Calculating distribution...");
summary.table <- return.dist.count(
	sv.data = sv,
	positions = size[size$chr %in% chr,1:3],
	type = 'full'
	);

# write counts to file
write.table(
	summary.table,
	file = paste(file.prefix, "chr_distribution.txt", sep = "_"),
	sep = "\t",
	col.names = NA,
	quote = FALSE
	);

chr.dist.object <- create.distribution.plot(
	summary.data = summary.table,
	chr.summary = TRUE
	);
output.plot(
	trellis.objects = list(chr.dist.object),
	filename = paste(file.prefix, "chr_distribution.tiff", sep = "_"),
	width = 12
	);

summary.table.count <- apply(summary.table, 2, sum);
summary.table.count <- as.data.frame(t(summary.table.count));

if(sum(summary.table.count) > 0){
	overall.dist.object <- create.distribution.plot(
		summary.table.count,
		chr.summary = FALSE
		);
	output.plot(
		trellis.objects = list(overall.dist.object),
		filename = paste(file.prefix, "overall_distribution.tiff", sep = "_")
		);
	}


### CLUSTERING BREAKPOINTS #########################################################################
# set adjusted pvalue threshold
adjusted.p.threshold <- 0.1;

print("Clustering based on combined CNV and SV breakpoints");
goodness.fit.values.sv.cnv.results <- process.clustering(
	breakpoint.data = breaks,
	positions = size[size$chr %in% chr,1:3],
	type = 'both',
	threshold = adjusted.p.threshold,
	class = 'full',
	prefix = file.prefix
	);
goodness.fit.values.sv.cnv <- goodness.fit.values.sv.cnv.results$goodness.fit.value;
clustering.object.sv.cnv <- goodness.fit.values.sv.cnv.results$cluster.object;
output.plot(
	trellis.objects = list(clustering.object.sv.cnv),
	filename = paste(file.prefix, "clustering_sv_cnv.tiff", sep="_"),
	width = 12
	);

chr.pass.sv.cnv.clustering <- rownames(goodness.fit.values.sv.cnv[goodness.fit.values.sv.cnv$adjusted.pvalue <= adjusted.p.threshold,]);

# generates distribution plots for those that pass significance
chr.pass.cnv.all <- as.vector(rownames(significant.values.cnv.sv));

### JOIN ANALYSIS ##################################################################################
join.pvalue.threshold = 0.05;
join.analysis <- run.joins(sv, size[size$chr %in% chr,], file.prefix);
pass.pvalue <- join.analysis$pass.pvalue;
pass.pvalue <- pass.pvalue[c('sub.join.pvalue', 'others.join.pvalue')];
colnames(pass.pvalue) <- c('join.pvalue', 'others.pvalue');
trellis.object.list <- join.analysis$trellis.object.list;
chr.pass.join <- rownames(join.analysis$chr.pass.join);
others.dist.obj <- NULL;
chr.pass.others <- rownames(pass.pvalue[!is.na(pass.pvalue$others.pvalue) & pass.pvalue$others.pvalue < join.pvalue.threshold,]);

### OSCILLATING COUNTS #############################################################################
chr.pass.oscillating.values.results <- run.oscillating.test(cnv, all.breaks, size[size$chr %in% chr,], prefix = file.prefix);
chr.pass.oscillating.values <- chr.pass.oscillating.values.results$values;
chr.pass.oscillating <- chr.pass.oscillating.values.results$chr;
oscillating.object.list <- chr.pass.oscillating.values.results$oscillating.object.list;

write.table(
	chr.pass.oscillating.values,
	file = paste(file.prefix, "oscillating_test.txt", sep="_"),
	quote = FALSE,
	sep = "\t",
	col.names = NA
	);

chr.pass.oscillating <- unique(chr.pass.oscillating);

### INTERSPERSED LOH ###############################################################################
# read in AR files
AR <- read.table(
	opt$AR,
	sep = "\t",
	header = TRUE,
	as.is = TRUE
	);

#cn <- cnv;
loh.sign.p.threshold <- 0.1;

# test ar
ar.results <- return.ar.info(
	ar.data = AR,
	cn.data = cnv,
	chromosomes = size[size$chr %in% chr,],
	type = 'full'
	);

### RETURN SUMMARY #################################################################################
overall.results <- return.overall(
	sv.cnv.counts.regions = chr.pass.cnv.all,
	sv.cnv.counts = counts,
	hyper.test.results.sv.cnv = hyper.test.results.SV.CNV,
	clustering.sv.cnv.counts = goodness.fit.values.sv.cnv,
	join.values = pass.pvalue,
	oscillating.values = chr.pass.oscillating.values,
	chr.pass.oscillating = chr.pass.oscillating,
	ar.results = ar.results,
	adjusted.p.threshold = adjusted.p.threshold,
	join.p.value = join.pvalue.threshold,
	chromosomes = chr
	);

### TABLE OF COUNTS ################################################################################
print("Looking for evidence of chained events");
trans.count <- generate.trans.summary.count(
	input.data = sv,
	regions = size[size$chr %in% chr,],
	all = TRUE,
	type = 'full'
	);

inter.chained.summary <- return.chain.inter.reactions(
	sv.counts = trans.count,
	regions = size[size$chr %in% chr,]
	);
chained <- inter.chained.summary$chained.reaction;
inter <- inter.chained.summary$inter.reaction;

reordered.data <- return.trans.count.ordered(
	sv.data = trans.count,
	interchrom = chained);
reordered.trans.count <- reordered.data$sv.data;
clustered.info <- reordered.data$clustered.info;

write.table(
	x = reordered.trans.count,
	paste(file.prefix, "chromosome_counts.txt", sep="_"),
	sep = "\t",
	quote = FALSE,
	col.names = NA
	);

# convert all 0 to NA for plotting
reordered.trans.count[reordered.trans.count ==  0] <- NA;

overall.results$type <- NA;
overall.results$chain <- NA;
if(length(chained) > 0){
	for(i in 1:length(chained)){
		reaction <- chained[i];
		interactor.split <- strsplit(reaction,',')[[1]];
		overall.results[interactor.split,'type'] <- 'chain';
		overall.results[interactor.split,'chain'] <- reaction;
		}
	}

overall.inter <- verify.inter.chromosomal(
	inter.chromosomal.vector = inter,
	list.filters = list(
		chr.pass.cnv.all = chr.pass.cnv.all,
		chr.pass.sv.cnv.clustering = chr.pass.sv.cnv.clustering,
		chr.pass.join = chr.pass.join,
		chr.pass.oscillating.values = chr.pass.oscillating
		),
	chromosomes.pass.others = chr.pass.others
	);

overall.results$chromothripsis <- NA;

if(nrow(overall.inter) > 0){
	for(i in rownames(overall.inter)){
		overall.results[i,'type'] <- 'inter';
		if(!is.na(overall.inter[i,1])){
			overall.results[i,'chromothripsis'] <- overall.inter[i,1];
			}
		}
	}

if(length(chained) > 0){
	chromothripsis.class.chained <- verify.intra.chromosomal(
		chained = chained,
		list.filters = list(
			chr.pass.cnv.all = chr.pass.cnv.all,
			chr.pass.sv.cnv.clustering = chr.pass.sv.cnv.clustering,
			chr.pass.join = chr.pass.join,
			chr.pass.oscillating = chr.pass.oscillating
			),
		chromosomes.pass.others = chr.pass.others
		);

	if(nrow(chromothripsis.class.chained) > 0){
		for(i in rownames(chromothripsis.class.chained)){
			if(!is.na(chromothripsis.class.chained[i,1])){
				overall.results[i,'chromothripsis'] <- chromothripsis.class.chained[i,1];
				}
			}
		}
	}

chained <- unique(overall.results[!is.na(overall.results$chromothripsis) & (overall.results$chromothripsis == 'chromothripsis_chain' | overall.results$chromothripsis == 'chain'),]$chain);
overall.inter <- rownames(overall.results[!is.na(overall.results$chromothripsis) & (overall.results$chromothripsis == 'chromothripsis_inter'),]);
others.fail.inter <- rownames(overall.results[!is.na(overall.results$chromothripsis) & (overall.results$chromothripsis == 'others_fail_inter'),]) ;
others.fail.chain <- unique(overall.results[!is.na(overall.results$chromothripsis) & (overall.results$chromothripsis == 'others_failed_chain' | overall.results$chromothripsis == 'others_failed_chained'),]$chain);

write.table(
	x = overall.results,
	file = paste(file.prefix, "overall_results.txt", sep="_"),
	sep = "\t",
	quote = FALSE,
	col.names = NA
	);

chained.object <- create.chained.reaction.plot(
	sv.data = reordered.trans.count,
	intrachrom = chained,
	interchrom = overall.inter,
	others.fail.inter = others.fail.inter,
	others.fail.chain = others.fail.chain,
	clustered = clustered.info
	);

output.plot(
	trellis.objects = list(chained.object),
	filename = paste(file.prefix, "_chained.tiff", sep=""),
	height = 10,
	width = 10
	);

# create join plots based on chained reaction
if(length(chained) > 0){
	chained.plots.others <- chained.join.plots(chained, sv, chr.pass.join, size[size$chr %in% chr,], type = 'full', file.prefix);
	others.dist.obj <- chained.plots.others$others.dist.obj;
	}

# need to remove anything with X and Y for plotting purposes (no cnv data)
removed.chr <- remove.chr(c('chrX', 'chrY'), overall.inter, chained);
overall.inter <- removed.chr$inter;
chained <- removed.chr$chain;

removed.chr <- remove.chr(c('chrX', 'chrY'), others.fail.inter, others.fail.chain);
others.fail.inter <- removed.chr$inter;
others.fail.chain <- removed.chr$chain;

chained <- chained[!is.na(chained)];
others.fail.chain <- others.fail.chain[!is.na(others.fail.chain)];

### CREATE GAVIN'S PLOT ############################################################################
submitted.plots.results <- submit.plots(
	chained.reaction = chained,
	inter.chr.reaction = overall.inter,
	chain.reaction.failed.others = others.fail.chain, 
	inter.chr.reaction.failed.others = others.fail.inter,
	size = size[size$chr %in% chr,], 
	prefix = file.prefix,
	classification.file = paste(file.prefix, '_sv_classified.txt', sep = ""),
	segment.file = opt$seg,
	copynumber.file = opt$copynumber,
	AR.file = opt$AR
	);

chained.reaction.files <- submitted.plots.results$chained.reaction.files;
inter.files <- submitted.plots.results$inter.files;
chain.fail.files <- submitted.plots.results$chain.fail.files;
inter.fail.files <- submitted.plots.results$inter.fail.files;

### IDENTIFY WINDOWS ###############################################################################
windows <- return.chromothriptic.windows(sv, chromosome = chr);
windows <- combine.windows(windows, cnv,flanking.length = 0, min.count = 2, chrom = chr);
windows <- find.window.range(windows, cnv);
windows <- check.overlap.windows(windows, sv);

write.table(
	windows,
	file = paste(file.prefix, "windows_crossover_rate.txt", sep="_"),
	quote = FALSE,
	col.names = NA,
	sep = "\t"
	);

windows <- windows[is.na(windows$remove),];
windows <- windows[c('chr', 'start', 'end', 'event.count', 'type')];

# recombine for overlap
windows <- merging.windows(windows, chr);

windows$chr.tmp <- factor(windows$chr, levels = chr, ordered = TRUE);
windows <- windows[order(windows$chr.tmp, windows$start),];

windowed.results <- as.data.frame(
	matrix(
		data = 'FAIL',
		nrow = nrow(windows),
		ncol = 5
		),
	stringsAsFactors = FALSE
	);

colnames(windowed.results) <- c('sv.cnv.hyper.adjusted.p.windowed', 'sv.cnv.hyper.windowed', 'cluster.sv.cnv.windowed', 'cluster.sv.cnv.adjusted.p.windowed', 'sv.events.windowed');

if(nrow(windowed.results) > 0){
	rownames(windowed.results) <- paste(windows$chr, ":", windows$start, "-", windows$end, sep = "");
	windowed.results$sv.cnv.hyper.adjusted.p.windowed <- 1;
	windowed.results$sv.events.windowed <- windows$event.count;
	}

### RUN HYPERGEOMETIC TEST #########################################################################
if(nrow(windowed.results) > 0){
	# set adjusted pvalue threshold
	adjusted.p.threshold <- 0.1;

	windows$cnv.count <- 0;
	windows$sv.count <- 0;
	windows$sv.cnv.count <- 0;

	for(i in 1:nrow(windows)){
		windows[i,'sv.count'] <- nrow(all.breaks[all.breaks$chr == windows[i,'chr'] & all.breaks$type == 'sv' & all.breaks$pos >= windows[i,'start'] & all.breaks$pos <= windows[i,'end'],]);
		windows[i,'cnv.count'] <- nrow(all.breaks[all.breaks$chr == windows[i,'chr'] & all.breaks$type == 'cnv' & all.breaks$pos >= windows[i,'start'] & all.breaks$pos <= windows[i,'end'],]);
		windows[i,'sv.cnv.count'] <-nrow(breaks[breaks$chr == windows[i,'chr'] & breaks$pos >= windows[i,'start'] & breaks$pos <= windows[i,'end'],]);
		}

	hyper.test.results.SV.CNV <- run.hyper.test(
		count.data = windows[c('chr', 'start', 'end', 'sv.cnv.count')],
		type = 'windows'
		);

	write.table(
		hyper.test.results.SV.CNV,
		file = paste(file.prefix, "hyper_test_results_SV_CNV_windowed.txt", sep="_"),
		quote = FALSE,
		col.names = NA,
		sep = "\t"
		);

	hyper.test.results.SV.CNV.object.windowed <- create.clustering.plot(
		data = hyper.test.results.SV.CNV['neg.log10.adjusted.pvalue'],
		cutoff = adjusted.p.threshold,
		main = 'SV & CNV Count Distribution',
		xcex = 0.5
		);

	output.plot(
		trellis.objects = list(hyper.test.results.SV.CNV.object.windowed),
		filename = paste(file.prefix, "hyper_test_results_SV_CNV_windowed.tiff", sep="_"),
		width = 12
		);

	windowed.results$sv.cnv.hyper.adjusted.p.windowed <- hyper.test.results.SV.CNV$adjusted.p;
	windowed.results[rownames(hyper.test.results.SV.CNV[hyper.test.results.SV.CNV$adjusted.p <= adjusted.p.threshold,]),'sv.cnv.hyper.windowed'] <- 'PASS';
	windowed.results <- windowed.results[rownames(windowed.results) %in% paste(windows[,'chr'], ":", windows[,'start'], "-", windows[,'end'], sep=""),];
	}

### DISTRIBUTION PLOT ##############################################################################
if(nrow(windowed.results) > 0){
	print("Calculating distribution...");
	summary.table.windowed <- return.dist.count(
		sv.data = sv,
		positions = windows[1:3],
		type = 'windows'
		);

	# write counts to file
	write.table(
		x = summary.table.windowed,
		file = paste(file.prefix, "chr_distribution_windowed.txt", sep = "_"),
		sep = "\t",
		col.names = NA,
		quote = FALSE
		);

	chr.dist.object.windowed <- create.distribution.plot(
		summary.data = summary.table.windowed,
		xcex = 0.5
		);
	output.plot(
		trellis.objects = list(chr.dist.object.windowed),
		filename = paste(file.prefix, "chr_distribution.windowed.tiff", sep = "_"),
		width = 12
		);

	summary.table.count.windowed <- apply(summary.table.windowed, 2, sum);
	summary.table.count.windowed <- as.data.frame(t(summary.table.count.windowed));
	overall.dist.object.windowed <- create.distribution.plot(
		summary.data = summary.table.count.windowed,
		chr.summary = FALSE
		);
	output.plot(
		trellis.objects = list(overall.dist.object.windowed),
		filename = paste(file.prefix, "overall_distribution_windowed.tiff", sep = "_")
		);
	}

### CLUSTERING BREAKPOINTS #########################################################################
if(nrow(windowed.results) > 0){
	# set adjusted pvalue threshold
	adjusted.p.threshold <- 0.1;

	print("Clustering based on combined CNV and SV breakpoints");
	goodness.fit.values.sv.cnv.results.windowed <- process.clustering(
		breakpoint.data = breaks,
		positions = windows[1:3],
		type = 'both',
		threshold = adjusted.p.threshold,
		class = 'windows',
		prefix = file.prefix
		);
	goodness.fit.values.sv.cnv.windowed <- goodness.fit.values.sv.cnv.results.windowed$goodness.fit.value;
	clustering.object.sv.cnv.windowed <- goodness.fit.values.sv.cnv.results.windowed$cluster.object;
	output.plot(
		trellis.objects = list(clustering.object.sv.cnv.windowed),
		filename = paste(file.prefix, "clustering_sv_cnv_windows.tiff", sep="_"),
		width = 12
		);
	windowed.results$cluster.sv.cnv.adjusted.p.windowed <- goodness.fit.values.sv.cnv.windowed$adjusted.pvalue;
	if(sum(goodness.fit.values.sv.cnv.windowed$adjusted.pvalue <= adjusted.p.threshold) > 0){
		windowed.results[rownames(goodness.fit.values.sv.cnv.windowed[goodness.fit.values.sv.cnv.windowed$adjusted.pvalue <= adjusted.p.threshold,]),]$cluster.sv.cnv.windowed<- 'PASS';
		}

	}

### JOIN ANALYSIS ##################################################################################
if(nrow(windowed.results) > 0){

	window.level.joins <- run.joins(sv, windows, file.prefix, type = 'windows');
	pass.pvalue.windowed <- window.level.joins$pass.pvalue;
	colnames(pass.pvalue.windowed) <- c('chr', 'start', 'end', 'sub.join.pvalue.windowed', 'sub.join.distribution.windowed', 'others.join.pvalue.windowed', 'others.join.distribution.windowed');
	trellis.object.list.windowed <- window.level.joins$trellis.object.list;
	chr.pass.join.windowed <- rownames(window.level.joins$chr.pass.join);
	others.dist.obj.windowed <- NULL;
	chr.pass.others.windowed <- rownames(pass.pvalue.windowed[!is.na(pass.pvalue.windowed$others.join.pvalue.windowed) & pass.pvalue.windowed$others.join.pvalue.windowed < join.pvalue.threshold,]);

	### OSCILLATING WINDOWS #########################
	chr.pass.oscillating.values.windows.results <- run.oscillating.test(
		cnv.data = cnv,
		breakpoint.file = all.breaks,
		chromosome.summary = windows,
		type = 'windows',
		prefix = file.prefix
		);
	chr.pass.oscillating.values.windows <- chr.pass.oscillating.values.windows.results$values;
	chr.pass.oscillating.windows <- chr.pass.oscillating.values.windows.results$chr;
	oscillating.object.list.windows <- chr.pass.oscillating.values.windows.results$oscillating.object.list;

	write.table(
		x = chr.pass.oscillating.values.windows,
		file = paste(file.prefix, "oscillating_test_windows.txt", sep="_"),
		quote = FALSE,
		sep = "\t",
		col.names = NA
		);

	chr.pass.oscillating.windows <- unique(chr.pass.oscillating.windows);

	for(i in rownames(chr.pass.oscillating.values.windows)){
		pass.pvalue.windowed[i,'oscillating.all.rounded.median.value.windowed'] <- chr.pass.oscillating.values.windows[i,'all.rounded.median.value'];
		pass.pvalue.windowed[i,'oscillating.pass.rounded.median.value.windowed'] <- chr.pass.oscillating.values.windows[i,'rounded.median.values'];
		pass.pvalue.windowed[i,'oscillating.number.states.greater1.windowed'] <- chr.pass.oscillating.values.windows[i,'num.states.greater.1'];
		pass.pvalue.windowed[i,'oscillating.pass.count.windowed'] <- chr.pass.oscillating.values.windows[i,'total.count.pass.oscillating'];
		pass.pvalue.windowed[i,'pass.oscillating.windows'] <- chr.pass.oscillating.values.windows[i,'pass.fail'];
		}
	}

### AR ANALYSIS ####################################################################################
if(nrow(windowed.results) > 0){
	ar.results.windowed <- return.ar.info(
		ar.data = AR,
		cn.data = cnv,
		chromosomes = windows[c(1:3)],
		type = 'windows'
		);

	for(i in 1:nrow(ar.results.windowed)){
		pass.pvalue.windowed[rownames(ar.results.windowed[i,]),'ar.q.values.windowed'] <- ar.results.windowed[i,'ar.q.values'];
		pass.pvalue.windowed[rownames(ar.results.windowed[i,]),'ar.number.pass.q.windowed'] <- ar.results.windowed[i,'number.pass.q'];
		pass.pvalue.windowed[rownames(ar.results.windowed[i,]),'ar.total.number.windowed'] <- ar.results.windowed[i,'total.comparisons'];
		pass.pvalue.windowed[rownames(ar.results.windowed[i,]),'ar.pass.windowed'] <- ar.results.windowed[i,'ar.change'];
		}
	}

### TABLE OF COUNTS ################################################################################
if(nrow(windowed.results) > 0){
	print("Looking for evidence of chained events");
	trans.count.windows <- generate.trans.summary.count(
		input.data = sv,
		regions = windows[c('chr', 'start', 'end')],
		all = TRUE,
		type = 'windows'
		);

	inter.chained.summary.windows <- return.chain.inter.reactions(
		sv.counts = trans.count.windows,
		regions = windows[c('chr', 'start', 'end')],
		type = 'windows'
		);
	chained.windows <- inter.chained.summary.windows$chained.reaction;
	inter.windows <- inter.chained.summary.windows$inter.reaction;

	for(y in inter.windows){
		tmp.chr <- strsplit(strsplit(y, ":")[[1]], ";")[[1]][1];
		tmp.start <- as.numeric(strsplit(strsplit(y, ":")[[1]], "-")[[2]][1]);
		tmp.end <- as.numeric(strsplit(strsplit(y, ":")[[1]], "-")[[2]][2]);
		windows[windows$chr == tmp.chr & windows$start == tmp.start & windows$end == tmp.end,'type'] <- 'inter';
		}
	windows <- windows[c('chr', 'start', 'end', 'event.count', 'type')];
	windows <- check.overlap.windows.fraction(windows, sv, 0.61);
	windows$per.pass <- 1-windows$percent.fail;

	tmp.value <- windows;
	tmp.value$sample <- opt$output;
	write.table(
		tmp.value,
		file = paste(file.prefix, "windows_crossover.txt", sep="_"),
		quote = FALSE,
		row.names = FALSE,
		sep = "\t"
		);

	windows <- windows[is.na(windows$remove),];
	windows <- windows[c('chr', 'start', 'end', 'event.count', 'type')];

	new.inter.window <- vector();
	for(y in inter.windows){
		tmp.chr <- strsplit(strsplit(y, ":")[[1]], ";")[[1]][1];
		tmp.start <- as.numeric(strsplit(strsplit(y, ":")[[1]], "-")[[2]][1]);
		tmp.end <- as.numeric(strsplit(strsplit(y, ":")[[1]], "-")[[2]][2]);
		if(nrow(windows[windows$chr == tmp.chr & windows$start == tmp.start & windows$end == tmp.end,]) == 1){
			new.inter.window <- c(new.inter.window, paste(tmp.chr, ":", tmp.start, "-", tmp.end, sep = ""));
			}
		}

	removed.inter <- inter.windows[!inter.windows %in% new.inter.window ];
	inter.windows <- new.inter.window;

	if(nrow(trans.count.windows) > 1){
		reordered.data.windows <- return.trans.count.ordered(
			sv.data = trans.count.windows,
			interchrom = chained.windows
			);
		reordered.trans.count.windows <- reordered.data.windows$sv.data;
		clustered.info.windows <- reordered.data.windows$clustered.info;
	}else{
		reordered.trans.count.windows <- trans.count.windows;
		clustered.info.windows <- NULL;
		}

	write.table(
		x = reordered.trans.count.windows,
		file = paste(file.prefix, "region_counts.txt", sep="_"),
		sep = "\t",
		quote = FALSE,
		col.names = NA
		);

	reordered.trans.count.windows[reordered.trans.count.windows ==  0] <- NA;

	pass.pvalue.windowed$type.windows <- NA;
	pass.pvalue.windowed$chain.windows <- NA;
	if(length(chained.windows) > 0){
		for(i in 1:length(chained.windows)){
			reaction <- chained.windows[i];
			interactor.split <- strsplit(reaction,',')[[1]];
			pass.pvalue.windowed[interactor.split,'type.windows'] <- 'chain';
			pass.pvalue.windowed[interactor.split,'chain.windows'] <- reaction;
			}
		}

	overall.inter.windows <- verify.inter.chromosomal(
		inter.chromosomal.vector = inter.windows,
		list.filters = list(
			chr.pass.cnv.all = chr.pass.cnv.all,
			chr.pass.sv.cnv.clustering = chr.pass.sv.cnv.clustering,
			chr.pass.join = chr.pass.join.windowed,
			chr.pass.oscillating.values = chr.pass.oscillating.windows
			),
		chromosomes.pass.others = chr.pass.others.windowed
		);

	pass.pvalue.windowed$chromothripsis.windowed <- NA;
	if(nrow(overall.inter.windows) > 0){
		pass.pvalue.windowed[rownames(overall.inter.windows),'type.windows'] <- 'inter';
		for(i in rownames(overall.inter.windows)){
			if(!is.na(overall.inter.windows[i,1])){
				pass.pvalue.windowed[i,'chromothripsis.windowed'] <- overall.inter.windows[i,1];
				}
			}
		}

	if(length(chained.windows) > 0){
		chromothripsis.class.chained.windowed <- verify.intra.chromosomal(
			chained = chained.windows,
			list.filters = list(
				chr.pass.cnv.all = chr.pass.cnv.all,
				chr.pass.sv.cnv.clustering = chr.pass.sv.cnv.clustering,
				chr.pass.join = chr.pass.join.windowed,
				chr.pass.oscillating = chr.pass.oscillating.windows
				),
			chromosomes.pass.others = chr.pass.others.windowed
			);

		if(nrow(chromothripsis.class.chained.windowed) > 0){
			for(i in rownames(chromothripsis.class.chained.windowed)){
				pass.pvalue.windowed[i,'chromothripsis.windowed'] <- chromothripsis.class.chained.windowed[i,1];
				}
			}
		}

	chained.windowed <- unique(pass.pvalue.windowed[!is.na(pass.pvalue.windowed$chromothripsis.windowed) & (pass.pvalue.windowed$chromothripsis.windowed == 'chromothripsis_chain' | pass.pvalue.windowed$chromothripsis.windowed == 'chain'),]$chain.windows);
	overall.inter.windowed <- rownames(pass.pvalue.windowed[!is.na(pass.pvalue.windowed$chromothripsis.windowed) & (pass.pvalue.windowed$chromothripsis.windowed == 'chromothripsis_inter'),]);
	others.fail.inter.windowed <- rownames(pass.pvalue.windowed[!is.na(pass.pvalue.windowed$chromothripsis.windowed) & (pass.pvalue.windowed$chromothripsis.windowed == 'others_fail_inter'),]) ;
	others.fail.chain.windowed <- unique(pass.pvalue.windowed[!is.na(pass.pvalue.windowed$chromothripsis.windowed) & (pass.pvalue.windowed$chromothripsis.windowed == 'others_failed_chain' | pass.pvalue.windowed$chromothripsis.windowed == 'others_failed_chained'),]$chain.windows);

	windowed.results <- cbind(windowed.results, pass.pvalue.windowed);

	overall.summary.results <- combine.window.overall.results(
		windowed.data = windowed.results,
		nonwindowed.data = overall.results,
		sample = opt$output,
		chromosomes = chr
		);

	if(!is.null(clustered.info.windows)){
		chained.object.windowed <- create.chained.reaction.plot(
			sv.data = reordered.trans.count.windows,
			intrachrom = chained.windowed,
			interchrom = overall.inter.windowed,
			others.fail.inter = others.fail.inter.windowed,
			others.fail.chain = others.fail.chain.windowed,
			clustered = clustered.info.windows
			);

		output.plot(
			trellis.objects = list(chained.object.windowed),
			filename = paste(file.prefix, "_chained_windowed.tiff", sep=""),
			height = 10,
			width = 10
			);
		}else{
			chained.object.windowed <- NULL;
			}

	# create join plots based on chained reaction
	if(length(chained.windowed) > 0){
		chained.join.plots.windowed <- chained.join.plots(chains = chained.windows, sv.data = sv, chr.pass.join.values = chr.pass.join.windowed, type = 'windows', filename.prefix = file.prefix);
		others.dist.obj.windowed <- chained.join.plots.windowed$others.dist.obj;
		}

	if(nrow(windows) > 0){
		tmp.window <- data.frame(cbind(paste(windows[,1], ":", windows[,2], "-", windows[,3], sep = ""), windows$start, windows$end), stringsAsFactors = FALSE);
		tmp.window$X2 <- as.numeric(tmp.window$X2);
		tmp.window$X3 <- as.numeric(tmp.window$X3);
	}else{
		tmp.window <- as.data.frame(matrix(nrow = 0, ncol = 3));
		colnames(tmp.window) <- c('X1', 'X2', 'X3');
		}

	# need to remove anything with X and Y
	removed.chr <- remove.chr(c('chrX', 'chrY'), overall.inter.windowed, chained.windowed, type = 'windows');
	overall.inter.windowed <- removed.chr$inter;
	chained.windowed <- removed.chr$chain;

	removed.chr <- remove.chr(c('chrX', 'chrY'), others.fail.inter.windowed, others.fail.chain.windowed, type = 'windows');
	others.fail.inter.windowed <- removed.chr$inter;
	others.fail.chain.windowed <- removed.chr$chain;

	chained.windowed <- chained.windowed[!is.na(chained.windowed)];
	others.fail.chain.windowed <- others.fail.chain.windowed[!is.na(others.fail.chain.windowed)];

	chained.windowed  <- return.whole.chr(
		data = chained.windowed,
		type = 'chain'
		);

	overall.inter.windowed  <- return.whole.chr(
		data = overall.inter.windowed,
		type = 'inter'
		);

	others.fail.chain.windowed  <- return.whole.chr(
		data = others.fail.chain.windowed,
		type = 'chain'
		);

	others.fail.inter.windowed  <- return.whole.chr(
		data = others.fail.inter.windowed,
		type = 'inter'
		);

	monitor.jobs(file.prefix);
	submitted.plots.results.windows <- submit.plots(
		chained.reaction = chained.windowed,
		inter.chr.reaction = overall.inter.windowed,
		chain.reaction.failed.others = others.fail.chain.windowed,
		inter.chr.reaction.failed.others = others.fail.inter.windowed,
		size = rbind(tmp.window,setNames(size[size$chr %in% chr,], c('X1', 'X2', 'X3'))), 
		prefix = paste(file.prefix, "_windowed", sep=""),
		classification.file = paste(file.prefix, '_sv_classified.txt', sep = ""),
		segment.file = opt$seg,
		copynumber.file = opt$copynumber,
		AR.file = opt$AR
		);

	chained.reaction.files.windowed <- submitted.plots.results.windows$chained.reaction.files;
	inter.files.windowed <- submitted.plots.results.windows$inter.files;
	chain.fail.files.windowed <- submitted.plots.results.windows$chain.fail.files;
	inter.fail.files.windowed <- submitted.plots.results.windows$inter.fail.files;

}else{
	chained.windowed <- vector();
	overall.inter.windowed <- vector();
	others.fail.chain.windowed <- vector();
	others.fail.inter.windowed <- vector();

	windowed.results <- as.data.frame(
		matrix(
			nrow = 24,
			ncol = 24 
			)
		);

	colnames(windowed.results) <- c("sv.cnv.hyper.adjusted.p.windowed", "sv.cnv.hyper.windowed", "cluster.sv.cnv.windowed", "cluster.sv.cnv.adjusted.p.windowed", "sv.events.windowed", "chr", "start", "end", "sub.join.pvalue.windowed", "sub.join.distribution.windowed", "others.join.pvalue.windowed", "others.join.distribution.windowed", "oscillating.all.rounded.median.value.windowed", "oscillating.pass.rounded.median.value.windowed", "oscillating.number.states.greater1.windowed", "oscillating.pass.count.windowed", "pass.oscillating.windows", "ar.q.values.windowed", "ar.number.pass.q.windowed", "ar.total.number.windowed", "ar.pass.windowed", "type.windows", "chain.windows", "chromothripsis.windowed");

	windowed.results$chr <- chr;
	overall.summary.results <- combine.window.overall.results(windowed.data = windowed.results, nonwindowed.data = overall.results, opt$output);
	overall.summary.results$chr <- factor(overall.summary.results$chr, levels = chr, ordered = TRUE);
	overall.summary.results <- overall.summary.results[order(overall.summary.results$chr),];
	rownames(overall.summary.results) <-overall.summary.results$chr;

	removed.inter <- vector();
	}

write.table(
	overall.summary.results,
	file = paste(file.prefix, "_overall_summary_results.txt", sep=""),
	sep = "\t",
	quote = FALSE,
	row.names = FALSE,
	col.names = TRUE
	);

# quick summary
if(!is.null(opt$output)){
	sample <- rep(opt$output, nrow(overall.summary.results));
}else{
	sample <- rep("", nrow(overall.summary.results));
	}

tmp.summary <- cbind(sample, overall.summary.results[c('chr', 'start', 'end', 'sv.cnv.dist', 'sv.cnv.hyper', 'cluster.sv.cnv', 'sub.join.distribution', 'others.join.distribution', 'sv.cnv.hyper.windowed', 'sub.join.distribution.windowed', 'others.join.distribution.windowed', 'chr.pass.oscillating', 'pass.oscillating.windows', 'ar.pass', 'ar.pass.windowed', 'type', 'type.windows', 'chain', 'chain.windows', 'chromothripsis', 'chromothripsis.windowed')]);
write.table(
	tmp.summary,
	file = paste(file.prefix, "_overall_summary_results_brief.txt", sep=""),
	sep = "\t",
	quote = FALSE,
	row.names = FALSE,
	col.names = TRUE
	);

# write output
output.plot(
	trellis.objects = list(chained.object),
	filename = paste(file.prefix, '_chained.tiff', sep=""),
	width = 12,
	height = 10
	);

#summary total counts
overall.summary.results <- summarize.classification(overall.summary.results, removed.inter);

write.table(
	overall.summary.results,
	file = paste(file.prefix, "_overall_summary_results_complete.txt", sep=""),
	sep = "\t",
	quote = FALSE,
	row.names = FALSE,
	col.names = TRUE
	);

run.combine.plot <- overall.summary.results[!is.na(overall.summary.results$dist.bin.clust.join.oscil.count.window.join),];

# this if inter chr
run.combined.plot.inter <- run.combine.plot[run.combine.plot$type.windows == 'inter' | is.na(run.combine.plot$type.windows),];
dist.bin.clust.join.oscil.count.window.join.inter <- vector();
if(nrow(run.combined.plot.inter) > 0){
	for(i in 1:nrow(run.combined.plot.inter)){
		if(is.na(run.combined.plot.inter[i,'start'])){
			dist.bin.clust.join.oscil.count.window.join.inter <- c(dist.bin.clust.join.oscil.count.window.join.inter, as.character(run.combined.plot.inter[i,'chr']));
		}else{
			dist.bin.clust.join.oscil.count.window.join.inter <- c(dist.bin.clust.join.oscil.count.window.join.inter, paste(as.character(run.combined.plot.inter[i,'chr']), paste(run.combined.plot.inter[i,'start'], run.combined.plot.inter[i,'end'], sep = "-"), sep = ":"));
			}
		}
	}

dist.bin.clust.join.oscil.count.window.join.inter <- dist.bin.clust.join.oscil.count.window.join.inter[!dist.bin.clust.join.oscil.count.window.join.inter %in% others.fail.inter.windowed];
dist.bin.clust.join.oscil.count.window.join.inter <- dist.bin.clust.join.oscil.count.window.join.inter[!dist.bin.clust.join.oscil.count.window.join.inter %in% overall.inter.windowed];
dist.bin.clust.join.oscil.count.window.join.inter <- dist.bin.clust.join.oscil.count.window.join.inter[!is.na(dist.bin.clust.join.oscil.count.window.join.inter)];

run.combined.plot.intra <- run.combine.plot[run.combine.plot$type == 'chain',];
dist.bin.clust.join.oscil.count.window.join.chain <- vector();
if(nrow(run.combined.plot.intra) > 0){
	for(i in 1:nrow(run.combined.plot.intra)){
		if(!is.null(run.combined.plot.intra[i,'chain.windows'])){
			dist.bin.clust.join.oscil.count.window.join.chain <- c(dist.bin.clust.join.oscil.count.window.join.chain, run.combined.plot.intra[i,'chain.windows']);
		}else{
			dist.bin.clust.join.oscil.count.window.join.chain <- c(dist.bin.clust.join.oscil.count.window.join.chain, run.combined.plot.intra[i,'chain']);
			}
		}
	}

dist.bin.clust.join.oscil.count.window.join.chain <- unique(dist.bin.clust.join.oscil.count.window.join.chain);
dist.bin.clust.join.oscil.count.window.join.chain <- dist.bin.clust.join.oscil.count.window.join.chain[!dist.bin.clust.join.oscil.count.window.join.chain %in% others.fail.chain.windowed];
dist.bin.clust.join.oscil.count.window.join.chain <- dist.bin.clust.join.oscil.count.window.join.chain[!dist.bin.clust.join.oscil.count.window.join.chain %in% chained.windowed];
dist.bin.clust.join.oscil.count.window.join.chain <- dist.bin.clust.join.oscil.count.window.join.chain[!is.na(dist.bin.clust.join.oscil.count.window.join.chain)];

# need to remove anything with X and Y
removed.chr <- remove.chr(c('chrX', 'chrY'), dist.bin.clust.join.oscil.count.window.join.inter, dist.bin.clust.join.oscil.count.window.join.chain, type = 'windows');
dist.bin.clust.join.oscil.count.window.join.inter <- removed.chr$inter;
dist.bin.clust.join.oscil.count.window.join.chain <- removed.chr$chain;

dist.bin.clust.join.oscil.count.window.join.inter <- dist.bin.clust.join.oscil.count.window.join.inter[!dist.bin.clust.join.oscil.count.window.join.inter %in% others.fail.inter.windowed];
dist.bin.clust.join.oscil.count.window.join.inter <- dist.bin.clust.join.oscil.count.window.join.inter[!dist.bin.clust.join.oscil.count.window.join.inter %in% overall.inter.windowed];

dist.bin.clust.join.oscil.count.window.join.chain <- dist.bin.clust.join.oscil.count.window.join.chain[!is.na(dist.bin.clust.join.oscil.count.window.join.chain)];
dist.bin.clust.join.oscil.count.window.join.chain <- dist.bin.clust.join.oscil.count.window.join.chain[!dist.bin.clust.join.oscil.count.window.join.chain %in% others.fail.chain.windowed];
dist.bin.clust.join.oscil.count.window.join.chain <- dist.bin.clust.join.oscil.count.window.join.chain[!dist.bin.clust.join.oscil.count.window.join.chain %in% chained.windowed];

dist.bin.clust.join.oscil.count.window.join.chain <- return.whole.chr(
	data = dist.bin.clust.join.oscil.count.window.join.chain,
	type = 'chain'
	);

dist.bin.clust.join.oscil.count.window.join.inter <- return.whole.chr(
	data = dist.bin.clust.join.oscil.count.window.join.inter,
	type = 'inter'
	);

monitor.jobs(file.prefix);

submitted.plots.results.windows.rescued <- submit.plots(
	chained.reaction = NULL,
	inter.chr.reaction = NULL,
	chain.reaction.failed.others = dist.bin.clust.join.oscil.count.window.join.chain,
	inter.chr.reaction.failed.others = dist.bin.clust.join.oscil.count.window.join.inter,
	size = rbind(tmp.window,setNames(size[size$chr %in% chr,], c('X1', 'X2', 'X3'))),
	prefix = paste(file.prefix, "_windowed_rescued", sep = ""),
	classification.file = paste(file.prefix, '_sv_classified.txt', sep = ""),
	segment.file = opt$seg,
	copynumber.file = opt$copynumber,
	AR.file = opt$AR
	);
dist.bin.clust.join.oscil.count.window.join.files.inter <- submitted.plots.results.windows.rescued$inter.fail.files;
dist.bin.clust.join.oscil.count.window.join.files.intra <- submitted.plots.results.windows.rescued$chain.fail.files;

monitor.jobs(file.prefix);

### GENERATE OUTPUT PDF ############################################################################
out.results <- return.final.output.tables(chr.wide.summary = overall.results, windows.value = windows, window.wide.summary = overall.summary.results);
merged.tables <- out.results$chr.level.tables;
tmp.merged.tables <- out.results$window.level.tables;
final.output.list <- out.results$final.output.list;
chr.wide.summary <- out.results$chr.wide.summary;
tmp.output <- out.results$window.wide.summary;

trellis.object.list$scatterplot.sv.cnv.object <- scatterplot.sv.cnv.object;
trellis.object.list$chr.dist.object <- chr.dist.object;
trellis.object.list$clustering.object.sv.cnv <- clustering.object.sv.cnv;
trellis.object.list$hyper.test.results.SV.CNV.object <- hyper.test.results.SV.CNV.object;

if(!is.null(others.dist.obj)){
	trellis.object.list <- c(trellis.object.list, others.dist.obj);
	}

files <- c(inter.files, chained.reaction.files, inter.fail.files, chain.fail.files);

files1 <- NULL;
run.window = FALSE;
if(nrow(windows) > 0){
	trellis.object.list.windowed$hyper.test.results.SV.CNV.object.windowed <- hyper.test.results.SV.CNV.object.windowed;

	if(!is.null(others.dist.obj.windowed)){
		trellis.object.list.windowed <- c(trellis.object.list.windowed, others.dist.obj.windowed);
		}

	trellis.object.list.windowed.1 <- trellis.object.list.windowed;
	trellis.object.list.windowed <- c(trellis.object.list.windowed.1, oscillating.object.list.windows);

	files1 <- c(inter.files.windowed, chained.reaction.files.windowed, inter.fail.files.windowed, chain.fail.files.windowed, dist.bin.clust.join.oscil.count.window.join.files.inter, dist.bin.clust.join.oscil.count.window.join.files.intra);
	run.window = TRUE;
	}

# clear files to use if nothing in end if called
if(sum(!is.na(overall.summary.results$dist.bin.clust.join.oscil.count.window.join)) == 0){
	files1 <- NULL;
	}

# calculate number lines
new.lines <- final.output.list[[1]][,9];
new.lines <- new.lines[!is.na(new.lines)];
if(length(new.lines) > 0){
	lines.appearing <- 2 + (length(unlist(strsplit(new.lines, "\n"))) - length(new.lines)) + nrow(final.output.list[[1]]);
}else{
	lines.appearing <- 24;
	}


if(lines.appearing > 50){
	fontsize = 3.35;
	v.pad = 1.5;
}else{
	fontsize = 4;
	v.pad = 2;
	}

create.output.pdf(
	prefix = file.prefix,
	trellis.object.list = trellis.object.list,
	oscillating.object.list = oscillating.object.list,
	merged.tables = merged.tables,
	chr.wide.summary = chr.wide.summary,
	chained.object = chained.object,
	files = files,
	trellis.object.list.windowed = trellis.object.list.windowed,
	merged.tables.windowed = tmp.merged.tables,
	windowed.results = tmp.output,
	chained.object.windowed = chained.object.windowed,
	files.windowed = files1,
	final.output = final.output.list,
	windows = run.window
	);

create.output.pdf.simplified(
	prefix = file.prefix,
	files.windowed = files1,
	final.output = final.output.list,
	samplename = opt$output,
	windows = run.window
	);

create.output.pdf.simplified.reformatted(
	prefix = file.prefix,
	files.windowed = files1,
	final.output = final.output.list,
	samplename = opt$output,
	windows = run.window,
	fontsize = fontsize
	);


print("Processing complete");

### WRITE SESSION INFO #############################################################################
filename <- paste(file.prefix, "_session_info.txt", sep="");
sink(file = filename);
print(sessionInfo());
sink();

### run_celluloid_v0.11.R ##########################################################################
# code to run celluloid v0.11
# Usage:
# Rscript pipeline.r Data/tumour.wig Data/normal.wig Data/gc.wig Data/map.wig Data/AR.txt samplename

### PREAMBLE #######################################################################################
# load libr ry
library(getopt);

### FUNCTIONS ######################################################################################
# usage
# function to return usage
usage <- function(){
usage.text <-
	'\nrun_celluloid_0.11.R
		-h help
		-T tumor.wig (optional. if set requires normal.wig, gc.wig, map.wig, ar.file)
		-N normal.wig (optional. if set requires tumor.wig, gc.wig, map.wig, ar.file)
		-G gc.wig (optional. if set requires tumor.wig, normal.wig, map.wig, ar.file, defaults to /u/mchan/PCSI/michelle/references/autosome_gc.wig)
		-M map.wig (optional. if set requires tumor.wig, normal.wig, gc.wig, ar.file, defaults to /u/mchan/PCSI/michelle/references/autosome_map.wig)
		-A ar.file (optional. if set requires tumor.wig, normal.wig, gc.wig, map.wig)
		-O output.prefix (defaults to "")
		-n nc.rda (optional. if set requires n.seg.rda, tc.rda, t.seg.rda, t.seg.mask.rda, ar.rda, t.ar.seg.rda, copyAr.ra, cntr.rda)
		-r n.seg.rda (optional. if set requires nc.rda, tc.rda, t.seg.rda, t.seg.mask.rda, ar.rda, t.ar.seg.rda, copyAr.rda, cntr.rda)
		-t tc.rda (optional. if set requires nc.rda, n.seg.rda, t.seg.rda, t.seg.mask.rda, ar.rda, t.ar.seg.rda, copyAr.rda, cntr.rda)
		-u t.seg.rda (optional. if set requires nc.rda, n.seg.rda, tc.rda, t.seg.mask.rda, ar.rda, t.ar.seg.rda, copyAr.rda, cntr.rda)
		-m t.seg.mask.rda (optional. if set requires nc.rda, n.seg.rda, tc.rda, t.seg.rda, ar.rda, t.ar.seg.rda, copyAr.rda, cntr.rda)
		-a ar.rda (optional. if set requires nc.rda, n.seg.rda, tc.rda, t.seg.rda, t.seg.mask.rda, t.ar.seg.rda, copyAr.rda, cntr.rda)
		-s t.ar.seg.rda (optional. if set requires nc.rda, n.seg.rda, tc.rda, t.seg.rda, t.seg.mask.rda, ar.rda, copyAr.rda, cntr.rda)
		-c copyAr.rda (optional. if set requires nc.rda, n.seg.rda, tc.rda, t.seg.rda, t.seg.mask.rda, ar.rda, t.ar.seg.rda, cntr.rda)
		-C cntr.rda (optional. if set requires nc.rda, n.seg.rda, tc.rda, t.seg.rda, t.seg.mask.rda, ar.rda, t.ar.seg.rda, copyAr.rda)
	EXAMPLES:
		Rscript run_celluloid_0.11.R -T tumor.wig -N normal.wig -G gc.wig -M map.wig -A ar.file -O test\n
		Rscript run_celluloid_0.11.R -n nc.rda -r n.seg.rda -t tc.rda -u t.seg.rda -m t.seg.mask.rda -a ar.rda -s t.ar.seg.rda -c copyAR.rda -C cntr.rda -O test\n\n'
	return(usage.text);
	}

### OBTAIN COMMAND LINE ARGUMENTS ##################################################################
params <- matrix(
	c(
		'help', 'h', 0, "logical",
		'tumor', 'T', 2, "character",
		'normal', 'N', 2, "character",
		'gc', 'G', 2, "character",
		'map', 'M', 2, "character",
		'ar', 'A', 2, "character",
		'output', 'O', 2, "character",
		'nc.rda', 'n', 2, "character",
		'n.seg.rda', 'r', 2, "character",
		'tc.rda', 't', 2, "character",
		't.seg.rda', 'u', 2, "character",
		't.seg.mask.rda', 'm', 2, "character",
		'ar.rda', 'a', 2, "character",
		't.ar.seg.rda', 's', 2, "character",
		'copyAr.rda', 'c', 2, "character",
		'cntr.rda', 'C', 2, "character"
		),
	ncol = 4,
	byrow = TRUE
	);

opt = getopt(params);

if(is.null(opt$gc)){
	opt$gc <- '/u/mchan/PCSI/michelle/references/autosome_gc.wig';
	}

if(is.null(opt$map)){
	opt$map <- '/u/mchan/PCSI/michelle/references/autosome_map.wig';
	}

# verified required argments and sets defaults
if (!is.null(opt$help) || ((is.null(opt$ar) ||  is.null(opt$tumor) || is.null(opt$normal) || is.null(opt$gc) || is.null(opt$map)) && (is.null(opt$nc.rda) || is.null(opt$n.seg.rda) || is.null(opt$tc.rda) || is.null(opt$t.seg.rda) || is.null(opt$t.seg.mask.rda) || is.null(opt$ar.rda) || is.null(opt$t.ar.seg.rda) || is.null(opt$t.ar.seg.rda) || is.null(opt$copyAr.rda) || is.null(opt$cntr.rda)) ) ){
	cat(usage());
	q(status=1)
	}

if(is.null(opt$output)){
	opt$output <- "";
	}

### PREAMBLE #######################################################################################
library(celluloid);

# general parameters
date.value <- Sys.Date();

### PREPARE DATA ###################################################################################
if(is.null(opt$nc.rda)){
	# load and segment normal data
	n <- wigsToRangedData(
		opt$normal,
		opt$gc,
		opt$map
		);

	# gc correction
	nc <- gcCorrect(n);
	n.seg <- segmentSeqData(
		nc,
		k = 50,
		maskmap = 0.8
		);

	output.file = paste(date.value, "_", opt$output, "_n.seg.rda", sep = "");
	save(
		n.seg,
		file = output.file
		);

	output.file = paste(date.value, "_", opt$output, "_nc.rda", sep = "");
	save(
		nc,
		file = output.file
		);

}else{

	load(opt$n.seg.rda);
	load(opt$nc.rda);
	}


if(is.null(opt$tc.rda) ){
	# load and segment tumor data
	t <- wigsToRangedData(
	opt$tumor,
	opt$gc,
	opt$map
	);

	tc <- gcCorrect(
		t,
		sampletype = "tumor"
		);
	t.seg <- segmentSeqData(
		tc,
		k = 50
		);

	output.file = paste(date.value, "_", opt$output, "_t.seg.rda", sep = "");
	save(
		t.seg,
		file = output.file
		);

	output.file = paste(date.value, "_", opt$output, "_tc.rda", sep = "");
	save(
		tc,
		file = output.file
		);

}else{
	load(opt$t.seg.rda);
	load(opt$tc.rda);
	}


# identify segments not "normal" in normal
if(is.null(opt$t.seg.mask.rda) ){
	t.n.seg <- intersectSegments(
		t.seg,
		n.seg
		);
	t.n.seg <- t.n.seg[ !is.na( t.n.seg$mean) & !is.na( t.n.seg$mean.1),];
	sel <- n.seg$end.pos-n.seg$start.pos > 150000 & n.seg$meanmap>.8;

	bp <- boxplot(
		n.seg$mean[sel],
		range = 3,
		plot = F
		);

	nrange <- c(bp$stats[1,1], bp$stats[5,1] );
	mask <- t.n.seg$mean.1>nrange[2] | t.n.seg$mean.1<nrange[1];
	t.seg.mask <- t.n.seg[,1:8];
	t.seg.mask$mask <- mask;

	output.file = paste(date.value, "_", opt$output, "_t.seg.mask.rda", sep = "");
	save(
		t.seg.mask,
		file = output.file
		);

}else{
	load(opt$t.seg.mask.rda);
	}

if( is.null(opt$t.ar.seg.rda) ){
	# reading allelic ratio file and segmenting
	ar <- read.table(
		opt$ar,
		head =T
		);

	output.file = paste(date.value, "_", opt$output, "_ar.rda", sep = "");
	save(
		ar,
		file = output.file
		);

	ar.seg <- segmentAR(
		ar,
		tc
		);

	output.file = paste(date.value, "_", opt$output, "_ar.seg.rda", sep = "");
	save(
		ar.seg,
		file = output.file
		);

	t.ar.seg <- intersectSegments(
		t.seg.mask,
		ar.seg
		);
	t.ar.seg <- t.ar.seg[ !apply( is.na( t.ar.seg[,c("mean","meanar")] ), 1, any),];
	t.ar.seg <-arInSeg(
		t.ar.seg,
		ar,
		tumourrangedata = tc,
		minhet = 50
		);

	output.file = paste(date.value, "_", opt$output, "_t.ar.seg.rda", sep = "");
	save(
		t.ar.seg,
		file = output.file
		);

}else{
	load(opt$ar.rda);
	load(opt$t.ar.seg.rda);
	}

if( is.null(opt$copyAr.rda) ){
	# creates the object used to draw a contour plot.
	mask <- t.ar.seg$mask | is.na(t.ar.seg$mask);
	copyAr <- prepCopyAr(
		t.ar.seg[ !mask ,],
		ar,
		tc
		);

	output.file = paste(date.value, "_", opt$output, "_copyAr.rda", sep = "");
	save(
		copyAr,
		file = output.file
		);

}else{
	load(opt$copyAr.rda);
	}

if( is.null( opt$cntr.rda) ){
	cntr <- showTumourProfile(
		copyAr,
		flatten = .25,
		nlev = 20,
		noise = 0.01,
		maxPoints = 200000,
		plot = F
		);

	output.file = paste(date.value, "_", opt$output, "_cntr.rda", sep = "");
	save(
		cntr,
		file = output.file
		);

}else{
	load(opt$cntr.rda);
	}

### ANALYSIS #######################################################################################
# finding LOH curve
Sn <- estimateLOHcurve(t.ar.seg);

# grid search
sel <- t.ar.seg$mean < max( cntr$x ) & !is.na( t.ar.seg$mean ) & !is.na( t.ar.seg$p );
segmentsubset <-t.ar.seg[sel,];

# using optim

# upper bound for %normal cell set to Sn/.25, set to 1 if you think your sample's ploidy can 
# be above 9 (2/.25)

lm <- coverParamSpace(
	segments = segmentsubset,
	control = list(
		maxit = 1000
		),
	Sn = Sn,
	maxc = 12,
	optimFct = 1,
	nrep = 50,
	lowerF = c(0),
	upperF = c(Sn/.25),
	addToParamSpace = T
	);

output.file = paste(date.value, "_", opt$output, "_lm.rda", sep = "");
save(
	lm,
	file = output.file
	);

localSolutions <- getLocalSolutions(
	lm,
	max = TRUE
	);

# only solutions for which 10% or more of the genome was captured are in the output
localSolutions <- localSolutions[ localSolutions$value>0.10 , ];

pdf(
	paste(date.value, "_", opt$output, "contour_solutions_%1d.pdf", sep = ""),
	height = 12,
	width = 16
	);

par(
	mfrow = c(3,2)
	);

plot(
	paramSpace[,3],
	100*paramSpace[,1],
	pch = 19,
	ylim = c(0,100),
	xlab = "%normal",
	ylab = "percent captured by model"
	);

for( i in 1:nrow(localSolutions) ){

	r <- as.numeric(rownames(localSolutions)[i] );

	image(
		cntr,
		col = terrain.colors(50)
		);

	contour(
		cntr,
		nlev = 8,
		add = T
		);

	sel <- t.ar.seg$size>1000000 & !t.ar.seg$mask & t.ar.seg$meanmap>.9;
	le <- t.ar.seg$end.pos[sel]-t.ar.seg$start.pos[sel];

	cxcut <- as.integer( cut( le, c(100000,1000000,5000000,10000000,20000000,50000000,Inf) ) )/3;

	points(
		x <- t.ar.seg$mean[sel],
		y <- t.ar.seg$p[sel],
		pch = 21,
		col = "blue",
		lwd = 3,
		cex = cxcut
		);

	points(
		t.ar.seg$mean[sel],
		t.ar.seg$p[sel],
		pch = 19,
		col = "white",
		cex = cxcut  - .5
		);

	points(
		t.ar.seg$mean[sel],
		1-t.ar.seg$p[sel],
		pch = 21,
		col = "blue",
		lwd = 3,
		cex = cxcut
		);

	points(
		t.ar.seg$mean[sel],
		1-t.ar.seg$p[sel],
		pch = 19,
		col = "white",
		cex = cxcut -.5
		);

	xxx <- seq( Sn, 2, .01 );

	points(
		xxx,
		arloh <- ARloh( xxx , 1 , Sn ),
		type = 'l',
		lwd = 2
		);

	points(
		xxx,
		1-arloh,
		type = 'l',
		lwd = 2
		);

	epp <- plotModelPeaks(
		lm[[r]]$par,
		selectedPoints = NULL,
		cn = cn,
		epcol = "red",
		epcex = 1,
		eplwd = 3,
		addlabels = F
		);

	legend(
		0,
		1,
		floor(1000*( lm[[r]]$value) )/10,
		bg = "white"
		);
	}
dev.off()

pdf(
	paste(date.value, "_", opt$output, "final_solution.pdf", sep = ""),
	height = 12,
	width = 16
	);

	i <- 1;

	r <- as.numeric(rownames(localSolutions)[i] );

	image(
		cntr,
		col = terrain.colors(50)
		);

	contour(
		cntr,
		nlev = 8,
		add = T
		);

	sel <- t.ar.seg$size>1000000 & !t.ar.seg$mask & t.ar.seg$meanmap>.9;
	le <- t.ar.seg$end.pos[sel]-t.ar.seg$start.pos[sel];

	cxcut <- as.integer( cut( le, c(100000,1000000,5000000,10000000,20000000,50000000,Inf) ) )/3;

	xxx <- seq( Sn, 2, .01 );

	epp <- plotModelPeaks(
		lm[[r]]$par,
		selectedPoints = NULL,
		cn = cn,
		epcol = "red",
		epcex = 1,
		eplwd = 3,
		addlabels = F
		);

dev.off()

pdf(
	paste(date.value, "_", opt$output, "2nd_solution.pdf", sep = ""),
	height = 12,
	width = 16
	);

	i <- 2;

	r <- as.numeric(rownames(localSolutions)[i] );

	image(
		cntr,
		col = terrain.colors(50)
		);

	contour(
		cntr,
		nlev = 8,
		add = T
		);

	sel <- t.ar.seg$size>1000000 & !t.ar.seg$mask & t.ar.seg$meanmap>.9;
	le <- t.ar.seg$end.pos[sel]-t.ar.seg$start.pos[sel];

	cxcut <- as.integer( cut( le, c(100000,1000000,5000000,10000000,20000000,50000000,Inf) ) )/3;

	xxx <- seq( Sn, 2, .01 );

	epp <- plotModelPeaks(
		lm[[r]]$par,
		selectedPoints = NULL,
		cn = cn,
		epcol = "red",
		epcex = 1,
		eplwd = 3,
		addlabels = F
		);

dev.off()

## OUTPUTS RESULTS FOR 1st solution
prepCN(12,1,NULL);
ePP <- ePeakPos(
	S = localSolutions[1,'S'],
	t = c(localSolutions[1,'N'], localSolutions[1,'T1']),
	cn = cn
	);

#ePP <- ePeakPos( par=lm[[1]]$par  , cn=cn  )

tcs <- scaleReadCounts( tc , ePP );

tmp.seg <- subset(
	as.data.frame(tcs),
	select = c(space,start,end,icopy)
	 );
tmp.seg$sample <- opt$output;
tmp.seg <- tmp.seg[c('sample', 'space', 'start', 'end', 'icopy')];

# convert to median
tmp.seg$value <- log(tmp.seg$icopy/2,2);
tmp.seg[is.nan(tmp.seg$value),'value'] <- NA;

tmp.seg <- tmp.seg[c('sample', 'space', 'start', 'end', 'value')];

output.file = paste(date.value, "_", opt$output, ".seg", sep = "");
write.table(
	tmp.seg,
	output.file,
	sep = "\t",
	col.names = FALSE,
	row.names = FALSE,
	quote = FALSE
	);

segments <- scaleSegments(
	t.ar.seg,
	ePP
	);

segments <- annotateSegments(
	segments,
	ePP
	);

segmentsout <- segments; #[segments$mask == FALSE,];
rownames(segmentsout) <- 1:nrow(segmentsout);
segmentsout <- segmentsout[c('chrom', 'start.pos', 'end.pos', 'imean')];
segmentsout$median <- log(segmentsout$imean/2,2);
segmentsout <- segmentsout[c('chrom', 'start.pos', 'end.pos', 'median')];
colnames(segmentsout) <- c('chr', 'start', 'end', 'median');

output.file = paste(date.value, "_", opt$output, ".segments", sep = "");
write.table(
	segmentsout,
	output.file,
	sep = "\t",
	col.names = NA,
	row.names = TRUE,
	quote = FALSE
	);

plotSegment(
	tcs,
	segments,
	ar,
	file = paste(date.value, "segments_page%1d", sep = "_"),
	device = "png",
	width = 960,
	height = 1320,
	cex.axis = 2,
	cex.main = 2,
	cex.lab = 2,
	type = "cairo",
	tlwd = 8
	);

# XY chromosomes can't be annotated due to lack of ar
segmentsXY <- scaleSegments(
	t.seg,
	ePP
	);

plotSegment(
	tcs,
	segmentsXY,
	ar = NULL,
	file = paste(date.value, "segments_pageXY", sep = "_"),
	device = "png",
	width = 960,
	height = 1320,
	cex.axis = 2,
	cex.main = 2,
	cex.lab = 2,
	type = "cairo",
	chr = c("chrX","chrY"),
	tlwd = 8
	);

## OUTPUTS RESULTS FOR 2st solution
prepCN(12,1,NULL);
ePP <- ePeakPos(
	S = localSolutions[2,'S'],
	t = c(localSolutions[2,'N'], localSolutions[2,'T1']),
	cn = cn
	);

tcs <- scaleReadCounts( tc , ePP );

tmp.seg <- subset(
	as.data.frame(tcs),
	select = c(space,start,end,icopy)
	);
tmp.seg$sample <- opt$output;
tmp.seg <- tmp.seg[c('sample', 'space', 'start', 'end', 'icopy')];

# convert to median
tmp.seg$value <- log(tmp.seg$icopy/2,2);
tmp.seg[is.nan(tmp.seg$value),'value'] <- NA;

tmp.seg <- tmp.seg[c('sample', 'space', 'start', 'end', 'value')];

output.file = paste(date.value, "_", opt$output, "_2nd_solution.seg", sep = "");
write.table(
	tmp.seg,
	output.file,
	sep = "\t",
	col.names = FALSE,
	row.names = FALSE,
	quote = FALSE
	);

segments <- scaleSegments(
	t.ar.seg,
	ePP
	);

segments <- annotateSegments(
	segments,
	ePP
	);

segmentsout <- segments; #[segments$mask == FALSE,];
rownames(segmentsout) <- 1:nrow(segmentsout);
segmentsout <- segmentsout[c('chrom', 'start.pos', 'end.pos', 'imean')];
segmentsout$median <- log(segmentsout$imean/2,2);
segmentsout <- segmentsout[c('chrom', 'start.pos', 'end.pos', 'median')];
colnames(segmentsout) <- c('chr', 'start', 'end', 'median');

output.file = paste(date.value, "_", opt$output, "_2nd_solution.segments", sep = "");
write.table(
	segmentsout,
	output.file,
	sep = "\t",
	col.names = NA,
	row.names = TRUE,
	quote = FALSE
	);

plotSegment(
	tcs,
	segments,
	ar,
	file = paste(date.value, "2nd_solution", "segments_page%1d", sep = "_"),
	device = "png",
	width = 960,
	height = 1320,
	cex.axis = 2,
	cex.main = 2,
	cex.lab = 2,
	type = "cairo",
	tlwd = 8
	);

# XY chromosomes can't be annotated due to lack of ar
segmentsXY <- scaleSegments(
	t.seg,
	ePP
	);

plotSegment(
	tcs,
	segmentsXY,
	ar = NULL,
	file = paste(date.value, "2nd_solution", "segments_pageXY", sep = "_"),
	device = "png",
	width = 960,
	height = 1320,
	cex.axis = 2,
	cex.main = 2,
	cex.lab = 2,
	type = "cairo",
	chr = c("chrX","chrY"),
	tlwd = 8
	);

print("Processing complete");

### WRITE SESSION INFO #############################################################################
filename <- paste(date.value, "_celluloid_v0.11_session_info.txt", sep="");
sink(file = filename);
print(sessionInfo());
sink();

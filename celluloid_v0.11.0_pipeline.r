# Usage:
# Rscript celluloid_v0.11.0_pipeline.r  tumour.wig normal.wig gc.wig map.wig AR.txt samplename

today<- format( Sys.time() ,"%Y%m%d" )

.libPaths(new="/.mounts/labs/PCSI/R/3.3/")		# hardcoded PCSI library location
library(celluloid)

rdapipe<- paste( "RdaPipeline_",today,sep="") 

if( !file.exists( rdapipe  ) ){ 
  system( paste("mkdir", rdapipe ) )
}

########################################
#### saving the seed for reproducibility
######################################## 
invisible(runif(1)) 
save( .Random.seed, file=paste( rdapipe,"/Random.seed",format( Sys.time() ,"%Y%m%d%H%M%S" ) , sep="" ) ) 

#######

argv<-commandArgs( trailing=TRUE )

argc<-length(argv)
if( argc!=6 ){ 
  stop("Usage: Rscript pipeline.r tumourWigFile normalWigFile gcWigFile mapWigFile arFile samplename") }

##############
# prepare data
##############

tumourWigFile<-argv[1]
normalWigFile<-argv[2]
gcWigFile<- argv[3]
mapWigFile<- argv[4]
arFile<- argv[5]
name<-argv[6]
if( !file.exists(paste( rdapipe,"/nc.rda", sep="") ) ){
  # load and segment normal data
  n<- wigsToRangedData( normalWigFile, gcWigFile, mapWigFile )
  # gc correction
  nc<-gcCorrect(n) 
  n.seg <- segmentSeqData( nc, k=50 , maskmap = 0.8 )
  save( n.seg, file=paste( rdapipe,"/n.seg.rda", sep=""))
  save( nc, file=paste( rdapipe,"/nc.rda", sep=""))
} else {
  load(paste( rdapipe,"/n.seg.rda", sep=""))
  load(paste( rdapipe,"/nc.rda", sep=""))
}


if( !file.exists(paste( rdapipe,"/tc.rda", sep="")) ){
  # load and segment tumor data
  t <- wigsToRangedData( tumourWigFile, gcWigFile, mapWigFile )
  tc<-gcCorrect( t, sampletype="tumor" )
  t.seg <- segmentSeqData( tc , k=50  )
  save( t.seg, file=paste( rdapipe,"/t.seg.rda", sep=""))
  save( tc, file=paste( rdapipe,"/tc.rda", sep=""))
} else {
  load(paste( rdapipe,"/t.seg.rda", sep=""))
  load(paste( rdapipe,"/tc.rda", sep=""))
}


# identify segments not "normal" in normal
if(!file.exists(paste( rdapipe,"/t.seg.mask.rda", sep="")) ){
  t.n.seg<-intersectSegments( t.seg, n.seg )
  t.n.seg<-t.n.seg[ !is.na( t.n.seg$mean) & !is.na( t.n.seg$mean.1),]
  sel<-n.seg$end.pos-n.seg$start.pos > 150000 & n.seg$meanmap>.8 
  bp<-boxplot( n.seg$mean[sel], range=3, plot=F  )
  nrange<-c(bp$stats[1,1], bp$stats[5,1] )
  mask<- t.n.seg$mean.1>nrange[2] | t.n.seg$mean.1<nrange[1]
  t.seg.mask<-t.n.seg[,1:8]
  t.seg.mask$mask<- mask
  save(t.seg.mask, file=paste( rdapipe,"/t.seg.mask.rda", sep="")) 
} else {
  load(paste( rdapipe,"/t.seg.mask.rda", sep=""))
}

if( !file.exists(paste( rdapipe,"/t.ar.seg.rda" , sep="")) ){
  # reading allelic ratio file and segmenting
  ar<-read.table(arFile, head=T )
  save( ar, file=paste( rdapipe,"/ar.rda", sep="")) 
  ar.seg<- segmentAR( ar, tc ) 
  save( ar.seg, file=paste( rdapipe,"/ar.seg.rda", sep=""))
  t.ar.seg <- intersectSegments( t.seg.mask, ar.seg  )
  t.ar.seg<-t.ar.seg[ !apply( is.na( t.ar.seg[,c("mean","meanar")] ), 1, any),]
  t.ar.seg <-arInSeg( t.ar.seg, ar,  tumourrangedata=tc , minhet = 50 )
  save( t.ar.seg, file=paste( rdapipe,"/t.ar.seg.rda", sep=""))
} else {
  load(paste( rdapipe,"/ar.rda", sep=""))
  load(paste( rdapipe,"/t.ar.seg.rda", sep=""))
}


if( !file.exists(paste( rdapipe,"/copyAr.rda", sep="")) ){
  # creates the object used to draw a contour plot. 
  mask<- t.ar.seg$mask | is.na(t.ar.seg$mask)
  copyAr<-  prepCopyAr( t.ar.seg[ !mask ,], ar,  tc  )
  save(copyAr, file=paste( rdapipe,"/copyAr.rda", sep=""))
} else {
  load(paste( rdapipe,"/copyAr.rda", sep=""))
}

if( !file.exists( paste( rdapipe,"/cntr.rda", sep="")) ){
  cntr<-showTumourProfile( copyAr , flatten=.25 , nlev=20 , noise= 0.01 , 
                           maxPoints=200000 , plot=F  )
  save(cntr, file=paste( rdapipe,"/cntr.rda", sep=""))
} else {
  load( paste( rdapipe,"/cntr.rda", sep="") )
}


#################
# analysis starts
#################

# finding LOH curve; this may fail and require manual intervention...
Sn<- estimateLOHcurve(t.ar.seg)

if( !file.exists( paste( rdapipe,"/lm.rda", sep=""))  | !file.exists( paste(  rdapipe,"/paramSpace.rda", sep="" ) ) ){
  
  
  # grid search
  sel <- t.ar.seg$mean < max( cntr$x ) & !is.na( t.ar.seg$mean ) & !is.na( t.ar.seg$p )
  segmentsubset <-t.ar.seg[sel,] 
  
  # using optim
  
  # upper bound for %normal cell set to Sn/.25, set to 1 if you think your sample's ploidy can 
  # be above 8 (2/.25)
  
  lm<- coverParamSpace( segments=segmentsubset, control=list( maxit=1000 ), Sn=Sn ,
                        maxc=12,optimFct=1 , nrep=50, 
                        lowerF=c(0), upperF=c( Sn/.25 ), addToParamSpace=T  )
  
  save(lm, file=paste(  rdapipe,"/lm.rda", sep="" ) )
  save(paramSpace, file=paste(  rdapipe,"/paramSpace.rda", sep="" ) )
  
} else {
  load( paste(  rdapipe,"/lm.rda", sep="" ) )
  load( paste(  rdapipe,"/paramSpace.rda", sep="" ) )
  prepCN( 100,1,NULL )
}


localSolutions <-getLocalSolutions(lm, max=TRUE) 
# only solutions for which 10% or more of the genome was captured are in the output, and at most 5
localSolutions<-head(localSolutions[ localSolutions$value>0.10 , ],5)

pdf( paste( "contour_solutions_",name,"_", today,".pdf", sep="") , height=12, width=16 )
par( mfrow=c(3,2) )
plot(  paramSpace[,3] , 100*paramSpace[,1] , pch=19 , ylim=c(0,100), xlab="%normal", ylab="percent captured by model", main=name)
for( i in 1:nrow(localSolutions) ){
  r<-as.numeric(rownames(localSolutions)[i] )
  image( cntr, col=terrain.colors(50))
  contour(cntr, nlev=8, add=T )
  sel<-t.ar.seg$size>1000000 & !t.ar.seg$mask & t.ar.seg$meanmap>.9
  le<- t.ar.seg$end.pos[sel]-t.ar.seg$start.pos[sel]
  cxcut<- as.integer( cut( le, c(100000,1000000,5000000,10000000,20000000,50000000,Inf) ) )/3
  points( x<-t.ar.seg$mean[sel], y<-t.ar.seg$p[sel],  pch=21 , col="blue", lwd=3 , cex=cxcut  )
  points( t.ar.seg$mean[sel], t.ar.seg$p[sel],  pch=19 ,  col="white" , cex= cxcut  - .5 )
  points( t.ar.seg$mean[sel], 1-t.ar.seg$p[sel],  pch=21 ,  col="blue", lwd=3  , cex=cxcut  )
  points( t.ar.seg$mean[sel], 1-t.ar.seg$p[sel],  pch=19 , col="white" , cex=cxcut -.5 )
  xxx <- seq( Sn, 2, .01 ) 
  points( xxx , arloh<- ARloh( xxx , 1 , Sn ) , type='l' , lwd=2  )
  points( xxx , 1-arloh , type='l' , lwd=2  )
  
  epp<-  plotModelPeaks(lm[[r]]$par, selectedPoints=NULL,cn=cn, epcol="red",epcex=1,eplwd=3 , addlabels=F )
  legend( 0,1, paste(i,"/ ", floor(1000*( lm[[r]]$value) )/10, sep="" ), bg="white"  )
  
  
}
dev.off() 



# this creates a directory for each solutions, and outputs all segments and segments plots

localSolutions$Ploidy<- ( 2/localSolutions$S - 2*localSolutions$N )/(1-localSolutions$N )

for( i in 1:nrow(localSolutions) ){
  
  soldir<- paste(  "solution",i, sep="") 		# removed the date -RD
  if( ! file.exists( soldir ) ) 
    system(paste("mkdir", soldir ) )
  
  r<-as.numeric(rownames(localSolutions)[i] )
  prepCN(12,1,NULL)
  
  ePP<-ePeakPos( par= lm[[r]]$par , cn=cn  )
  tcs<- scaleReadCounts( tc , ePP )
  segments<-scaleSegments(t.ar.seg ,  ePP )
  segments<-annotateSegments(segments, ePP)
  plotSegment( tcs,segments, ar , file=paste( soldir, "/segmentPlots_",name,"_page%1d",sep="" ) ,device="png",
               width=960,height=1320, cex.axis=2, cex.main=2, cex.lab=2, 
               type="cairo" , tlwd=8 ) 
  
  # XY chromosomes can't be annotated due to lack of ar
  segmentsXY<-scaleSegments(t.seg ,  ePP )
  plotSegment( tcs,segmentsXY, ar=NULL , file=paste( soldir, "/segmentPlots_",name,"_pageXY", sep="") ,device="png",
               width=960,height=1320, cex.axis=2, cex.main=2, cex.lab=2, 
               type="cairo" , chr=c("chrX","chrY") , tlwd=8) 
  
  # adding empty columns
  selXY<- segmentsXY$chrom=="chrX" |   segmentsXY$chrom=="chrY" | segmentsXY$chrom=="chr23" | segmentsXY$chrom=="chr24" 
  segmentsXY <- segmentsXY[selXY,] 
  segmentsXY$size <- segmentsXY$end.pos - segmentsXY$start.pos 
  segmentsXY$mask <- NA
  segmentsXY$meanar <- NA
  segmentsXY$p <- NA
  segmentsXY$labels<- NA
  segmentsXY<-segmentsXY[,c("chrom","start.pos","end.pos","size","arm","sampleID","mean","meanmap","mask","meanar","p","imean", "labels" )]   # changed order of columns so the output is a valid bed file -RD
  
  segments<-segments[,c("chrom","start.pos","end.pos","size","arm","sampleID","mean","meanmap","mask","meanar","p","imean", "labels" )]
  segments<- rbind( segments, segmentsXY )
  segments$sampleID<-name 
  
  write.table( segments,  file=paste( soldir, "/segments_",name,".txt", sep="" ), quo=F, col=TRUE, row=F, sep="\t" )   # added tab delimiter -RD
  
  write.table( localSolutions[i,] ,  file=paste( soldir, "/parameters_",name,".txt", sep="" ), quo=F, col=TRUE, row=F )
  
  tmp.seg <- subset( as.data.frame(tcs), select = c(space,start,end,icopy))
  tmp.seg$sample<-   name
  tmp.seg <- tmp.seg[c('sample', 'space', 'start', 'end', 'icopy')];
  write.table(tmp.seg, file=paste( soldir, "/scaledReadCounts_",name,".seg", sep="" ), sep = "\t",
    col.names = FALSE,row.names = FALSE,quote = FALSE)
   
  # plot contour as png for each solution -RD
  png( filename=paste( soldir, "/contour_",name,".png", sep="") , height=400, width=800 , type="cairo")

  r<-as.numeric(rownames(localSolutions)[i] )
  image( cntr, col=terrain.colors(50))
  contour(cntr, nlev=8, add=T )
  sel<-t.ar.seg$size>1000000 & !t.ar.seg$mask & t.ar.seg$meanmap>.9
  le<- t.ar.seg$end.pos[sel]-t.ar.seg$start.pos[sel]
  cxcut<- as.integer( cut( le, c(100000,1000000,5000000,10000000,20000000,50000000,Inf) ) )/3
  points( x<-t.ar.seg$mean[sel], y<-t.ar.seg$p[sel],  pch=21 , col="blue", lwd=3 , cex=cxcut  )
  points( t.ar.seg$mean[sel], t.ar.seg$p[sel],  pch=19 ,  col="white" , cex= cxcut  - .5 )
  points( t.ar.seg$mean[sel], 1-t.ar.seg$p[sel],  pch=21 ,  col="blue", lwd=3  , cex=cxcut  )
  points( t.ar.seg$mean[sel], 1-t.ar.seg$p[sel],  pch=19 , col="white" , cex=cxcut -.5 )
  xxx <- seq( Sn, 2, .01 )
  points( xxx , arloh<- ARloh( xxx , 1 , Sn ) , type='l' , lwd=2  )
  points( xxx , 1-arloh , type='l' , lwd=2  )

  epp<-  plotModelPeaks(lm[[r]]$par, selectedPoints=NULL,cn=cn, epcol="red",epcex=1,eplwd=3 , addlabels=F )
  legend( 0,1, paste(i,"/ ", floor(1000*( lm[[r]]$value) )/10, sep="" ), bg="white"  )

  dev.off()


}


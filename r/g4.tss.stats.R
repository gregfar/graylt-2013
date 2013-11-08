args <- commandArgs(trailingOnly = T)

# check to see if the plyr package is installed.
if("plyr" %in% rownames(installed.packages()) == F) {
  install.packages("plyr")
} else {
  library(plyr)
}

peak.classes <- c("character","numeric","numeric","character","numeric","factor","numeric")
bed.names <- c("chr","start","end","name","score","strand")
bed.classes <- c("character","numeric","numeric","character","numeric","factor")

# read the TSS 1kb region table, XPB summit, XPD summit, and ORB expression tables
tss <- read.table(args[1],col.names=bed.names,colClasses=bed.classes)
xpb.peaks <- read.table(args[2],header=T,colClasses=peak.classes)
xpd.peaks <- read.table(args[3],header=T,colClasses=peak.classes)
xpb.regions <- read.table(args[4],header=T,colClasses=peak.classes)
xpd.regions <- read.table(args[5],header=T,colClasses=peak.classes)
xpb.summits <- read.table(args[6],col.names=bed.names,colClasses=bed.classes,skip=1)
xpd.summits <- read.table(args[7],col.names=bed.names,colClasses=bed.classes,skip=1)
g4 <- read.table(args[8],header=T,
                 colClasses=c("character","numeric","numeric","character","numeric","factor","numeric","numeric","numeric","numeric","numeric"))
# For G4 motifs, I need to do the merge in Perl, because there are too many regions for this
# to take a reasonable amount of time.
cat.g4 <- read.table(args[9],col.names=c("chr","start","end","name","score","strand"),
                 colClasses=c("character","numeric","numeric","character","numeric","factor"))
genome <- read.table(args[10],col.names=c("chr","start","end","score"),
                     colClasses=c("character","numeric","numeric","numeric"))
results <- args[11]

# this function will flag overlapping regions with a 1
bed.flag <- function(bed1,bed2) {
  
  result <- cbind(bed1,overlap=0)
  
  for (i in 1:nrow(bed1)) {
    
    sub2 <- bed2[bed2$chr == bed1$chr[i],]
    
    compare <- sub2[sub2$start <= bed1$end[i] & sub2$end >= bed1$start[i],]
    
    if(nrow(compare) > 0) {
      result$overlap[i] <- 1
    }
    
  }
  
  return(result)
  
}

# This function will merge overlapping bed regions into a single chr start end entry
# which is useful for calculating coverage.
bed.merge <- function(bed) {
    
  result <- data.frame(chr=character(0),start=numeric(0),end=numeric(0),name=character(0),score=numeric(0),strand=character(0))
    
    c.bed <- bed
    
    for (i in 1:nrow(bed)) {
        
        overlap <- c.bed[c.bed$chr == bed$chr[i] & c.bed$start <= bed$end[i] & c.bed$end >= bed$start[i],]
        
        if(nrow(overlap) > 0) {
            overlap$start[1] <- min(overlap$start)
            overlap$end[1] <- max(overlap$end)
            result <- rbind(result,overlap[1,])
        } else {
            result <- rbind(result,bed[i,])
        }
        
    }
    
    return(result)
    
}


# This function will return just the overlapping portions of bed1 that overlap bed2
bed.overlap <- function(bed1,bed2) {
  
  result <- data.frame(chr=character(0),start=numeric(0),end=numeric(0),name=character(0),score=numeric(0),strand=character(0))
  
  for (i in 1:nrow(bed1)) {
    
    sub2 <- bed2[bed2$chr == bed1$chr[i],]
    
    compare <- sub2[sub2$start <= bed1$end[i] & sub2$end >= bed1$start[i],]
    
    if(nrow(compare) > 0) {

      if(bed1$start[i] < min(compare$start)) {
          bed1$start[i] <- min(compare$start)
      }
      
      if(bed1$end[i] > max(compare$end)) {
          bed1$end[i] <- max(compare$end)
      }
      
      result <- rbind(result,bed1[i,])
      
    }
    
  }
  
  return(result)
  
}

bed.coverage <- function(bed) {
  
  result <- sum(bed$end - bed$start)
  
  return(result)
  
}

# Flag TSS for XPB peaks, regions, and G4s
xpb.summits <- bed.flag(xpb.summits,tss)
names(xpb.summits)[7] <- "tss"
xpb.peaks <- cbind(xpb.peaks,tss=0)
xpb.peaks$tss[xpb.peaks$name %in% xpb.summits$name[xpb.summits$tss == 1]] <- 1
write.table(xpb.peaks,paste(results,"/xpb_peaks_g4_tss.txt",sep=""),sep="\t",quote=F,row.names=F)

xpb.regions <- cbind(xpb.regions,tss=0)
xpb.regions$tss[xpb.regions$name %in% xpb.summits$name[xpb.summits$tss == 1]] <- 1
write.table(xpb.regions,paste(results,"/xpb_regions_g4_tss.txt",sep=""),sep="\t",quote=F,row.names=F)

xpd.summits <- bed.flag(xpd.summits,tss)
names(xpd.summits)[7] <- "tss"
xpd.peaks <- cbind(xpd.peaks,tss=0)
xpd.peaks$tss[xpd.peaks$name %in% xpd.summits$name[xpd.summits$tss == 1]] <- 1
write.table(xpd.peaks,paste(results,"/xpd_peaks_g4_tss.txt",sep=""),sep="\t",quote=F,row.names=F)

xpd.regions <- cbind(xpd.regions,tss=0)
xpd.regions$tss[xpd.regions$name %in% xpd.summits$name[xpd.summits$tss == 1]] <- 1
write.table(xpd.regions,paste(results,"/xpd_regions_g4_tss.txt",sep=""),sep="\t",quote=F,row.names=F)

# Flagging of G4s for TSS, XPB, and XPD was done previously with perl, because the G4 table is so large.

# Count the number of overlaps between peaks, tss, and g4
xpb.peak.counts <- count(xpb.peaks,c("tss","g4"))
xpb.peak.counts <- cbind(xpb.peak.counts,percent=0)
xpb.peak.counts$percent <- round((xpb.peak.counts$freq/sum(xpb.peak.counts$freq))*100,2)
write.table(xpb.peak.counts,paste(results,"/xpb_peaks_g4_tss_counts.txt",sep=""),sep="\t",quote=F,row.names=F)

xpd.peak.counts <- count(xpd.peaks,c("tss","g4"))
xpd.peak.counts <- cbind(xpd.peak.counts,percent=0)
xpd.peak.counts$percent <- round((xpd.peak.counts$freq/sum(xpd.peak.counts$freq))*100,2)
write.table(xpd.peak.counts,paste(results,"/xpd_peaks_g4_tss_counts.txt",sep=""),sep="\t",quote=F,row.names=F)

xpb.region.counts <- count(xpb.regions,c("tss","g4"))
xpb.region.counts <- cbind(xpb.region.counts,percent=0)
xpb.region.counts$percent <- round((xpb.region.counts$freq/sum(xpb.region.counts$freq))*100,2)
write.table(xpb.region.counts,paste(results,"/xpb_regions_g4_tss_counts.txt",sep=""),sep="\t",quote=F,row.names=F)

xpd.region.counts <- count(xpd.regions,c("tss","g4"))
xpd.region.counts <- cbind(xpd.region.counts,percent=0)
xpd.region.counts$percent <- round((xpd.region.counts$freq/sum(xpd.region.counts$freq))*100,2)
write.table(xpd.region.counts,paste(results,"/xpd_regions_g4_tss_counts.txt",sep=""),sep="\t",quote=F,row.names=F)

g4.counts <- count(g4,c("tss","xpb.peaks","xpd.peaks"))
g4.counts <- cbind(g4.counts,percent=0)
g4.counts$percent <- round((g4.counts$freq/sum(g4.counts$freq))*100,2)
write.table(g4.counts,paste(results,"/g4_xpb_xpd_tss_counts.txt",sep=""),sep="\t",quote=F,row.names=F)

# merge the overlapping tss regions to calculate coverage
cat.tss <- bed.merge(tss)

# get just the xpb, xpd, and g4 regions that overlap TSS for coverage calculations
xpb.peaks.tss <- bed.overlap(xpb.peaks[xpb.peaks$tss == 1,],cat.tss)
xpd.peaks.tss <- bed.overlap(xpd.peaks[xpd.peaks$tss == 1,],cat.tss)
xpb.regions.tss <- bed.overlap(xpb.regions[xpb.regions$tss == 1,],cat.tss)
xpd.regions.tss <- bed.overlap(xpd.regions[xpd.regions$tss == 1,],cat.tss)
g4.tss <- bed.overlap(g4[g4$tss == 1,],cat.tss)

# calculate coverage of regions
genome.cov <- bed.coverage(genome)
tss.cov <- bed.coverage(cat.tss)
not.tss.cov <- genome.cov - tss.cov

xpb.peak.cov <- bed.coverage(xpb.peaks)
xpb.peak.tss.cov <- bed.coverage(xpb.peaks.tss)
xpb.peak.not.tss.cov <- xpb.peak.cov - xpb.peak.tss.cov

xpd.peak.cov <- bed.coverage(xpd.peaks)
xpd.peak.tss.cov <- bed.coverage(xpd.peaks.tss)
xpd.peak.not.tss.cov <- xpd.peak.cov - xpd.peak.tss.cov

xpb.region.cov <- bed.coverage(xpb.regions)
xpb.region.tss.cov <- bed.coverage(xpb.regions.tss)
xpb.region.not.tss.cov <- xpb.region.cov - xpb.region.tss.cov

xpd.region.cov <- bed.coverage(xpd.regions)
xpd.region.tss.cov <- bed.coverage(xpd.regions.tss)
xpd.region.not.tss.cov <- xpd.region.cov - xpd.region.tss.cov

g4.cov <- bed.coverage(cat.g4)
g4.tss.cov <- bed.coverage(g4.tss)
g4.not.tss.cov <- g4.cov - g4.tss.cov

coverage.results <- data.frame(region=c("all","tss","not.tss"),genome=c(genome.cov,tss.cov,not.tss.cov),
                      xpb.peak=c(xpb.peak.cov,xpb.peak.tss.cov,xpb.peak.not.tss.cov),
                      xpd.peak=c(xpd.peak.cov,xpd.peak.tss.cov,xpd.peak.not.tss.cov),
                      xpb.region=c(xpb.region.cov,xpb.region.tss.cov,xpb.region.not.tss.cov),
                      xpd.region=c(xpd.region.cov,xpd.region.tss.cov,xpd.region.not.tss.cov),
                      g4=c(g4.cov,g4.tss.cov,g4.not.tss.cov))
write.table(coverage.results,paste(results,"/tss_coverage_results.txt",sep=""),sep="\t",quote=F,row.names=F)

percent.coverage.results <- coverage.results
percent.coverage.results[,2:7] <- percent.coverage.results[,2:7] / coverage.results$genome
percent.coverage.results[,2:7] <- round(percent.coverage.results[,2:7] * 100,2)
write.table(percent.coverage.results,paste(results,"/tss_percent_coverage_results.txt",sep=""),sep="\t",quote=F,row.names=F)

args <- commandArgs(trailingOnly = T)
i <- args[5]

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

xpb.summits <- xpb.peaks
xpb.summits$start <- xpb.peaks$start + round((xpb.peaks$end - xpb.peaks$start)/2)
xpb.summits$end <- xpb.summits$start + 1

xpd.peaks <- read.table(args[3],header=T,colClasses=peak.classes)

xpd.summits <- xpd.peaks
xpd.summits$start <- xpd.peaks$start + round((xpd.peaks$end - xpd.peaks$start)/2)
xpd.summits$end <- xpd.summits$start + 1

results <- args[4]

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


# Flag TSS for XPB peaks, regions, and G4s
xpb.summits <- bed.flag(xpb.summits,tss)
names(xpb.summits)[7] <- "tss"
xpb.peaks <- cbind(xpb.peaks,tss=xpb.summits$tss)
write.table(xpb.peaks,paste(results,"/xpb_peaks_g4_tss_",i,".txt",sep=""),sep="\t",quote=F,row.names=F)

xpd.summits <- bed.flag(xpd.summits,tss)
names(xpd.summits)[7] <- "tss"
xpd.peaks <- cbind(xpd.peaks,tss=xpd.summits$tss)
write.table(xpd.peaks,paste(results,"/xpd_peaks_g4_tss_",i,".txt",sep=""),sep="\t",quote=F,row.names=F)


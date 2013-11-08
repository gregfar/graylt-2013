# check to see if the plyr package is installed.
if("plyr" %in% rownames(installed.packages()) == F) {
  install.packages("plyr")
} else {
  library(plyr)
}

# check to see if the ggplot2 package is installed.
if("ggplot2" %in% rownames(installed.packages()) == F) {
  install.packages("ggplot2")
} else {
  library(ggplot2)
}

# read the TSS 1kb region table, XPB summit, XPD summit, and ORB expression tables
refseq.tss <- read.table("tss_1kb.bed",col.names=c("chr","start","end","name","score","strand"),
                         colClasses=c("character","numeric","numeric","character","numeric","factor"))
xpb <- read.table("xpb_summits.bed",col.names=c("chr","start","end","name","score","strand"),
                  colClasses=c("character","numeric","numeric","character","numeric","factor"),skip=1)
xpd <- read.table("xpd_summits.bed",col.names=c("chr","start","end","name","score","strand"),
                  colClasses=c("character","numeric","numeric","character","numeric","factor"),skip=1)
orb <- read.table("orb_gene_table.txt",header=T,sep="\t",
                  colClasses=c("character","character","character","numeric","numeric","numeric","numeric","numeric","character"))
xpb.g4 <- read.table("xpb_g4_peaks_summits.bed",col.names=c("chr","start","end","name","score","strand"),
                  colClasses=c("character","numeric","numeric","character","numeric","factor"))
xpd.g4 <- read.table("xpd_g4_peaks_summits.bed",col.names=c("chr","start","end","name","score","strand"),
                  colClasses=c("character","numeric","numeric","character","numeric","factor"))

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

# this function will flag overlapping regions with the expression category
bed.flag.ec <- function(bed1,bed2) {
  
  result <- cbind(bed1,exp.cat=0)
  
  for (i in 1:nrow(bed1)) {
    
    sub2 <- bed2[bed2$chr == bed1$chr[i],]
    
    compare <- sub2[sub2$start <= bed1$end[i] & sub2$end >= bed1$start[i],]

    if(nrow(compare) > 0) {
      result$exp.cat[i] <- min(compare$exp.cat)
    }
    
  }
  
  return(result)
  
}

# this function will join overlapping bed regions
bed.join <- function(bed1,bed2) {
  
  result <- cbind(bed1[1,],bed2[1,])
  
  for (i in 1:nrow(bed1)) {
    
    sub2 <- bed2[bed2$chr == bed1$chr[i],]
    
    compare <- sub2[sub2$start <= bed1$end[i] & sub2$end >= bed1$start[i],]
    
    if(nrow(compare) > 0) {
      for (j in 1:nrow(compare)) {
        result.row <- cbind(bed1[i,],compare[j,])
        result <- rbind(result,result.row)
      }
    }
    
  }
  
  result <- result[-1,]
  
  return(result)
  
}

xpb <- cbind(xpb,g4=0)
xpb$g4[xpb$name %in% xpb.g4$name] <- 1

xpd <- cbind(xpd,g4=0)
xpd$g4[xpd$name %in% xpd.g4$name] <- 1

## Figure 3B calculations
# get the mean expression scores from the ORB table, sort by score, and assign expression categories
orb.scores <- data.frame(name=orb$gene.symbol,exp=orb$mean.exp,stringsAsFactors=F)
orb.scores <- orb.scores[orb.scores$name != "---",]
orb.scores <- orb.scores[with(orb.scores,order(-exp)),]
exp.cat <- c(rep(1,2924),rep(2,2923),rep(3,2924),rep(4,2923),rep(5,2924))
orb.scores <- cbind(orb.scores,exp.cat=exp.cat)

# get unique rows of the refseq tss table
refseq.tss <- unique(refseq.tss)

# merge the ORB expression data with the refseq table by gene name
df <- merge(refseq.tss,orb.scores,all.x=T,by="name")
df$exp[is.na(df$exp)] <- 0
df$exp.cat[df$exp == 0] <- 6

# flag XPB and XPD peaks that overlap TSS regions
xpb <- bed.flag.ec(xpb,df)
xpd <- bed.flag.ec(xpd,df)

count.xpb <- count(xpb,c("g4","exp.cat"))
count.xpd <- count(xpd,c("g4","exp.cat"))

write.table(count.xpb,"xpb_g4_tss_binding.txt",sep="\t",row.names=F,quote=F)
write.table(count.xpd,"xpd_g4_tss_binding.txt",sep="\t",row.names=F,quote=F)

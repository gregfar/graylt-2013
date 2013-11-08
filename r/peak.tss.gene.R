args <- commandArgs(trailingOnly = T)
results <- args[5]

# check to see if the plyr package is installed.
if("plyr" %in% rownames(installed.packages()) == F) {
  install.packages("plyr")
} else {
  library(plyr)
}

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

bed.names <- c("chr","start","end","name","score","strand")
bed.classes <- c("character","numeric","numeric","character","numeric","character")

refseq <- read.table(args[1],sep="\t",
                     col.names=c("bin","acc","chr","strand","txn.start","txn.end","cds.start","cds.end","exon.count","exon.starts","exon.ends","score","name","cds.start.stat","cds.end.stat","exon.frames"),
                     colClasses=c("numeric","character","character","character","numeric","numeric","numeric","numeric","numeric","character","character","numeric","character","character","character","character"))

tss <- read.table(args[2],sep="\t",col.names=bed.names,colClasses=bed.classes)

plus <- data.frame(chr=tss$chr[tss$strand == "+"],start=tss$end[tss$strand == "+"],end=refseq$txn.end[refseq$strand == "+"],name=tss$name[tss$strand == "+"],score=tss$score[tss$strand == "+"],strand="+")
minus <- data.frame(chr=tss$chr[tss$strand == "-"],start=refseq$txn.start[refseq$strand == "-"],end=tss$start[tss$strand == "-"],name=tss$name[tss$strand == "-"],score=tss$score[tss$strand == "-"],strand="-")
no.tss <- rbind(plus,minus)

xpb.summits <- read.table(args[3],sep="\t",skip=1,col.names=bed.names,colClasses=bed.classes)
xpd.summits <- read.table(args[4],sep="\t",skip=1,col.names=bed.names,colClasses=bed.classes)

xpb.summits <- bed.flag(xpb.summits,tss)
names(xpb.summits)[7] <- "tss"
xpb.summits <- bed.flag(xpb.summits,no.tss)
names(xpb.summits)[8] <- "gene"

xpd.summits <- bed.flag(xpd.summits,tss)
names(xpd.summits)[7] <- "tss"
xpd.summits <- bed.flag(xpd.summits,no.tss)
names(xpd.summits)[8] <- "gene"

library("plyr")

xpb.counts <- count(xpb.summits,c("tss","gene"))
xpb.counts <- cbind(xpb.counts,percent=round(100*xpb.counts$freq/sum(xpb.counts$freq),2))
write.table(xpb.counts,paste(results,"/xpb_tss1kb_gene_counts.txt",sep=""),sep="\t",quote=F,row.names=F)

xpd.counts <- count(xpd.summits,c("tss","gene"))
xpd.counts <- cbind(xpd.counts,percent=round(100*xpd.counts$freq/sum(xpd.counts$freq),2))
write.table(xpd.counts,paste(results,"/xpd_tss1kb_gene_counts.txt",sep=""),sep="\t",quote=F,row.names=F)

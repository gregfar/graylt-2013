args <- c("03_peaks/XPB/xpb_peaks_gc.txt",
          "03_peaks/XPD/xpd_peaks_gc.txt",
          "10_results/xpb_peaks_g4_tss.txt",
          "10_results/xpd_peaks_g4_tss.txt",
          "10_results")

args <- commandArgs(trailingOnly = T)
results <- args[5]

b.gc <- read.table(args[1],sep="\t",col.names=c("name","gc"),stringsAsFactors=F)
d.gc <- read.table(args[2],sep="\t",col.names=c("name","gc"),stringsAsFactors=F)

bed.names <- c("chr","start","end","name","score","strand","g4","tss")

b.g4 <- read.table(args[3],sep="\t",col.names=bed.names,stringsAsFactors=F,skip=1)
d.g4 <- read.table(args[4],sep="\t",col.names=bed.names,stringsAsFactors=F,skip=1)

b.gc <- cbind(b.gc,g4=0,tss=0)
d.gc <- cbind(d.gc,g4=0,tss=0)

b.gc$g4[b.gc$name %in% b.g4$name[b.g4$g4 == 1] ] <- 1
b.gc$tss[b.gc$name %in% b.g4$name[b.g4$tss == 1] ] <- 1

d.gc$g4[d.gc$name %in% d.g4$name[d.g4$g4 == 1] ] <- 1
d.gc$tss[d.gc$name %in% d.g4$name[d.g4$tss == 1] ] <- 1

library(plyr)

b.result <- count(b.g4,c("g4","tss"))
d.result <- count(d.g4,c("g4","tss"))

calculate.gc <- function(rt,gc) {
  
  result <- cbind(rt,gc=0)

  for(i in 1:nrow(result)) {
    
    result$gc[i] <- mean(gc$gc[gc$g4 == rt$g4[i] & gc$tss == rt$tss[i] ])
    
  }
  
  return(result)
}

b.result <- calculate.gc(b.result,b.gc)
d.result <- calculate.gc(d.result,d.gc)

write.table(b.result,paste(results,"/xpb_g4_tss_gc_summary.txt",sep=""),row.names=F,quote=F,sep="\t")
write.table(d.result,paste(results,"/xpd_g4_tss_gc_summary.txt",sep=""),row.names=F,quote=F,sep="\t")

library(ggplot2)
source("XX_tools/r/theme_pub.R")
qplot(interaction(b.gc$g4,b.gc$tss),b.gc$gc,geom="boxplot") + theme_pub()
ggsave(paste(results,"/xpb_g4_tss.tiff",sep=""),height=3,width=3)
qplot(interaction(d.gc$g4,d.gc$tss),d.gc$gc,geom="boxplot") + theme_pub()
ggsave(paste(results,"/xpd_g4_tss.tiff",sep=""),height=3,width=3)

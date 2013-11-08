args <- commandArgs(trailingOnly = T)
results <- args[3]

gc <- read.table(args[1],sep="\t",col.names=c("name","gc"),stringsAsFactors=F)

bed.names <- c("chr","start","end","name","score","strand","g4","tss")

g4 <- read.table(args[2],sep="\t",col.names=bed.names,stringsAsFactors=F,skip=1)

gc <- cbind(gc,g4=0,tss=0)

gc$g4 <- g4$g4
gc$tss <- g4$tss

gc <- gc[gc$gc != "N",]
gc$gc <- as.numeric(gc$gc)

library(plyr)

result <- count(gc,c("g4","tss"))

calculate.gc <- function(rt,gc) {
  
  result <- cbind(rt,gc=0)

  for(i in 1:nrow(result)) {
    
    result$gc[i] <- mean(gc$gc[gc$g4 == rt$g4[i] & gc$tss == rt$tss[i] ])
    
  }
  
  return(result)
}

result <- calculate.gc(result,gc)

write.table(result,paste(results,"/random_gc_summary.txt",sep=""),row.names=F,quote=F,sep="\t")

#xpb.gc <- read.table(args[3],header=F,stringsAsFactors=F,col.names=c("name","gc"))
#xpd.gc <- read.table(args[4],header=F,stringsAsFactors=F,col.names=c("name","gc"))



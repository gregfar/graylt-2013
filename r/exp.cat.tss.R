args <- commandArgs(trailingOnly = T)
output <- args[3]

refseq.tss <- read.table(args[1],col.names=c("chr","start","end","name","score","strand"),
                         colClasses=c("character","numeric","numeric","character","numeric","factor"))
summary <- read.table(args[2],header=T,
                      colClasses=c("character","numeric","numeric","numeric","numeric","character"))

tss.summary <- merge(refseq.tss,summary,by="name")
tss.summary <- tss.summary[with(tss.summary,order(chr,start)),]

write.table(tss.summary[tss.summary$exp.cat == 1,][,c(2,3,4,1,5,6)],paste(output,"/tss_1kb_cat1.bed",sep=""),sep="\t",quote=F,col.names=F,row.names=F)
write.table(tss.summary[tss.summary$exp.cat == 2,][,c(2,3,4,1,5,6)],paste(output,"/tss_1kb_cat2.bed",sep=""),sep="\t",quote=F,col.names=F,row.names=F)
write.table(tss.summary[tss.summary$exp.cat == 3,][,c(2,3,4,1,5,6)],paste(output,"/tss_1kb_cat3.bed",sep=""),sep="\t",quote=F,col.names=F,row.names=F)
write.table(tss.summary[tss.summary$exp.cat == 4,][,c(2,3,4,1,5,6)],paste(output,"/tss_1kb_cat4.bed",sep=""),sep="\t",quote=F,col.names=F,row.names=F)
write.table(tss.summary[tss.summary$exp.cat == 5,][,c(2,3,4,1,5,6)],paste(output,"/tss_1kb_cat5.bed",sep=""),sep="\t",quote=F,col.names=F,row.names=F)
write.table(tss.summary[tss.summary$exp.cat == 6,][,c(2,3,4,1,5,6)],paste(output,"/tss_1kb_cat6.bed",sep=""),sep="\t",quote=F,col.names=F,row.names=F)
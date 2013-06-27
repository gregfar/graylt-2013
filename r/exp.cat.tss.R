refseq.tss <- read.table("tss_1kb.bed",col.names=c("chr","start","end","name","score","strand"),
                         colClasses=c("character","numeric","numeric","character","numeric","factor"))
summary <- read.table("gene_results_summary.txt",header=T,
                      colClasses=c("character","numeric","numeric","numeric","character","numeric","numeric"))

tss.summary <- merge(refseq.tss,summary,by="name")
tss.summary <- tss.summary[with(tss.summary,order(chr,start)),]

write.table(tss.summary[tss.summary$exp.cat == 1,][,c(2,3,4,1,5,6)],"tss_1kb_cat1.bed",sep="\t",quote=F,col.names=F,row.names=F)
write.table(tss.summary[tss.summary$exp.cat == 2,][,c(2,3,4,1,5,6)],"tss_1kb_cat2.bed",sep="\t",quote=F,col.names=F,row.names=F)
write.table(tss.summary[tss.summary$exp.cat == 3,][,c(2,3,4,1,5,6)],"tss_1kb_cat3.bed",sep="\t",quote=F,col.names=F,row.names=F)
write.table(tss.summary[tss.summary$exp.cat == 4,][,c(2,3,4,1,5,6)],"tss_1kb_cat4.bed",sep="\t",quote=F,col.names=F,row.names=F)
write.table(tss.summary[tss.summary$exp.cat == 5,][,c(2,3,4,1,5,6)],"tss_1kb_cat5.bed",sep="\t",quote=F,col.names=F,row.names=F)
write.table(tss.summary[tss.summary$exp.cat == 6,][,c(2,3,4,1,5,6)],"tss_1kb_cat6.bed",sep="\t",quote=F,col.names=F,row.names=F)
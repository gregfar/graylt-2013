tss.1kb <- read.table("tss_1kb.bed",col.names=c("chr","start","end","name","score","strand"))
tss.100bp <- read.table("tss_100bp.bed",col.names=c("chr","start","end","name","score","strand"))

tss.1kb <- unique(tss.1kb)
tss.100bp <- unique(tss.100bp)

tss.1kb <- tss.1kb[with(tss.1kb,order(chr,start)),]
tss.100bp <- tss.100bp[with(tss.100bp,order(chr,start)),]

write.table(tss.1kb,"tss_1kb.bed",col.names=F,sep="\t",quote=F,row.names=F)
write.table(tss.100bp,"tss_100bp.bed",col.names=F,sep="\t",quote=F,row.names=F)
args <- commandArgs(trailingOnly = T)
results <- args[7]

xpb <- read.table(args[1],col.names=c("chr","r.start","r.end","name","score","strand"))
xpd <- read.table(args[2],col.names=c("chr","r.start","r.end","name","score","strand"))

g4.xpb.regions <- read.table(args[3],col.names=c("chr","r.start","r.end","name","p","strand","peakid"))
g4.xpd.regions <- read.table(args[4],col.names=c("chr","r.start","r.end","name","p","strand","peakid"))

write.table(xpb[xpb$name %in% unique(g4.xpb.regions$peakid),],paste(results,"/xpb_g4_summits.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
write.table(xpd[xpd$name %in% unique(g4.xpd.regions$peakid),],paste(results,"/xpd_g4_summits.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

g4.xpb.peaks <- read.table(args[5],col.names=c("chr","r.start","r.end","name","p","strand","peakid"))
g4.xpd.peaks <- read.table(args[6],col.names=c("chr","r.start","r.end","name","p","strand","peakid"))

write.table(xpb[xpb$name %in% unique(g4.xpb.peaks$peakid),],paste(results,"/xpb_g4_peaks_summits.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
write.table(xpd[xpd$name %in% unique(g4.xpd.peaks$peakid),],paste(results,"/xpd_g4_peaks_summits.bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

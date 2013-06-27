# check to see if the plyr package is installed.
if("plyr" %in% rownames(installed.packages()) == F) {
  install.packages("plyr")
} else {
  library(plyr)
}

peaks.xpb <- read.table("xpb_1e-15_peaks.bed",sep="\t",header=F,col.names=c("chr","start","end","name","score","strand"),
                        colClasses=c("character","numeric","numeric","character","numeric","character"))

g4.xpb <- read.table("g4_xpb_peaks.txt",sep="\t",header=T,
                     colClasses=c("character","numeric","numeric","character","numeric","character","character"))

ap1.xpb <- read.table("ap1_xpb_peaks.txt",sep="\t",header=T,
                     colClasses=c("character","numeric","numeric","character","numeric","character","character"))
                     
ets.xpb <- read.table("ets_xpb_peaks.txt",sep="\t",header=T,
                     colClasses=c("character","numeric","numeric","character","numeric","character","character"))                     
                     
maz.xpb <- read.table("maz_xpb_peaks.txt",sep="\t",header=T,
                     colClasses=c("character","numeric","numeric","character","numeric","character","character"))

xpb.binary <- data.frame(id=peaks.xpb$name)
xpb.binary <- cbind(xpb.binary,g4=0)
xpb.binary$g4[xpb.binary$id %in% g4.xpb$peak_id] <- 1
xpb.binary <- cbind(xpb.binary,ap1=0)
xpb.binary$ap1[xpb.binary$id %in% ap1.xpb$peak_id] <- 1
xpb.binary <- cbind(xpb.binary,ets=0)
xpb.binary$ets[xpb.binary$id %in% ets.xpb$peak_id] <- 1
xpb.binary <- cbind(xpb.binary,maz=0)
xpb.binary$maz[xpb.binary$id %in% maz.xpb$peak_id] <- 1

xpb.counts <- count(xpb.binary,c("g4","ap1","ets","maz"))
xpb.counts <- cbind(xpb.counts,percent=xpb.counts$freq)
xpb.counts$percent <- round((xpb.counts$percent/sum(xpb.counts$freq))*100,2)
xpb.counts <- xpb.counts[with(xpb.counts,order(-percent)),]

write.table(xpb.counts,"xpb_motif_coincidence_peaks.txt",sep="\t",row.names=F,quote=F)

xpb.totals <- data.frame(motif=c("g4","ap1","ets","maz"),
                         count=c(sum(xpb.counts$freq[xpb.counts$g4 == 1]),
                                  sum(xpb.counts$freq[xpb.counts$ap1 == 1]),
                                  sum(xpb.counts$freq[xpb.counts$ets == 1]),
                                  sum(xpb.counts$freq[xpb.counts$maz == 1])),
                         percent=c(sum(xpb.counts$percent[xpb.counts$g4 == 1]),
                                  sum(xpb.counts$percent[xpb.counts$ap1 == 1]),
                                  sum(xpb.counts$percent[xpb.counts$ets == 1]),
                                  sum(xpb.counts$percent[xpb.counts$maz == 1])))

write.table(xpb.totals,"xpb_motif_count_totals.txt",sep="\t",row.names=F,quote=F)

peaks.xpd <- read.table("xpd_1e-15_peaks.bed",sep="\t",header=F,col.names=c("chr","start","end","name","score","strand"),
                        colClasses=c("character","numeric","numeric","character","numeric","character"))

g4.xpd <- read.table("g4_xpd_peaks.txt",sep="\t",header=T,
                     colClasses=c("character","numeric","numeric","character","numeric","character","character"))

ap1.xpd <- read.table("ap1_xpd_peaks.txt",sep="\t",header=T,
                     colClasses=c("character","numeric","numeric","character","numeric","character","character"))
                     
ets.xpd <- read.table("ets_xpd_peaks.txt",sep="\t",header=T,
                     colClasses=c("character","numeric","numeric","character","numeric","character","character"))                     
                     
maz.xpd <- read.table("maz_xpd_peaks.txt",sep="\t",header=T,
                     colClasses=c("character","numeric","numeric","character","numeric","character","character"))

xpd.binary <- data.frame(id=peaks.xpd$name)
xpd.binary <- cbind(xpd.binary,g4=0)
xpd.binary$g4[xpd.binary$id %in% g4.xpd$peak_id] <- 1
xpd.binary <- cbind(xpd.binary,ap1=0)
xpd.binary$ap1[xpd.binary$id %in% ap1.xpd$peak_id] <- 1
xpd.binary <- cbind(xpd.binary,ets=0)
xpd.binary$ets[xpd.binary$id %in% ets.xpd$peak_id] <- 1
xpd.binary <- cbind(xpd.binary,maz=0)
xpd.binary$maz[xpd.binary$id %in% maz.xpd$peak_id] <- 1

xpd.counts <- count(xpd.binary,c("g4","ap1","ets","maz"))
xpd.counts <- cbind(xpd.counts,percent=xpd.counts$freq)
xpd.counts$percent <- round((xpd.counts$percent/sum(xpd.counts$freq))*100,2)
xpd.counts <- xpd.counts[with(xpd.counts,order(-percent)),]

write.table(xpd.counts,"xpd_motif_coincidence_peaks.txt",sep="\t",row.names=F,quote=F)


xpd.totals <- data.frame(motif=c("g4","ap1","ets","maz"),
                         count=c(sum(xpd.counts$freq[xpd.counts$g4 == 1]),
                                  sum(xpd.counts$freq[xpd.counts$ap1 == 1]),
                                  sum(xpd.counts$freq[xpd.counts$ets == 1]),
                                  sum(xpd.counts$freq[xpd.counts$maz == 1])),
                         percent=c(sum(xpd.counts$percent[xpd.counts$g4 == 1]),
                                  sum(xpd.counts$percent[xpd.counts$ap1 == 1]),
                                  sum(xpd.counts$percent[xpd.counts$ets == 1]),
                                  sum(xpd.counts$percent[xpd.counts$maz == 1])))

write.table(xpd.totals,"xpd_motif_count_totals.txt",sep="\t",row.names=F,quote=F)
xpb.s <- read.table("xpb_summits.txt",header=T,
                    col.names=c("chr","p.start","p.end","name","score","strand","s.start","s.end","s.len","s.height","s.center","s.rel.center"),
                    colClasses=c("character","numeric","numeric","character","numeric","factor","numeric","numeric","numeric","numeric","numeric","numeric"))
xpd.s <- read.table("xpd_summits.txt",header=T,
                    col.names=c("chr","p.start","p.end","name","score","strand","s.start","s.end","s.len","s.height","s.center","s.rel.center"),
                    colClasses=c("character","numeric","numeric","character","numeric","factor","numeric","numeric","numeric","numeric","numeric","numeric"))

xpb <- read.table("xpb_regions.rmsk.bed",col.names=c("chr","r.start","r.end","name","score","strand"))
xpb <- merge(xpb,xpb.s)
xpb <- xpb[with(xpb, order(-s.height)),]
write.table(xpb[1:1000,c("chr","r.start","r.end","name","score","strand")],"xpb_regions.top.rmsk.bed",col.names=F,row.names=F,quote=F,sep="\t")

xpd <- read.table("xpd_regions.rmsk.bed",col.names=c("chr","r.start","r.end","name","score","strand"))
xpd <- merge(xpd,xpd.s)
xpd <- xpd[with(xpd, order(-s.height)),]
write.table(xpd[1:1000,c("chr","r.start","r.end","name","score","strand")],"xpd_regions.top.rmsk.bed",col.names=F,row.names=F,quote=F,sep="\t")

both <- read.table("both_regions.rmsk.bed",col.names=c("chr","r.start","r.end","name","score","strand"))
both <- merge(both,xpb.s)
both <- both[with(both, order(-s.height)),]
write.table(both[1:1000,c("chr","r.start","r.end","name","score","strand")],"both_regions.top.rmsk.bed",col.names=F,row.names=F,quote=F,sep="\t")

xpb.o <- read.table("xpb_only_regions.rmsk.bed",col.names=c("chr","r.start","r.end","name","score","strand"))
xpb.o <- merge(xpb.o,xpb.s)
xpb.o <- xpb.o[with(xpb.o,order(-s.height)),]
write.table(xpb.o[1:1000,c("chr","r.start","r.end","name","score","strand")],"xpb_only_regions.top.rmsk.bed",col.names=F,row.names=F,quote=F,sep="\t")

xpd.o <- read.table("xpd_only_regions.rmsk.bed",col.names=c("chr","r.start","r.end","name","score","strand"))
xpd.o <- merge(xpd.o,xpd.s)
xpd.o <- xpd.o[with(xpd.o,order(-s.height)),]
write.table(xpd.o[1:1000,c("chr","r.start","r.end","name","score","strand")],"xpd_only_regions.top.rmsk.bed",col.names=F,row.names=F,quote=F,sep="\t")
	
q()
n

	xpb <- read.table("xpb_summits.bed",col.names=c("chr","r.start","r.end","name","score","strand"))
	xpd <- read.table("xpd_summits.bed",col.names=c("chr","r.start","r.end","name","score","strand"))

	ap1.xpb.regions <- read.table("ap1_xpb.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid","q"))
	ets.xpb.regions <- read.table("ets_xpb.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid","q"))
	maz.xpb.regions <- read.table("maz_xpb.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid","q"))
	g4.xpb.regions <- read.table("g4_xpb.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid"))

	ap1.xpd.regions <- read.table("ap1_xpd.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid","q"))
	ets.xpd.regions <- read.table("ets_xpd.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid","q"))
	maz.xpd.regions <- read.table("maz_xpd.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid","q"))
	g4.xpd.regions <- read.table("g4_xpd.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid"))
	
	write.table(xpb[xpb$name %in% unique(ap1.xpb.regions$peakid),],"xpb_ap1_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")
	write.table(xpb[xpb$name %in% unique(ets.xpb.regions$peakid),],"xpb_ets_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")
	write.table(xpb[xpb$name %in% unique(maz.xpb.regions$peakid),],"xpb_maz_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")
	write.table(xpb[xpb$name %in% unique(g4.xpb.regions$peakid),],"xpb_g4_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")

	write.table(xpd[xpd$name %in% unique(ap1.xpd.regions$peakid),],"xpd_ap1_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")
	write.table(xpd[xpd$name %in% unique(ets.xpd.regions$peakid),],"xpd_ets_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")
	write.table(xpd[xpd$name %in% unique(maz.xpd.regions$peakid),],"xpd_maz_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")
	write.table(xpd[xpd$name %in% unique(g4.xpd.regions$peakid),],"xpd_g4_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")

	ap1.xpb.peaks <- read.table("ap1_xpb_peaks.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid","q"))
	ets.xpb.peaks <- read.table("ets_xpb_peaks.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid","q"))
	maz.xpb.peaks <- read.table("maz_xpb_peaks.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid","q"))
	g4.xpb.peaks <- read.table("g4_xpb_peaks.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid"))

	ap1.xpd.peaks <- read.table("ap1_xpd_peaks.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid","q"))
	ets.xpd.peaks <- read.table("ets_xpd_peaks.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid","q"))
	maz.xpd.peaks <- read.table("maz_xpd_peaks.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid","q"))
	g4.xpd.peaks <- read.table("g4_xpd_peaks.txt",col.names=c("chr","r.start","r.end","name","p","strand","peakid"))
	
	write.table(xpb[xpb$name %in% unique(ap1.xpb.peaks$peakid),],"xpb_ap1_peaks_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")
	write.table(xpb[xpb$name %in% unique(ets.xpb.peaks$peakid),],"xpb_ets_peaks_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")
	write.table(xpb[xpb$name %in% unique(maz.xpb.peaks$peakid),],"xpb_maz_peaks_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")
	write.table(xpb[xpb$name %in% unique(g4.xpb.peaks$peakid),],"xpb_g4_peaks_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")

	write.table(xpd[xpd$name %in% unique(ap1.xpd.peaks$peakid),],"xpd_ap1_peaks_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")
	write.table(xpd[xpd$name %in% unique(ets.xpd.peaks$peakid),],"xpd_ets_peaks_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")
	write.table(xpd[xpd$name %in% unique(maz.xpd.peaks$peakid),],"xpd_maz_peaks_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")
	write.table(xpd[xpd$name %in% unique(g4.xpd.peaks$peakid),],"xpd_g4_peaks_summits.bed",col.names=F,row.names=F,quote=F,sep="\t")

df <- read.table("~/2013_XPB_XPD/XX_tools/chromInfo.txt",sep="\t",stringsAsFactors=F)
df <- df[1:24,]
df$V2 <- df$V2/sum(as.numeric(df$V2))
write.table(df[,1:2],"~/2013_XPB_XPD/XX_tools/chromWeight.txt",sep="\t",col.names=F,row.names=F,quote=F)
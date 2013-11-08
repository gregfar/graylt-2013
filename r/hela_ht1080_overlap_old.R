#Load bioconductor packages
source("http://www.bioconductor.org/biocLite.R")
installed <- rownames(installed.packages())
if(!"Biobase" %in% installed) { biocLite("Biobase") }
if(!"GEOquery" %in% installed) { biocLite("GEOquery") }
if(!"limma" %in% installed) { biocLite("limma") }
library(limma)
library(plyr)
library(ggplot2)
library(Biobase)
library(GEOquery)
source("XX_tools/r/theme_pub.R")

#Retrieve the datasets
hela.1 <- getGEO("GSM472905",destdir=".",GSEMatrix=T)
hela.2 <- getGEO("GSM472916",destdir=".",GSEMatrix=T)
hela.3 <- getGEO("GSM472932",destdir=".",GSEMatrix=T)

#Retrieve the Platform annotations
meta <- Meta(hela.1)
platform <- meta$platform
arr <- getGEO(platform,destdir=".")
arr.data <- Table(arr)

#Convert the datasets to simpler tables and calculate a mean for each probe
hela.1.data <- Table(hela.1)
hela.2.data <- Table(hela.2)
hela.3.data <- Table(hela.3)
hela.vals <- data.frame(hela.1=hela.1.data$VALUE,hela.2=hela.2.data$VALUE,hela.3=hela.3.data$VALUE)
mean.data <- data.frame(ID=hela.1.data$ID_REF,VALUE=rowMeans(hela.vals))

#Filter each row to get just the core probesets for each gene
arr.core <- arr.data[arr.data$level == "core",]
core.data <- mean.data[mean.data$ID %in% arr.core$ID,]

#Average core probesets for each transcript_cluster_id, and transform to log2.
cluster.id <- data.frame(ID=arr.core$ID,transcript_cluster_id=arr.core$transcript_cluster_id)
core.data <- merge(core.data,cluster.id,by="ID")
cluster.data <- ddply(core.data,.(transcript_cluster_id),summarize,mean=mean(VALUE))
log2.cluster.data <- cluster.data
log2.cluster.data$mean <- log(cluster.data$mean,2)

#Remove negative values and get gene name annotations
log2.cluster.data <- log2.cluster.data[log2.cluster.data$mean > 0,]
cluster.genes <- data.frame(transcript_cluster_id=arr.core$transcript_cluster_id,gene_assignment=arr.core$gene_assignment)
cluster.genes$gene_assignment <- sub("^[\\S]+ [/]+ ","",cluster.genes$gene_assignment,perl=T)
cluster.genes$gene_assignment <- sub(" .+","",cluster.genes$gene_assignment,perl=T)
cluster.genes <- unique(cluster.genes)

log2.cluster.data <- merge(log2.cluster.data,cluster.genes,by="transcript_cluster_id")
names(log2.cluster.data)[2] <- "hela.mean"
names(log2.cluster.data)[3] <- "gene.symbol"

ht1080.data <- read.table("orb_gene_table.txt",header=T,stringsAsFactors=F)
ht1080.data <- ht1080.data[ht1080.data$gene.symbol != "---",]

ht1080.hela <- merge(ht1080.data,log2.cluster.data,by="gene.symbol")

all.plot <- ggplot(ht1080.hela,aes(x=mean.exp,y=hela.mean)) + geom_point()

mean.sd <- sd(ht1080.hela$mean.exp - ht1080.hela$hela.mean)

lt.1sd <- ht1080.hela[abs(ht1080.hela$mean.exp - ht1080.hela$hela.mean) < mean.sd,]

lt1sd.plot <- ggplot(lt.1sd,aes(x=mean.exp,y=hela.mean)) + geom_point()

he <- test$hela.mean
ht <- test$mean.exp
m <- lm(he~ht)
resid.sd <- sd(resid(m))

test <- cbind(ht1080.hela,include=1)
test$include[abs(resid(m)) > resid.sd] <- 0

ggplot(test) + 
  geom_point(aes(x=mean.exp,y=hela.mean,colour=as.factor(include)),size=0.4) + 
  scale_colour_manual(values=c("gray","black"),guide=F) + 
  stat_smooth(aes(x=mean.exp,y=hela.mean),method=lm,size=0.5) +
  ylim(0,15) + xlim(0,15) +
  theme_bw(base_size=6)

ggsave("10_results/Figure_S3.png",height=5,width=5)

write.table(test[test$include == 1,],"08_arrays/hela_ht1080_common.txt",row.names=F,sep="\t",quote=F)






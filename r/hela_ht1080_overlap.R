#Load bioconductor packages
source("http://www.bioconductor.org/biocLite.R")
installed <- rownames(installed.packages())
if(!"Biobase" %in% installed) { biocLite("Biobase") }
if(!"GEOquery" %in% installed) { biocLite("GEOquery") }
if(!"limma" %in% installed) { biocLite("limma") }
if(!"oligo" %in% installed) { biocLite("oligo") }

library(limma)
library(plyr)
library(ggplot2)
library(Biobase)
library(GEOquery)
library(oligo)
source("XX_tools/r/theme_pub.R")

#Retrieve the datasets
expect <- c("GSM472905_X_Hs_HELA_E_080117_03_DS7930_D.CEL","GSM472916_X_Hs_HELA_E_090616_01_DS8200_W.CEL","GSM472932_X_Hs_HELA_E_090723_05_DS11528_W.CEL")
cel.files <- list.files(".",pattern="CEL$")
if(sum(cel.files == expect) != 3) {
  getGEOSuppFiles("GSM472905",makeDirectory=F)
  getGEOSuppFiles("GSM472916",makeDirectory=F)
  getGEOSuppFiles("GSM472932",makeDirectory=F)

#Unzip the CEL files, and load them in as an expression set.
gz.files <- list.files(".",pattern="CEL.gz$")
sapply(gz.files, gunzip)
}
cel.files <- list.files(".",pattern="CEL$") 
raw.data <- read.celfiles(cel.files)

# RMA Normalize the data on the core gene cluster level.
norm.data <- rma(raw.data)

#Get the values, and calculate a mean across all 3 expression sets.
hela.vals <- as.data.frame(exprs(norm.data))
mean.data <- data.frame(transcript_cluster_id=row.names(hela.vals),mean=rowMeans(hela.vals))

#Get the annotation database from GEO for the Affy HuEx 1.0 ST v2 arrays
hela.1 <- getGEO("GSM472905",destdir=".",GSEMatrix=T)
arr <- getGEO(Meta(hela.1)$platform,destdir=".")
arr.data <- Table(arr)

#Get transcript_cluster_id to gene symbol mappings
cluster.genes <- data.frame(transcript_cluster_id=arr.data$transcript_cluster_id,gene_assignment=arr.data$gene_assignment)
cluster.genes$gene_assignment <- sub("^[\\S]+ [/]+ ","",cluster.genes$gene_assignment,perl=T)
cluster.genes$gene_assignment <- sub(" .+","",cluster.genes$gene_assignment,perl=T)
cluster.genes <- unique(cluster.genes)
cluster.genes <- cluster.genes[cluster.genes$gene_assignment != "---",]

# merge the gene symbols to the hela data
cluster.data <- merge(mean.data,cluster.genes,by="transcript_cluster_id")
names(cluster.data)[2] <- "hela.mean"
names(cluster.data)[3] <- "gene.symbol"

# read the HT1080 data
ht1080.data <- read.table("08_arrays/orb_gene_table.txt",header=T,stringsAsFactors=F)
ht1080.data <- ht1080.data[ht1080.data$gene.symbol != "---",]

# merge HT1080 and hela data by gene symbols
ht1080.hela <- merge(ht1080.data,cluster.data,by="gene.symbol")

# calculate a linear model for the HT1080 and HeLa data and a standard deviation of residuals
he <- ht1080.hela$hela.mean
ht <- ht1080.hela$mean.exp
m <- lm(he~ht)
resid.sd <- sd(resid(m))

# Add an "include" column, and set it to 0 for any genes with residuals > 1 sd
ht1080.hela <- cbind(ht1080.hela,include=1)
ht1080.hela$include[abs(resid(m)) > resid.sd] <- 0

# Plot relative expression levels for each gene. Included genes will be black, excluded will be gray.
ggplot(ht1080.hela) + 
  geom_point(aes(x=mean.exp,y=hela.mean,colour=as.factor(include)),size=0.4) + 
  scale_colour_manual(values=c("gray","black"),guide=F) + 
  stat_smooth(aes(x=mean.exp,y=hela.mean),method=lm,size=0.5) +
  ylim(0,15) + xlim(0,15) +
  theme_bw(base_size=6)

# Write the graph to a file
ggsave("10_results/Figure_S4.png",height=5,width=5)

# Write a table containing the genes that are in the include == 1 category.
write.table(ht1080.hela[ht1080.hela$include == 1,],"08_arrays/hela_ht1080_common.txt",row.names=F,sep="\t",quote=F)


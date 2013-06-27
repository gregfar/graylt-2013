# check to see if the plyr package is installed.
if("plyr" %in% rownames(installed.packages()) == F) {
  install.packages("plyr")
} else {
  library(plyr)
}

# check to see if the ggplot2 package is installed.
if("ggplot2" %in% rownames(installed.packages()) == F) {
  install.packages("ggplot2")
} else {
  library(ggplot2)
}

# read the TSS 1kb region table, XPB summit, XPD summit, and ORB expression tables
refseq.tss <- read.table("tss_1kb.bed",col.names=c("chr","start","end","name","score","strand"),
                         colClasses=c("character","numeric","numeric","character","numeric","factor"))
xpb <- read.table("xpb_summits.bed",col.names=c("chr","start","end","name","score","strand"),
                  colClasses=c("character","numeric","numeric","character","numeric","factor"),skip=1)
xpd <- read.table("xpd_summits.bed",col.names=c("chr","start","end","name","score","strand"),
                  colClasses=c("character","numeric","numeric","character","numeric","factor"),skip=1)
orb <- read.table("orb_gene_table.txt",header=T,sep="\t",
                  colClasses=c("character","character","character","numeric","numeric","numeric","numeric","numeric","character"))

# this function will flag overlapping regions with a 1
bed.flag <- function(bed1,bed2) {
  
  result <- cbind(bed1,overlap=0)
  
  for (i in 1:nrow(bed1)) {
    
    sub2 <- bed2[bed2$chr == bed1$chr[i],]
    
    compare <- sub2[sub2$start <= bed1$end[i] & sub2$end >= bed1$start[i],]
    
    if(nrow(compare) > 0) {
      result$overlap[i] <- 1
    }
    
  }
  
  return(result)
  
}

# this function will join overlapping bed regions
bed.join <- function(bed1,bed2) {
  
  result <- cbind(bed1[1,],bed2[1,])
  
  for (i in 1:nrow(bed1)) {
    
    sub2 <- bed2[bed2$chr == bed1$chr[i],]
    
    compare <- sub2[sub2$start <= bed1$end[i] & sub2$end >= bed1$start[i],]
    
    if(nrow(compare) > 0) {
      for (j in 1:nrow(compare)) {
        result.row <- cbind(bed1[i,],compare[j,])
        result <- rbind(result,result.row)
      }
    }
    
  }
  
  result <- result[-1,]
  
  return(result)
  
}


## Figure 3A calculations
tss.xpb <- bed.join(refseq.tss,xpb)
colnames(tss.xpb) <- c("tss.chr","tss.start","tss.end","tss.name","tss.score","tss.strand",
                       "sum.chr","sum.start","sum.end","sum.name","sum.score","sum.strand")
tss.xpd <- bed.join(refseq.tss,xpd)
colnames(tss.xpd) <- c("tss.chr","tss.start","tss.end","tss.name","tss.score","tss.strand",
                       "sum.chr","sum.start","sum.end","sum.name","sum.score","sum.strand")

tss.xpb <- cbind(tss.xpb,distance=tss.xpb$sum.start-(tss.xpb$tss.start+1000))
tss.xpd <- cbind(tss.xpd,distance=tss.xpd$sum.start-(tss.xpd$tss.start+1000))

u.tss.xpb <- ddply(tss.xpb,c("tss.name","tss.strand"),summarise,distance=min(abs(distance)))
u.tss.xpb$distance[u.tss.xpb$tss.strand == "-"] <- u.tss.xpb$distance[u.tss.xpb$tss.strand == "-"] * -1

u.tss.xpd <- ddply(tss.xpd,c("tss.name","tss.strand"),summarise,distance=min(abs(distance)))
u.tss.xpd$distance[u.tss.xpd$tss.strand == "-"] <- u.tss.xpd$distance[u.tss.xpd$tss.strand == "-"] * -1

theme_pub <- function (base_size = 6, base_family = "") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text = element_text(size = rel(1)), 
          axis.ticks = element_line(colour = "black"), 
          legend.key = element_rect(colour = "grey80"), 
          legend.position = "none",
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA, colour = "grey50"), 
          panel.grid.major = element_line(colour = "grey80", size = 1), 
          panel.grid.minor = element_line(colour = "grey95", size =1), 
          strip.background = element_rect(fill = "grey80", colour = "grey50"), 
          strip.background = element_rect(fill = "grey80", colour = "grey50"))
}

ggplot(data=u.tss.xpb, aes(x=distance)) +
  geom_histogram(aes(y = (..count..)/21555),binwidth=50,fill="mediumorchid4") +
  scale_x_continuous(limits=c(-1000,1000)) + ylim(0,.018) +
  theme_pub()
ggsave("Figure_3aL.png",width=3,height=3)

ggplot(data=u.tss.xpd, aes(x=distance)) +
  geom_histogram(aes(y = (..count..)/14570),binwidth=50,fill="darkgreen") +
  scale_x_continuous(limits=c(-1000,1000)) + ylim(0,.018) +
  theme_pub()
ggsave("Figure3aR.png",width=3,height=3)

## Figure 3B calculations
# get the mean expression scores from the ORB table, sort by score, and assign expression categories
orb.scores <- data.frame(name=orb$gene.symbol,exp=orb$mean.exp,stringsAsFactors=F)
orb.scores <- orb.scores[orb.scores$name != "---",]
orb.scores <- orb.scores[with(orb.scores,order(-exp)),]
exp.cat <- c(rep(1,2924),rep(2,2923),rep(3,2924),rep(4,2923),rep(5,2924))
orb.scores <- cbind(orb.scores,exp.cat=exp.cat)

# get unique rows of the refseq tss table
refseq.tss <- unique(refseq.tss)

# merge the ORB expression data with the refseq table by gene name
df <- merge(refseq.tss,orb.scores,all.x=T,by="name")
df$exp[is.na(df$exp)] <- 0
df$exp.cat[df$exp == 0] <- 6

# flag TSS that overlap XPB and XPD summits
df <- bed.flag(df,xpb)
colnames(df)[9] <- "xpb"
df <- bed.flag(df,xpd)
colnames(df)[10] <- "xpd"

# Use ddply to get a single line for each gene name. Genes with an xpb or xpd peak at any of
# their annotated TSS will be considered bound.
binding <- ddply(df,"name",summarise,xpb=max(xpb),xpd=max(xpd),exp=max(exp),exp.cat=min(exp.cat))

# binding combinations are converted to categories for summarizing
binding <- cbind(binding,bind.cat=0)
binding$bind.cat[binding$xpb == 1 & binding$xpd == 1] <- "both"
binding$bind.cat[binding$xpb == 1 & binding$xpd == 0] <- "xpb.only"
binding$bind.cat[binding$xpb == 0 & binding$xpd == 1] <- "xpd.only"
binding$bind.cat[binding$xpb == 0 & binding$xpd == 0] <- "unbound"

# count the number of genes with neither, both, xpb only, or xpd only.
counts <- count(binding,c("exp.cat","bind.cat"))
counts <- cbind(counts,percent=counts$freq)
counts <- ddply(counts,"exp.cat",transform,percent=percent/sum(percent))

# Plot the results
ggplot(data=counts[counts$bind.cat != "unbound",],aes(x=factor(exp.cat),y=percent,fill=factor(bind.cat))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("dodgerblue3","dodgerblue3","dodgerblue3")) +
  theme_pub()

# output the plot
ggsave("Figure_3b_onecolor.png",height=3,width=3)

# Plot the results
ggplot(data=counts[counts$bind.cat != "unbound",],aes(x=factor(exp.cat),y=percent,fill=factor(exp.cat))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#D55E00","#E69F00","gold2","#009E73","#56B4E9","grey40")) +
  theme_pub()

# output the plot
ggsave("Figure_3b_rainbow.png",height=3,width=3)
 
  
# output a summary table
write.table(binding,"tss_binding_summary.txt",sep="\t",row.names=F,quote=F)
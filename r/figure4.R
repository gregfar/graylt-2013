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

# Add the ggplot2 theme for publication
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

# read the exon table, and the tss binding results table
orb.exon <- read.table("orb_exon_table.txt",header=T,
                       colClasses=c("character","character","factor","numeric","numeric","character",
                                    "numeric","numeric","numeric","numeric","numeric"))
binding <- read.table("tss_binding_summary.txt",header=T,
                      colClasses=c("character","numeric","numeric","numeric","character"))

# remove unnamed genes from the exon table
orb.exon <- orb.exon[orb.exon$gene.symbol != "---",]

# separate the exon table into + and - strand
orb.exon.plus <- orb.exon[orb.exon$strand == "+",]
orb.exon.minus <- orb.exon[orb.exon$strand == "-",]

# calculate the ratio of exon 1 to downstream exons for each set of exon data
orb.exon.plus.ratio <- ddply(orb.exon.plus,"gene.symbol",summarise,ratio=mean.exp[start == min(start)]/mean(mean.exp[start != min(start)]))
orb.exon.minus.ratio <- ddply(orb.exon.minus,"gene.symbol",summarise,ratio=mean.exp[end == max(end)]/mean(mean.exp[end != max(end)]))

# bind the results together, and correct for the existence of multiple entries per gene name
# (some genes are annotated on both strands, for some reason)
orb.exon.ratio <- rbind(orb.exon.plus.ratio,orb.exon.minus.ratio)
orb.exon.ratio <- ddply(orb.exon.ratio,"gene.symbol",summarise,ratio=max(ratio))

# merge the new exon ratio results to the binding table
results <- merge(binding,orb.exon.ratio,by.x="name",by.y="gene.symbol",all.x=T)
results$ratio[is.na(results$ratio)] <- 0

# calculate the cutoff for non-processive genes (mean + sd of ratio)
mean(results$ratio[results$ratio > 0]) + sd(results$ratio[results$ratio > 0])

# add a flag column for non-processive genes, and flag 1 for each gene above the cutoff
results <- cbind(results,np=0)
results$np[results$ratio > mean(results$ratio[results$ratio > 0]) + sd(results$ratio[results$ratio > 0])] <- 1

# count the number of genes in each binding category that are non-processive, and
# calculate the percent of each category that are non-processive.
np.counts <- count(results[results$exp.cat < 6,],c("bind.cat","np"))
np.counts <- cbind(np.counts,percent=np.counts$freq)
np.counts <- ddply(np.counts,"bind.cat",transform,percent=percent/sum(freq))
np.counts$bind.cat <- ordered(np.counts$bind.cat, levels=c("unbound","xpb.only","xpd.only","both"))

# graph the percent of genes that are non-processive in each category.
ggplot(data=np.counts[np.counts$np == 1,],aes(x=factor(bind.cat),y=percent,fill=factor(bind.cat))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("gray50","mediumorchid4","darkgreen","dodgerblue3")) +
  theme_pub()

ggsave("Figure4b.png",width=3,height=3)

write.table(results,"gene_results_summary.txt",sep="\t",row.names=F,quote=F)

#both vs unbound
ggplot() + 
  geom_vline(aes(xintercept=1.43)) +
  geom_freqpoly(aes(x=results$ratio[results$exp.cat < 6 & results$ratio > 0 & results$bind.cat == "unbound"],y=(..count..)/sum(..count..)),binwidth=.05,color="gray40",size=1) +
  geom_freqpoly(aes(x=results$ratio[results$exp.cat < 6 & results$ratio > 0 & results$bind.cat == "both"],y=(..count..)/sum(..count..)),binwidth=.05,color="dodgerblue",size=1) +
  xlim(0,2) +
  ylim(0,.12) +
  theme_pub()
ggsave("Figure4cBoth.png",width=3,height=3)

#xpb only vs unbound
ggplot() + 
  geom_vline(aes(xintercept=1.43)) +
  geom_freqpoly(aes(x=results$ratio[results$exp.cat < 6 & results$ratio > 0  & results$bind.cat == "unbound"],y=(..count..)/sum(..count..)),binwidth=.05,color="gray40",size=1) +
  geom_freqpoly(aes(x=results$ratio[results$exp.cat < 6 & results$ratio > 0  & results$bind.cat == "xpb.only"],y=(..count..)/sum(..count..)),binwidth=.05,color="mediumorchid4",size=1) +
  xlim(0,2) +
  ylim(0,.12) +
  theme_pub()
ggsave("Figure4cXPB.png",width=3,height=3)

#xpd only vs unbound
ggplot() + 
  geom_vline(aes(xintercept=1.43)) +
  geom_freqpoly(aes(x=results$ratio[results$exp.cat < 6 & results$ratio > 0  & results$bind.cat == "unbound"],y=(..count..)/sum(..count..)),binwidth=.05,color="gray40",size=1) +
  geom_freqpoly(aes(x=results$ratio[results$exp.cat < 6 & results$ratio > 0  & results$bind.cat == "xpd.only"],y=(..count..)/sum(..count..)),binwidth=.05,color="darkgreen",size=1) +
  xlim(0,2) +
  ylim(0,.12) +
  theme_pub()
ggsave("Figure4cXPD.png",width=3,height=3)

#Calculate p-values and output figure 4a (data for all expressed genes)
fig4a <- ddply(results[results$exp.cat < 6 & results$ratio > 0,],"bind.cat",summarise,mean=mean(ratio),sd=sd(ratio))
temp <- c("all",mean(results$ratio[results$exp.cat < 6 & results$ratio > 0]),sd(results$ratio[results$exp.cat < 6 & results$ratio > 0]))
fig4a <- rbind(temp,fig4a)
fig4a <- fig4a[c(1,3,2,4,5),]
temp <- data.frame(p=c("---",
          "---",
          t.test(results$ratio[results$exp.cat < 6 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat < 6 & results$ratio > 0 & results$bind.cat == "both"])[[3]],
          t.test(results$ratio[results$exp.cat < 6 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat < 6 & results$ratio > 0 & results$bind.cat == "xpb.only"])[[3]],
          t.test(results$ratio[results$exp.cat < 6 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat < 6 & results$ratio > 0 & results$bind.cat == "xpd.only"])[[3]]
          ))
fig4a <- cbind(fig4a,temp)

write.table(fig4a,"Figure4a.txt",sep="\t",quote=F,row.names=F)

#Calculate p-values and output figure 4a (expression category 1 only)
fig4_cat1 <- ddply(results[results$exp.cat == 1 & results$ratio > 0,],"bind.cat",summarise,mean=mean(ratio),sd=sd(ratio))
temp <- c("all",mean(results$ratio[results$exp.cat == 1 & results$ratio > 0]),sd(results$ratio[results$exp.cat == 1 & results$ratio > 0]))
fig4_cat1 <- rbind(temp,fig4_cat1)
fig4_cat1 <- fig4_cat1[c(1,3,2,4,5),]
temp <- data.frame(p=c("---",
          "---",
          t.test(results$ratio[results$exp.cat == 1 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 1 & results$ratio > 0 & results$bind.cat == "both"])[[3]],
          t.test(results$ratio[results$exp.cat == 1 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 1 & results$ratio > 0 & results$bind.cat == "xpb.only"])[[3]],
          t.test(results$ratio[results$exp.cat == 1 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 1 & results$ratio > 0 & results$bind.cat == "xpd.only"])[[3]]
          ))
fig4_cat1 <- cbind(fig4_cat1,temp)

write.table(fig4_cat1,"Figure4_exp_cat_1.txt",sep="\t",quote=F,row.names=F)

#Calculate p-values and output figure 4a (expression category 2 only)
fig4_cat2 <- ddply(results[results$exp.cat == 2 & results$ratio > 0,],"bind.cat",summarise,mean=mean(ratio),sd=sd(ratio))
temp <- c("all",mean(results$ratio[results$exp.cat == 2 & results$ratio > 0]),sd(results$ratio[results$exp.cat == 2 & results$ratio > 0]))
fig4_cat2 <- rbind(temp,fig4_cat2)
fig4_cat2 <- fig4_cat2[c(1,3,2,4,5),]
temp <- data.frame(p=c("---",
          "---",
          t.test(results$ratio[results$exp.cat == 2 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 2 & results$ratio > 0 & results$bind.cat == "both"])[[3]],
          t.test(results$ratio[results$exp.cat == 2 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 2 & results$ratio > 0 & results$bind.cat == "xpb.only"])[[3]],
          t.test(results$ratio[results$exp.cat == 2 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 2 & results$ratio > 0 & results$bind.cat == "xpd.only"])[[3]]
          ))
fig4_cat2 <- cbind(fig4_cat2,temp)

write.table(fig4_cat2,"Figure4_exp_cat_2.txt",sep="\t",quote=F,row.names=F)

#Calculate p-values and output figure 4a (expression category 3 only)
fig4_cat3 <- ddply(results[results$exp.cat == 3 & results$ratio > 0,],"bind.cat",summarise,mean=mean(ratio),sd=sd(ratio))
temp <- c("all",mean(results$ratio[results$exp.cat == 3 & results$ratio > 0]),sd(results$ratio[results$exp.cat == 3 & results$ratio > 0]))
fig4_cat3 <- rbind(temp,fig4_cat3)
fig4_cat3 <- fig4_cat3[c(1,3,2,4,5),]
temp <- data.frame(p=c("---",
          "---",
          t.test(results$ratio[results$exp.cat == 3 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 3 & results$ratio > 0 & results$bind.cat == "both"])[[3]],
          t.test(results$ratio[results$exp.cat == 3 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 3 & results$ratio > 0 & results$bind.cat == "xpb.only"])[[3]],
          t.test(results$ratio[results$exp.cat == 3 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 3 & results$ratio > 0 & results$bind.cat == "xpd.only"])[[3]]
          ))
fig4_cat3 <- cbind(fig4_cat3,temp)

write.table(fig4_cat3,"Figure4_exp_cat_3.txt",sep="\t",quote=F,row.names=F)

#Calculate p-values and output figure 4a (expression category 4 only)
fig4_cat4 <- ddply(results[results$exp.cat == 4 & results$ratio > 0,],"bind.cat",summarise,mean=mean(ratio),sd=sd(ratio))
temp <- c("all",mean(results$ratio[results$exp.cat == 4 & results$ratio > 0]),sd(results$ratio[results$exp.cat == 4 & results$ratio > 0]))
fig4_cat4 <- rbind(temp,fig4_cat4)
fig4_cat4 <- fig4_cat4[c(1,3,2,4,5),]
temp <- data.frame(p=c("---",
          "---",
          t.test(results$ratio[results$exp.cat == 4 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 4 & results$ratio > 0 & results$bind.cat == "both"])[[3]],
          t.test(results$ratio[results$exp.cat == 4 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 4 & results$ratio > 0 & results$bind.cat == "xpb.only"])[[3]],
          t.test(results$ratio[results$exp.cat == 4 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 4 & results$ratio > 0 & results$bind.cat == "xpd.only"])[[3]]
          ))
fig4_cat4 <- cbind(fig4_cat4,temp)

write.table(fig4_cat4,"Figure4_exp_cat_4.txt",sep="\t",quote=F,row.names=F)

#Calculate p-values and output figure 4a (expression category 5 only)
fig4_cat5 <- ddply(results[results$exp.cat == 5 & results$ratio > 0,],"bind.cat",summarise,mean=mean(ratio),sd=sd(ratio))
temp <- c("all",mean(results$ratio[results$exp.cat == 5 & results$ratio > 0]),sd(results$ratio[results$exp.cat == 5 & results$ratio > 0]))
fig4_cat5 <- rbind(temp,fig4_cat5)
fig4_cat5 <- fig4_cat5[c(1,3,2,4,5),]
temp <- data.frame(p=c("---",
          "---",
          t.test(results$ratio[results$exp.cat == 5 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 5 & results$ratio > 0 & results$bind.cat == "both"])[[3]],
          t.test(results$ratio[results$exp.cat == 5 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 5 & results$ratio > 0 & results$bind.cat == "xpb.only"])[[3]],
          t.test(results$ratio[results$exp.cat == 5 & results$ratio > 0 & results$bind.cat == "unbound"],results$ratio[results$exp.cat == 5 & results$ratio > 0 & results$bind.cat == "xpd.only"])[[3]]
          ))
fig4_cat5 <- cbind(fig4_cat5,temp)

write.table(fig4_cat5,"Figure4_exp_cat_5.txt",sep="\t",quote=F,row.names=F)

temp <- data.frame(set1="unbound",set2="both",chisq.p=chisq.test(matrix(np.counts$freq[np.counts$bind.cat == "unbound" | np.counts$bind.cat == "both"],2,2))[[3]])
fig4b.chi <- temp
temp <- data.frame(set1="unbound",set2="xpb.only",chisq.p=chisq.test(matrix(np.counts$freq[np.counts$bind.cat == "unbound" | np.counts$bind.cat == "xpb.only"],2,2))[[3]])
fig4b.chi <- rbind(fig4b.chi,temp)
temp <- data.frame(set1="unbound",set2="xpd.only",chisq.p=chisq.test(matrix(np.counts$freq[np.counts$bind.cat == "unbound" | np.counts$bind.cat == "xpd.only"],2,2))[[3]])
fig4b.chi <- rbind(fig4b.chi,temp)

write.table(fig4b.chi,"Figure4b_chisq.txt",sep="\t",quote=F,row.names=F)

# check to see if the ggplot2 package is installed.
if("ggplot2" %in% rownames(installed.packages()) == F) {
  install.packages("ggplot2")
} else {
  library(ggplot2)
}

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

xpb.1 <- read.table("xpb_pileups/xpb_tss_1kb_cat1.pile",header=T,colClasses=c("numeric","numeric"))
xpb.2 <- read.table("xpb_pileups/xpb_tss_1kb_cat2.pile",header=T,colClasses=c("numeric","numeric"))
xpb.3 <- read.table("xpb_pileups/xpb_tss_1kb_cat3.pile",header=T,colClasses=c("numeric","numeric"))
xpb.4 <- read.table("xpb_pileups/xpb_tss_1kb_cat4.pile",header=T,colClasses=c("numeric","numeric"))
xpb.5 <- read.table("xpb_pileups/xpb_tss_1kb_cat5.pile",header=T,colClasses=c("numeric","numeric"))

xpb.all <- merge(xpb.1,xpb.2,by="pos")
colnames(xpb.all) <- c("pos","val.1","val.2")
xpb.all <- merge(xpb.all,xpb.3,by="pos")
colnames(xpb.all)[4] <- "val.3"
xpb.all <- merge(xpb.all,xpb.4,by="pos")
colnames(xpb.all)[5] <- "val.4"
xpb.all <- merge(xpb.all,xpb.5,by="pos")
colnames(xpb.all)[6] <- "val.5"

xpd.1 <- read.table("xpd_pileups/xpd_tss_1kb_cat1.pile",header=T,colClasses=c("numeric","numeric"))
xpd.2 <- read.table("xpd_pileups/xpd_tss_1kb_cat2.pile",header=T,colClasses=c("numeric","numeric"))
xpd.3 <- read.table("xpd_pileups/xpd_tss_1kb_cat3.pile",header=T,colClasses=c("numeric","numeric"))
xpd.4 <- read.table("xpd_pileups/xpd_tss_1kb_cat4.pile",header=T,colClasses=c("numeric","numeric"))
xpd.5 <- read.table("xpd_pileups/xpd_tss_1kb_cat5.pile",header=T,colClasses=c("numeric","numeric"))

xpd.all <- merge(xpd.1,xpd.2,by="pos")
colnames(xpd.all) <- c("pos","val.1","val.2")
xpd.all <- merge(xpd.all,xpd.3,by="pos")
colnames(xpd.all)[4] <- "val.3"
xpd.all <- merge(xpd.all,xpd.4,by="pos")
colnames(xpd.all)[5] <- "val.4"
xpd.all <- merge(xpd.all,xpd.5,by="pos")
colnames(xpd.all)[6] <- "val.5"

i.xpb.1 <- read.table("xpb_input_pileups/xpb_input_tss_1kb_cat1.pile",header=T,colClasses=c("numeric","numeric"))
i.xpb.2 <- read.table("xpb_input_pileups/xpb_input_tss_1kb_cat2.pile",header=T,colClasses=c("numeric","numeric"))
i.xpb.3 <- read.table("xpb_input_pileups/xpb_input_tss_1kb_cat3.pile",header=T,colClasses=c("numeric","numeric"))
i.xpb.4 <- read.table("xpb_input_pileups/xpb_input_tss_1kb_cat4.pile",header=T,colClasses=c("numeric","numeric"))
i.xpb.5 <- read.table("xpb_input_pileups/xpb_input_tss_1kb_cat5.pile",header=T,colClasses=c("numeric","numeric"))

i.xpb.all <- merge(i.xpb.1,i.xpb.2,by="pos")
colnames(i.xpb.all) <- c("pos","val.1","val.2")
i.xpb.all <- merge(i.xpb.all,i.xpb.3,by="pos")
colnames(i.xpb.all)[4] <- "val.3"
i.xpb.all <- merge(i.xpb.all,i.xpb.4,by="pos")
colnames(i.xpb.all)[5] <- "val.4"
i.xpb.all <- merge(i.xpb.all,i.xpb.5,by="pos")
colnames(i.xpb.all)[6] <- "val.5"

i.xpd.1 <- read.table("xpd_input_pileups/xpd_input_tss_1kb_cat1.pile",header=T,colClasses=c("numeric","numeric"))
i.xpd.2 <- read.table("xpd_input_pileups/xpd_input_tss_1kb_cat2.pile",header=T,colClasses=c("numeric","numeric"))
i.xpd.3 <- read.table("xpd_input_pileups/xpd_input_tss_1kb_cat3.pile",header=T,colClasses=c("numeric","numeric"))
i.xpd.4 <- read.table("xpd_input_pileups/xpd_input_tss_1kb_cat4.pile",header=T,colClasses=c("numeric","numeric"))
i.xpd.5 <- read.table("xpd_input_pileups/xpd_input_tss_1kb_cat5.pile",header=T,colClasses=c("numeric","numeric"))

i.xpd.all <- merge(i.xpd.1,i.xpd.2,by="pos")
colnames(i.xpd.all) <- c("pos","val.1","val.2")
i.xpd.all <- merge(i.xpd.all,i.xpd.3,by="pos")
colnames(i.xpd.all)[4] <- "val.3"
i.xpd.all <- merge(i.xpd.all,i.xpd.4,by="pos")
colnames(i.xpd.all)[5] <- "val.4"
i.xpd.all <- merge(i.xpd.all,i.xpd.5,by="pos")
colnames(i.xpd.all)[6] <- "val.5"

n.xpb <- cbind(pos=xpb.all[,1],xpb.all[,2:6]/i.xpb.all[,2:6])
n.xpd <- cbind(pos=xpd.all[,1],xpd.all[,2:6]/i.xpd.all[,2:6])

ggplot(data=n.xpb) + geom_vline(aes(xintercept=1000),size=.5) +
  geom_line(aes(pos,val.5,color="cat5"),size=1) +
  geom_line(aes(pos,val.4,color="cat4"),size=1) +
  geom_line(aes(pos,val.3,color="cat3"),size=1) +
  geom_line(aes(pos,val.2,color="cat2"),size=1) +
  geom_line(aes(pos,val.1,color="cat1"),size=1) + 
  scale_color_manual(name="Exp Lvl",breaks=c("cat5","cat4","cat3","cat2","cat1"),values=c("#D55E00","#E69F00","gold2","#009E73","#56B4E9")) +
  scale_y_continuous(name="XPB Enrichment (XPB ChIP/Input)",limits=c(0,12),breaks=seq(0,12,2)) +
  scale_x_continuous(name="Position Relative to TSS",limits=c(0,2000)) +
  guides(fill=F) +
  theme_pub()
ggsave("Figure3dXPB.png",width=3,height=3)

ggplot(data=n.xpd) + geom_vline(aes(xintercept=1000),size=.5) +
  geom_line(aes(pos,val.5,color="cat5"),size=1) +
  geom_line(aes(pos,val.4,color="cat4"),size=1) +
  geom_line(aes(pos,val.3,color="cat3"),size=1) +
  geom_line(aes(pos,val.2,color="cat2"),size=1) +
  geom_line(aes(pos,val.1,color="cat1"),size=1) + 
  scale_color_manual(name="Exp Lvl",breaks=c("cat5","cat4","cat3","cat2","cat1"),values=c("#D55E00","#E69F00","gold2","#009E73","#56B4E9")) +
  scale_y_continuous(name="XPD Enrichment (XPD ChIP/Input)",limits=c(0,12),breaks=seq(0,12,2)) +
  scale_x_continuous(name="Position Relative to TSS",limits=c(0,2000)) +
  theme_pub()
ggsave("Figure3dXPD.png",width=3,height=3)

args <- commandArgs(trailingOnly = T)

setwd("~/2013_XPB_XPD/")
args <- c("10_results/tss_binding_summary.txt",
          "09_treatments/halder_360a",
          "09_treatments/halder_phendc3",
          "09_treatments/halder_8979a",
          "09_treatments/gray_phendc3",
          "10_results")
results <- args[6]

summary <- read.table(args[1],sep="\t",header=T,stringsAsFactors=F)

h3.dn <- read.table(paste(args[2],"_dn.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)
h3.up <- read.table(paste(args[2],"_up.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)

p.dn <- read.table(paste(args[3],"_dn.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)
p.up <- read.table(paste(args[3],"_up.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)

c.dn <- read.table(paste(args[4],"_dn.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)
c.up <- read.table(paste(args[4],"_up.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)

g.dn <- read.table(paste(args[5],"_dn.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)
g.up <- read.table(paste(args[5],"_up.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)

df <- cbind(summary,h360a=0,h8979a=0,hphendc3=0,gphendc3=0)
df <- unique(df)

temp <- df[0,]
names <- unique(df$name)
for(i in 1:length(names)) {
  
  sub <- df[df$name == names[i],]
  
  if(nrow(sub) > 1) {
    
    new.xpb <- max(sub$xpb)
    new.xpd <- max(sub$xpd)
    
    sub <- sub[1,]
    sub$xpb <- new.xpb
    sub$xpd <- new.xpd
    
    sub$bind.cat[sub$xpb == 1 & sub$xpd == 1] <- "both"
    sub$bind.cat[sub$xpb == 1 & sub$xpd == 0] <- "xpb.only"
    sub$bind.cat[sub$xpb == 0 & sub$xpd == 1] <- "xpd.only"
    sub$bind.cat[sub$xpb == 0 & sub$xpd == 0] <- "unbound"
    
    temp <- rbind(temp,sub)
    
  } else {
    temp <- rbind(temp,sub)
  }
  
}
df <- temp

df$h360a[df$name %in% h3.dn$name] <- -1
df$h360a[df$name %in% h3.up$name] <- 1

df$h8979a[df$name %in% c.dn$name] <- -1
df$h8979a[df$name %in% c.up$name] <- 1

df$hphendc3[df$name %in% p.dn$name] <- -1
df$hphendc3[df$name %in% p.up$name] <- 1

df$gphendc3[df$name %in% g.dn$name] <- -1
df$gphendc3[df$name %in% g.up$name] <- 1

# All Gray PhenDC3 Genes
gphendc3.h360a <- phyper(sum(df$gphendc3 != 0 & df$h360a != 0),
                    sum(df$gphendc3 != 0),
                    sum(df$gphendc3 == 0),
                    sum(df$h360a != 0),
                    lower.tail=F)
gphendc3.hphendc3 <- phyper(sum(df$gphendc3 != 0 & df$hphendc3 != 0),
                       sum(df$gphendc3 != 0),
                       sum(df$gphendc3 == 0),
                       sum(df$hphendc3 != 0),
                       lower.tail=F)
gphendc3.h8979a <- phyper(sum(df$gphendc3 != 0 & df$h8979a != 0),
                     sum(df$gphendc3 != 0),
                     sum(df$gphendc3 == 0),
                     sum(df$h8979a != 0),
                     lower.tail=F)
gphendc3.gphendc3 <- phyper(sum(df$gphendc3 != 0 & df$gphendc3 != 0),
                       sum(df$gphendc3 != 0),
                       sum(df$gphendc3 == 0),
                       sum(df$gphendc3 != 0),
                       lower.tail=F)

gphendc3.h360a.up <- phyper(sum(df$gphendc3 != 0 & df$h360a == 1),
                       sum(df$gphendc3 != 0),
                       sum(df$gphendc3 == 0),
                       sum(df$h360a == 1),
                       lower.tail=F)
gphendc3.hphendc3.up <- phyper(sum(df$gphendc3 != 0 & df$hphendc3 == 1),
                          sum(df$gphendc3 != 0),
                          sum(df$gphendc3 == 0),
                          sum(df$hphendc3 == 1),
                          lower.tail=F)
gphendc3.h8979a.up <- phyper(sum(df$gphendc3 != 0 & df$h8979a == 1),
                        sum(df$gphendc3 != 0),
                        sum(df$gphendc3 == 0),
                        sum(df$h8979a == 1),
                        lower.tail=F)
gphendc3.gphendc3.up <- phyper(sum(df$gphendc3 != 0 & df$gphendc3 == 1),
                          sum(df$gphendc3 != 0),
                          sum(df$gphendc3 == 0),
                          sum(df$gphendc3 == 1),
                          lower.tail=F)

gphendc3.h360a.dn <- phyper(sum(df$gphendc3 != 0 & df$h360a == -1),
                       sum(df$gphendc3 != 0),
                       sum(df$gphendc3 == 0),
                       sum(df$h360a == -1),
                       lower.tail=F)
gphendc3.hphendc3.dn <- phyper(sum(df$gphendc3 != 0 & df$hphendc3 == -1),
                          sum(df$gphendc3 != 0),
                          sum(df$gphendc3 == 0),
                          sum(df$hphendc3 == -1),
                          lower.tail=F)
gphendc3.h8979a.dn <- phyper(sum(df$gphendc3 != 0 & df$h8979a == -1),
                        sum(df$gphendc3 != 0),
                        sum(df$gphendc3 == 0),
                        sum(df$h8979a == -1),
                        lower.tail=F)

#Up-regulated in Gray PhenDC3
up.gphendc3.h360a <- phyper(sum(df$gphendc3 == 1 & df$h360a != 0),
                         sum(df$gphendc3 == 1),
                         sum(df$gphendc3 == 0),
                         sum(df$h360a != 0),
                         lower.tail=F)
up.gphendc3.hphendc3 <- phyper(sum(df$gphendc3 == 1 & df$hphendc3 != 0),
                            sum(df$gphendc3 == 1),
                            sum(df$gphendc3 == 0),
                            sum(df$hphendc3 != 0),
                            lower.tail=F)
up.gphendc3.h8979a <- phyper(sum(df$gphendc3 == 1 & df$h8979a != 0),
                          sum(df$gphendc3 == 1),
                          sum(df$gphendc3 == 0),
                          sum(df$h8979a != 0),
                          lower.tail=F)
up.gphendc3.gphendc3 <- phyper(sum(df$gphendc3 == 1 & df$gphendc3 == 1),
                            sum(df$gphendc3 == 1),
                            sum(df$gphendc3 == 0),
                            sum(df$gphendc3 == 1),
                            lower.tail=F)

up.gphendc3.h360a.up <- phyper(sum(df$gphendc3 == 1 & df$h360a == 1),
                            sum(df$gphendc3 == 1),
                            sum(df$gphendc3 == 0),
                            sum(df$h360a == 1),
                            lower.tail=F)
up.gphendc3.hphendc3.up <- phyper(sum(df$gphendc3 == 1 & df$hphendc3 == 1),
                               sum(df$gphendc3 == 1),
                               sum(df$gphendc3 == 0),
                               sum(df$hphendc3 == 1),
                               lower.tail=F)
up.gphendc3.h8979a.up <- phyper(sum(df$gphendc3 == 1 & df$h8979a == 1),
                             sum(df$gphendc3 == 1),
                             sum(df$gphendc3 == 0),
                             sum(df$h8979a == 1),
                             lower.tail=F)
up.gphendc3.up.gphendc3.up <- phyper(sum(df$gphendc3 == 1 & df$gphendc3 == 1),
                               sum(df$gphendc3 == 1),
                               sum(df$gphendc3 == 0),
                               sum(df$gphendc3 == 1),
                               lower.tail=F)

up.gphendc3.h360a.dn <- phyper(sum(df$gphendc3 == 1 & df$h360a == -1),
                            sum(df$gphendc3 == 1),
                            sum(df$gphendc3 == 0),
                            sum(df$h360a == -1),
                            lower.tail=F)
up.gphendc3.hphendc3.dn <- phyper(sum(df$gphendc3 == 1 & df$hphendc3 == -1),
                               sum(df$gphendc3 == 1),
                               sum(df$gphendc3 == 0),
                               sum(df$hphendc3 == -1),
                               lower.tail=F)
up.gphendc3.h8979a.dn <- phyper(sum(df$gphendc3 == 1 & df$h8979a == -1),
                             sum(df$gphendc3 == 1),
                             sum(df$gphendc3 == 0),
                             sum(df$h8979a == -1),
                             lower.tail=F)

#Down-regulated in Gray PhenDC3
dn.gphendc3.h360a <- phyper(sum(df$gphendc3 == -1 & df$h360a != 0),
                            sum(df$gphendc3 == -1),
                            sum(df$gphendc3 == 0),
                            sum(df$h360a != 0),
                            lower.tail=F)
dn.gphendc3.hphendc3 <- phyper(sum(df$gphendc3 == -1 & df$hphendc3 != 0),
                               sum(df$gphendc3 == -1),
                               sum(df$gphendc3 == 0),
                               sum(df$hphendc3 != 0),
                               lower.tail=F)
dn.gphendc3.h8979a <- phyper(sum(df$gphendc3 == -1 & df$h8979a != 0),
                             sum(df$gphendc3 == -1),
                             sum(df$gphendc3 == 0),
                             sum(df$h8979a != 0),
                             lower.tail=F)
dn.gphendc3.gphendc3 <- phyper(sum(df$gphendc3 == -1 & df$gphendc3 == -1),
                               sum(df$gphendc3 == -1),
                               sum(df$gphendc3 == 0),
                               sum(df$gphendc3 == -1),
                               lower.tail=F)

dn.gphendc3.h360a.up <- phyper(sum(df$gphendc3 == -1 & df$h360a == 1),
                               sum(df$gphendc3 == -1),
                               sum(df$gphendc3 == 0),
                               sum(df$h360a == 1),
                               lower.tail=F)
dn.gphendc3.hphendc3.up <- phyper(sum(df$gphendc3 == -1 & df$hphendc3 == 1),
                                  sum(df$gphendc3 == -1),
                                  sum(df$gphendc3 == 0),
                                  sum(df$hphendc3 == 1),
                                  lower.tail=F)
dn.gphendc3.h8979a.up <- phyper(sum(df$gphendc3 == -1 & df$h8979a == 1),
                                sum(df$gphendc3 == -1),
                                sum(df$gphendc3 == 0),
                                sum(df$h8979a == 1),
                                lower.tail=F)
dn.gphendc3.dn.gphendc3.up <- phyper(sum(df$gphendc3 == -1 & df$gphendc3 == -1),
                                     sum(df$gphendc3 == -1),
                                     sum(df$gphendc3 == 0),
                                     sum(df$gphendc3 == -1),
                                     lower.tail=F)

dn.gphendc3.h360a.dn <- phyper(sum(df$gphendc3 == -1 & df$h360a == -1),
                               sum(df$gphendc3 == -1),
                               sum(df$gphendc3 == 0),
                               sum(df$h360a == -1),
                               lower.tail=F)
dn.gphendc3.hphendc3.dn <- phyper(sum(df$gphendc3 == -1 & df$hphendc3 == -1),
                                  sum(df$gphendc3 == -1),
                                  sum(df$gphendc3 == 0),
                                  sum(df$hphendc3 == -1),
                                  lower.tail=F)
dn.gphendc3.h8979a.dn <- phyper(sum(df$gphendc3 == -1 & df$h8979a == -1),
                                sum(df$gphendc3 == -1),
                                sum(df$gphendc3 == 0),
                                sum(df$h8979a == -1),
                                lower.tail=F)

gray.treatment.summary <- data.frame(Treatment=c("360a.all","360a.up","360a.dn","PhenDC3","PhenDC3.up","PhenDC3.dn","8979a","8979a.up","8979a.dn"),
                                All.Gray.PDC3=c(gphendc3.h360a,gphendc3.h360a.up,gphendc3.h360a.dn,gphendc3.hphendc3,gphendc3.hphendc3.up,gphendc3.hphendc3.dn,gphendc3.h8979a,gphendc3.h8979a.up,gphendc3.h8979a.dn),
                                Up.Gray.PDC3=c(up.gphendc3.h360a,up.gphendc3.h360a.up,up.gphendc3.h360a.dn,up.gphendc3.hphendc3,up.gphendc3.hphendc3.up,up.gphendc3.hphendc3.dn,up.gphendc3.h8979a,up.gphendc3.h8979a.up,up.gphendc3.h8979a.dn),
                                Dn.Gray.PDC3=c(dn.gphendc3.h360a,dn.gphendc3.h360a.up,dn.gphendc3.h360a.dn,dn.gphendc3.hphendc3,dn.gphendc3.hphendc3.up,dn.gphendc3.hphendc3.dn,dn.gphendc3.h8979a,dn.gphendc3.h8979a.up,dn.gphendc3.h8979a.dn))
gray.treatment.summary[,2:ncol(gray.treatment.summary)] <- signif(gray.treatment.summary[,2:ncol(gray.treatment.summary)],3)
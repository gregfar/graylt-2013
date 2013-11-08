args <- commandArgs(trailingOnly = T)
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

# All XPB Peaks
xpb.h360a <- phyper(sum(df$xpb == 1 & df$h360a != 0),
       sum(df$xpb == 1),
       sum(df$xpb == 0),
       sum(df$h360a != 0),
       lower.tail=F)
xpb.hphendc3 <- phyper(sum(df$xpb == 1 & df$hphendc3 != 0),
       sum(df$xpb == 1),
       sum(df$xpb == 0),
       sum(df$hphendc3 != 0),
       lower.tail=F)
xpb.h8979a <- phyper(sum(df$xpb == 1 & df$h8979a != 0),
       sum(df$xpb == 1),
       sum(df$xpb == 0),
       sum(df$h8979a != 0),
       lower.tail=F)
xpb.gphendc3 <- phyper(sum(df$xpb == 1 & df$gphendc3 != 0),
                     sum(df$xpb == 1),
                     sum(df$xpb == 0),
                     sum(df$gphendc3 != 0),
                     lower.tail=F)

xpb.h360a.up <- phyper(sum(df$xpb == 1 & df$h360a == 1),
                    sum(df$xpb == 1),
                    sum(df$xpb == 0),
                    sum(df$h360a == 1),
                    lower.tail=F)
xpb.hphendc3.up <- phyper(sum(df$xpb == 1 & df$hphendc3 == 1),
                       sum(df$xpb == 1),
                       sum(df$xpb == 0),
                       sum(df$hphendc3 == 1),
                       lower.tail=F)
xpb.h8979a.up <- phyper(sum(df$xpb == 1 & df$h8979a == 1),
                     sum(df$xpb == 1),
                     sum(df$xpb == 0),
                     sum(df$h8979a == 1),
                     lower.tail=F)
xpb.gphendc3.up <- phyper(sum(df$xpb == 1 & df$gphendc3 == 1),
                        sum(df$xpb == 1),
                        sum(df$xpb == 0),
                        sum(df$gphendc3 == 1),
                        lower.tail=F)

xpb.h360a.dn <- phyper(sum(df$xpb == 1 & df$h360a == -1),
                       sum(df$xpb == 1),
                       sum(df$xpb == 0),
                       sum(df$h360a == -1),
                       lower.tail=F)
xpb.hphendc3.dn <- phyper(sum(df$xpb == 1 & df$hphendc3 == -1),
                          sum(df$xpb == 1),
                          sum(df$xpb == 0),
                          sum(df$hphendc3 == -1),
                          lower.tail=F)
xpb.h8979a.dn <- phyper(sum(df$xpb == 1 & df$h8979a == -1),
                        sum(df$xpb == 1),
                        sum(df$xpb == 0),
                        sum(df$h8979a == -1),
                        lower.tail=F)
xpb.gphendc3.dn <- phyper(sum(df$xpb == 1 & df$gphendc3 == -1),
                        sum(df$xpb == 1),
                        sum(df$xpb == 0),
                        sum(df$gphendc3 == -1),
                        lower.tail=F)

# All XPD Peaks
xpd.h360a <- phyper(sum(df$xpd == 1 & df$h360a != 0),
                    sum(df$xpd == 1),
                    sum(df$xpd == 0),
                    sum(df$h360a != 0),
                    lower.tail=F)
xpd.hphendc3 <- phyper(sum(df$xpd == 1 & df$hphendc3 != 0),
                       sum(df$xpd == 1),
                       sum(df$xpd == 0),
                       sum(df$hphendc3 != 0),
                       lower.tail=F)
xpd.h8979a <- phyper(sum(df$xpd == 1 & df$h8979a != 0),
                     sum(df$xpd == 1),
                     sum(df$xpd == 0),
                     sum(df$h8979a != 0),
                     lower.tail=F)
xpd.gphendc3 <- phyper(sum(df$xpd == 1 & df$gphendc3 != 0),
                     sum(df$xpd == 1),
                     sum(df$xpd == 0),
                     sum(df$gphendc3 != 0),
                     lower.tail=F)

xpd.h360a.up <- phyper(sum(df$xpd == 1 & df$h360a == 1),
                       sum(df$xpd == 1),
                       sum(df$xpd == 0),
                       sum(df$h360a == 1),
                       lower.tail=F)
xpd.hphendc3.up <- phyper(sum(df$xpd == 1 & df$hphendc3 == 1),
                          sum(df$xpd == 1),
                          sum(df$xpd == 0),
                          sum(df$hphendc3 == 1),
                          lower.tail=F)
xpd.h8979a.up <- phyper(sum(df$xpd == 1 & df$h8979a == 1),
                        sum(df$xpd == 1),
                        sum(df$xpd == 0),
                        sum(df$h8979a == 1),
                        lower.tail=F)
xpd.gphendc3.up <- phyper(sum(df$xpd == 1 & df$gphendc3 == 1),
                        sum(df$xpd == 1),
                        sum(df$xpd == 0),
                        sum(df$gphendc3 == 1),
                        lower.tail=F)

xpd.h360a.dn <- phyper(sum(df$xpd == 1 & df$h360a == -1),
                       sum(df$xpd == 1),
                       sum(df$xpd == 0),
                       sum(df$h360a == -1),
                       lower.tail=F)
xpd.hphendc3.dn <- phyper(sum(df$xpd == 1 & df$hphendc3 == -1),
                          sum(df$xpd == 1),
                          sum(df$xpd == 0),
                          sum(df$hphendc3 == -1),
                          lower.tail=F)
xpd.h8979a.dn <- phyper(sum(df$xpd == 1 & df$h8979a == -1),
                        sum(df$xpd == 1),
                        sum(df$xpd == 0),
                        sum(df$h8979a == -1),
                        lower.tail=F)
xpd.gphendc3.dn <- phyper(sum(df$xpd == 1 & df$gphendc3 == -1),
                        sum(df$xpd == 1),
                        sum(df$xpd == 0),
                        sum(df$gphendc3 == -1),
                        lower.tail=F)

# Both XPB and XPD
both.h360a <- phyper(sum(df$bind.cat == "both" & df$h360a != 0),
                    sum(df$bind.cat == "both"),
                    sum(df$bind.cat != "both"),
                    sum(df$h360a != 0),
                    lower.tail=F)
both.hphendc3 <- phyper(sum(df$bind.cat == "both" & df$hphendc3 != 0),
                     sum(df$bind.cat == "both"),
                     sum(df$bind.cat != "both"),
                     sum(df$hphendc3 != 0),
                     lower.tail=F)
both.h8979a <- phyper(sum(df$bind.cat == "both" & df$h8979a != 0),
                        sum(df$bind.cat == "both"),
                        sum(df$bind.cat != "both"),
                        sum(df$h8979a != 0),
                        lower.tail=F)
both.gphendc3 <- phyper(sum(df$bind.cat == "both" & df$gphendc3 != 0),
                      sum(df$bind.cat == "both"),
                      sum(df$bind.cat != "both"),
                      sum(df$gphendc3 != 0),
                      lower.tail=F)

both.h360a.up <- phyper(sum(df$bind.cat == "both" & df$h360a == 1),
                     sum(df$bind.cat == "both"),
                     sum(df$bind.cat != "both"),
                     sum(df$h360a == 1),
                     lower.tail=F)
both.hphendc3.up <- phyper(sum(df$bind.cat == "both" & df$hphendc3 == 1),
                        sum(df$bind.cat == "both"),
                        sum(df$bind.cat != "both"),
                        sum(df$hphendc3 == 1),
                        lower.tail=F)
both.h8979a.up <- phyper(sum(df$bind.cat == "both" & df$h8979a == 1),
                      sum(df$bind.cat == "both"),
                      sum(df$bind.cat != "both"),
                      sum(df$h8979a == 1),
                      lower.tail=F)
both.gphendc3.up <- phyper(sum(df$bind.cat == "both" & df$gphendc3 == 1),
                         sum(df$bind.cat == "both"),
                         sum(df$bind.cat != "both"),
                         sum(df$gphendc3 == 1),
                         lower.tail=F)

both.h360a.dn <- phyper(sum(df$bind.cat == "both" & df$h360a == -1),
                        sum(df$bind.cat == "both"),
                        sum(df$bind.cat != "both"),
                        sum(df$h360a == -1),
                        lower.tail=F)
both.hphendc3.dn <- phyper(sum(df$bind.cat == "both" & df$hphendc3 == -1),
                           sum(df$bind.cat == "both"),
                           sum(df$bind.cat != "both"),
                           sum(df$hphendc3 == -1),
                           lower.tail=F)
both.h8979a.dn <- phyper(sum(df$bind.cat == "both" & df$h8979a == -1),
                         sum(df$bind.cat == "both"),
                         sum(df$bind.cat != "both"),
                         sum(df$h8979a == -1),
                         lower.tail=F)
both.gphendc3.dn <- phyper(sum(df$bind.cat == "both" & df$gphendc3 == -1),
                         sum(df$bind.cat == "both"),
                         sum(df$bind.cat != "both"),
                         sum(df$gphendc3 == -1),
                         lower.tail=F)


# XPB-Only peaks
xpb.only.h360a <- phyper(sum(df$bind.cat == "xpb.only" & df$h360a != 0),
                     sum(df$bind.cat == "xpb.only"),
                     sum(df$bind.cat != "xpb.only"),
                     sum(df$h360a != 0),
                     lower.tail=F)
xpb.only.hphendc3 <- phyper(sum(df$bind.cat == "xpb.only" & df$hphendc3 != 0),
                        sum(df$bind.cat == "xpb.only"),
                        sum(df$bind.cat != "xpb.only"),
                        sum(df$hphendc3 != 0),
                        lower.tail=F)
xpb.only.h8979a <- phyper(sum(df$bind.cat == "xpb.only" & df$h8979a != 0),
                      sum(df$bind.cat == "xpb.only"),
                      sum(df$bind.cat != "xpb.only"),
                      sum(df$h8979a != 0),
                      lower.tail=F)
xpb.only.gphendc3 <- phyper(sum(df$bind.cat == "xpb.only" & df$gphendc3 != 0),
                          sum(df$bind.cat == "xpb.only"),
                          sum(df$bind.cat != "xpb.only"),
                          sum(df$gphendc3 != 0),
                          lower.tail=F)

xpb.only.h360a.up <- phyper(sum(df$bind.cat == "xpb.only" & df$h360a == 1),
                         sum(df$bind.cat == "xpb.only"),
                         sum(df$bind.cat != "xpb.only"),
                         sum(df$h360a == 1),
                         lower.tail=F)
xpb.only.hphendc3.up <- phyper(sum(df$bind.cat == "xpb.only" & df$hphendc3 == 1),
                            sum(df$bind.cat == "xpb.only"),
                            sum(df$bind.cat != "xpb.only"),
                            sum(df$hphendc3 == 1),
                            lower.tail=F)
xpb.only.h8979a.up <- phyper(sum(df$bind.cat == "xpb.only" & df$h8979a == 1),
                          sum(df$bind.cat == "xpb.only"),
                          sum(df$bind.cat != "xpb.only"),
                          sum(df$h8979a == 1),
                          lower.tail=F)
xpb.only.gphendc3.up <- phyper(sum(df$bind.cat == "xpb.only" & df$gphendc3 == 1),
                             sum(df$bind.cat == "xpb.only"),
                             sum(df$bind.cat != "xpb.only"),
                             sum(df$gphendc3 == 1),
                             lower.tail=F)

xpb.only.h360a.dn <- phyper(sum(df$bind.cat == "xpb.only" & df$h360a == -1),
                            sum(df$bind.cat == "xpb.only"),
                            sum(df$bind.cat != "xpb.only"),
                            sum(df$h360a == -1),
                            lower.tail=F)
xpb.only.hphendc3.dn <- phyper(sum(df$bind.cat == "xpb.only" & df$hphendc3 == -1),
                               sum(df$bind.cat == "xpb.only"),
                               sum(df$bind.cat != "xpb.only"),
                               sum(df$hphendc3 == -1),
                               lower.tail=F)
xpb.only.h8979a.dn <- phyper(sum(df$bind.cat == "xpb.only" & df$h8979a == -1),
                             sum(df$bind.cat == "xpb.only"),
                             sum(df$bind.cat != "xpb.only"),
                             sum(df$h8979a == -1),
                             lower.tail=F)
xpb.only.gphendc3.dn <- phyper(sum(df$bind.cat == "xpb.only" & df$gphendc3 == -1),
                             sum(df$bind.cat == "xpb.only"),
                             sum(df$bind.cat != "xpb.only"),
                             sum(df$gphendc3 == -1),
                             lower.tail=F)

# XPD-Only peaks
xpd.only.h360a <- phyper(sum(df$bind.cat == "xpd.only" & df$h360a != 0),
                         sum(df$bind.cat == "xpd.only"),
                         sum(df$bind.cat != "xpd.only"),
                         sum(df$h360a != 0),
                         lower.tail=F)
xpd.only.hphendc3 <- phyper(sum(df$bind.cat == "xpd.only" & df$hphendc3 != 0),
                            sum(df$bind.cat == "xpd.only"),
                            sum(df$bind.cat != "xpd.only"),
                            sum(df$hphendc3 != 0),
                            lower.tail=F)
xpd.only.h8979a <- phyper(sum(df$bind.cat == "xpd.only" & df$h8979a != 0),
                          sum(df$bind.cat == "xpd.only"),
                          sum(df$bind.cat != "xpd.only"),
                          sum(df$h8979a != 0),
                          lower.tail=F)
xpd.only.gphendc3 <- phyper(sum(df$bind.cat == "xpd.only" & df$gphendc3 != 0),
                          sum(df$bind.cat == "xpd.only"),
                          sum(df$bind.cat != "xpd.only"),
                          sum(df$gphendc3 != 0),
                          lower.tail=F)

xpd.only.h360a.up <- phyper(sum(df$bind.cat == "xpd.only" & df$h360a == 1),
                            sum(df$bind.cat == "xpd.only"),
                            sum(df$bind.cat != "xpd.only"),
                            sum(df$h360a == 1),
                            lower.tail=F)
xpd.only.hphendc3.up <- phyper(sum(df$bind.cat == "xpd.only" & df$hphendc3 == 1),
                               sum(df$bind.cat == "xpd.only"),
                               sum(df$bind.cat != "xpd.only"),
                               sum(df$hphendc3 == 1),
                               lower.tail=F)
xpd.only.h8979a.up <- phyper(sum(df$bind.cat == "xpd.only" & df$h8979a == 1),
                             sum(df$bind.cat == "xpd.only"),
                             sum(df$bind.cat != "xpd.only"),
                             sum(df$h8979a == 1),
                             lower.tail=F)
xpd.only.gphendc3.up <- phyper(sum(df$bind.cat == "xpd.only" & df$gphendc3 == 1),
                             sum(df$bind.cat == "xpd.only"),
                             sum(df$bind.cat != "xpd.only"),
                             sum(df$gphendc3 == 1),
                             lower.tail=F)

xpd.only.h360a.dn <- phyper(sum(df$bind.cat == "xpd.only" & df$h360a == -1),
                            sum(df$bind.cat == "xpd.only"),
                            sum(df$bind.cat != "xpd.only"),
                            sum(df$h360a == -1),
                            lower.tail=F)
xpd.only.hphendc3.dn <- phyper(sum(df$bind.cat == "xpd.only" & df$hphendc3 == -1),
                               sum(df$bind.cat == "xpd.only"),
                               sum(df$bind.cat != "xpd.only"),
                               sum(df$hphendc3 == -1),
                               lower.tail=F)
xpd.only.h8979a.dn <- phyper(sum(df$bind.cat == "xpd.only" & df$h8979a == -1),
                             sum(df$bind.cat == "xpd.only"),
                             sum(df$bind.cat != "xpd.only"),
                             sum(df$h8979a == -1),
                             lower.tail=F)
xpd.only.gphendc3.dn <- phyper(sum(df$bind.cat == "xpd.only" & df$gphendc3 == -1),
                             sum(df$bind.cat == "xpd.only"),
                             sum(df$bind.cat != "xpd.only"),
                             sum(df$gphendc3 == -1),
                             lower.tail=F)

# Unbound genes
unbound.h360a <- phyper(sum(df$bind.cat == "unbound" & df$h360a != 0),
                         sum(df$bind.cat == "unbound"),
                         sum(df$bind.cat != "unbound"),
                         sum(df$h360a != 0),
                         lower.tail=F)
unbound.hphendc3 <- phyper(sum(df$bind.cat == "unbound" & df$hphendc3 != 0),
                            sum(df$bind.cat == "unbound"),
                            sum(df$bind.cat != "unbound"),
                            sum(df$hphendc3 != 0),
                            lower.tail=F)
unbound.h8979a <- phyper(sum(df$bind.cat == "unbound" & df$h8979a != 0),
                          sum(df$bind.cat == "unbound"),
                          sum(df$bind.cat != "unbound"),
                          sum(df$h8979a != 0),
                          lower.tail=F)
unbound.gphendc3 <- phyper(sum(df$bind.cat == "unbound" & df$gphendc3 != 0),
                         sum(df$bind.cat == "unbound"),
                         sum(df$bind.cat != "unbound"),
                         sum(df$gphendc3 != 0),
                         lower.tail=F)

unbound.h360a.up <- phyper(sum(df$bind.cat == "unbound" & df$h360a == 1),
                        sum(df$bind.cat == "unbound"),
                        sum(df$bind.cat != "unbound"),
                        sum(df$h360a == 1),
                        lower.tail=F)
unbound.hphendc3.up <- phyper(sum(df$bind.cat == "unbound" & df$hphendc3 == 1),
                           sum(df$bind.cat == "unbound"),
                           sum(df$bind.cat != "unbound"),
                           sum(df$hphendc3 == 1),
                           lower.tail=F)
unbound.h8979a.up <- phyper(sum(df$bind.cat == "unbound" & df$h8979a == 1),
                         sum(df$bind.cat == "unbound"),
                         sum(df$bind.cat != "unbound"),
                         sum(df$h8979a == 1),
                         lower.tail=F)
unbound.gphendc3.up <- phyper(sum(df$bind.cat == "unbound" & df$gphendc3 == 1),
                            sum(df$bind.cat == "unbound"),
                            sum(df$bind.cat != "unbound"),
                            sum(df$gphendc3 == 1),
                            lower.tail=F)

unbound.h360a.dn <- phyper(sum(df$bind.cat == "unbound" & df$h360a == -1),
                           sum(df$bind.cat == "unbound"),
                           sum(df$bind.cat != "unbound"),
                           sum(df$h360a == -1),
                           lower.tail=F)
unbound.hphendc3.dn <- phyper(sum(df$bind.cat == "unbound" & df$hphendc3 == -1),
                              sum(df$bind.cat == "unbound"),
                              sum(df$bind.cat != "unbound"),
                              sum(df$hphendc3 == -1),
                              lower.tail=F)
unbound.h8979a.dn <- phyper(sum(df$bind.cat == "unbound" & df$h8979a == -1),
                            sum(df$bind.cat == "unbound"),
                            sum(df$bind.cat != "unbound"),
                            sum(df$h8979a == -1),
                            lower.tail=F)
unbound.gphendc3.dn <- phyper(sum(df$bind.cat == "unbound" & df$gphendc3 == -1),
                            sum(df$bind.cat == "unbound"),
                            sum(df$bind.cat != "unbound"),
                            sum(df$gphendc3 == -1),
                            lower.tail=F)

treatment.summary <- data.frame(Treatment=c("360a.all","360a.up","360a.dn","PhenDC3","PhenDC3.up","PhenDC3.dn","8979a","8979a.up","8979a.dn","gPhenDC3","gPhenDC3.up","gPhenDC3.dn"),
                                XPB=c(xpb.h360a,xpb.h360a.up,xpb.h360a.dn,xpb.hphendc3,xpb.hphendc3.up,xpb.hphendc3.dn,xpb.h8979a,xpb.h8979a.up,xpb.h8979a.dn,
                                      xpb.gphendc3,xpb.gphendc3.up,xpb.gphendc3.dn),
                                XPD=c(xpd.h360a,xpd.h360a.up,xpd.h360a.dn,xpd.hphendc3,xpd.hphendc3.up,xpd.hphendc3.dn,xpd.h8979a,xpd.h8979a.up,xpd.h8979a.dn,
                                      xpd.gphendc3,xpd.gphendc3.up,xpd.gphendc3.dn),
                                Both=c(both.h360a,both.h360a.up,both.h360a.dn,both.hphendc3,both.hphendc3.up,both.hphendc3.dn,both.h8979a,both.h8979a.up,both.h8979a.dn,
                                       both.gphendc3,both.gphendc3.up,both.gphendc3.dn),
                                XPB.Only=c(xpb.only.h360a,xpb.only.h360a.up,xpb.only.h360a.dn,xpb.only.hphendc3,xpb.only.hphendc3.up,xpb.only.hphendc3.dn,xpb.only.h8979a,xpb.only.h8979a.up,xpb.only.h8979a.dn,
                                           xpb.only.gphendc3,xpb.only.gphendc3.up,xpb.only.gphendc3.dn),
                                XPD.Only=c(xpd.only.h360a,xpd.only.h360a.up,xpd.only.h360a.dn,xpd.only.hphendc3,xpd.only.hphendc3.up,xpd.only.hphendc3.dn,xpd.only.h8979a,xpd.only.h8979a.up,xpd.only.h8979a.dn,
                                           xpd.only.gphendc3,xpd.only.gphendc3.up,xpd.only.gphendc3.dn),
                                Unbound=c(unbound.h360a,unbound.h360a.up,unbound.h360a.dn,unbound.hphendc3,unbound.hphendc3.up,unbound.hphendc3.dn,unbound.h8979a,unbound.h8979a.up,unbound.h8979a.dn,
                                          unbound.gphendc3,unbound.gphendc3.up,unbound.gphendc3.dn))
treatment.summary[,2:7] <- signif(treatment.summary[,2:7],3)

write.table(treatment.summary,paste(results,"/treatment_hypergeometric.txt",sep=""),sep="\t",quote=F,row.names=F)

xpb.genes <- list(all=df$name[df$xpb == 1],
                  all.360a=df$name[df$xpb == 1 & df$h360a != 0],
                  up.360a=df$name[df$xpb == 1 & df$h360a == 1],
                  dn.360a=df$name[df$xpb == 1 & df$h360a == -1],
                  all.phendc3=df$name[df$xpb == 1 & df$hphendc3 != 0],
                  up.phendc3=df$name[df$xpb == 1 & df$hphendc3 == 1],
                  dn.phendc3=df$name[df$xpb == 1 & df$hphendc3 == -1],
                  all.8979a=df$name[df$xpb == 1 & df$h8979a != 0],
                  up.8979a=df$name[df$xpb == 1 & df$h8979a == 1],
                  dn.8979a=df$name[df$xpb == 1 & df$h8979a == -1])

xpd.genes <- list(all=df$name[df$xpd == 1],
                  all.360a=df$name[df$xpd == 1 & df$h360a != 0],
                  up.360a=df$name[df$xpd == 1 & df$h360a == 1],
                  dn.360a=df$name[df$xpd == 1 & df$h360a == -1],
                  all.phendc3=df$name[df$xpd == 1 & df$hphendc3 != 0],
                  up.phendc3=df$name[df$xpd == 1 & df$hphendc3 == 1],
                  dn.phendc3=df$name[df$xpd == 1 & df$hphendc3 == -1],
                  all.8979a=df$name[df$xpd == 1 & df$h8979a != 0],
                  up.8979a=df$name[df$xpd == 1 & df$h8979a == 1],
                  dn.8979a=df$name[df$xpd == 1 & df$h8979a == -1])

both.genes <- list(all=df$name[df$bind.cat == "both"],
                  all.360a=df$name[df$bind.cat == "both" & df$h360a != 0],
                  up.360a=df$name[df$bind.cat == "both" & df$h360a == 1],
                  dn.360a=df$name[df$bind.cat == "both" & df$h360a == -1],
                  all.phendc3=df$name[df$bind.cat == "both" & df$hphendc3 != 0],
                  up.phendc3=df$name[df$bind.cat == "both" & df$hphendc3 == 1],
                  dn.phendc3=df$name[df$bind.cat == "both" & df$hphendc3 == -1],
                  all.8979a=df$name[df$bind.cat == "both" & df$h8979a != 0],
                  up.8979a=df$name[df$bind.cat == "both" & df$h8979a == 1],
                  dn.8979a=df$name[df$bind.cat == "both" & df$h8979a == -1])

xpb.only.genes <- list(all=df$name[df$bind.cat == "xpb.only"],
                   all.360a=df$name[df$bind.cat == "xpb.only" & df$h360a != 0],
                   up.360a=df$name[df$bind.cat == "xpb.only" & df$h360a == 1],
                   dn.360a=df$name[df$bind.cat == "xpb.only" & df$h360a == -1],
                   all.phendc3=df$name[df$bind.cat == "xpb.only" & df$hphendc3 != 0],
                   up.phendc3=df$name[df$bind.cat == "xpb.only" & df$hphendc3 == 1],
                   dn.phendc3=df$name[df$bind.cat == "xpb.only" & df$hphendc3 == -1],
                   all.8979a=df$name[df$bind.cat == "xpb.only" & df$h8979a != 0],
                   up.8979a=df$name[df$bind.cat == "xpb.only" & df$h8979a == 1],
                   dn.8979a=df$name[df$bind.cat == "xpb.only" & df$h8979a == -1])

xpd.only.genes <- list(all=df$name[df$bind.cat == "xpd.only"],
                       all.360a=df$name[df$bind.cat == "xpd.only" & df$h360a != 0],
                       up.360a=df$name[df$bind.cat == "xpd.only" & df$h360a == 1],
                       dn.360a=df$name[df$bind.cat == "xpd.only" & df$h360a == -1],
                       all.phendc3=df$name[df$bind.cat == "xpd.only" & df$hphendc3 != 0],
                       up.phendc3=df$name[df$bind.cat == "xpd.only" & df$hphendc3 == 1],
                       dn.phendc3=df$name[df$bind.cat == "xpd.only" & df$hphendc3 == -1],
                       all.8979a=df$name[df$bind.cat == "xpd.only" & df$h8979a != 0],
                       up.8979a=df$name[df$bind.cat == "xpd.only" & df$h8979a == 1],
                       dn.8979a=df$name[df$bind.cat == "xpd.only" & df$h8979a == -1])

any.genes  <- list(all=df$name[df$bind.cat != "unbound"],
                   all.360a=df$name[df$bind.cat != "unbound" & df$h360a != 0],
                   up.360a=df$name[df$bind.cat != "unbound" & df$h360a == 1],
                   dn.360a=df$name[df$bind.cat != "unbound" & df$h360a == -1],
                   all.phendc3=df$name[df$bind.cat != "unbound" & df$hphendc3 != 0],
                   up.phendc3=df$name[df$bind.cat != "unbound" & df$hphendc3 == 1],
                   dn.phendc3=df$name[df$bind.cat != "unbound" & df$hphendc3 == -1],
                   all.8979a=df$name[df$bind.cat != "unbound" & df$h8979a != 0],
                   up.8979a=df$name[df$bind.cat != "unbound" & df$h8979a == 1],
                   dn.8979a=df$name[df$bind.cat != "unbound" & df$h8979a == -1])
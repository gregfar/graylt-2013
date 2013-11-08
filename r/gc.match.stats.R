args <- commandArgs(trailingOnly = T)
results <- args[7]

xpb <- read.table(args[1],col.names=c("name","gc"),stringsAsFactors=F)
xpb.g4 <- read.table(args[2],header=T,stringsAsFactors=F)
xpb <- cbind(xpb,g4=xpb.g4$g4,tss=xpb.g4$tss)

xpd <- read.table(args[3],col.names=c("name","gc"),stringsAsFactors=F)
xpd.g4 <- read.table(args[4],header=T,stringsAsFactors=F)
xpd <- cbind(xpd,g4=xpd.g4$g4,tss=xpd.g4$tss)

rand <- read.table(args[5],col.names=c("name","gc"),stringsAsFactors=F)
rand.g4 <- read.table(args[6],header=T,stringsAsFactors=F)

rand <- cbind(rand,g4=rand.g4$g4,tss=rand.g4$tss)

rand <- rand[rand$gc != "N",]

rand$gc <- as.numeric(rand$gc)

cats <- data.frame(g4=c(0,0,1,1),tss=c(0,1,0,1))

t.results <- cbind(cats,xpb=0,xpd=0,rand=0,xpb.rand=0,xpd.rand=0,xpb.xpd=0)
for (i in 1:nrow(cats)) {
  t.results[i,"xpb.rand"] <- t.test(xpb$gc[xpb$g4 == cats$g4[i] & xpb$tss == cats$tss[i]],rand$gc[rand$g4 == cats$g4[i] & rand$tss == cats$tss[i]])[[3]]
  t.results[i,"xpb"] <- mean(xpb$gc[xpb$g4 == cats$g4[i] & xpb$tss == cats$tss[i]])
}
for (i in 1:nrow(cats)) {
  t.results[i,"xpd.rand"] <- t.test(xpd$gc[xpd$g4 == cats$g4[i] & xpd$tss == cats$tss[i]],rand$gc[rand$g4 == cats$g4[i] & rand$tss == cats$tss[i]])[[3]]
  t.results[i,"xpd"] <- mean(xpd$gc[xpd$g4 == cats$g4[i] & xpd$tss == cats$tss[i]])
}
for (i in 1:nrow(cats)) {
  t.results[i,"xpb.xpd"] <- t.test(xpb$gc[xpb$g4 == cats$g4[i] & xpb$tss == cats$tss[i]],xpd$gc[xpd$g4 == cats$g4[i] & xpd$tss == cats$tss[i]])[[3]]
  t.results[i,"rand"] <- mean(rand$gc[rand$g4 == cats$g4[i] & rand$tss == cats$tss[i]])
}
t.results[,3:5] <- t.results[,3:5] * 100
t.results[,3:8] <- signif(t.results[,3:8],3)

write.table(t.results,paste(results,"/gc_percent_t-tests.txt",sep=""),sep="\t",row.names=F,quote=F)

get.gc.match <- function(random,regions) {
  
  ra.df <- random
  ra.df$gc <- round(ra.df$gc,2)
    
  re.df <- regions
  re.df$gc <- round(re.df$gc,2)
  
  tab <- table(re.df$gc)

  result <- data.frame(name=character(),gc=numeric(),g4=integer(),tss=integer())
  
  for(i in 1:length(tab)) {
  
    perc <- as.numeric(names(tab[i]))
    num <- tab[i]
    r.match <- ra.df[ra.df$gc == perc,]
    
    pull <- r.match[sample(nrow(r.match),num),]
    
    result <- rbind(result, pull)
  }

  return(result)
}

rand.xpb.match <- get.gc.match(rand,xpb)
rand.xpd.match <- get.gc.match(rand,xpd)
library(plyr)
match.results <- cbind(cats,xpb=count(xpb,c("g4","tss"))$freq,
                            xpd=count(xpd,c("g4","tss"))$freq,
                            rand.xpb=count(rand.xpb.match,c("g4","tss"))$freq,
                            rand.xpd=count(rand.xpd.match,c("g4","tss"))$freq)

match.results.perc <- match.results
match.results.perc$xpb <- round(match.results$xpb/sum(match.results$xpb)*100,2)
match.results.perc$xpd <- round(match.results$xpd/sum(match.results$xpd)*100,2)
match.results.perc$rand.xpb <- round(match.results$rand.xpb/sum(match.results$rand.xpb)*100,2)
match.results.perc$rand.xpd <- round(match.results$rand.xpd/sum(match.results$rand.xpd)*100,2)

xpb.chisq <- chisq.test(cbind(match.results$xpb,match.results$rand.xpb))
xpd.chisq <- chisq.test(cbind(match.results$xpd,match.results$rand.xpd))

match.results.all <- list(counts=match.results,percents=match.results.perc,xpb.rand.chisq=xpb.chisq,xpd.rand.chisq=xpd.chisq)

write.table(match.results.all$counts,paste(results,"/gc_percent_match_counts.txt",sep=""),sep="\t",quote=F,row.names=F)
write.table(match.results.all$percents,paste(results,"/gc_percent_match_percents.txt",sep=""),sep="\t",quote=F,row.names=F)


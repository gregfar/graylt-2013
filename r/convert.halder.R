library(xlsx)

n <- c("id","name")

p.dn <- read.xlsx("10_treatments/halder_tableS4.xlsx",startRow=3,colIndex=c(1:2),sheetIndex=1,stringsAsFactors=F)
names(p.dn) <- n
p.dn <- p.dn[complete.cases(p.dn),]
write.table(p.dn,"10_treatments/halder_phendc3_dn.txt",sep="\t",row.names=F,quote=F)

p.up <- read.xlsx("10_treatments/halder_tableS4.xlsx",startRow=3,colIndex=c(4:5),sheetIndex=1,stringsAsFactors=F)
names(p.up) <- n
p.up <- p.up[complete.cases(p.up),]
write.table(p.up,"10_treatments/halder_phendc3_up.txt",sep="\t",row.names=F,quote=F)


a.dn <- read.xlsx("10_treatments/halder_tableS4.xlsx",startRow=3,colIndex=c(7:8),sheetIndex=1,stringsAsFactors=F)
names(a.dn) <- n
a.dn <- a.dn[complete.cases(a.dn),]
write.table(a.dn,"10_treatments/halder_360a_dn.txt",sep="\t",row.names=F,quote=F)

a.up <- read.xlsx("10_treatments/halder_tableS4.xlsx",startRow=3,colIndex=c(10:11),sheetIndex=1,stringsAsFactors=F)
names(a.up) <- n
a.up <- a.up[complete.cases(a.up),]
write.table(a.up,"10_treatments/halder_360a_up.txt",sep="\t",row.names=F,quote=F)

c.dn <- read.xlsx("10_treatments/halder_tableS4.xlsx",startRow=3,colIndex=c(13:14),sheetIndex=1,stringsAsFactors=F)
names(c.dn) <- n
c.dn <- c.dn[complete.cases(c.dn),]
write.table(c.dn,"10_treatments/halder_8979a_dn.txt",sep="\t",row.names=F,quote=F)

c.up <- read.xlsx("10_treatments/halder_tableS4.xlsx",startRow=3,colIndex=c(16:17),sheetIndex=1,stringsAsFactors=F)
names(c.up) <- n
c.up <- c.up[complete.cases(c.up),]
write.table(c.up,"10_treatments/halder_8979a_up.txt",sep="\t",row.names=F,quote=F)

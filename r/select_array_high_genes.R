setwd("~/2013_XPB_XPD/")

n <- c("array.id","acc","name","long.name","l.fold","m.fold","h.fold","un","l","m","h","mean","acc2","acc.source","desc")
orb.gene <- read.delim("08_arrays/Gray_010112_Human_ExonST_Results_Gene.txt",header=F,col.names=n,stringsAsFactors=F)

abs <- orb.gene[abs(orb.gene$h.fold) > 1.5, ]
up <- orb.gene[orb.gene$h.fold > 1.5,]
up <- up[!up$name == "---",]
up$name <- sub(" ","",up$name)
up$acc <- sub(" ","",up$acc)

dn <- orb.gene[orb.gene$h.fold < (1/1.5),]
dn <- dn[!dn$name == "---",]
dn$name <- sub(" ","",dn$name)
dn$acc <- sub(" ","",dn$acc)

write.table(up[2:3],"09_treatments/gray_phendc3_up.txt",sep="\t",row.names=F,quote=F)
write.table(dn[2:3],"09_treatments/gray_phendc3_dn.txt",sep="\t",row.names=F,quote=F)
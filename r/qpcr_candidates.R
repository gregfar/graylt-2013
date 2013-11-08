setwd("~/2013_XPB_XPD/")

results <- read.table("10_results/tss_binding_summary.txt",header=T,stringsAsFactors=F)

tss <- read.table("10_results/tss_1kb_g4.txt",header=T,stringsAsFactors=F)

pdc3.up <- read.delim("09_treatments/halder_phendc3_up.txt",header=T,stringsAsFactors=F)
pdc3.up <- pdc3.up[pdc3.up$name != "",]
pdc3.dn <- read.delim("09_treatments/halder_phendc3_dn.txt",header=T,stringsAsFactors=F)
pdc3.dn <- pdc3.dn[pdc3.dn$name != "",]

rod <- read.table("09_treatments/rodriguez.txt",header=T,stringsAsFactors=F)

tss.full <- merge(tss,results,by="name")
tss.full <- cbind(tss.full,pdc3.up=0,pdc3.dn=0,rod=0)

tss.full$pdc3.up[tss.full$name %in% pdc3.up$name] <- 1
tss.full$pdc3.dn[tss.full$name %in% pdc3.dn$name] <- 1
tss.full$rod[tss.full$name %in% rod$name] <- 1

pdc3.rod.candidates <- tss.full[(tss.full$pdc3.up == 1 | tss.full$pdc3.dn == 1) & tss.full$rod == 1,]

candidates.dn <- tss.full[tss.full$bind.cat == "both" & tss.full$nt.overlap > 1 & tss.full$exp.cat < 6 & tss.full$pdc3.dn == 1,]

candidates.up <- tss.full[tss.full$bind.cat == "both" & tss.full$nt.overlap > 1 & tss.full$exp.cat < 6 & tss.full$pdc3.up == 1,]

go.positive <- read.delim("XX_tools/ontologies/GO_0045787_positive_regulation_of_cell_cycle.txt",stringsAsFactors=F,header=F,sep="\t")
names(go.positive)[3] <- "name"
go.negative <- read.delim("XX_tools/ontologies/GO_0045786_negative_regulation_of_cell_cycle.txt",stringsAsFactors=F,header=F,sep="\t")
names(go.negative)[3] <- "name"

candidates.up.go.pos <- candidates.up[candidates.up$name %in% go.positive$name,]
candidates.up.go.neg <- candidates.up[candidates.up$name %in% go.negative$name,]

candidates.dn.go.pos <- candidates.dn[candidates.dn$name %in% go.positive$name,]
candidates.dn.go.neg <- candidates.dn[candidates.dn$name %in% go.negative$name,]

chrx.candidates <- tss.full[tss.full$chr == "chrX" & tss.full$pdc3.up + tss.full$pdc3.dn > 0 & !tss.full$bind.cat == "unbound",]
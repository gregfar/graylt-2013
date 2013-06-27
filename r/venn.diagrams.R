# check to see if the VennDiagram package is installed.
if("VennDiagram" %in% rownames(installed.packages()) == F) {
  install.packages("VennDiagram")
} else {
  library(VennDiagram)
}

xpb.motifs <- read.table("xpb_motif_coincidence_peaks.txt",sep="\t",header=T,
                         colClasses=c(rep("numeric",6)))
xpd.motifs <- read.table("xpd_motif_coincidence_peaks.txt",sep="\t",header=T,
                         colClasses=c(rep("numeric",6)))

par(mar=c(0,0,0,0),oma=c(0,0,0,0))

venn.plot <- draw.quad.venn(
  area1 = sum(xpb.motifs$freq[xpb.motifs$g4 == 1]),
  area2 = sum(xpb.motifs$freq[xpb.motifs$maz == 1]),
  area3 = sum(xpb.motifs$freq[xpb.motifs$ap1 == 1]),
  area4 = sum(xpb.motifs$freq[xpb.motifs$ets == 1]),
  n12 = sum(xpb.motifs$freq[xpb.motifs$g4 == 1 & xpb.motifs$maz == 1]),
  n13 = sum(xpb.motifs$freq[xpb.motifs$g4 == 1 & xpb.motifs$ap1 == 1]),
  n14 = sum(xpb.motifs$freq[xpb.motifs$g4 == 1 & xpb.motifs$ets == 1]),
  n23 = sum(xpb.motifs$freq[xpb.motifs$maz == 1 & xpb.motifs$ap1 == 1]),
  n24 = sum(xpb.motifs$freq[xpb.motifs$maz == 1 & xpb.motifs$ets == 1]),
  n34 = sum(xpb.motifs$freq[xpb.motifs$ap1 == 1 & xpb.motifs$ets == 1]),
  n123 = sum(xpb.motifs$freq[xpb.motifs$g4 == 1 & xpb.motifs$maz == 1 & xpb.motifs$ap1 == 1]),
  n124 = sum(xpb.motifs$freq[xpb.motifs$g4 == 1 & xpb.motifs$maz == 1 & xpb.motifs$ets == 1]),
  n134 = sum(xpb.motifs$freq[xpb.motifs$g4 == 1 & xpb.motifs$ap1 == 1 & xpb.motifs$ets == 1]),
  n234 = sum(xpb.motifs$freq[xpb.motifs$maz == 1 & xpb.motifs$ap1 == 1 & xpb.motifs$ets == 1]),
  n1234 = sum(xpb.motifs$freq[xpb.motifs$g4 == 1 & xpb.motifs$maz == 1 & xpb.motifs$ap1 == 1 & xpb.motifs$ets == 1]),
  category = c("G4", "MAZ", "AP1", "ETS"),
  fill = c("orange", "red", "khaki", "cadetblue1"),
  lty = "solid",
  lwd = "5",
  cex=2.3,
  cat.cex=2.3,
  cat.col = c("orange2", "red2", "khaki", "cadetblue1"),
  fontface = "plain",
  fontfamily = "Helvetica",
  cat.fontfamily = "Helvetica"
)
tiff(filename = "xpb_motif_counts_venn.tiff", width=1200, height = 1200, compression = "lzw", pointsize=24)
plot.new()
grid.draw(venn.plot)
text(0.88,0,paste("None:",xpb.motifs$freq[xpb.motifs$g4 == 0 & xpb.motifs$ap1 == 0 & xpb.motifs$ets == 0 & xpb.motifs$maz == 0],sep=" "),cex=2.3,family="Helvetica")
dev.off()

venn.plot <- draw.quad.venn(
  area1 = sum(xpb.motifs$percent[xpb.motifs$g4 == 1]),
  area2 = sum(xpb.motifs$percent[xpb.motifs$maz == 1]),
  area3 = sum(xpb.motifs$percent[xpb.motifs$ap1 == 1]),
  area4 = sum(xpb.motifs$percent[xpb.motifs$ets == 1]),
  n12 = sum(xpb.motifs$percent[xpb.motifs$g4 == 1 & xpb.motifs$maz == 1]),
  n13 = sum(xpb.motifs$percent[xpb.motifs$g4 == 1 & xpb.motifs$ap1 == 1]),
  n14 = sum(xpb.motifs$percent[xpb.motifs$g4 == 1 & xpb.motifs$ets == 1]),
  n23 = sum(xpb.motifs$percent[xpb.motifs$maz == 1 & xpb.motifs$ap1 == 1]),
  n24 = sum(xpb.motifs$percent[xpb.motifs$maz == 1 & xpb.motifs$ets == 1]),
  n34 = sum(xpb.motifs$percent[xpb.motifs$ap1 == 1 & xpb.motifs$ets == 1]),
  n123 = sum(xpb.motifs$percent[xpb.motifs$g4 == 1 & xpb.motifs$maz == 1 & xpb.motifs$ap1 == 1]),
  n124 = sum(xpb.motifs$percent[xpb.motifs$g4 == 1 & xpb.motifs$maz == 1 & xpb.motifs$ets == 1]),
  n134 = sum(xpb.motifs$percent[xpb.motifs$g4 == 1 & xpb.motifs$ap1 == 1 & xpb.motifs$ets == 1]),
  n234 = sum(xpb.motifs$percent[xpb.motifs$maz == 1 & xpb.motifs$ap1 == 1 & xpb.motifs$ets == 1]),
  n1234 = sum(xpb.motifs$percent[xpb.motifs$g4 == 1 & xpb.motifs$maz == 1 & xpb.motifs$ap1 == 1 & xpb.motifs$ets == 1]),
  category = c("G4", "MAZ", "AP1", "ETS"),
  fill = c("orange", "red", "khaki", "cadetblue1"),
  lty = "solid",
  lwd = "5",
  cex=2.3,
  cat.cex=2.3,
  cat.col = c("orange2", "red2", "khaki", "cadetblue1"),
  fontface = "plain",
  fontfamily = "Helvetica",
  cat.fontfamily = "Helvetica"
)
tiff(filename = "xpb_motif_percents_venn.tiff", width=1200, height = 1200, compression = "lzw", pointsize=24)
plot.new()
grid.draw(venn.plot)
text(0.88,0,paste("None:",xpb.motifs$percent[xpb.motifs$g4 == 0 & xpb.motifs$ap1 == 0 & xpb.motifs$ets == 0 & xpb.motifs$maz == 0],sep=" "),cex=2.3,family="Helvetica")
dev.off()

venn.plot <- draw.quad.venn(
  area1 = sum(xpd.motifs$freq[xpd.motifs$g4 == 1]),
  area2 = sum(xpd.motifs$freq[xpd.motifs$maz == 1]),
  area3 = sum(xpd.motifs$freq[xpd.motifs$ap1 == 1]),
  area4 = sum(xpd.motifs$freq[xpd.motifs$ets == 1]),
  n12 = sum(xpd.motifs$freq[xpd.motifs$g4 == 1 & xpd.motifs$maz == 1]),
  n13 = sum(xpd.motifs$freq[xpd.motifs$g4 == 1 & xpd.motifs$ap1 == 1]),
  n14 = sum(xpd.motifs$freq[xpd.motifs$g4 == 1 & xpd.motifs$ets == 1]),
  n23 = sum(xpd.motifs$freq[xpd.motifs$maz == 1 & xpd.motifs$ap1 == 1]),
  n24 = sum(xpd.motifs$freq[xpd.motifs$maz == 1 & xpd.motifs$ets == 1]),
  n34 = sum(xpd.motifs$freq[xpd.motifs$ap1 == 1 & xpd.motifs$ets == 1]),
  n123 = sum(xpd.motifs$freq[xpd.motifs$g4 == 1 & xpd.motifs$maz == 1 & xpd.motifs$ap1 == 1]),
  n124 = sum(xpd.motifs$freq[xpd.motifs$g4 == 1 & xpd.motifs$maz == 1 & xpd.motifs$ets == 1]),
  n134 = sum(xpd.motifs$freq[xpd.motifs$g4 == 1 & xpd.motifs$ap1 == 1 & xpd.motifs$ets == 1]),
  n234 = sum(xpd.motifs$freq[xpd.motifs$maz == 1 & xpd.motifs$ap1 == 1 & xpd.motifs$ets == 1]),
  n1234 = sum(xpd.motifs$freq[xpd.motifs$g4 == 1 & xpd.motifs$maz == 1 & xpd.motifs$ap1 == 1 & xpd.motifs$ets == 1]),
  category = c("G4", "MAZ", "AP1", "ETS"),
  fill = c("orange", "red", "khaki", "cadetblue1"),
  lty = "solid",
  lwd = "5",
  cex=2.3,
  cat.cex=2.3,
  cat.col = c("orange2", "red2", "khaki", "cadetblue1"),
  fontface = "plain",
  fontfamily = "Helvetica",
  cat.fontfamily = "Helvetica"
)
tiff(filename = "xpd_motif_counts_venn.tiff", width=1200, height = 1200, compression = "lzw", pointsize=24)
plot.new()
grid.draw(venn.plot)
text(0.88,0,paste("None:",xpd.motifs$freq[xpd.motifs$g4 == 0 & xpd.motifs$ap1 == 0 & xpd.motifs$ets == 0 & xpd.motifs$maz == 0],sep=" "),cex=2.3,family="Helvetica")
dev.off()

venn.plot <- draw.quad.venn(
  area1 = round(sum(xpd.motifs$percent[xpd.motifs$g4 == 1]),2),
  area2 = round(sum(xpd.motifs$percent[xpd.motifs$maz == 1]),2),
  area3 = round(sum(xpd.motifs$percent[xpd.motifs$ap1 == 1]),2),
  area4 = round(sum(xpd.motifs$percent[xpd.motifs$ets == 1]),2),
  n12 = round(sum(xpd.motifs$percent[xpd.motifs$g4 == 1 & xpd.motifs$maz == 1]),2),
  n13 = round(sum(xpd.motifs$percent[xpd.motifs$g4 == 1 & xpd.motifs$ap1 == 1]),2),
  n14 = round(sum(xpd.motifs$percent[xpd.motifs$g4 == 1 & xpd.motifs$ets == 1]),2),
  n23 = round(sum(xpd.motifs$percent[xpd.motifs$maz == 1 & xpd.motifs$ap1 == 1]),2),
  n24 = round(sum(xpd.motifs$percent[xpd.motifs$maz == 1 & xpd.motifs$ets == 1]),2),
  n34 = round(sum(xpd.motifs$percent[xpd.motifs$ap1 == 1 & xpd.motifs$ets == 1]),2),
  n123 = round(sum(xpd.motifs$percent[xpd.motifs$g4 == 1 & xpd.motifs$maz == 1 & xpd.motifs$ap1 == 1]),2),
  n124 = round(sum(xpd.motifs$percent[xpd.motifs$g4 == 1 & xpd.motifs$maz == 1 & xpd.motifs$ets == 1]),2),
  n134 = round(sum(xpd.motifs$percent[xpd.motifs$g4 == 1 & xpd.motifs$ap1 == 1 & xpd.motifs$ets == 1]),2),
  n234 = round(sum(xpd.motifs$percent[xpd.motifs$maz == 1 & xpd.motifs$ap1 == 1 & xpd.motifs$ets == 1]),2),
  n1234 = round(sum(xpd.motifs$percent[xpd.motifs$g4 == 1 & xpd.motifs$maz == 1 & xpd.motifs$ap1 == 1 & xpd.motifs$ets == 1]),2),
  category = c("G4", "MAZ", "AP1", "ETS"),
  fill = c("orange", "red", "khaki", "cadetblue1"),
  lty = "solid",
  lwd = "5",
  cex=2.3,
  cat.cex=2.3,
  cat.col = c("orange2", "red2", "khaki", "cadetblue1"),
  fontface = "plain",
  fontfamily = "Helvetica",
  cat.fontfamily = "Helvetica"
)
tiff(filename = "xpd_motif_percents_venn.tiff", width=1200, height = 1200, compression = "lzw", pointsize=24)
plot.new()
grid.draw(venn.plot)
text(0.88,0,paste("None:",xpd.motifs$percent[xpd.motifs$g4 == 0 & xpd.motifs$ap1 == 0 & xpd.motifs$ets == 0 & xpd.motifs$maz == 0],sep=" "),cex=2.3,family="Helvetica")
dev.off()
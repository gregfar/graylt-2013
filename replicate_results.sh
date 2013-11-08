# These variables define the locations of the directories for scripts and the various
# programs required to replicate this data analysis.

# R binary location
R=Rscript
# Bowtie directory location
BOWTIEDIR=~/tools/bowtie-0.12.7
# IGVTools binary location
IGVTOOLS=~/tools/IGVTools/igvtools
# MACS binary location (v1.40beta required to exactly reproduce results)
MACS=macs14
# CEAS binary location
CEAS=ceas
# UCSC Kent Source Tree is required for use of the bedGraphToBigWig script.
# UCSC bedGraphToBigWig binary location
BGTOBW=~/bin/bedGraphToBigWig
# UCSC twoBitToFa
TBTOFA=~/bin/twoBitToFa
# QuadParser binary location
QP=~/quadparser/quadparser
# Great calculateBinomialP tool binary location
GREAT=~/tools/greatTools/calculateBinomialP

#Directory structure
TOOLS=tools
mkdir $TOOLS

# Perl script directory
# BioPerl Bio::DB::BigWig module is required for some perl scripts
PERLDIR=perl
# R script directory
RDIR=r

READS=01_reads
mkdir $READS
MAPPED=02_mapped
mkdir $MAPPED
OVERLAPS=03_overlaps
mkdir $OVERLAPS
PEAKS=04_peaks
mkdir $PEAKS
mkdir $PEAKS/XPB
mkdir $PEAKS/XPD

XPB_PEAKS=$PEAKS/XPB/xpb_1e-15_peaks.bed
XPB_SUMMITS=$PEAKS/XPB/xpb_summits.bed
XPB_REGIONS=$PEAKS/XPB/xpb_regions.bed

XPD_PEAKS=$PEAKS/XPD/xpd_1e-15_peaks.bed
XPD_SUMMITS=$PEAKS/XPD/xpd_summits.bed
XPD_REGIONS=$PEAKS/XPD/xpd_regions.bed

CEAS_DIR=05_ceas
mkdir $CEAS_DIR

QP_DIR=06_quadparser
mkdir $QP_DIR

G4=$QP_DIR/g4-12.bed

ARRAYS=07_arrays
mkdir $ARRAYS

PILEUPS=08_pileups
mkdir $PILEUPS

ARRAYS=09_arrays
mkdir $ARRAYS

TREATMENTS=10_treatments
mkdir $TREATMENTS

RESULTS=00_results
mkdir $RESULTS

# Original files
# Input: input_1.fastq and input_2.fastq
# XPB: xpb_chip_1.fastq and xpb_chip_2.fastq
# XPD: xpd_chip_1.fastq and xpd_chip_2.fastq


## Read quality control
# Trim to 60 reads
$PERLDIR/fastq_filters/trim.pl $READS/input_1.fastq $READS/input_1.trim.fastq 60
$PERLDIR/fastq_filters/trim.pl $READS/input_2.fastq $READS/input_2.trim.fastq 60
$PERLDIR/fastq_filters/trim.pl $READS/xpb_chip_1.fastq $READS/xpb_chip_1.trim.fastq 60
$PERLDIR/fastq_filters/trim.pl $READS/xpb_chip_2.fastq $READS/xpb_chip_2.trim.fastq 60
$PERLDIR/fastq_filters/trim.pl $READS/xpd_chip_1.fastq $READS/xpd_chip_1.trim.fastq 60
$PERLDIR/fastq_filters/trim.pl $READS/xpd_chip_2.fastq $READS/xpd_chip_2.trim.fastq 60

# Delete the raw files to save space
rm $READS/input_1.fastq
rm $READS/input_2.fastq
rm $READS/xpb_chip_1.fastq
rm $READS/xpb_chip_2.fastq
rm $READS/xpd_chip_1.fastq
rm $READS/xpd_chip_2.fastq

# Remove pairs with 3 or more reads with a PHRED score below 20
$PERLDIR/fastq_filters/paired_min_qc_filter.pl $READS/input_1.trim.fastq $READS/input_1.filter.trim.fastq $READS/input_2.trim.fastq $READS/input_2.filter.trim.fastq 20 3
$PERLDIR/fastq_filters/paired_min_qc_filter.pl $READS/xpb_chip_1.trim.fastq $READS/xpb_chip_1.filter.trim.fastq $READS/xpb_chip_2.trim.fastq $READS/xpb_chip_2.filter.trim.fastq 20 3
$PERLDIR/fastq_filters/paired_min_qc_filter.pl $READS/xpd_chip_1.trim.fastq $READS/xpd_chip_1.filter.trim.fastq $READS/xpd_chip_2.trim.fastq $READS/xpd_chip_2.filter.trim.fastq 20 3

# Delete unfiltered files to save space
rm $READS/input_1.trim.fastq
rm $READS/input_2.trim.fastq
rm $READS/xpb_chip_1.trim.fastq
rm $READS/xpb_chip_2.trim.fastq
rm $READS/xpd_chip_1.trim.fastq
rm $READS/xpd_chip_2.trim.fastq

## Alignment of reads to the hg19 human genome using Bowtie v0.12.7
## Bowtie and the hg19 index were obtained from http://bowtie-bio.sourceforge.net/index.shtml
# Bowtie commands for SAM alignments (used to make bigWig files)
$BOWTIEDIR/bowtie -t -p 3 -n 0 -m 1 --best -S --chunkmbs 256 hg19 -1 $READS/input_1.filter.trim.fastq -2 $READS/input_2.filter.trim.fastq $MAPPED/input_hg19.sam
$BOWTIEDIR/bowtie -t -p 3 -n 0 -m 1 --best -S --chunkmbs 256 hg19 -1 $READS/xpb_chip_1.filter.trim.fastq -2 $READS/xpb_chip_2.filter.trim.fastq $MAPPED/xpb_hg19.sam
$BOWTIEDIR/bowtie -t -p 3 -n 0 -m 1 --best -S --chunkmbs 256 hg19 -1 $READS/xpd_chip_1.filter.trim.fastq -2 $READS/xpd_chip_2.filter.trim.fastq $MAPPED/xpd_hg19.sam

##Filter SAM files to remove unmapped reads or pairs with only 1 mapped mate
# hg19 alignments
$PERLDIR/sam_remove_unmapped_pairs.pl $MAPPED/input_hg19.sam $MAPPED/input_hg19.mapped.sam
$PERLDIR/sam_remove_unmapped_pairs.pl $MAPPED/xpb_hg19.sam $MAPPED/xpb_hg19.mapped.sam
$PERLDIR/sam_remove_unmapped_pairs.pl $MAPPED/xpd_hg19.sam $MAPPED/xpd_hg19.mapped.sam

# Delete files with unmapped reads to save space
rm $MAPPED/input_hg19.sam
rm $MAPPED/xpb_hg19.sam
rm $MAPPED/xpd_hg19.sam

# Select Input reads to match the number of mapped XPB and XPD reads
# Since alignments don't sort the read positions, the reads will be random with respect to genomic position
# and I can simply use head to get an equal number.
head -34399081 $MAPPED/input_hg19.mapped.sam > $MAPPED/xpb_input_hg19.mapped.sam
head -43952749 $MAPPED/input_hg19.mapped.sam > $MAPPED/xpd_input_hg19.mapped.sam

# Sort SAM files using IGVTools
# IGVTools are available from http://www.broadinstitute.org/igv/igvtools
$IGVTOOLS sort -m 2500000 $MAPPED/xpb_input_hg19.mapped.sam $MAPPED/xpb_input_hg19.sorted.mapped.sam
$IGVTOOLS sort -m 2500000 $MAPPED/xpd_input_hg19.mapped.sam $MAPPED/xpd_input_hg19.sorted.mapped.sam
$IGVTOOLS sort -m 2500000 $MAPPED/xpb_hg19.mapped.sam $MAPPED/xpb_hg19.sorted.mapped.sam
$IGVTOOLS sort -m 2500000 $MAPPED/xpd_hg19.mapped.sam $MAPPED/xpd_hg19.sorted.mapped.sam

# Delete unsorted files to save space
rm $MAPPED/xpb_input_hg19.mapped.sam
rm $MAPPED/xpd_input_hg19.mapped.sam
rm $MAPPED/xpb_hg19.mapped.sam
rm $MAPPED/xpd_hg19.mapped.sam

# Covert sorted SAM files to BEDGraph files
$PERLDIR/sam_sorted_to_bedgraph_faster.pl $MAPPED/xpb_input_hg19.sorted.mapped.sam $OVERLAPS/xpb_input.bedgraph 60 0
$PERLDIR/sam_sorted_to_bedgraph_faster.pl $MAPPED/xpd_input_hg19.sorted.mapped.sam $OVERLAPS/xpd_input.bedgraph 60 0
$PERLDIR/sam_sorted_to_bedgraph_faster.pl $MAPPED/xpb_hg19.sorted.mapped.sam $OVERLAPS/xpb.bedgraph 60 0
$PERLDIR/sam_sorted_to_bedgraph_faster.pl $MAPPED/xpd_hg19.sorted.mapped.sam $OVERLAPS/xpd.bedgraph 60 0

# Convert BEDGraph files to bigWig using UCSC Tools
# UCSC command line tools can be obtained from http://hgdownload.cse.ucsc.edu/admin/exe/
curl -L "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz" > $TOOLS/chromInfo.txt.gz
gzip -d $TOOLS/chromInfo.txt.gz

$BGTOBW $OVERLAPS/xpb_input.bedgraph $TOOLS/chromInfo.txt $OVERLAPS/xpb_input.bw
$BGTOBW $OVERLAPS/xpd_input.bedgraph $TOOLS/chromInfo.txt $OVERLAPS/xpd_input.bw
$BGTOBW $OVERLAPS/xpb.bedgraph $TOOLS/chromInfo.txt $OVERLAPS/xpb.bw
$BGTOBW $OVERLAPS/xpd.bedgraph $TOOLS/chromInfo.txt $OVERLAPS/xpd.bw

## Peak finding with MACS v1.4.2 using the Bowtie map files, above
# MACS can be obtained from http://liulab.dfci.harvard.edu/MACS/

# First, must convert the SAM files to Bowtie format, because MACS has a problem recognizing
# which strand reads in SAM format are on and throws an error.
$PERLDIR/sam_to_btm.pl $MAPPED/xpb_input_hg19.sorted.mapped.sam $MAPPED/xpb_input_hg19.btm
$PERLDIR/sam_to_btm.pl $MAPPED/xpd_input_hg19.sorted.mapped.sam $MAPPED/xpd_input_hg19.btm
$PERLDIR/sam_to_btm.pl $MAPPED/xpb_hg19.sorted.mapped.sam $MAPPED/xpb_hg19.btm
$PERLDIR/sam_to_btm.pl $MAPPED/xpd_hg19.sorted.mapped.sam $MAPPED/xpd_hg19.btm

# MACS peak finding
cd $PEAKS/XPB
$MACS -f BOWTIE -s 60 -t ../../$MAPPED/xpb_hg19.btm -c ../../$MAPPED/xpb_input_hg19.btm -n xpb_macs
cd ../../
cd $PEAKS/XPD
$MACS -f BOWTIE -s 60 -t ../../$MAPPED/xpd_hg19.btm -c ../../$MAPPED/xpd_input_hg19.btm -n xpd_macs
cd ../../

# Filter only peaks with a p-value < 1e-15
$PERLDIR/column_filter.pl $PEAKS/XPB/xpb_macs_peaks.bed $XPB_PEAKS 5 gt 150	
$PERLDIR/column_filter.pl $PEAKS/XPD/xpd_macs_peaks.bed $XPD_PEAKS 5 gt 150

# Sort peaks in chromosome and position order
$PERLDIR/bed_sort.pl $XPB_PEAKS temp.bed
mv temp.bed $XPB_PEAKS
$PERLDIR/bed_sort.pl $XPD_PEAKS temp.bed
mv temp.bed $XPD_PEAKS
rm temp.bed

# Correct MACS BED format to standard 6-column BED format and renumber peaks.
$PERLDIR/macs_to_bed.pl $XPB_PEAKS temp.bed xpb
mv temp.bed $XPB_PEAKS
$PERLDIR/macs_to_bed.pl $XPD_PEAKS temp.bed xpd
mv temp.bed $XPD_PEAKS

# Find peak summits using a Perl script that utilizes the Bio::DB::BigWig v1.01 BioPerl module
# Bio::DB::BigWig can be obtained at http://search.cpan.org/~lds/Bio-BigFile-1.01/lib/Bio/DB/BigWig.pm
$PERLDIR/bigwig_bed_summits.pl $OVERLAPS/xpb.bw $XPB_PEAKS $PEAKS/XPB/xpb_summits.txt
$PERLDIR/bigwig_bed_summits.pl $OVERLAPS/xpd.bw $XPD_PEAKS $PEAKS/XPD/xpd_summits.txt

# Make summit and summit region BED files
$PERLDIR/summit_to_bed.pl $PEAKS/XPB/xpb_summits.txt $XPB_SUMMITS 0 1
$PERLDIR/summit_to_bed.pl $PEAKS/XPD/xpd_summits.txt $XPD_SUMMITS 0 1
$PERLDIR/summit_to_bed.pl $PEAKS/XPB/xpb_summits.txt $XPB_REGIONS 50 50
$PERLDIR/summit_to_bed.pl $PEAKS/XPD/xpd_summits.txt $XPD_REGIONS 50 50

## Overlapping and non-overlapping XPB and XPD peaks
# XPB summit regions were used for "Both" summit regions.
$PERLDIR/bed_overlap.pl $XPB_PEAKS $XPD_PEAKS $PEAKS/both_peaks.bed
$PERLDIR/bed_subtract.pl $XPB_PEAKS $XPD_PEAKS $PEAKS/xpb_only_peaks.bed
$PERLDIR/bed_subtract.pl $XPD_PEAKS $XPB_PEAKS $PEAKS/xpd_only_peaks.bed
$PERLDIR/bed_overlap.pl $XPB_REGIONS $PEAKS/both_peaks.bed $PEAKS/both_regions.bed
$PERLDIR/bed_overlap.pl $XPB_REGIONS $PEAKS/xpb_only_peaks.bed $PEAKS/xpb_only_regions.bed
$PERLDIR/bed_overlap.pl $XPD_REGIONS $PEAKS/xpd_only_peaks.bed $PEAKS/xpd_only_regions.bed

## CEAS analysis of XPB and XPD peaks
# CEAS was obtained from http://liulab.dfci.harvard.edu/CEAS/
# CEAS requires a refGene table
curl -L "http://liulab.dfci.harvard.edu/CEAS/src/hg19.refGene.gz" > $TOOLS/hg19.refGene.gz
gzip -d $TOOLS/hg19.refGene.gz

#calculate CEAS results
cd $CEAS_DIR
$CEAS -b ../$XPB_PEAKS -g ../$TOOLS/hg19.refGene --name=xpb
$CEAS -b ../$XPD_PEAKS -g ../$TOOLS/hg19.refGene --name=xpd

# Generate the CEAS report PDFs using R
$R xpb.R
$R xpd.R
cd ../

## Quadparser searches in the hg19 genome
# Quadparser was downloaded from http://www-shankar.ch.cam.ac.uk/Applications/quadparser/
# The following shell script was then used to search each chromosome FASTA for G4 motifs,
# then compile and sort the results into a single BED-formatted file.

# Download the hg19 genome in FASTA format:
mkdir $QP_DIR/hg19_chromFa
curl -L "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz" > $QP_DIR/hg19_chromFa/hg19_chromFa.tar.gz
tar xvfz $QP_DIR/hg19_chromFa/hg19_chromFa.tar.gz

# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
# This resulted in one fasta file for each chromosome in a directory named hg19_chromFa.
# These files had to be converted to all caps using this shell script:
	
mkdir $QP_DIR/hg19_chromFa_allcaps

for f in $QP_DIR/hg19_chromFa/*.fa
do
	NOPATH=${f##*/}
	FILENAME=${NOPATH%.[^.]*}
	head -1 $f > $QP_DIR/hg19_chromFa_allcaps/$NOPATH
	tail -n +2 $f | tr '[a-z]' '[A-Z]' >> $QP_DIR/hg19_chromFa_allcaps/$NOPATH
done

mkdir $QP_DIR/hg19_qp
	
for f in $QP_DIR/hg19_chromFa_allcaps/*
do
	NOPATH=${f##*/}
	FILENAME=${NOPATH%.[^.]*}
	$QP -DAS $f GC 3 4 1 12 $QP_DIR/hg19_qp/$FILENAME.bed
	$PERLDIR/remove_lines.pl $QP_DIR/hg19_qp/$FILENAME.bed $QP_DIR/hg19_qp/temp.bed 1 track\ name=Quadruplexes
	$PERLDIR/convert_column.pl $QP_DIR/hg19_qp/temp.bed $QP_DIR/hg19_qp/$FILENAME.bed 1 $FILENAME
done
rm $QP_DIR/hg19_qp/temp.bed

cat $QP_DIR/hg19_qp/*.bed > $QP_DIR/hg19_qp/temp.bed
$PERLDIR/bed_sort.pl $QP_DIR/hg19_qp/temp.bed $G4
rm $QP_DIR/hg19_qp/temp.bed

#also need a concatenated version of the G4 file in order to calculate coverage later
$PERLDIR/cat_overlapping_regions.pl $G4 $QP_DIR/cat_g4-12.bed

# Convert g4-12.bed to bedgraph
$PERLDIR/bed_to_bedgraph.pl $G4 $QP_DIR/g4-12.bedgraph
$PERLDIR/cat_overlapping_bedgraph.pl $QP_DIR/g4-12.bedgraph $QP_DIR/temp.bedgraph
mv $QP_DIR/temp.bedgraph $QP_DIR/g4-12.bedgraph

# Make separate bedgraphs for + and - strands
$PERLDIR/get_lines.pl $G4 $QP_DIR/g4-12_plus.bed 6 +
$PERLDIR/bed_to_bedgraph.pl $QP_DIR/g4-12_plus.bed $QP_DIR/g4-12_plus.bedgraph

$PERLDIR/get_lines.pl $G4 $QP_DIR/g4-12_minus.bed 6 -
$PERLDIR/bed_to_bedgraph.pl $QP_DIR/g4-12_minus.bed $QP_DIR/g4-12_minus.bedgraph

# Remove chrM, chrUn and chr##_ non-standard chromosomes.
egrep -v "chr[0-9]+_.+$" $QP_DIR/g4-12.bedgraph > $QP_DIR/temp.bedgraph
mv $QP_DIR/temp.bedgraph $QP_DIR/g4-12.bedgraph
egrep -v "chr[M|U].+$" $QP_DIR/g4-12.bedgraph > $QP_DIR/temp.bedgraph
mv $QP_DIR/temp.bedgraph $QP_DIR/g4-12.bedgraph
egrep -v "chr[0-9]+_.+$" $QP_DIR/g4-12_plus.bedgraph > $QP_DIR/temp.bedgraph
mv $QP_DIR/temp.bedgraph $QP_DIR/g4-12_plus.bedgraph
egrep -v "chr[M|U].+$" $QP_DIR/g4-12_plus.bedgraph > $QP_DIR/temp.bedgraph
mv $QP_DIR/temp.bedgraph $QP_DIR/4-12_plus.bedgraph
egrep -v "chr[0-9]+_.+$" $QP_DIR/g4-12_minus.bedgraph > $QP_DIR/temp.bedgraph
mv $QP_DIR/temp.bedgraph $QP_DIR/g4-12_minus.bedgraph
egrep -v "chr[M|U].+$" $QP_DIR/g4-12_minus.bedgraph > $QP_DIR/temp.bedgraph
mv $QP_DIR/temp.bedgraph $QP_DIR/g4-12_minus.bedgraph

# Convert bedgraph files to bigwig for frequency maps
$BGTOBW $QP_DIR/g4-12.bedgraph $TOOLS/chromInfo.txt $QP_DIR/g4-12.bw
$BGTOBW $QP_DIR/g4-12_plus.bedgraph $TOOLS/chromInfo.txt $QP_DIR/g4-12_plus.bw
$BGTOBW $QP_DIR/g4-12_minus.bedgraph $TOOLS/chromInfo.txt $QP_DIR/g4-12_minus.bw

# Then, I used the GREAT calculateBinomialP tool to calculate significance of overlaps between
# G4 and XPB/XPD peaks.

$PERLDIR/bed_to_great.pl $XPB_PEAKS $PEAKS/XPB/xpb_peaks.rd
$PERLDIR/bed_to_great.pl $XPD_PEAKS $PEAKS/XPD/xpd_peaks.rd
$PERLDIR/bed_to_great.pl $XPB_REGIONS $PEAKS/XPB/xpb_regions.rd
$PERLDIR/bed_to_great.pl $XPD_REGIONS $PEAKS/XPD/xpd_regions.rd
	
echo -e "target_name\ttarget_count\tprobe_name\tprobe_count\toverlap_count\tp-value" > $RESULTS/g4_great_results.txt
STATS="$($PERLDIR/great_overlaps.pl $PEAKS/XPB/xpb_peaks.rd $G4)"
PVALUE="$($GREAT $PEAKS/XPB/xpb_peaks.rd $TOOLS/antigap.bed $STATS)"
echo -e "xpb_peaks\t26051\tg4-12\t$STATS\t$PVALUE" >> $RESULTS/g4_great_results.txt

STATS="$($PERLDIR/great_overlaps.pl $PEAKS/XPD/xpd_peaks.rd $G4)"
PVALUE="$($GREAT $PEAKS/XPD/xpd_peaks.rd $TOOLS/antigap.bed $STATS)"
echo -e "xpd_peaks\t16730\tg4-12\t$STATS\t$PVALUE" >> $RESULTS/g4_great_results.txt

STATS="$($PERLDIR/great_overlaps.pl $PEAKS/XPB/xpb_regions.rd $G4)"
PVALUE="$($GREAT $PEAKS/XPB/xpb_regions.rd $TOOLS/antigap.bed $STATS)"
echo -e "xpb_regions\t26051\tg4-12\t$STATS\t$PVALUE" >> $RESULTS/g4_great_results.txt

STATS="$($PERLDIR/great_overlaps.pl $PEAKS/XPD/xpd_regions.rd $G4)"
PVALUE="$($GREAT $PEAKS/XPD/xpd_regions.rd $TOOLS/antigap.bed $STATS)"
echo -e "xpd_regions\t16730\tg4-12\t$STATS\t$PVALUE" >> $RESULTS/g4_great_results.txt

## RefSeq TSS region selection
# The RefSeq gene table was obtained from the UCSC genome browser:
curl -L "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz" > $TOOLS/refGene.txt.gz
gzip -d $TOOLS/refGene.txt.gz

# Get TSS +/- 1kb and TSS +/- 100bp
$PERLDIR/refseq_tss_region.pl $TOOLS/refGene.txt $TOOLS/tss_1kb.bed 1000
$PERLDIR/refseq_tss_region.pl $TOOLS/refGene.txt $TOOLS/tss_100bp.bed 100

# Use R to get just unique lines of the TSS files:
$R $RDIR/unique.tss.R $TOOLS

# Calculate G4 overlaps with TSS regions, XPB and XPD peaks, and XPB and XPD regions	
$PERLDIR/bed_overlap_binary.pl $G4 $TOOLS/tss_1kb.bed temp.txt
mv temp.txt $RESULTS/g4-12_overlaps.txt
$PERLDIR/bed_overlap_binary.pl $RESULTS/g4-12_overlaps.txt $XPB_PEAKS temp.txt
mv temp.txt $RESULTS/g4-12_overlaps.txt
$PERLDIR/bed_overlap_binary.pl $RESULTS/g4-12_overlaps.txt $XPD_PEAKS temp.txt
mv temp.txt $RESULTS/g4-12_overlaps.txt
$PERLDIR/bed_overlap_binary.pl $RESULTS/g4-12_overlaps.txt $XPB_REGIONS temp.txt
mv temp.txt $RESULTS/g4-12_overlaps.txt
$PERLDIR/bed_overlap_binary.pl $RESULTS/g4-12_overlaps.txt $XPD_REGIONS temp.txt
head temp.txt
echo -e "chr\tstart\tend\tname\tscore\tstrand\ttss\txpb.peaks\txpd.peaks\txpb.regions\txpd.regions" > temp2.txt
cat temp2.txt temp.txt > $RESULTS/g4-12_overlaps.txt

rm temp.txt
rm temp2.txt

$PERLDIR/bed_overlap_binary.pl $XPB_PEAKS $G4 temp.txt
echo -e "chr\tstart\tend\tname\tscore\tstrand\tg4" > temp2.txt
cat temp2.txt temp.txt > $RESULTS/xpb_peaks_g4.txt
$PERLDIR/bed_overlap_binary.pl $XPD_PEAKS $G4 temp.txt
echo -e "chr\tstart\tend\tname\tscore\tstrand\tg4" > temp2.txt
cat temp2.txt temp.txt > $RESULTS/xpd_peaks_g4.txt
$PERLDIR/bed_overlap_binary.pl $XPB_REGIONS $G4 temp.txt
echo -e "chr\tstart\tend\tname\tscore\tstrand\tg4" > temp2.txt
cat temp2.txt temp.txt > $RESULTS/xpb_regions_g4.txt
$PERLDIR/bed_overlap_binary.pl $XPD_REGIONS $G4 temp.txt
echo -e "chr\tstart\tend\tname\tscore\tstrand\tg4" > temp2.txt
cat temp2.txt temp.txt > $RESULTS/xpd_regions_g4.txt

rm temp.txt
rm temp2.txt

$PERLDIR/bigwig_bed_overlap_stranded.pl $QP_DIR/g4-12_plus.bw $QP_DIR/g4-12_minus.bw $TOOLS/tss_1kb.bed $RESULTS/tss_1kb_g4.txt

# Calculate the frequency of TSS and G4 overlaps for XPB and XPD, as well as coverage calculations.
$R $RDIR/g4.tss.stats.R $TOOLS/tss_1kb.bed $RESULTS/xpb_peaks_g4.txt $RESULTS/xpd_peaks_g4.txt $RESULTS/xpb_regions_g4.txt $RESULTS/xpd_regions_g4.txt $XPB_SUMMITS $XPD_SUMMITS $RESULTS/g4-12_overlaps.txt $QP_DIR/cat_g4-12.bed $TOOLS/antigap.bed $RESULTS

# I generated tables of the locations of G4 motifs in XPB and XPD summit regions:
$PERLDIR/g4_results_table.pl $G4 $XPB_REGIONS $RESULTS/g4_xpb_regions.txt
$PERLDIR/g4_results_table.pl $G4 $XPD_REGIONS $RESULTS/g4_xpd_regions.txt
$PERLDIR/g4_results_table.pl $G4 $XPB_PEAKS $RESULTS/g4_xpb_peaks.txt
$PERLDIR/g4_results_table.pl $G4 $XPD_PEAKS $RESULTS/g4_xpd_peaks.txt
	
# XPB and XPD summit regions containing G4 motifs were selected using R for use with the online version of GREAT
$R $RDIR/motif.peaks.R $XPB_SUMMITS $XPD_SUMMITS $RESULTS/g4_xpb_regions.txt $RESULTS/g4_xpd_regions.txt $RESULTS/g4_xpb_peaks.txt $RESULTS/g4_xpb_peaks.txt $RESULTS

# Tabulate percent of peaks in TSS +/- 1kb regions, as well as peaks in RefSeq gene bodies:

$PERLDIR/refseq_tss_region.pl $TOOLS/refGene.txt $TOOLS/tss_1kb_all.bed 1000
$R $RDIR/peak.tss.gene.R $TOOLS/refGene.txt $TOOLS/tss_1kb_all.bed $XPB_SUMMITS $XPD_SUMMITS $RESULTS
	
# Calculate significance of enrichment of XPB and XPD summits near TSS:
	
$PERLDIR/bed_to_great.pl $TOOLS/tss_1kb.bed $TOOLS/tss_1kb.rd
$PERLDIR/bed_to_great.pl $TOOLS/tss_100bp.bed $TOOLS/tss_100bp.rd
	
echo -e "target_name\ttarget_count\tprobe_name\tprobe_count\toverlap_count\tp-value" > $RESULTS/tss_great_results.txt
STATS="$($PERLDIR/great_overlaps.pl $TOOLS/tss_1kb.rd $XPB_SUMMITS)"
PVALUE="$($GREAT $TOOLS/tss_1kb.rd $TOOLS/antigap.bed $STATS)"
echo -e "tss_1kb\t31069\txpb_summits\t$STATS\t$PVALUE" >> $RESULTS/tss_great_results.txt

STATS="$($PERLDIR/great_overlaps.pl $TOOLS/tss_1kb.rd $XPD_SUMMITS)"
PVALUE="$($GREAT $TOOLS/tss_1kb.rd $TOOLS/antigap.bed $STATS)"
echo -e "tss_1kb\t31069\txpd_summits\t$STATS\t$PVALUE" >> $RESULTS/tss_great_results.txt

STATS="$($PERLDIR/great_overlaps.pl $TOOLS/tss_100bp.rd $XPB_SUMMITS)"
PVALUE="$($GREAT $TOOLS/tss_100bp.rd $TOOLS/antigap.bed $STATS)"
echo -e "tss_100bp\t31069\txpb_summits\t$STATS\t$PVALUE" >> $RESULTS/tss_great_results.txt

STATS="$($PERLDIR/great_overlaps.pl $TOOLS/tss_100bp.rd $XPD_SUMMITS)"
PVALUE="$($GREAT $TOOLS/tss_100bp.rd $TOOLS/antigap.bed $STATS)"
echo -e "tss_100bp\t31069\txpd_summits\t$STATS\t$PVALUE" >> $RESULTS/tss_great_results.txt

# Calculate TSS binding and generate figures 3b, S3a, and S3b

$R $RDIR/figure3b_S1.R $TOOLS/tss_1kb.bed $XPB_SUMMITS $XPD_SUMMITS $ARRAYS/orb_gene_table.txt $RESULTS

# Get TSS locations for genes in each expression category

$R $RDIR/exp.cat.tss.R $TOOLS/tss_1kb.bed $RESULTS/gene_results_summary.txt $ARRAYS

# Calculate pileups for figure 3d
	
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpb.bw $ARRAYS/tss_1kb_cat1.bed 2000 $PILEUPS/xpb_tss_1kb_cat1.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpb.bw $ARRAYS/tss_1kb_cat2.bed 2000 $PILEUPS/xpb_tss_1kb_cat2.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpb.bw $ARRAYS/tss_1kb_cat3.bed 2000 $PILEUPS/xpb_tss_1kb_cat3.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpb.bw $ARRAYS/tss_1kb_cat4.bed 2000 $PILEUPS/xpb_tss_1kb_cat4.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpb.bw $ARRAYS/tss_1kb_cat5.bed 2000 $PILEUPS/xpb_tss_1kb_cat5.pile

$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpd.bw $ARRAYS/tss_1kb_cat1.bed 2000 $PILEUPS/xpd_tss_1kb_cat1.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpd.bw $ARRAYS/tss_1kb_cat2.bed 2000 $PILEUPS/xpd_tss_1kb_cat2.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpd.bw $ARRAYS/tss_1kb_cat3.bed 2000 $PILEUPS/xpd_tss_1kb_cat3.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpd.bw $ARRAYS/tss_1kb_cat4.bed 2000 $PILEUPS/xpd_tss_1kb_cat4.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpd.bw $ARRAYS/tss_1kb_cat5.bed 2000 $PILEUPS/xpd_tss_1kb_cat5.pile

$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpb_input.bw $ARRAYS/tss_1kb_cat1.bed 2000 $PILEUPS/xpb_input_tss_1kb_cat1.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpb_input.bw $ARRAYS/tss_1kb_cat2.bed 2000 $PILEUPS/xpb_input_tss_1kb_cat2.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpb_input.bw $ARRAYS/tss_1kb_cat3.bed 2000 $PILEUPS/xpb_input_tss_1kb_cat3.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpb_input.bw $ARRAYS/tss_1kb_cat4.bed 2000 $PILEUPS/xpb_input_tss_1kb_cat4.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpb_input.bw $ARRAYS/tss_1kb_cat5.bed 2000 $PILEUPS/xpb_input_tss_1kb_cat5.pile

$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpd_input.bw $ARRAYS/tss_1kb_cat1.bed 2000 $PILEUPS/xpd_input_tss_1kb_cat1.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpd_input.bw $ARRAYS/tss_1kb_cat2.bed 2000 $PILEUPS/xpd_input_tss_1kb_cat2.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpd_input.bw $ARRAYS/tss_1kb_cat3.bed 2000 $PILEUPS/xpd_input_tss_1kb_cat3.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpd_input.bw $ARRAYS/tss_1kb_cat4.bed 2000 $PILEUPS/xpd_input_tss_1kb_cat4.pile
$PERLDIR/bigwig_bed_pileup.pl $OVERLAPS/xpd_input.bw $ARRAYS/tss_1kb_cat5.bed 2000 $PILEUPS/xpd_input_tss_1kb_cat5.pile

# Generate graphs for figure 3c

$R $RDIR/tss.pileups.R $PILEUPS $RESULTS

# Calculate statistics for genes regulated in Halder, et al. (2012)
curl -L "http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3375199/bin/1756-0500-5-138-S2.XLSX" > $TREATMENTS/halder_tableS4.xlsx
$R $RDIR/convert.halder.R

$R $RDIR/treatment_comparisons.R $RESULTS/tss_binding_summary.txt $TREATMENTS/halder_360a $TREATMENTS/halder_phendc3 $TREATMENTS/halder_8979a $RESULTS

$R $RDIR/g4.tss.stats.R $TOOLS/tss_1kb.bed $RESULTS/xpb_peaks_g4.txt $RESULTS/xpd_peaks_g4.txt $RESULTS/xpb_regions_g4.txt $RESULTS/xpd_regions_g4.txt $XPB_SUMMITS $XPD_SUMMITS $RESULTS/g4-12_overlaps.txt $QP_DIR/cat_g4-12.bed $TOOLS/antigap.bed $RESULTS

#Calculate GC percents for peaks
curl -L "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/bigZips/hg19.2bit" > $TOOLS/hg19.2bit

$TBTOFA $TOOLS/hg19.2bit $XPB_PEAKS $PEAKS/XPB/xpb_peaks.fa
$TBTOFA $TOOLS/hg19.2bit $XPD_PEAKS $PEAKS/XPD/xpd_peaks.fa

$PERLDIR/fasta_gc_percent.pl $PEAKS/XPB/xpb-peaks.fa $PEAKS/XPB/xpb_peaks_gc.txt
$PERLDIR/fasta_gc_percent.pl $PEAKS/XPD/xpd-peaks.fa $PEAKS/XPD/xpd_peaks_gc.txt

#Calculate the GC percent of peaks
$R $RDIR/g4.tss.gc.stats.R $PEAKS/XPB/xpb_peaks_gc.txt $PEAKS/XPD/xpd_peaks_gc.txt $RESULTS/xpb_peaks_g4_tss.txt $RESULTS/xpd_peaks_g4_tss.txt $RESULTS

#Calculate the GC percent of randomly selected regions form comparison with peaks
mkdir $RESULTS/random_gc

echo -e "chr\tstart\tend\tname\tscore\tstrand\tg4\ttss" > header.txt
RAND_DIR=$RESULTS/random_gc

cp $XPB_PEAKS temp.bed

for i in {1..12}
do
	cat temp.bed $XPB_PEAKS > temp2.bed
	mv temp2.bed temp.bed
done
head -250000 temp.bed > $RAND_DIR/random_peaks.bed

$PERLDIR/bed_scramble.pl $RAND_DIR/random_peaks.bed $TOOLS/antigap.bed $TOOLS/chromWeight.txt temp.bed
mv temp.bed $RAND_DIR/random_peaks.bed

$PERLDIR/bed_sort.pl $RAND_DIR/random_peaks.bed temp.bed
mv temp.bed $RAND_DIR/random_peaks.bed

$PERLDIR/bed_center.pl $RAND_DIR/random_peaks.bed $RAND_DIR/random_summits.bed

$PERLDIR/convert_column.pl $RAND_DIR/random_peaks.bed temp.bed 5 0
mv temp.bed $RAND_DIR/random_peaks.bed

$TBTOFA $TOOLS/hg19.2bit $RAND_DIR/random_peaks.fa -bed=$RAND_DIR/random_peaks.bed
$PERLDIR/fasta_gc_percent.pl $RAND_DIR/random_peaks.fa $RAND_DIR/random_gc.txt

$PERLDIR/bed_overlap_binary.pl $RAND_DIR/random_summits.bed $G4 temp.txt
$PERLDIR/bed_overlap_binary.pl temp.txt $TOOLS/tss_1kb.bed temp2.txt
cat header.txt temp2.txt > $RAND_DIR/random_summits_g4_tss.txt

$R $RDIR/g4.tss.gc.stats.random.R $RAND_DIR/random_gc.txt $RAND_DIR/random_summits_g4_tss.txt $RESULTS

$R $RDIR/gc.match.stats.R $PEAKS/XPB/xpb_peaks_gc.txt $RESULTS/xpb_peaks_g4_tss.txt $PEAKS/XPD/xpd_peaks_gc.txt $RESULTS/xpd_peaks_g4_tss.txt $RESULTS/random_gc/random_gc.txt $RESULTS/random_gc/random_summits_g4_tss.txt $RESULTS 


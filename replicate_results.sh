# These variables define the locations of the directories for scripts and the various
# programs required to replicate this data analysis.
# Perl script directory
# BioPerl Bio::DB::BigWig module is required for some perl scripts
PERLDIR=~/perl
# R script directory
RDIR=~/Maizels/2012_reproducibility/r
# R binary location
R=R
# Bowtie binary location
BOWTIEDIR=~/chipseq/tools/bowtie-0.12.7
# IGVTools binary location
IGVTOOLS=~/chipseq/tools/IGVTools/igvtools
# MACS binary location
MACS=macs14
# CEAS binary location
CEAS=ceas
# UCSC Kent Source Tree is required for use of the bedGraphToBigWig script.
# UCSC bedGraphToBigWig binary location
BGTOBW=~/chipseq/tools/ucsc/bedGraphToBigWig
# MEME binary location
MEME=meme
# FIMO binary location
FIMO=fimo
# QuadParser binary location
QP=~/genomics/tools/quadparser/quadparser
# Great calculateBinomialP tool binary location
GREAT=~/chipseq/tools/greatTools/calculateBinomialP

# Original files
# Input: input_1.fastq and input_2.fastq
# XPB: xpb_chip_1.fastq and xpb_chip_2.fastq
# XPD: xpd_chip_1.fastq and xpd_chip_2.fastq

## Read quality control
# Trim to 60 reads
	$PERLDIR/fastq_filters/trim.pl input_1.fastq input_1.trim.fastq 60
	$PERLDIR/fastq_filters/trim.pl input_2.fastq input_2.trim.fastq 60
	$PERLDIR/fastq_filters/trim.pl xpb_chip_1.fastq xpb_chip_1.trim.fastq 60
	$PERLDIR/fastq_filters/trim.pl xpb_chip_2.fastq xpb_chip_2.trim.fastq 60
	$PERLDIR/fastq_filters/trim.pl xpd_chip_1.fastq xpd_chip_1.trim.fastq 60
	$PERLDIR/fastq_filters/trim.pl xpd_chip_2.fastq xpd_chip_2.trim.fastq 60

# Delete the raw files to save space
	rm input_1.fastq
	rm input_2.fastq
	rm xpb_chip_1.fastq
	rm xpb_chip_2.fastq
	rm xpd_chip_1.fastq
	rm xpd_chip_2.fastq

# Remove pairs with 3 or more reads with a PHRED score below 20
	$PERLDIR/fastq_filters/paired_min_qc_filter.pl input_1.trim.fastq input_1.filter.trim.fastq input_2.trim.fastq input_2.filter.trim.fastq 20 3
	$PERLDIR/fastq_filters/paired_min_qc_filter.pl xpb_chip_1.trim.fastq xpb_chip_1.filter.trim.fastq xpb_chip_2.trim.fastq xpb_chip_2.filter.trim.fastq 20 3
	$PERLDIR/fastq_filters/paired_min_qc_filter.pl xpd_chip_1.trim.fastq xpd_chip_1.filter.trim.fastq xpd_chip_2.trim.fastq xpd_chip_2.filter.trim.fastq 20 3

# Delete unfiltered files to save space
	rm input_1.trim.fastq
	rm input_2.trim.fastq
	rm xpb_chip_1.trim.fastq
	rm xpb_chip_2.trim.fastq
	rm xpd_chip_1.trim.fastq
	rm xpd_chip_2.trim.fastq

## Alignment of reads to the hg19 human genome using Bowtie v0.12.7
## Bowtie and the hg19 index were obtained from http://bowtie-bio.sourceforge.net/index.shtml
# Bowtie commands for SAM alignments (used to make bigWig files)
	$BOWTIEDIR/bowtie -t -p 3 -n 0 -m 1 --best -S --chunkmbs 256 hg19 -1 input_1.filter.trim.fastq -2 input_2.filter.trim.fastq input_hg19.sam
	$BOWTIEDIR/bowtie -t -p 3 -n 0 -m 1 --best -S --chunkmbs 256 hg19 -1 xpb_chip_1.filter.trim.fastq -2 xpb_chip_2.filter.trim.fastq xpb_hg19.sam
	$BOWTIEDIR/bowtie -t -p 3 -n 0 -m 1 --best -S --chunkmbs 256 hg19 -1 xpd_chip_1.filter.trim.fastq -2 xpd_chip_2.filter.trim.fastq xpd_hg19.sam

##Filter SAM files to remove unmapped reads or pairs with only 1 mapped mate
# hg19 alignments
	$PERLDIR/sam_remove_unmapped_pairs.pl input_hg19.sam input_hg19.mapped.sam
	$PERLDIR/sam_remove_unmapped_pairs.pl xpb_hg19.sam xpb_hg19.mapped.sam
	$PERLDIR/sam_remove_unmapped_pairs.pl xpd_hg19.sam xpd_hg19.mapped.sam

# Delete files with unmapped reads to save space
	rm input_hg19.sam
	rm xpb_hg19.sam
	rm xpd_hg19.sam

# Select Input reads to match the number of mapped XPB and XPD reads
# Since alignments don't sort the read positions, the reads will be random with respect to genomic position
# and I can simply use head to get an equal number.
	head -34399081 input_hg19.mapped.sam > xpb_input_hg19.mapped.sam
	head -43952749 input_hg19.mapped.sam > xpd_input_hg19.mapped.sam

# Sort SAM files using IGVTools
# IGVTools are available from http://www.broadinstitute.org/igv/igvtools
	$IGVTOOLS sort -m 2500000 xpb_input_hg19.mapped.sam xpb_input_hg19.sorted.mapped.sam
	$IGVTOOLS sort -m 2500000 xpd_input_hg19.mapped.sam xpd_input_hg19.sorted.mapped.sam
	$IGVTOOLS sort -m 2500000 xpb_hg19.mapped.sam xpb_hg19.sorted.mapped.sam
	$IGVTOOLS sort -m 2500000 xpd_hg19.mapped.sam xpd_hg19.sorted.mapped.sam

# Delete unsorted files to save space
	rm xpb_input_hg19.mapped.sam
	rm xpd_input_hg19.mapped.sam
	rm xpb_hg19.mapped.sam
	rm xpd_hg19.mapped.sam

# Covert sorted SAM files to BEDGraph files
	$PERLDIR/sam_sorted_to_bedgraph_faster.pl xpb_input_hg19.sorted.mapped.sam xpb_input.bedgraph 60 0
	$PERLDIR/sam_sorted_to_bedgraph_faster.pl xpd_input_hg19.sorted.mapped.sam xpd_input.bedgraph 60 0
	$PERLDIR/sam_sorted_to_bedgraph_faster.pl xpb_hg19.sorted.mapped.sam xpb.bedgraph 60 0
	$PERLDIR/sam_sorted_to_bedgraph_faster.pl xpd_hg19.sorted.mapped.sam xpd.bedgraph 60 0

# Convert BEDGraph files to bigWig using UCSC Tools
# UCSC command line tools can be obtained from http://hgdownload.cse.ucsc.edu/admin/exe/
	$BGTOBW xpb_input.bedgraph ~/chipseq/tools/ucsc/hg19.chrom.sizes xpb_input.bw
	$BGTOBW xpd_input.bedgraph ~/chipseq/tools/ucsc/hg19.chrom.sizes xpd_input.bw
	$BGTOBW xpb.bedgraph ~/chipseq/tools/ucsc/hg19.chrom.sizes xpb.bw
	$BGTOBW xpd.bedgraph ~/chipseq/tools/ucsc/hg19.chrom.sizes xpd.bw

## Peak finding with MACS v1.4.2 using the Bowtie map files, above
# MACS can be obtained from http://liulab.dfci.harvard.edu/MACS/

# First, must convert the SAM files to Bowtie format, because MACS has a problem recognizing
# which strand reads in SAM format are on and throws an error.
	$PERLDIR/sam_to_btm.pl xpb_input_hg19.sorted.mapped.sam xpb_input_hg19.btm
	$PERLDIR/sam_to_btm.pl xpd_input_hg19.sorted.mapped.sam xpd_input_hg19.btm
	$PERLDIR/sam_to_btm.pl xpb_hg19.sorted.mapped.sam xpb_hg19.btm
	$PERLDIR/sam_to_btm.pl xpd_hg19.sorted.mapped.sam xpd_hg19.btm

# MACS peak finding
	mkdir xpb_macs
	cd xpb_macs
	$MACS -f BOWTIE -s 60 -t ../xpb_hg19.btm -c ../xpb_input_hg19.btm -n xpb_macs
	cd ../
	mkdir xpd_macs
	cd xpd_macs
	$MACS -f BOWTIE -s 60 -t ../xpd_hg19.btm -c ../xpd_input_hg19.btm -n xpd_macs
	cd ../

# Filter only peaks with a p-value < 1e-15
	$PERLDIR/column_filter.pl xpb_macs/xpb_macs_peaks.bed xpb_1e-15_peaks.bed 5 gt 150
	$PERLDIR/column_filter.pl xpd_macs/xpd_macs_peaks.bed xpd_1e-15_peaks.bed 5 gt 150

# Sort peaks in chromosome and position order
	$PERLDIR/bed_sort.pl xpb_1e-15_peaks.bed temp.bed
	mv temp.bed xpb_1e-15_peaks.bed
	$PERLDIR/bed_sort.pl xpd_1e-15_peaks.bed temp.bed
	mv temp.bed xpd_1e-15_peaks.bed

# Correct MACS BED format to standard 6-column BED format and renumber peaks.
	$PERLDIR/macs_to_bed.pl xpb_1e-15_peaks.bed temp.bed xpb
	mv temp.bed xpb_1e-15_peaks.bed
	$PERLDIR/macs_to_bed.pl xpd_1e-15_peaks.bed temp.bed xpd
	mv temp.bed xpd_1e-15_peaks.bed

# Find peak summits using a Perl script that utilizes the Bio::DB::BigWig v1.01 BioPerl module
# Bio::DB::BigWig can be obtained at http://search.cpan.org/~lds/Bio-BigFile-1.01/lib/Bio/DB/BigWig.pm
	$PERLDIR/bigwig_bed_summits.pl xpb.bw xpb_1e-15_peaks.bed xpb_summits.txt
	$PERLDIR/bigwig_bed_summits.pl xpd.bw xpd_1e-15_peaks.bed xpd_summits.txt

# Make summit and summit region BED files
	$PERLDIR/summit_to_bed.pl xpb_summits.txt xpb_summits.bed 0 1
	$PERLDIR/summit_to_bed.pl xpd_summits.txt xpd_summits.bed 0 1
	$PERLDIR/summit_to_bed.pl xpb_summits.txt xpb_regions.bed 50 50
	$PERLDIR/summit_to_bed.pl xpd_summits.txt xpd_regions.bed 50 50

## Overlapping and non-overlapping XPB and XPD peaks
# XPB summit regions were used for "Both" summit regions.
	$PERLDIR/bed_overlap.pl xpb_1e-15_peaks.bed xpd_1e-15_peaks.bed both_peaks.bed
	$PERLDIR/bed_subtract.pl xpb_1e-15_peaks.bed xpd_1e-15_peaks.bed xpb_only_peaks.bed
	$PERLDIR/bed_subtract.pl xpd_1e-15_peaks.bed xpb_1e-15_peaks.bed xpd_only_peaks.bed
	$PERLDIR/bed_overlap.pl xpb_regions.bed both_peaks.bed both_regions.bed
	$PERLDIR/bed_overlap.pl xpb_regions.bed xpb_only_peaks.bed xpb_only_regions.bed
	$PERLDIR/bed_overlap.pl xpd_regions.bed xpd_only_peaks.bed xpd_only_regions.bed

## MEME searches for overrepresented motifs in XPB, XPD, and Both/Only region sets
# Regions were first filtered against the RepBase database to remove highly repetitive
# sequences that would yield false positives. The locations of RepBase elements were 
# retrieved from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
	curl -OL "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz"
	gzip -d rmsk.txt.gz
# The rmsk.txt file was first converted to BED format
	$PERLDIR/rmsk_to_bed.pl rmsk.txt rmsk.bed
# Then was used to filter all xpb and xpd region sets.
	$PERLDIR/bed_subtract.pl xpb_regions.bed rmsk.bed xpb_regions.rmsk.bed
	$PERLDIR/bed_subtract.pl xpd_regions.bed rmsk.bed xpd_regions.rmsk.bed
	$PERLDIR/bed_subtract.pl both_regions.bed rmsk.bed both_regions.rmsk.bed
	$PERLDIR/bed_subtract.pl xpb_only_regions.bed rmsk.bed xpb_only_regions.rmsk.bed
	$PERLDIR/bed_subtract.pl xpd_only_regions.bed rmsk.bed xpd_only_regions.rmsk.bed

# R was then used to retrieve summit height data for region files, sort these regions by
# summit height, and output the top 1000 regions for use with MEME.
	$R --vanilla $RDIR/sort.sumheight.1k.R

# The top peak files need to be sorted again based on position for retrieval:
	$PERLDIR/bed_sort.pl both_regions.top.rmsk.bed temp.bed
	mv temp.bed both_regions.top.rmsk.bed
	$PERLDIR/bed_sort.pl xpb_only_regions.top.rmsk.bed temp.bed
	mv temp.bed xpb_only_regions.top.rmsk.bed
	$PERLDIR/bed_sort.pl xpb_regions.top.rmsk.bed temp.bed
	mv temp.bed xpb_regions.top.rmsk.bed
	$PERLDIR/bed_sort.pl xpd_regions.top.rmsk.bed temp.bed
	mv temp.bed xpd_regions.top.rmsk.bed
	$PERLDIR/bed_sort.pl xpd_only_regions.top.rmsk.bed temp.bed
	mv temp.bed xpd_only_regions.top.rmsk.bed

# Download the hg19 genome in FASTA format:
	mkdir hg19_chromFa
	cd hg19_chromFa
	curl -L "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz" > hg19_chromFa.tar.gz
	tar xvfz hg19_chromFa.tar.gz
	cd ../

# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
# This resulted in one fasta file for each chromosome in a directory named hg19_chromFa.
# These files had to be converted to all caps using this shell script:
	
	mkdir hg19_chromFa_allcaps
	
	for f in hg19_chromFa/*.fa
	do
		NOPATH=${f##*/}
		FILENAME=${NOPATH%.[^.]*}
		head -1 $f > hg19_chromFa_allcaps/$NOPATH
		tail -n +2 $f | tr '[a-z]' '[A-Z]' >> hg19_chromFa_allcaps/$NOPATH
	done

# The FASTA sequences for regions defined in:
# xpb_regions.bed, xpd_regions.bed, xpb_regions.top.rmsk.bed, xpd_regions.top.rmsk.bed,  
# both_regions.top.rmsk.bed, xpb_only_regions.top.rmsk.bed, and xpd_only_regions.top.rmsk.bed
# need to be retrieved from the FASTA-formatted genome files.
	$PERLDIR/bed_retrieve_fasta.pl xpb_regions.top.rmsk.bed hg19_chromFa_allcaps xpb_regions.top.rmsk.fasta
	$PERLDIR/bed_retrieve_fasta.pl xpd_regions.top.rmsk.bed hg19_chromFa_allcaps xpd_regions.top.rmsk.fasta
	$PERLDIR/bed_retrieve_fasta.pl both_regions.top.rmsk.bed hg19_chromFa_allcaps both_regions.top.rmsk.fasta
	$PERLDIR/bed_retrieve_fasta.pl xpb_only_regions.top.rmsk.bed hg19_chromFa_allcaps xpb_only_regions.top.rmsk.fasta
	$PERLDIR/bed_retrieve_fasta.pl xpd_only_regions.top.rmsk.bed hg19_chromFa_allcaps xpd_only_regions.top.rmsk.fasta

# MEME was used to search these FASTA files. 
# MEME can be obtained from http://meme.nbcr.net/meme/meme-download.html
	$MEME xpb_regions.top.rmsk.fasta -mod zoops -nmotifs 5 -minw 6 -maxw 15 -revcomp -maxsize 100000 -p 3 -dna -oc meme_xpb_regions
	$MEME xpd_regions.top.rmsk.fasta -mod zoops -nmotifs 5 -minw 6 -maxw 15 -revcomp -maxsize 100000 -p 3 -dna -oc meme_xpd_regions
	$MEME both_regions.top.rmsk.fasta -mod zoops -nmotifs 5 -minw 6 -maxw 15 -revcomp -maxsize 100000 -p 3 -dna -oc meme_both_regions
	$MEME xpb_only_regions.top.rmsk.fasta -mod zoops -nmotifs 5 -minw 6 -maxw 15 -revcomp -maxsize 100000 -p 3 -dna -oc meme_xpb_only_regions
	$MEME xpd_only_regions.top.rmsk.fasta -mod zoops -nmotifs 5 -minw 6 -maxw 15 -revcomp -maxsize 100000 -p 3 -dna -oc meme_xpd_only_regions

# MEME-formatted motif files were generated for motifs of interest
	$PERLDIR/makemotif.pl TGACTCA AP1 ap1.motif
	$PERLDIR/makemotif.pl RYTTCCKG ETS ets.motif
	$PERLDIR/makemotif.pl AGRRARRR MAZ maz.motif

# The unfiltered summit region FASTA sequences were obtained for xpb and xpd:
	$PERLDIR/bed_retrieve_fasta.pl xpb_regions.bed hg19_chromFa_allcaps xpb_regions.fasta
	$PERLDIR/bed_retrieve_fasta.pl xpd_regions.bed hg19_chromFa_allcaps xpd_regions.fasta
	$PERLDIR/bed_retrieve_fasta.pl xpb_1e-15_peaks.bed hg19_chromFa_allcaps xpb_peaks.fasta
	$PERLDIR/bed_retrieve_fasta.pl xpd_1e-15_peaks.bed hg19_chromFa_allcaps xpd_peaks.fasta	

# All XPB and XPD summit regions were searched for these motifs using FIMO
	$FIMO --oc fimo_xpb_regions_ap1 ap1.motif xpb_regions.fasta
	$FIMO --oc fimo_xpb_regions_ets ets.motif xpb_regions.fasta
	$FIMO --oc fimo_xpb_regions_maz maz.motif xpb_regions.fasta
	$FIMO --oc fimo_xpd_regions_ap1 ap1.motif xpd_regions.fasta
	$FIMO --oc fimo_xpd_regions_ets ets.motif xpd_regions.fasta
	$FIMO --oc fimo_xpd_regions_maz maz.motif xpd_regions.fasta

	$FIMO --oc fimo_xpb_peaks_ap1 ap1.motif xpb_peaks.fasta
	$FIMO --oc fimo_xpb_peaks_ets ets.motif xpb_peaks.fasta
	$FIMO --oc fimo_xpb_peaks_maz maz.motif xpb_peaks.fasta
	$FIMO --oc fimo_xpd_peaks_ap1 ap1.motif xpd_peaks.fasta
	$FIMO --oc fimo_xpd_peaks_ets ets.motif xpd_peaks.fasta
	$FIMO --oc fimo_xpd_peaks_maz maz.motif xpd_peaks.fasta

# Scripts were used to generate results tables for these motif searches
	$PERLDIR/fimo_results_table_new.pl fimo_xpb_regions_ap1/fimo.txt xpb_regions.bed ap1_xpb_regions.txt
	$PERLDIR/fimo_results_table_new.pl fimo_xpb_regions_ets/fimo.txt xpb_regions.bed ets_xpb_regions.txt
	$PERLDIR/fimo_results_table_new.pl fimo_xpb_regions_maz/fimo.txt xpb_regions.bed maz_xpb_regions.txt
	$PERLDIR/fimo_results_table_new.pl fimo_xpd_regions_ap1/fimo.txt xpd_regions.bed ap1_xpd_regions.txt
	$PERLDIR/fimo_results_table_new.pl fimo_xpd_regions_ets/fimo.txt xpd_regions.bed ets_xpd_regions.txt
	$PERLDIR/fimo_results_table_new.pl fimo_xpd_regions_maz/fimo.txt xpd_regions.bed maz_xpd_regions.txt
	
	$PERLDIR/fimo_results_table_new.pl fimo_xpb_peaks_ap1/fimo.txt xpb_1e-15_peaks.bed ap1_xpb_peaks.txt
	$PERLDIR/fimo_results_table_new.pl fimo_xpb_peaks_ets/fimo.txt xpb_1e-15_peaks.bed ets_xpb_peaks.txt
	$PERLDIR/fimo_results_table_new.pl fimo_xpb_peaks_maz/fimo.txt xpb_1e-15_peaks.bed maz_xpb_peaks.txt
	$PERLDIR/fimo_results_table_new.pl fimo_xpd_peaks_ap1/fimo.txt xpd_1e-15_peaks.bed ap1_xpd_peaks.txt
	$PERLDIR/fimo_results_table_new.pl fimo_xpd_peaks_ets/fimo.txt xpd_1e-15_peaks.bed ets_xpd_peaks.txt
	$PERLDIR/fimo_results_table_new.pl fimo_xpd_peaks_maz/fimo.txt xpd_1e-15_peaks.bed maz_xpd_peaks.txt


## CEAS analysis of XPB and XPD peaks
# CEAS was obtained from http://liulab.dfci.harvard.edu/CEAS/
# CEAS requires a refGene table
	mkdir ceas
	cd ceas
	curl -L "http://liulab.dfci.harvard.edu/CEAS/src/hg19.refGene.gz" > hg19.refGene.gz
	gzip -d hg19.refGene.gz
	cd ceas
	$CEAS -b ../xpb_1e-15_peaks.bed -g hg19.refGene --name=xpb
	$CEAS -b ../xpd_1e-15_peaks.bed -g hg19.refGene --name=xpd

# Generate the CEAS report PDFs using R
	$R --vanilla < xpb.R
	$R --vanilla < xpd.R
	cd ../

## Quadparser searches in the hg19 genome
# Quadparser was downloaded from http://www-shankar.ch.cam.ac.uk/Applications/quadparser/
# The following shell script was then used to search each chromosome FASTA for G4 motifs,
# then compile and sort the results into a single BED-formatted file.

	mkdir hg19_qp
	
	for f in hg19_chromFa_allcaps/*
	do
		NOPATH=${f##*/}
		FILENAME=${NOPATH%.[^.]*}
		$QP -DAS $f GC 3 4 1 12 hg19_qp/$FILENAME.bed
		$PERLDIR/remove_lines.pl hg19_qp/$FILENAME.bed hg19_qp/temp.bed 1 track\ name=Quadruplexes
		$PERLDIR/convert_column.pl hg19_qp/temp.bed hg19_qp/$FILENAME.bed 1 $FILENAME
	done
	rm hg19_qp/temp.bed
	
	cat hg19_qp/*.bed > hg19_qp/temp.bed
	$PERLDIR/bed_sort.pl hg19_qp/temp.bed g4-12.bed
	rm hg19_qp/temp.bed
	
# Then, I used the GREAT calculateBinomialP tool to calculate significance of overlaps between
# G4 and XPB/XPD peaks.

	$PERLDIR/bed_to_great.pl xpb_1e-15_peaks.bed xpb_peaks.rd
	$PERLDIR/bed_to_great.pl xpd_1e-15_peaks.bed xpd_peaks.rd
	$PERLDIR/bed_to_great.pl xpb_regions.bed xpb_regions.rd
	$PERLDIR/bed_to_great.pl xpd_regions.bed xpd_regions.rd
	
	echo -e "target_name\ttarget_count\tprobe_name\tprobe_count\toverlap_count\tp-value" > g4_great_results.txt
	STATS="$($PERLDIR/great_overlaps.pl xpb_peaks.rd g4-12.bed)"
	PVALUE="$($GREAT xpb_peaks.rd antigap.bed $STATS)"
	echo -e "xpb_peaks\t26051\tg4-12\t$STATS\t$PVALUE" >> g4_great_results.txt

	STATS="$($PERLDIR/great_overlaps.pl xpd_peaks.rd g4-12.bed)"
	PVALUE="$($GREAT xpd_peaks.rd antigap.bed $STATS)"
	echo -e "xpd_peaks\t16730\tg4-12\t$STATS\t$PVALUE" >> g4_great_results.txt

	STATS="$($PERLDIR/great_overlaps.pl xpb_regions.rd g4-12.bed)"
	PVALUE="$($GREAT xpb_regions.rd antigap.bed $STATS)"
	echo -e "xpb_regions\t26051\tg4-12\t$STATS\t$PVALUE" >> g4_great_results.txt

	STATS="$($PERLDIR/great_overlaps.pl xpd_regions.rd g4-12.bed)"
	PVALUE="$($GREAT xpd_regions.rd antigap.bed $STATS)"
	echo -e "xpd_regions\t16730\tg4-12\t$STATS\t$PVALUE" >> g4_great_results.txt


# Calculate G4 overlaps with TSS regions, XPB and XPD peaks, and XPB and XPD regions
	
	$PERLDIR/bed_overlap_binary.pl g4-12.bed tss_1kb.bed temp.txt
	mv temp.txt g4-12_overlaps.txt
	$PERLDIR/bed_overlap_binary.pl g4-12_overlaps.txt xpb_1e-15_peaks.bed temp.txt
	mv temp.txt g4-12_overlaps.txt
	$PERLDIR/bed_overlap_binary.pl g4-12_overlaps.txt xpd_1e-15_peaks.bed temp.txt
	mv temp.txt g4-12_overlaps.txt
	$PERLDIR/bed_overlap_binary.pl g4-12_overlaps.txt xpb_regions.bed temp.txt
	mv temp.txt g4-12_overlaps.txt
	$PERLDIR/bed_overlap_binary.pl g4-12_overlaps.txt xpd_regions.bed temp.txt
	head temp.txt
	echo -e "chr\tstart\tend\tname\tscore\tstrand\ttss\txpb.peaks\txpd.peaks\txpb.regions\txpd.regions" > temp2.txt
	cat temp2.txt temp.txt > g4-12_overlaps.txt

	$PERLDIR/bed_overlap_binary.pl xpb_1e-15_peaks.bed g4-12.bed temp.txt
	echo -e "chr\tstart\tend\tname\tscore\tstrand\tg4" > temp2.txt
	cat temp2.txt temp.txt > xpb_peaks_g4.txt
	$PERLDIR/bed_overlap_binary.pl xpd_1e-15_peaks.bed g4-12.bed temp.txt
	echo -e "chr\tstart\tend\tname\tscore\tstrand\tg4" > temp2.txt
	cat temp2.txt temp.txt > xpd_peaks_g4.txt
	$PERLDIR/bed_overlap_binary.pl xpb_regions.bed g4-12.bed temp.txt
	echo -e "chr\tstart\tend\tname\tscore\tstrand\tg4" > temp2.txt
	cat temp2.txt temp.txt > xpb_regions_g4.txt
	$PERLDIR/bed_overlap_binary.pl xpd_regions.bed g4-12.bed temp.txt
	echo -e "chr\tstart\tend\tname\tscore\tstrand\tg4" > temp2.txt
	cat temp2.txt temp.txt > xpd_regions_g4.txt
	
	rm temp.txt
	rm temp2.txt
	
# Calculate the frequency of TSS and G4 overlaps for XPB and XPD, as well as coverage calculations.
	$R --vanilla < $RDIR/g4.tss.stats.R

# I generated tables of the locations of G4 motifs in XPB and XPD summit regions:
	$PERLDIR/g4_results_table.pl g4-12.bed xpb_regions.bed g4_xpb_regions.txt
	$PERLDIR/g4_results_table.pl g4-12.bed xpd_regions.bed g4_xpd_regions.txt
	$PERLDIR/g4_results_table.pl g4-12.bed xpb_1e-15_peaks.bed g4_xpb_peaks.txt
	$PERLDIR/g4_results_table.pl g4-12.bed xpd_1e-15_peaks.bed g4_xpd_peaks.txt
	
# XPB and XPD summit regions containing G4, AP1, ETS1, and MAZ motifs were selected using R for use with the online version of GREAT
	$R --vanilla < $RDIR/motif.peaks.R

# Calculate the frequency and coincidence of G4, AP1, ETS1, and MAZ motifs in whole peaks.

	$R --vanilla < $RDIR/motif.counts.R

# Generate Venn diagrams for these motif counts

	$R --vanilla < $RDIR/venn.diagrams.R

## RefSeq TSS region selection
# The RefSeq gene table was obtained from the UCSC genome browser:
	curl -L "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz" > refGene.txt.gz
	gzip -d refGene.txt.gz

# Get TSS +/- 1kb and TSS +/- 100bp

	$PERLDIR/refseq_tss_region.pl refGene.txt tss_1kb.bed 1000
	$PERLDIR/refseq_tss_region.pl refGene.txt tss_100bp.bed 100
	
# Use R to get just unique lines of the TSS files:
	
	$R --vanilla < $RDIR/unique.tss.R

# Tabulate percent of peaks in TSS +/- 1kb regions, as well as peaks in RefSeq gene bodies:

	$PERLDIR/refseq_tss_region.pl refGene.txt tss_1kb_all.bed 1000
	$R --vanilla < $RDIR/peak.tss.gene.R
	
# Calculate significance of enrichment of XPB and XPD summits near TSS:
	
	$PERLDIR/bed_to_great.pl tss_1kb.bed tss_1kb.rd
	$PERLDIR/bed_to_great.pl tss_100bp.bed tss_100bp.rd
	
	echo -e "target_name\ttarget_count\tprobe_name\tprobe_count\toverlap_count\tp-value" > tss_great_results.txt
	STATS="$($PERLDIR/great_overlaps.pl tss_1kb.rd xpb_summits.bed)"
	PVALUE="$($GREAT tss_1kb.rd antigap.bed $STATS)"
	echo -e "tss_1kb\t31069\txpb_summits\t$STATS\t$PVALUE" >> tss_great_results.txt

	STATS="$($PERLDIR/great_overlaps.pl tss_1kb.rd xpd_summits.bed)"
	PVALUE="$($GREAT tss_1kb.rd antigap.bed $STATS)"
	echo -e "tss_1kb\t31069\txpd_summits\t$STATS\t$PVALUE" >> tss_great_results.txt

	STATS="$($PERLDIR/great_overlaps.pl tss_100bp.rd xpb_summits.bed)"
	PVALUE="$($GREAT tss_100bp.rd antigap.bed $STATS)"
	echo -e "tss_100bp\t31069\txpb_summits\t$STATS\t$PVALUE" >> tss_great_results.txt

	STATS="$($PERLDIR/great_overlaps.pl tss_100bp.rd xpd_summits.bed)"
	PVALUE="$($GREAT tss_100bp.rd antigap.bed $STATS)"
	echo -e "tss_100bp\t31069\txpd_summits\t$STATS\t$PVALUE" >> tss_great_results.txt

# Calculate TSS binding and generate figures 3a and 3b

	$R --vanilla < $RDIR/figure3ab.R

# Get TSS locations for genes in each expression category

	$R --vanilla < $RDIR/exp.cat.tss.R

# Calculate pileups for figure 3d
	
	mkdir xpb_pileups
	cd xpb_pileups
	$PERLDIR/bigwig_bed_pileup.pl ../xpb.bw ../tss_1kb_cat1.bed 2000 xpb_tss_1kb_cat1.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpb.bw ../tss_1kb_cat2.bed 2000 xpb_tss_1kb_cat2.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpb.bw ../tss_1kb_cat3.bed 2000 xpb_tss_1kb_cat3.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpb.bw ../tss_1kb_cat4.bed 2000 xpb_tss_1kb_cat4.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpb.bw ../tss_1kb_cat5.bed 2000 xpb_tss_1kb_cat5.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpb.bw ../tss_1kb_cat6.bed 2000 xpb_tss_1kb_cat6.pile
	cd ../

	mkdir xpd_pileups
	cd xpd_pileups
	$PERLDIR/bigwig_bed_pileup.pl ../xpd.bw ../tss_1kb_cat1.bed 2000 xpd_tss_1kb_cat1.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpd.bw ../tss_1kb_cat2.bed 2000 xpd_tss_1kb_cat2.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpd.bw ../tss_1kb_cat3.bed 2000 xpd_tss_1kb_cat3.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpd.bw ../tss_1kb_cat4.bed 2000 xpd_tss_1kb_cat4.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpd.bw ../tss_1kb_cat5.bed 2000 xpd_tss_1kb_cat5.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpd.bw ../tss_1kb_cat6.bed 2000 xpd_tss_1kb_cat6.pile
	cd ../
	
	mkdir xpb_input_pileups
	cd xpb_input_pileups
	$PERLDIR/bigwig_bed_pileup.pl ../xpb_input.bw ../tss_1kb_cat1.bed 2000 xpb_input_tss_1kb_cat1.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpb_input.bw ../tss_1kb_cat2.bed 2000 xpb_input_tss_1kb_cat2.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpb_input.bw ../tss_1kb_cat3.bed 2000 xpb_input_tss_1kb_cat3.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpb_input.bw ../tss_1kb_cat4.bed 2000 xpb_input_tss_1kb_cat4.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpb_input.bw ../tss_1kb_cat5.bed 2000 xpb_input_tss_1kb_cat5.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpb_input.bw ../tss_1kb_cat6.bed 2000 xpb_input_tss_1kb_cat6.pile
	cd ../
	
	mkdir xpd_input_pileups
	cd xpd_input_pileups
	$PERLDIR/bigwig_bed_pileup.pl ../xpd_input.bw ../tss_1kb_cat1.bed 2000 xpd_input_tss_1kb_cat1.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpd_input.bw ../tss_1kb_cat2.bed 2000 xpd_input_tss_1kb_cat2.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpd_input.bw ../tss_1kb_cat3.bed 2000 xpd_input_tss_1kb_cat3.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpd_input.bw ../tss_1kb_cat4.bed 2000 xpd_input_tss_1kb_cat4.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpd_input.bw ../tss_1kb_cat5.bed 2000 xpd_input_tss_1kb_cat5.pile
	$PERLDIR/bigwig_bed_pileup.pl ../xpd_input.bw ../tss_1kb_cat6.bed 2000 xpd_input_tss_1kb_cat6.pile	
	cd ../
		
# Generate graphs for figure 3d
	
	$R --vanilla < $RDIR/tss.pileups.R
	
# Use R to calculate and output results for Figure 4 (processivity of genes bound by XPB and XPD)
	
	$R --vanilla < $RDIR/figure4.R
	
	
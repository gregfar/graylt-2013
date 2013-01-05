#!/usr/bin/perl
# refseq_tss_regions.pl - generates a BED file with locations around each TSS in the UCSC RefSeq Genes table.
# usage: refseq_tss_regions.pl in_table out_bed distance

if(@ARGV != 3) {

	print "usage: refseq_tss_regions.pl in_table out_file distance\n";

} else {
	
	#use arguments to define which files to read
	$in_file = @ARGV[0];
	$out_file = @ARGV[1];
	$distance = @ARGV[2];
	
	#set which exon and intron to retrieve for the subroutines.
	$distance = 1000;
	$exon = 1;
	$intron = 1;

	#subroutines for retrieving TSS regions, first exon, and first intron:
	sub get_tss_region {
	
		if(@line_split[3] eq "+") {
			
			$start = @exon_starts[0] - $distance;
			$end = $start + (2 * $distance);
			
			print OUTPUT $start . "\t" . $end;
			
		} elsif (@line_split[3] eq "-") {
			
			$start = @exon_ends[-1] - $distance;
			$end = $start + (2 * $distance);
			
			print OUTPUT $start . "\t" . $end;
			
		}
			
	}

	#open the input and output files for reading and writing
	open(INPUT, "<$in_file");
	open(OUTPUT, ">$out_file");
	
	while(<INPUT>) {
		
		#for each line of the input, remove the line break, then split the line by tabs
		
		$line = $_;
		chomp($line);
			
		@line_split = split(/\t/,$line);
		
		#split up the exon start and end locations for retrieval by the subroutines.
		@exon_starts = split(/,/,@line_split[9]);
		@exon_ends = split(/,/,@line_split[10]);
		
		#print the chromosome to the beginning of a new line
		print OUTPUT @line_split[2] . "\t";
		
		#then get a region around the TSS defined by the $distance variable at the beginning of the script.
		&get_tss_region;
		
		#next, print out the gene name, a score (0 as a placeholder), and the strand. This makes the TSS
		#region look like a BED formatted file, which will retain compatibility with my other scripts.
		print OUTPUT "\t" . @line_split[12] . "\t0\t" . @line_split[3] . "\t";	
		
	}
	
	close(INPUT);
	close(OUTPUT);

}
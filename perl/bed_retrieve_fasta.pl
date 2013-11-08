#!/usr/bin/perl
#bed_retrieve_fasta.pl - Reads a folder of FASTA files and retrieves the
#corresponding sequences from a BED file. 
#CAVEATS:
#The BED file must be sorted first.
#BED regions cannot overlap.
#This script will also have problems if the BED regions are very close together.

if(@ARGV != 4) {
	print "usage: bed_retrieve_fasta.pl in_bed in_fasta_dir out_fasta out_scores\n";
} else {
	
	$in_bed = @ARGV[0];
	$in_fasta_dir =  @ARGV[1];
	$out_fasta = @ARGV[2];
	$out_scores = @ARGV[3];
	
	# Read the BED file, split the file into separate arrays for each chromosome
	# and make an array containing all of the chromosome names.
	
	open(INPUT, "<$in_bed");
	
	$cur_chr = "nomatch";
	
	while(<INPUT>) {
		
		$line = $_;
		chomp($line);
		
		@bed_split = split(/\t/, $line);
		
		$chr = @bed_split[0];
		push(@{"$chr"},$line);
		
		if($chr ne $cur_chr) {
			push(@chr_array,$chr);
			$cur_chr = $chr;
		}
		
	}
	
	# for each chromosome, read the FASTA file, and if the position falls within the
	# region specified, write the sequence to the output file.
	
	for(@chr_array) {
		$chr = $_;
		$chr_array_len = @{"$chr"};
		
		$bed_pos = 0;
		$bed_line = @{"$chr"}[$bed_pos];
		@bed_split = split(/\t/,$bed_line);
		$bed_start = @bed_split[1];
		$bed_end = @bed_split[2];
		$out_line = @bed_split[3] . "\t" . @bed_split[5] . "\t";
		
		$fasta = $in_fasta_dir . "/" . $chr . ".fa";
				
		open(FASTA, "<$fasta");
		
		$fasta_line = 0;
		$fasta_len = 0;
		$fasta_pos = 0;
		
		while(<FASTA>) {
			
			$line = $_;
			chomp($line);
			
			@line_split = split(//, $line);
			
			# find out how many characters are on each line
			if($fasta_line == 1) {
				$fasta_len = length($line);
				$fasta_start = 0;
				$fasta_end = $fasta_len - 1;
			}
			
			# if this isn't the header, check to see if the fasta position falls
			# within a bed target site.
			
			if($fasta_line > 0) {
				
				if($bed_end < $fasta_start) {
					# if the fasta region is further along the chromosome than the bed region,
					# advance the bed file.
					$bed_pos++;
					push(@out_array,$out_line);
					if($bed_pos > $chr_array_len) {
						# if there aren't any more bed lines for this chromosome, stop
						# looking through the fasta file.
						last;
					} else {
						# otherwise, get the info for the new bed region.
						$bed_line = @{"$chr"}[$bed_pos];
						@bed_split = split(/\t/,$bed_line);
						$bed_start = @bed_split[1];
						$bed_end = @bed_split[2];
						$out_line = @bed_split[3] . "\t" . @bed_split[5] . "\t";
					}
				} elsif (($bed_start <= $fasta_start) & ($bed_end > $fasta_end)) {
					# if the fasta line is from the middle of the bed region, add
					# the whole fasta line to the output line.
					$out_line = $out_line . $line;
#					print "W" . "\n";
				} elsif (($bed_start >= $fasta_start) & ($bed_end <= $fasta_end)) {
					# if the fasta line contains the whole bed region, output just the
					# part of the line that matches the region.
					$sub_start = $bed_start - $fasta_start;
					$sub_length = $bed_end - $bed_start;
#					print "A" . "\n";
					$sub_fa = substr($line,$sub_start,$sub_length);
					
					$out_line = $out_line . $sub_fa;
				} elsif (($bed_start >= $fasta_start) & ($bed_start <= $fasta_end) & ($bed_end > $fasta_end)) {
					# if the fasta line has only the beginning of the bed region,
					# get just the relevant part of the line.
					$sub_start = $bed_start - $fasta_start;
					$sub_length = $fasta_len - $sub_start;
#					print "B" . "\n";
#					print $sub_start . "\t" . $sub_length . "\n";
					$sub_fa = substr($line,$sub_start,$sub_length);
					
					$out_line = $out_line . $sub_fa;
				} elsif (($bed_start < $fasta_start) & ($bed_end >= $fasta_start) & ($bed_end <= $fasta_end)) {
					# if the fasta line has only the end of the bed region,
					# get just the relevant part of the line.
					$sub_start = 0;
					$sub_length = $bed_end - $fasta_start;
#					print "C" . "\n";
					
					$sub_fa = substr($line,$sub_start,$sub_length);
					
					$out_line = $out_line . $sub_fa;
				}
				
				$fasta_start = $fasta_end + 1;
				$fasta_end = $fasta_end + $fasta_len;
			}
			
			$fasta_line++;
			
		}


	}
	
	open(OUTPUT, ">$out_fasta");
	open(OUT_SCORES, ">$out_scores");

	foreach(@out_array) {
		$out_line = $_;
		($out_name,$out_strand,$out_seq) = split(/\t/,$out_line);
		$l = length($out_seq);
		

		if($l > 0) {

			$c = grep /[GgCc]/, split(//,$out_seq);
			$out_percent = $c / $l;

			$n = grep /[Nn]/, split(//,$out_seq);

			if($n > 0) {
				print OUT_SCORES $out_name . "\tN\n"
			} else {
				print OUT_SCORES $out_name . "\t" . $out_percent . "\n";
			}

			if($out_strand eq "+") {
				print OUTPUT ">" . $out_name . "\n" . $out_seq . "\n";
			} elsif ($out_strand eq "-") {
				
				$revcomp = reverse($out_seq);
				$revcomp =~ tr/ACGTacgt/TGCAtgca/;
				print OUTPUT ">" . $out_name . "\n" . $revcomp . "\n";
				
			}

		}
		
	}
	
}
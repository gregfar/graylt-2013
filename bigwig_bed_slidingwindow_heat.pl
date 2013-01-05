#!/usr/bin/perl
# bigwig_bed_gc_slidingwindow_chr.pl - Uses the BioPerl Bio-BigFile scripts to directly access 
# bigWig scores in each region in a BED file. This will calculate the average value in a 
# sliding window centered on each base in the BED file. Results will be
# converted to a bedgraph-like format for graphing in R using the
# ggplot2 geom_rect tool. This format consists of an ID (from the name column of the bed
# file), the relative start and end positions, the score for the region from the bigwig,
# and a rank based on the order of the bed file.
# 
# This is designed to work with my GC bw files that are separate for each chromosome.
#
# usage: bigwig_bed_gc_slidingwindow_chr.pl in_bigwig in_bed out_heat window_size

if(@ARGV != 4) {
	print "usage: bigwig_bed_heatmap.pl in_bigwig_dir in_bed out_heat window_size\n";
} else {
	
	#load the BioPerl Bio-BigFile BigWig methods
	use Bio::DB::BigWig 'binMean','binStdev';
	
	$in_bw_dir = @ARGV[0];
	$in_bed = @ARGV[1];
	$out_file = @ARGV[2];
	$window_size = @ARGV[3];
	$half_window = int($window_size/2);
	
	open(INPUT, "<$in_bed");
	open(OUTPUT, ">$out_file");
	
	# put a header in the output for easier importing into R
	print OUTPUT "id\tstart\tend\tval\trank\n";
	
	$rank = 0;
	
	while(<INPUT>) {
		
		$line = $_;
		chomp($line);
		@bed_split = split(/\t/,$line);
		
		$chr = @bed_split[0];
		$start = @bed_split[1];
		$end = @bed_split[2];
		$name = @bed_split[3];
		$length = $end - $start;
		$strand = @bed_split[5];
		
		#DB::BigWig has a bug that makes it so that the bins retrieved are each an extra
		#base long. So, I'll have to break down each region into $bin_size-length regions
		#myself, then read each one out of the bigwig file.
		
		# tell the BioPerl methods where the bigwig file is
		$bw_file = $in_bw_dir . "/" . $chr . ".bw";
		my $bw = Bio::DB::BigWig->new(-bigwig=>$bw_file);
		
		$out_mean = 0;

		for($i = 0; $i < $length; $i++) {
			
			$window_start = ($start - $half_window) + $i;
			$window_end = $window_start + $window_size;

			#retrieve the values in the bigWig file
			my @features = $bw->features(-seq_id=>$chr,-start=>$window_start,-end=>$window_end);
		
			#calculate the mean for this position in the region
			$total = 0;
			
			for my $f (@features) {
				$f_start = $f->start;
				$f_end = $f->end;
				$score = $f->score;
				
				$total += $score * (($f_end - $f_start) + 1);
			}
			
			$mean = $total/$window_size;
			
			#adjust the position based on the strand
			if( $strand eq '+' ) {
				$rel_start = $i;
				$rel_end = $rel_start + 1;
			} elsif ( $strand eq '-' ) {
				$rel_start = ($length - $i) - 1;
				$rel_end = $rel_start + 1;
			}
			
			# Then, check to see if the mean is different from the last position.
			
			if($mean != $out_mean) {
			
				if($i > 0) {
					#If the mean is different, output the region with the last mean value.
					print OUTPUT $name . "\t" . $out_start . "\t" . $out_end . "\t" . $out_mean . "\t" . $rank . "\n";
				}
				
				#Set up a new region with the current mean value.
				$out_mean = $mean;
				$out_start = $rel_start;
				$out_end = $rel_end;
			
			} elsif ($mean == $out_mean) {
				#If the mean is the same, just extend the output region
				if( $strand eq '+' ) {
					$out_end = $rel_end;
				} elsif( $strand eq '-' ) {
					$out_start = $rel_start;
				}
			}

		}
		
		#Must do one last output for the end of the region.
		print OUTPUT $name . "\t" . $out_start . "\t" . $out_end . "\t" . $out_mean . "\t" . $rank . "\n";

		#increment the rank
		$rank++;
		
	}
	
	close(OUTPUT);
	close(INPUT);
	
}
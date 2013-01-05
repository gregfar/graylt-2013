#!/usr/bin/perl
# bigwig_bed_gc_slidingwindow_chr_pileup.pl - Uses the BioPerl Bio-BigFile scripts to directly access 
# bigWig scores in each region in a BED file. This will calculate the average value in a 
# sliding window centered on each base in the BED file. An average is calculated for each base
# in the region size provided.
# 
# This is designed to work with my GC bw files that are separate for each chromosome.
#
# usage: bigwig_bed_gc_slidingwindow_chr_pileup.pl in_bigwig in_bed out_heat window_size

if(@ARGV != 5) {
	print "usage: bigwig_bed_heatmap.pl in_bigwig_dir in_bed out_heat window_size positions\n";
} else {
	
	#load the BioPerl Bio-BigFile BigWig methods
	use Bio::DB::BigWig 'binMean','binStdev';
	
	$in_bw_dir = @ARGV[0];
	$in_bed = @ARGV[1];
	$out_file = @ARGV[2];
	$window_size = @ARGV[3];
	$positions = @ARGV[4];
	$length = $positions - 1;
	$half_window = int($window_size/2);

	#Initiate the pileup map
	@pileup = ();
	
	foreach(0..$length) {
		push(@pileup, 0);
	}
	
	open(INPUT, "<$in_bed");
	
	$bed_count = 0;
	while(<INPUT>) {
		
		$line = $_;
		chomp($line);
		@bed_split = split(/\t/,$line);
		
		$chr = @bed_split[0];
		$start = @bed_split[1];
		$end = @bed_split[2];
		$name = @bed_split[3];
		$strand = @bed_split[5];
		
		#DB::BigWig has a bug that makes it so that the bins retrieved are each an extra
		#base long. So, I'll have to break down each region into $bin_size-length regions
		#myself, then read each one out of the bigwig file.
		
		# tell the BioPerl methods where the bigwig file is
		$bw_file = $in_bw_dir . "/" . $chr . ".bw";
		my $bw = Bio::DB::BigWig->new(-bigwig=>$bw_file);
		
		for($i = 0; $i <= $length; $i++) {
			
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
			} elsif ( $strand eq '-' ) {
				$rel_start = $length - $i;
			}
			
			# Then, update the pileup map
			
			$new_val = @pileup[$rel_start] + $mean;
			$pileup[$rel_start] = $new_val;
			
			
		}
		
		#increment the count
		$bed_count++;
		
	}

	open(OUTPUT, ">$out_file");
	
	# put a header in the output for easier importing into R
	print OUTPUT "pos\tval\n";

	$pos = 0;
	
	foreach(@pileup) {
		$val = $_ / $bed_count;
		print OUTPUT $pos . "\t" . $val . "\n";
		$pos++;
	}
	
	close(OUTPUT);
	close(INPUT);
	
}
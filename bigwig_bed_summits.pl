#!/usr/bin/perl
# bigwig_bed_summits.pl - Uses the BioPerl Bio-BigFile scripts to directly access 
# bigWig scores in each region in a BED file. After retrieval, the region with the highest
# score is identified, and the BED is output with additional information about the summit.
#
# Currently, this is only set up to work with my bedGraph files that have been converted
# to bigWig format using UCSC tools. Some features of this script will have problems if
# a variable-step or single base step wiggle file were used (particularly figuring out 
# what the center of the summit region is - this will just grab the 5'-most base in that
# case).
#
# usage: bigwig_bed_summits.pl in_bigwig in_bed out_bed

if(@ARGV != 3) {
	print "usage: bigwig_bed_summits.pl in_bigwig in_bed out_bed\n";
} else {
	
	#load the BioPerl Bio-BigFile BigWig methods
	use Bio::DB::BigWig 'binMean','binStdev';
	
	$in_bw = @ARGV[0];
	$in_bed = @ARGV[1];
	$out_file = @ARGV[2];
	
	# tell the BioPerl methods where the bigwig file is
	my $bw = Bio::DB::BigWig->new(-bigwig=>$in_bw);
	
	open(INPUT, "<$in_bed");
	open(OUTPUT, ">$out_file");
	
	# put a header in the output, since it's hard to keep track of all the different values
	print OUTPUT "chr\tstart\tend\tname\tscore\tstrand\tmax_start\tmax_end\tmax_length\tmax_score\tmax_center\tmax_center_rel\n";
	
	while(<INPUT>) {
		
		$line = $_;
		chomp($line);
		@bed_split = split(/\t/,$line);
		
		$chr = @bed_split[0];
		$start = @bed_split[1];
		$end = @bed_split[2];
		$strand = @bed_split[5];
		
		#retrieve the values in the bigWig file
		my @features = $bw->features(-seq_id=>$chr,-start=>$start,-end=>$end);
		
		#determine which feature has the highest score
		
		$max_score = 0;
		for my $f (@features) {
			my $score = $f->score;
			
			if($score > $max_score) {
				$max_start = $f->start;
				$max_end = $f->end;
				$max_score = $score;
			}
						
		}
		
		#Now, I have the max_start, max_end, and max_score from the bigwig. I also want to
		#calculate the center of this region, the position of the center relative to
		#the start and end of the BED region (with respect to strand), and the length.
		
		$max_length = $max_end - $max_start;
		
		$max_center = int(($max_start + $max_end)/2);
		
		if( $strand eq '+' ) {
			$max_center_rel = $max_center - $start;
		} elsif ( $strand eq '-' ) {
			$max_center_rel = $end - $max_center;
		}
		
		#Then, print the original BED file plus new columns for the max data
		
		print OUTPUT $line . "\t" . $max_start . "\t" . $max_end . "\t" . $max_length . "\t" . $max_score . "\t" . $max_center . "\t" . $max_center_rel . "\n";
		
	}
	
	close(OUTPUT);
	close(INPUT);
	
}
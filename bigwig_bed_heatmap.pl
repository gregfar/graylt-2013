#!/usr/bin/perl
# bigwig_bed_heatmap.pl - Uses the BioPerl Bio-BigFile scripts to directly access 
# bigWig scores in each region in a BED file. After retrieval, the locations of each 
# bigWig region are converted to a bedgraph-like format for graphing in R using the
# ggplot2 geom_rect tool. This format consists of an ID (from the name column of the bed
# file), the relative start and end positions, the score for the region from the bigwig,
# and a rank based on the order of the bed file.
#
# usage: bigwig_bed_heatmap.pl in_bigwig in_bed out_heat

if(@ARGV != 3) {
	print "usage: bigwig_bed_heatmap.pl in_bigwig in_bed out_heat\n";
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
		
		#retrieve the values in the bigWig file
		my @features = $bw->features(-seq_id=>$chr,-start=>$start,-end=>$end);
		
		#retrieve the start, end and score, then adjust the positions based on the bed region
		
		$max_score = 0;
		for my $f (@features) {
			my $bw_start = $f->start;
			my $bw_end = $f->end;
			my $bw_score = $f->score;
			
			if( $strand eq '+' ) {
				$rel_start = $bw_start - $start;
				$rel_end = $bw_end - $start + 1;
			} elsif ( $strand eq '-' ) {
				$rel_start = $end - $bw_end;
				$rel_end = $end - $bw_start + 1;
			}
			
			# Then, output the adjusted region
			
			print OUTPUT $name . "\t" . $rel_start . "\t" . $rel_end . "\t" . $bw_score . "\t" . $rank . "\n";
			
		}
		
		#increment the rank
		$rank++;
		
	}
	
	close(OUTPUT);
	close(INPUT);
	
}
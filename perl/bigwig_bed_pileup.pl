#!/usr/bin/perl
# bigwig_bed_pileup.pl - Uses the BioPerl Bio-BigFile scripts to directly access 
# bigWig scores in each region in a BED file. After retrieval, the locations of each 
# bigWig region are adjusted based on strand, and then added to a pileup map. The values
# in the map are then divided by the number of regions in the BED file to give an average
# score.
# The regions in the BED file should all be the same size, given as the 'positions' argument.
#
# usage: bigwig_bed_heatmap.pl in_bigwig in_bed positions out_heat

if(@ARGV != 4) {
	print "usage: bigwig_bed_heatmap.pl in_bigwig in_bed positions out_heat\n";
} else {
	
	#load the BioPerl Bio-BigFile BigWig methods
	use Bio::DB::BigWig 'binMean','binStdev';
	
	$in_bw = @ARGV[0];
	$in_bed = @ARGV[1];
	$positions = @ARGV[2] - 1;
	$out_file = @ARGV[3];
	
	# tell the BioPerl methods where the bigwig file is
	my $bw = Bio::DB::BigWig->new(-bigwig=>$in_bw);
	
	open(INPUT, "<$in_bed");
	
	@pileup = ();
	
	foreach(0..$positions) {
		push(@pileup, 0);
	}
	
	$bed_count = 0;
	
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
				$rel_end = $bw_end - $start;
			} elsif ( $strand eq '-' ) {
				$rel_start = $end - $bw_end;
				$rel_end = $end - $bw_start;
			}
			
			# Then, update the pileup map
			
			foreach($rel_start..$rel_end) {
				
				$new_val = @pileup[$_] + $bw_score;
				
				$pileup[$_] = $new_val;
				
			}
						
		}
		
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
#!/usr/bin/perl
# cat_overlapping_bedgraph.pl - takes overlapping regions in bedGraph format and outputs
# single entries based on the outer boundaries of overlapping regions.
# usage: cat_overlapping_regions.pl in_bedgraph out_bedgraph

if (@ARGV != 2 ) {
	
	print "cat_overlapping_regions.pl in_bedgraph out_bedgraph";
	
} else {

	$in_file = @ARGV[0];
	$out_file = @ARGV[1];
	
	open(INPUT, "<$in_file");
	open(OUTPUT, ">$out_file");
	

	
	while(<INPUT>) {
	
		$line = $_;
		chomp($line);
		
		@line_split = split(/\t/,$line);
		
		@region_split = split(/\t/,$region_line);
		
		if ((@line_split[0] eq @region_split[0]) && (@line_split[1] <= @region_split[2])) {

			if(@line_split[1] < @region_split[1]) { $start = @line_split[1]; } else { $start = @region_split[1]; }
			if(@line_split[2] > @region_split[2]) { $end = @line_split[2]; } else { $end = @region_split[2]; }
			if(@line_split[3] > $region_score) { $region_score = @line_split[3]; }

			$region_line = @region_split[0] . "\t" . $start . "\t" . $end . "\t" . $region_score;
			
		} else {
			
			if($notfirstline == 1) {
				print OUTPUT $region_line . "\n";
				$region_score = 0;
			}
			
			$region_line = $line;

			
		}
		
		$notfirstline = 1;
		
	}
	
	print OUTPUT $region_line;
	
	
	close(INPUT);
	close(OUTPUT);

}
#!/usr/bin/perl
# cat_overlapping_regions.pl - takes overlapping regions in BED format and outputs
# single entries based on the outer boundaries of overlapping regions.
# This won't take any strandedness information into account.
# usage: cat_overlapping_regions.pl in_bed out_bed

if (@ARGV != 2 ) {
	
	print "cat_overlapping_regions.pl in_bed out_bed";
	
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
			$region_id = @region_split[3] . "," . @line_split[3];
			$region_line = @region_split[0] . "\t" . $start . "\t" . $end . "\t" . $region_id . "\t" . @region_split[4] . "\t" . @region_split[5];
			
		} else {
			
			if($notfirstline == 1) {
				print OUTPUT $region_line . "\n";
			}
			
			$region_line = $line;

			
		}
		
		$notfirstline = 1;
		
	}
	
	print OUTPUT $region_line;
	
	
	close(INPUT);
	close(OUTPUT);

}
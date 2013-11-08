#!/usr/bin/perl
# bed_to_bedgraph.pl - converts a bed file directly to a simple bedgraph file. Not made for
# bed files with overlapping regions. Use cat_overlapping_regions.pl first to flatten regions.
# usage: bed_to_bedgraph.pl bed_file bedgraph_file

if(@ARGV != 2) {
	print "usage: bed_coverage.pl bed_file bedgraph_file\n";
} else {

	$bed_file = @ARGV[0];
	$bedgraph_file = @ARGV[1];
	
	open(INPUT, "<$bed_file");
	open(OUTPUT, ">$bedgraph_file");
	
	while(<INPUT>) {
		$line = $_;
		chomp($line);
		
		@line_split = split(/\t/,$line);
		
		$coverage = @line_split[2] - @line_split[1];
		
		print OUTPUT @line_split[0] . "\t" . @line_split[1] . "\t" . @line_split[2] . "\t" . 1 . "\n";
		
		
	}
	
	close(INPUT);
	close(OUTPUT);
	

}
#!/usr/bin/perl
# bed_center.pl - outputs the center of each region in a BED file
# usage: bed_center.pl in_file out_file

if(@ARGV != 2) {
	print "usage: bed_center.pl in_file out_file\n";
} else {
	
	$in_file = @ARGV[0];
	$out_file = @ARGV[1];
	
	open(INFILE, "<$in_file");
	open(OUTFILE, ">$out_file");
	
	while(<INFILE>) {
		
		$line = $_;
		chomp($line);
		
		@line_split = split(/\t/,$line);
		
		$center_start = int((@line_split[2] + @line_split[1])/2);
		$center_end = $center_start + 1;
		
		print OUTFILE @line_split[0] . "\t" . $center_start . "\t" . $center_end . "\t" . @line_split[3] . "\t" . @line_split[4] . "\t" . @line_split[5] . "\n";
		
	}
	
	close(OUTFILE);
	close(INFILE);
		
}
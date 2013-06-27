#!/usr/bin/perl
# bed_to_great.pl - converts a bed file directly to a regulatory domain file
# for use with GREAT.
# usage: bed_coverage.pl bed_file rd_file

if(@ARGV != 2) {
	print "usage: bed_coverage.pl bed_file rd_file\n";
} else {

	$bed_file = @ARGV[0];
	$rd_file = @ARGV[1];
	
	open(INPUT, "<$bed_file");
	open(OUTPUT, ">$rd_file");
	
	while(<INPUT>) {
		$line = $_;
		chomp($line);
		
		@line_split = split(/\t/,$line);
		
		$coverage = @line_split[2] - @line_split[1];
		
		print OUTPUT @line_split[0] . "\t" . @line_split[1] . "\t" . @line_split[2] . "\t" . @line_split[3] . "\t" . $coverage . "\t" . @line_split[5] . "\n";
		
		
	}
	

}
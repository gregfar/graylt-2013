#!/usr/bin/perl

if(@ARGV != 5) {
	print "usage: column_filter.pl in_file out_file column [gt|lt] value\n";
} else {

	$input = @ARGV[0];
	$out_file = @ARGV[1];
	$column = @ARGV[2] - 1;
	$filter_type = @ARGV[3];
	$value = @ARGV[4];
		
	open(INPUT, "<$input");
	open(OUTPUT, ">$out_file");
	
	$sum = 0;

	while(<INPUT>) {
		
		$line = $_;
		chomp($line);
		
		@line_split = split(/\t/,$line);
		
		if($filter_type == "gt") {
			
			if(@line_split[$column] > $value) {
				print OUTPUT $line . "\n";
			}
			
		} elsif($filter_type == "lt") {
			
			if(@line_split[$column] < $value) {
				print OUTPUT $line . "\n";
			}
			
		}
		
	}

	close(INPUT);
	close(OUTPUT);

}
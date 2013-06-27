#!/usr/bin/perl
# convert_column.pl - changes the values in one column to the specified value
# usage: convert_column.pl in_file out_file column value

$in_file = @ARGV[0];
$out_file = @ARGV[1];
$column = @ARGV[2] - 1;
$value = @ARGV[3];

open(INPUT, "<$in_file");
open(OUTPUT, ">$out_file");

while(<INPUT>) {
	
	$line = $_;
	chomp($line);
	
	@line_split = split(/\t/, $line);
	
	for($i = 0; $i < @line_split; $i++) {
		if($i == $column) {
			print OUTPUT $value;
		} else {
			print OUTPUT @line_split[$i];
		}
		
		if($i == @line_split - 1) {
			print OUTPUT "\n";
		} else {
			print OUTPUT "\t";
		}
	}
	
}
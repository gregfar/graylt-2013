#!/usr/bin/perl
#remove_lines.pl - looks for instances of a string in the specified column
#outputs a file with only lines that do not contain the string in that column
#usage: remove_lines.pl in_file out_file column string

if(@ARGV != 4) {
	print "usage: remove_lines.pl in_file out_file column string\n";
}

$in_file = @ARGV[0];
$out_file = @ARGV[1];
$column = @ARGV[2] - 1;
$target = @ARGV[3];

open(IN, "<$in_file");
open(OUTPUT, ">$out_file");


while(<IN>) {

	$line = $_;
	chomp($line);
	@columns = split(/\t/, $line);

	if (@columns[$column] ne $target) {

		print OUTPUT $line . "\n";
	
	}

}

close(IN);
close(OUTPUT);
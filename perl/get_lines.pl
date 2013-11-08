#!/usr/bin/perl
#get_lines.pl - looks for instances of a string in the specified column
#outputs a file with only lines that contain the string in that column
#usage: get_lines.pl in_file out_file column string

if(@ARGV != 4) {
	print "usage: get_lines.pl in_file out_file column string\n";
} else {

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

	if (@columns[$column] eq $target) {
		print OUTPUT $line . "\n";
	}

}

close(IN);
close(OUTPUT);

}
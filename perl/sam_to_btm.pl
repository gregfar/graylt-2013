#!/usr/bin/perl
# sam_to_btm.pl - converts paired-end SAM files to bowtie map files.
# usage: sam_to_btm.pl in_sam out_btm

if(@ARGV != 2) {
	print "usage: sam_to_btm.pl in_sam out_btm\n";
} else {
	
$in_file = @ARGV[0];
$out_file = @ARGV[1];

open(INPUT, "<$in_file");
open(OUTPUT, ">$out_file");

while(<INPUT>) {
	
	$line = $_;
	chomp($line);
	
	@line_split = split(/\t/,$line);
	
	if((@line_split[1] == 163) || (@line_split[1] == 99)) {
		print OUTPUT @line_split[0] . "/1\t+\t" . @line_split[2] . "\t" . @line_split[3] . "\t" . @line_split[9] . "\t" . @line_split[10] . "\t0\n";
	} elsif ((@line_split[1] == 83) || (@line_split[1] == 147)) {
		print OUTPUT @line_split[0] . "/2\t-\t" . @line_split[2] . "\t" . @line_split[3] . "\t" . @line_split[9] . "\t" . @line_split[10] . "\t0\n";	
	}
	
}
	
close(INPUT);
close(OUTPUT);
	
}
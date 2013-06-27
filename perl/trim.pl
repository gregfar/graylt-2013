#!/usr/bin/perl
# trim.pl trims FASTQ files to the given number of cycles.
# usage: trim.pl infile outfile cycles

if(@ARGV != 3) {
	print "usage: trim.pl infile outfile cycles\n"
}

$in_file = @ARGV[0];
$out_file = @ARGV[1];
$cycles = @ARGV[2];

open(INPUT, "<$in_file");

open(OUTPUT, ">$out_file");

$x = 0;

while(<INPUT>) {
	
	$line = $_;
	chomp($line);
	
	if($x == 1) {
		
		$trimmed = substr($line, 0, $cycles);
		
		print OUTPUT $trimmed . "\n";
	
		$x = 0;
		
	} else {
		print OUTPUT $line . "\n";
		$x = 1;
	}
	
}

close(INPUT);
close(OUTPUT);
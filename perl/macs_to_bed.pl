#!/usr/bin/perl
# macs_to_bed.pl - converts MACS peak results to BED format by adding a strand column and
# renumbers peaks.
# usage: macs_to_bed.pl in_file out_file name

if(@ARGV != 3) {
	print "macs_to_bed.pl in_file out_file name\n";
} else {
	
	$in_file = @ARGV[0];
	$out_file = @ARGV[1];
	$name = @ARGV[2];
	
	open(INPUT, "<$in_file");
	open(OUTPUT, ">$out_file");
	
	$i = 1;
	
	while(<INPUT>) {
		
		$line = $_;
		chomp($line);
		
		@s = split(/\t/,$line);
		
		print OUTPUT @s[0] . "\t" . @s[1] . "\t" . @s[2] . "\t" . $name . "_" . $i . "\t" . @s[4] . "\t+\n";
		
		$i++;
		
	}
	
}
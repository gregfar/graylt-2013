#!/usr/bin/perl
# macs_to_bed.pl - converts peak summit results to BED format with 
# the specified 5' and 3' padding
# usage: macs_to_bed.pl in_file out_file 5_pad 3_pad

if(@ARGV != 4) {
	print "macs_to_bed.pl in_file out_file 5_pad 3_pad\n";
} else {
	
	$in_file = @ARGV[0];
	$out_file = @ARGV[1];
	$pad5 = @ARGV[2];
	$pad3 = @ARGV[3];
	
	open(INPUT, "<$in_file");
	open(OUTPUT, ">$out_file");
		
	while(<INPUT>) {
		
		$line = $_;
		chomp($line);
		
		@s = split(/\t/,$line);
		
		$start = @s[10] - $pad5;
		$end = @s[10] + $pad3;
		
		print OUTPUT @s[0] . "\t" . $start . "\t" . $end . "\t" . @s[3] . "\t" . @s[4] . "\t" . @s[5] . "\n";
				
	}
	
}
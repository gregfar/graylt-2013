#!/usr/bin/perl
# rmsk_to_bed.pl - converts the UCSC RepeatMasker table to BED format
# usage: rmsk_to_bed.pl in_file out_file

if(@ARGV != 2) {
	print "rmsk_to_bed.pl in_file out_file\n";
} else {
	
	$in_file = @ARGV[0];
	$out_file = @ARGV[1];
	
	open(INPUT, "<$in_file");
	open(OUTPUT, ">$out_file");
		
	while(<INPUT>) {
		
		$line = $_;
		chomp($line);
		
		@s = split(/\t/,$line);
		
		print OUTPUT @s[5] . "\t" . @s[6] . "\t" . @s[7] . "\t" . @s[10] . "\t" . 0 . "\t" . @s[9] . "\n";
				
	}
	
}
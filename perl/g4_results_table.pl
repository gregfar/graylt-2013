#!/usr/bin/perl
# compares two bed files to build a results table that should match GREAT results.


if(@ARGV != 3) {
	print "g4_results_table.pl g4.bed target.bed out_table.txt\n";
} else {
	
	$g4 = @ARGV[0];
	$target = @ARGV[1];
	$output = @ARGV[2];
	
	open(TARGET, "<$target");
	
	while(<TARGET>) {
		
		$line = $_;
		chomp($line);
		
		@line_split = split(/\t/,$line);
		
		push(@{"@line_split[0]"},$line);
		
	}
	
	close(TARGET);
	
	open(G4, "<$g4");
	
	open(OUTPUT, ">$output");
	print OUTPUT "chr\tstart\tend\tsequence\tlength\tstrand\tpeak_id\n";
	
	while(<G4>) {
		
		$line = $_;
		chomp($line);
		
		@g4_split = split(/\t/,$line);
		
		@seq_split = split(/-/,@g4_split[3]);
		
		$chr = @g4_split[0];
		$start = @g4_split[1];
		$end = @g4_split[2];
		$length = $end - $start;
		$center = $start + int(($end - $start) / 2);
		$strand = @g4_split[4];
				
		for($i = 0; $i < @{"$chr"};$i++) {

			@bed_split = split(/\t/,@{"$chr"}[$i]);
			
			if($center < @bed_split[1]) {
			
				last;
			
			} elsif((@bed_split[1] <= $center) && (@bed_split[2] >= $center)) {
			
				$strand = @g4_split[5];
				$id = @seq_split[2];

				if($strand eq "-") {
					
					$sequence = reverse(@seq_split[3]);
					$sequence =~ tr/ACGT/TGCA/;
					
				} else {
					
					$sequence = @seq_split[3];
					
				}
						
				$start = @g4_split[1];
				$end = @g4_split[2];
				$peak_name = @bed_split[3];
				
				print OUTPUT "$chr\t$start\t$end\t$sequence\t$length\t$strand\t$peak_name\n";
				last;
			}
			
		}
		
	}
	
	close(OUTPUT);
	close(G4);
	
}
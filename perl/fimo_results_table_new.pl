#!/usr/bin/perl

if(@ARGV != 3) {
	print "fimo_results_table.pl fimo.txt input.bed out_table.txt\n";
} else {
	
	$fimo = @ARGV[0];
	$bed = @ARGV[1];
	$output = @ARGV[2];
	
	open(BED, "<$bed");
	
	while(<BED>) {
		
		$line = $_;
		chomp($line);
		
		@line_split = split(/\t/,$line);
				
		push(@bed,$line);
		
	}
	
	close(BED);
	
	open(FIMO, "<$fimo");
	
	open(OUTPUT, ">$output");
	print OUTPUT "chr\tstart\tend\tsequence\tp-value\tstrand\tpeak_id\tq-value\n";
	
	while(<FIMO>) {
		
		$line = $_;
		chomp($line);
		
		@fimo_split = split(/\t/,$line);
				
		if(@fimo_split[2] < @fimo_split[3]) {
			$fimo_strand = "+";
			$fimo_start = @fimo_split[2];
			$fimo_end = @fimo_split[3];
		} else {
			$fimo_strand = "-";
			$fimo_start = @fimo_split[3];
			$fimo_end = @fimo_split[2];
		}
		
		for($i = 0; $i < @bed;$i++) {
			
			@bed_split = split(/\t/,@bed[$i]);

			if(@bed_split[3] eq @fimo_split[1]) {
				
				$chr = @bed_split[0];
				$start = @bed_split[1] + $fimo_start;
				$end = @bed_split[1] + $fimo_end;
				$name = @fimo_split[7];
				$p = @fimo_split[5];
				$strand = $fimo_strand;
				$peak_name = @bed_split[3];
				$q = @fimo_split[6];
				
				print OUTPUT "$chr\t$start\t$end\t$name\t$p\t$strand\t$peak_name\t$q\n";
				
				last;
			}
			
		}
		
	}
	
	close(OUTPUT);
	close(FIMO);
	
}
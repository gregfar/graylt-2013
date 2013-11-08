#!/usr/bin/perl
# bigwig_bed_overlap_stranded.pl - Uses the BioPerl Bio-BigFile scripts to directly access 
# bigWig scores using separate files for each strand of the genome in each region in a BED file.
# Overlap counts for each strand and a combined binary overlap count will be reported as additional
# columns added to the end of the BED file.
#
# usage: bigwig_bed_overlap_stranded.pl in_bigwig in_bed out_overlaps

if(@ARGV != 4) {
	print "usage: bigwig_bed_overlap_stranded.pl plus_bw minus_bw in_bed out_overlaps\n";
} else {
	
	#load the BioPerl Bio-BigFile BigWig methods
	use Bio::DB::BigWig 'binMean','binStdev';
	
	$plus_bw = @ARGV[0];
	$minus_bw = @ARGV[1];
	$in_bed = @ARGV[2];
	$out_file = @ARGV[3];
	
	# tell the BioPerl methods where the bigwig files are
	my $pbw = Bio::DB::BigWig->new(-bigwig=>$plus_bw);
	my $mbw = Bio::DB::BigWig->new(-bigwig=>$minus_bw);

	open(INPUT, "<$in_bed");
	open(OUTPUT, ">$out_file");
	
	$first = 1;
	
	while(<INPUT>) {
		
		if($first == 1) {
			
			print OUTPUT "chr\tstart\tend\tname\tscore\tstrand\tt.overlap\tnt.overlap\tsum.overlap\n";
			$first = 0;
			
		} else {
		
			$line = $_;
			chomp($line);
			@bed_split = split(/\t/,$line);
			
			$chr = @bed_split[0];
			$start = @bed_split[1];
			$end = @bed_split[2];
			$strand = @bed_split[5];
			
			$hit_nt = 0;
			$hit_t = 0;
			$hit_sum = 0;
			
			#retrieve the values in the bigWig files
			my @pfa = $pbw->features(-seq_id=>$chr,-start=>$start,-end=>$end);
			my @mfa = $mbw->features(-seq_id=>$chr,-start=>$start,-end=>$end);
			
			#retrieve the start, end and score, then adjust the positions based on the bed region
			# for the positive strand.
			
			for my $pf (@pfa) {
				my $bw_start = $pf->start;
				my $bw_end = $pf->end;
				my $bw_score = $pf->score;
				
				if( $strand eq '+' ) {
					$hit_nt++;
				} elsif ( $strand eq '-' ) {
					$hit_t++;
				}
				
			}
			
			#Do the same for the negative strand
			
			for my $mf (@mfa) {
				my $bw_start = $mf->start;
				my $bw_end = $mf->end;
				my $bw_score = $mf->score;
				
				if( $strand eq '+' ) {
					$hit_t++;
				} elsif ( $strand eq '-' ) {
					$hit_nt++;	
				}
				
			}
			
			$hit_sum = $hit_t + $hit_nt;
			
			print OUTPUT $line . "\t" . $hit_t . "\t" . $hit_nt . "\t" . $hit_sum . "\n";
		
		}	
		
	}
	
	close(INPUT);
	close(OUTPUT);

	
}
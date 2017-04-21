#!/usr/local/bin/perl -w

#perl SAMtoBED.pl myFile.SAM

$sam = $ARGV[0];

if (! $ARGV[0])	{	
print STDERR "\nConvert a SAM file to a BED file\n";
print STDERR "USAGE: $0 in.sam > out.bed\n\n";
exit;
}

open (SAM, "$sam") || die "Major problem: cannot open $sam for reading: $!";

while (<SAM>) {
	if (! /^@/){
	chomp($_);					
	@data = split /\t/, $_;
	$start = $data[3] - 1;
	$length = length $data[9];
 	$end = $start + $length;
	$bedRow = "$data[2]\t$start\t$end";
		if ($data[1] == 0) {
			$strand = "+";
		}
		elsif ($data[1] == 16) {
			$strand = "-";
		}
		else {
			print STDERR "Strand field ($data[1]) is not recognized\n";
		}
		# Get the number of mismatches (0 - 2)
		$numMismatches = substr ($data[13], -1, 1);
		# Create a score that's 2 minus the number of mismatches
		$score = 2 - $numMismatches;
		$bedRow .= "\t$.\t$score\t$strand";
		print "$bedRow\n";
	}
}
close (SAM);

########

sub samExtractor {
my ($refSeq, $bam, $chr, $st, $end)=@_;
 use Bio::DB::Sam;
 my $sam = Bio::DB::Sam->new(-bam  =>"$bam",
                             -fasta=>"$refSeq",
                             );

 my @targets    = $sam->seq_ids;
 my @alignments = $sam->get_features_by_location(-seq_id => "$chr",
                                                 -start  => $st,
                                                 -end    => $end);
 for my $a (@alignments) {

    # where does the alignment start in the reference sequence
    my $seqid  = $a->seq_id;
    my $start  = $a->start;
    my $end    = $a->end;
    my $strand = $a->strand;
    my $cigar  = $a->cigar_str;
    #my $paired = $a->get_tag_values('PAIRED');

    # where does the alignment start in the query sequence
    my $query_start = $a->query->start;     
    my $query_end   = $a->query->end;

    my $ref_dna   = $a->dna;        # reference sequence bases
    my $query_dna = $a->query->dna; # query sequence bases

    my @scores    = $a->qscore;     # per-base quality scores
    my $match_qual= $a->qual;       # quality of the match
 }

}

##########

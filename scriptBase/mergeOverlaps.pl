#!/usr/bin/perl
use strict;
use warnings;
use 5.010;

# Filter out the range overlaps from tab seperated alignment file. (filterlap overlaps script outfile)
# USAGE: perl mergeOverlaps.pl infile strand/blocks > outfile

open my $fh, '<', $ARGV[0];

my @terms;
while (<$fh>) {
    chomp;
    push @terms, [split /\t/];
}

my $biggest = 0;
my $id = '';
my @ar; my @art;
my $strand;
my $cnt=0;
my $signal=''; my $ofSize;
my @allids;
my $refName="RefName";
my $tarName="TarName";
my $type = "scaffolds";


#Strand fundtion here ... not fully funtional 

#chr1	3000701	3000741	WIGTC-HISEQ:5:1208:16926:179883#ACAGTG/1;1	255	+	3000701	3000741	255,0,0	1	40	0
#725478	NODE_2_length_7674_cov_46.7841_ID_3	+	7674	0	7674	NODE_2_length_7674_cov_46.7841_ID_3	+	7674	0	7674	7674/7674	100.0%	7674/7674	100.0%

if ($ARGV[1] eq "strand") {
for my $term (sort sorter @terms) {
	if ((($id ne $term->[1]) or ($strand ne $term->[2])) and ($cnt ne 0)) {
         	my @ar2 = sort {$a <=> $b} (@ar); my $nextStrand='';
		# Note: Overlappingfragments have different orientation in tar.
		if ($term->[4] < $ar[-1]) { $signal="OverlapFrag !!"; $ofSize=$term->[4]-$ar[-1];} else {$signal="NoOverlapFrag"; $ofSize=$term->[4]-$ar[-1]}
		if ($id ne $term->[1]) {$signal="END"; $ofSize="NA"; $nextStrand="NA"} else {$nextStrand=$term->[2]}
		my $fids = join( ",", (uniq(@allids)));		
		print "$refName\t$term->[1]\t$ar2[0]\t$ar[-1]\t$fids\tNA\tNA\t$strand\t$tarName\t$type\t$signal\t$ofSize\t$nextStrand **\n";
		undef @ar; $biggest = 0; undef @allids;
		}
    if ($term->[4] >= $biggest) {
        say join "\t", @$term;
	push @ar, ($term->[4], $term->[5]);
	push @art, ($term->[9], $term->[10]);
	push @allids, $term->[6]; 
        $biggest = $term->[4];
    }
    $id = $term->[1];
    $strand =  $term->[2];
    if ($cnt == $#terms) {
	my @ar2 = sort {$a <=> $b} (@ar);
	my $fids = join( ",", (uniq(@allids))); 	
	print "$refName\t$term->[1]\t$ar2[0]\t$ar[-1]\t$fids\tNA\tNA\t$strand\t$tarName\t$type\tEND\tEND\tNA ***\n"; }
    $cnt++;
} }

#General fundtion here ... funtional ... Just reformat 

#chr1	3000701	3000741	WIGTC-HISEQ:5:1208:16926:179883#ACAGTG/1;1	255	+	3000701	3000741	255,0,0	1	40	0
#725478	NODE_2_length_7674_cov_46.7841_ID_3	+	7674	0	7674	NODE_2_length_7674_cov_46.7841_ID_3	+	7674	0	7674	7674/7674	100.0%	7674/7674	100.0%
#RefName	NODE_9_length_4027_cov_38.3482_ID_17	0	4027	NODE_9_length_4027_cov_38.3482_ID_17	NA	NA	+	TarName	scaffolds	END	END	NA ***
elsif ($ARGV[1] eq "general") {
for my $term (sort sorter @terms) {

	print "$refName\t$term->[1]\t$term->[4]\t$term->[5]\t$term->[6]\t$term->[9]\t$term->[10]\t$term->[7]\tTarName\tscaffolds\t$term->[12]\t$term->[14]\n"; 

}}


sub sorter {
	$a->[1] cmp $b->[1] ||
     $a->[4] <=> $b->[4]
  || $b->[5] <=> $a->[5]
  || $b->[6] cmp $a->[6]
  || $b->[9] <=> $a->[9]
  || $b->[10] <=> $a->[10]
}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}



__END__

#!/usr/bin/perl
use strict;
use warnings;
use 5.010;

# Filter out the range overlaps from tab seperated alignment file. (filterlap overlaps script outfile)
# USAGE: perl mergeOverlaps.pl infile strand/blocks > outfile

open my $fh, '<', $ARGV[0];

my @terms;
while (<$fh>) {
    chomp;
    push @terms, [split /\t/];
}

my $biggest = 0;
my $id = '';
my $id2 = '';
my @ar; my @art;
my $strand;
my $cnt=0;
my $signal=''; my $ofSize;
my @allids;
my $refName="RefName";
my $tarName="TarName";
my $type = "scaffolds";

if ($ARGV[1] eq "strand") {
for my $term (sort sorter @terms) {
	if ((($id ne $term->[4]) or ($strand ne $term->[7])) and ($cnt ne 0)) {
         	my @ar2 = sort {$a <=> $b} (@ar); my $nextStrand='';
		#Note: Overlappingfragments have different orientation in tar.
		if ($term->[4] < $ar[-1]) { $signal="OverlapFrag !!"; $ofSize=$term->[4]-$ar[-1];} else {$signal="NoOverlapFrag"; $ofSize=$term->[4]-$ar[-1]}
		if ($id ne $term->[4]) {$signal="END"; $ofSize="NA"; $nextStrand="NA"} else {$nextStrand=$term->[7]}
		my $fids = join( ",", (uniq(@allids)));		
		print "$refName\t$term->[0]\t$ar2[0]\t$ar[-1]\t$fids\tNA\tNA\t$strand\t$tarName\t$type\t$signal\t$ofSize\t$nextStrand\n";
		undef @ar; $biggest = 0; undef @allids;
		}
    if ($term->[4] >= $biggest) {
        say join "\t", @$term;
	push @ar, ($term->[4], $term->[5]);
	push @art, ($term->[9], $term->[10]);
	push @allids, $term->[6]; 
        $biggest = $term->[4];
    }
    $id = $term->[4];
    $strand =  $term->[7];
    if ($cnt == $#terms) {
	my @ar2 = sort {$a <=> $b} (@ar);
	my $fids = join( ",", (uniq(@allids))); 	
	print "$refName\t$term->[0]\t$ar2[0]\t$ar[-1]\t$fids\tNA\tNA\t$strand\t$tarName\t$type\tEND\tEND\tNA\n"; }
    $cnt++;
} }


sub sorter {
	$a->[1] cmp $b->[1] ||
     $a->[4] <=> $b->[4]
  || $b->[5] <=> $a->[5]
  || $b->[6] cmp $a->[6]
  || $b->[9] <=> $a->[9]
  || $b->[10] <=> $a->[10]
}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}



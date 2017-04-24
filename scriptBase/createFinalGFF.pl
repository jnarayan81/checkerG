#!/usr/bin/perl
use strict;
use warnings;

#Create the final GFF file for visualization in IGV
#Jitendra

my $filename = $ARGV[0];
my $exSize=$ARGV[3];
my $avgCov=$ARGV[5];

open(my $fh, '<:encoding(UTF-8)', $filename) or die "Could not open file '$filename' $!";

while (my $row = <$fh>) {
  chomp $row;
  #refname	contig_10_modified	tarname	62721	62823	Break	contig_10	0	1	contig_10	64019	64020
  my @data = split /\t/, $row;
  next if $data[5] eq 'PseudoBreak'; #Ignoring all the Pseudo-breaks
  my $breakST=$data[3]-$exSize; my $breakED=$data[4]+$exSize;
  my $brkSize=$data[4]-$data[3]; 
  my $repeats='NA'; my $color='#FF9900'; my $palDecision='';
	#check the repeats in blocks then
	my $exSeq = uc (extractSeq($ARGV[4],$data[1],$breakST,$breakED));
	my $palRes = checkPal($exSeq);
	if ( $palRes ) { $palDecision="$palRes";} else { $palDecision='NotPalindromic'; }


	my $trfR=checkTRF ($ARGV[1], $data[1], $breakST, $breakED); #Check the repeats
	my @trfR=@$trfR; my $trfStr = join ',', @trfR;
	if (@trfR) {$repeats="REPEATS in breakpoint present"; $color='#181009';} else { $repeats="NO REPEATS in breakpoint region";}
	my $coverage=checkCov ($ARGV[2], $data[1], $breakST, $breakED, $avgCov); 

  print "$data[1]\tcheckerG\tBreakpoint[$data[3]-$data[4]]\t$breakST\t$breakED\t.\t+\t.\tcolor=$color; PalRes:$palDecision; GlobalAvarageCoverage:$avgCov; Coverage=$coverage; Real_breakSize=$brkSize; BreakStrand=+/-; BREAKS=Break size by extending the coordinates +/- $exSize; Repeats:$repeats : $trfStr\n";
}



# Checks if a provided two coordinates overlaps or not it return 1 if overlaps
sub checkCorOverlaps {
my ($x1, $x2, $y1, $y2)=@_;
return $x1 <= $y2 && $y1 <= $x2;
}

sub checkTRF {
my ($file, $name, $cor1, $cor2)=@_;
my @trfData;
open(my $fh, '<:encoding(UTF-8)', $file) or die "Could not open file '$file' $!";
while (my $row = <$fh>) {
  	chomp $row;
	my @data = split /\t/, $row;
	next if $name ne $data[0];
	my $res=checkCorOverlaps ($data[1], $data[2], $cor1, $cor2);
	if ($res) {push @trfData, $data[6]}

}
return \@trfData;
}

sub checkCov {
my ($file, $name, $cor1, $cor2, $avgCov)=@_;
my @covData; my $sum=0; my $avCov=0;
open(my $fh, '<:encoding(UTF-8)', $file) or die "Could not open file '$filename' $!";
while (my $row = <$fh>) {
  	chomp $row;
	#Contig_0	0	1	99
	my @data = split /\t/, lc ($row);
	next if $name ne $data[0];
	my $res=checkCorOverlaps ($data[1], $data[2], $cor1, $cor2);
	if ($res) {push @covData, $data[3]}

}
$sum += $_ for @covData;
$avCov=$sum/scalar(@covData);

return $avCov;
}

sub checkPal {
my ($seq)=@_;
my $pp = qr/(?: (\w) (?1) \g{-1} | \w? )/ix;
    while ($seq =~ /(?=($pp))/g) {
        return "$-[0] - $1" if length($1) > 10; #Palindrome of minimum size 50
    }
}

#Note: Case sensistive contig name match
sub extractSeq {
my ($file, $chr, $st, $ed)=@_;
use Bio::DB::Fasta;
my $db = Bio::DB::Fasta->new($file);
my $seq = $db->seq($chr, $st => $ed);
return $seq;
}


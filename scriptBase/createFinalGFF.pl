#!/usr/bin/perl
use strict;
use warnings;

my $filename = $ARGV[0];
my $exSize=$ARGV[3];

open(my $fh, '<:encoding(UTF-8)', $filename) or die "Could not open file '$filename' $!";
 
while (my $row = <$fh>) {
  chomp $row;
  my @data = split /\t/, $row;
  next if $data[5] eq 'PseudoBreak';
  my $breakST=$data[3]-$exSize; my $breakED=$data[4]+$exSize;
  my $blkSize=$data[4]-$data[3]; 
  my $repeats='NA'; my $color='#FF9900'; my $palDecision='';
  
	#check the repeats in blocks then
	my $exSeq = extractSeq($ARGV[4],$data[1], $data[3],$data[4]);
	my $palRes = checkPal($exSeq);
	if ( $palRes ) { $palDecision="Palindromic:$palRes";} else { $palDecision='NotPalindromic'; }
	$color='#181009';


	my $trfR=checkTRF ($ARGV[1], $data[1], $breakST, $breakED);
	my @trfR=@$trfR; my $trfStr = join ',', @trfR;
	if (@trfR) {$repeats="REPEATS in breakpoint present";}
	my $coverage=checkCov ($ARGV[2], $data[1], $breakST, $breakED);

  print "$data[1]\tCheckER\tBrekapoint, Blocks[$data[2]-$data[3]] have repeats:$breakST\t$breakED\tcolor=$color; $palDecision; Coverage=$coverage; blockSize=$blkSize; BreakStrand=+/-; BREAKS=Break size by extending the coordinates +/- $exSize; $repeats; $trfStr\n";
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
my ($file, $name, $cor1, $cor2)=@_;
my @covData; my $sum=0; my $avCov=0;
open(my $fh, '<:encoding(UTF-8)', $file) or die "Could not open file '$filename' $!";
while (my $row = <$fh>) {
  	chomp $row;
	my @data = split /\t/, $row;
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
        return "$-[0] - $1" if length($1) > 10;
    }
}


sub extractSeq {
my ($file, $chr, $st, $ed)=@_;
use Bio::DB::Fasta;
my $db = Bio::DB::Fasta->new($file);
my $seq = $db->seq($chr, $st => $ed);
return $seq;
}

#!/usr/bin/perl
use strict;
use warnings;
use 5.010;

# Filter out the exact/direct overlaps from tab seperated alignment file. (lastz format=general- ready)
# Do not inclide header in lastz outfile 

# USAGE: perl filterOverlaps.pl infile > outfile

open my $fh, '<', $ARGV[0];

my @terms;
while (<$fh>) {
    chomp;
    push @terms, [split /\t/];
}

my $biggest = 0;
my $id = '';

for my $term (sort sorter @terms) {
	$biggest = 0 if $id ne $term->[1];
    if ($term->[5] > $biggest) {
        say join "\t", @$term;
        $biggest = $term->[5];
    }
    $id = $term->[1];    
}

sub sorter {
	$a->[1] cmp $b->[1] ||
     $a->[4] <=> $b->[4]
  || $b->[5] <=> $a->[5]
}


#!/usr/bin/perl
use strict;
use warnings;

my $filename = $ARGV[0];
my $averageCov=$ARGV[1];

open(my $fh, '<:encoding(UTF-8)', $filename) or die "Could not open file '$filename' $!";
 
while (my $row = <$fh>) {
  chomp $row;
  my @data = split /\t/, $row;
  if ($data[3] >= ($averageCov * 2)) {
	print "$row\n";
}
close $fh;

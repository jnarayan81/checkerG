#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
use File::Basename;
use lib dirname (__FILE__); # check for local folder for modules

use checkerMOD;

open my $InFile, '<', $ARGV[0];
open my $OutFile, '>', $ARGV[1]; #Write the brk result in this file
my $chrFile= $ARGV[2];
my $increase= $ARGV[3];

my @array;
while (<$InFile>) {
    chomp;
    push @array, $_;
}


#chicken:300K	1	671505	1155740	129	19828	735441	+	Pygoscelis_adeliae	Scaffolds
#RefName	NODE_2_length_7674_cov_46.7841_ID_3	0	7674	NODE_2_length_7674_cov_46.7841_ID_3	0	7674	+	TarName	scaffolds	100.0%	100.0%

my @sorted_array = sort { (split "\t", $a)[1] cmp (split "\t", $b)[1] || (split "\t", $a)[2] <=> (split "\t", $b)[2] && (split "\t", $a)[3] <=> (split "\t", $b)[3] } @array;
push @sorted_array, "==="; ## End line to terminate;

checkerMOD::createBreaks(\@sorted_array, $increase, $InFile, $OutFile, $chrFile);

close $InFile;
close $OutFile;



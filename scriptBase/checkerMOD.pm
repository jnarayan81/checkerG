# $Id$ checkerMOD
# Perl module for checkerG;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2017 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

checkerMOD  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

=head1 CONTACT

Jitendra <jnarayan81@gmail.com>

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##

package checkerMOD;
use strict;
use warnings;

use Exporter;

our @EXPORT_OK = "checkerMOD";

sub createBreaks {

my ($sorted_array_ref, $increase, $InFile, $OutFile, $chrFile)= @_;
my @sorted_array=@$sorted_array_ref;

my $count; my $chr_num;
my (@ref_org, @chr_ref, @brk_cor1, @brk_cor2, @tar_org, @chrId, @chrStart, @chrEnd, @sign);
my $telomere1 = 0;
foreach my $line (@sorted_array) {
	chomp $line; 
	$line=lc($line);	
	if ($line =~ /^\s*#/) { next; }
	if ($line =~ /^$/) { next;} ## if a blank line
	$count++;
	my @tmp = split/\t/, lc($line);
	for (@tmp)  { s/^\s+//; s/\s+$//; } #replace one or more spaces 
	if ($count == 1) { $chr_num = $tmp[1];}

	##--------------------- increase or decrease the breakpoints size start -----------
 
	if ($line =~ /^\=/) { $tmp[1]=0; $tmp[2]=0; $tmp[3]=0; }         ## change here in increments .. need to check????
	my $tmp3cor1=$tmp[3]-$increase;      
	my $tmp2cor2=$tmp[2]+$increase;
	
	##-------------------- increase or decrease the breakpoints size end ----------- 

		if ($tmp[1] eq $chr_num) {
			
			if ($telomere1 == 0) { print $OutFile "$tmp[0]\t$tmp[1]\t$tmp[8]\t1\t$tmp2cor2\tPseudoBreak\t0\t0\t0\t0\t0\t0\n"; $telomere1 = 1}
			push @ref_org, $tmp[0];
			push @chr_ref, $tmp[1];
			push @brk_cor1, $tmp3cor1;
			push @brk_cor2, $tmp2cor2;

			push @tar_org, $tmp[8];
			push @chrId, $tmp[4];
			push @chrStart, $tmp[6];
			push @chrEnd, $tmp[5];
			push @sign, $tmp[7];
			}

		else    {
			for my $x(0.. $#ref_org) {
				if ($x != $#ref_org) {
					my $targetSpeciesCoordiPahala=$chrStart[$x]+1; ## added for getting target species coordinates ... If increasing the size we need to think of it !!!!!
					my $targetSpeciesCoordiDusara=$chrEnd[$x+1]-1;
					print $OutFile "$ref_org[$x]\t$chr_ref[$x]\t$tar_org[$x]\t$brk_cor1[$x]\t$brk_cor2[$x+1]\tBreak\t$chrId[$x]\t$chrStart[$x]\t$targetSpeciesCoordiPahala\t$chrId[$x+1]\t$targetSpeciesCoordiDusara\t$chrEnd[$x+1]\n";
					# print  "$ref_org[$x]\t$chr_ref[$x]\t$tar_org[$x]\t$brk_cor1[$x]\t$brk_cor2[$x+1]\n";
						if ($brk_cor1[$x] >= $brk_cor2[$x+1]) { 
							print "Overlappings $chr_ref[$x],$brk_cor1[$x],$brk_cor2[$x+1]\n";
						}
					undef $targetSpeciesCoordiPahala; undef $targetSpeciesCoordiDusara;
					}
				}
			my $maxChrSize=findChrSize($chr_ref[-1], $chrFile);
			if (!$maxChrSize) { print "Fail to open $chr_ref[-1]\n";}
			print $OutFile "$ref_org[-1]\t$chr_ref[-1]\t$tar_org[-1]\t$brk_cor1[-1]\t$maxChrSize\tPseudoBreak\t0\t0\t0\t0\t0\t0\n";
			$chr_num=$tmp[1];

			undef @ref_org; undef @chr_ref; undef @brk_cor1; undef @brk_cor2; undef @tar_org; undef @sign, undef @chrId, undef @chrStart, undef @chrEnd;

			push @ref_org, $tmp[0];
			push @chr_ref, $tmp[1];
			push @brk_cor1, $tmp3cor1;
			push @brk_cor2, $tmp2cor2;

			push @tar_org, $tmp[8];
			push @chrId, $tmp[4];
			push @chrStart, $tmp[5];
			push @chrEnd, $tmp[6];
			push @sign, $tmp[7];
			
		if ($tmp[1]) { print $OutFile "$tmp[0]\t$tmp[1]\t$tmp[8]\t1\t$tmp2cor2\tPseudoBreak\t0\t0\t0\t0\t0\t0\n";}
			}
}
undef @ref_org; undef @chr_ref; undef @brk_cor1; undef @brk_cor2; undef @tar_org; undef @sign, undef @chrId, undef @chrStart, undef @chrEnd;
} ##create breaks ends here 

##Find the maximum chromosome size

sub findChrSize {
my $chr=shift;
my $chrFile=shift;
my %hash;
open(CHRFILE, "$chrFile") || warn "fail to open chr file";
while (<CHRFILE>) { chomp; $_=trim($_); my ($key, $val) = split /\t/, lc($_); $hash{$key} = $val;} ### We can read and store it ... !!!!
close CHRFILE or die "cant close";
#foreach my $key ( sort {$a <=> $b} keys %hash){ push (@arrayChr, $key); } ## Store the chromosome data

return $hash{$chr};
}

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^[\t\s]+//;
	$string =~ s/[\t\s]+$//;
	$string =~ s/[\r\n]+$//; ## remove odd or bad newline ...
	return $string;
}

1;

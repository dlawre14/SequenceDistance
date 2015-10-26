#!/usr/bin/perl

# usage:
#
# perl ./JSD.pl -f ffp1.txt -s ffp2.txt -o JSD.txt
#
# written by: Se-Ran Jun

use warnings;
use strict;
use Getopt::Std;

my %args;

# read arguments ...
getopt('f:s:o', \%args);

my $ffp1_file = $args{f}; # ffp profile
my $ffp2_file = $args{s}; # ffp profile
my $output = $args{o}; # output

die "usage: ./JSD.pl -f <first ffp profile> -s <second ffp profile> -o <Jensen-Shannon Divergence>\n"
unless defined ($ffp1_file) && defined ($ffp2_file) && defined ($output);

my $log2 = log(2);

# read first ffp profile
my %ffp1 = (); my $sum = 0.0;
open (IN, "<$ffp1_file");
while(<IN>){

	chomp $_;
	$_ = trimwhitespace($_);
	(my $ft, my $fq) = split(' ',$_);

	$ffp1{$ft} = $fq;
	$sum = $sum + $fq;
}
close IN;

foreach my $ft ( keys %ffp1 ) {
	$ffp1{$ft} = $ffp1{$ft}/$sum; 
}

# read second ffp profile
my %ffp2 = (); $sum = 0.0;
open (IN, "<$ffp2_file");
while(<IN>){

	chomp $_;
	$_ = trimwhitespace($_);
	(my $ft, my $fq) = split(' ',$_);

	$ffp2{$ft} = $fq;
	$sum = $sum + $fq;
}
close IN;

foreach my $ft ( keys %ffp2 ) {
	$ffp2{$ft} = $ffp2{$ft}/$sum; 
}

# Jensen-Shannon divergence
my %affp=();
foreach my $ft ( keys %ffp1 ) {
 	if (defined $ffp2{$ft} ) {
    		$affp{$ft} = ( $ffp1{$ft} + $ffp2{$ft} )/2;
	}else{
    		$affp{$ft} = $ffp1{$ft}/2;
	}
}

foreach my $ft (keys %ffp2 ) {
	if ( !(defined $ffp1{$ft}) ) {
    		$affp{$ft} = $ffp2{$ft}/2;
	}
}

my $JSD = 0.0;
foreach my $ft ( keys %affp ) {
    	$JSD = $JSD - $affp{$ft} * log($affp{$ft})/$log2;
    	if ( defined $ffp1{$ft}) {
        	$JSD = $JSD + 0.5 * $ffp1{$ft} * log($ffp1{$ft})/$log2;
	}
    	if ( defined $ffp2{$ft} ) {
        	$JSD = $JSD + 0.5 * $ffp2{$ft} * log($ffp2{$ft})/$log2;
	}
}

open (OUT, ">$output");
printf( OUT "%17.15f\n", $JSD );
close OUT;

exit;

###############################################################################################
# SUBROUTINE
###############################################################################################

# trimwhitespace : remove whitespace from the start and end of the string
sub trimwhitespace {

        use strict;
        use warnings;

        my $string=shift;

        $string=~s/^\s+//;
        $string=~s/\s+$//;

        return $string;
}

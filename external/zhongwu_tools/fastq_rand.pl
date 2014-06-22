#!/usr/bin/perl -w

use Getopt::Std;
use strict;

our ( $opt_p );

getopts('p:');

my $p = $opt_p ? $opt_p : 0.10; # the probability
my %hash;
my ($R1, $R2) = @ARGV;
my $n = 0;
my $c = 0;

my ($O1, $O2) = ("$R1.$p", "$R2.$p");
$R1 = "gunzip -c $R1 |" if ( $R1 =~ /.gz$/ );
open( R1, $R1 );
open( O1, ">$O1" );
my $seq = "";
while(<R1>) {
    $n++;
    $_ = "+\n" if ( $n%4 == 3 );
    $seq .= $_;
    if ( $n%4 == 0 ) {
	$c++;
	if ( rand() < $p ) {
	    $hash{ $c } = 1;
	    print O1 $seq;
	}
	$seq = "";
    }
}
close(O1);
close(R1);
if ( $R2 ) {
    $R2 = "gunzip -c $R2 |" if ( $R2 =~ /.gz$/ );
    open( R2, $R2 );
    open( O2, ">$O2" );
    $n = 0;
    $c = 0;
    $seq = "";
    while(<R2>) {
        $n++;
	$_ = "+\n" if ( $n%4 == 3 );
	$seq .= $_;
	if ( $n%4 == 0 ) {
	    $c++;
	    if ( $hash{ $c } ) {
		print O2 $seq;
	    }
	    $seq = "";
	}
    }
    close(R2);
    close(O2);
}

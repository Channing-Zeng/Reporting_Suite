#!/usr/bin/perl -w

use Getopt::Std;
use strict;

our ($opt_n);

getopts( 'n:' );
my %hash;
while( <> ) {
    chomp;
    my $s = $1 if ( /([^\/]+)_L0/ );
    $s = $1 if ( $opt_n && /$opt_n/ );
    s/\*$//;
    s/\@$//;
    push( @{ $hash{ $s }}, $_ );
}
while(my ($s, $r) = each %hash ) {
    print join("\t", $s, @$r), "\n";
}

#!/usr/bin/perl -w

# used when there're more than one fastq files per sample.
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
    push( @{ $hash{ $s }->{ R1 }}, $_ ) if ( /_R1/ );
    push( @{ $hash{ $s }->{ R2 }}, $_ ) if ( /_R2/ );
}
while(my ($s, $r) = each %hash ) {
    print join("\t", $s, join(",", @{$r->{ R1 }}), join(",", @{$r->{ R2 } })), "\n";
}

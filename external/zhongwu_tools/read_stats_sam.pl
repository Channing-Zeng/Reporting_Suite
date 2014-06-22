#!/usr/bin/perl -w

# The program will generate the read stats from sam file

use Getopt::Std;
use strict;

our ($opt_s);
getopts('s:');
my %count;
my %hash;
my $total = 0;
my $pairs = 0;
my $reads = 0;
my $mapped = 0;
my $unmapped = 0;
if ( $opt_s ) {
    keys %hash = $opt_s;
}
while( <> ) {
    my @a = split();
    my $r = $a[0];
    next if ( $r =~ /^@/ );
    $total++;
    $pairs++ if ( $a[6] ne "*" );
    if ( $hash{ $r } ) {
        my $p = $hash{ $r };
	$count{ $p }--;
	$count{ $p + 1 }++;
	$hash{ $r }++;
    } else {
        ## Assuming unmapped reads are output only once
	if ( $a[2] ne "*" ) {
	    $hash{ $r } = 1;
	    $count{ 1 }++;
	    $mapped++
	} else {
	    $unmapped++;
	}
	$reads++;
    }
}

print "Total Reads: $reads\n";
print "Mapped Reads: $mapped\n";
print "Unmapped Reads: $unmapped\n";
print "Reads are mapped to $total locations.\n";
print "Mapped pairs: $pairs\n";
foreach( sort { $a <=> $b } keys %count ) {
    print "Reads mapped $_ times: $count{ $_ }\n" if ( $count{ $_ } > 0 );
}

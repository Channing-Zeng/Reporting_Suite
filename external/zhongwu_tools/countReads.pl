#!/usr/bin/perl -w

# Count the total reads in a fastq

use Getopt::Std;
use strict;

our ($opt_N, $opt_n, $opt_i);
getopts('in:N:');
my $file = $ARGV[0];
my $sample;
$sample = $1 if ( $ARGV[0] && $opt_n && $ARGV[0] =~ /$opt_n/ );
$sample = $opt_N if ( $opt_N );
my $count = 0;
if ( $file && $file =~ /gz$/ ) {
    $count = `gunzip -c $file | wc -l`;
} elsif ($file) {
    $count = `wc -l $file`;
} else {
    while( <> ) {
        $count++;
    }
}
chomp $count;
$count = $opt_i ? $count/4 : $count/2;
print "$sample\t$count\n";

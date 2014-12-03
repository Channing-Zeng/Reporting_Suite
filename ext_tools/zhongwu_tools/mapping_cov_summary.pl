#!/usr/bin/perl -w

# the program will produce the mapping coverage summary

use Getopt::Std;
use Stat::Basic;
use strict;

our ($opt_d, $opt_s, $opt_c);
getopts('d:s:c:');
my @depths = $opt_d ? (split(/:/, $opt_d)) : (1, 5, 10, 25, 50, 100, 500, 1000, 5000, 10000, 50000);
my $stat = new Stat::Basic;
my %hash;
my %target;
my %cov;
my $col_s = $opt_s ? $opt_s - 1 : 6; # amplicon size
my $col_c = $opt_c ? $opt_c - 1 : 7; # mean amplicon coverage
while( <> ) {
    chomp;
    next unless(/Whole/);
    next if ( /^Sample/ );
    my @a = split(/\t/);
    $target{ $a[0] } += $a[$col_s];
    $cov{ $a[0] } += $a[$col_s]*$a[$col_c];
    for(my $i = $col_c+1; $i < @a; $i++) {
	$hash{ $a[0] }->{ $depths[$i-$col_c-1] } += $a[$col_s] * $a[$i];
	#print $depths[$i-5], "\t", $hash{ $a[0] }->{ $depths[$i-5] }, "$_\n";
    }
}
print join("\t", "Sample", "MeanDepth", (map { $_ . "X"; } @depths)), "\n";
my @samples = sort(keys %hash);
foreach my $s (@samples) {
    print join("\t", $s, int($cov{$s}/$target{$s}), (map { sprintf("%.4f", $hash{ $s }->{ $_ }/$target{$s}); } @depths)), "\n";
}

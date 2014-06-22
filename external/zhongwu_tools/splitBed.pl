#!/usr/bin/perl -w

# Split bed file
use strict;

my $bed = shift;
my $size = shift;

my ($lines, $t) = split(/\s+/, `wc -l $bed`);
my $seg = int($lines/$size);
$seg++ if ( $lines % $size > 0 && $lines % $size <= $seg );

my $n = 0;
my $N = 0;
my $out;
open(BED, $bed);
my $base = `basename $bed`; chomp $base;
while(<BED>) {
    if ( $n % $seg == 0 ) {
        close($out) if ( $out );
	$N++;
	open( $out, "> $base.$N" );
    }
    print $out $_;
    $n++;
}
close( $out );
close( BED );

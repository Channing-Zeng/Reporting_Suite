#!/usr/bin/perl -w

use Getopt::Std;
use strict;

our ($opt_s);
getopts('s:');

my $sample = $opt_s ? $opt_s : "SAMPLE";
my %hash;
my %mates;
my %single;
my $file = shift;
@ARGV = ("-") unless( $file );
my $IN;
if ( $file && $file =~ /\.bam$/i ) {
    open( $IN, "samtools view $file |" );
} else {
    $IN = \*ARGV;
}
while( <$IN> ) {
    next if ( /^@/ );
    my @a = split(/\t/);
    my $k = $1 if ( /NH:i:(\d+)/ );
    my $dir = $a[1] & 0x40 ? 1 : 2;
    $hash{ $k }->{ $dir }++;
    #$mates{ $k }++ unless( $a[6] eq "*" );
    if ( $a[1] & 0x8 ) {
        $single{ $k }->{ $dir }++;
    } else {
	$mates{ $k }->{ $dir }++;
    }
}
close( $IN );

my ($total1, $pairs1, $singletons1) = (0, 0, 0);
my ($total2, $pairs2, $singletons2) = (0, 0, 0);
while(my ($k, $v) = each %hash) {
    print join("\t", $sample, $k, $v->{1}, $v->{2}, ($v->{1} + $v->{2}), ($v->{1} + $v->{2})/$k, "All"), "\n";
    $total1 += $v->{1}/$k;
    $total2 += $v->{2}/$k;
}
while(my ($k, $v) = each %mates) {
    print join("\t", $sample, $k, $v->{1}, $v->{2}, ($v->{1} + $v->{2}), ($v->{1} + $v->{2})/$k, "Pairs"), "\n";
    $pairs1 += $v->{1}/$k;
    $pairs2 += $v->{2}/$k;
}
while(my ($k, $v) = each %single) {
    print join("\t", $sample, $k, $v->{1}, $v->{2}, ($v->{1} + $v->{2}), ($v->{1} + $v->{2})/$k, "Singletons"), "\n";
    $singletons1 += $v->{1}/$k;
    $singletons2 += $v->{2}/$k;
}
print "Total mapped:\t$total1\t$total2\t", $total1 + $total2, "\n";
print "Pairs mapped:\t$pairs1\t$pairs2\t", $pairs1 + $pairs2, "\n";
print "Singltons mapped:\t$singletons1\t$singletons2\t", $singletons1 + $singletons2, "\n";

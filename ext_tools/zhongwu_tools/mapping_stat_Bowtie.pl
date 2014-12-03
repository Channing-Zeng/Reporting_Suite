#!/usr/bin/perl -w

use Getopt::Std;
use strict;

our ($opt_s);
getopts('s:');

my $sample = $opt_s ? $opt_s : "SAMPLE";
my %hash;
my %mates;
my %single;
my $file = $ARGV[0];
@ARGV = ("-") unless( $file );
my $IN;
if ( $file && $file =~ /\.bam$/i ) {
    open( $IN, "samtools view $file |" );
} else {
    $IN = \*ARGV;
}
my ($total1, $pairs1, $singletons1) = (0, 0, 0);
my ($total2, $pairs2, $singletons2) = (0, 0, 0);
my %unmap = (1 => 0, 2 => 0);
while( <$IN> ) {
    next if ( /^@/ );
    my @a = split(/\t/);
    my $k = $1 if ( /XM:i:(\S+)/ );
    my $dir = $a[1] & 0x40 ? 1 : 2;
    if ( $a[1] & 0x4 ) {
        $unmap{ $dir }++;
        next;
    }
    $dir eq 1 ? $total1++ : $total2++;
    $hash{ $k }->{ $dir }++;
    #$mates{ $k }++ unless( $a[6] eq "*" );
    if ( $a[1] & 0x8 ) {
        $single{ $k }->{ $dir }++;
    } else {
	$mates{ $k }->{ $dir }++;
    }
}
close( $IN );

while(my ($k, $v) = each %hash) {
    my $v1 = $v->{ 1 } ? $v->{ 1 } : 0;
    my $v2 = $v->{ 2 } ? $v->{ 2 } : 0;
    print join("\t", $sample, $k, $v1, $v2, ($v1 + $v2), "All"), "\n";
}
while(my ($k, $v) = each %mates) {
    my $v1 = $v->{ 1 } ? $v->{ 1 } : 0;
    my $v2 = $v->{ 2 } ? $v->{ 2 } : 0;
    print join("\t", $sample, $k, $v1, $v2, ($v1 + $v2), "Pairs"), "\n";
    $pairs1 += $v1;
    $pairs2 += $v2;
}
while(my ($k, $v) = each %single) {
    my $v1 = $v->{ 1 } ? $v->{ 1 } : 0;
    my $v2 = $v->{ 2 } ? $v->{ 2 } : 0;
    print join("\t", $sample, $k, $v1, $v2, ($v1 + $v2), "Singletons"), "\n";
    $singletons1 += $v1;
    $singletons2 += $v2;
}
print "Total mapped:\t$total1\t$total2\t", $total1 + $total2, "\t$sample\n";
print "Pairs mapped:\t$pairs1\t$pairs2\t", $pairs1 + $pairs2, "\t$sample\n";
print "Singltons mapped:\t$singletons1\t$singletons2\t", $singletons1 + $singletons2, "\t$sample\n";
print "Unmapped:\t$unmap{ 1 }\t$unmap{2}\t", $unmap{ 1 } + $unmap{ 2 }, "\t$sample\n";

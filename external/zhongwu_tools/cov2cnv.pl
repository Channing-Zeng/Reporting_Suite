#!/usr/bin/perl -w

# Normalize the coverage from targeted sequencing to CNV log2 ratio.  The algorithm assumes the medium 
# is diploid, thus not suitable for homogeneous samples (e.g. parent-child).

use Stats::Basic;
use Getopt::Std;
use strict;

our ($opt_c);

getopts( 'c:' );

my $CNT = shift;
my %cnt;
open( CNT, $CNT );
my @cnt;
while(<CNT>) {
    chomp;
    next if ( /Undetermined/ );
    next if ( /Sample/ || /Total/ );
    my @a = split(/\t/);
    $cnt{$a[0]} = $a[1];
    push(@cnt, $a[1]);
}
close(CNT);
my $stat = new Stats::Basic;
my $meanreads = $stat->mean(\@cnt);
my %factor; # to adjust sequencing coverage
while(my ($k, $v) = each %cnt) {
    $factor{ $k } = $meanreads/$v;
}
my %raw; # key: gene value: hash of (key: sample; value: raw depth)
my %norm1; # key: gene value: hash of (key: sample; value: normalized depth by sequencing distribution)
my %norm1b; # key: gene value: hash of (key: sample; value: normalized depth by gene)
my %norm2; # key: gene value: hash of (key: sample; value: normalized by gene median scaling log2)
my %norm3; # key: gene value: hash of (key: sample; value: normalized by sample median scaling log2)
my %samples;
my @depth;
while( <> ) {
    next unless( /Whole-/ );
    next if ( /_CONTROL_/ );
    my ($sample, $gene, $chr, $s, $e, $t, $len, $depth) = split(/\t/);
    next if ( $sample eq "Undetermined" );
    my $k = join("\t", $gene, $chr, $s, $e, $len);
    $norm1{ $k }->{ $sample } = sprintf("%.3f", $depth*$factor{ $sample });
    $raw{ $k }->{ $sample } = $depth;
    $samples{ $sample } = 1;
    push(@depth, $depth*$factor{ $sample });
    #print join("\t", $sample, $k, $depth*$factor{ $sample }), "\n";
}
my $meddepth = $stat->median(\@depth);
my %factor2;  # Gene factor
while( my ($k, $v) = each %norm1) {
    my @t = values %$v;
    $factor2{ $k } = $stat->median(\@t) ? $meddepth/$stat->median(\@t) : 0;
}
my %samplemedian;
my @samples = keys %samples;
foreach my $s (@samples) {
    my @tmp = ();
    while( my ($k, $v) = each %norm1 ) {
        push( @tmp, $v->{ $s } );
    }
    $samplemedian{ $s } = $stat->median( \@tmp );
    #print STDERR "$s\t", $stat->median( \@tmp ), "\n";
}

while( my ($k, $v) = each %norm1) {
    foreach my $s (@samples) {
	    $norm1b{ $k }->{ $s } = $v->{$s} * $factor2{ $k }+0.1;
        print "$k : $v->{$s} : $factor2{ $k } : $norm1b{ $k }->{ $s }";
        $norm2{ $k }->{ $s } = sprintf("%.3f", log(($v->{$s} * $factor2{ $k }+0.1)/$meddepth)/log(2));
        $norm3{ $k }->{ $s } = sprintf("%.3f", log(($v->{$s} * $factor2{ $k }+0.1)/$samplemedian{ $s })/log(2));
	    #$v->{$s} = log($v->{$s}/$meddepth)/log(2);
    }
}

print join("\t", qw(Sample Gene Chr Start Stop Length MeanDepth MeanDepth_Norm1 MeanDepth_Norm2 log2Ratio_norm1 log2Ratio_norm2));
print "\tlog2Ratio_normContr" if ( $opt_c );
print "\n";
while( my ($k, $v) = each %norm2) {
    while( my ($s, $d) = each %$v ) {
        print join("\t", $s, $k, $raw{ $k }->{ $s }, $norm1{ $k }->{ $s }, $norm1b{ $k }->{ $s }, $d, $norm3{ $k }->{ $s });
        if ( $opt_c ) {
            my @controls = split(/:/, $opt_c);
            my @tcntl = map { $norm1b{ $k }->{ $_ } } @controls;
            my $cntl_ave = $stat->mean( \@tcntl );
            print "\t", sprintf("%.3f", log($norm1b{ $k }->{ $s }/$cntl_ave)/log(2));
            #print "\t", sprintf("%.3f", log($norm1b{ $k }->{ $s }/$norm1b{ $k }->{ $opt_c })/log(2));
        }
        print "\n";
    }
}

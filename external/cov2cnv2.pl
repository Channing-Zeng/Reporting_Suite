#!/usr/bin/perl -w

# Normalize the coverage from targeted sequencing to CNV log2 ratio.  The algorithm assumes the medium
# is diploid, thus not suitable for homogeneous samples (e.g. parent-child).

use Getopt::Std;
use strict;

our ($opt_c, $opt_a, $opt_i);

getopts( 'aic:' );

sub median {
    my @sorted = sort @_;
    @sorted[($#sorted / 2)];
}

sub mean {
    my @arr = @_;
    my $sum = 0;
    for ( @arr ) {
        $sum += $_;
    }
    $sum / $#arr;
}

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
my $meanreads = mean(@cnt);
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
my %data;
my %loc;
while( <> ) {
    next if ( /Whole-/ );
    next if ( /_CONTROL_/ );
    my ($sample, $gene, $chr, $s, $e, $t, $len, $depth) = split(/\t/);
    next if ( $sample eq "Undetermined" );
    my $k = $opt_a ? join("\t", $gene, $chr, $s, $e, $len) : $gene;
    $data{ $k }->{ $sample }->{ len } += $len;
    $loc{ $k }->{ chr } = $chr;
    $loc{ $k }->{ start } = $s unless( $loc{ $k }->{ start } && $loc{ $k }->{ start } < $s );
    $loc{ $k }->{ end } = $e unless( $loc{ $k }->{ end } && $loc{ $k }->{ end } > $e );
    $data{ $k }->{ $sample }->{ cov } += $len*$depth;
}

while(my($k, $v) = each %data) {
    while( my($sample, $dv) = each %$v ) {
        $raw{ $k }->{ $sample } = sprintf("%.2f", $dv->{cov}/$dv->{len});
        $norm1{ $k }->{ $sample } = sprintf("%.2f", $raw{ $k }->{ $sample }*$factor{ $sample });
        $samples{ $sample } = 1;
        push(@depth, $norm1{ $k }->{ $sample } );
        #print join("\t", $sample, $k, $depth*$factor{ $sample }), "\n";
    }
}

my $meddepth = median(\@depth);
my %factor2;  # Gene factor
while( my ($k, $v) = each %norm1) {
    my @t = values %$v;
    $factor2{ $k } = median(\@t) ? $meddepth/median(\@t) : 0;
}
my %samplemedian;
my @samples = keys %samples;
foreach my $s (@samples) {
    my @tmp = ();
    while( my ($k, $v) = each %norm1 ) {
        push( @tmp, $v->{ $s } );
    }
    $samplemedian{ $s } = median( \@tmp );
}

while( my ($k, $v) = each %norm1) {
    foreach my $s (@samples) {
        $norm1b{ $k }->{ $s } = sprintf("%.2f", $v->{$s} * $factor2{ $k }+0.1);
        $norm2{ $k }->{ $s } = sprintf("%.2f", log(($v->{$s} * $factor2{ $k }+0.1)/$meddepth)/log(2));
        $norm3{ $k }->{ $s } = sprintf("%.2f", log(($v->{$s} * $factor2{ $k }+0.1)/$samplemedian{ $s })/log(2));
    }
}

print join("\t", qw(Sample Gene Chr Start Stop Length MeanDepth MeanDepth_Norm1 MeanDepth_Norm2 log2Ratio_norm1 log2Ratio_norm2));
print "\tlog2Ratio_normContr" if ( $opt_c );
print "\n";
while( my ($k, $v) = each %norm2) {
    while( my ($s, $d) = each %$v ) {
        next if ( $opt_i && $factor2{ $k } == 0);
        my $t = $opt_a ? $k : "$k\t$loc{$k}->{chr}\t$loc{$k}->{start}\t$loc{$k}->{end}\t$data{$k}->{$s}->{len}";
        print join("\t", $s, $t, $raw{ $k }->{ $s }, $norm1{ $k }->{ $s }, $norm1b{ $k }->{ $s }, $d, $norm3{ $k }->{ $s });
        if ( $opt_c ) {
            my @controls = split(/:/, $opt_c);
            my @tcntl = map { $norm1b{ $k }->{ $_ } } @controls;
            my $cntl_ave = mean( \@tcntl );
            print "\t", sprintf("%.3f", log($norm1b{ $k }->{ $s }/$cntl_ave)/log(2));
            #print "\t", sprintf("%.3f", log($norm1b{ $k }->{ $s }/$norm1b{ $k }->{ $opt_c })/log(2));
        }
        print "\n";
    }
}

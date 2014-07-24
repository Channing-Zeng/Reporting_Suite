#!/usr/bin/perl -w

# Normalize the coverage from targeted sequencing to CNV log2 ratio.  The algorithm assumes the medium
# is diploid, thus not suitable for homogeneous samples (e.g. parent-child).

use Stat::Basic;
use Statistics::TTest;
use Getopt::Std;
use strict;

our ($opt_c, $opt_a, $opt_P, $opt_F);

getopts( 'aPc:F:' );

my $FAILEDFACTOR = defined($opt_F) ? $opt_F : 0.1;

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
my $stat = new Stat::Basic;
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
my %data;
my %loc;
while( <> ) {
    next if ( /Whole-/ );
    next if ( /_CONTROL_/ );
    next if ( /^Sample/ );
    my ($sample, $gene, $chr, $s, $e, $t, $len, $depth) = split(/\t/);
    next if ( $sample eq "Undetermined" );
    my $k = $opt_a ? join("\t", $gene, $chr, $s, $e, $len) : $gene;
    $loc{ $k }->{ chr } = $chr;
    $loc{ $k }->{ start } = $s unless( $data{ $k }->{ start } && $data{ $k }->{ start } < $s );
    $loc{ $k }->{ end } = $s unless( $data{ $k }->{ end } && $data{ $k }->{ end } > $e );
    $data{ $k }->{ $sample }->{ len } += $len;
    $data{ $k }->{ $sample }->{ cov } += $depth;
}

while(my($k, $v) = each %data) {
    while( my($sample, $dv) = each %$v ) {
        $raw{ $k }->{ $sample } = $dv->{cov};
        $norm1{ $k }->{ $sample } = sprintf("%.2f", $raw{ $k }->{ $sample }*$factor{ $sample });
        $samples{ $sample } = 1;
        push(@depth, $norm1{ $k }->{ $sample } );
    }
}

my @samples = keys %samples;
my $meddepth = $stat->median(\@depth);

# remove genes/amplicons that failed
my %bad;
my @gooddepth = ();
while(my($k, $v) = each %data) {
    my @tmp = map { $norm1{ $k }->{ $_} } @samples;
    my $kmax = $stat->max(\@tmp);
    if ( $kmax < $meddepth*$FAILEDFACTOR ) {
        $bad{ $k } = 1;
    } else {
        push(@gooddepth, @tmp);
    }
}
$meddepth = $stat->median(\@gooddepth); # re-adjust median depth using only those from good amplicons/genes

my %factor2;  # Gene/amplicon factor
while( my ($k, $v) = each %norm1) {
    my @t = values %$v;
    $factor2{ $k } = $stat->median(\@t) ? $meddepth/$stat->median(\@t) : 0;
}

my %samplemedian;
foreach my $s (@samples) {
    my @tmp = ();
    while( my ($k, $v) = each %norm1 ) {
        next if ( $bad{ $k } ); # ignore failed genes/amplicons
        push( @tmp, $v->{ $s } );
    }
    $samplemedian{ $s } = $stat->median( \@tmp );
    #print STDERR "$s\t", $stat->median( \@tmp ), "\n";
}

while( my ($k, $v) = each %norm1) {
    next if ( $bad{ $k } );
    foreach my $s (@samples) {
        $norm1b{ $k }->{ $s } = sprintf("%.2f", $v->{$s} * $factor2{ $k }+0.1);
        $norm2{ $k }->{ $s } = sprintf("%.2f", log(($v->{$s} * $factor2{ $k }+0.1)/$meddepth)/log(2));
        $norm3{ $k }->{ $s } = sprintf("%.2f", log(($v->{$s} * $factor2{ $k }+0.1)/$samplemedian{ $s })/log(2));
    }
}

#print $opt_a ? join("\t", qw(Sample Gene Chr Start Stop Length MeanDepth MeanDepth_Norm1 MeanDepth_Norm2 log2Ratio_norm1 log2Ratio_norm2)) : join("\t", qw(Sample Gene Length MeanDepth MeanDepth_Norm1 MeanDepth_Norm2 log2Ratio_norm1 log2Ratio_norm2));
print join("\t", qw(Sample Gene Chr Start Stop Length MeanDepth MeanDepth_Norm1 MeanDepth_Norm2 log2Ratio_norm1 log2Ratio_norm2));
print "\tlog2Ratio_normContr" if ( $opt_c );
print "\n";
my %g2amp;
while( my ($k, $v) = each %norm2) {
    next if ($bad{ $k });
    while( my ($s, $d) = each %$v ) {
        my $t = $opt_a ? $k : "$k\t$loc{$k}->{chr}\t$loc{$k}->{start}\t$loc{$k}->{end}\t$data{$k}->{$s}->{len}";
        my $str = join("\t", $s, $t, $raw{ $k }->{ $s }, $norm1{ $k }->{ $s }, $norm1b{ $k }->{ $s }, $d, $norm3{ $k }->{ $s });
        if ( $opt_c ) {
            my @controls = split(/:/, $opt_c);
            my @tcntl = map { $norm1b{ $k }->{ $_ } } @controls;
            my $cntl_ave = $stat->mean( \@tcntl );
            $str .= "\t" .  sprintf("%.3f", log($norm1b{ $k }->{ $s }/$cntl_ave)/log(2));
            #print "\t", sprintf("%.3f", log($norm1b{ $k }->{ $s }/$norm1b{ $k }->{ $opt_c })/log(2));
        }
        #my @a = split(/\t/, $str);
        #push(@{ $g2amp{ $s }->{ $a[1] } }, \@a);
        print "$str\n";
    }
}

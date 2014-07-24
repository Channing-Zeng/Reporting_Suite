#!/usr/bin/perl -w

# Normalize the coverage from targeted sequencing to CNV log2 ratio.  The algorithm assumes the medium
# is diploid, thus not suitable for homogeneous samples (e.g. parent-child).

use Stat::Basic;
use Statistics::TTest;
use Getopt::Std;
use strict;

our ($opt_c, $opt_a, $opt_s);

getopts( 'ac:s:' );

my $MINSEGS = $opt_s ? $opt_s : 1;
print join("\t", qw(Sample Gene Chr Start Stop Length MeanDepth MeanDepth_Norm1 MeanDepth_Norm2 log2Ratio_norm1 log2Ratio_norm2));
print "\tlog2Ratio_normContr" if ( $opt_c );
print "\n";
my %g2amp;
my $stat = new Stat::Basic;
my $ttest = new Statistics::TTest;
while( <> ) {
    s/\r//g;
    chomp;
    next if ( /^Sample/ );
    my @a = split(/\t/);
    my ($s, $g) = @a[0, 1];
    push(@{ $g2amp{ $s }->{ $a[1] } }, \@a);
}

while(my ($s, $v) = each %g2amp) {
    while(my ($g, $vv) = each %$v) {
        my @lr = map { $opt_c ? $_->[11] : $_->[10]; } @$vv;
        my $lr = @lr > 1 ? $stat->median(\@lr) : $lr[0];
        my ($sig, $bp, $type, $affected, $total, $siglr) = checkBP($vv);
        unless($sig ne "-1") {
            if ( $lr > 1 ) {
                ($sig, $bp, $type, $affected, $total, $siglr) = ("0", "Whole", "Amp", @lr+0, @lr+0, $lr);
            } elsif ( $lr < -1 ) {
                ($sig, $bp, $type, $affected, $total, $siglr) = ("0", "Whole", "Del", @lr+0, @lr+0, $lr);
            }
        }
        print join("\t", $s, $g, $lr, $sig, $bp, $type, $affected, $total, $siglr), "\n"; # if ( $sig ne "-1" );
    }
}

sub checkBP {
    my $ref = shift;
    return (-1, "", "", "", @$ref+0, "") if ( @$ref < 4 );
    my @a = map { $opt_c ? [$_->[3], $_->[11]] : [$_->[3], $_->[10]]; } @$ref;
    my @lr = map { $opt_c ? $_->[11] : $_->[10]; } @$ref;
    @a = sort { $a->[0] <=> $b->[0]; } @a;
    for(my $i = 0; $i < @a; $i++) {
        $a[$i]->[2] = $i+1;
    }
    my $max = $stat->max(\@lr);
    my $min = $stat->min(\@lr);
    my $mid = ($max + $min)/2;
    my @up = ();
    my @bm = ();
    my @lrup = ();
    my @lrbm = ();
    foreach my $a (@a) {
        $a->[1] > $mid ? push(@up, $a) : push(@bm, $a);
        $a->[1] > $mid ? push(@lrup, $a->[1]) : push(@lrbm, $a->[1]);
    }
    #print STDERR "UP: ", join("\t", (map { $_->[1]; } @up)), "\n";
    #print STDERR "BM: ", join("\t", (map { $_->[1]; } @bm)), "\n";
    my $lrupm = $stat->median(\@lrup);
    my $lrbmm = $stat->median(\@lrbm);
    my $cn = $lrbmm < -0.5 && abs($lrbmm) > abs($lrupm) ? "Del" : ($lrupm > 0.5 && abs($lrbmm) > abs($lrupm) ? "Amp" : "NA");
    if ( isConsecutive(\@up) && isConsecutive(\@bm) ) {
        my $sig = isSig(\@up, \@bm);
        return ($sig, "BP", $cn, $cn eq "Del" ? (@bm+0, @a+0, $lrbmm) : ($cn eq "Amp" ? (@up+0, @a+0, $lrupm) : (@bm+0 > @up+0 ? (@up+0, @a+0, $lrupm) : (@bm+0, @a+0, $lrbmm)))) if ($sig != -1);
    } elsif ( isConsecutive(\@up) ) {
        my $sig = isSig(\@up, \@bm);
        return ($sig, "BP", $cn, $cn eq "Del" ? (@bm+0, @a+0, $lrbmm) : ($cn eq "Amp" ? (@up+0, @a+0, $lrupm) : (@bm+0 > @up+0 ? (@up+0, @a+0, $lrupm) : (@bm+0, @a+0, $lrbmm)))) if ($sig != -1);
    } elsif ( isConsecutive(\@bm) ) {
        my $sig = isSig(\@bm, \@up);
        return ($sig, "BP", $cn, $cn eq "Del" ? (@bm+0, @a+0, $lrbmm) : ($cn eq "Amp" ? (@up+0, @a+0, $lrupm) : (@bm+0 > @up+0 ? (@up+0, @a+0, $lrupm) : (@bm+0, @a+0, $lrbmm)))) if ($sig != -1);
    }
    return ("-1", "", "", "", @a+0, "");
}

sub isSig {
    my ($a, $b) = @_;
    my @x = map { $_->[1]; } @$a;
    my @y = map { $_->[1]; } @$b;
    if (@$a >= 3 && @$b >= 3) {
        $ttest->load_data(\@x, \@y);
        my $p = $ttest->{ t_prob };
        #print STDERR "p: $p\n";
        return sprintf("%.5f", $p) if ( $p < 0.01 && $ttest->mean_difference() >= 0.75 );
    } elsif( @$a >= $MINSEGS && @$b >= 3 ) {
        my $med = $stat->median(\@y);
        my $mad = $stat->mad(\@y, 1);
        my @t = map { ($_->[1]-$med+0.1)/($mad+0.1); } @$a;
        my $mean = $stat->mean(\@t);
        return sprintf("%.2f", abs($mean)) if ( abs($mean) > 25 && abs($stat->mean(\@x)-$stat->mean(\@y)) > 1 );
    }
    return -1;  # Either too few to tell or not sig
}

sub isConsecutive {
    my $ref = shift;
    for(my $i = 1; $i < @$ref; $i++) {
        return 0 if ($ref->[$i]->[2] - $ref->[$i-1]->[2] > 1);
    }
    return 1;
}

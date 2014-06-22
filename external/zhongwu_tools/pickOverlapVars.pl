#!/usr/bin/perl -w

# Pick overlapping variants for a given bed file
use Getopt::Std;
use strict;

our ($opt_h, $opt_c, $opt_s, $opt_e, $opt_x, $opt_v, $opt_H);
getopts('hHvc:s:e:x:');
USAGE() if ( $opt_H );

my $col_c = $opt_c ? $opt_c - 1 : 4; # the column for chromosome in var file
my $col_s = $opt_s ? $opt_s - 1 : 5; # the column for start in var file
my $col_e = $opt_e ? $opt_e - 1 : 6; # the column for end in var file

my $bed = shift;
my %bed;
open( BED, $bed );
while( <BED> ) {
    chomp;
    next if (/^@/);
    next if (/^#/);
    my @a = split(/\t/);
    next if ( /^browser/ || /^track/ );
    $a[0] = "chr$a[0]" unless( $a[0] =~ /^chr/ );
    my ($s, $e) = @a[1,2];
    ($s, $e) = ($s-$opt_x, $e+$opt_x) if ( $opt_x );
    push( @{ $bed{ $a[0] } }, [$s, $e] );
}
close( BED );

my $N = 1;
while( <> ) {
    if ( $N && $opt_h ) {
        print;
	$N = 0;
	next;
    }
    my @a = split(/\t/);
    my ($chr, $s, $e) = @a[$col_c, $col_s, $col_e];
    $chr = "chr$chr" unless( $chr =~ /^chr/ );
    $a[$col_c] = $chr;
    if (checkOverlap( $chr, $s, $e)) {
	print join("\t", @a) unless( $opt_v ); 
    } else {
	print join("\t", @a) if( $opt_v ); 
    }
}

sub checkOverlap {
    my ($chr, $s, $e) = @_;
    foreach my $r (@{ $bed{ $chr } }) {
        return 1 if ( $e >= $r->[0] && $s <= $r->[1] );
    }
    return 0;
}

sub USAGE {
    print <<USAGE;
    Usage: $0 [-hv] [-c chr_col] [-s start_col] [-e end_col] bed varfile

    The program will pick entries in varfile whose genomic locations are within bed file.

    Files:
    bed	  BED file
    varfile  Any file containing genomic locations

    Options:
    -h Indicate there's a header and should be printed
    -v Print non-overlapping entries instead.
    -c chr_col
        The column containing chromosome
    -s start_col
        The column containing start position
    -e end_col
        The column containing end position
    -x num_bp
        Number of basepair to extend for BED file.  For capture based technology, sequence
	outside of target region can also be captured.
USAGE
    exit(0);
}

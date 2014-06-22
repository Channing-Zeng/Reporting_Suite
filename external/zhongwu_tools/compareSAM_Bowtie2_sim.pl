#!/usr/bin/perl -w

# To compare to SAM files from aligned using Bowtie2
# The score is determined by the sum of alignment scores of both pairs

use Getopt::Std;
use Stat::Basic;
use strict;

our ($opt_s);
getopts('s:');
my $suffix = $opt_s ? $opt_s : "";

my $bam1 = shift;
my $bam2 = shift;

my $comd1 = $bam1;
$comd1 = "samtools view $bam1 | " if ( $bam1 =~ /bam$/ );
open( SAM1, $comd1);
my %sam1;
keys(%sam1) = 6000000;
my $n = 0;
my $np = 0;
print STDERR "Step 1: ", `date`;
my $stat = new Stat::Basic;
my $OW1 = 0;
while( <SAM1> ) {
    my @a = split(/\t/);
    $n++;
    #my $dir = $a[1] & 0x40 ? "1" : "2";
    my $read = "$a[0]";
    my $score = $1 if ( /AS:i:(\d+)/ );
    $score += $1 if ( /YS:i:(\d+)/ );
    #my $mate = $a[1] & 0x8 ? 0 : 1; # Mate is unmapped
    # overwrite if new mapping is better
    if ( ! $sam1{ $read } ) {
	$sam1{ $read } = $score;
	$np++;
    } elsif ( $sam1{ $read } < $score ) {
        $sam1{ $read } = $score;
	$OW1++;
    }
    #push( @{ $sam1{ $read } }, { S => $score, Q => $a[4], M => $mate } ); #, S => $_ } );
    print STDERR "$n mapping processed in sample 1\n" if ( $n%1000000 == 0 );
}
print STDERR "$n total mapping processed in sample 1\n";
print STDERR "$np total fragments mapped in sample 1\n";
print STDERR "$OW1 pairs got overwritten by better score in sample 1\n";
close( SAM1 );

open( SAM1R, ">reads_sample1$suffix.txt" );
open( SAM2R, ">reads_sample2$suffix.txt" );
open( AMBR, ">reads_ambigous$suffix.txt" );

print STDERR "Step 2: ", `date`;
my $comd2 = $bam2;
$comd2 = "samtools view $bam2 | " if ( $bam2 =~ /bam$/ );
open( SAM2, $comd2 );
my %sam2;
my %sam2u;
keys( %sam2 ) = 6000000;
$n = 0;
$np = 0;
my $N = 0;
my $OW2 = 0;
while( <SAM2> ) {
    my @a = split(/\t/);
    my $read = "$a[0]";
    $N++;
    print STDERR "$N total mapping processed in sample 2\n" if ( $N%1000000 == 0 );
    unless( $sam1{ $read } ) {
	unless( $sam2u{ $read } ) {
	    print SAM2R "$read\tSample2 only\n";
	    $np++;
	}
	$sam2u{ $read } = 1;
	next;
    }
    $n++;
    #my $dir = $a[1] & 0x40 ? "1" : "2";
    my $score = $1 if ( /AS:i:(\d+)/ );
    $score += $1 if ( /YS:i:(\d+)/ );
    if ( ! $sam2{ $read } ) {
        $sam2{ $read } = $score;
	$np++;
    } elsif ( $sam2{ $read } < $score ) {
        $sam2{ $read } = $score;
	$OW2++;
    }
    #push( @{ $sam2{ $read } }, { S => $score, Q => $a[4], M => $mate } ); #, S => $_ } );
    print STDERR "$n overalpping mapping processed in sample 2\n" if ( $n%1000000 == 0 );
}
print STDERR "$n total mapping for overlapping reads processed in sample 2\n";
print STDERR "$N total mapping in sample 2\n";
print STDERR "$np total fragments mapped in sample 2\n";
print STDERR "$OW2 pairs got overwritten by better score in sample 2\n";

print STDERR "Step 3: ", `date`;
#my $resolved_by_better_mate = 0;
#my $resolved_by_mate = 0;
#my $resolved_by_quality = 0;
#my $resolved_conflict = 0;
#my $unresolved = 0;
#my $unresolved_conflict = 0;
#my $resolved_2_by_better_mate = 0;
#my $unresolved_with_same_mate = 0;
#my $resolved_2_uniq_mate = 0;
#my $resolved_2_conflict = 0;
#my $better2 = 0;
my $new = 0;
my %hash;
while( my ($read, $score1) = each %sam1 ) {
    unless( $sam2{ $read } ) {
        $new++;
	print SAM1R "$read\tSample1_Only\n";
	#print STDERR "New S2R:\t", $_->{ S } foreach(@$map);
	next;
    }
    if ( $score1 > $sam2{ $read } ) {
        print SAM1R "$read\tBetter\n";
    } elsif ( $score1 < $sam2{ $read } ) {
        print SAM2R "$read\tBetter\n";
    } else {
        print AMBR "$read\tEqual\n";
    }
}

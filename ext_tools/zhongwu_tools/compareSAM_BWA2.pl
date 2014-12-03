#!/usr/bin/perl -w

# To compare to SAM files from explant sequencing

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
print STDERR "Step 1: ", `date`;
my $stat = new Stat::Basic;
while( <SAM1> ) {
    my @a = split(/\t/);
    $n++;
    my $dir = $a[1] & 0x40 ? "1" : "2";
    my $read = "$a[0]";
    my $nm = $1 if ( /NM:i:(\d+)/ );
    my $nh = /X0:i:(\d+)/ ? $1 : 1;
    my @S = $a[5] =~ m/(\d+)S/g;
    my $S = @S > 0 ? $stat->sum(\@S) : 0;
    my $gp = /XO:i:(\d+)/ ? $1 : 0;
    my $score = $nm + $nh + $S + $gp;
    my $mate = $a[1] & 0x8 ? 0 : 1; # Mate is unmapped
    $sam1{ $read }->{$dir} = { S => -$score, Q => $a[4], M => $mate } unless( $sam1{ $read }->{$dir} );
    # overwrite if new mapping is better
    if ( $sam1{ $read }->{$dir} && $sam1{ $read }->{$dir}->{ S } > -$score ) {
	#print STDERR "Better in 1: $_";
        $sam1{ $read }->{$dir} = { S => -$score, Q => $a[4], M => $mate };
    }
    #push( @{ $sam1{ $read } }, { S => $score, Q => $a[4], M => $mate } ); #, S => $_ } );
    print STDERR "$n mapping processed in sample 1\n" if ( $n%1000000 == 0 );
}
print STDERR "$n total mapping processed in sample 1\n";
close( SAM1 );

print STDERR "Step 2: ", `date`;
my $comd2 = $bam2;
$comd2 = "samtools view $bam2 | " if ( $bam2 =~ /bam$/ );
open( SAM2, $comd2 );
my %sam2;
keys( %sam2 ) = 6000000;
$n = 0;
my $N = 0;
while( <SAM2> ) {
    my @a = split(/\t/);
    my $read = "$a[0]";
    $N++;
    print STDERR "$N total mapping processed in sample 2\n" if ( $N%1000000 == 0 );
    next unless( $sam1{ $read } );
    $n++;
    my $dir = $a[1] & 0x40 ? "1" : "2";
    my $nm = $1 if ( /NM:i:(\d+)/ );
    my $nh = /X0:i:(\d+)/ ? $1 : 1;
    my @S = $a[5] =~ m/(\d+)S/g;
    my $S = @S > 0 ? $stat->sum(\@S) : 0;
    my $gp = /XO:i:(\d+)/ ? $1 : 0;
    my $score = $nm + $nh + $S + $gp;
    my $mate = $a[1] & 0x8 ? 0 : 1; # Mate is unmapped
    $sam2{ $read }->{ $dir } = { S => -$score, Q => $a[4], M => $mate } unless( $sam2{ $read }->{ $dir } );
    # overwrite if new mapping is better
    if ( $sam2{ $read }->{ $dir } && $sam2{ $read }->{ $dir }->{ S } > -$score ) {
	#print STDERR "Better in 2: $_";
        $sam2{ $read }->{ $dir } = { S => -$score, Q => $a[4], M => $mate };
    }
    #push( @{ $sam2{ $read } }, { S => $score, Q => $a[4], M => $mate } ); #, S => $_ } );
    print STDERR "$n overalpping mapping processed in sample 2\n" if ( $n%1000000 == 0 );
}
print STDERR "$n total mapping for overlapping reads processed in sample 2\n";
print STDERR "$N total mapping in sample 2\n";

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
open( SAM1R, ">reads_sample1$suffix.txt" );
open( SAM2R, ">reads_sample2$suffix.txt" );
open( AMBR, ">reads_ambigous$suffix.txt" );
my %hash;
while( my ($read, $map) = each %sam1 ) {
    unless( $sam2{ $read } ) {
        $new++;
	print SAM1R "$read\tSample1_Only\n";
	#print STDERR "New S2R:\t", $_->{ S } foreach(@$map);
	next;
    }
    my ($s1a, $s1b) = ($sam1{ $read }->{ 1 } ? $sam1{ $read }->{ 1 }->{ S } : undef, $sam1{ $read }->{ 2 } ? $sam1{ $read }->{ 2 }->{ S } : undef);
    my ($s2a, $s2b) = ($sam2{ $read }->{ 1 } ? $sam2{ $read }->{ 1 }->{ S } : undef, $sam2{ $read }->{ 2 } ? $sam2{ $read }->{ 2 }->{ S } : undef);
    my @info = ();
    push( @info, [$s1a, 1] ) if ( $s1a );
    push( @info, [$s1b, 1] ) if ( $s1b );
    push( @info, [$s2a, 2] ) if ( $s2a );
    push( @info, [$s2b, 2] ) if ( $s2b );
    @info = sort { $b->[0] <=> $a->[0] } @info;
    my $to = 0;
    my $cat = "";
    if ( @info == 4 ) {
	my $flag = 0;
        $flag = 1 if ( ($s1a > $s2a) && ($s1b < $s2b) );
        $flag = 1 if ( ($s1a < $s2a) && ($s1b > $s2b) );
	if ( $flag ) {
	    if ( $info[0]->[0] == $info[1]->[0] && $info[0]->[1] != $info[1]->[1] ) {
	        if ( $s1a + $s1b > $s2a + $s2b) {
		    $to = 1;
	        } elsif ( $s1a + $s1b < $s2a + $s2b) {
		    $to = 2;
		}
	    } else {
	        $to = $info[0]->[0];
	    }
	    $cat = $to ? "Resolved_Conflict" : "Unresolved_Conflict";
	} else {
	    if ( $info[0]->[0] == $info[1]->[0] && $info[0]->[1] != $info[1]->[1] ) {
		if ( $info[2]->[0] > $info[3]->[0] ) {
		    $to = $info[2]->[1];
		}
		$cat = $to ? "Better_Mapped_Mate" : "Unresolved_Both_Pair_Same_Score";
	    } elsif ( $info[0]->[0] == $info[1]->[0] && $info[0]->[1] == $info[1]->[1] ) {
	        $to = $info[0]->[1];
		$cat = "Better_Score_Both_Pairs";
	    } elsif ( $info[0]->[0] > $info[1]->[0] ) {
	        $to = $info[0]->[1];
		$cat = "Better_Score_Single_End";
	    } else {
		print STDERR "Should never happen3.\n";
	    }
	}
    } elsif ( @info == 3 ) {
	if ( $info[0]->[0] > $info[1]->[0] ) {
	    if ( $info[0]->[1] == $info[1]->[1] ) {
		$to = $info[0]->[1];
		$cat = "Better_and_unique_mate";
	    } elsif ( $info[0]->[1] == $info[2]->[1] ) {
		$to = $info[0]->[1];
		$cat = "Better_and_unique_mate";
	    } else {
		$to = $info[0]->[1];
		$cat = "Better_but_less_mate";
	    }
	} else {
	    $to = $info[2]->[1];
	    $cat = "Unique_Mapped_Mate";
	}
    } elsif( @info == 2 ) {
	if ( $info[0]->[0] > $info[1]->[0] ) {
	    $to = $info[0]->[1];
	    $cat = "Better_Score_No_Mate";
	} else {
	    $cat = "Unresolved_Same_Score_No_Mate";
	}
    } else { 
        ($to, $cat) = ($info[0]->[1], "Only_Mapping");
    }
    if ( $to == 1 ) {
        print SAM1R "$read\t$cat\n";
	$hash{ $to }->{ $cat }++;
    } elsif( $to == 2 ) {
        print SAM2R "$read\t$cat\n";
	$hash{ $to }->{ $cat }++;
    } else {
        print AMBR "$read\t$cat\n";
	$hash{ UN }->{ $cat }++;
    }
}
close( SAM2 );
close( SAM1R );
close( SAM2R );
close( AMBR );

#print "Resolved to sample 1 better mapping quality in sample 1:\t$resolved_by_quality\n";
#print "Resolved to sample 1 by better mate in sample 1:\t$resolved_by_better_mate\n";
#print "Resovled to sample 1 by sample 1 only mate:\t$resolved_by_mate\n";
print "Resolved to sample 1 by sample 1 only mapping:\t$new\n";
#print "Resolved to sample 1 by sample 1 best conflict read:\t$resolved_conflict\n";
#print "RESOLVED to sample 2 by better mapping quality in sample 2:\t$better2\n";
#print "RESOLVED to sample 2 by sample 2 only mate:\t$resolved_2_uniq_mate\n";
#print "RESOLVED to sample 2 by better mate in sample 2:\t$resolved_2_by_better_mate\n";
#print "RESOLVED to sample 2 by sample 2 best conflict read:\t$resolved_2_conflict\n";
#print "Unresolved (both reads have the same mapping quality):\t$unresolved_with_same_mate\n";
#print "Unresolved (mate is not mapped in both samples):\t$unresolved\n";
#print "Unresolved (conflicts):\t$unresolved_conflict\n";
foreach my $d (1, 2, "UN") {
    while( my ($k, $v) = each %{ $hash{ $d } } ) {
        print "$d $k:\t$v\n";
    }
}
print STDERR "Done: ", `date`, "\n";

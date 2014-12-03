#!/usr/bin/perl -w

use Stat::Basic;
# To compare to SAM files from explant sequencing

use Getopt::Std;
use strict;

our ($opt_s);
getopts('s:');
my $suffix = $opt_s ? $opt_s : "";

my $bam1 = shift;
my $bam2 = shift;

my $RLEN = 90; # the read length
my $stat = new Stat::Basic;

my $comd1 = $bam1;
$comd1 = "samtools view $bam1 | " if ( $bam1 =~ /bam$/ );
open( SAM1, $comd1);
my %sam1;
keys(%sam1) = 6000000;
my $n = 0;
print STDERR "Step 1: ", `date`;
while( <SAM1> ) {
    my @a = split(/\t/);
    $n++;
    my $dir = $a[1] & 0x40 ? "1" : "2";
    my $read = "$a[0]$dir";
    my $nm = $1 if ( /NM:i:(\d+)/ );
    my $nh = /X0:i:(\d+)/ ? $1 : 1;
    my @S = $a[5] =~ m/(\d+)S/g;
    my $S = @S > 0 ? $stat->sum(\@S) : 0;
    my $gp = /XO:i:(\d+)/ ? $1 : 0;  # Additional penalty for gap opens
    my $score = $nm + $nh + $gp + $S;
    my $mate = $a[1] & 0x8 ? 0 : 1; # Mate is unmapped
    $sam1{ $read } = { S => -$score, Q => $a[4], M => $mate } unless( $sam1{ $read } );
    # overwrite if new mapping is better
    if ( $sam1{ $read } && $sam1{ $read }->{ S } > -$score ) {
	#print STDERR "Better in 1: $_";
        $sam1{ $read } = { S => -$score, Q => $a[4], M => $mate };
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
    my $t1 = "$a[0]1";
    my $t2 = "$a[0]2";
    $N++;
    print STDERR "$N total mapping processed in sample 2\n" if ( $N%1000000 == 0 );
    next unless( $sam1{ $t1 } || $sam1{ $t2 } );
    $n++;
    my $dir = $a[1] & 0x40 ? "1" : "2";
    my $read = "$a[0]$dir";
    my $nm = $1 if ( /NM:i:(\d+)/ );
    my $nh = /X0:i:(\d+)/ ? $1 : 1;
    my @S = $a[5] =~ m/(\d+)S/g;
    my $S = @S > 0 ? $stat->sum(\@S) : 0;
    my $gp = /XO:i:(\d+)/ ? $1 : 0;  # Additional penalty for gap opens
    my $score = $nm + $nh + $gp + $S;
    my $mate = $a[1] & 0x8 ? 0 : 1; # Mate is unmapped
    $sam2{ $read } = { S => -$score, Q => $a[4], M => $mate } unless( $sam2{ $read } );
    # overwrite if new mapping is better
    if ( $sam2{ $read } && $sam2{ $read }->{ S } > -$score ) {
	#print STDERR "Better in 2: $_";
        $sam2{ $read } = { S => -$score, Q => $a[4], M => $mate };
    }
    #push( @{ $sam2{ $read } }, { S => $score, Q => $a[4], M => $mate } ); #, S => $_ } );
    print STDERR "$n overalpping mapping processed in sample 2\n" if ( $n%1000000 == 0 );
}
print STDERR "$n total mapping for overlapping reads processed in sample 2\n";
print STDERR "$N total mapping in sample 2\n";

print STDERR "Step 3: ", `date`;
my $resolved_by_better_mate = 0;
my $resolved_by_mate = 0;
my $resolved_by_quality = 0;
my $resolved_conflict = 0;
my $unresolved = 0;
my $unresolved_conflict = 0;
my $resolved_2_by_better_mate = 0;
my $unresolved_with_same_mate = 0;
my $resolved_2_uniq_mate = 0;
my $resolved_2_conflict = 0;
my $better2 = 0;
my $new = 0;
open( SAM1R, ">reads_sample1$suffix.txt" );
open( SAM2R, ">reads_sample2$suffix.txt" );
open( AMBR, ">reads_ambigous$suffix.txt" );
while( my ($read, $map) = each %sam1 ) {
    my $READ = $read; $READ =~ s/.$//;  # The original read name
    my $m = $read;
    $read =~ /1$/ ? ($m =~ s/1$/2/) : ($m =~ s/2$/1/);
    unless( $sam2{ $read } ) {
        $new++;
	print SAM1R "$READ\tSample1_Only\n";
	#print STDERR "New S2R:\t", $_->{ S } foreach(@$map);
	next;
    }
    my ($s1a, $s1b) = ($sam1{ $read }->{ S }, $sam1{ $m } ? $sam1{ $m }->{ S } : undef);
    my ($s2a, $s2b) = ($sam2{ $read }->{ S }, $sam2{ $m } ? $sam2{ $m }->{ S } : undef);
    if ( $s1a && $s1b && $s2a && $s2b ) {
	my $flag = 0;
        $flag = 1 if ( ($s1a > $s2a) && ($s1b < $s2b) );
        $flag = 1 if ( ($s1a < $s2a) && ($s1b > $s2b) );
	if ( $flag ) {
	    #print STDERR join("', '", $READ, $s1a, $s1b, $s2a, $s2b), "\n";
	    my @tmp = sort { $b->[0] <=> $a->[0] } ([$s1a, 1], [$s1b, 1], [$s2a, 2], [$s2b, 2]);
	    my $to = 0;
	    if ( $tmp[0]->[0] == $tmp[1]->[0] && $tmp[0]->[1] != $tmp[1]->[1] ) {
	        if ( $s1a + $s1b > $s2a + $s2b) {
		    $to = 1;
	        } elsif ( $s1a + $s1b < $s2a + $s2b) {
		    $to = 2;
		}
	    } elsif ( $tmp[0]->[1] == 1 ) {
		$to = 1;
	    } else {
	        $to = 2;
	    }
	    if ( $to == 1 ) {
		$resolved_conflict++;
	        print SAM1R "$READ\tConflict but 1\n";
		#print STDERR "$read has conflicting pair scores but to 1.\n" if ( $flag );
	    } elsif ( $to == 2 ) {
		$resolved_2_conflict++;
	        print SAM2R "$READ\tConflict but 2\n";
		#print STDERR "$read has conflicting pair scores but to 2.\n" if ( $flag );
	    } else {
		$unresolved_conflict++;
		print AMBR "$READ\tUnresolved conflicts\n";
		#print STDERR "$read is unresolved conflicting pair scores but to 2.\n" if ( $flag );
	    }
	    next;
	}
    }
    if ( $s1a == $s2a ) {
        if ( $sam2{ $m } && $sam1{ $m } ) {  # the mate is also mapped
	    #my $msum1 = $sam1{ $m }->[0]->{ NM } + $sam1{ $m }->[0]->{ NH };
	    #my $msum2 = $sam2{ $m }->[0]->{ NM } + $sam2{ $m }->[0]->{ NH };
	    if ( $s1b > $s2b ) {
	        $resolved_by_better_mate++;
		print SAM1R "$READ\tBetter_Mapped_Mate\n";
	    } elsif ( $s1b < $s2b ) {
	        $resolved_2_by_better_mate++;
		print SAM2R "$READ\tBetter_Mapped_Mate\n";
	    } else {
	        $unresolved_with_same_mate++;
		print AMBR "$READ\tBoth_Pair_Same_Score\n";
		#print STDERR "Unresolved_with_mate: $read\n";
		#print STDERR "UN Mate S2R:\t", $_->{ S } foreach(@$map);
		#print STDERR "UN Mate S2M:\t", $_->{ S } foreach(@{ $sam2{ $m } });
		#print STDERR "UN Mate S1R:\t", $_->{ S } foreach(@{ $sam1{ $read } });
		#print STDERR "UN Mate S1M:\t", $_->{ S } foreach(@{ $sam1{ $m } });
	    }
	} elsif ( $sam2{ $m } ) { # mate is mapped in SAM2, but not SAM1
	    $resolved_2_uniq_mate++;
	    print SAM2R "$READ\tUnique_Mapped_Mate\n";
	    #print STDERR "Unresolved_with_mate_sam2: $read\n";
	    #print STDERR "UN Mate sam2 S2R:\t", $_->{ S } foreach(@$map);
	    #print STDERR "UN Mate sam2 S2M:\t", $_->{ S } foreach(@{ $sam2{ $m } });
	    #print STDERR "UN Mate sam2 S1R:\t", $_->{ S } foreach(@{ $sam1{ $read } });
	} elsif( $sam1{ $m } ) {  # mate is mapped in SAM1, but not SAM2
	    $resolved_by_mate++;  # so it belong to SAM1
	    print SAM1R "$READ\tUnique_Mapped_Mate\n";
	} else { # Mate is mapped in neither SAM1 nor SAM2
	    $unresolved++;
	    print AMBR "$READ\tUnmapped_mate\n";
	    #print STDERR "Unresolved: $read\n";
	    #print STDERR "UN S2R:\t", $_->{ S } foreach(@$map);
	    #print STDERR "UN S1R:\t", $_->{ S } foreach(@{ $sam1{ $read } });
	}
    } elsif ( $s1a > $s2a ) {
        $resolved_by_quality++; # So it belong to SAM1
	print SAM1R "$READ\tUnique_Mapped_Mate\n";
    } else {
        $better2++;
	print SAM2R "$READ\tBetter_Score\n";
	#print STDERR "Better2: $read\n";
	#print STDERR "Better2 S2R:\t", $_->{ S } foreach(@$map);
	#print STDERR "Better2 S1R:\t", $_->{ S } foreach(@{ $sam1{ $read } });
    }
}
close( SAM2 );
close( SAM1R );
close( SAM2R );
close( AMBR );

print "Resolved to sample 1 by sample 1 only mapping:\t$new\n";
print "Resolved to sample 1 better mapping quality in sample 1:\t$resolved_by_quality\n";
print "Resolved to sample 1 by better mate in sample 1:\t$resolved_by_better_mate\n";
print "Resovled to sample 1 by sample 1 only mate:\t$resolved_by_mate\n";
print "Resolved to sample 1 by sample 1 best conflict read:\t$resolved_conflict\n";
print "RESOLVED to sample 2 by better mapping quality in sample 2:\t$better2\n";
print "RESOLVED to sample 2 by better mate in sample 2:\t$resolved_2_by_better_mate\n";
print "RESOLVED to sample 2 by sample 2 only mate:\t$resolved_2_uniq_mate\n";
print "RESOLVED to sample 2 by sample 2 best conflict read:\t$resolved_2_conflict\n";
print "Unresolved (both reads have the same mapping quality):\t$unresolved_with_same_mate\n";
print "Unresolved (mate is not mapped in both samples):\t$unresolved\n";
print "Unresolved (conflicts):\t$unresolved_conflict\n";
print STDERR "Done: ", `date`, "\n";

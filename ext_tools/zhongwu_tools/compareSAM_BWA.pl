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
while( <SAM1> ) {
    my @a = split(/\t/);
    next if ( $a[1] & 0x4 );
    $n++;
    my $dir = $a[1] & 0x40 ? "1" : "2";
    #$a[0] =~ s/^[^:]+://;
    my $read = "$a[0]$dir";
    #print STDERR "$read\n";
    my $nm = $1 if ( /NM:i:(\d+)/ ); # Edit distance
    my $nh = /X0:i:(\d+)/ ? $1 : 1;
    my @S = $a[5] =~ m/(\d+)S/g;
    my $S = @S > 0 ? $stat->sum(\@S) : 0;
    my $mate = $a[1] & 0x8 ? 0 : 1;
    push( @{ $sam1{ $read } }, { NM => $nm, NH => $nh, s => $S, Q => $a[4], M => $mate } ); #, S => $_ } );
    print STDERR "$n reads processed in sample 1\n" if ( $n%1000000 == 0 );
}
close( SAM1 );

my $comd2 = $bam2;
$comd2 = "samtools view $bam2 | " if ( $bam2 =~ /bam$/ );
open( SAM2, $comd2 );
my %sam2;
keys( %sam2 ) = 6000000;
$n = 0;
while( <SAM2> ) {
    my @a = split(/\t/);
    next if ( $a[1] & 0x4 );
    my $dir = $a[1] & 0x40 ? "1" : "2";
    #$a[0] =~ s/^[^:]+://;
    my $t1 = "$a[0]1";
    my $t2 = "$a[0]2";
    next unless( $sam1{ $t1 } || $sam1{ $t2 } );
    $n++;
    my $read = "$a[0]$dir";
    my $nm = $1 if ( /NM:i:(\d+)/ );
    my $nh = /X0:i:(\d+)/ ? $1 : 1;
    my @S = $a[5] =~ m/(\d+)S/g;
    my $S = @S > 0 ? $stat->sum(\@S) : 0;
    my $mate = $a[1] & 0x8 ? 0 : 1;
    push( @{ $sam2{ $read } }, { NM => $nm, NH => $nh, s => $S, Q => $a[4], M => $mate } ); #, S => $_ } );
    print STDERR "$n reads processed in sample 2\n" if ( $n%1000000 == 0 );
}
close( SAM2 );

print STDERR "Step 3: ", time(), "\n";
my $resolved_by_better_mate = 0;
my $resolved_by_mate = 0;
my $resolved_by_quality = 0;
my $unresolved = 0;
my $resolved_2_by_better_mate = 0;
my $unresolved_with_same_mate = 0;
my $resolved_2_uniq_mate = 0;
my $better2 = 0;
my $new = 0;
open( SAM1R, ">reads_sample1.txt" );
open( SAM2R, ">reads_sample2.txt" );
open( AMBR, ">reads_ambigous.txt" );
while( my ($read, $map) = each %sam2 ) {
    my $READ = $read; $READ =~ s/.$//;  # The original read name
    unless( $sam1{ $read } ) { # Only mapped in sample 2
        $new++;
	print SAM2R "$READ\tSample2_Only\n";
	#print STDERR "New S2R:\t", $_->{ S } foreach(@$map);
	next;
    }
    my $m = $read;
    $read =~ /1$/ ? ($m =~ s/1$/2/) : ($m =~ s/2$/1/);
    my $sum1 = $sam1{ $read }->[0]->{ NM } + $sam1{ $read }->[0]->{ NH } + $sam1{ $read }->[0]->{ s };
    my $sum2 = $map->[0]->{ NM } + $map->[0]->{ NH } + $map->[0]->{ s };
    #my $sum1 = $sam1{ $read }->[0]->{ Q }; # the mapping quality
    #my $sum2 = $map->[0]->{ Q };
    if ( $sum1 == $sum2 ) {
        if ( $sam2{ $m } && $sam1{ $m } ) {  # the mate is also mapped
	    #my $msum1 = $sam1{ $m }->[0]->{ Q };
	    #my $msum2 = $sam2{ $m }->[0]->{ Q };
	    my $msum1 = $sam1{ $m }->[0]->{ NM } + $sam1{ $m }->[0]->{ NH } + $sam1{ $m }->[0]->{ s };
	    my $msum2 = $sam2{ $m }->[0]->{ NM } + $sam2{ $m }->[0]->{ NH } + $sam2{ $m }->[0]->{ s };
	    if ( $msum2 > $msum1 ) {
	        $resolved_by_better_mate++;
		print SAM1R "$READ\tBetter_Mapped_Mate\n";
	    } elsif ( $msum2 < $msum1 ) {
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
    } elsif ( $sum1 < $sum2 ) {
        $resolved_by_quality++;
	print SAM1R "$READ\tBetter_Score\n";
    } else {
        $better2++;
	print SAM2R "$READ\tBetter_Score\n";
	#print STDERR "Better2: $read\n";
	#print STDERR "Better2 S2R:\t", $_->{ S } foreach(@$map);
	#print STDERR "Better2 S1R:\t", $_->{ S } foreach(@{ $sam1{ $read } });
    }
}
close( SAM1R );
close( SAM2R );
close( AMBR );

print "Resolved to sample 1 better mapping quality in sample 1:\t$resolved_by_quality\n";
print "Resolved to sample 1 by better mate in sample 1:\t$resolved_by_better_mate\n";
print "Resovled to sample 1 by sample 1 only mate:\t$resolved_by_mate\n";
print "RESOLVED to sample 2 by sample 2 only mapping:\t$new\n";
print "RESOLVED to sample 2 by better mapping quality in sample 2:\t$better2\n";
print "RESOLVED to sample 2 by sample 2 only mate:\t$resolved_2_uniq_mate\n";
print "RESOLVED to sample 2 by better mate in sample 2:\t$resolved_2_by_better_mate\n";
print "Unresolved (both reads have the same mapping quality):\t$unresolved_with_same_mate\n";
print "Unresolved (mate is not mapped in both samples):\t$unresolved\n";

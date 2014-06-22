#!/usr/bin/perl -w

# The script will take a mutation list, and output from VarDict.pl and -D option
# and use the mutation genotype to extract info from the output.
use Getopt::Std;
use strict;

our ($opt_f, $opt_d);
getopts('f:d:');

$opt_f = $opt_f ? $opt_f : 0.001;
my $MINDP = defined($opt_d) ? $opt_d : 5; # the minimum depeth

my %mut;
my $MUT = shift;
open(MUT, $MUT);
#chr7    55249063        G       A       EGFR    Q787Q
#chr7    55242464        AGGAATTAAGAGAAGC        A       EGFR    KELREA745K
while(<MUT>) {
    chomp;
    next if ( /^#/ );
    my @a = split(/\t/);
    my $nt = substr($a[2], 0, 1);
    if ( length($a[2]) > 1 ) {
	my $p = $a[1];
	if ( length($a[3]) == 1 ) {
	    $p++;
	    $nt = substr($a[2], 1, 1);
	    my $m = "-" . (length($a[2]) - 1);
	    $mut{ "$a[0]\t$p\t$nt\t$m" } = "-" . (length($a[2]) - 1);
	} elsif ( length($a[2]) == length($a[3]) ) {
	    $mut{ "$a[0]\t$p\t$nt\t$a[3]" } = $a[3];
	} else {
	    $mut{ "$a[0]\t$p\t$nt\t$a[3]" } = $a[3];
	}
    } elsif ( length($a[3]) > 1 ) {
	my $m = "I+" . substr($a[3], 1);
        $mut{ "$a[0]\t$a[1]\t$nt\t$m" } = "I+" . substr($a[3], 1);
    } else {
	$mut{ "$a[0]\t$a[1]\t$nt\t$a[3]" } = $a[3];
    }
}
close( MUT );
my %gt;
my %samples;
my %cov;
while(<>) {
    chomp;
    my @a = split(/\t/);
    next unless( /:F-/ );
    my @b = split(/ & /, $a[$#a] );
    my $p = $a[3];
    my $nt = substr($a[5], 0, 1);
    my $amp = $a[$#a-2]; # The amplicon
    #next unless( $a[7] > $MINDP );
    $cov{ $a[0] }->{ "$a[2]\t$p\t$nt" }->{ $amp } = join("\t", @a[7,9,10]);
    if ( length($a[5]) > 1 && length($a[6]) == 1 && $a[5] =~ /^$a[6]/ ) {
        $p++; $nt = substr($a[5], 1, 1);
    }
    my $k = "$a[2]\t$p\t$nt";
    $cov{ $a[0] }->{ $k }->{ $amp } = join("\t", @a[7,9,10]);
    $samples{$a[0]} = 1;
    foreach my $b (@b) {
	#T:3282:F-1638:R-1644:0.99636:2:41.2:9.45:37.3:4.6
	last if ( $b =~ /^Ref/ || $b =~ /^Alt/ );
        my ($NT, $dp, $f, $r, $af, $bias, $pmean, $pstd, $qmean, $qstd, $hifreq, $mq, $sn) = split(/:/, $b);
	next if ( ($pmean < 5) ); #&& $af < 0.1;
	my $K = "$k\t$NT";
	#$af = 0 unless( $af > $opt_f );
	#$hifreq = 0 unless( $hifreq > $opt_f );
	if ($mut{$K} ) {
	    #$gt{ $a[0] }->{ $K } = join("\t", @a[7,9,10], $dp, $af, sprintf("%.5f", $hifreq), $f, $r, $bias, $pmean, $qmean, $mq, $sn);
	    $gt{ $a[0] }->{ $K }->{ $amp } = join("\t", @a[7,9,10], $dp, $af, sprintf("%.5f", $hifreq), $f, $r, $bias, $pmean, $qmean, $mq, $sn, $a[$#a-1]);
	}
    }
}

foreach my $sam (keys %samples) {
    foreach my $k (keys %mut) {
	my $tk = $k; $tk =~ s/\t\S+$//;
	unless( $gt{ $sam }->{ $k } ) { # mutation not observed
	    if ($cov{ $sam }->{ $tk }) {
	        while(my ($amp, $cv) = each %{ $cov{ $sam }->{ $tk } }) {
		    print join("\t", $sam, $k, $cv, 0, 0, 0, "F-0", "R-0", "", 0, 0, 0, 0, $amp), "\n";
		}
	    } else {  # No Coverage
		print join("\t", $sam, $k, $cov{ $sam }->{ $tk } ? $cov{ $sam }->{ $tk } : "0\t0\t0", 0, "", "", "F-0", "R-0", "", 0, 0, 0, 0, 0), "\n";
	    }
	    next;
	}
	while(my ($amp, $cv) = each %{ $cov{ $sam }->{ $tk } }) {
	    if ( $gt{ $sam }->{ $k }->{ $amp } ) {  # mutations is observed in the amplicon
	        print join("\t", $sam, $k, $gt{ $sam }->{ $k }->{ $amp }), "\n";
	    } else {  # mutations is NOT observed in the amplicon
		print join("\t", $sam, $k, $cv, 0, 0, 0, "F-0", "R-0", "", 0, 0, 0, 0, $amp), "\n";
	    }
	}
    }
}
#    unless( $flag ) {
#	print join("\t", @a[0..1], $k, $mut{ $k }, $a[7], 0, 0, "", 0, 0, 0, 0), "\n" if ( $mut{$k} );
    #}
#50      TP53    chr17   7577120 7577120 C       T       105641  35366   24503   45671   12499   22867   C/T     0.33478 2;2     23.7    1       35.1    1       59.9    23.637  0.33933 0       0       0       GGCACAAACA      GCACCTCAAA      A:81:F-20:R-61:0.00077:2:23.4:1:13.8:1:9.00045002250112e-05:59.1:0.124 & T:35366:F-12499:R-22867:0.33478:2:23.7:1:35.1:1:0.339326966348317:59.9:23.637 & C:70174:F-24503:R-45671:0.66427:2:23.6:1:34.5:1:0.660573028651433:59.9:16.031 & G:20:F-7:R-13:0.00019:2:24.4:1:9.8:1:1.00005000250013e-05:55.0:0.051 & Ref|23.6|1|34.5|1|0.660573028651433|59.9|16.031 & Alt|23.7|1|35.1|1|0.339326966348317|59.9|23.637


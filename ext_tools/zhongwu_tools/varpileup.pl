#!/usr/bin/perl

# The program will produce all variants, including noise, from the output of checkVar.pl with -D option
use strict;

print join("\t", qw(Sample Gene Chr Start End Ref Alt TotalDepth AltDepth RefFwd RefRev AltFwd AltRev AF Bias Pmean Pstd Qmean Qstd Type)), "\n";
while( <> ) {
    chomp;
    my @a = split(/\t/);
    my @b = split(/ & /, $a[22]);
    #pop(@b); pop(@b);
    my @c = ();
    my $rb = $1 if ( $a[15] =~ /^(\d)/ );
    foreach $b (@b) {
	next if ( $b =~ /^Ref/ || $b =~ /^Alt/ );
        my @t = split(/:/, $b);
	shift(@t) if ($b =~ /^I:/ );
	$t[2] =~ s/F-//;
	$t[3] =~ s/R-//;
	$t[5] = "$rb;$t[5]";
	push(@c, [$t[0], @t[1..9], $t[1]*$t[8]]);
    }
    @c = sort {$b->[1] <=> $a->[1]} @c;
    my $n = 1;
    foreach my $c (@c) {
	my $type = $c->[0] eq $a[5] ? "Ref" : $n++;
        print join("\t", @a[0..5], $c->[0], $a[7], $c->[1], @a[9,10], @$c[2..9], $type), "\n";
    }
}
#NCI-H1975     PIK3CA  chr3    178952085       178952085       A       G       14049   12      8903    5112    7       5       A/G     0.00085 2;2     36.1    27.13   32.1    10.0    AATGATGCAC      TCATGGTGGC        A:14015:F-8903:R-5112:0.99758:2:28.0:27.51:34.6:5.3 & T:8:F-4:R-4:0.00057:2:50.5:12.06:38.4:1.1 & C:4:F-3:R-1:0.00028:2:42.8:26.99:15.2:0.5 & N:10:F-6:R-4:0.00071:2:31.2:26.67:2.0:0.0 & G:12:F-7:R-5:0.00085:2:36.1:27.13:32.1:10.0 & I++T:1:F-1:R-0:0.00007:0:4.0:0.00:34.0:0.0 & Ref|28.0|27.51|34.6|5.3 & Alt|36.1|27.13|32.1|10.0
#PromegaMale     AKT1    chr14   105239394       105239394       G       C       10981   2       5449    5500    2       0       0.00018 2;0     55.5    6.36    15.0    1.4     3
#PromegaMale     AKT1    chr14   105239394       105239394       G       Ref|56.9|9.81|35.7|6.3  10981           5449    5500                            2;                                      4
#PromegaMale     AKT1    chr14   105239394       105239394       G       Alt|55.5|11.50|17.6|5.6 10981           5449    5500                            2;                                      5


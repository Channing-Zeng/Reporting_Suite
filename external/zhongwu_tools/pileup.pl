#!/usr/bin/perl -w

# Calculate the pileup of a given region

use strict;

my %hash;
my $chr;
while( <> ) {
    my @a = split(/\t/);
    my $start = $a[3];
    $chr = $a[2];
    my $n = 0;
    while( $a[5] =~ /(\d+)([A-Z])/g ) {
        my $m = $1;
	my $C = $2;
	if ( $C eq "N" ) {
	    $start += $m;
	    next;
	} elsif ( $C eq "S" ) {
	    $n += $m;
	    next;
	} elsif ( $C eq "I" ) {
	    my $s = substr($a[9], $n, $m);
	    $hash{ $start - 1 }->{ I }->{ $s }++;
	    $n += $m;
	    next;
	} elsif ( $C eq "D" ) {
	    my $s = "-$m";
	    $hash{ $start - 1 }->{ $s }++;
	    $start += $m;
	    next;
	}
	for(my $i = 0; $i < $m; $i++) {
	    my $s = substr($a[9], $n, 1);
	    $s = "-$m" if ( $C eq "D" );
	    my $t = $C eq "D" ? $start - 1 : $start;
	    $hash{ $t }->{ $s }++;
	    $start++ unless( $C eq "I" );
	    $n++ unless( $C eq "D" );
	}
    }
}

while( my ($p, $v) = each %hash ) {
    my @tmp = ();
    while( my ($n, $cnt) = each %$v ) {
        push( @tmp, "$n:$cnt" ) unless( $n eq "I");
    }
    if ( $v->{ I } ) {
        while( my ($n, $cnt) = each %{ $v->{ I } } ) {
	    push( @tmp, "I:$n:$cnt" );
	}
    }
    print join("\t", $chr, $p, @tmp), "\n";
}
#FCB02N4ACXX:3:1106:11473:57062#GATGGTTC 163     chr3    38181980        50      80M188N10M      =       38182274        667     CTGAAGTTGTGTGTGTCTGACCGCGATGTCCTGCCTGGCACCTGTGTCTGGTGTATTGCTAGTGAGCTCATCGAAAAGAGGTGCCGCCGG ___\ceacgbe^cghghhhhhhfhifdhfhhhhefheb_agfe^MWWaegfa_MW\_S\Z\c^ddgece]acbURZ^b_```[bc^ca[a AS:i:-4 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:52C37      YT:Z:UU XS:A:+  NH:i:1  RG:Z:15
#FCB02N4ACXX:3:2206:20108:2526#GATGGTTC  163     chr3    38181981        50      79M188N11M      =       38182275        667     TGAAGTTGTGTGTGTCTGACCGCGATGTCCTGCCTGGCACCTGTGTCTGGTCTATTGCTAGTGAGCTCATCGTAAAGAGGTGCCGCCGGG \YY`c`\ZQPJ`e`b]e_Sbabc[^Ybfaega_^cafhR[U^ee[ec][R\Z\__ZZbZ\_\`Z`d^`Zb]bBBBBBBBBBBBBBBBBBB AS:i:-8 XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:72A16A0    YT:Z:UU XS:A:+  NH:i:1  RG:Z:15


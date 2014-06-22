#!/usr/bin/perl -w

use strict;

my %hash;

my $unmapped = 0;
while( <> ) {
    chomp;
    my @a = split(/\t/);
    my ($r, $f, $s, $q) = @a[0,1, 9,10];
    if ( $f & 0x10 ) {
        $s = reverse($s);
	$s =~ tr/ATGC/TACG/;
        $q = reverse($q);
    }
    $unmapped++ if ( $f & 0x4 );
    $r =~ s/\/\d$//;
    my $g = $f & 0x40 ? 1 : 2;
    #print STDERR "Conflict of flat $g $f for $r\n" if ( $f & 0x80 && $g == 1 );
    $r =~ s/\/.$//;
    $hash{ $r }->{ $g } = { s => $s, q => $q, um => $f & 0x4 };
}

open( FQ1, ">1.fastq");
open( FQ2, ">2.fastq");
open( SIN, ">singletons.fastq");
open( STAT, ">SAM_stat.txt" );
my ($paired, $unpaired) = (0, 0);
my ($ums, $ms) = (0, 0);
while( my ($read, $ref) = each %hash ) {
    unless( $ref->{ 1 } && $ref->{ 2 } ) {
        $unpaired++;
	my $t = $ref->{ 1 } ? $ref->{ 1 } : $ref->{ 2 };
	if ( $ref->{ 1 } ) {
	    print FQ1 "\@$read/1\n$t->{ s }\n+\n$t->{ q }\n";
	} else {
	    print FQ2 "\@$read/2\n$t->{ s }\n+\n$t->{ q }\n";
	}
	#print SIN "\@$read/2\n$t->{ s }\n+\n$t->{ q }\n";
	if ( $t->{ um } ) {
	    $ums++;
	} else {
	    $ms++;
	}
	next;
    }
    $paired++;
    print FQ1 "\@$read/1\n$ref->{ 1 }->{ s }\n+\n$ref->{ 1 }->{ q }\n";
    print FQ2 "\@$read/2\n$ref->{ 2 }->{ s }\n+\n$ref->{ 2 }->{ q }\n";
}
close( FQ1 );
close( FQ2 );
close( SIN );
print STAT "Paired reads: $paired pairs.\n";
print STAT "Unpaired reads: $unpaired singletons.\n";
print STAT "Unmapped reads: $unmapped reads.\n";
print STAT "Unmapped singletons: $ums reads.\n";
print STAT "Mapped singletons: $ms reads.\n";
close( STAT );

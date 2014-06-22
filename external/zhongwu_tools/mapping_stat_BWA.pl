#!/usr/bin/perl -w

use strict;

my $sample = shift;
my %hash;
my %XT;
my %mates;
my %unmap; $unmap{ 1 } = 0; $unmap{ 2 } = 0;
my %singleton; $singleton{ 1 } = 0; $singleton{ 2 } = 0;
my %total; $total{ 1 } = 0; $total{ 2 } = 0;
my %nokey; $nokey{ 1 } = 0; $nokey{ 2 } = 0;
while( <> ) {
    next if ( /^@/ );
    my @a = split(/\t/);
    my $dir = $a[1] & 0x40 ? 1 : 2;
    if ( $a[1] & 0x4 ) {
        $unmap{ $dir }++;
	next;
    }
    $total{ $dir }++;
    my $k = $1 if ( /X0:i:(\d+)/ );
    my $xt = $1 if ( /XT:A:(\S+)/ );
    $k = $xt unless ( $k );
    unless( $k ) {
        $k = "NM\t$1" if ( /NM:i:(\d+)/ );
    }
    unless( $k ) {
        #print STDERR;
	$nokey{ $dir }++;
	next;
    }
    $XT{ $xt }{ $dir }++ if ( $xt );
    $hash{ $k }{ $dir }++;
    if( $a[1] & 0x8 ) {
        $singleton{ $dir }++;
    } else {
	$mates{ $k }->{ $dir }++ 
    }
}

while(my ($k, $v) = each %hash) {
    my $v1 = $v->{ 1 } ? $v->{ 1 } : 0;
    my $v2 = $v->{ 2 } ? $v->{ 2 } : 0;
    print join("\t", $sample, $k, $v1, $v2, $k =~ /^[\d\.]+$/ ? ($v1 + $v2)/$k : $v1 + $v2), "\n";
}
while(my ($k, $v) = each %mates) {
    my $v1 = $v->{ 1 } ? $v->{ 1 } : 0;
    my $v2 = $v->{ 2 } ? $v->{ 2 } : 0;
    print join("\t", $sample, $k, $v1, $v2, $v1 + $v2, "Pairs"), "\n";
}
while(my ($k, $v) = each %XT) {
    my $v1 = $v->{ 1 } ? $v->{ 1 } : 0;
    my $v2 = $v->{ 2 } ? $v->{ 2 } : 0;
    print join("\t", $sample, $k, $v1, $v2, $v1 + $v2, "XT"), "\n";
}

print join("\t", "Total:", $total{ 1 }, $total{ 2 }, $total{ 1 } + $total{ 2 }, $sample), "\n";
print join("\t", "Unmapped:", $unmap{ 1 }, $unmap{ 2 }, $unmap{ 1 } + $unmap{ 2 }, $sample), "\n";
print join("\t", "Singleton:", $singleton{ 1 }, $singleton{ 2 }, $singleton{ 1 } + $singleton{ 2 }, $sample), "\n";
print join("\t", "Nokey:", $nokey{ 1 }, $nokey{ 2 }, $nokey{ 1 } + $nokey{ 2 }, $sample), "\n";

#!/usr/bin/perl -w

# Process MSI locations and produce histogram for all alelles

use strict;

print join("\t", "Sample", "Gene", "Chr", "Pos", "MSI-n", "Read_Cnt"), "\n";
my %hash; # key:  sample-gene-chr-pos Value: hash of n, cnt
while( <> ) {
    my @a = split(/\t/);
    my ($sample, $gene, $chr, $pos) = @a[0, 1, 2, 3];
    my @b = split(/ & /, $a[$#a]); # the alleles
    my $nt = $a[6];
    if ( length($a[5]) > 1 ) { # Deletion
        $pos++;
	#$a[5] = substr($a[5], 1, 1);
	$nt = substr($a[5], 1, 1);
	#$nt = $a[5];
    } elsif ( length($a[6]) > 1 ) { # Insertion
        $pos++;
	$nt = substr($a[6], 1, 1);
    }
    my $k = join("\t", $sample, $gene, $chr, $pos);
    next unless($b[0] =~ /:/);
    for(my $i = 0; $i < @b; $i++) {
	next if ( $b[$i] =~ /^Ref/ || $b[$i] =~ /Alt/ );
	my @c = split( /:/, $b[$i] );
	my ($n, $cnt) = (0, 0);
        if ( $b[$i] =~ /^I\+($nt+):/ && length($a[6]) > 1 && length($a[5]) == 1 ) { # Insertion
	    ($n, $cnt) = (length($1), $c[1]);
	} elsif (length($a[5]) > 1 && length($a[6]) ==1) {  # Deletion
	    if ( $c[0] =~ /^(-\d+)$/ ) {
	        ($n, $cnt) = ($1, $c[1]);
	    } elsif ( $c[0] eq $nt ) {
	        ($n, $cnt) = (0, $c[1]);
		next if ( length($a[6]) > 1 ); # Not use when it's insertion
	    } else {
	        # Ignore non-ref nucleotides
		#print STDERR "$c[1] $c[0] vs $a[5] ignored\n";
		#print join("\t", $sample, $gene, $chr, $pos, -0.5, $c[1]), "\n" if ( $c[1] > 0 );
		#$hash{ $k }->{ 0.5 } = $c[1];
		next;
	    }
	}
	#print STDERR "$n\t$cnt\n";
	$hash{ $k }->{ $n } = $cnt if ( $cnt > 0 );
	$hash{ $k }->{ I } += $cnt if ( $cnt > 0 && $n >= 1 );
	#print join("\t", $sample, $gene, $chr, $pos, $n, $cnt), "\n" if ( $cnt > 0 );
    }
}
while(my ($k, $v) = each %hash) {
    while(my ($n, $cnt) = each %$v) {
        next if ( $n eq "I" );
	$cnt -= $v->{ I } if ( $n == 0 && $v->{ I } );
	print join("\t", $k, $n, $cnt), "\n";
    }
}

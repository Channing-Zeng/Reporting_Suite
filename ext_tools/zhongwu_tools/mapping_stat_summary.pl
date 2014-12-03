#!/usr/bin/perl -w

use Getopt::Std;
use strict;
# Generate the mapping stats

our ($opt_s);
getopts('s:');

my %total;
my $all = 0;
if ( $opt_s ) {  # The read count summary
    open( SUM, $opt_s );
    while( <SUM> ) {
        chomp;
	my @a = split(/\t/);
	$total{ $a[0] } = $a[1];
	$all += $a[1];
    }
    close( SUM );
}

my %hash;
my $sample;
while( <> ) {
    $sample = $1 if ( /(.*)_sam\t/ );
    my @a = split( /\t|\n/ );
    $sample = $a[4] if ( @a > 4 && /:/);
    next unless( $sample );
    $sample =~ s/_sam$//;
    $hash{ $sample }->{ mapped } = $a[3] if ( /Total/ );
    $hash{ $sample }->{ forward } = $a[1] if ( /Total/ );
    $hash{ $sample }->{ reverse } = $a[2] if ( /Total/ );
    $hash{ $sample }->{ unmap } = $a[3] if ( /Unmap/ );
    if ( /Total/ || /Unmap/ ) {
        $all += $a[3] unless( $opt_s );
    }
}

print join("\t", "Sample", "Mapped", "Unmapped", "Total", '%Mapped', '%total'), "\n";
my @samples;
while( my ($s, $r) = each %hash ) {
    push(@samples, $s) unless( $s eq "Undetermined" );
}
@samples = sort (@samples);
my ($tmap, $tunmap) = (0,0);
push(@samples, "Undetermined") if ( $hash{ "Undetermined" } );

foreach my $s (@samples) {
    my $r = $hash{ $s };
    #print join("\t", $s, $r->{ mapped }, $r->{forward}, $r->{ reverse } ? $r->{ reverse } : 0, $r->{unmap}, $r->{ mapped }+$r->{unmap}, $r->{ mapped }/($r->{ mapped }+$r->{unmap})), "\n";
    if ( $opt_s ) {
	print join("\t", $s, $total{ $s } - $r->{ unmap }, $r->{unmap}, $total{ $s }, sprintf("%.4f", 1 - $r->{ unmap }/$total{ $s }), sprintf("%.4f", $total{ $s }/$all) ), "\n";
	$tmap += $total{ $s } - $r->{ unmap };
	$tunmap += $r->{unmap};
    } else {
	print join("\t", $s, $r->{ mapped }, $r->{unmap}, $r->{ mapped }+$r->{unmap}, sprintf("%.4f", $r->{ mapped }/($r->{ mapped }+$r->{unmap})), sprintf("%.4f", ($r->{ mapped }+$r->{unmap})/$all) ), "\n";
	$tmap += $r->{ mapped };
	$tunmap += $r->{unmap};
    }
}
print join("\t", "Total", $tmap, $tunmap, $all, sprintf("%.4f", $tmap/$all), 1.00), "\n";

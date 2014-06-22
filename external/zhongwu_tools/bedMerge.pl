#!/usr/bin/perl -w

# Merge overlapping segments 

use strict;

my ($pg, $pc, $ps, $pe);
my @amp = ();

while( <> ) {
    my @a = split(/\t/);
    my ($g, $c, $s, $e) = @a[6, 2, 3, 4];
    if ( $pg && $pg eq $g && $pc eq $c && $s < $pe ) {
        $pe = $e; 
	push(@amp, $a[1]);
    } elsif ( ($pe && $pg eq $g && $s > $pe) || ($pe && $pg ne $g ) ) {
        print join("\t", $pg, $pc, $ps, $pe, join(";", @amp)), "\n";
	@amp = ($a[1]);
	($pg, $pc, $ps, $pe) = ($g, $c, $s, $e);
    } else {
        ($pg, $pc, $ps, $pe) = ($g, $c, $s, $e);
	push(@amp, $a[1]);
    }
}
print join("\t", $pg, $pc, $ps, $pe, join(";", @amp)), "\n";

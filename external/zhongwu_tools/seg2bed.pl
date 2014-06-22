#!/usr/bin/perl -w

# extract the design region from amplicon designs.  Assuming the coordinates are ordered.

#52      CARD11_ex1_1a   chr7    2946213 2946384 Fluidigm_1269AAP-1501AAP_union  CARD11  AAA0007039
#34      CARD11_ex1_2    chr7    2946314 2946491 Fluidigm_1269AAP-1501AAP_union  CARD11  AAA0007016
#22      CARD11_ex1_3    chr7    2946416 2946606 Fluidigm_1269AAP-1501AAP_union  CARD11  AAA0006996
#24      CARD11_ex2_1a   chr7    2949575 2949774 Fluidigm_1269AAP-1501AAP_union  CARD11  AAA0007000
#40      CARD11_ex2_2a   chr7    2949675 2949865 Fluidigm_1269AAP-1501AAP_union  CARD11  AAA0007024

use Getopt::Std;
use strict;

our ($opt_g, $opt_c, $opt_s, $opt_e);
getopts("g:c:s:e:");

my $col_g = $opt_g ? $opt_g - 1 : 3;
my $col_c = $opt_c ? $opt_c - 1 : 0;
my $col_s = $opt_s ? $opt_s - 1 : 1;
my $col_e = $opt_e ? $opt_e - 1 : 2;
my ($pg, $pc, $ps, $pe);
my @amp = ();

while( <> ) {
    chomp;
    my @a = split(/\t/);
    my ($g, $c, $s, $e) = @a[$col_g, $col_c, $col_s, $col_e];
    if ( $pg && $pg eq $g && $pc eq $c && $s < $pe ) {
        $pe = $e; 
    } elsif ( ($pe && $pg eq $g && $s > $pe) || ($pe && $pg ne $g ) ) {
        print join("\t", $pc, $ps, $pe, $pg), "\n";
	($pg, $pc, $ps, $pe) = ($g, $c, $s, $e);
    } else {
        ($pg, $pc, $ps, $pe) = ($g, $c, $s, $e);
    }
}
print join("\t", $pc, $ps, $pe, $pg), "\n";

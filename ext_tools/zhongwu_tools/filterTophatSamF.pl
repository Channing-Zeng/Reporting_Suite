#!/usr/bin/perl -w

# Filter any read maps with F (Fusion) in Cigar in SAM output by tophat-fusion, so that it can be recognized by samtools
# and converted into BAM for visualization

use strict;

while( <> ) {
    my @a = split( /\t/ );
    next if ( $a[5] && $a[5] =~ /\S+F\S+/ );
    print;
}

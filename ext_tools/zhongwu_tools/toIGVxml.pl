#!/usr/bin/perl -w

use Getopt::Std;
use strict;

# Convert a list of BAM files to xml for IGV data server

our ($opt_n, $opt_t, $opt_c, $opt_h, $opt_C);
getopts('t:n:c:h:C:');

if ( $opt_t ) {
    print <<HEADER;
<?xml version="1.0" encoding="UTF-8"?>
<Global name="$opt_t"  infolink="http://www.broadinstitute.org/igv/" version="1">
HEADER
}
my $host = $opt_h ? $opt_h : "http://luna899.usbod.astrazeneca.net/~kdld047/";
my $indent = "";
if ( $opt_C ) {
    $indent = "    ";
    print qq|    <Category name="$opt_C">\n|;
}
print qq|$indent    <Category name="$opt_c">\n| if ( $opt_c );
while(<>) {
    my $sample = "";
    if ( $opt_n ) {
        /$opt_n/ && ($sample = $1);
    } else {
        /([^\/]+).bam/ && ($sample = $1);
    }
    my $link = $1 if ( /(igv.*bam)/ );
    #$link = "http://luna899.usbod.astrazeneca.net/~kdld047/" . $link;
    $link = $host . $link;
    print <<LINK;
$indent    <Resource name="$sample"
$indent            path="$link" />
LINK
}
print "$indent    </Category>\n" if ( $opt_c );
print "    </Category>\n" if ( $opt_C );
print "</Global>\n" if ( $opt_t );

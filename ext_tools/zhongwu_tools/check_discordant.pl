#!/usr/bin/perl -w

while( <> ) {
    @a = split(/\t/); 
    next if ( $a[6] eq "*" ); 
    print join("\t", @a[2,3,6, 7,8]), "\n" if ( abs($a[8]) > 100000 && $a[6] eq "="); 
    print join("\t", @a[2,3,6, 7,8]), "\n" if ( $a[6] ne "=" && $a[2] ne $a[6]);
}

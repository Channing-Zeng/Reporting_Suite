#!/usr/bin/perl -w

# The program will add queues to qsub for big NGS jobs.

use strict;

my @queues = ("-l huge_ram=1 -l h=chara", "-l huge_ram=1 -l h=rask", "-l huge_ram=1 -l h=chara", "-l huge_ram=1 -l h=rask", "-l huge_ram=1 -l h=chara", "-l huge_ram=1 -l h=rask", "-l huge_ram=1 -l h=rask");

while( <> ) {
    if ( /^qsub/ ) {
        my $q = $queues[int(rand()*7)];
	s/qsub/qsub $q/;
    }
    print;
}

#!/usr/bin/perl -w

# The program will add queues to qsub for big NGS jobs.

use Getopt::Std;
use strict;

our ($opt_o);
getopts("o:");

my $options = $opt_o ? $opt_o : "";

my @queues = ("-l huge_ram=1 -l h=chara", "-l huge_ram=1 -l h=rask", "-l h=orr", "-l huge_ram=1 -l h=chara", "-l h=espo", "-l huge_ram=1 -l h=rask", "-l huge_ram=1 -l h=chara", "-l huge_ram=1 -l h=rask", "-l huge_ram=1 -l h=rask");

while( <> ) {
    if ( /^qsub/ ) {
        my $q = $queues[int(rand()*9)];
	s/qsub/qsub $q $options/;
    }
    print;
}

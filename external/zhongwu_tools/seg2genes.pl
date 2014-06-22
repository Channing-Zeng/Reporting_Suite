#!/usr/bin/perl -w

use Getopt::Std;
use strict;

our ($opt_H, $opt_h, $opt_S, $opt_c, $opt_s, $opt_e, $opt_d);

getopts('ShHc:s:e:d:') || USAGE();
my $chr_c = $opt_c ? $opt_c -1 : 0;
my $start_c = $opt_s ? $opt_s -1 : 1;
my $end_c = $opt_e ? $opt_e -1 : 2;
my $DEL = $opt_d ? qr/$opt_d/ : qr/\t/;
USAGE() if ( $opt_h || $opt_H );

## Process gene location information
my %bin2gene;
open(GLOC, "/users/kdld047/igv/genomes/hg19/refGene.txt"); # IGV's reference gene file
my %gloc;
while( <GLOC> ) {
    chomp;
    s/\r//g;
    my @a = split(/\t/);
    my ($g, $chr, $s, $e) = @a[12,2,4,5];  # from IGV Broad's mapping file refGene.txt
    next unless( $chr =~ /chr[0-9XY]+$/ );
    #next if ( $g =~ /^MIR/ || /NR_\d+\t/ ); # Ignore microRNA, pseudo, or non-coding RNAs
    if ( $gloc{ $g } ) {
        my $flag = 0;
        while( my ($chrt, $v) = each %{ $gloc{ $g } } ) {
	    foreach my $vv (@$v) {
		my ($st, $et) = @$vv;
		if ( $chr eq $chrt && $s < $et && $st < $e ) {
		    $s = $st < $s ? $st : $s;
		    $e = $et > $e ? $et : $e;
		    $vv->[0] = $s;
		    $vv->[1] = $e;
		    $flag = 1;
		}
	    }
        }
        push(@{$gloc{ $g }->{ $chr }}, [$s, $e]) unless ($flag);
    } else {
        push( @{ $gloc{ $g }->{ $chr } }, [$s, $e] );
    }
}

my $BIN = 500000000;
while(my ($g, $v) = each %gloc) {
    while(my ($chr, $vv) = each %$v) {
        foreach my $t (@$vv) {
	    my ($s, $e) = @$t;
	    my $b1 = int($s/$BIN);
	    my $b2 = int($e/$BIN);
	    for(my $i = $b1; $i <= $b2; $i++) {
		push(@{ $bin2gene{ $chr }->{ $i } }, [$g, $chr, $s, $e]);
	    }
	}
    }
}
close( GLOC );
###################

my $chr = "";
my %result;
while( <> ) {
    chomp;
    s/\r//g;
    s/^\s+//;
    my @a = split(/$DEL/);
    my ($c, $s, $e) = @a[$chr_c, $start_c, $end_c];
    my $r = getGenes($c, $s, $e);
    my @genes = ();
    while(my($g, $v) = each %$r) {
        push(@genes, "$g-$v");
	if ( $opt_S ) {
	    print join("\t", $_, $g, $v), "\n";
	}
    }
    unless ( $opt_S ) {
	print join("\t", $_, join(":", @genes)), "\n";
    #} else {
	#print join("\t", $_, ""), "\n" unless( @genes > 0 );
    }
}

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text
}

sub getGenes {
    my ($chr, $start, $end) = @_;
    my $b1 = int($start/$BIN);
    my $b2 = int($end/$BIN);
    my %rs = ();
    for(my $i = $b1; $i <= $b2; $i++) {
        foreach my $r ( @{ $bin2gene{ $chr }->{ $i } } ) {
	    my ($g, $c, $s, $e) = @$r;
	    if ( $s <= $end && $start <= $e ) {
		my @t = sort { $a <=> $b } ($s, $start, $end, $e);
		my $cov = ($t[2]-$t[1])/($e-$s);
	        if ( $rs{ $g } ) {
		    #print STDERR "Two same genes on one segments: $g, $cov, $rs{$g}, $c, $s, $e, $chr, $start, $end\n";
		    $rs{ $g } = $cov;
		} else {
		    $rs{ $g } = $cov;
		}
	    }
	}
    }
    return \%rs;
}

sub USAGE {
    print <<USAGE;
    $0 [-S] [-c chr] [-s start] [-e end] segment_file

    The program will extract the genes from a segment file and append the genes, and fraction that the gene is covered by the segment.
    
    Options are (Columns are numbered from 1):

    -S If set, the segment data will be replicated for each gene it contains.  By default, it'll append the gene list
       concatinated by : to the end of the line.
    -c int
        The column has the chromosome, in the format of chr*.  Default 1, or first column.
    -s int
        The column has the segment start.  Default 2, or 2nd column.
    -e int
        The column has the segment end.  Default 3, or 3rd column.
    -d char
        The delimiter.  Default to "\\t"
USAGE
    exit(0);
}

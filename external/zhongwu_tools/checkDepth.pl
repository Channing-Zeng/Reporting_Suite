#!/usr/bin/perl -w
# Parse a list of refseq and check CDS coverage

use lib '/users/kdld047/lib';
use Getopt::Std;
use Stat::Basic;
use strict;

our ($opt_h, $opt_b, $opt_s, $opt_c, $opt_S, $opt_E, $opt_n, $opt_e, $opt_g, $opt_x, $opt_o, $opt_d, $opt_z, $opt_C, $opt_N);
USAGE() unless( getopts( 'hzCb:s:e:S:E:n:c:g:x:o:d:N:' ) );
USAGE() if ( $opt_h );

my $BAM = $opt_b; # the bam file
my $sample = $1 if ( $BAM =~ /([^\/.]+)[\/]*.bam/ );
if ( $opt_n ) {
    $sample = $1 if ( $BAM =~ /$opt_n/ );
}
$sample = $opt_N if ( $opt_N );
my $c_col = $opt_c ? $opt_c - 1 : 2;
my $S_col = $opt_S ? $opt_S - 1 : 6;
my $E_col = $opt_E ? $opt_E - 1 : 7;
my $s_col = $opt_s ? $opt_s - 1 : 9;
my $e_col = $opt_e ? $opt_e - 1 : 10;
my $g_col = $opt_g ? $opt_g - 1 : 12;
my $o_col = $opt_o ? $opt_o - 1 : 3;

my $SPLICE = $opt_x ? $opt_x : 0;
my $stat = new Stat::Basic;
my %regions;

while( <> ) {
    next if ( /^track/i );
    next if ( /^browser/i );
    next if ( /^#/ );
    chomp;
    my @A = split(/\t/);
    my ($chr, $ori, $cdss, $cdse, $gene) = @A[$c_col,$o_col,$S_col,$E_col,$g_col];
    my @starts = split(/,/, $A[$s_col]);
    my @ends = split(/,/, $A[$e_col]);
    my @CDS = ();
    for(my $i = 0; $i < @starts; $i++) {
        my ($s, $e) = ($starts[$i], $ends[$i]);
	next if ( $cdss > $e ); # not a coding exon
	last if ( $cdse < $s ); # No more coding exon
	$s = $cdss if ( $s < $cdss );
	$e = $cdse if ( $e > $cdse );
	$s -= $SPLICE unless ( $s == $cdss );
	$e += $SPLICE unless ( $e == $cdse );
	$s += 1 if ( $opt_z );
	push(@CDS, [$s, $e]);
    }
    push( @{ $regions{ $gene }->{ CDS } }, @CDS);
    $regions{ $gene }->{ chr } = $chr;
    $regions{ $gene }->{ ori } = $ori;
}

my @depths = $opt_d ? split(/:/, $opt_d) : (1, 5, 10, 25, 50, 100, 250, 500, 1000); # The depth criteria
print join("\t", "Sample", "Gene", "Chr", "Start", "End", "Tag", "Length", "MeanDepth", @depths), "\n";
while( my ($gene, $r) = each %regions ) {
    my $cov = 0;
    my @exoncov;
    my $CDS = $r->{ CDS };
    my $ori = $r->{ ori };
    my $chr = $r->{ chr };
    my $len = 0;
    my ($gene_start, $gene_end) = (1000000000, 0);
    for(my $i = 0; $i < @{ $CDS }; $i++) {
	my $EXN = ($ori eq "+" || $ori eq "1") ? $i + 1 : @{ $CDS } - $i; # Exon number
	my ($START, $END) = @{$CDS->[$i]};
	my $tchr = $chr;
	$tchr =~ s/chr// if ( $opt_C );
	my $reads = `samtools view $BAM $tchr:$START-$END | wc -l`;
	$exoncov[$EXN] = $reads;
	$EXN = "0$EXN" if ( $EXN < 10 );
	print join("\t", $sample, $gene, $chr, $START, $END, "Amplicon", $END-$START+1, sprintf("%.2f", $reads/($END-$START+1))), "\n";
	$len += $END-$START+1;
	$cov += $reads;
	$gene_start = $START if ( $START < $gene_start );
	$gene_end = $END if ( $END > $gene_end );
    }
    print join("\t", $sample, $gene, $chr, $gene_start, $gene_end, "Whole-Gene", $len, sprintf("%.2f", $cov/$len)), "\n";
}

sub USAGE {
    print STDERR <<USAGE;
    $0 [-hz] [-n name_reg] [-b bam] [-c chr] [-S start] [-E end] [-s seg_starts] [-e seg_ends] [-x #_nu] [-g gene] [-o ori] [-d depth] region_info

    The program will calculate candidate variance for a given region(s) in an indexed BAM file.  The default
    input is IGV's one or more entries in refGene.txt, but can be regions

    -h Print this help
    -C Indicate whether to trim chr from chromosome name.  Set if the BAM file chrs are 1, 2... 
       instead of chr1, chr2...
    -n name_reg
       The regular expression to extract sample name from bam filename
    -N name
       Mutual exclusive to -n.  Set the sample name to name
    -b bam
       The indexed BAM file
    -c chr
       The column for chr
    -S start
       The column for region start, e.g. gene start
    -E end
       The column for region end, e.g. gene end
    -s seg_starts
       The column for segment starts in the region, e.g. exon starts
    -e seg_ends
       The column for segment ends in the region, e.g. exon ends
    -g gene
       The column for gene name
    -o orientation
       The column for orientation, as "+" or "-"
    -x num
       The number of nucleotide to extend for each segment, default: 2
    -d depths
       Desired depths, default to 1:5:10:25:50:100:250:500:1000
    -z 
       Indicate whether it's zero based numbering, default is 1-based
USAGE
   exit(0);
}

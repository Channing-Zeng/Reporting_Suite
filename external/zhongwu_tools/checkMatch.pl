#!/usr/bin/perl -w
# Parse a list of refseq and check CDS coverage

use lib "/users/kdld047/lib/perl5";
use lib "/users/kdld047/aris/lib";
use Getopt::Std;
use Stat::Basic;
use Fasta;
use strict;

our ($opt_h, $opt_b, $opt_D, $opt_d, $opt_s, $opt_c, $opt_S, $opt_E, $opt_n, $opt_N, $opt_e, $opt_g, $opt_x, $opt_f, $opt_r, $opt_B, $opt_z, $opt_v, $opt_p, $opt_F, $opt_C);
unless( getopts( 'hvzpDCFd:b:s:e:S:E:n:c:g:x:f:r:B:N:' )) {
    USAGE();
}
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

$s_col = $S_col if ( $opt_S && (!$opt_s) );
$e_col = $E_col if ( $opt_E && (!$opt_e) );

my $fasta = new Fasta( -fasta => "/users/kdld047/work/NGS/NGS/genomes/hg19_complete/hg19_complete.fa");
#899	NM_007300	chr17	-	41196311	41277500	41197694	41276113	24	41196311,41199659,41201137,41203079,41209068,41215349,41215890,41219624,41222944,41226347,41228504,41231350,41234420,41242960,41243451,41247862,41249260,41251791,41256138,41256884,41258472,41267742,41276033,41277287,	41197819,41199720,41201211,41203134,41209152,41215390,41215968,41219712,41223255,41226538,41228628,41231416,41234592,41243049,41246877,41247939,41249306,41251897,41256278,41256973,41258550,41267796,41276132,41277500,	0	BRCA1	cmpl	cmpl	1,0,1,0,0,1,1,0,1,2,1,1,0,1,1,2,1,0,1,2,2,2,0,-1,
my $SPLICE = defined($opt_x) ? $opt_x : 2;
my $FREQ = $opt_f ? $opt_f : 0.05;
my $BIAS = 0.1; # The cutoff to decide whether a positin has read strand bias
my $MINB = $opt_B ? $opt_B : 3; # The minimum reads for bias calculation
my $MINR = $opt_r ? $opt_r : 3; # The minimum reads for variance allele
if ( $opt_p ) {
    $FREQ = -1;
    $MINR = 0;
}
my $stat = new Stat::Basic;
$opt_d = "\t" unless( $opt_d );
my %hash;
while( <> ) {
    chomp;
    next if ( /^#/ );
    next if ( /^browser/i );
    next if ( /^track/i );
    my @A = split(/$opt_d/);
    my ($chr, $cdss, $cdse, $gene) = @A[$c_col, $S_col, $E_col, $g_col];
    my @starts = split(/,/, $A[$s_col]);
    my @ends = split(/,/, $A[$e_col]);
    my @CDS = ();
    my @EXONSEQ = ();
    my %REF;
    $chr = "chr$chr" unless ($chr =~ /^chr/ );
    for(my $i = 0; $i < @starts; $i++) {
        my ($s, $e) = ($starts[$i], $ends[$i]);
	next if ( $cdss > $e ); # not a coding exon
	last if ( $cdse < $s ); # No more coding exon
	$s = $cdss if ( $s < $cdss );
	$e = $cdse if ( $e > $cdse );
	$s -= $SPLICE unless ( $s == $cdss );
	$e += $SPLICE unless ( $e == $cdse );
	$s++ if ( $opt_z );
	push(@CDS, [$s, $e]);
	my $s_start = $s - $SPLICE - 100 < 1 ? 1 : $s - $SPLICE - 100;
	my $exon = $fasta->getSeq( -id => $chr, -ori => 1, -start => $s_start, -end => $e + $SPLICE + 100);
	for(my $i = $s_start; $i <= $e + $SPLICE + 100; $i++) {
	    $REF{ $i } = uc(substr( $exon, $i - ($s_start), 1 ));
	}
	push(@EXONSEQ, $exon);
    }
    my %cov;
    my %var;
    for(my $i = 0; $i < @CDS; $i++) {
	my ($START, $END) = @{$CDS[$i]};
	$chr =~ s/^chr// if ( $opt_C );
	open(SAM, "samtools view $BAM $chr:$START-$END |");
	while( <SAM> ) {
	    my @a = split(/\t/);
	    my $start = $a[3];
	    my $n = 0;
	    my $p = 0;
	    my $dir = $a[1] & 0x10 ? "-" : "+";
	    my @segs = $a[5] =~ /(\d+)[MI]/g;
	    my $rlen = $stat->sum(\@segs); # The read length
	    while( $a[5] =~ /(\d+)([A-Z])/g ) {
		my $m = $1;
		my $C = $2;
		if ( $C eq "N" ) {
		    $start += $m;
		    next;
		} elsif ( $C eq "S" ) {
		    $n += $m;
		    next;
		} elsif ( $C eq "I" ) {
		    $n += $m;
		    $p += $m;
		    next;
		} elsif ( $C eq "D" ) {
		    $start += $m;
		    next;
		}
		for(my $i = 0; $i < $m; $i++) {
		    my $s = substr($a[9], $n, 1);
		    my $q = substr($a[10], $n, 1);
		    if ( $start >= $START && $start <= $END ) {
			if ( $s eq $REF{ $start } ) {
			    $hash{ $n + 1 }->{ match }++;
			} else {
			    $hash{ $n + 1 }->{ mis }++;
			}
		    }
		    $start++ unless( $C eq "I" );
		    $n++ unless( $C eq "D" );
		    $p++ unless( $C eq "D" );
		}
	    }
	}
	close( SAM );
    }
}
my @pos = sort { $a <=> $b } ( keys %hash );
foreach my $p ( @pos ) {
    my $v = $hash{ $p };
    my $mis = $v->{ mis } ? $v->{ mis } : 0;
    print join("\t", $p, $v->{ match }, $mis, sprintf("%.5f", $mis/($v->{match}+$mis))), "\n";
}
#AZ01    chr6    106536253       G       38535   2       G/A
#FCB02N4ACXX:3:1106:11473:57062#GATGGTTC 163     chr3    38181980        50      80M188N10M      =       38182274        667     CTGAAGTTGTGTGTGTCTGACCGCGATGTCCTGCCTGGCACCTGTGTCTGGTGTATTGCTAGTGAGCTCATCGAAAAGAGGTGCCGCCGG ___\ceacgbe^cghghhhhhhfhifdhfhhhhefheb_agfe^MWWaegfa_MW\_S\Z\c^ddgece]acbURZ^b_```[bc^ca[a AS:i:-4 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:52C37      YT:Z:UU XS:A:+  NH:i:1  RG:Z:15
#FCB02N4ACXX:3:2206:20108:2526#GATGGTTC  163     chr3    38181981        50      79M188N11M      =       38182275        667     TGAAGTTGTGTGTGTCTGACCGCGATGTCCTGCCTGGCACCTGTGTCTGGTCTATTGCTAGTGAGCTCATCGTAAAGAGGTGCCGCCGGG \YY`c`\ZQPJ`e`b]e_Sbabc[^Ybfaega_^cafhR[U^ee[ec][R\Z\__ZZbZ\_\`Z`d^`Zb]bBBBBBBBBBBBBBBBBBB AS:i:-8 XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:72A16A0    YT:Z:UU XS:A:+  NH:i:1  RG:Z:15
sub USAGE {
    print STDERR <<USAGE;
    $0 [-n name_reg] [-b bam] [-c chr] [-S start] [-E end] [-s seg_starts] [-e seg_ends] [-x #_nu] [-g gene] [-f freq] [-r #_reads] [-B #_reads] region_info

    The program will calculate candidate variance for a given region(s) in an indexed BAM file.  The default
    input is IGV's one or more entries in refGene.txt, but can be any regions

    -h Print this help page
    -z Indicate wehther is zero-based cooridates, as IGV does.
    -v VCF format output
    -p Do pileup regarless the frequency
    -C Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
    -D Debug mode.  Will print some error messages and append full genotype at the end.
    -d delimiter
       The delimiter for split region_info, default to tab "\t"
    -n name_reg
       The regular expression to extract sample name from bam filenames
    -N sample name
       The sample name to be used directly.  Will overwrite -n
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
    -x num
       The number of nucleotide to extend for each segment, default: 2
    -f frequency
       The threshold for allele frequency, default: 0.05
    -F Indicate to calculate the frequency as the sum of all non-reference variants, 
       instead of just the most frequent allele, which is the default
    -r minimum reads
       The minimum # of variance reads
    -B minimum reads
       The minimum # of reads to determine strand bias
USAGE
   exit(0);
}

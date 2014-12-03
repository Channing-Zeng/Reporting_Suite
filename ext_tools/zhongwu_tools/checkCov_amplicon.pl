#!/usr/bin/perl -w
# Parse a list of refseq and check CDS coverage

use Getopt::Std;
use Stat::Basic;
use Fasta;
use strict;

our ($opt_h, $opt_b, $opt_s, $opt_c, $opt_S, $opt_E, $opt_n, $opt_e, $opt_g, $opt_x, $opt_o, $opt_z, $opt_t, $opt_N, $opt_a, $opt_p);
USAGE() unless( getopts( 'hztb:s:e:S:E:n:c:g:x:o:N:a:p:' ) );
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
my $PRIMERLEN = $opt_p ? $opt_p : 0;
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
    my @pos = ();
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
	for(my $i = $s; $i <= $e; $i++) {
	    push(@pos, $i);
	}
    }
    push( @{ $regions{ $gene }->{ CDS } }, @CDS);
    push( @{ $regions{ $gene }->{ POS } }, @pos);
    $regions{ $gene }->{ chr } = $chr;
    $regions{ $gene }->{ ori } = $ori;
}

print join("\t", "Sample", "Gene", "Chr", "Start", "End", "Tag", "Length", "#Reads"), "\n";
while( my ($gene, $r) = each %regions ) {
    my %cov;
    my @exoncov;
    my $CDS = $r->{ CDS };
    my $pos = $r->{ POS };
    my $ori = $r->{ ori };
    my $chr = $r->{ chr };
    my $len = @{ $pos };
    for(my $i = 0; $i < @{ $CDS }; $i++) {
	my $EXN = ($ori eq "+" || $ori eq "1") ? $i + 1 : @{ $CDS } - $i; # Exon number
	my ($START, $END) = @{$CDS->[$i]};
	my $tchr = $chr;
	$tchr =~ s/chr// if ( $opt_t );
	open(SAM, "samtools view $BAM $tchr:$START-$END |");
	$exoncov[$EXN] = 0;
	while( <SAM> ) {
	    my @a = split(/\t/);
	    my $start = $a[3];
	    my $n = 0;
	    my $dir = $a[1] & 0x10 ? "-" : "+";
	    if ( $opt_a ) {
		my @segs = $a[5] =~ /(\d+)[MI]/g; # Only match and insertion counts toward read length
		my $rlen = $stat->sum(\@segs); # The read length for matched bases
		my ($fstart, $fend) = ($start, $start + $rlen); # fragment start and end positions
		unless ( $a[1] & 0x8 && $a[6] eq "=" ) { # the mate is mapped
		    $fstart = $a[7] if ( $a[7] < $fstart );
		    $fend = $a[7] + $rlen if ( $fend < $a[7] + $rlen ); # assuming the mate has the same read length
		}
		my @t = sort { $a <=> $b } ($START-$PRIMERLEN, $END+$PRIMERLEN, $fstart, $fend);
		if ( ($t[2]-$t[1]+1)/($fend-$fstart+1) < $opt_a ) {
		    next;
		}
		next unless ($fstart <= $START-$PRIMERLEN || $END+$PRIMERLEN <= $fend );
		    #print STDERR "SKIPPED: $START $END @t $rlen $_";
		#print STDERR "$a[0]/$dir\n";
	    }

	    $exoncov[$EXN]++;
	}
	close( SAM );
    }
    my $total = 0;
    my @covs;
    my $gene_length = 0;
    my ($gene_start, $gene_end) = (500000000, 0);
    for(my $i = 0; $i < @{ $CDS }; $i++) {
	my $EXN = ($ori eq "+" || $ori eq "1")? $i + 1 : @{ $CDS } - $i; # Exon number
	my ($START, $END) = @{$CDS->[$i]};
	$gene_length += $END - $START + 1;
	$gene_start = $START if ( $START < $gene_start );
	$gene_end = $END if ( $END > $gene_end );
	$total += $exoncov[$EXN];
	#print STDERR "Coverage should be equal $etotal and $exoncov[$EXN] for $chr $START $END\n" unless( $etotal == $exoncov[$EXN] );
	#$EXN = "0$EXN" if ( $EXN < 10 );
	#print join("\t", $sample, $gene, $chr, $START, $END, "Seg-$EXN:$START-$END", $END-$START+1, map { sprintf("%.3f", $_/($END-$START+1)); } ($etotal, @ecovs)), "\n";
	#print join("\t", $sample, $gene, $chr, $START, $END, "Amplicon", $END-$START+1, map { sprintf("%.3f", $_/($END-$START+1)); } ($etotal, @ecovs)), "\n";
	print join("\t", $sample, $gene, $chr, $START, $END, "Amplicon", $END-$START+1, $exoncov[$EXN]), "\n";
    }
    #print STDERR "Length should be equal $gene_length and $len\n" unless ( $gene_length == $len );
    #print join("\t", $sample, $gene, $chr, $gene_start, $gene_end, "Whole-Gene", $gene_length, map { sprintf("%.3f", $_/$len); } ($total, @covs)), "\n";
    print join("\t", $sample, $gene, $chr, $gene_start, $gene_end, "Whole-Gene", $gene_length, $total), "\n";
}

sub USAGE {
    print STDERR <<USAGE;
    $0 [-hz] [-n name_reg] [-b bam] [-c chr] [-S start] [-E end] [-s seg_starts] [-e seg_ends] [-x #_nu] [-g gene] [-o ori] [-d depth] [-a overlap] region_info

    The program will calculate candidate variance for a given region(s) in an indexed BAM file.  The default
    input is IGV's one or more entries in refGene.txt, but can be regions.

    -h Print this help
    -t Indicate whether to trim chr from chromosome name.  Set if the BAM file chrs are 1, 2... 
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
       Desired depths, default to 1:10:20:50:100:250:500:1000
    -z 
       Indicate whether it's zero based numbering, default is 1-based
    -a num
       Indicate that it's in amplicon mode.  Reads/inserts overlaping with amplicon greater than [num] is considered
       belonging to amplicon.  recommend: 0.6 or 60%
    -p int
       The primer length.  Used with -a option and the BED file doesn't include primer.
USAGE
   exit(0);
}

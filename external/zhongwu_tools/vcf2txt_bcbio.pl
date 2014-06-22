#!/usr/bin/perl -w

use lib '/users/kdld047/lib/perl5';
use Getopt::Std;
use Stat::Basic;
use strict;

our ($opt_F, $opt_R, $opt_r, $opt_n, $opt_N, $opt_H, $opt_q, $opt_p, $opt_b, $opt_f, $opt_c, $opt_u, $opt_D, $opt_Q, $opt_P);
getopts('uHbR:F:f:n:r:p:q:c:D:P:Q:') || USAGE();
USAGE() if ( $opt_H );
my $FRACTION = $opt_r ? $opt_r : 0.4;
my $MAXRATIO = $opt_R ? $opt_R : 0.75;
my $CNT = $opt_n ? $opt_n : 10;
my $AVEFREQ = $opt_f ? $opt_f : 0.15;
my $MINPMEAN = $opt_p ? $opt_p : 5;
my $MINQMEAN = $opt_q ? $opt_q : 25;
my $FILPMEAN = $opt_P ? $opt_P : 0; # will be filtered on the first place
my $FILQMEAN = $opt_Q ? $opt_Q : 0; # will be filtered on the first place
my $FILDEPTH = $opt_D ? $opt_D : 0; # will be filtered on the first place
my $MINFREQ = $opt_F ? $opt_F : 0.05;
my $MAF = 0.0025;
my $control = $opt_c;


# SBF: Strand Bias Fisher Exact test
#my @columns = qw(CDS AA END DP AF BIAS PMEAN PSTD QUAL QSTD SBF GMAF dbSNPBuildID);
my @columns = qw(DP AF QUAL SS CDS AA CNT GMAF dbSNPBuildID);
#print join("\t", qw(Sample Chr Start ID Ref Alt Type Effect Functional_Class Codon_Change Amino_Acid_Change Amino_Acid_Length Gene Transcript_bioType Gene_Coding Transcript Exon COSMIC_CDS_Change COSMIC_AA_Change End Depth AlleleFreq Bias Pmean Pstd Qual Qstd SBF GMAF dbSNPBuildID N_samples N_Var Pcnt_sample Ave_AF PASS Var_Class)), "\n";
print join("\t", qw(Sample Chr Start ID Ref Alt Type Effect Functional_Class Codon_Change Amino_Acid_Change Amino_Acid_Length Gene Transcript_bioType Gene_Coding Transcript Exon Depth AlleleFreq Qual Status COSMIC_CDS_Change COSMIC_AA_Change COSMIC_CNT GMAF dbSNPBuildID N_samples N_Var Pcnt_sample Ave_AF PASS Var_Class)), "\n";
my @data;
my %sample;
my %var;
my $stat = new Stat::Basic;
my %CONTROL;
my $SAMPLE = "";
while( <> ) {
    chomp;
    my @a = split(/\t/);
    if ( /^#CHROM/ ) {
        $SAMPLE = $a[9];
	$SAMPLE = $a[10] if ( $a[9] eq "none" );
	next;
    }
    next if ( /^#/ );
    $a[7] .= ";";
    my %d;
    while( $a[7] =~ /([^=]+)=([^=]+);/g ) {
	my ($k, $v) = ($1, $2);
	$k = "EFF" if ( $k =~ /;EFF$/ );
	$d{ $k } = $v;
    }
    my @fmt = split(/:/, $a[8]);
    my @fd = split(/:/, ($a[9] =~ /0.00:0$/ ? $a[10] : $a[9]));
    for(my $i = 0; $i < @fmt; $i++) {
        $d{ $fmt[$i] } = $fd[$i];
    }
    $d{ QUAL } = $d{ BQ };
    $d{ AF } = $d{ FA };
    my @effs = split(/,/, $d{ EFF });
    my $vark = join(":", @a[0,1,3,4]); # Chr Pos Ref Alt
    next if ( $FILDEPTH && $d{ DP } < $FILDEPTH );
    #next if ( $FILPMEAN && $d{ PMEAN } < $FILPMEAN );
    next if ( $FILQMEAN && $d{ QUAL } < $FILQMEAN );
    if ( $control && $control eq $SAMPLE ) {
	my ($pmean, $qmean) = ($d{ PMEAN }, $d{ QUAL });
	my $pass = "TRUE";
	#$pass = "FALSE" unless ( $d{PSTD} > 0 );
	$pass = "FALSE" if ($qmean < $MINQMEAN );
	#$pass = "FALSE" if ($pmean < $MINPMEAN );
	$pass = "FALSE" if ( $d{AF} < $MINFREQ );
	my $class = $a[2] =~ /COSM/ ? "COSMIC" : ($a[2] =~ /^rs/ ? "dbSNP" : "Novel");
        $CONTROL{ $vark } = 1 if ( $pass eq "TRUE" && $class ne "dbSNP" );
        #$CONTROL{ $vark } = 1 if ( $class eq "Novel" );
    }
    unless( $opt_u && $SAMPLE =~ /Undetermined/i ) { # Undetermined won't count toward samples
	$sample{ $SAMPLE } = 1;
	push( @{ $var{ $vark } }, $d{ AF } );
    }
    foreach my $eff (@effs) {
        $eff =~ s/\)$//;
	my @e = split(/\|/, $eff, -1);
	my ($type, $effect) = split(/\(/, $e[0]);
	my @tmp = map { defined($d{ $_ }) ? $d{ $_ } : ""; } @columns;
	my @tmp2= map { defined($_) ? $_ : ""; } @e[1..9];
	push(@data, [$SAMPLE, @a[0..4], $type, $effect, @tmp2, @tmp]);
    }
}

my @samples = keys %sample;
my $sam_n = @samples + 0;
foreach my $d (@data) {
    my $vark = join(":", @$d[1, 2, 4, 5]); # Chr Pos Ref Alt
    next unless( $var{ $vark } ); # Likely just in Undetermined.
    #my ($pmean, $qmean) = @$d[23,25];
    my ($pmean, $qmean) = @$d[17,17];
    my $varn = @{ $var{ $vark } } + 0;
    my $ave_af = $stat->mean( $var{ $vark } );
    my $pass = ($varn/$sam_n > $FRACTION && $varn >= $CNT && $ave_af < $AVEFREQ && $d->[3] eq "." && $d->[18] < 0.25 ) ? "FALSE" : "TRUE"; # novel and present in $MAXRATIO samples
    #$pass = "FALSE" unless ( $d->[24] > 0 ); # all variants from one position in reads
    #$pass = "TRUE" if ( $d->[24] ==  0 && $d->[22] !~ /1$/ && $d->[22] !~ /0$/ ); # all variants from one position in reads
    $pass = "FALSE" if ($d->[3] eq "." && $varn/$sam_n > $MAXRATIO ); # novel and present in $MAXRATIO samples, regardless of frequency
    $pass = "FALSE" if ($qmean < $MINQMEAN );
    #$pass = "FALSE" if ($pmean < $MINPMEAN );
    $pass = "FALSE" if ( $d->[18] < $MINFREQ );
    my $class = $d->[3] =~ /COSM/ ? "COSMIC" : ($d->[3] =~ /^rs/ ? "dbSNP" : "Novel");
    $class = "dbSNP" if ( $d->[24] && $d->[24] > $MAF ); # if there's MAF with frequency, it'll be considered dbSNP regardless of COSMIC
    $pass = "FALSE" if ( $control && $CONTROL{ $vark } );
    #$pass = "FALSE" if ( $opt_b && $class eq "Novel" && ($d->[22] eq "2;1" || $d->[22] eq "2;0") ); # Filter novel variants with strand bias.
    print join("\t", @$d, $sam_n, $varn, sprintf("%.3f", $varn/$sam_n), $ave_af, $pass, $class), "\n";
}

sub USAGE {
print <<USAGE;
    The program will convert an annotated vcf files by snfEFF using dbSNP and COSMIC back to txt format.  It also checks for quality
    and add "PASS" column.  It will not perform any filtering.
    
    Usage: $0 [-H] [-F var_fraction] [-n sample_cnt] [-f freq] [-p pos] [-q quality] vcf_files
    
    The program accepts more than one vcf files.

    Options:
    -H Print this help page
    -b Novel variants with strand bias "2;1" will be considered as false positive
    -u Undeteremined won't be counted for the sample count.
    -r DOUBLE
        When a novel variant is present in more than [fraction] of samples and mean allele frequency is less than -f, it's 
	considered as likely false positive. Default 0.4.
    
    -R DOUBLE
        When a novel variant is present in more than [fraction] of samples, regardless allele frequency, it's considered as 
	likely false positive. Default 0.75.

    -f DOUBLE
        When the ave allele frequency is also below the [freq], the variant is considered likely false positive.  Default 0.15.

    -F DOUBLE
        When indivisual allele frequency < feq for variants, it was considered likely false poitives. Default: 0.05 or 5%

    -n INT
        When the variant is detected in greater or equal [sample_cnt] samples, the variant is considered likely false positive.  Default 10.

    -p INT
        The minimum mean position in reads for variants.  Default: 5bp

    -q DOUBLE
        The minimum mean base quality phred score for variants.  Default: 25

    -P INT
        The filtering mean position in reads for variants.  The raw variant will be filtered on first place if the mean 
	posisiton is less then INT.  Default: 0bp

    -Q DOUBLE
        The filtering mean base quality phred score for variants.  The raw variant will be filtered on first place 
	if the mean quality is less then DOUBLE.  Default: 0

    -D INT
        The filtering depth for variants.  The raw variant will be filtered on first place if the depth is less then INT.  Default: 0

    -c Control
        The control sample name.  Any novel or COSMIC variants passing all above filters but also detected in Control sample will be deemed considered
	false positive.  Use only when there's control sample.

    A novel variant (non-dbSNP, non-COSMIC) is considered false positive if all three conditions (-r -f -n) are met. Any variant meeting the -p
    or -q conditions are also considered likely false positive.  False positive variants are annotated "FALSE" in column PASS, "TRUE" otherwise.
USAGE
    exit(0);
}

#!/usr/bin/perl -w

use strict;

print "##fileformat=VCFv4.1\n#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO\n";
#chr9    93606272        .       AG      A       .       PASS    DP=47549;END=93606273;AF=0.014;BIAS=2;2;SM=AZ47;PMEAN=33.0;PSTD=9.55
#chr9    93606299        .       T       G       .       PASS    DP=51402;END=93606299;AF=0.013;BIAS=2;1;SM=AZ47;PMEAN=1.2;PSTD=1.58
#chr9    93606397        .       A       C       .       PASS    DP=55670;END=93606397;AF=0.011;BIAS=2;1;SM=AZ47;PMEAN=14.3;PSTD=3.52

my %hash;
while(<>) {
    chomp;
    my @a = split(/\t/);
    push( @{ $hash{ $a[2] }->{ $a[3] } }, $_ );
}
my @chrs = map { "chr$_"; } (1..22);
push(@chrs, "chrX", "chrY", "chrM");
foreach my $chr (@chrs) {
    my @pos = sort { $a <=> $b } (keys %{ $hash{ $chr } });
    foreach my $p (@pos) {
	foreach my $d (@{ $hash{ $chr }->{ $p } }) {
	    my @a = split(/\t/, $d);
	    my $oddratio = $a[21];
	    if ( $oddratio eq "Inf" ) {
	        $oddratio = 0;
	    } elsif ( $oddratio < 1 && $oddratio > 0 ) {
	        $oddratio = 1/$oddratio;
	    }
	    print  join("\t", $a[2], $a[3], ".", @a[5,6], ".", "PASS", "SAMPLE=$a[0];DP=$a[7];END=$a[4];VP=$a[8];AF=$a[14];BIAS=$a[15];REFBIAS=$a[9]:$a[10];VARBIAS=$a[11]:$a[12];PMEAN=$a[16];PSTD=$a[17];QUAL=$a[18];QSTD=$a[19];SBF=$a[20];ODDRATIO=$oddratio;MQ=$a[22];QRATIO=$a[23];LSEQ=$a[24];RSEQ=$a[25]"), "\n";
	}
    }
}
#AZ01	EZH2	chr7	148504716	148504717	AG	A	20852	17250	3495	0	17249	1	-1/G	0.827	1;1	41.5	2.44	CAGAGG	GGGGGA	A:28:F-28:R-0	T:12:F-12:R-0	C:17:F-17:R-0	-2:50:F-50:R-0	G:3495:F-3495:R-0	-1:17250:F-17249:R-1
#AZ01	EZH2	chr7	148506396	148506396	A	C	17774	15940	1801	1	15940	0	C/A	0.897	1;1	34.1	2.31	AAAGGT	CCTACC	A:1802:F-1801:R-1	T:17:F-17:R-0	C:15940:F-15940:R-0	G:9:F-9:R-0	-1:6:F-6:R-0
##fileformat=VCFv4.1
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
#chr9    93606272        .       AG      A       .       PASS    DP=47549;END=93606273;AF=0.014;BIAS=2;2;SM=AZ47;PMEAN=33.0;PSTD=9.55
#chr9    93606299        .       T       G       .       PASS    DP=51402;END=93606299;AF=0.013;BIAS=2;1;SM=AZ47;PMEAN=1.2;PSTD=1.58
#chr9    93606397        .       A       C       .       PASS    DP=55670;END=93606397;AF=0.011;BIAS=2;1;SM=AZ47;PMEAN=14.3;PSTD=3.52


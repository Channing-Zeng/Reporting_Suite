#!/usr/bin/env perl

use strict;
use warnings;
use POSIX;

my $usage = "\n\nusage: $0 file.vcf [interval_length=5000] [count_homozygous_variants=0]\n\n";

my $vcf_file = $ARGV[0] or die $usage;
my $interval_length = $ARGV[1] || 5000;
my $count_homozygous_variants = $ARGV[2] || 0;

my $MIN_COVERAGE = 1;

main: {


    my %interval_to_count;
    
    my $counter = 0;

    print STDERR "-processing VCF file: $vcf_file\n";
	open (my $fh, $vcf_file) or die "Error, cannot open file $vcf_file";
	while (<$fh>) {
		chomp;
		if (/^\#/) { 
			next;
		}
		
        $counter++;
        if ($counter % 1000 == 0) {
            print STDERR "\r[$counter]  ";
        }
        
		my @x = split(/\t/);

		my $contig = $x[0];
        my $pos = $x[1];
        
		my $base = $x[3];
		unless ($base =~ /[GATC]/i) { next; }
		

        
        my $alt = $x[4];
        unless (length($base) == length($alt) || $alt =~ /,/) { next; } ## skip the indel lines 
        

        my $snp_info_line = $x[7];
        
        my $depth_of_cov = 0;

        if ($snp_info_line =~ /DP=(\d+)/) {
            $depth_of_cov = $1;
        }

        
        my $interval_id = $contig . "-I" . sprintf("%05i", ceil($pos / $interval_length));
        
        
        if ($depth_of_cov < $MIN_COVERAGE) {
            next;
        }
        
        # count base as having sufficient read support.
        $interval_to_count{$contig}->{$interval_id}->{cov}++;
        
        # increment depth of coverage for feature.
        $interval_to_count{$contig}->{$interval_id}->{depth} += $depth_of_cov;
        
        if ($x[6] =~ /GATKStandard|LowQual/) {
            ## ignore low quality snps, but still count the coverage info.
            next; 
        }
        
        if ($alt ne '.') {
            
            my $is_homozygous_variant = 0;

            if (/AF=([^;]+)/) {
                my $allele_freq = $1;
                if ($allele_freq == 1) {
                    $is_homozygous_variant = 1;
                }
            }
            
            unless ($is_homozygous_variant && ! $count_homozygous_variants) {
            
                # track the counts of SNPs for feature
                $interval_to_count{$contig}->{$interval_id}->{mut}++;
            }
        }
        
    }
    
    foreach my $contig (sort keys %interval_to_count) {
        print STDERR "// processing genes on $contig\n";
        
        foreach my $interval_id (sort keys %{$interval_to_count{$contig}} ) {
            
            my $covered_bases = $interval_to_count{$contig}->{$interval_id}->{cov} || 0;
            my $SNP_bases = $interval_to_count{$contig}->{$interval_id}->{mut} || 0;
            my $depth_of_coverage = $interval_to_count{$contig}->{$interval_id}->{depth};
            
            print "$interval_id\t$covered_bases\t$SNP_bases" 
                . sprintf("\t%.2f", $SNP_bases/$covered_bases*100) 
                . sprintf("\t%.1f", $depth_of_coverage/$covered_bases)
                . "\n";
        }
        
    }
    
    
	exit(0);
}


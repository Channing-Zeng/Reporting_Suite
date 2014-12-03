#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "\n\nusage: $0 vcf.annotated\n\n";

my $vcf_file = $ARGV[0] or die $usage;

my $MIN_COVERAGE = 1;

main: {


    my %trans_to_SNP_type_count;
    
    my $counter = 0;

    my %base_feature_counts;

    print STDERR "-processing VCF file: $vcf_file\n";
	open (my $fh, $vcf_file) or die "Error, cannot open file $vcf_file";
	while (<$fh>) {
		chomp;
		if (/^\#/) { 
			next;
		}
		unless (/\w/) { next; }

        $counter++;
        if ($counter % 1000 == 0) {
            print STDERR "\r[$counter]  ";
        }
        
		my @x = split(/\t/);

        unless (scalar(@x) >= 10) { next; }
        
		my $contig = $x[0];
        
		my $base = $x[3];
		
        unless ( defined ($base)) {
            die "Error, line doesn't parse: $_";
        }

        unless ($base =~ /[GATC]/i) { next; }
		
        my $alt = $x[4];


		my @annots = @x[10..$#x];
		
        
        		
        my $snp_info_line = $x[7];
        
        my $depth_of_cov = 0;

        if ($snp_info_line =~ /DP=(\d+)/) {
            $depth_of_cov = $1;
        }

                        
		foreach my $annot (@annots) {
			unless ($annot =~ /\w/) { 
				next;
			}
        			
			my ($type, @rest) = split(/\s+|,/, $annot);
            my ($gene_id, $model_id, $orient, $coding_pos, $frame) = @rest;
            
            $base_feature_counts{$type}->{total}++;


            if ($type eq "intergenic") {
                my $prev_gene_coords = "start";
                if ($annot =~ /prev_gene\(\d+-(\d+)/) {
                    $prev_gene_coords = $1+1;
                }
                my $next_gene_coords = "end";
                if ($annot =~ /next_gene\((\d+)-\d+/) {
                    $next_gene_coords = $1 - 1;
                }
                $model_id = "intergenic:$contig:$prev_gene_coords-$next_gene_coords";

            }
            
            if ($depth_of_cov < $MIN_COVERAGE) {
                next;
            }

            # count base as having sufficient read support.
            $trans_to_SNP_type_count{$contig}->{$model_id}->{$type}->{cov}++;
            $base_feature_counts{$type}->{cov}++;
            

            # increment depth of coverage for feature.
            $trans_to_SNP_type_count{$contig}->{$model_id}->{$type}->{depth} += $depth_of_cov;
            


            if ($x[6] =~ /GATKStandard|LowQual/) {
                ## ignore low quality snps, but still count the coverage info.
                next; 
            }
            
            if ($alt ne '.') {
                # track the counts of SNPs for feature
                $trans_to_SNP_type_count{$contig}->{$model_id}->{$type}->{mut}++;
            }
        }
    }
    
    
    my @types = qw(CDS intron p3UTR p5UTR gene intergenic);
    
    my %type_to_fh;
    foreach my $type (@types) {
        open (my $ofh, ">$vcf_file.$type") or die $!;
        $type_to_fh{$type} = $ofh;
    }
    
    foreach my $contig (keys %trans_to_SNP_type_count) {
        print STDERR "// processing genes on $contig\n";
        
        
        foreach my $model_id (keys %{$trans_to_SNP_type_count{$contig}} ) {
            
            if ($model_id =~ /^intergenic/) {
                
                my $ofh = $type_to_fh{intergenic};
                my $covered_bases = $trans_to_SNP_type_count{$contig}->{$model_id}->{intergenic}->{cov} || 0;
                my $SNP_bases = $trans_to_SNP_type_count{$contig}->{$model_id}->{intergenic}->{mut} || 0;
                my $depth_of_coverage = $trans_to_SNP_type_count{$contig}->{$model_id}->{intergenic}->{depth};
                
                print $ofh "$model_id\t$covered_bases\t$SNP_bases" 
                    . sprintf("\t%.2f", $SNP_bases/$covered_bases*100) 
                    . sprintf("\t%.1f", $depth_of_coverage/$covered_bases)
                    . "\n";
                
                next;
            }
            
            my $type_counts_href = $trans_to_SNP_type_count{$contig}->{$model_id};

            
            my $cds_length = $type_counts_href->{CDS}->{cov} || 0;
            my $cds_snp_count = $type_counts_href->{CDS}->{mut} || 0;
            my $cds_depth_of_cov = $type_counts_href->{CDS}->{depth} || 0;
            if ($cds_length) {
                my $ofh = $type_to_fh{CDS};
                print $ofh "$model_id\t$cds_length\t$cds_snp_count" 
                    . sprintf("\t%.2f", $cds_snp_count/$cds_length*100) 
                    . sprintf("\t%.1f", $cds_depth_of_cov/$cds_length) 
                    . "\n";
            }
            
            my $intron_length = $type_counts_href->{intron}->{cov} || 0;
            my $intron_snp_count = $type_counts_href->{intron}->{mut} || 0;
            my $intron_depth_of_cov = $type_counts_href->{intron}->{depth} || 0;
            
            if ($intron_length) {
                my $ofh = $type_to_fh{intron};
                print $ofh "$model_id\t$intron_length\t$intron_snp_count" 
                    . sprintf("\t%.2f", $intron_snp_count/$intron_length*100) 
                    . sprintf("\t%.1f", $intron_depth_of_cov/$intron_length) 
                    . "\n";
            }
            
            my $p5UTR_length = $type_counts_href->{'p5UTR'}->{cov} || 0;
            my $p5UTR_snp_count = $type_counts_href->{'p5UTR'}->{mut} || 0;
            my $p5UTR_depth_of_cov = $type_counts_href->{'p5UTR'}->{depth} || 0;

            if ($p5UTR_length) {
                my $ofh = $type_to_fh{'p5UTR'};
                print $ofh "$model_id\t$p5UTR_length\t$p5UTR_snp_count" 
                    . sprintf("\t%.2f", $p5UTR_snp_count/$p5UTR_length*100) 
                    . sprintf("\t%.1f", $p5UTR_depth_of_cov/$p5UTR_length) 
                    . "\n";
            }
            
            
            my $p3UTR_length = $type_counts_href->{'p3UTR'}->{cov} || 0;
            my $p3UTR_snp_count = $type_counts_href->{'p3UTR'}->{mut} || 0;
            my $p3UTR_depth_of_cov = $type_counts_href->{'p3UTR'}->{depth} || 0;
            
            if ($p3UTR_length) {
                my $ofh = $type_to_fh{'p3UTR'};
                print $ofh "$model_id\t$p3UTR_length\t$p3UTR_snp_count" 
                    . sprintf("\t%.2f", $p3UTR_snp_count/$p3UTR_length*100) 
                    . sprintf("\t%.1f", $p3UTR_depth_of_cov/$p3UTR_length) 
                    . "\n";
            }
            
            my $gene_length = $cds_length + $intron_length + $p5UTR_length + $p3UTR_length;
            my $gene_snps = $cds_snp_count + $intron_snp_count + $p5UTR_snp_count + $p3UTR_snp_count;
            my $gene_depth_of_cov = $cds_depth_of_cov + $intron_depth_of_cov + $p5UTR_depth_of_cov + $p3UTR_depth_of_cov;
            
            if ($gene_length) {
                my $ofh = $type_to_fh{'gene'};
                print $ofh "$model_id\t$gene_length\t$gene_snps" 
                    . sprintf("\t%.2f", $gene_snps/$gene_length*100) 
                    . sprintf("\t%.1f", $gene_depth_of_cov/$gene_length) 
                    . "\n";
            }
            
            
        }
        
    }
    
    foreach my $fh (values %type_to_fh) {
        close $fh;
    }

    ## write summary of coverage by feature type.
    {
        open (my $ofh, ">$vcf_file.summary") or die $!;
        foreach my $type (sort keys %base_feature_counts) {
            my $count = $base_feature_counts{$type}->{cov};
            my $total = $base_feature_counts{$type}->{total};

            print $ofh "$type\t$count\t$total\t" . sprintf("%.2f", $count/$total*100) . "\n";
        }

        close $ofh;
    }
    



    print STDERR "Done.\n";
    
	exit(0);
}


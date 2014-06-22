#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 feature_file_prefix max_cov_val\n\n";

my $prefix = $ARGV[0] or die $usage;
my $max_cov_val = $ARGV[1] or die $usage;

my @features = qw(intergenic intron gene p5UTR CDS p3UTR);

main: {

    print "max_y = 0\n";

    foreach my $feature_type (@features) {
        
        my $file = "$prefix.$feature_type";
        
        print "$feature_type = read.table(\"$file\")\n";
        print "$feature_type.cov = $feature_type\[,5]\n";
        print "$feature_type.cov_adj = $feature_type.cov[$feature_type.cov <= $max_cov_val]\n";
        print "$feature_type.hist = hist($feature_type.cov_adj, br=100, plot=F)\n";
        print "max_y = max(max_y, $feature_type.hist\$density)\n";


        print "\n";
    }

    print "line_colors = rainbow(" . scalar(@features) . ")\n";
    
    for (my $i = 0; $i <= $#features; $i++) {
        
        my $feature_type = $features[$i];
        
        if ($i == 0) {
            print "plot($feature_type.hist\$mids, $feature_type.hist\$density, col=line_colors[" . ($i+1) . "], t='l', ylim=c(0,max_y))\n";
        }
        else {
            print "points($feature_type.hist\$mids, $feature_type.hist\$density, col=line_colors[" . ($i+1) . "], t='l')\n";
        }

        
    }

    print "legend('topright', legend=c('" . join("','", @features) . "\'), pch=15, col=line_colors)\n";
    


    exit(0);
}

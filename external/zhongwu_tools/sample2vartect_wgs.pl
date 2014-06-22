#!/usr/bin/perl -w

# From a sample list to shell script for running variant detection using checkVar.pl

use Getopt::Std;
use strict;

our ($opt_p, $opt_b, $opt_H, $opt_a, $opt_A, $opt_o);
getopts("Hp:b:a:A:o:") || USAGE();
$opt_H && USAGE();

my $option = $opt_o ? $opt_o : "";

my $prog = $opt_p ? $opt_p : "/group/cancer_informatics/tools_resources/NGS/bin/vartect_hg19_wgs.sh";
my $bed = $opt_b ? $opt_b : "/group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/human_hg19_5k_segs_";
while( <> ) {
    chomp;
    my ($sample, $bam) = split(/\t/);
    $bam =~ s/\@$//;
    $sample .= $opt_a if ( $opt_a );
    $sample = $opt_A . $sample if ( $opt_A );
    print "mkdir $sample\n";
    print "cd $sample\n";
    for(my $i = 1; $i <= 40; $i++) {
	print "qsub $option -cwd -V -S /bin/bash -N VT_${i}_$sample $prog $bam $sample $bed$i.txt $i\n";
    }
    print "cd ..\n\n";
}

sub USAGE {
    print <<USAGE;
    Usage: $0 -H -p alignment_program -b bed_file sample2bam

    The program will generate shell scripts for submitting vartect jobs to clusters.  When submitting big jobs, pipe the output
    to addQueue.pl, which will submit jobs only to chara, rask, espo, and orr

    -H Print this help page.
    -p alignment_program
       The ABSOLUTE path to the alignment script.  Default: /group/cancer_informatics/tools_resources/NGS/bin/vartect_hg19.sh,
       which will perform alignment using BWA 0.7.4 and call variants using checkVar.pl
    -b bed_file
       The ABSOLUTE path to the bed file.  If not supplied, it'll make variant calls to whole gene using 1Mb segments.
       Default: /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/human_hg19_1M_segs.txt
    -o options
       Options to be appended to qsub, such as "-l huge_ram=1", "-l mem_free=4G", etc.
    -a string
       A string to append to the sample name, e.g. -mm10 when running against mouse genome
    -A string
       A string to suffix to the sample name, e.g. dis- when running after mouse disambiguation

USAGE
    exit(0);
}

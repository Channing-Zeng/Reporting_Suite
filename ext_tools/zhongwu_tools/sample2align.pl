#!/usr/bin/perl -w

# From a sample list to shell script for running alignments

use Getopt::Std;
use strict;

our ($opt_p, $opt_b, $opt_H, $opt_a, $opt_A, $opt_f, $opt_q);
getopts("Hp:b:a:A:f:q:") || USAGE();
$opt_H && USAGE();

my $prog = $opt_p ? $opt_p : "/group/cancer_informatics/tools_resources/NGS/bin/runBWA074_hg19.sh";
my $bed = defined($opt_b) ? $opt_b : "/group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/human_hg19_1M_segs.txt";
my $freq = defined($opt_f) ? $opt_f : 0.01;
my $queue = $opt_q ? $opt_q : "batch.q";
while( <> ) {
    chomp;
    my ($sample, $read1, $read2) = split(/\t/);
    $read1 =~ s/\@//;
    $read2 = "" unless( $read2 );
    $read2 =~ s/\@//;
    $sample .= $opt_a if ( $opt_a );
    $sample = $opt_A . $sample if ( $opt_A );
    print "if [ ! -d \"$sample\" ]; then\n";
    print "    mkdir $sample\n";
    print "fi\n";
    print "cd $sample\n";
    #print "qsub -l mem_free=24G -cwd -V -pe smp 8 -N runBWA_$sample -S /bin/bash $prog $read1 $read2 $sample $bed\n";
    #print "qsub -l huge_ram=1 -l mem_free=24G -cwd -V -pe smp 8 -N ${sample}_BWA -S /bin/bash $prog $read1 $read2 $sample $bed $freq\n";
    print "qsub -l mem_free=16G -q $queue -cwd -V -pe smp 8 -N ${sample}_BWA -S /bin/bash $prog $read1 $read2 $sample $bed $freq\n";
    print "cd ..\n\n";
}

sub USAGE {
    print <<USAGE;
    Usage: $0 -H -p alignment_program -b bed_file sample2fastq

    The program will generate shell scripts for submitting aligment jobs to clusters.  When submitting big jobs, pipe the output
    to addQueue.pl, which will submit jobs only to chara, rask, espo, and orr

    -H Print this help page.
    -p alignment_program
       The ABSOLUTE path to the alignment script.  Default: /group/cancer_informatics/tools_resources/NGS/bin/runBWA074_hg19.sh,
       which will perform alignment using BWA 0.7.4 and call variants using checkVar.pl
    -b bed_file
       The ABSOLUTE path to the bed file.  If not supplied, it'll make variant calls to whole gene using 1Mb segments.
       Default: /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Annotation/human_hg19_1M_segs.txt
    -a string
       A string to append to the sample name, e.g. -mm10 when running against mouse genome
    -A string
       A string to suffix to the sample name, e.g. mm10- when running against mouse genome
    -f allele frequency
       The lowest allele frequency.  Default to 1%, or 0.01
    -q queue
       The queue to use.  Default to batch.q

USAGE
    exit(0);
}

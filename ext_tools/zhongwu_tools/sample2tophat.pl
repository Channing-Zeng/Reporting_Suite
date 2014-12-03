#!/usr/bin/perl -w

# From a sample list to shell script for running alignments

use Getopt::Std;
use strict;

our ($opt_o, $opt_p, $opt_H, $opt_a, $opt_A, $opt_D);
getopts("Hp:a:A:D:o:") || USAGE();
$opt_H && USAGE();

my $prog = $opt_p ? $opt_p : "/group/cancer_informatics/tools_resources/NGS/bin/runTophat2_hg19.sh";
while( <> ) {
    chomp;
    my ($sample, $read1, $read2) = split(/\t/);
    $read1 =~ s/\@//;
    $read2 =~ s/\@//;
    $sample .= $opt_a if ( $opt_a );
    $sample = $opt_A . $sample if ( $opt_A );
    my $dir = $sample;
    $dir = $1 if ( $opt_D and $sample =~ /$opt_D/ );
    print "if [ ! -d \"$dir\" ]; then\n";
    print "    mkdir $dir\nfi\n";
    print "cd $dir\n";
    #print "qsub -l mem_free=24G -cwd -V -pe smp 8 -N runBWA_$sample -S /bin/bash $prog $read1 $read2 $sample $bed\n";
    my $options = $opt_o ? $opt_o : "-q batch.q";
    print "qsub $options -cwd -V -pe smp 8 -N ${sample}_TH -S /bin/bash $prog $sample $read1 $read2\n";
    print "cd ..\n\n";
}

sub USAGE {
    print <<USAGE;
    Usage: $0 -H -p alignment_program -b bed_file sample2fastq

    The program will generate shell scripts for submitting aligment jobs to clusters.  When submitting big jobs, pipe the output
    to addQueue.pl, which will submit jobs only to chara, rask, espo, and orr

    -H Print this help page.
    -o options
       Options for qsub, default to "-l huge_ram=1 -l mem_free=24G"
    -p alignment_program
       The ABSOLUTE path to the alignment script.  Default: /group/cancer_informatics/tools_resources/NGS/bin/runBWA074_hg19.sh,
       which will perform alignment using BWA 0.7.4 and call variants using checkVar.pl
    -a string
       A string to append to the sample name, e.g. -mm10 when running against mouse genome
    -A string
       A string to suffix to the sample name, e.g. mm10- when running against mouse genome
    -D regx
       A regular expression to capture the short sample name.  Used if a sample has more than one pairs of fastq files

USAGE
    exit(0);
}

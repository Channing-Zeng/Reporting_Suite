#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

import sys

from source.logger import err, critical
from source.file_utils import file_exists, verify_file, file_transaction


def main():
    if len(sys.argv) <= 2:
        critical('Usage: ' + __file__ + ' path_to_.fa genome_build_name')

    seq_fpath = sys.argv[1]
    genome_build = sys.argv[2]
    get_chr_lengths(seq_fpath, genome_build)


def get_chr_lengths(seq_fpath, genome_build, silence=False):
    chr_lengths = []

    verify_file(seq_fpath, is_critical=True)

    if verify_file(seq_fpath + '.fai'):
        err('Reading genome index file (.fai) to get chromosome lengths')
        with open(seq_fpath + '.fai', 'r') as handle:
            for line in handle:
                line = line.strip()
                if line:
                    chrom, length = line.split()[0], line.split()[1]
                    chr_lengths.append((chrom, length))
    else:
        err('Reading genome sequence (.fa) to get chromosome lengths')
        with open(seq_fpath, 'r') as handle:
            from Bio import SeqIO
            reference_records = SeqIO.parse(handle, 'fasta')
            for record in reference_records:
                chrom = record.id
                chr_lengths.append((chrom, len(record.seq)))

    #chr_lengths = sorted(chr_lengths, key=lambda (k, l): k.get_key())
    if silence:
        return chr_lengths
    for c, l in chr_lengths:
        sys.stdout.write(c + '\t' + str(l) + '\n')


if __name__ == '__main__':
    main()

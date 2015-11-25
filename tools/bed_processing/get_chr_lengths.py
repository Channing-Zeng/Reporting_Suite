import sys

from source.logger import info, critical
from source.file_utils import file_exists, verify_file, file_transaction
from source.targetcov.Region import SortableByChrom


def main():
    if len(sys.argv) > 1 :
        seq_fpath = sys.argv[1]
    else:
        seq_fpath = '/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa'

    get_chr_lengths(seq_fpath)


def get_chr_lengths(seq_fpath):
    chr_lengths = []

    verify_file(seq_fpath, is_critical=True)

    if verify_file(seq_fpath + '.fai'):
        info('Reading genome index file (.fai) to get chromosome lengths')
        with open(seq_fpath + '.fai', 'r') as handle:
            for line in handle:
                line = line.strip()
                if line:
                    chrom, length = line.split()[0], line.split()[1]
                    chr_lengths.append([SortableByChrom(chrom), length])
    else:
        info('Reading genome sequence (.fa) to get chromosome lengths')
        with open(seq_fpath, 'r') as handle:
            from Bio import SeqIO
            reference_records = SeqIO.parse(handle, 'fasta')
            for record in reference_records:
                chrom = record.id
                chr_lengths.append([SortableByChrom(chrom), len(record.seq)])

    chr_lengths = sorted(chr_lengths, key=lambda (k, l): (k.get_key(), l))

    for c, l in chr_lengths:
        sys.stdout.write(c.chrom + '\t' + str(l) + '\n')


if __name__ == '__main__':
    main()
